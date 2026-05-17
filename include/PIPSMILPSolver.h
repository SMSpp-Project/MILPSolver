/*--------------------------------------------------------------------------*/
/*--------------------------- File PIPSMILPSolver.h ------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the PIPSMILPSolver class.
 *
 * PIPSMILPSolver is an SMS++ MILPSolver interface for PIPS-IPM++.
 * The solver builds a two-level PIPS DistributedInputTree from an SMS++ Block
 * hierarchy and provides all matrix/vector data to PIPS through callbacks.
 */
/*--------------------------------------------------------------------------*/

#ifndef __PIPSMILPSOLVER_H
 #define __PIPSMILPSOLVER_H

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include <algorithm>
#include <stdexcept>
#include <unordered_map>
#include <utility>
#include <vector>

#include "mpi.h"

#include "DistributedInputTree.h"
#include "PIPSIPMppInterface.hpp"
#include "MILPSolver.h"
#include "PIPSIPMppOptions.h"

/*--------------------------------------------------------------------------*/
/*-------------------------- PIPS CALLBACK TYPES ---------------------------*/
/*--------------------------------------------------------------------------*/
/**
 * PIPS-IPM++ asks the user code for dimensions, matrices and vectors through
 * plain C-style callbacks. The first argument, user_data, is an opaque pointer
 * that this solver passes as `this`; every callback casts it back to
 * PIPSMILPSolver* and accesses the already-built node data.
 */
extern "C" {

/// Callback returning a dimension or a number of nonzeros.
typedef int (*FNNZ)( void * user_data , int id , int * nnz );

/// Callback returning a sparse matrix in CSR row-major format.
typedef int (*FMAT)( void * user_data , int id , int * krowM ,
                    int * jcolM , double * M );

/// Callback returning a dense vector.
typedef int (*FVEC)( void * user_data , int id , double * vec , int len );

} // extern "C"

/*--------------------------------------------------------------------------*/
/*----------------------------- NAMESPACE ----------------------------------*/
/*--------------------------------------------------------------------------*/

namespace SMSpp_di_unipi_it {

class LinearFunction;
class FRowConstraint;
class ColVariable;

/*--------------------------------------------------------------------------*/
/*----------------------- CLASS PIPSMILPSolver -----------------------------*/
/*--------------------------------------------------------------------------*/

class PIPSMILPSolver : public MILPSolver {

 public:

/*--------------------------------------------------------------------------*/
/*---------------------------- PUBLIC TYPES --------------------------------*/
/*--------------------------------------------------------------------------*/

 using Subset = Block::Subset;
 using c_Subset = Block::c_Subset;
 using Range = Block::Range;

 using Index = Block::Index;  // "import" Index from Block

/*--------------------------------------------------------------------------*/
/*--------------------- CONSTRUCTOR AND DESTRUCTOR -------------------------*/
/*--------------------------------------------------------------------------*/

 PIPSMILPSolver( void );
 ~PIPSMILPSolver() override;

/*--------------------------------------------------------------------------*/
/*--------------------- DERIVED METHODS OF BASE CLASS ----------------------*/
/*--------------------------------------------------------------------------*/

 /// Sets the SMS++ Block to be solved.
 void set_Block( Block * block ) override;

 /// Loads the SMS++ Block and builds the PIPS DistributedInputTree.
 void load_problem( void ) override;

 /// Solves the currently loaded problem with PIPS-IPM++.
 int compute( bool changedvars = false ) override;

 /// returns a valid lower bound on the optimal objective function value
 OFValue get_lb( void ) override;

 /// returns a valid upper bound on the optimal objective function value
 OFValue get_ub( void ) override;

 /// Clears the current PIPS tree/interface and the base MILPSolver data.
 void clear_problem( unsigned int what ) override;

 protected:

/*--------------------------------------------------------------------------*/
/*------------------------- PIPS MATRIX HELPERS ----------------------------*/
/*--------------------------------------------------------------------------*/

 /// Sparse matrix in CSR row format, as required by PIPS matrix callbacks.
 struct CSRMatrix {
  std::vector< int > krow;    ///< row pointers, size = rows + 1
  std::vector< int > jcol;    ///< local column indices
  std::vector< double > val;  ///< coefficient values

  void clear() {
   krow.clear();
   jcol.clear();
   val.clear();
  }

  [[nodiscard]] Index nnz() const {
   return( static_cast< Index >( val.size() ) );
  }
 };

 /// Extracts a CSR matrix directly from SMS++ linear constraints and columns.
 CSRMatrix extract_block_matrix(
  const std::vector< const FRowConstraint * > & rows ,
  const std::vector< const ColVariable * > & cols
  ) const;

 /// Extracts a CSR submatrix from MILPSolver's global column-wise matrix.
 CSRMatrix extractSubmatrixToCRS(
  const std::vector< int > & selectedRows ,
  const int nRows ,
  const std::vector< int > & selectedCols ,
  const int nCols );

 /// Copies a selected matrix block into the arrays supplied by PIPS.
 int ExtractMatrix( int id , int * krowM , int * jcolM , double * M ,
                    std::vector< const FRowConstraint * > node_cons ,
                    const int nCons ,
                    std::vector< const ColVariable * > vars ,
                    const int nVars );

 /// Counts the nonzeros of a selected PIPS matrix block.
 int EvaluateNnz( int id , int * nnz ,
                  std::vector< const FRowConstraint * > node_cons ,
                  const int nCons ,
                  std::vector< const ColVariable * > vars ,
                  const int nVars );

 /// Computes global MILPSolver row indices for a set of constraints.
 std::vector< int > compute_cons_global_idxs(
  std::vector< const FRowConstraint * > cons , const int nCons );

 /// Computes global MILPSolver column indices for a set of variables.
 std::vector< int > compute_vars_global_idxs(
  std::vector< const ColVariable * > vars , const int nVars );

/*--------------------------------------------------------------------------*/
/*-------------------------- TREE/SCAN HELPERS -----------------------------*/
/*--------------------------------------------------------------------------*/

 /// Recursively collects a Block subtree and maps each Block to its PIPS leaf.
 int collect_subtree( Block * block , std::vector< Block * > & subtree ,
                        Index parent_leaf );

 /// Scans an SMS++ variable/constraint group and dispatches to the right scan.
 template< typename T >
 void scan_group( const boost::any & gr , Block * qb , Index num_node ,
                  un_any_type< T > );

 /// Scans a simple SMS++ variable/constraint group.
 template< typename T >
 void scan_simple_group( const boost::any & gr , Block * qb , Index num_node ,
                         un_any_type< T > );

 template< typename T >
 void scan_multiarray_group( const boost::any & gr , Block * qb , 
                                Index num_node , un_any_type< T > );

 /// Classifies one constraint as node-local or global-linking.
 void scan_constraint( const FRowConstraint & con , Index num_node );

/*--------------------------------------------------------------------------*/
/*----------------------- VECTOR EXTRACTION HELPERS ------------------------*/
/*--------------------------------------------------------------------------*/

 /// Extracts equality RHS or inequality upper-bound values for selected rows.
 int ExtractRhsVector( int id , double * vec , int len ,
                       std::vector< const FRowConstraint * > node_cons ,
                       const int nCons ,
                       std::vector< double > rhs ,
                       std::vector< char > sense ,
                       std::vector< double > ranges );

 /// Extracts inequality lower-bound values for selected rows.
 int ExtractLhsVector( int id , double * vec , int len ,
                       std::vector< const FRowConstraint * > node_cons ,
                       const int nCons ,
                       std::vector< double > rhs ,
                       std::vector< char > sense ,
                       std::vector< double > ranges );

 /// Extracts active flags for selected row upper bounds.
 int ExtractRhsActiveFlag( int id , double * vec , int len ,
                           std::vector< const FRowConstraint * > node_cons ,
                           const int nCons ,
                           std::vector< double > rhs ,
                           std::vector< char > sense ,
                           std::vector< double > ranges );

 /// Extracts active flags for selected row lower bounds.
 int ExtractLhsActiveFlag( int id , double * vec , int len ,
                           std::vector< const FRowConstraint * > node_cons ,
                           const int nCons ,
                           std::vector< double > rhs ,
                           std::vector< char > sense ,
                           std::vector< double > ranges );

 /// Extracts variable bounds for a node.
 int ExtractVarBounds( int id , double * vec , int len ,
                       std::vector< const ColVariable * > node_vars ,
                       const int nVars ,
                       std::vector< double > bounds );

 /// Extracts active flags for variable bounds.
 int ExtractFlagVarBounds( int id , double * vec , int len ,
                            std::vector< const ColVariable * > node_vars ,
                            const int nVars ,
                            std::vector< double > bounds );

 /// Extracts variable objective coefficients.
 int ExtractObj( int id , double* vec , int len ,
                                std::vector< const ColVariable * > node_vars ,
                                const int nVars ,
                                std::vector< double > obj_value ,
                                int objsense );

/*--------------------------------------------------------------------------*/
/*--------------------------- PIPS CALLBACKS -------------------------------*/
/*--------------------------------------------------------------------------*/

 /** @name Dimension callbacks
  * These callbacks tell PIPS how many variables/rows each node owns before
  * matrix and vector data are requested.
  * @{ */

 /// Number of variables in node id.
 static int No_VarinNode( void * user_data , int id , int * nnz );

 /// Number of node-local equality rows in node id.
 static int No_EqConsinNode( void * user_data , int id , int * nnz );

 /// Number of node-local inequality rows in node id.
 static int No_InEqConsinNode( void * user_data , int id , int * nnz );

 /// Number of global linking equality rows shared by all nodes.
 static int No_LinkEqCons( void * user_data , int id , int * nnz );

 /// Number of global linking inequality rows shared by all nodes.
 static int No_LinkInEqCons( void * user_data , int id , int * nnz );

 /** @} */

 /** @name Nonzero-count callbacks
  * Each callback returns the number of nonzeros of the corresponding matrix
  * callback, allowing PIPS to allocate arrays of the right size.
  * @{ */

 /// Nonzeros of the diagonal/local equality block: A for root, D_i for child.
 static int nnzEqConsDiag( void * user_data , int id , int * nnz );

 /// Nonzeros of the vertical equality block: C_i, root-variable part of child rows.
 static int nnzEqConsVert( void * user_data , int id , int * nnz );

 /// Nonzeros of the diagonal/local inequality block.
 static int nnzInEqConsDiag( void * user_data , int id , int * nnz );

 /// Nonzeros of the vertical inequality block.
 static int nnzInEqConsVert( void * user_data , int id , int * nnz );

 /// Nonzeros of the global linking equality contribution of node id: B_i.
 static int nnzLinkEqCons( void * user_data , int id , int * nnz );

 /// Nonzeros of the global linking inequality contribution of node id.
 static int nnzLinkInEqCons( void * user_data , int id , int * nnz );

 /// Nonzeros of the quadratic objective matrix Q; zero for LPs.
 static int nnzAllZero( void * user_data , int id , int * nnz );

 /** @} */

 /** @name Matrix callbacks
  * PIPS expects all matrices in CSR row-major format: krowM, jcolM, M.
  * @{ */

 /// Diagonal/local equality matrix: A for root, D_i for child node i.
 static int MatEqConsDiag( void * user_data , int id , int * krowM ,
                           int * jcolM , double * M );

 /// Vertical equality matrix: C_i, coefficients of root variables in child rows.
 static int MatEqConsVert( void * user_data , int id , int * krowM ,
                           int * jcolM , double * M );

 /// Diagonal/local inequality matrix.
 static int MatInEqConsDiag( void * user_data , int id , int * krowM ,
                             int * jcolM , double * M );

 /// Vertical inequality matrix: root-variable part of child inequalities.
 static int MatInEqConsVert( void * user_data , int id , int * krowM ,
                             int * jcolM , double * M );

 /// Linking equality matrix contribution of node id to the global linking rows.
 static int MatLinkEqCons( void * user_data , int id , int * krowM ,
                           int * jcolM , double * M );

 /// Linking inequality matrix contribution of node id to global linking rows.
 static int MatLinkInEqCons( void * user_data , int id , int * krowM ,
                             int * jcolM , double * M );

 /// Empty matrix callback used for zero quadratic objective Q in LPs.
 static int matAllZero( void * user_data , int id , int * krowM ,
                        int * jcolM , double * M );

 /** @} */

 /** @name Vector callbacks
  * These callbacks provide objective coefficients, RHS values, row bounds,
  * variable bounds, and their active flags.
  * @{ */

 /// Linear objective coefficients for node variables.
 static int ObjVars( void * user_data , int id , double * vec , int len );

 /// RHS of node-local equality rows.
 static int RhsEqCons( void * user_data , int id , double * vec , int len );

 /// Upper bounds of node-local inequality rows.
 static int RhsInEqCons( void * user_data , int id , double * vec , int len );

 /// Lower bounds of node-local inequality rows.
 static int LhsInEqCons( void * user_data , int id , double * vec , int len );

 /// RHS of global linking equality rows.
 static int RhsLinkEqCons( void * user_data , int id , double * vec , int len );

 /// Upper bounds of global linking inequality rows.
 static int RhsLinkInEqCons( void * user_data , int id , double * vec , int len );

 /// Lower bounds of global linking inequality rows.
 static int LhsLinkInEqCons( void * user_data , int id , double * vec , int len );

 /// Active flags for upper bounds of node-local inequality rows.
 static int FlagRhsInEqCons( void * user_data , int id , double * vec , int len );

 /// Active flags for lower bounds of node-local inequality rows.
 static int FlagLhsInEqCons( void * user_data , int id , double * vec , int len );

 /// Active flags for upper bounds of linking inequality rows.
 static int FlagRhsLinkInEqCons( void * user_data , int id , double * vec , int len );

 /// Active flags for lower bounds of linking inequality rows.
 static int FlagLhsLinkInEqCons( void * user_data , int id , double * vec , int len );

 /// Upper bounds of node variables.
 static int UBVars( void * user_data , int id , double * vec , int len );

 /// Lower bounds of node variables.
 static int LBVars( void * user_data , int id , double * vec , int len );

 /// Active flags for variable upper bounds.
 static int FlagUBVars( void * user_data , int id , double * vec , int len );

 /// Active flags for variable lower bounds.
 static int FlagLBVars( void * user_data , int id , double * vec , int len );

 /** @} */

/*--------------------------------------------------------------------------*/
/*------------------------- PIPS INTERNAL DATA -----------------------------*/
/*--------------------------------------------------------------------------*/

 /// PIPS distributed problem tree owned by this solver.
 DistributedInputTree * pips_tree;

 /// PIPS solver interface owned by this solver.
 PIPSIPMppInterface * pips_interface;

 /// Number of PIPS nodes: 0 is root, 1..n_nodes-1 are leaves.
 Index n_nodes = 0;

 /// Block subtrees grouped by PIPS node.
 std::vector< std::vector< Block * > > nodes_subtrees;

 /// Map from each SMS++ Block to the PIPS leaf containing it.
 std::unordered_map< const Block * , std::size_t > blockToleaf;

 /// Number of variables in each node.
 std::vector< int > n_varNode;

 /// Variables assigned to each node.
 std::vector< std::vector< const ColVariable * > > varNode;

 /// Number of node-local equality constraints in each node.
 std::vector< int > n_EqConsNode;

 /// Node-local equality constraints.
 std::vector< std::vector< const FRowConstraint * > > EqConsNode;

 /// Number of node-local inequality constraints in each node.
 std::vector< int > n_InEqConsNode;

 /// Node-local inequality constraints.
 std::vector< std::vector< const FRowConstraint * > > InEqConsNode;

 /// Number of global linking equality constraints.
 int n_LinkEqCons = 0;

 /// Global linking equality constraints.
 std::vector< const FRowConstraint * > LinkEqCons;

 /// Number of global linking inequality constraints.
 int n_LinkInEqCons = 0;

 /// Global linking inequality constraints.
 std::vector< const FRowConstraint * > LinkInEqCons;

 bool mpi_initialized_by_this_solver = false;

/*--------------------------------------------------------------------------*/
/*----------------------------- OTHER DATA ---------------------------------*/
/*--------------------------------------------------------------------------*/

 double UpCutOff;  ///< externally set upper cutoff
 double LwCutOff;  ///< externally set lower cutoff

 /** vector containing the filenames used to load of the Configuration of
  * the "Configuration DB" */
 std::vector< std::string > ConfigDBFName;

 /// the "Configuration DB" istself
 std::vector< Configuration * > v_ConfigDB;

 /// the mutex to ensure that SCIP threads do not overstep in the callback
 std::mutex f_callback_mutex;

/*--------------------------------------------------------------------------*/
/*---------------------- PRIVATE PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

 private:

 /// Converts PIPS-IPM++ status codes into SMS++ solver status codes.
 static int decode_pips_status( TerminationStatus status );

 SMSpp_insert_in_factory_h;

};  // end( class( PIPSMILPSolver ) )

/*--------------------------------------------------------------------------*/

}  // end( namespace SMSpp_di_unipi_it )

/*--------------------------------------------------------------------------*/

#endif  /* __PIPSMILPSOLVER_H */

/*--------------------------------------------------------------------------*/
/*------------------------ End File PIPSMILPSolver.h -----------------------*/
/*--------------------------------------------------------------------------*/