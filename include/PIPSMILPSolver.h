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
#include "PIPSIPMppOptions.h"
#include "MILPSolver.h"

// Include the proper PIPS parameter mapping
#include "PIPS_defs.h"

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

 /// enum for integer parameters
 enum int_par_type_PIPS {
  // note: intCutSepPar has moved to MILPSolver base (enum int_par_type_MILP)
  intFirstPIPSPar = intLastAlgParMILP ,  ///< first PIPS int/long parameter
  /// first allowed new int parameter for derived classes
  intLastAlgParPIPS = intFirstPIPSPar + PIPS_NUM_INT_PARS
  };

 /// enum for double parameters
 enum dbl_par_type_PIPS {
  /// first PIPS double parameter
  dblFirstPIPSPar = dblLastAlgParMILP,
  /// first allowed new double parameter for derived classes
  dblLastAlgParPIPS = dblFirstPIPSPar + PIPS_NUM_DBL_PARS
  };

 /// enum for string parameters
 enum str_par_type_PIPS {
  /// first PIPS string parameter
  strFirstPIPSPar = strLastAlgParMILP,
  /// first allowed new string parameter for derived classes
  strLastAlgParPIPS = strFirstPIPSPar + PIPS_NUM_STR_PARS
  };

/*--------------------------------------------------------------------------*/
 // "importing" a few types from Block

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

 // note: the public compute() entry point is inherited from MILPSolver;
 // PIPSMILPSolver implements only the PIPS-specific solve in
 // guts_of_compute() below (protected)

 /// returns a valid lower bound on the optimal objective function value
 OFValue get_lb( void ) override;

 /// returns a valid upper bound on the optimal objective function value
 OFValue get_ub( void ) override;

 /// writes the current solution in the Block
 void get_var_solution( Configuration * solc = nullptr ) override;

 /// Clears the current PIPS tree/interface and the base MILPSolver data.
 void clear_problem( unsigned int what ) override;

 /** @} ---------------------------------------------------------------------*/
/*------------------- METHODS FOR HANDLING THE PARAMETERS ------------------*/
/*--------------------------------------------------------------------------*/
/** @name Methods for handling parameters
 * @{ */

 /// sets an integer parameter with the given value

 void set_par( idx_type par , int value ) override;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// sets a double parameter with the given value
 void set_par( idx_type par , double value ) override;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// sets a string parameter with the given value
 void set_par( idx_type par , std::string && value ) override;

/*--------------------------------------------------------------------------*/
 /// returns the number of integer parameters
 [[nodiscard]] idx_type get_num_int_par( void ) const override;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// returns the number of double parameters
 [[nodiscard]] idx_type get_num_dbl_par( void ) const override;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// returns the number of string parameters
 [[nodiscard]] idx_type get_num_str_par( void ) const override;

/*--------------------------------------------------------------------------*/
 /// returns the default value of the specified integer parameter
 [[nodiscard]] int get_dflt_int_par( idx_type par ) const override;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// returns the default value of the specified double parameter
 [[nodiscard]] double get_dflt_dbl_par( idx_type par ) const override;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /** Returns the default value of the specified string parameter
  * @note
  * Due to a limit in the implementation, the string referenced by
  * the return value is *overwritten* each time the method is called with
  * par as a PIPS parameter. */
 [[nodiscard]] const std::string & get_dflt_str_par( idx_type par )
  const override;

/*--------------------------------------------------------------------------*/
 /// returns the value of the specified integer parameter
 [[nodiscard]] int get_int_par( idx_type par ) const override;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// returns the value of the specified double parameter
 [[nodiscard]] double get_dbl_par( idx_type par ) const override;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /** returns the value of the specified string parameter
  * @note
  * Due to a limit in the implementation, the string referenced by
  * the return value is *overwritten* each time the method is called with
  * par as a PIPS parameter. */
 [[nodiscard]] const std::string & get_str_par( idx_type par ) const override;
 
/*--------------------------------------------------------------------------*/
 /// returns the index of the int parameter with the specified name
 [[nodiscard]] idx_type int_par_str2idx( const std::string & name )
  const override;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /** returns the name of the int parameter with the specified index
  * @note
  * Due to a limit in the implementation, the string referenced by
  * the return value is *overwritten* each time the method is called with
  * par as a PIPS parameter. */
 [[nodiscard]] const std::string & int_par_idx2str( idx_type idx )
  const override;

/*--------------------------------------------------------------------------*/
 /// Returns the index of the double parameter with the specified name
 [[nodiscard]] idx_type dbl_par_str2idx( const std::string & name )
  const override;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /** returns the name of the double parameter with the specified index
  * @note
  * Due to a limit in the implementation, the string referenced by
  * the return value is *overwritten* each time the method is called with
  * par as a PIPS parameter. */
 [[nodiscard]] const std::string & dbl_par_idx2str( idx_type idx )
  const override;

/*--------------------------------------------------------------------------*/
 /// returns the index of the string parameter with the specified name
 [[nodiscard]] idx_type
  str_par_str2idx( const std::string & name ) const override;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /** Returns the name of the string parameter with the specified index
  * @note
  * Due to a limit in the implementation, the string referenced by
  * the return value is *overwritten* each time the method is called with
  * par as a PIPS parameter. */
 [[nodiscard]] const std::string & str_par_idx2str( idx_type idx )
  const override;

 protected:

/*--------------------------------------------------------------------------*/
/*-------------------- PROTECTED METHODS OF THE CLASS ----------------------*/
/*--------------------------------------------------------------------------*/

/// PIPS-IPM++ back-end solve, called by MILPSolver::compute()
 /** Performs the actual PIPS-IPM++ optimisation. Locking, Modification
  * processing and the LP cut-separation loop (intRelaxIntVars == 2) are
  * all handled by MILPSolver::compute(). */

 int guts_of_compute( void ) override;

/** @} ---------------------------------------------------------------------*/
 /// maps a Solver integer parameter into a PIPS one
 /** Maps the Solver integer parameter \p par into a PIPS one;
  * returns a string with the PIPS parameter name */
 std::string pips_int_par_map( idx_type par ) const;
 
/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// maps a Solver double parameter into a PIPS one
 std::string pips_dbl_par_map( idx_type par ) const;

/** @name Handling of PIPS parameters
  *
  * The following maps are used to keep a relationship between SMS++ parameter
  * system and PIPS parameters. This allows us to use PIPS parameters
  * (See https://pips-ipmpp.gitlab.io/GMSPIPS.html#OPTIONS) as they were SMS++
  * parameters with the same names, for example in configuration files.
  *
  * Note: since SMS++ does not support bool parameters, both int and
  *       bool PIPS parameters are handled as SMS++ int parameters.
  * @{ */

 const static std::array< std::string , PIPS_NUM_INT_PARS > SMSpp_to_PIPS_int_pars;
 const static std::array< std::string , PIPS_NUM_DBL_PARS > SMSpp_to_PIPS_dbl_pars;
 const static std::array< std::string , PIPS_NUM_STR_PARS > SMSpp_to_PIPS_str_pars;

 const static std::array< std::pair< std::string , int > , PIPS_NUM_INT_PARS >
  PIPS_to_SMSpp_int_pars;
 const static std::array< std::pair< std::string , int > , PIPS_NUM_DBL_PARS >
  PIPS_to_SMSpp_dbl_pars;
 const static std::array< std::pair< std::string , int > , PIPS_NUM_STR_PARS >
  PIPS_to_SMSpp_str_pars;

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
  const int nCols ) const;

 /// Extracts and validates a matrix block while building the callback cache.
 CSRMatrix build_cached_matrix(
  int id , const char * name ,
  const std::vector< const FRowConstraint * > & rows ,
  const std::vector< const ColVariable * > & cols ) const;

 /// Builds every matrix block requested by PIPS callbacks.
 void build_matrix_cache();

 /// Copies a cached matrix block into the arrays supplied by PIPS.
 static int copy_cached_matrix( const CSRMatrix & matrix , int * krowM ,
                                int * jcolM , double * M );

 /// Computes global MILPSolver row indices for a set of constraints.
 std::vector< int > compute_cons_global_idxs(
  const std::vector< const FRowConstraint * > & cons ) const;

 /// Computes global MILPSolver column indices for a set of variables.
 std::vector< int > compute_vars_global_idxs(
  const std::vector< const ColVariable * > & vars ) const;

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
                       const std::vector< const FRowConstraint * > & node_cons ,
                       const std::vector< double > & rhs ,
                       const std::vector< char > & sense ,
                       const std::vector< double > & ranges );

 /// Extracts inequality lower-bound values for selected rows.
 int ExtractLhsVector( int id , double * vec , int len ,
                       const std::vector< const FRowConstraint * > & node_cons ,
                       const std::vector< double > & rhs ,
                       const std::vector< char > & sense ,
                       const std::vector< double > & ranges );

 /// Extracts active flags for selected row upper bounds.
 int ExtractRhsActiveFlag( int id , double * vec , int len ,
                       const std::vector< const FRowConstraint * > & node_cons ,
                       const std::vector< double > & rhs ,
                       const std::vector< char > & sense ,
                       const std::vector< double > & ranges );

 /// Extracts active flags for selected row lower bounds.
 int ExtractLhsActiveFlag( int id , double * vec , int len ,
                       const std::vector< const FRowConstraint * > & node_cons ,
                       const std::vector< double > & rhs ,
                       const std::vector< char > & sense ,
                       const std::vector< double > & ranges );

 /// Extracts variable bounds for a node.
 int ExtractVarBounds( int id , double * vec , int len ,
                       const std::vector< const ColVariable * > & node_vars ,
                       const std::vector< double > & bounds );

 /// Extracts active flags for variable bounds.
 int ExtractFlagVarBounds( int id , double * vec , int len ,
                       const std::vector< const ColVariable * > & node_vars ,
                       const std::vector< double > & bounds );

 /// Extracts variable objective coefficients.
 int ExtractObj( int id , double* vec , int len ,
                  const std::vector< const ColVariable * > & node_vars ,
                  const std::vector< double > & obj_value , int objsense );

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
 pipsipmpp::DistributedInputTree * pips_tree;

 /// PIPS solver interface owned by this solver.
 pipsipmpp::PIPSIPMppInterface * pips_interface;

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

 /// CSR blocks prepared once per load and copied by the PIPS callbacks.
 struct NodeMatrixCache {
  CSRMatrix eq_diag;
  CSRMatrix eq_vert;
  CSRMatrix ineq_diag;
  CSRMatrix ineq_vert;
  CSRMatrix link_eq;
  CSRMatrix link_ineq;
 };

 std::vector< NodeMatrixCache > matrix_cache;

 bool mpi_initialized_by_this_solver = false;

 std::unordered_map< const ColVariable *, int > var_to_node;

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

 /// Clears all PIPS-owned objects and data derived from the current Block.
 void reset_pips_data();

 /// Converts PIPS-IPM++ status codes into SMS++ solver status codes.
 static int decode_pips_status( pipsipmpp::TerminationStatus status );

 SMSpp_insert_in_factory_h;

};  // end( class( PIPSMILPSolver ) )

/*--------------------------------------------------------------------------*/

}  // end( namespace SMSpp_di_unipi_it )

/*--------------------------------------------------------------------------*/

#endif  /* __PIPSMILPSOLVER_H */

/*--------------------------------------------------------------------------*/
/*------------------------ End File PIPSMILPSolver.h -----------------------*/
/*--------------------------------------------------------------------------*/
