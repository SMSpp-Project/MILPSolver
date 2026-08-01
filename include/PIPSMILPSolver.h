/*--------------------------------------------------------------------------*/
/*-------------------------- File PIPSMILPSolver.h -------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the PIPSMILPSolver class, implementing the MILPSolver
 * interface for LP problems using the PIPS-IPM++ parallel interior-point
 * solver.
 *
 * \author Enrico Calandrini \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \author Antonio Frangioni \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \copyright &copy; by Enrico Calandrini, Antonio Frangioni
 */
/*--------------------------------------------------------------------------*/
/*----------------------------- DEFINITIONS --------------------------------*/
/*--------------------------------------------------------------------------*/

#ifndef __PIPSMILPSOLVER_H
 #define __PIPSMILPSOLVER_H
                      /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include <array>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include <mpi.h>

#include <DistributedInputTree.h>
#include <PIPSIPMppInterface.hpp>
#include <PIPSIPMppOptions.h>

#include "MILPSolver.h"

#include "PIPS_defs.h"

/*--------------------------------------------------------------------------*/
/*--------------------------- NAMESPACE ------------------------------------*/
/*--------------------------------------------------------------------------*/

/// namespace for the Structured Modeling System++ (SMS++)
namespace SMSpp_di_unipi_it
{
 class LinearFunction;  // forward declaration of LinearFunction

/*--------------------------------------------------------------------------*/
/*------------------------ CLASS PIPSMILPSolver ----------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------- GENERAL NOTES --------------------------------*/
/*--------------------------------------------------------------------------*/
/// class for solving LP Blocks via the PIPS-IPM++ parallel IP solver
/** The PIPSMILPSolver implements the MILPSolver interface using PIPS-IPM++,
 * an open-source, MPI-parallel interior-point solver for block-structured
 * linear programs. It is registered in the Solver factory, hence it can be
 * selected by name ("PIPSMILPSolver") in a BlockSolverConfig.
 *
 * PIPS-IPM++ solves problems with a dual block-angular structure. The Block
 * hierarchy is mapped onto a two-level PIPS tree as follows: the root node
 * holds the Variable and Constraint of the Block itself, and each nested
 * Block (together with its whole sub-tree of Block) becomes one leaf. The
 * Constraint of a Block are automatically classified as node-local or
 * linking according to the Variable they involve; all matrix and vector
 * data is provided to PIPS through callbacks that extract it from the
 * standard MILPSolver representation.
 *
 * The following restrictions apply:
 *
 * - only LP models are supported: models with a quadratic objective or
 *   quadratic constraints cause load_problem() to throw, and integer
 *   variables are only accepted if relaxed via intRelaxIntVars;
 *
 * - the Block needs at least one nested Block (i.e., the PIPS tree needs
 *   at least one leaf);
 *
 * - PIPS-IPM++ does not support changing the loaded model: whenever any
 *   Modification is received, the model (and the PIPS tree with it) is
 *   cleared and reconstructed from scratch at the next compute();
 *
 * - the solver runs on MPI_COMM_WORLD. MPI is initialised on demand at the
 *   first set_Block() (with the funneled threading level) if the
 *   application has not done it already, in which case it is finalised at
 *   program exit; it is never finalised by the destructor. In an SPMD run
 *   the primal and dual solutions are written into the Block on rank 0
 *   only;
 *
 * - no unbounded ray or infeasibility certificate is provided, i.e., the
 *   has_var_direction() / has_dual_direction() interface of CDASolver is
 *   not implemented;
 *
 * - the PIPS options are process-global, hence shared by all the
 *   PIPSMILPSolver instances in the same process; see the parameter
 *   handling methods for details. */

class PIPSMILPSolver : public MILPSolver
{
/*--------------------------------------------------------------------------*/
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

 public:

/*--------------------------------------------------------------------------*/
/*---------------------------- PUBLIC TYPES --------------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Public types of PIPSMILPSolver
 * @{ */

 /// public enum for the integer algorithmic parameters
 /** Public enum describing the different algorithmic parameters of int type
  * that PIPSMILPSolver has in addition to those of MILPSolver. The value
  * intLastAlgParPIPS is provided so that the list can be easily further
  * extended by derived classes. */
 enum int_par_type_PIPS {
  intFirstPIPSPar = intLastAlgParMILP ,  ///< first PIPS int parameter
  /// first allowed new int parameter for derived classes
  intLastAlgParPIPS = intFirstPIPSPar + PIPS_NUM_INT_PARS
  };

 /// public enum for the double algorithmic parameters
 /** Public enum describing the different algorithmic parameters of double
  * type that PIPSMILPSolver has in addition to those of MILPSolver. The
  * value dblLastAlgParPIPS is provided so that the list can be easily
  * further extended by derived classes. */
 enum dbl_par_type_PIPS {
  dblFirstPIPSPar = dblLastAlgParMILP ,  ///< first PIPS double parameter
  /// first allowed new double parameter for derived classes
  dblLastAlgParPIPS = dblFirstPIPSPar + PIPS_NUM_DBL_PARS
  };

 /// public enum for the string algorithmic parameters
 /** Public enum describing the different algorithmic parameters of string
  * type that PIPSMILPSolver has in addition to those of MILPSolver. The
  * value strLastAlgParPIPS is provided so that the list can be easily
  * further extended by derived classes. */
 enum str_par_type_PIPS {
  strFirstPIPSPar = strLastAlgParMILP ,  ///< first PIPS string parameter
  /// first allowed new string parameter for derived classes
  strLastAlgParPIPS = strFirstPIPSPar + PIPS_NUM_STR_PARS
  };

/*--------------------------------------------------------------------------*/
 // "importing" a few types from Block

 using Subset = Block::Subset;
 using c_Subset = Block::c_Subset;
 using Range = Block::Range;

 using Index = Block::Index;

/** @} ---------------------------------------------------------------------*/
/*------------------- CONSTRUCTING AND DESTRUCTING -------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Constructing and destructing PIPSMILPSolver
 * @{ */

 /// void constructor
 PIPSMILPSolver( void );

/*--------------------------------------------------------------------------*/
 /// destructor: releases the PIPS tree and interface
 /** Destructor. Note that it does *not* finalise MPI even if it was this
  * class that initialised it, since other instances may still be created
  * afterwards; MPI is finalised at program exit instead. */
 ~PIPSMILPSolver() override;

/** @} ---------------------------------------------------------------------*/
/*--------------------- DERIVED METHODS OF BASE CLASS ----------------------*/
/*--------------------------------------------------------------------------*/
/** @name Public methods derived from base classes
 * @{ */

 /// sets the Block that the Solver has to solve
 /** Sets the Block that the Solver has to solve; besides the standard
  * MILPSolver processing, MPI is initialised here if the application has
  * not done it already (see the general notes). */
 void set_Block( Block * block ) override;

/*--------------------------------------------------------------------------*/
 /// loads the problem from the Block and builds the PIPS input tree
 /** Loads the LP problem from the Block into the standard MILPSolver
  * representation, then builds the PIPS DistributedInputTree out of it (see
  * the general notes for the mapping). It throws if the model has a
  * quadratic objective or quadratic constraints, if (non-relaxed) integer
  * variables are present, or if the Block has no nested Block. */
 void load_problem( void ) override;

 // note: the public compute() entry point is inherited from MILPSolver;
 // PIPSMILPSolver only implements the PIPS-specific part in the protected
 // guts_of_compute()

/*--------------------------------------------------------------------------*/
 /// returns a valid lower bound on the optimal objective function value
 OFValue get_lb( void ) override;

 /// returns a valid upper bound on the optimal objective function value
 OFValue get_ub( void ) override;

/*--------------------------------------------------------------------------*/
 /// writes the current solution in the Block
 /** Writes the current primal solution in the Block. In an SPMD run the
  * solution is gathered on rank 0 and written into the Block there only;
  * on the other ranks the Variable of the Block are left untouched. */
 void get_var_solution( Configuration * solc = nullptr ) override;

 /// returns true if a dual solution is available
 bool has_dual_solution( void ) override;

 /// returns true if the current dual solution is dual feasible
 bool is_dual_feasible( void ) override;

 /// writes the current dual solution in the Block
 /** Writes the current dual solution in the Block; the same rank 0
  * restriction as in get_var_solution() applies. */
 void get_dual_solution( Configuration * solc = nullptr ) override;

/*--------------------------------------------------------------------------*/
 /// clears the PIPS tree/interface together with the base MILPSolver data
 void clear_problem( unsigned int what ) override;

/** @} ---------------------------------------------------------------------*/
/*------------------- METHODS FOR HANDLING THE PARAMETERS ------------------*/
/*--------------------------------------------------------------------------*/
/** @name Methods for handling the parameters of the PIPSMILPSolver
 *
 * PIPSMILPSolver supports, besides the parameters of MILPSolver, all the
 * options of PIPS-IPM++ under their native names (as found in the PIPS
 * options documentation), for instance in Configuration files; the mapping
 * is generated by the pips_pars tool. Since SMS++ has no bool parameters,
 * bool PIPS options are handled as int parameters.
 *
 * Note that the PIPS options are process-global: they are shared by all
 * the PIPSMILPSolver instances in the same process, and setting them is
 * not thread safe.
 * @{ */

 /// sets an integer parameter with the given value
 /** Sets an integer parameter with the given value. Besides the PIPS
  * options, the following MILPSolver / Solver parameters are mapped into
  * PIPS options:
  *
  * - intMaxIter is mapped into "IPM_MAX_ITER";
  *
  * - intLogVerb is mapped into "SILENT" (with inverted value, since the
  *   PIPS option has the opposite meaning).
  *
  * Any other parameter is handled by the base class; in particular,
  * parameters that have no PIPS counterpart (say, dblUpCutOff and
  * dblLwCutOff) are silently ignored. */
 void set_par( idx_type par , int value ) override;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// sets a double parameter with the given value
 /** Sets a double parameter with the given value; dblMaxTime is mapped
  * into "IPM_TIMELIMIT" and dblRelAcc into "OUTER_BICG_TOL", any other
  * non-PIPS parameter is handled by the base class. */
 void set_par( idx_type par , double value ) override;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// sets a string parameter with the given value
 /** Sets a string parameter with the given value; strOutputFile also sets
  * the PIPS option "WRITE_ORIGINAL_PROBLEM_TO_LP", any non-PIPS parameter
  * is handled by the base class. */
 void set_par( idx_type par , std::string && value ) override;

/*--------------------------------------------------------------------------*/
 /// returns the number of integer parameters
 [[nodiscard]] idx_type get_num_int_par( void ) const override;

 /// returns the number of double parameters
 [[nodiscard]] idx_type get_num_dbl_par( void ) const override;

 /// returns the number of string parameters
 [[nodiscard]] idx_type get_num_str_par( void ) const override;

/*--------------------------------------------------------------------------*/
 /// returns the default value of the specified integer parameter
 /** Returns the default value of the specified integer parameter; for
  * parameters mapped into PIPS options the current value of the
  * (process-global) option is returned, since PIPS does not expose the
  * default values of its options. */
 [[nodiscard]] int get_dflt_int_par( idx_type par ) const override;

 /// returns the default value of the specified double parameter
 /** Returns the default value of the specified double parameter; the same
  * note as in get_dflt_int_par() applies. */
 [[nodiscard]] double get_dflt_dbl_par( idx_type par ) const override;

 /// returns the default value of the specified string parameter
 /** Returns the default value of the specified string parameter; the same
  * note as in get_dflt_int_par() applies.
  *
  * @note the string referenced by the return value is *overwritten* each
  *       time the method is called with \p par a PIPS parameter. */
 [[nodiscard]] const std::string & get_dflt_str_par( idx_type par )
  const override;

/*--------------------------------------------------------------------------*/
 /// returns the value of the specified integer parameter
 [[nodiscard]] int get_int_par( idx_type par ) const override;

 /// returns the value of the specified double parameter
 [[nodiscard]] double get_dbl_par( idx_type par ) const override;

 /// returns the value of the specified string parameter
 /** Returns the value of the specified string parameter.
  *
  * @note the string referenced by the return value is *overwritten* each
  *       time the method is called with \p par a PIPS parameter. */
 [[nodiscard]] const std::string & get_str_par( idx_type par )
  const override;

/*--------------------------------------------------------------------------*/
 /// returns the index of the int parameter with the specified name
 [[nodiscard]] idx_type int_par_str2idx( const std::string & name )
  const override;

 /// returns the name of the int parameter with the specified index
 /** Returns the name of the int parameter with the specified index.
  *
  * @note the string referenced by the return value is *overwritten* each
  *       time the method is called with \p idx a PIPS parameter. */
 [[nodiscard]] const std::string & int_par_idx2str( idx_type idx )
  const override;

/*--------------------------------------------------------------------------*/
 /// returns the index of the double parameter with the specified name
 [[nodiscard]] idx_type dbl_par_str2idx( const std::string & name )
  const override;

 /// returns the name of the double parameter with the specified index
 /** Returns the name of the double parameter with the specified index.
  *
  * @note the string referenced by the return value is *overwritten* each
  *       time the method is called with \p idx a PIPS parameter. */
 [[nodiscard]] const std::string & dbl_par_idx2str( idx_type idx )
  const override;

/*--------------------------------------------------------------------------*/
 /// returns the index of the string parameter with the specified name
 [[nodiscard]] idx_type str_par_str2idx( const std::string & name )
  const override;

 /// returns the name of the string parameter with the specified index
 /** Returns the name of the string parameter with the specified index.
  *
  * @note the string referenced by the return value is *overwritten* each
  *       time the method is called with \p idx a PIPS parameter. */
 [[nodiscard]] const std::string & str_par_idx2str( idx_type idx )
  const override;

/** @} ---------------------------------------------------------------------*/
/*--------------------- PROTECTED PART OF THE CLASS ------------------------*/
/*--------------------------------------------------------------------------*/

 protected:

/*--------------------------------------------------------------------------*/
/*-------------------- PROTECTED METHODS OF THE CLASS ----------------------*/
/*--------------------------------------------------------------------------*/
/** @name Solving the problem
 * @{ */

 /// PIPS-IPM++ back-end solve, called by MILPSolver::compute()
 /** Performs the actual PIPS-IPM++ optimisation; locking, Modification
  * processing and the LP cut-separation loop (intRelaxIntVars == 2) are
  * all handled by MILPSolver::compute(). */
 int guts_of_compute( void ) override;

/** @} ---------------------------------------------------------------------*/
/*------------------- METHODS FOR MODIFYING THE PROBLEM --------------------*/
/*--------------------------------------------------------------------------*/
/** @name Methods for modifying the problem
 *
 * PIPS-IPM++ does not support changing the loaded model, hence every
 * Modification marks the model for a full rebuild: the base MILPSolver
 * representation and the PIPS tree are reconstructed from scratch (by
 * load_problem()) before the next solve, and all the other Modification
 * in the queue are discarded.
 * @{ */

 /// handles a Variable Modification by marking the model for a rebuild
 void var_modification( const VariableMod * mod ) override;

 /// handles an Objective Modification by marking the model for a rebuild
 void objective_modification( const ObjectiveMod * mod ) override;

 /// handles a Constraint Modification by marking the model for a rebuild
 void const_modification( const ConstraintMod * mod ) override;

 /// handles a bound Modification by marking the model for a rebuild
 void bound_modification( const OneVarConstraintMod * mod ) override;

 /// handles an Objective Function Modification by marking a rebuild
 void objective_function_modification( const FunctionMod * mod ) override;

 /// handles a Constraint Function Modification by marking a rebuild
 void constraint_function_modification( const FunctionMod * mod ) override;

 /// handles an Objective FunctionModVars by marking the model for a rebuild
 void objective_fvars_modification( const FunctionModVars * mod ) override;

 /// handles a Constraint FunctionModVars by marking the model for a rebuild
 void constraint_fvars_modification( const FunctionModVars * mod ) override;

 /// handles a dynamic Modification by marking the model for a rebuild
 void dynamic_modification( const BlockModAD * mod ) override;

/** @} ---------------------------------------------------------------------*/
/*------------------- HANDLING OF THE PIPS PARAMETERS ----------------------*/
/*--------------------------------------------------------------------------*/
/** @name Handling of the PIPS parameters
 *
 * The following methods and maps keep the relationship between the SMS++
 * parameter system and the PIPS options, allowing the latter to be used
 * under their native names as if they were SMS++ parameters (say, in
 * Configuration files); since SMS++ does not support bool parameters, both
 * int and bool PIPS options are handled as SMS++ int parameters.
 * @{ */

 /// maps a Solver integer parameter into the corresponding PIPS option
 /** Maps the Solver integer parameter \p par into the corresponding PIPS
  * option, returning its name, or the empty string if there is none. */
 std::string pips_int_par_map( idx_type par ) const;

 /// maps a Solver double parameter into the corresponding PIPS option
 std::string pips_dbl_par_map( idx_type par ) const;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 // parameter maps generated by the pips_pars tool

 const static std::array< std::string , PIPS_NUM_INT_PARS >
  SMSpp_to_PIPS_int_pars;
 const static std::array< std::string , PIPS_NUM_DBL_PARS >
  SMSpp_to_PIPS_dbl_pars;
 const static std::array< std::string , PIPS_NUM_STR_PARS >
  SMSpp_to_PIPS_str_pars;

 const static std::array< std::pair< std::string , int > ,
			  PIPS_NUM_INT_PARS > PIPS_to_SMSpp_int_pars;
 const static std::array< std::pair< std::string , int > ,
			  PIPS_NUM_DBL_PARS > PIPS_to_SMSpp_dbl_pars;
 const static std::array< std::pair< std::string , int > ,
			  PIPS_NUM_STR_PARS > PIPS_to_SMSpp_str_pars;

/** @} ---------------------------------------------------------------------*/
/*------------------------- PIPS MATRIX HELPERS ----------------------------*/
/*--------------------------------------------------------------------------*/
/** @name PIPS matrix helpers
 * @{ */

 /// sparse matrix in CSR row format, as required by PIPS matrix callbacks
 struct CSRMatrix {
  std::vector< int > krow;    ///< row pointers, size = rows + 1
  std::vector< int > jcol;    ///< local column indices
  std::vector< double > val;  ///< coefficient values

  void clear( void ) {
   krow.clear();
   jcol.clear();
   val.clear();
   }

  [[nodiscard]] Index nnz( void ) const {
   return( static_cast< Index >( val.size() ) );
   }
  };

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// extracts a CSR matrix directly from linear constraints and columns
 CSRMatrix extract_block_matrix(
		    const std::vector< const FRowConstraint * > & rows ,
		    const std::vector< const ColVariable * > & cols ) const;

 /// extracts a CSR submatrix from the MILPSolver column-wise matrix
 CSRMatrix extract_submatrix_to_csr(
			      const std::vector< int > & selected_rows ,
			      int n_rows ,
			      const std::vector< int > & selected_cols ,
			      int n_cols ) const;

 /// copies a selected matrix block into the arrays supplied by PIPS
 int extract_matrix( int id , int * krowM , int * jcolM , double * M ,
		     const std::vector< const FRowConstraint * > & node_cons ,
		     const std::vector< const ColVariable * > & vars ) const;

 /// counts the nonzeros of a selected PIPS matrix block
 int evaluate_nnz( int id , int * nnz ,
		   const std::vector< const FRowConstraint * > & node_cons ,
		   const std::vector< const ColVariable * > & vars ) const;

 /// computes the global MILPSolver row indices of a set of constraints
 std::vector< int > compute_cons_global_idxs(
		 const std::vector< const FRowConstraint * > & cons ) const;

 /// computes the global MILPSolver column indices of a set of variables
 std::vector< int > compute_vars_global_idxs(
		    const std::vector< const ColVariable * > & vars ) const;

/** @} ---------------------------------------------------------------------*/
/*-------------------------- TREE/SCAN HELPERS -----------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Tree and scan helpers
 * @{ */

 /// recursively collects a Block subtree and maps each Block to its leaf
 Index collect_subtree( Block * block , std::vector< Block * > & subtree ,
			Index parent_leaf );

 /// scans a group of Variable / Constraint, dispatching to the right scan
 template< typename T >
 void scan_group( const boost::any & gr , Block * qb , Index num_node ,
		  un_any_type< T > );

 /// scans a simple (non multi-array) group of Variable / Constraint
 template< typename T >
 void scan_simple_group( const boost::any & gr , Block * qb ,
			 Index num_node , un_any_type< T > );

 /// scans a multi-array group of Variable / Constraint
 template< typename T >
 void scan_multiarray_group( const boost::any & gr , Block * qb ,
			     Index num_node , un_any_type< T > );

 /// classifies one constraint as node-local or global-linking
 void scan_constraint( const FRowConstraint & con , Index num_node );

/** @} ---------------------------------------------------------------------*/
/*----------------------- VECTOR EXTRACTION HELPERS ------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Vector extraction helpers
 * @{ */

 /// extracts equality rhs or inequality upper bounds for selected rows
 int extract_rhs_vector( int id , double * vec , int len ,
		    const std::vector< const FRowConstraint * > & node_cons ,
		    const std::vector< double > & rhs ,
		    const std::vector< char > & sense ,
		    const std::vector< double > & ranges );

 /// extracts inequality lower bounds for selected rows
 int extract_lhs_vector( int id , double * vec , int len ,
		    const std::vector< const FRowConstraint * > & node_cons ,
		    const std::vector< double > & rhs ,
		    const std::vector< char > & sense ,
		    const std::vector< double > & ranges );

 /// extracts the active flags of the upper bounds of selected rows
 int extract_rhs_active_flag( int id , double * vec , int len ,
		    const std::vector< const FRowConstraint * > & node_cons ,
		    const std::vector< double > & rhs ,
		    const std::vector< char > & sense ,
		    const std::vector< double > & ranges );

 /// extracts the active flags of the lower bounds of selected rows
 int extract_lhs_active_flag( int id , double * vec , int len ,
		    const std::vector< const FRowConstraint * > & node_cons ,
		    const std::vector< double > & rhs ,
		    const std::vector< char > & sense ,
		    const std::vector< double > & ranges );

 /// extracts the variable bounds of a node
 int extract_var_bounds( int id , double * vec , int len ,
		    const std::vector< const ColVariable * > & node_vars ,
		    const std::vector< double > & bounds );

 /// extracts the active flags of the variable bounds of a node
 int extract_flag_var_bounds( int id , double * vec , int len ,
		    const std::vector< const ColVariable * > & node_vars ,
		    const std::vector< double > & bounds );

 /// extracts the objective coefficients of the variables of a node
 int extract_obj( int id , double * vec , int len ,
		  const std::vector< const ColVariable * > & node_vars ,
		  const std::vector< double > & obj_value , int objsense );

/** @} ---------------------------------------------------------------------*/
/*--------------------------- PIPS CALLBACKS -------------------------------*/
/*--------------------------------------------------------------------------*/
/** @name PIPS callbacks
 *
 * PIPS-IPM++ asks the user code for dimensions, matrices and vectors
 * through plain C-style callbacks (see the PipsFNNZ, PipsFMAT and PipsFVEC
 * types of DistributedInputTree); the first argument is an opaque pointer
 * that this class passes as `this`, and every callback casts it back to
 * PIPSMILPSolver * to access the node data built by load_problem().
 * @{ */

 /// number of variables in node id
 static int no_var_in_node( void * user_data , int id , int * nnz );

 /// number of node-local equality rows in node id
 static int no_eq_cons_in_node( void * user_data , int id , int * nnz );

 /// number of node-local inequality rows in node id
 static int no_ineq_cons_in_node( void * user_data , int id , int * nnz );

 /// number of global linking equality rows shared by all nodes
 static int no_link_eq_cons( void * user_data , int id , int * nnz );

 /// number of global linking inequality rows shared by all nodes
 static int no_link_ineq_cons( void * user_data , int id , int * nnz );

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 // nonzero-count callbacks, one for each matrix callback below

 /// nonzeros of the diagonal equality block: A for root, D_i for a leaf
 static int nnz_eq_cons_diag( void * user_data , int id , int * nnz );

 /// nonzeros of the vertical equality block of a leaf
 static int nnz_eq_cons_vert( void * user_data , int id , int * nnz );

 /// nonzeros of the diagonal inequality block
 static int nnz_ineq_cons_diag( void * user_data , int id , int * nnz );

 /// nonzeros of the vertical inequality block
 static int nnz_ineq_cons_vert( void * user_data , int id , int * nnz );

 /// nonzeros of the linking equality contribution of node id
 static int nnz_link_eq_cons( void * user_data , int id , int * nnz );

 /// nonzeros of the linking inequality contribution of node id
 static int nnz_link_ineq_cons( void * user_data , int id , int * nnz );

 /// identically zero count, used for the (absent) quadratic objective
 static int nnz_all_zero( void * user_data , int id , int * nnz );

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 // matrix callbacks; all matrices are in CSR row-major format

 /// diagonal equality matrix: A for the root, D_i for leaf i
 static int mat_eq_cons_diag( void * user_data , int id , int * krowM ,
			      int * jcolM , double * M );

 /// vertical equality matrix: root-variable part of the leaf rows
 static int mat_eq_cons_vert( void * user_data , int id , int * krowM ,
			      int * jcolM , double * M );

 /// diagonal inequality matrix
 static int mat_ineq_cons_diag( void * user_data , int id , int * krowM ,
				int * jcolM , double * M );

 /// vertical inequality matrix
 static int mat_ineq_cons_vert( void * user_data , int id , int * krowM ,
				int * jcolM , double * M );

 /// linking equality contribution of node id to the global linking rows
 static int mat_link_eq_cons( void * user_data , int id , int * krowM ,
			      int * jcolM , double * M );

 /// linking inequality contribution of node id to the global linking rows
 static int mat_link_ineq_cons( void * user_data , int id , int * krowM ,
				int * jcolM , double * M );

 /// identically empty matrix, used for the (absent) quadratic objective
 static int mat_all_zero( void * user_data , int id , int * krowM ,
			  int * jcolM , double * M );

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 // vector callbacks: objective coefficients, rhs / row bounds, variable
 // bounds, and their active flags

 /// linear objective coefficients of the variables of node id
 static int obj_vars( void * user_data , int id , double * vec , int len );

 /// rhs of the node-local equality rows
 static int rhs_eq_cons( void * user_data , int id , double * vec ,
			 int len );

 /// upper bounds of the node-local inequality rows
 static int rhs_ineq_cons( void * user_data , int id , double * vec ,
			   int len );

 /// lower bounds of the node-local inequality rows
 static int lhs_ineq_cons( void * user_data , int id , double * vec ,
			   int len );

 /// rhs of the global linking equality rows
 static int rhs_link_eq_cons( void * user_data , int id , double * vec ,
			      int len );

 /// upper bounds of the global linking inequality rows
 static int rhs_link_ineq_cons( void * user_data , int id , double * vec ,
				int len );

 /// lower bounds of the global linking inequality rows
 static int lhs_link_ineq_cons( void * user_data , int id , double * vec ,
				int len );

 /// active flags of the upper bounds of the node-local inequality rows
 static int flag_rhs_ineq_cons( void * user_data , int id , double * vec ,
				int len );

 /// active flags of the lower bounds of the node-local inequality rows
 static int flag_lhs_ineq_cons( void * user_data , int id , double * vec ,
				int len );

 /// active flags of the upper bounds of the linking inequality rows
 static int flag_rhs_link_ineq_cons( void * user_data , int id ,
				     double * vec , int len );

 /// active flags of the lower bounds of the linking inequality rows
 static int flag_lhs_link_ineq_cons( void * user_data , int id ,
				     double * vec , int len );

 /// upper bounds of the variables of node id
 static int ub_vars( void * user_data , int id , double * vec , int len );

 /// lower bounds of the variables of node id
 static int lb_vars( void * user_data , int id , double * vec , int len );

 /// active flags of the variable upper bounds
 static int flag_ub_vars( void * user_data , int id , double * vec ,
			  int len );

 /// active flags of the variable lower bounds
 static int flag_lb_vars( void * user_data , int id , double * vec ,
			  int len );

/** @} ---------------------------------------------------------------------*/
/*--------------------- PROTECTED FIELDS OF THE CLASS ----------------------*/
/*--------------------------------------------------------------------------*/

 /// the PIPS distributed problem tree owned by this solver
 pipsipmpp::DistributedInputTree * pips_tree;

 /// the PIPS solver interface owned by this solver
 pipsipmpp::PIPSIPMppInterface * pips_interface;

 /// number of PIPS nodes: 0 is the root, 1 .. n_nodes - 1 are the leaves
 Index n_nodes = 0;

 /// the Block subtrees grouped by PIPS node
 std::vector< std::vector< Block * > > nodes_subtrees;

 /// map from each Block to the PIPS leaf containing it
 std::unordered_map< const Block * , std::size_t > block_to_leaf;

 /// number of variables in each node
 std::vector< int > n_var_node;

 /// the variables assigned to each node
 std::vector< std::vector< const ColVariable * > > var_node;

 /// number of node-local equality constraints in each node
 std::vector< int > n_eq_cons_node;

 /// the node-local equality constraints
 std::vector< std::vector< const FRowConstraint * > > eq_cons_node;

 /// number of node-local inequality constraints in each node
 std::vector< int > n_ineq_cons_node;

 /// the node-local inequality constraints
 std::vector< std::vector< const FRowConstraint * > > ineq_cons_node;

 /// number of global linking equality constraints
 int n_link_eq_cons = 0;

 /// the global linking equality constraints
 std::vector< const FRowConstraint * > link_eq_cons;

 /// number of global linking inequality constraints
 int n_link_ineq_cons = 0;

 /// the global linking inequality constraints
 std::vector< const FRowConstraint * > link_ineq_cons;

 /// map from each variable to the PIPS node it is assigned to
 std::unordered_map< const ColVariable * , int > var_to_node;

/*--------------------------------------------------------------------------*/
/*---------------------- PRIVATE PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

 private:

/*--------------------------------------------------------------------------*/
/*-------------------- PRIVATE METHODS OF THE CLASS ------------------------*/
/*--------------------------------------------------------------------------*/

 /// clears all PIPS-owned objects and data derived from the current Block
 void reset_pips_data( void );

 /// converts PIPS-IPM++ termination statuses into Solver status codes
 static int decode_pips_status( pipsipmpp::TerminationStatus status );

/*--------------------------------------------------------------------------*/
 // definition of the private name() method required by the factory

 SMSpp_insert_in_factory_h;

/*--------------------------------------------------------------------------*/

 };  // end( class( PIPSMILPSolver ) )

/*--------------------------------------------------------------------------*/

}  // end( namespace SMSpp_di_unipi_it )

/*--------------------------------------------------------------------------*/

#endif  /* PIPSMILPSolver.h included */

/*--------------------------------------------------------------------------*/
/*----------------------- End File PIPSMILPSolver.h ------------------------*/
/*--------------------------------------------------------------------------*/
