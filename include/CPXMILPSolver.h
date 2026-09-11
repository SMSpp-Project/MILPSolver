/*--------------------------------------------------------------------------*/
/*--------------------------- File CPXMILPSolver.h -------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the CPXMILPSolver class.
 *
 * CPXMILPSolver implements a general purpose solver that is able to tackle a
 * MI-QCQP problem (the objective can be nonconvex but all quadratic
 * constraints must be convex) expressed by a Block using IBM CLPEX.
 *
 * \author Enrico Calandrini \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \author Antonio Frangioni \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \author Niccolo' Iardella \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \copyright &copy; by Enrico Calandrini, Antonio Frangioni,
 *                   Niccolo' Iardella
 */
/*--------------------------------------------------------------------------*/
/*----------------------------- DEFINITIONS --------------------------------*/
/*--------------------------------------------------------------------------*/

#ifndef __CPXMILPSOLVER_H
 #define __CPXMILPSOLVER_H
                      /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include <ilcplex/cplex.h>

#include <unordered_set>

#include "MILPSolver.h"

// Include the proper CPLEX parameter mapping
#include <boost/preprocessor/cat.hpp>
#include <boost/preprocessor/stringize.hpp>
#include BOOST_PP_STRINGIZE( BOOST_PP_CAT( BOOST_PP_CAT( CPX, CPX_VERSION ), _defs.h ) )

/*--------------------------------------------------------------------------*/
/*----------------------------- NAMESPACE ----------------------------------*/
/*--------------------------------------------------------------------------*/

/// namespace for the Structured Modeling System++ (SMS++)
namespace SMSpp_di_unipi_it {

 class LinearFunction;  // forward declaration of LinearFunction

/*--------------------------------------------------------------------------*/
/*----------------------- CLASS CPXMILPSolver ------------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------- GENERAL NOTES --------------------------------*/
/*--------------------------------------------------------------------------*/
/// class for solving MILP problems via CPLEX
/** The CPXMILPSolver class derives from MILPSolver and extends the
 * base class to solve MILP problems using CPLEX.
 *
 * The CPXMILPSolver can be registered to any kind of Block (assuming that
 * it contains a MILP formulation) and it uses the base class functionalities
 * to build a matricial representation of a MILP problem, then it solves it
 * using CPLEX through its Callable Library. Moreover, it implements the
 * interface that MILPSolver provides for processing modifications.
 *
 * The main logic is in compute(). This method copies the vectors that describe
 * the LP problem into a CPLEX environment, it processes the modifications and
 * then solves the problem.
 * get_var_solution() retrieves the values of the variables from CPLEX, saves
 * them into the Block variables and evaluates the objective function.
 *
 * Besides the configuration parameters already present in MILPSolver,
 * this class adds one int parameter that specifies if an exception must be
 * thrown if there is inconsistency when storing a reduced cost
 * (See set_par( idx_type, int )).
 * Moreover, the user can include in the configuration all the parameters
 * supported by CPXsetintparam(), CPXsetdblparam() and CPXsetstrparam()
 * (See the CPLEX Callable Library reference manual for all of them). */

class CPXMILPSolver : public MILPSolver {

/*--------------------------------------------------------------------------*/
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

 public:

/*--------------------------------------------------------------------------*/
/*---------------------------- PUBLIC TYPES --------------------------------*/
/*--------------------------------------------------------------------------*/
/** @name public types of CPXMILPSolver
 * @{ */

 /// enum for integer parameters
 enum int_par_type_CPXS {
  // note: intCutSepPar has moved to MILPSolver base (enum int_par_type_MILP)
  intFirstCPLEXPar = intLastAlgParMILP ,  ///< first CPLEX int/long parameter
  /// first allowed new int parameter for derived classes
  intLastAlgParCPXS = intFirstCPLEXPar + CPX_NUM_INT_PARS
  };

 /// enum for double parameters
 enum dbl_par_type_CPXS {
  /// first CPLEX double parameter
  dblFirstCPLEXPar = dblLastAlgParMILP,
  /// first allowed new double parameter for derived classes
  dblLastAlgParCPXS = dblFirstCPLEXPar + CPX_NUM_DBL_PARS
  };

 /// enum for string parameters
 enum str_par_type_CPXS {
  /// first CPLEX string parameter
  strFirstCPLEXPar = strLastAlgParMILP,
  /// first allowed new string parameter for derived classes
  strLastAlgParCPXS = strFirstCPLEXPar + CPX_NUM_STR_PARS
  };

 /// enum for vector-of-int parameters
 enum vint_par_type_CPXS {
  /// indices of separation Configurations in the "Configuration DB"
  vintCutSepCfgInd = vintLastAlgParMILP ,
  /// first allowed new vector-of-iint parameter for derived classes
  vintLastAlgParCPXS
  };

 /// enum for vector-of-string parameters
 enum cstr_par_type_CPXS {
 /// filenames to define the "Configuration DB"
  vstrConfigDBFName = vstrLastAlgParMILP ,
  /// first allowed new vector-of-string parameter for derived classes
  vstrLastAlgParCPXS
  };

/*--------------------------------------------------------------------------*/
 // "importing" a few types from Block

 using Subset = Block::Subset;
 
 using c_Subset = Block::c_Subset;

 using Range = Block::Range;

/** @} ---------------------------------------------------------------------*/
/*--------------------- CONSTRUCTOR AND DESTRUCTOR -------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Constructor and Destructor
 * @{ */

 CPXMILPSolver( void );

 ~CPXMILPSolver() override;

/** @} ---------------------------------------------------------------------*/
/*--------------------- DERIVED METHODS OF BASE CLASS ----------------------*/
/*--------------------------------------------------------------------------*/
/** @name Public Methods derived from base classes
 * @{ */

 /// sets the Block that the Solver has to solve and initializes CPLEX
 void set_Block( Block * block ) override;

 // note: the public compute() entry point is inherited from MILPSolver;
 // CPXMILPSolver implements only the CPLEX-specific solve in
 // guts_of_compute() below (protected)

 /// returns a valid lower bound on the optimal objective function value
 OFValue get_lb( void ) override;

 /// returns a valid upper bound on the optimal objective function value
 OFValue get_ub( void ) override;

 /// returns the value of the current solution, if any
 OFValue get_var_value( void ) override;

 /// tells whether a solution is available
 bool has_var_solution( void ) override;

 /// tells whether the current solution is feasible
 bool is_var_feasible( void ) override;

 /// writes the current solution in the Block
 void get_var_solution( Configuration * solc = nullptr ) override;

 /// writes a given solution vector in the Block
 /** Implementation of get_var_solution(9 taking the values of the solution
  * to be written in the Block out of a std::vector< double > at least as
  * long as there are columns (no checks performed). */
 void get_var_solution( const std::vector< double > & x );

 /// tells whether an unbounded direction is available
 bool has_var_direction( void ) override;

 /// writes the current unbounded direction in the Block
 void get_var_direction( Configuration * dirc = nullptr ) override;

 /// tells whether a dual solution is available
 bool has_dual_solution( void ) override;

 /// tells whether the current dual solution is feasible
 bool is_dual_feasible( void ) override;

 /// writes the current dual solution in the Block
 void get_dual_solution( Configuration * solc = nullptr ) override;

 /// tells whether a dual unbounded direction is available
 bool has_dual_direction( void ) override;

 /// writes the current dual unbounded direction in the Block
 void get_dual_direction( Configuration * dirc = nullptr ) override;

 /// writes the LP on the specified file
 void write_lp( const std::string & filename ) override;

 /// returns the number of nodes used to solve a MIP
 [[nodiscard]] int get_nodes( void ) const override;

 /// clears the CPLEX environment
 void clear_problem( unsigned int what ) override;

 /// loads the problem into CPLEX
 void load_problem( void ) override;

 /// returns the number of nodes used to solve a MIP
 [[nodiscard]] int get_explored_nodes( void ) const override;

 /// returns the estimated number of nodes left
 [[nodiscard]] long get_left_nodes( void ) const override;

 /// Returns a true value if a feasible solution is known, 
 //  false otherwise.
 [[nodiscard]] bool has_feasible_sol( void ) override;

 /// Returns elapsed solver runtime (in second).
 [[nodiscard]] double get_runtime( void ) const override;

 /// Returns a unique identifier for the node currently being explored  
 //  in the branch-and-bound algorithm for a MIP problem.  
 //  
 /// NOTE: This method should only be called during the callback process  
 //  and in specific situations (e.g., when a new incumbent solution is found,  
 //  and you need to identify the node from which it originates).  
 [[nodiscard]] long get_id_node( void ) const override;

/** 
 * Adds multiple MIP starts to a MIP problem. This function allows the solver 
 * to receive multiple sets of starting values by providing vectors of variable 
 * indices and corresponding values for each start.
 * 
 * NOTE: Partial solutions are allowed. In such cases, the solver will attempt 
 * to infer values for the unspecified variables.
 */
 void add_mip_starts( 
  std::vector< std::vector<int> > varidxs, 
  std::vector< std::vector<double> > varvalues ) override;

 #ifdef MILPSolver_DEBUG
  /// check the dictionaries for inconsistencies
  void check_status( void ) override;
 #endif
 
/** @} ---------------------------------------------------------------------*/
/*------------------- METHODS FOR HANDLING THE PARAMETERS ------------------*/
/*--------------------------------------------------------------------------*/
/** @name Methods for handling parameters
 * @{ */

 /// sets an integer parameter with the given value
 /** Set the "int" parameters specific of CPXMILPSolver, together with the
  * parameters of MILPSolver that CPXMILPSolver actually "listens to" and all
  * parameters supported by Cplex:
  *
  * - intCutSepPar [0]: coded bit-wise, indicate if and when separation of
  *                     either user cuts or lazy constraints is performed:
  *
  *   bit 0 : 1 (+1) if separation of user cuts is performed at the root
  *           node only
  *
  *   bit 1 : 1 (+2) if separation of user cuts is performed at every other
  *           node except the root one
  *
  *   bit 2 : 1 (+4) if separation of lazy constraints is performed each time
  *           a feasible solution is generated
  *
  *   bit 3-4: encode how user cuts are added to Cplex
  *            0 (+0) as CPX_USECUT_FILTER, i.e., "The cut is treated exactly
  *                   as cuts generated by CPLEX; that is, CPLEX applies its
  *                   filtering process and can possibly not even add the cut
  *                   to the relaxation, for example, if CPLEX deems other
  *                   cuts more effective, or if the cut is too dense."
  *            1 (+8) as CPX_USECUT_PURGE, i.e., "The cut is added to the
  *                   relaxation but can be purged later on if CPLEX deems
  *                   the cut ineffective."
  *            2 (+16) as CPX_USECUT_FORCE, i.e., "The cut is added to the
  *                    relaxation and stays there"
  *
  *   See vintCutSepCfgInd for properly setting Configurations for the
  *   corresponding calls to generate_dynamic_constraint(). */

 void set_par( idx_type par , int value ) override;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// sets a double parameter with the given value
 void set_par( idx_type par , double value ) override;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// sets a string parameter with the given value
 void set_par( idx_type par , std::string && value ) override;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// sets a vector-of-int parameter with the given value
 /** Set the vector-of-int parameters specific of CPXMILPSolver (note that
  * Cplex itself does not have any):
  *
  * - vintCutSepCfgInd [empty]: sets the Configuration for the various user
  *                             cuts / lazy constraints separations (see
  *   intCutSepPar) in terms of their indices in the "Configuration DataBase"
  *   (see vstrConfigDBFName). In particular:
  *
  *   = the 1st element sets the Configuration to be passed to
  *     generate_dynamic_constraint() when user cuts are to be separated at
  *     the root node
  *
  *   = the 2nd element sets the Configuration to be passed to
  *     generate_dynamic_constraint() when user cuts are to be separated at
  *     any other node except the root
  *
  *   = the 3rd element sets the Configuration to be passed to
  *     generate_dynamic_constraint() when user lazy constraints are to be
  *     separated for any feasible solution
  *
  *   If the passed vector is shorter than 3 elements, any missing ones are
  *   treated as "pass no Configuration" (nullptr). Similarly, if one entry
  *   is either negative or >= the size of the "Configuration DataBase", then
  *   "pass no Configuration" is assumed. */

 void set_par( idx_type par , std::vector< int > && value ) override;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// sets a vector-of-string parameter with the given value
 /** Set the vector-of-string parameters specific of CPXMILPSolver (note that
  * Cplex itself does not have any):
  *
  * - vstrConfigDBFName [empty]: provides file names used to construct the
  *                              "Configuration DataBase" that can be used
  *   to configure some operations on the underlying Block (e.g., user cuts
  *   or lazy constraints separation). Each entry in the vector is used as
  *   a filename out of which load a Configuration object that is then
  *   stored. This Configuration object is then "named" with the index that
  *   the filename has in this vector of string, so that it can be used for
  *   possibly multiple tasks. Note that it is assumed that using the
  *   Configuration objects does not change them. Note that the file names
  *   can actually be empty or "wrong", in which case nullptr is used. */

 void set_par( idx_type par , std::vector< std::string > && value ) override;

/*--------------------------------------------------------------------------*/
 /// returns the number of integer parameters
 [[nodiscard]] idx_type get_num_int_par( void ) const override;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// returns the number of double parameters
 [[nodiscard]] idx_type get_num_dbl_par( void ) const override;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// returns the number of string parameters
 [[nodiscard]] idx_type get_num_str_par( void ) const override;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// returns the number of vector-of-int parameters
 [[nodiscard]] idx_type get_num_vint_par( void ) const override;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// returns the number of vector-of-string parameters
 [[nodiscard]] idx_type get_num_vstr_par( void ) const override;

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
  * par as a CPLEX parameter. */
 [[nodiscard]] const std::string & get_dflt_str_par( idx_type par )
  const override;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// returns the default value of the specified vector-of-int parameter
 [[nodiscard]] const std::vector< int > & get_dflt_vint_par( idx_type par )
  const override;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// returns the default value of the specified vector-of-string parameter
 [[nodiscard]] const std::vector< std::string > & get_dflt_vstr_par(
					       idx_type par ) const override;

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
  * par as a CPLEX parameter. */
 [[nodiscard]] const std::string & get_str_par( idx_type par ) const override;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// returns the value of the specified vector-of-int parameter
 [[nodiscard]] const std::vector< int > & get_vint_par( idx_type par )
  const override;
 
/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// returns the value of the specified vector-of-string parameter
 [[nodiscard]] const std::vector< std::string > & get_vstr_par( idx_type par )
  const override;
 
/*--------------------------------------------------------------------------*/
 /// returns the index of the int parameter with the specified name
 [[nodiscard]] idx_type int_par_str2idx( const std::string & name )
  const override;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /** returns the name of the int parameter with the specified index
  * @note
  * Due to a limit in the implementation, the string referenced by
  * the return value is *overwritten* each time the method is called with
  * par as a CPLEX parameter. */
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
  * par as a CPLEX parameter. */
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
  * par as a CPLEX parameter. */
 [[nodiscard]] const std::string & str_par_idx2str( idx_type idx )
  const override;

/*--------------------------------------------------------------------------*/
 /// returns the index of the vector-of-int parameter with the specified name
 [[nodiscard]] idx_type vint_par_str2idx( const std::string & name )
  const override;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// returns the name of the vector-of-int parameter with the specified index
 [[nodiscard]] const std::string & vint_par_idx2str( idx_type idx )
  const override;

/*--------------------------------------------------------------------------*/
 /// returns the index of the vector-of-string parameter with the given name
 [[nodiscard]] idx_type vstr_par_str2idx( const std::string & name )
  const override;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// returns the name of the vector-of-string parameter with the given index
 [[nodiscard]] const std::string & vstr_par_idx2str( idx_type idx )
  const override;

/*--------------------------------------------------------------------------*/

 double up_cut_off( void ) const { return( UpCutOff ); }

 double lw_cut_off( void ) const { return( LwCutOff ); }

/*--------------------------------------------------------------------------*/
 /// callback implemented as a method of the class
 /** The implementation of CPLEX "generic" callback, which is used to check
  * for having reached prescribed upper/lower bounds and for user cuts / lazy
  * constraint separation, just calls this method.
  *
  * IMPORTANT NOTE: CPLEX has a different stance than SMS++ on dynamic
  *                 Constraint, in the sense that those that are added inside
  * a callback are not permanently added to the formulation and may be
  * discarded whole. In contrast, for SMS++ dynamic Constraint are
  * first-class citizens of the formulation. To reconcile this two different
  * viewpoints,
  *
  *     THE Modification ADDING DYNAMIC Constraint ARE *NOT* REMOVED FROM
  *     THE QUEUE OF ACTIVE Modification
  *
  * As a result, when Cplex terminates and gets re-solved (if ever), the
  * dynamic Constraint will be properly added to the formulation. This is
  * consistent with the view that Modification happening when the Solver is
  * running must not *necessarily* be immediately acted upon by changing the
  * model that the Solver is solving. */

 int callback( CPXCALLBACKCONTEXTptr context , CPXLONG contextid );
 
/** @} ---------------------------------------------------------------------*/
/*--------------------- PROTECTED PART OF THE CLASS ------------------------*/
/*--------------------------------------------------------------------------*/

 protected:

/*--------------------------------------------------------------------------*/
/*-------------------- PROTECTED METHODS OF THE CLASS ----------------------*/
/*--------------------------------------------------------------------------*/

 /// CPLEX back-end solve, called by MILPSolver::compute()
 /** Performs the actual CPLEX optimisation. Locking, Modification
  * processing and the LP cut-separation loop (intRelaxIntVars == 2) are
  * all handled by MILPSolver::compute(). */

 int guts_of_compute( void ) override;

/*--------------------------------------------------------------------------*/

 /** @name Get variable bounds for the problem
  *
  * The following two methods retrieve the upper and lower bound for the
  * given variable considering both the Variable bounds and all the active
  * OneVarConstraints active for that Variable. */

 /// gets the LB for the given variable in the problem
 double get_problem_lb( const ColVariable & var ) const override;

 /// gets the UB for the given variable in the problem
 double get_problem_ub( const ColVariable & var ) const override;

 /// gets both bounds for the given variable in the problem
 std::array< double , 2 > get_problem_bounds( const ColVariable & var )
  const override;

/** @} ---------------------------------------------------------------------*/
/*-------------------- METHODS FOR MODIFYING THE PROBLEM -------------------*/
/*--------------------------------------------------------------------------*/
/** @name Methods for modifying the constructed CPLEX problem
 *
 *  These methods implement the "Modification interface" of MILPSolver, so
 *  that all Modification are applied to the CPLEX formulation.
 *  @{ */

 /// handles a Variable Modification
 void var_modification( const VariableMod * mod ) override;

 /// handles an Objective Modification
 void objective_modification( const ObjectiveMod * mod ) override;

 /// handles a Constraint Modification
 void const_modification( const ConstraintMod * mod ) override;

 /// handles a bound (OneVarConstraint) Modification
 void bound_modification( const OneVarConstraintMod * mod ) override;

 /// handles a Function Modification applied to the Objective
 void objective_function_modification( const FunctionMod * mod ) override;

 /// handles a Function Modification applied to a Constraint
 void constraint_function_modification( const FunctionMod * mod ) override;

 /// handles a Function Variable Modification applied to the Objective
 void objective_fvars_modification( const FunctionModVars * mod )
  override;

 /// handles a Function Variable Modification applied to a Constraint
 void constraint_fvars_modification( const FunctionModVars * mod )
  override;

 // handles a dynamic Modification
 // no point in defining it, just calls the base class method
 // void dynamic_modification( const BlockModAD * mod ) override;

 /// adds a single new dynamic FRowConstraint
 void add_dynamic_constraint( const FRowConstraint * con ) override;

 /// adds a single new dynamic bound (OneVarConstraint)
 void add_dynamic_bound( const OneVarConstraint * con ) override;

 /// adds a single new dynamic ColVariable
 void add_dynamic_variable( const ColVariable * var ) override;

 /// removes a single dynamic FRowConstraint
 void remove_dynamic_constraint( const FRowConstraint * con ) override;

 /// removes a single dynamic ColVariable
 void remove_dynamic_variable( const ColVariable * var ) override;

 /// removes a single dynamic bound (OneVarConstraint)
 void remove_dynamic_bound( const OneVarConstraint * con ) override;

/** @} ---------------------------------------------------------------------*/
/*----------------------- METHODS FOR CUT SEPARATION -----------------------*/
/*--------------------------------------------------------------------------*/
/** @name Methods for separating user cuts / lazy constraints
 *  @{ */

 /** From within the callback, run the cut separation invoking
  * generate_dynamic_constraints() with the given Configuration and then
  * examining the list of Modification to see if some dynamic Constraint have
  * been added; if so they are reported back under the form needed to be
  * added as user cuts or lazy constraints (which is the same).
  *
  * Note that all vectors are supposed to be empty at the beginning of the
  * call, and they will still be empty if no cuts are found. */

 void perform_separation( Configuration * cfg ,
			  std::vector< int > & rmatbeg ,
			  std::vector< int > & rmatind ,
			  std::vector< double > & rmatval ,
			  std::vector< double > & rhs , 
			  std::vector< char > & sense );

/** @} ---------------------------------------------------------------------*/
 /// maps a Solver integer parameter into a Cplex one
 /** Maps the Solver integer parameter \p par into a Cplex one;
  * returns a positive number of it is an int parameter and a negative
  * number if it is a long one. Returns 0 if not a Cplex parameter. */
 int cpx_int_par_map( idx_type par ) const;
 
/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// maps a Solver double parameter into a Cplex one (or 0)
 int cpx_dbl_par_map( idx_type par ) const;

/*--------------------------------------------------------------------------*/
/*-------------------- PROTECTED FIELDS OF THE CLASS -----------------------*/
/*--------------------------------------------------------------------------*/

 CPXENVptr env; ///< CPLEX environment
 CPXLPptr lp;   ///< CPLEX LP problem

 bool f_callback_set;  // true if the callback has been set

 // note: CutSepPar (the bitwise-encoded cut/lazy separation parameter)
 // is now an inherited member of MILPSolver base

 /** vector containing the indices of the Configuration for the various
  * user cuts / lazy constraints separations in the "Configuration DB" */
 std::vector< int > CutSepCfgInd;

 /** vector containing the filenames used to load of the Configuration of
  * the "Configuration DB" */
 std::vector< std::string > ConfigDBFName;

 /// the "Configuration DB" istself
 std::vector< Configuration * > v_ConfigDB;

 /** pointer used to keep track of the current context of the callback */
 CPXCALLBACKCONTEXTptr current_Cntx;

 /** an integer value specifying the context in which the callback is invoked. */
 CPXLONG current_Cntx_id;

 /** double storing the timestamp when optimization begins. */
 double starting_time;

 /// the mutex to ensure that CPLEX threads do not overstep in the callback
 /** Since CPLEX is multi-threaded, lock()-ing the Block with the f_id of
  * CPXMILPSolver is not enough to prevent concurrent access to it. This is
  * an issue in che callback(), in particular when user cuts / lazy
  * constraints separation is required, and therefore 1) a solution has to
  * be written in the Variable, 2) generate_dynamic_constraints() has to be
  * called, which may cause the addition of new dynamic Constraint to the
  * Block. Thus, CPXMILPSolver will use this mutex to ensure mutual exclusion
  * of the CPLEX threads for the critical sections of the callback(). */
 std::mutex f_callback_mutex;

 /* In CPXMILPSolver we handle quadratic constraints like 
  * q x + x^T Q x <= q_0 by considering different scenarios:
  *
  *  - if q is null, then we simply add the constraint x^T Q x <= q_0
  *
  *  - otherwise, we build two separate constraint: q x + v <= q_0 
  *    and v >= x^T Q x, with v being an auxiliary variable. This is 
  *    because CPLEX does not allow to directly modify quadratic 
  *    constraints. Thus, we will need to store for each quadratic constraint 
  *    the CPLEX index of relative auxiliary variable and constraint being 
  *    built. To achieve this goal we will use two auxiliary vectors 
  *    cpx_quad_var_aux and cpx_quad_con_aux, with length equal to the 
  *    number of rows and value -1 for linear constraint. In the vector
  *    cpx_idx_aux_qvar we will simply keep track of the indices of 
  *    auxiliary variables built for this pourpose.
  *
  * NOTE: The set of indices of quadratic and linear rows are disjoint. 
  * For this reason, if the n-th constraint is quadratic, we will store 
  * in cpx_quad_con_aux[n] the index of the quadratic constraint in the
  * relative set. */
  std::vector< int > cpx_quad_var_aux;
  std::vector< int > cpx_quad_con_aux;
  std::vector< int > cpx_idx_aux_qvar; // Need to be sorted

 // function to retrieve actual idx of variable considering auxiliary ones
 int cpx_index_of_variable( const ColVariable * var ) const;

 // function to retrieve actual idx of dynamic variable considering auxiliary ones
 int cpx_index_of_dynamic_variable( const ColVariable * var ) const;

 // function to retrieve actual idx of constraint. In CPLEX indices of linear and
 // quadratic constraint are disjoint, so we need to retrieve the actual index 
 // based on the type of constraint.
 int cpx_index_of_linear_constraint( const FRowConstraint * con ) const;
 
 /** @name Handling of CPLEX parameters
  *
  * The following maps are used to keep a relationship between SMS++ parameter
  * system and CPLEX parameters. This allows us to use CPLEX parameters
  * (See CPLEX Parameters Reference Manual from IBM) as they were SMS++
  * parameters with the same names, for example in configuration files.
  *
  * Note: since SMS++ does not support long parameters, both int and
  *       long CPLEX parameters are handled as SMS++ int parameters.
  * @{ */

 const static std::array< int , CPX_NUM_INT_PARS > SMSpp_to_CPLEX_int_pars;
 const static std::array< int , CPX_NUM_DBL_PARS > SMSpp_to_CPLEX_dbl_pars;
 const static std::array< int , CPX_NUM_STR_PARS > SMSpp_to_CPLEX_str_pars;

 const static std::array< std::pair< int , int > , CPX_NUM_INT_PARS >
  CPLEX_to_SMSpp_int_pars;
 const static std::array< std::pair< int , int > , CPX_NUM_DBL_PARS >
  CPLEX_to_SMSpp_dbl_pars;
 const static std::array< std::pair< int , int > , CPX_NUM_STR_PARS >
  CPLEX_to_SMSpp_str_pars;

/** @} ---------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

 double UpCutOff;  ///< externally set upper cutoff to terminate
 double LwCutOff;  ///< externally set lower cutoff to terminate
 
/*--------------------------------------------------------------------------*/
/*---------------------- PRIVATE PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

 private:

 /** Returns the SMS++ status corresponding to the given
  * CPLEX status returned by CPXgetstat() in case of a LP/QP,
  * or by CPXgetsubstat() in case of a subproblem of a MIP. */
 static int decode_lqp_status( int status );

 /** Returns the SMS++ status corresponding to the given
  * CPLEX status returned by CPXgetstat() in case of a MIP. */
 static int decode_mip_status( int status );

 /** Returns the SMS++ status corresponding to the given
  * CPLEX error returned by CPXmipopt(), CPXlpopt() or CPXqpopt(). */
 static int decode_cpx_error( int error );

 /** Reloads a constraint.
  * To be used as fallback method for constraint FunctionMods. */
 // void reload_constraint( const LinearFunction * lf );

 /** Reloads the objective.
  * To be used as fallback method for objective FunctionMods. */
 // void reload_objective( Function * f );

 /// update problem type: false for a linear one, true for a quadratic one
 void update_problem_type( bool quad );

 // get the right Configuration for ci = 0, 1, 2
 Configuration * get_cfg( Index ci ) const;

 /** Create the structures used to provide the quadratic objective matrix 
  * to CPLEX with the function CPXcopyquad(). */
 void generate_qobj_matrix( std::vector< int > & qmatbeg ,
			  std::vector< int > & qmatcnt ,
			  std::vector< int > & qmatind ,
			  std::vector< double > & qmatval );

 /** Create the structures used to provide the quadratic matrix for the
 * constraint of index row to CPLEX. */
 void generate_qcon_matrix( std::vector< int > & qidx1 ,
			  std::vector< int > & qidx2 ,
			  std::vector< double > & qcoeff ,
        Index row ,
        bool lin_null );

 /** Evaluate the gradient of a specific quadratic constraint
  * in the optimum find by CPLEX. */
 double evaluate_dual_qcon( Index row ,
        std::vector< double > x_sol );

 /// FRowConstraints with \f$ lhs > rhs \f$ in their current state
 /** CPLEX encodes a ranged row with \f$ \mathit{rngval} =
  *  \mathit{rhs} - \mathit{lhs} \f$. When \f$ \mathit{rngval} < 0 \f$
  *  (i.e. \f$ \mathit{lhs} > \mathit{rhs} \f$), CPLEX silently re-interprets
  *  the row as the swapped interval
  *  \f$ [\mathit{rhs}, \mathit{lhs}] \f$, masking the intended
  *  infeasibility. To preserve the SMS++ semantics we detect this case in
  *  const_modification() / add_dynamic_constraint(), encode the row in
  *  CPLEX as a feasible equality at \f$ \mathit{rhs} \f$, and record the
  *  affected FRowConstraint here. compute() short-circuits to
  *  kInfeasible if this set is non-empty. */
 std::unordered_set< const FRowConstraint * > f_inverted_rows;

/*--------------------------------------------------------------------------*/

 SMSpp_insert_in_factory_h;

/*--------------------------------------------------------------------------*/

 };  // end( class( CPXMILPSolver ) )

/*--------------------------------------------------------------------------*/

}  // end( namespace SMSpp_di_unipi_it )

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

#endif  /* CPXMILPSolver.h included */

/*--------------------------------------------------------------------------*/
/*------------------------ End File CPXMILPSolver.h ------------------------*/
/*--------------------------------------------------------------------------*/
