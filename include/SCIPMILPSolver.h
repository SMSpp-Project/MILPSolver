/*--------------------------------------------------------------------------*/
/*--------------------------- File SCIPMILPSolver.h ------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the SCIPMILPSolver class.
 *
 * SCIPMILPSolver derives from MILPSolver and it uses the facilities
 * provided by the base class to implements a general purpose MILP solver
 * using calls to the ZIB SCIP API.
 *
 * \author Antonio Frangioni \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \author Niccolo' Iardella \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 * 
 * \author Enrico Calandrini \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \copyright &copy; by Antonio Frangioni, Niccolo' Iardella
 */
/*--------------------------------------------------------------------------*/
/*----------------------------- DEFINITIONS --------------------------------*/
/*--------------------------------------------------------------------------*/

#ifndef __SCIPMILPSOLVER_H
 #define __SCIPMILPSOLVER_H
                      /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "MILPSolver.h"

#include <scip/scip.h>
#include <objscip/objscip.h>

// Include the proper SCIP parameter mapping
#include <boost/preprocessor/cat.hpp>
#include <boost/preprocessor/stringize.hpp>
#include BOOST_PP_STRINGIZE( BOOST_PP_CAT( BOOST_PP_CAT( SCIP, SCIP_VERSION ), _defs.h ) )

/*--------------------------------------------------------------------------*/
/*----------------------------- NAMESPACE ----------------------------------*/
/*--------------------------------------------------------------------------*/

/// namespace for the Structured Modeling System++ (SMS++)
namespace SMSpp_di_unipi_it
{
/*--------------------------------------------------------------------------*/
/*----------------------- CLASS SCIPMILPSolver -----------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------- GENERAL NOTES --------------------------------*/
/*--------------------------------------------------------------------------*/
/// class for solving MILP problems via SCIP.
/** The SCIPMILPSolver class derives from MILPSolver and extends the
 * base class to solve MILP problems using SCIP.
 *
 * The SCIPMILPSolver can be registered to any kind of Block (assuming that
 * it contains a MILP formulation) and it uses the base class functionalities
 * to build a matricial representation of a MILP problem, then it solves it
 * using SCIP. Moreover, it implements the interface that MILPSolver provides
 * for processing modifications.
 *
 * The main logic is in compute(). This method copies the vectors that
 * decribe the MILP problem into a SCIP environment, it processes the
 * modifications and then solves the problem. get_var_solution() retrieves
 * the values of the variables from SCIP, saves them into the Block
 * variables and evaluates the objective function.
 *
 * Besides the configuration parameters already present in MILPSolver,
 * the user can include in the configuration all the parameters
 * supported by SCIP (See https://www.scipopt.org/doc/html/PARAMETERS.php).
 */

class SCIPMILPSolver : public MILPSolver
{
/*--------------------------------------------------------------------------*/
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

 public:

/*--------------------------------------------------------------------------*/
/*---------------------------- PUBLIC TYPES --------------------------------*/
/*--------------------------------------------------------------------------*/

 /// Types of integer parameters
 enum int_par_type_SCPS {
  /// throws exception if there is inconsistency when storing a reduced cost
  intThrowReducedCostException = intLastAlgParMILP ,
  intCutSepPar ,  ///< parameter for deciding if/when cut separation is done
  /// First SCIP int/long parameter
  intFirstSCIPPar ,
  /// First allowed new int parameter for derived classes
  intLastAlgParSCPS = intFirstSCIPPar + SCIP_NUM_INT_PARS
  };

 /// Types of double parameters
 enum dbl_par_type_SCPS {
  /// First SCIP double parameter
  dblFirstSCIPPar = dblLastAlgParMILP ,
  /// First allowed new double parameter for derived classes
  dblLastAlgParSCPS = dblFirstSCIPPar + SCIP_NUM_DBL_PARS
  };

 /// Types of string parameters
 enum str_par_type_SCPS {
  /// First SCIP string parameter
  strFirstSCIPPar = strLastAlgParMILP ,
  /// First allowed new string parameter for derived classes
  strLastAlgParSCPS = strFirstSCIPPar + SCIP_NUM_STR_PARS
  };

 /// enum for vector-of-int parameters
 enum vint_par_type_SCPS {
  /// indices of separation Configurations in the "Configuration DB"
  vintCutSepCfgInd = vintLastAlgParMILP ,
  /// first allowed new vector-of-iint parameter for derived classes
  vintLastAlgParCPXS
  };

 /// enum for vector-of-string parameters
 enum cstr_par_type_SCPS {
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

/*--------------------------------------------------------------------------*/
/*--------------------- CONSTRUCTOR AND DESTRUCTOR -------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Constructor and Destructor
 *  @{ */

 SCIPMILPSolver();

 ~SCIPMILPSolver() override;

/** @} ---------------------------------------------------------------------*/
/*--------------------- DERIVED METHODS OF BASE CLASS ----------------------*/
/*--------------------------------------------------------------------------*/
/** @name Public Methods derived from base classes
 *  @{ */

 /// sets the Block that the Solver has to solve and initializes CPLEX.
 void set_Block( Block * block ) override;

 /// optimizes the problem with SCIP
 int compute( bool changedvars = false ) override;

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

 /// clears the SCIP environment and the matrix representation
 void clear_problem( unsigned int what ) override;

 /// loads the problem into SCIP
 void load_problem( void ) override;

 // get the right Configuration for ci = 0, 1, 2
 Configuration * get_cfg( Index ci ) const;

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
			  std::vector< double > & lhs );

  /* From within the class SCIPMILPSolver_Conhdlr it is not possible to set
   * some protected field of the class useful to avoid collision between threads 
   * when performing separation. Thus, the two following public functions allows 
   * us to obtain this results from external class. */
  
  void set_f_cb_mutex( void );

  void unset_f_cb_mutex( void );

  /// get the actual SCIP var used
  std::vector< SCIP_VAR * > get_SCIP_var( void );

  #ifdef MILPSOLVER_DEBUG
  /// check the dictionaries for inconsistencies
   void check_status( void ) override;
  #endif

/** @} ---------------------------------------------------------------------*/
/*------------------- METHODS FOR HANDLING THE PARAMETERS ------------------*/
/*--------------------------------------------------------------------------*/
/** @name Methods for handling parameters
 *  @{ */

 /// sets an integer parameter with the given value
 void set_par( idx_type par , int value ) override;

 /// sets a double parameter with the given value
 void set_par( idx_type par , double value ) override;

 /// sets a string parameter with the given value
 void set_par( idx_type par , std::string && value ) override;

 /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// sets a vector-of-int parameter with the given value
 /** Set the vector-of-int parameters specific of SCIPMILPSolver (note that
  * SCIP itself does not have any):
  * 
  * - intThrowReducedCostException [0]: it indicates whether an exception must
  *                                     be thrown if there is an inconsistency
  *   when a reduced cost is being stored during a call to get_dual_solution()
  *   or get_dual_direction(). The reduced cost of a Variable is stored in at
  *   most one OneVarConstraint on that Variable. It may happen that a
  *   Variable has no OneVarConstraint, in which case its reduced cost will
  *   not be stored and will be lost. Usually, the reduced cost of a Variable
  *   is of interest if the Variable has a finite nonzero lower or upper
  *   bound. In this case, if a OneVariableConstraint for that Variable is not
  *   found, an exception is thrown. More specifically, there are two cases in
  *   which an exception is thrown:
  *
  *   1) The Variable is fixed to a finite nonzero value and there is no
  *      OneVarConstraint on that Variable whose lower and upper bounds are
  *      both equal to the value of that Variable.
  *
  *   2) The Variable is not fixed, it has a finite nonzero lower or upper
  *      bound and there is no OneVarConstraint on that Variable whose lower
  *      or upper bound match the bounds of the Variable.
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
 /** Set the vector-of-string parameters specific of SCIPMILPSolver (note that
  * SCIP itself does not have any):
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

 /// gets the number of integer parameters
 [[nodiscard]] idx_type get_num_int_par( void ) const override;

 /// gets the number of double parameters
 [[nodiscard]] idx_type get_num_dbl_par( void ) const override;

 /// gets the number of string parameters
 [[nodiscard]] idx_type get_num_str_par( void ) const override;

 /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// returns the number of vector-of-int parameters
 [[nodiscard]] idx_type get_num_vint_par( void ) const override;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// returns the number of vector-of-string parameters
 [[nodiscard]] idx_type get_num_vstr_par( void ) const override;

 /// gets the default value of the specified integer parameter
 [[nodiscard]] int get_dflt_int_par( idx_type par ) const override;

 /// gets the default value of the specified double parameter
 [[nodiscard]] double get_dflt_dbl_par( idx_type par ) const override;

 /** Gets the default value of the specified string parameter
  * @note
  * Due to a limit in the implementation, the string referenced by
  * the return value is *overwritten* each time the method is called with
  * par as a SCIP parameter. */

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

 /// gets the value of the specified integer parameter
 [[nodiscard]] int get_int_par( idx_type par ) const override;

 /// gets the value of the specified double parameter
 [[nodiscard]] double get_dbl_par( idx_type par ) const override;

 /** Gets the value of the specified string parameter
  * @note
  * Due to a limit in the implementation, the string referenced by
  * the return value is *overwritten* each time the method is called with
  * par as a SCIP parameter. */
 
 [[nodiscard]] const std::string & get_str_par( idx_type par ) const override;

 /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// returns the value of the specified vector-of-int parameter
 [[nodiscard]] const std::vector< int > & get_vint_par( idx_type par )
  const override;
 
/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// returns the value of the specified vector-of-string parameter
 [[nodiscard]] const std::vector< std::string > & get_vstr_par( idx_type par )
  const override;

 /// returns the index of the int parameter with the specified name
 [[nodiscard]] idx_type int_par_str2idx( const std::string & name )
  const override;

 /// returns the name of the int parameter with the specified index
 [[nodiscard]] const std::string & int_par_idx2str( idx_type idx )
  const override;

 /// returns the index of the double parameter with the specified name
 [[nodiscard]] idx_type dbl_par_str2idx( const std::string & name )
  const override;

 /// returns the name of the double parameter with the specified index
 [[nodiscard]] const std::string & dbl_par_idx2str( idx_type idx )
  const override;

 /// returns the index of the string parameter with the specified name
 [[nodiscard]] idx_type str_par_str2idx( const std::string & name )
  const override;

 /// returns the name of the string parameter with the specified index
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

/** @} ---------------------------------------------------------------------*/
/*--------------------- PROTECTED PART OF THE CLASS ------------------------*/
/*--------------------------------------------------------------------------*/

 protected:

/*--------------------------------------------------------------------------*/
/*-------------------- PROTECTED FIELDS OF THE CLASS -----------------------*/
/*--------------------------------------------------------------------------*/

 bool f_callback_set;  // true if the callback has been set

  /** This variable indicates whether an exception must be thrown if there is
  * an inconsistency when a reduced cost is being stored during a call to
  * get_dual_solution() or get_dual_direction(). */
 bool throw_reduced_cost_exception;

 /** bitwise-encoded parameter for deciding if and when separation of user
  * cuts and lazy constraints is performed */
 unsigned char CutSepPar;

 /** vector containing the indices of the Configuration for the various
  * user cuts / lazy constraints separations in the "Configuration DB" */
 std::vector< int > CutSepCfgInd;

 /** vector containing the filenames used to load of the Configuration of
  * the "Configuration DB" */
 std::vector< std::string > ConfigDBFName;

 /// the "Configuration DB" istself
 std::vector< Configuration * > v_ConfigDB;

 /// the mutex to ensure that SCIP threads do not overstep in the callback
 std::mutex f_callback_mutex;

 /// SCIP environment
 SCIP * scip{};

 std::vector< SCIP_VAR * > vars;   ///< SCIP variables
 std::vector< SCIP_CONS * > cons;  ///< SCIP constraints

 /// SCIP auxiliary variables for QPs
 std::vector< SCIP_VAR * > aux_vars;
 /// SCIP auxiliary constraints for QPs
 std::vector< SCIP_CONS * > aux_cons;

 double UpCutOff;  ///< externally set upper cutoff to terminate
 double LwCutOff;  ///< externally set lower cutoff to terminate

/*--------------------------------------------------------------------------*/
/*------------------- PROTECTED METHODS OF THE CLASS -----------------------*/
/*--------------------------------------------------------------------------*/
 /** @name Get variable bounds for the problem
  *
  * The following two methods retrieve the upper and lower bound for the
  * given variable considering both the Variable bounds and all the active
  * OneVarConstraints active for that Variable.
  * @{ */

 /// gets the LB fot the given variable in the problem
 double get_problem_lb( const ColVariable & var ) const override;

 /// gets the UB fot the given variable in the problem
 double get_problem_ub( const ColVariable & var ) const override;

 /// gets both bounds for the given variable in the problem
 std::array< double , 2 > get_problem_bounds( const ColVariable & var )
  const override;

/** @} ---------------------------------------------------------------------*/
/*-------------------- METHODS FOR MODIFYING THE PROBLEM -------------------*/
/*--------------------------------------------------------------------------*/
/** @name Methods for modifying the constructed SCIP problem
 *  @{ */

 /// handles a variable modification
 void var_modification( const VariableMod * mod ) override;

 /// handles an objective modification
 void objective_modification( const ObjectiveMod * mod ) override;

 /// handles a constraint modification
 void const_modification( const ConstraintMod * mod ) override;

 /// handles a bound modification
 void bound_modification( const OneVarConstraintMod * mod ) override;

 /// handles a function modification applied to the objective
 void objective_function_modification( const FunctionMod * mod ) override;

 /// handles a function modification applied to a constraint
 void constraint_function_modification( const FunctionMod * mod ) override;

 /// handles a function vars modification to the objective
 void objective_fvars_modification( const FunctionModVars * mod ) override;

 /// handles a function vars modification to a constraint
 void constraint_fvars_modification( const FunctionModVars * mod ) override;

 /// handles a dynamic modification
 // no point in defining it, just calls the base class method
 // void dynamic_modification( const BlockModAD * mod ) override;

 /// adds a single new dynamic constraint
 void add_dynamic_constraint( const FRowConstraint * con ) override;

 /// adds a single new dynamic bound
 void add_dynamic_bound( const OneVarConstraint * con ) override;

 /// adds a single new dynamic variable
 void add_dynamic_variable( const ColVariable * var ) override;

 /// removes a single dynamic constraint
 void remove_dynamic_constraint( const FRowConstraint * con ) override;

 /// removes a single dynamic variable
 void remove_dynamic_variable( const ColVariable * var ) override;

 /// removes a single dynamic bound
 void remove_dynamic_bound( const OneVarConstraint * con ) override;

/** @} ---------------------------------------------------------------------*/
/*--------------------- PRIVATE FIELDS OF THE CLASS ------------------------*/
/*--------------------------------------------------------------------------*/

 private:
 
/*--------------------------------------------------------------------------*/
 /** @name Handling of SCIP parameters
 *
 * The following maps are used to keep a relationship between SMS++ parameter
 * system and SCIP parameters. This allows us to use SCIP parameters
 * (See https://www.scipopt.org/doc/html/PARAMETERS.php) as they were SMS++
 * parameters with the same names, for example in configuration files.
 *
 * Bool, int and long SCIP parameters are handled as SMS++ int parameters.
 * Real SCIP parameters are handled as SMS++ double parameters.
 * Char and string SCIP parameters are handled as SMS++ string parameters.
 * @{ */

 const static std::array< std::string , SCIP_NUM_INT_PARS >
  SMSpp_to_SCIP_int_pars;

 const static std::array< std::string , SCIP_NUM_DBL_PARS >
  SMSpp_to_SCIP_dbl_pars;
 
 const static std::array< std::string , SCIP_NUM_STR_PARS >
  SMSpp_to_SCIP_str_pars;

 const static std::array< std::pair< std::string , int > , SCIP_NUM_INT_PARS >
  SCIP_to_SMSpp_int_pars;

 const static std::array< std::pair< std::string , int > , SCIP_NUM_DBL_PARS >
  SCIP_to_SMSpp_dbl_pars;

 const static std::array< std::pair< std::string , int > , SCIP_NUM_STR_PARS >
  SCIP_to_SMSpp_str_pars;

/** @} ---------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

 SMSpp_insert_in_factory_h;

/*--------------------------------------------------------------------------*/

 };  // end( class SCIPMILPSolver )

/*--------------------------------------------------------------------------*/
/*--------------------- Class SCIPMILPSolver_Conhdlr -----------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*----------------------------- DEFINITIONS --------------------------------*/
/*--------------------------------------------------------------------------*/

#ifndef __SCIPMILPSOLVER_CONHDLR_H
 #define __SCIPMILPSOLVER_CONHDLR_H
 
/*--------------------------- GENERAL NOTES --------------------------------*/
/*--------------------------------------------------------------------------*/
/*
 * SCIPMILPSolver_Conhdlr is a tool developed for SCIPMILPSolver, in order 
 * to generate valid inequalities or even facets of the polyhedron described 
 * by a single constraint or a subset of the constraints of a single 
 * constraint class. It is essentially used to for user cuts / lazy

* The SCIPMILPSolver_Conhdlr class derives from scip::ObjConshdlr and it
  * creates specific user cut/ or lazy constraint to be added within a
  * SMS++ model handled by SCIPMILPSolver.  */

class SCIPMILPSolver_Conhdlr : public scip::ObjConshdlr
{
/*--------------------------------------------------------------------------*/
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

 public:

/*--------------------------------------------------------------------------*/
/*--------------------- CONSTRUCTOR AND DESTRUCTOR -------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Constructor and Destructor
 *  @{ */

 SCIPMILPSolver_Conhdlr( SCIP* scip, /**< SCIP data structure */
    SMSpp_di_unipi_it::SCIPMILPSolver* scipmilpsolver, /**< "parent" 
                                    * SCIPMILPSolver from which the Constraint 
                                    * handler has been called. */
    unsigned char SeparationPar    /**< separation decider parameter */
    );

 ~SCIPMILPSolver_Conhdlr() override;

/** @} ---------------------------------------------------------------------*/
/*--------------------- FUNDAMENTAL CALLBACK METHODS -----------------------*/
/*--------------------------------------------------------------------------*/

/** constraint enforcing method of constraint handler for LP solutions
*
*  The method is called at the end of the node processing loop for a node 
*  where the LP was solved. The LP solution has to be checked for 
*  feasibility.
*
*  In this function we add new lazy constraints (lc) in the model.
*  NOTE: we are sure that a new lc can be added, because first a SCIP_CHECK
*  method has declared that separation is possible. 
*
*  Possible return values for *result:
*  - SCIP_CUTOFF     : the node is infeasible in the variable's bounds and 
*                      can be cut off
*  - SCIP_SEPARATED  : a cutting plane was generated
*  - SCIP_INFEASIBLE : at least one constraint is infeasible, but it was not 
*                      resolved
*  - SCIP_FEASIBLE   : all constraints of the handler are feasible
*/
 virtual SCIP_DECL_CONSENFOLP(scip_enfolp) override;

/** separation method of constraint handler for LP solution
*
*  Separates all constraints of the constraint handler. The method is called in 
*  the LP solution loop, which means that a valid LP solution exists.
*
*  In this function we add new user cut (uc) in the model.
*  NOTE: we are sure that a new uc can be added, because first a SCIP_CHECK
*  method has declared that separation is possible. 
*
*  possible return values for *result (if more than one applies, the first in 
*  the list should be used):
*  - SCIP_CUTOFF     : the node is infeasible in the variable's bounds and 
*                      can be cut off
*  - SCIP_SEPARATED  : a cutting plane was generated
*  - SCIP_INFEASIBLE : at least one constraint is infeasible, but it was not 
*                      resolved
*  - SCIP_FEASIBLE   : all constraints of the handler are feasible
*/
   virtual SCIP_DECL_CONSSEPALP(scip_sepalp) override;

/** constraint enforcing method of constraint handler for pseudo solutions
*
*  The method is called at the end of the node processing loop for a node 
*  where the LP was not solved. The pseudo solution has to be checked for 
*  feasibility. If possible, an infeasibility should be resolved by
*  branching, reducing a variable's domain to exclude the solution or adding 
*  an additional constraint. Separation is not possible, since the LP is 
*  not processed at the current node. All LP informations like
*  LP solution, slack values, or reduced costs are invalid and must not 
*  be accessed.
*
*  NOTE: At the moment in the solving loop of the algorithm SMS++ can't
*  separate pseudo-solution. For this reason this function is not yet
*  implemented
*
* Possible return values for *result:
*  - SCIP_DIDNOTRUN  : the enforcement was skipped 
*/
   virtual SCIP_DECL_CONSENFOPS(scip_enfops) override;

/** feasibility check method of constraint handler for primal solutions
*
*  The given solution has to be checked for feasibility.
*
*  In this method we check if possible new user cut / lazy constraint
*  can be added to the model.
* 
*  Possible return values for *result:
*  - SCIP_INFEASIBLE : at least one constraint of the handler is infeasible
*                      (i.e., a new cut can be added)
*  - SCIP_FEASIBLE   : all constraints of the handler are feasible
*/
   virtual SCIP_DECL_CONSCHECK(scip_check) override;

/** variable rounding lock method of constraint handler
*
*  This method is called, after a constraint is added or removed from 
*  the transformed problem. It should update the rounding locks of all 
*  associated variables with calls to SCIPaddVarLocksType(),
*  depending on the way, the variable is involved in the constraint:
*  - If the constraint may get violated by decreasing the value of a 
*    variable, it should call SCIPaddVarLocksType(scip, var, 
*    SCIP_LOCKTYPE_MODEL, nlockspos, nlocksneg), saying that rounding 
*    down is potentially rendering the (positive) constraint infeasible 
*    and rounding up is potentially rendering the negation of the constraint 
*    infeasible.
*  - If the constraint may get violated by increasing the value of a variable, 
*    it should call SCIPaddVarLocksType(scip, var, SCIP_LOCKTYPE_MODEL, 
*    nlocksneg, nlockspos), saying that rounding up is potentially rendering 
*    the constraint's negation infeasible and rounding up is potentially 
*    rendering the constraint itself infeasible.
*  - If the constraint may get violated by changing the variable in any direction,
*    it should call SCIPaddVarLocksType(scip, var, SCIP_LOCKTYPE_MODEL, 
*    nlockspos + nlocksneg, nlockspos + nlocksneg).
*
*    NOTE: when creating the caallback method, we don't know which variables 
*    may get involved in future cuts or lazy constraints and in which direction.
*    Thus, the last method will be always called on all the problem variables,
*    in order to avoid SCIP from fixing them. 
*/
   virtual SCIP_DECL_CONSLOCK(scip_lock) override;

/** @} ---------------------------------------------------------------------*/
/*--------------------- ADDITIONAL CALLBACK METHODS ------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Additional Callback Methods
 *  @{ */

/** transforms constraint data into data belonging to the transformed problem */
   virtual SCIP_DECL_CONSTRANS(scip_trans) override;

/** frees specific constraint data */
   virtual SCIP_DECL_CONSDELETE(scip_delete) override;

 protected:

/*--------------------------------------------------------------------------*/
/*-------------------- PROTECTED METHODS OF THE CLASS ----------------------*/
/*--------------------------------------------------------------------------*/

 /** bitwise-encoded parameter for deciding if and when separation of user
  * cuts and lazy constraints is performed */
 unsigned char CutSepPar;

 /* parent *milpsolver from which the separator* is called */
 SMSpp_di_unipi_it::SCIPMILPSolver* parent_scipmilpsolver;

/*--------------------------------------------------------------------------*/

 };  // end( class SCIPMILPSolver_Conhdlr )

/** @} ---------------------------------------------------------------------*/
/*---------------------- METHODS FOR ADDING A CALLBACK ---------------------*/
/*--------------------------------------------------------------------------*/
/** @name Methods for add a new constraint which will eventually be 
 *  enforced or separated producing new cuts or lazy constraints.
 *  @{ */

/** creates and captures a constraint used which will be used as a separator */
 SCIP_RETCODE SCIPcreateSCIPMILPSolver_cb(
   SCIP*        scip,               /**< SCIP data structure */
   SCIP_CONS**  cons,               /**< pointer to hold the created 
                                         constraint */
   const char*  name,               /**< name of constraint */
   std::vector< SCIP_VAR * > vars,  /**< SCIP vars */
   SCIP_Bool    initial,            /**< should the LP relaxation of 
                                         constraint be in the initial LP? */
   SCIP_Bool    separate,           /**< should the constraint be 
                                         separated during LP processing? */
   SCIP_Bool    enforce,            /**< should the constraint be enforced 
                                         during node processing? */
   SCIP_Bool    check,              /**< should the constraint be checked 
                                         for feasibility? */
   SCIP_Bool    propagate,          /**< should the constraint be propagated 
                                         during node processing? */
   SCIP_Bool    local,              /**< is constraint only valid locally? */
   SCIP_Bool    modifiable,         /**< is constraint modifiable (subject 
                                         to column generation)? */
   SCIP_Bool    dynamic,            /**< is constraint dynamic? */
   SCIP_Bool    removable           /**< should the constraint be removed 
                                         from the LP due to aging or cleanup? */
   );

/** creates and captures a a constraint which will be used as a separator
 *  with all its constraint flags set to their default values */
SCIP_RETCODE SCIPcreateSCIPMILPSolver_basiccb(
   SCIP*        scip,               /**< SCIP data structure */
   SCIP_CONS**  cons,               /**< pointer to hold the created constraint */
   const char*  name,               /**< name of constraint */
   std::vector< SCIP_VAR * > vars   /**< SCIP vars */
   );

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

#endif  /* SCIPMILPSolver_Conhdlr.h included */

/*--------------------------------------------------------------------------*/
/*------------ End methods of Class SCIPMILPSolver_Conhdlr -----------------*/
/*--------------------------------------------------------------------------*/


/*--------------------------------------------------------------------------*/

}  // end( namespace SMSpp_di_unipi_it )

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

#endif  /* SCIPMILPSolver.h included */

/*--------------------------------------------------------------------------*/
/*----------------------- End File SCIPMILPSolver.h ------------------------*/
/*--------------------------------------------------------------------------*/
