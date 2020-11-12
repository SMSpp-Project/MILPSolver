/*--------------------------------------------------------------------------*/
/*--------------------------- File CPXMILPSolver.h -------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the CPXMILPSolver class.
 *
 * CPXMILPSolver implements a general purpose solver that is able to tackle a
 * MILP problem expressed by a Block using IBM CLPEX.
 *
 * \author Antonio Frangioni \n
 *         Operations Research Group \n
 *         Dipartimento di Informatica \n
 *         Università di Pisa \n
 *
 * \author Niccolò Iardella \n
 *         Operations Research Group \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \copyright &copy; Antonio Frangioni, Niccolò Iardella
 */
/*--------------------------------------------------------------------------*/
/*----------------------------- DEFINITIONS --------------------------------*/
/*--------------------------------------------------------------------------*/

#ifndef __CPXMILPSOLVER_H
#define __CPXMILPSOLVER_H

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include <ilcplex/cplex.h>

#include <SMSTypedefs.h>
#include <Observer.h>
#include <Block.h>
#include <Solver.h>
#include <ColVariable.h>
#include <FRealObjective.h>
#include <FRowConstraint.h>
#include <OneVarConstraint.h>

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

/*--------------------------------------------------------------------------*/
/*----------------------- CLASS CPXMILPSolver ------------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------- GENERAL NOTES --------------------------------*/
/*--------------------------------------------------------------------------*/

/// Class for solving MILP problems via CPLEX.
/**
 * The CPXMILPSolver class derives from MILPSolver and extends the
 * base class to solve MILP problems using CPLEX.
 *
 * The CPXMILPSolver can be registered to any kind of Block (assuming that
 * it contains a MILP formulation) and it uses the base class functionalities
 * to build a matricial representation of a MILP problem, then it solves it
 * using CPLEX through its Callable Library. Moreover, it implements the
 * interface that MILPSolver provides for processing modifications.
 *
 * The main logic is in compute(). This method copies the vectors that decribe
 * the LP problem into a CPLEX environment, it processes the modifications and
 * then solves the problem.
 * get_var_solution() retrieves the values of the variables from CPLEX, saves
 * them into the Block variables and evaluates the objective function.
 *
 * Besides the configuration parameters already present in Solver, this class
 * adds two string parameters that allow to specify the CPLEX name of the
 * problem and an output file for CPLEX to write out the problem.
 * Moreover, the user can include in the configuration all the parameters
 * supported by CPXsetintparam(), CPXsetdblparam() and CPXsetstrparam().
 * (See the CPLEX Callable Library reference manual for all of them)
 */
class CPXMILPSolver : public MILPSolver {

/*--------------------------------------------------------------------------*/
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

 public:

/*--------------------------------------------------------------------------*/
/*---------------------------- PUBLIC TYPES --------------------------------*/
/*--------------------------------------------------------------------------*/

 /// Types of integer parameters
 enum int_par_type_CPXS {
  /// First CPLEX int/long parameter
  intFirstCPLEXPar = intLastAlgParMILP,
  /// First allowed new int parameter for derived classes
  intLastAlgParCPXS = intFirstCPLEXPar + CPX_NUM_INT_PARS
 };

 /// Types of double parameters
 enum dbl_par_type_CPXS {
  /// First CPLEX double parameter
  dblFirstCPLEXPar = dblLastAlgParMILP,
  /// First allowed new double parameter for derived classes
  dblLastAlgParCPXS = dblFirstCPLEXPar + CPX_NUM_DBL_PARS
 };

 /// Types of string parameters
 enum str_par_type_CPXS {
  /// First CPLEX string parameter
  strFirstCPLEXPar = strLastAlgParMILP,
  /// First allowed new string parameter for derived classes
  strLastAlgParCPXS = strFirstCPLEXPar + CPX_NUM_STR_PARS
 };

/*--------------------------------------------------------------------------*/
/*--------------------- CONSTRUCTOR AND DESTRUCTOR -------------------------*/
/*--------------------------------------------------------------------------*/

/**
 * @name Constructor and Destructor
 * @{
 */

 CPXMILPSolver();

 ~CPXMILPSolver() override;
 /// @}

/*--------------------------------------------------------------------------*/
/*--------------------- DERIVED METHODS OF BASE CLASS ----------------------*/
/*--------------------------------------------------------------------------*/

/**
 * @name Public Methods derived from base classes
 * @{
 */

 /// It sets the Block that the Solver has to solve and initializes CPLEX
 void set_Block( Block * block ) override;

 /// Optimizes the problem with CPLEX
 int compute( bool changedvars = false ) override;

 /// Returns a valid lower bound on the optimal objective function value
 OFValue get_lb() override;

 /// Returns a valid upper bound on the optimal objective function value
 OFValue get_ub() override;

 /// Returns the value of the current solution, if any
 OFValue get_var_value() override;

 /// Tells whether a solution is available
 bool has_var_solution() override;

 /// Tells whether the current solution is feasible
 bool is_var_feasible() override;

 /// Writes the current solution in the Block
 void get_var_solution( Configuration * solc = nullptr ) override;

 /// Tells whether a dual solution is available
 bool has_dual_solution() override;

 /// Tells whether the current dual solution is feasible
 bool is_dual_feasible() override;

 /// Writes the current dual solution in the Block
 void get_dual_solution( Configuration * solc = nullptr ) override;

 /// Tells whether a dual unbounded direction is available
 bool has_dual_direction() override;

 /// Writes the current dual unbounded direction in the Block
 void get_dual_direction( Configuration * dirc = nullptr ) override;

 /// Writes the LP on the specified file
 void write_lp( const std::string & filename ) override;

 /// Returns the number of nodes used to solve a MIP
 [[nodiscard]] int get_nodes() const override;

 /// It clears the CPLEX environment
 void clear_problem( unsigned int what ) override;

 /// It loads the problem into CPLEX
 void load_problem() override;

#ifdef MILPSOLVER_DEBUG

 /// Check the dictionaries for inconsistencies
 void check_status() override;

#endif
 /// @}

/*--------------------------------------------------------------------------*/
/*------------------- METHODS FOR HANDLING THE PARAMETERS ------------------*/
/*--------------------------------------------------------------------------*/
/**
 * @name Methods for handling parameters
 * @{
 */

 /// Sets an integer parameter with the given value
 void set_par( idx_type par, int value ) override;

 /// Sets a double parameter with the given value
 void set_par( idx_type par, double value ) override;

 /// Sets a string parameter with the given value
 void set_par( idx_type par, const std::string & value ) override;

 /// Gets the number of integer parameters
 [[nodiscard]] idx_type get_num_int_par() const override;

 /// Gets the number of double parameters
 [[nodiscard]] idx_type get_num_dbl_par() const override;

 /// Gets the number of string parameters
 [[nodiscard]] idx_type get_num_str_par() const override;

 /// Gets the default value of the specified integer parameter
 [[nodiscard]] int get_dflt_int_par( idx_type par ) const override;

 /// Gets the default value of the specified double parameter
 [[nodiscard]] double get_dflt_dbl_par( idx_type par ) const override;

 /// Gets the default value of the specified string parameter
 [[nodiscard]] const std::string &
 get_dflt_str_par( idx_type par ) const override;

 /// Gets the value of the specified integer parameter
 [[nodiscard]] int get_int_par( idx_type par ) const override;

 /// Gets the value of the specified double parameter
 [[nodiscard]] double get_dbl_par( idx_type par ) const override;

 /// Gets the value of the specified string parameter
 [[nodiscard]] const std::string & get_str_par( idx_type par ) const override;

 /// Returns the index of the int parameter with the specified name
 [[nodiscard]] idx_type
 int_par_str2idx( const std::string & name ) const override;

 /// Returns the name of the int parameter with the specified index
 [[nodiscard]] const std::string &
 int_par_idx2str( idx_type idx ) const override;

 /// Returns the index of the double parameter with the specified name
 [[nodiscard]] idx_type
 dbl_par_str2idx( const std::string & name ) const override;

 /// Returns the name of the double parameter with the specified index
 [[nodiscard]] const std::string &
 dbl_par_idx2str( idx_type idx ) const override;

 /// Returns the index of the string parameter with the specified name
 [[nodiscard]] idx_type
 str_par_str2idx( const std::string & name ) const override;

 /// Returns the name of the string parameter with the specified index
 [[nodiscard]] const std::string &
 str_par_idx2str( idx_type idx ) const override;
 /// @}

/*--------------------------------------------------------------------------*/
/*-------------------- PROTECTED FIELDS OF THE CLASS -----------------------*/
/*--------------------------------------------------------------------------*/

 protected:

 CPXENVptr env; ///< CPLEX environment
 CPXLPptr lp;   ///< CPLEX LP problem

 /** @name Get variable bounds for the problem
  *
  * The following two methods retrieve the upper and lower bound for the
  * given variable considering both the Variable bounds and all the active
  * OneVarConstraints active for that Variable.
  */

 /// Gets the LB fot the given variable in the problem
 double get_problem_lb( const ColVariable & var ) override;

 /// Gets the UB fot the given variable in the problem
 double get_problem_ub( const ColVariable & var ) override;
 /// @}

 /** @name Handling of CPLEX parameters
  *
  * The following maps are used to keep a relationship between SMS++ parameter
  * system and CPLEX parameters. This allows us to use CPLEX parameters
  * (See CPLEX Parameters Reference Manual from IBM) as they were SMS++
  * parameters with the same names, for example in configuration files.
  *
  * Note: since SMS++ does not support long parameters, both int and
  * long CPLEX parameters are handled as SMS++ int parameters.
  *
  * @{
  */

 const static std::array< int, CPX_NUM_INT_PARS > SMSpp_to_CPLEX_int_pars;
 const static std::array< int, CPX_NUM_DBL_PARS > SMSpp_to_CPLEX_dbl_pars;
 const static std::array< int, CPX_NUM_STR_PARS > SMSpp_to_CPLEX_str_pars;

 const static std::array< std::pair< int, int >, CPX_NUM_INT_PARS > CPLEX_to_SMSpp_int_pars;
 const static std::array< std::pair< int, int >, CPX_NUM_DBL_PARS > CPLEX_to_SMSpp_dbl_pars;
 const static std::array< std::pair< int, int >, CPX_NUM_STR_PARS > CPLEX_to_SMSpp_str_pars;

 /// @}

/*--------------------------------------------------------------------------*/
/*-------------------- METHODS FOR MODIFYING THE PROBLEM -------------------*/
/*--------------------------------------------------------------------------*/

/**
 * @name Methods for modifying the constructed CPLEX problem
 * @{
 */

 /// It handles a variable modification
 void var_modification( VariableMod * mod ) override;

 /// It handles an objective modification
 void objective_modification( ObjectiveMod * mod ) override;

 /// It handles a constraint modification
 void const_modification( ConstraintMod * mod ) override;

 /// It handles a bound modification
 void bound_modification( OneVarConstraintMod * mod ) override;

 /// It handles a function modification applied to the objective
 void objective_function_modification( FunctionMod * mod ) override;

 /// It handles a function modification applied to a constraint
 void constraint_function_modification( FunctionMod * mod ) override;

 /// It handles a function vars modification to the objective
 void objective_fvars_modification( FunctionModVars * mod ) override;

 /// It handles a function vars modification to a constraint
 void constraint_fvars_modification( FunctionModVars * mod ) override;

 /// It handles a dynamic modification
 void dynamic_modification( BlockModAD * mod ) override;

 /// It adds a single new dynamic constraint
 void add_dynamic_constraint( FRowConstraint * p_const ) override;

 /// It adds a single new dynamic bound
 void add_dynamic_bound( OneVarConstraint * p_bound ) override;

 /// It adds a single new dynamic variable
 void add_dynamic_variable( ColVariable * p_var ) override;

 /// It removes a single dynamic constraint
 void remove_dynamic_constraint( const FRowConstraint * p_const ) override;

 /// It removes a single dynamic variable
 void remove_dynamic_variable( const ColVariable * p_var ) override;

 /// It removes a single dynamic bound
 void remove_dynamic_bound( const OneVarConstraint * p_bound ) override;
 /// @}

/*--------------------------------------------------------------------------*/
/*--------------------- PRIVATE FIELDS OF THE CLASS ------------------------*/
/*--------------------------------------------------------------------------*/
 private:

/**
  * Gets the number of CPLEX integer vars.
  * Optionally, the read CPLEX variable types can be saved in a vector.
  * Such vector is automatically resized and filled.
  *
  * @param ctype A vector where CPLEX variable types are saved
  */
 int CPXgetintvars( std::vector< char > * ctype = nullptr );

 /**
  * Returns the SMS++ status corresponding to the given
  * CPLEX status returned by CPXgetstat() in case of a LP/QP,
  * or by CPXgetsubstat() in case of a subproblem of a MIP.
  */
 static int decode_lqp_status( int status );

 /**
  * Returns the SMS++ status corresponding to the given
  * CPLEX status returned by CPXgetstat() in case of a MIP.
  */
 static int decode_mip_status( int status );

 /**
  * Returns the SMS++ status corresponding to the given
  * CPLEX error returned by CPXmipopt(), CPXlpopt() or CPXqpopt().
  */
 static int decode_cpx_error( int error );

 SMSpp_insert_in_factory_h;
};

}

#endif
