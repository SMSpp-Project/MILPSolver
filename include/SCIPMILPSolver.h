/*--------------------------------------------------------------------------*/
/*--------------------------- File SCIPMILPSolver.h -------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the SCIPMILPSolver class.
 *
 * SCIPMILPSolver implements a general purpose solver that is able to tackle a
 * MILP problem expressed by a Block using ZIB SCIP.
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

#ifndef __SCIPMILPSOLVER_H
#define __SCIPMILPSOLVER_H

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include <scip/scip.h>

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
#include BOOST_PP_STRINGIZE( BOOST_PP_CAT( BOOST_PP_CAT( SCIP, SCIP_VERSION ), _defs.h ) )

/*--------------------------------------------------------------------------*/
/*----------------------------- NAMESPACE ----------------------------------*/
/*--------------------------------------------------------------------------*/

/// namespace for the Structured Modeling System++ (SMS++)
namespace SMSpp_di_unipi_it {

/*--------------------------------------------------------------------------*/
/*----------------------- CLASS SCIPMILPSolver ------------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------- GENERAL NOTES --------------------------------*/
/*--------------------------------------------------------------------------*/

/// Class for solving MILP problems via SCIP.
/**
 * The SCIPMILPSolver class derives from MILPSolver and extends the
 * base class to solve MILP problems using SCIP.
 *
 * The SCIPMILPSolver can be registered to any kind of Block (assuming that
 * it contains a MILP formulation) and it uses the base class functionalities
 * to build a matricial representation of a MILP problem, then it solves it
 * using SCIP. Moreover, it implements the interface that MILPSolver provides
 * for processing modifications.
 *
 * The main logic is in compute(). This method copies the vectors that decribe
 * the LP problem into a SCIP environment, it processes the modifications and
 * then solves the problem.
 * get_var_solution() retrieves the values of the variables from SCIP, saves
 * them into the Block variables and evaluates the objective function.
 */
class SCIPMILPSolver : public MILPSolver {

/*--------------------------------------------------------------------------*/
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

 public:

/*--------------------------------------------------------------------------*/
/*---------------------------- PUBLIC TYPES --------------------------------*/
/*--------------------------------------------------------------------------*/

 /// Types of integer parameters
 enum int_par_type_SCPS {
  /// First SCIP int/long parameter
  intFirstSCIPPar = intLastAlgParMILP,
  /// First allowed new int parameter for derived classes
  intLastAlgParSCPS = intFirstSCIPPar + SCIP_NUM_INT_PARS
 };

 /// Types of double parameters
 enum dbl_par_type_SCPS {
  /// First SCIP double parameter
  dblFirstSCIPPar = dblLastAlgParMILP,
  /// First allowed new double parameter for derived classes
  dblLastAlgParSCPS = dblFirstSCIPPar + SCIP_NUM_DBL_PARS
 };

 /// Types of string parameters
 enum str_par_type_SCPS {
  /// First SCIP string parameter
  strFirstSCIPPar = strLastAlgParMILP,
  /// First allowed new string parameter for derived classes
  strLastAlgParSCPS = strFirstSCIPPar + SCIP_NUM_STR_PARS
 };

/*--------------------------------------------------------------------------*/
/*--------------------- CONSTRUCTOR AND DESTRUCTOR -------------------------*/
/*--------------------------------------------------------------------------*/

/**
 * @name Constructor and Destructor
 * @{
 */

 SCIPMILPSolver();

 ~SCIPMILPSolver() override;
 /// @}

/*--------------------------------------------------------------------------*/
/*--------------------- DERIVED METHODS OF BASE CLASS ----------------------*/
/*--------------------------------------------------------------------------*/

/**
 * @name Public Methods derived from base classes
 * @{
 */

 /// It sets the Block that the Solver has to solve and initializes CPLEX.
 void set_Block( Block * block ) override;

 /// Optimizes the problem with SCIP
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

 /// It clears the SCIP environment
 void clear_problem( unsigned int what ) override;

 /// It loads the problem into SCIP
 void load_problem() override;
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
 /// SCIP environment
 SCIP * scip{};

 std::vector< SCIP_VAR * > vars;   ///< SCIP variables
 std::vector< SCIP_CONS * > cons;  ///< SCIP constraints

 /// SCIP auxiliary variables for QPs
 std::vector< SCIP_VAR * > aux_vars;
 /// SCIP auxiliary constraints for QPs
 std::vector< SCIP_CONS * > aux_cons;

 /** @name Get variable bounds for the problem
  *
  * The following two methods retrieve the upper and lower bound for the
  * given variable considering both the Variable bounds and all the active
  * OneVarConstraints active for that Variable.
  *
  * @{
  */

 /// Gets the lower bound for the given variable in the problem
 /**
  * The method scans the active bounds of the given variable
  * as well as the lower bound implied by the variable type, returning
  * the value of the one which is most stringent.
  * This method simply calls the overridden one from the MILPSolver class
  * and replaces Inf() values with SCIPinfinity().
  *
  * @param var a variable
  * @param con a pointer where the valid bound is saved. Can be null.
  * @return the valid lower bound
  */
 double
 get_problem_lb( const ColVariable & var, OneVarConstraint * con ) override;

 /// Gets the upper bound for the given variable in the problem
 /**
  * The method scans the active bounds of the given variable
  * as well as the upper bound implied by the variable type, returning
  * the value of the one which is most stringent.
  * This method simply calls the overridden one from the MILPSolver class
  * and replaces Inf() values with SCIPinfinity().
  *
  * @param var a variable
  * @param con a pointer where the valid bound is saved. Can be null.
  * @return the valid upper bound
  */
 double
 get_problem_ub( const ColVariable & var, OneVarConstraint * con ) override;
 /// @}

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
 *
 * @{
 */

 const static std::array< std::string, SCIP_NUM_INT_PARS > SMSpp_to_SCIP_int_pars;
 const static std::array< std::string, SCIP_NUM_DBL_PARS > SMSpp_to_SCIP_dbl_pars;
 const static std::array< std::string, SCIP_NUM_STR_PARS > SMSpp_to_SCIP_str_pars;

 const static std::array< std::pair< std::string, int >, SCIP_NUM_INT_PARS > SCIP_to_SMSpp_int_pars;
 const static std::array< std::pair< std::string, int >, SCIP_NUM_DBL_PARS > SCIP_to_SMSpp_dbl_pars;
 const static std::array< std::pair< std::string, int >, SCIP_NUM_STR_PARS > SCIP_to_SMSpp_str_pars;

 /// @}

/*--------------------------------------------------------------------------*/
/*-------------------- METHODS FOR MODIFYING THE PROBLEM -------------------*/
/*--------------------------------------------------------------------------*/

/**
 * @name Methods for modifying the constructed SCIP problem
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
 void add_dynamic_constraint( FRowConstraint * con ) override;

 /// It adds a single new dynamic bound
 void add_dynamic_bound( OneVarConstraint * con ) override;

 /// It adds a single new dynamic variable
 void add_dynamic_variable( ColVariable * var ) override;

 /// It removes a single dynamic constraint
 void remove_dynamic_constraint( const FRowConstraint * con ) override;

 /// It removes a single dynamic variable
 void remove_dynamic_variable( const ColVariable * var ) override;

 /// It removes a single dynamic bound
 void remove_dynamic_bound( const OneVarConstraint * con ) override;
 /// @}

/*--------------------------------------------------------------------------*/
/*--------------------- PRIVATE FIELDS OF THE CLASS ------------------------*/
/*--------------------------------------------------------------------------*/
 private:

 SMSpp_insert_in_factory_h;
};

}

#endif
