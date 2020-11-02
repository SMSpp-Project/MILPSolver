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
 /** Types of integer parameters.
  * Public enum describing the different types of algorithmic parameters
  * of "int" type that the SCIPMILPSolver has.
  */
 enum int_par_type_SCPS {
  intUseCustomNames = intLastAlgParMILP, ///< Use custom names for rows/columns
  intLastAlgParSCPS
 };

 /** Types of double parameters.
* Public enum describing the different types of algorithmic parameters
* of "double" type that the SCIPMILPSolver has.
*/
 enum dbl_par_type_CPXS {
  dblFirstCPLEXPar = dblLastAlgParMILP, ///< First double parameter
  dblLastAlgParSCPS
 };

 /** Types of string parameters.
  * Public enum describing the different types of algorithmic parameters
  * of "string" type that the SCIPMILPSolver has.
  */
 enum str_par_type_SCPS {
  strProblemName = strLastAlgParMILP, ///< Problem name
  strOutputFile,                      ///< Output .lp file
  strLastAlgParSCPS
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

 /// Optimizes the problem with CPLEX.
 int compute( bool changedvars = true ) override;

 OFValue get_lb() override;

 OFValue get_ub() override;

 bool has_var_solution() override;

 bool is_var_feasible() override;

 OFValue get_var_value() override;

 void get_var_solution( Configuration * solc = nullptr ) override;

 bool has_dual_solution() override;

 void get_dual_solution( Configuration * solc = nullptr ) override;

 void write_lp( const std::string & filename ) override;
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

 /// Gets the number of string parameters
 [[nodiscard]] idx_type get_num_str_par() const override;

 /// Gets the number of double parameters
 [[nodiscard]] idx_type get_num_dbl_par() const override;

 /// Gets the default value of the specified integer parameter
 [[nodiscard]] int get_dflt_int_par( idx_type par ) const override;

 /// Gets the default value of the specified integer parameter
 [[nodiscard]] double get_dflt_dbl_par( idx_type par ) const override;

 /// Gets the default value of the specified integer parameter
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
 SCIP * scip{};

 std::vector< SCIP_VAR * > vars;   ///< SCIP variables
 std::vector< SCIP_CONS * > conss; ///< SCIP constraints

 std::string prob_name;   ///< SCIP problem name
 std::string output_file; ///< Output file

 /// If Variable/Constraint names should be used
 bool use_custom_names = true;

 /** @name Clear and load the problem
  *
  * The following two methods include the main logic of the class.
  * They are called in the set_Block() method to build the problem when the
  * Solver is [re]registered to a Block, but also when a NBModification is
  * processed.
  * Note: These two override the base class methods but call them after running
  * the specific operations for this class.
  *
  * @{
  */

 void clear_problem() override; ///< It clears all the LP vectors

 void load_problem() override;  ///< It loads all the LP vectors
 /// @}

 /** @name Get variable bounds for the problem
 *
 * The following two methods retrieve the upper and lower bound for the
 * given variable considering both the Variable bounds and all the active
 * OneVarConstraints active for that Variable.
 */

 /// Gets the LB fot the given variable in the problem
 double get_problem_lb(const ColVariable & var) override;

 /// Gets the UB fot the given variable in the problem
 double get_problem_ub(const ColVariable & var) override;
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
 void of_modification( ObjectiveMod * mod ) override;

 /// It handles a constraint modification
 void const_modification( ConstraintMod * mod ) override;

 /// It handles a bound modification
 void bound_modification( OneVarConstraintMod * mod ) override;

 /// It handles a function modification
 void function_modification( FunctionMod * mod ) override;

 /// It handles a function vars modification
 void function_vars_modification( FunctionModVars * mod ) override;

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
  * It sets the value of a Colvariable taking it from x[i], then increments i.
  * This method is meant to be used inside get_var_solution(), in conjunction
  * with un_any_const_static() or un_any_const_dynamic().
  *
  * @param lvar The ColVariable to set
  * @param x The array containing the values
  * @param i The position in the array
  */
 void set_var_value( ColVariable & lvar, double * x, int & i );

 /**
  * It sets the dual value of a FRowConstraint taking it from pi[i],
  * then increments i.
  * This method is meant to be used inside get_dual_solution(), in conjunction
  * with un_any_const_static() or un_any_const_dynamic().
  *
  * @param lconst The FRowConstraint to set
  * @param pi The array containing the values
  * @param i The position in the array
  */
 void set_dual_value( FRowConstraint & lconst, double * pi, int & i );

 // void fix_integer_vars();

 SMSpp_insert_in_factory_h;
};

}

#endif
