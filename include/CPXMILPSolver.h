/*--------------------------------------------------------------------------*/
/*--------------------------- File CPXMILPSolver.h -------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the CPXMILPSolver class.
 *
 * CPXMILPSolver implements a general purpose solver that is able to tackle a
 * MILP problem expressed by a Block using IBM CLPEX.
 *
 * \version 0.90
 *
 * \date 14 - 06 - 2019
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
 * \author Kostas Tavlaridis-Gyparakis \n
 *         Operations Research Group \n
 *         Dipartimento di Informatica \n
 *         Università di Pisa \n
 *
 * \copyright &copy; Antonio Frangioni, Kostas Tavlaridis-Gyparakis, Niccolò Iardella
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
 * The CPXMILPSolver can be registered to any kind of Block (9assuming that
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

 /** Types of string parameters.
  * Public enum describing the different types of algorithmic parameters
  * of "string" type that the CPXMILPSolver has.
  */
 enum str_par_type_CPXS {
  strProblemName = strLastAlgParMILP, ///< Problem name
  strOutputFile,                      ///< Output .lp file
  strLastAlgParCPXS
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

 /// It sets the Block that the Solver has to solve and initializes CPLEX.
 void set_Block( Block * block ) override;

 /// Optimizes the problem with CPLEX.
 int compute( bool changedvars ) override;

 OFValue get_lb() override;

 OFValue get_ub() override;

 bool has_var_solution() override;

 bool is_var_feasible() override;

 OFValue get_var_value() override;

 void get_var_solution( Configuration * solc ) override;

 bool has_dual_solution() override;

 void get_dual_solution( Configuration * solc ) override;

 void set_par( idx_type par, int value ) override;

 void set_par( idx_type par, double value ) override;

 void set_par( idx_type par, const std::string & value ) override;

 void write_lp(const std::string & filename) override;
 /// @}

/*--------------------------------------------------------------------------*/
/*------------------- METHODS FOR HANDLING THE PARAMETERS ------------------*/
/*--------------------------------------------------------------------------*/
/**
 * @name Methods for handling parameters
 * @{
 */

 idx_type get_num_str_par() const override;

 const std::string & get_str_par( idx_type par ) const override;

 idx_type str_par_str2idx( const std::string & name ) const override;

 const std::string & dbl_par_idx2str( idx_type idx ) const override;
 /// @}

/*--------------------------------------------------------------------------*/
/*-------------------- PROTECTED FIELDS OF THE CLASS -----------------------*/
/*--------------------------------------------------------------------------*/

 protected:

 CPXENVptr env; /// CPLEX environment
 CPXLPptr milp; /// CPLEX LP problem

 std::string prob_name;   /// CPLEX problem name
 std::string output_file; /// Output file for CPXwriteprob

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
 void of_modification( ObjectiveMod * mod ) override;

 /// It handles a constraint modification
 void const_modification( ConstraintMod * mod ) override;

 /// It handles a bound modification
 void bound_modification( OneVarConstraintMod * mod ) override;

 /// It handles a function modification
 void function_modification( FunctionMod * mod ) override;

 /// It handles a dynamic modification
 void dynamic_modification( BlockModAD * mod ) override;

 /// It adds a single new dynamic constraint
 void add_dynamic_constraint( FRowConstraint * p_const ) override;

 /// It adds a single new dynamic variable
 void add_dynamic_variable( ColVariable * p_var ) override;

 /// It removes a single dynamic constraint
 void remove_dynamic_constraint( FRowConstraint * p_const ) override;

 /// It removes a single dynamic variable
 void remove_dynamic_variable( ColVariable * p_var ) override;
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
