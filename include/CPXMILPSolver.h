/*--------------------------------------------------------------------------*/
/*--------------------------- File CPXMILPSolver.h -------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the CPLEX MILP Solver class, that defines a General Purpose
 * Solver that is able to tackle a MILP Problem expressed by a Block via the
 * use of IMB CLPEX.
 *
 * \version 0.20
 *
 * \date 31 - 05 - 2019
 *
 * \author Antonio Frangioni \n
 *         Operations Research Group \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \author Niccolò Iardella \n
 *         Operations Research Group \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \author Kostas Tavlaridis-Gyparakis \n
 *         Operations Research Group \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * Copyright &copy by Antonio Frangioni, Kostas Tavlaridis-Gyparakis, Niccolò Iardella
 */
/*--------------------------------------------------------------------------*/
/*----------------------------- DEFINITIONS --------------------------------*/
/*--------------------------------------------------------------------------*/

#ifndef __CPXMILPSOLVER_H
#define __CPXMILPSOLVER_H

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include <cplex.h>

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

/// derived class for solving MILP Problems via CPLEX
/**
 * The CPXMILPSolver class derives from MILPSolver [see MILPSolver.h] and
 * extends the base class to solve MILP Problems using CPLEX.
 *
 * The CPXMILPSolver can be registered to any kind of Block, assuming that
 * it contains a MILP formulation, and it uses the base class functionalities
 * to build a matricial representation of a MILP, then it solves the problem
 * using CPLEX through its Callable Library.
 * Moreover, it implements the interface that MILPSolver provides for
 * processing modifications.
 *
 * The main logic is in compute(). This method copies the LP vectors of the
 * class into a CPLEX environment, it processes the modifications and then
 * it solves the problem.
 * get_var_solution() retrieves the values of the variables from CPLEX, saves
 * them into the Block variables and evaluates the objective function. 
 */
class CPXMILPSolver : public MILPSolver {

/*--------------------------------------------------------------------------*/
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

 public:

/*--------------------------------------------------------------------------*/
/*---------------------------- PUBLIC TYPES --------------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Public Types
 *  @{ */

/**
 * Public enum describing the different types of algorithmic parameters
 * of "string" type that the CPXMILPSolver has.
 */
 enum str_par_type_CPXS {
  strProblemName = strLastAlgParMILP, ///< Problem name
  strOutputFile,                     ///< Output .lp file
  strLastAlgParCPXS
 };

/*@}------------------------------------------------------------------------*/
/*--------------------- CONSTRUCTOR AND DESTRUCTOR -------------------------*/
/*--------------------------------------------------------------------------*/

/** @name Constructor and Destructor
 *  @{ */

 CPXMILPSolver();

 ~CPXMILPSolver() override;

/*@}------------------------------------------------------------------------*/
/*--------------------- DERIVED METHODS OF BASE CLASS ----------------------*/
/*--------------------------------------------------------------------------*/

/** @name Public Methods derived of the Base Class
 *  @{ */

 void set_Block(Block* block) override;

 int compute(bool changedvars) override;

 OFValue get_lb() override;

 OFValue get_ub() override;

 void get_var_solution(Configuration* solc) override;

 void get_dual_solution(Configuration* solc) override {}

 void set_par(idx_type par, int value) override;

 void set_par(idx_type par, double value) override;

 void set_par(idx_type par, const std::string & value) override;


/*@}------------------------------------------------------------------------*/
/*------------------- METHODS FOR HANDLING THE PARAMETERS ------------------*/
/*--------------------------------------------------------------------------*/
/** @name Handling the parameters of the CPXMILPSolver
 *  @{ */

 idx_type get_num_str_par() const override;

 const std::string& get_str_par(idx_type par) const override;

 idx_type str_par_str2idx(const std::string& name) const override;

 const std::string& dbl_par_idx2str(idx_type idx) const override;

/*@}------------------------------------------------------------------------*/
/*-------------------- PROTECTED FIELDS OF THE CLASS -----------------------*/
/*--------------------------------------------------------------------------*/

/** @name Protected Fields of the MILPSolver
 *  @{ */

 protected:

 CPXENVptr env; /// CPLEX environment
 CPXLPptr milp; /// CPLEX LP problem

 std::string prob_name;   /// CPLEX problem name
 std::string output_file; /// Output file for CPXwriteprob

 /**
  * It removes the problem from the CPLEX environment
  */
 void clear_problem() override;

 /**
  * It loads the problem into CPLEX environment
  */
 void load_problem() override;

/*@}------------------------------------------------------------------------*/
/*-------------------- METHODS FOR MODIFYING THE PROBLEM -------------------*/
/*--------------------------------------------------------------------------*/

 /** @name Methods for modifying the constructed CPLEX Problem
  *  @{ */

 /// It handles a Variable Modification
 void var_modification(VariableMod* mod) override;

 /// It handles an Objective Modification
 void of_modification(ObjectiveMod* mod) override;

 /// It handles a RowConstraint Modification
 void const_modification(ConstraintMod* mod) override;

 /// It handles a OneVarConstraint Modification
 void bound_modification(OneVarConstraintMod* mod) override;

 /// It handles a Function Modification
 void function_modification(FunctionMod* mod) override;

 /// It handles a dynamic Modification
 void dynamic_modification(BlockModAD* mod) override;

 /// It adds a single new dynamic constraint
 void add_dynamic_constraint(FRowConstraint* p_const) override;

 /// It adds a single new dynamic variable
 void add_dynamic_variable(ColVariable* p_var) override;

 /// It removes a single dynamic constraint
 void remove_dynamic_constraint(FRowConstraint* p_const) override;

 /// It removes a single dynamic variable
 void remove_dynamic_variable(ColVariable* p_var) override;

/*@}------------------------------------------------------------------------*/
/*--------------------- PRIVATE FIELDS OF THE CLASS ------------------------*/
/*--------------------------------------------------------------------------*/
 private:

/**
 * @name Private methods
 *  @{
 */

 /**
 * It sets the value of a Colvariable taking it from tmpx[i], then increments i.
 * This method is meant to be used inside get_var_solution(), in conjunction
 * with un_any_const_static() or un_any_const_dynamic().
 *
 * @param lvar The ColVariable to set
 * @param tmpx The array containing the values
 * @param i The position in the array
 */
 void set_var_value(ColVariable& lvar, double* tmpx, int& i);
/*@}------------------------------------------------------------------------*/
 SMSpp_insert_in_factory_h;
};

}

#endif
