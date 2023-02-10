/*--------------------------------------------------------------------------*/
/*--------------------------- File test_dual.cpp ---------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * This file contains the implementation of a series of tests to verify the
 * dual solution returned by CPXMILPSolver. These tests are intended to ensure
 * that the convention established by RowConstraint regarding dual solutions
 * is respected by the Solver.
 *
 * \author Rafael Durbano Lobato \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * Copyright &copy; by Rafael Durbano Lobato
 */

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "AbstractBlock.h"
#include "CPXMILPSolver.h"
#include "FRealObjective.h"
#include "FRowConstraint.h"
#include "LinearFunction.h"

/*--------------------------------------------------------------------------*/
/*--------------------------- NAMESPACE ------------------------------------*/
/*--------------------------------------------------------------------------*/

// namespace for the tests of the Structured Modeling System++ (SMS++)
using namespace SMSpp_di_unipi_it;

/*--------------------------------------------------------------------------*/
/*------------------------------- FUNCTIONS --------------------------------*/
/*--------------------------------------------------------------------------*/

CPXMILPSolver * build_solver() {
 auto solver = new CPXMILPSolver;
 solver->set_par( CPXPARAM_Preprocessing_Presolve , 0 );
 solver->set_par( CPXPARAM_Preprocessing_QCPDuals , CPX_QCPDUALS_FORCE );
 solver->set_par( CPXMILPSolver::intThrowReducedCostException , 1 );
 return( solver );
}

/*--------------------------------------------------------------------------*/

double get_obj_sign( Objective::of_type sense ) {
 switch( sense ) {
  case( Objective::eMin ): return( -1 );
  case( Objective::eMax ): return( 1 );
  default: return( 0 );
 }
}

/*--------------------------------------------------------------------------*/

void test_lower( Objective::of_type sense , const bool frow_constraint ) {

 /* min/max [+/-] x
  * s.t. x >= 1
  */

 auto lp = new AbstractBlock;

 // Variable

 auto x = new ColVariable;
 lp->add_static_variable( *x , "x" );

 // Objective

 const auto obj_sign = get_obj_sign( sense );

 auto objective_function = new LinearFunction();
 objective_function->add_variable( x , - obj_sign );

 auto objective = new FRealObjective( lp , objective_function );
 objective->set_sense( sense );

 lp->set_objective( objective );

 // Constraint

 RowConstraint * constraint{};

 if( frow_constraint ) {
  auto constraint_ = new FRowConstraint;
  auto function = new LinearFunction();
  function->add_variable( x , 1.0 );
  constraint_->set_function( function );
  constraint_->set_lhs( 1.0 );
  constraint_->set_rhs( Inf<double>() );
  lp->add_static_constraint( *constraint_ , "constraint" );
  constraint = constraint_;
 }
 else {
  auto constraint_ = new LBConstraint;
  constraint_->set_variable( x );
  constraint_->set_lhs( 1.0 );
  constraint_->set_rhs( Inf<double>() );
  lp->add_static_constraint( *constraint_ , "constraint" );
  constraint = constraint_;
 }

 // Solver

 auto solver = build_solver();

 lp->register_Solver( solver );

 // Solve and check solution

 auto status = solver->compute();

 assert( status == Solver::kOK );

 assert( solver->has_var_solution() );
 solver->get_var_solution();

 assert( solver->has_dual_solution() );
 solver->get_dual_solution();

 const auto dual = constraint->get_dual();

 assert( obj_sign * dual == 1 );

 delete lp;
}

/*--------------------------------------------------------------------------*/

void test_upper( Objective::of_type sense , const bool frow_constraint ) {

 /* min/max [-/+] x
  * s.t. x <= 1
  */

 auto lp = new AbstractBlock;

 // Variable

 auto x = new ColVariable;
 lp->add_static_variable( *x , "x" );

 // Objective

 const auto obj_sign = get_obj_sign( sense );

 auto objective_function = new LinearFunction();
 objective_function->add_variable( x , obj_sign );

 auto objective = new FRealObjective( lp , objective_function );
 objective->set_sense( sense );

 lp->set_objective( objective );

 // Constraint

 RowConstraint * constraint{};

 if( frow_constraint ) {
  auto constraint_ = new FRowConstraint;
  auto function = new LinearFunction();
  function->add_variable( x , 1.0 );
  constraint_->set_function( function );
  constraint_->set_lhs( -Inf<double>() );
  constraint_->set_rhs( 1.0 );
  lp->add_static_constraint( *constraint_ , "constraint" );
  constraint = constraint_;
 }
 else {
  auto constraint_ = new UBConstraint;
  constraint_->set_variable( x );
  constraint_->set_lhs( -Inf<double>() );
  constraint_->set_rhs( 1.0 );
  lp->add_static_constraint( *constraint_ , "constraint" );
  constraint = constraint_;
 }

 // Solver

 auto solver = build_solver();

 lp->register_Solver( solver );

 // Solve and check solution

 auto status = solver->compute();

 assert( status == Solver::kOK );

 assert( solver->has_var_solution() );
 solver->get_var_solution();

 assert( solver->has_dual_solution() );
 solver->get_dual_solution();

 const auto dual = constraint->get_dual();

 assert( obj_sign * dual == - 1 );

 delete lp;
}

/*--------------------------------------------------------------------------*/

void test_equality( Objective::of_type sense , const bool frow_constraint ) {

 /* min/max x
  * s.t. x = 1
  */

 auto lp = new AbstractBlock;

 // Variable

 auto x = new ColVariable;
 lp->add_static_variable( *x , "x" );

 // Objective

 auto objective_function = new LinearFunction();
 objective_function->add_variable( x , 1.0 );

 auto objective = new FRealObjective( lp , objective_function );
 objective->set_sense( sense );

 lp->set_objective( objective );

 // Constraint

 RowConstraint * constraint{};

 if( frow_constraint ) {
  auto constraint_ = new FRowConstraint;
  auto function = new LinearFunction();
  function->add_variable( x , 1.0 );
  constraint_->set_function( function );
  constraint_->set_both( 1.0 );
  lp->add_static_constraint( *constraint_ , "constraint" );
  constraint = constraint_;
 }
 else {
  auto constraint_ = new BoxConstraint;
  constraint_->set_variable( x );
  constraint_->set_both( 1.0 );
  lp->add_static_constraint( *constraint_ , "constraint" );
  constraint = constraint_;
 }

 // Solver

 auto solver = build_solver();

 lp->register_Solver( solver );

 // Solve and check solution

 auto status = solver->compute();

 assert( status == Solver::kOK );

 assert( solver->has_var_solution() );
 solver->get_var_solution();

 assert( solver->has_dual_solution() );
 solver->get_dual_solution();

 const auto dual = constraint->get_dual();

 assert( dual == -1 );

 delete lp;
}

/*--------------------------------------------------------------------------*/

void run() {

 // Test dual solution

 for( auto frow_constraint : { true , false } ) {
  for( auto sense : { Objective::eMin , Objective::eMax } ) {
   test_lower( sense , frow_constraint );
   test_upper( sense , frow_constraint );
   test_equality( sense , frow_constraint );
  }
 }
}

/*--------------------------------------------------------------------------*/

int main() {
 run();
 return( 0 );
}

/*--------------------------------------------------------------------------*/
/*------------------------ End File test_dual.cpp --------------------------*/
/*--------------------------------------------------------------------------*/
