/*--------------------------------------------------------------------------*/
/*--------------------------- File test_dual.cpp ---------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * This file contains the implementation of a series of tests to verify the
 * dual solution returned by *MILPSolver. These tests are intended to ensure
 * that the convention established by RowConstraint regarding dual solutions
 * is respected by the Solver.
 *
 * \author Rafael Durbano Lobato \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 * 
 * \author Enrico Calandrini \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \copyright &copy; by Rafael Durbano Lobato
 */

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include <cmath>
#include <iostream>
#include <list>

#include "AbstractBlock.h"

#include "BlockSolverConfig.h"

#include "FRealObjective.h"

#include "FRowConstraint.h"

#include "LinearFunction.h"

#include "CDASolver.h"

#include "OneVarConstraint.h"
#include "Solver.h"

/*--------------------------------------------------------------------------*/
/*--------------------------- NAMESPACE ------------------------------------*/
/*--------------------------------------------------------------------------*/

// namespace for the tests of the Structured Modeling System++ (SMS++)
using namespace SMSpp_di_unipi_it;

/*--------------------------------------------------------------------------*/
/*------------------------------- FUNCTIONS --------------------------------*/
/*--------------------------------------------------------------------------*/
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
  constraint_->set_rhs( Inf< double >() );
  lp->add_static_constraint( *constraint_ , "constraint" );
  constraint = constraint_;
 }
 else {
  auto constraint_ = new LBConstraint;
  constraint_->set_variable( x );
  constraint_->set_lhs( 1.0 );
  constraint_->set_rhs( Inf< double >() );
  lp->add_static_constraint( *constraint_ , "constraint" );
  constraint = constraint_;
 }

 // Solver

 auto lpbsc = dynamic_cast< BlockSolverConfig * >(
		     Configuration::deserialize( "LPPar-dual.txt" ) );
 if( ! lpbsc ) {
  std::cerr << "Error: configuration file not a BlockSolverConfig" << std::endl;
  exit( 1 );    
  }

 lpbsc->apply( lp );
 lpbsc->clear();  // keep the clear()-ed BlockSolverConfig for final cleanup

 // Solve and check solution

 Solver * solver = ( lp->get_registered_solvers() ).front();

 auto CDASp = dynamic_cast< CDASolver * >( solver );
 
 auto status = CDASp->compute();

 assert( status == Solver::kOK );

 assert( CDASp->has_var_solution() );
 solver->get_var_solution();

 assert( CDASp->has_dual_solution() );
 CDASp->get_dual_solution();

 const auto dual = constraint->get_dual();

 assert( obj_sign * dual == 1 );
 
 lpbsc->apply( lp );

 delete lpbsc;
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
  constraint_->set_lhs( -Inf< double >() );
  constraint_->set_rhs( 1.0 );
  lp->add_static_constraint( *constraint_ , "constraint" );
  constraint = constraint_;
 }
 else {
  auto constraint_ = new UBConstraint;
  constraint_->set_variable( x );
  constraint_->set_lhs( -Inf< double >() );
  constraint_->set_rhs( 1.0 );
  lp->add_static_constraint( *constraint_ , "constraint" );
  constraint = constraint_;
 }

 // Solver

 auto lpbsc = dynamic_cast< BlockSolverConfig * >(
		     Configuration::deserialize( "LPPar-dual.txt" ) );
 if( ! lpbsc ) {
  std::cerr << "Error: configuration file not a BlockSolverConfig" << std::endl;
  exit( 1 );    
  }

 lpbsc->apply( lp );
 lpbsc->clear();  // keep the clear()-ed BlockSolverConfig for final cleanup

 // Solve and check solution

 Solver * solver = ( lp->get_registered_solvers() ).front();
 
 auto CDASp = dynamic_cast< CDASolver * >( solver );
 
 auto status = CDASp->compute();

 assert( CDASp->has_var_solution() );
 CDASp->get_var_solution();

 assert( CDASp->has_dual_solution() );
 CDASp->get_dual_solution();

 const auto dual = constraint->get_dual();

 assert( obj_sign * dual == - 1 );

 lpbsc->apply( lp );

 delete lpbsc;
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

 auto lpbsc = dynamic_cast< BlockSolverConfig * >(
		     Configuration::deserialize( "LPPar-dual.txt" ) );
 if( ! lpbsc ) {
  std::cerr << "Error: configuration file not a BlockSolverConfig" << std::endl;
  exit( 1 );    
  }

 lpbsc->apply( lp );
 lpbsc->clear();  // keep the clear()-ed BlockSolverConfig for final cleanup

 // Solve and check solution

 Solver * solver = ( lp->get_registered_solvers() ).front();
 
 auto CDASp = dynamic_cast< CDASolver * >( solver );
 
 auto status = CDASp->compute();

 assert( CDASp->has_var_solution() );
 CDASp->get_var_solution();

 assert( CDASp->has_dual_solution() );
 CDASp->get_dual_solution();

 const auto dual = constraint->get_dual();

 assert( dual == -1 );

 lpbsc->apply( lp );

 delete lpbsc;
 delete lp;
}

/*--------------------------------------------------------------------------*/

int test_dynamic_in_wrapper( void ) {
 /* The coupling row of a Benders subproblem is a dynamic Constraint of a
  * Block whose sub-Blocks carry the model:
  *
  *    wrapper                     y[ 1 ] <= 1   (dynamic)
  *     +-- inner
  *          +-- sub 0    min - y[ 0 ] , y[ 0 ] <= 2 , 0 <= y[ 0 ] <= 10
  *          +-- sub 1    min - 2 y[ 1 ] , y[ 1 ] <= 3 , 0 <= y[ 1 ] <= 10
  *
  * whose optimum is - 2 - 2 = - 4 and whose derivative in the right-hand
  * side of the dynamic row is - 2, that row being the one that binds
  * y[ 1 ]. The dynamic rows reach the Solver after all the static ones,
  * whatever Block they belong to, so a dual written by walking the groups
  * of each Block in turn lands on the wrong row. The checks are explicit
  * rather than assert()-ed, an optimized build having no assert(). */

 int failures = 0;

 for( const std::string name :
       { "CPXMILPSolver" , "GRBMILPSolver" , "HiGHSMILPSolver" } ) {
  Solver * solver;
  try {
   solver = Solver::new_Solver( name );
   }
  catch( ... ) {
   std::cout << name << ": not in the factory, skipped" << std::endl;
   continue;
   }
  solver->set_par( Solver::intLogVerb , 0 );

  auto wrapper = new AbstractBlock;
  auto inner = new AbstractBlock( wrapper );
  wrapper->add_nested_Block( inner );

  ColVariable * y[ 2 ];
  FRowConstraint * cap[ 2 ];

  for( int k = 0 ; k < 2 ; ++k ) {
   auto sub = new AbstractBlock( inner );

   y[ k ] = new ColVariable;
   y[ k ]->set_type( ColVariable::kContinuous );
   sub->add_static_variable( * y[ k ] , "y" );

   auto box = new BoxConstraint;
   box->set_variable( y[ k ] );
   box->set_lhs( 0 );
   box->set_rhs( 10 );
   sub->add_static_constraint( * box , "box" );

   cap[ k ] = new FRowConstraint;
   auto cf = new LinearFunction();
   cf->add_variable( y[ k ] , 1.0 );
   cap[ k ]->set_function( cf );
   cap[ k ]->set_lhs( - Inf< double >() );
   cap[ k ]->set_rhs( 2.0 + k );
   sub->add_static_constraint( * cap[ k ] , "cap" );

   auto of = new LinearFunction();
   of->add_variable( y[ k ] , - 1.0 - k );
   auto objective = new FRealObjective( sub , of );
   objective->set_sense( Objective::eMin );
   sub->set_objective( objective );

   inner->add_nested_Block( sub );
   }

  auto wobj = new FRealObjective( wrapper , new LinearFunction() );
  wobj->set_sense( Objective::eMin );
  wrapper->set_objective( wobj );

  auto link = new std::list< FRowConstraint >( 1 );
  auto lf = new LinearFunction();
  lf->add_variable( y[ 1 ] , 1.0 );
  link->front().set_function( lf );
  link->front().set_lhs( - Inf< double >() );
  link->front().set_rhs( 1.0 );
  wrapper->add_dynamic_constraint( * link , "link" );

  wrapper->register_Solver( solver );

  const auto status = solver->compute();
  auto CDASp = dynamic_cast< CDASolver * >( solver );

  bool ok = ( status == Solver::kOK ) && CDASp && CDASp->has_dual_solution();
  if( ok ) {
   solver->get_var_solution();
   CDASp->get_dual_solution();

   ok = ( std::abs( solver->get_var_value() + 4.0 ) <= 1e-9 ) &&
        ( std::abs( link->front().get_dual() - 2.0 ) <= 1e-9 ) &&
        ( std::abs( cap[ 0 ]->get_dual() - 1.0 ) <= 1e-9 ) &&
        ( std::abs( cap[ 1 ]->get_dual() ) <= 1e-9 );
   }

  std::cout << name << ": value = " << solver->get_var_value()
            << ", dual of the dynamic row = " << link->front().get_dual()
            << " (2 expected), of the static ones = "
            << cap[ 0 ]->get_dual() << " (1) and " << cap[ 1 ]->get_dual()
            << " (0)" << ( ok ? " -> OK" : " -> KO" ) << std::endl;

  if( ! ok )
   ++failures;

  wrapper->unregister_Solver( solver );
  delete solver;
  delete wrapper;
  }

 return( failures );
 }

/*--------------------------------------------------------------------------*/

int run() {

 // Test dual solution

 for( auto frow_constraint : { true , false } ) {
  for( auto sense : { Objective::eMin , Objective::eMax } ) {
   test_lower( sense , frow_constraint );
   test_upper( sense , frow_constraint );
   test_equality( sense , frow_constraint );
  }
 }

 // the dual of a dynamic row of a Block whose sub-Blocks carry the model
 return( test_dynamic_in_wrapper() );
}

/*--------------------------------------------------------------------------*/

int main( int argc , char ** argv )
{
 return( run() );
}

/*--------------------------------------------------------------------------*/
/*------------------------ End File test_dual.cpp --------------------------*/
/*--------------------------------------------------------------------------*/
