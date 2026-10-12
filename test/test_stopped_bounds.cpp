/*--------------------------------------------------------------------------*/
/*---------------------- File test_stopped_bounds.cpp ----------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Checks the bounds that a :MILPSolver gives on a continuous problem whose
 * solve does not end optimal: the problem
 *
 *     min ( max )  +- ( x_1^2 + x_2^2 + x_3^2 ) - x_1 - 2 x_2 - 3 x_3
 *     s.t.         -10 <= x_i <= 10
 *
 * is solved with dblMaxTime = 0, and the bound on the side of the optimum
 * that the Solver has not proved (get_lb() when minimizing, get_ub() when
 * maximizing) has to be the infinite one, and not a large finite number
 * that a caller would take for a bound; the problem is then solved without
 * the time limit by a new Solver, and both bounds have to be the optimal
 * value. The Solver is GRBMILPSolver, whose bound on a continuous problem
 * comes from GUROBI's OBJBOUND unless NonConvex is 0, and which publishes
 * none there; the test is skipped (exit status 77) if GRBMILPSolver is not
 * in the build, and it checks nothing on the side of a solve that ends
 * within the time limit anyway.
 *
 * \author Donato Meoli \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 */

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include <cmath>
#include <iostream>

#include "AbstractBlock.h"

#include "FRealObjective.h"

#include "OneVarConstraint.h"

#include "QuadFunction.h"

#include "Solver.h"

/*--------------------------------------------------------------------------*/

using namespace SMSpp_di_unipi_it;

/*--------------------------------------------------------------------------*/

static AbstractBlock * make_problem( bool min )
{
 auto block = new AbstractBlock;

 std::vector< ColVariable > & x = * new std::vector< ColVariable >( 3 );
 block->add_static_variable( x , "x" );

 auto & box = * new std::vector< BoxConstraint >( 3 );
 for( int i = 0 ; i < 3 ; ++i ) {
  box[ i ].set_variable( & x[ i ] );
  box[ i ].set_lhs( -10 );
  box[ i ].set_rhs( 10 );
  }
 block->add_static_constraint( box , "box" );

 QuadFunction::v_coeff_triple diag;
 for( int i = 0 ; i < 3 ; ++i )
  diag.emplace_back( & x[ i ] , - double( i + 1 ) , min ? 1.0 : -1.0 );

 auto obj = new FRealObjective( block ,
                                new QuadFunction( std::move( diag ) ) );
 obj->set_sense( min ? Objective::eMin : Objective::eMax , eNoMod );
 block->set_objective( obj , eNoMod );

 return( block );
 }

/*--------------------------------------------------------------------------*/

int main( void )
{
 const std::string name = "GRBMILPSolver";
 if( ! Solver::has_Solver( name ) ) {
  std::cout << name << " is not in this build: nothing to do" << std::endl;
  return( 77 );
  }

 int failures = 0;

 for( const bool min : { true , false } ) {
  const std::string what = name + ( min ? " (min)" : " (max)" );
  bool ok = true;

  // stopped by the time limit before it has proved anything
  auto s = Solver::new_Solver( name );
  s->set_par( Solver::intLogVerb , 0 );
  s->set_par( Solver::dblMaxTime , 0.0 );
  auto b = make_problem( min );
  b->register_Solver( s );
  const auto st = s->compute();
  if( ( st == Solver::kOK ) || ( st >= Solver::kError ) )
   std::cout << what << ": status " << st << " within the time limit, "
             << "nothing to check" << std::endl;
  else {
   const auto bound = min ? s->get_lb() : s->get_ub();
   ok = min ? ( bound == - Inf< double >() ) : ( bound == Inf< double >() );
   std::cout << what << ": stopped (status " << st << "), "
             << ( min ? "lb" : "ub" ) << " = " << bound
             << ( ok ? "" : " is not infinite" ) << std::endl;
   }
  b->unregister_Solver( s );
  delete s;

  // solved: both bounds are the optimum, 7 / 2 away from 0 on either side
  s = Solver::new_Solver( name );
  s->set_par( Solver::intLogVerb , 0 );
  b->register_Solver( s );
  if( s->compute() != Solver::kOK ) {
   std::cout << what << ": not solved to optimality" << std::endl;
   ok = false;
   }
  else {
   const double opt = min ? -3.5 : 3.5;
   const double lb = s->get_lb() , ub = s->get_ub();
   std::cout << what << ": optimal, bounds [ " << lb << " , " << ub << " ]"
             << std::endl;
   if( ( std::abs( lb - opt ) > 1e-6 ) || ( std::abs( ub - opt ) > 1e-6 ) )
    ok = false;
   }
  b->unregister_Solver( s );
  delete s;

  std::cout << what << ( ok ? " -> OK" : " -> KO" ) << std::endl;
  if( ! ok )
   ++failures;
  }

 return( failures ? 1 : 0 );
 }

/*--------------------------------------------------------------------------*/
/*-------------------- End File test_stopped_bounds.cpp --------------------*/
/*--------------------------------------------------------------------------*/
