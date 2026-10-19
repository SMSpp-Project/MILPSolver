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
 * is solved with dblMaxTime = 0, and the bound on the side of the optimum that
 * the Solver has not proved (get_lb() when minimizing, get_ub() when
 * maximizing) has to be either infinite or a valid bound, i.e., not beyond the
 * optimum, and in no case a large finite number that a caller would take for a
 * bound; then the time limit is lifted (dblMaxTime = INF) on the same Solver,
 * which has to solve the problem to optimality, with both bounds equal to the
 * optimal value. Each :MILPSolver among CPXMILPSolver, GRBMILPSolver,
 * MOSEKMILPSolver, HiGHSMILPSolver, SCIPMILPSolver, DAQPMILPSolver,
 * QPOASESMILPSolver, ClarabelMILPSolver, ProxQPMILPSolver, PIQPMILPSolver,
 * OSQPMILPSolver, QPMADMILPSolver and QuadProgMILPSolver that is not in the
 * factory is skipped, and the test is skipped (exit status 77) if none of
 * them is; nothing is checked on the first solve of a Solver that ends
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
 int failures = 0 , tried = 0;

 for( const std::string name : { "CPXMILPSolver" , "GRBMILPSolver" ,
                                 "MOSEKMILPSolver" ,
                                 "HiGHSMILPSolver" , "SCIPMILPSolver" ,
                                 "DAQPMILPSolver" , "QPOASESMILPSolver" ,
                                 "ClarabelMILPSolver" , "ProxQPMILPSolver" ,
                                 "PIQPMILPSolver" , "OSQPMILPSolver" ,
                                 "QPMADMILPSolver" ,
                                 "QuadProgMILPSolver" } ) {
  if( ! Solver::has_Solver( name ) ) {
   std::cout << name << " is not in this build, skipped" << std::endl;
   continue;
   }
  ++tried;

  for( const bool min : { true , false } ) {
   const std::string what = name + ( min ? " (min)" : " (max)" );
   // the optimum is 7 / 2 away from 0, on either side
   const double opt = min ? -3.5 : 3.5;
   bool ok = true;

   auto s = Solver::new_Solver( name );
   s->set_par( Solver::intLogVerb , 0 );
   auto b = make_problem( min );
   b->register_Solver( s );

   // stopped by the time limit before it has proved anything
   s->set_par( Solver::dblMaxTime , 0.0 );
   const auto st = s->compute();
   if( ( st == Solver::kOK ) || ( st >= Solver::kError ) )
    std::cout << what << ": status " << st << " within the time limit, "
              << "nothing to check" << std::endl;
   else {
    const auto bound = min ? s->get_lb() : s->get_ub();
    const bool infinite = min ? ( bound == - Inf< double >() )
                              : ( bound == Inf< double >() );
    const bool valid = ( std::abs( bound ) < 1e+15 ) &&
                       ( min ? ( bound <= opt + 1e-6 )
                             : ( bound >= opt - 1e-6 ) );
    ok = infinite || valid;
    std::cout << what << ": stopped (status " << st << "), "
              << ( min ? "lb" : "ub" ) << " = " << bound
              << ( ok ? "" : " is no bound" ) << std::endl;
    }

   // the time limit lifted on the same Solver
   s->set_par( Solver::dblMaxTime , Inf< double >() );
   if( s->compute() != Solver::kOK ) {
    std::cout << what << ": not solved to optimality once the time limit "
              << "is lifted" << std::endl;
    ok = false;
    }
   else {
    const double lb = s->get_lb() , ub = s->get_ub();
    std::cout << what << ": optimal once the time limit is lifted, bounds [ "
              << lb << " , " << ub << " ]" << std::endl;
    if( ( std::abs( lb - opt ) > 1e-6 ) || ( std::abs( ub - opt ) > 1e-6 ) )
     ok = false;
    }

   b->unregister_Solver( s );
   delete s;

   std::cout << what << ( ok ? " -> OK" : " -> KO" ) << std::endl;
   if( ! ok )
    ++failures;
   }
  }

 if( ! tried ) {
  std::cout << "no :MILPSolver in this build: nothing to do" << std::endl;
  return( 77 );
  }

 return( failures ? 1 : 0 );
 }

/*--------------------------------------------------------------------------*/
/*-------------------- End File test_stopped_bounds.cpp --------------------*/
/*--------------------------------------------------------------------------*/
