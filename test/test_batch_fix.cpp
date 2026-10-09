/*--------------------------------------------------------------------------*/
/*------------------------- File test_batch_fix.cpp ------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Main for testing the :MILPSolver on a batch of changes that names the same
 * column more than once
 *
 * The changes issued on a channel reach a :MILPSolver as one
 * GroupModification, whose fixings and changes of the bounds are executed in
 * one operation of the back-end each [see MILPSolver::change_variables() and
 * MILPSolver::change_bounds()]. A column can be named several times in the
 * same batch, e.g., a ColVariable fixed and unfixed again, and what the
 * model has to end up with is the last state, exactly as when the same
 * changes arrive one at a time.
 *
 * The program is
 *
 *     min  - x - y  s.t.  0 <= x <= 4 ,  0 <= y <= 4
 *
 * the bounds being two BoxConstraint. Each step below issues its changes on
 * one channel, then the Solver computes and the optimal value and the value
 * of x are checked:
 *
 * - x fixed at 1 and unfixed: x = 4;
 *
 * - x fixed at 1, unfixed and fixed again: x = 1;
 *
 * - x (fixed) unfixed and fixed again: x = 1;
 *
 * - x unfixed: x = 4;
 *
 * - the upper bound of x set to 3 and then to 2: x = 2;
 *
 * - the upper bound of x set to 3, x fixed at 1 and unfixed, the upper bound
 *   set to 2.5: x = 2.5;
 *
 * - x fixed at 1 and its upper bound set to 2: x = 1;
 *
 * - x unfixed: x = 2.
 *
 * Every registered :MILPSolver goes through all the steps.
 *
 * \author Donato Meoli \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \copyright &copy; by Donato Meoli
 */
/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include <cmath>
#include <functional>
#include <iostream>
#include <string>
#include <vector>

#include "AbstractBlock.h"

#include "FRealObjective.h"

#include "LinearFunction.h"

#include "MILPSolver.h"

#include "OneVarConstraint.h"

/*--------------------------------------------------------------------------*/
/*-------------------------------- USING -----------------------------------*/
/*--------------------------------------------------------------------------*/

using namespace SMSpp_di_unipi_it;

/*--------------------------------------------------------------------------*/
/*----------------------------- CONSTANTS ----------------------------------*/
/*--------------------------------------------------------------------------*/

/// the :MILPSolver to try, each skipped if it is not in the Solver factory

static const std::vector< std::string > SolverNames =
 { "CPXMILPSolver" , "GRBMILPSolver" , "SCIPMILPSolver" , "HiGHSMILPSolver" };

/// the tolerance of every comparison

static constexpr double Eps = 1e-7;

/*--------------------------------------------------------------------------*/
/*------------------------------ FUNCTIONS ---------------------------------*/
/*--------------------------------------------------------------------------*/

/// goes through the steps of the file comment with the given Solver
/** Returns true if every step gives the right x and optimal value, printing
 * what is wrong otherwise. */

static bool run( const std::string & solver_name )
{
 auto ab = new AbstractBlock();

 auto v = new std::vector< ColVariable >( 2 );
 auto b = new std::vector< BoxConstraint >( 2 );
 ab->add_static_variable( *v , "v" );
 ab->add_static_constraint( *b , "b" );

 ColVariable * x = & (*v)[ 0 ];
 BoxConstraint * bx = & (*b)[ 0 ];
 for( int i = 0 ; i < 2 ; ++i ) {
  (*v)[ i ].set_type( ColVariable::kContinuous , eNoMod );
  (*b)[ i ].set_variable( & (*v)[ i ] , eNoMod );
  (*b)[ i ].set_lhs( 0 , eNoMod );
  (*b)[ i ].set_rhs( 4 , eNoMod );
  }

 auto obj = new FRealObjective( ab , new LinearFunction(
  { std::make_pair( & (*v)[ 0 ] , -1.0 ) ,
    std::make_pair( & (*v)[ 1 ] , -1.0 ) } ) );
 obj->set_sense( Objective::eMin , eNoMod );
 ab->set_objective( obj , eNoMod );

 auto solver = static_cast< CDASolver * >( Solver::new_Solver( solver_name ) );
 solver->set_par( MILPSolver::intLogVerb , 0 );
 ab->register_Solver( solver );

 bool ok = true;

 // issues on one channel what f does, computes and checks that x is xv
 auto step = [ & ]( const char * what , double xv ,
                    const std::function< void( ModParam ) > & f ) {
  auto chnl = ab->open_channel();
  f( Observer::make_par( eModBlck , chnl ) );
  ab->close_channel( chnl );

  std::cout << solver_name << ", " << what << ": ";
  const auto status = solver->compute( false );
  if( status != Solver::kOK ) {
   std::cout << "status " << status << " -> Error" << std::endl;
   ok = false;
   return;
   }

  solver->get_var_solution();
  bool sok = true;
  if( std::abs( x->get_value() - xv ) > Eps ) {
   std::cout << "x = " << x->get_value() << " != " << xv << " ";
   sok = false;
   }
  if( std::abs( solver->get_var_value() + xv + 4 ) > Eps ) {
   std::cout << "value " << solver->get_var_value() << " != "
             << - xv - 4 << " ";
   sok = false;
   }
  std::cout << ( sok ? "-> OK" : "-> Error" ) << std::endl;
  ok = ok && sok;
  };

 if( solver->compute( false ) != Solver::kOK ) {
  std::cout << solver_name << ": first solve failed -> Error" << std::endl;
  ok = false;
  }
 else {
  step( "fix , unfix" , 4 , [ & ]( ModParam p ) {
   x->set_value( 1 );
   x->is_fixed( true , p );
   x->is_fixed( false , p );
   } );

  step( "fix , unfix , fix" , 1 , [ & ]( ModParam p ) {
   x->set_value( 1 );
   x->is_fixed( true , p );
   x->is_fixed( false , p );
   x->is_fixed( true , p );
   } );

  step( "unfix , fix" , 1 , [ & ]( ModParam p ) {
   x->is_fixed( false , p );
   x->is_fixed( true , p );
   } );

  step( "unfix" , 4 , [ & ]( ModParam p ) {
   x->is_fixed( false , p );
   } );

  step( "ub 3 , ub 2" , 2 , [ & ]( ModParam p ) {
   bx->set_rhs( 3 , p );
   bx->set_rhs( 2 , p );
   } );

  step( "ub 3 , fix , unfix , ub 2.5" , 2.5 , [ & ]( ModParam p ) {
   bx->set_rhs( 3 , p );
   x->set_value( 1 );
   x->is_fixed( true , p );
   x->is_fixed( false , p );
   bx->set_rhs( 2.5 , p );
   } );

  step( "fix , ub 2" , 1 , [ & ]( ModParam p ) {
   x->set_value( 1 );
   x->is_fixed( true , p );
   bx->set_rhs( 2 , p );
   } );

  step( "unfix" , 2 , [ & ]( ModParam p ) {
   x->is_fixed( false , p );
   } );
  }

 ab->unregister_Solver( solver );
 delete solver;
 delete ab;

 return( ok );
 }

/*--------------------------------------------------------------------------*/
/*-------------------------------- main() ----------------------------------*/
/*--------------------------------------------------------------------------*/

int main( void )
{
 bool all_passed = true;
 bool any_solver = false;

 for( const auto & name : SolverNames ) {
  // skip the :MILPSolver that are not in the build
  if( ! Solver::has_Solver( name ) )
   continue;
  any_solver = true;

  if( ! run( name ) )
   all_passed = false;
  }

 if( ! any_solver ) {
  std::cout << "no :MILPSolver in this build, nothing to check" << std::endl;
  return( 0 );
  }

 if( all_passed )
  std::cout << "All tests passed!!" << std::endl;
 else
  std::cout << "Shit happened!!" << std::endl;

 return( all_passed ? 0 : 1 );
 }

/*--------------------------------------------------------------------------*/
/*------------------------- File test_batch_fix.cpp ------------------------*/
/*--------------------------------------------------------------------------*/
