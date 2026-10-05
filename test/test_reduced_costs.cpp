/*--------------------------------------------------------------------------*/
/*---------------------- File test_reduced_costs.cpp -----------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Main for testing that a :MILPSolver gives the reduced costs to the bounds
 * of the right column when the model has dynamic columns
 *
 * A :MILPSolver builds first the columns of all the static Variable, Block by
 * Block in breadth-first order, and then those of the dynamic ones, in the
 * order they arrive; the reduced cost of a column becomes the dual of the
 * OneVarConstraint bounding its Variable. The program is written so that a
 * dynamic Variable belongs to a Block that comes before the one of the static
 * Variable: the root Block has the dynamic y, and its only sub-Block has the
 * static x0 and x1, so that the columns are ( x0 , x1 , y ) while walking the
 * Variable of each Block in turn gives ( y , x0 , x1 ). The program is
 *
 *     min  3 y - x0 - 2 x1  s.t.  x0 + x1 <= 1.5 ,
 *                                 0 <= x0 <= 1 ,  0 <= x1 <= 1 ,  y >= 1
 *
 * the row and the bounds of x0 and x1 in the sub-Block, the bound of y in the
 * root Block, with the objective split in the same way. Its optimum is 0.5,
 * at ( x0 , x1 , y ) = ( 0.5 , 1 , 1 ): the dual of the row is 1, the reduced
 * costs are 0 for x0, which is strictly inside its bounds, -1 for x1 and 3 for
 * y, and these are what the duals of the three bounds must be (up to the sign
 * convention of MILPSolver::write_dual_solution(), which is why their absolute
 * value is compared).
 *
 * Every registered :MILPSolver among CPXMILPSolver, GRBMILPSolver and
 * HiGHSMILPSolver solves the program, and the optimal value and the duals are
 * checked. SCIPMILPSolver is left out because it gives the duals only when
 * intComputeDuals is set, its presolve otherwise solving this program
 * without an LP; the reduced costs are handed to the bounds by
 * MILPSolver::write_dual_solution(), which is the same for all of them.
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
#include <iostream>
#include <list>
#include <vector>

#include "AbstractBlock.h"

#include "BlockSolverConfig.h"

#include "FRealObjective.h"

#include "FRowConstraint.h"

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
 { "CPXMILPSolver" , "GRBMILPSolver" , "HiGHSMILPSolver" };

/// the tolerance of every comparison

static constexpr double Eps = 1e-7;

/*--------------------------------------------------------------------------*/
/*------------------------------ FUNCTIONS ---------------------------------*/
/*--------------------------------------------------------------------------*/

/// solves the program of the file comment with the given Solver
/** Returns true if the optimal value and the duals are right, printing what
 * is wrong otherwise. */

static bool solve( const std::string & solver_name )
{
 auto ab = new AbstractBlock();
 auto sb = new AbstractBlock( ab );

 // the root Block: the dynamic y, its bound and its part of the objective
 auto y = new std::list< ColVariable >( 1 );
 ab->add_dynamic_variable( *y , "y" );
 auto & yv = y->front();

 auto yb = new std::vector< BoxConstraint >( 1 );
 (*yb)[ 0 ].set_variable( & yv , eNoMod );
 (*yb)[ 0 ].set_lhs( 1 , eNoMod );
 (*yb)[ 0 ].set_rhs( Inf< double >() , eNoMod );
 ab->add_static_constraint( *yb , "y_bound" );

 auto aobj = new FRealObjective( ab , new LinearFunction(
                                          { std::make_pair( & yv , 3.0 ) } ) );
 aobj->set_sense( Objective::eMin , eNoMod );
 ab->set_objective( aobj , eNoMod );

 // the sub-Block: the static x0 and x1, their bounds, the row and their
 // part of the objective
 auto x = new std::vector< ColVariable >( 2 );
 sb->add_static_variable( *x , "x" );

 auto xb = new std::vector< BoxConstraint >( 2 );
 for( int i = 0 ; i < 2 ; ++i ) {
  (*xb)[ i ].set_variable( & (*x)[ i ] , eNoMod );
  (*xb)[ i ].set_lhs( 0 , eNoMod );
  (*xb)[ i ].set_rhs( 1 , eNoMod );
  }
 sb->add_static_constraint( *xb , "x_bound" );

 auto r = new std::vector< FRowConstraint >( 1 );
 (*r)[ 0 ].set_function( new LinearFunction(
  { std::make_pair( & (*x)[ 0 ] , 1.0 ) ,
    std::make_pair( & (*x)[ 1 ] , 1.0 ) } ) , eNoMod );
 (*r)[ 0 ].set_lhs( -Inf< double >() , eNoMod );
 (*r)[ 0 ].set_rhs( 1.5 , eNoMod );
 sb->add_static_constraint( *r , "r" );

 auto sobj = new FRealObjective( sb , new LinearFunction(
  { std::make_pair( & (*x)[ 0 ] , -1.0 ) ,
    std::make_pair( & (*x)[ 1 ] , -2.0 ) } ) );
 sobj->set_sense( Objective::eMin , eNoMod );
 sb->set_objective( sobj , eNoMod );

 ab->add_nested_Block( sb );

 auto bsc = new BlockSolverConfig( 1 );
 bsc->add_ComputeConfig( std::string( solver_name ) , nullptr );
 bsc->apply( ab );
 auto solver = static_cast< CDASolver * >(
                                    ab->get_registered_solvers().front() );
 solver->set_par( MILPSolver::intLogVerb , 0 );

 const auto status = solver->compute( false );
 bool ok = true;

 std::cout << solver_name << ": ";
 if( status != Solver::kOK ) {
  std::cout << "status " << status << " -> Error" << std::endl;
  ok = false;
  }
 else {
  solver->get_var_solution();
  solver->get_dual_solution();

  if( std::abs( solver->get_var_value() - 0.5 ) > Eps ) {
   std::cout << "value " << solver->get_var_value() << " != 0.5 ";
   ok = false;
   }

  auto check = [ & ok ]( const char * name , double dual , double right ) {
   if( std::abs( std::abs( dual ) - right ) > Eps ) {
    std::cout << "dual of " << name << " " << dual << " != " << right << " ";
    ok = false;
    }
   };

  check( "the row" , (*r)[ 0 ].get_dual() , 1 );
  check( "the bound of x0" , (*xb)[ 0 ].get_dual() , 0 );
  check( "the bound of x1" , (*xb)[ 1 ].get_dual() , 1 );
  check( "the bound of y" , (*yb)[ 0 ].get_dual() , 3 );

  std::cout << ( ok ? "-> OK" : "-> Error" ) << std::endl;
  }

 bsc->clear();
 bsc->apply( ab );
 delete bsc;
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

  if( ! solve( name ) )
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
/*---------------------- End File test_reduced_costs.cpp -------------------*/
/*--------------------------------------------------------------------------*/
