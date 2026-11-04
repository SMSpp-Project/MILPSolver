/*--------------------------------------------------------------------------*/
/*------------------------- File test_qcp_dual.cpp -------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Checks the duals that each :MILPSolver of the build able to solve a
 * convex problem with quadratic rows writes in them, and in the linear rows
 * that follow them. The problem is
 *
 *     min  - v + x + y - 2 u - z - 3 t
 *     s.t. v <= 4                        (r0 , linear)
 *          x^2 + y^2 <= 2                (r1 , quadratic , no linear part)
 *          u <= 5                        (r2 , linear)
 *          z + w^2 <= 1                  (r3 , quadratic with a linear part)
 *          t <= 3                        (r4 , linear)
 *
 * with every Variable in [ -10 , 10 ]. Its optimum is v = 4, x = y = -1,
 * u = 5, z = 1, w = 0 and t = 3, of value -26, every row is active there and
 * its multipliers are 1, 1 / 2, 2, 1 and 3, those of the two quadratic rows
 * being what makes the gradient of the objective vanish together with the
 * gradient of the row, e.g., 1 + \mu_1 2 x = 0 at x = -1. The duals have to
 * be these, all with the sign that the Solver gives to r0, whose convention
 * test_dual checks: a dual of r1 or r3 with the wrong value, or a dual of r2
 * or r4 taken from another row (the linear rows of the back-end not being
 * those of the Block when a quadratic row has no linear part) is an error.
 *
 * The :MILPSolver are those of LPPar-qcp.txt that the build has, i.e., the
 * ones that take quadratic rows; one that gives no duals of them is said so
 * and skipped, except CPXMILPSolver and GRBMILPSolver, which have to solve
 * the problem. The test fails if one of them does not, or if any gives a
 * wrong dual.
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
#include <set>
#include <string>
#include <vector>

#include "AbstractBlock.h"

#include "CDASolver.h"

#include "DQuadFunction.h"

#include "FRealObjective.h"

#include "FRowConstraint.h"

#include "LinearFunction.h"

#include "OneVarConstraint.h"

#include "test_common.h"

/*--------------------------------------------------------------------------*/

using namespace SMSpp_di_unipi_it;

/*--------------------------------------------------------------------------*/

/// the Variable v , x , y , u , z , w , t and the rows r0 , ... , r4
struct Problem {
 AbstractBlock * block;
 std::vector< ColVariable > * vars;
 std::vector< FRowConstraint > * rows;
 };

static const double XStar[ 7 ] = { 4 , -1 , -1 , 5 , 1 , 0 , 3 };
static const double MuStar[ 5 ] = { 1 , 0.5 , 2 , 1 , 3 };
static constexpr double ValStar = -26;

/*--------------------------------------------------------------------------*/

static Problem make_problem( void )
{
 Problem pb;
 pb.block = new AbstractBlock;
 pb.vars = new std::vector< ColVariable >( 7 );
 auto & x = * pb.vars;
 for( auto & xi : x )
  xi.set_type( ColVariable::kContinuous );
 pb.block->add_static_variable( x , "x" );

 auto bounds = new std::vector< BoxConstraint >( 7 );
 for( int i = 0 ; i < 7 ; ++i ) {
  ( * bounds )[ i ].set_variable( & x[ i ] );
  ( * bounds )[ i ].set_lhs( -10 );
  ( * bounds )[ i ].set_rhs( 10 );
  }
 pb.block->add_static_constraint( * bounds , "bounds" );

 pb.rows = new std::vector< FRowConstraint >( 5 );
 auto & r = * pb.rows;
 auto lin = []( ColVariable * v ) {
  return( new LinearFunction( LinearFunction::v_coeff_pair{ { v , 1.0 } } ) );
  };
 r[ 0 ].set_function( lin( & x[ 0 ] ) );
 r[ 1 ].set_function( new DQuadFunction( DQuadFunction::v_coeff_triple{
		     { & x[ 1 ] , 0.0 , 1.0 } , { & x[ 2 ] , 0.0 , 1.0 } } ) );
 r[ 2 ].set_function( lin( & x[ 3 ] ) );
 r[ 3 ].set_function( new DQuadFunction( DQuadFunction::v_coeff_triple{
		     { & x[ 4 ] , 1.0 , 0.0 } , { & x[ 5 ] , 0.0 , 1.0 } } ) );
 r[ 4 ].set_function( lin( & x[ 6 ] ) );
 const double rhs[ 5 ] = { 4 , 2 , 5 , 1 , 3 };
 for( int i = 0 ; i < 5 ; ++i ) {
  r[ i ].set_lhs( - Inf< double >() );
  r[ i ].set_rhs( rhs[ i ] );
  }
 pb.block->add_static_constraint( r , "rows" );

 auto of = new FRealObjective( pb.block , new LinearFunction(
  LinearFunction::v_coeff_pair{ { & x[ 0 ] , -1.0 } , { & x[ 1 ] , 1.0 } ,
                                { & x[ 2 ] , 1.0 } , { & x[ 3 ] , -2.0 } ,
                                { & x[ 4 ] , -1.0 } , { & x[ 6 ] , -3.0 } }
  ) );
 of->set_sense( Objective::eMin );
 pb.block->set_objective( of );

 return( pb );
 }

/*--------------------------------------------------------------------------*/

static bool close( double a , double b , double eps )
{
 return( std::abs( a - b ) <= eps * ( 1 + std::abs( b ) ) );
 }

/*--------------------------------------------------------------------------*/

int main( void )
{
 auto bsc = dynamic_cast< BlockSolverConfig * >(
		    Configuration::deserialize( "LPPar-qcp.txt" ) );
 Block::Index there = 0;
 if( bsc )
  for( Block::Index i = 0 ; i < bsc->num_ComputeConfig() ; ++i )
   if( Solver::has_Solver( bsc->get_SolverName( i ) ) )
    ++there;
 if( bsc && ( ! there ) ) {
  std::cout << "no :MILPSolver in this build" << std::endl;
  return( 1 );
  }
 keep_available_Solvers( bsc , "LPPar-qcp.txt" , there );

 std::vector< std::string > names;
 for( Block::Index i = 0 ; i < bsc->num_ComputeConfig() ; ++i )
  names.push_back( bsc->get_SolverName( i ) );

 const std::set< std::string > must = { "CPXMILPSolver" , "GRBMILPSolver" };

 auto pb = make_problem();
 bsc->apply( pb.block );
 bsc->clear();  // kept, clear()-ed, for the final cleanup

 int failures = 0;
 std::size_t i = 0;
 for( auto s : pb.block->get_registered_solvers() ) {
  const auto & name = names[ i++ ];
  auto solver = static_cast< CDASolver * >( s );

  int status;
  try {
   status = solver->compute();
   }
  catch( std::exception & ) {
   status = Solver::kError;
   }

  if( ( status != Solver::kOK ) || ( ! solver->has_dual_solution() ) ) {
   if( must.count( name ) ) {
    std::cout << name << ": no optimal solution with duals -> KO"
              << std::endl;
    ++failures;
    }
   else
    std::cout << name << ": does not solve it, skipped" << std::endl;
   continue;
   }

  solver->get_var_solution();
  solver->get_dual_solution();

  bool ok = close( solver->get_var_value() , ValStar , 1e-5 );
  for( int j = 0 ; j < 7 ; ++j )
   if( ! close( ( * pb.vars )[ j ].get_value() , XStar[ j ] , 1e-4 ) )
    ok = false;

  const double sign = ( ( * pb.rows )[ 0 ].get_dual() < 0 ) ? -1 : 1;
  std::cout << name << ": duals";
  for( int j = 0 ; j < 5 ; ++j ) {
   const double d = ( * pb.rows )[ j ].get_dual();
   std::cout << " " << d;
   if( ! close( d , sign * MuStar[ j ] , 1e-4 ) )
    ok = false;
   }
  std::cout << " (" << sign << " times 1 0.5 2 1 3)"
            << ( ok ? " -> OK" : " -> KO" ) << std::endl;
  if( ! ok )
   ++failures;
  }

 bsc->apply( pb.block );
 delete bsc;
 delete pb.block;

 return( failures ? 1 : 0 );
 }

/*--------------------------------------------------------------------------*/
/*------------------------ End File test_qcp_dual.cpp ----------------------*/
/*--------------------------------------------------------------------------*/
