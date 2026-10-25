/*--------------------------------------------------------------------------*/
/*---------------------- File test_free_workspace.cpp ----------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Checks that a GRBMILPSolver freeing the workspace of GUROBI after each
 * solve (intFreeWorkspace 1) gives what one keeping it gives, across warm
 * re-solves. The LP is a dispatch over 4 instants with demand d_t:
 *
 *     min   sum_t ( c g_t + 1000 y_t )
 *     s.t.  g_t + y_t + p_t = d_t                   (balance)
 *           v_t - v_{t-1} + p_t = 0 , v_{-1} = 0    (storage)
 *           0 <= g_t <= u , 0 <= y_t <= 1000 , 0 <= v_t <= 50
 *
 * with d = ( 40 , 40 , 160 , 160 ), u = 100 and c = 1. The two Solver of
 * LPPar-ws.txt solve it after each of these changes, which are those a
 * decomposition makes to a subproblem and the edge cases of the kept
 * solution:
 *
 * - none, the first solve;
 * - the demand, i.e., the sides of the balance rows;
 * - the bound of g_2, i.e., a bound;
 * - the cost of g, i.e., the Objective;
 * - a demand above u + 1000 at t = 3, which makes the LP infeasible, whose
 *   solve the workspace is not freed after;
 * - the demand back, a solve after an infeasible one;
 * - none again, a solve of a model with no change since the kept solution.
 *
 * After each one, the status, the value, the value of every Variable and the
 * dual of every row and of every bound of the two have to be the same.
 *
 * The test is skipped (exit status 77) if GRBMILPSolver is not in the build.
 *
 * \author Donato Meoli \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 */

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include <cmath>
#include <functional>
#include <iostream>
#include <tuple>

#include "AbstractBlock.h"

#include "BlockSolverConfig.h"

#include "CDASolver.h"

#include "FRealObjective.h"

#include "FRowConstraint.h"

#include "LinearFunction.h"

#include "OneVarConstraint.h"

/*--------------------------------------------------------------------------*/

using namespace SMSpp_di_unipi_it;

/*--------------------------------------------------------------------------*/

static constexpr int T = 4;

/// the LP and the parts of it that the test changes or reads
struct Problem {
 AbstractBlock * block;
 std::vector< ColVariable > * vars;
 std::vector< FRowConstraint > * rows;
 std::vector< BoxConstraint > * bounds;
 LinearFunction * obj;
 };

/*--------------------------------------------------------------------------*/

static Problem make_problem( void )
{
 Problem pb;
 pb.block = new AbstractBlock;
 auto b = pb.block;

 // g_t , y_t , p_t , v_t in this order
 pb.vars = new std::vector< ColVariable >( 4 * T );
 auto g = [ & ]( int t ) { return( & ( * pb.vars )[ t ] ); };
 auto y = [ & ]( int t ) { return( & ( * pb.vars )[ T + t ] ); };
 auto p = [ & ]( int t ) { return( & ( * pb.vars )[ 2 * T + t ] ); };
 auto v = [ & ]( int t ) { return( & ( * pb.vars )[ 3 * T + t ] ); };
 b->add_static_variable( * pb.vars , "x" );

 // the bounds of g_t , y_t , v_t in this order
 pb.bounds = new std::vector< BoxConstraint >( 3 * T );
 for( int t = 0 ; t < T ; ++t ) {
  ( * pb.bounds )[ t ].set_variable( g( t ) );
  ( * pb.bounds )[ t ].set_lhs( 0 );
  ( * pb.bounds )[ t ].set_rhs( 100 );
  ( * pb.bounds )[ T + t ].set_variable( y( t ) );
  ( * pb.bounds )[ T + t ].set_lhs( 0 );
  ( * pb.bounds )[ T + t ].set_rhs( 1000 );
  ( * pb.bounds )[ 2 * T + t ].set_variable( v( t ) );
  ( * pb.bounds )[ 2 * T + t ].set_lhs( 0 );
  ( * pb.bounds )[ 2 * T + t ].set_rhs( 50 );
  }
 b->add_static_constraint( * pb.bounds , "bounds" );

 // the balance rows , then the storage ones
 const double d[ T ] = { 40 , 40 , 160 , 160 };
 pb.rows = new std::vector< FRowConstraint >( 2 * T );
 for( int t = 0 ; t < T ; ++t ) {
  auto & bal = ( * pb.rows )[ t ];
  bal.set_function( new LinearFunction( LinearFunction::v_coeff_pair{
   { g( t ) , 1 } , { y( t ) , 1 } , { p( t ) , 1 } } ) );
  bal.set_both( d[ t ] );

  LinearFunction::v_coeff_pair sv = { { v( t ) , 1 } , { p( t ) , 1 } };
  if( t > 0 )
   sv.push_back( { v( t - 1 ) , -1 } );
  ( * pb.rows )[ T + t ].set_function(
                                    new LinearFunction( std::move( sv ) ) );
  ( * pb.rows )[ T + t ].set_both( 0 );
  }
 b->add_static_constraint( * pb.rows , "rows" );

 LinearFunction::v_coeff_pair obj;
 for( int t = 0 ; t < T ; ++t ) {
  obj.push_back( { g( t ) , 1 } );
  obj.push_back( { y( t ) , 1000 } );
  }
 pb.obj = new LinearFunction( std::move( obj ) );
 auto o = new FRealObjective( b , pb.obj );
 o->set_sense( Objective::eMin , eNoMod );
 b->set_objective( o , eNoMod );

 return( pb );
 }

/*--------------------------------------------------------------------------*/
/// what a Solver says after a solve: status, value, primal and dual values

struct Outcome {
 int status;
 double value = std::nan( "" );
 std::vector< double > x;
 std::vector< double > duals;
 };

static Outcome solve( Solver * s , const Problem & pb )
{
 Outcome out;
 out.status = s->compute();
 if( out.status != Solver::kOK )
  return( out );

 out.value = s->get_var_value();
 s->get_var_solution();
 for( auto & x : * pb.vars )
  out.x.push_back( x.get_value() );

 static_cast< CDASolver * >( s )->get_dual_solution();
 for( auto & r : * pb.rows )
  out.duals.push_back( r.get_dual() );
 for( auto & r : * pb.bounds )
  out.duals.push_back( r.get_dual() );

 return( out );
 }

/*--------------------------------------------------------------------------*/

static bool same( const std::vector< double > & a ,
                  const std::vector< double > & b )
{
 if( a.size() != b.size() )
  return( false );
 for( std::size_t i = 0 ; i < a.size() ; ++i )
  if( std::abs( a[ i ] - b[ i ] ) > 1e-6 * ( 1 + std::abs( a[ i ] ) ) )
   return( false );
 return( true );
 }

/*--------------------------------------------------------------------------*/

int main( void )
{
 if( ! Solver::has_Solver( "GRBMILPSolver" ) ) {
  std::cout << "GRBMILPSolver is not in this build: nothing to do"
            << std::endl;
  return( 77 );
  }

 auto pb = make_problem();

 auto bsc = dynamic_cast< BlockSolverConfig * >(
                          Configuration::deserialize( "LPPar-ws.txt" ) );
 if( ! bsc ) {
  std::cerr << "Error: LPPar-ws.txt is not a BlockSolverConfig" << std::endl;
  return( 1 );
  }
 bsc->apply( pb.block );
 bsc->clear();  // kept, clear()-ed, for the final cleanup

 auto & solvers = pb.block->get_registered_solvers();
 auto keep = solvers.front();
 auto frees = solvers.back();

 auto & rows = * pb.rows;
 auto & bounds = * pb.bounds;

 // the change and its description; the first and the last change nothing
 const std::vector< std::tuple< std::string , std::function< void() > > >
  steps = {
   { "first solve" , [](){} } ,
   { "the demand" , [ & ]() {
     rows[ 0 ].set_both( 60 );
     rows[ 2 ].set_both( 120 );
     } } ,
   { "the bound of g_2" , [ & ]() { bounds[ 2 ].set_rhs( 80 ); } } ,
   { "the cost of g" , [ & ]() {
     for( int t = 0 ; t < T ; ++t )
      pb.obj->modify_coefficient(
                        pb.obj->is_active( & ( * pb.vars )[ t ] ) , 3 );
     } } ,
   { "an infeasible demand" , [ & ]() { rows[ 3 ].set_both( 2000 ); } } ,
   { "the demand back" , [ & ]() { rows[ 3 ].set_both( 160 ); } } ,
   { "no change" , [](){} } };

 int failures = 0;
 for( const auto & [ what , change ] : steps ) {
  change();
  const auto a = solve( keep , pb );
  const auto b = solve( frees , pb );
  const bool ok = ( a.status == b.status ) &&
                  ( ( a.status != Solver::kOK ) ||
                    ( same( { a.value } , { b.value } ) &&
                      same( a.x , b.x ) && same( a.duals , b.duals ) ) );
  std::cout << what << ": status " << a.status << " and " << b.status
            << ", value " << a.value << " and " << b.value
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
/*------------------- End File test_free_workspace.cpp ---------------------*/
/*--------------------------------------------------------------------------*/
