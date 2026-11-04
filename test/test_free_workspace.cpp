/*--------------------------------------------------------------------------*/
/*---------------------- File test_free_workspace.cpp ----------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Checks that every :MILPSolver of LPPar-all.txt that the build has gives
 * the same solution as the others across warm re-solves, after each of
 * which a Solver may free the memory that its back-end keeps and that the
 * next solve does not need (the workspace of GUROBI, for instance). The LP
 * is a dispatch over 4 instants with demand d_t:
 *
 *     min   sum_t ( c g_t + 1000 y_t )
 *     s.t.  g_t + y_t = d_t                         (balance)
 *           0 <= g_t <= u_t , 0 <= y_t <= 1000
 *
 * with d = ( 40 , 40 , 160 , 160 ), u = 100 and c = 1, whose optimum, primal
 * and dual, is unique after each of these changes, which are those a
 * decomposition makes to a subproblem and the edge cases of a solution
 * kept by the Solver:
 *
 * - none, the first solve;
 * - the demand, i.e., the sides of the balance rows;
 * - the bound of g_2, i.e., a bound;
 * - the cost of g, i.e., the Objective;
 * - a demand above u + 1000 at t = 3, which makes the LP infeasible;
 * - the demand back, a solve after an infeasible one;
 * - none again, a solve of a model with no change since the last one.
 *
 * After each one, the status, the value, the value of every Variable and the
 * dual of every row and of every bound of each Solver have to be those of
 * the first one. The test is skipped (exit status 77) if the build has fewer
 * than two of them.
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

#include "test_common.h"

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

 // g_t , y_t in this order
 pb.vars = new std::vector< ColVariable >( 2 * T );
 auto g = [ & ]( int t ) { return( & ( * pb.vars )[ t ] ); };
 auto y = [ & ]( int t ) { return( & ( * pb.vars )[ T + t ] ); };
 b->add_static_variable( * pb.vars , "x" );

 // the bounds of g_t , y_t in this order
 pb.bounds = new std::vector< BoxConstraint >( 2 * T );
 for( int t = 0 ; t < T ; ++t ) {
  ( * pb.bounds )[ t ].set_variable( g( t ) );
  ( * pb.bounds )[ t ].set_lhs( 0 );
  ( * pb.bounds )[ t ].set_rhs( 100 );
  ( * pb.bounds )[ T + t ].set_variable( y( t ) );
  ( * pb.bounds )[ T + t ].set_lhs( 0 );
  ( * pb.bounds )[ T + t ].set_rhs( 1000 );
  }
 b->add_static_constraint( * pb.bounds , "bounds" );

 // the balance rows
 const double d[ T ] = { 40 , 40 , 160 , 160 };
 pb.rows = new std::vector< FRowConstraint >( T );
 for( int t = 0 ; t < T ; ++t ) {
  auto & bal = ( * pb.rows )[ t ];
  bal.set_function( new LinearFunction( LinearFunction::v_coeff_pair{
   { g( t ) , 1 } , { y( t ) , 1 } } ) );
  bal.set_both( d[ t ] );
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
  if( std::abs( a[ i ] - b[ i ] ) > 1e-5 * ( 1 + std::abs( a[ i ] ) ) )
   return( false );
 return( true );
 }

/*--------------------------------------------------------------------------*/

int main( void )
{
 const std::string fn = "LPPar-all.txt";
 auto bsc = dynamic_cast< BlockSolverConfig * >(
                                          Configuration::deserialize( fn ) );
 Block::Index there = 0;
 if( bsc )
  for( Block::Index i = 0 ; i < bsc->num_ComputeConfig() ; ++i )
   if( Solver::has_Solver( bsc->get_SolverName( i ) ) )
    ++there;
 keep_available_Solvers( bsc , fn , std::max( there , Block::Index( 2 ) ) );

 auto pb = make_problem();
 bsc->apply( pb.block );
 bsc->clear();  // kept, clear()-ed, for the final cleanup

 auto & solvers = pb.block->get_registered_solvers();
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
  Outcome ref;
  bool first = true;
  for( auto s : solvers ) {
   const auto o = solve( s , pb );
   if( first ) {
    ref = o;
    first = false;
    }
   const bool ok = ( o.status == ref.status ) &&
                   ( ( o.status != Solver::kOK ) ||
                     ( same( { o.value } , { ref.value } ) &&
                       same( o.x , ref.x ) && same( o.duals , ref.duals ) ) );
   std::cout << what << ", " << s->classname() << ": status " << o.status
             << ", value " << o.value << ( ok ? " -> OK" : " -> KO" )
             << std::endl;
   if( ( ! ok ) && ( o.status == Solver::kOK ) ) {
    for( std::size_t i = 0 ; i < o.x.size() ; ++i )
     if( ! same( { o.x[ i ] } , { ref.x[ i ] } ) )
      std::cout << "  x[ " << i << " ] = " << o.x[ i ] << " instead of "
                << ref.x[ i ] << std::endl;
    for( std::size_t i = 0 ; i < o.duals.size() ; ++i )
     if( ! same( { o.duals[ i ] } , { ref.duals[ i ] } ) )
      std::cout << "  dual[ " << i << " ] = " << o.duals[ i ]
                << " instead of " << ref.duals[ i ] << std::endl;
    }
   if( ! ok )
    ++failures;
   }
  }

 bsc->apply( pb.block );
 delete bsc;
 delete pb.block;

 return( failures ? 1 : 0 );
 }

/*--------------------------------------------------------------------------*/
/*------------------- End File test_free_workspace.cpp ---------------------*/
/*--------------------------------------------------------------------------*/
