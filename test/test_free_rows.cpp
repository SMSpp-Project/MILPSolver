/*--------------------------------------------------------------------------*/
/*------------------------ File test_free_rows.cpp -------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Checks the rows with no finite side, and the relaxed ones, of a
 * :MILPSolver across warm re-solves. The LP is a dispatch over 4 instants
 * with demand d = ( 40 , 40 , 160 , 160 ):
 *
 *     min   1000 sum_t y_t
 *     s.t.  s g_t + y_t + p_t = d_t                 (balance)
 *           v_t - v_{t-1} + p_t = 0 , v_{-1} = 0    (storage)
 *           - inf <= p_t <= + inf                   (free rows, two per t)
 *           0 <= g_t <= 100 , 0 <= y_t <= 1000 , 0 <= v_t <= 50
 *
 * with s = 1, whose optimum is 70000 (50 of the surplus 120 of the first
 * two instants is stored); the free rows are those of a battery with no
 * limit on its power. The same Solver then solves it again after each of
 * these changes, and its value has to be the expected one and that of a
 * Solver of the same type loaded from scratch:
 *
 * - s = 0.5, a change of the coefficients: 200000;
 * - the storage rows relaxed: 0;
 * - the storage rows enforced again: 200000;
 * - p_2 <= 5 and p_3 <= 5 on the first free row of each: 210000;
 * - the free rows free again: 200000.
 *
 * On the first change, a free row written as p_t <= 1e100 (or <= 1e33)
 * makes GUROBI 13 end the warm re-solve in numerical trouble (status 12),
 * while <= 1e32 does not.
 *
 * Each :MILPSolver among CPXMILPSolver, GRBMILPSolver, HiGHSMILPSolver,
 * SCIPMILPSolver, DAQPMILPSolver, QPOASESMILPSolver and ClarabelMILPSolver
 * that is not in the factory is skipped, and the test is skipped (exit status
 * 77) if none of them is.
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
#include <tuple>

#include <boost/multi_array.hpp>
#include <iostream>

#include "AbstractBlock.h"

#include "FRealObjective.h"

#include "FRowConstraint.h"

#include "LinearFunction.h"

#include "OneVarConstraint.h"

#include "Solver.h"

/*--------------------------------------------------------------------------*/

using namespace SMSpp_di_unipi_it;

/*--------------------------------------------------------------------------*/

static constexpr int T = 4;

/// the LP and its rows, the Block owning them
struct Problem {
 AbstractBlock * block;
 std::vector< ColVariable > * g;
 std::vector< FRowConstraint > * balance;
 std::vector< FRowConstraint > * storage;
 boost::multi_array< FRowConstraint , 2 > * open;
 };

/*--------------------------------------------------------------------------*/

static Problem make_problem( void )
{
 Problem pb;
 pb.block = new AbstractBlock;
 auto b = pb.block;

 pb.g = new std::vector< ColVariable >( T );
 auto y = new std::vector< ColVariable >( T );
 auto p = new std::vector< ColVariable >( T );
 auto v = new std::vector< ColVariable >( T );
 b->add_static_variable( * pb.g , "g" );
 b->add_static_variable( * y , "y" );
 b->add_static_variable( * p , "p" );
 b->add_static_variable( * v , "v" );

 auto gb = new std::vector< BoxConstraint >( T );
 auto yb = new std::vector< BoxConstraint >( T );
 auto vb = new std::vector< BoxConstraint >( T );
 for( int t = 0 ; t < T ; ++t ) {
  ( * gb )[ t ].set_variable( & ( * pb.g )[ t ] );
  ( * gb )[ t ].set_lhs( 0 );
  ( * gb )[ t ].set_rhs( 100 );
  ( * yb )[ t ].set_variable( & ( * y )[ t ] );
  ( * yb )[ t ].set_lhs( 0 );
  ( * yb )[ t ].set_rhs( 1000 );
  ( * vb )[ t ].set_variable( & ( * v )[ t ] );
  ( * vb )[ t ].set_lhs( 0 );
  ( * vb )[ t ].set_rhs( 50 );
  }
 b->add_static_constraint( * gb , "g_bounds" );
 b->add_static_constraint( * yb , "y_bounds" );
 b->add_static_constraint( * vb , "v_bounds" );

 const double d[ T ] = { 40 , 40 , 160 , 160 };
 pb.balance = new std::vector< FRowConstraint >( T );
 pb.storage = new std::vector< FRowConstraint >( T );
 pb.open = new boost::multi_array< FRowConstraint , 2 >(
                    boost::multi_array< FRowConstraint , 2 >::extent_gen()
                    [ 2 ][ T ] );

 for( int t = 0 ; t < T ; ++t ) {
  auto & bal = ( * pb.balance )[ t ];
  bal.set_function( new LinearFunction( LinearFunction::v_coeff_pair{
   { & ( * pb.g )[ t ] , 1 } , { & ( * y )[ t ] , 1 } ,
   { & ( * p )[ t ] , 1 } } ) );
  bal.set_both( d[ t ] );

  LinearFunction::v_coeff_pair sv = { { & ( * v )[ t ] , 1 } ,
                                      { & ( * p )[ t ] , 1 } };
  if( t > 0 )
   sv.push_back( { & ( * v )[ t - 1 ] , -1 } );
  ( * pb.storage )[ t ].set_function(
                                    new LinearFunction( std::move( sv ) ) );
  ( * pb.storage )[ t ].set_both( 0 );

  for( int i = 0 ; i < 2 ; ++i ) {
   ( * pb.open )[ i ][ t ].set_function( new LinearFunction(
    LinearFunction::v_coeff_pair{ { & ( * p )[ t ] , 1 } } ) );
   ( * pb.open )[ i ][ t ].set_lhs( - Inf< double >() );
   ( * pb.open )[ i ][ t ].set_rhs( Inf< double >() );
   }
  }
 b->add_static_constraint( * pb.balance , "balance" );
 b->add_static_constraint( * pb.storage , "storage" );
 b->add_static_constraint( * pb.open , "free" );

 LinearFunction::v_coeff_pair obj;
 for( int t = 0 ; t < T ; ++t )
  obj.push_back( { & ( * y )[ t ] , 1000 } );
 auto o = new FRealObjective( b , new LinearFunction( std::move( obj ) ) );
 o->set_sense( Objective::eMin , eNoMod );
 b->set_objective( o , eNoMod );

 return( pb );
 }

/*--------------------------------------------------------------------------*/
/// the value that a Solver name loaded from scratch finds on b, NaN if it
/// does not solve it

static double from_scratch( Block * b , const std::string & name )
{
 auto s = Solver::new_Solver( name );
 s->set_par( Solver::intLogVerb , 0 );
 b->register_Solver( s );
 const double v = ( s->compute() == Solver::kOK ) ? s->get_var_value()
                                                  : std::nan( "" );
 b->unregister_Solver( s );
 delete s;
 return( v );
 }

/*--------------------------------------------------------------------------*/

int main( void )
{
 int failures = 0 , tried = 0;

 for( const std::string name : { "CPXMILPSolver" , "GRBMILPSolver" ,
                                 "HiGHSMILPSolver" , "SCIPMILPSolver" ,
                                 "DAQPMILPSolver" , "QPOASESMILPSolver" ,
                                 "ClarabelMILPSolver" } ) {
  if( ! Solver::has_Solver( name ) ) {
   std::cout << name << " is not in this build, skipped" << std::endl;
   continue;
   }
  ++tried;

  auto pb = make_problem();
  auto s = Solver::new_Solver( name );
  s->set_par( Solver::intLogVerb , 0 );
  pb.block->register_Solver( s );

  // the change, its description and the expected value; the first one
  // changes nothing
  const std::vector< std::tuple< std::string , std::function< void() > ,
                                 double > > steps = {
   { "first solve" , [](){} , 70000 } ,
   { "s = 0.5" , [ & ]() {
     for( int t = 0 ; t < T ; ++t ) {
      auto f = static_cast< LinearFunction * >(
                                      ( * pb.balance )[ t ].get_function() );
      f->modify_coefficient( f->is_active( & ( * pb.g )[ t ] ) , 0.5 );
      }
     } , 200000 } ,
   { "storage rows relaxed" , [ & ]() {
     for( auto & r : * pb.storage )
      r.relax( true );
     } , 0 } ,
   { "storage rows enforced" , [ & ]() {
     for( auto & r : * pb.storage )
      r.relax( false );
     } , 200000 } ,
   { "p_2 , p_3 <= 5" , [ & ]() {
     ( * pb.open )[ 0 ][ 2 ].set_rhs( 5 );
     ( * pb.open )[ 0 ][ 3 ].set_rhs( 5 );
     } , 210000 } ,
   { "free rows free again" , [ & ]() {
     ( * pb.open )[ 0 ][ 2 ].set_rhs( Inf< double >() );
     ( * pb.open )[ 0 ][ 3 ].set_rhs( Inf< double >() );
     } , 200000 } };

  for( const auto & [ what , change , expected ] : steps ) {
   change();
   const auto st = s->compute();
   const double v = ( st == Solver::kOK ) ? s->get_var_value()
                                          : std::nan( "" );
   const double fresh = from_scratch( pb.block , name );
   const bool ok = ( std::abs( v - expected ) <= 1e-6 * ( 1 + expected ) )
                   && ( std::abs( fresh - expected ) <=
                        1e-6 * ( 1 + expected ) );
   std::cout << name << ", " << what << ": " << v << " (status " << st
             << "), from scratch " << fresh << ", expected " << expected
             << ( ok ? " -> OK" : " -> KO" ) << std::endl;
   if( ! ok )
    ++failures;
   }

  pb.block->unregister_Solver( s );
  delete s;
  delete pb.block;
  }

 if( ! tried ) {
  std::cout << "no :MILPSolver in this build: nothing to do" << std::endl;
  return( 77 );
  }

 return( failures ? 1 : 0 );
 }

/*--------------------------------------------------------------------------*/
/*---------------------- End File test_free_rows.cpp -----------------------*/
/*--------------------------------------------------------------------------*/
