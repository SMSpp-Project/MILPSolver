/*--------------------------------------------------------------------------*/
/*---------------------- File test_ipm_free_rows.cpp -----------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Checks the interior point of each :MILPSolver, without presolve, on an LP
 * with a row with no finite side followed by a row with two finite sides:
 *
 *     min   x + y
 *     s.t.  - inf <= x - y <= + inf      (row 0, free)
 *           1     <= x + y <= 4          (row 1)
 *           0.5   <= x - y <= + inf      (row 2)
 *           0 <= x , y <= 10
 *
 * whose optimum is 1, with dual values 0 on row 0 and -1 on row 1 (in the
 * sign convention of RowConstraint, see test_dual). Up to
 * HiGHS 1.15.1, the LP that HiGHS builds for its interior point (IPX) drops
 * the free rows and renumbers the others, but puts the slack of a row with
 * two finite sides in the row with its original index (ERGO-Code/HiGHS#3319):
 * here the slack of row 1 goes to row 2, and IPX says the LP infeasible. The
 * same Solver then solves the LP again after each of these changes, and its
 * value has to be the expected one:
 *
 * - row 0 at least 2: 2;
 * - row 0 free again: 1;
 * - row 1 relaxed: 0.5;
 * - row 2 at most 3 as well, a row with two finite sides after two free
 *   ones: 0.5;
 * - row 1 enforced again: 1.
 *
 * The Solver are those of LPPar-ipm.txt that are in the factory, all of them
 * attached at once, and the test is skipped (exit status 77) if none is.
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

int main( void )
{
 const std::string fn = "LPPar-ipm.txt";
 auto bsc = dynamic_cast< BlockSolverConfig * >(
                                           Configuration::deserialize( fn ) );
 if( ! bsc ) {
  std::cerr << "Error: " << fn << " is not a BlockSolverConfig" << std::endl;
  return( 1 );
  }

 Block::Index there = 0;
 for( Block::Index i = 0 ; i < bsc->num_ComputeConfig() ; ++i )
  if( Solver::has_Solver( bsc->get_SolverName( i ) ) )
   ++there;
 keep_available_Solvers( bsc , fn , std::max( there , Block::Index( 1 ) ) );

 // the LP- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 auto block = new AbstractBlock;
 auto & x = * new std::vector< ColVariable >( 2 );
 block->add_static_variable( x , "x" );

 auto & box = * new std::vector< BoxConstraint >( 2 );
 for( int j = 0 ; j < 2 ; ++j ) {
  box[ j ].set_variable( & x[ j ] );
  box[ j ].set_lhs( 0 );
  box[ j ].set_rhs( 10 );
  }
 block->add_static_constraint( box , "box" );

 auto & row = * new std::vector< FRowConstraint >( 3 );
 const double a[ 3 ][ 2 ] = { { 1 , -1 } , { 1 , 1 } , { 1 , -1 } };
 for( int i = 0 ; i < 3 ; ++i )
  row[ i ].set_function( new LinearFunction( LinearFunction::v_coeff_pair{
   { & x[ 0 ] , a[ i ][ 0 ] } , { & x[ 1 ] , a[ i ][ 1 ] } } ) );
 row[ 0 ].set_lhs( - Inf< double >() );
 row[ 0 ].set_rhs( Inf< double >() );
 row[ 1 ].set_lhs( 1 );
 row[ 1 ].set_rhs( 4 );
 row[ 2 ].set_lhs( 0.5 );
 row[ 2 ].set_rhs( Inf< double >() );
 block->add_static_constraint( row , "row" );

 auto obj = new FRealObjective( block , new LinearFunction(
  LinearFunction::v_coeff_pair{ { & x[ 0 ] , 1 } , { & x[ 1 ] , 1 } } ) );
 obj->set_sense( Objective::eMin , eNoMod );
 block->set_objective( obj , eNoMod );

 bsc->apply( block );
 bsc->clear();  // the clear()-ed BlockSolverConfig is for the final cleanup

 // the changes - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // the change, its description and the expected value; the first one
 // changes nothing
 const std::vector< std::tuple< std::string , std::function< void() > ,
                                double > > steps = {
  { "first solve" , [](){} , 1 } ,
  { "row 0 >= 2" , [ & ]() { row[ 0 ].set_lhs( 2 ); } , 2 } ,
  { "row 0 free again" , [ & ]() {
    row[ 0 ].set_lhs( - Inf< double >() ); } , 1 } ,
  { "row 1 relaxed" , [ & ]() { row[ 1 ].relax( true ); } , 0.5 } ,
  { "row 2 <= 3" , [ & ]() { row[ 2 ].set_rhs( 3 ); } , 0.5 } ,
  { "row 1 enforced" , [ & ]() { row[ 1 ].relax( false ); } , 1 } };

 int failures = 0;
 for( const auto & [ what , change , expected ] : steps ) {
  change();
  for( auto s : block->get_registered_solvers() ) {
   const auto name = s->classname();
   const auto st = s->compute();
   const double v = ( st == Solver::kOK ) ? s->get_var_value()
                                          : std::nan( "" );
   bool ok = ( std::abs( v - expected ) <= 1e-6 * ( 1 + expected ) );
   std::cout << name << ", " << what << ": " << v << " (status " << st
             << "), expected " << expected;

   // on the first solve, the dual values of the free row and of row 1
   if( ok && ( what == "first solve" ) ) {
    auto cda = dynamic_cast< CDASolver * >( s );
    if( cda && cda->has_dual_solution() ) {
     cda->get_dual_solution();
     const double d0 = row[ 0 ].get_dual() , d1 = row[ 1 ].get_dual();
     std::cout << ", duals " << d0 << " " << d1;
     ok = ( std::abs( d0 ) <= 1e-6 ) && ( std::abs( d1 + 1 ) <= 1e-6 );
     }
    else {
     std::cout << ", no dual solution";
     ok = false;
     }
    }

   std::cout << ( ok ? " -> OK" : " -> KO" ) << std::endl;
   if( ! ok )
    ++failures;
   }
  }

 bsc->apply( block );
 delete bsc;
 delete block;

 return( failures ? 1 : 0 );
 }

/*--------------------------------------------------------------------------*/
/*-------------------- End File test_ipm_free_rows.cpp ---------------------*/
/*--------------------------------------------------------------------------*/
