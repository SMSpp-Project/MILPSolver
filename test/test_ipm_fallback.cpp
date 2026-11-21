/*--------------------------------------------------------------------------*/
/*---------------------- File test_ipm_fallback.cpp ------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Main for testing that a :MILPSolver whose simplex is stopped by
 * intIPMFallback solves the LP by the interior point method, and that the
 * next solve, started from the basis that the crossover leaves, is right
 *
 * The LP is a transportation problem with Supplies sources and Demands
 * destinations, whose rows are equalities (one of them redundant, hence the
 * LP is degenerate) and whose columns are nonnegative; it is solved Steps
 * times, each time after changing the cost of a few arcs, so that every
 * solve but the first starts from the basis of the previous one.
 *
 * Every :MILPSolver of LPPar-ipm-fallback.txt that is in the build solves
 * the sequence twice, with the simplex alone and with intIPMFallback = 1,
 * i.e., the interior point method after the first iteration of the
 * simplex, and with dblKKTTol set, so that the optimality conditions of
 * every solution are checked by MILPSolver::compute(); the optimal values
 * must be the same, up to Eps relative, for all of them.
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
#include <iomanip>
#include <iostream>
#include <memory>
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

/// the BlockSolverConfig with the :MILPSolver to compare

static const std::string Config = "LPPar-ipm-fallback.txt";

static constexpr int Supplies = 4;   ///< number of sources

static constexpr int Demands = 6;    ///< number of destinations

static constexpr int Steps = 6;      ///< number of solves

static constexpr int Changes = 5;    ///< arcs whose cost changes at each step

/// the relative tolerance of the comparison of the optimal values

static constexpr double Eps = 1e-7;

/*--------------------------------------------------------------------------*/
/*------------------------------ FUNCTIONS ---------------------------------*/
/*--------------------------------------------------------------------------*/

/// a deterministic sequence of integers in [ 1 , 20 ]

static double next_cost( unsigned int & seed )
{
 seed = seed * 1103515245u + 12345u;
 return( 1 + ( ( seed >> 16 ) % 20 ) );
 }

/*--------------------------------------------------------------------------*/
/// the transportation problem of the file comment; x gets its columns

static AbstractBlock * build( std::vector< ColVariable > * & x ,
                              LinearFunction * & obj )
{
 static const double supply[ Supplies ] = { 30 , 25 , 35 , 20 };
 static const double demand[ Demands ] = { 15 , 20 , 10 , 25 , 30 , 10 };

 auto ab = new AbstractBlock();

 x = new std::vector< ColVariable >( Supplies * Demands );
 ab->add_static_variable( *x , "x" );

 auto xb = new std::vector< BoxConstraint >( x->size() );
 for( std::size_t k = 0 ; k < x->size() ; ++k ) {
  (*xb)[ k ].set_variable( & (*x)[ k ] , eNoMod );
  (*xb)[ k ].set_lhs( 0 , eNoMod );
  (*xb)[ k ].set_rhs( Inf< double >() , eNoMod );
  }
 ab->add_static_constraint( *xb , "x_bound" );

 auto rows = new std::vector< FRowConstraint >( Supplies + Demands );
 for( int i = 0 ; i < Supplies ; ++i ) {
  LinearFunction::v_coeff_pair terms;
  for( int j = 0 ; j < Demands ; ++j )
   terms.push_back( std::make_pair( & (*x)[ i * Demands + j ] , 1.0 ) );
  (*rows)[ i ].set_function( new LinearFunction( std::move( terms ) ) ,
                             eNoMod );
  (*rows)[ i ].set_both( supply[ i ] , eNoMod );
  }
 for( int j = 0 ; j < Demands ; ++j ) {
  LinearFunction::v_coeff_pair terms;
  for( int i = 0 ; i < Supplies ; ++i )
   terms.push_back( std::make_pair( & (*x)[ i * Demands + j ] , 1.0 ) );
  (*rows)[ Supplies + j ].set_function(
                           new LinearFunction( std::move( terms ) ) , eNoMod );
  (*rows)[ Supplies + j ].set_both( demand[ j ] , eNoMod );
  }
 ab->add_static_constraint( *rows , "rows" );

 unsigned int seed = 1;
 LinearFunction::v_coeff_pair terms;
 for( auto & v : *x )
  terms.push_back( std::make_pair( & v , next_cost( seed ) ) );
 obj = new LinearFunction( std::move( terms ) );
 auto fobj = new FRealObjective( ab , obj );
 fobj->set_sense( Objective::eMin , eNoMod );
 ab->set_objective( fobj , eNoMod );

 return( ab );
 }

/*--------------------------------------------------------------------------*/
/// solves the sequence with solver; returns false if a solve is not kOK

static bool solve( Solver * solver , std::vector< double > & values )
{
 std::vector< ColVariable > * x;
 LinearFunction * obj;
 std::unique_ptr< AbstractBlock > ab( build( x , obj ) );
 ab->register_Solver( solver );

 unsigned int seed = 7;
 bool ok = true;
 values.clear();
 for( int s = 0 ; s < Steps ; ++s ) {
  if( s > 0 )
   for( int c = 0 ; c < Changes ; ++c ) {
    const auto k = Block::Index( next_cost( seed ) ) % x->size();
    obj->modify_coefficient( k , next_cost( seed ) );
    }

  const auto status = solver->compute( false );
  if( status != Solver::kOK ) {
   std::cout << " step " << s << ": status " << status;
   ok = false;
   values.push_back( Inf< double >() );
   }
  else
   values.push_back( static_cast< CDASolver * >( solver )->get_var_value() );
  }

 ab->unregister_Solver( solver );
 return( ok );
 }

/*--------------------------------------------------------------------------*/
/*-------------------------------- main() ----------------------------------*/
/*--------------------------------------------------------------------------*/

int main( void )
{
 auto c = Configuration::deserialize( Config );
 std::unique_ptr< BlockSolverConfig > bsc(
                                    dynamic_cast< BlockSolverConfig * >( c ) );
 if( ! bsc ) {
  delete c;
  std::cout << Config << " is not a BlockSolverConfig" << std::endl;
  return( 1 );
  }

 bool all_passed = true;
 std::vector< double > reference;
 std::string reference_name;

 for( Block::Index i = 0 ; i < bsc->num_ComputeConfig() ; ++i ) {
  const auto sname = bsc->get_SolverName( i );
  if( ! Solver::has_Solver( sname ) )
   continue;                  // this one is not in the build, nothing to say

  std::unique_ptr< Solver > solver( Solver::new_Solver( sname ) );
  if( auto cfg = bsc->get_SolverConfig( i ) )
   solver->set_ComputeConfig( cfg );
  const bool fallback = solver->get_int_par(
                     solver->int_par_str2idx( "intIPMFallback" ) ) > 0;

  std::cout << std::left << std::setw( 16 ) << sname
            << ( fallback ? " barrier: " : " simplex: " );

  std::vector< double > values;
  bool ok = solve( solver.get() , values );

  if( ok && reference.empty() ) {
   reference = values;
   reference_name = sname;
   }
  else
   if( ok )
    for( int s = 0 ; s < Steps ; ++s )
     if( std::abs( values[ s ] - reference[ s ] ) >
         Eps * std::max( 1.0 , std::abs( reference[ s ] ) ) ) {
      std::cout << " step " << s << ": " << values[ s ] << " != "
                << reference[ s ] << " of " << reference_name;
      ok = false;
      }

  std::cout << ( ok ? "-> OK" : " -> Error" ) << std::endl;
  if( ! ok )
   all_passed = false;
  }

 if( reference.empty() ) {
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
/*---------------------- End File test_ipm_fallback.cpp --------------------*/
/*--------------------------------------------------------------------------*/
