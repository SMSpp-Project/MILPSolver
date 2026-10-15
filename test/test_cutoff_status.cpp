/*--------------------------------------------------------------------------*/
/*---------------------- File test_cutoff_status.cpp -----------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Checks the status of each :MILPSolver stopped by its upper cutoff on a
 * maximization MIP, the knapsack
 *
 *     max   sum_i p_i x_i
 *     s.t.  sum_i w_i x_i <= 157 , x_i in { 0 , 1 } , i = 0 , ... , 29
 *
 * with w_i = 10 + ( 37 i mod 23 ) and p_i = w_i + ( 11 i mod 7 ), whose
 * optimum is 201. Each Solver solves it, which loads the problem and sets
 * the cutoffs back to infinite, and then solves it again with dblUpCutOff =
 * 1 [see MILPCfg-cutoff.txt]: the first integer solution reaches the
 * cutoff, and the Solver may stop there. A Solver that stops is not to
 * throw: HiGHS ends such a run with the model status Interrupt, which
 * HiGHSMILPSolver did not know. The status has to be kOK, with the optimal
 * value, or kError, that of an aborted run.
 *
 * The Solver are those of MILPPar-cutoff.txt that are in the factory, all
 * of them attached at once, and the test is skipped (exit status 77) if none
 * is; MILPCfg-cutoff.txt is the ComputeConfig of each of them.
 *
 * \author Donato Meoli \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 */

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include <cmath>
#include <exception>
#include <iostream>

#include "AbstractBlock.h"

#include "BlockSolverConfig.h"

#include "FRealObjective.h"

#include "FRowConstraint.h"

#include "LinearFunction.h"

#include "test_common.h"

/*--------------------------------------------------------------------------*/

using namespace SMSpp_di_unipi_it;

/*--------------------------------------------------------------------------*/

static constexpr int N = 30;

/*--------------------------------------------------------------------------*/

int main( void )
{
 const std::string fn = "MILPPar-cutoff.txt";
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

 // the knapsack- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 auto block = new AbstractBlock;
 auto & x = * new std::vector< ColVariable >( N );
 for( auto & xi : x )
  xi.set_type( ColVariable::kBinary , eNoMod );
 block->add_static_variable( x , "x" );

 LinearFunction::v_coeff_pair wx , px;
 for( int i = 0 ; i < N ; ++i ) {
  const double w = 10 + ( 37 * i ) % 23;
  wx.push_back( { & x[ i ] , w } );
  px.push_back( { & x[ i ] , w + ( 11 * i ) % 7 } );
  }

 auto & cap = * new FRowConstraint;
 cap.set_function( new LinearFunction( std::move( wx ) ) );
 cap.set_lhs( - Inf< double >() );
 cap.set_rhs( 157 );
 block->add_static_constraint( cap , "capacity" );

 auto obj = new FRealObjective( block ,
                                new LinearFunction( std::move( px ) ) );
 obj->set_sense( Objective::eMax , eNoMod );
 block->set_objective( obj , eNoMod );

 bsc->apply( block );
 bsc->clear();  // the clear()-ed BlockSolverConfig is for the final cleanup

 // the solves- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 auto cfg = dynamic_cast< ComputeConfig * >(
                     Configuration::deserialize( "MILPCfg-cutoff.txt" ) );
 if( ! cfg ) {
  std::cerr << "Error: MILPCfg-cutoff.txt is not a ComputeConfig"
            << std::endl;
  return( 1 );
  }

 int failures = 0;
 for( const bool cutoff : { false , true } )
  for( auto s : block->get_registered_solvers() ) {
   if( cutoff )
    s->set_ComputeConfig( cfg );
   const auto name = s->classname();
   std::cout << name << ( cutoff ? ", with" : ", without" ) << " cutoff: ";
   bool ok;
   try {
    const auto st = s->compute();
    std::cout << "status " << st;
    if( st == Solver::kOK ) {
     const double v = s->get_var_value();
     std::cout << ", value " << v;
     ok = ( std::abs( v - 201 ) <= 1e-6 * 201 );
     }
    else
     ok = cutoff && ( st == Solver::kError );
    }
   catch( const std::exception & e ) {
    std::cout << "threw \"" << e.what() << "\"";
    ok = false;
    }

   std::cout << ( ok ? " -> OK" : " -> KO" ) << std::endl;
   if( ! ok )
    ++failures;
   }

 delete cfg;

 bsc->apply( block );
 delete bsc;
 delete block;

 return( failures ? 1 : 0 );
 }

/*--------------------------------------------------------------------------*/
/*-------------------- End File test_cutoff_status.cpp ---------------------*/
/*--------------------------------------------------------------------------*/
