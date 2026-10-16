/*--------------------------------------------------------------------------*/
/*---------------------- File test_cutoff_status.cpp -----------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Checks the cutoffs of each :MILPSolver [see dblUpCutOff and dblLwCutOff in
 * Solver.h] on the knapsack, with a constant in the objective,
 *
 *     max   1000 + sum_i p_i x_i
 *     s.t.  sum_i w_i x_i <= 157 , x_i in { 0 , 1 } , i = 0 , ... , 29
 *
 * with w_i = 10 + ( 37 i mod 23 ) and p_i = w_i + ( 11 i mod 7 ), whose
 * optimum is 1201, and on the same problem as a minimization, the objective
 * and the constant with the opposite sign, whose optimum is -1201. Any
 * solution, x = 0 included, is worth at least 1000 (at most -1000), which
 * a Solver that leaves out the constant does not see.
 *
 * Each Solver of MILPPar-cutoff.txt that is in the factory is attached to
 * either problem with its ComputeConfig, whose cutoffs it has to keep when
 * it loads the problem, and the first solve is to stop at a cutoff that the
 * first solution reaches (kCutOff, with bounds that show it reached and stay
 * valid) and to solve the problem otherwise (kOK and the optimal value). The
 * cutoffs of each ComputeConfig are read in the BlockSolverConfig, and are
 * those of
 *
 * - MILPCfg-cutoff.txt, an upper cutoff of 1000: the maximization stops at
 *   its first solution, the minimization is solved;
 *
 * - MILPCfg-cutoff-lw.txt, a lower cutoff of -1000: the minimization stops
 *   at its first solution, the maximization is solved;
 *
 * - MILPCfg-cutoff-bnd.txt, an upper cutoff of 10000 and a lower one of
 *   1100: the minimization stops at its first solution, while the
 *   maximization is solved, since its upper bound is never at most 1100
 *   even if its first solutions are.
 *
 * Each Solver then solves its problem again with the ComputeConfig of
 * MILPCfg-cutoff.txt, which gives a cutoff to a Solver that has solved its
 * problem with none it could reach: it may stop at the cutoff (kCutOff) or
 * give the optimal value, but it is not to throw or to fail. With kCutOff,
 * the value of the solution, seen as a maximization, has to be at least the
 * upper cutoff that the Solver has and at most the optimal value, and the
 * bound on the other side at least the optimal value. The test is skipped
 * (exit status 77) if no Solver is in the factory.
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
#include <vector>

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

static constexpr double C = 1000;     // the constant of the objective

static constexpr double Opt = 1201;   // the optimal value, maximizing

static constexpr double First = 1000; // the worst solution, x = 0

/*--------------------------------------------------------------------------*/
// the knapsack, maximizing or minimizing

static AbstractBlock * knapsack( bool max )
{
 auto block = new AbstractBlock;
 auto & x = * new std::vector< ColVariable >( N );
 for( auto & xi : x )
  xi.set_type( ColVariable::kBinary , eNoMod );
 block->add_static_variable( x , "x" );

 const double s = max ? 1 : -1;
 LinearFunction::v_coeff_pair wx , px;
 for( int i = 0 ; i < N ; ++i ) {
  const double w = 10 + ( 37 * i ) % 23;
  wx.push_back( { & x[ i ] , w } );
  px.push_back( { & x[ i ] , s * ( w + ( 11 * i ) % 7 ) } );
  }

 auto & cap = * new FRowConstraint;
 cap.set_function( new LinearFunction( std::move( wx ) ) );
 cap.set_lhs( - Inf< double >() );
 cap.set_rhs( 157 );
 block->add_static_constraint( cap , "capacity" );

 auto obj = new FRealObjective( block ,
                                new LinearFunction( std::move( px ) ,
                                                    s * C ) );
 obj->set_sense( max ? Objective::eMax : Objective::eMin , eNoMod );
 block->set_objective( obj , eNoMod );

 return( block );
 }

/*--------------------------------------------------------------------------*/
// the value of the double parameter name in cfg, or dflt if it is not there

static double par( const ComputeConfig * cfg , const std::string & name ,
                   double dflt )
{
 for( const auto & p : cfg->dbl_pars )
  if( p.first == name )
   return( p.second );
 return( dflt );
 }

/*--------------------------------------------------------------------------*/

int main( void )
{
 const std::string fn = "MILPPar-cutoff.txt";

 auto again = dynamic_cast< ComputeConfig * >(
                     Configuration::deserialize( "MILPCfg-cutoff.txt" ) );
 if( ! again ) {
  std::cerr << "Error: MILPCfg-cutoff.txt is not a ComputeConfig"
            << std::endl;
  return( 1 );
  }

 int failures = 0;
 for( const bool max : { true , false } ) {
  auto bsc = dynamic_cast< BlockSolverConfig * >(
                                           Configuration::deserialize( fn ) );
  if( ! bsc ) {
   std::cerr << "Error: " << fn << " is not a BlockSolverConfig"
             << std::endl;
   return( 1 );
   }

  Block::Index there = 0;
  for( Block::Index i = 0 ; i < bsc->num_ComputeConfig() ; ++i )
   if( Solver::has_Solver( bsc->get_SolverName( i ) ) )
    ++there;
  keep_available_Solvers( bsc , fn , std::max( there , Block::Index( 1 ) ) );

  // the cutoffs of each Solver, in the order in which they are attached
  std::vector< std::pair< double , double > > cutoffs;
  for( Block::Index i = 0 ; i < bsc->num_ComputeConfig() ; ++i ) {
   const auto cfg = bsc->get_SolverConfig( i );
   cutoffs.emplace_back( par( cfg , "dblUpCutOff" , Inf< double >() ) ,
                         par( cfg , "dblLwCutOff" , - Inf< double >() ) );
   }

  auto block = knapsack( max );
  bsc->apply( block );
  bsc->clear();  // the clear()-ed BlockSolverConfig is for the cleanup

  const double opt = max ? Opt : - Opt;
  int i = 0;
  for( auto s : block->get_registered_solvers() ) {
   // the cutoffs as seen by the maximization: the minimization of - f
   // stops where the maximization of f with the opposite cutoffs does
   const double up = max ? cutoffs[ i ].first : - cutoffs[ i ].second;
   const double lw = max ? cutoffs[ i ].second : - cutoffs[ i ].first;
   const bool stop = ( up <= First );
   if( ( ! stop ) && ( ( up <= Opt ) || ( lw >= Opt ) ) ) {
    std::cerr << "Error: the cutoffs of Solver " << i << " in " << fn
              << " neither stop at the first solution nor let the "
              << "problem be solved" << std::endl;
    return( 1 );
    }

   for( const bool first : { true , false } ) {
    if( ! first )
     s->set_ComputeConfig( again );
    std::cout << s->classname() << ( max ? ", max" : ", min" )
              << ", cutoffs " << cutoffs[ i ].first << " "
              << cutoffs[ i ].second
              << ( first ? ", first solve: " : ", second solve: " );
    bool ok;
    try {
     const auto st = s->compute();
     std::cout << "status " << st;
     if( st == Solver::kOK ) {
      const double v = s->get_var_value();
      std::cout << ", value " << v;
      ok = ( std::abs( v - opt ) <= 1e-6 * Opt ) && ! ( first && stop );
      }
     else
      if( st == Solver::kCutOff ) {
       // the bounds as seen by the maximization: the solution, which has
       // reached the upper cutoff, and the bound of the search
       const double sol = max ? s->get_lb() : - s->get_ub();
       const double bnd = max ? s->get_ub() : - s->get_lb();
       const double tol = 1e-6 * Opt;
       std::cout << ", solution " << sol << ", bound " << bnd;
       // the upper cutoff that the Solver has, as seen by the
       // maximization: the ComputeConfig of the second solve only adds to
       // the cutoffs of the first
       const double s_up = max ? s->get_dbl_par( Solver::dblUpCutOff )
                               : - s->get_dbl_par( Solver::dblLwCutOff );
       ok = ( sol >= s_up - tol ) && ( sol <= Opt + tol ) &&
            ( bnd >= Opt - tol ) && ! ( first && ( ! stop ) );
       }
      else
       ok = false;
     }
    catch( const std::exception & e ) {
     std::cout << "threw \"" << e.what() << "\"";
     ok = false;
     }

    if( first )
     std::cout << ( stop ? ", to stop" : ", to solve" );
    std::cout << ( ok ? " -> OK" : " -> KO" ) << std::endl;
    if( ! ok )
     ++failures;
    }
   ++i;
   }

  bsc->apply( block );
  delete bsc;
  delete block;
  }

 delete again;

 return( failures ? 1 : 0 );
 }

/*--------------------------------------------------------------------------*/
/*-------------------- End File test_cutoff_status.cpp ---------------------*/
/*--------------------------------------------------------------------------*/
