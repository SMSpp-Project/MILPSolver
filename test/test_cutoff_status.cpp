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
 * bound on the other side at least the optimal value.
 *
 * Then a continuous problem, the LP
 *
 *     min   1000 + sum_j c_j x_j
 *     s.t.  sum_j a_ij x_j >= b_i , i = 0 , ... , 39 , x >= 0
 *
 * with 60 columns, c_j = 5 + ( 7 j mod 11 ), b_i = 10 + ( 3 i mod 13 ) and
 * a_ij = 1 + ( ( 5 i + 3 j ) mod 9 ) if ( i + 2 j ) mod 5 = 0, else 0, whose
 * optimum is 1138.8188034, and the same as a maximization, with the
 * objective and the constant with the opposite sign. Each Solver of
 * LPPar-cutoff.txt is given the cutoffs of LPCfg-cutoff.txt, 1069.4 and
 * -1069.4, which the dual simplex passes before the optimum in either
 * sense: CPXMILPSolver and GRBMILPSolver in either sense, and
 * HiGHSMILPSolver (whose objective bound only acts when minimizing) in the
 * minimization, have to stop there (kCutOff) with a certified bound between
 * the cutoff and the optimal value, no bound on the other side short of the
 * optimal value (the solution they stop at is not primal feasible) and a
 * dual solution that is feasible, i.e., one that gives a valid Benders cut;
 * the others may solve the LP.
 * Solved again with the cutoffs of LPCfg-cutoff-far.txt, which nothing
 * reaches, each has to give the optimal value, which tells that the
 * parameters changed for the cutoff are given back. The test is skipped
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

#include "CDASolver.h"

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
// the LP, minimizing or maximizing, with its rows in rows

static constexpr int LN = 60;              // the columns of the LP
static constexpr int LM = 40;              // the rows of the LP
static constexpr double LOpt = 1138.8188034;  // its optimal value, minimizing

static double lp_c( int j ) { return( 5 + ( 7 * j ) % 11 ); }
static double lp_b( int i ) { return( 10 + ( 3 * i ) % 13 ); }
static double lp_a( int i , int j ) {
 return( ( ( i + 2 * j ) % 5 == 0 ) ? 1 + ( 5 * i + 3 * j ) % 9 : 0 ); }

static AbstractBlock * lp( bool max , std::vector< FRowConstraint > * & rows )
{
 auto block = new AbstractBlock;
 auto & x = * new std::vector< ColVariable >( LN );
 for( auto & xj : x )
  xj.is_positive( true , eNoMod );
 block->add_static_variable( x , "x" );

 rows = new std::vector< FRowConstraint >( LM );
 for( int i = 0 ; i < LM ; ++i ) {
  LinearFunction::v_coeff_pair ax;
  for( int j = 0 ; j < LN ; ++j )
   if( lp_a( i , j ) != 0 )
    ax.push_back( { & x[ j ] , lp_a( i , j ) } );
  ( * rows )[ i ].set_function( new LinearFunction( std::move( ax ) ) );
  ( * rows )[ i ].set_lhs( lp_b( i ) );
  ( * rows )[ i ].set_rhs( Inf< double >() );
  }
 block->add_static_constraint( * rows , "rows" );

 const double s = max ? -1 : 1;
 LinearFunction::v_coeff_pair cx;
 for( int j = 0 ; j < LN ; ++j )
  cx.push_back( { & x[ j ] , s * lp_c( j ) } );
 auto obj = new FRealObjective( block ,
                                new LinearFunction( std::move( cx ) ,
                                                    s * C ) );
 obj->set_sense( max ? Objective::eMax : Objective::eMin , eNoMod );
 block->set_objective( obj , eNoMod );

 return( block );
 }

/*--------------------------------------------------------------------------*/
// whether the duals of the rows, with one sign or the other, are a feasible
// dual solution of the minimization whose value is at least least, i.e.,
// y >= 0 , A^T y <= c and b^T y >= least: one that gives a valid cut

static bool dual_feasible( const std::vector< FRowConstraint > & rows ,
                           double least )
{
 for( const double sgn : { 1.0 , -1.0 } ) {
  bool ok = true;
  double by = 0;
  std::vector< double > ay( LN , 0 );
  for( int i = 0 ; i < LM ; ++i ) {
   const double y = sgn * rows[ i ].get_dual();
   if( y < -1e-9 )
    ok = false;
   by += lp_b( i ) * y;
   for( int j = 0 ; j < LN ; ++j )
    ay[ j ] += lp_a( i , j ) * y;
   }
  for( int j = 0 ; j < LN ; ++j )
   if( ay[ j ] > lp_c( j ) + 1e-7 )
    ok = false;
  if( ok && ( by >= least - 1e-6 ) )
   return( true );
  }
 return( false );
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

 // the LP - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 const std::string lfn = "LPPar-cutoff.txt";
 auto far = dynamic_cast< ComputeConfig * >(
                    Configuration::deserialize( "LPCfg-cutoff-far.txt" ) );
 if( ! far ) {
  std::cerr << "Error: LPCfg-cutoff-far.txt is not a ComputeConfig"
            << std::endl;
  return( 1 );
  }
 const double cut = 1069.4;   // the cutoffs of LPCfg-cutoff.txt, +-cut
 const double ltol = 1e-6 * LOpt;

 for( const bool max : { false , true } ) {
  auto bsc = dynamic_cast< BlockSolverConfig * >(
                                          Configuration::deserialize( lfn ) );
  if( ! bsc ) {
   std::cerr << "Error: " << lfn << " is not a BlockSolverConfig"
             << std::endl;
   return( 1 );
   }
  Block::Index there = 0;
  for( Block::Index i = 0 ; i < bsc->num_ComputeConfig() ; ++i )
   if( Solver::has_Solver( bsc->get_SolverName( i ) ) )
    ++there;
  keep_available_Solvers( bsc , lfn , std::max( there , Block::Index( 1 ) ) );

  std::vector< FRowConstraint > * rows = nullptr;
  auto block = lp( max , rows );
  bsc->apply( block );
  bsc->clear();

  for( auto s : block->get_registered_solvers() ) {
   const std::string name = s->classname();
   const bool must = ( name == "CPXMILPSolver" ) ||
                     ( name == "GRBMILPSolver" ) ||
                     ( ( name == "HiGHSMILPSolver" ) && ! max );
   for( const bool first : { true , false } ) {
    if( ! first )
     s->set_ComputeConfig( far );
    std::cout << name << ", LP" << ( max ? ", max" : ", min" )
              << ( first ? ", cutoff: " : ", no cutoff reached: " );
    bool ok;
    try {
     const auto st = s->compute();
     std::cout << "status " << st;
     // the value and the bound as seen by the minimization
     if( st == Solver::kOK ) {
      const double v = max ? - s->get_var_value() : s->get_var_value();
      std::cout << ", value " << v;
      ok = ( std::abs( v - LOpt ) <= ltol ) && ! ( first && must );
      }
     else
      if( ( st == Solver::kCutOff ) && first ) {
       // the bound of the dual simplex, and that on the other side, which
       // the solution it stopped at, not primal feasible, does not give
       const double bnd = max ? - s->get_ub() : s->get_lb();
       const double other = max ? - s->get_lb() : s->get_ub();
       auto cds = static_cast< CDASolver * >( s );
       bool duals = cds->has_dual_solution();
       if( duals ) {
        cds->get_dual_solution();
        duals = dual_feasible( * rows , cut - C );
        }
       std::cout << ", bound " << bnd << ", other bound " << other
                 << ( duals ? ", feasible duals" : ", no feasible duals" );
       ok = ( bnd >= cut - ltol ) && ( bnd <= LOpt + ltol ) && duals &&
            ( other >= LOpt - ltol );
       }
      else
       ok = false;
     }
    catch( const std::exception & e ) {
     std::cout << "threw \"" << e.what() << "\"";
     ok = false;
     }
    std::cout << ( ok ? " -> OK" : " -> KO" ) << std::endl;
    if( ! ok )
     ++failures;
    }
   }

  bsc->apply( block );
  delete bsc;
  delete block;
  }

 delete far;

 return( failures ? 1 : 0 );
 }

/*--------------------------------------------------------------------------*/
/*-------------------- End File test_cutoff_status.cpp ---------------------*/
/*--------------------------------------------------------------------------*/
