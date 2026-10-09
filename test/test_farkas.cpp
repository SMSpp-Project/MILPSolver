/*--------------------------------------------------------------------------*/
/*-------------------------- File test_farkas.cpp --------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Main for testing the infeasibility certificate of the :MILPSolver
 *
 * When a Block is infeasible a Solver may provide a Farkas certificate, the
 * dual ray proving that no feasible point exists, which is what a Benders
 * scheme turns into a feasibility cut. Each :MILPSolver reads it from its own
 * solver, whose sign convention is its own: Gurobi hands out FARKASDUAL with
 * the sign opposite to that of its optimal duals, SCIP likewise, CPLEX and
 * HiGHS are meant to already agree. Whoever consumes the certificate,
 * however, reads it off the Constraint with the very same code that reads the
 * optimal duals, and therefore needs the two to agree: a certificate that
 * arrives with the wrong sign is a cut in the wrong direction, which is worse
 * than no cut at all.
 *
 * This is what is checked here, on a linear program small enough that its
 * certificate is unique up to a positive factor. The same rows are given
 * twice, once with a right-hand side that makes the program feasible and once
 * with one that makes it infeasible; every registered :MILPSolver is asked
 * for the optimal duals of the former and for the certificate of the latter,
 * and the two must have the same sign on the same row. The check thus needs
 * to know no solver's convention, only that a solver keeps one.
 *
 * The program is
 *
 *     min  x + y  s.t.  2x +  y >= 3 ,  x + 2y >= 3 ,  x + y <= u ,  x , y >= 0
 *
 * with u = 3 in the feasible instance, whose optimum is x = y = 1, and u = 7/5
 * in the infeasible one. It is written this way on purpose: a certificate is
 * a ray of the dual simplex, and so it exists only where the infeasibility is
 * proved by the simplex. Were the two sides of a single variable to conflict,
 * every solver would see the empty domain by bound propagation alone and stop
 * before the first LP iteration, with no dual basis to build the ray from.
 * Here the activity bounds conclude nothing, u bounding each variable by 7/5
 * and the largest activity of the first row being 21/5, above its 3; the
 * infeasibility only shows as a combination, a third of the first row plus a
 * third of the second giving x + y >= 2 against x + y <= 7/5.
 *
 * The solvers that state they have no certificate to give are reported and
 * skipped, that being a legitimate answer: the infeasibility may have been
 * proved by the presolve, with no dual basis to build a ray from.
 *
 * Two more properties are checked on a second program, the same one with the
 * third row made two-sided and the two variables given an explicit box, so
 * that the certificate has to live on rows and bounds alike.
 *
 * The first is that the certificate is a certificate, i.e. that its value is
 * positive: summing \f$ - \pi_i b_i \f$ over rows and bounds, each taken on
 * the side its multiplier points at (they reach the Constraint negated, see
 * MILPSolver::write_dual_solution()), must give a positive number, that being
 * what proves the program infeasible and what a feasibility cut cuts away.
 * The one-sided program cannot see this: with one side infinite the side rule
 * is never exercised, so a consumer picking the wrong one is never caught.
 *
 * The second is that the certificate does not depend on the Objective. An
 * infeasibility certificate is a statement about the constraint system alone,
 * so giving the same program two different Objectives must leave the ray
 * unchanged up to the positive factor each solver is free to scale it by.
 * This is checked only where the Solver knows intHomogeneousDirection, the
 * parameter that asks for the multipliers of the ray, \f$ - A' y \f$, rather
 * than the reduced costs of an optimal dual point, \f$ c - A' y \f$.
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
#include <vector>

#include "AbstractBlock.h"

#include "FRealObjective.h"

#include "FRowConstraint.h"

#include "LinearFunction.h"

#include "OneVarConstraint.h"

#include "MILPSolver.h"

/*--------------------------------------------------------------------------*/
/*-------------------------------- USING -----------------------------------*/
/*--------------------------------------------------------------------------*/

using namespace SMSpp_di_unipi_it;

/*--------------------------------------------------------------------------*/
/*----------------------------- CONSTANTS ----------------------------------*/
/*--------------------------------------------------------------------------*/

/// the :MILPSolver to try, each skipped if it is not in the Solver factory

static const std::vector< std::string > SolverNames =
 { "CPXMILPSolver" , "GRBMILPSolver" , "SCIPMILPSolver" , "HiGHSMILPSolver" };

/// a dual value below this in absolute terms is taken to be zero

static constexpr double Eps = 1e-9;

/*--------------------------------------------------------------------------*/
/*------------------------------ FUNCTIONS ---------------------------------*/
/*--------------------------------------------------------------------------*/

/// builds the program of the file comment with the given right-hand side
/** Every Farkas certificate of the infeasible instance is a positive multiple
 * of ( 1/3 , 1/3 , 1 ): whichever way a solver scales it, the sign of each
 * component is determined, and so the comparison against the optimal duals of
 * the feasible instance is well posed. */

static AbstractBlock * build( double u ,
                              std::vector< FRowConstraint > * & rows )
{
 auto ab = new AbstractBlock();

 auto x = new std::vector< ColVariable >( 2 );
 for( auto & v : *x ) {
  v.set_type( ColVariable::kContinuous );
  v.is_positive( true );
  }
 ab->add_static_variable( *x , "x" );

 rows = new std::vector< FRowConstraint >( 3 );

 // 2x + y >= 3
 (*rows)[ 0 ].set_function( new LinearFunction(
  { std::make_pair( &(*x)[ 0 ] , 2.0 ) ,
    std::make_pair( &(*x)[ 1 ] , 1.0 ) } ) );
 (*rows)[ 0 ].set_lhs( 3 );
 (*rows)[ 0 ].set_rhs( Inf< double >() );

 // x + 2y >= 3
 (*rows)[ 1 ].set_function( new LinearFunction(
  { std::make_pair( &(*x)[ 0 ] , 1.0 ) ,
    std::make_pair( &(*x)[ 1 ] , 2.0 ) } ) );
 (*rows)[ 1 ].set_lhs( 3 );
 (*rows)[ 1 ].set_rhs( Inf< double >() );

 // x + y <= u
 (*rows)[ 2 ].set_function( new LinearFunction(
  { std::make_pair( &(*x)[ 0 ] , 1.0 ) ,
    std::make_pair( &(*x)[ 1 ] , 1.0 ) } ) );
 (*rows)[ 2 ].set_lhs( -Inf< double >() );
 (*rows)[ 2 ].set_rhs( u );

 for( auto & c : *rows )
  c.set_Block( ab );

 ab->add_static_constraint( *rows , "rows" );

 auto obj = new FRealObjective( ab , new LinearFunction(
  { std::make_pair( &(*x)[ 0 ] , 1.0 ) ,
    std::make_pair( &(*x)[ 1 ] , 1.0 ) } ) );
 obj->set_sense( Objective::eMin );
 ab->set_objective( obj );

 ab->generate_abstract_variables();
 ab->generate_abstract_constraints();
 ab->generate_objective();

 return( ab );
 }

/*--------------------------------------------------------------------------*/

/// builds the same program with the third row two-sided and a box on x and y
/** The rows are those of build(), save that the third is stated as
 * \f$ -u \leq x + y \leq u \f$, and each variable carries a BoxConstraint
 * \f$ 0 \leq x_j \leq 2 \f$ rather than being declared non-negative: the
 * certificate then has to be read off two-sided rows and bounds, each on the
 * side its multiplier points at, which is what the one-sided program cannot
 * exercise.
 *
 * \p cost is the coefficient the Objective gives to both variables, and is
 * there only to be changed: no certificate may depend on it. */

static AbstractBlock * build_boxed( double u , double cost ,
                                    std::vector< FRowConstraint > * & rows ,
                                    std::vector< BoxConstraint > * & bounds )
{
 auto ab = new AbstractBlock();

 auto x = new std::vector< ColVariable >( 2 );
 for( auto & v : *x )
  v.set_type( ColVariable::kContinuous );
 ab->add_static_variable( *x , "x" );

 bounds = new std::vector< BoxConstraint >( 2 );
 for( unsigned int j = 0 ; j < 2 ; ++j ) {
  (*bounds)[ j ].set_variable( &(*x)[ j ] );
  (*bounds)[ j ].set_lhs( 0 );
  (*bounds)[ j ].set_rhs( 2 );
  (*bounds)[ j ].set_Block( ab );
  }
 ab->add_static_constraint( *bounds , "bounds" );

 rows = new std::vector< FRowConstraint >( 3 );

 // 2x + y >= 3
 (*rows)[ 0 ].set_function( new LinearFunction(
  { std::make_pair( &(*x)[ 0 ] , 2.0 ) ,
    std::make_pair( &(*x)[ 1 ] , 1.0 ) } ) );
 (*rows)[ 0 ].set_lhs( 3 );
 (*rows)[ 0 ].set_rhs( Inf< double >() );

 // x + 2y >= 3
 (*rows)[ 1 ].set_function( new LinearFunction(
  { std::make_pair( &(*x)[ 0 ] , 1.0 ) ,
    std::make_pair( &(*x)[ 1 ] , 2.0 ) } ) );
 (*rows)[ 1 ].set_lhs( 3 );
 (*rows)[ 1 ].set_rhs( Inf< double >() );

 // -u <= x + y <= u, the two-sided one
 (*rows)[ 2 ].set_function( new LinearFunction(
  { std::make_pair( &(*x)[ 0 ] , 1.0 ) ,
    std::make_pair( &(*x)[ 1 ] , 1.0 ) } ) );
 (*rows)[ 2 ].set_lhs( -u );
 (*rows)[ 2 ].set_rhs( u );

 for( auto & c : *rows )
  c.set_Block( ab );

 ab->add_static_constraint( *rows , "rows" );

 auto obj = new FRealObjective( ab , new LinearFunction(
  { std::make_pair( &(*x)[ 0 ] , cost ) ,
    std::make_pair( &(*x)[ 1 ] , cost ) } ) );
 obj->set_sense( Objective::eMin );
 ab->set_objective( obj );

 ab->generate_abstract_variables();
 ab->generate_abstract_constraints();
 ab->generate_objective();

 return( ab );
 }

/*--------------------------------------------------------------------------*/

/// the value of the certificate written in the given Constraints
/** Sums \f$ - \pi_i b_i \f$ over the given rows and bounds, each taken on
 * the side the sign of its multiplier points at: a multiplier reaches the
 * Constraint negated, so a non-positive one belongs to the left-hand side and
 * a positive one to the right-hand side. An infinite side carries no
 * information and is skipped, the multiplier of a side that constrains
 * nothing being zero. */

template< typename T >
static double certificate_value( const std::vector< T > & constraints )
{
 double value = 0;

 for( const auto & c : constraints ) {
  const auto dual = c.get_dual();
  if( dual == 0 )
   continue;

  const auto b = ( dual > 0 ) ? c.get_rhs() : c.get_lhs();
  if( ( b <= -Inf< double >() ) || ( b >= Inf< double >() ) )
   continue;

  value -= dual * b;
  }

 return( value );
 }

/*--------------------------------------------------------------------------*/

/// the sign of \p v, zero when it is below the threshold

static int sgn( double v )
{
 if( v >  Eps ) return(  1 );
 if( v < -Eps ) return( -1 );
 return( 0 );
 }

/*--------------------------------------------------------------------------*/
/*-------------------------------- MAIN ------------------------------------*/
/*--------------------------------------------------------------------------*/

int main( void )
{
 bool all_passed = true;
 bool any_solver = false;

 // the sign the first row's multiplier has in the certificate, as agreed on
 // by the solvers seen so far; Inf marks "nobody has spoken yet"
 int agreed_sign = 0;
 std::string agreed_by;

 for( const auto & name : SolverNames ) {

  if( ! Solver::has_Solver( name ) )
   continue;                 // this one is not in the build, nothing to say
  auto solver = Solver::new_Solver( name );
  any_solver = true;

  // silence the solver's own log, and ask for the duals: SCIP does not
  // compute them unless told to, and refuses to hand them out otherwise
  solver->set_par( MILPSolver::intLogVerb , 0 );

  // a certificate is a ray of the dual simplex, so it exists only if the
  // infeasibility is proved by the simplex: on a program this small every
  // presolve gets there first, and then there is no dual basis to build it
  // from. Each solver spells the switch its own way, and the ones it does
  // not know are skipped
  auto try_int = [ &solver ]( const char * name , int value ) {
   const auto par = solver->int_par_str2idx( name );
   if( par < Inf< Solver::idx_type >() )
    solver->set_par( par , value );
   };
  auto try_str = [ &solver ]( const char * name , const char * value ) {
   const auto par = solver->str_par_str2idx( name );
   if( par < Inf< Solver::idx_type >() )
    solver->set_par( par , std::string( value ) );
   };

  try_int( "intComputeDuals" , 1 );                     // SCIP
  try_int( "Presolve" , 0 );                            // Gurobi
  try_int( "InfUnbdInfo" , 1 );                         // Gurobi
  try_int( "CPXPARAM_Preprocessing_Presolve" , 0 );     // CPLEX
  try_int( "presolving/maxrounds" , 0 );                // SCIP
  try_str( "presolve" , "off" );                        // HiGHS
  try_str( "solver" , "simplex" );                      // HiGHS

  std::cout << std::left << std::setw( 16 ) << name;

  // the feasible instance, whose optimum is x = y = 1 - - - - - - - - - - -

  std::vector< FRowConstraint > * frows;
  auto feasible = build( 3 , frows );
  feasible->register_Solver( solver );

  auto status = solver->compute();
  int solution_sign = 0;

  if( status != Solver::kOK )
   std::cout << " the feasible instance is not solved (status "
             << status << ")";
  else {
   auto cda = dynamic_cast< CDASolver * >( solver );
   try {
    if( cda && cda->has_dual_solution() ) {
     cda->get_dual_solution();
     solution_sign = sgn( (*frows)[ 0 ].get_dual() );
     }
    }
   catch( const std::exception & e ) {
    std::cout << " the optimal duals are refused: " << e.what();
    }
   }

  feasible->unregister_Solver( solver );
  delete feasible;

  if( ! solution_sign ) {
   std::cout << " no optimal dual on the binding row, nothing to compare"
             << std::endl;
   all_passed = false;
   delete solver;
   continue;
   }

  // the infeasible instance - - - - - - - - - - - - - - - - - - - - - - - -

  std::vector< FRowConstraint > * irows;
  auto infeasible = build( 1.4 , irows );
  infeasible->register_Solver( solver );

  status = solver->compute();

  if( status != Solver::kInfeasible ) {
   std::cout << " does not report the instance as infeasible (status "
             << status << ")" << std::endl;
   all_passed = false;
   }
  else {
   auto cda = dynamic_cast< CDASolver * >( solver );

   bool has = false;
   try { has = cda && cda->has_dual_direction(); }
   catch( const std::exception & e ) {
    std::cout << " has_dual_direction() throws: " << e.what() << std::endl;
    all_passed = false;
    }

   if( ! has )
    // a legitimate answer: without a dual basis there is no ray to give
    std::cout << " has no certificate to give" << std::endl;
   else {
    cda->get_dual_direction();

    const auto d0 = (*irows)[ 0 ].get_dual();
    const auto d1 = (*irows)[ 1 ].get_dual();
    const auto d2 = (*irows)[ 2 ].get_dual();

    if( ( ! sgn( d0 ) ) && ( ! sgn( d1 ) ) && ( ! sgn( d2 ) ) ) {
     std::cout << " gives an all-zero certificate" << std::endl;
     all_passed = false;
     }
    else {
     std::cout << " certificate ( " << std::showpos << std::setw( 9 ) << d0
               << " , " << std::setw( 9 ) << d1
               << " , " << std::setw( 9 ) << d2 << std::noshowpos << " )";

     // the certificate must have, on the same row, the sign the optimal dual
     // has: whoever reads it does so with the same code
     if( sgn( d0 ) != solution_sign ) {
      std::cout << " -> KO, the sign is not that of the optimal dual";
      all_passed = false;
      }
     else if( agreed_sign && ( sgn( d0 ) != agreed_sign ) ) {
      std::cout << " -> KO, it disagrees with " << agreed_by;
      all_passed = false;
      }
     else {
      if( ! agreed_sign ) { agreed_sign = sgn( d0 ); agreed_by = name; }
      std::cout << " -> OK";
      }

     std::cout << std::endl;
     }
    }
   }

  infeasible->unregister_Solver( solver );
  delete infeasible;

  // the two-sided program: is the certificate a certificate? - - - - - - - -

  auto ray_of = [ & ]( double cost , std::vector< double > & ray ,
		       double & value ) -> bool {
   std::vector< FRowConstraint > * brows;
   std::vector< BoxConstraint > * bbnds;
   auto boxed = build_boxed( 1.4 , cost , brows , bbnds );
   boxed->register_Solver( solver );

   bool got = false;
   if( solver->compute() == Solver::kInfeasible ) {
    auto cda = dynamic_cast< CDASolver * >( solver );
    bool has = false;
    try { has = cda && cda->has_dual_direction(); }
    catch( const std::exception & ) {}

    if( has ) {
     cda->get_dual_direction();

     ray.clear();
     for( const auto & c : *brows )
      ray.push_back( c.get_dual() );
     for( const auto & c : *bbnds )
      ray.push_back( c.get_dual() );

     value = certificate_value( *brows ) + certificate_value( *bbnds );
     got = true;
     }
    }

   boxed->unregister_Solver( solver );
   delete boxed;
   return( got );
   };

  std::cout << std::left << std::setw( 16 ) << name << " two-sided:";

  std::vector< double > ray;
  double value = 0;

  if( ! ray_of( 1.0 , ray , value ) )
   std::cout << " no certificate to give" << std::endl;
  else {
   std::cout << " value " << value;

   if( value <= Eps ) {
    std::cout << " -> KO, a certificate has to be positive";
    all_passed = false;
    }
   else
    std::cout << " -> OK";

   /* The multipliers of a ray are tied to one another: with the homogeneous
    * ones asked for, the multiplier of the bound of each column is minus the
    * combination of the rows through that column. In the sign convention the
    * Constraints carry (they hold the multipliers negated) this reads
    *
    *     b_j + sum_i a_ij d_i = 0
    *
    * which holds whichever of the rays the solver happens to return, and is
    * what the reduced costs of an optimal dual point, c - A' y, do not
    * satisfy: the Objective has no part in a certificate. */
   const auto par = solver->int_par_str2idx( "intHomogeneousDirection" );
   if( par >= Inf< Solver::idx_type >() )
    std::cout << ", the columns are not checked, this Solver has no"
	      << " intHomogeneousDirection";
   else {
    // a Solver that cannot produce them says so rather than taking the
    // parameter and going on, which is a legitimate answer and not a failure
    bool asked = true;
    try { solver->set_par( par , 1 ); }
    catch( const std::exception & ) {
     asked = false;
     std::cout << ", the columns are not checked, this Solver refuses the"
	       << " homogeneous direction";
     }

    if( asked ) {
     std::vector< double > hray;
     double hvalue = 0;

     if( ! ray_of( 1.0 , hray , hvalue ) )
      std::cout << ", the columns are not checked, no certificate with it";
     else {
      // the coefficient of column j in row i, as build_boxed() states them
      static const double a[ 3 ][ 2 ] = { { 2 , 1 } , { 1 , 2 } , { 1 , 1 } };

      double worst = 0;
      for( std::size_t j = 0 ; j < 2 ; ++j ) {
       double sum = hray[ 3 + j ];
       for( std::size_t i = 0 ; i < 3 ; ++i )
	sum += a[ i ][ j ] * hray[ i ];
       worst = std::max( worst , std::abs( sum ) );
       }

      if( worst > 1e-6 ) {
       std::cout << ", the columns are off by " << worst << " -> KO";
       all_passed = false;
       }
      else
       std::cout << ", the columns are the ray's -> OK";
      }

     solver->set_par( par , 0 );
     }
    }

   std::cout << std::endl;
   }

  delete solver;
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
/*-------------------------- File test_farkas.cpp --------------------------*/
/*--------------------------------------------------------------------------*/
