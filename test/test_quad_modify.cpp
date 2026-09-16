/*--------------------------------------------------------------------------*/
/*----------------------- File test_quad_modify.cpp ------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Checks that a change of an off-diagonal coefficient of a QuadFunction in
 * the Objective reaches a :MILPSolver right: the problem
 *
 *     min  x_1^2 + x_2^2 + x_3^2 + a x_1 x_2 + b x_2 x_3 - x_1 - 2 x_2 - 3 x_3
 *     s.t. -10 <= x_i <= 10
 *
 * is loaded with ( a0 , b0 ), a and b are changed by
 * QuadFunction::modify_term() and the problem is solved again, and the value
 * is compared with that of a fresh Solver of the same type loaded with the
 * new coefficients. Each :MILPSolver among CPXMILPSolver, GRBMILPSolver and
 * HiGHSMILPSolver that is not in the factory is skipped.
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

#include "AbstractBlock.h"

#include "FRealObjective.h"

#include "OneVarConstraint.h"

#include "QuadFunction.h"

#include "Solver.h"

/*--------------------------------------------------------------------------*/

using namespace SMSpp_di_unipi_it;

/*--------------------------------------------------------------------------*/

struct Problem {
 AbstractBlock * block;
 QuadFunction * f;
 };

static Problem make_problem( double a , double b )
{
 auto block = new AbstractBlock;

 std::vector< ColVariable > & x = * new std::vector< ColVariable >( 3 );
 block->add_static_variable( x , "x" );

 auto & box = * new std::vector< BoxConstraint >( 3 );
 for( int i = 0 ; i < 3 ; ++i ) {
  box[ i ].set_variable( & x[ i ] );
  box[ i ].set_lhs( -10 );
  box[ i ].set_rhs( 10 );
  }
 block->add_static_constraint( box , "box" );

 QuadFunction::v_coeff_triple diag;
 for( int i = 0 ; i < 3 ; ++i )
  diag.emplace_back( & x[ i ] , - double( i + 1 ) , 1.0 );

 QuadFunction::v_off_diag_term offdiag = { { 1 , 0 , a } , { 2 , 1 , b } };

 auto f = new QuadFunction( std::move( diag ) , std::move( offdiag ) );
 auto obj = new FRealObjective( block , f );
 obj->set_sense( Objective::eMin , eNoMod );
 block->set_objective( obj , eNoMod );

 return( Problem{ block , f } );
 }

/*--------------------------------------------------------------------------*/

static bool solve( Solver * s , double & value )
{
 const auto status = s->compute();
 if( ( status < Solver::kOK ) || ( status >= Solver::kError ) )
  return( false );
 value = s->get_var_value();
 return( true );
 }

/*--------------------------------------------------------------------------*/

int main( void )
{
 const double a0 = 0.3 , b0 = -0.2 , a1 = -0.7 , b1 = 0.6;
 int failures = 0;

 for( const std::string name :
       { "CPXMILPSolver" , "GRBMILPSolver" , "HiGHSMILPSolver" } ) {
  Solver * s;
  try {
   s = Solver::new_Solver( name );
   }
  catch( ... ) {
   std::cout << name << ": not in the factory, skipped" << std::endl;
   continue;
   }

  s->set_par( Solver::intLogVerb , 0 );

  auto p = make_problem( a0 , b0 );
  p.block->register_Solver( s );

  double before , after , fresh;
  if( ! solve( s , before ) ) {
   std::cout << name << ": the first solve failed" << std::endl;
   ++failures;
   continue;
   }

  p.f->modify_term( 1 , 0 , a1 );
  p.f->modify_term( 2 , 1 , b1 );
  const bool ok_after = solve( s , after );

  auto q = make_problem( a1 , b1 );
  auto s2 = Solver::new_Solver( name );
  s2->set_par( Solver::intLogVerb , 0 );
  q.block->register_Solver( s2 );
  const bool ok_fresh = solve( s2 , fresh );

  const bool ok = ok_after && ok_fresh &&
   ( std::abs( after - fresh ) <= 1e-6 * std::max( 1.0 , std::abs( fresh ) ) );

  std::cout << name << ": before = " << before << ", after the change = "
            << after << ", loaded anew = " << fresh
            << ( ok ? " -> OK" : " -> KO" ) << std::endl;
  if( ! ok )
   ++failures;

  p.block->unregister_Solver( s );
  q.block->unregister_Solver( s2 );
  delete s;
  delete s2;
  }

 return( failures ? 1 : 0 );
 }

/*--------------------------------------------------------------------------*/
/*--------------------- End File test_quad_modify.cpp ----------------------*/
/*--------------------------------------------------------------------------*/
