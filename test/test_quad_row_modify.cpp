/*--------------------------------------------------------------------------*/
/*--------------------- File test_quad_row_modify.cpp ----------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Checks that a change of a row with a DQuadFunction reaches a :MILPSolver
 * right: the problem
 *
 *     min  - x_1 - 2 x_2
 *     s.t. q_1 x_1^2 + q_2 x_2^2 + l_1 x_1 <= r
 *          -10 <= x_i <= 10
 *
 * is loaded and solved, then (i) the coefficients q_1 and l_1 are changed by
 * DQuadFunction::modify_term(), (ii) the right-hand side r is changed, and
 * (iii) a Variable is removed from the DQuadFunction, and after each change
 * the problem is solved again and the value is compared with that of a fresh
 * Solver of the same type loaded with the changed data. Each :MILPSolver
 * among CPXMILPSolver, GRBMILPSolver, MOSEKMILPSolver and SCIPMILPSolver,
 * the ones that take quadratic rows, that is not in the factory is skipped.
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

#include "DQuadFunction.h"

#include "FRealObjective.h"

#include "FRowConstraint.h"

#include "LinearFunction.h"

#include "OneVarConstraint.h"

#include "Solver.h"

/*--------------------------------------------------------------------------*/

using namespace SMSpp_di_unipi_it;

/*--------------------------------------------------------------------------*/

struct Data {
 double q1 , q2 , l1 , r;
 bool both;  // whether x_2 is in the row
 };

struct Problem {
 AbstractBlock * block;
 std::vector< ColVariable > * x;
 FRowConstraint * row;
 DQuadFunction * f;
 };

static Problem make_problem( const Data & d )
{
 auto block = new AbstractBlock;

 auto x = new std::vector< ColVariable >( 2 );
 block->add_static_variable( * x , "x" );

 auto & box = * new std::vector< BoxConstraint >( 2 );
 for( int i = 0 ; i < 2 ; ++i ) {
  box[ i ].set_variable( & ( * x )[ i ] );
  box[ i ].set_lhs( -10 );
  box[ i ].set_rhs( 10 );
  }
 block->add_static_constraint( box , "box" );

 DQuadFunction::v_coeff_triple terms = { { & ( * x )[ 0 ] , d.l1 , d.q1 } };
 if( d.both )
  terms.emplace_back( & ( * x )[ 1 ] , 0.0 , d.q2 );

 auto & row = * new FRowConstraint;
 auto f = new DQuadFunction( std::move( terms ) );
 row.set_function( f );
 row.set_lhs( - Inf< double >() );
 row.set_rhs( d.r );
 block->add_static_constraint( row , "row" );

 auto obj = new FRealObjective( block , new LinearFunction(
		     { { & ( * x )[ 0 ] , -1.0 } , { & ( * x )[ 1 ] , -2.0 } } ) );
 obj->set_sense( Objective::eMin , eNoMod );
 block->set_objective( obj , eNoMod );

 return( Problem{ block , x , & row , f } );
 }

/*--------------------------------------------------------------------------*/

static bool solve( Solver * s , double & value )
{
 int status;
 try {
  status = s->compute();
  }
 catch( std::exception & e ) {
  std::cout << "  " << e.what() << std::endl;
  return( false );
  }
 if( ( status < Solver::kOK ) || ( status >= Solver::kError ) )
  return( false );
 value = s->get_var_value();
 return( true );
 }

/*--------------------------------------------------------------------------*/

static bool fresh_value( const std::string & name , const Data & d ,
			 double & value )
{
 auto q = make_problem( d );
 auto s = Solver::new_Solver( name );
 s->set_par( Solver::intLogVerb , 0 );
 q.block->register_Solver( s );
 const bool ok = solve( s , value );
 q.block->unregister_Solver( s );
 delete s;
 return( ok );
 }

/*--------------------------------------------------------------------------*/

int main( void )
{
 int failures = 0;

 for( const std::string name :
      { "CPXMILPSolver" , "GRBMILPSolver" , "MOSEKMILPSolver" ,
	"SCIPMILPSolver" } ) {
  if( ! Solver::has_Solver( name ) ) {
   std::cout << name << ": not in the factory, skipped" << std::endl;
   continue;
   }

  Data d{ 1 , 1 , 0 , 4 , true };
  auto p = make_problem( d );
  auto s = Solver::new_Solver( name );
  s->set_par( Solver::intLogVerb , 0 );
  p.block->register_Solver( s );

  double value , fresh;
  if( ! solve( s , value ) ) {
   std::cout << name << ": the first solve failed -> KO" << std::endl;
   ++failures;
   }
  else
   for( const std::string change : { "coefficients" , "side" , "variable" } ) {
    if( change == "coefficients" ) {
     d.q1 = 4;
     d.l1 = 1;
     p.f->modify_term( 0 , d.l1 , d.q1 );
     }
    else
     if( change == "side" ) {
      d.r = 9;
      p.row->set_rhs( d.r );
      }
     else {
      d.both = false;
      p.f->remove_variable( 1 );
      }

    const bool ok_after = solve( s , value );
    const bool ok_fresh = fresh_value( name , d , fresh );
    const bool ok = ok_after && ok_fresh &&
     ( std::abs( value - fresh ) <=
       1e-6 * std::max( 1.0 , std::abs( fresh ) ) );

    std::cout << name << ", change of the " << change << ": "
	      << value << " against " << fresh << " loaded anew"
	      << ( ok ? " -> OK" : " -> KO" ) << std::endl;
    if( ! ok )
     ++failures;
    }

  p.block->unregister_Solver( s );
  delete s;
  }

 return( failures ? 1 : 0 );
 }

/*--------------------------------------------------------------------------*/
/*-------------------- End File test_quad_row_modify.cpp -------------------*/
/*--------------------------------------------------------------------------*/
