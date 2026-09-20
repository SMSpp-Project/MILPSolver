/*--------------------------------------------------------------------------*/
/*-------------------------- File test_groups.cpp --------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Main for testing the :MILPSolver on groups of Variable and Constraint made
 * of several arrays
 *
 * A group of static Variable or Constraint of a Block need not be a single
 * array: it can be a std::vector of std::vector or a boost::multi_array of
 * std::vector (and a group of dynamic ones a std::vector or a
 * boost::multi_array of std::list), whose elements are not contiguous in
 * memory as a whole even when those of each inner container are. A :MILPSolver
 * maps the static ones to columns and rows by address, one run of contiguous
 * elements at a time, and the dynamic ones one element at a time; what is
 * checked here is that the model it builds is the same whatever the shape of
 * the groups.
 *
 * The program is
 *
 *     min  x0 + 2 x1 + 3 x2  s.t.  x0 + x1 >= 1 ,  x1 + x2 >= 2 ,
 *                                  x0 + x2 >= 3 ,  x0 , x1 , x2 >= 0
 *
 * whose optimum is 7, at x = ( 2 , 0 , 1 ), with duals ( 0 , 2 , 1 ) (up to
 * the sign convention of MILPSolver::write_dual_solution(), which is why their
 * absolute value is compared). It is written in four shapes, the same three
 * Variable and three Constraint being always split as [ [ 0 , 1 ] , [ 2 ] ]:
 *
 * - flat: a std::vector of three ColVariable and one of three FRowConstraint;
 *
 * - irregular: the ColVariable in a boost::multi_array of std::vector, the
 *   FRowConstraint in a std::vector of std::vector, both static;
 *
 * - vector of list: both in a std::vector of std::list, dynamic;
 *
 * - multi_array of list: both in a boost::multi_array of std::list, dynamic.
 *
 * Every registered :MILPSolver solves each shape, and the optimal value, the
 * duals and the feasibility of the primal solution are checked. A
 * boost::multi_array is always given two dimensions, since a group that is a
 * one-dimensional boost::multi_array is not recognised (a std::vector is the
 * one-dimensional group).
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
#include <iostream>
#include <list>
#include <vector>

#include "AbstractBlock.h"

#include "BlockSolverConfig.h"

#include "FRealObjective.h"

#include "FRowConstraint.h"

#include "LinearFunction.h"

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

/// the shapes of the groups [see the file comment]

enum Shape { eFlat = 0 , eIrregular , eVectorOfList , eMultiArrayOfList };

static const char * const ShapeNames[] =
 { "flat" , "irregular" , "vector of list" , "multi_array of list" };

/// the tolerance of every comparison

static constexpr double Eps = 1e-7;

/*--------------------------------------------------------------------------*/
/*------------------------------ FUNCTIONS ---------------------------------*/
/*--------------------------------------------------------------------------*/

/// solves the program of the file comment with the given shape and Solver
/** Returns true if the optimal value, the duals and the primal solution are
 * right, printing what is wrong otherwise. */

static bool solve( Shape shape , const std::string & solver_name )
{
 auto ab = new AbstractBlock();

 // the Variable and Constraint, split as [ [ 0 , 1 ] , [ 2 ] ] in a container
 // of the given shape, which the AbstractBlock takes ownership of
 ColVariable * x[ 3 ];
 FRowConstraint * r[ 3 ];

 switch( shape ) {
  case( eFlat ): {
   auto v = new std::vector< ColVariable >( 3 );
   auto c = new std::vector< FRowConstraint >( 3 );
   for( int i = 0 ; i < 3 ; ++i ) {
    x[ i ] = & (*v)[ i ];
    r[ i ] = & (*c)[ i ];
    }
   ab->add_static_variable( *v , "x" );
   ab->add_static_constraint( *c , "r" );
   break;
   }
  case( eIrregular ): {
   auto v = new boost::multi_array< std::vector< ColVariable > , 2 >(
                                                     boost::extents[ 1 ][ 2 ] );
   (*v)[ 0 ][ 0 ].resize( 2 );
   (*v)[ 0 ][ 1 ].resize( 1 );
   auto c = new std::vector< std::vector< FRowConstraint > >( 2 );
   (*c)[ 0 ].resize( 2 );
   (*c)[ 1 ].resize( 1 );
   x[ 0 ] = & (*v)[ 0 ][ 0 ][ 0 ];
   x[ 1 ] = & (*v)[ 0 ][ 0 ][ 1 ];
   x[ 2 ] = & (*v)[ 0 ][ 1 ][ 0 ];
   r[ 0 ] = & (*c)[ 0 ][ 0 ];
   r[ 1 ] = & (*c)[ 0 ][ 1 ];
   r[ 2 ] = & (*c)[ 1 ][ 0 ];
   ab->add_static_variable( *v , "x" );
   ab->add_static_constraint( *c , "r" );
   break;
   }
  case( eVectorOfList ): {
   auto v = new std::vector< std::list< ColVariable > >( 2 );
   (*v)[ 0 ].resize( 2 );
   (*v)[ 1 ].resize( 1 );
   auto c = new std::vector< std::list< FRowConstraint > >( 2 );
   (*c)[ 0 ].resize( 2 );
   (*c)[ 1 ].resize( 1 );
   x[ 0 ] = & (*v)[ 0 ].front();
   x[ 1 ] = & (*v)[ 0 ].back();
   x[ 2 ] = & (*v)[ 1 ].front();
   r[ 0 ] = & (*c)[ 0 ].front();
   r[ 1 ] = & (*c)[ 0 ].back();
   r[ 2 ] = & (*c)[ 1 ].front();
   ab->add_dynamic_variable( *v , "x" );
   ab->add_dynamic_constraint( *c , "r" );
   break;
   }
  default: {
   auto v = new boost::multi_array< std::list< ColVariable > , 2 >(
                                                     boost::extents[ 1 ][ 2 ] );
   (*v)[ 0 ][ 0 ].resize( 2 );
   (*v)[ 0 ][ 1 ].resize( 1 );
   auto c = new boost::multi_array< std::list< FRowConstraint > , 2 >(
                                                     boost::extents[ 2 ][ 1 ] );
   (*c)[ 0 ][ 0 ].resize( 2 );
   (*c)[ 1 ][ 0 ].resize( 1 );
   x[ 0 ] = & (*v)[ 0 ][ 0 ].front();
   x[ 1 ] = & (*v)[ 0 ][ 0 ].back();
   x[ 2 ] = & (*v)[ 0 ][ 1 ].front();
   r[ 0 ] = & (*c)[ 0 ][ 0 ].front();
   r[ 1 ] = & (*c)[ 0 ][ 0 ].back();
   r[ 2 ] = & (*c)[ 1 ][ 0 ].front();
   ab->add_dynamic_variable( *v , "x" );
   ab->add_dynamic_constraint( *c , "r" );
   }
  }

 for( int i = 0 ; i < 3 ; ++i )
  x[ i ]->set_type( ColVariable::kNonNegative , eNoMod );

 // the two Variable of each row, and its left-hand side
 const int rows[ 3 ][ 2 ] = { { 0 , 1 } , { 1 , 2 } , { 0 , 2 } };
 const double lhs[ 3 ] = { 1 , 2 , 3 };
 for( int i = 0 ; i < 3 ; ++i ) {
  r[ i ]->set_function( new LinearFunction(
   { std::make_pair( x[ rows[ i ][ 0 ] ] , 1.0 ) ,
     std::make_pair( x[ rows[ i ][ 1 ] ] , 1.0 ) } ) , eNoMod );
  r[ i ]->set_lhs( lhs[ i ] , eNoMod );
  r[ i ]->set_rhs( Inf< double >() , eNoMod );
  }

 auto obj = new FRealObjective( ab , new LinearFunction(
  { std::make_pair( x[ 0 ] , 1.0 ) , std::make_pair( x[ 1 ] , 2.0 ) ,
    std::make_pair( x[ 2 ] , 3.0 ) } ) );
 obj->set_sense( Objective::eMin , eNoMod );
 ab->set_objective( obj , eNoMod );

 auto bsc = new BlockSolverConfig( 1 );
 bsc->add_ComputeConfig( std::string( solver_name ) , nullptr );
 bsc->apply( ab );
 auto solver = static_cast< CDASolver * >(
                                    ab->get_registered_solvers().front() );
 solver->set_par( MILPSolver::intLogVerb , 0 );

 const auto status = solver->compute( false );
 bool ok = true;

 std::cout << solver_name << ", " << ShapeNames[ shape ] << ": ";
 if( status != Solver::kOK ) {
  std::cout << "status " << status << " -> Error" << std::endl;
  ok = false;
  }
 else {
  solver->get_var_solution();
  solver->get_dual_solution();

  if( std::abs( solver->get_var_value() - 7 ) > Eps ) {
   std::cout << "value " << solver->get_var_value() << " != 7 ";
   ok = false;
   }

  const double duals[ 3 ] = { 0 , 2 , 1 };
  for( int i = 0 ; i < 3 ; ++i )
   if( std::abs( std::abs( r[ i ]->get_dual() ) - duals[ i ] ) > Eps ) {
    std::cout << "dual of row " << i << " " << r[ i ]->get_dual()
              << " != " << duals[ i ] << " ";
    ok = false;
    }

  for( int i = 0 ; i < 3 ; ++i )
   if( x[ rows[ i ][ 0 ] ]->get_value() + x[ rows[ i ][ 1 ] ]->get_value()
       < lhs[ i ] - Eps ) {
    std::cout << "row " << i << " violated ";
    ok = false;
    }

  std::cout << ( ok ? "-> OK" : "-> Error" ) << std::endl;
  }

 bsc->clear();
 bsc->apply( ab );
 delete bsc;
 delete ab;

 return( ok );
 }

/*--------------------------------------------------------------------------*/
/*-------------------------------- main() ----------------------------------*/
/*--------------------------------------------------------------------------*/

int main( void )
{
 bool all_passed = true;
 bool any_solver = false;

 for( const auto & name : SolverNames ) {
  // skip the :MILPSolver that are not in the build
  auto probe = Solver::new_Solver( name );
  if( ! probe )
   continue;
  delete probe;
  any_solver = true;

  for( int shape = eFlat ; shape <= eMultiArrayOfList ; ++shape )
   if( ! solve( Shape( shape ) , name ) )
    all_passed = false;
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
/*-------------------------- File test_groups.cpp --------------------------*/
/*--------------------------------------------------------------------------*/
