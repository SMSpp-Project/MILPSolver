/*--------------------------------------------------------------------------*/
/*--------------------- File test_partial_solution.cpp ---------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Checks that get_var_solution() and get_dual_solution() of each :MILPSolver
 * of the build write the part of the solution that a
 * SimpleConfiguration< std::vector< int > > asks for, and nothing else
 * [see MILPSolver::write_var_solution()]. The Block is
 *
 *     R    z , z <= 4
 *     |- A    a , a <= 3 , and the dynamic row a + a1 <= 4
 *     |  `- A1    a1 , a1 <= 2
 *     `- B    b , b <= 1 , and the dynamic column d in [ 0 , 6 ]
 *
 * with every other Variable in [ 0 , 10 ] by a BoxConstraint of its Block
 * and the Objective of each Block maximizing its Variable (d with
 * coefficient 2), so that every row and every bound has a dual value of its
 * own. After a first solve whose whole
 * solution is the reference, every value and every dual is overwritten with
 * a sentinel and then asked for again with:
 *
 * - { 0 }, sub-Block A and its sub-Block A1, the dynamic row included;
 * - { 1 }, sub-Block B, the dynamic column included;
 * - { -1 }, what R has itself, without its sub-Blocks;
 * - { 1 , -1 }, the two together;
 * - { }, nothing at all;
 * - a SimpleConfiguration< int >, which does not select and gives all;
 *
 * and each time what is asked for has to be the reference and what is not
 * has to be still the sentinel. Finally { 2 } and { -2 } have to throw
 * std::invalid_argument, there being no third sub-Block of R.
 *
 * The :MILPSolver are those of LPPar-all.txt that the build has, all
 * attached to the Block and checked one after the other; the test fails if
 * the build has none, or if any of them gives a wrong part.
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
#include <list>
#include <string>
#include <vector>

#include "AbstractBlock.h"

#include "CDASolver.h"

#include "test_common.h"

#include "FRealObjective.h"

#include "FRowConstraint.h"

#include "LinearFunction.h"

#include "OneVarConstraint.h"

/*--------------------------------------------------------------------------*/

using namespace SMSpp_di_unipi_it;

/*--------------------------------------------------------------------------*/

static constexpr double SENTINEL = -77;

/// the Block and what the test reads of it, each item with the Block it is in
struct Problem {
 AbstractBlock * root;
 std::vector< ColVariable * > vars;          ///< z , a , a1 , b , d
 std::vector< Constraint * > cons;           ///< rows and bounds
 std::vector< char > var_of , con_of;        ///< 'R' , 'A' , 'B'
 };

/*--------------------------------------------------------------------------*/
/// adds to b a Variable in [ 0 , 10 ], its row x <= cap and max x

static ColVariable * add_item( Problem & pb , AbstractBlock * b , char who ,
                               double cap , const std::string & name )
{
 auto x = new ColVariable;
 x->set_type( ColVariable::kContinuous );
 b->add_static_variable( * x , name.c_str() );

 auto box = new BoxConstraint;
 box->set_variable( x );
 box->set_lhs( 0 );
 box->set_rhs( 10 );
 b->add_static_constraint( * box , ( name + "_box" ).c_str() );

 auto row = new FRowConstraint;
 row->set_function( new LinearFunction(
                       LinearFunction::v_coeff_pair{ { x , 1.0 } } ) );
 row->set_lhs( - Inf< double >() );
 row->set_rhs( cap );
 b->add_static_constraint( * row , ( name + "_cap" ).c_str() );

 auto of = new FRealObjective( b , new LinearFunction(
                       LinearFunction::v_coeff_pair{ { x , 1.0 } } ) );
 of->set_sense( Objective::eMax );
 b->set_objective( of );

 pb.vars.push_back( x );
 pb.var_of.push_back( who );
 pb.cons.push_back( row );
 pb.cons.push_back( box );
 pb.con_of.push_back( who );
 pb.con_of.push_back( who );
 return( x );
 }

/*--------------------------------------------------------------------------*/

static Problem make_problem( void )
{
 Problem pb;
 pb.root = new AbstractBlock;
 auto A = new AbstractBlock( pb.root );
 auto A1 = new AbstractBlock( A );
 auto B = new AbstractBlock( pb.root );
 pb.root->add_nested_Block( A );
 pb.root->add_nested_Block( B );
 A->add_nested_Block( A1 );

 add_item( pb , pb.root , 'R' , 4 , "z" );
 auto a = add_item( pb , A , 'A' , 3 , "a" );
 auto a1 = add_item( pb , A1 , 'A' , 2 , "a1" );
 add_item( pb , B , 'B' , 1 , "b" );

 // the dynamic row of A , which binds: a + a1 <= 4 < 3 + 2
 auto rows = new std::list< FRowConstraint >( 1 );
 rows->front().set_function(
		       new LinearFunction( LinearFunction::v_coeff_pair{
                                             { a , 1.0 } , { a1 , 1.0 } } ) );
 rows->front().set_lhs( - Inf< double >() );
 rows->front().set_rhs( 4 );
 A->add_dynamic_constraint( * rows , "link" );
 pb.cons.push_back( & rows->front() );
 pb.con_of.push_back( 'A' );

 // the dynamic column of B , in the Objective of B with its own bound
 auto cols = new std::list< ColVariable >( 1 );
 auto d = & cols->front();
 d->set_type( ColVariable::kContinuous );
 B->add_dynamic_variable( * cols , "d" );
 auto boxes = new std::list< BoxConstraint >( 1 );
 auto dbox = & boxes->front();
 dbox->set_variable( d );
 dbox->set_lhs( 0 );
 dbox->set_rhs( 6 );
 B->add_dynamic_constraint( * boxes , "d_box" );
 static_cast< LinearFunction * >( static_cast< FRealObjective * >(
     B->get_objective() )->get_function() )->add_variable( d , 2.0 );
 pb.vars.push_back( d );
 pb.var_of.push_back( 'B' );
 pb.cons.push_back( dbox );
 pb.con_of.push_back( 'B' );

 return( pb );
 }

/*--------------------------------------------------------------------------*/

static double dual( Constraint * c )
{
 if( auto r = dynamic_cast< RowConstraint * >( c ) )
  return( r->get_dual() );
 return( std::nan( "" ) );
 }

/*--------------------------------------------------------------------------*/

static void overwrite( Problem & pb )
{
 for( auto v : pb.vars )
  v->set_value( SENTINEL );
 for( auto c : pb.cons )
  if( auto r = dynamic_cast< RowConstraint * >( c ) )
   r->set_dual( SENTINEL );
 }

/*--------------------------------------------------------------------------*/
/// true if what `asked` says is the reference, and all the rest the sentinel

static bool check( const Problem & pb , const std::string & asked ,
                   const std::vector< double > & x ,
                   const std::vector< double > & y )
{
 auto close = []( double u , double v ) {
  return( std::abs( u - v ) <= 1e-7 * ( 1 + std::abs( v ) ) );
  };

 bool ok = true;
 for( std::size_t i = 0 ; i < pb.vars.size() ; ++i ) {
  const bool in = asked.find( pb.var_of[ i ] ) != std::string::npos;
  if( ! close( pb.vars[ i ]->get_value() , in ? x[ i ] : SENTINEL ) )
   ok = false;
  }
 for( std::size_t i = 0 ; i < pb.cons.size() ; ++i ) {
  const bool in = asked.find( pb.con_of[ i ] ) != std::string::npos;
  if( ! close( dual( pb.cons[ i ] ) , in ? y[ i ] : SENTINEL ) )
   ok = false;
  }
 return( ok );
 }

/*--------------------------------------------------------------------------*/

int main( void )
{
 auto bsc = dynamic_cast< BlockSolverConfig * >(
		    Configuration::deserialize( "LPPar-all.txt" ) );
 Block::Index there = 0;
 if( bsc )
  for( Block::Index i = 0 ; i < bsc->num_ComputeConfig() ; ++i )
   if( Solver::has_Solver( bsc->get_SolverName( i ) ) )
    ++there;
 if( bsc && ( ! there ) ) {
  std::cout << "no :MILPSolver in this build" << std::endl;
  return( 1 );
  }
 keep_available_Solvers( bsc , "LPPar-all.txt" , there );

 std::vector< std::string > names;
 for( Block::Index i = 0 ; i < bsc->num_ComputeConfig() ; ++i )
  names.push_back( bsc->get_SolverName( i ) );

 auto pb = make_problem();
 bsc->apply( pb.root );
 bsc->clear();  // kept, clear()-ed, for the final cleanup

 int failures = 0;
 std::size_t i = 0;
 for( auto s : pb.root->get_registered_solvers() ) {
  const auto & name = names[ i++ ];
  auto solver = static_cast< CDASolver * >( s );

  bool ok = ( solver->compute() == Solver::kOK ) &&
            solver->has_var_solution() && solver->has_dual_solution();

  std::vector< double > x , y;
  if( ok ) {
   solver->get_var_solution();
   solver->get_dual_solution();
   for( auto v : pb.vars )
    x.push_back( v->get_value() );
   for( auto c : pb.cons )
    y.push_back( dual( c ) );
   }
  else
   std::cout << name << ": no solution to start from" << std::endl;

  const std::vector< std::pair< std::vector< int > , std::string > > parts =
   { { { 0 } , "A" } , { { 1 } , "B" } , { { -1 } , "R" } ,
     { { 1 , -1 } , "BR" } , { { } , "" } };

  for( const auto & [ idx , asked ] : parts ) {
   if( ! ok )
    break;
   SimpleConfiguration< std::vector< int > > solc( idx );
   overwrite( pb );
   solver->get_var_solution( & solc );
   solver->get_dual_solution( & solc );
   if( ! check( pb , asked , x , y ) ) {
    std::cout << name << ": wrong part for \"" << asked << "\""
              << std::endl;
    ok = false;
    }
   }

  if( ok ) {  // a Configuration that does not select gives all
   SimpleConfiguration< int > all( 1 );
   overwrite( pb );
   solver->get_var_solution( & all );
   solver->get_dual_solution( & all );
   if( ! check( pb , "RAB" , x , y ) ) {
    std::cout << name << ": not all with a SimpleConfiguration< int >"
              << std::endl;
    ok = false;
    }
   }

  for( int k : { 2 , -2 } ) {
   if( ! ok )
    break;
   SimpleConfiguration< std::vector< int > > solc( { k } );
   bool var_threw = false , dual_threw = false;
   try { solver->get_var_solution( & solc ); }
   catch( std::invalid_argument & ) { var_threw = true; }
   try { solver->get_dual_solution( & solc ); }
   catch( std::invalid_argument & ) { dual_threw = true; }
   if( ! ( var_threw && dual_threw ) ) {
    std::cout << name << ": index " << k << " does not throw" << std::endl;
    ok = false;
    }
   }

  std::cout << name << ( ok ? ": OK" : ": KO" ) << std::endl;
  if( ! ok )
   ++failures;
  }

 bsc->apply( pb.root );
 delete bsc;
 delete pb.root;

 return( failures ? 1 : 0 );
 }

/*--------------------------------------------------------------------------*/
/*-------------------- End File test_partial_solution.cpp ------------------*/
/*--------------------------------------------------------------------------*/
