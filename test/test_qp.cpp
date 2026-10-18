/*--------------------------------------------------------------------------*/
/*---------------------------- File test_qp.cpp ----------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Checks what each :MILPSolver does with a quadratic model, on small
 * problems built by hand whose optimum is known:
 *
 * - QP: the convex QP
 *
 *       min  sum_i a_i x_i^2 + sum_{i > j} a_ij x_i x_j + c^T x + 7
 *       s.t. x_0 + x_1 + x_4 <= -2              ( mu_0 = -2   )
 *            x_1 + x_2 >= -1.5                  ( mu_1 =  1   )
 *            x_0 + x_2 + x_3 = 4.5              ( mu_2 =  0.5 )
 *            -10 <= x_1 - x_3 + x_4 <= -6       ( mu_3 = -1   )
 *            -5 <= x_0 - x_3 <= 5               ( inactive    )
 *            x_2 + x_4 + x_5 <= 10              ( inactive    )
 *            1 <= x_0 <= 5 ( nu_0 = 3 ) ,  -1 <= x_2 <= 2 ,  x_3 >= 0 ,
 *            x_4 <= 0 ,  x_1 and x_5 free
 *
 *   with a = ( 1.5 , 1.5 , 1 , 1 , 1 , 1 ), built backwards from the
 *   solution x* = ( 1 , -2 , 0.5 , 3 , -1 , 2 ) and the multipliers above
 *   (the derivatives of the optimal value in the side of the row or of the
 *   bound), c being - H x* + A^T mu + nu. It is solved with the off-diagonal
 *   terms a_10 = 1, a_21 = -1, a_32 = 0.5, a_43 = 1 and a_50 = 1 in a
 *   QuadFunction, where c = ( -1.5 , 3.5 , -3 , -3.75 , -4 , -5 ) and the
 *   value is 0.25, and without them in a DQuadFunction, where
 *   c = ( -1.5 , 4 , 0.5 , -4.5 , -1 , -4 ) and the value is -1, each also
 *   as the maximization of the opposite (concave) objective. The value and
 *   x* are checked, and so are the duals of the rows and of the bounds when
 *   the Solver gives them: by the convention of
 *   MILPSolver::write_dual_solution() the dual of a row or of a bound is
 *   - mu or - nu in a minimization and mu or nu in the maximization of the
 *   opposite objective, i.e., grad f( x* ) + A^T y + d = 0 for the objective
 *   f as it is given, and this stationarity is checked as well.
 *
 * - MIQP: min x_0^2 + x_1^2 + y^2 + 0.5 x_0 x_1 - 0.8 x_0 - 3.2 x_1
 *   - 0.6 y + 1 s.t. x_0 + x_1 + y <= 1.5, x_0 and x_1 integer in
 *   [ -3 , 3 ], -5 <= y <= 5, whose optimum is -1.29 at ( 0 , 1 , 0.3 )
 *   while that of its continuous relaxation is about -1.5885.
 *
 * - Modifications: after a first solve of the QP with the QuadFunction a
 *   diagonal coefficient and then a linear coefficient of the objective are
 *   changed, the side of an inactive row is changed so that it binds, a row
 *   is added and then removed, and a variable is fixed; after each of them
 *   the Solver solves again, and its value is compared with that of a fresh
 *   Solver of the same type loaded with the modified problem.
 *
 * Each model is in the only sub-Block of a root Block with an empty
 * objective, as PIPSMILPSolver needs at least one nested Block.
 *
 * Each :MILPSolver among CPXMILPSolver, GRBMILPSolver, SCIPMILPSolver,
 * HiGHSMILPSolver, DAQPMILPSolver, QPOASESMILPSolver, ClarabelMILPSolver and
 * PIPSMILPSolver that is not in the factory is skipped, and the test exits
 * with 77 if none of them is; a name on the command line restricts the run to
 * that Solver. A Solver is run on each case it is meant to handle, and on each
 * of the others it has to refuse the model cleanly, by an exception or by a
 * status not below kError, rather than give a wrong answer. A Solver that does
 * not give the duals of a case is not a failure, but one whose duals are wrong
 * is. SCIPMILPSolver solves a QP through the outer approximation of an
 * auxiliary constraint, which leaves the value exact to its feasibility
 * tolerance but the solution only to about the square root of it, hence its
 * solutions are compared to 1e-3 rather than 1e-5.
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
#include <functional>
#include <iostream>
#include <list>
#include <sstream>
#include <vector>

#include "AbstractBlock.h"

#include "FRealObjective.h"

#include "FRowConstraint.h"

#include "LinearFunction.h"

#include "MILPSolver.h"

#include "OneVarConstraint.h"

#include "QuadFunction.h"

/*--------------------------------------------------------------------------*/
/*-------------------------------- USING -----------------------------------*/
/*--------------------------------------------------------------------------*/

using namespace SMSpp_di_unipi_it;

/*--------------------------------------------------------------------------*/
/*----------------------------- CONSTANTS ----------------------------------*/
/*--------------------------------------------------------------------------*/

/// what each :MILPSolver is meant to solve
struct Ability {
 std::string name;  ///< the name in the factory
 bool qp;           ///< a convex QP with linear rows
 bool miqp;         ///< a convex QP with integer variables
 double xtol;       ///< the tolerance on the solution
 };

static const std::vector< Ability > Solvers = {
 { "CPXMILPSolver"   , true  , true  , 1e-5 } ,
 { "GRBMILPSolver"   , true  , true  , 1e-5 } ,
 { "SCIPMILPSolver"  , true  , true  , 1e-3 } ,
 { "HiGHSMILPSolver" , true  , false , 1e-5 } ,
 { "DAQPMILPSolver"  , true  , false , 1e-5 } ,
 { "QPOASESMILPSolver" , true , false , 1e-5 } ,
 { "ClarabelMILPSolver" , true , false , 1e-5 } ,
 { "PIPSMILPSolver"  , true  , false , 1e-5 } };

/// the tolerance on the values and on the solutions

static constexpr double Eps = 1e-5;

/// the data of the QP of the file comment

static constexpr int NV = 6;
static constexpr int NR = 6;

static const double QuadC[ NV ] = { 1.5 , 1.5 , 1 , 1 , 1 , 1 };
static const double LinQ[ NV ] = { -1.5 , 3.5 , -3 , -3.75 , -4 , -5 };
static const double LinD[ NV ] = { -1.5 , 4 , 0.5 , -4.5 , -1 , -4 };
static const QuadFunction::v_off_diag_term OffDiag =
 { { 1 , 0 , 1 } , { 2 , 1 , -1 } , { 3 , 2 , 0.5 } , { 4 , 3 , 1 } ,
   { 5 , 0 , 1 } };
static const double Const = 7;

static const double XStar[ NV ] = { 1 , -2 , 0.5 , 3 , -1 , 2 };
static const double ValQ = 0.25;
static const double ValD = -1;

static const double A[ NR ][ NV ] = { { 1 , 1 , 0 ,  0 , 1 , 0 } ,
				      { 0 , 1 , 1 ,  0 , 0 , 0 } ,
				      { 1 , 0 , 1 ,  1 , 0 , 0 } ,
				      { 0 , 1 , 0 , -1 , 1 , 0 } ,
				      { 1 , 0 , 0 , -1 , 0 , 0 } ,
				      { 0 , 0 , 1 ,  0 , 1 , 1 } };
static const double Lhs[ NR ] = { -Inf< double >() , -1.5 , 4.5 , -10 , -5 ,
				  -Inf< double >() };
static const double Rhs[ NR ] = { -2 , Inf< double >() , 4.5 , -6 , 5 , 10 };
static const double Mu[ NR ] = { -2 , 1 , 0.5 , -1 , 0 , 0 };
static const double Nu[ NV ] = { 3 , 0 , 0 , 0 , 0 , 0 };

/*--------------------------------------------------------------------------*/
/*------------------------------- TYPES ------------------------------------*/
/*--------------------------------------------------------------------------*/

/// the QP of the file comment and the pieces of it that are looked at

struct QP {
 AbstractBlock * block;                   ///< the (empty) root Block
 AbstractBlock * sub;                     ///< its sub-Block, with the model
 std::vector< ColVariable > * x;
 std::vector< FRowConstraint > * rows;
 OneVarConstraint * bound[ NV ];          ///< nullptr for a free variable
 DQuadFunction * f;
 std::list< FRowConstraint > * dyn;       ///< the dynamic rows
 bool nonsep;
 double sign;                             ///< 1 for min, -1 for max
 };

/*--------------------------------------------------------------------------*/
/*------------------------------ FUNCTIONS ---------------------------------*/
/*--------------------------------------------------------------------------*/

static void make_row( std::vector< ColVariable > & x ,
				const double * a , double lhs , double rhs ,
				FRowConstraint & row )
{
 LinearFunction::v_coeff_pair coeffs;
 for( int j = 0 ; j < NV ; ++j )
  if( a[ j ] != 0 )
   coeffs.emplace_back( & x[ j ] , a[ j ] );
 row.set_function( new LinearFunction( std::move( coeffs ) ) , eNoMod );
 row.set_lhs( lhs , eNoMod );
 row.set_rhs( rhs , eNoMod );
 }

/*--------------------------------------------------------------------------*/

static QP make_qp( bool nonsep , Objective::of_type sense )
{
 QP p;
 p.nonsep = nonsep;
 p.sign = ( sense == Objective::eMin ) ? 1 : -1;
 p.block = new AbstractBlock;
 p.sub = new AbstractBlock( p.block );

 p.x = new std::vector< ColVariable >( NV );
 p.sub->add_static_variable( * p.x , "x" );

 p.rows = new std::vector< FRowConstraint >( NR );
 for( int i = 0 ; i < NR ; ++i )
  make_row( * p.x , A[ i ] , Lhs[ i ] , Rhs[ i ] , ( * p.rows )[ i ] );
 p.sub->add_static_constraint( * p.rows , "rows" );

 // finite and infinite bounds, by each kind of OneVarConstraint
 auto box = new std::vector< BoxConstraint >( 2 );
 ( * box )[ 0 ].set_variable( & ( * p.x )[ 0 ] , eNoMod );
 ( * box )[ 0 ].set_lhs( 1 , eNoMod );
 ( * box )[ 0 ].set_rhs( 5 , eNoMod );
 ( * box )[ 1 ].set_variable( & ( * p.x )[ 2 ] , eNoMod );
 ( * box )[ 1 ].set_lhs( -1 , eNoMod );
 ( * box )[ 1 ].set_rhs( 2 , eNoMod );
 p.sub->add_static_constraint( * box , "box" );

 auto lb = new LBConstraint;
 lb->set_variable( & ( * p.x )[ 3 ] , eNoMod );
 lb->set_lhs( 0 , eNoMod );
 p.sub->add_static_constraint( * lb , "lb" );

 auto ub = new UBConstraint;
 ub->set_variable( & ( * p.x )[ 4 ] , eNoMod );
 ub->set_rhs( 0 , eNoMod );
 p.sub->add_static_constraint( * ub , "ub" );

 p.bound[ 0 ] = & ( * box )[ 0 ];
 p.bound[ 1 ] = nullptr;
 p.bound[ 2 ] = & ( * box )[ 1 ];
 p.bound[ 3 ] = lb;
 p.bound[ 4 ] = ub;
 p.bound[ 5 ] = nullptr;

 p.dyn = new std::list< FRowConstraint >;
 p.sub->add_dynamic_constraint( * p.dyn , "dyn" );

 const double * lin = nonsep ? LinQ : LinD;
 QuadFunction::v_coeff_triple diag;
 for( int i = 0 ; i < NV ; ++i )
  diag.emplace_back( & ( * p.x )[ i ] , p.sign * lin[ i ] ,
		     p.sign * QuadC[ i ] );

 if( nonsep ) {
  auto offdiag = OffDiag;
  for( auto & t : offdiag )
   std::get< 2 >( t ) *= p.sign;
  p.f = new QuadFunction( std::move( diag ) , std::move( offdiag ) ,
			  p.sign * Const );
  }
 else
  p.f = new DQuadFunction( std::move( diag ) , p.sign * Const );

 auto obj = new FRealObjective( p.sub , p.f );
 obj->set_sense( sense , eNoMod );
 p.sub->set_objective( obj , eNoMod );

 p.block->add_nested_Block( p.sub );
 auto robj = new FRealObjective( p.block , new LinearFunction() );
 robj->set_sense( sense , eNoMod );
 p.block->set_objective( robj , eNoMod );

 return( p );
 }

/*--------------------------------------------------------------------------*/

static Solver * new_solver( const std::string & name , bool duals = false )
{
 auto s = Solver::new_Solver( name );
 s->set_par( Solver::intLogVerb , 0 );
 s->set_par( MILPSolver::intRelaxIntVars , 0 );
 // SCIPMILPSolver gives the duals only if asked beforehand
 if( duals && ( name == "SCIPMILPSolver" ) )
  (void) s->set_par( std::string( "intComputeDuals" ) , 1 );
 return( s );
 }

/*--------------------------------------------------------------------------*/

/// registers (if asked) and computes: false and the reason if it fails
/** An exception, a status out of [ kOK , kError ) and a value that is not
 * finite are all failures, \p why saying which. */

static bool run( Solver * s , Block * b , bool reg , double & value ,
		 std::string & why )
{
 try {
  if( reg )
   b->register_Solver( s );
  const auto status = s->compute();
  if( ( status < Solver::kOK ) || ( status >= Solver::kError ) ) {
   why = "status " + std::to_string( status );
   return( false );
   }
  value = s->get_var_value();
  if( ! std::isfinite( value ) ) {
   why = "value " + std::to_string( value );
   return( false );
   }
  s->get_var_solution();
  return( true );
  }
 catch( std::exception & e ) {
  why = std::string( "exception: " ) + e.what();
  return( false );
  }
 }

/*--------------------------------------------------------------------------*/

static bool close( double a , double b )
{
 return( std::abs( a - b ) <= Eps * std::max( 1.0 , std::abs( b ) ) );
 }

/*--------------------------------------------------------------------------*/

/// the case that a Solver is not meant to handle: has it refused it?
/** It has if run() failed (an exception or a status from kError on);
 * the answer is wrong otherwise, whatever it is. */

static bool refused( const std::string & name , const std::string & what ,
		     Solver * s , Block * b )
{
 double value;
 std::string why;
 if( run( s , b , true , value , why ) ) {
  std::cout << name << ", " << what << ": not refused, value = " << value
	    << " -> KO" << std::endl;
  return( false );
  }
 std::cout << name << ", " << what << ": refused (" << why << ") -> OK"
	   << std::endl;
 return( true );
 }

/*--------------------------------------------------------------------------*/

/// QP: value, solution, duals and stationarity

static int test_qp( const Ability & ab , bool nonsep ,
		    Objective::of_type sense )
{
 std::ostringstream what;
 what << "QP " << ( nonsep ? "QuadFunction" : "DQuadFunction" )
      << ( sense == Objective::eMin ? " min" : " max" );

 auto p = make_qp( nonsep , sense );
 auto s = new_solver( ab.name , true );
 int failures = 0;

 if( ! ab.qp )
  failures += ! refused( ab.name , what.str() , s , p.block );
 else {
  double value;
  std::string why;
  std::ostringstream out;
  bool ok = run( s , p.block , true , value , why );
  if( ! ok )
   out << why;
  else {
   const double exp_val = p.sign * ( nonsep ? ValQ : ValD );
   if( ! close( value , exp_val ) ) {
    out << "value " << value << " != " << exp_val << " ";
    ok = false;
    }
   for( int j = 0 ; j < NV ; ++j )
    if( std::abs( ( * p.x )[ j ].get_value() - XStar[ j ] ) > ab.xtol ) {
     out << "x_" << j << " = " << ( * p.x )[ j ].get_value() << " != "
	 << XStar[ j ] << " ";
     ok = false;
     }
   if( ok )
    out << "value = " << value;
   }

  std::cout << ab.name << ", " << what.str() << ": " << out.str()
	    << ( ok ? " -> OK" : " -> KO" ) << std::endl;
  failures += ! ok;

  // the duals, if the Solver gives them
  if( ok ) {
   auto cda = dynamic_cast< CDASolver * >( s );
   bool has = false;
   try {
    has = cda && cda->has_dual_solution();
    if( has )
     cda->get_dual_solution();
    }
   catch( std::exception & e ) {
    std::cout << ab.name << ", " << what.str() << ", duals: exception: "
	      << e.what() << " -> KO" << std::endl;
    ++failures;
    has = false;
    ok = false;
    }

   if( ! has ) {
    if( ok )
     std::cout << ab.name << ", " << what.str()
	       << ", duals: not given" << std::endl;
    }
   else {
    std::ostringstream dout;
    bool dok = true;

    // the expected duals: - mu in a min, mu in the max of - f
    for( int i = 0 ; i < NR ; ++i ) {
     const double y = ( * p.rows )[ i ].get_dual();
     if( std::abs( y + p.sign * Mu[ i ] ) > Eps ) {
      dout << "y_" << i << " = " << y << " != " << - p.sign * Mu[ i ] << " ";
      dok = false;
      }
     }
    for( int j = 0 ; j < NV ; ++j ) {
     const double d = p.bound[ j ] ? p.bound[ j ]->get_dual() : 0;
     if( std::abs( d + p.sign * Nu[ j ] ) > Eps ) {
      dout << "d_" << j << " = " << d << " != " << - p.sign * Nu[ j ] << " ";
      dok = false;
      }
     }

    // stationarity: grad f( x ) + A^T y + d = 0 at the given solution
    double x[ NV ] , g[ NV ];
    for( int j = 0 ; j < NV ; ++j )
     x[ j ] = ( * p.x )[ j ].get_value();
    const double * lin = nonsep ? LinQ : LinD;
    for( int j = 0 ; j < NV ; ++j )
     g[ j ] = p.sign * ( 2 * QuadC[ j ] * x[ j ] + lin[ j ] );
    if( nonsep )
     for( const auto & [ i , j , v ] : OffDiag ) {
      g[ i ] += p.sign * v * x[ j ];
      g[ j ] += p.sign * v * x[ i ];
      }
    double res = 0;
    for( int j = 0 ; j < NV ; ++j ) {
     double r = g[ j ] + ( p.bound[ j ] ? p.bound[ j ]->get_dual() : 0 );
     for( int i = 0 ; i < NR ; ++i )
      r += A[ i ][ j ] * ( * p.rows )[ i ].get_dual();
     res = std::max( res , std::abs( r ) );
     }
    if( res > Eps ) {
     dout << "stationarity residual " << res << " ";
     dok = false;
     }
    if( dok )
     dout << "KKT residual " << res;

    std::cout << ab.name << ", " << what.str() << ", duals: " << dout.str()
	      << ( dok ? " -> OK" : " -> KO" ) << std::endl;
    failures += ! dok;
    }
   }
  }

 try {
  p.block->unregister_Solver( s );
  }
 catch( ... ) {}
 delete s;
 delete p.block;
 return( failures );
 }

/*--------------------------------------------------------------------------*/

/// MIQP: the value against the known optimum

static int test_miqp( const Ability & ab )
{
 auto block = new AbstractBlock;
 auto sub = new AbstractBlock( block );
 auto x = new std::vector< ColVariable >( 3 );
 sub->add_static_variable( * x , "x" );
 ( * x )[ 0 ].set_type( ColVariable::kInteger , eNoMod );
 ( * x )[ 1 ].set_type( ColVariable::kInteger , eNoMod );

 auto box = new std::vector< BoxConstraint >( 3 );
 for( int i = 0 ; i < 3 ; ++i ) {
  ( * box )[ i ].set_variable( & ( * x )[ i ] , eNoMod );
  ( * box )[ i ].set_lhs( i < 2 ? -3 : -5 , eNoMod );
  ( * box )[ i ].set_rhs( i < 2 ? 3 : 5 , eNoMod );
  }
 sub->add_static_constraint( * box , "box" );

 auto row = new FRowConstraint;
 const double a[ NV ] = { 1 , 1 , 1 , 0 , 0 , 0 };
 make_row( * x , a , -Inf< double >() , 1.5 , * row );
 sub->add_static_constraint( * row , "row" );

 QuadFunction::v_coeff_triple diag = { { & ( * x )[ 0 ] , -0.8 , 1 } ,
				       { & ( * x )[ 1 ] , -3.2 , 1 } ,
				       { & ( * x )[ 2 ] , -0.6 , 1 } };
 QuadFunction::v_off_diag_term offdiag = { { 1 , 0 , 0.5 } };
 auto f = new QuadFunction( std::move( diag ) , std::move( offdiag ) , 1 );
 auto obj = new FRealObjective( sub , f );
 obj->set_sense( Objective::eMin , eNoMod );
 sub->set_objective( obj , eNoMod );

 block->add_nested_Block( sub );
 auto robj = new FRealObjective( block , new LinearFunction() );
 robj->set_sense( Objective::eMin , eNoMod );
 block->set_objective( robj , eNoMod );

 auto s = new_solver( ab.name );
 int failures = 0;

 if( ! ab.miqp )
  failures += ! refused( ab.name , "MIQP" , s , block );
 else {
  double value;
  std::string why;
  bool ok = run( s , block , true , value , why );
  if( ok )
   ok = close( value , -1.29 ) &&
        ( std::abs( ( * x )[ 0 ].get_value() ) <= ab.xtol ) &&
        ( std::abs( ( * x )[ 1 ].get_value() - 1 ) <= ab.xtol ) &&
        ( std::abs( ( * x )[ 2 ].get_value() - 0.3 ) <= ab.xtol );
  std::cout << ab.name << ", MIQP: ";
  if( why.empty() )
   std::cout << "value = " << value << " (-1.29 expected), x = ( "
	     << ( * x )[ 0 ].get_value() << " , " << ( * x )[ 1 ].get_value()
	     << " , " << ( * x )[ 2 ].get_value() << " )";
  else
   std::cout << why;
  std::cout << ( ok ? " -> OK" : " -> KO" ) << std::endl;
  failures += ! ok;
  }

 try {
  block->unregister_Solver( s );
  }
 catch( ... ) {}
 delete s;
 delete block;
 return( failures );
 }

/*--------------------------------------------------------------------------*/

/// Modifications: each one solved again and compared with a fresh Solver

static int test_mods( const Ability & ab )
{
 if( ! ab.qp )
  return( 0 );

 auto p = make_qp( true , Objective::eMin );
 auto s = new_solver( ab.name );
 int failures = 0;

 double value;
 std::string why;
 if( ! run( s , p.block , true , value , why ) ) {
  std::cout << ab.name << ", Modifications: the first solve failed ("
	    << why << ") -> KO" << std::endl;
  delete s;
  delete p.block;
  return( 1 );
  }

 auto & x = * p.x;

 const std::vector< std::pair< std::string , std::function< void() > > >
  steps = {
  { "diagonal coefficient of x_5 from 1 to 3" ,
    [ & ]() { p.f->modify_term( 5 , LinQ[ 5 ] , 3 ); } } ,
  { "linear coefficient of x_2 from -3 to 2" ,
    [ & ]() { p.f->modify_linear_coefficient( 2 , 2 ); } } ,
  { "rhs of row 5 from 10 to -2" ,
    [ & ]() { ( * p.rows )[ 5 ].set_rhs( -2 ); } } ,
  { "row x_5 >= 1 added" ,
    [ & ]() {
     std::list< FRowConstraint > nr( 1 );
     const double a[ NV ] = { 0 , 0 , 0 , 0 , 0 , 1 };
     make_row( x , a , 1 , Inf< double >() , nr.front() );
     p.sub->add_dynamic_constraints( * p.dyn , nr );
     } } ,
  { "row x_5 >= 1 removed" ,
    [ & ]() {
     p.sub->remove_dynamic_constraint( * p.dyn , p.dyn->begin() );
     } } ,
  { "x_3 fixed to 2.5" ,
    [ & ]() { x[ 3 ].set_value( 2.5 ); x[ 3 ].is_fixed( true ); } } };

 for( const auto & [ what , step ] : steps ) {
  double after = 0 , fresh = 0;
  std::string why_a , why_f;
  bool ok_a , ok_f;
  try {
   step();
   ok_a = run( s , p.block , false , after , why_a );
   }
  catch( std::exception & e ) {
   ok_a = false;
   why_a = std::string( "exception: " ) + e.what();
   }

  auto s2 = new_solver( ab.name );
  ok_f = run( s2 , p.block , true , fresh , why_f );
  try {
   p.block->unregister_Solver( s2 );
   }
  catch( ... ) {}
  delete s2;

  const bool ok = ok_a && ok_f && close( after , fresh );
  std::cout << ab.name << ", Modifications, " << what << ": ";
  if( ok_a )
   std::cout << "value = " << after;
  else
   std::cout << why_a;
  std::cout << ", loaded anew = ";
  if( ok_f )
   std::cout << fresh;
  else
   std::cout << why_f;
  std::cout << ( ok ? " -> OK" : " -> KO" ) << std::endl;
  failures += ! ok;
  }

 try {
  p.block->unregister_Solver( s );
  }
 catch( ... ) {}
 delete s;
 delete p.block;
 return( failures );
 }

/*--------------------------------------------------------------------------*/
/*-------------------------------- main() ----------------------------------*/
/*--------------------------------------------------------------------------*/

int main( int argc , char ** argv )
{
 int failures = 0;
 bool any = false;

 for( const auto & ab : Solvers ) {
  // a name on the command line restricts the run to that Solver
  if( ( argc > 1 ) && ( ab.name != argv[ 1 ] ) )
   continue;
  if( ! Solver::has_Solver( ab.name ) ) {
   std::cout << ab.name << ": not in the factory, skipped" << std::endl;
   continue;
   }
  any = true;

  for( bool nonsep : { true , false } )
   for( auto sense : { Objective::eMin , Objective::eMax } )
    failures += test_qp( ab , nonsep , sense );

  failures += test_miqp( ab );

  failures += test_mods( ab );
  }

 if( ! any ) {
  std::cout << "no :MILPSolver in this build: nothing to do" << std::endl;
  return( 77 );
  }

 std::cout << ( failures ? "Shit happened!!" : "All tests passed!!" )
	   << std::endl;
 return( failures ? 1 : 0 );
 }

/*--------------------------------------------------------------------------*/
/*------------------------- End File test_qp.cpp ---------------------------*/
/*--------------------------------------------------------------------------*/
