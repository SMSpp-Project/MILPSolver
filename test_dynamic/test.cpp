/*--------------------------------------------------------------------------*/
/*-------------------------- File test.cpp ---------------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Main for testing Linear Programs
 *
 * A "random" Linear Program is constructed that represents a polyhedral
 * function, either convex or concave (the Objective contains a single
 * variable with coefficient 1 that is minimised in the convex case and
 * maximised in the concave one, plus either an upper or a lower bound in
 * the convex/concave case), so that it is always feasible (but it may be
 * unbounded below/above, represented in terms of linear inequalities in
 * an otherwise "empty" AbstractBlock. The AbstractBlock is solved by any
 * number of different *MILPSolver (you can decide how many solver attach
 * to the Block by changing the LPPar.txt BlockSolverConfig file) and the
 * results are compared.  The Block is then repeatedly randomly modified
 * "in all possible ways", and re-solved several times.
 *
 * \author Antonio Frangioni \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \author Enrico Calandrini \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \copyright &copy; by Antonio Frangioni, Enrico Calandrini
 */
/*--------------------------------------------------------------------------*/
/*-------------------------------- MACROS ----------------------------------*/
/*--------------------------------------------------------------------------*/

#define LOG_LEVEL 0
// 0 = only pass/fail
// 1 = result of each test
// 2 = + solver log
// 3 = + save LP file
// 4 = + print data
//
// note: to always save the LP file with the same name it would be enough to
//       directly set strOutputFile in the configuration file, but the
//       tester rather saves the LP file of each iteration i in a different
//       LPBlock-<i>.lp file, which cannot be done with just the config file

#if( LOG_LEVEL >= 1 )
 #define LOG1( x ) cout << x
 #define CLOG1( y , x ) if( y ) cout << x

 #if( LOG_LEVEL >= 2 )
  #define LOG_ON_COUT 1
 #endif
#else
 #define LOG1( x )
 #define CLOG1( y , x )
#endif

/*--------------------------------------------------------------------------*/

// if nonzero, we avoid that bounds on variable initialized with rhs (or
// equally lhs) infinite become ranged (i.e., both rhs and lhs finite).
// This is beecause some *MILPSolver (e.g. GRBMILPSolver) could have 
// restrictions on the use of ranged constraints and in this way we make 
// sure that no constraint changes from non ranged to ranged one.
#define CONTROL_RANGED 0

/*--------------------------------------------------------------------------*/

// if nonzero, we are considering only variables with finite bound.
// This is because some *MILPSolver (e.g. SCIPMILPSolver) could have 
// some problems with interior point method in the case of unbounded variables.
#define BOUND_FINITE 1

/*--------------------------------------------------------------------------*/

// if nonzero, the Solver attached to the LPBlock is detached and re-attached
// to it at all iterations

#define DETACH_LP 0

/*--------------------------------------------------------------------------*/
// if nonzero, the Block is not solved at every round of changes, but
// only every SKIP_BEAT + 1 rounds. this allows changes to accumulate, and
// therefore puts more pressure on the Modification handling of the Solver
// (in case this tries to do "smart" things rather than dumbly processing
// each one in turn)
//
// note that the number of rounds of changes is them multiplied by
// SKIP_BEAT + 1, so that the input parameter still dictates the number of
// Block solutions

#define SKIP_BEAT 2

/*--------------------------------------------------------------------------*/

#define PANICMSG { cout << endl << "something very bad happened!" << endl; \
		   exit( 1 ); \
                   }

#define PANIC( x ) if( ! ( x ) ) PANICMSG

#define USECOLORS 1
#if( USECOLORS )
 #define RED( x ) "\x1B[31m" #x "\033[0m"
 #define GREEN( x ) "\x1B[32m" #x "\033[0m"
#else
 #define RED( x ) #x
 #define GREEN( x ) #x
#endif

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include <fstream>
#include <sstream>
#include <iomanip>

#include <random>

#include "AbstractBlock.h"

#include "BlockSolverConfig.h"

#include "FRealObjective.h"

#include "FRowConstraint.h"

#if( LOG_LEVEL >= 3 )
 #include "MILPSolver.h"
#endif

#include "LinearFunction.h"

#include "OneVarConstraint.h"

/*--------------------------------------------------------------------------*/
/*-------------------------------- USING -----------------------------------*/
/*--------------------------------------------------------------------------*/

using namespace std;

using namespace SMSpp_di_unipi_it;

/*--------------------------------------------------------------------------*/
/*-------------------------------- TYPES -----------------------------------*/
/*--------------------------------------------------------------------------*/

using Index = Block::Index;
using c_Index = Block::c_Index;

using Range = Block::Range;
using c_Range = Block::c_Range;

using Subset = Block::Subset;
using c_Subset = Block::c_Subset;

using FunctionValue = Function::FunctionValue;
using c_FunctionValue = Function::c_FunctionValue;

using RealVector = std::vector< FunctionValue >;
///< a real n-vector, useful for both the rows of A and b

using c_RealVector = const RealVector;   ///< a const RealVector

using MultiVector = std::vector< RealVector >;

using p_LF = LinearFunction *;

/*--------------------------------------------------------------------------*/
/*------------------------------- CONSTANTS --------------------------------*/
/*--------------------------------------------------------------------------*/

const double scale = 10;
const char *const logF = "log.bn";

const FunctionValue INF = SMSpp_di_unipi_it::Inf< FunctionValue >();

/*--------------------------------------------------------------------------*/
/*------------------------------- GLOBALS ----------------------------------*/
/*--------------------------------------------------------------------------*/

AbstractBlock * LPBlock;   // the problem expressed as an LP

bool convex = true;        // true if the polyhedral function is convex

double bound = 1000;       // a tentative bound to detect unbounded instances

FunctionValue BND;         // the bound in the PolyhedralFunction (if any)

Index nvar = 10;           // number of variables
Index nsvar;              // number of static variables
Index ndvar;              // number of dynamic variables

Index m;                   // number of rows
Index nranged;             // number of ranged constraint in the lpblock

std::mt19937 rg;           // base random generator
std::uniform_real_distribution<> dis( 0.0 , 1.0 );

MultiVector A;
RealVector b;

ColVariable * vLP;                 // pointer to v LP variable

std::vector< ColVariable > * xLP;  // pointer to (static) x LP variables

std::list< ColVariable > * xLPd;  // pointer to (dynamic) x LP variables

std::list< FRowConstraint > * LPbnd;  // FRowConstrait for LPBlock

#if CONTROL_RANGED
  // vector to store information about variable bound: 
  // if bound_ranged[i] == true, then the i-th variable has been initialized
  // with ranged bound.
  std::vector< bool > * bound_ranged;
#endif

/*--------------------------------------------------------------------------*/
/*------------------------------ FUNCTIONS ---------------------------------*/
/*--------------------------------------------------------------------------*/

// convex ==> minimize ==> negative numbers

static double rs( double x ) { return( convex ? -x : x ); }

/*--------------------------------------------------------------------------*/

template< class T >
static void Str2Sthg( const char* const str , T &sthg )
{
 istringstream( str ) >> sthg;
 }

/*--------------------------------------------------------------------------*/

static double rndfctr( void )
{
 // return a random number between 0.5 and 2, with 50% probability of being
 // < 1
 double fctr = dis( rg ) - 0.5;
 return( fctr < 0 ? - fctr : fctr * 4 );
 }

/*--------------------------------------------------------------------------*/

static void GenerateA( Index nr , Index nc )
{
 A.resize( nr );

 for( auto & Ai : A ) {
  Ai.resize( nc );
  for( auto & aij : Ai )
   aij = scale * ( 2 * dis( rg ) - 1 );
  }
 }

/*--------------------------------------------------------------------------*/

static void Generateb( Index nr )
{
 b.resize( nr );

 for( auto & bj : b )
  bj = scale * nvar * ( 2 * dis( rg ) - 1 ) / 4;
 }

/*--------------------------------------------------------------------------*/

static void GenerateAb( Index nr , Index nc )
{
 // rationale: the solution x^* will be more or less the solution of some
 // square sub-system A_B x = b_B. We want x^* to be "well scaled", i.e.,
 // the entries to be ~= 1 (in absolute value). The average of each row A_i
 // is 0, the maximum (and minimum) expected value is something like
 // scale * nvar / 2. So we take each b_j in +- scale * nvar / 4

 GenerateA( nr , nc );
 Generateb( nr );
 }

/*--------------------------------------------------------------------------*/

static void GenerateBND( void )
{
 // rationale: we expect the solution x^* to have entries ~= 1 (in absolute
 // value, and the coefficients of A are <= scale (in absolute value), so
 // the LHS should be at most around - scale * nvar; the RHS can add it
 // a further - scale * nvar / 4, so we expect - (5/4) * scale * nvar to
 // be a "natural" LB. We therefore set the LB to a mean of 1/2 of that
 // (tight) 33% of the time, a mean of 2 times that (loose) 33% of the time,
 // and -INF the rest

 if( dis( rg ) <= 0.333 ) {   // "tight" bound
  BND = rs( dis( rg ) * 5 * scale * nvar / 4 );
  return;
  }

 if( dis( rg ) <= 0.333 ) {  // "loose" bound
  BND = rs( dis( rg ) * 5 * scale * nvar );
  return;
  }

 BND = INF;
 }

/*--------------------------------------------------------------------------*/

static Subset GenerateRand( Index m , Index k )
{
 // generate a sorted random k-vector of unique integers in 0 ... m - 1

 Subset rnd( m );
 std::iota( rnd.begin() , rnd.end() , 0 );
 std::shuffle( rnd.begin() , rnd.end() , rg );
 rnd.resize( k );
 sort( rnd.begin() , rnd.end() );

 return( std::move( rnd ) );
 }

 /*--------------------------------------------------------------------------*/

static void ConstructLPConstraint( Index i , FRowConstraint & ci ,
				   bool setblock = true )
{
 // construct constraint ci out of A[ i ] and b[ i ]:
 //
 // in the convex case, the constraint is
 //
 //          b[ i ] <= vLP1 - \sum_j Ai[ j ] * xLP[ j ] <= INF
 //
 // in the concave case, the constraint is
 //
 //          -INF <= vLP1 - \sum_j Ai[ j ] * xLP[ j ] <= b[ i ]
 //
 // note: constraints are constructed dense (elements == 0, which are
 //       anyway quite unlikely, are ignored) to make things simpler
 //
 // note: variable x[ i ] is given index i + 1, variable v has index 0

 if( convex ) {
  ci.set_lhs( b[ i ] );
  ci.set_rhs( INF );
  }
 else {
  ci.set_lhs( -INF );
  ci.set_rhs( b[ i ] );
  }
 LinearFunction::v_coeff_pair vars( nvar + 1 );

 Index j = 0;

 // first, v
 vars[ j ] = std::make_pair( vLP , 1 );

 // then, static x
 for( ; j < nsvar ; ++j )
  vars[ j + 1 ] = std::make_pair( &((*xLP)[ j ] ) , - A[ i ][ j ] );

 // finally, dynamic x
 auto xLPdit = xLPd->begin();
 for( ; j < nvar ; ++j , ++xLPdit )
   vars[ j + 1 ] = std::make_pair( &(*xLPdit) , - A[ i ][ j ] );

 ci.set_function( new LinearFunction( std::move( vars ) ) );
 if( setblock )
  ci.set_Block( LPBlock );
 }

/*--------------------------------------------------------------------------*/

static void ChangeLPConstraint( Index i , FRowConstraint & ci , ModParam iAM )
{
 // change the constant == LHS or RHS of the constraint (depending on convex)
 if( convex )
  ci.set_lhs( b[ i ] , iAM );
 else
  ci.set_rhs( b[ i ] , iAM );

 // now change the coefficients, except that of v that is always 1
 LinearFunction::Vec_FunctionValue coeffs( nvar );

 for( Index j = 0 ; j < nvar ; ++j )
  coeffs[ j ] = - A[ i ][ j ];

 auto f = static_cast< p_LF >( ci.get_function() );
 f->modify_coefficients( std::move( coeffs ) , Range( 1 , nvar + 1 ) , iAM );
 }

/*--------------------------------------------------------------------------*/

static std::pair< double , double > Generate_lhs_rhs( double const p ) {
  double lhs , rhs;
  #if BOUNS_FINITE == 1
    auto p2 = dis( rg );
    lhs = p2 < 0.5 ? p2 : 0;
    rhs = p2 < 0.5 ? 1 : p2;
    return { lhs , rhs };
  #endif
  if( p < 0.333 ) { // lhs finite, rhs INF
      lhs = dis( rg );
      rhs = INF;
    }
    else if( p >= 0.333 && p < 0.666) { // both lhs and rhs finite
      auto p2 = dis( rg );
      lhs = p2 < 0.5 ? p2 : 0;
      rhs = p2 < 0.5 ? 1 : p2;
    }
    else{ // lhs -INF, rhs finite
      lhs = -INF;
      rhs = dis( rg );
    }
  return { lhs , rhs };
}

/*--------------------------------------------------------------------------*/

static inline void SetFRow( ColVariable & LPxi )
{
  LPbnd->resize( LPbnd->size() + 1 );
  LinearFunction::v_coeff_pair vars_LP( 1 );
  vars_LP[ 0 ] = std::make_pair( & LPxi , 1 );
  LPbnd->back().set_function( new LinearFunction( std::move( vars_LP ) ) );
  auto p = dis( rg );
  std::pair< double , double > bounds = Generate_lhs_rhs( p );
  LPbnd->back().set_lhs( bounds.first , eNoMod );
  LPbnd->back().set_rhs( bounds.second , eNoMod );
  ++nranged;
  #if CONTROL_RANGED
    if( bounds.first == -INF || bounds.second == INF )
      bound_ranged->push_back( false );
    else
      bound_ranged->push_back( true );
  #endif
 }

/*--------------------------------------------------------------------------*/

static void RemoveFRow( AbstractBlock & AB , Range rng )
{
 // the dynamic variable from the "xd" group in the Range are removed: if
 // anything is "active" in those is a FRowConstraint from the "xbnd" group
 // that has to be removed as well

 auto xd = AB.get_dynamic_variable< ColVariable >( "xd" );
 auto itxd = std::next( xd->begin() , rng.first );
 #if CONTROL_RANGED
  auto itcontrol = std::next( bound_ranged->begin() , rng.first + nsvar );
 #endif
 for( Index i = rng.first ; i < rng.second ; ++i , ++itxd ) {
  if( ! itxd->get_num_active() )
   continue;
  --nranged;
  int numbox = itxd->get_num_active();
  for( int j = 0 ; j < numbox ; ++j ) {
   std::vector< typename std::list< FRowConstraint >::iterator > rmvd;
   auto & frow = *(AB.get_dynamic_constraint< FRowConstraint >( "xbnd" ));
   auto rc = dynamic_cast< FRowConstraint * >( itxd->get_active( 0 ) );
  if( ! rc ) {
   cout << "Unexpected stuff active in to-be-deleted Variable" << endl;
   exit( 1 );
   }
   auto to_remove = std::find_if( frow.begin() , frow.end() ,
			  [ rc ]( FRowConstraint & x ) {
			   return( & x == rc );
			   } );
   if( to_remove == frow.end() ) {
    cout << "FRowConstraint not found" << endl;
    exit( 1 );
    }
   rmvd.push_back( to_remove );
   AB.remove_dynamic_constraints( frow , rmvd ); 
   }
  #if CONTROL_RANGED
    bound_ranged->erase( itcontrol );
  #endif
  }
 }

/*--------------------------------------------------------------------------*/

static void RemoveFRow( AbstractBlock & AB , const Subset & sbst )
{
 // the dynamic variable from the "xd" group in the (ordered) Subset are
 // removed: if anything is "active" in those is a FRowConstraint from the
 // "xbnd" group that has to be removed as well

 auto xd = AB.get_dynamic_variable< ColVariable >( "xd" );
 Index prev = 0;
 auto itxd = xd->begin();
 int n_removed = 0;
 for( auto ind : sbst ) {
  itxd = std::next( itxd , ind - prev );
  #if CONTROL_RANGED
    auto itcontrol = std::next( bound_ranged->begin() , nsvar + ind - n_removed );
    bound_ranged->erase( itcontrol );
  #endif
  ++n_removed;
  --nranged;
  prev = ind;
  if( ! itxd->get_num_active() )
   continue;
  int numbox = itxd->get_num_active();
  for( int j = 0 ; j < numbox ; ++j ) {
   std::vector< typename std::list< FRowConstraint >::iterator > rmvd;
   auto & frow = *(AB.get_dynamic_constraint< FRowConstraint >( "xbnd" ));
   auto rc = dynamic_cast< FRowConstraint * >( itxd->get_active( 0 ) );
   if( ! rc ) {
    cout << "Unexpected stuff active in to-be-deleted Variable" << endl;
    exit( 1 );
    }
   auto to_remove = std::find_if( frow.begin() , frow.end() ,
			  [ rc ]( FRowConstraint & x ) {
			   return( & x == rc );
			   } );
   if( to_remove == frow.end() ) {
    cout << "FRowConstraint not found" << endl;
    exit( 1 );
    }
   rmvd.push_back( to_remove );
   AB.remove_dynamic_constraints( frow , rmvd ); 
   }
  }
 }

/*--------------------------------------------------------------------------*/

static void ChangeFRow( AbstractBlock & AB , const Subset & sbst , 
                const bool control_rep )
{
 auto frow = AB.get_dynamic_constraint< FRowConstraint >( "xbnd" );
 Index prev = 0;
 auto frowit = frow->begin();

 #if CONTROL_RANGED
  auto itcontrol = bound_ranged->begin();
 #endif

 for( auto ind : sbst ) {
  #if CONTROL_RANGED
    itcontrol = std::next( itcontrol , ind - prev );
    frowit = std::next( frowit , ind - prev );
    prev = ind;
    double lhs, rhs;
    if( *itcontrol == false ) {
      // we have to check that the bound doesn't become ranged
      auto p = dis( rg );
      lhs = p < 0.5 ? p : -INF;
      rhs = p < 0.5 ? INF : p;
      }
    else{ // any situation could be reproduced
      auto p = dis( rg );
      std::pair< double , double > bounds = Generate_lhs_rhs( p );
      lhs = bounds.first;
      rhs = bounds.second;
      if( control_rep == true ){
        if( lhs == -INF )
         lhs = 0;
        else if( rhs == INF ){
         rhs = 1;
        }
       }
      }
    (*frowit).set_lhs( lhs );
    (*frowit).set_rhs( rhs );
  #else
    frowit = std::next( frowit , ind - prev );
    prev = ind;
    auto p = dis( rg );
    std::pair< double , double > bounds = Generate_lhs_rhs( p );
    (*frowit).set_lhs( bounds.first );
    (*frowit).set_rhs( bounds.second );
  #endif
 }
}

/*--------------------------------------------------------------------------*/

static void ChangeFRow( AbstractBlock & AB , Range rng ,
                  const bool control_rep )
{
 auto frow = AB.get_dynamic_constraint< FRowConstraint >( "xbnd" );
 auto frowit = std::next( frow->begin() , rng.first );
 #if CONTROL_RANGED
  auto itcontrol = std::next( bound_ranged->begin() , rng.first );
 #endif
 for( Index i = rng.first ; i < rng.second ; ++i , ++frowit ) {
  #if CONTROL_RANGED
    double lhs;
    double rhs;
    if( *itcontrol == false ) {
      // we have to check that the bound doesn't become ranged
      auto p = dis( rg );
      lhs = p < 0.5 ? p : -INF;
      rhs = p < 0.5 ? INF : p;
      }
    else{ // any situation could be reproduced
      auto p = dis( rg );
      std::pair< double , double > bounds = Generate_lhs_rhs( p );
      lhs = bounds.first;
      rhs = bounds.second;
      if( control_rep == true ){
        if( lhs == -INF )
         lhs = 0;
        else if( rhs == INF ){
         rhs = 1;    
        }
       }
      }
    (*frowit).set_lhs( lhs );
    (*frowit).set_rhs( rhs );
    ++itcontrol;
  #else
    auto p = dis( rg );
    std::pair< double , double > bounds = Generate_lhs_rhs( p );
    (*frowit).set_lhs( bounds.first );
    (*frowit).set_rhs( bounds.second );
  #endif
 }
}

/*--------------------------------------------------------------------------*/

static void printAb( const MultiVector & tA , const RealVector & tb ,
		     double bound )
{
 PANIC( ( tA.size() == tb.size() ) || ( tA.size() + 1 == tb.size() ) );
 PANIC( tA.size() == m );
 for( auto & tai : tA )
  PANIC( tai.size() == nvar );

 cout << "n = " << nvar << ", m = " << m;
 if( std::abs( bound ) == INF )
  cout << " (no bound)" << endl;
 else
  cout << ", bound = " << bound << endl;

 for( Index i = 0 ; i < m ; ++i ) {
  cout << "A[ " << i << " ] = [ ";
  for( Index j = 0 ; j < nvar ; ++j )
   cout << tA[ i ][ j ] << " ";
   cout << "], b[ " << i << " ] = " << tb[ i ] << endl;
  }
 }

/*--------------------------------------------------------------------------*/

// Some functions used to check results of solvers in loop

bool CompareSolution( double const d1 , double const d2 )
{
 return( abs( d1 - d2 ) >= 2e-7 *
	 max( double( 1 ) , abs( max( d1 , d2 ) ) ) );
 }

bool allEqual( std::vector< double > const & v )
{
 return( std::adjacent_find( v.begin() , v.end() , CompareSolution )
	 == v.end() );
 }

bool allTrue( std::vector< bool > const & v )
{
 return( std::all_of( v.begin() , v.end() ,
		      []( bool i ) { return( i == true ); } ) );
 }

bool allInfeasible( std::vector< int > const & v )
{
 return( std::all_of( v.begin() , v.end() ,
		      []( int i ) { return( i == Solver::kInfeasible ); } ) );
 }

bool allUnbounded( std::vector< int > const  & v )
{
 return( std::all_of( v.begin() , v.end() ,
		      []( int i ) { return( i == Solver::kUnbounded ); } ) );
 }

/*--------------------------------------------------------------------------*/

static bool SolveAll( void ) 
{
 try {
  auto slvr_list = LPBlock->get_registered_solvers();
  auto itslvr = slvr_list.begin();
  int num_slvr = slvr_list.size();
  std::vector< int > rtrnLP( num_slvr );
  std::vector< bool > hsLP( num_slvr );
  std::vector< double > foLP( num_slvr );

  for( int j = 0 ; j < num_slvr ; ++j , ++itslvr ) {
    // solve the LPBlock with the j^th solvers- - - - - - - - - - - - - - - - - - 

    Solver * slvrLP = *itslvr;
    #if DETACH_LP
     LPBlock->unregister_Solver( slvrLP );
     LPBlock->register_Solver( slvrLP , true );  // push it to the front
    #endif
    rtrnLP[ j ] = slvrLP->compute( false );
    hsLP[ j ] = ( ( rtrnLP[j] >= Solver::kOK ) && ( rtrnLP[ j ] < Solver::kError ) )
                || ( rtrnLP[ j ] == Solver::kLowPrecision );
    foLP[ j ] = hsLP[ j ] ? ( convex ? slvrLP->get_ub() : slvrLP->get_lb() )
                      : ( convex ? INF : -INF );
    }

  if( allTrue( hsLP ) && ( allEqual( foLP ) ) ) {
   LOG1( "OK(f)" << endl );
   return( true );
   }

  if( allInfeasible( rtrnLP ) ) {
    LOG1( "OK(?e?)" << endl );
    return( true );
    }

  if( allUnbounded( rtrnLP ) ) {
   LOG1( "OK(u)" << endl );
   return( true );
   }

  #if( LOG_LEVEL >= 1 )
   for( int j = 0 ; j < num_slvr ; ++j ) {
    cout << "Solver" << j <<  " = ";
    if( hsLP[ j ] ){
     cout << foLP[ j ] << " -- ";
     }
    else
     if( rtrnLP[ j ] == Solver::kInfeasible )
      cout << " Unfeas(?) -- ";
    else
     if( rtrnLP[ j ] == Solver::kUnbounded )
      cout << " Unbounded -- ";
     else
      cout << " Error! -- ";
    }
   cout << endl;
  #endif

  return( false );
  }
 catch( exception &e ) {
  cerr << e.what() << endl;
  exit( 1 );
  }
 catch(...) {
  cerr << "Error: unknown exception thrown" << endl;
  exit( 1 );
  }
 }

/*--------------------------------------------------------------------------*/

int main( int argc , char **argv )
{
 // reading command line parameters - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 assert( SKIP_BEAT >= 0 );

 long int seed = 0;
 Index wchg = 127;
 double dens = 4;  
 double p_change = 0.5;
 Index n_change = 10;
 Index n_repeat = 40;

 switch( argc ) {
  case( 8 ): Str2Sthg( argv[ 7 ] , p_change );
  case( 7 ): Str2Sthg( argv[ 6 ] , n_change );
  case( 6 ): Str2Sthg( argv[ 5 ] , n_repeat );
  case( 5 ): Str2Sthg( argv[ 4 ] , dens );
  case( 4 ): Str2Sthg( argv[ 3 ] , nvar );
  case( 3 ): Str2Sthg( argv[ 2 ] , wchg );
  case( 2 ): Str2Sthg( argv[ 1 ] , seed );
             break;
  default: cerr << "Usage: " << argv[ 0 ] <<
	   " seed [wchg nvar dens #rounds #chng %chng]"
 		<< endl <<
           "       wchg: what to change, coded bit-wise [255]"
		<< endl <<
           "             0 = add rows, 1 = delete rows "
		<< endl <<
           "             2 = modify rows, 3 = modify constants"
		<< endl <<
           "             4 = change global lower/upper bound"
		<< endl <<
           "             5 = add variables, 6 = delete variables"
    << endl <<
           "             7 = change variables bounds"
	        << endl <<
           "       nvar: number of variables [10]"
	        << endl <<
           "       dens: rows / variables [4]"
	        << endl <<
           "       #rounds: how many iterations [40]"
	        << endl <<
           "       #chng: number changes [10]"
	        << endl <<
           "       %chng: probability of changing [0.5]"
	        << endl;
	   return( 1 );
  }

 if( nvar < 1 ) {
  cout << "error: nvar too small";
  exit( 1 );
  }

 nsvar = nvar / 2;      // half of the variables are dynamic
 ndvar = nvar - nsvar;  // the other half are static
 nranged = 0;

 m = nvar * dens;
 if( m < 1 ) {
  cout << "error: dens too small";
  exit( 1 );
  }

 rg.seed( seed );  // seed the pseudo-random number generator

 // constructing the data of the problem- - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // choosing whether convex or concave: toss a(n unbiased, two-sided) coin
 convex = ( dis( rg ) < 0.5 );

 // construct the matrix m x nvar matrix A and the m-vector b

 GenerateAb( m , nvar );
 GenerateBND();

 cout.setf( ios::scientific, ios::floatfield );
 cout << setprecision( 10 );

 #if( LOG_LEVEL >= 4 )
  printAb( A , b , BND );
 #endif

 // construction and loading of the objects - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 // construct the LP- - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 {
  // ensure all original pointers go out of scope immediately after that
  // the construction has finished

  LPBlock = new AbstractBlock();

  // construct the Variable
  xLP = new std::vector< ColVariable >( nsvar );
  xLPd = new std::list< ColVariable >( ndvar );

  vLP = new ColVariable;
  vLP->set_Block( LPBlock );

  // construct the m dynamic Constraint
  auto ALP = new std::list< FRowConstraint >( m );
  auto ALPit = ALP->begin();
  for( Index i = 0 ; i < m ; )
   ConstructLPConstraint( i++ , *(ALPit++) );

  // construct the static lower bound Constraint
  auto LBc = new BoxConstraint( LPBlock , vLP , -INF , INF );
  if( BND != INF ) {
   if( convex )
    LBc->set_lhs( -BND );
   else
    LBc->set_rhs( BND );
   }

  // construct the Objective
  auto objLP = new FRealObjective();
  objLP->set_function( new LinearFunction( { std::make_pair( vLP , 1 ) } ) );
  objLP->set_sense( convex ? Objective::eMin : Objective::eMax , eNoMod );
  
  // now set the Variable, Constraint and Objective in the AbstractBlock
  LPBlock->add_static_variable( *vLP , "v" );
  LPBlock->add_static_variable( *xLP , "x" );
  LPBlock->add_dynamic_variable( *xLPd , "xd" );
  LPBlock->add_dynamic_constraint( *ALP , "cuts" );
  LPBlock->add_static_constraint( *LBc , "vbnd" );
  LPBlock->set_objective( objLP );
  }

 // define bound constraints- - - - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 LPbnd = new std::list< FRowConstraint >;
 #if CONTROL_RANGED
  bound_ranged = new std::vector< bool >;
 #endif
 auto & LPx = *(LPBlock->get_static_variable_v< ColVariable >( "x" ));
 for( Index i = 0 ; i < nsvar ; ++i )
  SetFRow( LPx[ i ] );

 auto LPxd = LPBlock->get_dynamic_variable< ColVariable >( "xd" )->begin();
 for( Index i = 0 ; i < ndvar ; ++i )
  SetFRow( *(LPxd++) );

 // note: the list may be empty, but it is intentionally added anyway
 LPBlock->add_dynamic_constraint( *LPbnd , "xbnd" );
 
 // attach (at least) two Solver to the LPBlock - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // do it by using a single a BlockSolverConfig, read from file
 
 auto lpbsc = dynamic_cast< BlockSolverConfig * >(
		     Configuration::deserialize( "LPPar.txt" ) );
 if( ! lpbsc ) {
  cerr << "Error: configuration file not a BlockSolverConfig" << endl;
  exit( 1 );    
  }

 lpbsc->apply( LPBlock );
 lpbsc->clear();  // keep the clear()-ed BlockSolverConfig for final cleanup

 // check Solvers - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( LPBlock->get_registered_solvers().empty() ) {
  cerr << "Error: BlockSolverConfig did not register any Solver" << endl;
  exit( 1 );    
  }

 // open log-file - - - - - - - - - - -  - - - - - - - - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 #if( LOG_LEVEL >= 3 )
   auto slvr_list = LPBlock->get_registered_solvers();
   auto itslvr = slvr_list.begin();
   for( int j = 0 ; j < slvr_list.size() ; ++j )
    (*(itslvr++))->set_par( MILPSolver::strOutputFile , 
			    "Solver" + std::to_string( j ) + ".lp" );
 #endif

 // first solver call - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 LOG1( "First call: " );

 bool AllPassed = SolveAll();
 
 // main loop - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // now, for n_repeat times in p_change% of the cases
 // - up to n_change rows are added
 // - up to n_change rows are deleted
 // - up to n_change rows are modified
 // - up to n_change rows are modified
 // - the bound is modified
 // - min( 1 , rnd( nsvar / 4 ) ) variables are added
 // - up to ndvar variables are removed
 //
 // then the two problems are re-solved

 for( Index rep = 0 ; rep < n_repeat * ( SKIP_BEAT + 1 ) ; ) {
  if( ! AllPassed )
   break;

  LOG1( rep << ": ");

  // add rows - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  if( ( wchg & 1 ) && ( dis( rg ) <= p_change ) )
   if( Index tochange = Index( dis( rg ) * n_change ) ) {
    LOG1( "added " << tochange << " rows - " );

    GenerateAb( tochange , nvar );

    // add them to the first LP
    vLP = LPBlock->get_static_variable< ColVariable >( "v" );
    xLP = LPBlock->get_static_variable_v< ColVariable >( "x" );
    xLPd = LPBlock->get_dynamic_variable< ColVariable >( "xd" );

    std::list< FRowConstraint > nc( tochange );
    auto ncit = nc.begin();
    for( Index i = 0 ; i < tochange ; )
     ConstructLPConstraint( i++ , *(ncit++) );
    auto cnst = LPBlock->get_dynamic_constraint< FRowConstraint >( "cuts" );
    LPBlock->add_dynamic_constraints( *cnst , nc );

    // update m
    m += tochange;

    // sanity checks
    PANIC( m == cnst->size() );
    }

  // delete rows- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  if( ( wchg & 2 ) && ( dis( rg ) <= p_change ) )
   if( Index tochange = min( m - 1 , Index( dis( rg ) * n_change ) ) ) {
    LOG1( "deleted " << tochange << " rows" );

    auto cnst = LPBlock->get_dynamic_constraint< FRowConstraint >( "cuts" );
    
    if( dis( rg ) <= 0.5 ) {  // in 50% of the cases do a ranged change
     LOG1( "(r) - " );

     Index strt = dis( rg ) * ( m - tochange );
     Index stp = strt + tochange;

     // remove them from the first LP
     LPBlock->remove_dynamic_constraints( *cnst , Range( strt , stp ) );
     }
    else {  // in the other 50% of the cases, do a sparse change
     LOG1( "(s) - " );
     Subset nms( GenerateRand( m , tochange ) );

     // remove them from the LP
     if( tochange == 1 )
      LPBlock->remove_dynamic_constraint( *cnst , std::next( cnst->begin() ,
							     nms[ 0 ] ) );
     else
      LPBlock->remove_dynamic_constraints( *cnst , Subset( nms ) , true );
     }

    // update m
    m -= tochange;

    // sanity checks
    PANIC( m == cnst->size() );
    }

  // modify rows- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  if( ( wchg & 4 ) && ( dis( rg ) <= p_change ) )
   if( Index tochange = std::min( m , Index( dis( rg ) * n_change ) ) ) {
    LOG1( "modified " << tochange << " rows" );

    GenerateAb( tochange , nvar );

    vLP = LPBlock->get_static_variable< ColVariable >( "v" );
    xLP = LPBlock->get_static_variable_v< ColVariable >( "x" );
    xLPd = LPBlock->get_dynamic_variable< ColVariable >( "xd" );
    auto cnst = LPBlock->get_dynamic_constraint< FRowConstraint >( "cuts" );

    if( dis( rg ) <= 0.5 ) {  // in 50% of the cases do a ranged change
     LOG1( "(r) - " );

     Index strt = dis( rg ) * ( m - tochange );
     Index stp = strt + tochange;

     // send all the Modification to the same channel
     Observer::ChnlName chnl = LPBlock->open_channel();
     const auto iAM = Observer::make_par( eModBlck , chnl );

     // modify them in the LP
     auto cit = std::next( cnst->begin() , strt );
     for( Index i = 0 ; i < tochange ; ++i )
      ChangeLPConstraint( i , *(cit++) , iAM );

     LPBlock->close_channel( chnl );  // close the channel
     }
    else {  // in the other 50% of the cases, do a sparse change
     LOG1( "(s) - " );
     Subset nms( GenerateRand( m , tochange ) );

     // send all the Modification to the same channel
     Observer::ChnlName chnl = LPBlock->open_channel();
     const auto iAM = Observer::make_par( eModBlck , chnl );

     // modify them in the LP
     Index prev = 0;
     auto cit = cnst->begin();
     for( Index i = 0 ; i < tochange ; ++i ) {
      cit = std::next( cit , nms[ i ] - prev );
      prev = nms[ i ];
      ChangeLPConstraint( i , *cit , iAM );
      }

     LPBlock->close_channel( chnl );  // close the channel
    }
   }

  // modify constants - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  if( ( wchg & 8 ) && ( dis( rg ) <= p_change ) )
   if( Index tochange = std::min( m , Index( dis( rg ) * n_change ) ) ) {
    LOG1( "modified " << tochange << " constants" );

    Generateb( tochange );
     
    auto cnst = LPBlock->get_dynamic_constraint< FRowConstraint >( "cuts" );

    if( dis( rg ) <= 0.5 ) {  // in 50% of the cases do a ranged change
     LOG1( "(r) - " );

     Index strt = dis( rg ) * ( m - tochange );
     Index stp = strt + tochange;

     // change them in the LP
     auto cit = std::next( cnst->begin() , strt );
     if( convex )
      for( Index i = 0 ; i < tochange ; )
       (*(cit++)).set_lhs( b[ i++ ] );
     else
      for( Index i = 0 ; i < tochange ; )
       (*(cit++)).set_rhs( b[ i++ ] );
     }
    else {  // in the other 50% of the cases, do a sparse change
     LOG1( "(s) - " );
     Subset nms( GenerateRand( m , tochange ) );

     // change them in the LP
     Index prev = 0;
     auto cit = cnst->begin();
     if( convex )
      for( Index i = 0 ; i < tochange ; ) {
       cit = std::next( cit , nms[ i ] - prev );
       prev = nms[ i ];
       (*cit).set_lhs( b[ i++ ] );
       }
     else
      for( Index i = 0 ; i < tochange ; ) {
       cit = std::next( cit , nms[ i ] - prev );
       prev = nms[ i ];
       (*cit).set_rhs( b[ i++ ] );
       }
     }
    }

  // modify bound - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  if( ( wchg & 16 ) && ( dis( rg ) <= p_change ) ) {
   LOG1( "modified bound - " );

   GenerateBND();

   // change it in the LP
   auto cnst = LPBlock->get_static_constraint< BoxConstraint >( "vbnd" );
   if( convex )
    cnst->set_lhs( -BND );
   else
    cnst->set_rhs( BND );
   }

 // add variables- - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 Index added_ndvar = 0; // number of dynamic variables added in current loop

  if( ( wchg & 32 ) && ( dis( rg ) <= p_change ) ) {
   Index tochange = std::max( Index( 1 ) , Index( dis( rg ) * nsvar / 4 ) );
   added_ndvar = tochange;
   LOG1( "added " << tochange << " variables - " );

   GenerateA( m , tochange );

   // add them in the LP
   std::list< ColVariable > nxLPd( tochange );
   std::vector< ColVariable * > nxp( tochange );
   auto nxlpit = nxLPd.begin();
   for( Index i = 0 ; i < tochange ; )
    nxp[ i++ ] = &(*(nxlpit++));

   LPBlock->add_dynamic_variables(
	   *(LPBlock->get_dynamic_variable< ColVariable >( "xd" )) , nxLPd );

   auto cnst_it =
        LPBlock->get_dynamic_constraint< FRowConstraint >( "cuts" )->begin();
   if( tochange == 1 )
    for( Index i = 0 ; i < m ; ++i ) {
     auto fi = static_cast< p_LF >( (cnst_it++)->get_function() );
     fi->add_variable( nxp[ 0 ] , - A[ i ][ 0 ] );
     }
   else
    for( Index i = 0 ; i < m ; ++i ) {
     auto fi = static_cast< p_LF >( (cnst_it++)->get_function() );
     LinearFunction::v_coeff_pair ncp( tochange );
     for( Index j = 0 ; j < ncp.size() ; ++j ) {
      ncp[ j ].first = nxp[ j ];
      ncp[ j ].second = - A[ i ][ j ];
      }
     fi->add_variables( std::move( ncp) );
     }

   // generate bound constraints
   auto LPxd = LPBlock->get_dynamic_variable< ColVariable >( 0 );
   auto LPxd_it = LPxd->begin();
   LPxd_it = std::next( LPxd_it , ndvar );
   
   LPbnd = new std::list< FRowConstraint >;

   for( ; LPxd_it != LPxd->end() ; )
    SetFRow( *(LPxd_it++) );

   if( ! LPbnd->empty() )
    LPBlock->add_dynamic_constraints(
	   *(LPBlock->get_dynamic_constraint< FRowConstraint >( "xbnd" )) , *LPbnd );

   // update nvar and ndvar
   nvar += tochange;
   ndvar += tochange;

   // sanity checks
   PANIC( ndvar ==
	      LPBlock->get_dynamic_variable< ColVariable >( 0 )->size() );
   for( auto & ci :
 	    *(LPBlock->get_dynamic_constraint< FRowConstraint >( "cuts" )) )
    PANIC( nvar + 1 == ci.get_num_active_var() );
   }

  // remove variables - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  if( ( wchg & 64 ) && ( dis( rg ) <= p_change ) )
   if( Index tochange = Index( dis( rg ) * ndvar ) ) {
    LOG1( "removed " << tochange << " variables" );

    if( dis( rg ) <= 0.5 ) {  // in 50% of the cases do a ranged removal
     LOG1( "(r) - " );

     Index strt = dis( rg ) * ( ndvar - tochange );
     Index stp = strt + tochange;

     // remove them from the LP
     auto xLPd = LPBlock->get_dynamic_variable< ColVariable >( 0 );
     auto cnst_it =
             LPBlock->get_dynamic_constraint< FRowConstraint >( 0 )->begin();
     if( tochange == 1 )
      for( Index i = 0 ; i < m ; ++i ) {
       auto fi = static_cast< p_LF >( (cnst_it++)->get_function() );
       fi->remove_variable( strt + nsvar + 1 );
       }
     else
      for( Index i = 0 ; i < m ; ++i ) {
       auto fi = static_cast< p_LF >( (cnst_it++)->get_function() );
       fi->remove_variables( Range( strt + nsvar + 1 , stp + nsvar + 1 ) );
       }

     // the variables can now only be active in the associated frow
     // constraint, if any: exploit this to identify the frow constraint
     // and remove it
     RemoveFRow( *LPBlock , Range( strt , stp ) );
     
     LPBlock->remove_dynamic_variables( *xLPd , Range( strt , stp ) );
     }
    else {  // in the other 50% of the cases, do a sparse change
     LOG1( "(s) - " );
     Subset nms( GenerateRand( ndvar , tochange ) );

     // remove them from the LP
     auto xLPd = LPBlock->get_dynamic_variable< ColVariable >( 0 );
     auto cnst_it =
             LPBlock->get_dynamic_constraint< FRowConstraint >( 0 )->begin();
     if( tochange == 1 ) {
      for( Index i = 0 ; i < m ; ++i ) {
       auto fi = static_cast< p_LF >( (cnst_it++)->get_function() );
       fi->remove_variable( nms[ 0 ] + nsvar + 1 );
      }
      
      // the variables can now only be active in the associated frow
      // constraint, if any: exploit this to identify the frow constraint
      // and remove it
      RemoveFRow( *LPBlock , Range( nms[ 0 ] , nms[ 0 ] + 1 ) );

      auto vp = std::next( xLPd->begin() , nms[ 0 ] );
      LPBlock->remove_dynamic_variable( *xLPd , vp );
      }
     else {
      for( Index i = 0 ; i < m ; ++i ) {
       auto fi = static_cast< p_LF >( (cnst_it++)->get_function() );

       Subset nms1( nms );
       for( auto & n1i : nms1 )
	      n1i = n1i + nsvar + 1;

       fi->remove_variables( std::move( nms1 ) , true );
       }

      // the variables can now only be active in the associated frow
      // constraint, if any: exploit this to identify the frow constraint
      // and remove it
      RemoveFRow( *LPBlock , nms );
      
      LPBlock->remove_dynamic_variables( *xLPd , Subset( nms ) );
     }
    }

    // update ndvar and nvar
    ndvar -= tochange;
    nvar -= tochange;

    // sanity checks
    PANIC( ndvar ==
	         LPBlock->get_dynamic_variable< ColVariable >( 0 )->size() );
    for( auto & ci :
	          *(LPBlock->get_dynamic_constraint< FRowConstraint >( 0 )) )
     PANIC( nvar + 1 == ci.get_num_active_var() );
    }

  // change ranged constraint rhs and lhs - - - - - - - - - - - - - - - - - - - -

  // in order to avoid strange behaviour, we don't want that the bounds on variables
  // added in the current loop changes
  Index n_oldranged = nranged - added_ndvar;
  if( ( wchg & 128 ) && ( dis( rg ) <= p_change ) )
   if( Index tochange = Index( dis( rg ) * n_oldranged ) ) {
    LOG1( "changed " << tochange << " rng limit" );

    if( dis( rg ) <= 0.5 ) {  // in 50% of the cases do a ranged removal
     LOG1( "(r) - " );

     Index strt = dis( rg ) * ( n_oldranged - tochange );
     Index stp = strt + tochange;

     // the variables can now only be active in the associated frow
     // constraint, if any: exploit this to identify the frow constraint
     // and remove it
     ChangeFRow( *LPBlock , Range( strt , stp ) , ( rep % ( SKIP_BEAT + 1 ) != 0 ) );
     }
    else {  // in the other 50% of the cases, do a sparse change
     LOG1( "(s) - " );
     Subset nms( GenerateRand( n_oldranged , tochange ) );

     ChangeFRow( *LPBlock , nms , ( rep % ( SKIP_BEAT + 1 ) != 0 )  );
     }
   }

  // if verbose, print out stuff- - - - - - - - - - - - - - - - - - - - - - -

  #if( LOG_LEVEL >= 3 )
   auto slvr_list = LPBlock->get_registered_solvers();
   auto itslvr = slvr_list.begin();
   for( int j = 0 ; j < slvr_list.size() ; ++j ) {
    (*itslvr)->set_par( MILPSolver::strOutputFile , 
      "Solver" + std::to_string( j ) + "-LPBlock-" + 
        std::to_string( rep ) + ".lp" );
    ++itslvr;
    }
  #endif

  // finally, re-solve the problems- - - - - - - - - - - - - - - - - - - - -
  // ... every SKIP_BEAT + 1 rounds

  if( ! ( ++rep % ( SKIP_BEAT + 1 ) ) )
   AllPassed &= SolveAll();
  #if( LOG_LEVEL >= 1 )
  else
   cout << endl;
  #endif

  }  // end( main loop )- - - - - - - - - - - - - - - - - - - - - - - - - - -
     // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( AllPassed )
  cout << GREEN( All tests passed!! ) << endl;
 else
  cout << RED( Shit happened!! ) << endl;
 
 // destroy the Block - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 lpbsc->apply( LPBlock );

 delete lpbsc;

 // delete the Blocks
 delete( LPBlock );

 // terminate - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 return( AllPassed ? 0 : 1 );

 }  // end( main )

/*--------------------------------------------------------------------------*/
/*------------------------ End File test.cpp -------------------------------*/
/*--------------------------------------------------------------------------*/
