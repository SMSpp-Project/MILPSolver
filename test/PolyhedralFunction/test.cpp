/*--------------------------------------------------------------------------*/
/*-------------------------- File test.cpp -----------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Main for testing PolyhedralFunction
 *
 * A "random" PolyhedralFunction is constructed and represented in terms of 
 * linear inequalities for two otherwise "empty" Block. The two Block are 
 * solved by two different *MILPSolver and the results are compared. 
 * The two Block are then repeatedly randomly modified "in the same way", 
 * and re-solved several times.
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
// if HAVE_CONSTRAINTS == 1, then about 50% of the variables will have a
// non-negativity constraint implemented via ColVariable::is_positive()
// if HAVE_CONSTRAINTS == 2, then about 50% of the variables will have
// bound constraints; of these, 33% will only have 0 lower bound, 33% will
// only have random upper bound, and the rest will have both. of the
// remaining 50% of the variables, another 50%  will have a
// non-negativity constraint implemented via ColVariable::is_positive()
// if HAVE_CONSTRAINT == 3, then the same situation described in the case 2 
// will be reproduced, but the bound constraint will be FRowConstraint.

#define HAVE_CONSTRAINTS 2

/*--------------------------------------------------------------------------*/

// if nonzero, the Solver attached to the LPBlock is detached and re-attached
// to it at all iterations

#define DETACH_LP 0

/*--------------------------------------------------------------------------*/
// if nonzero, the two Block are not solved at every round of changes, but
// only every SKIP_BEAT + 1 rounds. this allows changes to accumulate, and
// therefore puts more pressure on the Modification handling of the Solver
// (in case this tries to do "smart" things rather than dumbly processing
// each one in turn)
//
// note that the number of rounds of changes is them multiplied by
// SKIP_BEAT + 1, so that the input parameter still dictates the number of
// Block solutions

#define SKIP_BEAT 0

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

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

#define DYNAMIC_VARS 0
// if 1, half of the variables are dynamic

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

#include "PolyhedralFunction.h"

#include "CPXMILPSolver.h"

#include "GRBMILPSolver.h"

//#include "SCIPMILPSolver.h"

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

using MultiVector = PolyhedralFunction::MultiVector;
using RealVector = PolyhedralFunction::RealVector;

using p_LF = LinearFunction *;
using p_PF = PolyhedralFunction *;

/*--------------------------------------------------------------------------*/
/*------------------------------- CONSTANTS --------------------------------*/
/*--------------------------------------------------------------------------*/

const double scale = 10;
const char *const logF = "log.bn";

const FunctionValue INF = SMSpp_di_unipi_it::Inf< FunctionValue >();

/*--------------------------------------------------------------------------*/
/*------------------------------- GLOBALS ----------------------------------*/
/*--------------------------------------------------------------------------*/

AbstractBlock * LPBlock1;   // the problem expressed as an LP ( with the 1st 
                            // solver attached )

AbstractBlock * LPBlock2;   // the problem expressed as an LP ( with the 2nd 
                            // solver attached )

bool convex = true;        // true if the PolyhedralFunction is convex

double bound = 1000;       // a tentative bound to detect unbounded instances

FunctionValue BND;         // the bound in the PolyhedralFunction (if any)

Index nvar = 10;           // number of variables
#if DYNAMIC_VARS > 0
 Index nsvar;              // number of static variables
 Index ndvar;              // number of dynamic variables
#else
 #define nsvar nvar        // all variables are static
#endif

Index m;                   // number of rows

std::mt19937 rg;           // base random generator
std::uniform_real_distribution<> dis( 0.0 , 1.0 );

MultiVector A;
RealVector b;

ColVariable * vLP1;                 // pointer to v LP1 variable
ColVariable * vLP2;                 // pointer to v LP2 variable

std::vector< ColVariable > * xLP1;  // pointer to (static) x LP1 variables
std::vector< ColVariable > * xLP2;  // pointer to (static) x LP2 variables
#if DYNAMIC_VARS > 0
 std::list< ColVariable > * xLP1d;  // pointer to (dynamic) x LP1 variables
 std::list< ColVariable > * xLP2d;  // pointer to (dynamic) x LP2 variables
#endif

#if HAVE_CONSTRAINTS == 2
 std::list< BoxConstraint > * LPbnd1;   // BoxConstraint for first LPBlock
 std::list< BoxConstraint > * LPbnd2;  // BoxConstraint for second LPBlock
#endif
#if HAVE_CONSTRAINTS == 3
 std::list< FRowConstraint > * LPbnd1;  // FRowConstrait for first LPBlock
 std::list< FRowConstraint > * LPbnd2;  // FRowConstrait for second LPBlock
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

static void ConstructfirstLPConstraint( Index i , FRowConstraint & ci ,
				   bool setblock = true )
{
 // construct constraint ci out of A[ i ] and b[ i ]:
 //
 // in the convex case, the constraint is
 //
 //          b[ i ] <= vLP1 - \sum_j Ai[ j ] * xLP1[ j ] <= INF
 //
 // in the concave case, the constraint is
 //
 //          -INF <= vLP1 - \sum_j Ai[ j ] * xLP1[ j ] <= b[ i ]
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
 vars[ j ] = std::make_pair( vLP1 , 1 );

 // then, static x
 for( ; j < nsvar ; ++j )
  vars[ j + 1 ] = std::make_pair( &((*xLP1)[ j ] ) , - A[ i ][ j ] );

 #if DYNAMIC_VARS > 0
  // finally, dynamic x
  auto xLPdit = xLP1d->begin();
  for( ; j < nvar ; ++j , ++xLPdit )
   vars[ j + 1 ] = std::make_pair( &(*xLPdit) , - A[ i ][ j ] );
 #endif

 ci.set_function( new LinearFunction( std::move( vars ) ) );
 if( setblock )
  ci.set_Block( LPBlock1 );
 }

/*--------------------------------------------------------------------------*/

static void ConstructsecondLPConstraint( Index i , FRowConstraint & ci ,
				   bool setblock = true )
{
 // construct constraint ci out of A[ i ] and b[ i ]:
 //
 // in the convex case, the constraint is
 //
 //          b[ i ] <= vLP2 - \sum_j Ai[ j ] * xLP[ j ] <= INF
 //
 // in the concave case, the constraint is
 //
 //          -INF <= vLP2 - \sum_j Ai[ j ] * xLP[ j ] <= b[ i ]
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
 vars[ j ] = std::make_pair( vLP2 , 1 );

 // then, static x
 for( ; j < nsvar ; ++j )
  vars[ j + 1 ] = std::make_pair( &((*xLP2)[ j ] ) , - A[ i ][ j ] );

 #if DYNAMIC_VARS > 0
  // finally, dynamic x
  auto xLPdit = xLP2d->begin();
  for( ; j < nvar ; ++j , ++xLPdit )
   vars[ j + 1 ] = std::make_pair( &(*xLPdit) , - A[ i ][ j ] );
 #endif

 ci.set_function( new LinearFunction( std::move( vars ) ) );
 if( setblock )
  ci.set_Block( LPBlock2 );
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

#if HAVE_CONSTRAINTS > 0

static inline void SetNN( ColVariable & LPxi1 , ColVariable & LPxi2 )
{
 if( dis( rg ) < 0.5 ) {
  LPxi1.is_positive( true , eNoMod );
  LPxi2.is_positive( true , eNoMod );
  }
 }

/*--------------------------------------------------------------------------*/

#if DYNAMIC_VARS > 0

static void RemoveBox( AbstractBlock & AB , Range rng )
{
 // the dynamic variable from the "xd" group in the Range are removed: if
 // anything is "active" in those are BoxConstraints from the "xbnd" group
 // that has to be removed as well

 auto xd = AB.get_dynamic_variable< ColVariable >( 0 );
 auto it = std::next( xd->begin() , rng.first );
 for( Index i = rng.first ; i < rng.second ; ++i , ++it ) {
  if( ! it->get_num_active() )
   continue;
  int numbox = it->get_num_active();
  for( int j = 0 ; j < numbox ; ++j ){
   std::vector< typename std::list< BoxConstraint >::iterator > rmvd;
   auto & box = *(AB.get_dynamic_constraint< BoxConstraint >( "xbnd" ));
   auto bc = dynamic_cast< BoxConstraint * >( it->get_active( 0 ) );
   if( ! bc ) {
   cout << "Unexpected stuff active in to-be-deleted Variable" << endl;
   exit( 1 );
   }
   auto to_remove = std::find_if( box.begin() , box.end() ,
			  [ bc ]( BoxConstraint & x ) {
			   return( & x == bc );
			   } );
  if( to_remove == box.end() ) {
   cout << "BoxConstraint not found" << endl;
   exit( 1 );
   }
  rmvd.push_back( to_remove );
  AB.remove_dynamic_constraints( box , rmvd ); 
  }
 }
 }

/*--------------------------------------------------------------------------*/

static void RemoveBox( AbstractBlock & AB , const Subset & sbst )
{
 // the dynamic variable from the "xd" group in the (ordered) Subset are
 // removed: if anything is "active" in those is a BoxConstraint from the
 // "xbnd" group that has to be removed as well

 auto xd = AB.get_dynamic_variable< ColVariable >( "xd" );
 Index prev = 0;
 auto it = xd->begin();
 for( auto ind : sbst ) {
  it = std::next( it , ind - prev );
  prev = ind;
  if( ! it->get_num_active() )
   continue;
  int numbox = it->get_num_active();
  for( int j = 0 ; j < numbox ; ++j ){
   std::vector< typename std::list< BoxConstraint >::iterator > rmvd;
   auto & box = *(AB.get_dynamic_constraint< BoxConstraint >( "xbnd" ));
   auto bc = dynamic_cast< BoxConstraint * >( it->get_active( 0 ) );
   if( ! bc ) {
    cout << "Unexpected stuff active in to-be-deleted Variable" << endl;
    exit( 1 );
    }
   auto to_remove = std::find_if( box.begin() , box.end() ,
			  [ bc ]( BoxConstraint & x ) {
			   return( & x == bc );
			   } );
   if( to_remove == box.end() ) {
    cout << "BoxConstraint not found" << endl;
    exit( 1 );
    }
   rmvd.push_back( to_remove );
   AB.remove_dynamic_constraints( box , rmvd ); 
   }
  }
 }

/*--------------------------------------------------------------------------*/

#endif // DYNAMIC_VARS > 0

#if HAVE_CONSTRAINTS == 2

static inline void SetBox( ColVariable & LPxi1 , ColVariable & LPxi2 )
{
 if( dis( rg ) < 0.5 ) {
  LPbnd1->resize( LPbnd1->size() + 1 );
  LPbnd2->resize( LPbnd2->size() + 1 );
  LPbnd1->back().set_variable( & LPxi1 );
  LPbnd2->back().set_variable( & LPxi2 );
  auto p = dis( rg );
  auto lhs = p < 0.666 ? 0 : -INF;
  auto rhs = p > 0.333 ? dis( rg ) : INF;
  LPbnd1->back().set_lhs( lhs , eNoMod );
  LPbnd2->back().set_lhs( lhs , eNoMod );
  LPbnd1->back().set_rhs( rhs , eNoMod );
  LPbnd2->back().set_rhs( rhs , eNoMod );
  }
 else
  SetNN( LPxi1 , LPxi2 );
 }

/*--------------------------------------------------------------------------*/

#endif // HAVE CONSTRAINT == 2

#if HAVE_CONSTRAINTS == 3

static inline void SetFRow( ColVariable & LPxi1 , ColVariable & LPxi2 )
{
 if( dis( rg ) < 0.5 ) {
  LPbnd1->resize( LPbnd1->size() + 1 );
  LPbnd2->resize( LPbnd2->size() + 1 );
  LinearFunction::v_coeff_pair vars_LP1( 1 );
  LinearFunction::v_coeff_pair vars_LP2( 1 );
  vars_LP1[ 0 ] = std::make_pair( & LPxi1 , 1 );
  vars_LP2[ 0 ] = std::make_pair( & LPxi2 , 1 );
  LPbnd1->back().set_function( new LinearFunction( std::move( vars_LP1 ) ) );
  LPbnd2->back().set_function( new LinearFunction( std::move( vars_LP2 ) ) );
  auto p = dis( rg );
  auto lhs = p < 0.666 ? 0 : -INF;
  auto rhs = p > 0.333 ? dis( rg ) : INF;
  LPbnd1->back().set_lhs( lhs , eNoMod );
  LPbnd2->back().set_lhs( lhs , eNoMod );
  LPbnd1->back().set_rhs( rhs , eNoMod );
  LPbnd2->back().set_rhs( rhs , eNoMod );
  }
 else
  SetNN( LPxi1 , LPxi2 );
 }

/*--------------------------------------------------------------------------*/

#if DYNAMIC_VARS > 0

static void RemoveFRow( AbstractBlock & AB , Range rng )
{
 // the dynamic variable from the "xd" group in the Range are removed: if
 // anything is "active" in those is a FRowConstraint from the "xbnd" group
 // that has to be removed as well

 auto xd = AB.get_dynamic_variable< ColVariable >( "xd" );
 auto it = std::next( xd->begin() , rng.first );
 for( Index i = rng.first ; i < rng.second ; ++i , ++it ) {
  if( ! it->get_num_active() )
   continue;
  int numbox = it->get_num_active();
  for( int j = 0 ; j < numbox ; ++j ){
   std::vector< typename std::list< FRowConstraint >::iterator > rmvd;
   auto & frow = *(AB.get_dynamic_constraint< FRowConstraint >( "xbnd" ));
   auto rc = dynamic_cast< FRowConstraint * >( it->get_active( 0 ) );
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

static void RemoveFRow( AbstractBlock & AB , const Subset & sbst )
{
 // the dynamic variable from the "xd" group in the (ordered) Subset are
 // removed: if anything is "active" in those is a FRowConstraint from the
 // "xbnd" group that has to be removed as well

 auto xd = AB.get_dynamic_variable< ColVariable >( "xd" );
 auto & frow = *(AB.get_dynamic_constraint< FRowConstraint >( "xbnd" ));
 std::vector< typename std::list< FRowConstraint >::iterator > rmvd;
 Index prev = 0;
 auto it = xd->begin();
 for( auto ind : sbst ) {
  it = std::next( it , ind - prev );
  prev = ind;
  if( ! it->get_num_active() )
   continue;
  int numbox = it->get_num_active();
  for( int j = 0 ; j < numbox ; ++j ){
   std::vector< typename std::list< FRowConstraint >::iterator > rmvd;
   auto & frow = *(AB.get_dynamic_constraint< FRowConstraint >( "xbnd" ));
   auto rc = dynamic_cast< FRowConstraint * >( it->get_active( 0 ) );
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

#endif // DYNAMIC_VARS > 0

#endif // HAVE_CONSTRAINT == 3

#endif // HAVE_CONSTRAINT > 0

/*--------------------------------------------------------------------------*/

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

static bool SolveBoth( void ) 
{
 try {
  // solve the first LPBlock- - - - - - - - - - - - - - - - - - - - - - - - -
  Solver * slvrLP1 = (LPBlock1->get_registered_solvers()).front();
  #if DETACH_LP
   LPBlock1->unregister_Solver( slvrLP1 );
   LPBlock1->register_Solver( slvrLP1 , true );  // push it to the front
  #endif
  int rtrnLP1 = slvrLP1->compute( false );
  bool hsLP1 = ( ( rtrnLP1 >= Solver::kOK ) && ( rtrnLP1 < Solver::kError ) )
              || ( rtrnLP1 == Solver::kLowPrecision );
  double foLP1 = hsLP1 ? ( convex ? slvrLP1->get_ub() : slvrLP1->get_lb() )
                     : ( convex ? INF : -INF );

// solve the second LPBlock- - - - - - - - - - - - - - - - - - - - - - - - -
  Solver * slvrLP2 = (LPBlock2->get_registered_solvers()).front();
  #if DETACH_LP
   LPBlock2->unregister_Solver( slvrLP2 );
   LPBlock2->register_Solver( slvrLP2 , true );  // push it to the front
  #endif
  int rtrnLP2 = slvrLP2->compute( false );
  bool hsLP2 = ( ( rtrnLP2 >= Solver::kOK ) && ( rtrnLP2 < Solver::kError ) )
              || ( rtrnLP2 == Solver::kLowPrecision );
  double foLP2 = hsLP2 ? ( convex ? slvrLP2->get_ub() : slvrLP2->get_lb() )
                     : ( convex ? INF : -INF );

  if( hsLP1 && hsLP2 && ( abs( foLP1 - foLP2 ) <= 2e-7 *
			 max( double( 1 ) , abs( max( foLP1 , foLP2 ) ) ) ) ) {
   LOG1( "OK(f)" << endl );
   return( true );
   }

  if( ( rtrnLP1 == Solver::kInfeasible ) &&
      ( rtrnLP2== Solver::kInfeasible ) ) {
    LOG1( "OK(?e?)" << endl );
    return( true );
    }

  if( ( rtrnLP1 == Solver::kUnbounded ) &&
      ( rtrnLP2 == Solver::kUnbounded ) ) {
   LOG1( "OK(u)" << endl );
   return( true );
   }

  #if( LOG_LEVEL >= 1 )
   cout << "LPBlock1 = ";
   if( hsLP1 )
    cout << foLP1;
   else
    if( rtrnLP1 == Solver::kInfeasible )
     cout << "    Unfeas(?)";
    else
     if( rtrnLP1 == Solver::kUnbounded )
      cout << "      Unbounded";
     else
      cout << "      Error!";

   cout << "LPBlock2 = ";
   if( hsLP2 )
    cout << foLP2;
   else
    if( rtrnLP2 == Solver::kInfeasible )
     cout << "    Unfeas(?)";
    else
     if( rtrnLP2 == Solver::kUnbounded )
      cout << "      Unbounded";
     else
      cout << "      Error!";
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
           "       wchg: what to change, coded bit-wise [127]"
		<< endl <<
           "             0 = add rows, 1 = delete rows "
		<< endl <<
           "             2 = modify rows, 3 = modify constants"
		<< endl <<
           "             4 = change global lower/upper bound"
          #if DYNAMIC_VARS > 0
		<< endl <<
           "             5 = add variables, 6 = delete variables"
	  #endif
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

 #if DYNAMIC_VARS > 0
  nsvar = nvar / 2;      // half of the variables are dynamic
  ndvar = nvar - nsvar;  // the other half are static
 #endif

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

 // construct the first LP- - - - - - - - - - - - - - - - - - - - - - - - - -
 {
  // ensure all original pointers go out of scope immediately after that
  // the construction has finished

  LPBlock1 = new AbstractBlock();

  // construct the Variable
  xLP1 = new std::vector< ColVariable >( nsvar );
  #if DYNAMIC_VARS > 0
   xLP1d = new std::list< ColVariable >( ndvar );
  #endif

  vLP1 = new ColVariable;
  vLP1->set_Block( LPBlock1 );

  // construct the m dynamic Constraint
  auto ALP1 = new std::list< FRowConstraint >( m );
  auto ALP1it = ALP1->begin();
  for( Index i = 0 ; i < m ; )
   ConstructfirstLPConstraint( i++ , *(ALP1it++) );

  // construct the static lower bound Constraint
  auto LB1c = new BoxConstraint( LPBlock1 , vLP1 , -INF , INF );
  if( BND != INF ) {
   if( convex )
    LB1c->set_lhs( -BND );
   else
    LB1c->set_rhs( BND );
   }

  // construct the Objective
  auto objLP1 = new FRealObjective();
  objLP1->set_function( new LinearFunction( { std::make_pair( vLP1 , 1 ) } ) );
  objLP1->set_sense( convex ? Objective::eMin : Objective::eMax , eNoMod );
  
  // now set the Variable, Constraint and Objective in the AbstractBlock
  LPBlock1->add_static_variable( *vLP1 , "v" );
  LPBlock1->add_static_variable( *xLP1 , "x" );
  #if DYNAMIC_VARS > 0
   LPBlock1->add_dynamic_variable( *xLP1d , "xd" );
  #endif
  LPBlock1->add_dynamic_constraint( *ALP1 , "cuts" );
  LPBlock1->add_static_constraint( *LB1c , "vbnd" );
  LPBlock1->set_objective( objLP1 );
  }

 // construct the second LP- - - - - - - - - - - - - - - - - - - - - - - - - -
 {
  // ensure all original pointers go out of scope immediately after that
  // the construction has finished

  LPBlock2 = new AbstractBlock();

  // construct the Variable
  xLP2 = new std::vector< ColVariable >( nsvar );
  #if DYNAMIC_VARS > 0
   xLP2d = new std::list< ColVariable >( ndvar );
  #endif

  vLP2 = new ColVariable;
  vLP2->set_Block( LPBlock2 );

  // construct the m dynamic Constraint
  auto ALP2 = new std::list< FRowConstraint >( m );
  auto ALP2it = ALP2->begin();
  for( Index i = 0 ; i < m ; )
   ConstructsecondLPConstraint( i++ , *(ALP2it++) );

  // construct the static lower bound Constraint
  auto LB2c = new BoxConstraint( LPBlock2 , vLP2 , -INF , INF );
  if( BND != INF ) {
   if( convex )
    LB2c->set_lhs( -BND );
   else
    LB2c->set_rhs( BND );
   }

  // construct the Objective
  auto objLP2 = new FRealObjective();
  objLP2->set_function( new LinearFunction( { std::make_pair( vLP2 , 1 ) } ) );
  objLP2->set_sense( convex ? Objective::eMin : Objective::eMax , eNoMod );
  
  // now set the Variable, Constraint and Objective in the AbstractBlock
  LPBlock2->add_static_variable( *vLP2 , "v" );
  LPBlock2->add_static_variable( *xLP2 , "x" );
  #if DYNAMIC_VARS > 0
   LPBlock2->add_dynamic_variable( *xLP2d , "xd" );
  #endif
  LPBlock2->add_dynamic_constraint( *ALP2 , "cuts" );
  LPBlock2->add_static_constraint( *LB2c , "vbnd" );
  LPBlock2->set_objective( objLP2 );
  }

 // define bound constraints- - - - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 #if HAVE_CONSTRAINTS == 1
 {
  auto & LP1x = *(LPBlock1->get_static_variable_v< ColVariable >( "x" ));
  auto & LP2x = *(LPBlock2->get_static_variable_v< ColVariable >( "x" ));
  for( Index i = 0 ; i < nsvar ; ++i )
   SetNN( LP1x[ i ] , LP2x[ i ] );
  #if DYNAMIC_VARS > 0
   auto LP1xd = LPBlock1->get_dynamic_variable< ColVariable >( "xd" )->begin();
   auto LP2xd = LPBlock2->get_dynamic_variable< ColVariable >( "xd")->begin();
   for( Index i = 0 ; i < ndvar ; ++i )
    SetNN( *(LP1xd++) , *(LP2xd++) );
  #endif
  }
 #endif
 #if HAVE_CONSTRAINTS == 2
 {
  LPbnd1 = new std::list< BoxConstraint >;
  LPbnd2 = new std::list< BoxConstraint >;
  auto & LP1x = *(LPBlock1->get_static_variable_v< ColVariable >( "x" ));
  auto & LP2x = *(LPBlock2->get_static_variable_v< ColVariable >( "x" ));
  for( Index i = 0 ; i < nsvar ; ++i )
   SetBox( LP1x[ i ] , LP2x[ i ] );
  #if DYNAMIC_VARS > 0
   auto LP1xd = LPBlock1->get_dynamic_variable< ColVariable >( "xd" )->begin();
   auto LP2xd = LPBlock2->get_dynamic_variable< ColVariable >( "xd" )->begin();
   for( Index i = 0 ; i < ndvar ; ++i )
    SetBox( *(LP1xd++) , *(LP2xd++) );
  #endif

  // note: the list may be empty, but it is intentionally added anyway
  LPBlock1->add_dynamic_constraint( *LPbnd1 , "xbnd" );
  LPBlock2->add_dynamic_constraint( *LPbnd2 , "xbnd" );
  }
 #endif
 #if HAVE_CONSTRAINTS == 3
 {
  LPbnd1 = new std::list< FRowConstraint >;
  LPbnd2 = new std::list< FRowConstraint >;
  auto & LP1x = *(LPBlock1->get_static_variable_v< ColVariable >( "x" ));
  auto & LP2x = *(LPBlock2->get_static_variable_v< ColVariable >( "x" ));
  for( Index i = 0 ; i < nsvar ; ++i )
   SetFRow( LP1x[ i ] , LP2x[ i ] );
  #if DYNAMIC_VARS > 0
   auto LP1xd = LPBlock1->get_dynamic_variable< ColVariable >( "xd" )->begin();
   auto LP2xd = LPBlock2->get_dynamic_variable< ColVariable >( "xd" )->begin();
   for( Index i = 0 ; i < ndvar ; ++i )
    SetFRow( *(LP1xd++) , *(LP2xd++) );
  #endif

  // note: the list may be empty, but it is intentionally added anyway
  LPBlock1->add_dynamic_constraint( *LPbnd1 , "xbnd" );
  LPBlock2->add_dynamic_constraint( *LPbnd2 , "xbnd" );
  }
 #endif
 
 // attach the Solver to the Block- - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // Here you can decide which solver attach to specific block - - - - - - - -

 Solver * solver1 = new CPXMILPSolver();
 LPBlock1->register_Solver( solver1 );

 //Solver * solver1 = new GRBMILPSolver();
 //LPBlock1->register_Solver( solver1 );

 //Solver * solver1 = new SCIPMILPSolver();
 //LPBlock1->register_Solver( solver1 );

 Solver * solver2 = new GRBMILPSolver();
 LPBlock2->register_Solver( solver2 );

 //Solver * solver2 = new CPXMILPSolver();
 //LPBlock2->register_Solver( solver2 );

 //Solver * solver2 = new SCIPMILPSolver();
 //LPBlock2->register_Solver( solver2 );

 // open log-file - - - - - - - - - - -  - - - - - - - - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 #if( LOG_LEVEL >= 3 )
  ((LPBlock1->get_registered_solvers()).front())->set_par(
                            MILPSolver::strOutputFile , "LP1Block.lp" );
  ((LPBlock2->get_registered_solvers()).front())->set_par(
                            MILPSolver::strOutputFile , "LP2Block.lp" );
 #endif

 // first solver call - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 LOG1( "First call: " );

 bool AllPassed = SolveBoth();
 
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
    vLP1 = LPBlock1->get_static_variable< ColVariable >( "v" );
    xLP1 = LPBlock1->get_static_variable_v< ColVariable >( "x" );
    #if DYNAMIC_VARS > 0
     xLP1d = LPBlock1->get_dynamic_variable< ColVariable >( "xd" );
    #endif

    std::list< FRowConstraint > nc1( tochange );
    auto nc1it = nc1.begin();
    for( Index i = 0 ; i < tochange ; )
     ConstructfirstLPConstraint( i++ , *(nc1it++) );
    auto cnst1 = LPBlock1->get_dynamic_constraint< FRowConstraint >( "cuts" );
    LPBlock1->add_dynamic_constraints( *cnst1 , nc1 );

    // add them to the second LP
    vLP2 = LPBlock2->get_static_variable< ColVariable >( "v" );
    xLP2 = LPBlock2->get_static_variable_v< ColVariable >( "x" );
    #if DYNAMIC_VARS > 0
     xLP2d = LPBlock2->get_dynamic_variable< ColVariable >( "xd" );
    #endif

    std::list< FRowConstraint > nc2( tochange );
    auto nc2it = nc2.begin();
    for( Index i = 0 ; i < tochange ; )
     ConstructsecondLPConstraint( i++ , *(nc2it++) );
    auto cnst2 = LPBlock2->get_dynamic_constraint< FRowConstraint >( "cuts" );
    LPBlock2->add_dynamic_constraints( *cnst2 , nc2 );


    // update m
    m += tochange;

    // sanity checks
    PANIC( m == cnst1->size() );
    PANIC( m == cnst2->size() );
    }

  // delete rows- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  if( ( wchg & 2 ) && ( dis( rg ) <= p_change ) )
   if( Index tochange = min( m - 1 , Index( dis( rg ) * n_change ) ) ) {
    LOG1( "deleted " << tochange << " rows" );

    auto cnst1 = LPBlock1->get_dynamic_constraint< FRowConstraint >( "cuts" );
    auto cnst2 = LPBlock2->get_dynamic_constraint< FRowConstraint >( "cuts" );
    
    if( dis( rg ) <= 0.5 ) {  // in 50% of the cases do a ranged change
     LOG1( "(r) - " );

     Index strt = dis( rg ) * ( m - tochange );
     Index stp = strt + tochange;

     // remove them from the first LP
     LPBlock1->remove_dynamic_constraints( *cnst1 , Range( strt , stp ) );

     // remove them from the second LP
     LPBlock2->remove_dynamic_constraints( *cnst2 , Range( strt , stp ) );
     }
    else {  // in the other 50% of the cases, do a sparse change
     LOG1( "(s) - " );
     Subset nms( GenerateRand( m , tochange ) );

     // remove them from the first LP
     if( tochange == 1 )
      LPBlock1->remove_dynamic_constraint( *cnst1 , std::next( cnst1->begin() ,
							     nms[ 0 ] ) );
     else
      LPBlock1->remove_dynamic_constraints( *cnst1 , Subset( nms ) , true );
    
     // remove them from the second LP
     if( tochange == 1 )
      LPBlock2->remove_dynamic_constraint( *cnst2 , std::next( cnst2->begin() ,
							     nms[ 0 ] ) );
     else
      LPBlock2->remove_dynamic_constraints( *cnst2 , Subset( nms ) , true );
     }

    // update m
    m -= tochange;

    // sanity checks
    PANIC( m == cnst1->size() );
    PANIC( m == cnst2->size() );
    }

  // modify rows- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  if( ( wchg & 4 ) && ( dis( rg ) <= p_change ) )
   if( Index tochange = std::min( m , Index( dis( rg ) * n_change ) ) ) {
    LOG1( "modified " << tochange << " rows" );

    GenerateAb( tochange , nvar );

    vLP1 = LPBlock1->get_static_variable< ColVariable >( "v" );
    xLP1 = LPBlock1->get_static_variable_v< ColVariable >( "x" );
    #if DYNAMIC_VARS > 0
     xLP1d = LPBlock1->get_dynamic_variable< ColVariable >( "xd" );
    #endif
    auto cnst1 = LPBlock1->get_dynamic_constraint< FRowConstraint >( "cuts" );

    vLP2 = LPBlock2->get_static_variable< ColVariable >( "v" );
    xLP2 = LPBlock2->get_static_variable_v< ColVariable >( "x" );
    #if DYNAMIC_VARS > 0
     xLP2d = LPBlock2->get_dynamic_variable< ColVariable >( "xd" );
    #endif
    auto cnst2 = LPBlock2->get_dynamic_constraint< FRowConstraint >( "cuts" );

    if( dis( rg ) <= 0.5 ) {  // in 50% of the cases do a ranged change
     LOG1( "(r) - " );

     Index strt = dis( rg ) * ( m - tochange );
     Index stp = strt + tochange;

     // send all the Modification to the same channel
     Observer::ChnlName chnl1 = LPBlock1->open_channel();
     const auto iAM1 = Observer::make_par( eModBlck , chnl1 );

     // modify them in the first LP
     auto cit1 = std::next( cnst1->begin() , strt );
     for( Index i = 0 ; i < tochange ; ++i )
      ChangeLPConstraint( i , *(cit1++) , iAM1 );

     LPBlock1->close_channel( chnl1 );  // close the channel

     // send all the Modification to the same channel
     Observer::ChnlName chnl2 = LPBlock2->open_channel();
     const auto iAM2 = Observer::make_par( eModBlck , chnl2 );

     // modify them in the second LP
     auto cit2 = std::next( cnst2->begin() , strt );
     for( Index i = 0 ; i < tochange ; ++i )
      ChangeLPConstraint( i , *(cit2++) , iAM2 );

     LPBlock2->close_channel( chnl2 );  // close the channel
     }
    else {  // in the other 50% of the cases, do a sparse change
     LOG1( "(s) - " );
     Subset nms( GenerateRand( m , tochange ) );

     // send all the Modification to the same channel
     Observer::ChnlName chnl1 = LPBlock1->open_channel();
     const auto iAM1 = Observer::make_par( eModBlck , chnl1 );

     // modify them in the first LP
     Index prev = 0;
     auto cit1 = cnst1->begin();
     for( Index i = 0 ; i < tochange ; ++i ) {
      cit1 = std::next( cit1 , nms[ i ] - prev );
      prev = nms[ i ];
      ChangeLPConstraint( i , *cit1 , iAM1 );
      }

     LPBlock1->close_channel( chnl1 );  // close the channel

     // send all the Modification to the same channel
     Observer::ChnlName chnl2 = LPBlock2->open_channel();
     const auto iAM2 = Observer::make_par( eModBlck , chnl2 );

     // modify them in the first LP
     prev = 0;
     auto cit2 = cnst2->begin();
     for( Index i = 0 ; i < tochange ; ++i ) {
      cit2 = std::next( cit2 , nms[ i ] - prev );
      prev = nms[ i ];
      ChangeLPConstraint( i , *cit2 , iAM2 );
      }

     LPBlock2->close_channel( chnl2 );  // close the channel
    }
   }

  // modify constants - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  if( ( wchg & 8 ) && ( dis( rg ) <= p_change ) )
   if( Index tochange = std::min( m , Index( dis( rg ) * n_change ) ) ) {
    LOG1( "modified " << tochange << " constants" );

    Generateb( tochange );
     
    auto cnst1 = LPBlock1->get_dynamic_constraint< FRowConstraint >( "cuts" );
    auto cnst2 = LPBlock2->get_dynamic_constraint< FRowConstraint >( "cuts" );


    if( dis( rg ) <= 0.5 ) {  // in 50% of the cases do a ranged change
     LOG1( "(r) - " );

     Index strt = dis( rg ) * ( m - tochange );
     Index stp = strt + tochange;

     // change them in the first LP
     auto cit1 = std::next( cnst1->begin() , strt );
     if( convex )
      for( Index i = 0 ; i < tochange ; )
       (*(cit1++)).set_lhs( b[ i++ ] );
     else
      for( Index i = 0 ; i < tochange ; )
       (*(cit1++)).set_rhs( b[ i++ ] );
 
     // change them in the second LP
     auto cit2 = std::next( cnst2->begin() , strt );
     if( convex )
      for( Index i = 0 ; i < tochange ; )
       (*(cit2++)).set_lhs( b[ i++ ] );
     else
      for( Index i = 0 ; i < tochange ; )
       (*(cit2++)).set_rhs( b[ i++ ] );
     }
    else {  // in the other 50% of the cases, do a sparse change
     LOG1( "(s) - " );
     Subset nms( GenerateRand( m , tochange ) );

     // change them in the first LP
     Index prev = 0;
     auto cit1 = cnst1->begin();
     if( convex )
      for( Index i = 0 ; i < tochange ; ) {
       cit1 = std::next( cit1 , nms[ i ] - prev );
       prev = nms[ i ];
       (*cit1).set_lhs( b[ i++ ] );
       }
     else
      for( Index i = 0 ; i < tochange ; ) {
       cit1 = std::next( cit1 , nms[ i ] - prev );
       prev = nms[ i ];
       (*cit1).set_rhs( b[ i++ ] );
       }

     // change them in the second LP
     prev = 0;
     auto cit2 = cnst2->begin();
     if( convex )
      for( Index i = 0 ; i < tochange ; ) {
       cit2 = std::next( cit2 , nms[ i ] - prev );
       prev = nms[ i ];
       (*cit2).set_lhs( b[ i++ ] );
       }
     else
      for( Index i = 0 ; i < tochange ; ) {
       cit2 = std::next( cit2 , nms[ i ] - prev );
       prev = nms[ i ];
       (*cit2).set_rhs( b[ i++ ] );
       }
     }
    }

  // modify bound - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  if( ( wchg & 16 ) && ( dis( rg ) <= p_change ) ) {
   LOG1( "modified bound - " );

   GenerateBND();

   // change it in the first LP
   auto cnst1 = LPBlock1->get_static_constraint< BoxConstraint >( "vbnd" );
   if( convex )
    cnst1->set_lhs( -BND );
   else
    cnst1->set_rhs( BND );

   // change it in the second LP
   auto cnst2 = LPBlock2->get_static_constraint< BoxConstraint >( "vbnd" );
   if( convex )
    cnst2->set_lhs( -BND );
   else
    cnst2->set_rhs( BND );
   }

 // add variables- - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  #if DYNAMIC_VARS > 0
  if( ( wchg & 32 ) && ( dis( rg ) <= p_change ) ) {
   Index tochange = std::max( Index( 1 ) , Index( dis( rg ) * nsvar / 4 ) );
   LOG1( "added " << tochange << " variables - " );

   GenerateA( m , tochange );

   // add them in the first LP
   std::list< ColVariable > nxLP1d( tochange );
   std::vector< ColVariable * > nxp1( tochange );
   auto nxlp1it = nxLP1d.begin();
   for( Index i = 0 ; i < tochange ; )
    nxp1[ i++ ] = &(*(nxlp1it++));

   LPBlock1->add_dynamic_variables(
	   *(LPBlock1->get_dynamic_variable< ColVariable >( "xd" )) , nxLP1d );

   auto cnst1_it =
        LPBlock1->get_dynamic_constraint< FRowConstraint >( "cuts" )->begin();
   if( tochange == 1 )
    for( Index i = 0 ; i < m ; ++i ) {
     auto fi1 = static_cast< p_LF >( (cnst1_it++)->get_function() );
     fi1->add_variable( nxp1[ 0 ] , - A[ i ][ 0 ] );
     }
   else
    for( Index i = 0 ; i < m ; ++i ) {
     auto fi1 = static_cast< p_LF >( (cnst1_it++)->get_function() );
     LinearFunction::v_coeff_pair ncp1( tochange );
     for( Index j = 0 ; j < ncp1.size() ; ++j ) {
      ncp1[ j ].first = nxp1[ j ];
      ncp1[ j ].second = - A[ i ][ j ];
      }
     fi1->add_variables( std::move( ncp1 ) );
     }

   // add them in the second LP
   std::list< ColVariable > nxLP2d( tochange );
   std::vector< ColVariable * > nxp2( tochange );
   auto nxlp2it = nxLP2d.begin();
   for( Index i = 0 ; i < tochange ; )
    nxp2[ i++ ] = &(*(nxlp2it++));

   LPBlock2->add_dynamic_variables(
	   *(LPBlock2->get_dynamic_variable< ColVariable >( "xd" )) , nxLP2d );

   auto cnst2_it =
        LPBlock2->get_dynamic_constraint< FRowConstraint >( "cuts" )->begin();
   if( tochange == 1 )
    for( Index i = 0 ; i < m ; ++i ) {
     auto fi2 = static_cast< p_LF >( (cnst2_it++)->get_function() );
     fi2->add_variable( nxp2[ 0 ] , - A[ i ][ 0 ] );
     }
   else
    for( Index i = 0 ; i < m ; ++i ) {
     auto fi2 = static_cast< p_LF >( (cnst2_it++)->get_function() );
     LinearFunction::v_coeff_pair ncp2( tochange );
     for( Index j = 0 ; j < ncp2.size() ; ++j ) {
      ncp2[ j ].first = nxp2[ j ];
      ncp2[ j ].second = - A[ i ][ j ];
      }
     fi2->add_variables( std::move( ncp2 ) );
     }

   // generate bound constraints
   #if HAVE_CONSTRAINTS > 0
    auto & LP1xd = *(LPBlock1->get_dynamic_variable< ColVariable >( "xd" ));
    auto LP1xd_it = LP1xd.begin();
    auto & LP2xd = *(LPBlock2->get_dynamic_variable< ColVariable >( "xd" ));
    auto LP2xd_it = LP2xd.begin();
    std::next( LP1xd_it , ndvar );
    std::next( LP2xd_it , ndvar );

    #if HAVE_CONSTRAINTS == 1
     for( ; LP1xd_it != LP1xd.end() ; )
      SetNN( *(LP1xd_it++) , *(LP2xd_it++) );
    #endif
    #if HAVE_CONSTRAINTS == 2
     LPbnd1 = new std::list< BoxConstraint >;
     LPbnd2 = new std::list< BoxConstraint >;

     for( ; LP1xd_it != LP1xd.end() ; )
      SetBox( *(LP1xd_it++) , *(LP2xd_it++) );

     if( ! LPbnd1->empty() ) {
      LPBlock1->add_dynamic_constraints(
	 *(LPBlock1->get_dynamic_constraint< BoxConstraint >( "xbnd" )) ,
	 *LPbnd1 );
      LPBlock2->add_dynamic_constraints(
	 *(LPBlock2->get_dynamic_constraint< BoxConstraint >( "xbnd" )) ,
	 *LPbnd2 );
      }
    #endif
    #if HAVE_CONSTRAINTS == 3
     LPbnd1 = new std::list< FRowConstraint >;
     LPbnd2 = new std::list< FRowConstraint >;

     for( ; LP1xd_it != LP1xd.end() ; )
      SetFRow( *(LP1xd_it++) , *(LP2xd_it++) );

     if( ! LPbnd1->empty() ) {
      LPBlock1->add_dynamic_constraints(
	 *(LPBlock1->get_dynamic_constraint< FRowConstraint >( "xbnd" )) ,
	 *LPbnd1 );
      LPBlock2->add_dynamic_constraints(
	 *(LPBlock2->get_dynamic_constraint< FRowConstraint >( "xbnd" )) ,
	 *LPbnd2 );
      }
    #endif
   #endif

   // update nvar and ndvar
   nvar += tochange;
   ndvar += tochange;

   // sanity checks
   PANIC( ndvar ==
	      LPBlock1->get_dynamic_variable< ColVariable >( "xd" )->size() );
   for( auto & ci :
 	    *(LPBlock1->get_dynamic_constraint< FRowConstraint >( "cuts" )) )
    PANIC( nvar + 1 == ci.get_num_active_var() );

    PANIC( ndvar ==
	      LPBlock2->get_dynamic_variable< ColVariable >( "xd" )->size() );
   for( auto & ci :
 	    *(LPBlock2->get_dynamic_constraint< FRowConstraint >( "cuts" )) )
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

     // remove them from the first LP
     auto xLP1d = LPBlock1->get_dynamic_variable< ColVariable >( 0 );
     auto cnst1_it =
             LPBlock1->get_dynamic_constraint< FRowConstraint >( 0 )->begin();
     if( tochange == 1 )
      for( Index i = 0 ; i < m ; ++i ) {
       auto fi1 = static_cast< p_LF >( (cnst1_it++)->get_function() );
       fi1->remove_variable( strt + nsvar + 1 );
       }
     else
      for( Index i = 0 ; i < m ; ++i ) {
       auto fi1 = static_cast< p_LF >( (cnst1_it++)->get_function() );
       fi1->remove_variables( Range( strt + nsvar + 1 , stp + nsvar + 1 ) );
       }

     // remove them from the second LP
     auto xLP2d = LPBlock2->get_dynamic_variable< ColVariable >( 0 );
     auto cnst2_it =
             LPBlock2->get_dynamic_constraint< FRowConstraint >( 0 )->begin();
     if( tochange == 1 )
      for( Index i = 0 ; i < m ; ++i ) {
       auto fi2 = static_cast< p_LF >( (cnst2_it++)->get_function() );
       fi2->remove_variable( strt + nsvar + 1 );
       }
     else
      for( Index i = 0 ; i < m ; ++i ) {
       auto fi2 = static_cast< p_LF >( (cnst2_it++)->get_function() );
       fi2->remove_variables( Range( strt + nsvar + 1 , stp + nsvar + 1 ) );
      }
    
     #if HAVE_CONSTRAINTS == 2
      // the variables can now only be active in the associated box
      // constraint, if any: exploit this to identify the box constraint
      // and remove it
      RemoveBox( *LPBlock1 , Range( strt , stp ) );
      RemoveBox( *LPBlock2 , Range( strt , stp ) );
     #endif
     #if HAVE_CONSTRAINTS == 3
      // the variables can now only be active in the associated frow
      // constraint, if any: exploit this to identify the frow constraint
      // and remove it
      RemoveFRow( *LPBlock1 , Range( strt , stp ) );
      RemoveFRow( *LPBlock2 , Range( strt , stp ) );
     #endif
     
     LPBlock1->remove_dynamic_variables( *xLP1d , Range( strt , stp ) );
     LPBlock2->remove_dynamic_variables( *xLP2d , Range( strt , stp ) );
     }
    else {  // in the other 50% of the cases, do a sparse change
     LOG1( "(s) - " );
     Subset nms( GenerateRand( ndvar , tochange ) );

     // remove them from the first LP
     auto xLP1d = LPBlock1->get_dynamic_variable< ColVariable >( 0 );
     auto cnst1_it =
             LPBlock1->get_dynamic_constraint< FRowConstraint >( 0 )->begin();
     auto xLP2d = LPBlock2->get_dynamic_variable< ColVariable >( 0 );
     auto cnst2_it =
             LPBlock2->get_dynamic_constraint< FRowConstraint >( 0 )->begin();
     if( tochange == 1 ) {
      for( Index i = 0 ; i < m ; ++i ) {
       auto fi1 = static_cast< p_LF >( (cnst1_it++)->get_function() );
       fi1->remove_variable( nms[ 0 ] + nsvar + 1 );
       
       auto fi2 = static_cast< p_LF >( (cnst2_it++)->get_function() );
       fi2->remove_variable( nms[ 0 ] + nsvar + 1 );
       }

      #if HAVE_CONSTRAINTS == 2
       // the variables can now only be active in the associated box
       // constraint, if any: exploit this to identify the box constraint
       // and remove it
       RemoveBox( *LPBlock1 , Range( nms[ 0 ] , nms[ 0 ] + 1 ) );
       RemoveBox( *LPBlock2 , Range( nms[ 0 ] , nms[ 0 ] + 1 ) );
      #endif
      #if HAVE_CONSTRAINTS == 3
       // the variables can now only be active in the associated frow
       // constraint, if any: exploit this to identify the frow constraint
       // and remove it
       RemoveFRow( *LPBlock1 , Range( nms[ 0 ] , nms[ 0 ] + 1 ) );
       RemoveFRow( *LPBlock2 , Range( nms[ 0 ] , nms[ 0 ] + 1 ) );
      #endif

      auto vp1 = std::next( xLP1d->begin() , nms[ 0 ] );
      LPBlock1->remove_dynamic_variable( *xLP1d , vp1 );
      auto vp2 = std::next( xLP2d->begin() , nms[ 0 ] );
      LPBlock2->remove_dynamic_variable( *xLP2d , vp2 );
      }
     else {
      for( Index i = 0 ; i < m ; ++i ) {
       auto fi1 = static_cast< p_LF >( (cnst1_it++)->get_function() );
       auto fi2 = static_cast< p_LF >( (cnst2_it++)->get_function() );

       Subset nms1( nms );
       Subset nms2( nms );
       for( auto & n1i : nms1 )
	      n1i = n1i + nsvar + 1;
       for( auto & n2i : nms2 )
	      n2i = n2i + nsvar + 1;

       fi1->remove_variables( std::move( nms1 ) , true );
       fi2->remove_variables( std::move( nms2 ) , true );
       }

      #if HAVE_CONSTRAINTS == 2
       // the variables can now only be active in the associated box
       // constraint, if any: exploit this to identify the box constraint
       // and remove it
       RemoveBox( *LPBlock1 , nms );
       RemoveBox( *LPBlock2 , nms );
      #endif
      #if HAVE_CONSTRAINTS == 3
       // the variables can now only be active in the associated frow
       // constraint, if any: exploit this to identify the frow constraint
       // and remove it
       RemoveFRow( *LPBlock1 , nms );
       RemoveFRow( *LPBlock2 , nms );
      #endif
      
      LPBlock1->remove_dynamic_variables( *xLP1d , Subset( nms ) );
      LPBlock2->remove_dynamic_variables( *xLP2d , Subset( nms ) );
     }
    }

    // update ndvar and nvar
    ndvar -= tochange;
    nvar -= tochange;

    // sanity checks
    PANIC( ndvar ==
	         LPBlock1->get_dynamic_variable< ColVariable >( 0 )->size() );
    PANIC( ndvar ==
	         LPBlock2->get_dynamic_variable< ColVariable >( 0 )->size() );
    for( auto & ci :
	          *(LPBlock1->get_dynamic_constraint< FRowConstraint >( 0 )) )
     PANIC( nvar + 1 == ci.get_num_active_var() );
    for( auto & ci :
	          *(LPBlock2->get_dynamic_constraint< FRowConstraint >( 0 )) )
     PANIC( nvar + 1 == ci.get_num_active_var() );
    }
  #endif

  // if verbose, print out stuff- - - - - - - - - - - - - - - - - - - - - - -

  #if( LOG_LEVEL >= 3 )
   ((LPBlock1->get_registered_solvers()).front())->set_par(
		                     MILPSolver::strOutputFile , "LP1Block-" +
		                     std::to_string( rep ) + ".lp" );
    ((LPBlock2->get_registered_solvers()).front())->set_par(
		                     MILPSolver::strOutputFile , "LP2Block-" +
		                     std::to_string( rep ) + ".lp" );
  #endif

  // finally, re-solve the problems- - - - - - - - - - - - - - - - - - - - -
  // ... every SKIP_BEAT + 1 rounds

  if( ! ( ++rep % ( SKIP_BEAT + 1 ) ) )
   AllPassed &= SolveBoth();
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

 // delete the Blocks
 delete( LPBlock1 );
 delete( LPBlock2 );


 // terminate - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 return( AllPassed ? 0 : 1 );

 }  // end( main )

/*--------------------------------------------------------------------------*/
/*------------------------ End File test.cpp -------------------------------*/
/*--------------------------------------------------------------------------*/