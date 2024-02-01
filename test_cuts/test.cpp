/*--------------------------------------------------------------------------*/
/*-------------------------- File test.cpp ---------------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Main for testing lazy constraints separation in :MILPSolver
 *
 * Given one input parameter n, a NCoCubeBlock with that size is defined.
 * The NCoCubeBlock is a Block that encodes for the N-Co-Cube, i.e., the
 * unitary ball in the L_1 norm, with a simple linear Objective. Optimizing
 * over this is trivial, since the N-Co-Cube only has 2n vertices of the
 * form [ 0 , 0 , ... \pm 1 , ... , 0 ], and therefore also has the
 * integrality property. However, NCoCubeBlock implements the "crazy"
 * formulation in the original variable space which has 2^n constraints of
 * the form
 *
 *   \pm 1 x_1 + \pm 1 x_2 + ... + \pm 1 x_n <= 1
 *
 * Of course these are implemented as dynamic constraint, whose separation
 * (trivial and linear) is performed in generate_dynamic_constraint().
 *
 * The problem is repeatedly solved a number of times randomly changing the
 * Objective (which is the only thing that can change), and the optimal
 * value is compared with that obtained by the trivial optimization.
 *
 * \author Antonio Frangioni \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \copyright &copy; by Antonio Frangioni
 */
/*--------------------------------------------------------------------------*/
/*-------------------------------- MACROS ----------------------------------*/
/*--------------------------------------------------------------------------*/

#define LOG_LEVEL 0
// 0 = only pass/fail
// 1 = result of each test
// 2 = + solver log
// 3 = + save LP file

#if( LOG_LEVEL >= 1 )
 #define LOG1( x ) cout << x
 #define CLOG1( y , x ) if( y ) cout << x

 #if( LOG_LEVEL >= 2 )
  #define LOG_ON_COUT 1
  // if nonzero, the NDO Solver log is sent on cout rather than on a file
 #endif
#else
 #define LOG1( x )
 #define CLOG1( y , x )
#endif

/*--------------------------------------------------------------------------*/
// if nonzero, the :MILPSolver attached to the NCoCubeBlock is detached and
// re-attached to it at all iterations

#define DETACH_LP 0

/*--------------------------------------------------------------------------*/

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

// #include "AbstractBlock.h"

#include "BlockSolverConfig.h"

#include "LinearFunction.h"

#include "FRowConstraint.h"

#include "FRealObjective.h"

#if( LOG_LEVEL >= 3 )
 #include "MILPSolver.h"
#endif

/*--------------------------------------------------------------------------*/
/*-------------------------------- USING -----------------------------------*/
/*--------------------------------------------------------------------------*/

using namespace std;

using namespace SMSpp_di_unipi_it;

/*--------------------------------------------------------------------------*/
/*---------------------------- SIMPLE TYPES --------------------------------*/
/*--------------------------------------------------------------------------*/

using Index = Block::Index;

using Range = Block::Range;

using Subset = Block::Subset;

using FunctionValue = Function::FunctionValue;
using c_FunctionValue = Function::c_FunctionValue;

using v_coeff_pair = LinearFunction::v_coeff_pair;

// using p_AB = AbstractBlock *;
using p_LF = LinearFunction *;

/*--------------------------------------------------------------------------*/
/*------------------------- CLASS NCoCubeBlock -----------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------- GENERAL NOTES --------------------------------*/
/*--------------------------------------------------------------------------*/
/// an N-Co-Cube, i.e., the unitary ball in the L_1 norm
/** NCoCubeBlock is a Block that encodes for the N-Co-Cube, i.e., the unitary
 * ball in the L_1 norm, with a simple linear Objective. Optimizing over this
 * is trivial, since the N-Co-Cube only has 2n vertices of the form
 * [ 0 , 0 , ... \pm 1 , ... , 0 ], and therefore also has the integrality
 * property. However, NCoCubeBlock implements the "crazy" formulation in the
 * original variable space which has 2^n constraints of  the form
 *
 *   \pm 1 x_1 + \pm 1 x_2 + ... + \pm 1 x_n <= 1
 *
 * Of course these are implemented as dynamic constraint, whose separation
 * (trivial and linear) is performed in generate_dynamic_constraint(). */

class NCoCubeBlock : public Block {

/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/

public:

/*--------------------- PUBLIC METHODS OF THE CLASS ------------------------*/
/*--------------------- CONSTRUCTOR AND DESTRUCTOR -------------------------*/

 /// constructor of NCoCubeBlock, taking a pointer to the father Block

 explicit NCoCubeBlock( Block *father = nullptr ) : Block( father ) {}

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// destructor of NCoCubeBlock

 virtual ~NCoCubeBlock() {
  Constraint::clear( v_cuts );   // clear the cuts
  f_obj.clear();                 // clear the Objective
  v_x.clear();                   // delete Variable
  // explicitly reset all Constraint and Variable
  reset_static_variables();
  reset_dynamic_constraints();
  reset_objective();
  }

/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
 /// load the data of the problem, i.e., the cost vector, moving it
 /** The only data of the problem is the cost vector; this also provides the
  * size of the space (number of variables) by its size(). */

 void load( std::vector< FunctionValue > && c ) {
  v_cost = std::move( c );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// load the data of the problem, i.e., the cost vector, copying it

 void load( const std::vector< FunctionValue > & c ) {
  load( std::vector< FunctionValue >( c ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// placeholder for load( std::istream & input  )

 void load( std::istream & input , char frmt = 0 ) override {}
 
/*--------------------------------------------------------------------------*/
 /// extends Block::deserialize( netCDF::NcGroup )

 void deserialize( const netCDF::NcGroup & group ) override {
  throw( std::logic_error( "NCoCubeBlock::deserialize not implemented yet"
			   ) );
  }

/*--------------------------------------------------------------------------*/
 /// generate the abstract variables of the NCoCubeBlock
 /** Method that generates the abstract Variable of the NCoCubeBlock, i.e.,
  * a (static) std::vector of ColVariable with the same size as the costs,
  * bounded to be integers between -1 and 1.
  * It is assumed that load() has been called already. Also, there should be
  * checks that this does not happen twice, but since we are using this
  * class in a very controlled way we forego them. */

 void generate_abstract_variables( Configuration *stvv = nullptr ) override {
  v_x.resize( v_cost.size() );
  for( auto & var : v_x ) {
   var.is_integer( true , eNoBlck );
   var.is_unitary( true , eNoBlck );
   }

  add_static_variable( v_x , "x" );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// generate the static constraint of the NCoCubeBlock
 /** Method that generates the abstract Constraint of the NCoCubeBlock, i.e.,
  * a (dynamic) std::list of FRowConstraint, initially empty. Hence, this
  * method may be called even if load() and generate_abstract_variables()
  * have not yet. Also, there should be checks that this does not happen
  * twice, but since we are using this class in a very controlled way we
  * forego them. */

 void generate_abstract_constraints( Configuration * stcc = nullptr ) override
 {
  v_cuts.clear();  // should not be necessary
  add_dynamic_constraint( v_cuts , "cuts" );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// generate the objective of the NCoCubeBlock
 /** Method that generates the objective of the NCoCubeBlock. It is assumed
  * that load() and generate_abstract_variables() have been called already.
  * Also, there should be checks that this does not happen twice, but since
  * we are using this class in a very controlled way we forego them. */

 void generate_objective( Configuration *objc = nullptr ) override {
  v_coeff_pair p( v_cost.size() );
  auto cit = v_cost.begin();
  auto vit = v_x.begin();
  for( auto & el : p )
   el = std::make_pair( & (*(vit++)) , *(cit++) );

  f_obj.set_sense( Objective::eMin );
  f_obj.set_function( new LinearFunction( std::move( p ) ) , eNoMod );
  set_objective( & f_obj , eNoMod );
  }

/*--------------------------------------------------------------------------*/
 /// separate the lazy constraints
 /** The lazy constraints have the form
  *
  *   \pm 1 x_1 + \pm 1 x_2 + ... + \pm 1 x_n <= 1
  *
  * To see if any is violated we maximise the LSH, i.e., if x_i > 0 then we
  * take it (and we put +1 in the coefficient), otherwise we take - x_i
  * (and we put -1 in the coefficient). If the corresponding sum is > 1
  * then the constraint is violated and it is added, otherwise none of the
  * lazy constraints are violated. */

 void generate_dynamic_constraints( Configuration * dycc = nullptr ) override
 {
  v_coeff_pair cut( v_cost.size() );
  double viol = -1;
  for( Index i = 0 ; i < v_cost.size() ; ++i ) {
   cut[ i ].first = & v_x[ i ];
   auto vxi = v_x[ i ].get_value();
   if( vxi >= 0 ) {
    viol += vxi;
    cut[ i ].second = 1;
    }
   else {
    viol -= vxi;
    cut[ i ].second = -1;
    }
   }

  if( viol >= 1e-8 ) {  // the 1e-8 should not be necessary, v_x is integer
   std::list< FRowConstraint > newcut( 1 );
   newcut.front().set_function( new LinearFunction( std::move( cut ) ) ,
				eNoMod );
   newcut.front().set_rhs( 1 );
   newcut.front().set_lhs( -Inf< FunctionValue >() );
   add_dynamic_constraints( v_cuts , newcut );
   }
  }

/*-------------- METHODS FOR READING TE DATA OF NCoCubeBlock ---------------*/

 const std::vector< FunctionValue > & get_costs( void ) { return( v_cost ); }
 
/*------------- METHODS FOR ADDING / REMOVING / CHANGING DATA --------------*/

 void chg_costs( std::vector< FunctionValue > && c ) {
  if( c.size() != v_cost.size() )
   throw( std::invalid_argument( "NCoCubeBlock::chg_cost: wrong size" ) );
  load( std::move( c ) );
  auto lf = static_cast< p_LF >( f_obj.get_function() );
  lf->modify_coefficients( std::vector< FunctionValue >( v_cost ) );
  }

/*-------------------- PROTECTED PART OF THE CLASS -------------------------*/

 protected:

/*--------------------------- PROTECTED FIELDS  ----------------------------*/

 std::vector< FunctionValue > v_cost;
 
 std::vector< ColVariable > v_x;

 std::list< FRowConstraint > v_cuts;

 FRealObjective f_obj;
 
/*--------------------------------------------------------------------------*/

 SMSpp_insert_in_factory_h;  // insert NCoCubeBlock in the Block factory

/*--------------------------------------------------------------------------*/

 };  // end( class NCoCubeBlock )

/*--------------------------------------------------------------------------*/
/*-------------------------- FACTORY MANAGEMENT ----------------------------*/
/*--------------------------------------------------------------------------*/

// register NCoCubeBlock to the Block factory

SMSpp_insert_in_factory_cpp_1( NCoCubeBlock );

/*--------------------------------------------------------------------------*/
/*------------------------------- CONSTANTS --------------------------------*/
/*--------------------------------------------------------------------------*/

static constexpr double scale = 10;

c_FunctionValue INF = SMSpp_di_unipi_it::Inf< FunctionValue >();

/*--------------------------------------------------------------------------*/
/*------------------------------- GLOBALS ----------------------------------*/
/*--------------------------------------------------------------------------*/

std::mt19937 rg;           // base random generator
std::uniform_real_distribution<> dis( - 1.0 , 1.0 );

Index nvar = 10;           // number of variables

NCoCubeBlock NCCB;         // the NCoCubeBlock

/*--------------------------------------------------------------------------*/
/*------------------------------ FUNCTIONS ---------------------------------*/
/*--------------------------------------------------------------------------*/

template< class T >
static void Str2Sthg( const char * const str , T & sthg )
{
 istringstream( str ) >> sthg;
 }

/*--------------------------------------------------------------------------*/

std::vector< FunctionValue > generate_costs( void )
{
 std::vector< FunctionValue > retval( nvar );
 for( auto & el : retval )
  el = scale * dis( rg );

 return( retval );
 }

/*--------------------------------------------------------------------------*/

static bool SolveBoth( void ) 
{
 try {
  // solve the NCoCubeBlock - - - - - - - - - - - - - - - - - - - - - - - - -
  Solver * slvr = NCCB.get_registered_solvers().front();
  #if DETACH_LP
   NCCB.unregister_Solver( slvr );
   NCCB.register_Solver( slvr , true );  // push to front
  #endif
  int rtrn = slvr->compute( false );
  bool hs = ( ( rtrn >= Solver::kOK ) && ( rtrn < Solver::kError ) )
            || ( rtrn == Solver::kLowPrecision );

  if( ! hs ) {
   if( rtrn == Solver::kInfeasible )
    cout << "Unfeas(?)" << endl;
   else
    if( rtrn == Solver::kUnbounded )
     cout << "Unbounded(?)" << endl;
    else
     cout << "Error!" << endl;
   return( false );
   }

  double fo = hs ? slvr->get_ub() : INF;  // get Solver optimal value

  // manually compute the optimal value
  double opt = INF;
  for( auto ci : NCCB.get_costs() )
   if( - std::abs( ci ) < opt )
    opt = - std::abs( ci );

  if( abs( fo - opt ) <= 1e-7 * max( double( 1 ) , abs( opt ) ) ) {
   LOG1( "OK" << endl );
   return( true );
   }

  LOG1( "Error: fo = " << fo << ", opt = " << opt << endl );
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
 }  // end( SolveBoth )

/*--------------------------------------------------------------------------*/

int main( int argc , char **argv )
{
 // reading command line parameters - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 long int seed = 0;
 Index n_repeat = 10;

 switch( argc ) {
  case( 4 ): Str2Sthg( argv[ 3 ] , n_repeat );
  case( 3 ): Str2Sthg( argv[ 2 ] , nvar );
  case( 2 ): Str2Sthg( argv[ 1 ] , seed );
             break;
  default: cerr << "Usage: " << argv[ 0 ] << " seed [nvar n_repeat]"
	        << endl << "       nvar: number of variables [10]"
                << endl << "       n_repeat: number of repetitions [10]"
	        << endl;
	   return( 1 );
  }

 if( nvar < 1 ) {
  cout << "error: nvar too small";
  exit( 1 );
  }

 rg.seed( seed );  // seed the pseudo-random number generator

 cout.setf( ios::scientific, ios::floatfield );
 cout << setprecision( 10 );

 // load the problem- - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 NCCB.load( generate_costs() );

 // do NOT generate the abstract representation, since this is done by
 // CPXMILPSolver when it is registered, and since there is no check in
 // NCoCubeBlock it would end up being generated twice
 // NCCB.generate_abstract_variables();
 // NCCB.generate_abstract_constraints();
 // NCCB.generate_objective();

 // attach the Solver to the Block- - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // do it by using a single a BlockSolverConfig, read from file
 
 auto bsc = dynamic_cast< BlockSolverConfig * >(
		               Configuration::deserialize( "MILPPar.txt" ) );
 if( ! bsc ) {
  cerr << "Error: configuration file not a BlockSolverConfig" << endl;
  exit( 1 );    
  }

 bsc->apply( & NCCB );
 bsc->clear();  // keep the clear()-ed BlockSolverConfig for final cleanup

 // check Solvers - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( NCCB.get_registered_solvers().empty() ) {
  cerr << "Error: BlockSolverConfig did not register any Solver" << endl;
  exit( 1 );    
  }

 // open log-file - - - - - - - - - - -  - - - - - - - - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 #if( LOG_LEVEL >= 2 )
  #if( LOG_ON_COUT )
   NCCB.get_registered_solvers().front()->set_log( & cout );
  #else
   ofstream LOGFile( "log.txt" , ofstream::out );
   if( ! LOGFile.is_open() )
    cerr << "Warning: cannot open log file" << endl;
   else {
    LOGFile.setf( ios::scientific , ios::floatfield );
    LOGFile << setprecision( 10 );
    NCCB.get_registered_solvers().front()->set_log( & LOGFile );
    }
  #endif
 #endif

 // open log-file - - - - - - - - - - -  - - - - - - - - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 #if( LOG_LEVEL >= 3 )
   NCCB.get_registered_solvers().front()->set_par(
		                     MILPSolver::strOutputFile , "NCCBlock.lp" );
 #endif

 // first solver call - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 LOG1( "0: " );

 bool AllPassed = SolveBoth();
 
 // main loop - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // now, for n_repeat times, all costs are changed then the problem is
// re-solved

 for( Index rep = 0 ; rep < n_repeat ; ++rep ) {
  LOG1( rep << ": ");

  NCCB.chg_costs( generate_costs() );

  #if( LOG_LEVEL >= 3 )
   NCCB.get_registered_solvers().front()->set_par(
		                     MILPSolver::strOutputFile , "NCCBlock-" +
                         std::to_string( rep ) + ".lp" );
  #endif

  AllPassed &= SolveBoth();

  }  // end( main loop )- - - - - - - - - - - - - - - - - - - - - - - - - - -
     // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( AllPassed )
  cout << GREEN( All test passed!! ) << endl;
 else
  cout << RED( Shit happened!! ) << endl;
 
 // final cleanup - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

 bsc->apply( & NCCB );  // remove the Solver by apply()-ing the cleared bsc

 delete bsc;            // delete the BlockSolverConfig

 // all done- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

 return( AllPassed ? 0 : 1 );

 }  // end( main )

/*--------------------------------------------------------------------------*/
/*------------------------ End File test.cpp -------------------------------*/
/*--------------------------------------------------------------------------*/
