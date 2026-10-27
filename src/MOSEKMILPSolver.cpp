/*--------------------------------------------------------------------------*/
/*------------------------ File MOSEKMILPSolver.cpp ------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Implementation of the MOSEKMILPSolver class, which implements the
 * MILPSolver interface for Mixed-Integer Linear and convex Quadratic
 * Programs using the Optimizer API for C of MOSEK.
 *
 * \author Donato Meoli \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \copyright &copy; by Donato Meoli
 */
/*--------------------------------------------------------------------------*/
/*---------------------------- IMPLEMENTATION ------------------------------*/
/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include <algorithm>
#include <cmath>
#include <iostream>
#include <map>
#include <numeric>
#include <stdexcept>
#include <vector>

#include <LinearFunction.h>

#include <QuadFunction.h>

#include "MOSEKMILPSolver.h"

#ifdef MILPSolver_DEBUG
 #define DEBUG_LOG( stuff ) std::cout << "[MILPSolver DEBUG] " << stuff
#else
 #define DEBUG_LOG( stuff )
#endif

#include <boost/preprocessor/cat.hpp>
#include <boost/preprocessor/stringize.hpp>
#include BOOST_PP_STRINGIZE( BOOST_PP_CAT( BOOST_PP_CAT( BOOST_PP_CAT( BOOST_PP_CAT( MOSEK , MSK_VERSION_MAJOR ) , MSK_VERSION_MINOR ) , MSK_VERSION_REVISION ) , _maps.h ) )

/*--------------------------------------------------------------------------*/
/*------------------------- NAMESPACE AND USING ----------------------------*/
/*--------------------------------------------------------------------------*/

using namespace SMSpp_di_unipi_it;

/*--------------------------------------------------------------------------*/
/*-------------------------- FACTORY MANAGEMENT ----------------------------*/
/*--------------------------------------------------------------------------*/

SMSpp_insert_in_factory_cpp_0( MOSEKMILPSolver );

/*--------------------------------------------------------------------------*/
/*----------------------------- FUNCTIONS ----------------------------------*/
/*--------------------------------------------------------------------------*/

namespace {

// the callback of MOSEK, which defers to the class method

MSKint32t MSKAPI MOSEKMILPSolver_callback( MSKtask_t , MSKuserhandle_t h ,
					   MSKcallbackcodee caller ,
					   const MSKrealt * douinf ,
					   const MSKint32t * intinf ,
					   const MSKint64t * )
{
 return( static_cast< MOSEKMILPSolver * >( h )->callback( caller , douinf ,
							  intinf ) );
 }

/*--------------------------------------------------------------------------*/

// the stream of MOSEK printing the log on std::cout

void MSKAPI MOSEKMILPSolver_print( MSKuserhandle_t , const char * str )
{
 std::cout << str << std::flush;
 }

/*--------------------------------------------------------------------------*/

// the text of a response code of MOSEK

std::string code_text( MSKrescodee r )
{
 std::array< char , MSK_MAX_STR_LEN > sym = {} , str = {};
 if( MSK_getcodedesc( r , sym.data() , str.data() ) != MSK_RES_OK )
  return( "MOSEK response code " + std::to_string( int( r ) ) );
 return( std::string( sym.data() ) + " (" + str.data() + ")" );
 }

/*--------------------------------------------------------------------------*/

// true if r says that the license could not be checked out

bool is_license_error( MSKrescodee r )
{
 return( ( r >= MSK_RES_ERR_LICENSE ) &&
	 ( r <= MSK_RES_ERR_LICENSE_NO_SERVER_LINE ) );
 }

/*--------------------------------------------------------------------------*/

// true if MOSEK reads the side v as an infinite one

bool is_inf( double v ) { return( std::abs( v ) >= MSK_INFINITY ); }

/*--------------------------------------------------------------------------*/

// a side as MOSEK is given it: an infinite one, which the bound key ignores
// anyway, is +- MSK_INFINITY rather than the infinity of the double

double msk_side( double v )
{
 return( std::clamp( v , - MSK_INFINITY , MSK_INFINITY ) );
 }

}  // end( anonymous namespace )

/*--------------------------------------------------------------------------*/
/*------------------------- THE SHARED ENVIRONMENT -------------------------*/
/*--------------------------------------------------------------------------*/

unsigned int MOSEKMILPSolver::f_env_users = 0;
MSKenv_t MOSEKMILPSolver::f_shared_env = nullptr;
std::mutex MOSEKMILPSolver::f_env_mutex;

/*--------------------------------------------------------------------------*/
/*--------------------- CONSTRUCTOR AND DESTRUCTOR -------------------------*/
/*--------------------------------------------------------------------------*/

MOSEKMILPSolver::MOSEKMILPSolver( void ) :
 MILPSolver() , env( nullptr ) , task( nullptr ) ,
 f_cb_bound( - Inf< double >() ) , f_cb_incumbent( Inf< double >() ) ,
 UpCutOff( Inf< double >() ) , LwCutOff( - Inf< double >() )
{
 // the environment is one for the whole process: creating it does not
 // check out the license, which the first solve does, and the license then
 // stays with the environment rather than being checked out by each task
 {
  std::lock_guard< std::mutex > lock( f_env_mutex );

  if( ! f_shared_env ) {
   auto r = MSK_makeenv( & f_shared_env , nullptr );
   if( r != MSK_RES_OK ) {
    f_shared_env = nullptr;
    throw( std::runtime_error( "MOSEKMILPSolver::MOSEKMILPSolver: "
			       "MSK_makeenv failed with " + code_text( r ) ) );
    }
   }

  env = f_shared_env;
  ++f_env_users;
  }

 auto r = MSK_maketask( env , 0 , 0 , & task );
 if( r != MSK_RES_OK ) {
  std::lock_guard< std::mutex > lock( f_env_mutex );
  if( ( --f_env_users == 0 ) && f_shared_env ) {
   MSK_deleteenv( & f_shared_env );
   f_shared_env = nullptr;
   }
  throw( std::runtime_error( "MOSEKMILPSolver::MOSEKMILPSolver: "
			     "MSK_maketask failed with " + code_text( r ) ) );
  }

 // the log is off until intLogVerb says otherwise
 MSK_putintparam( task , MSK_IPAR_LOG , 0 );
 }

/*--------------------------------------------------------------------------*/

MOSEKMILPSolver::~MOSEKMILPSolver()
{
 if( task )
  MSK_deletetask( & task );

 // the environment goes when the last of us goes
 std::lock_guard< std::mutex > lock( f_env_mutex );
 if( ( f_env_users > 0 ) && ( --f_env_users == 0 ) && f_shared_env ) {
  MSK_deleteenv( & f_shared_env );
  f_shared_env = nullptr;
  }
 env = nullptr;
 }

/*--------------------------------------------------------------------------*/
/*--------------------- DERIVED METHODS OF BASE CLASS ----------------------*/
/*--------------------------------------------------------------------------*/

void MOSEKMILPSolver::clear_problem( unsigned int what )
{
 MILPSolver::clear_problem( 0 );

 // the task is emptied rather than deleted, so that it keeps the parameters
 // it has been given
 MSKint32t n = 0;
 MSK_getnumvar( task , & n );
 if( n > 0 ) {
  std::vector< MSKint32t > idx( n );
  std::iota( idx.begin() , idx.end() , 0 );
  check( MSK_removevars( task , n , idx.data() ) , "clear_problem" );
  }

 MSK_getnumcon( task , & n );
 if( n > 0 ) {
  std::vector< MSKint32t > idx( n );
  std::iota( idx.begin() , idx.end() , 0 );
  check( MSK_removecons( task , n , idx.data() ) , "clear_problem" );
  }

 check( MSK_putqobj( task , 0 , nullptr , nullptr , nullptr ) ,
	"clear_problem" );

 q_hessian.clear();
 q_hessian_nnz = 0;
 f_hessian_changed = false;
 f_soltype = MSK_SOL_ITR;
 }

/*--------------------------------------------------------------------------*/

void MOSEKMILPSolver::load_problem( void )
{
 MILPSolver::load_problem();

 // the problems MOSEKMILPSolver does not take are refused before anything
 // reaches the task: with quadratic rows MILPSolver stores the matrix by
 // rows, which would be misread below
 if( numquadrows > 0 )
  throw( std::invalid_argument( "MOSEKMILPSolver::load_problem: "
				"quadratic constraints are not supported" ) );

 for( int j = 0 ; j < numcols ; ++j )
  if( ( xctype[ j ] == 'S' ) || ( xctype[ j ] == 'N' ) )
   throw( std::invalid_argument(
     "MOSEKMILPSolver::load_problem: column " + std::to_string( j ) +
     " is semi-continuous or semi-integer, which MOSEK does not support" ) );

 // the task is emptied first, its parameters staying as they are
 {
  MSKint32t n = 0;
  MSK_getnumvar( task , & n );
  if( n > 0 ) {
   std::vector< MSKint32t > idx( n );
   std::iota( idx.begin() , idx.end() , 0 );
   check( MSK_removevars( task , n , idx.data() ) , "load_problem" );
   }
  MSK_getnumcon( task , & n );
  if( n > 0 ) {
   std::vector< MSKint32t > idx( n );
   std::iota( idx.begin() , idx.end() , 0 );
   check( MSK_removecons( task , n , idx.data() ) , "load_problem" );
   }
  }

 check( MSK_putobjsense( task , objsense == -1 ? MSK_OBJECTIVE_SENSE_MAXIMIZE
				               : MSK_OBJECTIVE_SENSE_MINIMIZE ) ,
	"load_problem" );

 if( numrows > 0 )
  check( MSK_appendcons( task , numrows ) , "load_problem" );
 if( numcols > 0 )
  check( MSK_appendvars( task , numcols ) , "load_problem" );

 if( numcols > 0 ) {
  check( MSK_putcslice( task , 0 , numcols , objective.data() ) ,
	 "load_problem" );

  // the bounds of the columns, each with its key
  std::vector< MSKboundkeye > bk( numcols );
  std::vector< double > bl( numcols ) , bu( numcols );
  for( int j = 0 ; j < numcols ; ++j ) {
   bk[ j ] = bound_key( lb[ j ] , ub[ j ] );
   bl[ j ] = msk_side( lb[ j ] );
   bu[ j ] = msk_side( ub[ j ] );
   }
  check( MSK_putvarboundslice( task , 0 , numcols , bk.data() , bl.data() ,
			       bu.data() ) , "load_problem" );

  // the matrix, by columns, as MILPSolver keeps it
  std::vector< MSKint32t > ptre( numcols );
  for( int j = 0 ; j < numcols ; ++j )
   ptre[ j ] = matbeg[ j ] + matcnt[ j ];
  if( ! matval.empty() )
   check( MSK_putacolslice( task , 0 , numcols , matbeg.data() ,
			    ptre.data() , matind.data() , matval.data() ) ,
	  "load_problem" );

  // the type of the columns
  std::vector< MSKint32t > cols( numcols );
  std::iota( cols.begin() , cols.end() , 0 );
  std::vector< MSKvariabletypee > vt( numcols );
  for( int j = 0 ; j < numcols ; ++j )
   vt[ j ] = ( ( xctype[ j ] == 'B' ) || ( xctype[ j ] == 'I' ) ) ?
             MSK_VAR_TYPE_INT : MSK_VAR_TYPE_CONT;
  check( MSK_putvartypelist( task , numcols , cols.data() , vt.data() ) ,
	 "load_problem" );
  }

 // the sides of the rows, out of the sense, the right-hand side and the
 // range of MILPSolver
 if( numrows > 0 ) {
  std::vector< MSKboundkeye > bk( numrows );
  std::vector< double > bl( numrows ) , bu( numrows );
  for( int i = 0 ; i < numrows ; ++i ) {
   switch( sense[ i ] ) {
    case( 'L' ): bl[ i ] = - Inf< double >(); bu[ i ] = rhs[ i ]; break;
    case( 'G' ): bl[ i ] = rhs[ i ]; bu[ i ] = Inf< double >(); break;
    case( 'E' ): bl[ i ] = bu[ i ] = rhs[ i ]; break;
    case( 'R' ):
     if( rngval[ i ] > 0 ) {
      bl[ i ] = rhs[ i ];
      bu[ i ] = rhs[ i ] + rngval[ i ];
      }
     else {
      bl[ i ] = rhs[ i ] + rngval[ i ];
      bu[ i ] = rhs[ i ];
      }
     break;
    default:
     throw( std::invalid_argument( "MOSEKMILPSolver::load_problem: row " +
				   std::to_string( i ) + " has sense " +
				   std::string( 1 , sense[ i ] ) ) );
    }
   bk[ i ] = bound_key( bl[ i ] , bu[ i ] );
   bl[ i ] = msk_side( bl[ i ] );
   bu[ i ] = msk_side( bu[ i ] );
   }
  check( MSK_putconboundslice( task , 0 , numrows , bk.data() , bl.data() ,
			       bu.data() ) , "load_problem" );
  }

 // the quadratic part of the objective
 generate_qobj_hessian();
 pass_hessian();

 if( use_custom_names ) {
  for( int j = 0 ; j < numcols ; ++j )
   if( colname[ j ] )
    MSK_putvarname( task , j , colname[ j ] );
  for( int i = 0 ; i < numrows ; ++i )
   if( rowname[ i ] )
    MSK_putconname( task , i , rowname[ i ] );
  }

 f_soltype = MSK_SOL_ITR;

 // the base representation isn't needed anymore
 MILPSolver::clear_problem( 15 );

 }  // end( MOSEKMILPSolver::load_problem )

/*--------------------------------------------------------------------------*/
/*-------------------- METHODS FOR SOLVING THE MODEL -----------------------*/
/*--------------------------------------------------------------------------*/

int MOSEKMILPSolver::guts_of_compute( void )
{
 // Note: locking, process_modifications() and the LP cut separation loop
 // (when intRelaxIntVars == 2) are all handled by MILPSolver::compute().
 // This method is only responsible for the actual MOSEK call.

 // the changes to the quadratic objective go to MOSEK in one go
 pass_hessian();

 if( ! output_file.empty() )
  write_lp( output_file );

 // the MIP path is taken only if the model has integer columns, i.e., with
 // intRelaxIntVars == 0: with 1 or 2 the integer Variable are given to MOSEK
 // as continuous ones and int_vars is 0
 const bool mip = ( int_vars > 0 );

 if( mip && ( CutSepPar & 7 ) )
  throw( std::invalid_argument(
     "MOSEKMILPSolver::guts_of_compute: MOSEK has no callback for user cuts "
     "or lazy constraints, so intCutSepPar cannot ask for separation in a "
     "MIP solve (intRelaxIntVars == 2 separates on the LP relaxation)" ) );

 // the callback is only there for the cutoffs of a MIP
 const bool cutoffs = mip && ( ( UpCutOff < Inf< double >() ) ||
			       ( LwCutOff > - Inf< double >() ) );
 if( cutoffs ) {
  f_cb_bound = - Inf< double >();
  f_cb_incumbent = Inf< double >();
  if( objsense == -1 )
   std::swap( f_cb_bound , f_cb_incumbent );
  check( MSK_putcallbackfunc( task , & MOSEKMILPSolver_callback , this ) ,
	 "guts_of_compute" );
  f_callback_set = true;
  }
 else
  if( f_callback_set ) {
   MSK_putcallbackfunc( task , nullptr , nullptr );
   f_callback_set = false;
   }

 MSKrescodee trm = MSK_RES_OK;
 const MSKrescodee r = MSK_optimizetrm( task , & trm );

 // the solution to read: the integer one of a MIP, the basic one if the
 // solve has produced one, which the simplex algorithms and the basis
 // identification after the interior point do, the interior point one
 // otherwise
 MSKbooleant def = 0;
 if( mip )
  f_soltype = MSK_SOL_ITG;
 else
  if( ( MSK_solutiondef( task , MSK_SOL_BAS , & def ) == MSK_RES_OK ) &&
      def ) {
   MSKsolstae sta = MSK_SOL_STA_UNKNOWN;
   MSK_getsolsta( task , MSK_SOL_BAS , & sta );
   MSKbooleant itr = 0;
   MSK_solutiondef( task , MSK_SOL_ITR , & itr );
   f_soltype = ( ( sta != MSK_SOL_STA_UNKNOWN ) || ( ! itr ) ) ?
               MSK_SOL_BAS : MSK_SOL_ITR;
   }
  else
   f_soltype = MSK_SOL_ITR;

 sol_status = decode_status( r , trm );
 return( sol_status );

 }  // end( MOSEKMILPSolver::guts_of_compute )

/*--------------------------------------------------------------------------*/

// what a status that is a failure of the run says, if the log is on, before
// it becomes kError: the caller reads a failure and not an answer about the
// model, and kError on its own tells neither which failure it was nor that
// the model is fine

static int failed( const std::string & why , bool say )
{
 if( say )
  std::cerr << "MOSEKMILPSolver::guts_of_compute: warning: MOSEK stopped on "
	    << why << ", the model was not solved" << std::endl;
 return( Solver::kError );
 }

/*--------------------------------------------------------------------------*/

int MOSEKMILPSolver::decode_status( MSKrescodee r , MSKrescodee trm )
{
 DEBUG_LOG( "MSK_optimizetrm returned " << r << ", termination " << trm
	    << std::endl );

 // an error of the call, rather than of the model: the license, which has
 // nothing to do with the model, is said out loud and is no exception, a
 // computation that may have been running for hours being not worth taking
 // down for it; a non-convex objective is a model MOSEK cannot solve; the
 // others are errors of this class or of MOSEK
 if( r != MSK_RES_OK ) {
  if( is_license_error( r ) ) {
   std::cerr << "MOSEKMILPSolver::guts_of_compute: warning: the license "
                "could not be checked out [" << code_text( r )
	     << "], the model was not solved" << std::endl;
   return( kError );
   }

  if( ( r == MSK_RES_ERR_OBJ_Q_NOT_PSD ) ||
      ( r == MSK_RES_ERR_OBJ_Q_NOT_NSD ) )
   return( failed( "an objective that is not convex [" + code_text( r ) +
		   "]" , log_verbosity > 0 ) );

  // warnings are below 1000 and are no failure
  if( r >= 1000 )
   throw( std::runtime_error( "MOSEKMILPSolver::guts_of_compute: "
			      "MSK_optimizetrm failed with " +
			      code_text( r ) ) );
  }

 MSKprostae prosta = MSK_PRO_STA_UNKNOWN;
 MSKsolstae solsta = MSK_SOL_STA_UNKNOWN;
 MSKbooleant def = 0;
 if( ( MSK_solutiondef( task , f_soltype , & def ) == MSK_RES_OK ) && def ) {
  MSK_getprosta( task , f_soltype , & prosta );
  MSK_getsolsta( task , f_soltype , & solsta );
  }

 switch( trm ) {
  case( MSK_RES_OK ): break;
  case( MSK_RES_TRM_MAX_TIME ): return( kStopTime );
  case( MSK_RES_TRM_MAX_ITERATIONS ):
  case( MSK_RES_TRM_MIO_NUM_RELAXS ):
  case( MSK_RES_TRM_MIO_NUM_BRANCHES ): return( kStopIter );
  case( MSK_RES_TRM_NUM_MAX_NUM_INT_SOLUTIONS ):
  case( MSK_RES_TRM_OBJECTIVE_RANGE ): return( kOK );
  case( MSK_RES_TRM_USER_CALLBACK ):
   // stopped by the callback at a cutoff: MILPSolver::compute() turns
   // this into kCutOff
   if( f_cutoff_stop )
    return( kError );
   return( failed( "an interruption of the callback" , log_verbosity > 0 ) );
  case( MSK_RES_TRM_STALL ):
  case( MSK_RES_TRM_NUMERICAL_PROBLEM ):
  case( MSK_RES_TRM_MAX_NUM_SETBACKS ):
   // the interior point stops on slow progress also when what it has is
   // already an answer, typically the certificate of an infeasible or of
   // an unbounded problem, which MOSEK then gives with its status: that is
   // used as it is, while a solution with no accuracy promise attached is
   // what kLowPrecision says
   switch( solsta ) {
    case( MSK_SOL_STA_OPTIMAL ):
    case( MSK_SOL_STA_INTEGER_OPTIMAL ): return( kOK );
    case( MSK_SOL_STA_PRIM_INFEAS_CER ): return( kInfeasible );
    case( MSK_SOL_STA_DUAL_INFEAS_CER ): return( kUnbounded );
    default: break;
    }
   if( has_primal() )
    return( kLowPrecision );
   return( failed( code_text( trm ) , log_verbosity > 0 ) );
  default:
   return( failed( code_text( trm ) , log_verbosity > 0 ) );
  }

 switch( solsta ) {
  case( MSK_SOL_STA_OPTIMAL ):
  case( MSK_SOL_STA_INTEGER_OPTIMAL ): return( kOK );
  case( MSK_SOL_STA_PRIM_INFEAS_CER ): return( kInfeasible );
  case( MSK_SOL_STA_DUAL_INFEAS_CER ): return( kUnbounded );
  default: break;
  }

 switch( prosta ) {
  case( MSK_PRO_STA_PRIM_INFEAS ):
  case( MSK_PRO_STA_PRIM_AND_DUAL_INFEAS ): return( kInfeasible );
  case( MSK_PRO_STA_DUAL_INFEAS ):
  // as GRB_INF_OR_UNBD for GRBMILPSolver, which MOSEK only gives on a MIP
  case( MSK_PRO_STA_PRIM_INFEAS_OR_UNBOUNDED ): return( kUnbounded );
  default: break;
  }

 // a MIP whose integer solution has not been found, nor the problem
 // proved anything, has nothing to give
 if( has_primal() )
  return( kLowPrecision );

 return( failed( "solution status " + std::to_string( int( solsta ) ) +
		 ", problem status " + std::to_string( int( prosta ) ) ,
		 log_verbosity > 0 ) );

 }  // end( MOSEKMILPSolver::decode_status )

/*--------------------------------------------------------------------------*/

int MOSEKMILPSolver::callback( MSKcallbackcodee caller ,
			       const MSKrealt * douinf ,
			       const MSKint32t * intinf )
{
 if( ( int_vars == 0 ) || ( ! douinf ) || ( ! intinf ) )
  return( 0 );

 // the incumbent and the bound are in the sense of the problem, but lack
 // the constant of the objective, which the cutoffs have; the bound only
 // means something after the first relaxation, the incumbent after the
 // first integer solution
 const bool min = ( objsense != -1 );
 double solv = min ? Inf< double >() : - Inf< double >();
 double bndv = - solv;

 if( intinf[ MSK_IINF_MIO_NUM_INT_SOLUTIONS ] > 0 ) {
  solv = douinf[ MSK_DINF_MIO_OBJ_INT ];
  if( is_inf( solv ) )
   solv = min ? Inf< double >() : - Inf< double >();
  else
   solv += constant_value;
  }

 if( intinf[ MSK_IINF_MIO_NUM_RELAX ] > 0 ) {
  bndv = douinf[ MSK_DINF_MIO_OBJ_BOUND ];
  if( is_inf( bndv ) )
   bndv = min ? - Inf< double >() : Inf< double >();
  else
   bndv += constant_value;
  }

 f_cb_bound = bndv;
 f_cb_incumbent = solv;

 if( min ) {  // solv is an upper bound and bndv a lower bound
  if( ( bndv >= up_cut_off() ) || ( solv <= lw_cut_off() ) ) {
   f_cutoff_stop = true;
   return( 1 );
   }
  }
 else         // solv is a lower bound and bndv an upper bound
  if( ( solv >= up_cut_off() ) || ( bndv <= lw_cut_off() ) ) {
   f_cutoff_stop = true;
   return( 1 );
   }

 return( 0 );

 }  // end( MOSEKMILPSolver::callback )

/*--------------------------------------------------------------------------*/
/*---------------------- METHODS FOR READING RESULTS -----------------------*/
/*--------------------------------------------------------------------------*/

MSKsolstae MOSEKMILPSolver::solution_status( void ) const
{
 MSKbooleant def = 0;
 MSKsolstae sta = MSK_SOL_STA_UNKNOWN;
 if( ( MSK_solutiondef( task , f_soltype , & def ) == MSK_RES_OK ) && def )
  MSK_getsolsta( task , f_soltype , & sta );
 return( sta );
 }

/*--------------------------------------------------------------------------*/

bool MOSEKMILPSolver::has_primal( void ) const
{
 switch( solution_status() ) {
  case( MSK_SOL_STA_OPTIMAL ):
  case( MSK_SOL_STA_INTEGER_OPTIMAL ):
  case( MSK_SOL_STA_PRIM_FEAS ):
  case( MSK_SOL_STA_PRIM_AND_DUAL_FEAS ): return( true );
  default: return( false );
  }
 }

/*--------------------------------------------------------------------------*/

Solver::OFValue MOSEKMILPSolver::primal_value( void ) const
{
 double v = 0;
 MSK_getprimalobj( task , f_soltype , & v );
 return( v + constant_value );
 }

/*--------------------------------------------------------------------------*/

Solver::OFValue MOSEKMILPSolver::stopped_dual_bound( OFValue none ) const
{
 if( int_vars == 0 )
  return( none );

 MSKint32t nrelax = 0;
 MSK_getintinf( task , MSK_IINF_MIO_NUM_RELAX , & nrelax );
 if( nrelax <= 0 )
  return( none );

 double bound = 0;
 if( ( MSK_getdouinf( task , MSK_DINF_MIO_OBJ_BOUND , & bound ) !=
       MSK_RES_OK ) || is_inf( bound ) )
  return( none );
 return( bound + constant_value );
 }

/*--------------------------------------------------------------------------*/

Solver::OFValue MOSEKMILPSolver::get_lb( void )
{
 if( objsense != -1 ) {  // minimization problem- - - - - - - - - - - - - -
  switch( sol_status ) {
   case( kUnbounded ):  return( - Inf< OFValue >() );
   case( kInfeasible ): return( Inf< OFValue >() );
   case( kOK ):
    // the bound of the search for a MIP, which the gaps leave below the
    // incumbent, the optimal value otherwise
    if( int_vars > 0 )
     return( stopped_dual_bound( primal_value() ) );
    return( primal_value() );
   // a stopped solve, or one ended without accuracy, has no lower bound
   // but the one of the search of a MIP
   case( kLowPrecision ):
   case( kStopIter ):
   case( kCutOff ):
   case( kStopTime ):   return( stopped_dual_bound( - Inf< OFValue >() ) );
   case( kUnEval ):     return( f_cb_bound );
   default:             return( - Inf< OFValue >() );
   }
  }

 // maximization problem- - - - - - - - - - - - - - - - - - - - - - - - - - -
 switch( sol_status ) {
  case( kUnbounded ):  return( Inf< OFValue >() );
  case( kInfeasible ): return( - Inf< OFValue >() );
  case( kOK ):
  case( kLowPrecision ):
  case( kStopIter ):
  case( kCutOff ):
  case( kStopTime ):
   return( has_primal() ? primal_value() : - Inf< OFValue >() );
  case( kUnEval ):     return( f_cb_incumbent );
  default:             return( - Inf< OFValue >() );
  }
 }

/*--------------------------------------------------------------------------*/

Solver::OFValue MOSEKMILPSolver::get_ub( void )
{
 if( objsense != -1 ) {  // minimization problem- - - - - - - - - - - - - -
  switch( sol_status ) {
   case( kUnbounded ):  return( - Inf< OFValue >() );
   case( kInfeasible ): return( Inf< OFValue >() );
   case( kOK ):
   case( kLowPrecision ):
   case( kStopIter ):
   case( kCutOff ):
   case( kStopTime ):
    return( has_primal() ? primal_value() : Inf< OFValue >() );
   case( kUnEval ):     return( f_cb_incumbent );
   default:             return( Inf< OFValue >() );
   }
  }

 // maximization problem- - - - - - - - - - - - - - - - - - - - - - - - - - -
 switch( sol_status ) {
  case( kUnbounded ):  return( Inf< OFValue >() );
  case( kInfeasible ): return( - Inf< OFValue >() );
  case( kOK ):
   if( int_vars > 0 )
    return( stopped_dual_bound( primal_value() ) );
   return( primal_value() );
  case( kLowPrecision ):
  case( kStopIter ):
  case( kCutOff ):
  case( kStopTime ):   return( stopped_dual_bound( Inf< OFValue >() ) );
  case( kUnEval ):     return( f_cb_bound );
  default:             return( Inf< OFValue >() );
  }
 }

/*--------------------------------------------------------------------------*/

Solver::OFValue MOSEKMILPSolver::get_var_value( void )
{
 return( objsense == -1 ? get_lb() : get_ub() );
 }

/*--------------------------------------------------------------------------*/

bool MOSEKMILPSolver::has_var_solution( void )
{
 return( has_primal() );
 }

/*--------------------------------------------------------------------------*/

void MOSEKMILPSolver::get_var_solution( Configuration * solc )
{
 std::vector< double > x( numcols , 0 );
 if( numcols > 0 )
  check( MSK_getxx( task , f_soltype , x.data() ) , "get_var_solution" );

 MILPSolver::write_var_solution( x , solc );
 }

/*--------------------------------------------------------------------------*/

bool MOSEKMILPSolver::has_var_direction( void )
{
 return( ( int_vars == 0 ) &&
	 ( solution_status() == MSK_SOL_STA_DUAL_INFEAS_CER ) );
 }

/*--------------------------------------------------------------------------*/

void MOSEKMILPSolver::get_var_direction( Configuration * dirc )
{
 if( ! has_var_direction() )
  throw( std::logic_error( "MOSEKMILPSolver::get_var_direction: "
			   "no unbounded direction available" ) );

 // the primal certificate of dual infeasibility is in the primal part of
 // the solution
 std::vector< double > x( numcols , 0 );
 if( numcols > 0 )
  check( MSK_getxx( task , f_soltype , x.data() ) , "get_var_direction" );

 MILPSolver::write_var_solution( x , dirc );
 }

/*--------------------------------------------------------------------------*/

bool MOSEKMILPSolver::has_dual_solution( void )
{
 if( int_vars > 0 )
  return( false );

 switch( solution_status() ) {
  case( MSK_SOL_STA_OPTIMAL ):
  case( MSK_SOL_STA_DUAL_FEAS ):
  case( MSK_SOL_STA_PRIM_AND_DUAL_FEAS ): return( true );
  default: return( false );
  }
 }

/*--------------------------------------------------------------------------*/

bool MOSEKMILPSolver::is_dual_feasible( void )
{
 // MOSEK says whether the dual part of its solution is feasible
 return( has_dual_solution() );
 }

/*--------------------------------------------------------------------------*/

void MOSEKMILPSolver::get_dual_solution( Configuration * solc )
{
 if( int_vars > 0 )
  throw( std::logic_error( "MOSEKMILPSolver::get_dual_solution: "
			   "no dual solution for a MIP" ) );

 std::vector< double > pi( numrows , 0 );
 if( numrows > 0 )
  check( MSK_gety( task , f_soltype , pi.data() ) , "get_dual_solution" );

 // y is slc - suc, and the reduced cost slx - sux: A' y + slx - sux = c
 // (+ Q x) for both senses of the objective
 std::vector< double > rc( numcols , 0 );
 if( numcols > 0 ) {
  std::vector< double > sux( numcols , 0 );
  check( MSK_getslx( task , f_soltype , rc.data() ) , "get_dual_solution" );
  check( MSK_getsux( task , f_soltype , sux.data() ) , "get_dual_solution" );
  for( int j = 0 ; j < numcols ; ++j )
   rc[ j ] -= sux[ j ];
  }

 MILPSolver::write_dual_solution( pi , rc , solc );
 }

/*--------------------------------------------------------------------------*/

bool MOSEKMILPSolver::has_dual_direction( void )
{
 return( ( int_vars == 0 ) &&
	 ( solution_status() == MSK_SOL_STA_PRIM_INFEAS_CER ) );
 }

/*--------------------------------------------------------------------------*/

void MOSEKMILPSolver::get_dual_direction( Configuration * dirc )
{
 if( ! has_dual_direction() )
  throw( std::logic_error( "MOSEKMILPSolver::get_dual_direction: "
			   "no dual unbounded direction available" ) );

 // the certificate of primal infeasibility is in the dual part of the
 // solution: y = slc - suc, with y_i >= 0 on a row whose lower side is used
 // and y_i <= 0 on one whose upper side is, which is the sign of the Farkas
 // certificates of the other :MILPSolver
 std::vector< double > y( numrows , 0 );
 std::vector< double > slc( numrows , 0 ) , suc( numrows , 0 );
 if( numrows > 0 ) {
  check( MSK_gety( task , f_soltype , y.data() ) , "get_dual_direction" );
  check( MSK_getslc( task , f_soltype , slc.data() ) , "get_dual_direction" );
  check( MSK_getsuc( task , f_soltype , suc.data() ) , "get_dual_direction" );
  }

 std::vector< double > slx( numcols , 0 ) , sux( numcols , 0 );
 if( numcols > 0 ) {
  check( MSK_getslx( task , f_soltype , slx.data() ) , "get_dual_direction" );
  check( MSK_getsux( task , f_soltype , sux.data() ) , "get_dual_direction" );
  }

 // the value of the certificate, over the finite sides
 double value = 0;
 if( numrows > 0 ) {
  std::vector< MSKboundkeye > bk( numrows );
  std::vector< double > bl( numrows ) , bu( numrows );
  check( MSK_getconboundslice( task , 0 , numrows , bk.data() , bl.data() ,
			       bu.data() ) , "get_dual_direction" );
  for( int i = 0 ; i < numrows ; ++i ) {
   if( ( bk[ i ] == MSK_BK_LO ) || ( bk[ i ] == MSK_BK_RA ) ||
       ( bk[ i ] == MSK_BK_FX ) )
    value += bl[ i ] * slc[ i ];
   if( ( bk[ i ] == MSK_BK_UP ) || ( bk[ i ] == MSK_BK_RA ) ||
       ( bk[ i ] == MSK_BK_FX ) )
    value -= bu[ i ] * suc[ i ];
   }
  }

 // the multipliers of the columns, c - A' y as CPXdjfrompi() gives them,
 // or - A' y with intHomogeneousDirection, computed out of the matrix
 // rather than taken as slx - sux, so that they are those of y exactly
 std::vector< double > dj( numcols , 0 );
 if( numcols > 0 ) {
  std::vector< MSKboundkeye > bk( numcols );
  std::vector< double > bl( numcols ) , bu( numcols );
  check( MSK_getvarboundslice( task , 0 , numcols , bk.data() , bl.data() ,
			       bu.data() ) , "get_dual_direction" );
  for( int j = 0 ; j < numcols ; ++j ) {
   if( ( bk[ j ] == MSK_BK_LO ) || ( bk[ j ] == MSK_BK_RA ) ||
       ( bk[ j ] == MSK_BK_FX ) )
    value += bl[ j ] * slx[ j ];
   if( ( bk[ j ] == MSK_BK_UP ) || ( bk[ j ] == MSK_BK_RA ) ||
       ( bk[ j ] == MSK_BK_FX ) )
    value -= bu[ j ] * sux[ j ];
   }

  if( ! homogeneous_direction )
   check( MSK_getcslice( task , 0 , numcols , dj.data() ) ,
	  "get_dual_direction" );

  MSKint64t nnz = 0;
  check( MSK_getacolslicenumnz64( task , 0 , numcols , & nnz ) ,
	 "get_dual_direction" );
  if( nnz > 0 ) {
   std::vector< MSKint64t > ptrb( numcols ) , ptre( numcols );
   std::vector< MSKint32t > sub( nnz );
   std::vector< double > val( nnz );
   check( MSK_getacolslice64( task , 0 , numcols , nnz , ptrb.data() ,
			      ptre.data() , sub.data() , val.data() ) ,
	  "get_dual_direction" );
   for( int j = 0 ; j < numcols ; ++j )
    for( auto k = ptrb[ j ] ; k < ptre[ j ] ; ++k )
     dj[ j ] -= val[ k ] * y[ sub[ k ] ];
   }
  }

 f_dual_direction_value = value;

 MILPSolver::write_dual_solution( y , dj , dirc );
 }

/*--------------------------------------------------------------------------*/

void MOSEKMILPSolver::write_lp( const std::string & filename )
{
 pass_hessian();

 // the format is chosen from the extension of the file, so a name that has
 // none is written as an LP
 std::string name = filename;
 if( name.find_last_of( '.' ) == std::string::npos )
  name += ".lp";

 // names the LP format does not take (say, with brackets) are replaced by
 // the generic ones MOSEK makes up
 const auto write = [ & name ]( MSKtask_t t ) {
  auto res = MSK_writedata( t , name.c_str() );
  if( ( res >= MSK_RES_ERR_WRITE_LP_INVALID_VAR_NAMES ) &&
      ( res <= MSK_RES_ERR_WRITE_LP_DUPLICATE_CON_NAMES ) ) {
   MSKint32t generic = 0;
   MSK_getintparam( t , MSK_IPAR_WRITE_GENERIC_NAMES , & generic );
   MSK_putintparam( t , MSK_IPAR_WRITE_GENERIC_NAMES , MSK_ON );
   res = MSK_writedata( t , name.c_str() );
   MSK_putintparam( t , MSK_IPAR_WRITE_GENERIC_NAMES , generic );
   }
  return( res );
  };

 auto r = write( task );

 // a format that has neither free nor ranged rows, as the LP one, gets a
 // copy of the task with no free row, which binds nothing, and with each
 // ranged row lhs <= a'x <= rhs written as the two a'x >= lhs and
 // a'x <= rhs, the latter after all the rows: the model is the same
 if( ( r == MSK_RES_ERR_INVALID_FILE_FORMAT_FOR_FREE_CONSTRAINTS ) ||
     ( r == MSK_RES_ERR_INVALID_FILE_FORMAT_FOR_RANGED_CONSTRAINTS ) ) {
  MSKtask_t copy = nullptr;
  check( MSK_clonetask( task , & copy ) , "write_lp" );

  MSKint32t m = 0 , n = 0;
  MSK_getnumcon( copy , & m );
  MSK_getnumvar( copy , & n );
  std::vector< MSKint32t > free_rows , sub( n );
  std::vector< double > val( n );
  r = MSK_RES_OK;
  for( MSKint32t i = 0 ; ( i < m ) && ( r == MSK_RES_OK ) ; ++i ) {
   MSKboundkeye bk;
   double bl , bu;
   MSK_getconbound( copy , i , & bk , & bl , & bu );
   if( bk == MSK_BK_FR )
    free_rows.push_back( i );
   if( bk != MSK_BK_RA )
    continue;
   MSKint32t nz = 0;
   r = MSK_getarow( copy , i , & nz , sub.data() , val.data() );
   if( r == MSK_RES_OK )
    r = MSK_putconbound( copy , i , MSK_BK_LO , bl , bu );
   if( r == MSK_RES_OK )
    r = MSK_appendcons( copy , 1 );
   MSKint32t last = 0;
   MSK_getnumcon( copy , & last );
   if( r == MSK_RES_OK )
    r = MSK_putarow( copy , last - 1 , nz , sub.data() , val.data() );
   if( r == MSK_RES_OK )
    r = MSK_putconbound( copy , last - 1 , MSK_BK_UP , bl , bu );
   }

  if( ( r == MSK_RES_OK ) && ( ! free_rows.empty() ) )
   r = MSK_removecons( copy , free_rows.size() , free_rows.data() );
  if( r == MSK_RES_OK )
   r = write( copy );
  MSK_deletetask( & copy );
  }

 check( r , "write_lp" );
 }

/*--------------------------------------------------------------------------*/

int MOSEKMILPSolver::get_nodes( void ) const
{
 MSKint32t nodes = 0;
 MSK_getintinf( task , MSK_IINF_MIO_NUM_SOLVED_NODES , & nodes );
 return( nodes );
 }

/*--------------------------------------------------------------------------*/

int MOSEKMILPSolver::get_explored_nodes( void ) const
{
 return( get_nodes() );
 }

/*--------------------------------------------------------------------------*/

long MOSEKMILPSolver::get_left_nodes( void ) const
{
 MSKint32t nodes = 0;
 MSK_getintinf( task , MSK_IINF_MIO_NUM_ACTIVE_NODES , & nodes );
 return( nodes );
 }

/*--------------------------------------------------------------------------*/

bool MOSEKMILPSolver::has_feasible_sol( void )
{
 return( has_primal() );
 }

/*--------------------------------------------------------------------------*/

double MOSEKMILPSolver::get_runtime( void ) const
{
 double t = 0;
 MSK_getdouinf( task , MSK_DINF_OPTIMIZER_TIME , & t );
 return( t );
 }

/*--------------------------------------------------------------------------*/

void MOSEKMILPSolver::add_mip_starts(
				std::vector< std::vector< int > > varidxs ,
				std::vector< std::vector< double > > varvalues )
{
 if( varidxs.empty() || ( numcols == 0 ) )
  return;

 if( varidxs.size() > 1 )
  std::cout << "WARNING: MOSEK supports only one starting solution. "
	       "Only the first one will be used." << std::endl;

 // the values not given are the ones of the last integer solution, if any
 std::vector< double > x( numcols , 0 );
 MSKbooleant def = 0;
 if( ( MSK_solutiondef( task , MSK_SOL_ITG , & def ) == MSK_RES_OK ) && def )
  MSK_getxx( task , MSK_SOL_ITG , x.data() );

 const auto & idx = varidxs[ 0 ];
 const auto & val = varvalues[ 0 ];
 for( std::size_t k = 0 ; k < idx.size() && k < val.size() ; ++k )
  if( ( idx[ k ] >= 0 ) && ( idx[ k ] < numcols ) )
   x[ idx[ k ] ] = val[ k ];

 check( MSK_putxx( task , MSK_SOL_ITG , x.data() ) , "add_mip_starts" );
 check( MSK_putintparam( task , MSK_IPAR_MIO_CONSTRUCT_SOL , MSK_ON ) ,
	"add_mip_starts" );
 }

/*--------------------------------------------------------------------------*/
/*------------------ METHODS FOR MODIFYING THE PROBLEM ---------------------*/
/*--------------------------------------------------------------------------*/

void MOSEKMILPSolver::var_modification( const VariableMod * mod )
{
 // call the method of MILPSolver to update dictionaries (and int_vars)
 MILPSolver::var_modification( mod );

 auto var = static_cast< const ColVariable * >( mod->variable() );
 auto idx = index_of_variable( var );
 if( idx == Inf< int >() )  // the Variable is not (yet) there (?)
  return;                   // nothing to do

 // react to changes in the integrality - - - - - - - - - - - - - - - - - - -
 if( ColVariable::is_integer( mod->old_state() ) !=
     ColVariable::is_integer( mod->new_state() ) )
  check( MSK_putvartype( task , idx , var_type( var ) ) ,
	 "var_modification" );

 // react to fix / unfix- - - - - - - - - - - - - - - - - - - - - - - - - - -
 if( Variable::is_fixed( mod->old_state() ) !=
     Variable::is_fixed( mod->new_state() ) ) {
  if( Variable::is_fixed( mod->new_state() ) )  // fix the variable
   put_var_bound( idx , { double( var->get_value() ) ,
			  double( var->get_value() ) } );
  else                                          // un fix the variable
   put_var_bound( idx , get_problem_bounds( *var ) );
  }
 }  // end( MOSEKMILPSolver::var_modification )

/*--------------------------------------------------------------------------*/

void MOSEKMILPSolver::objective_modification( const ObjectiveMod * mod )
{
 // call the method of MILPSolver to update sense direction
 MILPSolver::objective_modification( mod );

 switch( mod->type() ) {
  case( ObjectiveMod::eSetMin ):
   check( MSK_putobjsense( task , MSK_OBJECTIVE_SENSE_MINIMIZE ) ,
	  "objective_modification" );
   break;
  case( ObjectiveMod::eSetMax ):
   check( MSK_putobjsense( task , MSK_OBJECTIVE_SENSE_MAXIMIZE ) ,
	  "objective_modification" );
   break;
  default:
   throw( std::invalid_argument( "MOSEKMILPSolver::objective_modification: "
				 "invalid type of ObjectiveMod" ) );
  }
 }

/*--------------------------------------------------------------------------*/

void MOSEKMILPSolver::const_modification( const ConstraintMod * mod )
{
 auto * con = dynamic_cast< FRowConstraint * >( mod->constraint() );
 if( ! con )  // this should not happen
  return;     // but in case, nothing to do

 int index = index_of_constraint( con );
 if( index == Inf< int >() )  // the FRowConstraint is not (yet?) there
  return;                     // nothing to do

 std::array< double , 2 > sd;

 switch( mod->type() ) {
  case( ConstraintMod::eRelaxConst ):
   // a relaxed row is a free one, which MOSEK has
   sd = { - Inf< double >() , Inf< double >() };
   break;

  case( ConstraintMod::eEnforceConst ):
  case( RowConstraintMod::eChgLHS ):
  case( RowConstraintMod::eChgRHS ):
  case( RowConstraintMod::eChgBTS ):
   sd = row_sides( con );
   break;

  default:
   throw( std::invalid_argument( "MOSEKMILPSolver::const_modification: "
				 "invalid type of ConstraintMod" ) );
  }

 check( MSK_putconbound( task , index , bound_key( sd[ 0 ] , sd[ 1 ] ) ,
			 msk_side( sd[ 0 ] ) , msk_side( sd[ 1 ] ) ) ,
	"const_modification" );

 }  // end( MOSEKMILPSolver::const_modification )

/*--------------------------------------------------------------------------*/

void MOSEKMILPSolver::bound_modification( const OneVarConstraintMod * mod )
{
 auto con = static_cast< OneVarConstraint * >( mod->constraint() );
 auto var = static_cast< const ColVariable * >( con->get_active_var( 0 ) );
 if( ! var )  // this should never happen
  return;     // but in case, there is nothing to do

 // a fixed Variable has its bounds used to fix it, hence a change of the
 // bounds is ignored here, as in the other :MILPSolver
 if( var->is_fixed() )
  return;

 auto vi = index_of_variable( var );
 if( vi == Inf< int >() )  // the ColVariable has been removed
  return;                  // is strange, but there is nothing to do

 switch( mod->type() ) {
  case( RowConstraintMod::eChgLHS ):
  case( RowConstraintMod::eChgRHS ):
  case( RowConstraintMod::eChgBTS ):
   put_var_bound( vi , get_problem_bounds( *var ) );
   break;
  default:
   throw( std::invalid_argument( "MOSEKMILPSolver::bound_modification: "
				 "invalid type of OneVarConstraintMod" ) );
  }
 }  // end( MOSEKMILPSolver::bound_modification )

/*--------------------------------------------------------------------------*/

bool MOSEKMILPSolver::add_columns( const std::vector< Variable * > & vars ,
                                   const GroupModification * gmod )
{
 // the entries are read before anything is written, so that a group that
 // turns out not to be a plain column can still be handed back to the
 // one-by-one path [see MILPSolver::read_column_group()]
 std::vector< const FRowConstraint * > econs;
 std::vector< Index > evars;
 std::vector< double > evals;
 std::vector< double > cost;

 if( ! read_column_group( vars , gmod , econs , evars , evals , cost ) )
  return( false );

 // the initial objective coefficients carried by the additions are used
 // even if a later queued Modification has removed a Variable from its
 // Function
 std::fill( cost.begin() , cost.end() , 0.0 );
 for( const auto & submod : gmod->sub_Modifications() ) {
  auto modl = dynamic_cast< const LinearFunctionModVarsAddd * >(
                                                             submod.get() );
  if( ( ! modl ) ||
      ( ! dynamic_cast< const Objective * >(
                                      modl->function()->get_Observer() ) ) )
   continue;

  for( Index i = 0 ; i < modl->vars().size() ; ++i ) {
   auto it = std::find( vars.begin() , vars.end() , modl->vars()[ i ] );
   if( it == vars.end() )
    return( false );
   cost[ std::distance( vars.begin() , it ) ] += modl->coeff()[ i ];
   }
  }

 // the columns themselves, each born empty with its bounds and its type,
 // exactly as it is when the Modification arrive one by one
 std::vector< int > vidx( vars.size() );
 for( Index i = 0 ; i < vars.size() ; ++i ) {
  auto var = static_cast< const ColVariable * >( vars[ i ] );
  add_dynamic_variable( var );
  vidx[ i ] = index_of_variable( var );
  }

 // the entries, all in one call
 std::vector< MSKint32t > subi , subj;
 std::vector< double > valij;
 for( Index e = 0 ; e < econs.size() ; ++e ) {
  auto row = index_of_constraint( econs[ e ] );
  if( row == Inf< int >() )  // the Constraint is not (yet?) there, which is
   continue;                  // what the handler does with it one by one
  subi.push_back( row );
  subj.push_back( vidx[ evars[ e ] ] );
  valij.push_back( evals[ e ] );
  }

 if( ! subi.empty() )
  check( MSK_putaijlist( task , subi.size() , subi.data() , subj.data() ,
			 valij.data() ) , "add_columns" );

 std::vector< MSKint32t > oidx;
 std::vector< double > oval;
 for( Index i = 0 ; i < vars.size() ; ++i )
  if( cost[ i ] != 0 ) {  // the column is born with a zero cost
   oidx.push_back( vidx[ i ] );
   oval.push_back( cost[ i ] );
   }

 if( ! oidx.empty() )
  check( MSK_putclist( task , oidx.size() , oidx.data() , oval.data() ) ,
	 "add_columns" );

 return( true );

 }  // end( MOSEKMILPSolver::add_columns )

/*--------------------------------------------------------------------------*/

bool MOSEKMILPSolver::remove_columns( const std::vector< Variable * > & vars ,
                                      const GroupModification * gmod )
{
 // deleting the column takes its coefficients away with it, so not one of
 // the rows it appears in is touched
 for( auto v : vars )
  remove_dynamic_variable( static_cast< const ColVariable * >( v ) );

 return( true );

 }  // end( MOSEKMILPSolver::remove_columns )

/*--------------------------------------------------------------------------*/

bool MOSEKMILPSolver::change_coefficients(
                          const std::vector< const FunctionMod * > & mods )
{
 // nothing is written until the whole group has been read, so that a group
 // that turns out not to be a plain change of coefficients can still be
 // handed back to the one-by-one path
 std::vector< MSKint32t > rows , cols;
 std::vector< double > vals;

 // the costs, where a change is a delta on what the column has
 std::vector< int > ocols;
 std::vector< double > odeltas;

 for( auto mod : mods ) {
  auto modl = dynamic_cast< const C05FunctionModLin * >( mod );
  if( ! modl )
   return( false );

  auto lf = dynamic_cast< const LinearFunction * >( mod->function() );
  if( ! lf )
   return( false );

  auto modlr = dynamic_cast< const C05FunctionModLinRngd * >( modl );
  auto modls = dynamic_cast< const C05FunctionModLinSbst * >( modl );
  if( ( ! modlr ) && ( ! modls ) )
   return( false );

  auto obs = lf->get_Observer();

  if( auto con = dynamic_cast< const FRowConstraint * >( obs ) ) {
   auto idxs = modlr ? lf->map_index( modl->vars() , modlr->range() )
                     : lf->map_index( modl->vars() , modls->subset() );
   auto row = index_of_constraint( con );
   if( row == Inf< int >() )  // the Constraint is not (yet?) there, which is
    continue;                 // what the handler does with it one by one

   auto & cp = lf->get_v_var();
   for( Block::Index i = 0 ; i < modl->vars().size() ; ++i )
    if( auto idx = idxs[ i ] ; idx < Inf< Index >() ) {
     auto var = static_cast< const ColVariable * >( modl->vars()[ i ] );
     rows.push_back( row );
     cols.push_back( index_of_variable( var ) );
     vals.push_back( cp[ idx ].second );
     }

   continue;
   }

  if( dynamic_cast< const Objective * >( obs ) ) {
   // the deltas describe queued changes, not the live Function state
   for( Block::Index i = 0 ; i < modl->vars().size() ; ++i ) {
    auto var = static_cast< const ColVariable * >( modl->vars()[ i ] );
    if( auto idx = index_of_variable( var ) ; idx < Inf< int >() ) {
     ocols.push_back( idx );
     odeltas.push_back( modl->delta()[ i ] );
     }
    }

   continue;
   }

  return( false );  // a Function of something this Solver does not have
  }

 // an entry named twice takes the last value, as one at a time
 if( ! rows.empty() )
  check( MSK_putaijlist( task , rows.size() , rows.data() , cols.data() ,
			 vals.data() ) , "change_coefficients" );

 if( ! ocols.empty() )
  change_obj_costs( ocols , odeltas , true );

 return( true );

 }  // end( MOSEKMILPSolver::change_coefficients )

/*--------------------------------------------------------------------------*/

bool MOSEKMILPSolver::change_variables(
                          const std::vector< const VariableMod * > & mods )
{
 // a change of integrality is one call per column anyway, and it has to be
 // executed in the order it comes with the fixings: the batch is refused
 for( auto mod : mods )
  if( ColVariable::is_integer( mod->old_state() ) !=
      ColVariable::is_integer( mod->new_state() ) )
   return( false );

 // a column named twice keeps the last state, exactly as it does when the
 // Modification arrive one by one
 std::map< int , std::array< double , 2 > > bnds;

 for( auto mod : mods ) {
  // the bookkeeping of the base class is the same it does one by one
  MILPSolver::var_modification( mod );

  if( Variable::is_fixed( mod->old_state() ) ==
      Variable::is_fixed( mod->new_state() ) )
   continue;  // nothing that reaches the model

  auto var = static_cast< const ColVariable * >( mod->variable() );
  auto vi = index_of_variable( var );
  if( vi == Inf< int >() )  // the Variable is not (yet) there
   continue;

  if( Variable::is_fixed( mod->new_state() ) )  // fix it at its value
   bnds[ vi ] = { double( var->get_value() ) , double( var->get_value() ) };
  else                                          // give its bounds back
   bnds[ vi ] = get_problem_bounds( *var );
  }

 std::vector< int > idx;
 std::vector< std::array< double , 2 > > bd;
 for( const auto & el : bnds ) {
  idx.push_back( el.first );
  bd.push_back( el.second );
  }
 put_var_bounds( idx , bd );

 return( true );

 }  // end( MOSEKMILPSolver::change_variables )

/*--------------------------------------------------------------------------*/

bool MOSEKMILPSolver::change_bounds(
                   const std::vector< const OneVarConstraintMod * > & mods )
{
 // the last word on a column the group touches more than once is the one
 // that goes out
 std::map< int , std::array< double , 2 > > bnds;

 for( auto mod : mods ) {
  switch( mod->type() ) {
   case( RowConstraintMod::eChgLHS ):
   case( RowConstraintMod::eChgRHS ):
   case( RowConstraintMod::eChgBTS ):
    break;
   default:  // nothing has been written yet, so the group can still be
    return( false );  // taken apart and each of them throw on its own
   }

  auto con = static_cast< OneVarConstraint * >( mod->constraint() );
  auto var = static_cast< const ColVariable * >( con->get_active_var( 0 ) );
  if( ! var )  // this should never happen
   continue;   // but in case, there is nothing to do

  // a fixed Variable has its bounds used to fix it [see
  // bound_modification()]
  if( var->is_fixed() )
   continue;

  auto vi = index_of_variable( var );
  if( vi == Inf< int >() )  // the ColVariable has been removed
   continue;

  bnds[ vi ] = get_problem_bounds( *var );
  }

 std::vector< int > idx;
 std::vector< std::array< double , 2 > > bd;
 for( const auto & el : bnds ) {
  idx.push_back( el.first );
  bd.push_back( el.second );
  }
 put_var_bounds( idx , bd );

 return( true );

 }  // end( MOSEKMILPSolver::change_bounds )

/*--------------------------------------------------------------------------*/

bool MOSEKMILPSolver::change_sides(
                      const std::vector< const RowConstraintMod * > & mods )
{
 std::map< int , std::array< double , 2 > > sides;

 for( auto mod : mods ) {
  switch( mod->type() ) {
   case( RowConstraintMod::eChgLHS ):
   case( RowConstraintMod::eChgRHS ):
   case( RowConstraintMod::eChgBTS ):
    break;
   default:  // relaxing and enforcing a Constraint is not a side, and
    return( false );  // nothing has been written yet
   }

  auto con = dynamic_cast< FRowConstraint * >( mod->constraint() );
  if( ! con )  // this should not happen
   return( false );

  auto index = index_of_constraint( con );
  if( index == Inf< int >() )  // the FRowConstraint is not (yet?) there,
   continue;                   // which is what the handler does one by one

  sides[ index ] = row_sides( con );
  }

 std::vector< int > idx;
 std::vector< std::array< double , 2 > > sd;
 for( const auto & el : sides ) {
  idx.push_back( el.first );
  sd.push_back( el.second );
  }
 put_con_bounds( idx , sd );

 return( true );

 }  // end( MOSEKMILPSolver::change_sides )

/*--------------------------------------------------------------------------*/

void MOSEKMILPSolver::objective_function_modification(
						     const FunctionMod * mod )
{
 auto f = mod->function();

 // C05FunctionModLin - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 if( auto modl = dynamic_cast< const C05FunctionModLin * >( mod ) ) {

  if( dynamic_cast< const LinearFunction * >( f ) ) {
   if( ( ! dynamic_cast< const C05FunctionModLinRngd * >( modl ) ) &&
       ( ! dynamic_cast< const C05FunctionModLinSbst * >( modl ) ) )
    throw( std::logic_error(
               "MOSEKMILPSolver: unknown type of C05FunctionModLin" ) );

   // the delta() of the Modification is applied to the cost, which is
   // (potentially) a sum of terms; queued deltas are applied even if a
   // later Modification has removed the Variable from the live Function
   std::vector< int > cidx;
   std::vector< double > delta;
   for( Block::Index i = 0 ; i < modl->vars().size() ; ++i ) {
    auto var = static_cast< const ColVariable * >( modl->vars()[ i ] );
    if( auto idx = index_of_variable( var ) ; idx < Inf< int >() ) {
     cidx.push_back( idx );
     delta.push_back( modl->delta()[ i ] );
     }
    }

   change_obj_costs( cidx , delta , true );
   return;
   }

  if( auto qf = dynamic_cast< const DQuadFunction * >( f ) ) {
   Subset idxs;
   if( auto modlr = dynamic_cast< const C05FunctionModLinRngd * >( modl ) )
    idxs = qf->map_index( modl->vars() , modlr->range() );
   else
    if( auto modls = dynamic_cast< const C05FunctionModLinSbst * >( modl ) )
     idxs = qf->map_index( modl->vars() , modls->subset() );
    else
     throw( std::logic_error(
                "MOSEKMILPSolver: unknown type of C05FunctionModLin" ) );

   std::vector< int > cidx;
   std::vector< double > delta;
   auto idxit = idxs.begin();
   for( Block::Index i = 0 ; i < modl->vars().size() ; ++i )
    if( *( idxit++ ) < Inf< Index >() ) {
     cidx.push_back( index_of_variable(
                   static_cast< const ColVariable * >( modl->vars()[ i ] ) ) );
     delta.push_back( modl->delta()[ i ] );
     }

   change_obj_costs( cidx , delta , true );
   return;
   }

  throw( std::invalid_argument(
             "MOSEKMILPSolver: unknown type of Objective Function" ) );
  }

 // C05FunctionMod- - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 if( auto modl = dynamic_cast< const C05FunctionMod * >( mod ) ) {

  if( modl->type() == C05FunctionMod::NothingChanged ) {
   const auto shift = modl->shift();

   if( ( shift == FunctionMod::INFshift ) ||
       ( shift == - FunctionMod::INFshift ) || std::isnan( shift ) )
    throw( std::logic_error(
               "MOSEKMILPSolver::objective_function_modification: "
               "unexpected *C05FunctionMod* from Objective Function" ) );

   constant_value += shift;
   return;
   }

  auto qf = dynamic_cast< const QuadFunction * >( f );
  auto dqf = dynamic_cast< const DQuadFunction * >( f );
  if( ( ! qf ) && ( ! dqf ) )
   throw( std::logic_error(
              "MOSEKMILPSolver::objective_function_modification: "
              "unexpected *C05FunctionMod* from Linear Objective" ) );

  auto fqf = qf ? qf : dqf;

  if( auto modld = dynamic_cast< const DQuadFunctionModRngd * >( modl ) ;
      modld || dynamic_cast< const DQuadFunctionModSbst * >( modl ) ) {
   // the delta() gives the difference between the new and the old value of
   // both the linear and the quadratic coefficient
   Subset idxs;
   c_Vec_p_Var * vars;
   c_v_coeff_pair * delta_coeff;
   if( modld ) {
    idxs = fqf->map_index( modld->vars() , modld->range() );
    vars = & modld->vars();
    delta_coeff = & modld->delta();
    }
   else {
    auto modls = static_cast< const DQuadFunctionModSbst * >( modl );
    idxs = fqf->map_index( modls->vars() , modls->subset() );
    vars = & modls->vars();
    delta_coeff = & modls->delta();
    }

   std::vector< int > cidx;
   std::vector< double > delta;
   auto idxit = idxs.begin();
   auto dcoeffit = delta_coeff->begin();

   for( auto v : *vars ) {
    const auto i = *( idxit++ );
    if( i < Inf< Index >() ) {
     auto c = index_of_variable( static_cast< const ColVariable * >( v ) );
     cidx.push_back( c );
     delta.push_back( std::get< 0 >( *dcoeffit ) );
     // the diagonal takes the value that the Function has, not the delta:
     // a sequence of changes that brings a coefficient back would otherwise
     // leave a rounding error on it, which on a singular Hessian is enough
     // for MOSEK to refuse the objective as not convex
     set_hessian( c , c , 2 * fqf->get_quadratic_coefficient( i ) );
     }
    ++dcoeffit;
    }

   change_obj_costs( cidx , delta , true );
   return;
   }

  if( auto modlq = dynamic_cast< const QuadFunctionModSbst * >( modl ) ) {
   // one coefficient at a time for each Modification
   c_Vec_p_Var vars = modlq->vars();
   if( vars.size() != 2 )
    throw( std::logic_error(
               "MOSEKMILPSolver::objective_function_modification: "
               "expected single coefficient Modification in "
               "QuadFunctionModSbst" ) );

   int idx1 = index_of_variable( static_cast< ColVariable * >( vars[ 0 ] ) );
   int idx2 = index_of_variable( static_cast< ColVariable * >( vars[ 1 ] ) );
   add_to_hessian( idx1 , idx2 , modlq->delta() );
   return;
   }

  throw( std::logic_error(
             "MOSEKMILPSolver: unknown type of *QuadFunctionMod*" ) );
  }

 }  // end( MOSEKMILPSolver::objective_function_modification )

/*--------------------------------------------------------------------------*/

void MOSEKMILPSolver::constraint_function_modification(
						     const FunctionMod * mod )
{
 auto lf = dynamic_cast< const LinearFunction * >( mod->function() );
 if( ! lf )
  return;

 auto con = dynamic_cast< const FRowConstraint * >( lf->get_Observer() );
 if( ! con )
  return;

 auto row = index_of_constraint( con );
 if( row == Inf< int >() )  // the constraint is not (yet?) there
  return;

 auto modl = dynamic_cast< const C05FunctionModLin * >( mod );
 if( ! modl )
  throw( std::logic_error(
             "MOSEKMILPSolver::constraint_function_modification: "
             "unexpected *C05FunctionModLin* from FRowConstraint" ) );

 Subset idxs;
 if( auto modlr = dynamic_cast< const C05FunctionModLinRngd * >( modl ) )
  idxs = lf->map_index( modl->vars() , modlr->range() );
 else
  if( auto modls = dynamic_cast< const C05FunctionModLinSbst * >( modl ) )
   idxs = lf->map_index( modl->vars() , modls->subset() );
  else
   throw( std::logic_error(
              "MOSEKMILPSolver::constraint_function_modification: "
              "unknown type of C05FunctionModLin" ) );

 auto idxit = idxs.begin();
 auto & cp = lf->get_v_var();
 std::vector< MSKint32t > rows , cols;
 std::vector< double > vals;

 for( auto v : modl->vars() )
  if( auto idx = *( idxit++ ) ; idx < Inf< Index >() ) {
   rows.push_back( row );
   cols.push_back( index_of_variable(
                               static_cast< const ColVariable * >( v ) ) );
   vals.push_back( cp[ idx ].second );
   }

 if( ! rows.empty() )
  check( MSK_putaijlist( task , rows.size() , rows.data() , cols.data() ,
			 vals.data() ) , "constraint_function_modification" );

 }  // end( MOSEKMILPSolver::constraint_function_modification )

/*--------------------------------------------------------------------------*/

void MOSEKMILPSolver::objective_fvars_modification(
						 const FunctionModVars * mod )
{
 // this is called in response to Variable being added to / removed from the
 // Objective; all Variable are supposed to exist at the time this is called,
 // so index_of_variable() is always correct

 if( mod->vars().empty() )
  return;

 if( ( ! dynamic_cast< const LinearFunctionModVarsAddd * >( mod ) ) &&
     ( ! dynamic_cast< const DQuadFunctionModVarsAddd * >( mod ) ) &&
     ( ! dynamic_cast< const QuadFunctionModVarsAddd * >( mod ) ) &&
     ( ! dynamic_cast< const C05FunctionModVarsRngd * >( mod ) ) &&
     ( ! dynamic_cast< const C05FunctionModVarsSbst * >( mod ) ) )
  throw( std::invalid_argument(
             "MOSEKMILPSolver: this type of FunctionModVars is not handled"
				) );

 std::vector< int > indices;
 std::vector< double > values;

 auto f = mod->function();

 if( dynamic_cast< const LinearFunction * >( f ) ) {
  // the coefficients of an addition and the deltas of a removal refer to
  // the issuing Function only, and they are applied to the total cost in
  // arrival order
  const Function::Vec_FunctionValue * delta = nullptr;
  if( mod->added() ) {
   auto modl = dynamic_cast< const LinearFunctionModVarsAddd * >( mod );
   if( ! modl )
    throw( std::invalid_argument(
               "MOSEKMILPSolver: invalid LinearFunction addition" ) );
   delta = & modl->coeff();
   }
  else
   if( auto modr = dynamic_cast< const C05FunctionModVarsRngd * >( mod ) ) {
    if( modr->has_linear_delta() )
     delta = & modr->linear_delta();
    }
   else
    if( auto mods = dynamic_cast< const C05FunctionModVarsSbst * >( mod ) )
     if( mods->has_linear_delta() )
      delta = & mods->linear_delta();

  for( Block::Index i = 0 ; i < mod->vars().size() ; ++i ) {
   auto var = static_cast< const ColVariable * >( mod->vars()[ i ] );
   if( auto idx = index_of_variable( var ) ; idx < Inf< int >() ) {
    indices.push_back( idx );
    // removals without a delta keep zeroing the cost
    values.push_back( delta ? ( *delta )[ i ] : 0.0 );
    }
   }

  change_obj_costs( indices , values , delta != nullptr , ! mod->added() );
  return;
  }

 if( auto qf = dynamic_cast< const QuadFunction * >( f ) ) {
  if( ! mod->added() ) {  // the Variable leave the Objective, and so does
   // every term they appear in
   for( auto v : mod->vars() )
    if( auto idx = index_of_variable( static_cast< const ColVariable * >( v ) )
        ; idx < Inf< int >() ) {
     for( int j = 0 ; j <= idx ; ++j )
      set_hessian( idx , j , 0 );
     for( int i = idx + 1 ; i < int( q_hessian.size() ) ; ++i )
      set_hessian( i , idx , 0 );
     indices.push_back( idx );
     values.push_back( 0 );
     }

   change_obj_costs( indices , values , false );
   return;
   }

  auto modq = dynamic_cast< const QuadFunctionModVarsAddd * >( mod );
  if( ! modq )
   throw( std::invalid_argument(
              "MOSEKMILPSolver: unexpected type of Objective Function "
              "Modification" ) );

  // the off-diagonal terms; the diagonal ones and the linear coefficients
  // follow below, a QuadFunction being a DQuadFunction
  for( auto t : modq->od_terms() ) {
   auto var1 = qf->get_active_var( std::get< 0 >( t ) );
   auto var2 = qf->get_active_var( std::get< 1 >( t ) );
   int idx1 = index_of_variable( static_cast< const ColVariable * >( var1 ) );
   int idx2 = index_of_variable( static_cast< const ColVariable * >( var2 ) );
   if( ( idx1 < Inf< int >() ) && ( idx2 < Inf< int >() ) )
    add_to_hessian( idx1 , idx2 , std::get< 2 >( t ) );
   }
  }

 if( dynamic_cast< const DQuadFunction * >( f ) ) {
  if( ! mod->added() ) {  // the Variable leave the Objective
   for( auto v : mod->vars() )
    if( auto idx = index_of_variable( static_cast< const ColVariable * >( v ) )
        ; idx < Inf< int >() ) {
     set_hessian( idx , idx , 0 );
     indices.push_back( idx );
     values.push_back( 0 );
     }

   change_obj_costs( indices , values , false );
   return;
   }

  auto modq = dynamic_cast< const DQuadFunctionModVarsAddd * >( mod );
  if( ! modq )
   throw( std::invalid_argument(
              "MOSEKMILPSolver: unexpected type of Objective Function "
              "Modification" ) );

  for( Block::Index i = 0 ; i < mod->vars().size() ; ++i ) {
   auto var = static_cast< const ColVariable * >( mod->vars()[ i ] );
   if( auto idx = index_of_variable( var ) ; idx < Inf< int >() ) {
    indices.push_back( idx );
    values.push_back( modq->coeff()[ i ].first );
    add_to_hessian( idx , idx , 2 * modq->coeff()[ i ].second );
    }
   }

  change_obj_costs( indices , values , true );
  return;
  }

 throw( std::invalid_argument(
            "MOSEKMILPSolver::objective_fvars_modification: "
            "unknown type of Objective Function" ) );

 }  // end( MOSEKMILPSolver::objective_fvars_modification )

/*--------------------------------------------------------------------------*/

void MOSEKMILPSolver::constraint_fvars_modification(
						const FunctionModVars * mod )
{
 auto lf = dynamic_cast< const LinearFunction * >( mod->function() );
 if( ! lf )
  return;

 if( mod->vars().empty() )
  return;

 if( ( ! dynamic_cast< const C05FunctionModVarsAddd * >( mod ) ) &&
     ( ! dynamic_cast< const C05FunctionModVarsRngd * >( mod ) ) &&
     ( ! dynamic_cast< const C05FunctionModVarsSbst * >( mod ) ) )
  throw( std::invalid_argument(
             "MOSEKMILPSolver: this type of FunctionModVars is not handled"
				) );

 auto con = dynamic_cast< const FRowConstraint * >( lf->get_Observer() );
 if( ! con )  // non-FRowConstraint Mods are silently ignored
  return;

 auto cidx = index_of_constraint( con );
 if( cidx == Inf< int >() )
  return;

 auto nav = lf->get_num_active_var();

 // the Modification are managed asynchronously with the model changes: the
 // added Variable may no longer be active in the LinearFunction, in which
 // case their coefficient is 0, as that of the removed ones
 std::vector< MSKint32t > rows , cols;
 std::vector< double > vals;
 for( auto v : mod->vars() ) {
  auto var = static_cast< const ColVariable * >( v );
  if( auto vidx = index_of_variable( var ) ; vidx < Inf< int >() ) {
   double value = 0;
   if( mod->added() ) {
    auto idx = lf->is_active( var );
    value = idx < nav ? lf->get_coefficient( idx ) : 0;
    }
   rows.push_back( cidx );
   cols.push_back( vidx );
   vals.push_back( value );
   }
  }

 if( ! rows.empty() )
  check( MSK_putaijlist( task , rows.size() , rows.data() , cols.data() ,
			 vals.data() ) , "constraint_fvars_modification" );

 }  // end( MOSEKMILPSolver::constraint_fvars_modification )

/*--------------------------------------------------------------------------*/

void MOSEKMILPSolver::add_dynamic_constraint( const FRowConstraint * con )
{
 add_dynamic_constraints( { con } );
 }

/*--------------------------------------------------------------------------*/

void MOSEKMILPSolver::add_dynamic_constraints(
          const std::vector< const FRowConstraint * > & cons )
{
 if( cons.empty() )
  return;

 // the rows are appended at the end of the task, in the order of cons,
 // which is the one the dictionaries give them
 std::vector< MSKint64t > ptrb , ptre;
 std::vector< MSKint32t > index;
 std::vector< double > value;
 std::vector< std::array< double , 2 > > sd;

 const auto idx_of = [ this ]( const ColVariable * v ) {
                      return( index_of_variable( v ) );
                      };

 std::vector< const FRowConstraint * > rows;
 for( auto * c : cons ) {
  if( ! c )
   continue;

  auto lf = dynamic_cast< const LinearFunction * >( c->get_function() );
  if( ! lf )
   throw( std::invalid_argument(
             "MOSEKMILPSolver::add_dynamic_constraints: "
             "the FRowConstraint is not linear" ) );

  rows.push_back( c );
  ptrb.push_back( index.size() );
  scatter_lf_to_csr( lf , idx_of , index , value );
  ptre.push_back( index.size() );
  sd.push_back( row_sides( c ) );
  }

 if( rows.empty() )
  return;

 // update the dictionaries (only) for every row
 for( auto * c : rows )
  MILPSolver::add_dynamic_constraint( c );

 MSKint32t first = 0;
 MSK_getnumcon( task , & first );
 const MSKint32t n = rows.size();

 check( MSK_appendcons( task , n ) , "add_dynamic_constraints" );
 if( ! index.empty() )
  check( MSK_putarowslice64( task , first , first + n , ptrb.data() ,
			     ptre.data() , index.data() , value.data() ) ,
	 "add_dynamic_constraints" );

 std::vector< int > idx( n );
 std::iota( idx.begin() , idx.end() , first );
 put_con_bounds( idx , sd );

 }  // end( MOSEKMILPSolver::add_dynamic_constraints )

/*--------------------------------------------------------------------------*/

void MOSEKMILPSolver::add_dynamic_variable( const ColVariable * var )
{
 // call the method of MILPSolver to update dictionaries (and int_vars)
 MILPSolver::add_dynamic_variable( var );

 auto bd = get_problem_bounds( *var );
 fix_bounds( *var , bd );  // a Variable may come already fixed

 // the column is born with no entry and a zero cost
 check( MSK_appendvars( task , 1 ) , "add_dynamic_variable" );
 const int j = numcols - 1;
 put_var_bound( j , bd );
 check( MSK_putvartype( task , j , var_type( var ) ) ,
	"add_dynamic_variable" );

 // the column is born with no quadratic term, but a quadratic part that is
 // there is given again with one more column
 q_hessian.resize( numcols );

 }  // end( MOSEKMILPSolver::add_dynamic_variable )

/*--------------------------------------------------------------------------*/

void MOSEKMILPSolver::add_dynamic_bound( const OneVarConstraint * con )
{
 auto var = static_cast< const ColVariable * >( con->get_active_var( 0 ) );
 if( ! var )
  throw( std::logic_error( "MOSEKMILPSolver::add_dynamic_bound: "
			   "added a bound on no Variable" ) );

 auto idx = index_of_variable( var );
 if( idx == Inf< int >() )
  throw( std::logic_error( "MOSEKMILPSolver::add_dynamic_bound: "
			   "added a bound on unknown Variable" ) );

 if( var->is_fixed() )  // the bounds of a fixed Variable are its value
  return;

 put_var_bound( idx , get_problem_bounds( *var ) );
 }

/*--------------------------------------------------------------------------*/

void MOSEKMILPSolver::remove_dynamic_constraint( const FRowConstraint * con )
{
 int index = index_of_dynamic_constraint( con );
 if( index == Inf< int >() )
  throw( std::runtime_error( "MOSEKMILPSolver::remove_dynamic_constraint: "
			     "dynamic constraint not found" ) );

 MSKint32t sub = index;
 check( MSK_removecons( task , 1 , & sub ) , "remove_dynamic_constraint" );

 // call the method of MILPSolver to update the dictionaries (only)
 MILPSolver::remove_dynamic_constraint( con );
 }

/*--------------------------------------------------------------------------*/

void MOSEKMILPSolver::remove_dynamic_variable( const ColVariable * var )
{
 int index = index_of_dynamic_variable( var );
 if( index == Inf< int >() )
  throw( std::runtime_error( "MOSEKMILPSolver::remove_dynamic_variable: "
			     "dynamic variable not found" ) );

 MSKint32t sub = index;
 check( MSK_removevars( task , 1 , & sub ) , "remove_dynamic_variable" );

 // MOSEK takes the column away from its quadratic part, the copy follows
 remove_hessian_column( index );

 // call the method of MILPSolver to update the dictionaries (only)
 MILPSolver::remove_dynamic_variable( var );
 }

/*--------------------------------------------------------------------------*/

void MOSEKMILPSolver::remove_dynamic_bound( const OneVarConstraint * con )
{
 // this only works because the removed OneVarConstraint is not cleared,
 // hence the ColVariable it was about can be found
 auto var = static_cast< const ColVariable * >( con->get_active_var( 0 ) );
 if( ! var )  // this should never happen
  return;     // but in case, there is nothing to do

 int idx = index_of_variable( var );
 if( idx == Inf< int >() )  // the ColVariable has been removed
  return;                   // is strange, but there is nothing to do

 if( var->is_fixed() )  // the bounds of a fixed Variable are its value
  return;

 put_var_bound( idx , get_problem_bounds( *var ) );
 }

/*--------------------------------------------------------------------------*/
/*------------------- METHODS FOR HANDLING THE PARAMETERS ------------------*/
/*--------------------------------------------------------------------------*/

std::string MOSEKMILPSolver::msk_int_par_map( idx_type par ) const
{
 switch( par ) {
  case( intMaxIter ):   return( "MSK_IPAR_MIO_MAX_NUM_BRANCHES" );
  case( intMaxSol ):    return( "MSK_IPAR_MIO_MAX_NUM_SOLUTIONS" );
  case( intMaxThread ): return( "MSK_IPAR_NUM_THREADS" );
  }

 if( ( par >= intFirstMOSEKPar ) && ( par < intLastAlgParMOSEKS ) )
  return( SMSpp_to_MOSEK_int_pars[ par - intFirstMOSEKPar ] );

 return( "" );
 }

/*--------------------------------------------------------------------------*/

std::vector< std::string > MOSEKMILPSolver::msk_dbl_par_map( idx_type par )
 const
{
 using v_str = std::vector< std::string >;

 switch( par ) {
  case( dblMaxTime ): return( v_str{ "MSK_DPAR_OPTIMIZER_MAX_TIME" } );
  case( dblRelAcc ):  return( v_str{ "MSK_DPAR_MIO_TOL_REL_GAP" } );
  case( dblAbsAcc ):  return( v_str{ "MSK_DPAR_MIO_TOL_ABS_GAP" } );
  case( dblFAccSol ):
   return( v_str{ "MSK_DPAR_INTPNT_TOL_PFEAS" ,
		  "MSK_DPAR_INTPNT_QO_TOL_PFEAS" , "MSK_DPAR_BASIS_TOL_X" ,
		  "MSK_DPAR_MIO_TOL_FEAS" } );
  }

 if( ( par >= dblFirstMOSEKPar ) && ( par < dblLastAlgParMOSEKS ) )
  return( v_str{ SMSpp_to_MOSEK_dbl_pars[ par - dblFirstMOSEKPar ] } );

 return( v_str() );
 }

/*--------------------------------------------------------------------------*/

void MOSEKMILPSolver::set_par( idx_type par , int value )
{
 if( par == intLogVerb ) {
  log_verbosity = value;
  check( MSK_putintparam( task , MSK_IPAR_LOG , std::max( value , 0 ) ) ,
	 "set_par" );
  link_log();
  return;
  }

 const auto name = msk_int_par_map( par );
 if( ! name.empty() ) {
  auto r = MSK_putnaintparam( task , name.c_str() , value );
  if( r != MSK_RES_OK )
   throw( std::invalid_argument( "MOSEKMILPSolver::set_par: " + name +
				 " = " + std::to_string( value ) +
				 " refused with " + code_text( r ) ) );
  return;
  }

 MILPSolver::set_par( par , value );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

void MOSEKMILPSolver::set_par( idx_type par , double value )
{
 switch( par ) {
  case( dblUpCutOff ): UpCutOff = value; return;
  case( dblLwCutOff ): LwCutOff = value; return;
  case( dblMaxTime ):
   // MOSEK has no time limit when it is negative
   if( value >= Inf< double >() )
    value = -1;
   else
    value = std::max( value , 0.0 );
   break;
  }

 const auto names = msk_dbl_par_map( par );
 if( ! names.empty() ) {
  for( const auto & name : names ) {
   auto r = MSK_putnadouparam( task , name.c_str() , value );
   if( r != MSK_RES_OK )
    throw( std::invalid_argument( "MOSEKMILPSolver::set_par: " + name +
				  " = " + std::to_string( value ) +
				  " refused with " + code_text( r ) ) );
   }
  return;
  }

 MILPSolver::set_par( par , value );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

void MOSEKMILPSolver::set_par( idx_type par , std::string && value )
{
 if( par == strLogFileName ) {
  f_log_file = value;
  link_log();
  return;
  }

 if( ( par >= strFirstMOSEKPar ) && ( par < strLastAlgParMOSEKS ) ) {
  const auto & name = SMSpp_to_MOSEK_str_pars[ par - strFirstMOSEKPar ];
  auto r = MSK_putnastrparam( task , name.c_str() , value.c_str() );
  if( r != MSK_RES_OK )
   throw( std::invalid_argument( "MOSEKMILPSolver::set_par: " + name +
				 " = " + value + " refused with " +
				 code_text( r ) ) );
  return;
  }

 MILPSolver::set_par( par , std::move( value ) );
 }

/*--------------------------------------------------------------------------*/

Solver::idx_type MOSEKMILPSolver::get_num_int_par( void ) const {
 return( MILPSolver::get_num_int_par()
	 + intLastAlgParMOSEKS - intLastAlgParMILP );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

Solver::idx_type MOSEKMILPSolver::get_num_dbl_par( void ) const {
 return( MILPSolver::get_num_dbl_par()
	 + dblLastAlgParMOSEKS - dblLastAlgParMILP );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

Solver::idx_type MOSEKMILPSolver::get_num_str_par( void ) const {
 return( MILPSolver::get_num_str_par()
	 + strLastAlgParMOSEKS - strLastAlgParMILP );
 }

/*--------------------------------------------------------------------------*/

// the defaults of MOSEK are those of a task that nobody has touched, which
// is made once and kept for the whole process

static MSKtask_t default_task( MSKenv_t env )
{
 static std::mutex m;
 static MSKtask_t dt = nullptr;
 std::lock_guard< std::mutex > lock( m );
 if( ! dt )
  if( MSK_maketask( env , 0 , 0 , & dt ) != MSK_RES_OK )
   dt = nullptr;
 return( dt );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

int MOSEKMILPSolver::get_dflt_int_par( idx_type par ) const
{
 if( par == intLogVerb )
  return( 0 );

 const auto name = msk_int_par_map( par );
 if( ! name.empty() )
  if( auto dt = default_task( env ) ) {
   MSKint32t value = 0;
   if( MSK_getnaintparam( dt , name.c_str() , & value ) == MSK_RES_OK )
    return( value );
   }

 return( MILPSolver::get_dflt_int_par( par ) );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

double MOSEKMILPSolver::get_dflt_dbl_par( idx_type par ) const
{
 switch( par ) {
  case( dblUpCutOff ): return( Inf< double >() );
  case( dblLwCutOff ): return( - Inf< double >() );
  case( dblMaxTime ):  return( Inf< double >() );
  }

 const auto names = msk_dbl_par_map( par );
 if( ! names.empty() )
  if( auto dt = default_task( env ) ) {
   double value = 0;
   if( MSK_getnadouparam( dt , names[ 0 ].c_str() , & value ) == MSK_RES_OK )
    return( value );
   }

 return( MILPSolver::get_dflt_dbl_par( par ) );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

const std::string & MOSEKMILPSolver::get_dflt_str_par( idx_type par ) const
{
 static std::string value;

 if( ( par >= strFirstMOSEKPar ) && ( par < strLastAlgParMOSEKS ) ) {
  value.clear();
  if( auto dt = default_task( env ) ) {
   std::array< char , MSK_MAX_STR_LEN > buf = {};
   MSKint32t len = 0;
   if( MSK_getnastrparam( dt , SMSpp_to_MOSEK_str_pars[
                            par - strFirstMOSEKPar ].c_str() ,
                          buf.size() , & len , buf.data() ) == MSK_RES_OK )
    value = buf.data();
   }
  return( value );
  }

 return( MILPSolver::get_dflt_str_par( par ) );
 }

/*--------------------------------------------------------------------------*/

int MOSEKMILPSolver::get_int_par( idx_type par ) const
{
 if( par == intLogVerb )
  return( log_verbosity );

 const auto name = msk_int_par_map( par );
 if( ! name.empty() ) {
  MSKint32t value = 0;
  if( MSK_getnaintparam( task , name.c_str() , & value ) == MSK_RES_OK )
   return( value );
  }

 return( MILPSolver::get_int_par( par ) );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

double MOSEKMILPSolver::get_dbl_par( idx_type par ) const
{
 switch( par ) {
  case( dblUpCutOff ): return( UpCutOff );
  case( dblLwCutOff ): return( LwCutOff );
  }

 const auto names = msk_dbl_par_map( par );
 if( ! names.empty() ) {
  double value = 0;
  if( MSK_getnadouparam( task , names[ 0 ].c_str() , & value ) ==
      MSK_RES_OK ) {
   if( ( par == dblMaxTime ) && ( value < 0 ) )
    return( Inf< double >() );
   return( value );
   }
  }

 return( MILPSolver::get_dbl_par( par ) );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

const std::string & MOSEKMILPSolver::get_str_par( idx_type par ) const
{
 static std::string value;

 if( par == strLogFileName )
  return( f_log_file );

 if( ( par >= strFirstMOSEKPar ) && ( par < strLastAlgParMOSEKS ) ) {
  std::array< char , MSK_MAX_STR_LEN > buf = {};
  MSKint32t len = 0;
  value.clear();
  if( MSK_getnastrparam( task , SMSpp_to_MOSEK_str_pars[
                           par - strFirstMOSEKPar ].c_str() ,
                         buf.size() , & len , buf.data() ) == MSK_RES_OK )
   value = buf.data();
  return( value );
  }

 return( MILPSolver::get_str_par( par ) );
 }

/*--------------------------------------------------------------------------*/

Solver::idx_type MOSEKMILPSolver::int_par_str2idx(
					     const std::string & name ) const
{
 idx_type idx = MILPSolver::int_par_str2idx( name );
 if( idx < Inf< idx_type >() )
  return( idx );

 auto it = std::find( SMSpp_to_MOSEK_int_pars.begin() ,
		      SMSpp_to_MOSEK_int_pars.end() , name );
 if( it != SMSpp_to_MOSEK_int_pars.end() )
  return( MOSEK_to_SMSpp_int_pars[ std::distance(
                             SMSpp_to_MOSEK_int_pars.begin() , it ) ].second );

 return( Inf< idx_type >() );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

const std::string & MOSEKMILPSolver::int_par_idx2str( idx_type idx ) const
{
 if( ( idx >= intFirstMOSEKPar ) && ( idx < intLastAlgParMOSEKS ) )
  return( SMSpp_to_MOSEK_int_pars[ idx - intFirstMOSEKPar ] );

 return( MILPSolver::int_par_idx2str( idx ) );
 }

/*--------------------------------------------------------------------------*/

Solver::idx_type MOSEKMILPSolver::dbl_par_str2idx(
					     const std::string & name ) const
{
 idx_type idx = MILPSolver::dbl_par_str2idx( name );
 if( idx < Inf< idx_type >() )
  return( idx );

 auto it = std::find( SMSpp_to_MOSEK_dbl_pars.begin() ,
		      SMSpp_to_MOSEK_dbl_pars.end() , name );
 if( it != SMSpp_to_MOSEK_dbl_pars.end() )
  return( MOSEK_to_SMSpp_dbl_pars[ std::distance(
                             SMSpp_to_MOSEK_dbl_pars.begin() , it ) ].second );

 return( Inf< idx_type >() );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

const std::string & MOSEKMILPSolver::dbl_par_idx2str( idx_type idx ) const
{
 if( ( idx >= dblFirstMOSEKPar ) && ( idx < dblLastAlgParMOSEKS ) )
  return( SMSpp_to_MOSEK_dbl_pars[ idx - dblFirstMOSEKPar ] );

 return( MILPSolver::dbl_par_idx2str( idx ) );
 }

/*--------------------------------------------------------------------------*/

Solver::idx_type MOSEKMILPSolver::str_par_str2idx(
					     const std::string & name ) const
{
 idx_type idx = MILPSolver::str_par_str2idx( name );
 if( idx < Inf< idx_type >() )
  return( idx );

 auto it = std::find( SMSpp_to_MOSEK_str_pars.begin() ,
		      SMSpp_to_MOSEK_str_pars.end() , name );
 if( it != SMSpp_to_MOSEK_str_pars.end() )
  return( MOSEK_to_SMSpp_str_pars[ std::distance(
                             SMSpp_to_MOSEK_str_pars.begin() , it ) ].second );

 return( Inf< idx_type >() );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

const std::string & MOSEKMILPSolver::str_par_idx2str( idx_type idx ) const
{
 if( ( idx >= strFirstMOSEKPar ) && ( idx < strLastAlgParMOSEKS ) )
  return( SMSpp_to_MOSEK_str_pars[ idx - strFirstMOSEKPar ] );

 return( MILPSolver::str_par_idx2str( idx ) );
 }

/*--------------------------------------------------------------------------*/

#ifdef MILPSolver_DEBUG

void MOSEKMILPSolver::check_status( void )
{
 MSKint32t n = 0;
 MSK_getnumvar( task , & n );
 if( numcols != n )
  DEBUG_LOG( "numcols is " << numcols << " but MSK_getnumvar returns "
	     << n << std::endl );

 MSK_getnumcon( task , & n );
 if( numrows != n )
  DEBUG_LOG( "numrows is " << numrows << " but MSK_getnumcon returns "
	     << n << std::endl );

 MILPSolver::check_status();
 }

#endif

/*--------------------------------------------------------------------------*/
/*-------------------- PRIVATE METHODS OF THE CLASS ------------------------*/
/*--------------------------------------------------------------------------*/

void MOSEKMILPSolver::check( MSKrescodee r , const char * where ) const
{
 // the warnings, below 1000, are no error
 if( ( r == MSK_RES_OK ) || ( r < 1000 ) )
  return;

 throw( std::runtime_error( std::string( "MOSEKMILPSolver::" ) + where +
			    ": MOSEK failed with " + code_text( r ) ) );
 }

/*--------------------------------------------------------------------------*/

MSKboundkeye MOSEKMILPSolver::bound_key( double lo , double up )
{
 const bool fl = ( lo > - Inf< double >() ) && ( ! is_inf( lo ) );
 const bool fu = ( up < Inf< double >() ) && ( ! is_inf( up ) );

 if( fl && fu )
  return( lo == up ? MSK_BK_FX : MSK_BK_RA );
 if( fl )
  return( MSK_BK_LO );
 if( fu )
  return( MSK_BK_UP );
 return( MSK_BK_FR );
 }

/*--------------------------------------------------------------------------*/

std::array< double , 2 > MOSEKMILPSolver::row_sides(
					     const FRowConstraint * con )
{
 return( std::array< double , 2 >{ double( con->get_lhs() ) ,
				   double( con->get_rhs() ) } );
 }

/*--------------------------------------------------------------------------*/

void MOSEKMILPSolver::put_var_bound( int j ,
				     const std::array< double , 2 > & bd )
{
 check( MSK_putvarbound( task , j , bound_key( bd[ 0 ] , bd[ 1 ] ) ,
			 msk_side( bd[ 0 ] ) , msk_side( bd[ 1 ] ) ) ,
	"put_var_bound" );
 }

/*--------------------------------------------------------------------------*/

void MOSEKMILPSolver::put_var_bounds(
			 const std::vector< int > & idx ,
			 const std::vector< std::array< double , 2 > > & bd )
{
 if( idx.empty() )
  return;

 std::vector< MSKint32t > sub( idx.begin() , idx.end() );
 std::vector< MSKboundkeye > bk( idx.size() );
 std::vector< double > bl( idx.size() ) , bu( idx.size() );
 for( std::size_t k = 0 ; k < idx.size() ; ++k ) {
  bl[ k ] = bd[ k ][ 0 ];
  bu[ k ] = bd[ k ][ 1 ];
  bk[ k ] = bound_key( bl[ k ] , bu[ k ] );
  bl[ k ] = msk_side( bl[ k ] );
  bu[ k ] = msk_side( bu[ k ] );
  }

 check( MSK_putvarboundlist( task , sub.size() , sub.data() , bk.data() ,
			     bl.data() , bu.data() ) , "put_var_bounds" );
 }

/*--------------------------------------------------------------------------*/

void MOSEKMILPSolver::put_con_bounds(
			 const std::vector< int > & idx ,
			 const std::vector< std::array< double , 2 > > & sd )
{
 if( idx.empty() )
  return;

 std::vector< MSKint32t > sub( idx.begin() , idx.end() );
 std::vector< MSKboundkeye > bk( idx.size() );
 std::vector< double > bl( idx.size() ) , bu( idx.size() );
 for( std::size_t k = 0 ; k < idx.size() ; ++k ) {
  bl[ k ] = sd[ k ][ 0 ];
  bu[ k ] = sd[ k ][ 1 ];
  bk[ k ] = bound_key( bl[ k ] , bu[ k ] );
  bl[ k ] = msk_side( bl[ k ] );
  bu[ k ] = msk_side( bu[ k ] );
  }

 check( MSK_putconboundlist( task , sub.size() , sub.data() , bk.data() ,
			     bl.data() , bu.data() ) , "put_con_bounds" );
 }

/*--------------------------------------------------------------------------*/

MSKvariabletypee MOSEKMILPSolver::var_type( const ColVariable * var ) const
{
 return( ( var->is_integer() && ( ! relax_int_vars ) ) ? MSK_VAR_TYPE_INT
	                                               : MSK_VAR_TYPE_CONT );
 }

/*--------------------------------------------------------------------------*/

void MOSEKMILPSolver::link_log( void )
{
 MSK_unlinkfuncfromtaskstream( task , MSK_STREAM_LOG );

 if( ! f_log_file.empty() )
  check( MSK_linkfiletotaskstream( task , MSK_STREAM_LOG ,
				   f_log_file.c_str() , 0 ) , "link_log" );
 else
  if( log_verbosity > 0 )
   check( MSK_linkfunctotaskstream( task , MSK_STREAM_LOG , nullptr ,
				    & MOSEKMILPSolver_print ) , "link_log" );
 }

/*--------------------------------------------------------------------------*/

void MOSEKMILPSolver::generate_qobj_hessian( void )
{
 q_hessian.assign( numcols , {} );
 q_hessian_nnz = 0;

 // the diagonal first, so that each column is born ordered, then the terms
 // x_i x_j, i != j, that MILPSolver keeps in the three parallel vectors;
 // the objective c' x + sum_j q_j x_j^2 + sum_k ndq_k x_i x_j is
 // c' x + 1/2 x' Q x with Q_jj = 2 q_j and Q_ij = Q_ji = ndq_k
 for( int j = 0 ; j < numcols ; ++j )
  if( ( j < int( q_objective.size() ) ) && ( q_objective[ j ] != 0 ) ) {
   q_hessian[ j ].emplace_back( j , 2 * q_objective[ j ] );
   ++q_hessian_nnz;
   }

 for( Index k = 0 ; k < ndq_objective.size() ; ++k )
  add_to_hessian( ndq_rowind[ k ] , ndq_colind[ k ] , ndq_objective[ k ] );

 f_hessian_changed = true;

 }  // end( MOSEKMILPSolver::generate_qobj_hessian )

/*--------------------------------------------------------------------------*/

void MOSEKMILPSolver::add_to_hessian( int i , int j , double delta )
{
 if( delta == 0 )
  return;

 if( i < j )
  std::swap( i , j );

 auto & col = q_hessian[ j ];
 auto it = std::lower_bound( col.begin() , col.end() , i ,
                             []( const auto & e , int r ) {
                              return( e.first < r );
                              } );

 if( ( it != col.end() ) && ( it->first == i ) ) {
  it->second += delta;
  if( it->second == 0 ) {  // the term is gone
   col.erase( it );
   --q_hessian_nnz;
   }
  }
 else {
  col.emplace( it , i , delta );
  ++q_hessian_nnz;
  }

 f_hessian_changed = true;

 }  // end( MOSEKMILPSolver::add_to_hessian )

/*--------------------------------------------------------------------------*/

void MOSEKMILPSolver::set_hessian( int i , int j , double value )
{
 if( i < j )
  std::swap( i , j );

 auto & col = q_hessian[ j ];
 auto it = std::lower_bound( col.begin() , col.end() , i ,
                             []( const auto & e , int r ) {
                              return( e.first < r );
                              } );

 if( ( it != col.end() ) && ( it->first == i ) ) {
  if( value == 0 ) {
   col.erase( it );
   --q_hessian_nnz;
   }
  else
   it->second = value;
  }
 else {
  if( value == 0 )
   return;
  col.emplace( it , i , value );
  ++q_hessian_nnz;
  }

 f_hessian_changed = true;

 }  // end( MOSEKMILPSolver::set_hessian )

/*--------------------------------------------------------------------------*/

void MOSEKMILPSolver::remove_hessian_column( int c )
{
 if( c >= int( q_hessian.size() ) )
  return;

 if( q_hessian_nnz == 0 ) {  // nothing to renumber
  q_hessian.erase( q_hessian.begin() + c );
  return;
  }

 // the terms of the column go with it, the ones of the row are found in the
 // columns before it, and every index past c moves back by one; MOSEK does
 // the same on its own copy, which therefore needs not be given again
 q_hessian_nnz -= q_hessian[ c ].size();
 q_hessian.erase( q_hessian.begin() + c );

 for( int j = 0 ; j < int( q_hessian.size() ) ; ++j ) {
  auto & col = q_hessian[ j ];
  if( col.empty() || ( col.back().first < c ) )
   continue;
  auto it = std::lower_bound( col.begin() , col.end() , c ,
                              []( const auto & e , int r ) {
                               return( e.first < r );
                               } );
  if( ( it != col.end() ) && ( it->first == c ) ) {
   it = col.erase( it );
   --q_hessian_nnz;
   }
  for( ; it != col.end() ; ++it )
   --( it->first );
  }

 }  // end( MOSEKMILPSolver::remove_hessian_column )

/*--------------------------------------------------------------------------*/

void MOSEKMILPSolver::pass_hessian( void )
{
 if( ! f_hessian_changed )
  return;

 f_hessian_changed = false;

 // MOSEK wants the triplets of the lower triangle, which is what q_hessian
 // holds column by column
 std::vector< MSKint32t > subi , subj;
 std::vector< double > val;
 subi.reserve( q_hessian_nnz );
 subj.reserve( q_hessian_nnz );
 val.reserve( q_hessian_nnz );

 for( int j = 0 ; j < int( q_hessian.size() ) ; ++j )
  for( const auto & e : q_hessian[ j ] ) {
   subi.push_back( e.first );
   subj.push_back( j );
   val.push_back( e.second );
   }

 check( MSK_putqobj( task , val.size() , subi.data() , subj.data() ,
		     val.data() ) , "pass_hessian" );

 }  // end( MOSEKMILPSolver::pass_hessian )

/*--------------------------------------------------------------------------*/

void MOSEKMILPSolver::change_obj_costs( std::vector< int > & idx ,
                                        std::vector< double > & val ,
                                        bool add , bool removal )
{
 if( idx.empty() )
  return;

 // a column named more than once takes the sum of its deltas if they are
 // added to the old cost, and the last value otherwise
 std::vector< Index > perm( idx.size() );
 std::iota( perm.begin() , perm.end() , 0 );
 std::stable_sort( perm.begin() , perm.end() ,
                   [ & idx ]( Index a , Index b ) {
                    return( idx[ a ] < idx[ b ] );
                    } );

 std::vector< MSKint32t > set;
 std::vector< double > nval;
 set.reserve( idx.size() );
 nval.reserve( idx.size() );
 for( auto p : perm )
  if( set.empty() || ( set.back() != idx[ p ] ) ) {
   set.push_back( idx[ p ] );
   nval.push_back( val[ p ] );
   }
  else
   if( add )
    nval.back() += val[ p ];
   else
    nval.back() = val[ p ];

 if( add ) {  // the old costs are read once for the whole set
  std::vector< double > old( set.size() );
  check( MSK_getclist( task , set.size() , set.data() , old.data() ) ,
	 "change_obj_costs" );
  for( Index k = 0 ; k < set.size() ; ++k ) {
   const double change = nval[ k ];
   double newval = old[ k ] + change;
   if( removal ) {
    const double scale = std::max( { 1.0 , std::abs( old[ k ] ) ,
                                    std::abs( change ) } );
    if( std::isfinite( newval ) &&
        ( std::abs( newval ) <= obj_removal_tol * scale ) )
     newval = 0.0;
    }
   nval[ k ] = newval;
   }
  }

 check( MSK_putclist( task , set.size() , set.data() , nval.data() ) ,
	"change_obj_costs" );

 }  // end( MOSEKMILPSolver::change_obj_costs )

/*--------------------------------------------------------------------------*/
/*---------------------- End File MOSEKMILPSolver.cpp ----------------------*/
/*--------------------------------------------------------------------------*/
