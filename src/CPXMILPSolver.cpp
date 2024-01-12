/*--------------------------------------------------------------------------*/
/*------------------------- File CPXMILPSolver.cpp -------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Implementation of the CPXMILPSolver class.
 *
 * \author Antonio Frangioni \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \author Niccolo' Iardella \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \author Kostas Tavlaridis-Gyparakis \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \copyright &copy; by Antonio Frangioni, Niccolo' Iardella
 */
/*--------------------------------------------------------------------------*/
/*---------------------------- IMPLEMENTATION ------------------------------*/
/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include <queue>

#include <LinearFunction.h>

#include <DQuadFunction.h>

#include "CPXMILPSolver.h"

#ifdef MILPSOLVER_DEBUG
 #define DEBUG_LOG( stuff ) std::cout << "[MILPSolver DEBUG] " << stuff
#else
 #define DEBUG_LOG( stuff )
#endif

// include the proper CPLEX parameter mapping
#include <boost/preprocessor/cat.hpp>
#include <boost/preprocessor/stringize.hpp>
#include BOOST_PP_STRINGIZE( BOOST_PP_CAT( BOOST_PP_CAT( CPX , CPX_VERSION ) , _maps.h ) )

/*--------------------------------------------------------------------------*/
/*------------------------- NAMESPACE AND USING ----------------------------*/
/*--------------------------------------------------------------------------*/

using namespace SMSpp_di_unipi_it;

/*--------------------------------------------------------------------------*/
/*-------------------------- FACTORY MANAGEMENT ----------------------------*/
/*--------------------------------------------------------------------------*/

SMSpp_insert_in_factory_cpp_0( CPXMILPSolver );

/*--------------------------------------------------------------------------*/
/*----------------------------- FUNCTIONS ----------------------------------*/
/*--------------------------------------------------------------------------*/

int CPXMILPSolver_callback( CPXCALLBACKCONTEXTptr context ,
			    CPXLONG contextid , void * userhandle )
{
 // just defer to the class method
 return( static_cast< CPXMILPSolver * >( userhandle
					 )->callback( context , contextid ) );
 }

/*--------------------------------------------------------------------------*/
/*--------------------- CONSTRUCTOR AND DESTRUCTOR -------------------------*/
/*--------------------------------------------------------------------------*/

CPXMILPSolver::CPXMILPSolver( void ) :
 MILPSolver() , env( nullptr ) , lp( nullptr ) , f_callback_set( false ) ,
 throw_reduced_cost_exception( 0 ) , CutSepPar( 0 ) ,
 UpCutOff( Inf< double >() ) , LwCutOff( -Inf< double >() )
{
 int status = 0;
 env = CPXopenCPLEX( & status );
 if( ! env )
  throw( std::runtime_error( "CPXopenCPLEX returned with status " +
			     std::to_string( status ) ) );
 lp = nullptr;
 #ifdef MILPSOLVER_DEBUG
  CPXsetintparam( env , CPXPARAM_Read_DataCheck , CPX_DATACHECK_WARN );
 #endif
 }

/*--------------------------------------------------------------------------*/

CPXMILPSolver::~CPXMILPSolver()
{
 for( auto el : v_ConfigDB )
  delete el;

 if( lp )
  CPXfreeprob( env , & lp );

 CPXcloseCPLEX( & env );
 }

/*--------------------------------------------------------------------------*/
/*--------------------- DERIVED METHODS OF BASE CLASS ----------------------*/
/*--------------------------------------------------------------------------*/

void CPXMILPSolver::set_Block( Block * block )
{
 if( block == f_Block )
  return;

 MILPSolver::set_Block( block );
 UpCutOff = Inf< double >();
 LwCutOff = -Inf< double >();
 }

/*--------------------------------------------------------------------------*/

void CPXMILPSolver::clear_problem( unsigned int what )
{
 MILPSolver::clear_problem( 0 );

 if( lp ) {
  CPXfreeprob( env , & lp );
  lp = nullptr;
  }
 }

/*--------------------------------------------------------------------------*/

void CPXMILPSolver::load_problem( void )
{
 MILPSolver::load_problem();

 int status = 0;
 if( lp )
  CPXfreeprob( env , & lp );
 lp = CPXcreateprob( env , & status , prob_name.c_str() );

 std::vector< double > cpx_lb = lb;
 std::vector< double > cpx_ub = ub;
 std::vector< double > cpx_rhs = rhs;

 for( int i = 0 ; i < numcols ; ++i ) {
  if( cpx_lb[ i ] == -Inf< double >() )
   cpx_lb[ i ] = -CPX_INFBOUND;
  if( cpx_ub[ i ] == Inf< double >() )
   cpx_ub[ i ] = CPX_INFBOUND;
  }

 for( int i = 0 ; i < numrows ; ++i ) {
  if( cpx_rhs[ i ] == -Inf< double >() )
   cpx_rhs[ i ] = -CPX_INFBOUND;
  else
   if( cpx_rhs[ i ] == Inf< double >() )
    cpx_rhs[ i ] = CPX_INFBOUND;
  }

 if( use_custom_names )
  CPXcopylpwnames( env , lp , numcols , numrows , objsense ,
                   objective.data(), cpx_rhs.data() , sense.data() ,
                   matbeg.data() , matcnt.data() , matind.data() ,
                   matval.data() , cpx_lb.data() , cpx_ub.data() ,
                   rngval.data() , colname.data() , rowname.data() );
 else
  CPXcopylp( env , lp , numcols , numrows , objsense ,
             objective.data() , cpx_rhs.data() , sense.data() ,
             matbeg.data() , matcnt.data() , matind.data() ,
             matval.data() , cpx_lb.data() , cpx_ub.data() ,
             rngval.data() );

 bool is_qp = std::any_of( q_objective.begin() ,
                           q_objective.end() ,
                           []( double d ) { return( d != 0 ); } );
 if( is_qp ) {
  // CPLEX evaluates the corresponding objective with a factor
  // of 0.5 in front of the quadratic objective term.
  std::vector< double > double_q_obj = q_objective;
  for( auto & qi : double_q_obj )
   qi *= 2;

  // adding q_objective information automatically changes the problem type
  // from linear to quadratic
  CPXcopyqpsep( env , lp , double_q_obj.data() );
  }

 // adding ctype information automatically changes the problem type
 // from continuous to mixed integer
 if( int_vars > 0 )
  CPXcopyctype( env , lp , xctype.data() );

 // the base representation isn't needed anymore
 MILPSolver::clear_problem( 15 );

 UpCutOff = Inf< double >();
 LwCutOff = -Inf< double >();

 }  // end( CPXMILPSolver::load_problem )

/*--------------------------------------------------------------------------*/

double CPXMILPSolver::get_problem_lb( const ColVariable & var ) const
{
 double b = MILPSolver::get_problem_lb( var );
 if( b == -Inf< double >() )
  b = -CPX_INFBOUND;

 return( b );
 }

/*--------------------------------------------------------------------------*/

double CPXMILPSolver::get_problem_ub( const ColVariable & var ) const
{
 double b = MILPSolver::get_problem_ub( var );
 if( b == Inf< double >() )
  b = CPX_INFBOUND;

 return( b );
 }

/*--------------------------------------------------------------------------*/

std::array< double , 2 > CPXMILPSolver::get_problem_bounds(
					      const ColVariable & var ) const
{
 auto ret = MILPSolver::get_problem_bounds( var );
 if( ret[ 0 ] == -Inf< double >() )
  ret[ 0 ] = -CPX_INFBOUND;
 if( ret[ 1 ] == Inf< double >() )
  ret[ 1 ] = CPX_INFBOUND;
 return( ret );
 }

/*--------------------------------------------------------------------------*/

int CPXMILPSolver::compute( bool changedvars )
{
 lock();  // lock the mutex: this is done again inside MILPSolver::compute,
          // but that's OK since the mutex is recursive

 // process Modification: this is driven by MILPSolver- - - - - - - - - - - -
 if( MILPSolver::compute( changedvars ) != kOK )
  throw( std::runtime_error( "an error occurred in MILPSolver::compute()" ) );

 // if required, write the problem to file- - - - - - - - - - - - - - - - - -
 if( ! output_file.empty() )
  CPXwriteprob( env , lp , output_file.c_str() , NULL );

 // figure out which API function is to be called - - - - - - - - - - - - - -
 bool is_qp = false;
 switch( CPXgetprobtype( env , lp ) ) {
  case( CPXPROB_LP ):
   // DEBUG_LOG( "CPLEX problem type: LP" << std::endl );
   // break;
  case( CPXPROB_MILP ):
   // DEBUG_LOG( "CPLEX problem type: MILP" << std::endl );
   break;
  case( CPXPROB_FIXEDMILP ):
   // DEBUG_LOG( "CPLEX problem type: FIXEDMILP" << std::endl );
   break;
  case( CPXPROB_QP ):
   // DEBUG_LOG( "CPLEX problem type: QP" << std::endl );
   // is_qp = true;
   // break;
  case( CPXPROB_MIQP ):
   // DEBUG_LOG( "CPLEX problem type: MIQP" << std::endl );
   // is_qp = true;
   // break;
  case( CPXPROB_FIXEDMIQP ):
   // DEBUG_LOG( "CPLEX problem type: FIXEDMIQP" << std::endl );
   is_qp = true;
   break;
  case( CPXPROB_QCP ):
   // DEBUG_LOG( "CPLEX problem type: QCP" << std::endl );
   throw( std::runtime_error( "Unsupported CPLEX problem type" ) );
  case( CPXPROB_MIQCP ):
   // DEBUG_LOG( "CPLEX problem type: MIQCP" << std::endl );
   throw( std::runtime_error( "Unsupported CPLEX problem type" ) );
  default:
   throw( std::runtime_error( "Undefined CPLEX problem type" ) );
  }

 // the actual call to CPLEX- - - - - - - - - - - - - - - - - - - - - - - - -

 if( int_vars > 0 ) {  // the MIP case- - - - - - - - - - - - - - - - - - - -

  if( ( CutSepPar & 7 ) ||
      ( UpCutOff < Inf< double >() ) || ( LwCutOff > Inf< double >() ) ) {
   // the callback has to be set
   CPXLONG cntxt = 0;
   if( ( UpCutOff < Inf< double >() ) || ( LwCutOff > Inf< double >() ) )
    cntxt = CPX_CALLBACKCONTEXT_LOCAL_PROGRESS
          | CPX_CALLBACKCONTEXT_GLOBAL_PROGRESS;
   if( CutSepPar & 3 )
    cntxt |= CPX_CALLBACKCONTEXT_RELAXATION;
   if( CutSepPar & 4 )
    cntxt |= CPX_CALLBACKCONTEXT_CANDIDATE;

   CPXcallbacksetfunc( env , lp , cntxt , & CPXMILPSolver_callback , this );
   f_callback_set = true;
   }
  else
   if( f_callback_set ) {    // the callback was set
    CPXcallbacksetfunc( env , lp , 0 , nullptr , nullptr );  // un-set it
    f_callback_set = false;
    }

  if( int status = CPXmipopt( env , lp ) ) {  // error
   if( status == CPXERR_SUBPROB_SOLVE )
    sol_status = decode_lqp_status( CPXgetsubstat( env , lp ) );
   else
    sol_status = decode_cpx_error( status );

   goto Return_status;
   }

  sol_status = decode_mip_status( CPXgetstat( env , lp ) );
  goto Return_status;
  }

 // the continuous case - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( int status = is_qp ? CPXqpopt( env , lp ) : CPXlpopt( env , lp ) ) {
  sol_status = decode_cpx_error( status );  // error
  goto Return_status;
  }

 sol_status = decode_lqp_status( CPXgetstat( env , lp ) );

 Return_status:
 unlock();  // unlock the mutex
 return( sol_status );

 }  // end( CPXMILPSolver::compute )

/*--------------------------------------------------------------------------*/

int CPXMILPSolver::decode_mip_status( int status )
{
 DEBUG_LOG( "CPXgetstat() returned " << status << std::endl );

 /* The following are the symbols that may represent the status of
 * a CPLEX solution as returned by CPXgetstat() in case of a a MIP,
 * as listed in CPLEX Callable Library API manual.
 * Some cases are commented out as they should never occur in the
 * conditions posed by CPXMILPSolver. */

 switch( status ) {
  case( CPXMIP_ABORT_FEAS ):
   // Stopped, but an integer solution exists.
  case( CPXMIP_ABORT_INFEAS ):
   // Stopped; no integer solution.
  case( CPXMIP_ABORT_RELAXATION_UNBOUNDED ):
   // Could not bound convex relaxation of nonconvex (MI)QP.
   return( kError );
  case( CPXMIP_DETTIME_LIM_FEAS ):
   // Deterministic time limit exceeded, but integer solution exists.
  case( CPXMIP_DETTIME_LIM_INFEAS ):
   // Deterministic time limit exceeded; no integer solution.
   return( kStopTime );
  case( CPXMIP_FAIL_FEAS ):
   // Terminated because of an error, but integer solution exists.
  case( CPXMIP_FAIL_FEAS_NO_TREE ):
   // Out of memory, no tree available, integer solution exists.
  case( CPXMIP_FAIL_INFEAS ):
   // Terminated because of an error; no integer solution.
  case( CPXMIP_FAIL_INFEAS_NO_TREE ):
   // Out of memory, no tree available, no integer solution.
   return( kError );
  case( CPXMIP_INFEASIBLE ):
   // Solution is integer infeasible.
   return( kInfeasible );
  case( CPXMIP_INForUNBD ):
   // Problem has been proven either infeasible or unbounded.
   //!! note: this is typically given by the preprocessor when it finds an
   //!!       easy dual unfeasible ray: since a dual feasible solution has
   //!!       not been constructed yet it is not formally possible to declare
   //!!       the problem unbounded, but both an empty primal and an empty
   //!!       dual is a rare occurrence, so the most likely correct answer is
   //!!       that the problem is unbounded
   //!!   return( kInfeasible );
   return( kUnbounded );
  case( CPXMIP_MEM_LIM_FEAS ):
   // Limit on tree memory has been reached, but an integer solution exists.
  case( CPXMIP_MEM_LIM_INFEAS ):
   // Limit on tree memory has been reached; no integer solution.
   return( kError );
  case( CPXMIP_NODE_LIM_FEAS ):
   // Node limit has been exceeded but integer solution exists.
  case( CPXMIP_NODE_LIM_INFEAS ):
   // Node limit has been reached; no integer solution.
   return( kStopIter );
  case( CPXMIP_OPTIMAL ):
   // An optimal integer solution has been found.
  case( CPXMIP_OPTIMAL_INFEAS ):
   // Problem is optimal with unscaled infeasibilities.
  case( CPXMIP_OPTIMAL_TOL ):
   // An optimal solution within the tolerance defined by
   // the relative or absolute MIP gap has been found.
  case( CPXMIP_SOL_LIM ):
   // The limit on mixed integer solutions has been reached.
   return( kOK );
  case( CPXMIP_TIME_LIM_FEAS ):
   // Time limit exceeded, but integer solution exists.
  case( CPXMIP_TIME_LIM_INFEAS ):
   // Time limit exceeded; no integer solution.
   return( kStopTime );
  case( CPXMIP_UNBOUNDED ):
   // Problem has an unbounded ray.
   return( kUnbounded );
   // case( CPXMIP_ABORT_RELAXED ):
   // case( CPXMIP_FEASIBLE ):
   // case( CPXMIP_FEASIBLE_RELAXED_INF ):
   // case( CPXMIP_FEASIBLE_RELAXED_QUAD ):
   // case( CPXMIP_FEASIBLE_RELAXED_SUM ):
   // case( CPXMIP_POPULATESOL_LIM ):
   // case( CPXMIP_OPTIMAL_POPULATED ):
   // case( CPXMIP_OPTIMAL_POPULATED_TOL ):
   // case( CPXMIP_OPTIMAL_RELAXED_INF ):
   // case( CPXMIP_OPTIMAL_RELAXED_QUAD ):
   // case( CPXMIP_OPTIMAL_RELAXED_SUM ):
  default:;
  }

 throw( std::runtime_error( "CPXgetstat() returned unknown status " +
			    std::to_string( status ) ) );
 }

/*--------------------------------------------------------------------------*/

int CPXMILPSolver::decode_lqp_status( int status )
{
 DEBUG_LOG( "CPXgetstat() returned " << status << std::endl );

 /* The following are the symbols that may represent the status of
  * a CPLEX solution as returned by CPXgetstat() in case of a LP/QP,
  * or by CPXgetsubstat() in case of a subproblem of a MIP,
  * as listed in CPLEX Callable Library API manual.
  * Some cases are commented out as they should never occur in the
  * conditions posed by CPXMILPSolver. */

 switch( status ) {
  case( CPX_STAT_ABORT_DETTIME_LIM ):
   // Stopped due to a deterministic time limit.
   return( kStopTime );
  case( CPX_STAT_ABORT_DUAL_OBJ_LIM ):
   // Stopped due to a limit on the dual objective.
   return( kError );
  case( CPX_STAT_ABORT_IT_LIM ):
   // Stopped due to limit on number of iterations.
   return( kStopIter );
  case( CPX_STAT_ABORT_OBJ_LIM ):
   // Stopped due to an objective limit.
  case( CPX_STAT_ABORT_PRIM_OBJ_LIM ):
   // Stopped due to a limit on the primal objective.
   return( kError );
  case( CPX_STAT_ABORT_TIME_LIM ):
   // Stopped due to a time limit.
   return( kStopTime );
  case( CPX_STAT_ABORT_USER ):
   // Stopped due to a request from the user.
   return( kError );
  case( CPX_STAT_BENDERS_NUM_BEST ):
   // Solution is infeasible, but cannot be cut with
   // a Benders cut due to numerical difficulties.
  case( CPX_STAT_INFEASIBLE ):
   // Problem has been proven infeasible.
   return( kInfeasible );
  case( CPX_STAT_INForUNBD ):
   // Problem has been proven either infeasible or unbounded.
   //!! note: this is typically given by the preprocessor when it finds an
   //!!       easy dual unfeasible ray: since a dual feasible solution has
   //!!       not been constructed yet it is not formally possible to declare
   //!!       the problem unbounded, but both an empty primal and an empty
   //!!       dual is a rare occurrence, so the most likely correct answer is
   //!!       that the problem is unbounded
   //!! return( kInfeasible );
   return( kUnbounded );
  case( CPX_STAT_NUM_BEST ):
   // Solution is available, but not proved optimal,
   // due to numeric difficulties during optimization.
  case( CPX_STAT_OPTIMAL ):
   // Optimal solution is available.
  case( CPX_STAT_OPTIMAL_FACE_UNBOUNDED ):
   // Model has an unbounded optimal face.
  case( CPX_STAT_OPTIMAL_INFEAS ):
   // Optimal solution is available, but with infeasibilities after unscaling.
   return( kOK );
  case( CPX_STAT_UNBOUNDED ):
   // Problem has an unbounded ray.
   return( kUnbounded );
   // case( CPX_STAT_CONFLICT_ABORT_CONTRADICTION ):
   // case( CPX_STAT_CONFLICT_ABORT_DETTIME_LIM ):
   // case( CPX_STAT_CONFLICT_ABORT_IT_LIM ):
   // case( CPX_STAT_CONFLICT_ABORT_MEM_LIM ):
   // case( CPX_STAT_CONFLICT_ABORT_NODE_LIM ):
   // case( CPX_STAT_CONFLICT_ABORT_OBJ_LIM ):
   // case( CPX_STAT_CONFLICT_ABORT_TIME_LIM ):
   // case( CPX_STAT_CONFLICT_ABORT_USER ):
   // case( CPX_STAT_CONFLICT_FEASIBLE ):
   // case( CPX_STAT_CONFLICT_MINIMAL ):
   // case( CPX_STAT_FEASIBLE ):
   // case( CPX_STAT_FEASIBLE_RELAXED_INF ):
   // case( CPX_STAT_FEASIBLE_RELAXED_QUAD ):
   // case( CPX_STAT_FEASIBLE_RELAXED_SUM ):
   // case( CPX_STAT_FIRSTORDER ):
   // case( CPX_STAT_MULTIOBJ_INFEASIBLE ):
   // case( CPX_STAT_MULTIOBJ_INForUNBD ):
   // case( CPX_STAT_MULTIOBJ_NON_OPTIMAL ):
   // case( CPX_STAT_MULTIOBJ_OPTIMAL ):
   // case( CPX_STAT_MULTIOBJ_STOPPED ):
   // case( CPX_STAT_MULTIOBJ_UNBOUNDED ):
   // case( CPX_STAT_OPTIMAL_RELAXED_INF ):
   // case( CPX_STAT_OPTIMAL_RELAXED_QUAD ):
   // case( CPX_STAT_OPTIMAL_RELAXED_SUM ):
  default:;
  }

 throw( std::runtime_error( "CPXgetstat() returned unknown status " +
                            std::to_string( status ) ) );
 }

/*--------------------------------------------------------------------------*/

int CPXMILPSolver::decode_cpx_error( int error )
{
 DEBUG_LOG( "CPLEX returned " << error << std::endl );

 /* The following symbols represent error codes returned by CPLEX, for
  * example by CPXlpopt(), CPXqpopt() and CPXmipopt(). Not all codes are
  * returned by all the methods, and it's hard to know what returns what
  * without doing extensive tests. So we are implementing just the ones we
  * encounter as we go. */
 switch( error ) {
  // case( CPXERR_ABORT_STRONGBRANCH ):
  // case( CPXERR_ADJ_SIGN_QUAD ):
  // case( CPXERR_ADJ_SIGNS ):
  // case( CPXERR_ADJ_SIGN_SENSE ):
  // case( CPXERR_ARC_INDEX_RANGE ):
  // case( CPXERR_ARRAY_BAD_SOS_TYPE ):
  // case( CPXERR_ARRAY_NOT_ASCENDING ):
  // case( CPXERR_ARRAY_TOO_LONG ):
  // case( CPXERR_BAD_ARGUMENT ):
  // case( CPXERR_BAD_BOUND_SENSE ):
  // case( CPXERR_BAD_BOUND_TYPE ):
  // case( CPXERR_BAD_CHAR ):
  // case( CPXERR_BAD_CTYPE ):
  // case( CPXERR_BAD_DECOMPOSITION ):
  // case( CPXERR_BAD_DIRECTION ):
  // case( CPXERR_BAD_EXPONENT ):
  // case( CPXERR_BAD_EXPO_RANGE ):
  // case( CPXERR_BAD_FILETYPE ):
  // case( CPXERR_BAD_ID ):
  // case( CPXERR_BAD_INDCONSTR ):
  // case( CPXERR_BAD_INDICATOR ):
  // case( CPXERR_BAD_INDTYPE ):
  // case( CPXERR_BAD_LAZY_UCUT ):
  // case( CPXERR_BAD_LUB ):
  // case( CPXERR_BAD_METHOD ):
  // case( CPXERR_BAD_MULTIOBJ_ATTR ):
  // case( CPXERR_BAD_NUMBER ):
  // case( CPXERR_BAD_OBJ_SENSE ):
  // case( CPXERR_BAD_PARAM_NAME ):
  // case( CPXERR_BAD_PARAM_NUM ):
  // case( CPXERR_BAD_PIVOT ):
  // case( CPXERR_BAD_PRIORITY ):
  // case( CPXERR_BAD_PROB_TYPE ):
  // case( CPXERR_BAD_ROW_ID ):
  // case( CPXERR_BAD_SECTION_BOUNDS ):
  // case( CPXERR_BAD_SECTION_ENDATA ):
  // case( CPXERR_BAD_SECTION_QMATRIX ):
  // case( CPXERR_BAD_SENSE ):
  // case( CPXERR_BAD_SOS_TYPE ):
  // case( CPXERR_BAD_STATUS ):
  // case( CPXERR_BAS_FILE_SHORT ):
  // case( CPXERR_BAS_FILE_SIZE ):
  // case( CPXERR_BENDERS_MASTER_SOLVE ):
  // case( CPXERR_CALLBACK ):
  // case( CPXERR_CALLBACK_INCONSISTENT ):
  // case( CPXERR_CAND_NOT_POINT ):
  // case( CPXERR_CAND_NOT_RAY ):
  // case( CPXERR_CNTRL_IN_NAME ):
  // case( CPXERR_COL_INDEX_RANGE ):
  // case( CPXERR_COL_REPEAT_PRINT ):
  // case( CPXERR_COL_REPEATS ):
  // case( CPXERR_COL_ROW_REPEATS ):
  // case( CPXERR_COL_UNKNOWN ):
  // case( CPXERR_CONFLICT_UNSTABLE ):
  // case( CPXERR_COUNT_OVERLAP ):
  // case( CPXERR_COUNT_RANGE ):
  // case( CPXERR_CPUBINDING_FAILURE ):
  // case( CPXERR_DBL_MAX ):
  // case( CPXERR_DECOMPRESSION ):
  // case( CPXERR_DETTILIM_STRONGBRANCH ):
  // case( CPXERR_DUP_ENTRY ):
  // case( CPXERR_DYNFUNC ):
  // case( CPXERR_DYNLOAD ):
  // case( CPXERR_ENCODING_CONVERSION ):
  // case( CPXERR_EXTRA_BV_BOUND ):
  // case( CPXERR_EXTRA_FR_BOUND ):
  // case( CPXERR_EXTRA_FX_BOUND ):
  // case( CPXERR_EXTRA_INTEND ):
  // case( CPXERR_EXTRA_INTORG ):
  // case( CPXERR_EXTRA_SOSEND ):
  // case( CPXERR_EXTRA_SOSORG ):
  // case( CPXERR_FAIL_OPEN_READ ):
  // case( CPXERR_FAIL_OPEN_WRITE ):
  // case( CPXERR_FILE_ENTRIES ):
  // case( CPXERR_FILE_FORMAT ):
  // case( CPXERR_FILE_IO ):
  // case( CPXERR_FILTER_VARIABLE_TYPE ):
  // case( CPXERR_ILL_DEFINED_PWL ):
  // case( CPXERR_INDEX_NOT_BASIC ):
  // case( CPXERR_INDEX_RANGE ):
  // case( CPXERR_INDEX_RANGE_HIGH ):
  // case( CPXERR_INDEX_RANGE_LOW ):
  // case( CPXERR_IN_INFOCALLBACK ):
  // case( CPXERR_INT_TOO_BIG ):
  // case( CPXERR_INT_TOO_BIG_INPUT ):
  // case( CPXERR_INVALID_NUMBER ):
  // case( CPXERR_LIMITS_TOO_BIG ):
  // case( CPXERR_LINE_TOO_LONG ):
  // case( CPXERR_LO_BOUND_REPEATS ):
  // case( CPXERR_LOCK_CREATE ):
  // case( CPXERR_LP_NOT_IN_ENVIRONMENT ):
  // case( CPXERR_LP_PARSE ):
  // case( CPXERR_MASTER_SOLVE ):
  // case( CPXERR_MIPSEARCH_WITH_CALLBACKS ):
  // case( CPXERR_MISS_SOS_TYPE ):
  // case( CPXERR_MSG_NO_CHANNEL ):
  // case( CPXERR_MSG_NO_FILEPTR ):
  // case( CPXERR_MSG_NO_FUNCTION ):
  // case( CPXERR_MULTIOBJ_SUBPROB_SOLVE ):
  // case( CPXERR_MULTIPLE_PROBS_IN_REMOTE_ENVIRONMENT ):
  // case( CPXERR_NAME_CREATION ):
  // case( CPXERR_NAME_NOT_FOUND ):
  // case( CPXERR_NAME_TOO_LONG ):
  // case( CPXERR_NAN ):
  // case( CPXERR_NEED_OPT_SOLN ):
  // case( CPXERR_NEGATIVE_SURPLUS ):
  // case( CPXERR_NET_DATA ):
  // case( CPXERR_NET_FILE_SHORT ):
  // case( CPXERR_NO_BARRIER_SOLN ):
  // case( CPXERR_NO_BASIC_SOLN ):
  // case( CPXERR_NO_BASIS ):
  // case( CPXERR_NO_BOUND_SENSE ):
  // case( CPXERR_NO_BOUND_TYPE ):
  // case( CPXERR_NO_COLUMNS_SECTION ):
  // case( CPXERR_NO_CONFLICT ):
  // case( CPXERR_NO_DECOMPOSITION ):
  // case( CPXERR_NO_OBJ_NAME ):
  // case( CPXERR_NODE_INDEX_RANGE ):
  // case( CPXERR_NODE_ON_DISK ):
  // case( CPXERR_NO_DUAL_SOLN ):
  // case( CPXERR_NO_ENDATA ):
  // case( CPXERR_NO_ENVIRONMENT ):
  // case( CPXERR_NO_FILENAME ):
  // case( CPXERR_NO_ID ):
  // case( CPXERR_NO_ID_FIRST ):
  // case( CPXERR_NO_INT_X ):
  // case( CPXERR_NO_KAPPASTATS ):
  // case( CPXERR_NO_LU_FACTOR ):
  // case( CPXERR_NO_MEMORY ):
  // case( CPXERR_NO_MIPSTART ):
  // case( CPXERR_NO_NAMES ):
  // case( CPXERR_NO_NAME_SECTION ):
  // case( CPXERR_NO_NORMS ):
  // case( CPXERR_NO_NUMBER_BOUND ):
  // case( CPXERR_NO_NUMBER ):
  // case( CPXERR_NO_NUMBER_FIRST ):
  // case( CPXERR_NO_OBJECTIVE ):
  // case( CPXERR_NO_OBJ_SENSE ):
  // case( CPXERR_NO_OPERATOR ):
  // case( CPXERR_NO_OP_OR_SENSE ):
  // case( CPXERR_NO_ORDER ):
  // case( CPXERR_NO_PROBLEM ):
  // case( CPXERR_NO_QP_OPERATOR ):
  // case( CPXERR_NO_QUAD_EXP ):
  // case( CPXERR_NO_RHS_COEFF ):
  // case( CPXERR_NO_RHS_IN_OBJ ):
  // case( CPXERR_NO_ROW_NAME ):
  // case( CPXERR_NO_ROW_SENSE ):
  // case( CPXERR_NO_ROWS_SECTION ):
  // case( CPXERR_NO_SENSIT ):
  // case( CPXERR_NO_SOLN ):
  // case( CPXERR_NO_SOLNPOOL ):
  // case( CPXERR_NO_SOS ):
  // case( CPXERR_NO_TREE ):
  // case( CPXERR_NOT_DUAL_UNBOUNDED ):
  // case( CPXERR_NOT_FIXED ):
  // case( CPXERR_NOT_FOR_BENDERS ):
  // case( CPXERR_NOT_FOR_DISTMIP ):
  // case( CPXERR_NOT_FOR_MIP ):
  // case( CPXERR_NOT_FOR_MULTIOBJ ):
  // case( CPXERR_NOT_FOR_QCP ):
  // case( CPXERR_NOT_FOR_QP ):
  // case( CPXERR_NOT_MILPCLASS ):
  // case( CPXERR_NOT_MIN_COST_FLOW ):
  // case( CPXERR_NOT_MIP ):
  // case( CPXERR_NOT_MIQPCLASS ):
  // case( CPXERR_NOT_ONE_PROBLEM ):
  // case( CPXERR_NOT_QP ):
  // case( CPXERR_NOT_SAV_FILE ):
  // case( CPXERR_NOT_UNBOUNDED ):
  // case( CPXERR_NO_VECTOR_SOLN ):
  // case( CPXERR_NULL_POINTER ):
  // case( CPXERR_ORDER_BAD_DIRECTION ):
  // case( CPXERR_OVERFLOW ):
  // case( CPXERR_PARAM_INCOMPATIBLE ):
  // case( CPXERR_PARAM_TOO_BIG ):
  // case( CPXERR_PARAM_TOO_SMALL ):
  // case( CPXERR_PRESLV_ABORT ):
  // case( CPXERR_PRESLV_BAD_PARAM ):
  // case( CPXERR_PRESLV_BASIS_MEM ):
  // case( CPXERR_PRESLV_COPYORDER ):
  // case( CPXERR_PRESLV_COPYSOS ):
  // case( CPXERR_PRESLV_CRUSHFORM ):
  // case( CPXERR_PRESLV_DETTIME_LIM ):
  // case( CPXERR_PRESLV_DUAL ):
  // case( CPXERR_PRESLV_FAIL_BASIS ):
  // case( CPXERR_PRESLV_INF ):
  // case( CPXERR_PRESLV_INForUNBD ):
  // case( CPXERR_PRESLV_NO_BASIS ):
  // case( CPXERR_PRESLV_NO_PROB ):
  // case( CPXERR_PRESLV_SOLN_MIP ):
  // case( CPXERR_PRESLV_SOLN_QP ):
  // case( CPXERR_PRESLV_START_LP ):
  // case( CPXERR_PRESLV_TIME_LIM ):
  // case( CPXERR_PRESLV_UNBD ):
  // case( CPXERR_PRESLV_UNCRUSHFORM ):
  // case( CPXERR_PRIIND ):
  // case( CPXERR_PRM_DATA ):
  // case( CPXERR_PROTOCOL ):
  // case( CPXERR_QCP_SENSE ):
  // case( CPXERR_QCP_SENSE_FILE ):
  // case( CPXERR_Q_DIVISOR ):
  // case( CPXERR_Q_DUP_ENTRY ):
  // case( CPXERR_Q_NOT_INDEF ):
  // case( CPXERR_Q_NOT_POS_DEF ):
  // case( CPXERR_Q_NOT_SYMMETRIC ):
  // case( CPXERR_QUAD_EXP_NOT_2 ):
  // case( CPXERR_QUAD_IN_ROW ):
  // case( CPXERR_RANGE_SECTION_ORDER ):
  // case( CPXERR_RESTRICTED_VERSION ):
  // case( CPXERR_RHS_IN_OBJ ):
  // case( CPXERR_RIMNZ_REPEATS ):
  // case( CPXERR_RIM_REPEATS ):
  // case( CPXERR_RIM_ROW_REPEATS ):
  // case( CPXERR_ROW_INDEX_RANGE ):
  // case( CPXERR_ROW_REPEAT_PRINT ):
  // case( CPXERR_ROW_REPEATS ):
  // case( CPXERR_ROW_UNKNOWN ):
  // case( CPXERR_SAV_FILE_DATA ):
  // case( CPXERR_SAV_FILE_VALUE ):
  // case( CPXERR_SAV_FILE_WRITE ):
  // case( CPXERR_SBASE_ILLEGAL ):
  // case( CPXERR_SBASE_INCOMPAT ):
  // case( CPXERR_SINGULAR ):
  // case( CPXERR_STR_PARAM_TOO_LONG ):
  // case( CPXERR_SUBPROB_SOLVE ):
   // CPXmipopt failed to solve one of the subproblems in the
   // branch-and-cut tree. This failure can be due to a limit
   // (for example, an iteration limit) or due to numeric trouble.
  // case( CPXERR_SYNCPRIM_CREATE ):
  // case( CPXERR_SYSCALL ):
  // case( CPXERR_THREAD_FAILED ):
  // case( CPXERR_TILIM_CONDITION_NO ):
  // case( CPXERR_TILIM_STRONGBRANCH ):
  // case( CPXERR_TOO_MANY_COEFFS ):
  // case( CPXERR_TOO_MANY_COLS ):
  // case( CPXERR_TOO_MANY_RIMNZ ):
  // case( CPXERR_TOO_MANY_RIMS ):
  // case( CPXERR_TOO_MANY_ROWS ):
  // case( CPXERR_TOO_MANY_THREADS ):
  // case( CPXERR_TREE_MEMORY_LIMIT ):
  // case( CPXERR_TUNE_MIXED ):
  // case( CPXERR_UNIQUE_WEIGHTS ):
  // case( CPXERR_UNSUPPORTED_CONSTRAINT_TYPE ):
  // case( CPXERR_UNSUPPORTED_OPERATION ):
  // case( CPXERR_UP_BOUND_REPEATS ):
  // case( CPXERR_WORK_FILE_OPEN ):
  // case( CPXERR_WORK_FILE_READ ):
  // case( CPXERR_WORK_FILE_WRITE ):
  // case( CPXERR_XMLPARSE ):
  default:
   return( kError + error );
  }
 }

/*--------------------------------------------------------------------------*/

Solver::OFValue CPXMILPSolver::get_lb( void )
{
 OFValue lower_bound = 0;
 int probtype = CPXgetprobtype( env , lp );

 switch( CPXgetobjsen( env, lp ) ) {
  case( CPX_MIN ):  // Minimization problem - - - - - - - - - - - - - - - - - -
   switch( sol_status ) {
    case( kUnbounded ):  lower_bound = -Inf< OFValue >(); break;
    case( kInfeasible ): lower_bound = Inf< OFValue >();  break;
    case( kOK ):
    case( kStopIter ):
    case( kStopTime ):
     switch( probtype ) {
      case( CPXPROB_MILP ):
      case( CPXPROB_MIQP ):
      case( CPXPROB_FIXEDMILP ):
      case( CPXPROB_FIXEDMIQP ):
       CPXgetbestobjval( env , lp , & lower_bound );
       lower_bound += constant_value;
       break;
      default:
       // FIXME: It's unclear how to get a lb for a continuous problem here
       CPXgetobjval( env , lp , & lower_bound );
       lower_bound += constant_value;
      }
     break;

    default:
     // If Cplex does not state that an optimal solution has been found
     // then we do not have enough information to provide a "good" lower bound
     // (the problem may be unbounded and Cplex has not detected it yet).
     // Therefore, in this case, the lower bound should be -Inf.
     lower_bound = -Inf< OFValue >();
     break;
    }
   break;

  case( CPX_MAX ):  // Maximization problem - - - - - - - - - - - - - - - - - -
   switch( sol_status ) {
    case( kUnbounded ):  lower_bound = Inf< OFValue >(); break;
    case( kInfeasible ): lower_bound = -Inf< OFValue >(); break;

    // if the algorithm has been stopped, the bound only exists if a
    // feasible solution has been generated
    case( kStopIter ):
    case( kStopTime ):
     if( ! has_var_solution() ) {
      lower_bound = -Inf< OFValue >();
      break;
      }
     lower_bound += constant_value;

    case( kOK ):
     CPXgetobjval( env , lp , & lower_bound );
     lower_bound += constant_value;
     break;

    default:
     // Same as above
     lower_bound = -Inf< OFValue >();
     break;
    }
   break;

  default:
   throw( std::runtime_error( "Objective type not yet defined" ) );
   break;
  }

 return( lower_bound );
 }

/*--------------------------------------------------------------------------*/

Solver::OFValue CPXMILPSolver::get_ub( void )
{
 OFValue upper_bound = 0;
 int probtype = CPXgetprobtype( env , lp );

 switch( CPXgetobjsen( env , lp ) ) {
  case( CPX_MIN ):  // Minimization problem - - - - - - - - - - - - - - - - - -
   switch( sol_status ) {
    case( kUnbounded ):  upper_bound = -Inf< OFValue >(); break;
    case( kInfeasible ): upper_bound = Inf< OFValue >(); break;

    // if the algorithm has been stopped, the bound only exists if a
    // feasible solution has been generated
    case( kStopIter ):
    case( kStopTime ):
     if( ! has_var_solution() ) {
      upper_bound = Inf< OFValue >();
      break;
      }
     upper_bound += constant_value;

    case( kOK ):
     CPXgetobjval( env , lp , & upper_bound );
     upper_bound += constant_value;
     break;

    default:
     // If Cplex does not state that an optimal solution has been found
     // then we do not have enough information to provide a "good" upper bound
     // (the problem may be unbounded and Cplex has not detected it yet).
     // Therefore, in this case, the upper bound should be +Inf.
     upper_bound = Inf< OFValue >();
     break;
    }
   break;

  case( CPX_MAX ):  // Maximization problem - - - - - - - - - - - - - - - - - -
   switch( sol_status ) {
    case( kUnbounded ):  upper_bound = Inf< OFValue >(); break;
    case( kInfeasible ): upper_bound = -Inf< OFValue >(); break;

    case( kOK ):
    case( kStopIter ):
    case( kStopTime ):
     switch( probtype ) {
      case( CPXPROB_MILP ):
      case( CPXPROB_MIQP ):
      case( CPXPROB_FIXEDMILP ):
      case( CPXPROB_FIXEDMIQP ):
       CPXgetbestobjval( env , lp , & upper_bound );
       upper_bound += constant_value;
       break;
      default:
       // FIXME: It's unclear how to get a ub for a continuous problem here
       CPXgetobjval( env , lp , & upper_bound );
       upper_bound += constant_value;
      }
     break;

    default:
     // Same as above
     upper_bound = Inf< OFValue >();
     break;
    }
   break;

   // Sense not defined
 default: throw( std::runtime_error( "Objective type not yet defined" ) );
 }

 return( upper_bound );
 }

/*--------------------------------------------------------------------------*/

bool CPXMILPSolver::has_var_solution( void )
{
 int solnmethod, solntype, pfeasind, dfeasind;
 int status = CPXsolninfo( env , lp , & solnmethod , & solntype ,
                           & pfeasind , & dfeasind );
 if( status )
  throw( std::runtime_error( "An error occurred in CPXsolninfo()" ) );

 switch( solntype ) {
  case( CPX_BASIC_SOLN ):    // The problem has a simplex basis
  case( CPX_NONBASIC_SOLN ): // Primal and dual solution but no basis
  case( CPX_PRIMAL_SOLN ):   // Primal solution but no dual solution
   return( true );
   // case( CPX_NO_SOLN ):   No solution
  }

 return( false );
 }

/*--------------------------------------------------------------------------*/

bool CPXMILPSolver::is_var_feasible( void )
{
 int solnmethod , solntype , pfeasind , dfeasind;
 if( CPXsolninfo( env , lp , & solnmethod , & solntype ,
		  & pfeasind , & dfeasind ) )
  throw( std::runtime_error( "an error occurred in CPXsolninfo()" ) );

 return( bool( pfeasind ) );
 }

/*--------------------------------------------------------------------------*/

Solver::OFValue CPXMILPSolver::get_var_value( void )
{
 switch( CPXgetobjsen( env , lp ) ) {
  case( CPX_MIN ): return( get_ub() );
  case( CPX_MAX ): return( get_lb() );
  default: throw( std::runtime_error( "Objective type not yet defined" ) );
  }
 }

/*--------------------------------------------------------------------------*/

void CPXMILPSolver::get_var_solution( Configuration * solc )
{
 std::vector< double > x( numcols , 0 );
 if( CPXgetx( env , lp , x.data() , 0 , numcols - 1 ) )
  throw( std::runtime_error( "Unable to get the solution with CPXgetx()" ) );

 get_var_solution( x );
 }

/*--------------------------------------------------------------------------*/

void CPXMILPSolver::get_var_solution( const std::vector< double > & x )
{
 int col = 0;
 int dcol = static_vars;

 auto set = [ & x , & col ]( ColVariable & v ) {
  v.set_value( x[ col++ ] );
  };

 auto setd = [ & x , & dcol ]( ColVariable & v ) {
  v.set_value( x[ dcol++ ] );
  };

 for( auto qb : v_BFS ) {
  for( const auto & vi : qb->get_static_variables() )
   un_any_const_static( vi , set , un_any_type< ColVariable >() );

  for( const auto & vi : qb->get_dynamic_variables() )
   un_any_const_dynamic( vi , setd , un_any_type< ColVariable >() );
  }
 }

/*--------------------------------------------------------------------------*/

bool CPXMILPSolver::has_dual_solution( void )
{
 int solnmethod , solntype , pfeasind , dfeasind;
 if( CPXsolninfo( env , lp , & solnmethod , & solntype ,
		  & pfeasind , & dfeasind ) )
  throw( std::runtime_error( "An error occurred in CPXsolninfo()" ) );

 switch( solntype ) {
  case( CPX_BASIC_SOLN ):    // The problem has a simplex basis
  case( CPX_NONBASIC_SOLN ): // Primal and dual solution but no basis
   return( true );
   // case( CPX_PRIMAL_SOLN ):   Primal solution but no dual solution
   // case( CPX_NO_SOLN ):       No solution
  }

 return( false );
 }

/*--------------------------------------------------------------------------*/

bool CPXMILPSolver::is_dual_feasible( void )
{
 int solnmethod , solntype , pfeasind , dfeasind;
 if( CPXsolninfo( env , lp , & solnmethod , & solntype ,
		  & pfeasind , & dfeasind ) )
  throw( std::runtime_error( "an error occurred in CPXsolninfo()" ) );

 return( bool( dfeasind ) );
 }

/*--------------------------------------------------------------------------*/

void CPXMILPSolver::get_dual_solution( Configuration * solc )
{
 std::vector< double > pi( numrows , 0 );
 std::vector< double > dj( numcols , 0 );

 if( numrows > 0 )
  if( CPXgetpi( env , lp , pi.data() , 0 , numrows - 1 ) )
   throw( std::runtime_error( "Unable to get dual values with CPXgetpi()" ) );

 if( CPXgetdj( env , lp , dj.data() , 0 , numcols - 1 ) )
  throw( std::runtime_error( "Unable to get reduced costs with CPXgetdj()"
			     ) );
 int row = 0;
 int row_dynamic = static_cons;

 auto set = [ & pi , & row ]( FRowConstraint & c ) {
  c.set_dual( - pi[ row++ ] );
  };

 auto set_dynamic = [ & pi , & row_dynamic ]( FRowConstraint & c ) {
  c.set_dual( - pi[ row_dynamic++ ] );
  };

 for( auto qb : v_BFS ) {
  for( const auto & ci : qb->get_static_constraints() )
   un_any_const_static( ci , set , un_any_type< FRowConstraint >() );

  for( const auto & ci : qb->get_dynamic_constraints() )
   un_any_const_dynamic( ci, set_dynamic, un_any_type< FRowConstraint >() );
  }

 for( int i = 0 ; i < numcols ; ++i ) {

  auto var = variable_with_index( i );
  auto active_bounds = get_active_bounds( *var );

  // Bounds that will have the dual value set.
  OneVarConstraint * lhs_con = nullptr;
  OneVarConstraint * rhs_con = nullptr;

  auto var_lb = var->get_lb();
  auto var_ub = var->get_ub();

  const auto var_is_fixed = var->is_fixed();
  if( var_is_fixed ) {
   /* The Variable is fixed. There should be at least one OneVarConstraint
    * (for this Variable) whose lower and upper bounds are equal to the value
    * of this Variable. If such OneVarConstraint exists, the reduced cost of
    * this Variable will be dual of that OneVarConstraint. If there is no such
    * OneVarConstraint, the reduced cost of this variable will be lost. */
   var_lb = var->get_value();
   var_ub = var->get_value();

   for( auto b: active_bounds ) {
    b->set_dual( 0 );
    if( b->get_lhs() == var_lb && b->get_rhs() == var_lb ) {
     lhs_con = b;
     rhs_con = b;
     }
    }

   assert( lhs_con == rhs_con );
   }
  else {  // a non-fixed Variable
   for( auto b: active_bounds ) {
    b->set_dual( 0 );

    if( b->get_lhs() >= var_lb ) {
     var_lb = b->get_lhs();
     lhs_con = b;
     }

    if( b->get_rhs() <= var_ub ) {
     var_ub = b->get_rhs();
     rhs_con = b;
     }
    }
   }

  if( lhs_con && ( dj[ i ] >= 0 ) )
   lhs_con->set_dual( - dj[ i ] );
  else
   if( rhs_con && ( dj[ i ] <= 0 ) )
    rhs_con->set_dual( - dj[ i ] );
   else
    if( lhs_con || rhs_con )
     throw( std::logic_error(
	       "CPXMILPSolver::get_dual_solution: invalid dual value." ) );

  if( throw_reduced_cost_exception ) {
   if( var_is_fixed && ( ! lhs_con ) && ( var_lb != 0 ) ) {
    /* The Variable is fixed but it has no associated OneVarConstraint
     * with both bounds equal to the value of the Variable. */

    throw( std::logic_error(
     "CPXMILPSolver::get_dual_solution: variable with index " +
     std::to_string( i ) + " is fixed to " +
     std::to_string( var->get_value() ) + ", but it has no OneVarConstraint" +
     "with both bounds equal to the value of this variable." ) );
    }
   else
    if( ( ! var_is_fixed ) && ( ! lhs_con ) && ( ! rhs_con ) ) {
     /* The Variable is not fixed and it has no associated OneVarConstraint.
      * An exception is thrown if it has a finite nonzero bound. */

     if( ( ( var_lb != 0 ) && ( std::abs( var_lb ) < Inf< double >() ) ) ||
	 ( ( var_ub != 0 ) && ( std::abs( var_ub ) < Inf< double >() ) ) )
      throw( std::logic_error(
                "CPXMILPSolver::get_dual_solution: variable with index " +
		std::to_string( i ) + " has no OneVarConstraint." ) );
     }
   }
  }
 }  // end( CPXMILPSolver::get_dual_solution )

/*--------------------------------------------------------------------------*/

bool CPXMILPSolver::has_dual_direction( void )
{
 std::vector< double > y( numrows , 0 );
 double proof = 0;
 return( ! bool( CPXdualfarkas( env , lp , y.data() , & proof ) ) );
 }

/*--------------------------------------------------------------------------*/

void CPXMILPSolver::get_dual_direction( Configuration * dirc )
{
 std::vector< double > y( numrows , 0 );
 std::vector< double > dj( numcols , 0 );

 // CPXdualfarkas gives a Farkas certificate y so that:
 // y' * A * x >= y' * b
 //   If it is a <= constraint then y[i] <= 0 holds;
 //   If it is a >= constraint then y[i] >= 0 holds.

 double proof = 0;
 if( CPXdualfarkas( env , lp , y.data() , & proof ) )
  throw( std::runtime_error( "an error occurred in CPXdualfarkas()" ) );

 // CPXdjfrompi computes reduced costs from dual values
 // dj = c - A'y

 if( CPXdjfrompi( env , lp , y.data() , dj.data() ) )
  throw( std::runtime_error( "an error occurred in CPXdjfrompi()" ) );

 int row = 0;
 int row_dynamic = static_cons;
 auto set = [ & y , & row ]( FRowConstraint & c ) {
  c.set_dual( - y[ row++ ] );
  };

 auto set_dynamic = [ & y , & row_dynamic ]( FRowConstraint & c ) {
  c.set_dual( - y[ row_dynamic++ ] );
  };

 for( auto qb : v_BFS ) {
  for( const auto & ci : qb->get_static_constraints() )
   un_any_const_static( ci , set , un_any_type< FRowConstraint >() );

  for( const auto & ci : qb->get_dynamic_constraints() )
   un_any_const_dynamic( ci , set_dynamic , un_any_type< FRowConstraint >() );
  }

 for( int i = 0 ; i < numcols ; ++i ) {
  // Bounds that will have the dual value set.
  OneVarConstraint * lhs_con = nullptr;
  OneVarConstraint * rhs_con = nullptr;

  auto var = variable_with_index( i );
  auto active_bounds = get_active_bounds( *var );

  double var_lb = var->get_lb();
  double var_ub = var->get_ub();

  const auto var_is_fixed = var->is_fixed();
  if( var_is_fixed ) {
   /* The Variable is fixed. There should be at least one OneVarConstraint
    * (for this Variable) whose lower and upper bounds are equal to the value
    * of this Variable. If such a OneVarConstraint exists, the reduced cost of
    * this Variable will be dual of that OneVarConstraint. If there is no such
    * OneVarConstraint, the reduced cost of this variable will be lost. */
   var_lb = var->get_value();
   var_ub = var->get_value();

   for( auto b: active_bounds ) {
    b->set_dual( 0 );
    if( ( b->get_lhs() == var_lb ) && ( b->get_rhs() == var_lb ) ) {
     lhs_con = b;
     rhs_con = b;
     }
    }

   assert( lhs_con == rhs_con );
   }
  else {  // A non-fixed Variable
   for( auto b : active_bounds ) {
    b->set_dual( 0 );

    if( b->get_lhs() >= var_lb ) {
     var_lb = b->get_lhs();
     lhs_con = b;
     }

    if( b->get_rhs() <= var_ub ) {
     var_ub = b->get_rhs();
     rhs_con = b;
     }
    }
   }

  if( lhs_con && ( dj[ i ] >= 0 ) )
   lhs_con->set_dual( - dj[ i ] );
  else
   if( rhs_con && ( dj[ i ] <= 0 ) )
    rhs_con->set_dual( - dj[ i ] );
   else
    if( lhs_con || rhs_con )
     throw( std::logic_error(
	       "CPXMILPSolver::get_dual_direction: invalid dual value" ) );

  if( throw_reduced_cost_exception ) {
   if( var_is_fixed && ( ! lhs_con ) && ( var_lb != 0 ) ) {
    /* The Variable is fixed but it has no associated OneVarConstraint
     * with both bounds equal to the value of the Variable. */

    throw( std::logic_error(
     "CPXMILPSolver::get_dual_direction: variable with index " +
     std::to_string( i ) + " is fixed to " +
     std::to_string( var->get_value() ) + ", but it has no OneVarConstraint" +
     "with both bounds equal to the value of this variable." ) );
    }
   else
    if( ( ! var_is_fixed ) && ( ! lhs_con ) && ( ! rhs_con ) ) {
     /* The Variable is not fixed and it has no associated OneVarConstraint.
      * An exception is thrown if it has a finite nonzero bound. */

     if( ( ( var_lb != 0 ) && ( std::abs( var_lb ) < Inf< double >() ) ) ||
	 ( ( var_ub != 0 ) && ( std::abs( var_ub ) < Inf< double >() ) ) )
      throw( std::logic_error(
                "CPXMILPSolver::get_dual_direction: variable with index " +
		std::to_string( i ) + " has no OneVarConstraint." ) );
     }
   }
  }
 }  // end( CPXMILPSolver::get_dual_direction )

/*--------------------------------------------------------------------------*/

void CPXMILPSolver::write_lp( const std::string & filename )
{
 CPXwriteprob( env , lp , filename.c_str() , "LP" );
 }

/*--------------------------------------------------------------------------*/

int CPXMILPSolver::get_nodes( void ) const
{
 return( CPXgetnodecnt( env , lp ) );
 }

/*--------------------------------------------------------------------------*/
/*-------------------- PROTECTED FIELDS OF THE CLASS -----------------------*/
/*--------------------------------------------------------------------------*/

void CPXMILPSolver::var_modification( const VariableMod * mod )
{
 bool was_mip = int_vars > 0;  // was a MIP before the change

 // call the method of MILPSolver to update dictionaries (and int_vars)
 MILPSolver::var_modification( mod );

 bool is_mip = int_vars > 0;  // is a MIP after the change

 auto var = static_cast< const ColVariable * >( mod->variable() );
 auto idx = index_of_variable( var );
 if( idx == Inf< int >() )  // the Variable is not (yet) there (?)
  return;                   // nothing to do
 std::array< int , 2 > indices = { idx , idx };

 // react to changes in the integrality - - - - - - - - - - - - - - - - - - -
 if( ColVariable::is_integer( mod->old_state() ) !=
     ColVariable::is_integer( mod->new_state() ) ) {

  // construct new variable type
  char new_ctype;
  if( var->is_integer() && ( ! relax_int_vars ) ) {
   if( var->is_unitary() && var->is_positive() )
    new_ctype = 'B'; // Binary
   else
    new_ctype = 'I'; // Integer
   }
  else
   new_ctype = 'C';  // Continuous

  if( is_mip )    // it is a MIP
   if( was_mip )  // but it already was so, update only the one variable
    CPXchgctype( env , lp , 1 , indices.data() , & new_ctype );
   else {         // it was not a MIP
    // the first integer variable is added
    // all ctype values must be [re]added to the problem
    switch( CPXgetprobtype( env, lp ) ) {
     case( CPXPROB_LP ): CPXchgprobtype( env , lp , CPXPROB_MILP ); break;
     case( CPXPROB_QP ): CPXchgprobtype( env , lp , CPXPROB_MIQP ); break;
     default: throw( std::runtime_error( "Wrong CPLEX problem type" ) );
     }

    // all variables are continuous
    std::vector< char > ctype( numcols , 'C' );
    ctype[ idx ] = new_ctype;  // save the one that just became binary
    CPXcopyctype( env , lp , ctype.data() );
    }
  else              // it is not a MIP
   if( was_mip ) {  // but was so
    // the last integer variable was removed, remove all ctype values
    switch( CPXgetprobtype( env , lp ) ) {
     case( CPXPROB_LP ):
     case( CPXPROB_QP ):   break;
     case( CPXPROB_MILP ):
     case( CPXPROB_FIXEDMILP ): CPXchgprobtype( env , lp , CPXPROB_LP ); break;
     case( CPXPROB_MIQP ):
     case( CPXPROB_FIXEDMIQP ): CPXchgprobtype( env , lp , CPXPROB_QP ); break;
     default: throw( std::runtime_error( "Wrong CPLEX problem type" ) );
     }
    // else do nothing, the problem stays continuous; besides this should
    // never happen as the variable type has changed
    }
  }   // end( if( new integrality != old integrality ) )

 // react to fix / unfix- - - - - - - - - - - - - - - - - - - - - - - - - - -
 if( Variable::is_fixed( mod->old_state() ) !=
     Variable::is_fixed( mod->new_state() ) ) {
  if( Variable::is_fixed( mod->new_state() ) ) {  // fix the variable
    static std::array< char , 1 > lu = { 'B' };
    std::array< double , 1 > bd = { var->get_value() };
    CPXchgbds( env , lp , 1 , indices.data() , lu.data() , bd.data() );
    }
   else {                                         // un fix the variable
    static std::array< char , 2 > lu = { 'L' , 'U' };
    auto bd = CPXMILPSolver::get_problem_bounds( *var );
    CPXchgbds( env , lp , 2 , indices.data() , lu.data() , bd.data() );
    }
  }
 }  // end( CPXMILPSolver::var_modification )

/*--------------------------------------------------------------------------*/

void CPXMILPSolver::objective_modification( const ObjectiveMod * mod )
{
 // call the method of MILPSolver to update sense direction
 MILPSolver::objective_modification( mod );

 /* ObjectiveMod class does not include any modification types except
  * for eSetMin and eSetMax.
  * To change OF coefficents, a FunctionMod must be used. */

 switch( mod->type() ) {
  case( ObjectiveMod::eSetMin ): CPXchgobjsen( env , lp , CPX_MIN ); break;
  case( ObjectiveMod::eSetMax ): CPXchgobjsen( env , lp , CPX_MAX ); break;
  default: throw( std::invalid_argument( "Invalid type of ObjectiveMod" ) );
  }
 }

/*--------------------------------------------------------------------------*/

void CPXMILPSolver::const_modification( const ConstraintMod * mod )
{
 // no point in calling the method of MILPSolver, as it does nothing
 // MILPSolver::const_modification( mod );

 /* To change the coefficents, a FunctionMod must be used. */

 auto * con = dynamic_cast< FRowConstraint * >( mod->constraint() );
 if( ! con )  // this should not happen
  return;     // but in case, nothing to do

 int index = index_of_constraint( con );
 if( index == Inf< int >() )  // the FRowConstraint is not (yet?) there
  return;                     // nothing to do
 char sense;
 double rhs;
 double rngval;

 RowConstraint::RHSValue con_lhs = NAN;
 RowConstraint::RHSValue con_rhs = NAN;

 switch( mod->type() ) {
  case( ConstraintMod::eRelaxConst ):
   // In order to relax the constraint all we do is transform it
   // into an inequality (<=) with RHS equal to infinity

   sense = 'L';
   rhs = CPX_INFBOUND;
   CPXchgsense( env , lp , 1 , & index , & sense );
   CPXchgrhs( env , lp , 1 , & index , & rhs );
   break;

  case( ConstraintMod::eEnforceConst ):
  case( RowConstraintMod::eChgLHS ):
  case( RowConstraintMod::eChgRHS ):
  case( RowConstraintMod::eChgBTS ):
   // In order to enforce a relaxed constraint all we need to do is
   // reverse the process of relaxing it, by changing the sense and
   // the rhs back to the original form of the constraint
   // Moreover, for the way the LP vectors are built, handling
   // LHS/RHS/BTS cases separately is not worth it.

   con_lhs = con->get_lhs();
   con_rhs = con->get_rhs();

   if( con_lhs == con_rhs ) {
    sense = 'E';
    rhs = con_rhs;
    }
   else
    if( con_lhs == -Inf< double >() ) {
     sense = 'L';
     rhs = con_rhs;
     }
    else
     if( con_rhs == Inf< double >() ) {
      sense = 'G';
      rhs = con_lhs;
      }
     else {
      sense = 'R';
      rhs = con_lhs;
      rngval = con_rhs - con_lhs;
      }

   CPXchgrhs( env , lp , 1 , & index , & rhs );
   CPXchgsense( env , lp , 1 , & index , & sense );
   if( sense == 'R' )
    CPXchgrngval( env , lp , 1 , & index , & rngval );
   break;

  default:
   throw( std::invalid_argument( "Invalid type of ConstraintMod" ) );
  }
 }  // end( CPXMILPSolver::const_modification )

/*--------------------------------------------------------------------------*/

void CPXMILPSolver::bound_modification( const OneVarConstraintMod * mod )
{
 // no point in calling the method of MILPSolver, as it does nothing
 // MILPSolver::bound_modification( mod );

 /* One ColVariable can have more active OneVarConstraints, so each time we
  * modify one of them we have to check if LHS and RHS of the Variable
  * actually change (as they may not). */

 static std::array< char , 2 > lu = { 'L' , 'U' };

 auto con = static_cast< OneVarConstraint * >( mod->constraint() );
 auto var = static_cast< const ColVariable * >( con->get_active_var( 0 ) );
 if( ! var )  // this should never happen
  return;     // but in case, there is nothing to do

 // fixed variables are implemented in CPXMILPSolver by changing the bounds;
 // therefore, actual changes of the bounds are ignored here. note that we
 // are assuming the new bounds do not make the fixed value of the variable
 // unfeasible, as this will make the whole problem unfeasible
 // TODO: check this
 if( var->is_fixed() )
  return;

 auto vi = index_of_variable( var );
 if( vi == Inf< int >() )  // the ColVariable has been removed
  return;                  // is strange, but there is nothing to do

 std::array< int , 2 > ind = { vi , vi };

 switch( mod->type() ) {

  case( RowConstraintMod::eChgLHS ): {
   std::array< double , 1 > bd = { CPXMILPSolver::get_problem_lb( *var ) };
   CPXchgbds( env , lp , 1 , ind.data() , lu.data() , bd.data() );
   break;
   }

  case( RowConstraintMod::eChgRHS ): {
   std::array< double , 1 > bd = { CPXMILPSolver::get_problem_ub( *var ) };
   CPXchgbds( env , lp , 1 , ind.data() , lu.data() + 1 , bd.data() );
   break;
   }

  case( RowConstraintMod::eChgBTS ): {
   auto bd = CPXMILPSolver::get_problem_bounds( *var );
   CPXchgbds( env , lp , 2 , ind.data() , lu.data() , bd.data() );
   break;
   }

  default:
   throw( std::invalid_argument( "Invalid type of OneVarConstraintMod" ) );
  }
 }  // end( CPXMILPSolver::bound_modification )

/*--------------------------------------------------------------------------*/

void CPXMILPSolver::objective_function_modification( const FunctionMod * mod )
{
 // no point in calling the method of MILPSolver, as it does nothing
 // MILPSolver::objective_function_modification( mod );

 auto f = mod->function();

 // C05FunctionModLin - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 if( auto modl = dynamic_cast< const C05FunctionModLin * >( mod ) ) {

  if( auto lf = dynamic_cast< const LinearFunction * >( f ) ) {
   // Linear objective function

   Subset idxs;
   if( auto modlr = dynamic_cast< const C05FunctionModLinRngd * >( modl ) )
    idxs = lf->map_index( modl->vars() , modlr->range() );
   else
    if( auto modls = dynamic_cast< const C05FunctionModLinSbst * >( modl ) )
     idxs = lf->map_index( modl->vars() , modls->subset() );
    else
     throw( std::logic_error( "unknown type of C05FunctionModLinRngd" ) );

   std::vector< double > nval( idxs.size() );
   std::vector< int > cidx( idxs.size() );
   auto nvit = nval.begin();
   auto idxit = idxs.begin();
   auto cidxit = cidx.begin();
   auto & cp = lf->get_v_var();

   for( auto v :  modl->vars() )
    if( auto idx = *(idxit++) ; idx < Inf< Index >() ) {
     *(nvit++) = cp[ idx ].second;
     *(cidxit++) = index_of_variable( static_cast< const ColVariable * >( v ) );
     }

   auto nsz = std::distance( nval.begin() , nvit );
   cidx.resize( nsz );
   nval.resize( nsz );

   CPXchgobj( env , lp , cidx.size() , cidx.data() , nval.data() );
   return;
   }

  if( auto qf = dynamic_cast< const DQuadFunction * >( f ) ) {
   // quadratic objective function

   Subset idxs;
   if( auto modlr = dynamic_cast< const C05FunctionModLinRngd * >( modl ) )
    idxs = qf->map_index( modl->vars() , modlr->range() );
   else
    if( auto modls = dynamic_cast< const C05FunctionModLinSbst * >( modl ) )
     idxs = qf->map_index( modl->vars() , modls->subset() );
    else
     throw( std::logic_error( "unknown type of C05FunctionModLinRngd" ) );

   std::vector< double > nval( idxs.size() );
   std::vector< int > cidx( idxs.size() );
   auto nvit = nval.begin();
   auto idxit = idxs.begin();
   auto cidxit = cidx.begin();
   auto & cp = qf->get_v_var();

   for( auto v :  modl->vars() )
    if( auto idx = *(idxit++) ; idx < Inf< Index >() ) {
     *(nvit++) = std::get< 1 >( cp[ idx ] );
     *(cidxit++) = index_of_variable( static_cast< const ColVariable * >( v ) );
     }

   auto nsz = std::distance( nval.begin() , nvit );
   cidx.resize( nsz );
   nval.resize( nsz );

   CPXchgobj( env , lp , cidx.size() , cidx.data() , nval.data() );
   return;
   }

  // This should never happen
  throw( std::invalid_argument( "Unknown type of Objective Function" ) );
  }

 // C05FunctionMod- - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 if( auto modl = dynamic_cast< const C05FunctionMod * >( mod ) ) {

  if( modl->type() == C05FunctionMod::NothingChanged ) {

   const auto shift = modl->shift();

   if( ( shift == FunctionMod::INFshift ) ||
       ( shift == - FunctionMod::INFshift ) ||
       ( std::isnan( shift ) ) )
    throw( std::logic_error(
     "unexpected *C05FunctionMod* from Objective Function" ) );

   constant_value += shift;
   return;
   }

  auto qf = dynamic_cast< const DQuadFunction * >( f );
  if( ! qf )
   throw( std::logic_error(
		       "unexpected *C05FunctionMod* from Linear Objective" ) );

  Subset idxs;
  c_Vec_p_Var * vars;
  if( auto modlr = dynamic_cast< const C05FunctionModRngd * >( modl ) ) {
   idxs = qf->map_index( modlr->vars() , modlr->range() );
   vars = & modlr->vars();
   }
  else
   if( auto modls = dynamic_cast< const C05FunctionModSbst * >( modl ) ) {
    idxs = qf->map_index( modls->vars() , modls->subset() );
    vars = & modls->vars();
    }
   else
    throw( std::logic_error( "unknown type of C05FunctionModLinRngd" ) );

  std::vector< double > nval( idxs.size() );
  std::vector< int > cidx( idxs.size() );
  auto nvit = nval.begin();
  auto idxit = idxs.begin();
  auto cidxit = cidx.begin();
  auto & cp = qf->get_v_var();

  for( auto v : *vars )
   if( auto idx = *(idxit++) ; idx < Inf< Index >() ) {
    *(nvit++) = std::get< 1 >( cp[ idx ] );
    auto cidx = index_of_variable( static_cast< const ColVariable * >( v ) );
    *(cidxit++) = cidx;
    // quadratic coefficients need be changed one at a time
    CPXchgqpcoef( env , lp , cidx , cidx , 2 * std::get< 2 >( cp[ idx ] ) );
    }

  auto nsz = std::distance( nval.begin() , nvit );
  cidx.resize( nsz );
  nval.resize( nsz );

  CPXchgobj( env , lp , idxs.size() , cidx.data() , nval.data() );
  return;
  }

 // Fallback method - Update all costs
 // --------------------------------------------------------------------------
 // reload_objective( f );

 }  // end( CPXMILPSolver::objective_function_modification )

/*--------------------------------------------------------------------------*/

void CPXMILPSolver::constraint_function_modification( const FunctionMod *mod )
{
 // no point in calling the method of MILPSolver, as it does nothing
 // MILPSolver::constraint_function_modification( mod );

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
  throw( std::logic_error( "unexpected *C05FunctionModLin* from FRowConstraint" ) );

 Subset idxs;
 if( auto modlr = dynamic_cast< const C05FunctionModLinRngd * >( modl ) )
  idxs = lf->map_index( modl->vars() , modlr->range() );
 else
  if( auto modls = dynamic_cast< const C05FunctionModLinSbst * >( modl ) )
   idxs = lf->map_index( modl->vars() , modls->subset() );
  else
   throw( std::logic_error( "unknown type of C05FunctionModLinRngd" ) );

 std::vector< double > nval( idxs.size() );
 std::vector< int > cidx( idxs.size() );
 auto nvit = nval.begin();
 auto idxit = idxs.begin();
 auto cidxit = cidx.begin();
 auto & cp = lf->get_v_var();

 for( auto v :  modl->vars() )
  if( auto idx = *(idxit++) ; idx < Inf< Index >() ) {
   *(nvit++) = cp[ idx ].second;
   *(cidxit++) = index_of_variable( static_cast< const ColVariable * >( v ) );
   }

 auto nsz = std::distance( nval.begin() , nvit );
 cidx.resize( nsz );
 nval.resize( nsz );

 std::vector< int > rows( nsz , row );
 CPXchgcoeflist( env , lp , nsz , rows.data() , cidx.data() , nval.data() );

 // Fallback method - Reload all coefficients
 // --------------------------------------------------------------------------
 // reload_constraint( lf );

 }  // end( CPXMILPSolver::constraint_function_modification )

/*--------------------------------------------------------------------------*/

void CPXMILPSolver::objective_fvars_modification( const FunctionModVars *mod )
{
 // this is called in response to Variable being added to / removed from the
 // Objective; however, note that all Variable are supposed to exist at the
 // time this is called, so index_of_variable() is always correct

 // no point in calling the method of MILPSolver, as it does nothing
 // MILPSolver::objective_fvars_modification( mod );

 auto nv = mod->vars().size();
 if( ! nv )
  return;

 // check the modification type
 if( ( ! dynamic_cast< const C05FunctionModVarsAddd * >( mod ) ) &&
     ( ! dynamic_cast< const C05FunctionModVarsRngd * >( mod ) ) &&
     ( ! dynamic_cast< const C05FunctionModVarsSbst * >( mod ) ) )
  throw( std::invalid_argument( "This type of FunctionModVars is not handled"
				) );
 std::vector< int > indices;
 indices.reserve( nv );
 std::vector< double > values;
 values.reserve( nv );

 auto f = mod->function();
 auto nav = f->get_num_active_var();

 // while changing the coefficients, we have to be careful about the fact
 // that Modification are managed asynchronously with the model changes
 // although the added/removed Variable do exist in the internal data
 // structure of [CPX]MILPSolver since the Modification are managed
 // strictly in arrival order, they may no longer exist in the model;
 // more to the point, they may no longer be active in the LinearFunction

 if( auto lf = dynamic_cast< const LinearFunction * >( f ) ) {
  // Linear objective function

  for( auto v : mod->vars() ) {
   auto var = static_cast< const ColVariable * >( v );
   if( auto idx = index_of_variable( var ) ; idx < Inf< int >() ) {
    indices.push_back( idx );
    if( mod->added() ) {
     auto cidx = lf->is_active( var );
     values.push_back( cidx < nav ? lf->get_coefficient( cidx ) : 0 );
     }
    else
     values.push_back( 0 );
    }
   }

  CPXchgobj( env , lp , indices.size() , indices.data() , values.data() );
  return;
  }

 if( auto qf = dynamic_cast< const DQuadFunction * >( f ) ) {
  // Quadratic objective function

  for( auto v : mod->vars() ) {
   auto var = static_cast< const ColVariable * >( v );
   if( auto ind = index_of_variable( var ) ; ind < Inf< int >() ) {
    indices.push_back( ind );
    double value = 0;
    double q_value = 0;

    if( mod->added() )
     if( auto idx = qf->is_active( var ) ; idx < nav ) {
       value = qf->get_linear_coefficient( idx );
       q_value = qf->get_quadratic_coefficient( idx );
       }

    values.push_back( value );
    CPXchgqpcoef( env , lp , ind , ind , 2 * q_value );
    }
   }

  CPXchgobj( env , lp , indices.size() , indices.data() , values.data() );
  return;
  }

 // This should never happen
 throw( std::invalid_argument( "Unknown type of Objective Function" ) );

 }  // end( CPXMILPSolver::objective_fvars_modification )

/*--------------------------------------------------------------------------*/

void CPXMILPSolver::constraint_fvars_modification(
						const FunctionModVars * mod )
{
 // this is called in response to Variable being added to / removed from the
 // Constraint; however, note that all Variable are supposed to exist at the
 // time this is called, so index_of_variable() is always correct

 // no point in calling the method of MILPSolver, as it does nothing
 // MILPSolver::constraint_fvars_modification( mod );

 auto lf = dynamic_cast< const LinearFunction * >( mod->function() );
 if( ! lf )
  return;

 auto nv = mod->vars().size();
 if( ! nv )
  return;

 // check the modification type
 if( ( ! dynamic_cast< const C05FunctionModVarsAddd * >( mod ) ) &&
     ( ! dynamic_cast< const C05FunctionModVarsRngd * >( mod ) ) &&
     ( ! dynamic_cast< const C05FunctionModVarsSbst * >( mod ) ) )
  throw( std::invalid_argument( "This type of FunctionModVars is not handled"
				) );

 auto con = dynamic_cast< const FRowConstraint * >( lf->get_Observer() );
 if( ! con )  // TODO: Throw exception?
  return;

 auto cidx = index_of_constraint( con );
 if( cidx == Inf< int >() )
  return;

 std::vector< int > indices;
 indices.reserve( nv );
 std::vector< double > values;
 values.reserve( mod->vars().size() );
 std::vector< int > rows( nv , cidx );

 auto nav = lf->get_num_active_var();

 // while changing the coefficients, we have to be careful about the fact
 // that Modification are managed asynchronously with the model changes
 // although the added/removed Variable do exist in the internal data
 // structure of [CPX]MILPSolver since the Modification are managed
 // strictly in arrival order, they may no longer exist in the model;
 // more to the point, they may no longer be active in the LinearFunction

 // get indices and coefficients
 for( auto v : mod->vars() ) {
  auto var = static_cast< const ColVariable * >( v );
  if( auto vidx = index_of_variable( var ) ; vidx < Inf< int >() ) {
   indices.push_back( vidx );
   if( mod->added() ) {
    auto idx = lf->is_active( var );
    values.push_back( idx < nav ? lf->get_coefficient( idx ) : 0 );
    }
   else
    values.push_back( 0 );
   }
  }

 // update the coefficients
 CPXchgcoeflist( env , lp , indices.size() , rows.data() ,
		 indices.data() , values.data() );

 }  // end( CPXMILPSolver::constraint_fvars_modification )

/*----------------------------------------------------------------------------

void CPXMILPSolver::dynamic_modification( const BlockModAD * mod )
{
 MILPSolver::dynamic_modification( mod );
 }

----------------------------------------------------------------------------*/

void CPXMILPSolver::add_dynamic_constraint( const FRowConstraint * con )
{
 // call the method of MILPSolver to update the dictionaries (only)
 MILPSolver::add_dynamic_constraint( con );

 if( auto f = con->get_function() ) {
  auto lf = dynamic_cast< const LinearFunction * >( f );
  if( ! lf )
   throw( std::invalid_argument( "The Constraint is not linear" ) );

  int nzcnt = lf->get_num_active_var();

  std::array< int , 2 > rmatbeg = { 0 , nzcnt };
  std::vector< int > rmatind;
  rmatind.reserve( nzcnt );
  std::vector< double > rmatval;
  rmatval.reserve( nzcnt );

  // get the coefficients to fill the matrix
  for( auto & el : lf->get_v_var() )
   if( auto idx = index_of_variable( el.first ) ; idx < Inf< int >() ) {
    rmatind.push_back( idx );
    rmatval.push_back( el.second );
    }

  // get the bounds
  auto con_lhs = con->get_lhs();
  auto con_rhs = con->get_rhs();
  double rhs , rngval;
  char sense;

  if( con_lhs == con_rhs ) {
   sense = 'E';
   rhs = con_rhs;
   }
  else
   if( con_lhs == -Inf< double >() ) {
    sense = 'L';
    rhs = con_rhs;
    }
   else
    if( con_rhs == Inf< double >() ) {
     sense = 'G';
     rhs = con_lhs;
     }
    else {
     sense = 'R';
     rhs = con_lhs;
     rngval = con_rhs - con_lhs;
     }

  // update the CPLEX problem
  CPXaddrows( env , lp , 0 , 1 , rmatind.size() , & rhs , & sense ,
              rmatbeg.data() , rmatind.data() , rmatval.data() ,
              nullptr , nullptr );
  if( sense == 'R' ) {
   //!!  int index = index_of_dynamic_constraint( con );
   // the constraint has just been added at the end
   int index = numrows - 1;
   CPXchgrngval( env , lp , 1 , & index , & rngval );
   }
  }
 }  // end( CPXMILPSolver::add_dynamic_constraint )

/*--------------------------------------------------------------------------*/

void CPXMILPSolver::add_dynamic_variable( const ColVariable * var )
{
 bool was_mip = int_vars > 0;  // was a MIP before the change

 // call the method of MILPSolver to update dictionaries (and int_vars)
 MILPSolver::add_dynamic_variable( var );

 bool is_mip = int_vars > 0;  // is a MIP after the change

 // get the bounds
 auto bd = CPXMILPSolver::get_problem_bounds( *var );

 // update the CPLEX problem
 CPXaddcols( env , lp , 1 , 0 , nullptr , nullptr , nullptr , nullptr,
	     & bd[ 0 ] , & bd[ 1 ] , nullptr );

 // update problem type, if necessary
 if( is_mip ) {       // the problem is a MIP
  //!! int idx = index_of_dynamic_variable( var );
  // the variable has just been inserted at the end
  int idx = numcols - 1;

  char new_ctype;     // get the new variable type
  if( var->is_integer() && ( ! relax_int_vars ) ) {
   if( var->is_unitary() && var->is_positive() )
    new_ctype = 'B';  // Binary
   else
    new_ctype = 'I';  // Integer
   }
  else
   new_ctype = 'C';   // Continuous

  if( was_mip )       // and it was a MIP before
   // set the ctype for the one variable
   CPXchgctype( env , lp , 1 , & idx , & new_ctype );
  else {              // it was not a MIP before
   // the first integer variable is added
   // all ctype values must be [re]added to the problem
   switch( CPXgetprobtype( env, lp ) ) {
    case( CPXPROB_LP ): CPXchgprobtype( env , lp , CPXPROB_MILP ); break;
    case( CPXPROB_QP ): CPXchgprobtype( env , lp , CPXPROB_MIQP ); break;
    default: throw( std::runtime_error( "Wrong CPLEX problem type" ) );
    }

   // all variables are continuous
   std::vector< char > ctype( numcols , 'C' );
   ctype[ idx ] = new_ctype;  // save the one that just became binary
   CPXcopyctype( env , lp , ctype.data() );
   }
  }
 else                 // it is not a MIP
  if( was_mip ) {     // but was so
   // the last integer variable was removed, make the problem continuous
   switch( CPXgetprobtype( env , lp ) ) {
    case( CPXPROB_LP ):
    case( CPXPROB_QP ):   break;  // this should not happe, but ...
    case( CPXPROB_MILP ):
    case( CPXPROB_FIXEDMILP ): CPXchgprobtype( env , lp , CPXPROB_LP ); break;
    case( CPXPROB_MIQP ):
    case( CPXPROB_FIXEDMIQP ): CPXchgprobtype( env , lp , CPXPROB_QP ); break;
    default: throw( std::runtime_error( "Wrong CPLEX problem type" ) );
    }
   }
  // else do nothing, the problem stays continuous

 }  // end( CPXMILPSolver::add_dynamic_variable )

/*--------------------------------------------------------------------------*/

void CPXMILPSolver::add_dynamic_bound( const OneVarConstraint * con )
{
 // no point in calling the method of MILPSolver, as it does nothing
 // MILPSolver::add_dynamic_bound( con );

 auto var = static_cast< const ColVariable * >( con->get_active_var( 0 ) );
 if( ! var )
  throw( std::logic_error( "CPXMILPSolver: added a bound on no Variable" ) );

 auto idx = index_of_variable( var );
 if( idx == Inf< int >() )
  throw( std::logic_error( "CPXMILPSolver: added a bound on unknown Variable"
			   ) );

 std::array< int , 2 > indices = { idx , idx };
 static std::array< char , 2 > lu = { 'L' , 'U' };
 auto bd = CPXMILPSolver::get_problem_bounds( *var );

 CPXchgbds( env , lp , 2 , indices.data() , lu.data() , bd.data() );
 }

/*--------------------------------------------------------------------------*/

void CPXMILPSolver::remove_dynamic_constraint( const FRowConstraint * con )
{
 int index = index_of_dynamic_constraint( con );
 if( index == Inf< int >() )
  throw( std::runtime_error( "Dynamic constraint not found" ) );

 CPXdelrows( env , lp , index , index );

 // call the method of MILPSolver to update the dictionaries (only)
 MILPSolver::remove_dynamic_constraint( con );
 }

/*--------------------------------------------------------------------------*/

void CPXMILPSolver::remove_dynamic_variable( const ColVariable * var )
{
 int index = index_of_dynamic_variable( var );
 if( index == Inf< int >() )
  throw( std::runtime_error( "Dynamic variable not found" ) );

 CPXdelcols( env , lp , index , index );

 // call the method of MILPSolver to update the dictionaries (only)
 MILPSolver::remove_dynamic_variable( var );
 }

/*--------------------------------------------------------------------------*/

void CPXMILPSolver::remove_dynamic_bound( const OneVarConstraint * con )
{
 // no point in calling the method of MILPSolver, as it does nothing
 // MILPSolver::remove_dynamic_bound( con );

 // note: this only works because remove_dynamic_constraint[s]() do *not*
 //       clear the removed OneVarConstraint, and therefore we can easily
 //       reconstruct which ColVariable it was about
 auto var = static_cast< const ColVariable * >( con->get_active_var( 0 ) );
 if( ! var )  // this should never happen
  return;     // but in case, there is nothing to do

 int idx = index_of_variable( var );
 if( idx == Inf< int >() )  // the ColVariable has been removed
  return;                   // is strange, but there is nothing to do

 std::array< int , 2 > indices = { idx , idx };
 static std::array< char , 2 > lu = { 'L' , 'U' };
 auto bd = CPXMILPSolver::get_problem_bounds( *var );

 CPXchgbds( env , lp , 2 , indices.data() , lu.data() , bd.data() );
 }

/*--------------------------------------------------------------------------*/

int CPXMILPSolver::callback( CPXCALLBACKCONTEXTptr context ,
			     CPXLONG contextid )
{
 // main switch: depending on contextid - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 switch( contextid ) {
  case( CPX_CALLBACKCONTEXT_LOCAL_PROGRESS ):
  case( CPX_CALLBACKCONTEXT_GLOBAL_PROGRESS ): {
   // local or global progress- - - - - - - - - - - - - - - - - - - - - - - -
   // check upper / lower bounds and in case stop

   double solv , bndv;
   CPXcallbackgetinfodbl( context , CPXCALLBACKINFO_BEST_SOL , & solv );
   CPXcallbackgetinfodbl( context , CPXCALLBACKINFO_BEST_BND , & bndv );

   if( get_objsense() == 1 ) {
    // a minimization problem: solv is upper bound and bndv is lower bound
    if( solv >= 1e+75 )
     solv = Inf< double >();

    if( bndv <= - 1e+75 )
     bndv = -Inf< double >();

    if( ( bndv >= up_cut_off() ) || ( solv <= lw_cut_off() ) )
     CPXcallbackabort( context );
    }
   else {
    // a maximization problem: solv is lower bound and bndv is upper bound
    if( solv <= -1e+75 )
     solv = -Inf< double >();

    if( bndv >= 1e+75 )
     bndv = Inf< double >();

    if( ( solv >= up_cut_off() ) || ( bndv <= lw_cut_off() ) )
     CPXcallbackabort( context );
    }
   break;
   }

  case( CPX_CALLBACKCONTEXT_RELAXATION ): {
   // a relaxation has been solved- - - - - - - - - - - - - - - - - - - - - -
   if( ! ( CutSepPar & 3 ) )  // but we don't do user cut separation
    break;                    // nothing to do

   CPXLONG depth;             // find the depth of the current node
   #if CPX_VERSION <= 12090000
    if( CPXcallbackgetinfolong( context , CPXCALLBACKINFO_NODECOUNT ,
				& depth ) )
   #else
    if( CPXcallbackgetinfolong( context , CPXCALLBACKINFO_NODEDEPTH ,
				& depth ) )
   #endif
     throw( std::runtime_error(
                "Unable to get the depth with CPXcallbackgetinfolong()" ) );

   // if we are at a depth for which separation is not enabled
   if( ( ( ! depth ) && ( ! ( CutSepPar & 1 ) ) ) ||
       ( depth && ( ! ( CutSepPar & 2 ) ) ) )
    break;                    // nothing to do

   // this is a critical section where different CPLEX threads may compete
   // for access to the Block: ensure mutual exclusion
   f_callback_mutex.lock();

   // ensure no interference from other threads (except CPLEX ones) by also
   // lock()-ing the Block
   bool owned = f_Block->is_owned_by( f_id );
   if( ( ! owned ) && ( ! f_Block->lock( f_id ) ) )
    throw( std::runtime_error( "Unable to lock the Block" ) );

   // get the solution of the relaxation
   std::vector< double > x( numcols );
   if( CPXcallbackgetrelaxationpoint( context , x.data() , 0 , numcols - 1 ,
				      nullptr ) )
    throw( std::runtime_error(
       "Unable to get the solution with CPXcallbackgetrelaxationpoint()" ) );

   // write it in the Variable of the Block
   get_var_solution( x );

   // now perform the user cut separation with the right Configuration
   std::vector< int > rmatbeg;
   std::vector< int > rmatind;
   std::vector< double > rmatval;
   std::vector< double > rhs;
   std::vector< char > sense;
   perform_separation( get_cfg( depth ? 1 : 0 ) ,
		       rmatbeg , rmatind , rmatval , rhs , sense );
   if( ! owned )
    f_Block->unlock( f_id );  // unlock the Block

   // critical section ends here, release the mutex
   f_callback_mutex.unlock();

   // if any user cut was generated, add them
   if( ! rmatbeg.empty() ) {
    auto md = ( CutSepPar >> 3 ) & 3;
    std::vector< int > purgeable( rhs.size() ,
				  md == 0 ? CPX_USECUT_FILTER
				          : ( md == 1 ? CPX_USECUT_PURGE
				                      : CPX_USECUT_FORCE )
				  );
    std::vector< int > local( rhs.size() , 0 );
    if( CPXcallbackaddusercuts( context , rhs.size() , rmatind.size() ,
				rhs.data() , sense.data() , rmatbeg.data() ,
				rmatind.data() , rmatval.data() ,
				purgeable.data() , local.data() ) )
     throw( std::logic_error( "problem in CPXcallbackaddusercuts" ) );
    }

   break;
   }

  case( CPX_CALLBACKCONTEXT_CANDIDATE ): {
   // a feasible solution has been found- - - - - - - - - - - - - - - - - - -
   if( ! ( CutSepPar & 4 ) )  // but we don't do lazy constraint separation
    break;                    // nothing to do

   // this is a critical section where different CPLEX threads may compete
   // for access to the Block: ensure mutual exclusion
   f_callback_mutex.lock();

   // ensure no interference from other threads (except CPLEX ones) by also
   // lock()-ing the Block
   bool owned = f_Block->is_owned_by( f_id );
   if( ( ! owned ) && ( ! f_Block->lock( f_id ) ) )
    throw( std::runtime_error( "Unable to lock the Block" ) );

   // get the feasible solution
   std::vector< double > x( numcols );
   if( CPXcallbackgetcandidatepoint( context , x.data() , 0 , numcols - 1 ,
				     nullptr ) )
    throw( std::runtime_error(
       "Unable to get the solution with CPXcallbackgetcandidatepoint()" ) );

   // write it in the Variable of the Block
   get_var_solution( x );

   // now perform the lazy constraint separation with the right Configuration
   std::vector< int > rmatbeg;
   std::vector< int > rmatind;
   std::vector< double > rmatval;
   std::vector< double > rhs;
   std::vector< char > sense;
   perform_separation( get_cfg( 2 ) ,
		       rmatbeg , rmatind , rmatval , rhs , sense );
   if( ! owned )
    f_Block->unlock( f_id );  // unlock the Block

   // critical section ends here, release the mutex
   f_callback_mutex.unlock();

   // if any lazy constraint was generated, add them
   if( ! rmatbeg.empty() )
    if( CPXcallbackrejectcandidate( context , rhs.size() , rmatind.size() ,
				    rhs.data() , sense.data() ,
				    rmatbeg.data() , rmatind.data() ,
				    rmatval.data() ) )
     throw( std::logic_error( "problem in CPXcallbackrejectcandidate" ) );
   }
  }  // end( main switch )- - - - - - - - - - - - - - - - - - - - - - - - - -
     // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 return( 0 );
 }

/*--------------------------------------------------------------------------*/

void CPXMILPSolver::perform_separation( Configuration * cfg ,
					std::vector< int > & rmatbeg ,
					std::vector< int > & rmatind ,
					std::vector< double > & rmatval ,
					std::vector< double > & rhs ,
					std::vector< char > & sense )
{
 // note: we assume the Block to have been lock()-ed already and the solution
 //       (be it from the relaxation or feasible) to have been written in the
 //       Variable of the Block
 //
 // since the Block is lock()-ed we assume that we can freely work with the
 // Modification list as no-one has a reason to change it

 auto nM = v_mod.size();  // current number of Modification in the list
 auto it = v_mod.end();
 if( nM )                 // if the list is not empty
  it = prev( it );        // initialize an iterator to the last element

 // call generate_dynamic_constraint()
 f_Block->generate_dynamic_constraints( cfg );

 if( v_mod.size() <= nM )  // check if new Modification have been inserted
  return;                  // if not, nothing to do

 if( ! nM )                // if the list was empty at the beginning
  it = v_mod.begin();      // start from the beginning
 else                      // the list was nonempty
  ++it;                    // move to the first new element

 rmatbeg.push_back( 0 );   // first element of rmatbeg is fixed

 // main loop: check all new Modification for a Constraint addition
 for( ; it != v_mod.end() ; ++it ) {
  // check if the Modification indicates an added FRowConstraint
  auto tmod = dynamic_cast< const BlockModAdd< FRowConstraint > * >(
								it->get() );
  if( ! tmod )  // if not
   continue;    // next

  // add all the new constraint to the matrix, one by one
  for( auto con : tmod->added() ) {
   if( auto f = con->get_function() ) {
    auto * lf = dynamic_cast< const LinearFunction * >( f );
    if( ! lf )
     throw( std::invalid_argument( "The Constraint is not linear" ) );

    auto nzcnt = lf->get_num_active_var();
    auto sz = rmatind.size();
    rmatind.resize( sz + nzcnt );
    rmatval.resize( sz + nzcnt );

    // get the coefficients to fill the matrix
    auto iit = rmatind.begin() + sz;
    auto vit = rmatval.begin() + sz;
    for( auto & el : lf->get_v_var() ) {
     *(iit++) = index_of_variable( el.first );
     *(vit++) = el.second;
     }

    // get the bounds
    auto con_lhs = con->get_lhs();
    auto con_rhs = con->get_rhs();

    if( con_lhs == con_rhs ) {
     sense.push_back( 'E' );
     rhs.push_back( con_rhs );
     }
    else
     if( con_lhs == -Inf< double >() ) {
      sense.push_back( 'L' );
      rhs.push_back( con_rhs );
      }
     else
      if( con_rhs == Inf< double >() ) {
       sense.push_back( 'G' );
       rhs.push_back( con_lhs );
       }
      else {
       // kludge: the added constraint is ranged LHS <= lf( x ) <= RHS, but
       // CPLEX does not allow cuts to be ranged: hence, separately add
       // the two constraints lf( x ) >= LHS and lf( x ) <= RHS
       sense.push_back( 'G' );
       rhs.push_back( con_lhs );
       auto nsz = rmatind.size();
       rmatbeg.push_back( nsz );
       sense.push_back( 'L' );
       rhs.push_back( con_rhs );
       rmatind.resize( nsz + nzcnt );
       std::copy( rmatind.begin() + sz , rmatind.begin() + nsz ,
                  rmatind.begin() + nsz );
       rmatval.resize( nsz + nzcnt );
       std::copy( rmatval.begin() + sz , rmatval.begin() + nsz ,
                  rmatval.begin() + nsz );
       }

    rmatbeg.push_back( rmatind.size() );
    }
   }  // end( for each added FRowConstraint )
  }  // end( main loop )
 }  // end( CPXMILPSolver::perform_separation )

/*--------------------------------------------------------------------------*/

int CPXMILPSolver::cpx_int_par_map( idx_type par ) const
{
 switch( par ) {
  case( intMaxIter ): return( - CPXPARAM_MIP_Limits_Nodes );
  case( intMaxSol ):  return( CPXPARAM_MIP_Pool_Capacity );
  case( intLogVerb ): return( CPXPARAM_ScreenOutput );
  }

 // CPLEX parameters
 if( ( par >= intFirstCPLEXPar ) && ( par < intLastAlgParCPXS ) ) {
  int cplex_par = SMSpp_to_CPLEX_int_pars[ par - intFirstCPLEXPar ];

  // both int and long CPLEX parameters are handled as SMS++ int parameters
  int type;
  CPXgetparamtype( env , cplex_par , & type );

  if( type == CPX_PARAMTYPE_INT )
   return( cplex_par );
  else
   if( type == CPX_PARAMTYPE_LONG )
    return( - cplex_par );
   else
    throw( std::invalid_argument( std::to_string( par ) +
				  " is not a int/long Cplex parameter" ) );
  }

 return( 0 );
 }

/*--------------------------------------------------------------------------*/

int CPXMILPSolver::cpx_dbl_par_map( idx_type par ) const
{
 switch( par ) {
  case( dblMaxTime ): return( CPXPARAM_TimeLimit );
  case( dblRelAcc ):  return( CPXPARAM_MIP_Tolerances_MIPGap );
  case( dblAbsAcc ):  return( CPXPARAM_MIP_Tolerances_AbsMIPGap );
  case( dblRAccSol ): return( CPXPARAM_MIP_Pool_RelGap );
  case( dblAAccSol ): return( CPXPARAM_MIP_Pool_AbsGap );
  case( dblFAccSol ): return( CPXPARAM_Simplex_Tolerances_Feasibility );
  }

 if( ( par >= dblFirstCPLEXPar ) && ( par < dblLastAlgParCPXS ) )
  return( SMSpp_to_CPLEX_dbl_pars[ par - dblFirstCPLEXPar ] );

 return( 0 );
 }

/*--------------------------------------------------------------------------*/
/*------------------- METHODS FOR HANDLING THE PARAMETERS ------------------*/
/*--------------------------------------------------------------------------*/

void CPXMILPSolver::set_par( idx_type par , int value )
{
 if( par == intThrowReducedCostException ) {
  throw_reduced_cost_exception = bool( value );
  return;
  }

 if( par == intCutSepPar ) {
  CutSepPar = value;
  return;
  }

 if( int cp = cpx_int_par_map( par ) ) {
  if( cp > 0 )
   CPXsetintparam( env , cp , value );
  else
   CPXsetlongparam( env , - cp , value );
  return;
  }

 MILPSolver::set_par( par, value );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

void CPXMILPSolver::set_par( idx_type par , double value )
{
 // Solver parameters explicitly mapped in CPLEX
 switch( par ) {
  case( dblUpCutOff ): UpCutOff = value; return;
  case( dblLwCutOff ): LwCutOff = value; return;
  }

 if( int cp = cpx_dbl_par_map( par ) ) {
  CPXsetdblparam( env , cp , value );
  return;
  }

 MILPSolver::set_par( par , value );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

void CPXMILPSolver::set_par( idx_type par , std::string && value )
{
 // CPLEX parameters
 if( ( par >= strFirstCPLEXPar ) && ( par < strLastAlgParCPXS ) ) {
  int cplex_par = SMSpp_to_CPLEX_str_pars[ par - strFirstCPLEXPar ];
  CPXsetstrparam( env , cplex_par , value.c_str() );
  return;
  }

 MILPSolver::set_par( par, std::move( value ) );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

void CPXMILPSolver::set_par( idx_type par , std::vector< int > && value )
{
 if( par == vintCutSepCfgInd ) {
  CutSepCfgInd = std::move( value );
  return;
  }

 // MILPSolver and its ancestors have no set_par( std::vector< int > ),
 // so avoid calling it
 // MILPSolver::set_par( par, std::move( value ) );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

void CPXMILPSolver::set_par( idx_type par ,
			     std::vector< std::string > && value )
{
 if( par == vstrConfigDBFName ) {
  auto sz = value.size();
  // delete existing Configuration with index larger than the new size
  for( Index i = sz ; i < v_ConfigDB.size() ; ++i )
   delete v_ConfigDB[ i ];
  // resize the Configuration DB: if the new size is larger than the
  // old ones, the new Configuration defaut to nullptr
  v_ConfigDB.resize( sz , nullptr );
  // resize the configuration names: this makes the next step easier, as
  // any non-existing element will be an empty string and therefore not
  // equal to en existing one unless the existing is empty as well, but
  // this implies that the Configuration is nullptr so it works
  ConfigDBFName.resize( sz );
  // for each new configuration check if the filename is the same as
  // the existing one: if so leave the existing one, otherwise
  // substitute it with a newly loaded one
  for( Index i = 0 ; i < sz ; ++i )
   if( ConfigDBFName[ i ] != value[ i ] ) {
    delete v_ConfigDB[ i ];
    v_ConfigDB[ i ] = Configuration::deserialize( value[ i ] );
    }
  // finally store the new names in place of the existing ones
  ConfigDBFName = std::move( value );
  return;
  }

 // MILPSolver and its ancestors have no
 // set_par( std::vector< std::string > ), so avoid calling it
 // MILPSolver::set_par( par, std::move( value ) );
 }

/*--------------------------------------------------------------------------*/

Solver::idx_type CPXMILPSolver::get_num_int_par( void ) const {
 return( MILPSolver::get_num_int_par()
	 + intLastAlgParCPXS - intLastAlgParMILP );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

Solver::idx_type CPXMILPSolver::get_num_dbl_par( void ) const {
 return( MILPSolver::get_num_dbl_par()
	 + dblLastAlgParCPXS - dblLastAlgParMILP );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

Solver::idx_type CPXMILPSolver::get_num_str_par( void ) const {
 return( MILPSolver::get_num_str_par()
	 + strLastAlgParCPXS - strLastAlgParMILP );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

Solver::idx_type CPXMILPSolver::get_num_vint_par( void ) const {
 return( MILPSolver::get_num_vint_par()
	 + vintLastAlgParCPXS - vintLastAlgParMILP );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

Solver::idx_type CPXMILPSolver::get_num_vstr_par( void ) const {
 return( MILPSolver::get_num_vstr_par()
	 + vstrLastAlgParCPXS - vstrLastAlgParMILP );
 }

/*--------------------------------------------------------------------------*/

int CPXMILPSolver::get_dflt_int_par( idx_type par ) const
{
 if( ( par == intThrowReducedCostException ) || ( par == intCutSepPar ) )
  return( 0 );

 if( int cp = cpx_int_par_map( par ) ) {
  if( cp > 0 ) {
   int value;
   CPXinfointparam( env , cp , & value , nullptr , nullptr );
   return( value );
   }
  else {
   CPXLONG value;
   CPXinfolongparam( env , - cp , & value , nullptr , nullptr );
   return( ( int ) value );
   }
  }

 return( MILPSolver::get_dflt_int_par( par ) );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

double CPXMILPSolver::get_dflt_dbl_par( idx_type par ) const
{
 switch( par ) {
  case( dblUpCutOff ): return( Inf< double >() );
  case( dblLwCutOff ): return( -Inf< double >() );
  }

 if( int cp = cpx_dbl_par_map( par ) ) {
  double value;
  CPXinfodblparam( env , cp , & value , nullptr , nullptr );
  return( value );
  }

 return( MILPSolver::get_dflt_dbl_par( par ) );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

const std::string & CPXMILPSolver::get_dflt_str_par( idx_type par ) const
{
 // note: this implementation is not thread safe and it may lead to elements
 //       of value[] to be allocated more than once with some memory being
 //       lost, but the chances are too slim and the potential drawback too
 //       limted to warrant even a humble std::atomic_flag
 static std::vector< std::string > value( strLastAlgParCPXS -
					  strFirstCPLEXPar );

 if( ( par >= strFirstCPLEXPar ) && ( par < strLastAlgParCPXS ) ) {
  auto i = par - strFirstCPLEXPar;
  if( value[ i ].empty() ) {
   value[ i ].reserve( CPX_STR_PARAM_MAX );
   CPXinfostrparam( env , SMSpp_to_CPLEX_str_pars[ i ] , value[ i ].data() );
   }

  return( value[ i ] );
  }

 return( MILPSolver::get_dflt_str_par( par ) );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

const std::vector< int > & CPXMILPSolver::get_dflt_vint_par( idx_type par )
 const
{
 static std::vector< int > _empty;
 if( par == vintCutSepCfgInd )
  return( _empty );

 return( MILPSolver::get_dflt_vint_par( par ) );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

const std::vector< std::string > & CPXMILPSolver::get_dflt_vstr_par(
							idx_type par ) const
{
 static std::vector< std::string > _empty;
 if( par == vstrConfigDBFName )
  return( _empty );

 return( MILPSolver::get_dflt_vstr_par( par ) );
 }

/*--------------------------------------------------------------------------*/

int CPXMILPSolver::get_int_par( idx_type par ) const
{
 if( par == intThrowReducedCostException )
  return( throw_reduced_cost_exception );

 if( par == intCutSepPar )
  return( CutSepPar );

 if( int cp = cpx_int_par_map( par ) ) {
  if( cp > 0 ) {
   int value;
   CPXgetintparam( env , cp , & value );
   return( value );
   }
  else {
   CPXLONG value;
   CPXgetlongparam( env , - cp , & value );
   return( ( int ) value );
   }
  }

 return( MILPSolver::get_int_par( par ) );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

double CPXMILPSolver::get_dbl_par( idx_type par ) const
{
 switch( par ) {
  case( dblUpCutOff ): return( UpCutOff );
  case( dblLwCutOff ): return( LwCutOff );
  }

 if( int cp = cpx_dbl_par_map( par ) ) {
  double value;
  CPXgetdblparam( env , cp , & value );
  return( value );
  }

 return( MILPSolver::get_dbl_par( par ) );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

const std::string & CPXMILPSolver::get_str_par( idx_type par ) const
{
 static std::string value;

 if( ( par >= strFirstCPLEXPar ) && ( par < strLastAlgParCPXS ) ) {
  int cplex_par = SMSpp_to_CPLEX_str_pars[ par - strFirstCPLEXPar ];
  value.reserve( CPX_STR_PARAM_MAX );
  CPXgetstrparam( env , cplex_par , value.data() );
  return( value );
  }

 return( MILPSolver::get_str_par( par ) );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

const std::vector< int > & CPXMILPSolver::get_vint_par( idx_type par ) const
{
 if( par == vintCutSepCfgInd )
  return( CutSepCfgInd );

 return( MILPSolver::get_vint_par( par ) );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

const std::vector< std::string > & CPXMILPSolver::get_vstr_par( idx_type par )
 const
{
 if( par == vstrConfigDBFName )
  return( ConfigDBFName );

 return( MILPSolver::get_vstr_par( par ) );
 }

/*--------------------------------------------------------------------------*/

Solver::idx_type CPXMILPSolver::int_par_str2idx(
					     const std::string & name ) const
{
 if( name == "intThrowReducedCostException" )
  return( intThrowReducedCostException );

 if( name == "intCutSepPar" )
  return( intCutSepPar );

 /* In CPXMILPSolver::*_par_str2idx() methods we check with MILPSolver first
  * to hide the ugly warning that CPXgetparamnum() shows when a parameter name
  * is not found. */

 idx_type idx = MILPSolver::int_par_str2idx( name );
 if( idx < Inf< idx_type >() )
  return( idx );

 // CPLEX parameters
 int cplex_par;
 if( ! CPXgetparamnum( env , name.c_str() , & cplex_par ) ) {
  auto it = lower_bound( CPLEX_to_SMSpp_int_pars.begin() ,
                         CPLEX_to_SMSpp_int_pars.end() ,
                         std::make_pair( cplex_par , 0 ) );
  return( it->second );
  }

 return( Inf< idx_type >() );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

const std::string & CPXMILPSolver::int_par_idx2str( idx_type idx ) const
{
 static const std::array< std::string , 2 > _pars =
                     { "intThrowReducedCostException" , "intCutSepPar" };
 if( idx == intThrowReducedCostException )
  return( _pars[ 0 ] );

 if( idx == intCutSepPar )
  return( _pars[ 1 ] );

 // note: this implementation is not thread safe and it requires that the
 //       result is used immediately after the call (prior to any other call
 //       to int_par_idx2str()), this may have to be improved upon
 static std::string par_name;

 if( ( idx >= intFirstCPLEXPar ) && ( idx < intLastAlgParCPXS ) ) {
  int cplex_par = SMSpp_to_CPLEX_int_pars[ idx - intFirstCPLEXPar ];
  par_name.reserve( CPX_STR_PARAM_MAX );
  #if CPX_VERSION < 12090000
   CPXgetparamname( env , cplex_par, par_name.data() );
  #else
   CPXgetparamhiername( env , cplex_par, par_name.data() );
  #endif
  return( par_name );
  }

 return( MILPSolver::int_par_idx2str( idx ) );
 }

/*--------------------------------------------------------------------------*/

Solver::idx_type CPXMILPSolver::dbl_par_str2idx( const std::string & name )
 const
{
 /* In CPXMILPSolver::*_par_str2idx() methods we check with MILPSolver first
  * to hide the ugly warning that CPXgetparamnum() shows when a parameter name
  * is not found. */

 idx_type idx = MILPSolver::dbl_par_str2idx( name );
 if( idx < Inf< idx_type >() )
  return( idx );

 // CPLEX parameters
 int cplex_par;
 if( ! CPXgetparamnum( env , name.c_str() , & cplex_par ) ) {
  auto it = lower_bound( CPLEX_to_SMSpp_dbl_pars.begin() ,
                         CPLEX_to_SMSpp_dbl_pars.end() ,
                         std::make_pair( cplex_par , 0 ) );
  return( it->second );
  }

 return( Inf< idx_type >() );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

const std::string & CPXMILPSolver::dbl_par_idx2str( idx_type idx ) const
{
 // note: this implementation is not thread safe and it requires that the
 //       result is used immediately after the call (prior to any other call
 //       to int_par_idx2str()), this may have to be improved upon
 static std::string par_name;

 if( ( idx >= dblFirstCPLEXPar ) && ( idx < dblLastAlgParCPXS ) ) {
  int cplex_par = SMSpp_to_CPLEX_dbl_pars[ idx - dblFirstCPLEXPar ];
  par_name.reserve( CPX_STR_PARAM_MAX );
  #if CPX_VERSION < 12090000
   CPXgetparamname( env , cplex_par , par_name.data() );
  #else
   CPXgetparamhiername( env , cplex_par , par_name.data() );
  #endif
  return( par_name );
  }

 return( MILPSolver::dbl_par_idx2str( idx ) );
 }

/*--------------------------------------------------------------------------*/

Solver::idx_type CPXMILPSolver::str_par_str2idx( const std::string & name )
 const
{
 /* In CPXMILPSolver::*_par_str2idx() methods we check with MILPSolver first
  * to hide the ugly warning that CPXgetparamnum() shows when a parameter name
  * is not found. */

 idx_type idx = MILPSolver::str_par_str2idx( name );
 if( idx < Inf< idx_type >() )
  return( idx );

 // CPLEX parameters
 int cplex_par;
 if( ! CPXgetparamnum( env , name.c_str() , & cplex_par ) ) {
  auto it = lower_bound( CPLEX_to_SMSpp_str_pars.begin() ,
                         CPLEX_to_SMSpp_str_pars.end() ,
                         std::make_pair( cplex_par , 0 ) );
  return( it->second );
  }

 return( Inf< idx_type >() );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

const std::string & CPXMILPSolver::str_par_idx2str( idx_type idx ) const
{
 // note: this implementation is not thread safe and it requires that the
 //       result is used immediately after the call (prior to any other call
 //       to int_par_idx2str()), this may have to be improved upon
 static std::string par_name;

 if( ( idx >= strFirstCPLEXPar ) && ( idx < strLastAlgParCPXS ) ) {
  int cplex_par = SMSpp_to_CPLEX_str_pars[ idx - strFirstCPLEXPar ];
  par_name.reserve( CPX_STR_PARAM_MAX );
  #if CPX_VERSION < 12090000
   CPXgetparamname( env , cplex_par , par_name.data() );
  #else
   CPXgetparamhiername( env , cplex_par , par_name.data() );
  #endif
  return( par_name );
  }

 return( MILPSolver::str_par_idx2str( idx ) );
 }

/*--------------------------------------------------------------------------*/

Solver::idx_type CPXMILPSolver::vint_par_str2idx(
					     const std::string & name ) const
{
 if( name == "vintCutSepCfgInd" )
  return( vintCutSepCfgInd );

 return( MILPSolver::vint_par_str2idx( name ) );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

const std::string & CPXMILPSolver::vint_par_idx2str( idx_type idx ) const
{
 static const std::string _pars = "vintCutSepCfgInd";
 if( idx == vintCutSepCfgInd )
  return( _pars );

 return( MILPSolver::vint_par_idx2str( idx ) );
 }

/*--------------------------------------------------------------------------*/

Solver::idx_type CPXMILPSolver::vstr_par_str2idx(
					     const std::string & name ) const
{
 if( name == "vstrConfigDBFName" )
  return( vstrConfigDBFName );

 return( MILPSolver::vstr_par_str2idx( name ) );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

const std::string & CPXMILPSolver::vstr_par_idx2str( idx_type idx ) const
{
 static const std::string _pars = "vstrConfigDBFName";
 if( idx == vstrConfigDBFName )
  return( _pars );

 return( MILPSolver::vstr_par_idx2str( idx ) );
 }

/*--------------------------------------------------------------------------*/

#ifdef MILPSOLVER_DEBUG

void CPXMILPSolver::check_status( void )
{
 auto tmp = CPXgetnumcols( env , lp );
 if( numcols != tmp )
  DEBUG_LOG( "numcols is " << numcols << " but CPXgetnumcols() returns "
	     << tmp << std::endl );

 tmp = CPXgetnumrows( env , lp );
 if( numrows != tmp )
  DEBUG_LOG( "numrows is " << numrows << " but CPXgetnumrows() returns "
	     << tmp << std::endl );

 tmp = CPXgetnumbin( env , lp ) + CPXgetnumint( env , lp );
 if( int_vars != tmp )
  DEBUG_LOG( "int_vars is " << int_vars << " but CPLEX has actually "
	     << tmp << " integer variables" << std::endl );

 MILPSolver::check_status();
 }

#endif

/*--------------------------------------------------------------------------*/
/*-------------------- PRIVATE METHODS OF THE CLASS ------------------------*/
/*----------------------------------------------------------------------------

void CPXMILPSolver::reload_constraint( const LinearFunction * lf )
{
 // note: this is called in response to a FunctionMod, which means that
 // while the coefficients of the LinearFunction change, their number do
 // not; hence the set of nonzeros does not change. thus, IF WE WERE SURE
 // THAT NO Variable HAVE BEEN REMOVED TO THE CONSTRAINT, this would be
 // simple because we would just have to look at the Variable that are
 // currently active. unfortunately, this is not true because Modification
 // are managed asynchronously with the model changes, and therefore it is
 // possible that some Variable that are in the Constraint for CPXMILPSolver
 // (have a nonzero coefficient) may no longer be in the LinearFunction.
 // this may lead to some coefficients to mistakenly remain nonzero after
 // the call
 //
 // however, if this is the case then there is a Modification in the queue
 // following the one that caused this call where the Variable are removed
 // from the LinearFunction. during the management of that Modification the
 // erroneous nonzero coefficients will be zeroed
 
 auto nv = lf->get_num_active_var();
 if( ! nv )  // but this was an empty Constraint (?)
  return;    // so nothing changes

 if( nv > Index( numcols ) )  // check if some variables have been defined
  nv = numcols;               // that have not yet been added to the matrix

 std::vector< double > vals;
 vals.reserve( nv );
 std::vector< int > cols;
 cols.reserve( nv );
 auto row = index_of_constraint( static_cast< const FRowConstraint * >(
						      lf->get_Observer() ) );
 if( row == Inf< int >() )  // the constraint is not (yet?) there
  return;

 std::vector< int > rows( nv , row );

 for( auto & var : lf->get_v_var() )
  if( auto idx = index_of_variable( var.first ) ; idx < Inf< int >() ) {
   cols.push_back( idx );
   vals.push_back( var.second );
   }

 CPXchgcoeflist( env , lp , cols.size() , rows.data() , cols.data() ,
		 vals.data() );
 }

------------------------------------------------------------------------------

void CPXMILPSolver::reload_objective( Function * f )
{
 // note: this is called in response to a FunctionMod, which means that
 // while the coefficients of the LinearFunction change, their number do
 // not; hence the set of nonzeros does not change. thus, IF WE WERE SURE
 // THAT NO Variable HAVE BEEN REMOVED TO THE CONSTRAINT, this would be
 // simple because we would just have to look at the Variable that are
 // currently active. unfortunately, this is not true because Modification
 // are managed asynchronously with the model changes, and therefore it is
 // possible that some Variable that are in the Objective for CPXMILPSolver
 // (have a nonzero coefficient) may no longer be in the LinearFunction.
 // this may lead to some coefficients to mistakenly remain nonzero after
 // the call
 //
 // however, if this is the case then there is a Modification in the queue
 // following the one that caused this call where the Variable are removed
 // from the LinearFunction. during the management of that Modification the
 // erroneous nonzero coefficients will be zeroed

 auto nv = f->get_num_active_var();
 if( ! nv )  // but this was an empty Objective (?)
  return;    // so nothing changes

 if( nv > Index( numcols ) )  // check if some variables have been defined
  nv = numcols;               // that have not yet been added to the matrix

 std::vector< int > indices;
 indices.reserve( nv );
 std::vector< double > values;
 values.reserve( nv );

 if( auto lf = dynamic_cast< const LinearFunction * >( f ) ) {
  for( auto & el : lf->get_v_var() )
   if( auto idx = index_of_variable( el.first ) ; idx < Inf< int >() ) {
    indices.push_back( idx );
    values.push_back( el.second );
    }

  CPXchgobj( env , lp , indices.size() , indices.data() , values.data() );

  update_problem_type( false );
  return;
  }

 if( auto * qf = dynamic_cast< const DQuadFunction * >( f ) ) {
  for( auto & el : qf->get_v_var() )
   if( auto idx = index_of_variable( std::get< 0 >( el ) ) ;
       idx < Inf< int >() ) {
    // linear coefficients can be changed all at once with CPXchgobj
    indices.push_back( idx );
    values.push_back( std::get< 1 >( el ) );
    // quadratic coefficients need be changed one at a time
    CPXchgqpcoef( env , lp , idx , idx , 2 * std::get< 2 >( el ) );
    }

  CPXchgobj( env , lp , indices.size() , indices.data() , values.data() );
  
  update_problem_type( true );
  return;
  }

 // this should never happen
 throw( std::invalid_argument( "Unknown type of Objective Function" ) );

 }  // end( CPXMILPSolver::reload_objective )

----------------------------------------------------------------------------*/

void CPXMILPSolver::update_problem_type( bool quad )
{
 // TODO: I'm not really sure if this is done automatically by CPLEX, check

 if( ! quad ) {
  switch( CPXgetprobtype( env , lp ) ) {
   case( CPXPROB_LP ):
   case( CPXPROB_MILP ):
   case( CPXPROB_FIXEDMILP ): break;
   case( CPXPROB_QP ):        CPXchgprobtype( env , lp , CPXPROB_LP );
                              break;
   case( CPXPROB_MIQP ):      CPXchgprobtype( env , lp , CPXPROB_MILP );
                              break;
   case( CPXPROB_FIXEDMIQP ): CPXchgprobtype( env , lp , CPXPROB_FIXEDMILP );
                              break;
   default: throw( std::runtime_error( "Wrong CPLEX problem type" ) );
   }
  return;
  }

 switch( CPXgetprobtype( env , lp ) ) {
  case( CPXPROB_LP ):   CPXchgprobtype( env , lp , CPXPROB_QP ); break;
  case( CPXPROB_MILP ): CPXchgprobtype( env , lp , CPXPROB_MIQP ); break;
  case( CPXPROB_FIXEDMILP ): CPXchgprobtype( env , lp , CPXPROB_FIXEDMIQP );
                             break;
  case( CPXPROB_QP ):
  case( CPXPROB_MIQP ):
  case( CPXPROB_FIXEDMIQP ): break;
  default: throw( std::runtime_error( "Wrong CPLEX problem type" ) );
  }
 }

/*--------------------------------------------------------------------------*/

Configuration * CPXMILPSolver::get_cfg( Index ci ) const
{
 if( ci >= CutSepCfgInd.size() )
  return( nullptr );
 auto dbi = CutSepCfgInd[ ci ];
 if( ( dbi < 0 ) || ( Index( dbi ) >= v_ConfigDB.size() ) )
  return( nullptr );
 return( v_ConfigDB[ dbi ] );
 }

/*--------------------------------------------------------------------------*/
/*--------------------- End File CPXMILPSolver.cpp -------------------------*/
/*--------------------------------------------------------------------------*/
