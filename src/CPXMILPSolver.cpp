/*--------------------------------------------------------------------------*/
/*------------------------- File CPXMILPSolver.cpp -------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Implementation of the CPXMILPSolver class.
 *
 * \author Antonio Frangioni \n
 *         Operations Research Group \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \author Niccolò Iardella \n
 *         Operations Research Group \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \author Kostas Tavlaridis-Gyparakis \n
 *         Operations Research Group \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \copyright &copy by Antonio Frangioni, Kostas Tavlaridis-Gyparakis, Niccolò Iardella
 */
/*--------------------------------------------------------------------------*/
/*---------------------------- IMPLEMENTATION ------------------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include <functional>
#include <queue>
#include <iomanip>

#include <Block.h>
#include <MILPSolver.h>
#include <OneVarConstraint.h>
#include <LinearFunction.h>
#include <DQuadFunction.h>

#include "CPXMILPSolver.h"

// Include the proper CPLEX parameter mapping
#include <boost/preprocessor/cat.hpp>
#include <boost/preprocessor/stringize.hpp>
#include BOOST_PP_STRINGIZE( BOOST_PP_CAT( BOOST_PP_CAT( CPX, CPX_VERSION ), _maps.h ) )

/*--------------------------------------------------------------------------*/
/*------------------------- NAMESPACE AND USING ----------------------------*/
/*--------------------------------------------------------------------------*/

using namespace SMSpp_di_unipi_it;

SMSpp_insert_in_factory_cpp_0( CPXMILPSolver );

/*--------------------------------------------------------------------------*/
/*--------------------- CONSTRUCTOR AND DESTRUCTOR -------------------------*/
/*--------------------------------------------------------------------------*/

CPXMILPSolver::CPXMILPSolver() : MILPSolver() {
 int status = 0;
 env = CPXopenCPLEX( &status );
 if( env == nullptr ) {
  throw std::runtime_error( "CPXopenCPLEX returned with status " +
                            std::to_string( status ) );
 }
 lp = nullptr;
#ifdef MILPSOLVER_DEBUG
 CPXsetintparam( env, CPXPARAM_Read_DataCheck, CPX_DATACHECK_WARN );
#endif
}

CPXMILPSolver::~CPXMILPSolver() {
 if( lp ) {
  CPXfreeprob( env, &lp );
 }
 CPXcloseCPLEX( &env );
}

/*--------------------------------------------------------------------------*/
/*--------------------- DERIVED METHODS OF BASE CLASS ----------------------*/
/*--------------------------------------------------------------------------*/

void CPXMILPSolver::set_Block( Block * block ) {
 if( block == f_Block ) {
  return;
 }
 MILPSolver::set_Block( block );
}

/*--------------------------------------------------------------------------*/

void CPXMILPSolver::clear_problem( unsigned int what ) {
 MILPSolver::clear_problem( 0 );

 if( lp ) {
  CPXfreeprob( env, &lp );
 }
}

/*--------------------------------------------------------------------------*/

void CPXMILPSolver::load_problem() {
 MILPSolver::load_problem();

 int status = 0;
 lp = CPXcreateprob( env, &status, prob_name.c_str() );

 std::vector< double > cpx_lb = lb;
 std::vector< double > cpx_ub = ub;

 for( int i = 0; i < numcols; ++i ) {
  if( cpx_lb[ i ] == -Inf< double >() ) {
   cpx_lb[ i ] = -CPX_INFBOUND;
  }
  if( cpx_ub[ i ] == Inf< double >() ) {
   cpx_ub[ i ] = CPX_INFBOUND;
  }
 }

 if( use_custom_names ) {
  CPXcopylpwnames( env, lp,
                   numcols,
                   numrows,
                   objsense,
                   objective.data(),
                   rhs.data(),
                   sense.data(),
                   matbeg.data(),
                   matcnt.data(),
                   matind.data(),
                   matval.data(),
                   cpx_lb.data(),
                   cpx_ub.data(),
                   rngval.data(),
                   colname.data(),
                   rowname.data() );
 } else {
  CPXcopylp( env, lp,
             numcols,
             numrows,
             objsense,
             objective.data(),
             rhs.data(),
             sense.data(),
             matbeg.data(),
             matcnt.data(),
             matind.data(),
             matval.data(),
             cpx_lb.data(),
             cpx_ub.data(),
             rngval.data() );
 }

 bool is_qp = std::any_of( q_objective.begin(),
                           q_objective.end(),
                           []( double d ) { return d != 0; } );

 if( is_qp ) {
  // CPLEX evaluates the corresponding objective with a factor
  // of 0.5 in front of the quadratic objective term.
  std::vector< double > double_q_obj = q_objective;
  for( auto & i: double_q_obj ) {
   i = i * 2;
  }

  // Adding q_objective information automatically changes the problem type
  // from linear to quadratic
  CPXcopyqpsep( env, lp, double_q_obj.data() );
 }

 if( int_vars > 0 ) {
  // Adding ctype information automatically changes the problem type
  // from continuous to mixed integer
  CPXcopyctype( env, lp, xctype.data() );
 }

 // The base representation isn't needed anymore
 MILPSolver::clear_problem( 15 );
}

/*--------------------------------------------------------------------------*/

double CPXMILPSolver::get_problem_lb( const ColVariable & var ) {
 double b = MILPSolver::get_problem_lb( var );
 if( b == -Inf< double >() ) {
  b = -CPX_INFBOUND;
 }
 return b;
}

/*--------------------------------------------------------------------------*/

double CPXMILPSolver::get_problem_ub( const ColVariable & var ) {
 double b = MILPSolver::get_problem_ub( var );
 if( b == Inf< double >() ) {
  b = CPX_INFBOUND;
 }
 return b;
}

/*--------------------------------------------------------------------------*/

int CPXMILPSolver::compute( bool changedvars ) {
 if( MILPSolver::compute( changedvars ) != kOK ) {
  // This should never happen
  throw std::runtime_error( "An error occurred in MILPSolver::compute()" );
 }

 int status = 0;
 bool is_qp = false;

 if( !output_file.empty() ) {
  CPXwriteprob( env, lp, output_file.c_str(), "LP" );
 }

 int probtype = CPXgetprobtype( env, lp );
 switch( probtype ) {
  // TODO: remove asserts
  case CPXPROB_LP :
   DEBUG_LOG( "CPLEX problem type: LP" << std::endl );
   break;
  case CPXPROB_MILP :
   DEBUG_LOG( "CPLEX problem type: MILP" << std::endl );
   assert( int_vars > 0 );
   break;
  case CPXPROB_FIXEDMILP :
   DEBUG_LOG( "CPLEX problem type: FIXEDMILP" << std::endl );
   assert( int_vars > 0 );
   break;
  case CPXPROB_QP :
   DEBUG_LOG( "CPLEX problem type: QP" << std::endl );
   is_qp = true;
   break;
  case CPXPROB_MIQP :
   DEBUG_LOG( "CPLEX problem type: MIQP" << std::endl );
   assert( int_vars > 0 );
   is_qp = true;
   break;
  case CPXPROB_FIXEDMIQP :
   DEBUG_LOG( "CPLEX problem type: FIXEDMIQP" << std::endl );
   assert( int_vars > 0 );
   is_qp = true;
   break;
  case CPXPROB_QCP :
   DEBUG_LOG( "CPLEX problem type: QCP" << std::endl );
   throw std::runtime_error( "Unsupported CPLEX problem type" );
  case CPXPROB_MIQCP :
   DEBUG_LOG( "CPLEX problem type: MIQCP" << std::endl );
   throw std::runtime_error( "Unsupported CPLEX problem type" );
  default:
   throw std::runtime_error( "Undefined CPLEX problem type" );
 }

 if( int_vars > 0 ) {
  status = CPXmipopt( env, lp );
  if( status ) {
   // Error
   if( status == CPXERR_SUBPROB_SOLVE ) {
    int substatus = CPXgetsubstat( env, lp );
    sol_status = decode_lqp_status( substatus );
   } else {
    sol_status = decode_cpx_error( status );
   }
   return sol_status;
  }

  status = CPXgetstat( env, lp );
  sol_status = decode_mip_status( status );
  return sol_status;

 } else {
  if( is_qp ) {
   status = CPXqpopt( env, lp );
  } else {
   status = CPXlpopt( env, lp );
  }

  if( status ) {
   // Error
   sol_status = decode_cpx_error( status );
   return sol_status;
  }

  status = CPXgetstat( env, lp );
  sol_status = decode_lqp_status( status );
  return sol_status;
 }
}

/*--------------------------------------------------------------------------*/

int CPXMILPSolver::decode_mip_status( int status ) {
 DEBUG_LOG( "CPXgetstat() returned " << status << std::endl );

 /*
 * The following are the symbols that may represent the status of
 * a CPLEX solution as returned by CPXgetstat() in case of a a MIP,
 * as listed in CPLEX Callable Library API manual.
 * Some cases are commented out as they should never occur in the
 * conditions posed by CPXMILPSolver.
 */
 switch( status ) {
  case CPXMIP_ABORT_FEAS:
   // Stopped, but an integer solution exists.
  case CPXMIP_ABORT_INFEAS:
   // Stopped; no integer solution.
  case CPXMIP_ABORT_RELAXATION_UNBOUNDED:
   // Could not bound convex relaxation of nonconvex (MI)QP.
   return kError;
  case CPXMIP_DETTIME_LIM_FEAS:
   // Deterministic time limit exceeded, but integer solution exists.
  case CPXMIP_DETTIME_LIM_INFEAS:
   // Deterministic time limit exceeded; no integer solution.
   return kStopTime;
  case CPXMIP_FAIL_FEAS:
   // Terminated because of an error, but integer solution exists.
  case CPXMIP_FAIL_FEAS_NO_TREE:
   // Out of memory, no tree available, integer solution exists.
  case CPXMIP_FAIL_INFEAS:
   // Terminated because of an error; no integer solution.
  case CPXMIP_FAIL_INFEAS_NO_TREE:
   // Out of memory, no tree available, no integer solution.
   return kError;
  case CPXMIP_INFEASIBLE:
   // Solution is integer infeasible.
   return kInfeasible;
  case CPXMIP_INForUNBD:
   // Problem has been proven either infeasible or unbounded.
   //!! note: this is typically given by the preprocessor when it finds an
   //!!       easy dual unfeasible ray: since a dual feasible solution has
   //!!       not been constructed yet it is not formally possible to declare
   //!!       the problem unbounded, but both an empty primal and an empty
   //!!       dual is a rare occurrence, so the most likely correct answer is
   //!!       that the problem is unbounded
   //!!   return kInfeasible;
   return kUnbounded;
  case CPXMIP_MEM_LIM_FEAS:
   // Limit on tree memory has been reached, but an integer solution exists.
  case CPXMIP_MEM_LIM_INFEAS:
   // Limit on tree memory has been reached; no integer solution.
   return kError;
  case CPXMIP_NODE_LIM_FEAS:
   // Node limit has been exceeded but integer solution exists.
  case CPXMIP_NODE_LIM_INFEAS:
   // Node limit has been reached; no integer solution.
   return kStopIter;
  case CPXMIP_OPTIMAL:
   // An optimal integer solution has been found.
  case CPXMIP_OPTIMAL_INFEAS:
   // Problem is optimal with unscaled infeasibilities.
  case CPXMIP_OPTIMAL_TOL:
   // An optimal solution within the tolerance defined by
   // the relative or absolute MIP gap has been found.
  case CPXMIP_SOL_LIM:
   // The limit on mixed integer solutions has been reached.
   return kOK;
  case CPXMIP_TIME_LIM_FEAS:
   // Time limit exceeded, but integer solution exists.
  case CPXMIP_TIME_LIM_INFEAS:
   // Time limit exceeded; no integer solution.
   return kStopTime;
  case CPXMIP_UNBOUNDED:
   // Problem has an unbounded ray.
   return kUnbounded;
   // case CPXMIP_ABORT_RELAXED:
   // case CPXMIP_FEASIBLE:
   // case CPXMIP_FEASIBLE_RELAXED_INF:
   // case CPXMIP_FEASIBLE_RELAXED_QUAD:
   // case CPXMIP_FEASIBLE_RELAXED_SUM:
   // case CPXMIP_POPULATESOL_LIM:
   // case CPXMIP_OPTIMAL_POPULATED:
   // case CPXMIP_OPTIMAL_POPULATED_TOL:
   // case CPXMIP_OPTIMAL_RELAXED_INF:
   // case CPXMIP_OPTIMAL_RELAXED_QUAD:
   // case CPXMIP_OPTIMAL_RELAXED_SUM:
  default:;
 }

 throw std::runtime_error( "CPXgetstat() returned an unknown status: " +
                           std::to_string( status ) );
}

/*--------------------------------------------------------------------------*/

int CPXMILPSolver::decode_lqp_status( int status ) {
 DEBUG_LOG( "CPXgetstat() returned " << status << std::endl );

 /*
  * The following are the symbols that may represent the status of
  * a CPLEX solution as returned by CPXgetstat() in case of a LP/QP,
  * or by CPXgetsubstat() in case of a subproblem of a MIP,
  * as listed in CPLEX Callable Library API manual.
  * Some cases are commented out as they should never occur in the
  * conditions posed by CPXMILPSolver.
  */
 switch( status ) {
  case CPX_STAT_ABORT_DETTIME_LIM:
   // Stopped due to a deterministic time limit.
   return kStopTime;
  case CPX_STAT_ABORT_DUAL_OBJ_LIM:
   // Stopped due to a limit on the dual objective.
   return kError;
  case CPX_STAT_ABORT_IT_LIM:
   // Stopped due to limit on number of iterations.
   return kStopIter;
  case CPX_STAT_ABORT_OBJ_LIM:
   // Stopped due to an objective limit.
  case CPX_STAT_ABORT_PRIM_OBJ_LIM:
   // Stopped due to a limit on the primal objective.
   return kError;
  case CPX_STAT_ABORT_TIME_LIM:
   // Stopped due to a time limit.
   return kStopTime;
  case CPX_STAT_ABORT_USER:
   // Stopped due to a request from the user.
   return kError;
  case CPX_STAT_BENDERS_NUM_BEST:
   // Solution is infeasible, but cannot be cut with
   // a Benders cut due to numerical difficulties.
  case CPX_STAT_INFEASIBLE:
   // Problem has been proven infeasible.
   return kInfeasible;
  case CPX_STAT_INForUNBD:
   // Problem has been proven either infeasible or unbounded.
   //!! note: this is typically given by the preprocessor when it finds an
   //!!       easy dual unfeasible ray: since a dual feasible solution has
   //!!       not been constructed yet it is not formally possible to declare
   //!!       the problem unbounded, but both an empty primal and an empty
   //!!       dual is a rare occurrence, so the most likely correct answer is
   //!!       that the problem is unbounded
   //!! return kInfeasible;
   return kUnbounded;
  case CPX_STAT_NUM_BEST:
   // Solution is available, but not proved optimal,
   // due to numeric difficulties during optimization.
  case CPX_STAT_OPTIMAL:
   // Optimal solution is available.
  case CPX_STAT_OPTIMAL_FACE_UNBOUNDED:
   // Model has an unbounded optimal face.
  case CPX_STAT_OPTIMAL_INFEAS:
   // Optimal solution is available, but with infeasibilities after unscaling.
   return kOK;
  case CPX_STAT_UNBOUNDED:
   // Problem has an unbounded ray.
   return kUnbounded;
   // case CPX_STAT_CONFLICT_ABORT_CONTRADICTION:
   // case CPX_STAT_CONFLICT_ABORT_DETTIME_LIM:
   // case CPX_STAT_CONFLICT_ABORT_IT_LIM:
   // case CPX_STAT_CONFLICT_ABORT_MEM_LIM:
   // case CPX_STAT_CONFLICT_ABORT_NODE_LIM:
   // case CPX_STAT_CONFLICT_ABORT_OBJ_LIM:
   // case CPX_STAT_CONFLICT_ABORT_TIME_LIM:
   // case CPX_STAT_CONFLICT_ABORT_USER:
   // case CPX_STAT_CONFLICT_FEASIBLE:
   // case CPX_STAT_CONFLICT_MINIMAL:
   // case CPX_STAT_FEASIBLE:
   // case CPX_STAT_FEASIBLE_RELAXED_INF:
   // case CPX_STAT_FEASIBLE_RELAXED_QUAD:
   // case CPX_STAT_FEASIBLE_RELAXED_SUM:
   // case CPX_STAT_FIRSTORDER:
   // case CPX_STAT_MULTIOBJ_INFEASIBLE:
   // case CPX_STAT_MULTIOBJ_INForUNBD:
   // case CPX_STAT_MULTIOBJ_NON_OPTIMAL:
   // case CPX_STAT_MULTIOBJ_OPTIMAL:
   // case CPX_STAT_MULTIOBJ_STOPPED:
   // case CPX_STAT_MULTIOBJ_UNBOUNDED:
   // case CPX_STAT_OPTIMAL_RELAXED_INF:
   // case CPX_STAT_OPTIMAL_RELAXED_QUAD:
   // case CPX_STAT_OPTIMAL_RELAXED_SUM:
  default:;
 }

 throw std::runtime_error( "CPXgetstat() returned an unknown status: " +
                           std::to_string( status ) );
}

/*--------------------------------------------------------------------------*/

int CPXMILPSolver::decode_cpx_error( int error ) {
 DEBUG_LOG( "CPLEX returned " << error << std::endl );

 /*
  * The following symbols represent error codes returned by CPLEX,
  * for example by CPXlpopt(), CPXqpopt() and CPXmipopt().
  * Not all codes are returned by all the methods, and it's impossible
  * to know what returns what without doing extensive tests.
  * So we are implementing just the ones we encounter as we go.
 */
 switch( error ) {
  // case CPXERR_ABORT_STRONGBRANCH:
  // case CPXERR_ADJ_SIGN_QUAD:
  // case CPXERR_ADJ_SIGNS:
  // case CPXERR_ADJ_SIGN_SENSE:
  // case CPXERR_ARC_INDEX_RANGE:
  // case CPXERR_ARRAY_BAD_SOS_TYPE:
  // case CPXERR_ARRAY_NOT_ASCENDING:
  // case CPXERR_ARRAY_TOO_LONG:
  // case CPXERR_BAD_ARGUMENT:
  // case CPXERR_BAD_BOUND_SENSE:
  // case CPXERR_BAD_BOUND_TYPE:
  // case CPXERR_BAD_CHAR:
  // case CPXERR_BAD_CTYPE:
  // case CPXERR_BAD_DECOMPOSITION:
  // case CPXERR_BAD_DIRECTION:
  // case CPXERR_BAD_EXPONENT:
  // case CPXERR_BAD_EXPO_RANGE:
  // case CPXERR_BAD_FILETYPE:
  // case CPXERR_BAD_ID:
  // case CPXERR_BAD_INDCONSTR:
  // case CPXERR_BAD_INDICATOR:
  // case CPXERR_BAD_INDTYPE:
  // case CPXERR_BAD_LAZY_UCUT:
  // case CPXERR_BAD_LUB:
  // case CPXERR_BAD_METHOD:
  // case CPXERR_BAD_MULTIOBJ_ATTR:
  // case CPXERR_BAD_NUMBER:
  // case CPXERR_BAD_OBJ_SENSE:
  // case CPXERR_BAD_PARAM_NAME:
  // case CPXERR_BAD_PARAM_NUM:
  // case CPXERR_BAD_PIVOT:
  // case CPXERR_BAD_PRIORITY:
  // case CPXERR_BAD_PROB_TYPE:
  // case CPXERR_BAD_ROW_ID:
  // case CPXERR_BAD_SECTION_BOUNDS:
  // case CPXERR_BAD_SECTION_ENDATA:
  // case CPXERR_BAD_SECTION_QMATRIX:
  // case CPXERR_BAD_SENSE:
  // case CPXERR_BAD_SOS_TYPE:
  // case CPXERR_BAD_STATUS:
  // case CPXERR_BAS_FILE_SHORT:
  // case CPXERR_BAS_FILE_SIZE:
  // case CPXERR_BENDERS_MASTER_SOLVE:
  // case CPXERR_CALLBACK:
  // case CPXERR_CALLBACK_INCONSISTENT:
  // case CPXERR_CAND_NOT_POINT:
  // case CPXERR_CAND_NOT_RAY:
  // case CPXERR_CNTRL_IN_NAME:
  // case CPXERR_COL_INDEX_RANGE:
  // case CPXERR_COL_REPEAT_PRINT:
  // case CPXERR_COL_REPEATS:
  // case CPXERR_COL_ROW_REPEATS:
  // case CPXERR_COL_UNKNOWN:
  // case CPXERR_CONFLICT_UNSTABLE:
  // case CPXERR_COUNT_OVERLAP:
  // case CPXERR_COUNT_RANGE:
  // case CPXERR_CPUBINDING_FAILURE:
  // case CPXERR_DBL_MAX:
  // case CPXERR_DECOMPRESSION:
  // case CPXERR_DETTILIM_STRONGBRANCH:
  // case CPXERR_DUP_ENTRY:
  // case CPXERR_DYNFUNC:
  // case CPXERR_DYNLOAD:
  // case CPXERR_ENCODING_CONVERSION:
  // case CPXERR_EXTRA_BV_BOUND:
  // case CPXERR_EXTRA_FR_BOUND:
  // case CPXERR_EXTRA_FX_BOUND:
  // case CPXERR_EXTRA_INTEND:
  // case CPXERR_EXTRA_INTORG:
  // case CPXERR_EXTRA_SOSEND:
  // case CPXERR_EXTRA_SOSORG:
  // case CPXERR_FAIL_OPEN_READ:
  // case CPXERR_FAIL_OPEN_WRITE:
  // case CPXERR_FILE_ENTRIES:
  // case CPXERR_FILE_FORMAT:
  // case CPXERR_FILE_IO:
  // case CPXERR_FILTER_VARIABLE_TYPE:
  // case CPXERR_ILL_DEFINED_PWL:
  // case CPXERR_INDEX_NOT_BASIC:
  // case CPXERR_INDEX_RANGE:
  // case CPXERR_INDEX_RANGE_HIGH:
  // case CPXERR_INDEX_RANGE_LOW:
  // case CPXERR_IN_INFOCALLBACK:
  // case CPXERR_INT_TOO_BIG:
  // case CPXERR_INT_TOO_BIG_INPUT:
  // case CPXERR_INVALID_NUMBER:
  // case CPXERR_LIMITS_TOO_BIG:
  // case CPXERR_LINE_TOO_LONG:
  // case CPXERR_LO_BOUND_REPEATS:
  // case CPXERR_LOCK_CREATE:
  // case CPXERR_LP_NOT_IN_ENVIRONMENT:
  // case CPXERR_LP_PARSE:
  // case CPXERR_MASTER_SOLVE:
  // case CPXERR_MIPSEARCH_WITH_CALLBACKS:
  // case CPXERR_MISS_SOS_TYPE:
  // case CPXERR_MSG_NO_CHANNEL:
  // case CPXERR_MSG_NO_FILEPTR:
  // case CPXERR_MSG_NO_FUNCTION:
  // case CPXERR_MULTIOBJ_SUBPROB_SOLVE:
  // case CPXERR_MULTIPLE_PROBS_IN_REMOTE_ENVIRONMENT:
  // case CPXERR_NAME_CREATION:
  // case CPXERR_NAME_NOT_FOUND:
  // case CPXERR_NAME_TOO_LONG:
  // case CPXERR_NAN:
  // case CPXERR_NEED_OPT_SOLN:
  // case CPXERR_NEGATIVE_SURPLUS:
  // case CPXERR_NET_DATA:
  // case CPXERR_NET_FILE_SHORT:
  // case CPXERR_NO_BARRIER_SOLN:
  // case CPXERR_NO_BASIC_SOLN:
  // case CPXERR_NO_BASIS:
  // case CPXERR_NO_BOUND_SENSE:
  // case CPXERR_NO_BOUND_TYPE:
  // case CPXERR_NO_COLUMNS_SECTION:
  // case CPXERR_NO_CONFLICT:
  // case CPXERR_NO_DECOMPOSITION:
  // case CPXERR_NO_OBJ_NAME:
  // case CPXERR_NODE_INDEX_RANGE:
  // case CPXERR_NODE_ON_DISK:
  // case CPXERR_NO_DUAL_SOLN:
  // case CPXERR_NO_ENDATA:
  // case CPXERR_NO_ENVIRONMENT:
  // case CPXERR_NO_FILENAME:
  // case CPXERR_NO_ID:
  // case CPXERR_NO_ID_FIRST:
  // case CPXERR_NO_INT_X:
  // case CPXERR_NO_KAPPASTATS:
  // case CPXERR_NO_LU_FACTOR:
  // case CPXERR_NO_MEMORY:
  // case CPXERR_NO_MIPSTART:
  // case CPXERR_NO_NAMES:
  // case CPXERR_NO_NAME_SECTION:
  // case CPXERR_NO_NORMS:
  // case CPXERR_NO_NUMBER_BOUND:
  // case CPXERR_NO_NUMBER:
  // case CPXERR_NO_NUMBER_FIRST:
  // case CPXERR_NO_OBJECTIVE:
  // case CPXERR_NO_OBJ_SENSE:
  // case CPXERR_NO_OPERATOR:
  // case CPXERR_NO_OP_OR_SENSE:
  // case CPXERR_NO_ORDER:
  // case CPXERR_NO_PROBLEM:
  // case CPXERR_NO_QP_OPERATOR:
  // case CPXERR_NO_QUAD_EXP:
  // case CPXERR_NO_RHS_COEFF:
  // case CPXERR_NO_RHS_IN_OBJ:
  // case CPXERR_NO_ROW_NAME:
  // case CPXERR_NO_ROW_SENSE:
  // case CPXERR_NO_ROWS_SECTION:
  // case CPXERR_NO_SENSIT:
  // case CPXERR_NO_SOLN:
  // case CPXERR_NO_SOLNPOOL:
  // case CPXERR_NO_SOS:
  // case CPXERR_NO_TREE:
  // case CPXERR_NOT_DUAL_UNBOUNDED:
  // case CPXERR_NOT_FIXED:
  // case CPXERR_NOT_FOR_BENDERS:
  // case CPXERR_NOT_FOR_DISTMIP:
  // case CPXERR_NOT_FOR_MIP:
  // case CPXERR_NOT_FOR_MULTIOBJ:
  // case CPXERR_NOT_FOR_QCP:
  // case CPXERR_NOT_FOR_QP:
  // case CPXERR_NOT_MILPCLASS:
  // case CPXERR_NOT_MIN_COST_FLOW:
  // case CPXERR_NOT_MIP:
  // case CPXERR_NOT_MIQPCLASS:
  // case CPXERR_NOT_ONE_PROBLEM:
  // case CPXERR_NOT_QP:
  // case CPXERR_NOT_SAV_FILE:
  // case CPXERR_NOT_UNBOUNDED:
  // case CPXERR_NO_VECTOR_SOLN:
  // case CPXERR_NULL_POINTER:
  // case CPXERR_ORDER_BAD_DIRECTION:
  // case CPXERR_OVERFLOW:
  // case CPXERR_PARAM_INCOMPATIBLE:
  // case CPXERR_PARAM_TOO_BIG:
  // case CPXERR_PARAM_TOO_SMALL:
  // case CPXERR_PRESLV_ABORT:
  // case CPXERR_PRESLV_BAD_PARAM:
  // case CPXERR_PRESLV_BASIS_MEM:
  // case CPXERR_PRESLV_COPYORDER:
  // case CPXERR_PRESLV_COPYSOS:
  // case CPXERR_PRESLV_CRUSHFORM:
  // case CPXERR_PRESLV_DETTIME_LIM:
  // case CPXERR_PRESLV_DUAL:
  // case CPXERR_PRESLV_FAIL_BASIS:
  // case CPXERR_PRESLV_INF:
  // case CPXERR_PRESLV_INForUNBD:
  // case CPXERR_PRESLV_NO_BASIS:
  // case CPXERR_PRESLV_NO_PROB:
  // case CPXERR_PRESLV_SOLN_MIP:
  // case CPXERR_PRESLV_SOLN_QP:
  // case CPXERR_PRESLV_START_LP:
  // case CPXERR_PRESLV_TIME_LIM:
  // case CPXERR_PRESLV_UNBD:
  // case CPXERR_PRESLV_UNCRUSHFORM:
  // case CPXERR_PRIIND:
  // case CPXERR_PRM_DATA:
  // case CPXERR_PROTOCOL:
  // case CPXERR_QCP_SENSE:
  // case CPXERR_QCP_SENSE_FILE:
  // case CPXERR_Q_DIVISOR:
  // case CPXERR_Q_DUP_ENTRY:
  // case CPXERR_Q_NOT_INDEF:
  // case CPXERR_Q_NOT_POS_DEF:
  // case CPXERR_Q_NOT_SYMMETRIC:
  // case CPXERR_QUAD_EXP_NOT_2:
  // case CPXERR_QUAD_IN_ROW:
  // case CPXERR_RANGE_SECTION_ORDER:
  // case CPXERR_RESTRICTED_VERSION:
  // case CPXERR_RHS_IN_OBJ:
  // case CPXERR_RIMNZ_REPEATS:
  // case CPXERR_RIM_REPEATS:
  // case CPXERR_RIM_ROW_REPEATS:
  // case CPXERR_ROW_INDEX_RANGE:
  // case CPXERR_ROW_REPEAT_PRINT:
  // case CPXERR_ROW_REPEATS:
  // case CPXERR_ROW_UNKNOWN:
  // case CPXERR_SAV_FILE_DATA:
  // case CPXERR_SAV_FILE_VALUE:
  // case CPXERR_SAV_FILE_WRITE:
  // case CPXERR_SBASE_ILLEGAL:
  // case CPXERR_SBASE_INCOMPAT:
  // case CPXERR_SINGULAR:
  // case CPXERR_STR_PARAM_TOO_LONG:
  case CPXERR_SUBPROB_SOLVE:
   // CPXmipopt failed to solve one of the subproblems in the
   // branch-and-cut tree. This failure can be due to a limit
   // (for example, an iteration limit) or due to numeric trouble.
   return kError;
   // case CPXERR_SYNCPRIM_CREATE:
   // case CPXERR_SYSCALL:
   // case CPXERR_THREAD_FAILED:
   // case CPXERR_TILIM_CONDITION_NO:
   // case CPXERR_TILIM_STRONGBRANCH:
   // case CPXERR_TOO_MANY_COEFFS:
   // case CPXERR_TOO_MANY_COLS:
   // case CPXERR_TOO_MANY_RIMNZ:
   // case CPXERR_TOO_MANY_RIMS:
   // case CPXERR_TOO_MANY_ROWS:
   // case CPXERR_TOO_MANY_THREADS:
   // case CPXERR_TREE_MEMORY_LIMIT:
   // case CPXERR_TUNE_MIXED:
   // case CPXERR_UNIQUE_WEIGHTS:
   // case CPXERR_UNSUPPORTED_CONSTRAINT_TYPE:
   // case CPXERR_UNSUPPORTED_OPERATION:
   // case CPXERR_UP_BOUND_REPEATS:
   // case CPXERR_WORK_FILE_OPEN:
   // case CPXERR_WORK_FILE_READ:
   // case CPXERR_WORK_FILE_WRITE:
   // case CPXERR_XMLPARSE:
  default:;
 }

 throw std::runtime_error( "CPLEX returned an unmanaged error: " +
                           std::to_string( error ) );
}

/*--------------------------------------------------------------------------*/

Solver::OFValue CPXMILPSolver::get_lb() {

 OFValue lower_bound = 0;
 int probtype = CPXgetprobtype( env, lp );

 switch( objsense ) {

  // Minimization problem
  case CPX_MIN:

   switch( sol_status ) {

    case kUnbounded:
     lower_bound = -Inf< OFValue >();
     break;

    case kInfeasible:
     lower_bound = Inf< OFValue >();
     break;

    case kOK:
     switch( probtype ) {
      case CPXPROB_MILP:
      case CPXPROB_MIQP:
      case CPXPROB_FIXEDMILP:
      case CPXPROB_FIXEDMIQP:
       CPXgetbestobjval( env, lp, &lower_bound );
       break;
      default:
       // FIXME: It's unclear how to get a lb for a continuous problem here
       CPXgetobjval( env, lp, &lower_bound );
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

   // Maximization problem
  case CPX_MAX:

   switch( sol_status ) {

    case kUnbounded:
     lower_bound = Inf< OFValue >();
     break;

    case kInfeasible:
     lower_bound = -Inf< OFValue >();
     break;

    case kOK:
     CPXgetobjval( env, lp, &lower_bound );
     break;

    default:
     // Same as above
     lower_bound = -Inf< OFValue >();
     break;
   }
   break;

  default:
   throw std::runtime_error( "Objective type not yet defined" );
   break;
 }

 return lower_bound;
}

/*--------------------------------------------------------------------------*/

Solver::OFValue CPXMILPSolver::get_ub() {

 OFValue upper_bound = 0;
 int probtype = CPXgetprobtype( env, lp );

 switch( objsense ) {

  // Minimization problem
  case CPX_MIN:

   switch( sol_status ) {

    case kUnbounded:
     upper_bound = -Inf< OFValue >();
     break;

    case kInfeasible:
     upper_bound = Inf< OFValue >();
     break;

    case kOK:
     CPXgetobjval( env, lp, &upper_bound );
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

   // Maximization problem
  case CPX_MAX:

   switch( sol_status ) {

    case kUnbounded:
     upper_bound = Inf< OFValue >();
     break;

    case kInfeasible:
     upper_bound = -Inf< OFValue >();
     break;

    case kOK:
     switch( probtype ) {
      case CPXPROB_MILP:
      case CPXPROB_MIQP:
      case CPXPROB_FIXEDMILP:
      case CPXPROB_FIXEDMIQP:
       CPXgetbestobjval( env, lp, &upper_bound );
       break;
      default:
       // FIXME: It's unclear how to get a ub for a continuous problem here
       CPXgetobjval( env, lp, &upper_bound );
     }
     break;

    default:
     // Same as above
     upper_bound = Inf< OFValue >();
     break;
   }
   break;

   // Sense not defined
  default:
   throw std::runtime_error( "Objective type not yet defined" );
   break;
 }

 return upper_bound;
}

/*--------------------------------------------------------------------------*/

bool CPXMILPSolver::has_var_solution() {
 int solnmethod, solntype, pfeasind, dfeasind;
 int status = CPXsolninfo( env, lp, &solnmethod, &solntype,
                           &pfeasind, &dfeasind );
 if( status ) {
  throw std::runtime_error( "An error occurred in CPXsolninfo()" );
 }

 switch( solntype ) {
  case CPX_BASIC_SOLN:    // The problem has a simplex basis
  case CPX_NONBASIC_SOLN: // Primal and dual solution but no basis
  case CPX_PRIMAL_SOLN:   // Primal solution but no corresponding dual solution
   return true;
  case CPX_NO_SOLN:       // No solution
  default:
   return false;
 }
}

/*--------------------------------------------------------------------------*/

bool CPXMILPSolver::is_var_feasible() {
 int solnmethod, solntype, pfeasind, dfeasind;
 int status = CPXsolninfo( env, lp, &solnmethod, &solntype,
                           &pfeasind, &dfeasind );
 if( status ) {
  throw std::runtime_error( "An error occurred in CPXsolninfo()" );
 }

 return bool( pfeasind );
}

/*--------------------------------------------------------------------------*/

Solver::OFValue CPXMILPSolver::get_var_value() {
 switch( objsense ) {
  case 1: // Minimization problem
   return get_ub();
  case -1: // Maximization problem
   return get_lb();
  default:
   throw std::runtime_error( "Objective type not yet defined" );
 }
}

/*--------------------------------------------------------------------------*/

void CPXMILPSolver::get_var_solution( Configuration * solc ) {
 // TODO use vector
 auto * x = new double[numcols];
 int status = CPXgetx( env, lp, x, 0, numcols - 1 );
 if( status ) {
  delete[] x;
  throw std::runtime_error( "Unable to get the solution values with CPXgetx()" );
 }

 int col = 0;

 std::queue< Block * > Q;

 bool owned = f_Block->is_owned_by( f_id );
 if( !owned && !f_Block->lock( f_id ) ) {
  throw std::runtime_error( "Unable to lock the Block" );
 }

 Q.push( f_Block );

 while( !Q.empty() ) {
  Block * q_Block = Q.front();
  Q.pop();

  for( auto * i : q_Block->get_nested_Blocks() ) {
   Q.push( i );
  }

  auto set = [ &x, &col ]( ColVariable & v ) {
   v.set_value( x[ col++ ] );
  };

  for( const auto & i : q_Block->get_static_variables() ) {
   un_any_const_static( i, set, un_any_type< ColVariable >() );
  }

  for( const auto & i : q_Block->get_dynamic_variables() ) {
   un_any_const_dynamic( i, set, un_any_type< ColVariable >() );
  }
 }

 if( !owned ) {
  f_Block->unlock( f_id );
 }

 delete[]x;
}

/*--------------------------------------------------------------------------*/

bool CPXMILPSolver::has_dual_solution() {
 int solnmethod, solntype, pfeasind, dfeasind;
 int status = CPXsolninfo( env, lp, &solnmethod, &solntype,
                           &pfeasind, &dfeasind );
 if( status ) {
  throw std::runtime_error( "An error occurred in CPXsolninfo()" );
 }

 switch( solntype ) {
  case CPX_BASIC_SOLN:    // The problem has a simplex basis
  case CPX_NONBASIC_SOLN: // Primal and dual solution but no basis
   return true;
  case CPX_PRIMAL_SOLN:   // Primal solution but no corresponding dual solution
  case CPX_NO_SOLN:       // No solution
  default:
   return false;
 }
}

/*--------------------------------------------------------------------------*/

bool CPXMILPSolver::is_dual_feasible() {
 int solnmethod, solntype, pfeasind, dfeasind;
 int status = CPXsolninfo( env, lp, &solnmethod, &solntype,
                           &pfeasind, &dfeasind );
 if( status ) {
  throw std::runtime_error( "An error occurred in CPXsolninfo()" );
 }

 return bool( dfeasind );
}

/*--------------------------------------------------------------------------*/

void CPXMILPSolver::get_dual_solution( Configuration * solc ) {
 // TODO use vector
 auto * pi = new double[numrows];
 int status = CPXgetpi( env, lp, pi, 0, numrows - 1 );
 if( status ) {
  delete[] pi;
  throw std::runtime_error( "Unable to get the solution values with CPXgetpi()" );
 }

 int row = 0;

 std::queue< Block * > Q;

 bool owned = f_Block->is_owned_by( f_id );
 if( !owned && !f_Block->lock( f_id ) ) {
  throw std::runtime_error( "Unable to lock the Block" );
 }

 Q.push( f_Block );

 while( !Q.empty() ) {
  Block * q_Block = Q.front();
  Q.pop();

  for( auto * i : q_Block->get_nested_Blocks() ) {
   Q.push( i );
  }

  auto set = [ &pi, &row ]( FRowConstraint & c ) {
   c.set_dual( pi[ row++ ] );
  };

  for( const auto & i : q_Block->get_static_constraints() ) {
   un_any_const_static( i, set, un_any_type< FRowConstraint >() );
  }

  for( const auto & i : q_Block->get_dynamic_constraints() ) {
   un_any_const_dynamic( i, set, un_any_type< FRowConstraint >() );
  }
 }

 if( !owned ) {
  f_Block->unlock( f_id );
 }

 delete[]pi;
}

/*--------------------------------------------------------------------------*/

bool CPXMILPSolver::has_dual_direction() {
 // TODO use vector
 auto * y = new double[numrows];
 double proof = 0;
 int status = CPXdualfarkas( env, lp, y, &proof );
 delete[] y;
 return !bool( status );
}

/*--------------------------------------------------------------------------*/

void CPXMILPSolver::get_dual_direction( Configuration * dirc ) {
 // TODO use vector
 auto * y = new double[numrows];
 auto * v = new double[numcols];
 auto * w = new double[numcols];
 auto * dj = new double[numcols];
 double proof = 0;
 int status;

 // CPXdualfarkas gives a Farkas certificate y so that:
 // y' * A * x >= y' * b
 //   If it is a <= constraint then y[i] <= 0 holds;
 //   If it is a >= constraint then y[i] >= 0 holds.
 status = CPXdualfarkas( env, lp, y, &proof );

 if( status ) {
  delete[]y;
  delete[]v;
  delete[]dj;
  delete[]w;
  throw std::runtime_error( "An error occurred in CPXdualfarkas()" );
 }

 // CPXdjfrompi computes reduced costs from dual values
 // dj = c - A'y
 status = CPXdjfrompi( env, lp, y, dj );

 if( status ) {
  delete[]y;
  delete[]v;
  delete[]dj;
  delete[]w;
  throw std::runtime_error( "An error occurred in CPXdjfrompi()" );
 }

 // Dual multipliers for bounds
 for( int i = 0; i < numcols; ++i ) {
  if( dj[ i ] >= 0 ) { // <?
   v[ i ] = dj[ i ]; // - objective[i] ?
   w[ i ] = 0;
  } else {
   v[ i ] = 0;
   w[ i ] = dj[ i ]; // - objective[i] ?
  }
 }

 // ----------------
 int row = 0;

 std::queue< Block * > Q;

 bool owned = f_Block->is_owned_by( f_id );
 if( !owned && !f_Block->lock( f_id ) ) {
  throw std::runtime_error( "Unable to lock the Block" );
 }

 Q.push( f_Block );

 while( !Q.empty() ) {
  Block * q_Block = Q.front();
  Q.pop();

  for( auto * i : q_Block->get_nested_Blocks() ) {
   Q.push( i );
  }

  auto set = [ &y, &row ]( FRowConstraint & c ) {
   c.set_dual( y[ row++ ] );
  };

  for( const auto & i : q_Block->get_static_constraints() ) {
   un_any_const_static( i, set, un_any_type< FRowConstraint >() );
  }

  for( const auto & i : q_Block->get_dynamic_constraints() ) {
   un_any_const_dynamic( i, set, un_any_type< FRowConstraint >() );
  }
 }

 // TODO
 // for( int i = 0; i < numcols; ++i ) {
 //  used_bounds[ i ].first->set_dual( v[ i ] );
 //  used_bounds[ i ].second->set_dual( w[ i ] );
 // }

 if( !owned ) {
  f_Block->unlock( f_id );
 }

 delete[]y;
 delete[]v;
 delete[]dj;
 delete[]w;
}

/*--------------------------------------------------------------------------*/

void CPXMILPSolver::write_lp( const std::string & filename ) {
 CPXwriteprob( env, lp, filename.c_str(), "LP" );
}

int CPXMILPSolver::get_nodes() const {
 return CPXgetnodecnt( env, lp );
}
/*--------------------------------------------------------------------------*/
/*-------------------- PROTECTED FIELDS OF THE CLASS -----------------------*/
/*--------------------------------------------------------------------------*/

void CPXMILPSolver::var_modification( VariableMod * mod ) {
 try {
  MILPSolver::var_modification( mod );
 } catch( std::logic_error & e ) {}

 /*
  * VariableMod class does not include any modification types, so we "refresh"
  * all the CPLEX information for the variable. In particular we:
  *
  *  - reset the ctype (binary, integer, continuous)
  *  - if the variable is fixed, we fix it in CPLEX by doing LHS = RHS
  *  - if the variable is not fixed, we reset LHS and RHS values
  *
  * The CPLEX problem type must also be updated if we remove the last integer
  * variable of a MIP or add a integer variable to a LP/QP.
  */

 auto * var = dynamic_cast<ColVariable *>(mod->variable());
 int idx = index_of_variable( var );
 std::vector< int > indices( 2, idx );

 // Read old variable types
 std::vector< char > ctype;
 int is_mip = CPXgetintvars( &ctype );

 if( is_mip > 0 ) {
  if( ctype[ idx ] == 'B' || ctype[ idx ] == 'I' ) {
   // The variable to be changed was integer, decrease the number
   --is_mip;
  }
 }

 // Read new variable type
 char new_ctype;
 if( var->is_integer() ) {
  // The variable to be changed will be integer, increase the number
  ++is_mip;
  if( var->is_unitary() && var->is_positive() ) {
   new_ctype = 'B'; // Binary
  } else {
   new_ctype = 'I'; // Integer
  }
 } else {
  new_ctype = 'C';  // Continuous
 }

 // Update problem type
 if( is_mip == 0 ) {
  // The last integer variable was removed, or the problem stays continuous
  // This call removes all ctype values
  switch( CPXgetprobtype( env, lp ) ) {
   case CPXPROB_LP :
   case CPXPROB_QP :
    break;
   case CPXPROB_MILP :
   case CPXPROB_FIXEDMILP :
    CPXchgprobtype( env, lp, CPXPROB_LP );
    break;
   case CPXPROB_MIQP :
   case CPXPROB_FIXEDMIQP :
    CPXchgprobtype( env, lp, CPXPROB_QP );
    break;
   default:
    throw std::runtime_error( "Wrong CPLEX problem type" );
  }

 } else if( is_mip == 1 ) {
  // The first integer variable was added
  // All ctype values must be [re]added to the problem
  switch( CPXgetprobtype( env, lp ) ) {
   case CPXPROB_LP :
    CPXchgprobtype( env, lp, CPXPROB_MILP );
    break;
   case CPXPROB_QP :
    CPXchgprobtype( env, lp, CPXPROB_MIQP );
    break;
   default:
    throw std::runtime_error( "Wrong CPLEX problem type" );
  }

  ctype[ idx ] = new_ctype;
  CPXcopyctype( env, lp, ctype.data() );

 } else {
  // The problem stays a MIP, update only the one variable
  CPXchgctype( env, lp, 1, indices.data(), &new_ctype );
 }

 // Update bounds
 std::vector< char > lu;
 std::vector< double > bd;

 if( var->is_fixed() ) {
  lu.resize( 1 );
  bd.resize( 1 );
  lu[ 0 ] = 'B';
  bd[ 0 ] = var->get_value();
  CPXchgbds( env, lp, 1, indices.data(), lu.data(), bd.data() );

 } else {
  lu.resize( 2 );
  bd.resize( 2 );
  lu[ 0 ] = 'L';
  lu[ 1 ] = 'U';
  bd[ 0 ] = get_problem_lb( *var );
  bd[ 1 ] = get_problem_ub( *var );

  CPXchgbds( env, lp, 2, indices.data(), lu.data(), bd.data() );
 }
}

/*--------------------------------------------------------------------------*/

void CPXMILPSolver::objective_modification( ObjectiveMod * mod ) {
 try {
  MILPSolver::objective_modification( mod );
 } catch( std::logic_error & e ) {}

 /*
  * ObjectiveMod class does not include any modification types except
  * for eSetMin and eSetMax.
  * To change OF coefficents, a FunctionMod must be used.
  */

 switch( mod->type() ) {

  case ObjectiveMod::eSetMin:
   CPXchgobjsen( env, lp, CPX_MIN );
   break;

  case ObjectiveMod::eSetMax:
   CPXchgobjsen( env, lp, CPX_MAX );
   break;

  default:
   throw std::invalid_argument( "Invalid type of ObjectiveMod" );
 }
}

/*--------------------------------------------------------------------------*/

void CPXMILPSolver::const_modification( ConstraintMod * mod ) {
 try {
  MILPSolver::const_modification( mod );
 } catch( std::logic_error & e ) {}

 /*
  * To change the coefficents, a FunctionMod must be used.
  */

 auto * p_const = dynamic_cast<FRowConstraint *>(mod->constraint());

 const int cnt = 1;
 std::array< int, cnt > indices{};
 std::array< double, cnt > values{};
 std::array< char, cnt > sense{};
 std::array< double, cnt > rngval{};

 RowConstraint::RHSValue const_lhs = NAN;
 RowConstraint::RHSValue const_rhs = NAN;

 switch( mod->type() ) {

  case ConstraintMod::eRelaxConst:
   // In order to relax the constraint all we do is transform it
   // into an inequality with RHS equal to infinity

   sense[ 0 ] = 'G';
   values[ 0 ] = -CPX_INFBOUND;
   indices[ 0 ] = index_of_constraint( p_const );
   CPXchgrhs( env, lp, cnt, indices.data(), values.data() );
   CPXchgsense( env, lp, cnt, indices.data(), sense.data() );
   break;

  case ConstraintMod::eEnforceConst:
  case RowConstraintMod::eChgLHS:
  case RowConstraintMod::eChgRHS:
  case RowConstraintMod::eChgBTS:
   // In order to enforce a relaxed constraint all we need to do is
   // reverse the process of relaxing it, by changing the sense and
   // the rhs back to the original form of the constraint
   // Moreover, for the way the LP vectors are built, handling
   // LHS/RHS/BTS cases separately is not worth it.

   const_lhs = p_const->get_lhs();
   const_rhs = p_const->get_rhs();

   if( const_lhs == const_rhs ) {
    sense[ 0 ] = 'E';
    values[ 0 ] = const_rhs;
   } else if( const_lhs == -Inf< double >() ) {
    sense[ 0 ] = 'L';
    values[ 0 ] = const_rhs;
   } else if( const_rhs == Inf< double >() ) {
    sense[ 0 ] = 'G';
    values[ 0 ] = const_lhs;
   } else {
    sense[ 0 ] = 'R';
    values[ 0 ] = const_rhs;
    rngval[ 0 ] = const_rhs - const_lhs;
   }

   indices[ 0 ] = index_of_constraint( p_const );
   CPXchgrhs( env, lp, cnt, indices.data(), values.data() );
   CPXchgsense( env, lp, cnt, indices.data(), sense.data() );
   if( sense[ 0 ] == 'R' ) {
    CPXchgrngval( env, lp, cnt, indices.data(), rngval.data() );
   }
   break;

  default:
   throw std::invalid_argument( "Invalid type of ConstraintMod" );
 }
}

/*--------------------------------------------------------------------------*/

void CPXMILPSolver::bound_modification( OneVarConstraintMod * mod ) {
 try {
  MILPSolver::bound_modification( mod );
 } catch( std::logic_error & e ) {}

 /*
  * The same ColVariable can have more active OneVarConstraints,
  * so each time we modify one of them we have to check if LHS and RHS
  * of the Variable change.
  */

 auto * p_const = dynamic_cast<OneVarConstraint *>(mod->constraint());
 auto * p_var = dynamic_cast<ColVariable *>(p_const->get_active_var( 0 ));

 std::vector< int > indices( 2, index_of_variable( p_var ) );
 std::vector< char > lu;
 std::vector< double > bd;

 switch( mod->type() ) {

  case RowConstraintMod::eChgLHS:
   lu.resize( 1 );
   bd.resize( 1 );
   lu[ 0 ] = 'L';
   bd[ 0 ] = get_problem_lb( *p_var );

   CPXchgbds( env, lp, 1, indices.data(), lu.data(), bd.data() );
   break;

  case RowConstraintMod::eChgRHS:
   lu.resize( 1 );
   bd.resize( 1 );
   lu[ 0 ] = 'U';
   bd[ 0 ] = get_problem_ub( *p_var );

   CPXchgbds( env, lp, 1, indices.data(), lu.data(), bd.data() );
   break;

  case RowConstraintMod::eChgBTS:
   lu.resize( 2 );
   bd.resize( 2 );
   lu[ 0 ] = 'L';
   lu[ 1 ] = 'U';
   bd[ 0 ] = get_problem_lb( *p_var );
   bd[ 1 ] = get_problem_ub( *p_var );

   CPXchgbds( env, lp, 2, indices.data(), lu.data(), bd.data() );
   break;

  default:
   throw std::invalid_argument( "Invalid type of OneVarConstraintMod" );
 }
}

/*--------------------------------------------------------------------------*/

void CPXMILPSolver::objective_function_modification( FunctionMod * mod ) {
 try {
  MILPSolver::objective_function_modification( mod );
 } catch( std::logic_error & e ) {}

 auto * mod_f = mod->function();

 std::vector< int > indices;
 std::vector< double > values;
 std::vector< double > q_values;

 if( const auto * lf = dynamic_cast<const LinearFunction *> (mod_f) ) {
  // Linear objective function
  // TODO: Change only involved variables
  indices.reserve( lf->get_num_active_var() );
  values.reserve( lf->get_num_active_var() );

  for( auto el : lf->get_v_var() ) {
   indices.push_back( index_of_variable( el.first ) );
   values.push_back( el.second );
  }

  // FIXME: The following stuff doesn't work
  // if( auto * rngd = dynamic_cast<C05FunctionModRngd *>( mod ) ) {
  //  indices.reserve( rngd->vars().size() );
  //  values.reserve( rngd->vars().size() );
  //
  //  for( auto * it1 : rngd->vars() ) {
  //   for( auto it2: lf->get_v_var() ) {
  //    if( it1 == it2.first ) {
  //     indices.push_back( index_of_variable( it2.first ) );
  //     values.push_back( it2.second );
  //     break;
  //    }
  //   }
  //  }
  // }
  //
  // if( auto * sbst = dynamic_cast<C05FunctionModSbst *>( mod ) ) {
  //  indices.reserve( sbst->vars().size() );
  //  values.reserve( sbst->vars().size() );
  //
  //  for( auto * it1 : sbst->vars() ) {
  //   for( auto it2: lf->get_v_var() ) {
  //    if( it1 == it2.first ) {
  //     indices.push_back( index_of_variable( it2.first ) );
  //     values.push_back( it2.second );
  //     break;
  //    }
  //   }
  //  }
  // }

  if( !indices.empty() ) {
   CPXchgobj( env, lp, indices.size(), indices.data(), values.data() );
  }

  // Update problem type, if needed
  switch( CPXgetprobtype( env, lp ) ) {
   case CPXPROB_LP :
   case CPXPROB_MILP :
   case CPXPROB_FIXEDMILP :
    break;
   case CPXPROB_QP :
    CPXchgprobtype( env, lp, CPXPROB_LP );
    break;
   case CPXPROB_MIQP :
    CPXchgprobtype( env, lp, CPXPROB_MILP );
    break;
   case CPXPROB_FIXEDMIQP :
    CPXchgprobtype( env, lp, CPXPROB_FIXEDMILP );
    break;
   default:
    throw std::runtime_error( "Wrong CPLEX problem type" );
  }
  return;
 }

 if( const auto * qf = dynamic_cast<const DQuadFunction *> (mod_f) ) {
  // Quadratic objective function
  // TODO: Change only involved variables
  indices.reserve( qf->get_num_active_var() );
  values.reserve( qf->get_num_active_var() );
  q_values.reserve( qf->get_num_active_var() );

  for( auto el : qf->get_v_var() ) {
   // Linear coefficients can be changed all at once with CPXchgobj
   indices.push_back( index_of_variable( std::get< 0 >( el ) ) );
   values.push_back( std::get< 1 >( el ) );
   q_values.push_back( std::get< 2 >( el ) );

   // Quadratic coefficients can be changed one at a time
   CPXchgqpcoef( env, lp, indices.back(), indices.back(), q_values.back() );
  }

  // FIXME: The following stuff doesn't work
  // if( auto * rngd = dynamic_cast<C05FunctionModRngd *>( mod ) ) {
  //  indices.reserve( rngd->vars().size() );
  //  values.reserve( rngd->vars().size() );
  //  q_values.reserve( rngd->vars().size() );
  //
  //  for( auto * it1 : rngd->vars() ) {
  //   for( auto it2: qf->get_v_var() ) {
  //    if( it1 == std::get< 0 >( it2 ) ) {
  //     indices.push_back( index_of_variable( std::get< 0 >( it2 ) ) );
  //     values.push_back( std::get< 1 >( it2 ) );
  //     q_values.push_back( std::get< 2 >( it2 ) );
  //
  //     // Quadratic coefficients can be changed one at a time
  //     CPXchgqpcoef( env, lp, indices.back(), indices.back(), q_values.back() );
  //    }
  //   }
  //  }
  // }
  //
  // if( auto * sbst = dynamic_cast<C05FunctionModSbst *>( mod ) ) {
  //  indices.reserve( sbst->vars().size() );
  //  values.reserve( sbst->vars().size() );
  //  q_values.reserve( sbst->vars().size() );
  //
  //
  //  for( auto * it1 : sbst->vars() ) {
  //   for( auto it2: qf->get_v_var() ) {
  //    if( it1 == std::get< 0 >( it2 ) ) {
  //     indices.push_back( index_of_variable( std::get< 0 >( it2 ) ) );
  //     values.push_back( std::get< 1 >( it2 ) );
  //     q_values.push_back( std::get< 2 >( it2 ) );
  //
  //     // Quadratic coefficients can be changed one at a time
  //     CPXchgqpcoef( env, lp, indices.back(), indices.back(), q_values.back() );
  //    }
  //   }
  //  }
  // }

  if( !indices.empty() ) {
   CPXchgobj( env, lp, indices.size(), indices.data(), values.data() );
  }

  // Update problem type, if needed
  switch( CPXgetprobtype( env, lp ) ) {
   case CPXPROB_LP :
    CPXchgprobtype( env, lp, CPXPROB_QP );
    break;
   case CPXPROB_MILP :
    CPXchgprobtype( env, lp, CPXPROB_MIQP );
    break;
   case CPXPROB_FIXEDMILP :
    CPXchgprobtype( env, lp, CPXPROB_FIXEDMIQP );
    break;
   case CPXPROB_QP :
   case CPXPROB_MIQP :
   case CPXPROB_FIXEDMIQP :
    break;
   default:
    throw std::runtime_error( "Wrong CPLEX problem type" );
  }
  return;
 }

 // This should never happen
 throw std::invalid_argument( "Unknown type of Objective Function" );
}

/*--------------------------------------------------------------------------*/

void CPXMILPSolver::constraint_function_modification( FunctionMod * mod ) {
 try {
  MILPSolver::constraint_function_modification( mod );
 } catch( std::logic_error & e ) {}

 auto * mod_f = mod->function();
 const auto * lf = dynamic_cast<const LinearFunction *> (mod_f);

 if( lf == nullptr ) {
  return;
 }

 std::vector< int > indices;
 std::vector< double > values;
 std::vector< int > rows;

 auto * p_const = dynamic_cast<FRowConstraint *>(lf->get_Observer());

 if( auto * rngd = dynamic_cast<C05FunctionModRngd *>( mod ) ) {
  indices.reserve( rngd->vars().size() );
  values.reserve( rngd->vars().size() );
  rows.reserve( rngd->vars().size() );

  for( auto * it1 : rngd->vars() ) {
   for( auto it2: lf->get_v_var() ) {
    if( it1 == it2.first ) {
     indices.push_back( index_of_variable( it2.first ) );
     rows.push_back( index_of_constraint( p_const ) );
     values.push_back( it2.second );
    }
   }
  }
 }

 if( auto * sbst = dynamic_cast<C05FunctionModSbst *>( mod ) ) {
  indices.reserve( sbst->vars().size() );
  values.reserve( sbst->vars().size() );
  rows.reserve( sbst->vars().size() );

  for( auto * it1 : sbst->vars() ) {
   for( auto it2: lf->get_v_var() ) {
    if( it1 == it2.first ) {
     indices.push_back( index_of_variable( it2.first ) );
     rows.push_back( index_of_constraint( p_const ) );
     values.push_back( it2.second );
    }
   }
  }
 }

 if( !indices.empty() ) {
  CPXchgcoeflist( env, lp, indices.size(), rows.data(),
                  indices.data(), values.data() );
 }
}

/*--------------------------------------------------------------------------*/

void CPXMILPSolver::objective_fvars_modification( FunctionModVars * mod ) {
 try {
  MILPSolver::objective_fvars_modification( mod );
 } catch( std::logic_error & e ) {}

 auto * mod_f = mod->function();

 // Check the modification type
 // TODO: Remove this when debugging is done
 auto * add = dynamic_cast<C05FunctionModVarsAddd *>( mod );
 auto * rmvr = dynamic_cast<C05FunctionModVarsRngd *>( mod );
 auto * rmvs = dynamic_cast<C05FunctionModVarsSbst *>( mod );
 if( add == nullptr && rmvr == nullptr && rmvs == nullptr ) {
  throw std::invalid_argument( "This type of FunctionModVars is not handled" );
 }

 std::vector< int > indices;
 std::vector< double > values;
 std::vector< double > q_values;

 indices.reserve( mod->vars().size() );
 values.reserve( mod->vars().size() );

 // TODO: We should also check if new cols must be added/removed,
 //       but at this point it's already done by a dynamic modification.

 if( const auto * lf = dynamic_cast<const LinearFunction *> (mod_f) ) {
  // Linear objective function

  for( auto * it1 : mod->vars() ) {
   for( auto it2: lf->get_v_var() ) {
    if( it1 == it2.first ) {
     indices.push_back( index_of_variable( it2.first ) );
     if( mod->added() ) {
      values.push_back( it2.second );
     } else {
      values.push_back( 0 );
     }
     break;
    }
   }
  }

  if( !indices.empty() ) {
   CPXchgobj( env, lp, indices.size(), indices.data(), values.data() );
  }
  return;
 }

 if( const auto * qf = dynamic_cast<const DQuadFunction *> (mod_f) ) {
  // Quadratic objective function
  q_values.reserve( mod->vars().size() );

  for( auto * it1 : mod->vars() ) {
   for( auto it2: qf->get_v_var() ) {
    if( it1 == std::get< 0 >( it2 ) ) {
     indices.push_back( index_of_variable( std::get< 0 >( it2 ) ) );
     if( mod->added() ) {
      values.push_back( std::get< 1 >( it2 ) );
      q_values.push_back( std::get< 2 >( it2 ) );
     } else {
      values.push_back( 0 );
      q_values.push_back( 0 );
     }
     break;
    }
   }
   CPXchgqpcoef( env, lp, indices.back(), indices.back(), q_values.back() );
  }

  if( !indices.empty() ) {
   CPXchgobj( env, lp, indices.size(), indices.data(), values.data() );
  }
  return;
 }

 // This should never happen
 throw std::invalid_argument( "Unknown type of Objective Function" );
}

/*--------------------------------------------------------------------------*/

void CPXMILPSolver::constraint_fvars_modification( FunctionModVars * mod ) {
 try {
  MILPSolver::constraint_fvars_modification( mod );
 } catch( std::logic_error & e ) {}

 auto * mod_f = mod->function();
 const auto * lf = dynamic_cast<const LinearFunction *> (mod_f);

 if( lf == nullptr ) {
  return;
 }

 // Check the modification type
 // TODO: Remove this when debugging is done
 auto * add = dynamic_cast<C05FunctionModVarsAddd *>( mod );
 auto * rmvr = dynamic_cast<C05FunctionModVarsRngd *>( mod );
 auto * rmvs = dynamic_cast<C05FunctionModVarsSbst *>( mod );
 if( add == nullptr && rmvr == nullptr && rmvs == nullptr ) {
  throw std::invalid_argument( "This type of FunctionModVars is not handled" );
 }

 std::vector< int > indices;
 std::vector< double > values;
 std::vector< int > rows;

 indices.reserve( mod->vars().size() );
 values.reserve( mod->vars().size() );
 rows.reserve( mod->vars().size() );

 // Get indices and coefficients
 auto * p_const = dynamic_cast<FRowConstraint *>(lf->get_Observer());
 for( auto * it1 : mod->vars() ) {
  for( auto it2: lf->get_v_var() ) {
   if( it1 == it2.first ) {
    indices.push_back( index_of_variable( it2.first ) );
    rows.push_back( index_of_constraint( p_const ) );
    if( mod->added() ) {
     values.push_back( it2.second );
    } else {
     values.push_back( 0 );
    }
    break;
   }
  }
 }
 // Update the coefficients (all zeroes)
 if( !indices.empty() ) {
  CPXchgcoeflist( env, lp, indices.size(), rows.data(),
                  indices.data(), values.data() );
 }
}

/*--------------------------------------------------------------------------*/

void CPXMILPSolver::dynamic_modification( BlockModAD * mod ) {
 try {
  MILPSolver::dynamic_modification( mod );
 } catch( std::logic_error & e ) {}
}

/*--------------------------------------------------------------------------*/

void CPXMILPSolver::add_dynamic_constraint( FRowConstraint * p_const ) {
 try {
  MILPSolver::add_dynamic_constraint( p_const );
 } catch( std::logic_error & e ) {}

 const auto * p_fun =
  dynamic_cast<const LinearFunction *>(p_const->get_function());
 if( p_fun == nullptr ) {
  throw std::invalid_argument( "The Constraint is not linear" );
 }

 int nzcnt = p_const->get_num_active_var();

 std::array< int, 2 > rmatbeg = { 0, nzcnt };
 std::vector< int > rmatind( nzcnt );
 std::vector< double > rmatval( nzcnt );

 std::array< double, 1 > rhs{};
 std::array< double, 1 > rngval{};
 std::array< int, 1 > indices{};
 std::array< char, 1 > sense{};

 // Get the coefficients to fill the matrix
 for( int i = 0; i < nzcnt; ++i ) {
  auto * p_var = dynamic_cast<ColVariable *>(p_fun->get_active_var( i ));
  rmatind[ i ] = index_of_variable( p_var );
  rmatval[ i ] = p_fun->get_coefficient( i );
 }

 // Get the bounds
 int new_index = index_of_dynamic_constraint( p_const );
 auto const_lhs = p_const->get_lhs();
 auto const_rhs = p_const->get_rhs();

 if( const_lhs == const_rhs ) {
  sense[ 0 ] = 'E';
  rhs[ 0 ] = const_rhs;

 } else if( const_lhs == -Inf< double >() ) {
  sense[ 0 ] = 'L';
  rhs[ 0 ] = const_rhs;

 } else if( const_rhs == Inf< double >() ) {
  sense[ 0 ] = 'G';
  rhs[ 0 ] = const_lhs;

 } else {
  sense[ 0 ] = 'R';
  rhs[ 0 ] = const_lhs;
  rngval[ 0 ] = const_rhs - const_lhs;
  indices[ 0 ] = new_index;
 }

 // Update the CPLEX problem
 CPXaddrows( env, lp, 0, 1, nzcnt, rhs.data(),
             sense.data(),
             rmatbeg.data(),
             rmatind.data(),
             rmatval.data(),
             nullptr,
             nullptr );
 if( sense[ 0 ] == 'R' ) {
  CPXchgrngval( env, lp, 1, indices.data(), rngval.data() );
 }
}

/*--------------------------------------------------------------------------*/

void CPXMILPSolver::add_dynamic_variable( ColVariable * p_var ) {
 try {
  MILPSolver::add_dynamic_variable( p_var );
 } catch( std::logic_error & e ) {}

 // Build the coefficient matrix for the new variable
 std::vector< int > cmatind;
 std::vector< double > cmatval;

 auto active_constraints = get_active_constraints( *p_var );
 cmatind.reserve( active_constraints.size() );
 cmatval.reserve( active_constraints.size() );

 // Get the coefficients for this variable for each active constraint
 for( auto * p_const : active_constraints ) {

  const auto * p_fun =
   dynamic_cast<const LinearFunction *>(p_const->get_function());
  if( p_fun == nullptr ) {
   throw std::invalid_argument( "The Constraint is not linear" );
  }
  cmatind.push_back( index_of_constraint( p_const ) );
  auto it = find_if( p_fun->get_v_var().begin(),
                     p_fun->get_v_var().end(),
                     [ & ]( LinearFunction::coeff_pair pair ) {
                      return pair.first == p_var;
                     } );
  cmatval.push_back( it->second );
 }

 // Get the bounds
 double lb = get_problem_lb( *p_var );
 double ub = get_problem_ub( *p_var );

 // Update the CPLEX problem
 int nzcnt = cmatind.size();
 std::array< int, 2 > cmatbeg = { 0, nzcnt };
 CPXaddcols( env, lp, 1, nzcnt, nullptr, cmatbeg.data(),
             cmatind.data(), cmatval.data(), &lb, &ub, nullptr );

 //Get the new variable type
 char new_ctype;
 std::vector< char > old_ctype;
 int is_mip = CPXgetintvars( &old_ctype );

 if( p_var->is_integer() ) {
  ++is_mip;
  if( p_var->is_unitary() && p_var->is_positive() ) {
   new_ctype = 'B'; // Binary
  } else {
   new_ctype = 'I'; // Integer
  }
 } else {
  new_ctype = 'C';  // Continuous
 }

 // Update problem type, if necessary
 if( is_mip == 0 ) {
  // The problem wasn't and still isn't a MIP, nothing to do
 } else if( is_mip == 1 ) {
  // The first integer variable was added
  // All ctype values must be [re]added to the problem
  switch( CPXgetprobtype( env, lp ) ) {
   case CPXPROB_LP :
    CPXchgprobtype( env, lp, CPXPROB_MILP );
    break;
   case CPXPROB_QP :
    CPXchgprobtype( env, lp, CPXPROB_MIQP );
    break;
   default:
    throw std::runtime_error( "Wrong CPLEX problem type" );
  }

  old_ctype[ index_of_dynamic_variable( p_var ) ] = new_ctype;
  CPXcopyctype( env, lp, old_ctype.data() );
 } else {
  // The problem stays a MIP, update only the one variable
  std::array< int, 1 > indices = { index_of_dynamic_variable( p_var ) };
  CPXchgctype( env, lp, 1, indices.data(), &new_ctype );
 }
}

/*--------------------------------------------------------------------------*/

void
CPXMILPSolver::add_dynamic_bound( OneVarConstraint * p_bound ) {
 try {
  MILPSolver::add_dynamic_bound( p_bound );
 } catch( std::logic_error & e ) {}

 auto * p_var = dynamic_cast<ColVariable *>(p_bound->get_active_var( 0 ));
 if( p_var != nullptr ) {

  std::vector< int > indices( 2, index_of_variable( p_var ) );
  std::vector< char > lu( 2 );
  std::vector< double > bd( 2 );

  lu[ 0 ] = 'L';
  lu[ 1 ] = 'U';
  bd[ 0 ] = get_problem_lb( *p_var );
  bd[ 1 ] = get_problem_ub( *p_var );

  CPXchgbds( env, lp, 2, indices.data(), lu.data(), bd.data() );
 }
}

/*--------------------------------------------------------------------------*/

void
CPXMILPSolver::remove_dynamic_constraint( const FRowConstraint * p_const ) {

 int index = index_of_dynamic_constraint( p_const );
 if( index < Inf< int >() ) {
  CPXdelrows( env, lp, index, index );

  try {
   MILPSolver::remove_dynamic_constraint( p_const );
  } catch( std::logic_error & e ) {}
 } else {
  throw std::runtime_error( "Dynamic constraint not found" );
 }
}

/*--------------------------------------------------------------------------*/

void
CPXMILPSolver::remove_dynamic_variable( const ColVariable * p_var ) {

 int index = index_of_dynamic_variable( p_var );
 if( index < Inf< int >() ) {
  CPXdelcols( env, lp, index, index );

  try {
   MILPSolver::remove_dynamic_variable( p_var );
  } catch( std::logic_error & e ) {}
 } else {
  throw std::runtime_error( "Dynamic variable not found" );
 }
}

/*--------------------------------------------------------------------------*/

void
CPXMILPSolver::remove_dynamic_bound( const OneVarConstraint * p_bound ) {
 try {
  MILPSolver::remove_dynamic_bound( p_bound );
 } catch( std::logic_error & e ) {}

 auto * p_var = dynamic_cast<ColVariable *>(p_bound->get_active_var( 0 ));
 if( p_var != nullptr ) {

  std::vector< int > indices( 2, index_of_variable( p_var ) );
  std::vector< char > lu( 2 );
  std::vector< double > bd( 2 );

  lu[ 0 ] = 'L';
  lu[ 1 ] = 'U';
  bd[ 0 ] = get_problem_lb( *p_var );
  bd[ 1 ] = get_problem_ub( *p_var );

  CPXchgbds( env, lp, 2, indices.data(), lu.data(), bd.data() );
 }
}

/*--------------------------------------------------------------------------*/
/*------------------- METHODS FOR HANDLING THE PARAMETERS ------------------*/
/*--------------------------------------------------------------------------*/

void CPXMILPSolver::set_par( const idx_type par, const int value ) {
 // A switch-case here makes the code uglier

 // Solver parameters explicitly mapped in CPLEX
 if( par == intMaxIter ) {
  CPXsetlongparam( env, CPXPARAM_MIP_Limits_Nodes, value );
  return;
 }
 if( par == intMaxSol ) {
  CPXsetintparam( env, CPXPARAM_MIP_Pool_Capacity, value );
  return;
 }
 if( par == intLogVerb ) {
  CPXsetintparam( env, CPXPARAM_ScreenOutput, value );
  return;
 }

 // CPLEX parameters
 if( par >= intFirstCPLEXPar && par < intLastAlgParCPXS ) {
  int cplex_par = SMSpp_to_CPLEX_int_pars[ par - intFirstCPLEXPar ];

  // Both int and long CPLEX parameters are handled as SMS++ int parameters
  int type;
  CPXgetparamtype( env, cplex_par, &type );

  if( type == CPX_PARAMTYPE_INT ) {
   CPXsetintparam( env, cplex_par, value );
  }
  if( type == CPX_PARAMTYPE_LONG ) {
   CPXsetlongparam( env, cplex_par, value );
  }
  return;
 }

 MILPSolver::set_par( par, value );
}

/*--------------------------------------------------------------------------*/

void CPXMILPSolver::set_par( idx_type par, const double value ) {
 // A switch-case here makes the code uglier

 // Solver parameters explicitly mapped in CPLEX
 if( par == dblMaxTime ) {
  CPXsetdblparam( env, CPXPARAM_TimeLimit, value );
  return;
 }
 if( par == dblRelAcc ) {
  CPXsetdblparam( env, CPXPARAM_MIP_Tolerances_RelObjDifference, value );
  return;
 }
 if( par == dblAbsAcc ) {
  CPXsetdblparam( env, CPXPARAM_MIP_Tolerances_ObjDifference, value );
  return;
 }
 if( par == dblUpCutOff ) {
  CPXsetdblparam( env, CPXPARAM_MIP_Tolerances_UpperCutoff, value );
  return;
 }
 if( par == dblLwCutOff ) {
  CPXsetdblparam( env, CPXPARAM_MIP_Tolerances_LowerCutoff, value );
  return;
 }
 if( par == dblRAccSol ) {
  CPXsetdblparam( env, CPXPARAM_MIP_Tolerances_MIPGap, value );
  return;
 }
 if( par == dblAAccSol ) {
  CPXsetdblparam( env, CPXPARAM_MIP_Tolerances_AbsMIPGap, value );
  return;
 }
 if( par == dblFAccSol ) {
  CPXsetdblparam( env, CPXPARAM_Simplex_Tolerances_Feasibility, value );
  return;
 }

 // CPLEX parameters
 if( par >= dblFirstCPLEXPar && par < dblLastAlgParCPXS ) {
  int cplex_par = SMSpp_to_CPLEX_dbl_pars[ par - dblFirstCPLEXPar ];
  CPXsetdblparam( env, cplex_par, value );
  return;
 }

 MILPSolver::set_par( par, value );
}

/*--------------------------------------------------------------------------*/

void CPXMILPSolver::set_par( idx_type par, const std::string & value ) {

 // CPLEX parameters
 if( par >= strFirstCPLEXPar && par < strLastAlgParCPXS ) {
  int cplex_par = SMSpp_to_CPLEX_str_pars[ par - strFirstCPLEXPar ];
  CPXsetstrparam( env, cplex_par, value.c_str() );
  return;
 }

 MILPSolver::set_par( par, value );
}

/*--------------------------------------------------------------------------*/

ThinComputeInterface::idx_type CPXMILPSolver::get_num_int_par() const {
 return MILPSolver::get_num_int_par() + intLastAlgParCPXS - intLastAlgParMILP;
}

ThinComputeInterface::idx_type CPXMILPSolver::get_num_dbl_par() const {
 return MILPSolver::get_num_dbl_par() + dblLastAlgParCPXS - dblLastAlgParMILP;
}

ThinComputeInterface::idx_type CPXMILPSolver::get_num_str_par() const {
 return MILPSolver::get_num_str_par() + strLastAlgParCPXS - strLastAlgParMILP;
}

/*--------------------------------------------------------------------------*/

int CPXMILPSolver::get_int_par( idx_type par ) const {
 int value;
 CPXLONG long_value;

 // Solver parameters explicitly mapped in CPLEX
 if( par == intMaxIter ) {
  CPXgetlongparam( env, CPXPARAM_MIP_Limits_Nodes, &long_value );
  return ( int ) long_value;
 }
 if( par == intMaxSol ) {
  CPXgetintparam( env, CPXPARAM_MIP_Pool_Capacity, &value );
  return value;
 }
 if( par == intLogVerb ) {
  CPXgetintparam( env, CPXPARAM_ScreenOutput, &value );
  return value;
 }

 // CPLEX parameters
 if( par >= intFirstCPLEXPar && par < intLastAlgParCPXS ) {
  int cplex_par = SMSpp_to_CPLEX_int_pars[ par - intFirstCPLEXPar ];

  // Both int and long CPLEX parameters are handled as SMS++ int parameters
  int type;
  CPXgetparamtype( env, cplex_par, &type );

  switch( type ) {
   case CPX_PARAMTYPE_INT:
    CPXgetintparam( env, cplex_par, &value );
    return value;
   case CPX_PARAMTYPE_LONG:
    CPXgetlongparam( env, cplex_par, &long_value );
    return ( int ) long_value;
   default:
    break;
  }
 }

 return MILPSolver::get_int_par( par );
}

/*--------------------------------------------------------------------------*/

double CPXMILPSolver::get_dbl_par( idx_type par ) const {
 double value;

 // Solver parameters explicitly mapped in CPLEX
 if( par == dblMaxTime ) {
  CPXgetdblparam( env, CPXPARAM_TimeLimit, &value );
  return value;
 }
 if( par == dblRelAcc ) {
  CPXgetdblparam( env, CPXPARAM_MIP_Tolerances_RelObjDifference, &value );
  return value;
 }
 if( par == dblAbsAcc ) {
  CPXgetdblparam( env, CPXPARAM_MIP_Tolerances_ObjDifference, &value );
  return value;
 }
 if( par == dblUpCutOff ) {
  CPXgetdblparam( env, CPXPARAM_MIP_Tolerances_UpperCutoff, &value );
  return value;
 }
 if( par == dblLwCutOff ) {
  CPXgetdblparam( env, CPXPARAM_MIP_Tolerances_LowerCutoff, &value );
  return value;
 }
 if( par == dblRAccSol ) {
  CPXgetdblparam( env, CPXPARAM_MIP_Tolerances_MIPGap, &value );
  return value;
 }
 if( par == dblAAccSol ) {
  CPXgetdblparam( env, CPXPARAM_MIP_Tolerances_AbsMIPGap, &value );
  return value;
 }
 if( par == dblFAccSol ) {
  CPXgetdblparam( env, CPXPARAM_Simplex_Tolerances_Feasibility, &value );
  return value;
 }

 // CPLEX parameters
 if( par >= dblFirstCPLEXPar && par < dblLastAlgParCPXS ) {
  int cplex_par = SMSpp_to_CPLEX_dbl_pars[ par - dblFirstCPLEXPar ];
  CPXgetdblparam( env, cplex_par, &value );
 }

 return MILPSolver::get_dbl_par( par );
}

/*--------------------------------------------------------------------------*/

const std::string & CPXMILPSolver::get_str_par( const idx_type par ) const {
 // CPLEX parameters
 if( par >= strFirstCPLEXPar && par < strLastAlgParCPXS ) {
  int cplex_par = SMSpp_to_CPLEX_str_pars[ par - strFirstCPLEXPar ];
  char value[CPX_STR_PARAM_MAX];
  CPXgetstrparam( env, cplex_par, value );

  return std::move( std::string( value ) );
 }

 return MILPSolver::get_str_par( par );
}

/*--------------------------------------------------------------------------*/

int CPXMILPSolver::get_dflt_int_par( const idx_type par ) const {
 int value;
 CPXLONG long_value;

 // Solver parameters explicitly mapped in CPLEX
 if( par == intMaxIter ) {
  CPXinfolongparam( env, CPXPARAM_MIP_Limits_Nodes, &long_value,
                    nullptr, nullptr );
  return ( int ) long_value;
 }
 if( par == intMaxSol ) {
  CPXinfointparam( env, CPXPARAM_MIP_Pool_Capacity, &value, nullptr, nullptr );
  return value;
 }
 if( par == intLogVerb ) {
  CPXinfointparam( env, CPXPARAM_ScreenOutput, &value, nullptr, nullptr );
  return value;
 }

 // CPLEX parameters
 if( par >= intFirstCPLEXPar && par < intLastAlgParCPXS ) {
  int cplex_par = SMSpp_to_CPLEX_int_pars[ par - intFirstCPLEXPar ];

  // Both int and long CPLEX parameters are handled as SMS++ int parameters
  int type;
  CPXgetparamtype( env, cplex_par, &type );

  switch( type ) {
   case CPX_PARAMTYPE_INT:
    CPXinfointparam( env, cplex_par, &value, nullptr, nullptr );
    return value;
   case CPX_PARAMTYPE_LONG:
    CPXinfolongparam( env, cplex_par, &long_value, nullptr, nullptr );
    return ( int ) long_value;
   default:
    break;
  }
 }

 return MILPSolver::get_dflt_int_par( par );
}

/*--------------------------------------------------------------------------*/
double CPXMILPSolver::get_dflt_dbl_par( const idx_type par ) const {
 double value;

 // Solver parameters explicitly mapped in CPLEX
 if( par == dblMaxTime ) {
  CPXinfodblparam( env, CPXPARAM_TimeLimit, &value, nullptr, nullptr );
  return value;
 }
 if( par == dblRelAcc ) {
  CPXinfodblparam( env, CPXPARAM_MIP_Tolerances_RelObjDifference, &value,
                   nullptr, nullptr );
  return value;
 }
 if( par == dblAbsAcc ) {
  CPXinfodblparam( env, CPXPARAM_MIP_Tolerances_ObjDifference, &value,
                   nullptr, nullptr );
  return value;
 }
 if( par == dblUpCutOff ) {
  CPXinfodblparam( env, CPXPARAM_MIP_Tolerances_UpperCutoff, &value,
                   nullptr, nullptr );
  return value;
 }
 if( par == dblLwCutOff ) {
  CPXinfodblparam( env, CPXPARAM_MIP_Tolerances_LowerCutoff, &value,
                   nullptr, nullptr );
  return value;
 }
 if( par == dblRAccSol ) {
  CPXinfodblparam( env, CPXPARAM_MIP_Tolerances_MIPGap, &value,
                   nullptr, nullptr );
  return value;
 }
 if( par == dblAAccSol ) {
  CPXinfodblparam( env, CPXPARAM_MIP_Tolerances_AbsMIPGap, &value,
                   nullptr, nullptr );
  return value;
 }
 if( par == dblFAccSol ) {
  CPXinfodblparam( env, CPXPARAM_Simplex_Tolerances_Feasibility, &value,
                   nullptr, nullptr );
  return value;
 }

 // CPLEX parameters
 if( par >= dblFirstCPLEXPar && par < dblLastAlgParCPXS ) {
  int cplex_par = SMSpp_to_CPLEX_dbl_pars[ par - dblFirstCPLEXPar ];
  CPXinfodblparam( env, cplex_par, &value, nullptr, nullptr );
  return value;
 }

 return MILPSolver::get_dflt_dbl_par( par );
}

/*--------------------------------------------------------------------------*/

const std::string &
CPXMILPSolver::get_dflt_str_par( const idx_type par ) const {

 // CPLEX parameters
 if( par >= strFirstCPLEXPar && par < strLastAlgParCPXS ) {
  int cplex_par = SMSpp_to_CPLEX_str_pars[ par - strFirstCPLEXPar ];
  char value[CPX_STR_PARAM_MAX];
  CPXinfostrparam( env, cplex_par, value );

  return std::move( std::string( value ) );
 }

 return MILPSolver::get_dflt_str_par( par );
}

/*--------------------------------------------------------------------------*/

ThinComputeInterface::idx_type
CPXMILPSolver::int_par_str2idx( const std::string & name ) const {

 // CPLEX parameters
 int cplex_par;
 int status = CPXgetparamnum( env, name.c_str(), &cplex_par );
 if( status == 0 ) {
  auto it = lower_bound( CPLEX_to_SMSpp_int_pars.begin(),
                         CPLEX_to_SMSpp_int_pars.end(),
                         std::make_pair( cplex_par, 0 ) );
  return it->second;
 }

 return MILPSolver::int_par_str2idx( name );
}

/*--------------------------------------------------------------------------*/

const std::string &
CPXMILPSolver::int_par_idx2str( const idx_type idx ) const {

 // CPLEX parameters
 if( idx >= intFirstCPLEXPar && idx < intLastAlgParCPXS ) {
  int cplex_par = SMSpp_to_CPLEX_int_pars[ idx - intFirstCPLEXPar ];
  char par_name[CPX_STR_PARAM_MAX];
#if CPX_VERSION < 12090000
  int status = CPXgetparamname( env, cplex_par, par_name );
#else
  int status = CPXgetparamhiername( env, cplex_par, par_name );
#endif
  return std::move( std::string( par_name ) );
 }

 return MILPSolver::int_par_idx2str( idx );
}

/*--------------------------------------------------------------------------*/

ThinComputeInterface::idx_type
CPXMILPSolver::dbl_par_str2idx( const std::string & name ) const {

 // CPLEX parameters
 int cplex_par;
 int status = CPXgetparamnum( env, name.c_str(), &cplex_par );
 if( status == 0 ) {
  auto it = lower_bound( CPLEX_to_SMSpp_dbl_pars.begin(),
                         CPLEX_to_SMSpp_dbl_pars.end(),
                         std::make_pair( cplex_par, 0 ) );
  return it->second;
 }

 return MILPSolver::dbl_par_str2idx( name );
}

/*--------------------------------------------------------------------------*/

const std::string & CPXMILPSolver::dbl_par_idx2str( const idx_type idx ) const {
 // CPLEX parameters
 if( idx >= dblFirstCPLEXPar && idx < dblLastAlgParCPXS ) {
  int cplex_par = SMSpp_to_CPLEX_dbl_pars[ idx - dblFirstCPLEXPar ];
  char par_name[CPX_STR_PARAM_MAX];
#if CPX_VERSION < 12090000
  int status = CPXgetparamname( env, cplex_par, par_name );
#else
  int status = CPXgetparamhiername( env, cplex_par, par_name );
#endif
  return std::move( std::string( par_name ) );
 }

 return MILPSolver::dbl_par_idx2str( idx );
}

/*--------------------------------------------------------------------------*/

ThinComputeInterface::idx_type
CPXMILPSolver::str_par_str2idx( const std::string & name ) const {

 // CPLEX parameters
 int cplex_par;
 int status = CPXgetparamnum( env, name.c_str(), &cplex_par );
 if( status == 0 ) {
  auto it = lower_bound( CPLEX_to_SMSpp_str_pars.begin(),
                         CPLEX_to_SMSpp_str_pars.end(),
                         std::make_pair( cplex_par, 0 ) );
  return it->second;
 }

 return MILPSolver::str_par_str2idx( name );
}

/*--------------------------------------------------------------------------*/

const std::string &
CPXMILPSolver::str_par_idx2str( const idx_type idx ) const {

 // CPLEX parameters
 if( idx >= strFirstCPLEXPar && idx < strLastAlgParCPXS ) {
  int cplex_par = SMSpp_to_CPLEX_str_pars[ idx - strFirstCPLEXPar ];
  char par_name[CPX_STR_PARAM_MAX];
#if CPX_VERSION < 12090000
  int status = CPXgetparamname( env, cplex_par, par_name );
#else
  int status = CPXgetparamhiername( env, cplex_par, par_name );
#endif
  return std::move( std::string( par_name ) );
 }
 return MILPSolver::str_par_idx2str( idx );
}

/*--------------------------------------------------------------------------*/

#ifdef MILPSOLVER_DEBUG

void CPXMILPSolver::check_status() {
 DEBUG_LOG( "Checking CPXMILPSolver status" << std::endl );

 if( numcols != CPXgetnumcols( env, lp ) ) {
  DEBUG_LOG( "numcols is " << numcols
                           << "but CPXgetnumcols() returns "
                           << CPXgetnumcols( env, lp ) << std::endl );
 }

 if( numrows != CPXgetnumrows( env, lp ) ) {
  DEBUG_LOG( "numrows is " << numrows
                           << "but CPXgetnumrows() returns "
                           << CPXgetnumrows( env, lp ) << std::endl );
 }

 int iv = CPXgetintvars( nullptr );
 if( int_vars != iv ) {
  DEBUG_LOG( "int_vars is " << int_vars
                            << "but CPLEX has actually " << iv
                            << " integer variables" << std::endl );
 }
 MILPSolver::check_status();
}

#endif

/*--------------------------------------------------------------------------*/
/*--------------------- PRIVATE FIELDS OF THE CLASS ------------------------*/
/*--------------------------------------------------------------------------*/

int CPXMILPSolver::CPXgetintvars( std::vector< char > * ctype ) {
 int n;
 int current_cols = CPXgetnumcols( env, lp );
 std::vector< char > * good_ctype = nullptr;

 // Create a new vector if needed
 if( ctype == nullptr ) {
  good_ctype = new std::vector< char >( current_cols );
 } else {
  good_ctype = ctype;
  good_ctype->resize( current_cols, 'C' );
 }

 // Retrieve the ctype array
 int status = CPXgetctype( env, lp, good_ctype->data(), 0, current_cols - 1 );

 if( status == 0 ) {
  // Problem is MIP, get the number of integer variables
  n = std::count_if( good_ctype->begin(),
                     good_ctype->end(),
                     []( char c ) { return c != 'C'; } );
 } else if( status == CPXERR_NOT_MIP ) {
  // Problem is not MIP
  n = 0;
 } else {
  throw std::runtime_error( "CPXgetctype() returned " +
                            std::to_string( status ) );
 }

 if( ctype == nullptr ) {
  delete good_ctype;
 }
 return n;
}

/*--------------------------------------------------------------------------*/
/*--------------------- End File CPXMILPSolver.cpp -------------------------*/
/*--------------------------------------------------------------------------*/
