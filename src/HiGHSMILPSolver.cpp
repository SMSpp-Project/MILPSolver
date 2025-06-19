/*--------------------------------------------------------------------------*/
/*------------------------- File HiGHSMILPSolver.cpp -----------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Implementation of the HiGHSMILPSolver class.
 *
 * \author Antonio Frangioni \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \author Enrico Calandrini \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \copyright &copy by Antonio Frangioni, Enrico Calandrini
 */
/*--------------------------------------------------------------------------*/
/*---------------------------- IMPLEMENTATION ------------------------------*/
/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include <queue>

#include <LinearFunction.h>

#include <QuadFunction.h>

#include "HiGHSMILPSolver.h"

#ifdef MILPSolver_DEBUG
 #define DEBUG_LOG( stuff ) std::cout << "[MILPSolver DEBUG] " << stuff
#else
 #define DEBUG_LOG( stuff )
#endif

#include <boost/preprocessor/cat.hpp>
#include <boost/preprocessor/stringize.hpp>
#include BOOST_PP_STRINGIZE( BOOST_PP_CAT( BOOST_PP_CAT( BOOST_PP_CAT( BOOST_PP_CAT( HiGHS , HIGHS_VERSION_MAJOR ) , HIGHS_VERSION_MINOR ) , HIGHS_VERSION_PATCH ) , _maps.h ) )

/*--------------------------------------------------------------------------*/
/*------------------------- NAMESPACE AND USING ----------------------------*/
/*--------------------------------------------------------------------------*/

using namespace SMSpp_di_unipi_it;

/*--------------------------------------------------------------------------*/
/*-------------------------- FACTORY MANAGEMENT ----------------------------*/
/*--------------------------------------------------------------------------*/

SMSpp_insert_in_factory_cpp_0( HiGHSMILPSolver );

/*--------------------------------------------------------------------------*/
/*----------------------------- FUNCTIONS ----------------------------------*/
/*--------------------------------------------------------------------------*/

void HiGHSMILPSolver_callback( const int callback_type,
        const char* message,
        const HighsCallbackDataOut* data_out,
        HighsCallbackDataIn* data_in,
        void* user_callback_data )
{
 // just defer to the class method
 static_cast< HiGHSMILPSolver * >( user_callback_data
        )->callback( callback_type , data_out , data_in );

 return;
}

/*--------------------------------------------------------------------------*/
/*--------------------- CONSTRUCTOR AND DESTRUCTOR -------------------------*/
/*--------------------------------------------------------------------------*/

HiGHSMILPSolver::HiGHSMILPSolver( void ) :
 MILPSolver() , highs( nullptr ) , f_callback_set( false ) ,
 CutSepPar( 0 ) , UpCutOff( Inf< double >() ) , LwCutOff( - Inf< double >() )
{
 // Create a Highs instance
 highs = Highs_create();

 // Set default HiGHS log to 0
 Highs_setBoolOptionValue( highs , "output_flag" , 0 );
 }

/*--------------------------------------------------------------------------*/

HiGHSMILPSolver::~HiGHSMILPSolver()
{
 for( auto el : v_ConfigDB )
  delete el;

 Highs_destroy(highs);

 // Free auxiliary structures
 q_obj_begin.clear();
 q_obj_ind.clear();
 q_obj_val.clear();
 }

 /*--------------------------------------------------------------------------*/
/*--------------------- DERIVED METHODS OF BASE CLASS ----------------------*/
/*--------------------------------------------------------------------------*/

void HiGHSMILPSolver::set_Block( Block * block )
{
 if( block == f_Block )
  return;

 MILPSolver::set_Block( block );
 UpCutOff = Inf< double >();
 LwCutOff = - Inf< double >();
 }

/*--------------------------------------------------------------------------*/

void HiGHSMILPSolver::clear_problem( unsigned int what )
{
 MILPSolver::clear_problem( 0 );

 Highs_clearModel(highs);
 }

 /*--------------------------------------------------------------------------*/

void HiGHSMILPSolver::load_problem( void )
{
 MILPSolver::load_problem();

 q_obj_begin.clear();
 q_obj_ind.clear();
 q_obj_val.clear();

 int model_status = Highs_getModelStatus( highs );
 if( model_status == kHighsModelStatusModelEmpty )
  Highs_clearModel(highs);

 std::vector< double > highs_lb = lb;
 std::vector< double > highs_ub = ub;
 std::vector< double > highs_rhs = rhs;
 std::vector< int > highs_xctype( numcols );

 for( int i = 0 ; i < numcols ; ++i ) {
  // check LB/UB
  if( highs_lb[ i ] == -Inf< double >() )
   highs_lb[ i ] = -kHighsInf;
  if( highs_ub[ i ] == Inf< double >() )
   highs_ub[ i ] = kHighsInf;

  //check xctype
  switch( xctype[ i ] ) {
    case( 'C' ): highs_xctype[ i ] = kHighsVarTypeContinuous;
                break;
    case( 'B' ):
    case( 'I' ): highs_xctype[ i ] = kHighsVarTypeInteger;
                break;
    case( 'S' ): highs_xctype[ i ] = kHighsVarTypeSemiContinuous;
                break;
    case( 'N' ): highs_xctype[ i ] = kHighsVarTypeSemiInteger;
                break;
    default:
     throw( std::runtime_error( "xctype[" + std::to_string( i ) +
            "] not a valid type" ) );
     break;
   }
  }

 /** An array of length at least numrows containing the lefthand side value
  * for each constraint in the constraint matrix. */
 std::vector< double > highs_lhs( numrows );

 /** Due to HiGHS method (it doesn't use rowsense but always set lhs and rhs),
  *  we have to properly set rhs and lhs of each row before of passing them to
  *  the model. */
 for( int i = 0 ; i < numrows ; ++i) {
  switch( sense[ i ] ) {
    case( 'L' ): highs_lhs[ i ] = -kHighsInf;
              break;
    case( 'G' ): highs_lhs[ i ] = highs_rhs[ i ];
              highs_rhs[ i ] = kHighsInf;
              break;
    case( 'E' ): highs_lhs[ i ] = highs_rhs[ i ];
              break;
    case( 'R' ):
              if( rngval [ i ] > 0 ) {
                highs_lhs[ i ] = highs_rhs[ i ];
                highs_rhs[ i ] = highs_lhs[ i ] + rngval [ i ];
              }
              else
                highs_lhs[ i ] = highs_rhs[ i ] + rngval [ i ];
              break;
  }
 }
 
 bool is_mip = std::any_of( xctype.begin() ,
                           xctype.end() ,
                           []( char c ) { return( c == 'B' || c == 'I' || 
                                                c == 'N' ); } );

 // Separable quadratic problem
 bool is_sqp = std::any_of( q_objective.begin() ,
                           q_objective.end() ,
                           []( double d ) { return( d != 0 ); } );

 // General quadratic problem
 bool is_qp = std::any_of( ndq_objective.begin() ,
                           ndq_objective.end() ,
                           []( double d ) { return( d != 0 ); } );

 // Quadratic constrained problem
 bool is_qcp = ( numquadrows > 0 );

 // HiGHS uses different function to instantiate a model based on
 // his type
 int status;
 if( ! is_mip ) { // LP or QP problem
  status = Highs_passLp( highs , numcols , numrows ,
                        matval.size() , kHighsMatrixFormatColwise , objsense ,
                        0.0 , objective.data() , highs_lb.data() , 
                        highs_ub.data() , highs_lhs.data() , highs_rhs.data() ,
                        matbeg.data() , matind.data() , matval.data()
                        );

  if( status == kHighsStatusError )
    throw( std::runtime_error( "Highs_passLp returned with kHighsStatus " +
			      std::to_string( status ) ) );
  }
 else{ // MIP problem
  status = Highs_passMip( highs , numcols , numrows ,
                          matval.size() , kHighsMatrixFormatColwise , objsense ,
                          0.0 , objective.data() , highs_lb.data() , 
                          highs_ub.data() , highs_lhs.data() , highs_rhs.data() ,
                          matbeg.data() , matind.data() , matval.data() , 
                          highs_xctype.data() 
                          );

  if( status == kHighsStatusError )
    throw( std::runtime_error( "Highs_passMip returned with kHighsStatus " +
			      std::to_string( status ) ) );
 }
 
 if( is_sqp || is_qp ) { // QP problem, the Hessian matrix need to be added

  /* HiGHS read the Hessian matrix in sparse column form, so we have 
  * to prepare three different vector:
  * - q_obj_begin: An array of length [ numcols ] containing the starting index
  *   of each column in `index`;
  * - q_obj_ind: An array of length [ num_nz_q ] with indices of hessian matrix
  *   entries 
  * - q_obj_val: An array of length [ num_nz_q ] with values of hessian matrix
  *   entries. */

  // Call specific function to generate the structures required
  generate_qobj_hessian( );

  status = Highs_passHessian( highs, numcols , q_obj_val.size() ,
                      kHighsHessianFormatTriangular , q_obj_begin.data() ,
                      q_obj_ind.data() , q_obj_val.data() 
                      );
    
  if( status == kHighsStatusError )
   throw( std::runtime_error( "Highs_passHessian returned with kHighsStatus " +
        std::to_string( status ) ) );
 }

 if( is_qcp )
  throw( std::runtime_error( "HiGHS cannot solve QCP models" ) );

 // names must be added manually
 if( use_custom_names ) {
  for( int j = 0 ; j < numcols ; ++j )
    Highs_passColName( highs , j , colname[ j ] );

  for( int i = 0 ; i < numrows ; ++i )
    Highs_passRowName( highs , i , rowname[ i ] );
 }

 // the base representation isn't needed anymore
 MILPSolver::clear_problem( 15 );

 UpCutOff = Inf< double >();
 LwCutOff = - Inf< double >();

 }  // end( HiGHSMILPSolver::load_problem )

/*--------------------------------------------------------------------------*/

double HiGHSMILPSolver::get_problem_lb( const ColVariable & var ) const
{
 double b = MILPSolver::get_problem_lb( var );
 if( b == -Inf< double >() )
  b = -kHighsInf;

 return( b );
 }

/*--------------------------------------------------------------------------*/

double HiGHSMILPSolver::get_problem_ub( const ColVariable & var ) const
{
 double b = MILPSolver::get_problem_ub( var );
 if( b == Inf< double >() )
  b = kHighsInf;

 return( b );
 }

/*--------------------------------------------------------------------------*/

std::array< double , 2 > HiGHSMILPSolver::get_problem_bounds(
					      const ColVariable & var ) const
{
 auto ret = MILPSolver::get_problem_bounds( var );
 if( ret[ 0 ] == -Inf< double >() )
  ret[ 0 ] = -kHighsInf;
 if( ret[ 1 ] == Inf< double >() )
  ret[ 1 ] = kHighsInf;

 return( ret );
 }

/*--------------------------------------------------------------------------*/

int HiGHSMILPSolver::compute( bool changedvars )
{
 lock();  // lock the mutex: this is done again inside MILPSolver::compute,
          // but that's OK since the mutex is recursive

 // process Modification: this is driven by MILPSolver- - - - - - - - - - - -
 if( MILPSolver::compute( changedvars ) != kOK )
  throw( std::runtime_error( "an error occurred in MILPSolver::compute()" ) );

 // HiGHS doesn't actually support MIQP problem
 if( ( int_vars > 0 ) && ( q_obj_val.size() > 0 ) )
  if( ! relax_int_vars ) // we are not relaxing int variables
    throw( std::runtime_error( 
  "HiGHS cannot solve QP models where some of the variables must take integer values" ) );

 // if required, write the problem to file- - - - - - - - - - - - - - - - - -
 if( ! output_file.empty() ) {
  std::string output_file_lp;
  std::stringstream X(output_file);
  std::getline( X , output_file_lp , '.');
  output_file_lp = output_file_lp.append(".lp");
  Highs_writeModel( highs , output_file_lp.c_str() );
 }

 // the actual call to HiGHS- - - - - - - - - - - - - - - - - - - - - - - - -

 if( int_vars > 0 ) {  // the MIP case- - - - - - - - - - - - - - - - - - - -
  if( ( CutSepPar & 7 ) ||
      ( UpCutOff < Inf< double >() ) || ( LwCutOff > Inf< double >() ) ) {
   // the callback has to be set 
   Highs_setCallback( highs,  & HiGHSMILPSolver_callback, this );

   // Enable basic Callback type during MIP to check UB/LB
   Highs_startCallback( highs , kHighsCallbackMipLogging );
   Highs_startCallback( highs , kHighsCallbackMipInterrupt );
   f_callback_set = true;

   if( CutSepPar & 3 ) // we do user cut separation
    throw( std::runtime_error( 
      "HiGHS still doesn't support user cut separation" ) );
   
   if( CutSepPar & 4 )  // we do lazy constraint separation
    Highs_startCallback( highs , kHighsCallbackMipImprovingSolution );
   }
  else
   if( f_callback_set )  // the callback was set
    f_callback_set = false;
  }

 // Call HiGHS to solve the problem
 if( Highs_run( highs ) == -1 ) {
  // An error happened during the call of HiGHS_run. Notice that this is not
  // an error related to the model and so the model status could be incorrect.
  std::cerr << "WARNING: An unmanaged error occurred during the execution of  " <<
  "HiGHS_run" << std::endl;
 }

 int m_status;
 m_status = Highs_getModelStatus( highs );

 // If an LP is solved and then modified, the original optimal basis is used
 // to provide a starting basis for the modified LP.
 // Sometimes, HiGHS could produce an error when solving from the advanced 
 // basis, returning a model status Unknown. In such cases, we can try to 
 // re-solve the model from scratch.
 if( m_status == kHighsModelStatusUnknown ){
  Highs_clearSolver( highs );

  // Call HiGHS to solve the problem
  if( Highs_run( highs ) == -1 ) {
   // An error happened during the call of HiGHS_run. Notice that this is not
   // an error related to the model and so the model status could be incorrect.
   std::cerr << "WARNING: An unmanaged error occurred during the execution of  " <<
    "HiGHS_run" << std::endl;
   }
 
  m_status = Highs_getModelStatus( highs );
  }
  
 sol_status = decode_model_status( m_status );

 Return_status:
 unlock();  // unlock the mutex
 return( sol_status );

 }  // end( HiGHSMILPSolver::compute )

/*--------------------------------------------------------------------------*/


int HiGHSMILPSolver::decode_model_status( int status )
{
 DEBUG_LOG( "HiGHS_ModelStatus returned " << status << std::endl );

 /* The following are the symbols that may represent the status of
 * a HiGHS solution as returned by Highs_getModelStatus
 * as listed in the Enum section on HiGHS Documentation. */

 switch(status) {
  case( kHighsModelStatusNotset ):
   // The model status has not been set.
  case( kHighsModelStatusLoadError ):
   // There has been an error in the load of the model.
  case( kHighsModelStatusModelError ):
   // There is an error in the model.
  case( kHighsModelStatusPresolveError ):
   // There has been an error in the presolve phase.
  case( kHighsModelStatusSolveError ):
   // There has been an error when solving the model.
   return( kError );
  case( kHighsModelStatusModelEmpty ):
   // The model is empty.
   return( kError );
  case( kHighsModelStatusOptimal ):
   // The model has been solved to optimality.
   return( kOK );
  case( kHighsModelStatusInfeasible ):
   // The model is infeasible.
   return( kInfeasible );
  case( kHighsModelStatusUnboundedOrInfeasible ):
   // The model is unbounded or infeasible.
  case( kHighsModelStatusUnbounded ):
   // The model is unbounded.
   return( kUnbounded );
  case( kHighsModelStatusObjectiveBound ):
   // The bound on the model objective value has been reached.
  case( kHighsModelStatusObjectiveTarget ):
   // The target value for the model objective has been reached.
   return( kOK );
  case( kHighsModelStatusTimeLimit ):
   // The run time limit has been reached.
   return( kStopTime );
  case( kHighsModelStatusIterationLimit ):
   // The iteration limit has been reached.
   return( kStopIter );
  case( kHighsModelStatusUnknown ):
   //  The model status is unknown.
   return( kError );
  case( kHighsModelStatusSolutionLimit ):
   // The MIP solver has reached the limit on the number of LPs solved.
   return( kOK );
  default:;
 }

 throw( std::runtime_error( "HiGHS_ModelStatus returned unknown status " +
			    std::to_string( status ) ) );
 }

/*--------------------------------------------------------------------------*/

Solver::OFValue HiGHSMILPSolver::get_lb( void )
{
 OFValue lower_bound = 0;
 int sense;
 Highs_getObjectiveSense( highs , &sense );

 switch( sense ) {
  case( kHighsObjSenseMinimize ):  // Minimization problem - - - - - - - - - -
   switch( sol_status ) {
    case( kUnbounded ):  lower_bound = -Inf< OFValue >(); break;
    case( kInfeasible ): lower_bound = Inf< OFValue >();  break;
    case( kOK ):
    case( kStopIter ):
    case( kStopTime ):
    case( kUnEval ): // Sometimes it could be asked also during the computation

      // TODO: Here we should retrieve the bound
      lower_bound = Highs_getObjectiveValue( highs );
      lower_bound += constant_value;
      break;

    default:
     // If HiGHS does not state that an optimal solution has been found
     // then we do not have enough information to provide a "good" lower bound
     // (the problem may be unbounded and HiGHS has not detected it yet).
     // Therefore, in this case, the lower bound should be -Inf.
     lower_bound = -Inf< OFValue >();
    }
   break;

  case( kHighsObjSenseMaximize ):  // Maximization problem - - - - - - - - - -
   switch( sol_status ) {
    case( kUnbounded ):  lower_bound = Inf< OFValue >(); break;
    case( kInfeasible ): lower_bound = -Inf< OFValue >(); break;

    // if the algorithm has been stopped, the bound only exists if a
    // feasible solution has been generated
    case( kStopIter ):
    case( kStopTime ):
    case( kUnEval ): // Sometimes it could be asked also during the computation
     if( ! has_var_solution() ) {
      lower_bound = - Inf< OFValue >();
      break;
      }

    case( kOK ):
     lower_bound = Highs_getObjectiveValue( highs );
     lower_bound += constant_value;
     break;

    default:  // Same as above
     lower_bound = -Inf< OFValue >();
    }
   break;

  default: throw( std::runtime_error( "Objective type not yet defined" ) );
  }

 return( lower_bound );
 }

/*--------------------------------------------------------------------------*/

Solver::OFValue HiGHSMILPSolver::get_ub( void )
{
 OFValue upper_bound = 0;
 int sense;
 Highs_getObjectiveSense( highs , &sense );

 switch( sense ) {
  case( kHighsObjSenseMinimize ):  // Minimization problem - - - - - - - - - -
   switch( sol_status ) {
    case( kUnbounded ):  upper_bound = -Inf< OFValue >(); break;
    case( kInfeasible ): upper_bound = Inf< OFValue >(); break;

    // if the algorithm has been stopped, the bound only exists if a
    // feasible solution has been generated
    case( kStopIter ):
    case( kStopTime ):
    case( kUnEval ): // Sometimes it could be asked also during the computation
     if( ! has_var_solution() ) {
      upper_bound = Inf< OFValue >();
      break;
      }

    case( kOK ):
     upper_bound = Highs_getObjectiveValue( highs );
     upper_bound += constant_value;
     break;

    default:
     // If HiGHS does not state that an optimal solution has been found
     // then we do not have enough information to provide a "good" upper bound
     // (the problem may be unbounded and HiGHS has not detected it yet).
     // Therefore, in this case, the upper bound should be +Inf.
     upper_bound = Inf< OFValue >();
    }
   break;

  case( kHighsObjSenseMaximize ):  // Maximization problem - - - - - - - - - -
   switch( sol_status ) {
    case( kUnbounded ):  upper_bound = Inf< OFValue >(); break;
    case( kInfeasible ): upper_bound = -Inf< OFValue >(); break;

    case( kOK ):
    case( kStopIter ):
    case( kStopTime ):
    case( kUnEval ): // Sometimes it could be asked also during the computation

     // TODO: Here we should retrieve the bound
     upper_bound = Highs_getObjectiveValue( highs );
     upper_bound += constant_value;
     break;

    default:  // Same as above
     upper_bound = Inf< OFValue >();
    }
   break;

  // Sense not defined
  default: throw( std::runtime_error( "Objective type not yet defined" ) );
  }

 return( upper_bound );
 }

/*--------------------------------------------------------------------------*/

bool HiGHSMILPSolver::has_var_solution( void )
{
 int sol_status , status;
 status = Highs_getIntInfoValue( highs , "primal_solution_status" ,
				 &sol_status );

 if( status == kHighsStatusError )
  throw( std::runtime_error( 
  "An error occurred in getting primal_solution_status with Highs_getIntInfoValue" ) );

 if( sol_status == kHighsSolutionStatusFeasible )  // the solution is feasible
  return( true );
 else  // there is no solution information or the solution is not feasible
  return( false );
 }

/*--------------------------------------------------------------------------*/

Solver::OFValue HiGHSMILPSolver::get_var_value( void )
{
 int objsense;
 Highs_getObjectiveSense( highs , & objsense );

 switch( objsense ) {
  case( kHighsObjSenseMinimize ): return( get_ub() );
  case( kHighsObjSenseMaximize ): return( get_lb() );
  default: throw( std::runtime_error( "Objective type not yet defined" ) );
  }
 }

/*--------------------------------------------------------------------------*/

void HiGHSMILPSolver::get_var_solution( Configuration * solc )
{
 std::vector< double > col_value( numcols );
 std::vector< double > col_dual( numcols );
 std::vector< double > row_value( numrows );
 std::vector< double > row_dual( numrows );

 int status;
 status = Highs_getSolution( highs , col_value.data() , col_dual.data() ,
                            row_value.data() , row_dual.data() );

 if( status == kHighsStatusError )
  throw( std::runtime_error( "An error occurred in Highs_getSolution()" ) );

 MILPSolver::write_var_solution( col_value );
 }

/*--------------------------------------------------------------------------*/

bool HiGHSMILPSolver::has_var_direction( void )
{
 int has_primal_ray = 0;
 
 Highs_getPrimalRay( highs , & has_primal_ray , nullptr );
 return( bool( has_primal_ray ) );
}

/*--------------------------------------------------------------------------*/

void HiGHSMILPSolver::get_var_direction( Configuration * dirc )
{
 int has_primal_ray = 0;
 std::vector< double > col_value( numcols );
  
 int status;
 status = Highs_getPrimalRay( highs , & has_primal_ray , col_value.data() );

 if( status == kHighsStatusError )
  throw( std::runtime_error( "An error occurred in Highs_getSolution()" ) );

 MILPSolver::write_var_solution( col_value );
}

/*--------------------------------------------------------------------------*/

bool HiGHSMILPSolver::has_dual_solution( void )
{ 
 int dual_sol_status;
 Highs_getIntInfoValue( highs , "dual_solution_status" , & dual_sol_status );
 
 switch( dual_sol_status ) {
  case( kHighsSolutionStatusNone): // There is no solution information
    return( false );
  case( kHighsSolutionStatusInfeasible ): // The solution is not feasible.
   return( false );
  case( kHighsSolutionStatusFeasible ): // The solution is feasible.
   return( true );
  default:
   throw( std::runtime_error( "dual_solution_status not recognized" ) );
 }

 return( false );
 }

/*--------------------------------------------------------------------------*/

bool HiGHSMILPSolver::is_dual_feasible( void )
{
 int dual_sol_status;
 Highs_getIntInfoValue( highs , "dual_solution_status" , & dual_sol_status );
 
 switch( dual_sol_status ) {
  case( kHighsSolutionStatusNone): // There is no solution information
    return( false );
  case( kHighsSolutionStatusInfeasible ): // The solution is not feasible.
   return( false );
  case( kHighsSolutionStatusFeasible ): // The solution is feasible.
   return( true );
  default:
   throw( std::runtime_error( "dual_solution_status not recognized" ) );
 }

 return( false );
 }

/*--------------------------------------------------------------------------*/

void HiGHSMILPSolver::get_dual_solution( Configuration * solc )
{
 std::vector< double > col_value( numcols );
 std::vector< double > col_dual( numcols ); // == dj
 std::vector< double > row_value( numrows );
 std::vector< double > row_dual( numrows );  // == pi

 int status;
 status = Highs_getSolution( highs , col_value.data() , col_dual.data() ,
                            row_value.data() , row_dual.data() );

 if( status == kHighsStatusError )
  throw( std::runtime_error( "An error occurred in Highs_getSolution()" ) );

 // Call the method of the base class
 MILPSolver::write_dual_solution( row_dual , col_dual );
}

/*--------------------------------------------------------------------------*/

bool HiGHSMILPSolver::has_dual_direction( void )
{
 int has_dual_ray;
 Highs_getDualRay( highs , & has_dual_ray , nullptr );
 return( bool( has_dual_ray ) );
}

/*--------------------------------------------------------------------------*/

void HiGHSMILPSolver::get_dual_direction( Configuration * dirc )
{
 std::vector< double > y( numrows , 0 );
 std::vector< double > dj( numcols , 0 );

 int has_dual_ray;

 // We are searching a Farkas certificate y so that:
 // y' * A * x >= y' * b
 //   If it is a <= constraint then y[ i ] <= 0 holds;
 //   If it is a >= constraint then y[ i ] >= 0 holds.

 if( Highs_getDualRay( highs , & has_dual_ray , y.data() ) == kHighsStatusError )
  throw( std::runtime_error( "an error occurred in getting Farkas certificate" ) );

 // reverse the sign of y due to Gurobi approach
 //for( auto i = y.begin() ; i != y.end() ; ++i  )
  //*i = -*i;

 if( Highs_getSolution( highs , NULL , NULL , dj.data() , NULL ) == kHighsStatusError )
  throw( std::runtime_error( "Unable to get reduced costs with Highs_getSolution") );

 // Call the method of the base class
 MILPSolver::write_dual_solution( y , dj );
}

/*--------------------------------------------------------------------------*/

void HiGHSMILPSolver::write_lp( const std::string & filename )
{
 std::string output_file_lp;
 std::stringstream X(output_file);
 std::getline( X , output_file_lp , '.');
 output_file_lp = output_file_lp.append(".lp");
 Highs_writeModel( highs , output_file_lp.c_str() );
 }

/*--------------------------------------------------------------------------*/

int HiGHSMILPSolver::get_nodes( void ) const
{
 int64_t nodecnt;
 Highs_getInt64InfoValue( highs , "mip_node_count" , & nodecnt );
 return( (int)nodecnt );
 }

 /*--------------------------------------------------------------------------*/
/*-------------------- PROTECTED FIELDS OF THE CLASS -----------------------*/
/*--------------------------------------------------------------------------*/

void HiGHSMILPSolver::var_modification( const VariableMod * mod )
{
 bool was_mip = int_vars > 0;  // was a MIP before the change

 // call the method of MILPSolver to update dictionaries (and int_vars)
 MILPSolver::var_modification( mod );

 bool is_mip = int_vars > 0;  // is a MIP after the change

 auto var = static_cast< const ColVariable * >( mod->variable() );
 auto idx = index_of_variable( var );
 if( idx == Inf< int >() )  // the Variable is not (yet) there (?)
  return;                   // nothing to do
        
 // react to changes in the integrality - - - - - - - - - - - - - - - - - - -
 if( ColVariable::is_integer( mod->old_state() ) !=
     ColVariable::is_integer( mod->new_state() ) ) {

  // construct new variable type
  char new_ctype;
  double lb, ub; // if is binary, we have to set lhs and rhs to [0,1]
  if( var->is_integer() && ( ! relax_int_vars ) )
   // Integer or Binary
   new_ctype = kHighsVarTypeInteger;
  else
   new_ctype = kHighsVarTypeContinuous;  // Continuous

  Highs_changeColIntegrality( highs , idx , new_ctype );
  }   // end( if( new integrality != old integrality ) )

 // react to fix / unfix- - - - - - - - - - - - - - - - - - - - - - - - - - -
 if( Variable::is_fixed( mod->old_state() ) !=
     Variable::is_fixed( mod->new_state() ) ) {
  if( Variable::is_fixed( mod->new_state() ) ) {  // fix the variable
    std::array< double , 1 > bd = { var->get_value() };
    Highs_changeColBounds( highs , idx , bd[ 0 ] , bd[ 0 ] );
    }
   else {                                         // un fix the variable
    auto bd = HiGHSMILPSolver::get_problem_bounds( *var );
    Highs_changeColBounds( highs , idx , bd[ 0 ] , bd[ 1 ] );
    }
  }
 }  // end( HiGHSMILPSolver::var_modification )

/*--------------------------------------------------------------------------*/

void HiGHSMILPSolver::objective_modification( const ObjectiveMod * mod )
{
 // call the method of MILPSolver to update sense direction
 MILPSolver::objective_modification( mod );

 /* ObjectiveMod class does not include any modification types except
  * for eSetMin and eSetMax.
  * To change OF coefficients, a FunctionMod must be used. */

 switch( mod->type() ) {
  case( ObjectiveMod::eSetMin ):
   Highs_changeObjectiveSense( highs , kHighsObjSenseMinimize );
   break;
  case( ObjectiveMod::eSetMax ):
   Highs_changeObjectiveSense( highs , kHighsObjSenseMaximize );
   break;
  default: throw( std::invalid_argument( "Invalid type of ObjectiveMod" ) );
  }
 }

/*--------------------------------------------------------------------------*/

void HiGHSMILPSolver::const_modification( const ConstraintMod * mod )
{
 // no point in calling the method of MILPSolver, as it does nothing
 // MILPSolver::const_modification( mod );

 /* To change the coefficients, a FunctionMod must be used. */

 auto * con = dynamic_cast< FRowConstraint * >( mod->constraint() );
 if( ! con )  // this should not happen
  return;     // but in case, nothing to do

 int index = index_of_constraint( con );
 if( index == Inf< int >() )  // the FRowConstraint is not (yet?) there
  return;                     // nothing to do
 double rhs, lhs;
 double rngval;

 RowConstraint::RHSValue con_lhs = NAN;
 RowConstraint::RHSValue con_rhs = NAN;

 switch( mod->type() ) {
  case( ConstraintMod::eRelaxConst ):
    // In order to relax the constraint all we do is set both
    // lhs and rhs to infinity

    rhs = kHighsInf;
    lhs = -kHighsInf;
    Highs_changeRowBounds( highs , index , lhs , rhs );
    break;

  case( ConstraintMod::eEnforceConst ):
  case( RowConstraintMod::eChgLHS ):
  case( RowConstraintMod::eChgRHS ):
  case( RowConstraintMod::eChgBTS ):
   // In order to enforce a relaxed constraint all we need to do is
   // reverse the process of relaxing it, by changing the lhs and
   // the rhs back to the original form of the constraint.
   // Moreover, for the way the LP vectors are built, handling
   // LHS/RHS/BTS cases separately is not worth it.

   con_lhs = con->get_lhs();
   con_rhs = con->get_rhs();

   if( con_lhs == con_rhs ) {
    lhs = con_rhs;
    rhs = con_rhs;
    }
   else
    if( con_lhs == -Inf< double >() ) {
     lhs = -kHighsInf;
     rhs = con_rhs;
     }
    else
     if( con_rhs == Inf< double >() ) {
      lhs = con_lhs;
      rhs = kHighsInf;
      }
     else {
      rhs = con_rhs;
      lhs = con_lhs;
      }
   
   Highs_changeRowBounds( highs , index , lhs , rhs );
   
   break;

  default:
   throw( std::invalid_argument( "Invalid type of ConstraintMod" ) );
  }
 }  // end( HiGHSMILPSolver::const_modification )

/*--------------------------------------------------------------------------*/

void HiGHSMILPSolver::bound_modification( const OneVarConstraintMod * mod )
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

 // fixed variables are implemented in HiGHSMILPSolver by changing the bounds;
 // therefore, actual changes of the bounds are ignored here. note that we
 // are assuming the new bounds do not make the fixed value of the variable
 // unfeasible, as this will make the whole problem unfeasible
 // TODO: check this
 if( var->is_fixed() )
  return;

 auto vi = index_of_variable( var );
 if( vi == Inf< int >() )  // the ColVariable has been removed
  return;                  // is strange, but there is nothing to do

 // There is no need of switch based on the type of the modification,
 // because in HiGHSMILPSolver we have always to set both lhs and rhs.
 // However, we do this in order to check if the modification is an
 // actual recognized one.

 auto bd = HiGHSMILPSolver::get_problem_bounds( *var );
 switch( mod->type() ) {
  case( RowConstraintMod::eChgLHS ):
  case( RowConstraintMod::eChgRHS ): 
  case( RowConstraintMod::eChgBTS ): 
      Highs_changeColBounds( highs , vi , bd[ 0 ] , bd[ 1 ] );
      break;

  default:
   throw( std::invalid_argument( "Invalid type of OneVarConstraintMod" ) );
  }
 }  // end( HiGHSMILPSolver::bound_modification )

/*--------------------------------------------------------------------------*/

void HiGHSMILPSolver::objective_function_modification( const FunctionMod * mod )
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

   // we exploit the delta() vector of C05FunctionModLin, giving the difference
   // between the new and the old value of the linear coefficient, to update
   // the objective values without having to recompute them: since they are
   // (potentially) a sum of terms, recomputing them would require fetching
   // back all the terms, while the delta() can just be applied to the sum
   for( Block::Index i = 0 ; i < modl->vars().size() ; ++i ) {
    auto var = static_cast< const ColVariable * >( modl->vars()[ i ] );

    if( auto idx = *(idxit++) ; idx < Inf< Index >() ) {
     int vidx = index_of_variable( var );
     *(cidxit++) = vidx;
      
     // Retrieve old coefficient
     double oldval;
     int num_col, num_nz;
     Highs_getColsByRange( highs , vidx , vidx , &num_col, &oldval, NULL, NULL,
      &num_nz , NULL , NULL , NULL );

     // Update new coefficient
     *(nvit++) = oldval + modl->delta()[ i ];
     }
   }

   auto nsz = std::distance( nval.begin() , nvit );
   cidx.resize( nsz );
   nval.resize( nsz );

   Highs_changeColsCostBySet( highs , cidx.size() , cidx.data() , nval.data() );
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

   // we exploit the delta() vector of C05FunctionModLin, giving the difference
   // between the new and the old value of the linear coefficient, to update
   // the objective values without having to recompute them: since they are
   // (potentially) a sum of terms, recomputing them would require fetching
   // back all the terms, while the delta() can just be applied to the sum
   for( Block::Index i = 0 ; i < modl->vars().size() ; ++i ) {
    auto var = static_cast< const ColVariable * >( modl->vars()[ i ] );

    if( auto idx = *(idxit++) ; idx < Inf< Index >() ) {
     int vidx = index_of_variable( var );
     *(cidxit++) = vidx;
      
     // Retrieve old coefficient
     double oldval;
     int num_col, num_nz;
     Highs_getColsByRange( highs , vidx , vidx , &num_col, &oldval, NULL, NULL,
      &num_nz , NULL , NULL , NULL );

     // Update new coefficient
     *(nvit++) = oldval + modl->delta()[ i ];
     }
   }

   Highs_changeColsCostBySet( highs , cidx.size() , cidx.data() , nval.data() );
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

  auto qf = dynamic_cast< const QuadFunction * >( f );
  auto dqf = dynamic_cast< const DQuadFunction * >( f );
  if( ( ! qf ) && ( ! dqf ) )
    throw( std::logic_error(
		       "unexpected *C05FunctionMod* from Linear Objective" ) );

  // Select correct quadratic function
  auto fqf = ( qf ) ? qf : dqf;

  if( auto modlr = dynamic_cast< const DQuadFunctionModRngd * >( modl ) ) {
   // we exploit the delta() vector of DQuadFunctionModRngd, giving the difference
   // between the new and the old value of both linear and quadratic coefficient,
   // to update the objective values without having to recompute them: since they are
   // (potentially) a sum of terms, recomputing them would require fetching
   // back all the terms, while the delta() can just be applied to the sum
   Subset idxs = fqf->map_index( modlr->vars() , modlr->range() );
   c_Vec_p_Var * vars = & modlr->vars();
   c_v_coeff_pair * delta_coeff = & modlr->delta();

   std::vector< double > nval( idxs.size() );
   std::vector< int > cidx( idxs.size() );
   auto nvit = nval.begin();
   auto idxit = idxs.begin();
   auto cidxit = cidx.begin();
   auto dcoeffit = delta_coeff->begin();

   // In HiGHS we can pass quadratic coefficients to the model only by providing
   // the whole Hessian matrix. Thus, it is important to retrieve the old matrix 
   // and fill it with the new values.
   int nnz_old_hessian = q_obj_val.size();
   int nnz_new_hessian = nnz_old_hessian;

   for( auto v : *vars ){
    if( auto idx = *(idxit++) ; idx < Inf< Index >() ){
     auto cidx = index_of_variable( static_cast< const ColVariable * >( v ) );
     *(cidxit++) = cidx;
      
     // Retrieve old linear coefficient
     double oldlinval;
     int num_col, num_nz;
     Highs_getColsByRange( highs , cidx , cidx , &num_col, &oldlinval, NULL, NULL,
      &num_nz , NULL , NULL , NULL );
     // Update new linear coefficient
     *(nvit++) = oldlinval + std::get< 0 >( *dcoeffit );

     // Retrieve old quadratic coefficient.
     // Note: In HiGHS we store the Hessian matrix in sparse column form using 
     // the 3 vectors described in HiGHSMILPSolver.h. Moreover, the Hessian matrix
     // should be kept in upper triangular form. Thus, the diagonal term can be
     // found at the first position correponding to the coefficients related to
     // column cidx.
     if( std::get< 1 >( *dcoeffit ) != 0 ){
      if( q_obj_begin.empty() ){
        // This could happen when we started with a LinearFunction on the objective,
        // and now we are trying to add quadratic terms for already active variables.
        q_obj_begin.resize( numcols , 0 );
        q_obj_val.push_back( 2 * std::get< 1 >( *dcoeffit ) );
        q_obj_ind.push_back( cidx );
        ++nnz_new_hessian;
        
        // Update q_obj_begin for all the variables having a greater index
        auto it_q_begin = std::next( q_obj_begin.begin() , cidx + 1 );
        while( it_q_begin != q_obj_begin.end() ) {
        ++( *it_q_begin );
        ++it_q_begin;
        }
      }
      else{
        int var_coeff_begin = q_obj_begin[ cidx ];
        auto delta_q_coeff = 2 * std::get< 1 >( *dcoeffit ); 
        if( ( q_obj_ind[ var_coeff_begin ] != cidx ) && ( delta_q_coeff != 0 ) ) {
          // no quadratic coefficient was already set for the diagonal term and
          // the quadratic coefficient is nonzero
          ++nnz_new_hessian;
          auto it_q_begin = std::next( q_obj_begin.begin() , cidx + 1 );
          
          // Update q_obj_begin for all the variables having a greater index
          while( it_q_begin != q_obj_begin.end() ) {
          ++( *it_q_begin );
          ++it_q_begin;
          }

          // Insert new term in the other vectors
          q_obj_val.insert( std::next( q_obj_val.begin() , var_coeff_begin ) ,
                        delta_q_coeff );
          q_obj_ind.insert( std::next( q_obj_ind.begin() , var_coeff_begin ) , cidx );
        }
        else if( q_obj_ind[ var_coeff_begin ] == cidx ) {
        // there was already a value in the hessian diagonal for the variable cidx
          if( q_obj_val[ var_coeff_begin ] == - delta_q_coeff ) {
          // the delta value is the opposite of the old one (i.e. we are removing
          // the term from the hessian matrix)
          --nnz_new_hessian;
          auto it_q_begin = std::next( q_obj_begin.begin() , cidx + 1 );
          // Update q_obj_begin
          while( it_q_begin != q_obj_begin.end() ) {
            --( *it_q_begin );
            ++it_q_begin;
            }
            q_obj_val.erase( std::next( q_obj_val.begin() , var_coeff_begin ) );
            q_obj_ind.erase( std::next( q_obj_ind.begin() , var_coeff_begin ) );
          }
          else{ // just update the array q_obj_val with the new value
            q_obj_val[ var_coeff_begin ] += delta_q_coeff;
          }
        }
      }
     }
    }
    dcoeffit++;
   }

   // Update Hessian
   int status = 0;
   if( !q_obj_begin.empty() )
    status = Highs_passHessian( highs , numcols , nnz_new_hessian ,
                      kHighsHessianFormatTriangular , q_obj_begin.data() ,
                      q_obj_ind.data() , q_obj_val.data()
                      );
                  
   if( status == kHighsStatusError )
      throw( std::runtime_error( 
      "Redefinition of HiGHS hessian matrix returned with kHighsStatus " + 
          std::to_string( status ) ) );

   auto nsz = std::distance( nval.begin() , nvit );
   cidx.resize( nsz );
   nval.resize( nsz );

   // Update linear coefficients
   Highs_changeColsCostBySet( highs , cidx.size() , cidx.data() , nval.data() );
   return;
   }
  else if( auto modls = dynamic_cast< const DQuadFunctionModSbst * >( modl ) ){
   // we exploit the delta() vector of DQuadFunctionModSbst, giving the difference
   // between the new and the old value of both linear and quadratic coefficient,
   // to update the objective values without having to recompute them: since they are
   // (potentially) a sum of terms, recomputing them would require fetching
   // back all the terms, while the delta() can just be applied to the sum
   Subset idxs = fqf->map_index( modls->vars() , modls->subset() );
   c_Vec_p_Var * vars = & modls->vars();
   c_v_coeff_pair * delta_coeff = & modls->delta();

   std::vector< double > nval( idxs.size() );
   std::vector< int > cidx( idxs.size() );
   auto nvit = nval.begin();
   auto idxit = idxs.begin();
   auto cidxit = cidx.begin();
   auto dcoeffit = delta_coeff->begin();

   // In HiGHS we can pass quadratic coefficients to the model only by providing
   // the whole Hessian matrix. Thus, it is important to retrieve the old matrix 
   // and fill it with the new values.
   int nnz_old_hessian = q_obj_val.size();
   int nnz_new_hessian = nnz_old_hessian;

   for( auto v : *vars ){
    if( auto idx = *(idxit++) ; idx < Inf< Index >() ){
     auto cidx = index_of_variable( static_cast< const ColVariable * >( v ) );
     *(cidxit++) = cidx;
      
     // Retrieve old linear coefficient
     double oldlinval;
     int num_col, num_nz;
     Highs_getColsByRange( highs , cidx , cidx , &num_col, &oldlinval, NULL, NULL,
      &num_nz , NULL , NULL , NULL );
     // Update new linear coefficient
     *(nvit++) = oldlinval + std::get< 0 >( *dcoeffit );

     // Retrieve old quadratic coefficient.
     // Note: In HiGHS we store the Hessian matrix in sparse column form using 
     // the 3 vectors described in HiGHSMILPSolver.h. Moreover, the Hessian matrix
     // should be kept in upper triangular form. Thus, the diagonal term can be
     // found at the first position correponding to the coefficients related to
     // column cidx.
     if( std::get< 1 >( *dcoeffit ) != 0 ){
      if( q_obj_begin.empty() ){
        // This could happen when we started with a LinearFunction on the objective,
        // and now we are trying to add quadratic terms for already active variables.
        q_obj_begin.resize( numcols , 0 );
        q_obj_val.push_back( 2 * std::get< 1 >( *dcoeffit ) );
        q_obj_ind.push_back( cidx );
        ++nnz_new_hessian;
        
        // Update q_obj_begin for all the variables having a greater index
        auto it_q_begin = std::next( q_obj_begin.begin() , cidx + 1 );
        while( it_q_begin != q_obj_begin.end() ) {
        ++( *it_q_begin );
        ++it_q_begin;
        }
      }
      else{
        int var_coeff_begin = q_obj_begin[ cidx ];
        auto delta_q_coeff = 2 * std::get< 1 >( *dcoeffit );
        if( ( q_obj_ind[ var_coeff_begin ] != cidx ) && ( delta_q_coeff != 0 ) ) {
          // no quadratic coefficient was already set for the diagonal term and
          // the quadratic coefficient is nonzero
          ++nnz_new_hessian;
          auto it_q_begin = std::next( q_obj_begin.begin() , cidx + 1 );
          
          // Update q_obj_begin for all the variables having a greater index
          while( it_q_begin != q_obj_begin.end() ) {
          ++( *it_q_begin );
          ++it_q_begin;
          }

          // Insert new term in the other vectors
          q_obj_val.insert( std::next( q_obj_val.begin() , var_coeff_begin ) ,
                        delta_q_coeff );
          q_obj_ind.insert( std::next( q_obj_ind.begin() , var_coeff_begin ) , cidx );
        }
        else if( q_obj_ind[ var_coeff_begin ] == cidx ) {
        // there was already a value in the hessian diagonal for the variable cidx
          if( q_obj_val[ var_coeff_begin ] == - delta_q_coeff ) {
          // the delta value is the opposite of the old one (i.e. we are removing
          // the term from the hessian matrix)
          --nnz_new_hessian;
          auto it_q_begin = std::next( q_obj_begin.begin() , cidx + 1 );
          // Update q_obj_begin
          while( it_q_begin != q_obj_begin.end() ) {
            --( *it_q_begin );
            ++it_q_begin;
            }
            q_obj_val.erase( std::next( q_obj_val.begin() , var_coeff_begin ) );
            q_obj_ind.erase( std::next( q_obj_ind.begin() , var_coeff_begin ) );
          }
          else{ // just update the array q_obj_val with the new value
            q_obj_val[ var_coeff_begin ] += delta_q_coeff;
          }
        }
      }
     }
    }
    dcoeffit++;
   }
   
   // Update Hessian
   int status = 0;
   if( !q_obj_begin.empty() )
    status = Highs_passHessian( highs , numcols , nnz_new_hessian ,
                      kHighsHessianFormatTriangular , q_obj_begin.data() ,
                      q_obj_ind.data() , q_obj_val.data()
                      );
                  
   if( status == kHighsStatusError )
      throw( std::runtime_error( 
      "Redefinition of HiGHS hessian matrix returned with kHighsStatus " + 
          std::to_string( status ) ) );

   auto nsz = std::distance( nval.begin() , nvit );
   cidx.resize( nsz );
   nval.resize( nsz );

   // Update linear coefficients
   Highs_changeColsCostBySet( highs , cidx.size() , cidx.data() , nval.data() );
   return;
   }
  else if( auto modlq = dynamic_cast< const QuadFunctionModSbst * >( modl ) ) {
   // we exploit the delta() vector of QuadFunctionModSbst, giving the difference
   // between the new and the old value of both linear and quadratic coefficient,
   // to update the objective values without having to recompute them: since they are
   // (potentially) a sum of terms, recomputing them would require fetching
   // back all the terms, while the delta() can just be applied to the sum.
   // NOTE: in the actual version of QuadFunction, we expect to recieve one 
   // coefficient at time for each Modification.
   Subset idxs = fqf->map_index( modlq->vars() , modlq->subset() );
   c_Vec_p_Var vars = modlq->vars();
   Coefficient delta_coeff = 2 * modlq->delta();

   // In HiGHS we can pass quadratic coefficients to the model only by providing
   // the whole Hessian matrix. Thus, it is important to retrieve the old matrix 
   // and fill it with the new values.
   int nnz_old_hessian = q_obj_val.size();
   int nnz_new_hessian = nnz_old_hessian;

   if( idxs.size() != 2 )
    throw( std::logic_error(
		       "Expected single coefficient Modification in QuadFunctionModSbst" ) );

   int idx1 = index_of_variable( dynamic_cast< ColVariable * >( vars[ 0 ] ) );
   int idx2 = index_of_variable( dynamic_cast< ColVariable * >( vars[ 1 ] ) );

   // Retrieve old quadratic coefficient.
   // Note: In HiGHS we store the Hessian matrix in sparse column form using 
   // the 3 vectors described in HiGHSMILPSolver.h. 
   int var1_coeff_begin = q_obj_begin[ idx1 ];
   int var1_coeff_end = q_obj_begin[ idx1 + 1 ];
   // Iterator pointing to the right term
   auto qobj_indit = std::find( q_obj_ind.begin() + var1_coeff_begin , 
                                 q_obj_ind.begin() + var1_coeff_end , 
                                 idx2 );

   if( ( qobj_indit == q_obj_ind.begin() + var1_coeff_end ) &&
       ( delta_coeff != 0 ) ) {
    // no quadratic coefficient was already set for the term and
    // the quadratic coefficient is nonzero
    ++nnz_new_hessian;
    auto it_q_begin = std::next( q_obj_begin.begin() , idx1 + 1 );
      
    // Update q_obj_begin for all the variables having a greater index
    while( it_q_begin != q_obj_begin.end() ) {
     ++( *it_q_begin );
     ++it_q_begin;
     }

    // Find the right position to insert the term to keep the UT form
    auto itpos_to_insert = std::upper_bound( q_obj_ind.begin() + var1_coeff_begin , 
                                 q_obj_ind.begin() + var1_coeff_end , 
                                 idx2 );
    int pos_to_insert = std::distance( q_obj_ind.begin() , itpos_to_insert );

    // Insert new term in the other vectors
    q_obj_val.insert( std::next( q_obj_val.begin() , pos_to_insert ) ,
                     delta_coeff );
    q_obj_ind.insert( std::next( q_obj_ind.begin() , pos_to_insert ) , idx2 );
    }
   else if( qobj_indit != q_obj_ind.begin() + var1_coeff_end ) {
    // there was already a value in the hessian diagonal for the variable cidx
    int rpos = std::distance( q_obj_ind.begin() , qobj_indit );
    if( q_obj_val[ rpos ] == - delta_coeff ) {
      // the delta value is the opposite of the old one (i.e. we are removing
      // the term from the hessian matrix)
      --nnz_new_hessian;
      auto it_q_begin = std::next( q_obj_begin.begin() , idx1 + 1 );
      // Update q_obj_begin
      while( it_q_begin != q_obj_begin.end() ) {
       --( *it_q_begin );
       ++it_q_begin;
       }
      q_obj_val.erase( std::next( q_obj_val.begin() , rpos ) );
      q_obj_ind.erase( std::next( q_obj_ind.begin() , rpos ) );
     }
    else{ // just update the array q_obj_val with the new value
        q_obj_val[ rpos ] += delta_coeff;
      }
     }
   // Update Hessian
   int status = Highs_passHessian( highs , numcols , nnz_new_hessian ,
                      kHighsHessianFormatTriangular , q_obj_begin.data() ,
                      q_obj_ind.data() , q_obj_val.data()
                      );
                  
   if( status == kHighsStatusError )
      throw( std::runtime_error( 
      "Redefinition of HiGHS hessian matrix returned with kHighsStatus " + 
          std::to_string( status ) ) );
   return;
  }
  else
    throw( std::logic_error( "unknown type of *QuadFunctionMod*" ) );
 }

 // Fallback method - Update all costs
 // --------------------------------------------------------------------------
 // reload_objective( f );

 }  // end( HiGHSMILPSolver::objective_function_modification )

/*--------------------------------------------------------------------------*/

void HiGHSMILPSolver::constraint_function_modification( const FunctionMod *mod )
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

 auto idxit = idxs.begin();
 auto & cp = lf->get_v_var();

 for( auto v :  modl->vars() )
  if( auto idx = *(idxit++) ; idx < Inf< Index >() ) {
   double value = cp[ idx ].second;
   int var_idx = index_of_variable( static_cast< const ColVariable * >( v ) );
   Highs_changeCoeff( highs, row , var_idx , value );
   }

 // Fallback method - Reload all coefficients
 // --------------------------------------------------------------------------
 // reload_constraint( lf );

 }  // end( HiGHSMILPSolver::constraint_function_modification )

/*--------------------------------------------------------------------------*/

void HiGHSMILPSolver::objective_fvars_modification( const FunctionModVars *mod )
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
 if( ( ! dynamic_cast< const LinearFunctionModVarsAddd * >( mod ) ) &&
     ( ! dynamic_cast< const DQuadFunctionModVarsAddd * >( mod ) ) &&
     ( ! dynamic_cast< const QuadFunctionModVarsAddd * >( mod ) ) &&
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
 // structure of [ HiGHS ]MILPSolver since the Modification are managed
 // strictly in arrival order, they may no longer exist in the model;
 // more to the point, they may no longer be active in the LinearFunction

  if( auto lf = dynamic_cast< const LinearFunction * >( f ) ) {
  // Linear objective function modification
  
  // we exploit the coeff() vector of LinearFunctionModVarsAddd, giving the sum
  // between the new and the old value of the linear coefficient, to update
  // the objective values without having to recompute them: since they are
  // (potentially) a sum of terms, recomputing them would require fetching
  // back all the terms, while the coeff() can just be applied to the sum
  
  for( Block::Index i = 0 ; i < mod->vars().size() ; ++i ) {
    auto var = static_cast< const ColVariable * >( mod->vars()[ i ] );

    if( auto idx = index_of_variable( var ) ; idx < Inf< int >() ) {
      indices.push_back( idx );

      if( mod->added() ) {
        auto modl = dynamic_cast< const SMSpp_di_unipi_it::LinearFunctionModVarsAddd * >( mod );

        // Retrieve old coefficient
        double oldval;
        int num_col, num_nz;
        Highs_getColsByRange( highs , idx , idx , &num_col, &oldval, NULL, NULL,
          &num_nz , NULL , NULL , NULL );
        
        auto cidx = lf->is_active( var );
        values.push_back( cidx < nav ? oldval + modl->coeff()[ i ] : oldval );
      }
      else
        values.push_back( 0 );
     }
   }

  Highs_changeColsCostBySet( highs , indices.size() , indices.data() , values.data() );
  return;
 }

 if( auto qf = static_cast< const QuadFunction * >( f ) ) {
  // Quadratic objective function modification

  // In HiGHS we can pass quadratic coefficients to the model only by providing
  // the whole Hessian matrix. Thus, it is important to retrieve the old matrix 
  // and fill it with the new values.
  int nnz_old_hessian = Highs_getHessianNumNz( highs );
  int nnz_new_hessian = nnz_old_hessian;
  
  // Firstly check if we are simply removing variables
  if( ! mod->added() ) {
    for( Block::Index i = 0 ; i < mod->vars().size() ; ++i ) {
      auto var = static_cast< const ColVariable * >( mod->vars()[ i ] );

      if( auto idx = index_of_variable( var ) ; idx < Inf< int >() ) {
       // Retrieve old quadratic coefficient.
       // Note: In HiGHS we store the Hessian matrix in sparse column form using 
       // the 3 vectors described in HiGHSMILPSolver.h. Moreover, the Hessian matrix
       // should be kept in upper triangular form. Thus, the diagonal term can be
       // found at the first position correponding to the coefficients related to
       // column cidx.
       int var_coeff_begin = q_obj_begin[ idx ];
       if( q_obj_ind[ var_coeff_begin ] == idx ) {
        // There was a diagonal term set for variable of index idx
        --nnz_new_hessian;
        auto it_q_begin = std::next( q_obj_begin.begin() , idx + 1 );
        // Update q_obj_begin
        while( it_q_begin != q_obj_begin.end() ) {
         --( *it_q_begin );
         ++it_q_begin;
        }
        q_obj_val.erase( std::next( q_obj_val.begin() , var_coeff_begin ) );
        q_obj_ind.erase( std::next( q_obj_ind.begin() , var_coeff_begin ) );
       }

        indices.push_back( idx );
        values.push_back( 0 );
      }
      // TODO: BUilt a specific Modification to include non diagonal terms
      // to be set to 0.
    }
    // Update all linear coefficient at ones
    Highs_changeColsCostBySet( highs , indices.size() , indices.data() , values.data() );

    // Update all quadratic coefficient at ones if something changed
    if( nnz_old_hessian != nnz_new_hessian ) {
      int status = Highs_passHessian( highs , numcols , nnz_new_hessian ,
                        kHighsHessianFormatTriangular , q_obj_begin.data() ,
                        q_obj_ind.data() , q_obj_val.data()
                        );
                    
      if( status == kHighsStatusError )
        throw( std::runtime_error( 
          "Redefinition of HiGHS hessian matrix returned with kHighsStatus " + 
            std::to_string( status ) ) );
    }             
    return;
  }

  auto modq = static_cast< const SMSpp_di_unipi_it::QuadFunctionModVarsAddd * >( mod );
  if( ! modq )
    // This should never happen
    throw( std::invalid_argument( "Unexpected type of Objective Function Modification" ) );

  // we exploit the od_terms() vector of QuadFunctionModVarsAddd, giving the sum
  // between the new and the old value of the quadratic coefficient, to update 
  // the objective values without having to recompute them: since they are 
  // (potentially) a sum of terms, recomputing them would require fetching back 
  // all the terms, while the coeff() can just be applied to the sum.
  
  for( auto t : modq->od_terms() ) {
    int loc_idx1 = std::get<0>( t );
    int loc_idx2 = std::get<1>( t );

    auto var1 = qf->get_active_var( loc_idx1 );
    auto var2 = qf->get_active_var( loc_idx2 );

    int glob_idx1 = index_of_variable( dynamic_cast< ColVariable * >( var1 ) );
    int glob_idx2 = index_of_variable( dynamic_cast< ColVariable * >( var2 ) );

    // Retrieve old quadratic coefficient.
    // Note: In HiGHS we store the Hessian matrix in sparse column form using 
    // the 3 vectors described in HiGHSMILPSolver.h. 
    int var1_coeff_begin = q_obj_begin[ glob_idx1 ];
    int var1_coeff_end = q_obj_begin[ glob_idx1 + 1 ];
    // Iterator pointing to the right term
    auto qobj_indit = std::find( q_obj_ind.begin() + var1_coeff_begin , 
                                 q_obj_ind.begin() + var1_coeff_end , 
                                 glob_idx2 );

    if( qobj_indit == q_obj_ind.begin() + var1_coeff_end ) {
     // no quadratic coefficient was already set for the term
     ++nnz_new_hessian;
     auto it_q_begin = std::next( q_obj_begin.begin() , glob_idx1 + 1 );
      
     // Update q_obj_begin for all the variables having a greater index
     while( it_q_begin != q_obj_begin.end() ) {
      ++( *it_q_begin );
      ++it_q_begin;
     }

     // Find the right position to insert the term to keep the UT form
     auto itpos_to_insert = std::upper_bound( q_obj_ind.begin() + var1_coeff_begin , 
                                 q_obj_ind.begin() + var1_coeff_end , 
                                 glob_idx2 );
     int pos_to_insert = std::distance( q_obj_ind.begin() , itpos_to_insert );

     // Insert new term in the other vectors
     q_obj_val.insert( std::next( q_obj_val.begin() , pos_to_insert ) ,
                      2 * std::get<2>( t ) );
     q_obj_ind.insert( std::next( q_obj_ind.begin() , pos_to_insert ) , 
                      glob_idx2 );
    }
    else{ 
     // there was already a value in the hessian diagonal for the variable term
     int rpos = std::distance( q_obj_ind.begin() , qobj_indit );
     q_obj_val[ rpos ] += 2 * std::get<2>( t );
    }
  }
  // Update all quadratic coefficient at ones
  int status = Highs_passHessian( highs , numcols , nnz_new_hessian ,
                        kHighsHessianFormatTriangular , q_obj_begin.data() ,
                        q_obj_ind.data() , q_obj_val.data()
                        );
                    
  if( status == kHighsStatusError )
    throw( std::runtime_error( 
      "Redefinition of HiGHS hessian matrix returned with kHighsStatus " + 
            std::to_string( status ) ) );
  // Here we don't need any return, as we know that any QuadFunction
  // derives from a DQuadFunction
 }

 if( auto dqf = dynamic_cast< const DQuadFunction * >( f ) ) {
  // Separable quadratic objective function modification

  // In HiGHS we can pass quadratic coefficients to the model only by providing
  // the whole Hessian matrix. Thus, it is important to retrieve the old matrix 
  // and fill it with the new values.
  int nnz_old_hessian = Highs_getHessianNumNz( highs );
  int nnz_new_hessian = nnz_old_hessian;

  // Firstly check if we are simply removing variables
  if( ! mod->added() ) {
    for( Block::Index i = 0 ; i < mod->vars().size() ; ++i ) {
      auto var = static_cast< const ColVariable * >( mod->vars()[ i ] );

      if( auto idx = index_of_variable( var ) ; idx < Inf< int >() ) {
       // Retrieve old quadratic coefficient.
       // Note: In HiGHS we store the Hessian matrix in sparse column form using 
       // the 3 vectors described in HiGHSMILPSolver.h. Moreover, the Hessian matrix
       // should be kept in upper triangular form. Thus, the diagonal term can be
       // found at the first position correponding to the coefficients related to
       // column cidx.
       int var_coeff_begin = q_obj_begin[ idx ];
       if( q_obj_ind[ var_coeff_begin ] == idx ) {
        // There was a diagonal term set for variable of index idx
        --nnz_new_hessian;
        auto it_q_begin = std::next( q_obj_begin.begin() , idx + 1 );
        // Update q_obj_begin
        while( it_q_begin != q_obj_begin.end() ) {
         --( *it_q_begin );
         ++it_q_begin;
        }
        q_obj_val.erase( std::next( q_obj_val.begin() , var_coeff_begin ) );
        q_obj_ind.erase( std::next( q_obj_ind.begin() , var_coeff_begin ) );
       }

        indices.push_back( idx );
        values.push_back( 0 );
      }
    }
    // Update all linear coefficient at ones
    Highs_changeColsCostBySet( highs , indices.size() , indices.data() , values.data() );

    // Update all quadratic coefficient at ones if something changed
    if( nnz_old_hessian != nnz_new_hessian ) {
      int status = Highs_passHessian( highs , numcols , nnz_new_hessian ,
                        kHighsHessianFormatTriangular , q_obj_begin.data() ,
                        q_obj_ind.data() , q_obj_val.data()
                        );
                    
      if( status == kHighsStatusError )
        throw( std::runtime_error( 
         "Redefinition of HiGHS hessian matrix returned with kHighsStatus " + 
            std::to_string( status ) ) );
    }             
    return;
  }

  auto modq = dynamic_cast< const SMSpp_di_unipi_it::DQuadFunctionModVarsAddd * >( mod );
  if( ! modq )
    // This should never happen
    throw( std::invalid_argument( "Unexpected type of Objective Function Modification" ) );

  // we exploit the coeff() vector of DQuadFunctionModVarsAddd, giving the sum
  // between the new and the old value of both the linear and quadratic
  // coefficient, to update the objective values without having to recompute 
  // them: since they are (potentially) a sum of terms, recomputing them 
  // would require fetching back all the terms, while the coeff()
  // can just be applied to the sum
  
  for( Block::Index i = 0 ; i < mod->vars().size() ; ++i ) {
    auto var = static_cast< const ColVariable * >( mod->vars()[ i ] );
    
    if( auto idx = index_of_variable( var ) ; idx < Inf< int >() ) {
     indices.push_back( idx );

     // Retrieve old linear coefficient
     double lin_oldval;
     int num_col, num_nz;
     Highs_getColsByRange( highs , idx , idx , &num_col, &lin_oldval, NULL, NULL,
        &num_nz , NULL , NULL , NULL );
       
     // add new linear coefficient to final vector
     values.push_back( lin_oldval + modq->coeff()[ i ].first );

     if( modq->coeff()[ i ].second != 0 ) {
      // Search for the quadratic coefficient associated to idx
      // Note: In HiGHS we store the Hessian matrix in sparse column form using 
      // the 3 vectors described in HiGHSMILPSolver.h. Moreover, the Hessian matrix
      // should be kept in upper triangular form. Thus, the diagonal term can be
      // found at the first position correponding to the coefficients related to
      // column idx.
      int var_coeff_begin = q_obj_begin[ idx ];
      if( q_obj_ind[ var_coeff_begin ] == idx ) {
        // Simply update the coefficient
        q_obj_val[ var_coeff_begin ] += 2 * modq->coeff()[ i ].second;
      }
      else{
        // we have actually to add an entry in the Hessian matrix
        ++nnz_new_hessian;
        auto it_q_begin = std::next( q_obj_begin.begin() , idx );
        // Update q_obj_begin
        while( it_q_begin != q_obj_begin.end() ) {
          ++( *it_q_begin );
          ++it_q_begin;
          }
        q_obj_val.insert( std::next( q_obj_val.begin() , var_coeff_begin ) , 
                          2 * modq->coeff()[ i ].second );
        q_obj_ind.insert( std::next( q_obj_ind.begin() , var_coeff_begin ) , idx );
      }
     }
    }
  }
  // Update all linear coefficient at ones
  Highs_changeColsCostBySet( highs , indices.size() , indices.data() , values.data() );

  // Update all quadratic coefficient at ones
  int status = Highs_passHessian( highs , numcols , nnz_new_hessian ,
                        kHighsHessianFormatTriangular , q_obj_begin.data() ,
                        q_obj_ind.data() , q_obj_val.data()
                        );
                    
  if( status == kHighsStatusError )
    throw( std::runtime_error( 
      "Redefinition of HiGHS hessian matrix returned with kHighsStatus " + 
            std::to_string( status ) ) );
  return;
 }

 // This should never happen
 throw( std::invalid_argument( "Unknown type of Objective Function" ) );

 }  // end( HiGHSMILPSolver::objective_fvars_modification )

/*--------------------------------------------------------------------------*/

void HiGHSMILPSolver::constraint_fvars_modification(
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

 auto nav = lf->get_num_active_var();

 // while changing the coefficients, we have to be careful about the fact
 // that Modification are managed asynchronously with the model changes
 // although the added/removed Variable do exist in the internal data
 // structure of [ HiGHS ]MILPSolver since the Modification are managed
 // strictly in arrival order, they may no longer exist in the model;
 // more to the point, they may no longer be active in the LinearFunction

 // get indices and coefficients
 for( auto v : mod->vars() ) {
  auto var = static_cast< const ColVariable * >( v );
  if( auto vidx = index_of_variable( var ) ; vidx < Inf< int >() ) {
   double value = 0;

   if( mod->added() ) {
    auto idx = lf->is_active( var );
    value = idx < nav ? lf->get_coefficient( idx ) : 0;
    }
   Highs_changeCoeff( highs , cidx , vidx , value );
   }
  }
 }  // end( HiGHSMILPSolver::constraint_fvars_modification )

/*----------------------------------------------------------------------------

void HiGHSMILPSolver::dynamic_modification( const BlockModAD * mod )
{
 MILPSolver::dynamic_modification( mod );
 }

----------------------------------------------------------------------------*/

void HiGHSMILPSolver::add_dynamic_constraint( const FRowConstraint * con )
{
 // call the method of MILPSolver to update the dictionaries (only)
 MILPSolver::add_dynamic_constraint( con );

 auto lf = dynamic_cast< const LinearFunction * >( con->get_function() );
 if( ! lf )
  throw( std::invalid_argument( "the FRowConstraint is not linear" ) );

 int nzcnt = lf->get_num_active_var();

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
 double lhs = 0;
 double rhs = 0;
 char sense;

 if( con_lhs == con_rhs ) {
  lhs = rhs;
  rhs = con_rhs;
  }
 else
  if( con_lhs == -Inf< double >() ) {
   lhs = -kHighsInf;
   rhs = con_rhs;
   }
  else
   if( con_rhs == Inf< double >() ) {
    lhs = con_lhs;
    rhs = kHighsInf;
    }
   else {
    lhs = con_lhs;
    rhs = con_rhs;
    }

  Highs_addRow( highs , lhs , rhs , rmatind.size() ,
              rmatind.data() , rmatval.data() );

 }  // end( HiGHSMILPSolver::add_dynamic_constraint )

/*--------------------------------------------------------------------------*/

void HiGHSMILPSolver::add_dynamic_variable( const ColVariable * var )
{
 bool was_mip = int_vars > 0;  // was a MIP before the change

 // call the method of MILPSolver to update dictionaries (and int_vars)
 MILPSolver::add_dynamic_variable( var );

 bool is_mip = int_vars > 0;  // is a MIP after the change

 // get the bounds
 auto bd = HiGHSMILPSolver::get_problem_bounds( *var );

 char new_ctype;  // get the new variable type
 
 if( var->is_integer() && ( ! relax_int_vars ) )
  // Integer or Binary
  new_ctype = kHighsVarTypeInteger;
 else
  new_ctype = kHighsVarTypeContinuous;  // Continuous

 // update the HiGHS problem
 Highs_addVar( highs , bd[0] , bd[1] );
 Highs_changeColIntegrality( highs , numcols - 1 , new_ctype );

 }  // end( HiGHSMILPSolver::add_dynamic_variable )

/*--------------------------------------------------------------------------*/

void HiGHSMILPSolver::add_dynamic_bound( const OneVarConstraint * con )
{
 // no point in calling the method of MILPSolver, as it does nothing
 // MILPSolver::add_dynamic_bound( con );

 auto var = static_cast< const ColVariable * >( con->get_active_var( 0 ) );
 if( ! var )
  throw( std::logic_error( "HiGHSMILPSolver: added a bound on no Variable" ) );

 auto idx = index_of_variable( var );
 if( idx == Inf< int >() )
  throw( std::logic_error( "HiGHSMILPSolver: added a bound on unknown Variable"
			   ) );

 auto bd = HiGHSMILPSolver::get_problem_bounds( *var );

 Highs_changeColBounds( highs , idx , bd[ 0 ] , bd[ 1 ] );
 }

/*--------------------------------------------------------------------------*/

void HiGHSMILPSolver::remove_dynamic_constraint( const FRowConstraint * con )
{
 int index = index_of_dynamic_constraint( con );
 if( index == Inf< int >() )
  throw( std::runtime_error( "Dynamic constraint not found" ) );

 Highs_deleteRowsByRange( highs , index , index );

 // call the method of MILPSolver to update the dictionaries (only)
 MILPSolver::remove_dynamic_constraint( con );
 }

/*--------------------------------------------------------------------------*/

void HiGHSMILPSolver::remove_dynamic_variable( const ColVariable * var )
{
 int index = index_of_dynamic_variable( var );
 if( index == Inf< int >() )
  throw( std::runtime_error( "Dynamic variable not found" ) );

 Highs_deleteColsByRange( highs , index , index );

 // call the method of MILPSolver to update the dictionaries (only)
 MILPSolver::remove_dynamic_variable( var );
 }

/*--------------------------------------------------------------------------*/

void HiGHSMILPSolver::remove_dynamic_bound( const OneVarConstraint * con )
{
 // no point in calling the method of MILPSolver, as it does nothing
 // MILPSolver::remove_dynamic_bound( con );

 // note: this only works because remove_dynamic_constraint[ s ]() do *not*
 //       clear the removed OneVarConstraint, and therefore we can easily
 //       reconstruct which ColVariable it was about
 auto var = static_cast< const ColVariable * >( con->get_active_var( 0 ) );
 if( ! var )  // this should never happen
  return;     // but in case, there is nothing to do

 int idx = index_of_variable( var );
 if( idx == Inf< int >() )  // the ColVariable has been removed
  return;                   // is strange, but there is nothing to do

 auto bd = HiGHSMILPSolver::get_problem_bounds( *var );

 Highs_changeColBounds( highs , idx , bd[ 0 ] , bd[ 1 ] );
 }

/*--------------------------------------------------------------------------*/

int HiGHSMILPSolver::callback( const int callback_type,
  const HighsCallbackDataOut* data_out,
  HighsCallbackDataIn* data_in )
{
 // main switch: depending on callback_type - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 switch( callback_type ) {
  case( kHighsCallbackLogging ): break; /* Ignore Logging callback */
  case( kHighsCallbackSimplexInterrupt ): break; /* Ignore Simplex callback */
  case( kHighsCallbackIpmInterrupt ): break; /* Ignore IPM callback */
  case( kHighsCallbackMipLogging ): break; /* Ignore MIP Logging */
  case( kHighsCallbackMipInterrupt ):
   // Currently in MIP- - - - - - - - - - - - - - - - - - - - - - - -
   // check upper / lower bounds and in case stop

   double solv , bndv;
   solv = data_out->objective_function_value;
   bndv = data_out->mip_primal_bound;

   if( get_objsense() == 1 ) {
    // a minimization problem: solv is upper bound and bndv is lower bound
    if( solv >= 1e+75 )
     solv = Inf< double >();

    if( bndv <= - 1e+75 )
     bndv = - Inf< double >();

    if( ( bndv >= up_cut_off() ) || ( solv <= lw_cut_off() ) ) {
     // Here we should terminate the optimization process but HiGHS currently
     // does not provide such method
     }
    }
   else {
    // a maximization problem: solv is lower bound and bndv is upper bound
    if( solv <= -1e+75 )
     solv = - Inf< double >();

    if( bndv >= 1e+75 )
     bndv = Inf< double >();

    if( ( solv >= up_cut_off() ) || ( bndv <= lw_cut_off() ) ) {
      data_in->user_interrupt = 1; // Force interruption of Solver
      }
    }
   break;
  case( kHighsCallbackMipImprovingSolution ):
   // a feasible solution has been found- - - - - - - - - - - - - - - - - - -
   if( ! ( CutSepPar & 4 ) )  // but we don't do lazy constraint separation
    // strange: this callback type should not be enabled but nothing to do
    break;              

   // this is a critical section where different GUROBI threads may compete
   // for access to the Block: ensure mutual exclusion
   f_callback_mutex.lock();

   // ensure no interference from other threads (except GUROBI ones) by also
   // lock()-ing the Block
   bool owned = f_Block->is_owned_by( f_id );
   if( ( ! owned ) && ( ! f_Block->lock( f_id ) ) )
    throw( std::runtime_error( "Unable to lock the Block" ) );

   // get the feasible solution
   std::vector< double > x( numcols );
   x.assign(data_out->mip_solution, data_out->mip_solution + numcols);

   // write it in the Variable of the Block
   MILPSolver::write_var_solution( x );

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
   if( ! rmatbeg.empty() ) {
    // Here we should add the lazy constraints generated, but this feature
    // is currently under development in HiGHS. TBD
    std::cerr << "WARNING: Detected a new lazy constraint during HiGHS "
      "callback but this feature is not supported yet" << std::endl;
    }
   break;
  } // end( main switch )- - - - - - - - - - - - - - - - - - - - - - - - - -

 return( 0 );
}

/*--------------------------------------------------------------------------*/

void HiGHSMILPSolver::perform_separation( Configuration * cfg ,
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
 // Modification list as no-one has a reason tochange it

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
  auto * lf = dynamic_cast< const LinearFunction * >( con->get_function() );
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

  /* TBD: Modify this based on HiGHS version of lazy constraints.
  if( con_lhs == con_rhs ) {
   sense.push_back( GRB_EQUAL );
   rhs.push_back( con_rhs );
  }
  else if( con_lhs == -Inf< double >() ) {
   sense.push_back( GRB_LESS_EQUAL );
   rhs.push_back( con_rhs );
  }
  else if( con_rhs == Inf< double >() ) {
   sense.push_back( GRB_GREATER_EQUAL );
   rhs.push_back( con_lhs );
  }
  else {
   // kludge: the added constraint is ranged LHS <= lf( x ) <= RHS, but
   // GUROBI does not allow cuts to be ranged: hence, separately add
   // the two constraints lf( x ) >= LHS and lf( x ) <= RHS
   sense.push_back( GRB_GREATER_EQUAL );
   rhs.push_back( con_lhs );
   auto nsz = rmatind.size();
   rmatbeg.push_back( nsz );
   sense.push_back( GRB_LESS_EQUAL );
   rhs.push_back( con_rhs );
   rmatind.resize( nsz + nzcnt );
   std::copy( rmatind.begin() + sz , rmatind.begin() + nsz ,
                      rmatind.begin() + nsz );
   rmatval.resize( nsz + nzcnt );
   std::copy( rmatval.begin() + sz , rmatval.begin() + nsz ,
                      rmatval.begin() + nsz );
  }

  rmatbeg.push_back( rmatind.size() ); */

  }  // end( for each added FRowConstraint )
 }  // end( main loop )
}  // end( HiGHSMILPSolver::perform_separation )

/*--------------------------------------------------------------------------*/

std::string HiGHSMILPSolver::highs_int_par_map( idx_type par ) const
{
 switch( par ) {
  case( intMaxIter ): return( "mip_max_nodes" );
  case( intMaxSol ):  return( "mip_max_improving_sols" );
  case( intLogVerb ): return( "output_flag" );
  case( intMaxThread ): return( "threads" );
  }

 // HiGHS options
 if( ( par >= intFirstHiGHSPar ) && ( par < intLastAlgParHiGHS ) ) {
  return( SMSpp_to_HiGHS_int_pars[ par - intFirstHiGHSPar ] );

  }

 return( "" );
 }

/*--------------------------------------------------------------------------*/

std::string HiGHSMILPSolver::highs_dbl_par_map( idx_type par ) const
{
 switch( par ) {
  case( dblMaxTime ): return( "time_limit" );
  case( dblRelAcc ):  return( "mip_rel_gap" );
  case( dblAbsAcc ):  return( "mip_abs_gap" );
  case( dblRAccSol ): return( "NoPar" );
  case( dblAAccSol ): return( "NoPar" );
  case( dblFAccSol ): return( "primal_feasibility_tolerance" );
  }

 if( ( par >= dblFirstHiGHSPar ) && ( par < dblLastAlgParHiGHS ) )
  return( SMSpp_to_HiGHS_dbl_pars[ par - dblFirstHiGHSPar ] );

 return( "" );
 }

 /*--------------------------------------------------------------------------*/
/*------------------- METHODS FOR HANDLING THE PARAMETERS ------------------*/
/*--------------------------------------------------------------------------*/

void HiGHSMILPSolver::set_par( idx_type par , int value )
{
 if( par == intCutSepPar ) {
  CutSepPar = value;
  return;
  }

 std::string highs_opt = highs_int_par_map( par );
 if( highs_opt.size() > 0 ) {
  // NOTE: in SMS++ we treat both int and bool HiGHS options as int. Thus, it 
  // is important to retrieve the tybe before setting them
  int type;
  Highs_getOptionType( highs , highs_opt.data() , & type );
  switch( type ) {
    case( kHighsOptionTypeBool ): 
      Highs_setBoolOptionValue( highs , highs_opt.data() , value);
      break;
    case( kHighsOptionTypeInt ):
      Highs_setIntOptionValue( highs , highs_opt.data() , value);
      break;
    default:
      throw( std::logic_error( 
      "Option type not int or bool in set_par( idx_type par , int value )" ) );
   }
  return;
  }
 //else
  //throw( std::invalid_argument( "Parameter " + int_par_idx2str(par) + 
   //     " not correctly converted in HiGHS" ) );
  

 MILPSolver::set_par( par, value );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

void HiGHSMILPSolver::set_par( idx_type par , double value )
{
 // Solver parameters explicitly mapped in HiGHS
 switch( par ) {
  case( dblUpCutOff ): UpCutOff = value; return;
  case( dblLwCutOff ): LwCutOff = value; return;
  }

 std::string highs_opt = highs_dbl_par_map( par );

 if( highs_opt.size() > 0 ) {
  if( highs_opt != "NoPar" )
    Highs_setDoubleOptionValue( highs , highs_opt.data() , value );
  
  return;
  }

 MILPSolver::set_par( par , value );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

void HiGHSMILPSolver::set_par( idx_type par , std::string && value )
{
 // HiGHS option
 if( ( par >= strFirstHiGHSPar ) && ( par < strLastAlgParHiGHS ) ) {
  std::string highs_opt = SMSpp_to_HiGHS_str_pars[ par - strFirstHiGHSPar ];
  Highs_setStringOptionValue( highs, highs_opt.data() , value.c_str() );
  return;
  }

 MILPSolver::set_par( par, std::move( value ) );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

void HiGHSMILPSolver::set_par( idx_type par , std::vector< int > && value )
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

void HiGHSMILPSolver::set_par( idx_type par ,
			     std::vector< std::string > && value )
{
 if( par == vstrConfigDBFName ) {
  auto sz = value.size();
  // delete existing Configuration with index larger than the new size
  for( Index i = sz ; i < v_ConfigDB.size() ; ++i )
   delete v_ConfigDB[ i ];
  // resize the Configuration DB: if the new size is larger than the
  // old ones, the new Configuration default to nullptr
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

Solver::idx_type HiGHSMILPSolver::get_num_int_par( void ) const {
 return( MILPSolver::get_num_int_par()
	 + intLastAlgParHiGHS - intLastAlgParMILP );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

Solver::idx_type HiGHSMILPSolver::get_num_dbl_par( void ) const {
 return( MILPSolver::get_num_dbl_par()
	 + dblLastAlgParHiGHS - dblLastAlgParMILP );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

Solver::idx_type HiGHSMILPSolver::get_num_str_par( void ) const {
 return( MILPSolver::get_num_str_par()
	 + strLastAlgParHiGHS - strLastAlgParMILP );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

Solver::idx_type HiGHSMILPSolver::get_num_vint_par( void ) const {
 return( MILPSolver::get_num_vint_par()
	 + vintLastAlgParHiGHS - vintLastAlgParMILP );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

Solver::idx_type HiGHSMILPSolver::get_num_vstr_par( void ) const {
 return( MILPSolver::get_num_vstr_par()
	 + vstrLastAlgParHiGHS - vstrLastAlgParMILP );
 }

/*--------------------------------------------------------------------------*/

int HiGHSMILPSolver::get_dflt_int_par( idx_type par ) const
{
 if( par == intCutSepPar )
  return( 0 );

 std::string highs_opt = highs_int_par_map( par );
 if( highs_opt.size() > 0 ) {
   int value, default_value;
   int type;
   Highs_getOptionType( highs , highs_opt.data() , & type );
  switch( type ) {
    case( kHighsOptionTypeBool ): 
      Highs_getBoolOptionValues( highs , highs_opt.data() ,
                                   & value, & default_value);
      break;
    case( kHighsOptionTypeInt ):
      Highs_getIntOptionValues( highs , highs_opt.data() , & value, NULL ,
                                NULL , & default_value);
      break;
    default:
      throw( std::logic_error( 
    "Option type not int or bool in get_dflt_int_par( idx_type par )" ) );
   }
   return( default_value );
  }

 return( MILPSolver::get_dflt_int_par( par ) );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

double HiGHSMILPSolver::get_dflt_dbl_par( idx_type par ) const
{
 switch( par ) {
  case( dblUpCutOff ): return( Inf< double >() );
  case( dblLwCutOff ): return( - Inf< double >() );
  }

 std::string highs_opt = highs_dbl_par_map( par );
 if( highs_opt.size() > 0 ) {
  if( highs_opt == "NoPar" )
    return( 0 );
  else{
    double value, default_value;
    Highs_getDoubleOptionValues( highs , highs_opt.data() , & value, NULL ,
                              NULL , & default_value);
    return( default_value );
   }
 }

 return( MILPSolver::get_dflt_dbl_par( par ) );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

const std::string & HiGHSMILPSolver::get_dflt_str_par( idx_type par ) const
{
 // note: this implementation is not thread safe, and it may lead to elements
 //       of value[] to be allocated more than once with some memory being
 //       lost, but the chances are too slim and the potential drawback too
 //       limited to warrant even a humble std::atomic_flag
 std::string value;
 static std::vector< std::string > default_value( strLastAlgParHiGHS -
					  strFirstHiGHSPar );

 if( ( par >= strFirstHiGHSPar ) && ( par < strLastAlgParHiGHS ) ) {
  auto i = par - strFirstHiGHSPar;
  
  if( default_value[ i ].empty() ) {
    std::string str_option = SMSpp_to_HiGHS_str_pars[ i ];

    default_value[ i ].reserve( 512 );

    // List some of the default value for HiGHS options. 
    // NOTE: this should not be necessary, but currently Highs_getStringOptionValues
    // is not properly working.
    if( str_option == "presolve" )
      default_value[ i ] = "choose";
    else if( str_option == "solver" )
      default_value[ i ] = "choose";
    else if( str_option == "parallel" )
      default_value[ i ] = "choose";
    else if( str_option == "run_crossover" )
      default_value[ i ] = "on";
    else if( str_option == "ranging" )
      default_value[ i ] = "off";
    else
      Highs_getStringOptionValues( highs ,  SMSpp_to_HiGHS_str_pars[ i ].data() ,
                                      value.data() , default_value[ i ].data() );
   }

  return( default_value[ i ] );
  }

 return( MILPSolver::get_dflt_str_par( par ) );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

const std::vector< int > & HiGHSMILPSolver::get_dflt_vint_par( idx_type par )
 const
{
 static std::vector< int > _empty;
 if( par == vintCutSepCfgInd )
  return( _empty );

 return( MILPSolver::get_dflt_vint_par( par ) );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

const std::vector< std::string > & HiGHSMILPSolver::get_dflt_vstr_par(
							idx_type par ) const
{
 static std::vector< std::string > _empty;
 if( par == vstrConfigDBFName )
  return( _empty );

 return( MILPSolver::get_dflt_vstr_par( par ) );
 }

/*--------------------------------------------------------------------------*/

int HiGHSMILPSolver::get_int_par( idx_type par ) const
{
 if( par == intCutSepPar )
  return( CutSepPar );

 std::string highs_opt = highs_int_par_map( par );
  if( highs_opt.size() > 0 ) {

   int value, default_value;
   int type;
   Highs_getOptionType( highs , highs_opt.data() , & type );
   switch( type ) {
    case( kHighsOptionTypeBool ): 
      Highs_getBoolOptionValue( highs , highs_opt.data() , & value );
      break;
    case( kHighsOptionTypeInt ):
      Highs_getIntOptionValue( highs , highs_opt.data() , & value );
      break;
    default:
      throw( std::logic_error( 
        "Option type not int or bool in get_int_par( idx_type par )" ) );
   }

   return( value );
  }

 return( MILPSolver::get_int_par( par ) );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

double HiGHSMILPSolver::get_dbl_par( idx_type par ) const
{
 switch( par ) {
  case( dblUpCutOff ): return( UpCutOff );
  case( dblLwCutOff ): return( LwCutOff );
  }

 std::string highs_opt = highs_dbl_par_map( par );
 if( highs_opt.size() > 0 ) {
  double value;
  Highs_getDoubleOptionValue( highs , highs_opt.data() , & value );
  return( value );
  }

 return( MILPSolver::get_dbl_par( par ) );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

const std::string & HiGHSMILPSolver::get_str_par( idx_type par ) const
{
 static std::string value;

 if( ( par >= strFirstHiGHSPar ) && ( par < strLastAlgParHiGHS ) ) {
  std::string highs_opt = SMSpp_to_HiGHS_str_pars[ par - strFirstHiGHSPar ];
  value.reserve( 512 );
  Highs_getStringOptionValue( highs , highs_opt.data() , value.data() );
  return( value );
  }

 return( MILPSolver::get_str_par( par ) );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

const std::vector< int > & HiGHSMILPSolver::get_vint_par( idx_type par ) const
{
 if( par == vintCutSepCfgInd )
  return( CutSepCfgInd );

 return( MILPSolver::get_vint_par( par ) );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

const std::vector< std::string > & HiGHSMILPSolver::get_vstr_par( idx_type par )
 const
{
 if( par == vstrConfigDBFName )
  return( ConfigDBFName );

 return( MILPSolver::get_vstr_par( par ) );
 }

/*--------------------------------------------------------------------------*/

Solver::idx_type HiGHSMILPSolver::int_par_str2idx(
					     const std::string & name ) const
{
 if( name == "intCutSepPar" )
  return( intCutSepPar );

 /* In HiGHSMILPSolver::*_par_str2idx() methods we check with MILPSolver first */

 idx_type idx = MILPSolver::int_par_str2idx( name );
 if( idx < Inf< idx_type >() )
  return( idx );

 // HiGHS options
 std::string highs_opt = name;
 auto array_pos = std::find( SMSpp_to_HiGHS_int_pars.begin() ,
                        SMSpp_to_HiGHS_int_pars.end() ,
                        highs_opt );

 if( array_pos != SMSpp_to_HiGHS_int_pars.end() ) {
  int pos = std::distance( SMSpp_to_HiGHS_int_pars.begin(), array_pos );
  auto idx_par = HiGHS_to_SMSpp_int_pars[ pos ].second;
  return( idx_par );
  }

 return( Inf< idx_type >() );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

const std::string & HiGHSMILPSolver::int_par_idx2str( idx_type idx ) const
{
 static const std::array< std::string , 1 > _pars =
                     { "intCutSepPar" };
 if( idx == intCutSepPar )
  return( _pars[ 0 ] );

 // note: this implementation is not thread safe, and it requires that the
 //       result is used immediately after the call (prior to any other call
 //       to int_par_idx2str()), this may have to be improved upon
 static std::string par_name;
 par_name.reserve( 512 );

 if( ( idx >= intFirstHiGHSPar ) && ( idx < intLastAlgParHiGHS ) ) {
  par_name = SMSpp_to_HiGHS_int_pars[ idx - intFirstHiGHSPar ];
  return( par_name );
  }

 return( MILPSolver::int_par_idx2str( idx ) );
 }

/*--------------------------------------------------------------------------*/

Solver::idx_type HiGHSMILPSolver::dbl_par_str2idx( const std::string & name )
 const
{
 /* In HiGHSMILPSolver::*_par_str2idx() methods we check with MILPSolver first */

 idx_type idx = MILPSolver::dbl_par_str2idx( name );
 if( idx < Inf< idx_type >() )
  return( idx );

 // HiGHS options
 std::string highs_opt = name;
 auto array_pos = std::find( SMSpp_to_HiGHS_dbl_pars.begin() ,
                        SMSpp_to_HiGHS_dbl_pars.end() ,
                        highs_opt );

 if( array_pos != SMSpp_to_HiGHS_dbl_pars.end() ) {
  int pos = std::distance( SMSpp_to_HiGHS_dbl_pars.begin(), array_pos );
  auto idx_par = HiGHS_to_SMSpp_dbl_pars[ pos ].second;
  return( idx_par );
  }

 return( Inf< idx_type >() );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

const std::string & HiGHSMILPSolver::dbl_par_idx2str( idx_type idx ) const
{
 // note: this implementation is not thread safe, and it requires that the
 //       result is used immediately after the call (prior to any other call
 //       to int_par_idx2str()), this may have to be improved upon
 static std::string par_name;
 par_name.reserve( 512 );

 if( ( idx >= dblFirstHiGHSPar ) && ( idx < dblLastAlgParHiGHS ) ) {
  par_name = SMSpp_to_HiGHS_dbl_pars[ idx - dblFirstHiGHSPar ];
  return( par_name );
  }

 return( MILPSolver::dbl_par_idx2str( idx ) );
 }

/*--------------------------------------------------------------------------*/

Solver::idx_type HiGHSMILPSolver::str_par_str2idx( const std::string & name )
 const
{
 /* In HiGHSMILPSolver::*_par_str2idx() methods we check with MILPSolver first */

 idx_type idx = MILPSolver::str_par_str2idx( name );
 if( idx < Inf< idx_type >() )
  return( idx );

 // HiGHS options
 std::string highs_opt = name;
 auto array_pos = std::find( SMSpp_to_HiGHS_str_pars.begin() ,
                        SMSpp_to_HiGHS_str_pars.end() ,
                        highs_opt);

 if( array_pos != SMSpp_to_HiGHS_str_pars.end() ) {
  int pos = std::distance( SMSpp_to_HiGHS_str_pars.begin(), array_pos );
  auto idx_par = HiGHS_to_SMSpp_str_pars[ pos ].second;
  return( idx_par );
  }

 return( Inf< idx_type >() );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

const std::string & HiGHSMILPSolver::str_par_idx2str( idx_type idx ) const
{
 // note: this implementation is not thread safe, and it requires that the
 //       result is used immediately after the call (prior to any other call
 //       to int_par_idx2str()), this may have to be improved upon
 static std::string par_name;
 par_name.reserve( 512 );

 if( ( idx >= strFirstHiGHSPar ) && ( idx < strLastAlgParHiGHS ) ) {
  par_name = SMSpp_to_HiGHS_str_pars[ idx - strFirstHiGHSPar ];
  return( par_name );
  }

 return( MILPSolver::str_par_idx2str( idx ) );
 }

/*--------------------------------------------------------------------------*/

Solver::idx_type HiGHSMILPSolver::vint_par_str2idx(
					     const std::string & name ) const
{
 if( name == "vintCutSepCfgInd" )
  return( vintCutSepCfgInd );

 return( MILPSolver::vint_par_str2idx( name ) );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

const std::string & HiGHSMILPSolver::vint_par_idx2str( idx_type idx ) const
{
 static const std::string _pars = "vintCutSepCfgInd";
 if( idx == vintCutSepCfgInd )
  return( _pars );

 return( MILPSolver::vint_par_idx2str( idx ) );
 }

/*--------------------------------------------------------------------------*/

Solver::idx_type HiGHSMILPSolver::vstr_par_str2idx(
					     const std::string & name ) const
{
 if( name == "vstrConfigDBFName" )
  return( vstrConfigDBFName );

 return( MILPSolver::vstr_par_str2idx( name ) );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

const std::string & HiGHSMILPSolver::vstr_par_idx2str( idx_type idx ) const
{
 static const std::string _pars = "vstrConfigDBFName";
 if( idx == vstrConfigDBFName )
  return( _pars );

 return( MILPSolver::vstr_par_idx2str( idx ) );
 }

/*--------------------------------------------------------------------------*/

#ifdef MILPSolver_DEBUG

void HiGHSMILPSolver::check_status( void )
{
 int nvars;
 nvars = Highs_getNumCol( highs );
 if( numcols != nvars )
  DEBUG_LOG( "numcols is " << numcols << " but Highs_getNumCol returns "
	     << nvars << std::endl );

 int nconstr;
 nconstr = Highs_getNumRow( highs );
 if( numrows != nconstr )
  DEBUG_LOG( "numrows is " << numrows << " but Highs_getNumRow returns "
	     << nconstr << std::endl );

 int nint = 0;
 for( int i = 0 ; i < numcols ; ++i ) {
  int type;
  Highs_getColIntegrality( highs , i , & type );
  if( type == kHighsVarTypeInteger )
   ++nint;
 }

 if( int_vars != nint )
  DEBUG_LOG( "int_vars is " << int_vars << " but HiGHS has actually "
	     << nint << " integer variables" << std::endl );

 MILPSolver::check_status();
 }

#endif

/*--------------------------------------------------------------------------*/
/*-------------------- PRIVATE METHODS OF THE CLASS ------------------------*/
/*----------------------------------------------------------------------------

void HiGHSMILPSolver::reload_constraint( const LinearFunction * lf )
{
 // note: this is called in response to a FunctionMod, which means that
 // while the coefficients of the LinearFunction change, their number do
 // not; hence the set of nonzeros does not change. thus, IF WE WERE SURE
 // THAT NO Variable HAVE BEEN REMOVED TO THE CONSTRAINT, this would be
 // simple because we would just have to look at the Variable that are
 // currently active. unfortunately, this is not true because Modification
 // are managed asynchronously with the model changes, and therefore it is
 // possible that some Variable that are in the Constraint for HiGHSMILPSolver
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
  if( auto idx = grb_index_of_variable( var.first ) ; idx < Inf< int >() ) {
   cols.push_back( idx );
   vals.push_back( var.second );
   }

 CPXchgcoeflist( env , lp , cols.size() , rows.data() , cols.data() ,
		 vals.data() );
 }

------------------------------------------------------------------------------

void HiGHSMILPSolver::reload_objective( Function * f )
{
 // note: this is called in response to a FunctionMod, which means that
 // while the coefficients of the LinearFunction change, their number do
 // not; hence the set of nonzeros does not change. thus, IF WE WERE SURE
 // THAT NO Variable HAVE BEEN REMOVED TO THE CONSTRAINT, this would be
 // simple because we would just have to look at the Variable that are
 // currently active. unfortunately, this is not true because Modification
 // are managed asynchronously with the model changes, and therefore it is
 // possible that some Variable that are in the Objective for HiGHSMILPSolver
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

 }  // end( HiGHSMILPSolver::reload_objective )

----------------------------------------------------------------------------*/
/*
void HiGHSMILPSolver::update_problem_type( bool quad )
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
 }*/

/*--------------------------------------------------------------------------*/

Configuration * HiGHSMILPSolver::get_cfg( Index ci ) const
{
 if( ci >= CutSepCfgInd.size() )
  return( nullptr );
 auto dbi = CutSepCfgInd[ ci ];
 if( ( dbi < 0 ) || ( Index( dbi ) >= v_ConfigDB.size() ) )
  return( nullptr );
 return( v_ConfigDB[ dbi ] );
 }

/*--------------------------------------------------------------------------*/

void HiGHSMILPSolver::generate_qobj_hessian( void )
{
 q_obj_begin.resize( numcols , 0 );

 // In MILPSolver we are storing the strictly lower triangular part of the 
 // matrix. IN HiGHS we need to pass the upper one, so we have to "transpose" the
 // data.
 int count_nnz_col = 0;
 for( int j = 0 ; j < numcols ; j++ ) {
  // Update qmatbeg with the results found in the previous iteration
  if( j > 0)
    q_obj_begin[ j ] = q_obj_begin[ j - 1 ] + count_nnz_col;

  // Reinitialize counter of column nonzeros
  count_nnz_col = 0;

  // Firstly, search if the diagonal term is non-zero
  if( q_objective[ j ] != 0 ) {
    // Update count for variable j
    count_nnz_col++;

    q_obj_val.push_back( 2*q_objective[ j ] );
    q_obj_ind.push_back( j );
   }

  // Search for the first quadratic terms x_j*x_h 
  auto it_row = find( ndq_colind.begin() , ndq_colind.end() , j );

  while( it_row != ndq_colind.end() ) {
    int pos = it_row - ndq_colind.begin();
    int var2_ind = ndq_rowind[ pos ]; // collect h (h should be always lower than j)
    double q_coeff = ndq_objective[ pos ];
    
    // Update count for variable j
    count_nnz_col++;
    
    q_obj_val.push_back( q_coeff );
    q_obj_ind.push_back( var2_ind );

    // Search for the next occurrence of j
    it_row = find( it_row + 1 , ndq_colind.end() , j );
  }

  } // end( main loop )*/
 } // end( HiGHSMILPSolver::generate_qobj_hessian )

/*--------------------------------------------------------------------------*/
/*--------------------- End File HiGHSMILPSolver.cpp -------------------------*/
/*--------------------------------------------------------------------------*/
