/*--------------------------------------------------------------------------*/
/*------------------------- File GRBMILPSolver.cpp -------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Implementation of the GRBMILPSolver class.
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

#include <DQuadFunction.h>

#include "GRBMILPSolver.h"

#ifdef MILPSOLVER_DEBUG
 #define DEBUG_LOG( stuff ) std::cout << "[MILPSolver DEBUG] " << stuff
#else
 #define DEBUG_LOG( stuff )
#endif

// include the proper GUROBI parameter mapping
#include <boost/preprocessor/cat.hpp>
#include <boost/preprocessor/stringize.hpp>
#include BOOST_PP_STRINGIZE( BOOST_PP_CAT( BOOST_PP_CAT( BOOST_PP_CAT( BOOST_PP_CAT( GRB , GRB_VERSION_MAJOR ) , GRB_VERSION_MINOR ) , GRB_VERSION_TECHNICAL ) , _maps.h ) )

/*--------------------------------------------------------------------------*/
/*------------------------- NAMESPACE AND USING ----------------------------*/
/*--------------------------------------------------------------------------*/

using namespace SMSpp_di_unipi_it;

/*--------------------------------------------------------------------------*/
/*-------------------------- FACTORY MANAGEMENT ----------------------------*/
/*--------------------------------------------------------------------------*/

SMSpp_insert_in_factory_cpp_0( GRBMILPSolver );

/*--------------------------------------------------------------------------*/
/*----------------------------- FUNCTIONS ----------------------------------*/
/*--------------------------------------------------------------------------*/

int GRBMILPSolver_callback( GRBmodel * model , void * cbdata , int where ,
			    void * usrdata )
{
 // just defer to the class method
 return( static_cast< GRBMILPSolver * >( usrdata
					 )->callback( model , cbdata , where )
	 );
 }

/*--------------------------------------------------------------------------*/
/*--------------------- CONSTRUCTOR AND DESTRUCTOR -------------------------*/
/*--------------------------------------------------------------------------*/

GRBMILPSolver::GRBMILPSolver( void ) :
 MILPSolver() , env( nullptr ) , model( nullptr ) , f_callback_set( false ) ,
 last_static_rng_con( -1 ) , throw_reduced_cost_exception( 0 ) , CutSepPar( 0 ) ,
 UpCutOff( Inf< double >() ) , LwCutOff( - Inf< double >() )
{
 int status = 0;
 status = GRBemptyenv( &env );
 if( status != 0 )
  throw( std::runtime_error( "GRBemptyenv returned with status " +
			     std::to_string( status ) ) );

 GRBsetintparam( env , GRB_INT_PAR_LOGTOCONSOLE , 0 );
 // suppress Gurobi logging
 
 status = GRBstartenv( env );
 if( status != 0 )
  throw( std::runtime_error( "GRBstartenv returned with status " +
			     std::to_string( status ) ) );
 }

/*--------------------------------------------------------------------------*/

GRBMILPSolver::~GRBMILPSolver()
{
 for( auto el : v_ConfigDB )
  delete el;

 if( model )
  GRBfreemodel( model );

 GRBfreeenv( env );
 }

 /*--------------------------------------------------------------------------*/
/*--------------------- DERIVED METHODS OF BASE CLASS ----------------------*/
/*--------------------------------------------------------------------------*/

void GRBMILPSolver::set_Block( Block * block )
{
 if( block == f_Block )
  return;

 MILPSolver::set_Block( block );
 UpCutOff = Inf< double >();
 LwCutOff = - Inf< double >();
 }

/*--------------------------------------------------------------------------*/

void GRBMILPSolver::clear_problem( unsigned int what )
{
 MILPSolver::clear_problem( 0 );

 if( model ) {
  GRBfreemodel( model );
  model = nullptr;
  }
 }

 /*--------------------------------------------------------------------------*/

void GRBMILPSolver::load_problem( void )
{
 MILPSolver::load_problem();

 int status = 0;
 if( model ) {
  GRBfreemodel( model );
  model = nullptr;
 }

 std::vector< double > grb_lb = lb;
 std::vector< double > grb_ub = ub;
 std::vector< double > grb_rhs = rhs;

 for( int i = 0 ; i < numcols ; ++i ) {
  if( grb_lb[ i ] == -Inf< double >() )
   grb_lb[ i ] = -GRB_INFINITY;
  if( grb_ub[ i ] == Inf< double >() )
   grb_ub[ i ] = GRB_INFINITY;
  }

 for( int i = 0 ; i < numrows ; ++i ) {
  if( grb_rhs[ i ] == -Inf< double >() )
   grb_rhs[ i ] = -GRB_INFINITY;
  else
   if( grb_rhs[ i ] == Inf< double >() )
    grb_rhs[ i ] = GRB_INFINITY;
  }

 // creating model with only variables inside
 if( use_custom_names )
  status = GRBnewmodel( env , & model , prob_name.c_str() ,
  						numcols ,  objective.data() , grb_lb.data() , 
						grb_ub.data() , xctype.data() , colname.data() );
 else
  status = GRBnewmodel( env , & model , prob_name.c_str() ,
  						numcols ,  objective.data() , grb_lb.data() , 
						grb_ub.data() , xctype.data() , NULL );

 // setting model sense
 GRBsetintattr( model , GRB_INT_ATTR_MODELSENSE , objsense);

 bool is_qp = std::any_of( q_objective.begin() ,
                           q_objective.end() ,
                           []( double d ) { return( d != 0 ); } );
 if( is_qp ) {
  std::vector< double > double_q_obj = q_objective;
  //for( auto & qi : double_q_obj )
   //qi *= 0.5;
  
  std::vector< int > qp_indices;
  int n_qp = 0;

  // creating a vector containg only non-zero coefficients for quadratic terms and corresponding indices
  for( int i = 0 ; i < numcols ; ++i ) {
	 if( double_q_obj[n_qp] != 0 ) {
	  n_qp = n_qp + 1;
	  qp_indices.push_back( i );
	 }
	 else{
    auto iter = double_q_obj.begin();
	  double_q_obj.erase( std::next( iter , n_qp) );
    }
  }

  // adding q_objective information automatically changes the problem type
  // from linear to quadratic
  GRBaddqpterms( model , n_qp , qp_indices.data() , qp_indices.data() , double_q_obj.data() );
  }

 // transposing the coefficient matrix

 std::vector< int > n_nz_row( numrows, 0 );
 // retrieving number of nonzeros in each row
 for( size_t i = 0 ; i < matind.size() ; ++i ) {
  int row = matind[ i ];
  ++n_nz_row[ row ];
 }

 // constructing the transposed matrix and transforming sense in grb_sense
 std::vector< double > matval_t( matval.size() , 0.0 );
 std::vector< int > matind_t( matind.size() , 0);
 std::vector< int > matbeg_t( numrows , 0 );
 std::vector< char > grb_sense( numrows , 'R' );
 // filling matbeg_t
 for( int j = 0 ; j < numrows ; ++j ) {
  switch( sense[ j ] ) {
    case( 'L' ): grb_sense[ j ] = GRB_LESS_EQUAL;
            break;
    case( 'E' ): grb_sense[ j ] = GRB_EQUAL;
            break;
    case( 'G' ): grb_sense[ j ] = GRB_GREATER_EQUAL;
            break;
    }
  
  if( j > 0 )
    matbeg_t[ j ] = matbeg_t[ j - 1 ] + n_nz_row[ j - 1 ];
 }
 
 int z = 0;
 std::vector< int > inserted_el_row( numrows , 0 );
 // filling matind_t and matval_t
 for( int i = 0 ; i < matind.size() ; ++i ) {
  int row = matind[ i ];
  while( z != numcols - 1 && i == matbeg[ z + 1 ] ) // we stepped to the next column
    ++z;
  int pos = matbeg_t[ row ] + inserted_el_row[ row ]; // where we have to insert the new value
  ++inserted_el_row[ row ];
  matval_t[ pos ] = matval[ i ];
  matind_t[ pos ] = z;
 }
  
 // adding constraints (grouping non ranged and singularly ranged)
 int n_ranged_con = 0;
 for( int j = 0 ; j < numrows ; ) {

  int tmp = j;
  int tot_nnz = 0;
  int n_constrs = 0; // number of non ranged constraints in group
  while( grb_sense[ j ] != 'R' ) {
    ++n_constrs;
    tot_nnz = tot_nnz + n_nz_row[ j ];
    ++j;
    if( j == numrows )
      break;
  }

  if( sense[ tmp ] != 'R' ) { // not ranged case

    std::vector< int > matbeg_group_con( n_constrs , 0 );
  
    // filling matbeg for the group of constraint
    for( int i = 0 ; i < n_constrs ; ++i )
      matbeg_group_con[ i ] = matbeg_t[ tmp + i ] - matbeg_t[ tmp ];

    std::vector< int > matind_group_con( tot_nnz , 0 );
    std::vector< double > matval_group_con( tot_nnz , 0 );

    // filling matind and matval for the group of constraint
    for( int i = 0 ; i < tot_nnz ; ++i ) {
      matind_group_con[ i ] = matind_t[ matbeg_t[ tmp ] + i ];
      matval_group_con[ i ] = matval_t[ matbeg_t[ tmp ] + i ];
    }
    
    if( use_custom_names )
      GRBaddconstrs( model , n_constrs , tot_nnz , matbeg_group_con.data() , 
                    matind_group_con.data() , matval_group_con.data() , & grb_sense[ tmp ] ,
                    & grb_rhs[ tmp ] , & rowname[ tmp ] );
    else
      GRBaddconstrs( model , n_constrs , tot_nnz , & matbeg_t[ tmp ] , 
                    & matind_t[ tmp ] , & matval_t[ tmp ] , & grb_sense[ tmp ] ,
                    & grb_rhs[ tmp ] , NULL );
  }
  else { // ranged case
    char * name = use_custom_names ? rowname[ j ] : NULL; // retrieve constraint name
    
    // Filling map between ranged constraint and auxiliary variables built by Gurobi
    // See GRBMILPSolver.h for further informations
    map_rng_con_aux_var.push_back( { j , numcols + n_ranged_con } );

    if( j < static_cons)
      last_static_rng_con = n_ranged_con;

    if( rngval[j] > 0 )
      GRBaddrangeconstr( model , n_nz_row[ j ] , & matind_t[ matbeg_t[ j ] ] , 
                         & matval_t[ matbeg_t[ j ] ] , grb_rhs[ j ] , grb_rhs[ j ] + rngval[j] ,
                         name );
    else
      GRBaddrangeconstr( model , n_nz_row[ j ] , & matind_t[ matbeg_t[ j ] ] , 
                         & matval_t[ matbeg_t[ j ] ] , grb_rhs[ j ] + rngval[j] , grb_rhs[ j ] ,
                         name );
    ++n_ranged_con;
    ++j;
  }
 }

 status = GRBupdatemodel( model );
 if( status != 0 )
  throw( std::runtime_error( "GRBupdatemodel returned with status " +
			     std::to_string( status ) ) );

 // the base representation isn't needed anymore
 MILPSolver::clear_problem( 15 );

 UpCutOff = Inf< double >();
 LwCutOff = - Inf< double >();

 }  // end( GRBMILPSolver::load_problem )

/*--------------------------------------------------------------------------*/

double GRBMILPSolver::get_problem_lb( const ColVariable & var ) const
{
 double b = MILPSolver::get_problem_lb( var );
 if( b == -Inf< double >() )
  b = -GRB_INFINITY;

 return( b );
 }

/*--------------------------------------------------------------------------*/

double GRBMILPSolver::get_problem_ub( const ColVariable & var ) const
{
 double b = MILPSolver::get_problem_ub( var );
 if( b == Inf< double >() )
  b = GRB_INFINITY;

 return( b );
 }

/*--------------------------------------------------------------------------*/

std::array< double , 2 > GRBMILPSolver::get_problem_bounds(
					      const ColVariable & var ) const
{
 auto ret = MILPSolver::get_problem_bounds( var );
 if( ret[ 0 ] == -Inf< double >() )
  ret[ 0 ] = -GRB_INFINITY;
 if( ret[ 1 ] == Inf< double >() )
  ret[ 1 ] = GRB_INFINITY;

 return( ret );
 }

/*--------------------------------------------------------------------------*/

int GRBMILPSolver::compute( bool changedvars )
{
 lock();  // lock the mutex: this is done again inside MILPSolver::compute,
          // but that's OK since the mutex is recursive

 // process Modification: this is driven by MILPSolver- - - - - - - - - - - -
 if( MILPSolver::compute( changedvars ) != kOK )
  throw( std::runtime_error( "an error occurred in MILPSolver::compute()" ) );

 // if required, write the problem to file- - - - - - - - - - - - - - - - - -
 if( ! output_file.empty() ) {
  GRBwrite( model , output_file.c_str() );
 }

 // figure out which API function is to be called - - - - - - - - - - - - - -
 bool is_qp = false;
 int qp;
 GRBgetintattr( model , GRB_INT_ATTR_IS_QP , &qp );

 int qcp;
 GRBgetintattr( model , GRB_INT_ATTR_IS_QCP , &qcp );

 if( qp == 1 ) {
  DEBUG_LOG( "GUROBI problem type: QP" << std::endl );
  is_qp = true;
 }
 else if( qcp == 1 ) {
  DEBUG_LOG( "GUROBI problem type: QCP" << std::endl );
  is_qp = true;
 }
 else{
  DEBUG_LOG( "GUROBI problem type: MIP or LP" << std::endl );
 }

 // the actual call to GUROBI- - - - - - - - - - - - - - - - - - - - - - - - -

 if( int_vars > 0 ) {  // the MIP case- - - - - - - - - - - - - - - - - - - -

  if( ( CutSepPar & 7 ) ||
      ( UpCutOff < Inf< double >() ) || ( LwCutOff > Inf< double >() ) ) {
   // the callback has to be set
   GRBsetcallbackfunc( model , & GRBMILPSolver_callback , this );

   if( CutSepPar & 3 ) {  // we do user cut separation, thus we have to set the possibility in Gurobi
    GRBsetintparam( GRBgetenv( model ) , GRB_INT_PAR_PRECRUSH , 1 );
    auto md = ( CutSepPar >> 3 ) & 3;
    GRBsetintparam( GRBgetenv( model ) , GRB_INT_PAR_CUTS , md );
   }
   
   if( CutSepPar & 4 )  // we do lazy constraint separation, thus we have to set the possibility in Gurobi
    GRBsetintparam( GRBgetenv( model ) , GRB_INT_PAR_LAZYCONSTRAINTS , 1 );
   
   f_callback_set = true;
   }
  else {
   if( f_callback_set ) {    // the callback was set
    GRBsetcallbackfunc( model , NULL , nullptr );  // un-set it
    f_callback_set = false;
    }
  }

  if( int status = GRBoptimize( model ) ) { //error
   
   sol_status = decode_grb_error( status );
   goto Return_status;
   }

  int m_status;
  GRBgetintattr( model , GRB_INT_ATTR_STATUS , &m_status );

  sol_status = decode_model_status( m_status );
  goto Return_status;
  }

 // the continuous case - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( int status = GRBoptimize( model ) ) {
  sol_status = decode_grb_error( status );
  goto Return_status;
  }

 int m_status;
 GRBgetintattr( model , GRB_INT_ATTR_STATUS , &m_status );

 sol_status = decode_model_status( m_status );

 Return_status:
 unlock();  // unlock the mutex
 return( sol_status );

 }  // end( GRBMILPSolver::compute )

/*--------------------------------------------------------------------------*/


int GRBMILPSolver::decode_model_status( int status )
{
 DEBUG_LOG( "GRB_STATUS returned " << status << std::endl );

 /* The following are the symbols that may represent the status of
 * a GUROBI solution as returned by GRBgetintattr( model , GRB_INT_ATTR_STATUS , model_status ),
 * as listed in GUROBI Callable Library API manual.*/
 switch(status) {
  //case( GRB_LOADED )
  case( GRB_OPTIMAL ):
   // Model was solved to optimality (subject to tolerances), and
   // an optimal solution is available.
   return( kOK );
  case( GRB_INFEASIBLE ):
   // Problem was proven to be infeasible.
   return( kInfeasible );
  case( GRB_INF_OR_UNBD ):
   // Model was proven to be either infeasible or unbounded. To
   // obtain a more definitive conclusion, set the DualReductions
   // parameter to 0 and reoptimize
   //!!   return( kInfeasible );
   return( kUnbounded );
  case( GRB_UNBOUNDED ):
   // Model was proven to be unbounded. 
   //!! Important note: an unbounded status indicates the presence of an unbounded ray
   //!! that allows the objective to improve without limit. It says
   //!! nothing about whether the model has a feasible solution. If
   //!! you require information on feasibility, you should set the objective to zero and reoptimize
   return( kUnbounded );
  case( GRB_CUTOFF ):
   // Optimal objective for model was proven to be worse than
   // the value specified in the Cutoff parameter. No solution
   // information is available.
   return( kOK );
  case( GRB_ITERATION_LIMIT ):
   // Optimization terminated because the total number of simplex iterations performed exceeded
   // the value specified in the IterationLimit parameter, or because the total number of barrier
   // iterations exceeded the value specified in the BarIterLimit parameter
   return( kStopIter );
  case( GRB_NODE_LIMIT ):
   // Optimization terminated because the total number of branch-and-cut nodes
   // explored exceeded the value specified in the NodeLimit parameter
   return( kStopIter );
  case( GRB_TIME_LIMIT ):
   // Optimization terminated because the time expended exceeded the value 
   // specified in the TimeLimit parameter
   return( kStopTime );
  case( GRB_SOLUTION_LIMIT ):
   // Optimization terminated because the number of solutions
   // found reached the value specified in the SolutionLimit parameter
   return( kOK );
  case( GRB_INTERRUPTED ):
   //  Optimization was terminated by the user
   return( kError );
  case( GRB_NUMERIC ):
   // Optimization was terminated due to unrecoverable numerical
   // difficulties
   return( kError );
  case( GRB_SUBOPTIMAL ):
   // Unable to satisfy optimality tolerances; a sub-optimal solution is available
   return( kOK );
  case( GRB_INPROGRESS ):
   // An asynchronous optimization call was made, but the 
   // associated optimization run is not yet complete
   return( kError );
  case( GRB_USER_OBJ_LIMIT ):
   // User specified an objective limit (a bound on either the best
   // objective or the best bound), and that limit has been reached
   return( kOK );
  case( GRB_WORK_LIMIT ):
  case( GRB_MEM_LIMIT ):
   return( kError );
  default:;
 }

 throw( std::runtime_error( "GRB_STATUS returned unknown status " +
			    std::to_string( status ) ) );
 }

/*--------------------------------------------------------------------------*/

int GRBMILPSolver::decode_grb_error( int error )
{
 DEBUG_LOG( "GUROBI returned " << error << std::endl );

 /* The following symbols represent error codes returned by GUROBI, for
  * example by GRBoptimize(). Not all codes are
  * returned by all the methods, and it's hard to know what returns what
  * without doing extensive tests. So we are implementing just the ones we
  * encounter as we go. */
 switch( error ) {
  //case( GRB_ERROR_OUT_OF_MEMORY ):
  //case( GRB_ERROR_NULL_ARGUMENT ):
  //case( GRB_ERROR_INVALID_ARGUMENT ):
  //case( GRB_ERROR_UNKNOWN_ATTRIBUTE ):
  //case( GRB_ERROR_DATA_NOT_AVAILABLE ):
  //case( GRB_ERROR_INDEX_OUT_OF_RANGE ):
  //case( GRB_ERROR_UNKNOWN_PARAMETER ):
  //case( GRB_ERROR_VALUE_OUT_OF_RANGE ):
  //case( GRB_ERROR_NO_LICENSE ):
  //case( GRB_ERROR_SIZE_LIMIT_EXCEEDED ):
  //case( GRB_ERROR_CALLBACK ):
  //case( GRB_ERROR_FILE_READ ):
  //case( GRB_ERROR_FILE_WRITE ):
  //case( GRB_ERROR_NUMERIC ):
  //case( GRB_ERROR_IIS_NOT_INFEASIBLE ):
  //case( GRB_ERROR_NOT_FOR_MIP ):
  //case( GRB_ERROR_OPTIMIZATION_IN_PROGRESS ):
  //case( GRB_ERROR_DUPLICATES ):
  //case( GRB_ERROR_NODEFILE ):
  //case( GRB_ERROR_Q_NOT_PSD ):
  //case( GRB_ERROR_QCP_EQUALITY_CONSTRAINT ):
  //case( GRB_ERROR_NETWORK ):
  //case( GRB_ERROR_JOB_REJECTED ):
  //case( GRB_ERROR_NOT_SUPPORTED ):
  //case( GRB_ERROR_EXCEED_2B_NONZEROS ):
  //case( GRB_ERROR_INVALID_PIECEWISE_OBJ ):
  //case( GRB_ERROR_UPDATEMODE_CHANGE ):
  //case( GRB_ERROR_CLOUD ):
  //case( GRB_ERROR_MODEL_MODIFICATION ):
  //case( GRB_ERROR_CSWORKER ):
  //case( GRB_ERROR_TUNE_MODEL_TYPES ):
  //case( GRB_ERROR_SECURITY ):
  }

 throw( std::runtime_error( "GUROBI returned unmanaged error " +
			    std::to_string( error ) ) );
 }

/*--------------------------------------------------------------------------*/

Solver::OFValue GRBMILPSolver::get_lb( void )
{
 OFValue lower_bound = 0;
 int sense;
 GRBgetintattr( model , GRB_INT_ATTR_MODELSENSE , & sense );

 switch( sense ) {
  case( GRB_MINIMIZE ):  // Minimization problem- - - - - - - - - - - - - - - -
   switch( sol_status ) {
    case( kUnbounded ):  lower_bound = -Inf< OFValue >(); break;
    case( kInfeasible ): lower_bound = Inf< OFValue >();  break;
    case( kOK ):
    case( kStopIter ):
    case( kStopTime ):
      int m_status;
      GRBgetintattr( model , GRB_INT_ATTR_STATUS , &m_status );
      // when a gurobi model stop with cutoff status, 
      // no solution information is available
      if(m_status == GRB_CUTOFF)
        throw( std::runtime_error( "No solution information is available whit GRB_CUTOFF status" ) );

      GRBgetdblattr( model , GRB_DBL_ATTR_OBJVAL , & lower_bound );
      lower_bound += constant_value;
      break;

    default:
     // If Gurobi does not state that an optimal solution has been found
     // then we do not have enough information to provide a "good" lower bound
     // (the problem may be unbounded and Gurobi has not detected it yet).
     // Therefore, in this case, the lower bound should be -Inf.
     lower_bound = -Inf< OFValue >();
     break;
    }
   break;

  case( GRB_MAXIMIZE ):  // Maximization problem- - - - - - - - - - - - - - - -
   switch( sol_status ) {
    case( kUnbounded ):  lower_bound = Inf< OFValue >(); break;
    case( kInfeasible ): lower_bound = -Inf< OFValue >(); break;

    // if the algorithm has been stopped, the bound only exists if a
    // feasible solution has been generated
    case( kStopIter ):
    case( kStopTime ):
     if( ! has_var_solution() ) {
      lower_bound = - Inf< OFValue >();
      break;
      }
     lower_bound += constant_value;
     break;

    case( kOK ):
     int m_status;
     GRBgetintattr( model , GRB_INT_ATTR_STATUS , &m_status );
     // when a gurobi model stop with cutoff status, 
     // no solution information is available
     if(m_status == GRB_CUTOFF)
       throw( std::runtime_error( "No solution information is available whit GRB_CUTOFF status" ) );
     GRBgetdblattr( model , GRB_DBL_ATTR_OBJVAL , & lower_bound );
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

Solver::OFValue GRBMILPSolver::get_ub( void )
{
 OFValue upper_bound = 0;
 int sense;
 GRBgetintattr( model , GRB_INT_ATTR_MODELSENSE , & sense );

 switch( sense ) {
  case( GRB_MINIMIZE ):  // Minimization problem- - - - - - - - - - - - - - - -
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
     break;

    case( kOK ):
     int m_status;
     GRBgetintattr( model , GRB_INT_ATTR_STATUS , &m_status );
     // when a gurobi model stop with cutoff status, 
     // no solution information is available
     if(m_status == GRB_CUTOFF)
       throw( std::runtime_error( "No solution information is available whit GRB_CUTOFF status" ) );
     GRBgetdblattr( model , GRB_DBL_ATTR_OBJVAL , & upper_bound );
     upper_bound += constant_value;
     break;

    default:
     // If Gurobi does not state that an optimal solution has been found
     // then we do not have enough information to provide a "good" upper bound
     // (the problem may be unbounded and Gurobi has not detected it yet).
     // Therefore, in this case, the upper bound should be +Inf.
     upper_bound = Inf< OFValue >();
     break;
    }
   break;

  case( GRB_MAXIMIZE ):  // Maximization problem- - - - - - - - - - - - - - - -
   switch( sol_status ) {
    case( kUnbounded ):  upper_bound = Inf< OFValue >(); break;
    case( kInfeasible ): upper_bound = -Inf< OFValue >(); break;

    case( kOK ):
    case( kStopIter ):
    case( kStopTime ):
     int m_status;
     GRBgetintattr( model , GRB_INT_ATTR_STATUS , &m_status );
     // when a gurobi model stop with cutoff status, 
     // no solution information is available
     if(m_status == GRB_CUTOFF)
       throw( std::runtime_error( "No solution information is available whit GRB_CUTOFF status" ) );
     GRBgetdblattr( model , GRB_DBL_ATTR_OBJVAL , & upper_bound );
     upper_bound += constant_value;
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

bool GRBMILPSolver::has_var_solution( void )
{
 int n_solution;
 int status = GRBgetintattr( model , GRB_INT_ATTR_SOLCOUNT , & n_solution );

 if( status )
  throw( std::runtime_error( "An error occurred in getting GRB_SOLCOUNT" ) );

 if( n_solution > 0 ) // there is at least a solution
   return( true );
 else // no solution found
   return( false );
 }

/*--------------------------------------------------------------------------*/

Solver::OFValue GRBMILPSolver::get_var_value( void )
{
 int objsense;
 GRBgetintattr( model , GRB_INT_ATTR_MODELSENSE , & objsense );

 switch( objsense ) {
  case( GRB_MINIMIZE ): return( get_ub() );
  case( GRB_MAXIMIZE ): return( get_lb() );
  default: throw( std::runtime_error( "Objective type not yet defined" ) );
  }
 }

/*--------------------------------------------------------------------------*/

void GRBMILPSolver::get_var_solution( Configuration * solc )
{
 int n_ranged_con = map_rng_con_aux_var.size();
 std::vector< double > x( numcols , 0 );
 std::vector< double > x_grb( numcols + n_ranged_con, 0 );

 if( GRBgetdblattrarray( model , GRB_DBL_ATTR_X , 0 , numcols + n_ranged_con , x_grb.data() ) )
  throw( std::runtime_error( "Unable to get the solution with GRB_DBL_ATTR_X" ) );

 if( n_ranged_con == 0 ) // there are no ranged constraint. Thus, no aux var in Gurobi
  x = x_grb;
 else{
  int aux_counter = 0;
  for( int j = 0 ; j < numcols + n_ranged_con ; ++j ) {
    if( j != map_rng_con_aux_var[aux_counter].second ) // column j is not an auxiliary variable
      x[ j - aux_counter ] = x_grb[ j ];
    else
      ++aux_counter;
  }
 }

 get_var_solution( x );
 }

/*--------------------------------------------------------------------------*/

void GRBMILPSolver::get_var_solution( const std::vector< double > & x )
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

bool GRBMILPSolver::has_dual_solution( void )
{
 int isMIP = 1;
 if( GRBgetintattr( model , GRB_INT_ATTR_IS_MIP , &isMIP ) )
  throw( std::runtime_error( "An error occurred in getting GRB_INT_ATTR_IS_MIP" ) );
 
 int verbosity = 0;
 GRBgetintparam( env , GRB_INT_PAR_LOGTOCONSOLE , &verbosity );

 if( ( isMIP ) ){
 // The model is a MIP
    if( verbosity )
      DEBUG_LOG( "Dual solution for MIP model not available" << std::endl);
    return( false );  
  }

 int infunbd_info = 0;
 if( GRBgetintparam( env , GRB_INT_PAR_INFUNBDINFO , &infunbd_info ) )
  throw( std::runtime_error( "An error occurred in getting GRB_INT_PAR_INFUNBDINFO" ) );

 if( !infunbd_info ){
  if( verbosity )
    DEBUG_LOG( "In order to ask for the dual solution of"
                      "the model, the parameter GRB_INT_PAR_INFUNBDINFO" 
                      "should be set to 1" << std::endl);
  return( false );
 }

 return( true );
}

/*--------------------------------------------------------------------------*/

bool GRBMILPSolver::is_dual_feasible( void )
{
 return( false );  // TODO
 }

/*--------------------------------------------------------------------------*/

void GRBMILPSolver::get_dual_solution( Configuration * solc )
{
 int n_ranged_con = map_rng_con_aux_var.size();
 std::vector< double > pi( numrows , 0 );
 std::vector< double > dj( numcols , 0 );
 std::vector< double > dj_grb( numcols + n_ranged_con , 0 );

 if( numrows > 0 )
  if( GRBgetdblattrarray( model , GRB_DBL_ATTR_PI , 0 , numrows , pi.data() ) )
   throw( std::runtime_error( "Unable to get dual values querying the attribute GRB_PI" ) );

 if( GRBgetdblattrarray( model , GRB_DBL_ATTR_RC , 0 , numcols + n_ranged_con , dj_grb.data() ) )
  throw( std::runtime_error( "Unable to get reduced costs querying the attribute GRB_RC") );

 if( n_ranged_con == 0 ) // there are no ranged constraint. Thus, no aux var in Gurobi
  dj = dj_grb;
 else{
  int aux_counter = 0;
  for( int j = 0 ; j < numcols + n_ranged_con ; ++j ) {
    if( j != map_rng_con_aux_var[aux_counter].second ) // column j is not an auxiliary variable
      dj[ j - aux_counter ] = dj_grb[ j ];
    else
      ++aux_counter;
  }
 }

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
	       "GRBMILPSolver::get_dual_solution: invalid dual value." ) );

  if( throw_reduced_cost_exception ) {
   if( var_is_fixed && ( ! lhs_con ) && ( var_lb != 0 ) ) {
    /* The Variable is fixed but it has no associated OneVarConstraint
     * with both bounds equal to the value of the Variable. */

    throw( std::logic_error(
     "GRBMILPSolver::get_dual_solution: variable with index " +
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
                "GRBMILPSolver::get_dual_solution: variable with index " +
		std::to_string( i ) + " has no OneVarConstraint." ) );
     }
   }
  }
 }  // end( GRBMILPSolver::get_dual_solution )

/*--------------------------------------------------------------------------*/

bool GRBMILPSolver::has_dual_direction( void )
{
 double proof;
 std::vector< double > y( numrows , 0 );

 int model_status;
 GRBgetintattr( model , GRB_INT_ATTR_STATUS , & model_status );

 switch( model_status ) {
  case( GRB_INFEASIBLE ):
  case( GRB_INF_OR_UNBD ):
  case( GRB_UNBOUNDED ):  break;
  default: DEBUG_LOG( "Status of the Gurobi model not infeasible or "
      "unbounded" << std::endl );
    return( false );                       
 }

 int infunbd_info;
 GRBgetintparam( env , GRB_INT_PAR_INFUNBDINFO , & infunbd_info );
 if( !infunbd_info ){
  DEBUG_LOG( "In order to ask for the farkas proof of the model, the "
    "parameter GRB_INT_PAR_INFUNBDINFO should be set to 1" << std::endl );
  return( false );
 }
 
 return( true );
}

/*--------------------------------------------------------------------------*/

void GRBMILPSolver::get_dual_direction( Configuration * dirc )
{
 int n_ranged_con = map_rng_con_aux_var.size();
 std::vector< double > y( numrows , 0 );
 std::vector< double > dj( numcols , 0 );
 std::vector< double > dj_grb( numcols + n_ranged_con , 0 );

 double proof;

 // FARKASPROOF and FARKASDUAL gives a Farkas certificate y so that:
 // y' * A * x >= y' * b
 //   If it is a <= constraint then y[i] <= 0 holds;
 //   If it is a >= constraint then y[i] >= 0 holds.

 int status_proof = GRBgetdblattr( model , GRB_DBL_ATTR_FARKASPROOF , & proof );
 int status_y = GRBgetdblattrarray( model , GRB_DBL_ATTR_FARKASDUAL , 0 , numrows , y.data() );

 // reverse the sign of y due to Gurobi approach
 for( auto i = y.begin() ; i != y.end() ; ++i  )
  *i = -*i;

 if( status_proof != 0 || status_y != 0 )
  throw( std::runtime_error( "an error occurred in getting Farkas certificate" ) );

 if( GRBgetdblattrarray( model , GRB_DBL_ATTR_RC , 0 , numcols + n_ranged_con , dj_grb.data() ) )
  throw( std::runtime_error( "Unable to get reduced costs querying the attribute GBL_RC") );

 if( n_ranged_con == 0 ) // there are no ranged constraint. Thus, no aux var in Gurobi
  dj = dj_grb;
 else{
  int aux_counter = 0;
  for( int j = 0 ; j < numcols + n_ranged_con ; ++j ) {
    if( j != map_rng_con_aux_var[aux_counter].second ) // column j is not an auxiliary variable
      dj[ j - aux_counter ] = dj_grb[ j ];
    else
      ++aux_counter;
  }
 }

 int row = 0;
 int row_dynamic = static_cons;
 auto set = [ y , & row ]( FRowConstraint & c ) {
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
	       "GRBMILPSolver::get_dual_direction: invalid dual value" ) );

  if( throw_reduced_cost_exception ) {
   if( var_is_fixed && ( ! lhs_con ) && ( var_lb != 0 ) ) {
    /* The Variable is fixed but it has no associated OneVarConstraint
     * with both bounds equal to the value of the Variable. */

    throw( std::logic_error(
     "GRBMILPSolver::get_dual_direction: variable with index " +
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
                "GRBMILPSolver::get_dual_direction: variable with index " +
		std::to_string( i ) + " has no OneVarConstraint." ) );
     }
   }
  }
 }  // end( GRBMILPSolver::get_dual_direction )

/*--------------------------------------------------------------------------*/

void GRBMILPSolver::write_lp( const std::string & filename )
{
 std::string output_file_lp;
 std::stringstream X(output_file);
 std::getline( X , output_file_lp , '.');
 output_file_lp = output_file_lp.append(".lp");
 GRBwrite( model , output_file_lp.c_str() );
 }

/*--------------------------------------------------------------------------*/

int GRBMILPSolver::get_nodes( void ) const
{
 double nodecnt;
 GRBgetdblattr( model , GRB_DBL_ATTR_NODECOUNT , & nodecnt );
 return( nodecnt );
 }

/*--------------------------------------------------------------------------*/

int GRBMILPSolver::grb_index_of_variable( const ColVariable * var ) const
{
 auto idx = index_of_variable( var );
 if( idx == Inf< int >() )
  return( idx );

 int n_ranged_con = map_rng_con_aux_var.size();
 
 if( n_ranged_con != 0 ) {
  int tmp_count = 0;
    while( idx >= map_rng_con_aux_var[ tmp_count ].second  && tmp_count < n_ranged_con ) {
      ++tmp_count;
      ++idx;
    }
  }

  return( idx );
 }

 /*--------------------------------------------------------------------------*/

int GRBMILPSolver::grb_index_of_dynamic_variable( const ColVariable * var ) const
{
 auto idx = index_of_dynamic_variable( var );
 if( idx == Inf< int >() )
  return( idx );

 int n_ranged_con = map_rng_con_aux_var.size();
 
 if( n_ranged_con != 0 ) {
  int tmp_count = last_static_rng_con + 1;
    while( idx >= map_rng_con_aux_var[ tmp_count ].second  && tmp_count < n_ranged_con ) {
      ++tmp_count;
      ++idx;
    }
  }

  return( idx );
 }

 /*--------------------------------------------------------------------------*/
/*-------------------- PROTECTED FIELDS OF THE CLASS -----------------------*/
/*--------------------------------------------------------------------------*/

void GRBMILPSolver::var_modification( const VariableMod * mod )
{
 bool was_mip = int_vars > 0;  // was a MIP before the change

 // call the method of MILPSolver to update dictionaries (and int_vars)
 MILPSolver::var_modification( mod );

 bool is_mip = int_vars > 0;  // is a MIP after the change

 auto var = static_cast< const ColVariable * >( mod->variable() );
 auto idx = grb_index_of_variable( var );
 if( idx == Inf< int >() )  // the Variable is not (yet) there (?)
  return;                   // nothing to do
        
 // react to changes in the integrality - - - - - - - - - - - - - - - - - - -
 if( ColVariable::is_integer( mod->old_state() ) !=
     ColVariable::is_integer( mod->new_state() ) ) {

  // construct new variable type
  char new_ctype;
  if( var->is_integer() && ( ! relax_int_vars ) ) {
   if( var->is_unitary() && var->is_positive() )
    new_ctype = GRB_BINARY; // Binary
   else
    new_ctype = GRB_INTEGER; // Integer
   }
  else
   new_ctype = GRB_CONTINUOUS;  // Continuous

  GRBsetcharattrelement( model , GRB_CHAR_ATTR_VTYPE , idx , new_ctype );
  }   // end( if( new integrality != old integrality ) )

 // react to fix / unfix- - - - - - - - - - - - - - - - - - - - - - - - - - -
 if( Variable::is_fixed( mod->old_state() ) !=
     Variable::is_fixed( mod->new_state() ) ) {
  if( Variable::is_fixed( mod->new_state() ) ) {  // fix the variable
    std::array< double , 1 > bd = { var->get_value() };
    GRBsetdblattrelement( model , GRB_DBL_ATTR_LB , idx , bd[ 0 ] );
    GRBsetdblattrelement( model , GRB_DBL_ATTR_UB , idx , bd[ 0 ] );
    }
   else {                                         // un fix the variable
    auto bd = GRBMILPSolver::get_problem_bounds( *var );
    GRBsetdblattrelement( model , GRB_DBL_ATTR_LB , idx , bd[ 0 ] );
    GRBsetdblattrelement( model , GRB_DBL_ATTR_UB , idx , bd[ 1 ] );
    }
  }
 }  // end( GRBMILPSolver::var_modification )

/*--------------------------------------------------------------------------*/

void GRBMILPSolver::objective_modification( const ObjectiveMod * mod )
{
 // call the method of MILPSolver to update sense direction
 MILPSolver::objective_modification( mod );

 /* ObjectiveMod class does not include any modification types except
  * for eSetMin and eSetMax.
  * To change OF coefficients, a FunctionMod must be used. */

 switch( mod->type() ) {
  case( ObjectiveMod::eSetMin ):
   GRBsetintattr( model , GRB_INT_ATTR_MODELSENSE , GRB_MINIMIZE );
   break;
  case( ObjectiveMod::eSetMax ):
   GRBsetintattr( model , GRB_INT_ATTR_MODELSENSE , GRB_MAXIMIZE );
   break;
  default: throw( std::invalid_argument( "Invalid type of ObjectiveMod" ) );
  }
 }

/*--------------------------------------------------------------------------*/

void GRBMILPSolver::const_modification( const ConstraintMod * mod )
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
 char sense;
 double rhs;
 double rngval;

 RowConstraint::RHSValue con_lhs = NAN;
 RowConstraint::RHSValue con_rhs = NAN;

 // find if con is a ranged constraint
 auto it_rng = std::find_if( map_rng_con_aux_var.begin(), map_rng_con_aux_var.end(), 
    [&index]( std::pair< int , int > const& elem ) {
    return( elem.first == index );
  });
 bool is_rng = ( it_rng != map_rng_con_aux_var.end() ); // 0 isn't a ranged constraint

 switch( mod->type() ) {
  case( ConstraintMod::eRelaxConst ):
    // In order to relax the constraint all we do is transform it
    // into an inequality (<=) with RHS equal to infinity.
    // NOTE: for ranged constraint we can do the same BUT in the reverse 
    // process we have to remember to set the sense to ==

    sense = GRB_LESS_EQUAL;
    rhs = GRB_INFINITY;
    GRBsetcharattrelement( model , GRB_CHAR_ATTR_SENSE , index , sense );
    GRBsetdblattrelement( model , GRB_DBL_ATTR_RHS , index , rhs );
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
    sense = GRB_EQUAL;
    rhs = con_rhs;
    }
   else
    if( con_lhs == -Inf< double >() ) {
     sense = GRB_LESS_EQUAL;
     rhs = con_rhs;
     }
    else
     if( con_rhs == Inf< double >() ) {
      sense = GRB_GREATER_EQUAL;
      rhs = con_lhs;
      }
     else {
      sense = 'R';
      rhs = con_rhs;
      rngval = con_rhs - con_lhs;
      }
   
   if( sense != 'R' ) {
    GRBsetcharattrelement( model , GRB_CHAR_ATTR_SENSE , index , sense );
    GRBsetdblattrelement( model , GRB_DBL_ATTR_RHS , index , rhs );
    
    if( is_rng ) {
     // we are modifying a previously ranged constraint into a non ranged one.
     // We allow this to happen, but we have to set the bound on the auxiliary
     // variable to 0
     int idx_aux_var = ( *it_rng ).second;
     GRBsetdblattrelement( model , GRB_DBL_ATTR_UB , idx_aux_var , rngval );
    }
   }
   else{
    // GRBMILPSolver doesn't support change in linear constraint sense from
    // non ranged to ranged
    if( ! is_rng )
      throw( std::invalid_argument( "Tried to convert a non ranged constraint "
                                    "into a ranged one. GRBMILPSolver does not "
                                    "support this function." ) );
    
    int idx_aux_var = ( *it_rng ).second;
    GRBsetcharattrelement( model , GRB_CHAR_ATTR_SENSE , index , GRB_EQUAL );
    GRBsetdblattrelement( model , GRB_DBL_ATTR_UB , idx_aux_var , rngval );
    GRBsetdblattrelement( model , GRB_DBL_ATTR_RHS , index , rhs );
   }
   
   break;

  default:
   throw( std::invalid_argument( "Invalid type of ConstraintMod" ) );
  }
 }  // end( GRBMILPSolver::const_modification )

/*--------------------------------------------------------------------------*/

void GRBMILPSolver::bound_modification( const OneVarConstraintMod * mod )
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

 // fixed variables are implemented in GRBMILPSolver by changing the bounds;
 // therefore, actual changes of the bounds are ignored here. note that we
 // are assuming the new bounds do not make the fixed value of the variable
 // unfeasible, as this will make the whole problem unfeasible
 // TODO: check this
 if( var->is_fixed() )
  return;

 auto vi = grb_index_of_variable( var );
 if( vi == Inf< int >() )  // the ColVariable has been removed
  return;                  // is strange, but there is nothing to do

 switch( mod->type() ) {

  case( RowConstraintMod::eChgLHS ): {
   std::array< double , 1 > bd = { GRBMILPSolver::get_problem_lb( *var ) };
   GRBsetdblattrelement( model , GRB_DBL_ATTR_LB , vi , bd[0] );
   break;
   }

  case( RowConstraintMod::eChgRHS ): {
   std::array< double , 1 > bd = { GRBMILPSolver::get_problem_ub( *var ) };
   GRBsetdblattrelement( model , GRB_DBL_ATTR_UB , vi , bd[ 0 ] );
   break;
   }

  case( RowConstraintMod::eChgBTS ): {
   auto bd = GRBMILPSolver::get_problem_bounds( *var );
   GRBsetdblattrelement( model , GRB_DBL_ATTR_LB , vi , bd[0] );
   GRBsetdblattrelement( model , GRB_DBL_ATTR_UB , vi , bd[ 1 ] );
   break;
   }

  default:
   throw( std::invalid_argument( "Invalid type of OneVarConstraintMod" ) );
  }
 }  // end( GRBMILPSolver::bound_modification )

/*--------------------------------------------------------------------------*/

void GRBMILPSolver::objective_function_modification( const FunctionMod * mod )
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
     auto vi = grb_index_of_variable( static_cast< const ColVariable * >( v ) );

     *(cidxit++) = vi ;
    }

   auto nsz = std::distance( nval.begin() , nvit );
   cidx.resize( nsz );
   nval.resize( nsz );

   for( size_t i = 0 ; i < cidx.size() ; ++i )
    GRBsetdblattrelement( model , GRB_DBL_ATTR_OBJ , cidx[ i ] , nval[ i ]  );

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
     auto vi = grb_index_of_variable( static_cast< const ColVariable * >( v ) );

     *(cidxit++) = vi;
    }

   auto nsz = std::distance( nval.begin() , nvit );
   cidx.resize( nsz );
   nval.resize( nsz );

   for( size_t i = 0 ; i < cidx.size() ; ++i )
    GRBsetdblattrelement( model , GRB_DBL_ATTR_OBJ , cidx[ i ] , nval[ i ]  );

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

  // In Gurobi to change quadratic coefficients we need to retrieve all the old coeff.,
  // and then add the difference between the new and the old ones

  int nqz;
  GRBgetintattr( model , GRB_INT_ATTR_NUMQNZS  , & nqz );

  std::vector< int > oldind_row ( nqz );
  std::vector< int > oldind_col ( nqz );
  std::vector< double > oldval ( nqz );

  int status = GRBgetq( model, & nqz , oldind_row.data() , oldind_col.data() , oldval.data() );
  if( status != 0 )
   throw( std::runtime_error( "Error while querying quadratic coefficients with GRBgetq" ) );

  for( auto v : *vars )
   if( auto idx = *(idxit++) ; idx < Inf< Index >() ) {
    *(nvit++) = std::get< 1 >( cp[ idx ] );
    auto cidx = grb_index_of_variable( static_cast< const ColVariable * >( v ) );
     
    *(cidxit++) = cidx ;
    
    auto arr_idx_row = std::find(oldind_row.begin(), oldind_row.end(), cidx);
    auto arr_idx_col = std::find(oldind_col.begin(), oldind_col.end(), cidx);
    double old_var_q_coeff;

    if( arr_idx_row == oldind_row.end() ) // no quadratic coefficient was already set for the variable
      old_var_q_coeff = 0.0;
    else{
      if( *arr_idx_row != *arr_idx_col )
        throw( std::runtime_error( "Error while modifying quadratic coefficients" ) );

      old_var_q_coeff = oldval[ arr_idx_row - oldind_row.begin() ];
    }

    // quadratic coefficients need be changed one at a time and adding only 
    // the difference between the previous and the new value

    double q_diff = std::get< 2 >( cp[ idx ] ) - old_var_q_coeff;

    GRBaddqpterms( model, 1 , & cidx , & cidx , & q_diff );
    }

  auto nsz = std::distance( nval.begin() , nvit );
  cidx.resize( nsz );
  nval.resize( nsz );

  GRBsetdblattrlist( model , GRB_DBL_ATTR_OBJ , cidx.size() , cidx.data() , nval.data());

  GRBupdatemodel( model );
  return;
  }

 // Fallback method - Update all costs
 // --------------------------------------------------------------------------
 // reload_objective( f );

 }  // end( GRBMILPSolver::objective_function_modification )

/*--------------------------------------------------------------------------*/

void GRBMILPSolver::constraint_function_modification( const FunctionMod *mod )
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
   *(cidxit++) = grb_index_of_variable( static_cast< const ColVariable * >( v ) );
   }

 auto nsz = std::distance( nval.begin() , nvit );
 cidx.resize( nsz );
 nval.resize( nsz );

 std::vector< int > rows( nsz , row );
 GRBchgcoeffs( model , nsz , rows.data() , cidx.data() , nval.data() );
 GRBupdatemodel( model );

 // Fallback method - Reload all coefficients
 // --------------------------------------------------------------------------
 // reload_constraint( lf );

 }  // end( GRBMILPSolver::constraint_function_modification )

/*--------------------------------------------------------------------------*/

void GRBMILPSolver::objective_fvars_modification( const FunctionModVars *mod )
{
 // this is called in response to Variable being added to / removed from the
 // Objective; however, note that all Variable are supposed to exist at the
 // time this is called, so grb_index_of_variable() is always correct

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

 auto f = mod->function();
 auto nav = f->get_num_active_var();

 // while changing the coefficients, we have to be careful about the fact
 // that Modification are managed asynchronously with the model changes
 // although the added/removed Variable do exist in the internal data
 // structure of [GRB]MILPSolver since the Modification are managed
 // strictly in arrival order, they may no longer exist in the model;
 // more to the point, they may no longer be active in the LinearFunction

 if( auto lf = dynamic_cast< const LinearFunction * >( f ) ) {
  // Linear objective function

  for( auto v : mod->vars() ) {
   auto var = static_cast< const ColVariable * >( v );
   if( auto idx = grb_index_of_variable( var ) ; idx < Inf< int >() ) {
    double value = 0;
    
    if( mod->added() ) {
     auto cidx = lf->is_active( var );
     value = cidx < nav ? lf->get_coefficient( cidx ) : 0;
     }

    GRBsetdblattrelement( model , GRB_DBL_ATTR_OBJ , idx , value  );
    }
   }
  
  return;
  }

 if( auto qf = dynamic_cast< const DQuadFunction * >( f ) ) {
  // Quadratic objective function

  // In Gurobi to change quadratic coefficients we need to retrieve all the old coeff.,
  // and then add the difference between the new and the old ones. In this case there if we want to delete
  // a quadratic coefficient, we can just subtract its old value

  int nqz;
  GRBgetintattr( model , GRB_INT_ATTR_NUMQNZS , & nqz );

  std::vector< int > oldind_row ( nqz );
  std::vector< int > oldind_col ( nqz );
  std::vector< double > oldval ( nqz );

  int status = GRBgetq( model, & nqz , oldind_row.data() , oldind_col.data() , oldval.data() );
  if( status != 0 )
   throw( std::runtime_error( "Error while querying quadratic coefficients with GRBgetq" ) );

  for( auto v : mod->vars() ) {
   auto var = static_cast< const ColVariable * >( v );
   if( auto ind = grb_index_of_variable( var ) ; ind < Inf< int >() ) {
    double value = 0;
    double q_value = 0;

    if( mod->added() ) {
     if( auto idx = qf->is_active( var ) ; idx < nav ) {
       value = qf->get_linear_coefficient( idx );
       q_value = qf->get_quadratic_coefficient( idx );
       }
     }
    else { // removed variable
       auto arr_idx_row = std::find(oldind_row.begin(), oldind_row.end(), ind);
       auto arr_idx_col = std::find(oldind_col.begin(), oldind_col.end(), ind);

       if( *arr_idx_row != *arr_idx_col )
       throw( std::runtime_error( "Error while modifying quadratic coefficients" ) );
       
       q_value = - oldval[ *arr_idx_row ];
      }

    GRBsetdblattrelement( model , GRB_DBL_ATTR_OBJ , ind , value );
    GRBaddqpterms( model, 1 , & ind , & ind , & q_value );
    }
   }

  GRBupdatemodel( model );
  return;
  }

 // This should never happen
 throw( std::invalid_argument( "Unknown type of Objective Function" ) );

 }  // end( GRBMILPSolver::objective_fvars_modification )

/*--------------------------------------------------------------------------*/

void GRBMILPSolver::constraint_fvars_modification(
						const FunctionModVars * mod )
{
 // this is called in response to Variable being added to / removed from the
 // Constraint; however, note that all Variable are supposed to exist at the
 // time this is called, so grb_index_of_variable() is always correct

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
 // structure of [GRB]MILPSolver since the Modification are managed
 // strictly in arrival order, they may no longer exist in the model;
 // more to the point, they may no longer be active in the LinearFunction

 // get indices and coefficients
 for( auto v : mod->vars() ) {
  auto var = static_cast< const ColVariable * >( v );
  if( auto vidx = grb_index_of_variable( var ) ; vidx < Inf< int >() ) {

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
 GRBchgcoeffs( model , indices.size() , rows.data() , 
                indices.data() , values.data() );
 GRBupdatemodel( model );

 }  // end( GRBMILPSolver::constraint_fvars_modification )

/*----------------------------------------------------------------------------

void GRBMILPSolver::dynamic_modification( const BlockModAD * mod )
{
 MILPSolver::dynamic_modification( mod );
 }

----------------------------------------------------------------------------*/

void GRBMILPSolver::add_dynamic_constraint( const FRowConstraint * con )
{
 // call the method of MILPSolver to update the dictionaries (only)
 MILPSolver::add_dynamic_constraint( con );

 auto lf = dynamic_cast< const LinearFunction * >( con->get_function() );
 if( ! lf )
  throw( std::invalid_argument( "the FRowConstraint is not linear" ) );

 int nzcnt = lf->get_num_active_var();

 int n_ranged_con = map_rng_con_aux_var.size();
 std::array< int , 2 > rmatbeg = { 0 , nzcnt };
 std::vector< int > rmatind;
 rmatind.reserve( nzcnt );
 std::vector< double > rmatval;
 rmatval.reserve( nzcnt );

 // get the coefficients to fill the matrix
 for( auto & el : lf->get_v_var() )
  if( auto idx = grb_index_of_variable( el.first ) ; idx < Inf< int >() ) {
   rmatind.push_back( idx );
   rmatval.push_back( el.second );
   }

 // get the bounds
 auto con_lhs = con->get_lhs();
 auto con_rhs = con->get_rhs();
 double rhs , lhs , rngval;
 char sense;

 if( con_lhs == con_rhs ) {
  sense = GRB_EQUAL;
  rhs = con_rhs;
  }
 else
  if( con_lhs == -Inf< double >() ) {
   sense = GRB_LESS_EQUAL;
   rhs = con_rhs;
   }
  else
   if( con_rhs == Inf< double >() ) {
    sense = GRB_GREATER_EQUAL;
    rhs = con_lhs;
    }
   else {
    sense = 'R';
    lhs = con_lhs;
    rhs = con_rhs;
    }

 // update the GUROBI problem
 if( sense != 'R' )
  GRBaddconstr( model , rmatind.size() , rmatind.data() , 
                  rmatval.data() , sense , rhs ,
                  NULL );
 else{
   // Filling map between ranged constraint and auxiliary variables built by Gurobi
   // See GRBMILPSolver.h for further informations
   map_rng_con_aux_var.push_back( { numrows - 1 , numcols + n_ranged_con } );
   GRBaddrangeconstr( model , rmatind.size() , rmatind.data() , 
                         rmatval.data() , lhs , rhs ,
                         NULL );
  }

 GRBupdatemodel( model );

 }  // end( GRBMILPSolver::add_dynamic_constraint )

/*--------------------------------------------------------------------------*/

void GRBMILPSolver::add_dynamic_variable( const ColVariable * var )
{
 bool was_mip = int_vars > 0;  // was a MIP before the change

 // call the method of MILPSolver to update dictionaries (and int_vars)
 MILPSolver::add_dynamic_variable( var );

 bool is_mip = int_vars > 0;  // is a MIP after the change

 // get the bounds
 auto bd = GRBMILPSolver::get_problem_bounds( *var );

 char new_ctype;  // get the new variable type
 
 if( var->is_integer() && ( ! relax_int_vars ) ) {
   if( var->is_unitary() && var->is_positive() )
    new_ctype = GRB_BINARY;  // Binary
   else
    new_ctype = GRB_INTEGER;  // Integer
   }
  else
   new_ctype = GRB_CONTINUOUS;   // Continuous

 // update the GUROBI problem
 GRBaddvar( model , 0 , nullptr , nullptr , 0.0 , 
            bd[ 0 ] , bd[ 1 ] , new_ctype , nullptr );
 
 GRBupdatemodel( model );

 }  // end( GRBMILPSolver::add_dynamic_variable )

/*--------------------------------------------------------------------------*/

void GRBMILPSolver::add_dynamic_bound( const OneVarConstraint * con )
{
 // no point in calling the method of MILPSolver, as it does nothing
 // MILPSolver::add_dynamic_bound( con );

 auto var = static_cast< const ColVariable * >( con->get_active_var( 0 ) );
 if( ! var )
  throw( std::logic_error( "GRBMILPSolver: added a bound on no Variable" ) );

 auto idx = grb_index_of_variable( var );
 if( idx == Inf< int >() )
  throw( std::logic_error( "GRBMILPSolver: added a bound on unknown Variable"
			   ) );

 auto bd = GRBMILPSolver::get_problem_bounds( *var );

 GRBsetdblattrelement( model , GRB_DBL_ATTR_LB , idx , bd[ 0 ] );
 GRBsetdblattrelement( model , GRB_DBL_ATTR_UB , idx , bd[ 1 ] );
 }

/*--------------------------------------------------------------------------*/

void GRBMILPSolver::remove_dynamic_constraint( const FRowConstraint * con )
{
 int index = index_of_dynamic_constraint( con );
 if( index == Inf< int >() )
  throw( std::runtime_error( "Dynamic constraint not found" ) );

 int n_ranged_con = map_rng_con_aux_var.size();
 if( n_ranged_con != 0 ) {
  // find if con is a ranged constraint
  auto it_rng = std::find_if( map_rng_con_aux_var.begin() + last_static_rng_con + 1, 
      map_rng_con_aux_var.end(), [&index]( std::pair< int , int > const& elem ) {
      return( elem.first == index );
    });
  bool is_rng = ( it_rng != map_rng_con_aux_var.end() ); // 0 isn't a ranged constraint

  // The element ( index , aux_var ) has to be removed from the map and also the 
  // auxiliary variable has to be removed from the Gurobi model
  if( is_rng ) {
    int index_aux_var = (*it_rng).second;
    GRBdelvars( model , 1 , &index_aux_var );
    map_rng_con_aux_var.erase( it_rng );

    // Update map : find the first pair with idx aux var greater than index
    auto it_rng_var = std::find_if( map_rng_con_aux_var.begin(), map_rng_con_aux_var.end(), 
      [&index_aux_var]( std::pair< int , int > const& elem ) {
      return( elem.second > index_aux_var );
    });
    // Update map : decrease the idx of aux var
    while( it_rng_var != map_rng_con_aux_var.end() ) {
      --( *it_rng_var ).second;
      ++it_rng_var;
    }
  }

  // Update map : find the first pair with idx con greater than index
  auto it_rng_s = std::find_if( map_rng_con_aux_var.begin() + last_static_rng_con + 1,
      map_rng_con_aux_var.end(), [&index]( std::pair< int , int > const& elem ) {
      return( elem.first > index );
    });
  // Update map : decrease the idx of rng con
  while( it_rng_s != map_rng_con_aux_var.end() ) {
    --( *it_rng_s ).first;
    ++it_rng_s;
  }
 }

 GRBdelconstrs( model , 1 , &index );
 GRBupdatemodel( model );

 // call the method of MILPSolver to update the dictionaries (only)
 MILPSolver::remove_dynamic_constraint( con );
 }

/*--------------------------------------------------------------------------*/

void GRBMILPSolver::remove_dynamic_variable( const ColVariable * var )
{
 int index = grb_index_of_dynamic_variable( var );
 int n_ranged_con = map_rng_con_aux_var.size();

 if( n_ranged_con != 0 ) {
 // Update map : find the first pair with idx aux var greater than index
  auto it_rng = std::find_if( map_rng_con_aux_var.begin(), map_rng_con_aux_var.end(), 
      [&index]( std::pair< int , int > const& elem ) {
      return( elem.second > index );
    });
  // Update map : decrease the idx of aux var
  while( it_rng != map_rng_con_aux_var.end() ) {
    --( *it_rng ).second;
    ++it_rng;
  }
 }

 GRBdelvars( model , 1 , &index );
 GRBupdatemodel( model );

 // call the method of MILPSolver to update the dictionaries (only)
 MILPSolver::remove_dynamic_variable( var );
 }

/*--------------------------------------------------------------------------*/

void GRBMILPSolver::remove_dynamic_bound( const OneVarConstraint * con )
{
 // no point in calling the method of MILPSolver, as it does nothing
 // MILPSolver::remove_dynamic_bound( con );

 // note: this only works because remove_dynamic_constraint[s]() do *not*
 //       clear the removed OneVarConstraint, and therefore we can easily
 //       reconstruct which ColVariable it was about
 auto var = static_cast< const ColVariable * >( con->get_active_var( 0 ) );
 if( ! var )  // this should never happen
  return;     // but in case, there is nothing to do

 int idx = grb_index_of_variable( var );
 if( idx == Inf< int >() )  // the ColVariable has been removed
  return;                   // is strange, but there is nothing to do

 auto bd = GRBMILPSolver::get_problem_bounds( *var );

 GRBsetdblattrelement( model , GRB_DBL_ATTR_LB , idx , bd[ 0 ] );
 GRBsetdblattrelement( model , GRB_DBL_ATTR_UB , idx , bd[ 1 ] );
 }

/*--------------------------------------------------------------------------*/

int GRBMILPSolver::callback( GRBmodel *model,
           void *cbdata,
           int where )
{
 // main switch: depending on where - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 switch( where ) {
  case( GRB_CB_POLLING ): break; /* Ignore polling callback */
  case( GRB_CB_PRESOLVE ): break; /* Ignore presolve callback */
  case( GRB_CB_SIMPLEX ):
   // Currently in simplex- - - - - - - - - - - - - - - - - - - - - - - - 
   // check upper / lower bounds and in case stop ??
   break;
  case( GRB_CB_MIP ): {
   // Currently in MIP- - - - - - - - - - - - - - - - - - - - - - - -
   // check upper / lower bounds and in case stop

   double solv , bndv;
   GRBcbget( cbdata , where , GRB_CB_MIP_OBJBST , & solv );
   GRBcbget( cbdata , where , GRB_CB_MIP_OBJBND , & bndv );

   if( get_objsense() == 1 ) {
    // a minimization problem: solv is upper bound and bndv is lower bound
    if( solv >= 1e+75 )
     solv = Inf< double >();

    if( bndv <= - 1e+75 )
     bndv = - Inf< double >();

    if( ( bndv >= up_cut_off() ) || ( solv <= lw_cut_off() ) )
     GRBterminate( model );
    }
   else {
    // a maximization problem: solv is lower bound and bndv is upper bound
    if( solv <= -1e+75 )
     solv = - Inf< double >();

    if( bndv >= 1e+75 )
     bndv = Inf< double >();

    if( ( solv >= up_cut_off() ) || ( bndv <= lw_cut_off() ) )
     GRBterminate( model );
    }
   break;
   }

  case( GRB_CB_MIPNODE ): {
   // MIP node callback - - - - - - - - - - - - - - - - - - - - - -
   if( ! ( CutSepPar & 3 ) )  // but we don't do user cut separation
    break;                    // nothing to do

   double depth; // find the depth of the current node
   if( GRBcbget( cbdata , where , GRB_CB_MIPNODE_NODCNT , & depth ) )
    throw( std::runtime_error(
                "Unable to get the depth with GRB_CB_MIPNODE_NODCNT" ) );

   int status;
   GRBcbget( cbdata , where , GRB_CB_MIPNODE_STATUS , & status);
   if(status == GRB_OPTIMAL) {

    // if we are at a depth for which separation is not enabled
    if( ( ( ! depth ) && ( ! ( CutSepPar & 1 ) ) ) ||
        ( depth && ( ! ( CutSepPar & 2 ) ) ) )
     break;                    // nothing to do

    // this is a critical section where different GUROBI threads may compete
    // for access to the Block: ensure mutual exclusion
    f_callback_mutex.lock();

    // ensure no interference from other threads (except GUROBI ones) by also
    // lock()-ing the Block
    bool owned = f_Block->is_owned_by( f_id );
    if( ( ! owned ) && ( ! f_Block->lock( f_id ) ) )
     throw( std::runtime_error( "Unable to lock the Block" ) );

    // get the solution of the relaxation
    std::vector< double > x( numcols );
    if( GRBcbget( cbdata , where, GRB_CB_MIPNODE_REL , x.data() ) )
     throw( std::runtime_error(
        "Unable to get the solution with GRB_CB_MIPNODE_REL" ) );
    // GRBcbsolution( cbdata, x.data() , nullptr);

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
      for( size_t c = 0 ; c < rhs.size() ; ++c ) {
        int nnz; // number of nonzero coefficients in the actual cut
        int idx = rmatbeg[ c ];
        if( c < rhs.size() - 1)
          nnz = rmatbeg[ c + 1 ] - rmatbeg[ c ]; 
        else
          nnz = rmatind.size() - rmatbeg[ c ];

        if( GRBcbcut( cbdata , nnz , & rmatind[ idx ] ,
            & rmatval[ idx ] , sense[ c ] , rhs[ c ] ) )
          throw( std::logic_error( "problem in GRBcbcut" ) );

        }
      }
    }

   break;
   }
  case( GRB_CB_MIPSOL ): {
   // a feasible solution has been found- - - - - - - - - - - - - - - - - - -
   if( ! ( CutSepPar & 4 ) )  // but we don't do lazy constraint separation
    break;                    // nothing to do

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
   if( GRBcbget( cbdata , where , GRB_CB_MIPSOL_SOL, x.data() ) )
    throw( std::runtime_error(
       "Unable to get the solution with GRB_CB_MIPSOL_SOL" ) );

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
   if( ! rmatbeg.empty() ) {
    for( size_t c = 0 ; c < rhs.size() ; ++c ) {
      int nnz; // number of nonzero coefficients in the actual lazy constraint
      int idx = rmatbeg[ c ];
      if( c < rhs.size() - 1)
        nnz = rmatbeg[ c + 1 ] - rmatbeg[ c ]; 
      else
        nnz = rmatind.size() - rmatbeg[ c ];

      if( GRBcblazy( cbdata , nnz , & rmatind[ idx ] ,
                      & rmatval[ idx ] , sense[ c ] , rhs[ c ] ) )
       throw( std::logic_error( "problem in GRBcblazy" ) );
      }
    }
    
    break;
    }
  case( GRB_CB_MESSAGE ): break;
  case( GRB_CB_BARRIER ): break;
  case( GRB_CB_MULTIOBJ ): break;
  }  // end( main switch )- - - - - - - - - - - - - - - - - - - - - - - - - -
     // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 return( 0 );
 }

/*--------------------------------------------------------------------------*/

void GRBMILPSolver::perform_separation( Configuration * cfg ,
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
    *(iit++) = grb_index_of_variable( el.first );
    *(vit++) = el.second;
    }

   // get the bounds
   auto con_lhs = con->get_lhs();
   auto con_rhs = con->get_rhs();

   if( con_lhs == con_rhs ) {
    sense.push_back( GRB_EQUAL );
    rhs.push_back( con_rhs );
    }
   else
    if( con_lhs == -Inf< double >() ) {
     sense.push_back( GRB_LESS_EQUAL );
     rhs.push_back( con_rhs );
     }
    else
     if( con_rhs == Inf< double >() ) {
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

   rmatbeg.push_back( rmatind.size() );

   }  // end( for each added FRowConstraint )
  }  // end( main loop )
 }  // end( GRBMILPSolver::perform_separation )

/*--------------------------------------------------------------------------*/

std::string GRBMILPSolver::grb_int_par_map( idx_type par ) const
{
 switch( par ) {
  case( intMaxIter ): return( GRB_DBL_PAR_NODELIMIT );
  case( intMaxSol ):  return( GRB_INT_PAR_SOLUTIONLIMIT );
  case( intLogVerb ): return( GRB_INT_PAR_LOGTOCONSOLE );
  case( intMaxThread ): return( GRB_INT_PAR_THREADS );
  }

 // GUROBI parameters
 if( ( par >= intFirstGUROBIPar ) && ( par < intLastAlgParGRBS ) ) {
  return( SMSpp_to_GUROBI_int_pars[ par - intFirstGUROBIPar ] );

  }

 return( "" );
 }

/*--------------------------------------------------------------------------*/

std::string GRBMILPSolver::grb_dbl_par_map( idx_type par ) const
{
 switch( par ) {
  case( dblMaxTime ): return( GRB_DBL_PAR_TIMELIMIT );
  case( dblRelAcc ):  return( GRB_DBL_PAR_MIPGAP );
  case( dblAbsAcc ):  return( GRB_DBL_PAR_MIPGAPABS );
  case( dblRAccSol ): return( GRB_DBL_PAR_POOLGAP );
  case( dblAAccSol ): return( GRB_DBL_PAR_POOLGAPABS );
  case( dblFAccSol ): return( GRB_DBL_PAR_FEASIBILITYTOL );
  }

 if( ( par >= dblFirstGUROBIPar ) && ( par < dblLastAlgParGRBS ) )
  return( SMSpp_to_GUROBI_dbl_pars[ par - dblFirstGUROBIPar ] );

 return( "" );
 }

 /*--------------------------------------------------------------------------*/
/*------------------- METHODS FOR HANDLING THE PARAMETERS ------------------*/
/*--------------------------------------------------------------------------*/

void GRBMILPSolver::set_par( idx_type par , int value )
{
 if( par == intThrowReducedCostException ) {
  throw_reduced_cost_exception = bool( value );
  return;
  }

 if( par == intCutSepPar ) {
  CutSepPar = value;
  return;
  }

 if( par == intMaxIter ) { // intMaxIter is an int parameter in sms++ but a double in Gurobi
  set_par( par , (double)value );
  return;
 }

 std::string gp = grb_int_par_map( par );
 if( gp.size() > 0 ) {
  GRBsetintparam( env , gp.c_str() , value );
  return;
  }
 //else
  //throw( std::invalid_argument( "Parameter " + int_par_idx2str(par) + " not correctly converted in Gurobi" ) );
  

 MILPSolver::set_par( par, value );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

void GRBMILPSolver::set_par( idx_type par , double value )
{
 // Solver parameters explicitly mapped in GUROBI
 switch( par ) {
  case( dblUpCutOff ): UpCutOff = value; return;
  case( dblLwCutOff ): LwCutOff = value; return;
  }

 std::string gp;
 if( par == intMaxIter ) // intMaxIter is an int parameter in sms++ but a double in Gurobi
  gp = grb_int_par_map( par );
 else
  gp = grb_dbl_par_map( par );

 if( gp.size() > 0 ) {
  GRBsetdblparam( env , gp.data() , value );
  return;
  }

 MILPSolver::set_par( par , value );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

void GRBMILPSolver::set_par( idx_type par , std::string && value )
{
 // GUROBI parameters
 if( ( par >= strFirstGUROBIPar ) && ( par < strLastAlgParGRBS ) ) {
  std::string gurobi_par = SMSpp_to_GUROBI_str_pars[ par - strFirstGUROBIPar ];
  GRBsetstrparam( env , gurobi_par.data() , value.c_str() );
  return;
  }

 MILPSolver::set_par( par, std::move( value ) );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

void GRBMILPSolver::set_par( idx_type par , std::vector< int > && value )
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

void GRBMILPSolver::set_par( idx_type par ,
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

Solver::idx_type GRBMILPSolver::get_num_int_par( void ) const {
 return( MILPSolver::get_num_int_par()
	 + intLastAlgParGRBS - intLastAlgParMILP );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

Solver::idx_type GRBMILPSolver::get_num_dbl_par( void ) const {
 return( MILPSolver::get_num_dbl_par()
	 + dblLastAlgParGRBS - dblLastAlgParMILP );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

Solver::idx_type GRBMILPSolver::get_num_str_par( void ) const {
 return( MILPSolver::get_num_str_par()
	 + strLastAlgParGRBS - strLastAlgParMILP );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

Solver::idx_type GRBMILPSolver::get_num_vint_par( void ) const {
 return( MILPSolver::get_num_vint_par()
	 + vintLastAlgParGRBS - vintLastAlgParMILP );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

Solver::idx_type GRBMILPSolver::get_num_vstr_par( void ) const {
 return( MILPSolver::get_num_vstr_par()
	 + vstrLastAlgParGRBS - vstrLastAlgParMILP );
 }

/*--------------------------------------------------------------------------*/

int GRBMILPSolver::get_dflt_int_par( idx_type par ) const
{
 if( ( par == intThrowReducedCostException ) || ( par == intCutSepPar ) )
  return( 0 );

 std::string gp = grb_int_par_map( par );
 if( gp.size() > 0 ) {
   int value;
   GRBgetintparam( env , gp.data() , & value );
   return( value );
  }

 return( MILPSolver::get_dflt_int_par( par ) );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

double GRBMILPSolver::get_dflt_dbl_par( idx_type par ) const
{
 switch( par ) {
  case( dblUpCutOff ): return( Inf< double >() );
  case( dblLwCutOff ): return( - Inf< double >() );
  }

 std::string gp = grb_dbl_par_map( par );
 if( gp.size() > 0 ) {
  double value;
  GRBgetdblparam( env , gp.data() , & value );
  return( value );
  }

 return( MILPSolver::get_dflt_dbl_par( par ) );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

const std::string & GRBMILPSolver::get_dflt_str_par( idx_type par ) const
{
 // note: this implementation is not thread safe and it may lead to elements
 //       of value[] to be allocated more than once with some memory being
 //       lost, but the chances are too slim and the potential drawback too
 //       limited to warrant even a humble std::atomic_flag
 static std::vector< std::string > value( strLastAlgParGRBS -
					  strFirstGUROBIPar );

 if( ( par >= strFirstGUROBIPar ) && ( par < strLastAlgParGRBS ) ) {
  auto i = par - strFirstGUROBIPar;
  if( value[ i ].empty() ) {
   value[ i ].reserve( 512 );
   GRBgetstrparam( env , SMSpp_to_GUROBI_str_pars[ i ].data() , value[ i ].data() );
   }

  return( value[ i ] );
  }

 return( MILPSolver::get_dflt_str_par( par ) );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

const std::vector< int > & GRBMILPSolver::get_dflt_vint_par( idx_type par )
 const
{
 static std::vector< int > _empty;
 if( par == vintCutSepCfgInd )
  return( _empty );

 return( MILPSolver::get_dflt_vint_par( par ) );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

const std::vector< std::string > & GRBMILPSolver::get_dflt_vstr_par(
							idx_type par ) const
{
 static std::vector< std::string > _empty;
 if( par == vstrConfigDBFName )
  return( _empty );

 return( MILPSolver::get_dflt_vstr_par( par ) );
 }

/*--------------------------------------------------------------------------*/

int GRBMILPSolver::get_int_par( idx_type par ) const
{
 if( par == intThrowReducedCostException )
  return( throw_reduced_cost_exception );

 if( par == intCutSepPar )
  return( CutSepPar );

 std::string gp = grb_int_par_map( par );
  if( gp.size() > 0 ) {
   int value;
   GRBgetintparam( env , gp.data() , & value );
   return( value );
  }

 return( MILPSolver::get_int_par( par ) );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

double GRBMILPSolver::get_dbl_par( idx_type par ) const
{
 switch( par ) {
  case( dblUpCutOff ): return( UpCutOff );
  case( dblLwCutOff ): return( LwCutOff );
  }

 std::string gp = grb_dbl_par_map( par );
 if( gp.size() > 0 ) {
  double value;
  GRBgetdblparam( env , gp.data() , & value );
  return( value );
  }

 return( MILPSolver::get_dbl_par( par ) );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

const std::string & GRBMILPSolver::get_str_par( idx_type par ) const
{
 static std::string value;

 if( ( par >= strFirstGUROBIPar ) && ( par < strLastAlgParGRBS ) ) {
  std::string gurobi_par = SMSpp_to_GUROBI_str_pars[ par - strFirstGUROBIPar ];
  value.reserve( 512 );
  GRBgetstrparam( env , gurobi_par.data() , value.data() );
  return( value );
  }

 return( MILPSolver::get_str_par( par ) );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

const std::vector< int > & GRBMILPSolver::get_vint_par( idx_type par ) const
{
 if( par == vintCutSepCfgInd )
  return( CutSepCfgInd );

 return( MILPSolver::get_vint_par( par ) );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

const std::vector< std::string > & GRBMILPSolver::get_vstr_par( idx_type par )
 const
{
 if( par == vstrConfigDBFName )
  return( ConfigDBFName );

 return( MILPSolver::get_vstr_par( par ) );
 }

/*--------------------------------------------------------------------------*/

Solver::idx_type GRBMILPSolver::int_par_str2idx(
					     const std::string & name ) const
{
 if( name == "intThrowReducedCostException" )
  return( intThrowReducedCostException );

 if( name == "intCutSepPar" )
  return( intCutSepPar );

 /* In GRBMILPSolver::*_par_str2idx() methods we check with MILPSolver first */

 idx_type idx = MILPSolver::int_par_str2idx( name );
 if( idx < Inf< idx_type >() )
  return( idx );

 // GUROBI parameters
 std::string gurobi_par = name;
 auto array_pos = std::find( SMSpp_to_GUROBI_int_pars.begin() ,
                        SMSpp_to_GUROBI_int_pars.end() ,
                        gurobi_par);

 if( array_pos != SMSpp_to_GUROBI_int_pars.end() ) {
  int pos = std::distance( SMSpp_to_GUROBI_int_pars.begin(), array_pos );
  auto idx_par = GUROBI_to_SMSpp_int_pars[pos].second;
  return( idx_par );
  }

 return( Inf< idx_type >() );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

const std::string & GRBMILPSolver::int_par_idx2str( idx_type idx ) const
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
 par_name.reserve( 512 );

 if( ( idx >= intFirstGUROBIPar ) && ( idx < intLastAlgParGRBS ) ) {
  par_name = SMSpp_to_GUROBI_int_pars[ idx - intFirstGUROBIPar ];
  return( par_name );
  }

 return( MILPSolver::int_par_idx2str( idx ) );
 }

/*--------------------------------------------------------------------------*/

Solver::idx_type GRBMILPSolver::dbl_par_str2idx( const std::string & name )
 const
{
 /* In GRBMILPSolver::*_par_str2idx() methods we check with MILPSolver first */

 idx_type idx = MILPSolver::dbl_par_str2idx( name );
 if( idx < Inf< idx_type >() )
  return( idx );

 // GUROBI parameters
 std::string gurobi_par = name;
 auto array_pos = std::find( SMSpp_to_GUROBI_dbl_pars.begin() ,
                        SMSpp_to_GUROBI_dbl_pars.end() ,
                        gurobi_par);

 if( array_pos != SMSpp_to_GUROBI_dbl_pars.end() ) {
  int pos = std::distance( SMSpp_to_GUROBI_dbl_pars.begin(), array_pos );
  auto idx_par = GUROBI_to_SMSpp_dbl_pars[pos].second;
  return( idx_par );
  }

 return( Inf< idx_type >() );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

const std::string & GRBMILPSolver::dbl_par_idx2str( idx_type idx ) const
{
 // note: this implementation is not thread safe and it requires that the
 //       result is used immediately after the call (prior to any other call
 //       to int_par_idx2str()), this may have to be improved upon
 static std::string par_name;
 par_name.reserve( 512 );

 if( ( idx >= dblFirstGUROBIPar ) && ( idx < dblLastAlgParGRBS ) ) {
  par_name = SMSpp_to_GUROBI_dbl_pars[ idx - dblFirstGUROBIPar ];
  return( par_name );
  }

 return( MILPSolver::dbl_par_idx2str( idx ) );
 }

/*--------------------------------------------------------------------------*/

Solver::idx_type GRBMILPSolver::str_par_str2idx( const std::string & name )
 const
{
 /* In GRBMILPSolver::*_par_str2idx() methods we check with MILPSolver first */

 idx_type idx = MILPSolver::str_par_str2idx( name );
 if( idx < Inf< idx_type >() )
  return( idx );

 // GUROBI parameters
 std::string gurobi_par = name;
 auto array_pos = std::find( SMSpp_to_GUROBI_str_pars.begin() ,
                        SMSpp_to_GUROBI_str_pars.end() ,
                        gurobi_par);

 if( array_pos != SMSpp_to_GUROBI_str_pars.end() ) {
  int pos = std::distance( SMSpp_to_GUROBI_str_pars.begin(), array_pos );
  auto idx_par = GUROBI_to_SMSpp_str_pars[pos].second;
  return( idx_par );
  }

 return( Inf< idx_type >() );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

const std::string & GRBMILPSolver::str_par_idx2str( idx_type idx ) const
{
 // note: this implementation is not thread safe and it requires that the
 //       result is used immediately after the call (prior to any other call
 //       to int_par_idx2str()), this may have to be improved upon
 static std::string par_name;
 par_name.reserve( 512 );

 if( ( idx >= strFirstGUROBIPar ) && ( idx < strLastAlgParGRBS ) ) {
  par_name = SMSpp_to_GUROBI_str_pars[ idx - strFirstGUROBIPar ];
  return( par_name );
  }

 return( MILPSolver::str_par_idx2str( idx ) );
 }

/*--------------------------------------------------------------------------*/

Solver::idx_type GRBMILPSolver::vint_par_str2idx(
					     const std::string & name ) const
{
 if( name == "vintCutSepCfgInd" )
  return( vintCutSepCfgInd );

 return( MILPSolver::vint_par_str2idx( name ) );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

const std::string & GRBMILPSolver::vint_par_idx2str( idx_type idx ) const
{
 static const std::string _pars = "vintCutSepCfgInd";
 if( idx == vintCutSepCfgInd )
  return( _pars );

 return( MILPSolver::vint_par_idx2str( idx ) );
 }

/*--------------------------------------------------------------------------*/

Solver::idx_type GRBMILPSolver::vstr_par_str2idx(
					     const std::string & name ) const
{
 if( name == "vstrConfigDBFName" )
  return( vstrConfigDBFName );

 return( MILPSolver::vstr_par_str2idx( name ) );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

const std::string & GRBMILPSolver::vstr_par_idx2str( idx_type idx ) const
{
 static const std::string _pars = "vstrConfigDBFName";
 if( idx == vstrConfigDBFName )
  return( _pars );

 return( MILPSolver::vstr_par_idx2str( idx ) );
 }

/*--------------------------------------------------------------------------*/

#ifdef MILPSOLVER_DEBUG

void GRBMILPSolver::check_status( void )
{
 int nvars;
 GRBgetintattr( model , GRB_INT_ATTR_NUMVARS , &nvars );
 if( numcols != nvars )
  DEBUG_LOG( "numcols is " << numcols << " but GRB_INT_ATTR_NUMVARS returns "
	     << nvars << std::endl );

 int nconstr;
 GRBgetintattr( model , GRB_INT_ATTR_NUMCONSTRS , &nconstr );
 if( numrows != nconstr )
  DEBUG_LOG( "numrows is " << numrows << " but GRB_INT_ATTR_NUMCONSTRS returns "
	     << nconstr << std::endl );

 int nbin , nint;
 GRBgetintattr( model , GRB_INT_ATTR_NUMINTVARS , &nint );
 GRBgetintattr( model , GRB_INT_ATTR_NUMBINVARS , &nbin );
 if( int_vars != nint + nbin )
  DEBUG_LOG( "int_vars is " << int_vars << " but GUROBI has actually "
	     << nint + nbin << " integer variables" << std::endl );

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
  if( auto idx = grb_index_of_variable( var.first ) ; idx < Inf< int >() ) {
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
/*
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
 }*/

/*--------------------------------------------------------------------------*/

Configuration * GRBMILPSolver::get_cfg( Index ci ) const
{
 if( ci >= CutSepCfgInd.size() )
  return( nullptr );
 auto dbi = CutSepCfgInd[ ci ];
 if( ( dbi < 0 ) || ( Index( dbi ) >= v_ConfigDB.size() ) )
  return( nullptr );
 return( v_ConfigDB[ dbi ] );
 }

/*--------------------------------------------------------------------------*/
/*--------------------- End File GRBMILPSolver.cpp -------------------------*/
/*--------------------------------------------------------------------------*/
