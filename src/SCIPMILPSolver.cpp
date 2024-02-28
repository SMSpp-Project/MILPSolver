/*--------------------------------------------------------------------------*/
/*------------------------- File SCIPMILPSolver.cpp ------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Implementation of the SCIPMILPSolver class.
 *
 * \author Antonio Frangioni \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \author Niccolo' Iardella \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 * 
 * \author Enrico Calandrini \n
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

#include <functional>
#include <cstdio>
#include <cmath>

#include <LinearFunction.h>
#include <DQuadFunction.h>

#include "SCIPMILPSolver.h"

#include <scip/scipdefplugins.h>
#include <scip/cons_linear.h>

#ifdef MILPSOLVER_DEBUG
 #define DEBUG_LOG( stuff ) std::cout << "[MILPSolver DEBUG] " << stuff
#else
 #define DEBUG_LOG( stuff )
#endif

// include the proper SCIP parameter mapping
#include <boost/preprocessor/cat.hpp>
#include <boost/preprocessor/stringize.hpp>
#include BOOST_PP_STRINGIZE( BOOST_PP_CAT( BOOST_PP_CAT( SCIP, SCIP_VERSION ), _maps.h ) )

/*--------------------------------------------------------------------------*/
/*------------------------- NAMESPACE AND USING ----------------------------*/
/*--------------------------------------------------------------------------*/

using namespace SMSpp_di_unipi_it;

SMSpp_insert_in_factory_cpp_0( SCIPMILPSolver );

/*--------------------------------------------------------------------------*/
/*--------------------- CONSTRUCTOR AND DESTRUCTOR -------------------------*/
/*--------------------------------------------------------------------------*/

SCIPMILPSolver::SCIPMILPSolver() : 
MILPSolver() , f_callback_set( false ) , CutSepPar( 0 ) ,
throw_reduced_cost_exception( 0 ) , UpCutOff( Inf< double >() ) , 
LwCutOff( -Inf< double >() )
{
 SCIP_CALL_ABORT( SCIPcreate( & scip ) );
 SCIP_CALL_ABORT( SCIPincludeDefaultPlugins( scip ) );
 }

/*--------------------------------------------------------------------------*/

SCIPMILPSolver::~SCIPMILPSolver()
{
 for( auto el : v_ConfigDB )
  delete el;

 SCIPfree( & scip );
 }

/*--------------------------------------------------------------------------*/
/*--------------------- DERIVED METHODS OF BASE CLASS ----------------------*/
/*--------------------------------------------------------------------------*/

void SCIPMILPSolver::set_Block( Block * block )
{
 if( block == f_Block )
  return;

 MILPSolver::set_Block( block );
 UpCutOff = Inf< double >();
 LwCutOff = -Inf< double >();
 }

/*--------------------------------------------------------------------------*/

void SCIPMILPSolver::clear_problem( unsigned int what )
{
 MILPSolver::clear_problem( what );

 SCIP_CALL_ABORT( SCIPfreeProb( scip ) );
 }

/*--------------------------------------------------------------------------*/

void SCIPMILPSolver::load_problem( void )
{
 MILPSolver::load_problem();

 SCIP_CALL_ABORT( SCIPfreeProb( scip ) );
 SCIP_CALL_ABORT( SCIPcreateProbBasic( scip , prob_name.c_str() ) );

 // Set objective sense
 if( objsense == 1 )
  SCIP_CALL_ABORT( SCIPsetObjsense( scip , SCIP_OBJSENSE_MINIMIZE ) );
 else
  if( objsense == -1 )
   SCIP_CALL_ABORT( SCIPsetObjsense( scip , SCIP_OBJSENSE_MAXIMIZE ) );
  else
   SCIPABORT();

 // Add variables
 vars.resize( numcols );
 for( int i = 0 ; i < numcols ; ++i ) {

  SCIP_Real collb = ( lb[ i ] == -Inf< double >() ) ?
                    -SCIPinfinity( scip ) : lb[ i ];
  SCIP_Real colub = ( ub[ i ] == Inf< double >() ) ?
                    SCIPinfinity( scip ) : ub[ i ];

  SCIP_VARTYPE vartype;
  switch( xctype[ i ] ) {
   case( 'C' ): vartype = SCIP_VARTYPE_CONTINUOUS; break;
   case( 'B' ): vartype = SCIP_VARTYPE_BINARY;     break;
   case( 'I' ): vartype = SCIP_VARTYPE_INTEGER;    break;
   case( 'S' ):
   case( 'N' ):
   default:     SCIPABORT();
   }

  /* SCIP doesn't allows to declare binary variables with UB or LB
  *  not in [0,1] ( as it should be ... ). Thus, we have to check for
  *  this situation. */
 if( vartype == SCIP_VARTYPE_BINARY ){
  if( collb < 0 )
    collb = 0;
  if( collb > 1 )
    collb = 1;
  if( colub < 0 )
    colub = 0;
  if( colub > 1 )
    colub = 1;
  }

  SCIP_VAR * var = nullptr;
  char * name = use_custom_names ? colname[ i ] : nullptr;
  SCIP_CALL_ABORT( SCIPcreateVarBasic( scip , & var , name , collb ,  colub ,
                                       objective[ i ] , vartype ) );
  SCIP_CALL_ABORT( SCIPaddVar( scip , var ) );
  vars[ i ] = var;
  SCIP_CALL_ABORT( SCIPreleaseVar( scip , & var ) );
  }

 // add constraints
 cons.resize( numrows );
 for( int i = 0 ; i < numrows ; ++i ) {
  SCIP_Real con_lhs = NAN;
  SCIP_Real con_rhs = NAN;

  switch( sense[ i ] ) {
   case( 'L' ): con_lhs = -SCIPinfinity( scip );
                con_rhs = rhs[ i ];
	        break;
   case( 'E' ): con_lhs = con_rhs = rhs[ i ];
                break;
   case( 'G' ): con_lhs = rhs[ i ];
                con_rhs = SCIPinfinity( scip );
	        break;
   case( 'R' ): con_lhs = con_rhs = rhs[ i ];
                // TODO Check this
                if( rngval[ i ] > 0 )
		 con_rhs += rngval[ i ];
		else
		 con_lhs += rngval[ i ];
   }

  SCIP_CONS * con = nullptr;
  char * name = use_custom_names ? rowname[ i ] : nullptr;
  SCIP_CALL_ABORT( SCIPcreateConsBasicLinear( scip , & con , name , 0 ,
                                              nullptr , nullptr ,
                                              con_lhs , con_rhs ) );
  SCIP_CALL_ABORT( SCIPaddCons( scip , con ) );
  cons[ i ] = con;
  SCIP_CALL_ABORT( SCIPreleaseCons( scip , & con ) );
  }

 // add constraint coefficients
 for( int c = 0 ; c < numcols ; ++c )
  for( int i = matbeg[ c ] ; i < matbeg[ c + 1 ] ; ++i )
   SCIP_CALL_ABORT( SCIPaddCoefLinear( scip , cons[ matind[ i ] ] ,
                                       vars[ c ] , matval[ i ] ) );

 bool is_qp = std::any_of( q_objective.begin() , q_objective.end() ,
                           []( double d ) { return( d != 0 ); } );
 if( is_qp ) {
  /* Nonlinear objective functions are not supported by SCIP and must be 
  *  modeled as constraint function. Thus, for each variable v with nonzero
  *  quadratic coefficient we create a new aux_var z. In the objective z 
  *  goes with the q_objective value of v and a new quadratic
  *  aux_constraint is created such as: z - xˆ2 >= 0 */

  aux_vars.resize( numcols );
  aux_cons.resize( numcols );

  for( int i = 0 ; i < numcols ; ++i ) {
   if( q_objective[ i ] == 0 ) {
    aux_vars[ i ] = nullptr;
    aux_cons[ i ] = nullptr;
    continue;
    }

   // Add auxiliary variables
   SCIP_Real z_lb = - SCIPinfinity( scip );
   SCIP_Real z_ub = SCIPinfinity( scip ) ;

   SCIP_VARTYPE z_type = SCIP_VARTYPE_CONTINUOUS;

   SCIP_VAR * z = nullptr;
   SCIP_CALL_ABORT( SCIPcreateVarBasic( scip , & z , nullptr , 
                              z_lb  , z_ub ,  q_objective[ i ] , 
                              z_type ) );
   SCIP_CALL_ABORT( SCIPaddVar( scip , z ) );
   aux_vars[ i ] = z;
   SCIP_CALL_ABORT( SCIPreleaseVar( scip , & z ) );

   // Add auxiliary constraints z - xˆ2 >= 0
   SCIP_Real con_lhs = 0;
   SCIP_Real con_rhs = SCIPinfinity( scip );
   SCIP_Real lincoef = 1;
   SCIP_Real quadcoef = -1;

   SCIP_CONS * con = nullptr;
   SCIP_VAR * linvar = aux_vars[ i ];
   SCIP_VAR * quadvar = vars[ i ];
   const char * name = "aux_con";

   #if SCIP_VERSION < 800
    SCIP_CALL_ABORT( SCIPcreateConsBasicQuadratic( scip , & con , name ,
                     1 , & linvar , & lincoef , 1 , & quadvar , & quadvar ,
		     & quadcoef , con_lhs , con_rhs ) );
   #else
    SCIP_CALL_ABORT( SCIPcreateConsBasicQuadraticNonlinear( scip , & con ,
		     name , 1 , & linvar , & lincoef , 1 , & quadvar ,
		     & quadvar , & quadcoef , con_lhs , con_rhs ) );
   #endif

   SCIP_CALL_ABORT( SCIPaddCons( scip , con ) );
   aux_cons[ i ] = con;
   SCIP_CALL_ABORT( SCIPreleaseCons( scip , &con ) );

   // Add constraint coefficients
   // SCIP_CALL_ABORT( SCIPaddCoefLinear( scip,
   //                                     aux_cons[ i ],
   //                                     aux_vars[ i ],
   //                                     1 ) );
   // SCIP_CALL_ABORT( SCIPaddQuadVarQuadratic( scip,
   //                                           aux_cons[ i ],
   //                                           aux_vars[ i ],
   //                                           0,
   //                                           -1 ) );
   }
  }

 // the base representation isn't needed anymore
 MILPSolver::clear_problem( 15 );

 UpCutOff = Inf< double >();
 LwCutOff = -Inf< double >();

 }  // end( SCIPMILPSolver::load_problem )

/*--------------------------------------------------------------------------*/

double SCIPMILPSolver::get_problem_lb( const ColVariable & var ) const
{
 double b = MILPSolver::get_problem_lb( var );
 if( b == -Inf< double >() )
  b = -SCIPinfinity( scip );

 return( b );
 }

/*--------------------------------------------------------------------------*/

double SCIPMILPSolver::get_problem_ub( const ColVariable & var ) const
{
 double b = MILPSolver::get_problem_ub( var );
 if( b == Inf< double >() )
  b = SCIPinfinity( scip );

 return( b );
 }

/*--------------------------------------------------------------------------*/

std::array< double , 2 > SCIPMILPSolver::get_problem_bounds(
					      const ColVariable & var ) const
{
 auto ret = MILPSolver::get_problem_bounds( var );
 if( ret[ 0 ] == -Inf< double >() )
  ret[ 0 ] = -SCIPinfinity( scip );
 if( ret[ 1 ] == Inf< double >() )
  ret[ 1 ] = SCIPinfinity( scip );
 return( ret );
 }

/*--------------------------------------------------------------------------*/

int SCIPMILPSolver::compute( bool changedvars )
{
 lock();  // lock the mutex: this is done again inside MILPSolver::compute,
          // but that's OK since the mutex is recursive

 // process Modification: this is driven by MILPSolver- - - - - - - - - - - -
 if( MILPSolver::compute( changedvars ) != kOK )
  throw( std::runtime_error( "an error occurred in MILPSolver::compute()" ) );

 // if required, write the problem to file- - - - - - - - - - - - - - - - - -
 if( ! output_file.empty() )
  SCIP_CALL_ABORT( SCIPwriteOrigProblem( scip , output_file.c_str() , NULL ,
					 FALSE ) );

 // the actual call to SCIP - - - - - - - - - - - - - - - - - - - - - - - - -
 
 if( int_vars > 0 ) {  // the MIP case- - - - - - - - - - - - - - - - - - - -

 // SCIP constraint handler for adding cuts and lazy constraints

 if( ( CutSepPar & 7 ) ||
   ( UpCutOff < Inf< double >() ) || ( LwCutOff > Inf< double >() ) ) {
   if( f_callback_set == false ){
    // the callback has to be set
    SCIP_CALL_ABORT( SCIPincludeObjConshdlr( scip , 
                    new SCIPMILPSolver_Conhdlr( scip , this , CutSepPar ),
                    TRUE ) );
    
    SCIP_CONS* cons;
    SCIP_CALL_ABORT( SCIPcreateSCIPMILPSolver_basiccb(
        scip,                /**< SCIP data structure */
        &cons,               /**< pointer to hold the created constraint */
        "callback",          /**< name of constraint */
        vars                 /**< active set of SCIP variables */
      ) );

    SCIP_CALL( SCIPaddCons(scip, cons) );
    SCIP_CALL( SCIPreleaseCons(scip, &cons) );
    
    f_callback_set = true;
    }
    //else 
      // it was already set, nothing to do
   }
 else {
   if( f_callback_set ) {    // the callback was set
    f_callback_set = false;
    }
  }
 }

 auto st =  SCIPsolve( scip ) ;

 // decode SCIP exit status - - - - - - - - - - - - - - - - - - - - - - - - -
 switch( SCIPgetStatus( scip ) ) {
  case( SCIP_STATUS_OPTIMAL ):
  case( SCIP_STATUS_GAPLIMIT ):
  case( SCIP_STATUS_SOLLIMIT ):   sol_status = kOK; break;
  case( SCIP_STATUS_INFEASIBLE ): sol_status = kInfeasible; break;
  case( SCIP_STATUS_NODELIMIT ):  sol_status = kStopIter;   break;
  case( SCIP_STATUS_TIMELIMIT ):  sol_status = kStopTime;   break;
  case( SCIP_STATUS_INFORUNBD ):
  case( SCIP_STATUS_UNBOUNDED ):  sol_status = kUnbounded;  break;
  default:                        sol_status = kError + SCIPgetStatus( scip );
  }

 unlock();  // unlock the mutex
 return( sol_status );
 }

/*--------------------------------------------------------------------------*/

Solver::OFValue SCIPMILPSolver::get_lb( void )
{
 OFValue lower_bound = 0;

 switch( SCIPgetObjsense( scip ) ) {
  case( SCIP_OBJSENSE_MINIMIZE ):
   switch( sol_status ) {
    case( kUnbounded ):  lower_bound = -Inf< OFValue >(); break;
    case( kInfeasible ): lower_bound = Inf< OFValue >();  break;
    default:          lower_bound = SCIPgetDualbound( scip ) + constant_value;
    }
   break;
  case( SCIP_OBJSENSE_MAXIMIZE ):
   switch( sol_status ) {
    case( kUnbounded ):  lower_bound = Inf< OFValue >();  break;
    case( kInfeasible ): lower_bound = -Inf< OFValue >(); break;
    default:          lower_bound = SCIPgetPrimalbound( scip ) + constant_value;
    }
   break;
  default:
   throw( std::runtime_error( "Objective type not yet defined" ) );
  }

 return( lower_bound );
 }

/*--------------------------------------------------------------------------*/

Solver::OFValue SCIPMILPSolver::get_ub( void )
{
 OFValue upper_bound = 0;

 switch( SCIPgetObjsense( scip ) ) {
  case( SCIP_OBJSENSE_MINIMIZE ):
   switch( sol_status ) {
    case( kUnbounded ):  upper_bound = -Inf< OFValue >(); break;
    case( kInfeasible ): upper_bound = Inf< OFValue >();  break;
    default:          upper_bound = SCIPgetPrimalbound( scip ) + constant_value;
    }
   break;
  case( SCIP_OBJSENSE_MAXIMIZE ):
   switch( sol_status ) {
    case( kUnbounded ):  upper_bound = Inf< OFValue >();  break;
    case( kInfeasible ): upper_bound = -Inf< OFValue >(); break;
    default:          upper_bound = SCIPgetDualbound( scip ) + constant_value;
    }
   break;
  default:
   throw( std::runtime_error( "Objective type not yet defined" ) );
  }

 return( upper_bound );
 }

/*--------------------------------------------------------------------------*/

bool SCIPMILPSolver::has_var_solution( void )
{
 switch( sol_status ) {
  case( kOK ):
  case( kUnbounded ): return( true );
  }

 return( false );
 }

/*--------------------------------------------------------------------------*/

bool SCIPMILPSolver::is_var_feasible( void )
{
 return( sol_status != kInfeasible );
 }

/*--------------------------------------------------------------------------*/

Solver::OFValue SCIPMILPSolver::get_var_value( void )
{
 switch( SCIPgetObjsense( scip ) ) {
  case( SCIP_OBJSENSE_MINIMIZE ): return( get_ub() );
  case( SCIP_OBJSENSE_MAXIMIZE ): return( get_lb() );
  default:
   throw( std::runtime_error( "Objective type not yet defined" ) );
  }
 }

/*--------------------------------------------------------------------------*/

void SCIPMILPSolver::get_var_solution( Configuration * solc )
{
 SCIP_SOL * sol = SCIPgetBestSol( scip );
 if( ! sol )
  return;

 std::vector< double > x( numcols );
 SCIP_RETCODE status = SCIPgetSolVals( scip , sol , numcols ,
				       vars.data() , x.data() );
 if( status != SCIP_OKAY )
  throw( std::runtime_error( "Unable to get the SCIP solution values" ) );

 get_var_solution( x );
 } 

/*--------------------------------------------------------------------------*/

void SCIPMILPSolver::get_var_solution( const std::vector< double > & x )
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

bool SCIPMILPSolver::has_dual_solution( void )
{
 SCIP_Bool has_dual_solution;

 has_dual_solution = SCIPisDualSolAvailable( scip , FALSE );
 
 if( has_dual_solution )
   return( true );
 else
   return( false );
 }

/*--------------------------------------------------------------------------*/

bool SCIPMILPSolver::is_dual_feasible( void )
{
 return( false );  // TODO
 }

/*--------------------------------------------------------------------------*/

void SCIPMILPSolver::get_dual_solution( Configuration * solc )
{
 std::vector< SCIP_Real > pi( numrows , 0 );
 std::vector< SCIP_Real > dj( numcols , 0 );

 if( numrows > 0 )
  for( int i = 0 ; i < numrows ; ++i )
   SCIP_CALL_ABORT( SCIPgetDualSolVal( scip , cons[ i ] , & pi[ i ] , NULL) );

 for( int j = 0 ; j < numrows ; ++j )
  dj[ j ] = SCIPgetVarRedcost( scip , vars[ j ] );
  
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
	       "SCIPMILPSolver::get_dual_solution: invalid dual value." ) );

  if( throw_reduced_cost_exception ) {
   if( var_is_fixed && ( ! lhs_con ) && ( var_lb != 0 ) ) {
    /* The Variable is fixed but it has no associated OneVarConstraint
     * with both bounds equal to the value of the Variable. */

    throw( std::logic_error(
     "SCIPMILPSolver::get_dual_solution: variable with index " +
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
                "SCIPMILPSolver::get_dual_solution: variable with index " +
		std::to_string( i ) + " has no OneVarConstraint." ) );
     }
   }
  }
 }  // end( SCIPMILPSolver::get_dual_solution )

/*--------------------------------------------------------------------------*/

bool SCIPMILPSolver::has_dual_direction( void )
{
return(false);
 /*std::vector< double > y( numrows , 0 );

 for( int i = 0 ; i < numrows ; ++i )
   SCIP_CALL_ABORT( SCIPgetColFarkasCoef ( scip , )

 return( ! bool( CPXdualfarkas( env , lp , y.data() , & proof ) ) );
 }*/
 }
/*--------------------------------------------------------------------------*/

void SCIPMILPSolver::get_dual_direction( Configuration * dirc )
{
 SCIPABORT();  // TODO
 }

/*--------------------------------------------------------------------------*/

void SCIPMILPSolver::write_lp( const std::string & filename )
{
 SCIP_CALL_ABORT( SCIPwriteOrigProblem( scip , filename.c_str() ,
                                        "mps" , FALSE ) );
 }

/*--------------------------------------------------------------------------*/

int SCIPMILPSolver::get_nodes( void ) const
{
 return( ( int ) SCIPgetNTotalNodes( scip ) );
 }

/*--------------------------------------------------------------------------*/
/*-------------------- PROTECTED FIELDS OF THE CLASS -----------------------*/
/*--------------------------------------------------------------------------*/

void SCIPMILPSolver::var_modification( const VariableMod * mod )
{
 // call the method of MILPSolver to update dictionaries (and nothing else)
 MILPSolver::var_modification( mod );

 if( SCIPisTransformed( scip ) )
  SCIP_CALL_ABORT( SCIPfreeTransform( scip ) );

 auto var = static_cast< ColVariable * >( mod->variable() );
 int idx = index_of_variable( var );
 if( idx == Inf< int >() )  // the variable is no longer there (?)
  return;

 // read old variable type
 SCIP_VARTYPE oldtype = SCIPvarGetType( vars[ idx ] );
 SCIP_Bool infeas = FALSE;

  // react to changes in the integrality - - - - - - - - - - - - - - - - - - -
 if( ColVariable::is_integer( mod->old_state() ) !=
     ColVariable::is_integer( mod->new_state() ) ) {

  // construct new variable type
  SCIP_VARTYPE new_ctype;
  if( var->is_integer() && ( ! relax_int_vars ) ) {
   if( var->is_unitary() && var->is_positive() )
    new_ctype = SCIP_VARTYPE_BINARY; // Binary
   else
    new_ctype = SCIP_VARTYPE_INTEGER; // Integer
   }
  else
   new_ctype = SCIP_VARTYPE_CONTINUOUS;  // Continuous

  SCIP_CALL_ABORT( SCIPchgVarType( scip , vars[ idx ] ,
				                          new_ctype , & infeas ) );
  assert( ! infeas );
  }   // end( if( new integrality != old integrality ) )

 // react to fix / unfix- - - - - - - - - - - - - - - - - - - - - - - - - - -
 if( var->is_fixed() ) {  // fix the variable
  SCIP_Bool fixed = 0;
  SCIP_CALL_ABORT( SCIPfixVar( scip , vars[ idx ] , var->get_value() ,
			       & infeas , & fixed ) );
  assert( fixed );
  }
 else{
  auto bd = SCIPMILPSolver::get_problem_bounds( *var );

  if( SCIPvarGetLbOriginal( vars[ idx ] ) != bd[ 0 ] )
    SCIP_CALL_ABORT( SCIPchgVarLb( scip , vars[ idx ] , bd[ 0 ] ) );

  if( SCIPvarGetUbOriginal( vars[ idx ] ) != bd[ 1 ] )
    SCIP_CALL_ABORT( SCIPchgVarUb( scip , vars[ idx ] , bd[ 1 ] ) );
  }

 }  // end( SCIPMILPSolver::var_modification )

/*--------------------------------------------------------------------------*/

void SCIPMILPSolver::objective_modification( const ObjectiveMod * mod )
{
 // call the method of MILPSolver to update sense direction
 MILPSolver::objective_modification( mod );

 if( SCIPisTransformed( scip ) )
  SCIP_CALL_ABORT( SCIPfreeTransform( scip ) );

 /* ObjectiveMod class does not include any modification types except
  * for eSetMin and eSetMax.
  * To change OF coefficients, a FunctionMod must be used. */

 switch( mod->type() ) {
  case( ObjectiveMod::eSetMin ):
   SCIP_CALL_ABORT( SCIPsetObjsense( scip , SCIP_OBJSENSE_MINIMIZE ) );
   break;

  case( ObjectiveMod::eSetMax ):
   SCIP_CALL_ABORT( SCIPsetObjsense( scip , SCIP_OBJSENSE_MAXIMIZE ) );
   break;

  default:
   throw( std::invalid_argument( "Invalid type of ObjectiveMod" ) );
  }
 }

/*--------------------------------------------------------------------------*/

void SCIPMILPSolver::const_modification( const ConstraintMod * mod )
{
 // no point in calling the method of MILPSolver, as it does nothing
 // MILPSolver::const_modification( mod );

 if( SCIPisTransformed( scip ) )
  SCIP_CALL_ABORT( SCIPfreeTransform( scip ) );

 auto con = dynamic_cast< const FRowConstraint * >(mod->constraint());
 if( ! con )  // this should not happen
  return;     // but in case, nothing to do

 RowConstraint::RHSValue con_lhs = NAN;
 RowConstraint::RHSValue con_rhs = NAN;

 auto idx = index_of_constraint( con );
 if( idx == Inf< int >() )  // the Constraint is no longer there (?)
  return;                   // nothing to do
 
 SCIP_CONS * scip_con = cons[ idx ];

 switch( mod->type() ) {
  case( ConstraintMod::eRelaxConst ):
   // In order to relax the constraint all we do is transform it
   // into an inequality with RHS equal to infinity
   SCIP_CALL_ABORT( SCIPchgLhsLinear( scip , scip_con ,
				      -SCIPinfinity( scip ) ) );
   SCIP_CALL_ABORT( SCIPchgRhsLinear( scip , scip_con ,
				      SCIPinfinity( scip ) ) );
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

   con_lhs = con->get_lhs() == -Inf< double >() ?
             -SCIPinfinity( scip ) : con->get_lhs();
   con_rhs = con->get_rhs() == Inf< double >() ?
             SCIPinfinity( scip ) : con->get_rhs();

   SCIP_CALL_ABORT( SCIPchgLhsLinear( scip , scip_con , con_lhs ) );
   SCIP_CALL_ABORT( SCIPchgRhsLinear( scip , scip_con , con_rhs ) );
   break;

  default:
   throw( std::invalid_argument( "Invalid type of ObjectiveMod" ) );
  }
 }  // end( SCIPMILPSolver::const_modification )

/*--------------------------------------------------------------------------*/

void SCIPMILPSolver::bound_modification( const OneVarConstraintMod * mod )
{
 // no point in calling the method of MILPSolver, as it does nothing
 // MILPSolver::bound_modification( mod );

 /* The same ColVariable can have more active OneVarConstraints,
  * so each time we modify one of them we have to check if LHS and RHS
  * of the Variable change. */

 if( SCIPisTransformed( scip ) )
  SCIP_CALL_ABORT( SCIPfreeTransform( scip ) );

 auto con = static_cast< OneVarConstraint * >( mod->constraint() );
 auto var = static_cast< ColVariable * >( con->get_active_var( 0 ) );
 if( ! var )  // this should never happen
  return;     // but in case, there is nothing to do

 auto idx = index_of_variable( var );
 if( idx == Inf< int >() )  // the ColVariable has been removed
  return;                   // is strange, but there is nothing to do

 auto scip_var = vars[ idx ];

 switch( mod->type() ) {
  case( RowConstraintMod::eChgLHS ): {
   SCIP_Real lb = SCIPMILPSolver::get_problem_lb( *var );
   SCIP_VARTYPE type = SCIPvarGetType	(	scip_var );
   SCIP_Bool inf;
   if( type == SCIP_VARTYPE_BINARY && lb != 0 ){
    /* During presolving, an integer variable whose bound changes to {0,1} 
    *  is upgraded by SCIP to a binary variable. Thus, here we assume that 
    *  we are changing bound for a previously declared integer variable,
    *  which SCIP automatically changed into binary. */
    SCIP_VARTYPE newtype = SCIP_VARTYPE_INTEGER;
    SCIP_Bool inf;
    SCIP_CALL_ABORT( SCIPchgVarType( scip , scip_var , 
                        newtype , &inf ) );
    }

   SCIP_CALL_ABORT( SCIPchgVarLb( scip , scip_var , lb ) );
   break;
   }

  case( RowConstraintMod::eChgRHS ): {
   SCIP_Real ub = SCIPMILPSolver::get_problem_ub( *var );
   SCIP_VARTYPE type = SCIPvarGetType	(	scip_var );
   if( type == SCIP_VARTYPE_BINARY && ub != 1 ){
    /* During presolving, an integer variable whose bound changes to {0,1} 
    *  is upgraded by SCIP to a binary variable. Thus, here we assume that 
    *  we are changing bound for a previously declared integer variable,
    *  which SCIP automatically changed into binary. */
    SCIP_VARTYPE newtype = SCIP_VARTYPE_INTEGER;
    SCIP_Bool inf;
    SCIP_CALL_ABORT( SCIPchgVarType( scip , scip_var , 
                      newtype , &inf ) );
    }

   SCIP_CALL_ABORT( SCIPchgVarUb( scip , scip_var , ub ) );
   break;
   }

  case( RowConstraintMod::eChgBTS ): {
   auto bd = SCIPMILPSolver::get_problem_bounds( *var );
   SCIP_VARTYPE type = SCIPvarGetType	(	scip_var );
   if( type == SCIP_VARTYPE_BINARY && ( bd[ 0 ] != 0 || bd[ 1 ] != 1 ) ){
    /* During presolving, an integer variable whose bound changes to {0,1} 
    *  is upgraded by SCIP to a binary variable. Thus, here we assume that 
    *  we are changing bound for a previously declared integer variable,
    *  which SCIP automatically changed into binary. */
    SCIP_VARTYPE newtype = SCIP_VARTYPE_INTEGER;
    SCIP_Bool inf;
    SCIP_CALL_ABORT( SCIPchgVarType( scip , scip_var , 
                      newtype , &inf ) );
    }

   SCIP_CALL_ABORT( SCIPchgVarLb( scip , scip_var , bd[ 0 ] ) );
   SCIP_CALL_ABORT( SCIPchgVarUb( scip , scip_var , bd[ 1 ] ) );
   break;
   }

  default:
   throw( std::invalid_argument( "Invalid type of OneVarConstraintMod" ) );
  }
 }  // end( SCIPMILPSolver::bound_modification )

/*--------------------------------------------------------------------------*/

void SCIPMILPSolver::objective_function_modification( const FunctionMod * mod )
{
 // no point in calling the method of MILPSolver, as it does nothing
 // MILPSolver::objective_function_modification( mod );

 auto f = mod->function();

 if( SCIPisTransformed( scip ) )
  SCIP_CALL_ABORT( SCIPfreeTransform( scip ) );

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

   auto idxit = idxs.begin();
   auto & cp = lf->get_v_var();

   for( auto v :  modl->vars() )
    if( auto idx = *(idxit++) ; idx < Inf< Index >() ) {
     auto nval = cp[ idx ].second;
     auto vidx = index_of_variable( static_cast< const ColVariable * >( v ) );
     SCIP_CALL_ABORT( SCIPchgVarObj( scip , vars[ vidx ] , nval ) );
     }

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

   auto idxit = idxs.begin();
   auto & cp = qf->get_v_var();

   for( auto v :  modl->vars() )
    if( auto idx = *(idxit++) ; idx < Inf< Index >() ) {
     auto nval = std::get< 1 >( cp[ idx ] );
     auto vidx = index_of_variable( static_cast< const ColVariable * >( v ) );
     SCIP_CALL_ABORT( SCIPchgVarObj( scip , vars[ vidx ] , nval ) );
     }
     
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
  c_Vec_p_Var * milpvars;
  if( auto modlr = dynamic_cast< const C05FunctionModRngd * >( modl ) ) {
   idxs = qf->map_index( modlr->vars() , modlr->range() );
   milpvars = & modlr->vars();
   }
  else
   if( auto modls = dynamic_cast< const C05FunctionModSbst * >( modl ) ) {
    idxs = qf->map_index( modls->vars() , modls->subset() );
    milpvars = & modls->vars();
    }
   else
    throw( std::logic_error( "unknown type of C05FunctionModLinRngd" ) );

   auto idxit = idxs.begin();
   auto & cp = qf->get_v_var();

   for( auto v : *milpvars )
    if( auto idx = *(idxit++) ; idx < Inf< Index >() ) {
      auto nlinval = std::get< 1 >( cp[ idx ] );
      auto nquadval = std::get< 2 >( cp[ idx ] );
      auto vidx = index_of_variable( static_cast< const ColVariable * >( v ) );
      // change linear coefficient 
      SCIP_CALL_ABORT( SCIPchgVarObj( scip , vars[ vidx ] , nlinval ) );

      // Check if v already had a quad coeff associated
      if( aux_vars[ vidx ] != NULL ){
        
        if( nquadval != 0 )
          /* Variable v already had a quadratic coefficient and the new one is 
          *  nonzero. Thus, it is sufficient to change his linear coeff and the
          *  quad ones associated with aux var */

          // change quadratic coefficient by changing the lin coeff of aux var
          SCIP_CALL_ABORT( SCIPchgVarObj( scip , aux_vars[ vidx ] , nquadval ) );
        else{
          /* Variable v already had a quadratic coefficient and the new one is 
          *  zero. Thus, we have to remove the aux var and the aux_con created */
          SCIP_VAR * scip_aux_var = aux_vars[ vidx ];
          aux_vars[ vidx ] = nullptr;
          SCIP_CONS * scip_aux_con = aux_cons[ vidx ]; 
          aux_cons[ vidx ] = nullptr;
          SCIP_Bool deleted = 0;

          // Remove aux constraint
          SCIP_CALL_ABORT( SCIPdelCons (	scip, scip_aux_con ) );	

          // Remove aux variable
          SCIP_CALL_ABORT( SCIPdelVar( scip , scip_aux_var , & deleted ) );
          assert( deleted );
        }
      }
      else{
        /* Variable v didn't have a quadratic coefficient associated. Thus, it is 
        *  necessary to create a new aux_var and aux_con to handle it. 
        *  For more information, see SCIPMILPSolver::compute() */

        if( nquadval == 0 )
          continue; // Nothing to do

        SCIP_Real z_lb = - SCIPinfinity( scip );
        SCIP_Real z_ub = SCIPinfinity( scip ) ;

        SCIP_VARTYPE z_type = SCIP_VARTYPE_CONTINUOUS;

        SCIP_VAR * z = nullptr;
        SCIP_CALL_ABORT( SCIPcreateVarBasic( scip , & z , nullptr , 
                              z_lb  , z_ub ,  nquadval , 
                              z_type ) );
        SCIP_CALL_ABORT( SCIPaddVar( scip , z ) );
        aux_vars[ vidx ] = z;
        SCIP_CALL_ABORT( SCIPreleaseVar( scip , & z ) );

        // Add auxiliary constraints z - xˆ2 >= 0
        SCIP_Real con_lhs = 0;
        SCIP_Real con_rhs = SCIPinfinity( scip );
        SCIP_Real lincoef = 1;
        SCIP_Real quadcoef = -1;

        SCIP_CONS * con = nullptr;
        SCIP_VAR * linvar = aux_vars[ vidx ];
        SCIP_VAR * quadvar = vars[ vidx ];
        const char * name = "aux_con";

        #if SCIP_VERSION < 800
          SCIP_CALL_ABORT( SCIPcreateConsBasicQuadratic( scip , & con , name ,
                          1 , & linvar , & lincoef , 1 , & quadvar , & quadvar ,
              & quadcoef , con_lhs , con_rhs ) );
        #else
          SCIP_CALL_ABORT( SCIPcreateConsBasicQuadraticNonlinear( scip , & con ,
              name , 1 , & linvar , & lincoef , 1 , & quadvar ,
              & quadvar , & quadcoef , con_lhs , con_rhs ) );
        #endif

        SCIP_CALL_ABORT( SCIPaddCons( scip , con ) );
        aux_cons[ vidx ] = con;
        SCIP_CALL_ABORT( SCIPreleaseCons( scip , &con ) );
      }
    }
      
  return;
  }

 // Fallback method - Update all costs
 // --------------------------------------------------------------------------
 // reload_objective( f );

 }  // end( SCIPMILPSolver::objective_function_modification )

/*--------------------------------------------------------------------------*/

void SCIPMILPSolver::constraint_function_modification(
						    const FunctionMod * mod )
{
 // no point in calling the method of MILPSolver, as it does nothing
 // MILPSolver::constraint_function_modification( mod );

 if( SCIPisTransformed( scip ) )
  SCIP_CALL_ABORT( SCIPfreeTransform( scip ) );

 auto lf = dynamic_cast< const LinearFunction * >( mod->function() );
 if( ! lf )
  return;

 auto con = dynamic_cast< const FRowConstraint * >( lf->get_Observer() );
 if( ! con )
  return;

 auto cidx = index_of_constraint( con );
 if( cidx == Inf< int >() )  // the Constraint is no longer there (?)
  return;                    // nothing to do
 
 SCIP_CONS * scip_con = cons[ cidx ];

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
   SCIP_Real new_coeff = cp[ idx ].second;
   auto var_idx = index_of_variable( static_cast< const ColVariable * >( v ) );
   SCIP_CALL_ABORT( SCIPchgCoefLinear( scip , scip_con , vars[ var_idx ] ,
				       new_coeff ) );
   }

 // Fallback method - Reload all coefficients
 // --------------------------------------------------------------------------
 /*for( auto el : lf->get_v_var() )
  if( auto idx = index_of_variable( el.first ) ; idx < Inf< int >() )
   SCIP_CALL_ABORT( SCIPchgCoefLinear( scip , scip_con , vars[ idx ] ,
				       el.second ) );*/
  
 }  // end( SCIPMILPSolver::constraint_function_modification )

/*--------------------------------------------------------------------------*/

void SCIPMILPSolver::objective_fvars_modification(
						const FunctionModVars * mod )
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

 auto f = mod->function();
 auto nav = f->get_num_active_var();

 // while changing the coefficients, we have to be careful about the fact
 // that Modification are managed asynchronously with the model changes
 // although the added/removed Variable do exist in the internal data
 // structure of [SCIP]MILPSolver since the Modification are managed
 // strictly in arrival order, they may no longer exist in the model;
 // more to the point, they may no longer be active in the LinearFunction
 
 if( SCIPisTransformed( scip ) )
  SCIP_CALL_ABORT( SCIPfreeTransform( scip ) );

 if( auto lf = dynamic_cast< const LinearFunction * >( f ) ) {
  // Linear objective function

  for( auto v : mod->vars() ) {
   auto var = static_cast< const ColVariable * >( v );
   if( auto idx = index_of_variable( var ) ; idx < Inf< int >() ) {
    SCIP_VAR * scip_var = vars[ idx ];
    if( mod->added() ){
       auto cidx = lf->is_active( var );
       SCIP_Real new_coeff = cidx < nav ? 
                          lf->get_coefficient( cidx ) : 0;
       SCIP_CALL_ABORT( SCIPchgVarObj( scip, scip_var, new_coeff ) );
      }
    else
       SCIP_CALL_ABORT( SCIPchgVarObj( scip, scip_var, 0 ) );
    }
   }

  return;
  }

 if( auto * qf = dynamic_cast< const DQuadFunction * >( f ) ) {
  // Quadratic objective function
  
  for( auto v : mod->vars() ) {
   auto var = static_cast< const ColVariable * >( v );
   if( auto ind = index_of_variable( var ) ; ind < Inf< int >() ) {
    SCIP_VAR * scip_var = vars[ ind ];
    SCIP_Real value = 0;
    SCIP_Real q_value = 0;

    if( mod->added() ){

      /* We have to add new space in the auxiliary vectors. */
      SCIP_VAR * scip_aux_var = nullptr;
      SCIP_CONS * scip_aux_con = nullptr;
      aux_vars.push_back( scip_aux_var );
      aux_cons.push_back( scip_aux_con );
      SCIP_CALL_ABORT( SCIPreleaseVar( scip , & scip_aux_var ) );
      SCIP_CALL_ABORT( SCIPreleaseCons( scip , & scip_aux_con ) );

      if( auto idx = qf->is_active( var ) ; idx < nav ) {
       value = qf->get_linear_coefficient( idx );
       q_value = qf->get_quadratic_coefficient( idx );
       }
      
      // Change linear coeff
      SCIP_CALL_ABORT( SCIPchgVarObj( scip, scip_var, value ) );

      if( q_value != 0 ){
        /* A new variable has been added. Thus, if q_value is non zero, 
        *  it is necessary to create a new aux_var and aux_con to handle it. 
        *  For more information, see SCIPMILPSolver::compute() */

        SCIP_Real z_lb = - SCIPinfinity( scip );
        SCIP_Real z_ub = SCIPinfinity( scip ) ;

        SCIP_VARTYPE z_type = SCIP_VARTYPE_CONTINUOUS;

        SCIP_VAR * z = nullptr;
        SCIP_CALL_ABORT( SCIPcreateVarBasic( scip , & z , nullptr , 
                              z_lb  , z_ub ,  q_value , 
                              z_type ) );
        SCIP_CALL_ABORT( SCIPaddVar( scip , z ) );
        aux_vars[ ind ] = z;
        SCIP_CALL_ABORT( SCIPreleaseVar( scip , & z ) );

        // Add auxiliary constraints z - xˆ2 >= 0
        SCIP_Real con_lhs = 0;
        SCIP_Real con_rhs = SCIPinfinity( scip );
        SCIP_Real lincoef = 1;
        SCIP_Real quadcoef = -1;

        SCIP_CONS * con = nullptr;
        SCIP_VAR * linvar = aux_vars[ ind ];
        const char * name = "aux_con";

        #if SCIP_VERSION < 800
          SCIP_CALL_ABORT( SCIPcreateConsBasicQuadratic( scip , & con , name ,
                          1 , & linvar , & lincoef , 1 , & scip_var , & scip_var ,
              & quadcoef , con_lhs , con_rhs ) );
        #else
          SCIP_CALL_ABORT( SCIPcreateConsBasicQuadraticNonlinear( scip , & con ,
              name , 1 , & linvar , & lincoef , 1 , & scip_var ,
              & scip_var , & quadcoef , con_lhs , con_rhs ) );
        #endif

        SCIP_CALL_ABORT( SCIPaddCons( scip , con ) );
        aux_cons[ ind ] = con;
        SCIP_CALL_ABORT( SCIPreleaseCons( scip , &con ) );
       }
     }
    else{
      SCIP_CALL_ABORT( SCIPchgVarObj( scip, scip_var, 0 ) );

      /* We have to remove the aux var and the aux_con associated
      *  to scip_var */
      SCIP_VAR * scip_aux_var = aux_vars[ ind ];
      aux_vars[ ind ] = nullptr;
      SCIP_CONS * scip_aux_con = aux_cons[ ind ]; 
      aux_cons[ ind ] = nullptr;
      SCIP_Bool deleted = 0;

      // Remove aux constraint
      SCIP_CALL_ABORT( SCIPdelCons ( scip, scip_aux_con ) );	

      // Remove aux variable
      SCIP_CALL_ABORT( SCIPdelVar( scip , scip_aux_var , & deleted ) );
      assert( deleted );
      }
    }
  }
 }
 
 // This should never happen
 throw( std::invalid_argument( "Unknown type of Objective Function" ) );

 }  // end( SCIPMILPSolver::objective_fvars_modification )

/*--------------------------------------------------------------------------*/

void SCIPMILPSolver::constraint_fvars_modification(
						const FunctionModVars * mod )
{
 // this is called in response to Variable being added to / removed from the
 // Constraint; however, note that all Variable are supposed to exist at the
 // time this is called, so index_of_variable() is always correct

 // no point in calling the method of MILPSolver, as it does nothing
 // MILPSolver::constraint_fvars_modification( mod );

 auto * lf = dynamic_cast< const LinearFunction * >( mod->function() );
 if( ! lf )
  return;

 // check the modification type
 if( ( ! dynamic_cast< const C05FunctionModVarsAddd * >( mod ) ) &&
     ( ! dynamic_cast< const C05FunctionModVarsRngd * >( mod ) ) &&
     ( ! dynamic_cast< const C05FunctionModVarsSbst * >( mod ) ) )
  throw( std::invalid_argument( "This type of FunctionModVars is not handled"
				) );

 auto con = dynamic_cast< FRowConstraint * >( lf->get_Observer() );
 if( ! con )
  return;

 auto cidx = index_of_constraint( con );
 if( cidx == Inf< int >() )  // the Constraint is no longer there (?)
  return;                   // nothing to do

 SCIP_CONS * scip_con = cons[ cidx ];

 if( SCIPisTransformed( scip ) )
  SCIP_CALL_ABORT( SCIPfreeTransform( scip ) );

 auto nav = lf->get_num_active_var();
 for( auto v : mod->vars() ){
  auto var = static_cast< const ColVariable * >( v );
  if( auto vidx = index_of_variable( var ) ; vidx < Inf< int >() ) {
    SCIP_VAR * scip_var = vars[ vidx ];
    if( mod->added() ){
     auto idx = lf->is_active( var );
     SCIP_CALL_ABORT( SCIPaddCoefLinear( scip , scip_con , scip_var ,
					 idx < nav ? lf->get_coefficient( idx ) : 0 ) );
     }
    else
     SCIP_CALL_ABORT( SCIPdelCoefLinear( scip , scip_con , scip_var ) );
   }
  }
 }  // end( SCIPMILPSolver::constraint_fvars_modification )

/*----------------------------------------------------------------------------

void SCIPMILPSolver::dynamic_modification( const BlockModAD * mod )
{
 MILPSolver::dynamic_modification( mod );
 }

----------------------------------------------------------------------------*/

void SCIPMILPSolver::add_dynamic_constraint( const FRowConstraint * con )
{
 // call the method of MILPSolver to update the dictionaries (only)
 MILPSolver::add_dynamic_constraint( con );

 if( auto f = con->get_function() ) {
  auto lf = dynamic_cast< const LinearFunction * >( f );
  if( ! lf )
   throw( std::invalid_argument( "The Constraint is not linear" ) );

  if( SCIPisTransformed( scip ) )
   SCIP_CALL_ABORT( SCIPfreeTransform( scip ) );

  SCIP_CONS * scip_con = nullptr;

  SCIP_Real con_lhs = con->get_lhs() == -Inf< double >() ?
                      -SCIPinfinity( scip ) : con->get_lhs();
  SCIP_Real con_rhs = con->get_rhs() == Inf< double >() ?
                      SCIPinfinity( scip ) : con->get_rhs();

  char name[ 32 ];
  std::snprintf( name , sizeof( name ) , "%p" , ( void * ) con );
  SCIP_CALL_ABORT( SCIPcreateConsBasicLinear( scip , & scip_con , name , 0 ,
                                              nullptr , nullptr ,
                                              con_lhs , con_rhs ) );

  // get the coefficients to fill the matrix
  for( auto & el : lf->get_v_var() )
   if( auto idx = index_of_variable( el.first ) ; idx < Inf< int >() )
    SCIP_CALL_ABORT( SCIPaddCoefLinear( scip , scip_con , vars[ idx ] ,
                                        el.second ) );

  SCIP_CALL_ABORT( SCIPaddCons( scip , scip_con ) );
  cons.push_back( scip_con );
  SCIP_CALL_ABORT( SCIPreleaseCons( scip , & scip_con ) );
  }
 }  // end( SCIPMILPSolver::add_dynamic_constraint )

/*--------------------------------------------------------------------------*/

void SCIPMILPSolver::add_dynamic_variable( const ColVariable * var )
{
 // call the method of MILPSolver to update the dictionaries (only)
 MILPSolver::add_dynamic_variable( var );

 if( SCIPisTransformed( scip ) )
  SCIP_CALL_ABORT( SCIPfreeTransform( scip ) );

 auto bd = SCIPMILPSolver::get_problem_bounds( *var );
 SCIP_VARTYPE vartype = SCIP_VARTYPE_BINARY;

 // variable type
 if( var->is_integer() && ( ! relax_int_vars ) ) {
  if( var->is_unitary() && var->is_positive() )
   vartype = SCIP_VARTYPE_BINARY;
  else
   vartype = SCIP_VARTYPE_INTEGER;
  }
 else
  vartype = SCIP_VARTYPE_CONTINUOUS;

 SCIP_VAR * scip_var = nullptr;

 SCIP_CALL_ABORT( SCIPcreateVarBasic( scip , & scip_var , nullptr ,
                                      bd[ 0 ] , bd[ 1 ] , 0.0 , vartype ) );
 vars.push_back( scip_var );
 SCIP_CALL_ABORT( SCIPaddVar( scip, scip_var ) );

 SCIP_CALL_ABORT( SCIPreleaseVar( scip , & scip_var ) );

 }  // end( SCIPMILPSolver::add_dynamic_variable )

/*--------------------------------------------------------------------------*/

void SCIPMILPSolver::add_dynamic_bound( const OneVarConstraint * con )
{
 // no point in calling the method of MILPSolver, as it does nothing
 // MILPSolver::add_dynamic_bound( con );

 if( SCIPisTransformed( scip ) )
  SCIP_CALL_ABORT( SCIPfreeTransform( scip ) );

 auto var = static_cast< ColVariable * >( con->get_active_var( 0 ) );
 if( ! var )
  throw( std::logic_error( "SCIPMILPSolver: added a bound on no Variable" ) );

 auto idx = index_of_variable( var );
 if( idx == Inf< int >() )
  throw( std::logic_error( "SCIPMILPSolver: added a bound on unknown Variable"
			   ) );

 SCIP_VAR * scip_var = vars[ idx ];

 auto bd = SCIPMILPSolver::get_problem_bounds( *var );
 SCIP_CALL_ABORT( SCIPchgVarLb( scip , scip_var , bd[ 0 ] ) );
 SCIP_CALL_ABORT( SCIPchgVarUb( scip , scip_var , bd[ 1 ] ) );
 }

/*--------------------------------------------------------------------------*/

void SCIPMILPSolver::remove_dynamic_constraint( const FRowConstraint * con )
{
 if( SCIPisTransformed( scip ) )
  SCIP_CALL_ABORT( SCIPfreeTransform( scip ) );

 int index = index_of_dynamic_constraint( con );
 if( index < Inf< int >() ) {
  auto it = cons.begin() + index;
  SCIP_CONS * scip_con = *it;
  cons.erase( it );
  SCIP_CALL_ABORT( SCIPdelCons( scip , scip_con ) );

  // call the method of MILPSolver to update the dictionaries (only)
  MILPSolver::remove_dynamic_constraint( con );
  }
 else
  throw( std::runtime_error( "Dynamic constraint not found" ) );
 }

/*--------------------------------------------------------------------------*/

void SCIPMILPSolver::remove_dynamic_variable( const ColVariable * var )
{
 if( SCIPisTransformed( scip ) )
  SCIP_CALL_ABORT( SCIPfreeTransform( scip ) );

 int index = index_of_dynamic_variable( var );
 if( index < Inf< int >() ) {
  auto it = vars.begin() + index;
  SCIP_VAR * scip_var = *it;
  vars.erase( it );

  SCIP_Bool deleted = 0;
  SCIP_CALL_ABORT( SCIPdelVar( scip , scip_var , & deleted ) );
  assert( deleted );

  // call the method of MILPSolver to update the dictionaries (only)
  MILPSolver::remove_dynamic_variable( var );
  }
 else
  throw( std::runtime_error( "Dynamic variable not found" ) );
 }

/*--------------------------------------------------------------------------*/

void SCIPMILPSolver::remove_dynamic_bound( const OneVarConstraint * con )
{
 // no point in calling the method of MILPSolver, as it does nothing
 // MILPSolver::remove_dynamic_bound( con );

 if( SCIPisTransformed( scip ) )
  SCIP_CALL_ABORT( SCIPfreeTransform( scip ) );

 // note: this only works because remove_dynamic_constraint[s]() do *not*
 //       clear the removed OneVarConstraint, and therefore we can easily
 //       reconstruct which ColVariable it was about
 auto var = static_cast< ColVariable * >( con->get_active_var( 0 ) );
 if( ! var )  // this should never happen
  return;     // but in case, there is nothing to do

 int idx = index_of_variable( var );
 if( idx == Inf< int >() )  // the ColVariable has been removed
  return;                   // is strange, but there is nothing to do

 SCIP_VAR * scip_var = vars[ idx ];

 auto bd = SCIPMILPSolver::get_problem_bounds( *var );
 SCIP_CALL_ABORT( SCIPchgVarLb( scip , scip_var , bd[ 0 ] ) );
 SCIP_CALL_ABORT( SCIPchgVarUb( scip , scip_var , bd[ 1 ] ) );
 }

/*--------------------------------------------------------------------------*/

void SCIPMILPSolver::perform_separation( Configuration * cfg ,
					std::vector< int > & rmatbeg ,
					std::vector< int > & rmatind ,
					std::vector< double > & rmatval ,
					std::vector< double > & rhs ,
               std::vector< double > & lhs )
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

    rhs.push_back( con_rhs );
    lhs.push_back( con_lhs );

    rmatbeg.push_back( rmatind.size() );
    }
   }  // end( for each added FRowConstraint )
  }  // end( main loop )
 }  // end( SCIPMILPSolver::perform_separation )

/*--------------------------------------------------------------------------*/

void SCIPMILPSolver::set_f_cb_mutex( ){
 
 // we should be in critical section where different SCIP threads may compete
 // for access to the Block: ensure mutual exclusion
 f_callback_mutex.lock();

 // ensure no interference from other threads (except SCIP ones) by also
 // lock()-ing the Block
 bool owned = f_Block->is_owned_by( f_id );
 if( ( ! owned ) && ( ! f_Block->lock( f_id ) ) )
   throw( std::runtime_error( "Unable to lock the Block" ) );
}

/*--------------------------------------------------------------------------*/

void SCIPMILPSolver::unset_f_cb_mutex( ){
 
 bool owned = f_Block->is_owned_by( f_id );
 if( ! owned )
   f_Block->unlock( f_id );  // unlock the Block

 // critical section ends here, release the mutex
 f_callback_mutex.unlock();
}

/*--------------------------------------------------------------------------*/

std::vector< SCIP_VAR * > SCIPMILPSolver::get_SCIP_var( void ){
  // return the SCIP variable stored in the protected field of the class
  return vars;
}

/*--------------------------------------------------------------------------*/
/*------------------- METHODS FOR HANDLING THE PARAMETERS ------------------*/
/*--------------------------------------------------------------------------*/

void SCIPMILPSolver::set_par( idx_type par, int value )
{
 // Solver parameters explicitly mapped in SCIP
 switch( par ) {
  case( intThrowReducedCostException ):
   throw_reduced_cost_exception = bool( value );
   return;
  case( intCutSepPar ): 
   CutSepPar = value; 
   return;
  case( intMaxIter ):
   SCIP_CALL_ABORT( SCIPsetLongintParam( scip, "limits/nodes", value ) );
   return;
  case( intMaxSol ):
   SCIP_CALL_ABORT( SCIPsetIntParam( scip, "limits/solutions", value ) );
   return;
  case( intLogVerb ):
   SCIP_CALL_ABORT( SCIPsetIntParam( scip, "display/verblevel", value ) );
   return;
  }

 // SCIP parameters
 if( par >= intFirstSCIPPar && par < intLastAlgParSCPS ) {
  const std::string & scip_par =
   SMSpp_to_SCIP_int_pars[ par - intFirstSCIPPar ];

  // Bool, int and long SCIP parameters are handled as SMS++ int parameters
  SCIP_PARAM * param = SCIPgetParam( scip , scip_par.c_str() );
  SCIP_PARAMTYPE type = SCIPparamGetType( param );

  if( type == SCIP_PARAMTYPE_BOOL )
   SCIP_CALL_ABORT( SCIPsetBoolParam( scip , scip_par.c_str() , value ) );
  else
   if( type == SCIP_PARAMTYPE_INT )
    SCIP_CALL_ABORT( SCIPsetIntParam( scip , scip_par.c_str() , value ) );
   else
    if( type == SCIP_PARAMTYPE_LONGINT )
     SCIP_CALL_ABORT( SCIPsetLongintParam( scip , scip_par.c_str() , value ) );

  return;
  }

 MILPSolver::set_par( par, value );
 }

/*--------------------------------------------------------------------------*/

void SCIPMILPSolver::set_par( idx_type par , double value )
{
 // Solver parameters explicitly mapped in SCIP
 switch( par ) {
  case( dblMaxTime ):
   SCIP_CALL_ABORT( SCIPsetRealParam( scip , "limits/time" , value ) );
   return;
   // case( dblRelAcc ): // TODO
   //  return;
   // case( dblAbsAcc ): // TODO
   //  return;
  case( dblUpCutOff ): UpCutOff = value; return;
  case( dblLwCutOff ): LwCutOff = value; return;
  case( dblRAccSol ):
   SCIP_CALL_ABORT( SCIPsetRealParam( scip , "limits/gap" , value ) );
   return;
  case( dblAAccSol ):
   SCIP_CALL_ABORT( SCIPsetRealParam( scip , "limits/absgap" , value ) );
   return;
  case( dblFAccSol ):
   SCIP_CALL_ABORT( SCIPsetRealParam( scip , "numerics/feastol" , value ) );
   return;
  }

 // SCIP parameters
 if( par >= dblFirstSCIPPar && par < dblLastAlgParSCPS ) {
  const std::string & scip_par =
   SMSpp_to_SCIP_dbl_pars[ par - dblFirstSCIPPar ];
  SCIP_CALL_ABORT( SCIPsetRealParam( scip , scip_par.c_str() , value ) );
  return;
  }

 MILPSolver::set_par( par , value );
 }

/*--------------------------------------------------------------------------*/

void SCIPMILPSolver::set_par( idx_type par , std::string && value )
{
 // SCIP parameters
 if( par >= strFirstSCIPPar && par < strLastAlgParSCPS ) {
  const std::string & scip_par =
   SMSpp_to_SCIP_str_pars[ par - strFirstSCIPPar ];

  // char and string SCIP parameters are handled as SMS++ string parameters
  SCIP_PARAM * param = SCIPgetParam( scip , scip_par.c_str() );
  SCIP_PARAMTYPE type = SCIPparamGetType( param );

  if( type == SCIP_PARAMTYPE_CHAR )
   SCIP_CALL_ABORT( SCIPsetCharParam( scip , scip_par.c_str() ,
				      value[ 0 ] ) );
  else
   if( type == SCIP_PARAMTYPE_STRING )
    SCIP_CALL_ABORT( SCIPsetStringParam( scip , scip_par.c_str() ,
					 value.c_str() ) );
  return;
  }

 MILPSolver::set_par( par , std::move( value ) );
 }

/*--------------------------------------------------------------------------*/

void SCIPMILPSolver::set_par( idx_type par , std::vector< int > && value )
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

void SCIPMILPSolver::set_par( idx_type par ,
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

ThinComputeInterface::idx_type SCIPMILPSolver::get_num_int_par( void ) const
{
 return( MILPSolver::get_num_int_par() + intLastAlgParSCPS
	                               - intLastAlgParMILP );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

ThinComputeInterface::idx_type SCIPMILPSolver::get_num_str_par( void ) const
{
 return( MILPSolver::get_num_str_par() + strLastAlgParSCPS
	                               - strLastAlgParMILP );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

ThinComputeInterface::idx_type SCIPMILPSolver::get_num_dbl_par( void ) const
{
 return( MILPSolver::get_num_dbl_par() + dblLastAlgParSCPS
	                               - dblLastAlgParMILP );
 }

/*--------------------------------------------------------------------------*/

Solver::idx_type SCIPMILPSolver::get_num_vint_par( void ) const {
 return( MILPSolver::get_num_vint_par()
	 + vintLastAlgParCPXS - vintLastAlgParMILP );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

Solver::idx_type SCIPMILPSolver::get_num_vstr_par( void ) const {
 return( MILPSolver::get_num_vstr_par()
	 + vstrLastAlgParCPXS - vstrLastAlgParMILP );
 }

/*--------------------------------------------------------------------------*/

int SCIPMILPSolver::get_int_par( idx_type par ) const
{
 int value;
 SCIP_Bool bool_val;
 SCIP_Longint long_val;

 // solver parameters explicitly mapped in SCIP
 switch( par ) {
  case( intThrowReducedCostException ):
   return( throw_reduced_cost_exception );
  case( intCutSepPar ):
   return( CutSepPar );
  case( intMaxIter ):
   SCIP_CALL_ABORT( SCIPgetLongintParam( scip , "limits/nodes" , & long_val
					 ) );
   return( ( int ) long_val );
  case( intMaxSol ):
   SCIP_CALL_ABORT( SCIPgetIntParam( scip , "limits/solutions" , & value ) );
   return( value );
  case( intLogVerb ):
   SCIP_CALL_ABORT( SCIPgetIntParam( scip , "display/verblevel" , & value ) );
   return( value );
  }

 // SCIP parameters
 if( ( par >= intFirstSCIPPar ) && ( par < intLastAlgParSCPS ) ) {
  const std::string & scip_par =
   SMSpp_to_SCIP_int_pars[ par - intFirstSCIPPar ];

  // Bool, int and long SCIP parameters are handled as SMS++ int parameters
  SCIP_PARAM * param = SCIPgetParam( scip , scip_par.c_str() );
  SCIP_PARAMTYPE type = SCIPparamGetType( param );

  switch( type ) {
   case( SCIP_PARAMTYPE_BOOL ):
    SCIP_CALL_ABORT( SCIPgetBoolParam( scip , scip_par.c_str() , & bool_val
				       ) );
    return( ( int ) bool_val );
   case( SCIP_PARAMTYPE_INT ):
    SCIP_CALL_ABORT( SCIPgetIntParam( scip , scip_par.c_str() , & value ) );
    return( value );
   case( SCIP_PARAMTYPE_LONGINT ):
    SCIP_CALL_ABORT( SCIPgetLongintParam( scip , scip_par.c_str() ,
					  & long_val ) );
    return( ( int ) long_val );
   default:;  // here just to avoid a pesky warning
   }
  }

 return( MILPSolver::get_int_par( par ) );
 }

/*--------------------------------------------------------------------------*/

double SCIPMILPSolver::get_dbl_par( idx_type par ) const
{
 double value;

 // sxolver parameters explicitly mapped in SCIP
 switch( par ) {
  case( dblMaxTime ):
   SCIP_CALL_ABORT( SCIPgetRealParam( scip , "limits/time" , & value ) );
   return( value );
   // case( dblRelAcc ):   // TODO
   //  return( 1e-6 );
   // case( dblAbsAcc ):   // TODO
   //  return( Inf< OFValue >() );
  case( dblUpCutOff ): return( UpCutOff );
  case( dblLwCutOff ): return( LwCutOff );
  case( dblRAccSol ):
   SCIP_CALL_ABORT( SCIPgetRealParam( scip , "limits/gap" , & value ) );
   return( value );
  case( dblAAccSol ):
   SCIP_CALL_ABORT( SCIPgetRealParam( scip , "limits/absgap" , & value ) );
   return( value );
  case( dblFAccSol ):
   SCIP_CALL_ABORT( SCIPgetRealParam( scip , "numerics/feastol" , & value ) );
   return( value );
  }

 // SCIP parameters
 if( par >= dblFirstSCIPPar && par < dblLastAlgParSCPS ) {
  const std::string & scip_par =
   SMSpp_to_SCIP_dbl_pars[ par - dblFirstSCIPPar ];
  SCIP_CALL_ABORT( SCIPgetRealParam( scip, scip_par.c_str(), &value ) );
  return( value );
 }

 return( MILPSolver::get_dbl_par( par ) );
 }

/*--------------------------------------------------------------------------*/

const std::string & SCIPMILPSolver::get_str_par( idx_type par ) const
{
 static std::string value;
 char * str_val;

 // SCIP parameters
 if( ( par >= strFirstSCIPPar ) && ( par < strLastAlgParSCPS ) ) {
  const auto & scip_par = SMSpp_to_SCIP_str_pars[ par - strFirstSCIPPar ];

  // char and string SCIP parameters are handled as SMS++ string parameters
  SCIP_PARAM * param = SCIPgetParam( scip, scip_par.c_str() );
  SCIP_PARAMTYPE type = SCIPparamGetType( param );

  switch( type ) {
   case( SCIP_PARAMTYPE_CHAR ):
    value.resize( 1 );
    SCIP_CALL_ABORT( SCIPgetCharParam( scip , scip_par.c_str() ,
				       value.data() ) );
    return( value );
   case( SCIP_PARAMTYPE_STRING ):
    SCIP_CALL_ABORT( SCIPgetStringParam( scip , scip_par.c_str() ,
					 & str_val ) );
    value = str_val;
    return( value );
   default:;  // here just to avoid a pesky warning
   }
  }

 return( MILPSolver::get_str_par( par ) );
 }

/*--------------------------------------------------------------------------*/

const std::vector< int > & SCIPMILPSolver::get_vint_par( idx_type par ) const
{
 if( par == vintCutSepCfgInd )
  return( CutSepCfgInd );

 return( MILPSolver::get_vint_par( par ) );
 }

/*--------------------------------------------------------------------------*/

const std::vector< std::string > & SCIPMILPSolver::get_vstr_par( idx_type par )
 const
{
 if( par == vstrConfigDBFName )
  return( ConfigDBFName );

 return( MILPSolver::get_vstr_par( par ) );
 }

/*--------------------------------------------------------------------------*/

int SCIPMILPSolver::get_dflt_int_par( idx_type par ) const
{
 if( ( par == intThrowReducedCostException ) || ( par == intCutSepPar ) )
  return( 0 );

 int value;
 SCIP_Longint long_val;
 SCIP_PARAM * param;

 // Solver parameters explicitly mapped in SCIP
 switch( par ) {
  case( intMaxIter ):
   param = SCIPgetParam( scip , "limits/nodes" );
   return( ( int ) SCIPparamGetLongintDefault( param ) );
  case( intMaxSol ):
   param = SCIPgetParam( scip , "limits/solutions" );
   return( SCIPparamGetIntDefault( param ) );
  case( intLogVerb ):
   param = SCIPgetParam( scip, "display/verblevel" );
   return( SCIPparamGetIntDefault( param ) );
  }

 // SCIP parameters
 if( ( par >= intFirstSCIPPar ) && ( par < intLastAlgParSCPS ) ) {
  const std::string & scip_par =
   SMSpp_to_SCIP_int_pars[ par - intFirstSCIPPar ];

  // Bool, int and long SCIP parameters are handled as SMS++ int parameters
  param = SCIPgetParam( scip , scip_par.c_str() );

  switch( SCIPparamGetType( param ) ) {
   case( SCIP_PARAMTYPE_BOOL ):
    return( ( int ) SCIPparamGetBoolDefault( param ) );
   case( SCIP_PARAMTYPE_INT ):
    return( SCIPparamGetIntDefault( param ) );
   case( SCIP_PARAMTYPE_LONGINT ):
    return( ( int ) SCIPparamGetLongintDefault( param ) );
   default:;  // here just to avoid a pesky warning
   }
  }

 return( MILPSolver::get_dflt_int_par( par ) );
 }

/*--------------------------------------------------------------------------*/

double SCIPMILPSolver::get_dflt_dbl_par( idx_type par ) const
{
 double value;
 SCIP_PARAM * param;

 switch( par ) {
  case( dblMaxTime ):
   param = SCIPgetParam( scip , "limits/time" );
   return( SCIPparamGetRealDefault( param ) );
   // case( dblRelAcc ):   // TODO
   //  return( 1e-6 );
   // case( dblAbsAcc ):   // TODO
   //  return( Inf< OFValue >() );
  case( dblUpCutOff ): return( Inf< double >() );
  case( dblLwCutOff ): return( -Inf< double >() );
  case( dblRAccSol ):
   param = SCIPgetParam( scip , "limits/gap" );
   return( SCIPparamGetRealDefault( param ) );
  case( dblAAccSol ):
   param = SCIPgetParam( scip , "limits/absgap" );
   return( SCIPparamGetRealDefault( param ) );
  case( dblFAccSol ):
   param = SCIPgetParam( scip , "numerics/feastol" );
   return( SCIPparamGetRealDefault( param ) );
  }

 // SCIP parameters
 if( ( par >= dblFirstSCIPPar ) && ( par < dblLastAlgParSCPS ) ) {
  const std::string & scip_par =
   SMSpp_to_SCIP_dbl_pars[ par - dblFirstSCIPPar ];

  param = SCIPgetParam( scip , scip_par.c_str() );
  //?? SCIP_PARAMTYPE type = SCIPparamGetType( param );
  return( SCIPparamGetRealDefault( param ) );
  }

 return( MILPSolver::get_dflt_dbl_par( par ) );
 }

/*--------------------------------------------------------------------------*/

const std::string & SCIPMILPSolver::get_dflt_str_par( idx_type par ) const
{
 static std::string value;
 char * str_val;

 // SCIP parameters
 if( ( par >= strFirstSCIPPar ) && ( par < strLastAlgParSCPS ) ) {
  const std::string & scip_par =
   SMSpp_to_SCIP_str_pars[ par - strFirstSCIPPar ];

  // char and string SCIP parameters are handled as SMS++ string parameters
  SCIP_PARAM * param = SCIPgetParam( scip , scip_par.c_str() );

  switch( SCIPparamGetType( param ) ) {
   case( SCIP_PARAMTYPE_CHAR ):
    value.resize( 1 );
    value[ 0 ] = SCIPparamGetCharDefault( param );
    return( value );
   case( SCIP_PARAMTYPE_STRING ):
    str_val = SCIPparamGetStringDefault( param );
    value = str_val;
    return( value );
   default:;  // here just to avoid a pesky warning
   }
  }

 return( MILPSolver::get_dflt_str_par( par ) );
 }

/*--------------------------------------------------------------------------*/

const std::vector< int > & SCIPMILPSolver::get_dflt_vint_par( idx_type par )
 const
{
 static std::vector< int > _empty;
 if( par == vintCutSepCfgInd )
  return( _empty );

 return( MILPSolver::get_dflt_vint_par( par ) );
 }

/*--------------------------------------------------------------------------*/

const std::vector< std::string > & SCIPMILPSolver::get_dflt_vstr_par(
							idx_type par ) const
{
 static std::vector< std::string > _empty;
 if( par == vstrConfigDBFName )
  return( _empty );

 return( MILPSolver::get_dflt_vstr_par( par ) );
 }

/*--------------------------------------------------------------------------*/

Solver::idx_type SCIPMILPSolver::int_par_str2idx( const std::string & name )
 const
{
 if( name == "intThrowReducedCostException" )
  return( intThrowReducedCostException );

 if( name == "intCutSepPar" )
  return( intCutSepPar );

 // SCIP parameters
 auto it = find_if( SCIP_to_SMSpp_int_pars.begin(),
                    SCIP_to_SMSpp_int_pars.end(),
                    [ & ]( const auto & pair ) {
                     return( pair.first == name );
                     } );

 if( ( it != SCIP_to_SMSpp_int_pars.end() ) && ( it->first == name ) )
  return( it->second );

 return( MILPSolver::int_par_str2idx( name ) );
 }

/*--------------------------------------------------------------------------*/

const std::string & SCIPMILPSolver::int_par_idx2str( idx_type idx ) const
{
 static const std::array< std::string , 2 > _pars =
                     { "intThrowReducedCostException" , "intCutSepPar" };
 if( idx == intThrowReducedCostException )
  return( _pars[ 0 ] );

 if( idx == intCutSepPar )
  return( _pars[ 1 ] );

 // SCIP parameters
 if( ( idx >= intFirstSCIPPar ) && ( idx < intLastAlgParSCPS ) )
  return( SMSpp_to_SCIP_int_pars[ idx - intFirstSCIPPar ] );


 return( MILPSolver::int_par_idx2str( idx ) );
 }

/*--------------------------------------------------------------------------*/

Solver::idx_type SCIPMILPSolver::dbl_par_str2idx( const std::string & name )
 const
{
 // SCIP parameters
 auto it = find_if( SCIP_to_SMSpp_dbl_pars.begin() ,
                    SCIP_to_SMSpp_dbl_pars.end() ,
                    [ & ]( const auto & pair ) {
                     return( pair.first == name );
                     } );

 if( ( it != SCIP_to_SMSpp_dbl_pars.end() ) && ( it->first == name ) )
  return( it->second );

 return( MILPSolver::dbl_par_str2idx( name ) );
 }

/*--------------------------------------------------------------------------*/

const std::string & SCIPMILPSolver::dbl_par_idx2str( idx_type idx ) const
{
 // SCIP parameters
 if( ( idx >= dblFirstSCIPPar ) && ( idx < dblLastAlgParSCPS ) )
  return( SMSpp_to_SCIP_dbl_pars[ idx - dblFirstSCIPPar ] );

 return( MILPSolver::dbl_par_idx2str( idx ) );
 }

/*--------------------------------------------------------------------------*/

Solver::idx_type SCIPMILPSolver::str_par_str2idx( const std::string & name )
 const
{
 // SCIP parameters
 auto it = find_if( SCIP_to_SMSpp_str_pars.begin() ,
                    SCIP_to_SMSpp_str_pars.end() ,
                    [ & ]( const auto & pair ) {
		     return( pair.first == name );
                     } );

 if( ( it != SCIP_to_SMSpp_str_pars.end() ) && ( it->first == name ) )
  return( it->second );

 return( MILPSolver::str_par_str2idx( name ) );
}

/*--------------------------------------------------------------------------*/

const std::string & SCIPMILPSolver::str_par_idx2str( idx_type idx ) const
{
 // SCIP parameters
 if( ( idx >= strFirstSCIPPar ) && ( idx < strLastAlgParSCPS ) )
  return( SMSpp_to_SCIP_str_pars[ idx - strFirstSCIPPar ] );

 return( MILPSolver::str_par_idx2str( idx ) );
 }

/*--------------------------------------------------------------------------*/

Solver::idx_type SCIPMILPSolver::vint_par_str2idx(
					     const std::string & name ) const
{
 if( name == "vintCutSepCfgInd" )
  return( vintCutSepCfgInd );

 return( MILPSolver::vint_par_str2idx( name ) );
 }

/*--------------------------------------------------------------------------*/

const std::string & SCIPMILPSolver::vint_par_idx2str( idx_type idx ) const
{
 static const std::string _pars = "vintCutSepCfgInd";
 if( idx == vintCutSepCfgInd )
  return( _pars );

 return( MILPSolver::vint_par_idx2str( idx ) );
 }

/*--------------------------------------------------------------------------*/

Solver::idx_type SCIPMILPSolver::vstr_par_str2idx(
					     const std::string & name ) const
{
 if( name == "vstrConfigDBFName" )
  return( vstrConfigDBFName );

 return( MILPSolver::vstr_par_str2idx( name ) );
 }

/*--------------------------------------------------------------------------*/

const std::string & SCIPMILPSolver::vstr_par_idx2str( idx_type idx ) const
{
 static const std::string _pars = "vstrConfigDBFName";
 if( idx == vstrConfigDBFName )
  return( _pars );

 return( MILPSolver::vstr_par_idx2str( idx ) );
 }

/*--------------------------------------------------------------------------*/

#ifdef MILPSOLVER_DEBUG

bool NotIsNull( SCIP_VAR * v ) { return( ! ( v == NULL) ); }

void SCIPMILPSolver::check_status( void )
{
 auto n_aux_var = std::count_if( aux_vars.begin() , aux_vars.end() , NotIsNull );
 auto tmp = SCIPgetNVars( scip ) - n_aux_var;
 if( numcols != tmp )
  DEBUG_LOG( "numcols is " << numcols << " but SCIPgetNVars() returns "
	     << tmp << std::endl );

 tmp = SCIPgetNConss( scip ) - n_aux_var;
 if( numrows != tmp )
  DEBUG_LOG( "numrows is " << numrows << " but SCIPgetNConss() returns "
	     << tmp << std::endl );

 tmp = SCIPgetNBinVars( scip ) + SCIPgetNIntVars( scip );
 if( int_vars != tmp )
  DEBUG_LOG( "int_vars is " << int_vars << " but SCIP has actually "
	     << tmp << " integer variables" << std::endl );

 MILPSolver::check_status();
 }

#endif

/*--------------------------------------------------------------------------*/
/*-------------------- PRIVATE METHODS OF THE CLASS ------------------------*/
/*--------------------------------------------------------------------------*/

Configuration * SCIPMILPSolver::get_cfg( Index ci ) const
{
 if( ci >= CutSepCfgInd.size() )
  return( nullptr );
 auto dbi = CutSepCfgInd[ ci ];
 if( ( dbi < 0 ) || ( Index( dbi ) >= v_ConfigDB.size() ) )
  return( nullptr );
 return( v_ConfigDB[ dbi ] );
 }

/*--------------------------------------------------------------------------*/
/*--------- End Methods related to Class SCIPMILPSolver.cpp ----------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*------ Begin Methods related to Class SCIPMILPSolver_Conhdlr.cpp ---------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*---------------------------- DEFINITIONS ---------------------------------*/
/*--------------------------------------------------------------------------*/

/* fundamental constraint handler properties */
#define CONSHDLR_NAME          "SCIPMILPSolver_Conhdlr"
#define CONSHDLR_DESC          "Constraint handler for SCIPMILPSolver. This \
particular class will provide callback method inside a SCIPMILPSolver object, \
allowing it to generate user cuts or lazy constraints. "
#define CONSHDLR_ENFOPRIORITY  -2000000  // priority of the constraint handler for
                                    // constraint enforcing
#define CONSHDLR_CHECKPRIORITY -2000000  // priority of the constraint handler for 
                                    // checking feasibility
#define CONSHDLR_EAGERFREQ     -1   // frequency for using all instead of only 
                                    // the useful constraints in separation,
                                    //propagation and enforcement
#define CONSHDLR_NEEDSCONS     TRUE // should the constraint handler be skipped, 
                                    // if no constraints are available? 
#define CONSHDLR_SEPAPRIORITY  -2000000  // priority of the constraint handler for 
                                     // separation
#define CONSHDLR_SEPAFREQ      1.0   // frequency for separating cuts; zero 
                                     // means to separate only in the root node 
#define CONSHDLR_DELAYSEPA     FALSE // should separation method be delayed, if 
                                     // other separators found cuts?
#define CONSHDLR_PROPFREQ      -1    // frequency for propagating domains; 
                                     // zero means only preprocessing propagation
#define CONSHDLR_DELAYPROP     FALSE // should propagation method be delayed, 
                                     // if other propagators found reductions?
#define CONSHDLR_PROP_TIMING   SCIP_PROPTIMING_BEFORELP 
#define CONSHDLR_PRESOLTIMING  SCIP_PRESOLTIMING_MEDIUM 
#define CONSHDLR_MAXPREROUNDS  -1    // maximal number of presolving rounds the
                                     // constraint handler participates in 
                                     // (-1: no limit) */

/*--------------------------------------------------------------------------*/
/*------------------------------ DATA STRUCTURES ---------------------------*/
/*--------------------------------------------------------------------------*/

/** constraint data for lazy constraints */
struct SCIP_ConsData
{
   SCIP_Bool new_cut; // Parameter used to tell if the constraint handler found 
                      // with the CHECK function a possible cut
};

/*--------------------------------------------------------------------------*/
/*------------------------- CONSTRUCTOR AND DESTRUCTOR ---------------------*/
/*--------------------------------------------------------------------------*/

SCIPMILPSolver_Conhdlr::SCIPMILPSolver_Conhdlr( SCIP* scip,
    SMSpp_di_unipi_it::SCIPMILPSolver* scipmilpsolver,
     unsigned char SeparationPar
    ) : ObjConshdlr( scip , CONSHDLR_NAME , CONSHDLR_DESC, CONSHDLR_SEPAPRIORITY,
         CONSHDLR_ENFOPRIORITY, CONSHDLR_CHECKPRIORITY, CONSHDLR_SEPAFREQ, 
         CONSHDLR_PROPFREQ, CONSHDLR_EAGERFREQ , CONSHDLR_MAXPREROUNDS,
         CONSHDLR_DELAYSEPA, CONSHDLR_DELAYPROP, CONSHDLR_NEEDSCONS,
         CONSHDLR_PROP_TIMING, CONSHDLR_PRESOLTIMING )
   {
      parent_scipmilpsolver = scipmilpsolver; // set parent scipmilpsolver

      CutSepPar = SeparationPar; // set parameter for deciding if and when 
                                 // separation should be performed
   }

 SCIPMILPSolver_Conhdlr::~SCIPMILPSolver_Conhdlr(){
 }

/*--------------------------------------------------------------------------*/
/*--------------------------------- LOCAL METHODS --------------------------*/
/*--------------------------------------------------------------------------*/

/** local method used to perform the separation within the constraint handler.
*
*  The method is called from the verify_separation or scipmilpsolver_separation 
*  function, in order to create the structures to add a lazy constraint/user cut. 
*/
 static SCIP_RETCODE perform_separation( 
   SCIP*                 scip,                         /* SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,                     /* the constraint handler 
                                                         itself */
   SMSpp_di_unipi_it::SCIPMILPSolver* scipmilpsolver,  /* the milpsolver model from
                                                        * which the conhdlr has been 
                                                        * called */
   SCIP_CONS**           conss,                        /* array of constraints 
                                                        * to process */
   SCIP_SOL*             sol,                          /* primal solution that should 
                                                        * be separated */
   SCIP_Bool             enforce,                      /* whether we are in 
                                                        *enforcing */
   SCIP_RESULT*          result,                        /* pointer to store the result 
                                                        * of the separation call */
	std::vector< int >   & rmatbeg ,
	std::vector< int >   & rmatind ,
	std::vector< double > & rmatval ,
	std::vector< double > & rhs ,
 std::vector< double > & lhs
   )
 {
 std::vector< SCIP_VAR * > scip_vars;
 int nvars;

 // get variables data from the scipmilpsolver 
 /** NOTE: This is essential, because using SCIPMILPSolver we are adding auxiliary 
 * variables which don't have to be considered in the computation of the best 
 * solution. */
 scip_vars =  scipmilpsolver->get_SCIP_var();
 nvars = scip_vars.size();

 // this is a critical section where different SCIP threads may compete
 // for access to the Block: ensure mutual exclusion
 scipmilpsolver->set_f_cb_mutex();

 // get the feasible solution
 std::vector< double > x( nvars );
 for( int i = 0 ; i < nvars ; ++i ){
   SCIP_VAR* var = scip_vars[ i ];
   SCIP_Real best_sol = SCIPgetSolVal( scip , sol, var );
   x[i] = best_sol;
 }

 // write it in the Variable of the Block
 scipmilpsolver->get_var_solution( x );

 // now perform the lazy constraint separation with the right Configuration

 if( enforce == FALSE ){
  // Adding a user cut
  int depth;
  depth = SCIPgetSubscipDepth( scip ); // find the depth of the current node

  scipmilpsolver->perform_separation( 
   scipmilpsolver->get_cfg( depth ? 1 : 0 ) ,
   rmatbeg , rmatind , rmatval , rhs , lhs );
  }
 else{
  // Adding a lazy constraint
  scipmilpsolver->perform_separation( 
   scipmilpsolver->get_cfg( 2 ) ,
   rmatbeg , rmatind , rmatval , rhs , lhs );
  }

 // critical section ends here, release the mutex
 scipmilpsolver->unset_f_cb_mutex();

 return SCIP_OKAY;
}

/** local method used to check if a user cut/lazy constraint exists
*
*  The method is called from the scip_check function, in order to state if a 
*  lazy constraint/user cut is available and can be added. 
*/
 static SCIP_RETCODE verify_separation( 
   SCIP*                 scip,                         /* SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,                     /* the constraint handler 
                                                         itself */
   SMSpp_di_unipi_it::SCIPMILPSolver* scipmilpsolver,  /* the milpsolver model from
                                                        * which the conhdlr has been 
                                                        * called */
   SCIP_CONS**           conss,                        /* array of constraints 
                                                        * to process */
   SCIP_SOL*             sol,                          /* primal solution that should 
                                                        * be separated */
   SCIP_Bool             enforce,                      /* whether we are in 
                                                        *enforcing */
   SCIP_RESULT*          result                        /* pointer to store the result 
                                                        * of the separation call */
   )
 {
 assert(result != NULL);

 // get all required structures
 SCIP_CONSDATA* consdata;

 consdata = SCIPconsGetData(conss[0]);
 assert(consdata != NULL);

 std::vector< SCIP_VAR * > scip_vars;
 int nvars;

 // get variables data from the scipmilpsolver 
 /** NOTE: This is essential, because using SCIPMILPSolver we are adding auxiliary 
 * variables which don't have to be considered in the computation of the best 
 * solution. */
 scip_vars =  scipmilpsolver->get_SCIP_var();
 nvars = scip_vars.size();

 // now perform the lazy constraint separation with the right Configuration
 std::vector< int > rmatbeg;
 std::vector< int > rmatind;
 std::vector< double > rmatval;
 std::vector< double > rhs;
 std::vector< double > lhs;

 perform_separation( scip , conshdlr , scipmilpsolver , conss , sol ,
                     enforce , result , rmatbeg , rmatind , rmatval ,
                     rhs , lhs );

 // Understand if any user cut/ lazy constraint are available
 if( ! rmatbeg.empty() ){
    // at least acutting plane has been found
    *result = SCIP_INFEASIBLE;
    consdata->new_cut = TRUE; // inform constraint data that 
                              // a new cut has been found
  }
 else{
    // the solution is already feasible
    *result = SCIP_FEASIBLE;
    consdata->new_cut = TRUE; // inform constraint data that 
                              // a no cut is available
 }

 return SCIP_OKAY;
}


/** local method used to add a lazy constraint or a user cut
*
*  The method is called from a ENFOLP/SEPALP, respectively to add a new
*  lazy constraint or a user cut to the model. The type of constraint is
*  specified with the 'enforce' parameter.
*/
 static
 SCIP_RETCODE scipmilpsolver_separation( 
   SCIP*                 scip,                         /* SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,                     /* the constraint handler 
                                                         itself */
   SMSpp_di_unipi_it::SCIPMILPSolver* scipmilpsolver,  /* the milpsolver model from
                                                        * which the conhdlr has been 
                                                        * called */
   SCIP_CONS**           conss,                        /* array of constraints 
                                                        * to process */
   SCIP_SOL*             sol,                          /* primal solution that should 
                                                        * be separated */
   SCIP_Bool             enforce,                      /* whether we are in 
                                                        *enforcing */
   SCIP_RESULT*          result                        /* pointer to store the result 
                                                        * of the separation call */
   )
 {
 assert(result != NULL);

 // get all required structures
 SCIP_CONSDATA* consdata;

 consdata = SCIPconsGetData(conss[0]);
 assert(consdata != NULL);

 /* if a new cut is available, the constraint data must be already informed */
 if( !consdata->new_cut )
  // strange, but nothing to do
  return SCIP_OKAY;

 std::vector< SCIP_VAR * > scip_vars;
 int nvars;

 // get variables data from the scipmilpsolver 
 /** NOTE: This is essential, because using SCIPMILPSolver we are adding auxiliary 
 * variables which don't have to be considered in the computation of the best 
 * solution. */
 //scip_vars = consdata->vars;
 scip_vars =  scipmilpsolver->get_SCIP_var();
 nvars = scip_vars.size();

 // now perform the lazy constraint separation with the right Configuration
 std::vector< int > rmatbeg;
 std::vector< int > rmatind;
 std::vector< double > rmatval;
 std::vector< double > rhs;
 std::vector< double > lhs;

 perform_separation( scip , conshdlr , scipmilpsolver , conss , sol ,
                     enforce , result , rmatbeg , rmatind , rmatval ,
                     rhs , lhs );

 // if any lazy constraint/user cut was generated, add them
 if( ! rmatbeg.empty() )
   for( int c = 0 ; c < rhs.size() ; ++c ){
      SCIP_ROW* row;
      SCIP_CALL( SCIPcreateEmptyRowConshdlr( scip, &row, conshdlr, 
                     "scipmilpsolver_cut", lhs[ c ] , rhs[ c ], 
                     FALSE, FALSE, TRUE) );

      SCIP_CALL( SCIPcacheRowExtensions(scip, row) );

      int nnz; // number of nonzero coefficients in the actual lazy constraint
      int beg_idx = rmatbeg[ c ]; // idx from where new coefficients begin 
         
      // retrieve number of nonzero
      if( c < rhs.size() - 1)
         nnz = rmatbeg[ c + 1 ] - rmatbeg[ c ]; 
      else
         nnz = rmatind.size() - rmatbeg[ c ];

      for( int counter = 0 ; counter < nnz ; ++counter ){
         int var_idx = rmatind[ beg_idx + counter ];

         SCIP_CALL( SCIPaddVarToRow( scip , row , scip_vars[ var_idx ] , 
                        rmatval[ beg_idx + counter ] ) );
         }

      SCIP_CALL( SCIPflushRowExtensions(scip, row) );
      //SCIP_CALL( SCIPprintRow(scip, row, NULL));

      // Add violated cut. If we are enforcing, then this is enough to add 
      // the cut. Otherwise (we are separating), we check whether the
      // cut is efficacious.
      if( enforce || SCIPisCutEfficacious( scip , sol , row ) ){
         SCIP_Bool infeasible;
         SCIP_CALL( SCIPaddRow( scip , row , FALSE , &infeasible) );
         if ( infeasible )
            *result = SCIP_CUTOFF;
         else
            *result = SCIP_SEPARATED;
         }
      
      SCIP_CALL( SCIPreleaseRow(scip, &row) );
      }

   return SCIP_OKAY;
}

/*--------------------------------------------------------------------------*/
/*------------------- CALLBACK METHODS OF CONSTRAINT HANDLER ---------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/

/** constraint enforcing method of constraint handler for LP solutions */
SCIP_DECL_CONSENFOLP(SCIPMILPSolver_Conhdlr::scip_enfolp)
{  /*lint --e{715}*/

 assert( result != NULL );

 if( ! ( CutSepPar & 4 ) )  // but we don't do lazy constraint separation
   return SCIP_OKAY;         // nothing to do

 SCIP_CALL( scipmilpsolver_separation ( scip, conshdlr, parent_scipmilpsolver,
                conss , NULL , TRUE, result ) );

 return SCIP_OKAY;

}

/*--------------------------------------------------------------------------*/

/** separation method of constraint handler for LP solution */
SCIP_DECL_CONSSEPALP(SCIPMILPSolver_Conhdlr::scip_sepalp)
{
 assert( result != NULL );

 if( ! ( CutSepPar & 3 ) )  // but we don't do user cut separation
   return SCIP_OKAY;        // nothing to do

 int depth;
 depth = SCIPgetSubscipDepth( scip ); // find the depth of the current node

 // if we are at a depth for which separation is not enabled
 if( ( ( ! depth ) && ( ! ( CutSepPar & 1 ) ) ) ||
   ( depth && ( ! ( CutSepPar & 2 ) ) ) )
     return SCIP_OKAY;     // nothing to do

 SCIP_CALL( scipmilpsolver_separation ( scip, conshdlr , parent_scipmilpsolver,
                conss , NULL , FALSE, result) );

 return SCIP_OKAY;
}

/*--------------------------------------------------------------------------*/

/** constraint enforcing method of constraint handler for pseudo solutions */
SCIP_DECL_CONSENFOPS(SCIPMILPSolver_Conhdlr::scip_enfops)
{  /*lint --e{715}*/
   *result = SCIP_DIDNOTRUN;
   return SCIP_OKAY;
}

/** feasibility check method of constraint handler for primal solutions */
SCIP_DECL_CONSCHECK(SCIPMILPSolver_Conhdlr::scip_check)
{  /*lint --e{715}*/

   // Retrieve actual SCIP stage
   SCIP_STAGE current_stage = SCIPgetStage( scip );

   int sep = 0; // parameter to decide if separation is enabled

   sep = verify_separation( scip, conshdlr, parent_scipmilpsolver,
                conss , sol , TRUE, result);
   
   if( sep == 0 ) // No separation is available
    *result = SCIP_FEASIBLE;
   else // we can add some user cut/lazy constraint
    *result = SCIP_INFEASIBLE;

   return SCIP_OKAY;
}

/** variable rounding lock method of constraint handler */
SCIP_DECL_CONSLOCK(SCIPMILPSolver_Conhdlr::scip_lock)
{  /*lint --e{715}*/

 std::vector< SCIP_VAR * > scip_vars;
 int nvars;

 // get variables data from the scipmilpsolver 
 /** NOTE: This is essential, because using SCIPMILPSolver we are adding auxiliary 
 * variables which don't have to be considered in the computation of the best 
 * solution. */
 scip_vars =  parent_scipmilpsolver->get_SCIP_var();
 nvars = scip_vars.size();

 for( int i = 0; i < nvars; i++){
      SCIP_CALL( SCIPaddVarLocksType(scip, scip_vars[i], locktype, nlockspos + nlocksneg, nlockspos + nlocksneg) );
   }

 return SCIP_OKAY;
}

/** transforms constraint data into data belonging to the transformed problem */
SCIP_DECL_CONSTRANS(SCIPMILPSolver_Conhdlr::scip_trans){
   SCIP_CONSDATA* sourcedata;
   SCIP_CONSDATA* targetdata = NULL;

   sourcedata = SCIPconsGetData(sourcecons);
   assert( sourcedata != NULL );

   SCIP_CALL( SCIPallocBlockMemory(scip, &targetdata) );
   targetdata->new_cut = FALSE;

   /* create target constraint */
   SCIP_CALL( SCIPcreateCons(scip, targetcons, SCIPconsGetName(sourcecons), conshdlr, targetdata,
         SCIPconsIsInitial(sourcecons), SCIPconsIsSeparated(sourcecons), SCIPconsIsEnforced(sourcecons),
         SCIPconsIsChecked(sourcecons), SCIPconsIsPropagated(sourcecons),  SCIPconsIsLocal(sourcecons),
         SCIPconsIsModifiable(sourcecons), SCIPconsIsDynamic(sourcecons), SCIPconsIsRemovable(sourcecons),
         SCIPconsIsStickingAtNode(sourcecons)) );

   return SCIP_OKAY;
}

/** frees specific constraint data */
SCIP_DECL_CONSDELETE(SCIPMILPSolver_Conhdlr::scip_delete){  /*lint --e{715}*/
   
   assert(consdata != NULL);
   SCIPfreeBlockMemory(scip, consdata);

   return SCIP_OKAY;
}


/** creates and captures a lazy constraint with all its constraint 
 *  flags set to their default values */
SCIP_RETCODE SMSpp_di_unipi_it::SCIPcreateSCIPMILPSolver_basiccb(
   SCIP*        scip,               /**< SCIP data structure */
   SCIP_CONS**  cons,               /**< pointer to hold the created constraint */
   const char*  name,               /**< name of constraint */
   std::vector< SCIP_VAR * > vars   /**< SCIP vars */
   ){

   SCIP_CALL( SCIPcreateSCIPMILPSolver_cb(scip, cons, name, vars ,
         FALSE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, TRUE) );

   return SCIP_OKAY;
}

/** creates and captures a constraint used which will be used as a separator */
 SCIP_RETCODE SMSpp_di_unipi_it::SCIPcreateSCIPMILPSolver_cb(
   SCIP*        scip,               /**< SCIP data structure */
   SCIP_CONS**  cons,              /**< pointer to hold the created 
                                        constraint */
   const char*  name,               /**< name of constraint */
   std::vector< SCIP_VAR * > vars,   /**< SCIP vars */
   SCIP_Bool    initial,            /**< should the LP relaxation of 
                                        constraint be in the initial LP? */
   SCIP_Bool    separate,           /**< should the constraint be 
                                        separated during LP processing? */
   SCIP_Bool    enforce,            /**< should the constraint be enforced 
                                        during node processing? */
   SCIP_Bool    check,              /**< should the constraint be checked 
                                        for feasibility? */
   SCIP_Bool    propagate,          /**< should the constraint be propagated 
                                        during node processing? */
   SCIP_Bool    local,              /**< is constraint only valid locally? */
   SCIP_Bool    modifiable,         /**< is constraint modifiable (subject 
                                        to column generation)? */
   SCIP_Bool    dynamic,            /**< is constraint dynamic? */
   SCIP_Bool    removable           /**< should the constraint be removed 
                                        from the LP due to aging or cleanup? */
   ){
   
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSDATA* consdata = nullptr;
   int nvars;

   /* find the subtour constraint handler */
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   if( conshdlr == NULL )
   {
      SCIPerrorMessage("scipmilpsolver constraint handler not found\n");
      return SCIP_PLUGINNOTFOUND;
   }

   /* create constraint data */
   SCIP_CALL( SCIPallocBlockMemory(scip, &consdata) ); /*lint !e530*/

   consdata->new_cut = FALSE;

   /* create constraint */
   SCIP_CALL( SCIPcreateCons(scip, cons, name, conshdlr, consdata, initial, 
         separate, enforce, check, propagate, local, modifiable, dynamic, 
         removable, FALSE) );

   return SCIP_OKAY;
}

/*--------------------------------------------------------------------------*/
/*------- End Methods related to Class SCIPMILPSolver_Conhdlr.cpp ----------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*--------------------- End File SCIPMILPSolver.cpp -------------------------*/
/*--------------------------------------------------------------------------*/
