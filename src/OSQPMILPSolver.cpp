/*--------------------------------------------------------------------------*/
/*------------------------- File OSQPMILPSolver.cpp ------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Implementation of the OSQPMILPSolver class, which implements the
 * MILPSolver interface for convex Quadratic Programs with linear
 * constraints using the operator splitting solver OSQP.
 *
 * \author Donato Meoli \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \author Antonio Frangioni \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \copyright &copy; by Donato Meoli, Antonio Frangioni
 */
/*--------------------------------------------------------------------------*/
/*---------------------------- IMPLEMENTATION ------------------------------*/
/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include <algorithm>
#include <array>
#include <limits>
#include <stdexcept>
#include <tuple>

#include "OSQPMILPSolver.h"

/*--------------------------------------------------------------------------*/
/*-------------------------------- USING -----------------------------------*/
/*--------------------------------------------------------------------------*/

using namespace SMSpp_di_unipi_it;

/*--------------------------------------------------------------------------*/
/*-------------------------------- CONSTANTS -------------------------------*/
/*--------------------------------------------------------------------------*/

static const std::array< std::string , 6 > OSQP_int_pars_str = {
 "intOSQPWarmStart" , "intOSQPPolishing" , "intOSQPPolishRefineIter" ,
 "intOSQPScaling" , "intOSQPAdaptiveRho" , "intOSQPCheckTermination" };

static const std::array< std::string , 7 > OSQP_dbl_pars_str = {
 "dblOSQPEpsRel" , "dblOSQPEpsPrimInf" , "dblOSQPEpsDualInf" ,
 "dblOSQPRho" , "dblOSQPSigma" , "dblOSQPAlpha" , "dblOSQPDelta" };

// the settings of OSQP, which give the defaults of the parameters but those
// below

static const OSQPSettings OSQP_dflt = [] {
 OSQPSettings s;
 osqp_set_default_settings( &s );
 return( s );
 }();

// the default of dblFAccSol and of dblOSQPEpsRel, tighter than the 1e-3 of
// OSQP, which leaves errors of that order in the solution and, when the
// polishing fails, in the multipliers

static constexpr double OSQP_dflt_eps = 1e-6;

// the default of intOSQPAdaptiveRho, the adaptation of rho on the decrease
// of the error of the optimality conditions, with which the ADMM reaches
// the tolerances above on degenerate Linear Programs where the adaptation
// every some iterations of OSQP has rho oscillate without converging

static constexpr int OSQP_dflt_adaptive_rho =
 OSQP_ADAPTIVE_RHO_UPDATE_KKT_ERROR;

// the default of dblOSQPEpsPrimInf and of dblOSQPEpsDualInf, tighter than
// the 1e-4 of OSQP, with which a bounded problem having a large solution is
// declared unbounded (or a feasible one infeasible)

static constexpr double OSQP_dflt_eps_inf = 1e-8;

// the status of a solve whose set up, or the update of the data or of the
// settings, has failed

static constexpr int OSQP_setup_failed = -2;

/*--------------------------------------------------------------------------*/
/*------------------------- FACTORY MANAGEMENT -----------------------------*/
/*--------------------------------------------------------------------------*/

SMSpp_insert_in_factory_cpp_0( OSQPMILPSolver );

/*--------------------------------------------------------------------------*/
/*--------------------- CONSTRUCTOR AND DESTRUCTOR -------------------------*/
/*--------------------------------------------------------------------------*/

OSQPMILPSolver::OSQPMILPSolver( void ) : MILPSolver()
{
 d_settings = OSQP_dflt;
 d_settings.verbose = 0;
 d_settings.polishing = 1;
 d_settings.eps_abs = OSQP_dflt_eps;
 d_settings.eps_rel = OSQP_dflt_eps;
 d_settings.adaptive_rho = OSQP_dflt_adaptive_rho;
 d_settings.eps_prim_inf = OSQP_dflt_eps_inf;
 d_settings.eps_dual_inf = OSQP_dflt_eps_inf;
 }

/*--------------------------------------------------------------------------*/

OSQPMILPSolver::~OSQPMILPSolver()
{
 free_solver();
 }

/*--------------------------------------------------------------------------*/
/*--------------------- DERIVED METHODS OF BASE CLASS ----------------------*/
/*--------------------------------------------------------------------------*/

void OSQPMILPSolver::clear_problem( unsigned int what )
{
 MILPSolver::clear_problem( what );

 // the solver of OSQP and the last solution are kept while the problem is
 // loaded again, for the update of the data and the warm start
 if( ! d_reloading ) {
  free_solver();
  d_x.clear();
  d_y.clear();
  }
 d_n = d_m = 0;
 d_row_pos.clear();
 d_box_pos.clear();
 d_loaded = false;
 d_status = -1;
 }

/*--------------------------------------------------------------------------*/

void OSQPMILPSolver::load_problem( void )
{
 d_reloading = true;
 try {
  MILPSolver::load_problem();
  }
 catch( ... ) {
  d_reloading = false;
  throw;
  }
 d_reloading = false;

 if( numquadrows > 0 )
  throw( std::invalid_argument( "OSQPMILPSolver::load_problem: "
				"quadratic constraints are not supported" ) );

 if( ( int_vars > 0 ) && ( relax_int_vars == 0 ) )
  throw( std::invalid_argument( "OSQPMILPSolver::load_problem: integer "
				"variables are only supported if relaxed "
				"(intRelaxIntVars)" ) );

 const OSQPInt n = numcols;
 d_n = n;
 d_sign = ( objsense == -1 ) ? -1 : 1;  // OSQP minimizes

 // the objective: c' x + sum_j q_j x_j^2 + sum_k ndq_k x_i x_j is
 // 1/2 x' P x + q' x with P_jj = 2 q_j and P_ij = P_ji = ndq_k, of which
 // OSQP takes the upper triangle
 d_q.resize( n );
 for( OSQPInt j = 0 ; j < n ; ++j )
  d_q[ j ] = d_sign * objective[ j ];

 // ( column , row , value ) of the upper triangle of P
 std::vector< std::tuple< OSQPInt , OSQPInt , double > > trp;
 for( OSQPInt j = 0 ; j < n ; ++j )
  if( q_objective[ j ] != 0 )
   trp.emplace_back( j , j , d_sign * 2 * q_objective[ j ] );
 for( std::size_t k = 0 ; k < ndq_objective.size() ; ++k ) {
  const OSQPInt i = ndq_rowind[ k ];
  const OSQPInt j = ndq_colind[ k ];
  const double v = d_sign * ndq_objective[ k ];
  if( i == j )
   trp.emplace_back( i , i , 2 * v );
  else
   trp.emplace_back( std::max( i , j ) , std::min( i , j ) , v );
  }
 std::sort( trp.begin() , trp.end() ,
	    []( const auto & a , const auto & b ) {
	     return( std::tie( std::get< 0 >( a ) , std::get< 1 >( a ) ) <
		     std::tie( std::get< 0 >( b ) , std::get< 1 >( b ) ) );
	     } );

 d_Pp.assign( n + 1 , 0 );
 d_Pi.clear();
 d_Px.clear();
 for( std::size_t k = 0 ; k < trp.size() ; ) {
  const auto [ c , r , v0 ] = trp[ k ];
  double v = v0;
  while( ( ++k < trp.size() ) && ( std::get< 0 >( trp[ k ] ) == c ) &&
	 ( std::get< 1 >( trp[ k ] ) == r ) )
   v += std::get< 2 >( trp[ k ] );
  if( v != 0 ) {
   d_Pi.push_back( r );
   d_Px.push_back( v );
   ++d_Pp[ c + 1 ];
   }
  }
 for( OSQPInt j = 0 ; j < n ; ++j )
  d_Pp[ j + 1 ] += d_Pp[ j ];

 // the sides of the rows, with the infinite ones at OSQP_INFTY: those with
 // no finite side are left out, all the others are the first rows of A
 d_row_pos.resize( numrows );
 d_l.clear();
 d_u.clear();
 for( int i = 0 ; i < numrows ; ++i ) {
  double lo , hi;
  switch( sense[ i ] ) {
   case( 'L' ): lo = - Inf< double >();
                hi = rhs[ i ];
                break;
   case( 'G' ): lo = rhs[ i ];
                hi = Inf< double >();
                break;
   case( 'E' ): lo = hi = rhs[ i ];
                break;
   case( 'R' ): if( rngval[ i ] > 0 ) {
                 lo = rhs[ i ];
                 hi = rhs[ i ] + rngval[ i ];
                 }
                else {
                 lo = rhs[ i ] + rngval[ i ];
                 hi = rhs[ i ];
                 }
                break;
   default:
    throw( std::invalid_argument( "OSQPMILPSolver::load_problem: row " +
				  std::to_string( i ) + " has sense " +
				  std::string( 1 , sense[ i ] ) ) );
   }
  lo = std::max( lo , - OSQP_INFTY );
  hi = std::min( hi , OSQP_INFTY );
  if( ( lo <= - OSQP_INFTY ) && ( hi >= OSQP_INFTY ) )
   d_row_pos[ i ] = -1;
  else {
   d_row_pos[ i ] = d_l.size();
   d_l.push_back( lo );
   d_u.push_back( hi );
   }
  }

 // the bounds: a unit row for each variable with a finite one
 d_box_pos.assign( n , -1 );
 for( OSQPInt j = 0 ; j < n ; ++j ) {
  const double lo = std::max( double( lb[ j ] ) , - OSQP_INFTY );
  const double hi = std::min( double( ub[ j ] ) , OSQP_INFTY );
  if( ( lo > - OSQP_INFTY ) || ( hi < OSQP_INFTY ) ) {
   d_box_pos[ j ] = d_l.size();
   d_l.push_back( lo );
   d_u.push_back( hi );
   }
  }
 d_m = d_l.size();

 // A, by columns, from those of MILPSolver
 d_Ap.assign( n + 1 , 0 );
 d_Ai.clear();
 d_Ax.clear();
 std::vector< std::pair< OSQPInt , double > > col;
 for( OSQPInt j = 0 ; j < n ; ++j ) {
  col.clear();
  for( int k = matbeg[ j ] ; k < matbeg[ j ] + matcnt[ j ] ; ++k )
   if( ( d_row_pos[ matind[ k ] ] >= 0 ) && ( matval[ k ] != 0 ) )
    col.emplace_back( d_row_pos[ matind[ k ] ] , matval[ k ] );
  if( d_box_pos[ j ] >= 0 )
   col.emplace_back( d_box_pos[ j ] , 1.0 );
  std::sort( col.begin() , col.end() ,
	     []( const auto & a , const auto & b ) {
	      return( a.first < b.first );
	      } );
  for( std::size_t k = 0 ; k < col.size() ; ++k )
   if( ( ! d_Ai.empty() ) && ( d_Ap[ j + 1 ] > 0 ) &&
       ( d_Ai.back() == col[ k ].first ) )
    d_Ax.back() += col[ k ].second;
   else {
    d_Ai.push_back( col[ k ].first );
    d_Ax.push_back( col[ k ].second );
    ++d_Ap[ j + 1 ];
    }
  d_Ap[ j + 1 ] += d_Ap[ j ];
  }

 d_loaded = true;
 d_status = -1;

 }  // end( OSQPMILPSolver::load_problem )

/*--------------------------------------------------------------------------*/
/*-------------------- METHODS FOR SOLVING THE MODEL -----------------------*/
/*--------------------------------------------------------------------------*/

int OSQPMILPSolver::guts_of_compute( void )
{
 d_iter_limit = false;
 d_fval = 0;
 d_cert.clear();

 // a problem with no variables, or with crossing sides of a row or bounds
 // of a variable, is solved here, OSQP not accepting any of them
 bool feas = true;
 for( OSQPInt i = 0 ; i < d_m ; ++i )
  feas &= ( d_l[ i ] <= d_u[ i ] );
 if( d_n == 0 ) {
  for( OSQPInt i = 0 ; i < d_m ; ++i )
   feas &= ( d_l[ i ] <= 0 ) && ( d_u[ i ] >= 0 );
  if( feas ) {
   d_x.clear();
   d_y.assign( d_m , 0 );
   }
  }
 if( ( d_n == 0 ) || ( ! feas ) ) {
  d_status = feas ? OSQP_SOLVED : OSQP_PRIMAL_INFEASIBLE;
  return( sol_status = decode_status() );
  }

 // the settings that may change between two solves; an infinite intMaxIter
 // is the largest one that OSQP can count to
 d_settings.max_iter = ( d_max_iter < Inf< int >() ) ? d_max_iter
                       : std::numeric_limits< int >::max() - 1;
 d_settings.time_limit = ( d_max_time < Inf< double >() ) ?
  std::max( d_max_time , std::numeric_limits< double >::min() ) :
  OSQP_dflt.time_limit;

 const bool ws = d_warm_start && ( d_x.size() == std::size_t( d_n ) ) &&
                 ( d_y.size() == std::size_t( d_m ) );
 d_settings.warm_starting = ws;

 // the data are given to the solver already set up if the sizes and the
 // sparsity patterns are the same, which is set up again otherwise
 bool setup = d_setup || ( ! d_solver );
 if( d_loaded && ! setup ) {
  OSQPInt m , n;
  osqp_get_dimensions( d_solver , &m , &n );
  setup = ( m != d_m ) || ( n != d_n ) || ( d_Pp != d_sPp ) ||
          ( d_Pi != d_sPi ) || ( d_Ap != d_sAp ) || ( d_Ai != d_sAi );
  if( ! setup ) {
   setup = osqp_update_data_vec( d_solver , d_q.data() , d_l.data() ,
				 d_u.data() ) != 0;
   const bool newP = d_Px != d_sPx;
   const bool newA = d_Ax != d_sAx;
   if( ( ! setup ) && ( newP || newA ) )
    setup = osqp_update_data_mat( d_solver ,
				  newP ? d_Px.data() : OSQP_NULL ,
				  OSQP_NULL , d_Px.size() ,
				  newA ? d_Ax.data() : OSQP_NULL ,
				  OSQP_NULL , d_Ax.size() ) != 0;
   if( ! setup ) {
    d_sPx = d_Px;
    d_sAx = d_Ax;
    }
   }
  }

 if( setup ) {
  free_solver();

  // OSQP takes the arrays of the matrices as non-const, without changing
  // them, and needs them to be there also when they are empty
  d_Pi.reserve( 1 );
  d_Px.reserve( 1 );
  d_Ai.reserve( 1 );
  d_Ax.reserve( 1 );
  d_l.reserve( 1 );
  d_u.reserve( 1 );
  OSQPCscMatrix P , A;
  OSQPCscMatrix_set_data( &P , d_n , d_n , d_Px.size() , d_Px.data() ,
			  d_Pi.data() , d_Pp.data() );
  OSQPCscMatrix_set_data( &A , d_m , d_n , d_Ax.size() , d_Ax.data() ,
			  d_Ai.data() , d_Ap.data() );
  if( osqp_setup( &d_solver , &P , d_q.data() , &A , d_l.data() ,
		  d_u.data() , d_m , d_n , &d_settings ) != 0 ) {
   free_solver();
   d_loaded = true;
   d_status = OSQP_setup_failed;
   return( sol_status = decode_status() );
   }
  d_sPp = d_Pp;
  d_sPi = d_Pi;
  d_sPx = d_Px;
  d_sAp = d_Ap;
  d_sAi = d_Ai;
  d_sAx = d_Ax;
  d_setup = false;
  }
 else
  if( osqp_update_settings( d_solver , &d_settings ) != 0 ) {
   d_status = OSQP_setup_failed;
   return( sol_status = decode_status() );
   }

 d_loaded = false;

 // the start: the solution of the previous solve, if of the same size,
 // otherwise zero
 if( ws )
  osqp_warm_start( d_solver , d_x.data() , d_y.data() );

 osqp_solve( d_solver );

 d_status = d_solver->info->status_val;
 d_iter_limit = d_solver->info->iter >= d_settings.max_iter;

 // the certificate of an infeasibility, a ray of x for the dual one and a
 // combination of the rows for the primal one
 if( d_status == OSQP_DUAL_INFEASIBLE )
  d_cert.assign( d_solver->solution->dual_inf_cert ,
		 d_solver->solution->dual_inf_cert + d_n );
 if( d_status == OSQP_PRIMAL_INFEASIBLE )
  d_cert.assign( d_solver->solution->prim_inf_cert ,
		 d_solver->solution->prim_inf_cert + d_m );

 // the value is computed here, from the solution, rather than taken from
 // OSQP
 if( has_solution() ) {
  d_x.assign( d_solver->solution->x , d_solver->solution->x + d_n );
  // the bounds are rows of OSQP, satisfied up to the tolerances: x is
  // projected on them, so that a fixed variable gets exactly its value
  for( OSQPInt j = 0 ; j < d_n ; ++j )
   if( d_box_pos[ j ] >= 0 )
    d_x[ j ] = std::clamp( d_x[ j ] , d_l[ d_box_pos[ j ] ] ,
			   d_u[ d_box_pos[ j ] ] );
  d_y.assign( d_solver->solution->y , d_solver->solution->y + d_m );
  for( OSQPInt j = 0 ; j < d_n ; ++j ) {
   d_fval += d_q[ j ] * d_x[ j ];
   for( OSQPInt k = d_Pp[ j ] ; k < d_Pp[ j + 1 ] ; ++k )
    d_fval += ( d_Pi[ k ] == j ? 0.5 : 1.0 ) * d_Px[ k ] *
              d_x[ d_Pi[ k ] ] * d_x[ j ];
   }
  }

 return( sol_status = decode_status() );

 }  // end( OSQPMILPSolver::guts_of_compute )

/*--------------------------------------------------------------------------*/

int OSQPMILPSolver::decode_status( void ) const
{
 switch( d_status ) {
  case( OSQP_SOLVED ):                       return( kOK );
  case( OSQP_SOLVED_INACCURATE ):
   return( d_iter_limit ? kStopIter : kStopTime );
  case( OSQP_PRIMAL_INFEASIBLE ):
  case( OSQP_PRIMAL_INFEASIBLE_INACCURATE ): return( kInfeasible );
  case( OSQP_DUAL_INFEASIBLE ):
  case( OSQP_DUAL_INFEASIBLE_INACCURATE ):   return( kUnbounded );
  case( OSQP_MAX_ITER_REACHED ):             return( kStopIter );
  case( OSQP_TIME_LIMIT_REACHED ):           return( kStopTime );
  default:                                   return( kError );
  }
 }

/*--------------------------------------------------------------------------*/

void OSQPMILPSolver::free_solver( void )
{
 if( d_solver ) {
  osqp_cleanup( d_solver );
  d_solver = nullptr;
  }
 }

/*--------------------------------------------------------------------------*/

bool OSQPMILPSolver::solved( void ) const
{
 return( d_status == OSQP_SOLVED );
 }

/*--------------------------------------------------------------------------*/

bool OSQPMILPSolver::has_solution( void ) const
{
 return( ( d_status == OSQP_SOLVED ) ||
	 ( d_status == OSQP_SOLVED_INACCURATE ) );
 }

/*--------------------------------------------------------------------------*/
/*---------------------- METHODS FOR READING RESULTS -----------------------*/
/*--------------------------------------------------------------------------*/

Solver::OFValue OSQPMILPSolver::get_lb( void )
{
 if( solved() )
  return( d_sign * d_fval + constant_value );

 switch( sol_status ) {
  case( kInfeasible ): return( objsense == -1 ? - Inf< OFValue >()
			                       : Inf< OFValue >() );
  case( kUnbounded ):  return( objsense == -1 ? Inf< OFValue >()
			                       : - Inf< OFValue >() );
  default:             return( - Inf< OFValue >() );
  }
 }

/*--------------------------------------------------------------------------*/

Solver::OFValue OSQPMILPSolver::get_ub( void )
{
 if( solved() )
  return( d_sign * d_fval + constant_value );

 switch( sol_status ) {
  case( kInfeasible ): return( objsense == -1 ? - Inf< OFValue >()
			                       : Inf< OFValue >() );
  case( kUnbounded ):  return( objsense == -1 ? Inf< OFValue >()
			                       : - Inf< OFValue >() );
  default:             return( Inf< OFValue >() );
  }
 }

/*--------------------------------------------------------------------------*/

bool OSQPMILPSolver::has_var_solution( void )
{
 return( has_solution() );
 }

/*--------------------------------------------------------------------------*/

Solver::OFValue OSQPMILPSolver::get_var_value( void )
{
 if( has_solution() )
  return( d_sign * d_fval + constant_value );

 return( objsense == -1 ? get_lb() : get_ub() );
 }

/*--------------------------------------------------------------------------*/

void OSQPMILPSolver::get_var_solution( Configuration * solc )
{
 if( ! has_solution() )
  throw( std::logic_error( "OSQPMILPSolver::get_var_solution: "
			   "no solution available" ) );

 MILPSolver::write_var_solution(
  std::vector< double >( d_x.begin() , d_x.end() ) );
 }

/*--------------------------------------------------------------------------*/

bool OSQPMILPSolver::has_dual_solution( void )
{
 return( has_solution() );
 }

/*--------------------------------------------------------------------------*/

bool OSQPMILPSolver::is_dual_feasible( void )
{
 return( solved() );
 }

/*--------------------------------------------------------------------------*/

void OSQPMILPSolver::get_dual_solution( Configuration * solc )
{
 if( ! has_solution() )
  throw( std::logic_error( "OSQPMILPSolver::get_dual_solution: "
			   "no solution available" ) );

 // OSQP has P x + q + A' y = 0, with y > 0 on an active upper side and
 // y < 0 on an active lower one, for the problem it minimizes: the
 // multipliers of the problem of the Block are - y for a minimization and
 // y for a maximization, whose Objective OSQP has been given negated; that
 // of a row left out is zero
 std::vector< double > pi( d_row_pos.size() , 0 );
 for( std::size_t i = 0 ; i < d_row_pos.size() ; ++i )
  if( d_row_pos[ i ] >= 0 )
   pi[ i ] = - d_sign * d_y[ d_row_pos[ i ] ];

 // the multipliers of the bounds are those of the unit rows
 std::vector< double > rc( d_n , 0 );
 for( OSQPInt j = 0 ; j < d_n ; ++j )
  if( d_box_pos[ j ] >= 0 )
   rc[ j ] = - d_sign * d_y[ d_box_pos[ j ] ];

 MILPSolver::write_dual_solution( pi , rc );
 }

/*--------------------------------------------------------------------------*/

bool OSQPMILPSolver::has_var_direction( void )
{
 return( ( d_status == OSQP_DUAL_INFEASIBLE ) &&
	 ( d_cert.size() == std::size_t( d_n ) ) );
 }

/*--------------------------------------------------------------------------*/

void OSQPMILPSolver::get_var_direction( Configuration * dirc )
{
 if( ! has_var_direction() )
  throw( std::logic_error( "OSQPMILPSolver::get_var_direction: "
			   "no unbounded direction available" ) );

 // P x = 0 and q' x < 0 for the problem OSQP minimizes, whose q is that of
 // the Block negated for a maximization: x is a ray of the Block
 MILPSolver::write_var_solution(
		       std::vector< double >( d_cert.begin() , d_cert.end() ) );
 }

/*--------------------------------------------------------------------------*/

bool OSQPMILPSolver::has_dual_direction( void )
{
 return( ( d_status == OSQP_PRIMAL_INFEASIBLE ) &&
	 ( d_cert.size() == std::size_t( d_m ) ) );
 }

/*--------------------------------------------------------------------------*/

void OSQPMILPSolver::get_dual_direction( Configuration * dirc )
{
 if( ! has_dual_direction() )
  throw( std::logic_error( "OSQPMILPSolver::get_dual_direction: "
			   "no dual unbounded direction available" ) );

 // A' y = 0 with y > 0 on the upper side of a row and y < 0 on the lower
 // one: the certificates of the other :MILPSolver are nonnegative where
 // the lower side is used, hence the multiplier of a row is - y, that of a
 // row left out is zero, and the value of the certificate is
 // - u' max( y , 0 ) - l' min( y , 0 )
 double value = 0;
 for( OSQPInt r = 0 ; r < d_m ; ++r )
  if( d_cert[ r ] > 0 )
   value -= d_u[ r ] * d_cert[ r ];
  else
   if( d_cert[ r ] < 0 )
    value -= d_l[ r ] * d_cert[ r ];

 std::vector< double > pi( d_row_pos.size() , 0 );
 for( std::size_t i = 0 ; i < d_row_pos.size() ; ++i )
  if( d_row_pos[ i ] >= 0 )
   pi[ i ] = - d_cert[ d_row_pos[ i ] ];

 // the multipliers of the columns, c - A' y, or - A' y with
 // intHomogeneousDirection, computed out of the matrix so that they are
 // those of y exactly
 std::vector< double > dj( d_n , 0 );
 for( OSQPInt j = 0 ; j < d_n ; ++j ) {
  if( ! homogeneous_direction )
   dj[ j ] = objective[ j ];
  for( int k = matbeg[ j ] ; k < matbeg[ j ] + matcnt[ j ] ; ++k )
   dj[ j ] -= matval[ k ] * pi[ matind[ k ] ];
  }

 f_dual_direction_value = value;

 MILPSolver::write_dual_solution( pi , dj );
 }

/*--------------------------------------------------------------------------*/
/*------------------ METHODS FOR HANDLING THE PARAMETERS -------------------*/
/*--------------------------------------------------------------------------*/

void OSQPMILPSolver::set_par( idx_type par , int value )
{
 switch( par ) {
  case( intMaxIter ): d_max_iter = std::max( value , 1 );    return;
  case( intOSQPWarmStart ):  d_warm_start = bool( value );   return;
  case( intOSQPPolishing ):  d_settings.polishing = bool( value ); return;
  case( intOSQPPolishRefineIter ):
   d_settings.polish_refine_iter = std::max( value , 0 );
   return;
  case( intOSQPScaling ):      // the parameters of the set up, which is
  case( intOSQPAdaptiveRho ): {  // done again at the next solve
   OSQPInt & p = par == intOSQPScaling ? d_settings.scaling
                                       : d_settings.adaptive_rho;
   const OSQPInt v = par == intOSQPScaling ? std::max( value , 0 )
                                           : std::clamp( value , 0 , 3 );
   if( v != p ) {
    p = v;
    d_setup = true;
    }
   return;
   }
  case( intOSQPCheckTermination ):
   d_settings.check_termination = std::max( value , 0 );
   return;
  }

 MILPSolver::set_par( par , value );
 }

/*--------------------------------------------------------------------------*/

void OSQPMILPSolver::set_par( idx_type par , double value )
{
 switch( par ) {
  case( dblMaxTime ):        d_max_time = std::max( value , 0.0 ); return;
  case( dblFAccSol ):        d_settings.eps_abs = value;      return;
  case( dblOSQPEpsRel ):     d_settings.eps_rel = value;      return;
  case( dblOSQPEpsPrimInf ): d_settings.eps_prim_inf = value; return;
  case( dblOSQPEpsDualInf ): d_settings.eps_dual_inf = value; return;
  case( dblOSQPRho ):          // the parameters of the set up, which is
  case( dblOSQPSigma ): {      // done again at the next solve
   OSQPFloat & p = par == dblOSQPRho ? d_settings.rho : d_settings.sigma;
   if( value != p ) {
    p = value;
    d_setup = true;
    }
   return;
   }
  case( dblOSQPAlpha ):      d_settings.alpha = value;        return;
  case( dblOSQPDelta ):      d_settings.delta = value;        return;
  }

 MILPSolver::set_par( par , value );
 }

/*--------------------------------------------------------------------------*/

int OSQPMILPSolver::get_dflt_int_par( idx_type par ) const
{
 switch( par ) {
  case( intOSQPWarmStart ):        return( 1 );
  case( intOSQPPolishing ):        return( 1 );
  case( intOSQPPolishRefineIter ): return( OSQP_dflt.polish_refine_iter );
  case( intOSQPScaling ):          return( OSQP_dflt.scaling );
  case( intOSQPAdaptiveRho ):      return( OSQP_dflt_adaptive_rho );
  case( intOSQPCheckTermination ): return( OSQP_dflt.check_termination );
  }

 return( MILPSolver::get_dflt_int_par( par ) );
 }

/*--------------------------------------------------------------------------*/

double OSQPMILPSolver::get_dflt_dbl_par( idx_type par ) const
{
 switch( par ) {
  case( dblFAccSol ):        return( OSQP_dflt_eps );
  case( dblOSQPEpsRel ):     return( OSQP_dflt_eps );
  case( dblOSQPEpsPrimInf ): return( OSQP_dflt_eps_inf );
  case( dblOSQPEpsDualInf ): return( OSQP_dflt_eps_inf );
  case( dblOSQPRho ):        return( OSQP_dflt.rho );
  case( dblOSQPSigma ):      return( OSQP_dflt.sigma );
  case( dblOSQPAlpha ):      return( OSQP_dflt.alpha );
  case( dblOSQPDelta ):      return( OSQP_dflt.delta );
  }

 return( MILPSolver::get_dflt_dbl_par( par ) );
 }

/*--------------------------------------------------------------------------*/

int OSQPMILPSolver::get_int_par( idx_type par ) const
{
 switch( par ) {
  case( intMaxIter ):              return( d_max_iter );
  case( intOSQPWarmStart ):        return( d_warm_start );
  case( intOSQPPolishing ):        return( d_settings.polishing );
  case( intOSQPPolishRefineIter ): return( d_settings.polish_refine_iter );
  case( intOSQPScaling ):          return( d_settings.scaling );
  case( intOSQPAdaptiveRho ):      return( d_settings.adaptive_rho );
  case( intOSQPCheckTermination ): return( d_settings.check_termination );
  }

 return( MILPSolver::get_int_par( par ) );
 }

/*--------------------------------------------------------------------------*/

double OSQPMILPSolver::get_dbl_par( idx_type par ) const
{
 switch( par ) {
  case( dblMaxTime ):        return( d_max_time );
  case( dblFAccSol ):        return( d_settings.eps_abs );
  case( dblOSQPEpsRel ):     return( d_settings.eps_rel );
  case( dblOSQPEpsPrimInf ): return( d_settings.eps_prim_inf );
  case( dblOSQPEpsDualInf ): return( d_settings.eps_dual_inf );
  case( dblOSQPRho ):        return( d_settings.rho );
  case( dblOSQPSigma ):      return( d_settings.sigma );
  case( dblOSQPAlpha ):      return( d_settings.alpha );
  case( dblOSQPDelta ):      return( d_settings.delta );
  }

 return( MILPSolver::get_dbl_par( par ) );
 }

/*--------------------------------------------------------------------------*/

Solver::idx_type OSQPMILPSolver::int_par_str2idx( const std::string & name )
 const
{
 auto it = std::find( OSQP_int_pars_str.begin() ,
		      OSQP_int_pars_str.end() , name );
 if( it != OSQP_int_pars_str.end() )
  return( intLastAlgParMILP + ( it - OSQP_int_pars_str.begin() ) );

 return( MILPSolver::int_par_str2idx( name ) );
 }

/*--------------------------------------------------------------------------*/

const std::string & OSQPMILPSolver::int_par_idx2str( idx_type idx ) const
{
 if( ( idx >= intLastAlgParMILP ) && ( idx < intLastAlgParOSQP ) )
  return( OSQP_int_pars_str[ idx - intLastAlgParMILP ] );

 return( MILPSolver::int_par_idx2str( idx ) );
 }

/*--------------------------------------------------------------------------*/

Solver::idx_type OSQPMILPSolver::dbl_par_str2idx( const std::string & name )
 const
{
 auto it = std::find( OSQP_dbl_pars_str.begin() ,
		      OSQP_dbl_pars_str.end() , name );
 if( it != OSQP_dbl_pars_str.end() )
  return( dblLastAlgParMILP + ( it - OSQP_dbl_pars_str.begin() ) );

 return( MILPSolver::dbl_par_str2idx( name ) );
 }

/*--------------------------------------------------------------------------*/

const std::string & OSQPMILPSolver::dbl_par_idx2str( idx_type idx ) const
{
 if( ( idx >= dblLastAlgParMILP ) && ( idx < dblLastAlgParOSQP ) )
  return( OSQP_dbl_pars_str[ idx - dblLastAlgParMILP ] );

 return( MILPSolver::dbl_par_idx2str( idx ) );
 }

/*--------------------------------------------------------------------------*/
/*---------------------- End File OSQPMILPSolver.cpp -----------------------*/
/*--------------------------------------------------------------------------*/
