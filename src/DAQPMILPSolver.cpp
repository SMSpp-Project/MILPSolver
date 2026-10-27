/*--------------------------------------------------------------------------*/
/*------------------------- File DAQPMILPSolver.cpp ------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Implementation of the DAQPMILPSolver class, which implements the
 * MILPSolver interface for convex Quadratic Programs with linear
 * constraints using the dual active-set solver DAQP.
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
#include <cmath>
#include <limits>
#include <stdexcept>

#include "DAQPMILPSolver.h"

/*--------------------------------------------------------------------------*/
/*-------------------------------- USING -----------------------------------*/
/*--------------------------------------------------------------------------*/

using namespace SMSpp_di_unipi_it;

/*--------------------------------------------------------------------------*/
/*-------------------------------- CONSTANTS -------------------------------*/
/*--------------------------------------------------------------------------*/

static const std::array< std::string , 3 > DAQP_int_pars_str = {
 "intDAQPWarmStart" , "intDAQPEqReduction" , "intDAQPCycleTol" };

static const std::array< std::string , 5 > DAQP_dbl_pars_str = {
 "dblDAQPDualTol" , "dblDAQPZeroTol" , "dblDAQPPivotTol" ,
 "dblDAQPEpsProx" , "dblDAQPEtaProx" };

/*--------------------------------------------------------------------------*/
/*-------------------------------- FUNCTIONS -------------------------------*/
/*--------------------------------------------------------------------------*/

// a side of SMS++ in the format of DAQP, where an infinite one is DAQP_INF

static double daqp_side( double v )
{
 if( v >= DAQP_INF )
  return( DAQP_INF );
 if( v <= - DAQP_INF )
  return( - DAQP_INF );
 return( v );
 }

/*--------------------------------------------------------------------------*/
/*------------------------- FACTORY MANAGEMENT -----------------------------*/
/*--------------------------------------------------------------------------*/

SMSpp_insert_in_factory_cpp_0( DAQPMILPSolver );

/*--------------------------------------------------------------------------*/
/*--------------------- CONSTRUCTOR AND DESTRUCTOR -------------------------*/
/*--------------------------------------------------------------------------*/

DAQPMILPSolver::DAQPMILPSolver( void ) : MILPSolver()
{
 daqp_default_settings( & d_settings );
 }

/*--------------------------------------------------------------------------*/

DAQPMILPSolver::~DAQPMILPSolver()
{
 free_workspace();
 }

/*--------------------------------------------------------------------------*/
/*--------------------- DERIVED METHODS OF BASE CLASS ----------------------*/
/*--------------------------------------------------------------------------*/

void DAQPMILPSolver::clear_problem( unsigned int what )
{
 MILPSolver::clear_problem( what );
 free_workspace();
 d_ws.clear();
 d_ws_n = -1;
 d_n = d_m = 0;
 d_loaded = false;
 d_exitflag = 0;
 }

/*--------------------------------------------------------------------------*/

void DAQPMILPSolver::load_problem( void )
{
 // the workspace points to the old data, which are about to change; its
 // working set is kept across the clear_problem() of MILPSolver, so that
 // the problem loaded again whole after a Modification starts from it
 free_workspace();
 auto ws = std::move( d_ws );
 const int ws_n = d_ws_n;

 MILPSolver::load_problem();

 d_ws = std::move( ws );
 d_ws_n = ws_n;

 if( numquadrows > 0 )
  throw( std::invalid_argument( "DAQPMILPSolver::load_problem: "
				"quadratic constraints are not supported" ) );

 if( ( int_vars > 0 ) && ( relax_int_vars == 0 ) )
  throw( std::invalid_argument( "DAQPMILPSolver::load_problem: integer "
				"variables are only supported if relaxed "
				"(intRelaxIntVars)" ) );

 const int n = numcols;
 d_n = n;
 d_m = n + numrows;
 d_sign = ( objsense == -1 ) ? -1 : 1;  // DAQP minimizes

 // the objective: c' x + sum_j q_j x_j^2 + sum_k ndq_k x_i x_j is
 // 1/2 x' H x + f' x with H_jj = 2 q_j and H_ij = H_ji = ndq_k
 d_f.assign( n , 0 );
 for( int j = 0 ; j < n ; ++j )
  d_f[ j ] = d_sign * objective[ j ];

 d_H.assign( std::size_t( n ) * n , 0 );
 for( int j = 0 ; j < n ; ++j )
  d_H[ std::size_t( j ) * n + j ] = d_sign * 2 * q_objective[ j ];
 for( std::size_t k = 0 ; k < ndq_objective.size() ; ++k ) {
  const int i = ndq_rowind[ k ];
  const int j = ndq_colind[ k ];
  const double v = d_sign * ndq_objective[ k ];
  if( i == j )
   d_H[ std::size_t( i ) * n + i ] += 2 * v;
  else {
   d_H[ std::size_t( i ) * n + j ] += v;
   d_H[ std::size_t( j ) * n + i ] += v;
   }
  }

 // the rows, dense and row-major, from the columns of MILPSolver
 d_A.assign( std::size_t( numrows ) * n , 0 );
 for( int j = 0 ; j < n ; ++j )
  for( int k = matbeg[ j ] ; k < matbeg[ j ] + matcnt[ j ] ; ++k )
   d_A[ std::size_t( matind[ k ] ) * n + j ] += matval[ k ];

 // the sides: the bounds first, then the rows; a fixed variable and an
 // equality row are constraints that DAQP keeps active
 d_bu.resize( d_m );
 d_bl.resize( d_m );
 d_fixed.assign( d_m , 0 );
 for( int j = 0 ; j < n ; ++j ) {
  d_bl[ j ] = daqp_side( lb[ j ] );
  d_bu[ j ] = daqp_side( ub[ j ] );
  if( lb[ j ] == ub[ j ] )
   d_fixed[ j ] = DAQP_ACTIVE | DAQP_IMMUTABLE;
  }

 for( int i = 0 ; i < numrows ; ++i ) {
  const int r = n + i;
  switch( sense[ i ] ) {
   case( 'L' ): d_bl[ r ] = - DAQP_INF;
                d_bu[ r ] = daqp_side( rhs[ i ] );
                break;
   case( 'G' ): d_bl[ r ] = daqp_side( rhs[ i ] );
                d_bu[ r ] = DAQP_INF;
                break;
   case( 'E' ): d_bl[ r ] = d_bu[ r ] = rhs[ i ];
                d_fixed[ r ] = DAQP_ACTIVE | DAQP_IMMUTABLE;
                break;
   case( 'R' ): if( rngval[ i ] > 0 ) {
                 d_bl[ r ] = daqp_side( rhs[ i ] );
                 d_bu[ r ] = daqp_side( rhs[ i ] + rngval[ i ] );
                 }
                else {
                 d_bl[ r ] = daqp_side( rhs[ i ] + rngval[ i ] );
                 d_bu[ r ] = daqp_side( rhs[ i ] );
                 }
                if( d_bl[ r ] == d_bu[ r ] )
                 d_fixed[ r ] = DAQP_ACTIVE | DAQP_IMMUTABLE;
                break;
   default:
    throw( std::invalid_argument( "DAQPMILPSolver::load_problem: row " +
				  std::to_string( i ) + " has sense " +
				  std::string( 1 , sense[ i ] ) ) );
   }
  }

 d_x.assign( n , 0 );
 d_lam.assign( d_m , 0 );
 d_loaded = true;
 d_exitflag = 0;

 }  // end( DAQPMILPSolver::load_problem )

/*--------------------------------------------------------------------------*/
/*-------------------- METHODS FOR SOLVING THE MODEL -----------------------*/
/*--------------------------------------------------------------------------*/

int DAQPMILPSolver::guts_of_compute( void )
{
 if( d_loaded || ( ! d_has_work ) ) {
  free_workspace();

  // the constraints the problem fixes, and the working set of the previous
  // solve on the others, if any and of the same size; a constraint that the
  // previous problem fixed and this one does not is left out, since it was
  // active for being fixed rather than at the solution
  d_sense = d_fixed;
  if( d_warm_start && ( d_ws_n == d_n ) &&
      ( d_ws.size() == std::size_t( d_m ) ) )
   for( int i = 0 ; i < d_m ; ++i )
    if( ( ! d_fixed[ i ] ) && ( ! ( d_ws[ i ] & DAQP_IMMUTABLE ) ) )
     d_sense[ i ] |= d_ws[ i ] & ( DAQP_ACTIVE | 2 );

  d_qp = DAQPProblem{};
  d_qp.n = d_n;
  d_qp.m = d_m;
  d_qp.ms = d_n;
  d_qp.H = d_H.data();
  d_qp.f = d_f.data();
  d_qp.A = d_A.empty() ? nullptr : d_A.data();
  d_qp.bupper = d_bu.data();
  d_qp.blower = d_bl.data();
  d_qp.sense = d_sense.data();

  d_work = DAQPWorkspace{};
  d_work.settings = & d_settings;
  double setup_time = 0;
  const int flag = setup_daqp( & d_qp , & d_work , & setup_time );
  if( flag < 0 ) {
   d_work.settings = nullptr;
   free_daqp_workspace( & d_work );
   free_daqp_ldp( & d_work );
   d_exitflag = flag;
   d_loaded = false;
   return( sol_status = decode_exitflag( flag ) );
   }
  d_has_work = true;
  d_loaded = false;
  }

 // DAQP keeps its own iteration limit if none is set, and has no time
 // limit if its own is 0, so that a time limit of 0 is the smallest positive
 // one
 d_settings.iter_limit = ( d_max_iter < Inf< int >() ) ?
  d_max_iter : DAQP_DEFAULT_ITER_LIMIT;
 d_settings.time_limit = ( d_max_time < Inf< double >() ) ?
  std::max( d_max_time , std::numeric_limits< double >::min() ) : 0;

 d_res = DAQPResult{};
 d_res.x = d_x.data();
 d_res.lam = d_lam.data();
 daqp_solve( & d_res , & d_work );
 d_exitflag = d_res.exitflag;
 elapsed_iter += d_res.iter;

 // the value is computed here rather than taken from DAQP, whose value of a
 // problem solved by proximal-point iterations (e.g., a Linear Program) is
 // not that of the problem
 d_fval = 0;
 if( solved() )
  for( int i = 0 ; i < d_n ; ++i ) {
   double Hx = 0;
   for( int j = 0 ; j < d_n ; ++j )
    Hx += d_H[ std::size_t( i ) * d_n + j ] * d_x[ j ];
   d_fval += d_x[ i ] * ( d_f[ i ] + Hx / 2 );
   }

 return( sol_status = decode_exitflag( d_exitflag ) );

 }  // end( DAQPMILPSolver::guts_of_compute )

/*--------------------------------------------------------------------------*/

int DAQPMILPSolver::decode_exitflag( int exitflag )
{
 switch( exitflag ) {
  case( DAQP_EXIT_OPTIMAL ):
 #ifdef DAQP_EXIT_OPTIMAL_INEXACT
  case( DAQP_EXIT_OPTIMAL_INEXACT ):
 #endif
  case( DAQP_EXIT_SOFT_OPTIMAL ):
   return( kOK );
  case( DAQP_EXIT_INFEASIBLE ): return( kInfeasible );
  case( DAQP_EXIT_UNBOUNDED ):  return( kUnbounded );
  case( DAQP_EXIT_ITERLIMIT ):  return( kStopIter );
  case( DAQP_EXIT_TIMELIMIT ):  return( kStopTime );
  default:                      return( kError );
  }
 }

/*--------------------------------------------------------------------------*/

void DAQPMILPSolver::free_workspace( void )
{
 if( ! d_has_work )
  return;

 // the working set, for the next start
 d_ws.assign( d_work.sense , d_work.sense + d_work.m );
 d_ws_n = d_work.n;

 d_work.settings = nullptr;  // they are ours
 free_daqp_workspace( & d_work );
 free_daqp_ldp( & d_work );
 d_has_work = false;
 }

/*--------------------------------------------------------------------------*/

bool DAQPMILPSolver::solved( void ) const
{
 return( decode_exitflag( d_exitflag ) == kOK );
 }

/*--------------------------------------------------------------------------*/
/*---------------------- METHODS FOR READING RESULTS -----------------------*/
/*--------------------------------------------------------------------------*/

Solver::OFValue DAQPMILPSolver::get_lb( void )
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

Solver::OFValue DAQPMILPSolver::get_ub( void )
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

bool DAQPMILPSolver::has_var_solution( void )
{
 return( solved() );
 }

/*--------------------------------------------------------------------------*/

Solver::OFValue DAQPMILPSolver::get_var_value( void )
{
 return( objsense == -1 ? get_lb() : get_ub() );
 }

/*--------------------------------------------------------------------------*/

void DAQPMILPSolver::get_var_solution( Configuration * solc )
{
 if( ! solved() )
  throw( std::logic_error( "DAQPMILPSolver::get_var_solution: "
			   "no solution available" ) );

 MILPSolver::write_var_solution( d_x , solc );
 }

/*--------------------------------------------------------------------------*/

bool DAQPMILPSolver::has_dual_solution( void )
{
 return( solved() );
 }

/*--------------------------------------------------------------------------*/

bool DAQPMILPSolver::is_dual_feasible( void )
{
 return( solved() );
 }

/*--------------------------------------------------------------------------*/

void DAQPMILPSolver::get_dual_solution( Configuration * solc )
{
 if( ! solved() )
  throw( std::logic_error( "DAQPMILPSolver::get_dual_solution: "
			   "no solution available" ) );

 // DAQP has H x + f + A' lam = 0, with lam > 0 on an active upper side and
 // lam < 0 on an active lower one, for the problem it minimizes: the
 // multipliers of the problem of the Block are - lam for a minimization
 // and lam for a maximization, whose Objective DAQP has been given negated
 std::vector< double > rc( d_n );
 for( int j = 0 ; j < d_n ; ++j )
  rc[ j ] = - d_sign * d_lam[ j ];

 std::vector< double > pi( d_m - d_n );
 for( int i = d_n ; i < d_m ; ++i )
  pi[ i - d_n ] = - d_sign * d_lam[ i ];

 MILPSolver::write_dual_solution( pi , rc , solc );
 }

/*--------------------------------------------------------------------------*/
/*------------------ METHODS FOR HANDLING THE PARAMETERS -------------------*/
/*--------------------------------------------------------------------------*/

void DAQPMILPSolver::set_par( idx_type par , int value )
{
 switch( par ) {
  case( intMaxIter ): d_max_iter = std::max( value , 0 ); return;
  case( intDAQPWarmStart ):   d_warm_start = bool( value );       return;
  case( intDAQPEqReduction ):
 #ifdef DAQP_EQ_REDUCTION_AUTO
   d_settings.eq_reduction = value;
 #endif
   return;
  case( intDAQPCycleTol ):    d_settings.cycle_tol = value;       return;
  }

 MILPSolver::set_par( par , value );
 }

/*--------------------------------------------------------------------------*/

void DAQPMILPSolver::set_par( idx_type par , double value )
{
 switch( par ) {
  case( dblMaxTime ): d_max_time = std::max( value , 0.0 ); return;
  case( dblFAccSol ):     d_settings.primal_tol = value;  return;
  case( dblDAQPDualTol ): d_settings.dual_tol = value;    return;
  case( dblDAQPZeroTol ): d_settings.zero_tol = value;    return;
  case( dblDAQPPivotTol ): d_settings.pivot_tol = value;  return;
  case( dblDAQPEpsProx ): d_settings.eps_prox = value;    return;
  case( dblDAQPEtaProx ): d_settings.eta_prox = value;    return;
  }

 MILPSolver::set_par( par , value );
 }

/*--------------------------------------------------------------------------*/

int DAQPMILPSolver::get_dflt_int_par( idx_type par ) const
{
 switch( par ) {
  case( intDAQPWarmStart ):   return( 1 );
  case( intDAQPEqReduction ):
 #ifdef DAQP_EQ_REDUCTION_AUTO
   return( DAQP_EQ_REDUCTION_AUTO );
 #else
   return( 0 );
 #endif
  case( intDAQPCycleTol ):    return( DAQP_DEFAULT_CYCLE_TOL );
  }

 return( MILPSolver::get_dflt_int_par( par ) );
 }

/*--------------------------------------------------------------------------*/

double DAQPMILPSolver::get_dflt_dbl_par( idx_type par ) const
{
 switch( par ) {
  case( dblFAccSol ):      return( DAQP_DEFAULT_PRIM_TOL );
  case( dblDAQPDualTol ):  return( DAQP_DEFAULT_DUAL_TOL );
  case( dblDAQPZeroTol ):  return( DAQP_DEFAULT_ZERO_TOL );
  case( dblDAQPPivotTol ): return( DAQP_DEFAULT_PIVOT_TOL );
  case( dblDAQPEpsProx ):  return( DAQP_DEFAULT_EPS_PROX );
  case( dblDAQPEtaProx ):  return( DAQP_DEFAULT_ETA );
  }

 return( MILPSolver::get_dflt_dbl_par( par ) );
 }

/*--------------------------------------------------------------------------*/

int DAQPMILPSolver::get_int_par( idx_type par ) const
{
 switch( par ) {
  case( intMaxIter ): return( d_max_iter );
  case( intDAQPWarmStart ):   return( d_warm_start );
  case( intDAQPEqReduction ):
 #ifdef DAQP_EQ_REDUCTION_AUTO
   return( d_settings.eq_reduction );
 #else
   return( 0 );
 #endif
  case( intDAQPCycleTol ):    return( d_settings.cycle_tol );
  }

 return( MILPSolver::get_int_par( par ) );
 }

/*--------------------------------------------------------------------------*/

double DAQPMILPSolver::get_dbl_par( idx_type par ) const
{
 switch( par ) {
  case( dblMaxTime ): return( d_max_time );
  case( dblFAccSol ):      return( d_settings.primal_tol );
  case( dblDAQPDualTol ):  return( d_settings.dual_tol );
  case( dblDAQPZeroTol ):  return( d_settings.zero_tol );
  case( dblDAQPPivotTol ): return( d_settings.pivot_tol );
  case( dblDAQPEpsProx ):  return( d_settings.eps_prox );
  case( dblDAQPEtaProx ):  return( d_settings.eta_prox );
  }

 return( MILPSolver::get_dbl_par( par ) );
 }

/*--------------------------------------------------------------------------*/

Solver::idx_type DAQPMILPSolver::int_par_str2idx( const std::string & name )
 const
{
 auto it = std::find( DAQP_int_pars_str.begin() , DAQP_int_pars_str.end() ,
		      name );
 if( it != DAQP_int_pars_str.end() )
  return( intLastAlgParMILP + ( it - DAQP_int_pars_str.begin() ) );

 return( MILPSolver::int_par_str2idx( name ) );
 }

/*--------------------------------------------------------------------------*/

const std::string & DAQPMILPSolver::int_par_idx2str( idx_type idx ) const
{
 if( ( idx >= intLastAlgParMILP ) && ( idx < intLastAlgParDAQP ) )
  return( DAQP_int_pars_str[ idx - intLastAlgParMILP ] );

 return( MILPSolver::int_par_idx2str( idx ) );
 }

/*--------------------------------------------------------------------------*/

Solver::idx_type DAQPMILPSolver::dbl_par_str2idx( const std::string & name )
 const
{
 auto it = std::find( DAQP_dbl_pars_str.begin() , DAQP_dbl_pars_str.end() ,
		      name );
 if( it != DAQP_dbl_pars_str.end() )
  return( dblLastAlgParMILP + ( it - DAQP_dbl_pars_str.begin() ) );

 return( MILPSolver::dbl_par_str2idx( name ) );
 }

/*--------------------------------------------------------------------------*/

const std::string & DAQPMILPSolver::dbl_par_idx2str( idx_type idx ) const
{
 if( ( idx >= dblLastAlgParMILP ) && ( idx < dblLastAlgParDAQP ) )
  return( DAQP_dbl_pars_str[ idx - dblLastAlgParMILP ] );

 return( MILPSolver::dbl_par_idx2str( idx ) );
 }

/*--------------------------------------------------------------------------*/
/*---------------------- End File DAQPMILPSolver.cpp -----------------------*/
/*--------------------------------------------------------------------------*/
