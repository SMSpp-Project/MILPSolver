/*--------------------------------------------------------------------------*/
/*------------------------- File PIQPMILPSolver.cpp ------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Implementation of the PIQPMILPSolver class, which implements the
 * MILPSolver interface for convex Quadratic Programs with linear
 * constraints using the proximal interior point solver PIQP.
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

#include <piqp/piqp.hpp>

#include "PIQPMILPSolver.h"

/*--------------------------------------------------------------------------*/
/*-------------------------------- USING -----------------------------------*/
/*--------------------------------------------------------------------------*/

using namespace SMSpp_di_unipi_it;

/*--------------------------------------------------------------------------*/
/*-------------------------------- CONSTANTS -------------------------------*/
/*--------------------------------------------------------------------------*/

static const std::array< std::string , 5 > PIQP_int_pars_str = {
 "intPIQPBackend" , "intPIQPCheckDualityGap" , "intPIQPPrecondIter" ,
 "intPIQPPrecondScaleCost" , "intPIQPIterRefinement" };

static const std::array< std::string , 6 > PIQP_dbl_pars_str = {
 "dblPIQPEpsRel" , "dblPIQPEpsGapAbs" , "dblPIQPEpsGapRel" ,
 "dblPIQPInfeasThreshold" , "dblPIQPRhoInit" , "dblPIQPDeltaInit" };

// the settings of PIQP, which give the defaults of the parameters

static const piqp::Settings< double > PIQP_dflt;

// the default of dblFAccSol and of dblPIQPEpsGapAbs, tighter than the 1e-8
// of PIQP, which leaves errors of that order in the solution and its value

static constexpr double PIQP_dflt_eps_abs = 1e-9;

// the position of a row with no finite side, which is left out of PIQP

static constexpr int PIQP_free_row = std::numeric_limits< int >::min();

/*--------------------------------------------------------------------------*/
/*------------------------- THE SOLVERS OF PIQP ----------------------------*/
/*--------------------------------------------------------------------------*/

struct PIQPMILPSolver::Solvers {
 std::unique_ptr< piqp::DenseSolver< double > > dense;
 std::unique_ptr< piqp::SparseSolver< double , int > > sparse;
 };

/*--------------------------------------------------------------------------*/
/*------------------------- FACTORY MANAGEMENT -----------------------------*/
/*--------------------------------------------------------------------------*/

SMSpp_insert_in_factory_cpp_0( PIQPMILPSolver );

/*--------------------------------------------------------------------------*/
/*--------------------- CONSTRUCTOR AND DESTRUCTOR -------------------------*/
/*--------------------------------------------------------------------------*/

PIQPMILPSolver::PIQPMILPSolver( void ) : MILPSolver()
{
 d_check_gap = PIQP_dflt.check_duality_gap;
 d_precond_iter = int( PIQP_dflt.preconditioner_iter );
 d_precond_cost = PIQP_dflt.preconditioner_scale_cost;
 d_iter_ref = PIQP_dflt.iterative_refinement_always_enabled;
 d_eps_abs = PIQP_dflt_eps_abs;
 d_eps_rel = PIQP_dflt.eps_rel;
 d_eps_gap_abs = PIQP_dflt_eps_abs;
 d_eps_gap_rel = PIQP_dflt.eps_duality_gap_rel;
 d_infeas_thr = PIQP_dflt.infeasibility_threshold;
 d_rho_init = PIQP_dflt.rho_init;
 d_delta_init = PIQP_dflt.delta_init;
 }

/*--------------------------------------------------------------------------*/

PIQPMILPSolver::~PIQPMILPSolver()
{
 free_solver();
 }

/*--------------------------------------------------------------------------*/
/*--------------------- DERIVED METHODS OF BASE CLASS ----------------------*/
/*--------------------------------------------------------------------------*/

void PIQPMILPSolver::clear_problem( unsigned int what )
{
 MILPSolver::clear_problem( what );
 free_solver();
 d_n = d_neq = d_nin = 0;
 d_row_pos.clear();
 d_loaded = false;
 d_status = 0;
 }

/*--------------------------------------------------------------------------*/

void PIQPMILPSolver::load_problem( void )
{
 MILPSolver::load_problem();
 free_solver();

 if( numquadrows > 0 )
  throw( std::invalid_argument( "PIQPMILPSolver::load_problem: "
				"quadratic constraints are not supported" ) );

 if( ( int_vars > 0 ) && ( relax_int_vars == 0 ) )
  throw( std::invalid_argument( "PIQPMILPSolver::load_problem: integer "
				"variables are only supported if relaxed "
				"(intRelaxIntVars)" ) );

 const int n = numcols;
 d_n = n;
 d_sign = ( objsense == -1 ) ? -1 : 1;  // PIQP minimizes

 // the objective: c' x + sum_j q_j x_j^2 + sum_k ndq_k x_i x_j is
 // 1/2 x' P x + c' x with P_jj = 2 q_j and P_ij = P_ji = ndq_k
 d_c.resize( n );
 for( int j = 0 ; j < n ; ++j )
  d_c[ j ] = d_sign * objective[ j ];

 std::vector< Eigen::Triplet< double > > trp;
 for( int j = 0 ; j < n ; ++j )
  if( q_objective[ j ] != 0 )
   trp.emplace_back( j , j , d_sign * 2 * q_objective[ j ] );
 for( std::size_t k = 0 ; k < ndq_objective.size() ; ++k ) {
  const int i = ndq_rowind[ k ];
  const int j = ndq_colind[ k ];
  const double v = d_sign * ndq_objective[ k ];
  if( i == j )
   trp.emplace_back( i , i , 2 * v );
  else {
   trp.emplace_back( i , j , v );
   trp.emplace_back( j , i , v );
   }
  }
 d_P.resize( n , n );
 d_P.setFromTriplets( trp.begin() , trp.end() );
 d_P.prune( 0.0 );

 // the sides of the rows: those whose two sides coincide are equalities,
 // which go in A, those with no finite side (as PIQP sees it) are left out,
 // and all the others go in G
 std::vector< double > lo( numrows ) , hi( numrows );
 d_row_pos.resize( numrows );
 d_neq = d_nin = 0;
 for( int i = 0 ; i < numrows ; ++i ) {
  switch( sense[ i ] ) {
   case( 'L' ): lo[ i ] = - Inf< double >();
                hi[ i ] = rhs[ i ];
                break;
   case( 'G' ): lo[ i ] = rhs[ i ];
                hi[ i ] = Inf< double >();
                break;
   case( 'E' ): lo[ i ] = hi[ i ] = rhs[ i ];
                break;
   case( 'R' ): if( rngval[ i ] > 0 ) {
                 lo[ i ] = rhs[ i ];
                 hi[ i ] = rhs[ i ] + rngval[ i ];
                 }
                else {
                 lo[ i ] = rhs[ i ] + rngval[ i ];
                 hi[ i ] = rhs[ i ];
                 }
                break;
   default:
    throw( std::invalid_argument( "PIQPMILPSolver::load_problem: row " +
				  std::to_string( i ) + " has sense " +
				  std::string( 1 , sense[ i ] ) ) );
   }
  if( lo[ i ] == hi[ i ] )
   d_row_pos[ i ] = d_neq++;
  else
   if( ( lo[ i ] <= - PIQP_INF ) && ( hi[ i ] >= PIQP_INF ) )
    d_row_pos[ i ] = PIQP_free_row;
   else
    d_row_pos[ i ] = - ( ++d_nin );
  }

 d_b.resize( d_neq );
 d_hl.resize( d_nin );
 d_hu.resize( d_nin );
 for( int i = 0 ; i < numrows ; ++i )
  if( d_row_pos[ i ] >= 0 )
   d_b[ d_row_pos[ i ] ] = lo[ i ];
  else
   if( d_row_pos[ i ] != PIQP_free_row ) {
    d_hl[ - d_row_pos[ i ] - 1 ] = lo[ i ];
    d_hu[ - d_row_pos[ i ] - 1 ] = hi[ i ];
    }

 // the rows, from the columns of MILPSolver
 std::vector< Eigen::Triplet< double > > trpA , trpG;
 for( int j = 0 ; j < n ; ++j )
  for( int k = matbeg[ j ] ; k < matbeg[ j ] + matcnt[ j ] ; ++k ) {
   const int p = d_row_pos[ matind[ k ] ];
   if( p >= 0 )
    trpA.emplace_back( p , j , matval[ k ] );
   else
    if( p != PIQP_free_row )
     trpG.emplace_back( - p - 1 , j , matval[ k ] );
   }
 d_A.resize( d_neq , n );
 d_A.setFromTriplets( trpA.begin() , trpA.end() );
 d_G.resize( d_nin , n );
 d_G.setFromTriplets( trpG.begin() , trpG.end() );

 // the bounds, the box constraints of PIQP
 d_xl.resize( n );
 d_xu.resize( n );
 for( int j = 0 ; j < n ; ++j ) {
  d_xl[ j ] = lb[ j ];
  d_xu[ j ] = ub[ j ];
  }

 d_loaded = true;
 d_status = 0;

 }  // end( PIQPMILPSolver::load_problem )

/*--------------------------------------------------------------------------*/
/*-------------------- METHODS FOR SOLVING THE MODEL -----------------------*/
/*--------------------------------------------------------------------------*/

int PIQPMILPSolver::guts_of_compute( void )
{
 // the settings that are read by the setup of PIQP
 auto set_setup = [ this ]( piqp::Settings< double > & s ) {
  s.preconditioner_iter = d_precond_iter;
  s.preconditioner_scale_cost = d_precond_cost;
  s.verbose = false;
  s.compute_timings = false;
  };

 // the settings that are read by the solve
 auto set_solve = [ this ]( piqp::Settings< double > & s ) {
  s.eps_abs = d_eps_abs;
  s.eps_rel = d_eps_rel;
  s.check_duality_gap = d_check_gap;
  s.eps_duality_gap_abs = d_eps_gap_abs;
  s.eps_duality_gap_rel = d_eps_gap_rel;
  s.infeasibility_threshold = d_infeas_thr;
  s.rho_init = d_rho_init;
  s.delta_init = d_delta_init;
  s.iterative_refinement_always_enabled = d_iter_ref;
  s.max_iter = ( d_max_iter < Inf< int >() ) ? d_max_iter
                                             : PIQP_dflt.max_iter;
  };

 // a problem with no variables is solved here, PIQP needing some
 if( d_n == 0 ) {
  bool feas = true;
  for( int i = 0 ; i < d_neq ; ++i )
   feas &= ( d_b[ i ] == 0 );
  for( int i = 0 ; i < d_nin ; ++i )
   feas &= ( d_hl[ i ] <= 0 ) && ( d_hu[ i ] >= 0 );
  d_status = feas ? piqp::PIQP_SOLVED : piqp::PIQP_PRIMAL_INFEASIBLE;
  d_x.resize( 0 );
  d_y.setZero( d_neq );
  d_z.setZero( d_nin );
  d_zb.resize( 0 );
  d_backend = d_par_backend;
  d_fval = 0;
  d_loaded = false;
  return( sol_status = decode_status( d_status ) );
  }

 if( d_loaded || ( d_backend != d_par_backend ) || ( ! d_slv ) ) {
  free_solver();
  d_backend = d_par_backend;
  d_slv = std::make_unique< Solvers >();

  if( d_backend == 0 ) {
   auto & s = d_slv->dense;
   s = std::make_unique< piqp::DenseSolver< double > >();
   set_setup( s->settings() );
   set_solve( s->settings() );
   const Eigen::MatrixXd P( d_P ) , A( d_A ) , G( d_G );
   s->setup( P , d_c , A , d_b , G , d_hl , d_hu , d_xl , d_xu );
   }
  else {
   auto & s = d_slv->sparse;
   s = std::make_unique< piqp::SparseSolver< double , int > >();
   set_setup( s->settings() );
   set_solve( s->settings() );
   s->setup( d_P , d_c , d_A , d_b , d_G , d_hl , d_hu , d_xl , d_xu );
   }
  d_loaded = false;
  }

 piqp::Status st;
 const piqp::Result< double > * res;
 if( d_slv->dense ) {
  set_solve( d_slv->dense->settings() );
  st = d_slv->dense->solve();
  res = & d_slv->dense->result();
  }
 else {
  set_solve( d_slv->sparse->settings() );
  st = d_slv->sparse->solve();
  res = & d_slv->sparse->result();
  }

 d_status = int( st );

 // the value is computed here, from the solution, rather than taken from
 // PIQP
 d_fval = 0;
 if( solved() ) {
  d_x = res->x;
  d_y = res->y;
  d_z = res->z_u - res->z_l;
  d_zb = res->z_bu - res->z_bl;
  d_fval = d_x.dot( d_P * d_x ) / 2 + d_c.dot( d_x );
  }

 return( sol_status = decode_status( d_status ) );

 }  // end( PIQPMILPSolver::guts_of_compute )

/*--------------------------------------------------------------------------*/

int PIQPMILPSolver::decode_status( int status )
{
 switch( status ) {
  case( piqp::PIQP_SOLVED ):            return( kOK );
  case( piqp::PIQP_PRIMAL_INFEASIBLE ): return( kInfeasible );
  case( piqp::PIQP_DUAL_INFEASIBLE ):   return( kUnbounded );
  case( piqp::PIQP_MAX_ITER_REACHED ):  return( kStopIter );
  default:                              return( kError );
  }
 }

/*--------------------------------------------------------------------------*/

void PIQPMILPSolver::free_solver( void )
{
 d_slv.reset();
 }

/*--------------------------------------------------------------------------*/

bool PIQPMILPSolver::solved( void ) const
{
 return( d_status == piqp::PIQP_SOLVED );
 }

/*--------------------------------------------------------------------------*/
/*---------------------- METHODS FOR READING RESULTS -----------------------*/
/*--------------------------------------------------------------------------*/

Solver::OFValue PIQPMILPSolver::get_lb( void )
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

Solver::OFValue PIQPMILPSolver::get_ub( void )
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

bool PIQPMILPSolver::has_var_solution( void )
{
 return( solved() );
 }

/*--------------------------------------------------------------------------*/

Solver::OFValue PIQPMILPSolver::get_var_value( void )
{
 return( objsense == -1 ? get_lb() : get_ub() );
 }

/*--------------------------------------------------------------------------*/

void PIQPMILPSolver::get_var_solution( Configuration * solc )
{
 if( ! solved() )
  throw( std::logic_error( "PIQPMILPSolver::get_var_solution: "
			   "no solution available" ) );

 MILPSolver::write_var_solution(
  std::vector< double >( d_x.data() , d_x.data() + d_x.size() ) );
 }

/*--------------------------------------------------------------------------*/

bool PIQPMILPSolver::has_dual_solution( void )
{
 return( solved() );
 }

/*--------------------------------------------------------------------------*/

bool PIQPMILPSolver::is_dual_feasible( void )
{
 return( solved() );
 }

/*--------------------------------------------------------------------------*/

void PIQPMILPSolver::get_dual_solution( Configuration * solc )
{
 if( ! solved() )
  throw( std::logic_error( "PIQPMILPSolver::get_dual_solution: "
			   "no solution available" ) );

 // PIQP has P x + c + A' y + G' ( z_u - z_l ) + ( z_bu - z_bl ) = 0, with
 // z_l , z_u , z_bl , z_bu >= 0, for the problem it minimizes: the
 // multipliers of the problem of the Block are - y and - ( z_u - z_l ) for
 // a minimization and y and z_u - z_l for a maximization, whose Objective
 // PIQP has been given negated
 std::vector< double > pi( d_row_pos.size() , 0 );
 for( std::size_t i = 0 ; i < d_row_pos.size() ; ++i ) {
  const int p = d_row_pos[ i ];
  if( p != PIQP_free_row )
   pi[ i ] = - d_sign * ( p >= 0 ? d_y[ p ] : d_z[ - p - 1 ] );
  }

 std::vector< double > rc( d_n );
 for( int j = 0 ; j < d_n ; ++j )
  rc[ j ] = - d_sign * d_zb[ j ];

 MILPSolver::write_dual_solution( pi , rc );
 }

/*--------------------------------------------------------------------------*/
/*------------------ METHODS FOR HANDLING THE PARAMETERS -------------------*/
/*--------------------------------------------------------------------------*/

void PIQPMILPSolver::set_par( idx_type par , int value )
{
 switch( par ) {
  case( intMaxIter ):             d_max_iter = std::max( value , 1 ); return;
  case( intPIQPBackend ):         d_par_backend = value ? 1 : 0;      return;
  case( intPIQPCheckDualityGap ): d_check_gap = bool( value );        return;
  case( intPIQPIterRefinement ):  d_iter_ref = bool( value );         return;
  case( intPIQPPrecondIter ):  // the equilibration is computed when the
   value = std::max( value , 0 );  // solver is set up
   if( value != d_precond_iter ) {
    d_precond_iter = value;
    if( d_slv )
     d_loaded = true;
    }
   return;
  case( intPIQPPrecondScaleCost ):
   if( bool( value ) != d_precond_cost ) {
    d_precond_cost = bool( value );
    if( d_slv )
     d_loaded = true;
    }
   return;
  }

 MILPSolver::set_par( par , value );
 }

/*--------------------------------------------------------------------------*/

void PIQPMILPSolver::set_par( idx_type par , double value )
{
 switch( par ) {
  case( dblFAccSol ):             d_eps_abs = value;     return;
  case( dblPIQPEpsRel ):          d_eps_rel = value;     return;
  case( dblPIQPEpsGapAbs ):       d_eps_gap_abs = value; return;
  case( dblPIQPEpsGapRel ):       d_eps_gap_rel = value; return;
  case( dblPIQPInfeasThreshold ): d_infeas_thr = value;  return;
  case( dblPIQPRhoInit ):         d_rho_init = value;    return;
  case( dblPIQPDeltaInit ):       d_delta_init = value;  return;
  }

 MILPSolver::set_par( par , value );
 }

/*--------------------------------------------------------------------------*/

int PIQPMILPSolver::get_dflt_int_par( idx_type par ) const
{
 switch( par ) {
  case( intPIQPBackend ):          return( 0 );
  case( intPIQPCheckDualityGap ):  return( PIQP_dflt.check_duality_gap );
  case( intPIQPPrecondIter ):
   return( int( PIQP_dflt.preconditioner_iter ) );
  case( intPIQPPrecondScaleCost ):
   return( PIQP_dflt.preconditioner_scale_cost );
  case( intPIQPIterRefinement ):
   return( PIQP_dflt.iterative_refinement_always_enabled );
  }

 return( MILPSolver::get_dflt_int_par( par ) );
 }

/*--------------------------------------------------------------------------*/

double PIQPMILPSolver::get_dflt_dbl_par( idx_type par ) const
{
 switch( par ) {
  case( dblFAccSol ):             return( PIQP_dflt_eps_abs );
  case( dblPIQPEpsRel ):          return( PIQP_dflt.eps_rel );
  case( dblPIQPEpsGapAbs ):       return( PIQP_dflt_eps_abs );
  case( dblPIQPEpsGapRel ):       return( PIQP_dflt.eps_duality_gap_rel );
  case( dblPIQPInfeasThreshold ):
   return( PIQP_dflt.infeasibility_threshold );
  case( dblPIQPRhoInit ):         return( PIQP_dflt.rho_init );
  case( dblPIQPDeltaInit ):       return( PIQP_dflt.delta_init );
  }

 return( MILPSolver::get_dflt_dbl_par( par ) );
 }

/*--------------------------------------------------------------------------*/

int PIQPMILPSolver::get_int_par( idx_type par ) const
{
 switch( par ) {
  case( intMaxIter ):              return( d_max_iter );
  case( intPIQPBackend ):          return( d_par_backend );
  case( intPIQPCheckDualityGap ):  return( d_check_gap );
  case( intPIQPPrecondIter ):      return( d_precond_iter );
  case( intPIQPPrecondScaleCost ): return( d_precond_cost );
  case( intPIQPIterRefinement ):   return( d_iter_ref );
  }

 return( MILPSolver::get_int_par( par ) );
 }

/*--------------------------------------------------------------------------*/

double PIQPMILPSolver::get_dbl_par( idx_type par ) const
{
 switch( par ) {
  case( dblFAccSol ):             return( d_eps_abs );
  case( dblPIQPEpsRel ):          return( d_eps_rel );
  case( dblPIQPEpsGapAbs ):       return( d_eps_gap_abs );
  case( dblPIQPEpsGapRel ):       return( d_eps_gap_rel );
  case( dblPIQPInfeasThreshold ): return( d_infeas_thr );
  case( dblPIQPRhoInit ):         return( d_rho_init );
  case( dblPIQPDeltaInit ):       return( d_delta_init );
  }

 return( MILPSolver::get_dbl_par( par ) );
 }

/*--------------------------------------------------------------------------*/

Solver::idx_type PIQPMILPSolver::int_par_str2idx( const std::string & name )
 const
{
 auto it = std::find( PIQP_int_pars_str.begin() ,
		      PIQP_int_pars_str.end() , name );
 if( it != PIQP_int_pars_str.end() )
  return( intLastAlgParMILP + ( it - PIQP_int_pars_str.begin() ) );

 return( MILPSolver::int_par_str2idx( name ) );
 }

/*--------------------------------------------------------------------------*/

const std::string & PIQPMILPSolver::int_par_idx2str( idx_type idx ) const
{
 if( ( idx >= intLastAlgParMILP ) && ( idx < intLastAlgParPIQP ) )
  return( PIQP_int_pars_str[ idx - intLastAlgParMILP ] );

 return( MILPSolver::int_par_idx2str( idx ) );
 }

/*--------------------------------------------------------------------------*/

Solver::idx_type PIQPMILPSolver::dbl_par_str2idx( const std::string & name )
 const
{
 auto it = std::find( PIQP_dbl_pars_str.begin() ,
		      PIQP_dbl_pars_str.end() , name );
 if( it != PIQP_dbl_pars_str.end() )
  return( dblLastAlgParMILP + ( it - PIQP_dbl_pars_str.begin() ) );

 return( MILPSolver::dbl_par_str2idx( name ) );
 }

/*--------------------------------------------------------------------------*/

const std::string & PIQPMILPSolver::dbl_par_idx2str( idx_type idx ) const
{
 if( ( idx >= dblLastAlgParMILP ) && ( idx < dblLastAlgParPIQP ) )
  return( PIQP_dbl_pars_str[ idx - dblLastAlgParMILP ] );

 return( MILPSolver::dbl_par_idx2str( idx ) );
 }

/*--------------------------------------------------------------------------*/
/*---------------------- End File PIQPMILPSolver.cpp -----------------------*/
/*--------------------------------------------------------------------------*/
