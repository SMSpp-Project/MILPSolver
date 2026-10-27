/*--------------------------------------------------------------------------*/
/*------------------------ File ProxQPMILPSolver.cpp -----------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Implementation of the ProxQPMILPSolver class, which implements the
 * MILPSolver interface for convex Quadratic Programs with linear
 * constraints using the solver ProxQP of the ProxSuite library.
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
#include <stdexcept>

#include <proxsuite/proxqp/dense/dense.hpp>
#include <proxsuite/proxqp/sparse/sparse.hpp>

#include "ProxQPMILPSolver.h"

/*--------------------------------------------------------------------------*/
/*-------------------------------- USING -----------------------------------*/
/*--------------------------------------------------------------------------*/

using namespace SMSpp_di_unipi_it;

namespace pqp = proxsuite::proxqp;

/*--------------------------------------------------------------------------*/
/*-------------------------------- CONSTANTS -------------------------------*/
/*--------------------------------------------------------------------------*/

static const std::array< std::string , 5 > ProxQP_int_pars_str = {
 "intProxQPWarmStart" , "intProxQPInitialGuess" , "intProxQPBackend" ,
 "intProxQPDenseBackend" , "intProxQPCheckDualityGap" };

static const std::array< std::string , 6 > ProxQP_dbl_pars_str = {
 "dblProxQPEpsRel" , "dblProxQPEpsPrimalInf" , "dblProxQPEpsDualInf" ,
 "dblProxQPRho" , "dblProxQPMuEq" , "dblProxQPMuIn" };

// the settings of ProxQP for the PrimalDualLDLT factorization, which give
// the defaults of the parameters but those below

static const pqp::Settings< double > ProxQP_dflt;

// the default of dblFAccSol, tighter than the 1e-5 of ProxQP, which leaves
// errors of that order in the multipliers

static constexpr double ProxQP_dflt_eps_abs = 1e-9;

// the default of dblProxQPEpsPrimalInf and dblProxQPEpsDualInf, tighter than
// the 1e-4 of ProxQP, with which feasible and bounded LPs (e.g., one with
// costs of 1000) are declared primal infeasible

static constexpr double ProxQP_dflt_eps_inf = 1e-8;

/*--------------------------------------------------------------------------*/
/*------------------------- FACTORY MANAGEMENT -----------------------------*/
/*--------------------------------------------------------------------------*/

SMSpp_insert_in_factory_cpp_0( ProxQPMILPSolver );

/*--------------------------------------------------------------------------*/
/*--------------------- CONSTRUCTOR AND DESTRUCTOR -------------------------*/
/*--------------------------------------------------------------------------*/

ProxQPMILPSolver::ProxQPMILPSolver( void ) : MILPSolver()
{
 d_eps_abs = ProxQP_dflt_eps_abs;
 d_eps_rel = ProxQP_dflt.eps_rel;
 d_eps_primal_inf = ProxQP_dflt_eps_inf;
 d_eps_dual_inf = ProxQP_dflt_eps_inf;
 d_rho = ProxQP_dflt.default_rho;
 d_mu_eq = ProxQP_dflt.default_mu_eq;
 d_mu_in = ProxQP_dflt.default_mu_in;
 }

/*--------------------------------------------------------------------------*/

ProxQPMILPSolver::~ProxQPMILPSolver()
{
 free_solver();
 }

/*--------------------------------------------------------------------------*/
/*--------------------- DERIVED METHODS OF BASE CLASS ----------------------*/
/*--------------------------------------------------------------------------*/

void ProxQPMILPSolver::clear_problem( unsigned int what )
{
 MILPSolver::clear_problem( what );
 free_solver();
 d_x.resize( 0 );
 d_y.resize( 0 );
 d_z.resize( 0 );
 d_ws_backend = -1;
 d_n = d_neq = d_nin = 0;
 d_row_pos.clear();
 d_box_pos.clear();
 d_loaded = false;
 d_status = -1;
 }

/*--------------------------------------------------------------------------*/

void ProxQPMILPSolver::load_problem( void )
{
 MILPSolver::load_problem();
 free_solver();

 if( numquadrows > 0 )
  throw( std::invalid_argument( "ProxQPMILPSolver::load_problem: "
				"quadratic constraints are not supported" ) );

 if( ( int_vars > 0 ) && ( relax_int_vars == 0 ) )
  throw( std::invalid_argument( "ProxQPMILPSolver::load_problem: integer "
				"variables are only supported if relaxed "
				"(intRelaxIntVars)" ) );

 const int n = numcols;
 d_n = n;
 d_sign = ( objsense == -1 ) ? -1 : 1;  // ProxQP minimizes

 // the objective: c' x + sum_j q_j x_j^2 + sum_k ndq_k x_i x_j is
 // 1/2 x' H x + g' x with H_jj = 2 q_j and H_ij = H_ji = ndq_k
 d_g.resize( n );
 for( int j = 0 ; j < n ; ++j )
  d_g[ j ] = d_sign * objective[ j ];

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
 d_H.resize( n , n );
 d_H.setFromTriplets( trp.begin() , trp.end() );
 d_H.prune( 0.0 );

 // the sides of the rows: those whose two sides coincide are equalities,
 // which go in A, and all the others in C
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
    throw( std::invalid_argument( "ProxQPMILPSolver::load_problem: row " +
				  std::to_string( i ) + " has sense " +
				  std::string( 1 , sense[ i ] ) ) );
   }
  if( lo[ i ] == hi[ i ] )
   d_row_pos[ i ] = d_neq++;
  else
   d_row_pos[ i ] = - ( ++d_nin );
  }

 d_b.resize( d_neq );
 d_l.resize( d_nin );
 d_u.resize( d_nin );
 for( int i = 0 ; i < numrows ; ++i )
  if( d_row_pos[ i ] >= 0 )
   d_b[ d_row_pos[ i ] ] = lo[ i ];
  else {
   d_l[ - d_row_pos[ i ] - 1 ] = lo[ i ];
   d_u[ - d_row_pos[ i ] - 1 ] = hi[ i ];
   }

 // the rows, from the columns of MILPSolver
 std::vector< Eigen::Triplet< double > > trpA , trpC;
 for( int j = 0 ; j < n ; ++j )
  for( int k = matbeg[ j ] ; k < matbeg[ j ] + matcnt[ j ] ; ++k ) {
   const int p = d_row_pos[ matind[ k ] ];
   if( p >= 0 )
    trpA.emplace_back( p , j , matval[ k ] );
   else
    trpC.emplace_back( - p - 1 , j , matval[ k ] );
   }
 d_A.resize( d_neq , n );
 d_A.setFromTriplets( trpA.begin() , trpA.end() );
 d_C.resize( d_nin , n );
 d_C.setFromTriplets( trpC.begin() , trpC.end() );

 // the bounds
 d_lbox.resize( n );
 d_ubox.resize( n );
 for( int j = 0 ; j < n ; ++j ) {
  d_lbox[ j ] = lb[ j ];
  d_ubox[ j ] = ub[ j ];
  }

 d_loaded = true;
 d_status = -1;

 }  // end( ProxQPMILPSolver::load_problem )

/*--------------------------------------------------------------------------*/
/*-------------------- METHODS FOR SOLVING THE MODEL -----------------------*/
/*--------------------------------------------------------------------------*/

int ProxQPMILPSolver::guts_of_compute( void )
{
 // the settings that may change between two solves
 auto set = [ this ]( pqp::Settings< double > & s ) {
  s.eps_abs = d_eps_abs;
  s.eps_rel = d_eps_rel;
  s.max_iter = ( d_max_iter < Inf< int >() ) ? d_max_iter
                                             : ProxQP_dflt.max_iter;
  s.eps_primal_inf = d_eps_primal_inf;
  s.eps_dual_inf = d_eps_dual_inf;
  s.check_duality_gap = d_check_gap;
  s.verbose = false;
  };

 // a problem with no variables is solved here, ProxQP needing some
 if( d_n == 0 ) {
  bool feas = true;
  for( int i = 0 ; i < d_neq ; ++i )
   feas &= ( d_b[ i ] == 0 );
  for( int i = 0 ; i < d_nin ; ++i )
   feas &= ( d_l[ i ] <= 0 ) && ( d_u[ i ] >= 0 );
  d_status = int( feas ? pqp::QPSolverOutput::PROXQP_SOLVED
                       : pqp::QPSolverOutput::PROXQP_PRIMAL_INFEASIBLE );
  d_x.resize( 0 );
  d_y.setZero( d_neq );
  d_z.setZero( d_nin );
  d_ws_backend = d_backend = d_par_backend;
  d_fval = 0;
  d_loaded = false;
  return( sol_status = decode_status( d_status ) );
  }

 if( d_loaded || ( d_backend != d_par_backend ) ||
     ( ! ( d_dqp || d_sqp ) ) ) {
  free_solver();
  d_backend = d_par_backend;

  if( d_backend == 0 ) {  // dense, with the bounds as box constraints
   const auto ht = d_H.nonZeros() ? pqp::HessianType::Dense
                                  : pqp::HessianType::Zero;
   const auto db = d_dense_backend == 1 ? pqp::DenseBackend::PrimalDualLDLT
                 : d_dense_backend == 2 ? pqp::DenseBackend::PrimalLDLT
                                        : pqp::DenseBackend::Automatic;
   d_dqp = std::make_unique< DQP >( d_n , d_neq , d_nin , true , db , ht );
   set( d_dqp->settings );
   const Eigen::MatrixXd H( d_H ) , A( d_A ) , C( d_C );
   d_dqp->init( H , d_g , A , d_b , C , d_l , d_u , d_lbox , d_ubox ,
		true , d_rho , d_mu_eq , d_mu_in );
   }
  else {                  // sparse, with the finite bounds as rows of C
   d_box_pos.assign( d_n , -1 );
   std::vector< double > lo , hi;
   std::vector< Eigen::Triplet< double > > trp;
   for( int k = 0 ; k < d_C.outerSize() ; ++k )
    for( SpMat::InnerIterator it( d_C , k ) ; it ; ++it )
     trp.emplace_back( it.row() , it.col() , it.value() );
   int r = d_nin;
   for( int j = 0 ; j < d_n ; ++j )
    if( ( d_lbox[ j ] > - Inf< double >() ) ||
	( d_ubox[ j ] < Inf< double >() ) ) {
     trp.emplace_back( r , j , 1.0 );
     lo.push_back( d_lbox[ j ] );
     hi.push_back( d_ubox[ j ] );
     d_box_pos[ j ] = r++;
     }
   SpMat C( r , d_n );
   C.setFromTriplets( trp.begin() , trp.end() );
   Eigen::VectorXd l( r ) , u( r );
   l.head( d_nin ) = d_l;
   u.head( d_nin ) = d_u;
   for( int i = d_nin ; i < r ; ++i ) {
    l[ i ] = lo[ i - d_nin ];
    u[ i ] = hi[ i - d_nin ];
    }

   d_sqp = std::make_unique< SQP >( d_n , d_neq , r );
   set( d_sqp->settings );
   d_sqp->init( d_H , d_g , d_A , d_b , C , l , u ,
		true , d_rho , d_mu_eq , d_mu_in );
   }
  d_loaded = false;
  }

 pqp::Settings< double > & s = d_dqp ? d_dqp->settings : d_sqp->settings;
 pqp::Results< double > & res = d_dqp ? d_dqp->results : d_sqp->results;
 set( s );

 // the start: the solution of the previous solve, if of the same size and
 // of the same backend, otherwise the initial guess
 if( d_warm_start && ( d_ws_backend == d_backend ) &&
     ( d_x.size() == res.x.size() ) && ( d_y.size() == res.y.size() ) &&
     ( d_z.size() == res.z.size() ) ) {
  s.initial_guess = pqp::InitialGuessStatus::WARM_START;
  res.x = d_x;
  res.y = d_y;
  res.z = d_z;
  }
 else
  s.initial_guess = d_initial_guess ?
   pqp::InitialGuessStatus::EQUALITY_CONSTRAINED_INITIAL_GUESS :
   pqp::InitialGuessStatus::NO_INITIAL_GUESS;

 if( d_dqp )
  d_dqp->solve();
 else
  d_sqp->solve();

 d_status = int( res.info.status );

 // the ray of a dual infeasible problem, which ProxQP leaves in place of
 // the solution
 d_cx.resize( 0 );
 if( d_status == int( pqp::QPSolverOutput::PROXQP_DUAL_INFEASIBLE ) )
  d_cx = res.x;

 // the value is computed here, from the solution, rather than taken from
 // ProxQP
 d_fval = 0;
 if( solved() ) {
  d_x = res.x;
  d_y = res.y;
  d_z = res.z;
  d_ws_backend = d_backend;
  d_fval = d_x.dot( d_H * d_x ) / 2 + d_g.dot( d_x );
  }

 return( sol_status = decode_status( d_status ) );

 }  // end( ProxQPMILPSolver::guts_of_compute )

/*--------------------------------------------------------------------------*/

int ProxQPMILPSolver::decode_status( int status )
{
 switch( pqp::QPSolverOutput( status ) ) {
  case( pqp::QPSolverOutput::PROXQP_SOLVED ):            return( kOK );
  case( pqp::QPSolverOutput::PROXQP_PRIMAL_INFEASIBLE ):
  case( pqp::QPSolverOutput::PROXQP_SOLVED_CLOSEST_PRIMAL_FEASIBLE ):
                                                        return( kInfeasible );
  case( pqp::QPSolverOutput::PROXQP_DUAL_INFEASIBLE ):   return( kUnbounded );
  case( pqp::QPSolverOutput::PROXQP_MAX_ITER_REACHED ):  return( kStopIter );
  default:                                              return( kError );
  }
 }

/*--------------------------------------------------------------------------*/

void ProxQPMILPSolver::free_solver( void )
{
 d_dqp.reset();
 d_sqp.reset();
 }

/*--------------------------------------------------------------------------*/

bool ProxQPMILPSolver::solved( void ) const
{
 return( ( d_status >= 0 ) && ( decode_status( d_status ) == kOK ) );
 }

/*--------------------------------------------------------------------------*/
/*---------------------- METHODS FOR READING RESULTS -----------------------*/
/*--------------------------------------------------------------------------*/

Solver::OFValue ProxQPMILPSolver::get_lb( void )
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

Solver::OFValue ProxQPMILPSolver::get_ub( void )
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

bool ProxQPMILPSolver::has_var_solution( void )
{
 return( solved() );
 }

/*--------------------------------------------------------------------------*/

Solver::OFValue ProxQPMILPSolver::get_var_value( void )
{
 return( objsense == -1 ? get_lb() : get_ub() );
 }

/*--------------------------------------------------------------------------*/

void ProxQPMILPSolver::get_var_solution( Configuration * solc )
{
 if( ! solved() )
  throw( std::logic_error( "ProxQPMILPSolver::get_var_solution: "
			   "no solution available" ) );

 MILPSolver::write_var_solution(
  std::vector< double >( d_x.data() , d_x.data() + d_x.size() ) , solc );
 }

/*--------------------------------------------------------------------------*/

bool ProxQPMILPSolver::has_dual_solution( void )
{
 return( solved() );
 }

/*--------------------------------------------------------------------------*/

bool ProxQPMILPSolver::is_dual_feasible( void )
{
 return( solved() );
 }

/*--------------------------------------------------------------------------*/

void ProxQPMILPSolver::get_dual_solution( Configuration * solc )
{
 if( ! solved() )
  throw( std::logic_error( "ProxQPMILPSolver::get_dual_solution: "
			   "no solution available" ) );

 // ProxQP has H x + g + A' y + C' z = 0, with z > 0 on an active upper side
 // and z < 0 on an active lower one, for the problem it minimizes: the
 // multipliers of the problem of the Block are - y and - z for a
 // minimization and y and z for a maximization, whose Objective ProxQP has
 // been given negated
 std::vector< double > pi( d_row_pos.size() );
 for( std::size_t i = 0 ; i < d_row_pos.size() ; ++i ) {
  const int p = d_row_pos[ i ];
  pi[ i ] = - d_sign * ( p >= 0 ? d_y[ p ] : d_z[ - p - 1 ] );
  }

 // the multipliers of the bounds follow those of the rows in z: all of
 // them with the dense backend, those of the finite ones with the sparse
 std::vector< double > rc( d_n , 0 );
 for( int j = 0 ; j < d_n ; ++j )
  if( d_backend == 0 )
   rc[ j ] = - d_sign * d_z[ d_nin + j ];
  else
   if( d_box_pos[ j ] >= 0 )
    rc[ j ] = - d_sign * d_z[ d_box_pos[ j ] ];

 MILPSolver::write_dual_solution( pi , rc , solc );
 }

/*--------------------------------------------------------------------------*/

bool ProxQPMILPSolver::has_var_direction( void )
{
 return( ( d_status == int( pqp::QPSolverOutput::PROXQP_DUAL_INFEASIBLE ) )
	 && ( d_cx.size() == d_n ) );
 }

/*--------------------------------------------------------------------------*/

void ProxQPMILPSolver::get_var_direction( Configuration * dirc )
{
 if( ! has_var_direction() )
  throw( std::logic_error( "ProxQPMILPSolver::get_var_direction: "
			   "no unbounded direction available" ) );

 // H x = 0 and g' x < 0 for the problem ProxQP minimizes, whose g is that
 // of the Block negated for a maximization: x is a ray of the Block
 MILPSolver::write_var_solution(
  std::vector< double >( d_cx.begin() , d_cx.end() ) , dirc );
 }

/*--------------------------------------------------------------------------*/
/*------------------ METHODS FOR HANDLING THE PARAMETERS -------------------*/
/*--------------------------------------------------------------------------*/

void ProxQPMILPSolver::set_par( idx_type par , int value )
{
 switch( par ) {
  case( intMaxIter ): d_max_iter = std::max( value , 0 );   return;
  case( intProxQPWarmStart ):    d_warm_start = bool( value ); return;
  case( intProxQPInitialGuess ): d_initial_guess = bool( value ); return;
  case( intProxQPBackend ):      d_par_backend = value ? 1 : 0; return;
  case( intProxQPDenseBackend ):
   if( value != d_dense_backend ) {
    d_dense_backend = std::clamp( value , 0 , 2 );
    if( d_dqp )  // the factorization is chosen when the solver is set up
     d_loaded = true;
    }
   return;
  case( intProxQPCheckDualityGap ): d_check_gap = bool( value ); return;
  }

 MILPSolver::set_par( par , value );
 }

/*--------------------------------------------------------------------------*/

void ProxQPMILPSolver::set_par( idx_type par , double value )
{
 switch( par ) {
  case( dblFAccSol ):            d_eps_abs = value;        return;
  case( dblProxQPEpsRel ):       d_eps_rel = value;        return;
  case( dblProxQPEpsPrimalInf ): d_eps_primal_inf = value; return;
  case( dblProxQPEpsDualInf ):   d_eps_dual_inf = value;   return;
  case( dblProxQPRho ):          // the proximal parameters are given to
  case( dblProxQPMuEq ):         // the solver when it is set up
  case( dblProxQPMuIn ): {
   double & p = par == dblProxQPRho ? d_rho
              : par == dblProxQPMuEq ? d_mu_eq : d_mu_in;
   if( value != p ) {
    p = value;
    if( d_dqp || d_sqp )
     d_loaded = true;
    }
   return;
   }
  }

 MILPSolver::set_par( par , value );
 }

/*--------------------------------------------------------------------------*/

int ProxQPMILPSolver::get_dflt_int_par( idx_type par ) const
{
 switch( par ) {
  case( intProxQPWarmStart ):       return( 1 );
  case( intProxQPInitialGuess ):    return( 1 );
  case( intProxQPBackend ):         return( 0 );
  case( intProxQPDenseBackend ):    return( 0 );
  case( intProxQPCheckDualityGap ): return( ProxQP_dflt.check_duality_gap );
  }

 return( MILPSolver::get_dflt_int_par( par ) );
 }

/*--------------------------------------------------------------------------*/

double ProxQPMILPSolver::get_dflt_dbl_par( idx_type par ) const
{
 switch( par ) {
  case( dblFAccSol ):            return( ProxQP_dflt_eps_abs );
  case( dblProxQPEpsRel ):       return( ProxQP_dflt.eps_rel );
  case( dblProxQPEpsPrimalInf ): return( ProxQP_dflt_eps_inf );
  case( dblProxQPEpsDualInf ):   return( ProxQP_dflt_eps_inf );
  case( dblProxQPRho ):          return( ProxQP_dflt.default_rho );
  case( dblProxQPMuEq ):         return( ProxQP_dflt.default_mu_eq );
  case( dblProxQPMuIn ):         return( ProxQP_dflt.default_mu_in );
  }

 return( MILPSolver::get_dflt_dbl_par( par ) );
 }

/*--------------------------------------------------------------------------*/

int ProxQPMILPSolver::get_int_par( idx_type par ) const
{
 switch( par ) {
  case( intMaxIter ):               return( d_max_iter );
  case( intProxQPWarmStart ):       return( d_warm_start );
  case( intProxQPInitialGuess ):    return( d_initial_guess );
  case( intProxQPBackend ):         return( d_par_backend );
  case( intProxQPDenseBackend ):    return( d_dense_backend );
  case( intProxQPCheckDualityGap ): return( d_check_gap );
  }

 return( MILPSolver::get_int_par( par ) );
 }

/*--------------------------------------------------------------------------*/

double ProxQPMILPSolver::get_dbl_par( idx_type par ) const
{
 switch( par ) {
  case( dblFAccSol ):            return( d_eps_abs );
  case( dblProxQPEpsRel ):       return( d_eps_rel );
  case( dblProxQPEpsPrimalInf ): return( d_eps_primal_inf );
  case( dblProxQPEpsDualInf ):   return( d_eps_dual_inf );
  case( dblProxQPRho ):          return( d_rho );
  case( dblProxQPMuEq ):         return( d_mu_eq );
  case( dblProxQPMuIn ):         return( d_mu_in );
  }

 return( MILPSolver::get_dbl_par( par ) );
 }

/*--------------------------------------------------------------------------*/

Solver::idx_type ProxQPMILPSolver::int_par_str2idx( const std::string & name )
 const
{
 auto it = std::find( ProxQP_int_pars_str.begin() ,
		      ProxQP_int_pars_str.end() , name );
 if( it != ProxQP_int_pars_str.end() )
  return( intLastAlgParMILP + ( it - ProxQP_int_pars_str.begin() ) );

 return( MILPSolver::int_par_str2idx( name ) );
 }

/*--------------------------------------------------------------------------*/

const std::string & ProxQPMILPSolver::int_par_idx2str( idx_type idx ) const
{
 if( ( idx >= intLastAlgParMILP ) && ( idx < intLastAlgParProxQP ) )
  return( ProxQP_int_pars_str[ idx - intLastAlgParMILP ] );

 return( MILPSolver::int_par_idx2str( idx ) );
 }

/*--------------------------------------------------------------------------*/

Solver::idx_type ProxQPMILPSolver::dbl_par_str2idx( const std::string & name )
 const
{
 auto it = std::find( ProxQP_dbl_pars_str.begin() ,
		      ProxQP_dbl_pars_str.end() , name );
 if( it != ProxQP_dbl_pars_str.end() )
  return( dblLastAlgParMILP + ( it - ProxQP_dbl_pars_str.begin() ) );

 return( MILPSolver::dbl_par_str2idx( name ) );
 }

/*--------------------------------------------------------------------------*/

const std::string & ProxQPMILPSolver::dbl_par_idx2str( idx_type idx ) const
{
 if( ( idx >= dblLastAlgParMILP ) && ( idx < dblLastAlgParProxQP ) )
  return( ProxQP_dbl_pars_str[ idx - dblLastAlgParMILP ] );

 return( MILPSolver::dbl_par_idx2str( idx ) );
 }

/*--------------------------------------------------------------------------*/
/*--------------------- End File ProxQPMILPSolver.cpp ----------------------*/
/*--------------------------------------------------------------------------*/
