/*--------------------------------------------------------------------------*/
/*------------------------ File HPIPMMILPSolver.cpp ------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Implementation of the HPIPMMILPSolver class, which implements the
 * MILPSolver interface for convex Quadratic Programs with linear constraints
 * using the dense interior-point solver of HPIPM.
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
#include <array>
#include <cmath>
#include <stdexcept>

#include "HPIPMMILPSolver.h"

#include <hpipm_d_dense_qp.h>
#include <hpipm_d_dense_qp_dim.h>
#include <hpipm_d_dense_qp_ipm.h>
#include <hpipm_d_dense_qp_sol.h>

/*--------------------------------------------------------------------------*/
/*-------------------------------- USING -----------------------------------*/
/*--------------------------------------------------------------------------*/

using namespace SMSpp_di_unipi_it;

/*--------------------------------------------------------------------------*/
/*-------------------------------- CONSTANTS -------------------------------*/
/*--------------------------------------------------------------------------*/

static const std::array< std::string , 2 > HPIPM_int_pars_str = {
 "intHPIPMMode" , "intHPIPMWarmStart" };

static constexpr double HPIPM_dflt_tol = 1e-8;

/*--------------------------------------------------------------------------*/
/*----------------------------- HPIPMData ----------------------------------*/
/*--------------------------------------------------------------------------*/
// the structures of HPIPM live in memory that HPIPM is given; they are kept
// from a solve to the next one, so that a warm start finds the previous
// solution in sol, as long as the sizes have not changed

struct HPIPMMILPSolver::HPIPMData {
 d_dense_qp_dim dim;
 d_dense_qp qp;
 d_dense_qp_sol sol;
 d_dense_qp_ipm_arg arg;
 d_dense_qp_ipm_ws ws;

 std::vector< char > dim_mem , qp_mem , sol_mem , arg_mem , ws_mem;

 int nv = -1 , ne = -1 , nb = -1 , ng = -1;  ///< the sizes of the memory
 int mode = -1;                              ///< the mode of the memory
 bool has_sol = false;           ///< true if sol holds a solution

 // the problem in the form of HPIPM
 std::vector< double > A , b , C , lg , ug , lg_mask , ug_mask;
 std::vector< double > lb , ub , lb_mask , ub_mask;
 std::vector< int > idxb;
 std::vector< int > eq_of;       ///< variable or row of each equality
 std::vector< int > g_of;        ///< row of each general constraint

 /// (re)creates the structures of HPIPM for the given sizes and mode
 void create( int n_v , int n_e , int n_b , int n_g , int md ) {
  if( ( n_v == nv ) && ( n_e == ne ) && ( n_b == nb ) && ( n_g == ng ) &&
      ( md == mode ) )
   return;

  nv = n_v; ne = n_e; nb = n_b; ng = n_g; mode = md;

  dim_mem.resize( d_dense_qp_dim_memsize() );
  d_dense_qp_dim_create( &dim , dim_mem.data() );
  d_dense_qp_dim_set_nv( n_v , &dim );
  d_dense_qp_dim_set_ne( n_e , &dim );
  d_dense_qp_dim_set_nb( n_b , &dim );
  d_dense_qp_dim_set_ng( n_g , &dim );

  qp_mem.resize( d_dense_qp_memsize( &dim ) );
  d_dense_qp_create( &dim , &qp , qp_mem.data() );

  sol_mem.resize( d_dense_qp_sol_memsize( &dim ) );
  d_dense_qp_sol_create( &dim , &sol , sol_mem.data() );

  arg_mem.resize( d_dense_qp_ipm_arg_memsize( &dim ) );
  d_dense_qp_ipm_arg_create( &dim , &arg , arg_mem.data() );
  set_default( md );

  // the workspace is sized on the settings, hence it is created after them
  ws_mem.resize( d_dense_qp_ipm_ws_memsize( &dim , &arg ) );
  d_dense_qp_ipm_ws_create( &dim , &arg , &ws , ws_mem.data() );

  has_sol = false;
  }

 /// the settings of the mode, the removal of the dependent equalities of a
 /// Block (e.g., a row equal to the bound of a fixed variable) and, for a
 /// problem with no inequality, the null-space factorization of its KKT
 /// system: HPIPM then solves it once, and the Schur complement of H, its
 /// default, gives a solution that violates the equalities when H is
 /// singular, which it declares optimal
 void set_default( int md ) {
  d_dense_qp_ipm_arg_set_default( hpipm_mode( md ) , &arg );
  int rm = 1;
  d_dense_qp_ipm_arg_set_remove_lin_dep_eq( &rm , &arg );
  int kkt = ( nb + ng == 0 ) ? 0 : 1;
  d_dense_qp_ipm_arg_set_kkt_fact_alg( &kkt , &arg );
  }

 void reset( void ) { nv = ne = nb = ng = mode = -1; has_sol = false; }
 };

/*--------------------------------------------------------------------------*/
/*------------------------- FACTORY MANAGEMENT -----------------------------*/
/*--------------------------------------------------------------------------*/

SMSpp_insert_in_factory_cpp_0( HPIPMMILPSolver );

/*--------------------------------------------------------------------------*/
/*--------------------- CONSTRUCTOR AND DESTRUCTOR -------------------------*/
/*--------------------------------------------------------------------------*/

HPIPMMILPSolver::HPIPMMILPSolver( void ) : MILPSolver() ,
 d_hp( std::make_unique< HPIPMData >() ) , d_tol( HPIPM_dflt_tol ) {}

/*--------------------------------------------------------------------------*/

HPIPMMILPSolver::~HPIPMMILPSolver() = default;

/*--------------------------------------------------------------------------*/
/*--------------------- DERIVED METHODS OF BASE CLASS ----------------------*/
/*--------------------------------------------------------------------------*/

void HPIPMMILPSolver::clear_problem( unsigned int what )
{
 MILPSolver::clear_problem( what );
 d_hp->reset();
 d_n = d_nc = 0;
 d_status = kError;
 }

/*--------------------------------------------------------------------------*/

void HPIPMMILPSolver::load_problem( void )
{
 // the solution of the previous problem, which clear_problem() does not
 // keep in the memory of HPIPM, is moved to the columns that are still
 // there [see previous_column_index()], the new ones at 0, for the warm
 // start of the next solve
 auto px = std::move( d_guess_x );
 const bool had = d_has_guess;

 MILPSolver::load_problem();

 if( numquadrows > 0 )
  throw( std::invalid_argument( "HPIPMMILPSolver::load_problem: "
				"quadratic constraints are not supported" ) );

 if( ( int_vars > 0 ) && ( relax_int_vars == 0 ) )
  throw( std::invalid_argument( "HPIPMMILPSolver::load_problem: integer "
				"variables are only supported if relaxed "
				"(intRelaxIntVars)" ) );

 const int n = numcols;
 d_n = n;
 d_nc = numrows;
 d_sign = ( objsense == -1 ) ? -1 : 1;  // HPIPM minimizes

 // the objective: c' x + sum_j q_j x_j^2 + sum_k ndq_k x_i x_j is
 // 1/2 x' H x + h' x with H_jj = 2 q_j and H_ij = H_ji = ndq_k
 d_h.assign( n , 0 );
 for( int j = 0 ; j < n ; ++j )
  d_h[ j ] = d_sign * objective[ j ];

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
 d_A.assign( std::size_t( d_nc ) * n , 0 );
 for( int j = 0 ; j < n ; ++j )
  for( int k = matbeg[ j ] ; k < matbeg[ j ] + matcnt[ j ] ; ++k )
   d_A[ std::size_t( matind[ k ] ) * n + j ] += matval[ k ];

 // the bounds of the variables, then the sides of the rows
 d_lb.assign( lb.begin() , lb.begin() + n );
 d_ub.assign( ub.begin() , ub.begin() + n );

 d_lbA.resize( d_nc );
 d_ubA.resize( d_nc );
 for( int i = 0 ; i < d_nc ; ++i )
  switch( sense[ i ] ) {
   case( 'L' ): d_lbA[ i ] = - Inf< double >();
                d_ubA[ i ] = rhs[ i ];
                break;
   case( 'G' ): d_lbA[ i ] = rhs[ i ];
                d_ubA[ i ] = Inf< double >();
                break;
   case( 'E' ): d_lbA[ i ] = d_ubA[ i ] = rhs[ i ];
                break;
   case( 'R' ): if( rngval[ i ] > 0 ) {
                 d_lbA[ i ] = rhs[ i ];
                 d_ubA[ i ] = rhs[ i ] + rngval[ i ];
                 }
                else {
                 d_lbA[ i ] = rhs[ i ] + rngval[ i ];
                 d_ubA[ i ] = rhs[ i ];
                 }
                break;
   default:
    throw( std::invalid_argument( "HPIPMMILPSolver::load_problem: row " +
				  std::to_string( i ) + " has sense " +
				  std::string( 1 , sense[ i ] ) ) );
   }

 d_x.assign( n , 0 );

 d_has_guess = false;
 d_guess_x.clear();
 if( had && ( px.size() == std::size_t( previous_sizes().first ) ) ) {
  d_guess_x.assign( n , 0 );
  for( int j = 0 ; j < n ; ++j )
   if( const int pj = previous_column_index( j ) ; pj >= 0 )
    d_guess_x[ j ] = px[ pj ];
  d_has_guess = true;
  }
 d_y.assign( n + d_nc , 0 );
 d_status = kError;

 }  // end( HPIPMMILPSolver::load_problem )

/*--------------------------------------------------------------------------*/
/*-------------------- METHODS FOR SOLVING THE MODEL -----------------------*/
/*--------------------------------------------------------------------------*/

int HPIPMMILPSolver::guts_of_compute( void )
{
 d_fval = 0;

 // with no variable, the problem is feasible if and only if each row admits
 // 0
 if( d_n == 0 ) {
  d_status = kOK;
  for( int i = 0 ; i < d_nc ; ++i )
   if( ( d_lbA[ i ] > d_tol ) || ( d_ubA[ i ] < - d_tol ) )
    d_status = kInfeasible;
  d_y.assign( d_nc , 0 );
  return( sol_status = d_status );
  }

 auto & hp = *d_hp;
 const int n = d_n;

 // the form of HPIPM: the variables with equal bounds and the rows with
 // equal sides are equalities, the bounds of the other variables are box
 // constraints and the other rows general ones, an infinite side being
 // masked out (a row with both sides infinite is left out)
 hp.eq_of.clear();
 hp.idxb.clear();
 hp.g_of.clear();
 for( int j = 0 ; j < n ; ++j )
  if( d_lb[ j ] == d_ub[ j ] )
   hp.eq_of.push_back( j );
  else
   if( ( d_lb[ j ] > - Inf< double >() ) || ( d_ub[ j ] < Inf< double >() ) )
    hp.idxb.push_back( j );
 for( int i = 0 ; i < d_nc ; ++i )
  if( d_lbA[ i ] == d_ubA[ i ] )
   hp.eq_of.push_back( n + i );
  else
   if( ( d_lbA[ i ] > - Inf< double >() ) ||
       ( d_ubA[ i ] < Inf< double >() ) )
    hp.g_of.push_back( i );

 const int ne = hp.eq_of.size();
 const int nb = hp.idxb.size();
 const int ng = hp.g_of.size();

 // the matrices of HPIPM are column-major
 hp.A.assign( std::size_t( ne ) * n , 0 );
 hp.b.resize( ne );
 for( int k = 0 ; k < ne ; ++k )
  if( hp.eq_of[ k ] < n ) {
   hp.A[ std::size_t( hp.eq_of[ k ] ) * ne + k ] = 1;
   hp.b[ k ] = d_lb[ hp.eq_of[ k ] ];
   }
  else {
   const int i = hp.eq_of[ k ] - n;
   for( int j = 0 ; j < n ; ++j )
    hp.A[ std::size_t( j ) * ne + k ] = d_A[ std::size_t( i ) * n + j ];
   hp.b[ k ] = d_lbA[ i ];
   }

 auto side = [ & ]( double v , std::vector< double > & val ,
		    std::vector< double > & mask ) {
  const bool finite = std::abs( v ) < Inf< double >();
  val.push_back( finite ? v : 0 );
  mask.push_back( finite ? 1 : 0 );
  };

 hp.lb.clear(); hp.ub.clear(); hp.lb_mask.clear(); hp.ub_mask.clear();
 for( int k = 0 ; k < nb ; ++k ) {
  side( d_lb[ hp.idxb[ k ] ] , hp.lb , hp.lb_mask );
  side( d_ub[ hp.idxb[ k ] ] , hp.ub , hp.ub_mask );
  }

 hp.C.assign( std::size_t( ng ) * n , 0 );
 hp.lg.clear(); hp.ug.clear(); hp.lg_mask.clear(); hp.ug_mask.clear();
 for( int k = 0 ; k < ng ; ++k ) {
  const int i = hp.g_of[ k ];
  for( int j = 0 ; j < n ; ++j )
   hp.C[ std::size_t( j ) * ng + k ] = d_A[ std::size_t( i ) * n + j ];
  side( d_lbA[ i ] , hp.lg , hp.lg_mask );
  side( d_ubA[ i ] , hp.ug , hp.ug_mask );
  }

 // a warm start is from the solution of the previous solve, moved to the
 // current columns by load_problem(); the setters of HPIPM do not change
 // their arguments, but take them as non-const
 const bool warm = d_warm && d_has_guess &&
                   ( d_guess_x.size() == std::size_t( n ) );
 hp.create( n , ne , nb , ng , d_mode );
 if( warm )
  d_dense_qp_sol_set_v( d_guess_x.data() , &hp.sol );

 d_dense_qp_set_H( d_H.data() , &hp.qp );
 d_dense_qp_set_g( d_h.data() , &hp.qp );
 if( ne > 0 ) {
  d_dense_qp_set_A( hp.A.data() , &hp.qp );
  d_dense_qp_set_b( hp.b.data() , &hp.qp );
  }
 if( nb > 0 ) {
  d_dense_qp_set_idxb( hp.idxb.data() , &hp.qp );
  d_dense_qp_set_lb( hp.lb.data() , &hp.qp );
  d_dense_qp_set_ub( hp.ub.data() , &hp.qp );
  d_dense_qp_set_lb_mask( hp.lb_mask.data() , &hp.qp );
  d_dense_qp_set_ub_mask( hp.ub_mask.data() , &hp.qp );
  }
 if( ng > 0 ) {
  d_dense_qp_set_C( hp.C.data() , &hp.qp );
  d_dense_qp_set_lg( hp.lg.data() , &hp.qp );
  d_dense_qp_set_ug( hp.ug.data() , &hp.qp );
  d_dense_qp_set_lg_mask( hp.lg_mask.data() , &hp.qp );
  d_dense_qp_set_ug_mask( hp.ug_mask.data() , &hp.qp );
  }

 // the settings: those of the mode, then the tolerance, the iteration
 // limit and the warm start
 hp.set_default( d_mode );
 double tol = d_tol;
 d_dense_qp_ipm_arg_set_tol_stat( &tol , &hp.arg );
 d_dense_qp_ipm_arg_set_tol_eq( &tol , &hp.arg );
 d_dense_qp_ipm_arg_set_tol_ineq( &tol , &hp.arg );
 d_dense_qp_ipm_arg_set_tol_comp( &tol , &hp.arg );
 if( d_max_iter < Inf< int >() ) {
  int it = d_max_iter;
  d_dense_qp_ipm_arg_set_iter_max( &it , &hp.arg );
  }
 int ws = warm ? 1 : 0;
 d_dense_qp_ipm_arg_set_warm_start( &ws , &hp.arg );

 d_dense_qp_ipm_solve( &hp.qp , &hp.sol , &hp.arg , &hp.ws );

 int st , iter;
 d_dense_qp_ipm_get_status( &hp.ws , &st );
 d_dense_qp_ipm_get_iter( &hp.ws , &iter );
 elapsed_iter += iter;

 // a start on the boundary leaves the interior point little room, which
 // may keep it from ending: a warm-started solve that does not is repeated
 // from the initial point of HPIPM
 if( warm && ( st != SUCCESS ) ) {
  ws = 0;
  d_dense_qp_ipm_arg_set_warm_start( &ws , &hp.arg );
  d_dense_qp_ipm_solve( &hp.qp , &hp.sol , &hp.arg , &hp.ws );
  d_dense_qp_ipm_get_status( &hp.ws , &st );
  d_dense_qp_ipm_get_iter( &hp.ws , &iter );
  elapsed_iter += iter;
  }

 switch( st ) {
  case( SUCCESS ):   d_status = kOK;         break;
  case( MAX_ITER ):  d_status = kStopIter;   break;
  case( INCONS_EQ ): d_status = kInfeasible; break;
  default:           d_status = kError;  // MIN_STEP, NAN_SOL, FAILED_REG
  }

 hp.has_sol = ( d_status == kOK );
 if( d_status != kOK )
  return( sol_status = d_status );

 // the solution and the multipliers; HPIPM gives
 // H x + h = A' pi + lam_lb - lam_ub + C' ( lam_lg - lam_ug ), with
 // nonnegative lam, hence the multiplier of an equality is its pi and that
 // of an inequality the difference of the two lam, positive on an active
 // lower side and negative on an active upper one
 std::vector< double > pi( std::max( ne , 1 ) ) ,
  lam_lb( std::max( nb , 1 ) ) , lam_ub( std::max( nb , 1 ) ) ,
  lam_lg( std::max( ng , 1 ) ) , lam_ug( std::max( ng , 1 ) ) , none( 1 );
 d_dense_qp_sol_get_all( &hp.sol , d_x.data() , none.data() , none.data() ,
			 pi.data() , lam_lb.data() , lam_ub.data() ,
			 lam_lg.data() , lam_ug.data() , none.data() ,
			 none.data() );
 d_guess_x = d_x;
 d_has_guess = true;

 d_y.assign( n + d_nc , 0 );
 for( int k = 0 ; k < ne ; ++k )
  d_y[ hp.eq_of[ k ] ] = pi[ k ];
 for( int k = 0 ; k < nb ; ++k )
  d_y[ hp.idxb[ k ] ] = lam_lb[ k ] - lam_ub[ k ];
 for( int k = 0 ; k < ng ; ++k )
  d_y[ n + hp.g_of[ k ] ] = lam_lg[ k ] - lam_ug[ k ];

 // the value is computed here, on H as loaded
 for( int i = 0 ; i < n ; ++i ) {
  double Hx = 0;
  for( int j = 0 ; j < n ; ++j )
   Hx += d_H[ std::size_t( i ) * n + j ] * d_x[ j ];
  d_fval += d_x[ i ] * ( d_h[ i ] + Hx / 2 );
  }

 return( sol_status = d_status );

 }  // end( HPIPMMILPSolver::guts_of_compute )

/*--------------------------------------------------------------------------*/
/*---------------------- METHODS FOR READING RESULTS -----------------------*/
/*--------------------------------------------------------------------------*/

Solver::OFValue HPIPMMILPSolver::get_lb( void )
{
 if( solved() )
  return( d_sign * d_fval + constant_value );

 if( sol_status == kInfeasible )
  return( objsense == -1 ? - Inf< OFValue >() : Inf< OFValue >() );

 return( - Inf< OFValue >() );
 }

/*--------------------------------------------------------------------------*/

Solver::OFValue HPIPMMILPSolver::get_ub( void )
{
 if( solved() )
  return( d_sign * d_fval + constant_value );

 if( sol_status == kInfeasible )
  return( objsense == -1 ? - Inf< OFValue >() : Inf< OFValue >() );

 return( Inf< OFValue >() );
 }

/*--------------------------------------------------------------------------*/

bool HPIPMMILPSolver::has_var_solution( void )
{
 return( solved() );
 }

/*--------------------------------------------------------------------------*/

Solver::OFValue HPIPMMILPSolver::get_var_value( void )
{
 return( objsense == -1 ? get_lb() : get_ub() );
 }

/*--------------------------------------------------------------------------*/

void HPIPMMILPSolver::get_var_solution( Configuration * solc )
{
 if( ! solved() )
  throw( std::logic_error( "HPIPMMILPSolver::get_var_solution: "
			   "no solution available" ) );

 MILPSolver::write_var_solution( d_x );
 }

/*--------------------------------------------------------------------------*/

bool HPIPMMILPSolver::has_dual_solution( void )
{
 return( solved() );
 }

/*--------------------------------------------------------------------------*/

bool HPIPMMILPSolver::is_dual_feasible( void )
{
 return( solved() );
 }

/*--------------------------------------------------------------------------*/

void HPIPMMILPSolver::get_dual_solution( Configuration * solc )
{
 if( ! solved() )
  throw( std::logic_error( "HPIPMMILPSolver::get_dual_solution: "
			   "no solution available" ) );

 // the multipliers of the problem of the Block are y for a minimization and
 // - y for a maximization, whose Objective HPIPM has been given negated
 std::vector< double > rc( d_n );
 for( int j = 0 ; j < d_n ; ++j )
  rc[ j ] = d_sign * d_y[ j ];

 std::vector< double > pi( d_nc );
 for( int i = 0 ; i < d_nc ; ++i )
  pi[ i ] = d_sign * d_y[ d_n + i ];

 MILPSolver::write_dual_solution( pi , rc );
 }

/*--------------------------------------------------------------------------*/
/*------------------ METHODS FOR HANDLING THE PARAMETERS -------------------*/
/*--------------------------------------------------------------------------*/

void HPIPMMILPSolver::set_par( idx_type par , int value )
{
 switch( par ) {
  case( intMaxIter ):        d_max_iter = std::max( value , 1 );      return;
  case( intHPIPMMode ):      d_mode = std::clamp( value , 0 , 3 );    return;
  case( intHPIPMWarmStart ): d_warm = bool( value );                  return;
  }

 MILPSolver::set_par( par , value );
 }

/*--------------------------------------------------------------------------*/

void HPIPMMILPSolver::set_par( idx_type par , double value )
{
 switch( par ) {
  case( dblMaxTime ): d_max_time = std::max( value , 0.0 ); return;
  case( dblFAccSol ): d_tol = value;                        return;
  }

 MILPSolver::set_par( par , value );
 }

/*--------------------------------------------------------------------------*/

int HPIPMMILPSolver::get_dflt_int_par( idx_type par ) const
{
 switch( par ) {
  case( intHPIPMMode ):      return( 3 );
  case( intHPIPMWarmStart ): return( 0 );
  }

 return( MILPSolver::get_dflt_int_par( par ) );
 }

/*--------------------------------------------------------------------------*/

double HPIPMMILPSolver::get_dflt_dbl_par( idx_type par ) const
{
 if( par == dblFAccSol )
  return( HPIPM_dflt_tol );

 return( MILPSolver::get_dflt_dbl_par( par ) );
 }

/*--------------------------------------------------------------------------*/

int HPIPMMILPSolver::get_int_par( idx_type par ) const
{
 switch( par ) {
  case( intMaxIter ):        return( d_max_iter );
  case( intHPIPMMode ):      return( d_mode );
  case( intHPIPMWarmStart ): return( d_warm );
  }

 return( MILPSolver::get_int_par( par ) );
 }

/*--------------------------------------------------------------------------*/

double HPIPMMILPSolver::get_dbl_par( idx_type par ) const
{
 switch( par ) {
  case( dblMaxTime ): return( d_max_time );
  case( dblFAccSol ): return( d_tol );
  }

 return( MILPSolver::get_dbl_par( par ) );
 }

/*--------------------------------------------------------------------------*/

Solver::idx_type HPIPMMILPSolver::int_par_str2idx(
					      const std::string & name ) const
{
 auto it = std::find( HPIPM_int_pars_str.begin() ,
		      HPIPM_int_pars_str.end() , name );
 if( it != HPIPM_int_pars_str.end() )
  return( intLastAlgParMILP + ( it - HPIPM_int_pars_str.begin() ) );

 return( MILPSolver::int_par_str2idx( name ) );
 }

/*--------------------------------------------------------------------------*/

const std::string & HPIPMMILPSolver::int_par_idx2str( idx_type idx ) const
{
 if( ( idx >= intLastAlgParMILP ) && ( idx < intLastAlgParHPIPM ) )
  return( HPIPM_int_pars_str[ idx - intLastAlgParMILP ] );

 return( MILPSolver::int_par_idx2str( idx ) );
 }

/*--------------------------------------------------------------------------*/
/*--------------------- End File HPIPMMILPSolver.cpp -----------------------*/
/*--------------------------------------------------------------------------*/
