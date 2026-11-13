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

 // the working set of the previous problem, moved to the columns and rows
 // that are still there [see previous_column_index()], the new ones
 // inactive, so that a solve after columns or rows have been added or
 // removed starts from the previous one as well; it is kept on all the
 // columns and rows of MILPSolver, those that DAQP is not given included
 if( ( d_ws_n >= 0 ) && ( ( d_ws_n != numcols ) ||
			  ( d_ws.size() != std::size_t( numcols + numrows ) ) ) ) {
  const auto [ pn , pm ] = previous_sizes();
  if( ( d_ws_n == pn ) && ( d_ws.size() == std::size_t( pn + pm ) ) ) {
   std::vector< int > ws( numcols + numrows , 0 );
   for( int j = 0 ; j < numcols ; ++j )
    if( const int pj = previous_column_index( j ) ; pj >= 0 )
     ws[ j ] = d_ws[ pj ];
   for( int i = 0 ; i < numrows ; ++i )
    if( const int pi = previous_row_index( i ) ; pi >= 0 )
     ws[ numcols + i ] = d_ws[ pn + pi ];
   d_ws = std::move( ws );
   d_ws_n = numcols;
   }
  else {
   d_ws.clear();
   d_ws_n = -1;
   }
  }

 d_sign = ( objsense == -1 ) ? -1 : 1;  // DAQP minimizes

 // a fixed column is not given to DAQP, its value going to the sides of
 // the rows and to the linear term, unless all the columns are fixed
 d_xfix.assign( numcols , 0 );
 d_pos.assign( numcols , -1 );
 d_col.clear();
 for( int j = 0 ; j < numcols ; ++j )
  if( ( lb[ j ] == ub[ j ] ) && std::isfinite( lb[ j ] ) )
   d_xfix[ j ] = lb[ j ];
  else {
   d_pos[ j ] = d_col.size();
   d_col.push_back( j );
   }
 if( d_col.empty() )
  for( int j = 0 ; j < numcols ; ++j ) {
   d_pos[ j ] = j;
   d_col.push_back( j );
   }
 const int n = d_col.size();

 // the part of each row on the fixed columns, and whether the row has
 // entries on the others
 std::vector< double > shift( numrows , 0 );
 std::vector< char > used( numrows , 0 );
 for( int j = 0 ; j < numcols ; ++j )
  for( int k = matbeg[ j ] ; k < matbeg[ j ] + matcnt[ j ] ; ++k )
   if( d_pos[ j ] >= 0 )
    used[ matind[ k ] ] = 1;
   else
    shift[ matind[ k ] ] += matval[ k ] * d_xfix[ j ];

 // the sides of each row; a row with both sides infinite (e.g., a relaxed
 // one) is not given to DAQP, nor is one with no entries on the columns
 // that are given whose fixed part is within its sides
 std::vector< double > rbl( numrows ) , rbu( numrows );
 d_row.clear();
 d_rpos.assign( numrows , -1 );
 for( int i = 0 ; i < numrows ; ++i ) {
  switch( sense[ i ] ) {
   case( 'L' ): rbl[ i ] = - DAQP_INF;
                rbu[ i ] = daqp_side( rhs[ i ] );
                break;
   case( 'G' ): rbl[ i ] = daqp_side( rhs[ i ] );
                rbu[ i ] = DAQP_INF;
                break;
   case( 'E' ): rbl[ i ] = rbu[ i ] = rhs[ i ];
                break;
   case( 'R' ): if( rngval[ i ] > 0 ) {
                 rbl[ i ] = daqp_side( rhs[ i ] );
                 rbu[ i ] = daqp_side( rhs[ i ] + rngval[ i ] );
                 }
                else {
                 rbl[ i ] = daqp_side( rhs[ i ] + rngval[ i ] );
                 rbu[ i ] = daqp_side( rhs[ i ] );
                 }
                break;
   default:
    throw( std::invalid_argument( "DAQPMILPSolver::load_problem: row " +
				  std::to_string( i ) + " has sense " +
				  std::string( 1 , sense[ i ] ) ) );
   }
  if( ( rbl[ i ] <= - DAQP_INF ) && ( rbu[ i ] >= DAQP_INF ) )
   continue;
  if( ( ! used[ i ] ) && ( rbl[ i ] <= shift[ i ] ) &&
      ( shift[ i ] <= rbu[ i ] ) )
   continue;
  d_rpos[ i ] = d_row.size();
  d_row.push_back( i );
  }

 d_n = n;
 d_m = n + d_row.size();
 d_ncols = numcols;
 d_nrows = numrows;

 // the objective: c' x + sum_j q_j x_j^2 + sum_k ndq_k x_i x_j is
 // 1/2 x' H x + f' x with H_jj = 2 q_j and H_ij = H_ji = ndq_k; the term
 // H_ij x_i x_j / 2 of an entry on a fixed column goes to the linear term,
 // or to the constant if both columns are fixed
 d_f.assign( n , 0 );
 d_H.assign( std::size_t( n ) * n , 0 );
 d_fconst = 0;
 auto add_H = [ & ]( int i , int j , double h ) {
  const int pi = d_pos[ i ];
  const int pj = d_pos[ j ];
  if( ( pi >= 0 ) && ( pj >= 0 ) )
   d_H[ std::size_t( pi ) * n + pj ] += h;
  else if( pi >= 0 )
   d_f[ pi ] += h * d_xfix[ j ] / 2;
  else if( pj >= 0 )
   d_f[ pj ] += h * d_xfix[ i ] / 2;
  else
   d_fconst += h * d_xfix[ i ] * d_xfix[ j ] / 2;
  };

 for( int j = 0 ; j < numcols ; ++j ) {
  if( d_pos[ j ] >= 0 )
   d_f[ d_pos[ j ] ] += d_sign * objective[ j ];
  else
   d_fconst += d_sign * objective[ j ] * d_xfix[ j ];
  if( q_objective[ j ] != 0 )
   add_H( j , j , d_sign * 2 * q_objective[ j ] );
  }
 for( std::size_t k = 0 ; k < ndq_objective.size() ; ++k ) {
  const int i = ndq_rowind[ k ];
  const int j = ndq_colind[ k ];
  const double v = d_sign * ndq_objective[ k ];
  if( i == j )
   add_H( i , i , 2 * v );
  else {
   add_H( i , j , v );
   add_H( j , i , v );
   }
  }

 // the rows, dense and row-major, from the columns of MILPSolver
 d_A.assign( std::size_t( d_m - n ) * n , 0 );
 for( int k = 0 ; k < n ; ++k ) {
  const int j = d_col[ k ];
  for( int h = matbeg[ j ] ; h < matbeg[ j ] + matcnt[ j ] ; ++h )
   if( const int r = d_rpos[ matind[ h ] ] ; r >= 0 )
    d_A[ std::size_t( r ) * n + k ] += matval[ h ];
  }

 // the sides: the bounds first, then the rows less their fixed part; a
 // fixed variable and an equality row are constraints that DAQP keeps
 // active
 d_bu.resize( d_m );
 d_bl.resize( d_m );
 d_fixed.assign( d_m , 0 );
 for( int k = 0 ; k < n ; ++k ) {
  const int j = d_col[ k ];
  d_bl[ k ] = daqp_side( lb[ j ] );
  d_bu[ k ] = daqp_side( ub[ j ] );
  if( lb[ j ] == ub[ j ] )
   d_fixed[ k ] = DAQP_ACTIVE | DAQP_IMMUTABLE;
  }

 for( int r = n ; r < d_m ; ++r ) {
  const int i = d_row[ r - n ];
  d_bl[ r ] = ( rbl[ i ] <= - DAQP_INF ) ? - DAQP_INF : rbl[ i ] - shift[ i ];
  d_bu[ r ] = ( rbu[ i ] >= DAQP_INF ) ? DAQP_INF : rbu[ i ] - shift[ i ];
  if( d_bl[ r ] == d_bu[ r ] )
   d_fixed[ r ] = DAQP_ACTIVE | DAQP_IMMUTABLE;
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
 // a solve started from the working set of the previous one may end
 // unbounded or infeasible on a problem that a solve from scratch finds
 // optimal (DAQP 0.11 after one iteration, on a master problem of a
 // cutting-plane method whose components start empty): one more try with
 // no working set; likewise a solve that DAQP gives up as cycling, which
 // with a singular H it may do on a problem it then solves when allowed
 // more cycles (an SVM dual with more samples than features), is repeated
 // from scratch with 10 and then 100 times as many
 const int cycle_tol = d_settings.cycle_tol;
 int retry = 0;
 for( bool cold = false ; ; cold = true ) {
  bool warm = false;  // true if the solve starts from a working set
  if( d_loaded || ( ! d_has_work ) || cold ) {
   free_workspace();

   // the constraints the problem fixes, and the working set of the
   // previous solve on the others, if any and of the same size; a
   // constraint that the previous problem fixed and this one does not is
   // left out, since it was active for being fixed rather than at the
   // solution
   d_sense = d_fixed;
   if( d_warm_start && ( ! cold ) && ( d_ws_n == numcols ) &&
       ( d_ws.size() == std::size_t( numcols + numrows ) ) )
    for( int i = 0 ; i < d_m ; ++i ) {
     const int ws = d_ws[ ( i < d_n ) ? d_col[ i ]
                                      : numcols + d_row[ i - d_n ] ];
     if( ( ! d_fixed[ i ] ) && ( ! ( ws & DAQP_IMMUTABLE ) ) ) {
      d_sense[ i ] |= ws & ( DAQP_ACTIVE | 2 );
      if( ws & DAQP_ACTIVE )
       warm = true;
      }
     }

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
    d_settings.cycle_tol = cycle_tol;
    return( sol_status = decode_exitflag( flag ) );
    }
   d_has_work = true;
   d_loaded = false;
   }
  else
   warm = true;  // the workspace keeps the working set of its last solve

  // DAQP keeps its own iteration limit if none is set, and has no time
  // limit if its own is 0, so that a time limit of 0 is the smallest
  // positive one
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

  if( ( d_exitflag == DAQP_EXIT_CYCLE ) && ( retry++ < 2 ) ) {
   d_settings.cycle_tol = 10 * std::max( d_settings.cycle_tol , 1 );
   continue;
   }

  if( cold || ( ! warm ) || ( ( d_exitflag != DAQP_EXIT_UNBOUNDED ) &&
                              ( d_exitflag != DAQP_EXIT_INFEASIBLE ) ) )
   break;
  }
 d_settings.cycle_tol = cycle_tol;

 // with a singular H (the proximal-point iterations) DAQP may take a problem
 // unbounded below for an infeasible one: the rows are feasible if the same
 // problem with no linear term has a solution, and a convex QP with a
 // feasible point and no optimum is unbounded below
 if( ( d_exitflag == DAQP_EXIT_INFEASIBLE ) && rows_feasible() )
  d_exitflag = DAQP_EXIT_UNBOUNDED;

 // the value is computed here rather than taken from DAQP, whose value of a
 // problem solved by proximal-point iterations (e.g., a Linear Program) is
 // not that of the problem
 d_fval = d_fconst;
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

bool DAQPMILPSolver::rows_feasible( void )
{
 // copies of all that DAQP is given, which its setup may change in place
 std::vector< double > f0( d_n , 0 );
 std::vector< double > H( d_H ) , A( d_A ) , bu( d_bu ) , bl( d_bl );
 std::vector< int > sense( d_fixed );
 DAQPProblem qp = d_qp;
 qp.H = H.data();
 qp.f = f0.data();
 qp.A = A.empty() ? nullptr : A.data();
 qp.bupper = bu.data();
 qp.blower = bl.data();
 qp.sense = sense.data();
 DAQPSettings settings = d_settings;

 // daqp_quadprog() sets the problem up, solves it and frees what it made
 std::vector< double > x( d_n , 0 );
 std::vector< double > lam( d_m , 0 );
 DAQPResult res{};
 res.x = x.data();
 res.lam = lam.data();
 daqp_quadprog( &res , &qp , &settings );
 elapsed_iter += res.iter;
 const bool feasible = ( decode_exitflag( res.exitflag ) == kOK );
 return( feasible );
 }

/*--------------------------------------------------------------------------*/

void DAQPMILPSolver::free_workspace( void )
{
 if( ! d_has_work )
  return;

 // the working set, for the next start, on all the columns and rows of
 // MILPSolver (0 on those that DAQP has not been given)
 d_ws.assign( d_ncols + d_nrows , 0 );
 for( int k = 0 ; k < d_work.m ; ++k )
  d_ws[ ( k < d_work.n ) ? d_col[ k ]
                         : d_ncols + d_row[ k - d_work.n ] ] =
   d_work.sense[ k ];
 d_ws_n = d_ncols;

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

 std::vector< double > x( d_xfix );
 for( int k = 0 ; k < d_n ; ++k )
  x[ d_col[ k ] ] = d_x[ k ];

 MILPSolver::write_var_solution( x , solc );
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
 // and 0 on a row that DAQP has not been given; on a fixed column, the
 // multiplier of its bound is the one that makes its component of
 // H x + f + A' lam zero
 std::vector< double > pi( numrows , 0 );
 for( int r = d_n ; r < d_m ; ++r )
  pi[ d_row[ r - d_n ] ] = - d_sign * d_lam[ r ];

 std::vector< double > rc( numcols , 0 );
 for( int k = 0 ; k < d_n ; ++k )
  rc[ d_col[ k ] ] = - d_sign * d_lam[ k ];

 if( d_n < numcols ) {
  std::vector< double > x( d_xfix );
  for( int k = 0 ; k < d_n ; ++k )
   x[ d_col[ k ] ] = d_x[ k ];

  std::vector< double > g( numcols , 0 );  // H x + f + A' lam
  for( int j = 0 ; j < numcols ; ++j ) {
   g[ j ] = d_sign * ( objective[ j ] + 2 * q_objective[ j ] * x[ j ] );
   for( int h = matbeg[ j ] ; h < matbeg[ j ] + matcnt[ j ] ; ++h )
    if( const int r = d_rpos[ matind[ h ] ] ; r >= 0 )
     g[ j ] += matval[ h ] * d_lam[ d_n + r ];
   }
  for( std::size_t k = 0 ; k < ndq_objective.size() ; ++k ) {
   const int i = ndq_rowind[ k ];
   const int j = ndq_colind[ k ];
   const double v = d_sign * ndq_objective[ k ];
   if( i == j )
    g[ i ] += 2 * v * x[ i ];
   else {
    g[ i ] += v * x[ j ];
    g[ j ] += v * x[ i ];
    }
   }
  for( int j = 0 ; j < numcols ; ++j )
   if( d_pos[ j ] < 0 )
    rc[ j ] = d_sign * g[ j ];
  }

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
