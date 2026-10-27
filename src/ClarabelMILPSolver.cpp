/*--------------------------------------------------------------------------*/
/*----------------------- File ClarabelMILPSolver.cpp ----------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Implementation of the ClarabelMILPSolver class, which implements the
 * MILPSolver interface for convex Quadratic Programs with linear
 * constraints using the interior-point solver Clarabel.
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
#include <utility>

#include "ClarabelMILPSolver.h"

extern "C" {
#include <clarabel.h>
}

/*--------------------------------------------------------------------------*/
/*-------------------------------- USING -----------------------------------*/
/*--------------------------------------------------------------------------*/

using namespace SMSpp_di_unipi_it;

/*--------------------------------------------------------------------------*/
/*-------------------------------- CONSTANTS -------------------------------*/
/*--------------------------------------------------------------------------*/

static const std::array< std::string , 2 > Clarabel_dbl_pars_str = {
 "dblClarabelInfeasAbsTol" , "dblClarabelInfeasRelTol" };

// the value from which a side is infinite, that of the presolve of Clarabel

static constexpr double ClarabelInf = 1e20;

/*--------------------------------------------------------------------------*/
/*-------------------------------- FUNCTIONS -------------------------------*/
/*--------------------------------------------------------------------------*/

// the default settings of Clarabel

static const ClarabelDefaultSettings_f64 & clarabel_defaults( void )
{
 static const ClarabelDefaultSettings_f64 dflt =
  clarabel_DefaultSettings_f64_default();
 return( dflt );
 }

/*--------------------------------------------------------------------------*/

// the triplets ( column , row , value ) in compressed sparse column format
// with the given number of columns, sorted by row in each column and with
// the duplicates summed

static void triplets_to_csc( std::vector< std::pair< std::pair< int , int > ,
			                             double > > & t ,
			     int ncols , std::vector< std::uintptr_t > & beg ,
			     std::vector< std::uintptr_t > & ind ,
			     std::vector< double > & val )
{
 std::sort( t.begin() , t.end() ,
	    []( const auto & a , const auto & b ) {
	     return( a.first < b.first );
	     } );

 beg.assign( ncols + 1 , 0 );
 ind.clear();
 val.clear();
 ind.reserve( t.size() );
 val.reserve( t.size() );
 for( std::size_t k = 0 ; k < t.size() ; ) {
  const auto pos = t[ k ].first;
  double v = 0;
  for( ; ( k < t.size() ) && ( t[ k ].first == pos ) ; ++k )
   v += t[ k ].second;
  if( v == 0 )
   continue;
  ind.push_back( pos.second );
  val.push_back( v );
  ++beg[ pos.first + 1 ];
  }

 for( int j = 0 ; j < ncols ; ++j )
  beg[ j + 1 ] += beg[ j ];
 }

/*--------------------------------------------------------------------------*/

// a pointer that is never null, as Clarabel wants, to the data of a vector

template< class T >
static const T * nonnull( const std::vector< T > & v )
{
 static const T dummy{};
 return( v.empty() ? & dummy : v.data() );
 }

/*--------------------------------------------------------------------------*/
/*------------------------- FACTORY MANAGEMENT -----------------------------*/
/*--------------------------------------------------------------------------*/

SMSpp_insert_in_factory_cpp_0( ClarabelMILPSolver );

/*--------------------------------------------------------------------------*/
/*--------------------- CONSTRUCTOR AND DESTRUCTOR -------------------------*/
/*--------------------------------------------------------------------------*/

ClarabelMILPSolver::ClarabelMILPSolver( void ) : MILPSolver()
{
 const auto & dflt = clarabel_defaults();
 d_max_iter = dflt.max_iter;
 d_tol_gap_abs = dflt.tol_gap_abs;
 d_tol_gap_rel = dflt.tol_gap_rel;
 d_tol_feas = dflt.tol_feas;
 d_tol_infeas_abs = dflt.tol_infeas_abs;
 d_tol_infeas_rel = dflt.tol_infeas_rel;
 }

/*--------------------------------------------------------------------------*/
/*--------------------- DERIVED METHODS OF BASE CLASS ----------------------*/
/*--------------------------------------------------------------------------*/

void ClarabelMILPSolver::clear_problem( unsigned int what )
{
 MILPSolver::clear_problem( what );
 d_n = d_m = d_meq = 0;
 d_loaded = false;
 d_status = 0;
 }

/*--------------------------------------------------------------------------*/

void ClarabelMILPSolver::load_problem( void )
{
 MILPSolver::load_problem();

 if( numquadrows > 0 )
  throw( std::invalid_argument( "ClarabelMILPSolver::load_problem: "
				"quadratic constraints are not supported" ) );

 if( ( int_vars > 0 ) && ( relax_int_vars == 0 ) )
  throw( std::invalid_argument( "ClarabelMILPSolver::load_problem: integer "
				"variables are only supported if relaxed "
				"(intRelaxIntVars)" ) );

 const int n = numcols;
 d_n = n;
 d_sign = ( objsense == -1 ) ? -1 : 1;  // Clarabel minimizes

 // the objective: c' x + sum_j q_j x_j^2 + sum_k ndq_k x_i x_j is
 // 1/2 x' P x + q' x with P_jj = 2 q_j and P_ij = P_ji = ndq_k, of which
 // Clarabel only wants the upper triangle
 d_q.assign( n , 0 );
 for( int j = 0 ; j < n ; ++j )
  d_q[ j ] = d_sign * objective[ j ];

 std::vector< std::pair< std::pair< int , int > , double > > t;
 t.reserve( n + ndq_objective.size() );
 for( int j = 0 ; j < n ; ++j )
  if( q_objective[ j ] != 0 )
   t.push_back( { { j , j } , d_sign * 2 * q_objective[ j ] } );
 for( std::size_t k = 0 ; k < ndq_objective.size() ; ++k ) {
  const int i = ndq_rowind[ k ];
  const int j = ndq_colind[ k ];
  const double v = d_sign * ndq_objective[ k ];
  if( i == j )
   t.push_back( { { i , i } , 2 * v } );
  else
   t.push_back( { { std::max( i , j ) , std::min( i , j ) } , v } );
  }
 triplets_to_csc( t , n , d_Pbeg , d_Pind , d_Pval );

 // the sides of each constraint, the bounds first and then the rows, an
 // infinite one being none
 const int nc = n + numrows;
 std::vector< double > lo( nc ) , up( nc );
 for( int j = 0 ; j < n ; ++j ) {
  lo[ j ] = lb[ j ];
  up[ j ] = ub[ j ];
  }

 for( int i = 0 ; i < numrows ; ++i ) {
  const int r = n + i;
  switch( sense[ i ] ) {
   case( 'L' ): lo[ r ] = - Inf< double >();
                up[ r ] = rhs[ i ];
                break;
   case( 'G' ): lo[ r ] = rhs[ i ];
                up[ r ] = Inf< double >();
                break;
   case( 'E' ): lo[ r ] = up[ r ] = rhs[ i ];
                break;
   case( 'R' ): if( rngval[ i ] > 0 ) {
                 lo[ r ] = rhs[ i ];
                 up[ r ] = rhs[ i ] + rngval[ i ];
                 }
                else {
                 lo[ r ] = rhs[ i ] + rngval[ i ];
                 up[ r ] = rhs[ i ];
                 }
                break;
   default:
    throw( std::invalid_argument( "ClarabelMILPSolver::load_problem: row "
				  + std::to_string( i ) + " has sense " +
				  std::string( 1 , sense[ i ] ) ) );
   }
  }

 // the rows of A: a constraint with equal and finite sides is a row of the
 // zero cone, each other finite side a row of the nonnegative orthant, the
 // upper one as a' x + s = u and the lower one as - a' x + s = - l
 auto is_eq = [ & ]( int k ) {
  return( ( lo[ k ] == up[ k ] ) && ( std::abs( up[ k ] ) < ClarabelInf ) );
  };

 d_up.assign( nc , -1 );
 d_lo.assign( nc , -1 );
 d_b.clear();
 for( int k = 0 ; k < nc ; ++k )
  if( is_eq( k ) ) {
   d_up[ k ] = d_b.size();
   d_b.push_back( up[ k ] );
   }
 d_meq = d_b.size();

 for( int k = 0 ; k < nc ; ++k ) {
  if( is_eq( k ) )
   continue;
  if( up[ k ] < ClarabelInf ) {
   d_up[ k ] = d_b.size();
   d_b.push_back( up[ k ] );
   }
  if( lo[ k ] > - ClarabelInf ) {
   d_lo[ k ] = d_b.size();
   d_b.push_back( - lo[ k ] );
   }
  }
 d_m = d_b.size();

 // the columns of A, from those of MILPSolver and from the bounds
 t.clear();
 std::size_t nnz = n;
 for( int j = 0 ; j < n ; ++j )
  nnz += matcnt[ j ];
 t.reserve( 2 * nnz );
 auto add = [ & ]( int k , int j , double v ) {
  if( d_up[ k ] >= 0 )
   t.push_back( { { j , d_up[ k ] } , v } );
  if( d_lo[ k ] >= 0 )
   t.push_back( { { j , d_lo[ k ] } , - v } );
  };

 for( int j = 0 ; j < n ; ++j ) {
  add( j , j , 1 );
  for( int k = matbeg[ j ] ; k < matbeg[ j ] + matcnt[ j ] ; ++k )
   add( n + matind[ k ] , j , matval[ k ] );
  }
 triplets_to_csc( t , n , d_Abeg , d_Aind , d_Aval );

 d_x.assign( n , 0 );
 d_z.assign( d_m , 0 );
 d_loaded = true;
 d_status = 0;

 }  // end( ClarabelMILPSolver::load_problem )

/*--------------------------------------------------------------------------*/
/*-------------------- METHODS FOR SOLVING THE MODEL -----------------------*/
/*--------------------------------------------------------------------------*/

int ClarabelMILPSolver::guts_of_compute( void )
{
 d_status = ClarabelUnsolved;
 if( ! d_loaded )
  return( sol_status = kError );

 ClarabelCscMatrix_f64 P , A;
 clarabel_CscMatrix_f64_init( & P , d_n , d_n , d_Pbeg.data() ,
			      nonnull( d_Pind ) , nonnull( d_Pval ) );
 clarabel_CscMatrix_f64_init( & A , d_m , d_n , d_Abeg.data() ,
			      nonnull( d_Aind ) , nonnull( d_Aval ) );

 std::vector< ClarabelSupportedConeT_f64 > cones;
 if( d_meq > 0 ) {
  ClarabelSupportedConeT_f64 c{};
  c.tag = ClarabelZeroConeT_Tag;
  c.zero_cone_t = d_meq;
  cones.push_back( c );
  }
 if( d_m > d_meq ) {
  ClarabelSupportedConeT_f64 c{};
  c.tag = ClarabelNonnegativeConeT_Tag;
  c.nonnegative_cone_t = d_m - d_meq;
  cones.push_back( c );
  }

 // a time limit of 0 is the smallest positive one
 ClarabelDefaultSettings_f64 settings = clarabel_defaults();
 settings.max_iter = d_max_iter;
 settings.time_limit = ( d_max_time < Inf< double >() ) ?
  std::max( d_max_time , std::numeric_limits< double >::min() ) :
  Inf< double >();
 settings.verbose = ( log_verbosity > 0 );
 settings.tol_gap_abs = d_tol_gap_abs;
 settings.tol_gap_rel = d_tol_gap_rel;
 settings.tol_feas = d_tol_feas;
 settings.tol_infeas_abs = d_tol_infeas_abs;
 settings.tol_infeas_rel = d_tol_infeas_rel;

 auto solver = clarabel_DefaultSolver_f64_new( & P , nonnull( d_q ) , & A ,
					       nonnull( d_b ) , cones.size() ,
					       nonnull( cones ) ,
					       & settings );
 if( ! solver )  // the data are not accepted (Clarabel says why)
  return( sol_status = kError );

 clarabel_DefaultSolver_f64_solve( solver );
 const auto sol = clarabel_DefaultSolver_f64_solution( solver );
 d_status = sol.status;
 // x is also the ray of a dual infeasible problem, and z the certificate of
 // a primal infeasible one
 if( solved() || ( d_status == ClarabelDualInfeasible ) )
  std::copy( sol.x , sol.x + d_n , d_x.begin() );
 if( solved() || ( d_status == ClarabelPrimalInfeasible ) )
  std::copy( sol.z , sol.z + d_m , d_z.begin() );
 clarabel_DefaultSolver_f64_free( solver );

 // the value is computed here rather than taken from Clarabel, so that it
 // is that of the solution which is given out
 d_fval = 0;
 if( solved() )
  for( int j = 0 ; j < d_n ; ++j ) {
   d_fval += d_q[ j ] * d_x[ j ];
   for( auto k = d_Pbeg[ j ] ; k < d_Pbeg[ j + 1 ] ; ++k ) {
    const auto i = d_Pind[ k ];
    d_fval += ( i == csc_idx( j ) ? 0.5 : 1 ) * d_Pval[ k ] *
              d_x[ i ] * d_x[ j ];
    }
   }

 return( sol_status = decode_status( d_status ) );

 }  // end( ClarabelMILPSolver::guts_of_compute )

/*--------------------------------------------------------------------------*/

int ClarabelMILPSolver::decode_status( int status )
{
 switch( status ) {
  case( ClarabelSolved ):             return( kOK );
  case( ClarabelAlmostSolved ):       return( kLowPrecision );
  case( ClarabelPrimalInfeasible ):   return( kInfeasible );
  case( ClarabelDualInfeasible ):     return( kUnbounded );
  case( ClarabelMaxIterations ):      return( kStopIter );
  case( ClarabelMaxTime ):            return( kStopTime );
  default:                            return( kError );
  }
 }

/*--------------------------------------------------------------------------*/

bool ClarabelMILPSolver::solved( void ) const
{
 return( ( d_status == ClarabelSolved ) ||
	 ( d_status == ClarabelAlmostSolved ) );
 }

/*--------------------------------------------------------------------------*/
/*---------------------- METHODS FOR READING RESULTS -----------------------*/
/*--------------------------------------------------------------------------*/

Solver::OFValue ClarabelMILPSolver::get_lb( void )
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

Solver::OFValue ClarabelMILPSolver::get_ub( void )
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

bool ClarabelMILPSolver::has_var_solution( void )
{
 return( solved() );
 }

/*--------------------------------------------------------------------------*/

Solver::OFValue ClarabelMILPSolver::get_var_value( void )
{
 return( objsense == -1 ? get_lb() : get_ub() );
 }

/*--------------------------------------------------------------------------*/

void ClarabelMILPSolver::get_var_solution( Configuration * solc )
{
 if( ! solved() )
  throw( std::logic_error( "ClarabelMILPSolver::get_var_solution: "
			   "no solution available" ) );

 MILPSolver::write_var_solution( d_x , solc );
 }

/*--------------------------------------------------------------------------*/

bool ClarabelMILPSolver::has_dual_solution( void )
{
 return( solved() );
 }

/*--------------------------------------------------------------------------*/

bool ClarabelMILPSolver::is_dual_feasible( void )
{
 return( solved() );
 }

/*--------------------------------------------------------------------------*/

void ClarabelMILPSolver::get_dual_solution( Configuration * solc )
{
 if( ! solved() )
  throw( std::logic_error( "ClarabelMILPSolver::get_dual_solution: "
			   "no solution available" ) );

 // Clarabel has P x + q + A' z = 0, with z >= 0 on the nonnegative orthant:
 // a constraint whose upper side is the row a' x + s = u and whose lower
 // side is the row - a' x + s = - l has multiplier lam = z_u - z_l (z on
 // its row if it is an equality), positive on an active upper side and
 // negative on an active lower one, for the problem Clarabel minimizes;
 // the multipliers of the problem of the Block are - lam for a
 // minimization and lam for a maximization, whose Objective Clarabel has
 // been given negated
 auto lam = [ & ]( int k ) {
  double l = 0;
  if( d_up[ k ] >= 0 )
   l += d_z[ d_up[ k ] ];
  if( d_lo[ k ] >= 0 )
   l -= d_z[ d_lo[ k ] ];
  return( - d_sign * l );
  };

 std::vector< double > rc( d_n );
 for( int j = 0 ; j < d_n ; ++j )
  rc[ j ] = lam( j );

 std::vector< double > pi( numrows );
 for( int i = 0 ; i < numrows ; ++i )
  pi[ i ] = lam( d_n + i );

 MILPSolver::write_dual_solution( pi , rc , solc );
 }

/*--------------------------------------------------------------------------*/

bool ClarabelMILPSolver::has_var_direction( void )
{
 return( d_loaded && ( d_status == ClarabelDualInfeasible ) );
 }

/*--------------------------------------------------------------------------*/

void ClarabelMILPSolver::get_var_direction( Configuration * dirc )
{
 if( ! has_var_direction() )
  throw( std::logic_error( "ClarabelMILPSolver::get_var_direction: "
			   "no unbounded direction available" ) );

 // P x = 0 and q' x < 0 for the problem Clarabel minimizes, whose q is
 // that of the Block negated for a maximization: x is a ray of the Block
 MILPSolver::write_var_solution( d_x , dirc );
 }

/*--------------------------------------------------------------------------*/

bool ClarabelMILPSolver::has_dual_direction( void )
{
 return( d_loaded && ( d_status == ClarabelPrimalInfeasible ) );
 }

/*--------------------------------------------------------------------------*/

void ClarabelMILPSolver::get_dual_direction( Configuration * dirc )
{
 if( ! has_dual_direction() )
  throw( std::logic_error( "ClarabelMILPSolver::get_dual_direction: "
			   "no dual unbounded direction available" ) );

 // A' z = 0 and b' z < 0: the multiplier of a constraint is that of its
 // lower side minus that of its upper one, nonnegative where the lower side
 // is used, which is the sign of the Farkas certificates of the other
 // :MILPSolver, and the value of the certificate is - b' z
 auto y = [ & ]( int k ) {
  double v = 0;
  if( d_lo[ k ] >= 0 )
   v += d_z[ d_lo[ k ] ];
  if( d_up[ k ] >= 0 )
   v -= d_z[ d_up[ k ] ];
  return( v );
  };

 double value = 0;
 for( int r = 0 ; r < d_m ; ++r )
  value -= d_b[ r ] * d_z[ r ];

 std::vector< double > pi( numrows );
 for( int i = 0 ; i < numrows ; ++i )
  pi[ i ] = y( d_n + i );

 // the multipliers of the columns, c - A' y, or - A' y with
 // intHomogeneousDirection, computed out of the matrix so that they are
 // those of y exactly
 std::vector< double > dj( d_n , 0 );
 for( int j = 0 ; j < d_n ; ++j ) {
  if( ! homogeneous_direction )
   dj[ j ] = objective[ j ];
  for( int k = matbeg[ j ] ; k < matbeg[ j ] + matcnt[ j ] ; ++k )
   dj[ j ] -= matval[ k ] * pi[ matind[ k ] ];
  }

 f_dual_direction_value = value;

 MILPSolver::write_dual_solution( pi , dj , dirc );
 }

/*--------------------------------------------------------------------------*/
/*------------------ METHODS FOR HANDLING THE PARAMETERS -------------------*/
/*--------------------------------------------------------------------------*/

void ClarabelMILPSolver::set_par( idx_type par , int value )
{
 switch( par ) {
  case( intMaxIter ): d_max_iter = std::max( value , 0 ); return;
  case( intLogVerb ): log_verbosity = value;              return;
  }

 MILPSolver::set_par( par , value );
 }

/*--------------------------------------------------------------------------*/

void ClarabelMILPSolver::set_par( idx_type par , double value )
{
 switch( par ) {
  case( dblMaxTime ): d_max_time = std::max( value , 0.0 ); return;
  case( dblAbsAcc ):  d_tol_gap_abs = value;                return;
  case( dblRelAcc ):  d_tol_gap_rel = value;                return;
  case( dblFAccSol ): d_tol_feas = value;                   return;
  case( dblClarabelInfeasAbsTol ): d_tol_infeas_abs = value; return;
  case( dblClarabelInfeasRelTol ): d_tol_infeas_rel = value; return;
  }

 MILPSolver::set_par( par , value );
 }

/*--------------------------------------------------------------------------*/

int ClarabelMILPSolver::get_dflt_int_par( idx_type par ) const
{
 if( par == intMaxIter )
  return( clarabel_defaults().max_iter );

 return( MILPSolver::get_dflt_int_par( par ) );
 }

/*--------------------------------------------------------------------------*/

double ClarabelMILPSolver::get_dflt_dbl_par( idx_type par ) const
{
 const auto & dflt = clarabel_defaults();
 switch( par ) {
  case( dblAbsAcc ):               return( dflt.tol_gap_abs );
  case( dblRelAcc ):               return( dflt.tol_gap_rel );
  case( dblFAccSol ):              return( dflt.tol_feas );
  case( dblClarabelInfeasAbsTol ): return( dflt.tol_infeas_abs );
  case( dblClarabelInfeasRelTol ): return( dflt.tol_infeas_rel );
  }

 return( MILPSolver::get_dflt_dbl_par( par ) );
 }

/*--------------------------------------------------------------------------*/

int ClarabelMILPSolver::get_int_par( idx_type par ) const
{
 switch( par ) {
  case( intMaxIter ): return( d_max_iter );
  case( intLogVerb ): return( log_verbosity );
  }

 return( MILPSolver::get_int_par( par ) );
 }

/*--------------------------------------------------------------------------*/

double ClarabelMILPSolver::get_dbl_par( idx_type par ) const
{
 switch( par ) {
  case( dblMaxTime ):              return( d_max_time );
  case( dblAbsAcc ):               return( d_tol_gap_abs );
  case( dblRelAcc ):               return( d_tol_gap_rel );
  case( dblFAccSol ):              return( d_tol_feas );
  case( dblClarabelInfeasAbsTol ): return( d_tol_infeas_abs );
  case( dblClarabelInfeasRelTol ): return( d_tol_infeas_rel );
  }

 return( MILPSolver::get_dbl_par( par ) );
 }

/*--------------------------------------------------------------------------*/

Solver::idx_type ClarabelMILPSolver::dbl_par_str2idx(
					      const std::string & name ) const
{
 auto it = std::find( Clarabel_dbl_pars_str.begin() ,
		      Clarabel_dbl_pars_str.end() , name );
 if( it != Clarabel_dbl_pars_str.end() )
  return( dblLastAlgParMILP + ( it - Clarabel_dbl_pars_str.begin() ) );

 return( MILPSolver::dbl_par_str2idx( name ) );
 }

/*--------------------------------------------------------------------------*/

const std::string & ClarabelMILPSolver::dbl_par_idx2str( idx_type idx )
 const
{
 if( ( idx >= dblLastAlgParMILP ) && ( idx < dblLastAlgParClarabel ) )
  return( Clarabel_dbl_pars_str[ idx - dblLastAlgParMILP ] );

 return( MILPSolver::dbl_par_idx2str( idx ) );
 }

/*--------------------------------------------------------------------------*/
/*-------------------- End File ClarabelMILPSolver.cpp ---------------------*/
/*--------------------------------------------------------------------------*/
