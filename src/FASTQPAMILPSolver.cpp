/*--------------------------------------------------------------------------*/
/*----------------------- File FASTQPAMILPSolver.cpp -----------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Implementation of the FASTQPAMILPSolver class, which implements the
 * MILPSolver interface for strictly convex Mixed-Integer Quadratic Programs
 * with linear constraints by a depth-first branch-and-bound whose nodes are
 * bounded by the dual of their continuous relaxation, solved by the
 * feasible active-set method FAST-QPA.
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
#include <chrono>
#include <cmath>
#include <limits>
#include <stdexcept>

#include "FASTQPAMILPSolver.h"

/*--------------------------------------------------------------------------*/
/*-------------------------------- USING -----------------------------------*/
/*--------------------------------------------------------------------------*/

using namespace SMSpp_di_unipi_it;

/*--------------------------------------------------------------------------*/
/*-------------------------------- CONSTANTS -------------------------------*/
/*--------------------------------------------------------------------------*/

static const std::array< std::string , 5 > FQ_int_pars_str = {
 "intFQWarmStart" , "intFQEarlyPrune" , "intFQIntBoundsAsRows" ,
   "intFQMaxQPIter" ,
 "intFQPostSolve"
 };

static const std::array< std::string , 7 > FQ_dbl_pars_str = {
 "dblFQOptTol" , "dblFQPostOptTol" , "dblFQEpsFactor" , "dblFQDelta" ,
 "dblFQGamma" ,  "dblFQEta" ,        "dblFQPDTol"
 };

// a lower bound at least this large is taken as infinite: the dual of the
// node is unbounded, i.e., its relaxation is infeasible

static constexpr double FQ_BIG = 1e20;

// multipliers this large make the dual of a node be checked for a
// certificate of the infeasibility of its relaxation

static constexpr double FQ_LAM_BIG = 1e6;

// how far from an integer the bound of an integer variable can be

static constexpr double FQ_INT_TOL = 1e-9;

/*--------------------------------------------------------------------------*/
/*-------------------------------- FUNCTIONS -------------------------------*/
/*--------------------------------------------------------------------------*/

// true if v is a finite side or bound in the sense of SMS++

static bool is_finite( double v )
{
 return( ( v < Inf< double >() ) && ( v > -Inf< double >() ) );
 }

/*--------------------------------------------------------------------------*/
// the 2-norm of the projected gradient x - [ x - g ]^+ of min q , x >= 0

static double pg_norm( const double * x , const double * g , int m )
{
 double s = 0;
 for( int i = 0 ; i < m ; ++i ) {
  const double v = std::min( x[ i ] , g[ i ] );
  s += v * v;
  }
 return( std::sqrt( s ) );
 }

/*--------------------------------------------------------------------------*/
/*------------------------- FACTORY MANAGEMENT -----------------------------*/
/*--------------------------------------------------------------------------*/

SMSpp_insert_in_factory_cpp_0( FASTQPAMILPSolver );

/*--------------------------------------------------------------------------*/
/*--------------------- CONSTRUCTOR AND DESTRUCTOR -------------------------*/
/*--------------------------------------------------------------------------*/

FASTQPAMILPSolver::FASTQPAMILPSolver( void ) : MILPSolver()
{ }

/*--------------------------------------------------------------------------*/

FASTQPAMILPSolver::~FASTQPAMILPSolver() = default;

/*--------------------------------------------------------------------------*/
/*--------------------- DERIVED METHODS OF BASE CLASS ----------------------*/
/*--------------------------------------------------------------------------*/

void FASTQPAMILPSolver::clear_problem( unsigned int what )
{
 MILPSolver::clear_problem( what );
 f_dirty = true;
 f_prep = false;
 f_status = kUnEval;
 f_has_sol = false;
 }

/*--------------------------------------------------------------------------*/

void FASTQPAMILPSolver::load_problem( void )
{
 MILPSolver::load_problem();

 if( numquadrows > 0 )
  throw( std::invalid_argument( "FASTQPAMILPSolver::load_problem: "
				"quadratic constraints are not supported" ) );

 f_dirty = true;
 f_prep = false;
 f_status = kUnEval;
 f_has_sol = false;
 }

/*--------------------------------------------------------------------------*/
/*-------------------- METHODS FOR SOLVING THE MODEL -----------------------*/
/*--------------------------------------------------------------------------*/

void FASTQPAMILPSolver::convert( void )
{
 const int n = numcols;
 f_n = n;
 f_sign = ( objsense == -1 ) ? -1 : 1;
 f_bad_data = false;

 // the integer variables first, each group in the order of the columns
 f_perm.clear();
 f_perm.reserve( n );
 for( int j = 0 ; j < n ; ++j )
  if( ( std::size_t( j ) < xctype.size() ) && ( xctype[ j ] != 'C' ) )
   f_perm.push_back( j );
 f_n1 = f_perm.size();
 for( int j = 0 ; j < n ; ++j )
  if( ( std::size_t( j ) >= xctype.size() ) || ( xctype[ j ] == 'C' ) )
   f_perm.push_back( j );

 std::vector< int > ip( n );
 for( int k = 0 ; k < n ; ++k )
  ip[ f_perm[ k ] ] = k;

 // the objective: c' x + sum_j q_j x_j^2 + sum_k ndq_k x_i x_j of the Block
 // is x' Q x + c' x with Q_jj = q_j and Q_ij = Q_ji = ndq_k / 2, negated
 // for a maximization
 f_c.assign( n , 0 );
 f_Q.assign( std::size_t( n ) * n , 0 );
 for( int j = 0 ; j < n ; ++j ) {
  const int k = ip[ j ];
  f_c[ k ] = f_sign * objective[ j ];
  if( std::size_t( j ) < q_objective.size() )
   f_Q[ std::size_t( k ) * n + k ] = f_sign * q_objective[ j ];
  }
 for( std::size_t h = 0 ; h < ndq_objective.size() ; ++h ) {
  const int i = ip[ ndq_rowind[ h ] ];
  const int j = ip[ ndq_colind[ h ] ];
  const double v = f_sign * ndq_objective[ h ];
  if( i == j )
   f_Q[ std::size_t( i ) * n + i ] += v;
  else {
   f_Q[ std::size_t( i ) * n + j ] += v / 2;
   f_Q[ std::size_t( j ) * n + i ] += v / 2;
   }
  }

 // the windows of the integer variables
 f_ilo.assign( f_n1 , -Inf< double >() );
 f_ihi.assign( f_n1 , Inf< double >() );
 for( int k = 0 ; k < f_n1 ; ++k ) {
  const int j = f_perm[ k ];
  if( is_finite( lb[ j ] ) )
   f_ilo[ k ] = std::ceil( lb[ j ] - FQ_INT_TOL );
  if( is_finite( ub[ j ] ) )
   f_ihi[ k ] = std::floor( ub[ j ] + FQ_INT_TOL );
  if( f_ilo[ k ] > f_ihi[ k ] )
   f_bad_data = true;
  }

 // the rows of the Block, dense, in the order of the internal columns
 std::vector< double > R( std::size_t( numrows ) * n , 0 );
 for( int j = 0 ; j < n ; ++j )
  for( int h = matbeg[ j ] ; h < matbeg[ j ] + matcnt[ j ] ; ++h )
   R[ std::size_t( matind[ h ] ) * n + ip[ j ] ] += matval[ h ];

 f_A.clear();
 f_b.clear();
 f_m = 0;

 // adds the row s a' x <= s rhs, a being given by the n values at a (or
 // being e_k if a is nullptr); an empty row is only checked
 auto add_row = [ & ]( const double * a , int k , double s , double rhs ) {
  bool empty = true;
  if( a ) {
   for( int j = 0 ; j < n ; ++j )
    if( a[ j ] != 0 ) {
     empty = false;
     break;
     }
   }
  else
   empty = false;

  if( empty ) {
   if( s * rhs < -f_feas_tol * std::max( 1.0 , std::abs( rhs ) ) )
    f_bad_data = true;
   return;
   }

  const auto pos = f_A.size();
  f_A.resize( pos + n , 0 );
  if( a )
   for( int j = 0 ; j < n ; ++j )
    f_A[ pos + j ] = s * a[ j ];
  else
   f_A[ pos + k ] = s;
  f_b.push_back( s * rhs );
  ++f_m;
  };

 for( int i = 0 ; i < numrows ; ++i ) {
  const double * a = R.data() + std::size_t( i ) * n;
  switch( sense[ i ] ) {
   case( 'L' ):
    if( is_finite( rhs[ i ] ) )
     add_row( a , 0 , 1 , rhs[ i ] );
    break;
   case( 'G' ):
    if( is_finite( rhs[ i ] ) )
     add_row( a , 0 , -1 , rhs[ i ] );
    break;
   case( 'E' ):
    add_row( a , 0 , 1 , rhs[ i ] );
    add_row( a , 0 , -1 , rhs[ i ] );
    break;
   case( 'R' ): {
    const double lo = rhs[ i ] + std::min( rngval[ i ] , 0.0 );
    const double hi = rhs[ i ] + std::max( rngval[ i ] , 0.0 );
    if( is_finite( hi ) )
     add_row( a , 0 , 1 , hi );
    if( is_finite( lo ) )
     add_row( a , 0 , -1 , lo );
    break;
    }
   default:
    throw( std::invalid_argument( "FASTQPAMILPSolver::convert: row " +
				  std::to_string( i ) + " has sense " +
				  std::string( 1 , sense[ i ] ) ) );
   }
  }

 // the bounds of the continuous variables, and if so asked of the integer
 // ones (their windows, already rounded)
 for( int k = 0 ; k < n ; ++k ) {
  double l = lb[ f_perm[ k ] ];
  double u = ub[ f_perm[ k ] ];
  if( k < f_n1 ) {
   if( ! f_int_bounds_rows )
    continue;
   l = f_ilo[ k ];
   u = f_ihi[ k ];
   }
  if( is_finite( u ) )
   add_row( nullptr , k , 1 , u );
  if( is_finite( l ) )
   add_row( nullptr , k , -1 , l );
  }

 f_bsc.resize( f_m );
 for( int i = 0 ; i < f_m ; ++i )
  f_bsc[ i ] = std::max( 1.0 , std::abs( f_b[ i ] ) );

 }  // end( FASTQPAMILPSolver::convert )

/*--------------------------------------------------------------------------*/

bool FASTQPAMILPSolver::preprocess( void )
{
 const int n = f_n;
 const int m = f_m;
 const int n1 = f_n1;

 // the Cholesky factor L of Q, which must be positive definite
 double maxd = 0;
 for( int j = 0 ; j < n ; ++j )
  maxd = std::max( maxd , f_Q[ std::size_t( j ) * n + j ] );

 std::vector< double > L( std::size_t( n ) * n , 0 );
 for( int j = 0 ; j < n ; ++j ) {
  double d = f_Q[ std::size_t( j ) * n + j ];
  for( int k = 0 ; k < j ; ++k )
   d -= L[ std::size_t( j ) * n + k ] * L[ std::size_t( j ) * n + k ];
  if( ( ! ( d > f_pd_tol * maxd ) ) || ( ! ( d > 0 ) ) )
   return( false );
  const double ljj = std::sqrt( d );
  L[ std::size_t( j ) * n + j ] = ljj;
  for( int i = j + 1 ; i < n ; ++i ) {
   double v = f_Q[ std::size_t( i ) * n + j ];
   for( int k = 0 ; k < j ; ++k )
    v -= L[ std::size_t( i ) * n + k ] * L[ std::size_t( j ) * n + k ];
   L[ std::size_t( i ) * n + j ] = v / ljj;
   }
  }

 // Q^{-1} = L^{-T} L^{-1}, column by column
 std::vector< double > Qi( std::size_t( n ) * n , 0 );
 {
  std::vector< double > w( n );
  for( int c = 0 ; c < n ; ++c ) {
   for( int i = 0 ; i < n ; ++i ) {  // L w = e_c
    double v = ( i == c ) ? 1 : 0;
    for( int k = 0 ; k < i ; ++k )
     v -= L[ std::size_t( i ) * n + k ] * w[ k ];
    w[ i ] = v / L[ std::size_t( i ) * n + i ];
    }
   for( int i = n - 1 ; i >= 0 ; --i ) {  // L' w = w
    double v = w[ i ];
    for( int k = i + 1 ; k < n ; ++k )
     v -= L[ std::size_t( k ) * n + i ] * w[ k ];
    w[ i ] = v / L[ std::size_t( i ) * n + i ];
    }
   for( int i = 0 ; i < n ; ++i )
    Qi[ std::size_t( i ) * n + c ] = w[ i ];
   }
  }

 // the minimum with no row, - Q^{-1} c / 2, and its value c' x / 2
 f_xroot.assign( n , 0 );
 f_froot = 0;
 for( int i = 0 ; i < n ; ++i ) {
  double v = 0;
  for( int j = 0 ; j < n ; ++j )
   v += Qi[ std::size_t( i ) * n + j ] * f_c[ j ];
  f_xroot[ i ] = -v / 2;
  f_froot += f_c[ i ] * f_xroot[ i ] / 2;
  }

 // the Hessian of the dual of the root, A Q^{-1} A' / 4
 std::vector< double > Qt( std::size_t( m ) * m , 0 );
 {
  std::vector< double > T( std::size_t( n ) * m , 0 );  // Q^{-1} A'
  for( int i = 0 ; i < n ; ++i )
   for( int r = 0 ; r < m ; ++r ) {
    double v = 0;
    for( int j = 0 ; j < n ; ++j )
     v += Qi[ std::size_t( i ) * n + j ] * f_A[ std::size_t( r ) * n + j ];
    T[ std::size_t( i ) * m + r ] = v;
    }
  for( int r = 0 ; r < m ; ++r )
   for( int t = r ; t < m ; ++t ) {
    double v = 0;
    for( int i = 0 ; i < n ; ++i )
     v += f_A[ std::size_t( r ) * n + i ] * T[ std::size_t( i ) * m + t ];
    Qt[ std::size_t( r ) * m + t ] = Qt[ std::size_t( t ) * m + r ] = v / 4;
    }
  }

 // the last column of each row
 std::vector< int > last( m , -1 );
 for( int r = 0 ; r < m ; ++r )
  for( int j = n - 1 ; j >= 0 ; --j )
   if( f_A[ std::size_t( r ) * n + j ] != 0 ) {
    last[ r ] = j;
    break;
    }

 f_z.assign( n1 , {} );
 f_Az.assign( n1 , {} );
 f_s.assign( n1 , 0 );
 f_Qt.assign( n1 + 1 , {} );
 f_lmax.assign( n1 + 1 , 0 );
 f_eps.assign( n1 + 1 , 1 );

 // the depths l = 0 ... n1 - 1, Qi being Q_l^{-1}, of size k = n - l
 for( int l = 0 ; l <= n1 ; ++l ) {
  const int k = n - l;

  // a row on fixed variables alone has an exactly zero row of the Hessian,
  // which the downdates only leave tiny: the dual must see it unbounded
  // when the row is violated
  for( int r = 0 ; r < m ; ++r )
   if( last[ r ] < l )
    for( int t = 0 ; t < m ; ++t )
     Qt[ std::size_t( r ) * m + t ] = Qt[ std::size_t( t ) * m + r ] = 0;

  f_lmax[ l ] = max_eigenvalue( Qt , m );
  f_eps[ l ] = ( f_lmax[ l ] > 0 ) ? f_eps_factor / ( 2 * f_lmax[ l ] ) : 1;
  f_Qt[ l ] = Qt;

  if( l == n1 )
   break;

  // u = Q_l^{-1} e_1 and the direction z^l = u / u_1 along which the
  // minimum moves when the first free variable is fixed
  const double u0 = Qi[ 0 ];
  if( ! ( u0 > 0 ) )
   return( false );
  auto & z = f_z[ l ];
  z.resize( k );
  for( int j = 0 ; j < k ; ++j )
   z[ j ] = Qi[ std::size_t( j ) * k ] / u0;
  f_s[ l ] = 1 / u0;

  auto & Az = f_Az[ l ];
  Az.assign( m , 0 );
  for( int r = 0 ; r < m ; ++r ) {
   const double * a = f_A.data() + std::size_t( r ) * n + l;
   double v = 0;
   for( int j = 0 ; j < k ; ++j )
    v += a[ j ] * z[ j ];
   Az[ r ] = v;
   }

  // the next depth: Q_{l+1}^{-1} is the Schur complement of the first
  // pivot of Q_l^{-1}, and the Hessian of the dual loses ( A u )( A u )'
  // / ( 4 u_1 ) = u_1 A z^l ( A z^l )' / 4
  for( int r = 0 ; r < m ; ++r )
   for( int t = 0 ; t < m ; ++t )
    Qt[ std::size_t( r ) * m + t ] -= u0 * Az[ r ] * Az[ t ] / 4;

  std::vector< double > Qn( std::size_t( k - 1 ) * ( k - 1 ) );
  for( int i = 1 ; i < k ; ++i ) {
   const double ui = Qi[ std::size_t( i ) * k ];
   for( int j = 1 ; j < k ; ++j )
    Qn[ std::size_t( i - 1 ) * ( k - 1 ) + j - 1 ] =
      Qi[ std::size_t( i ) * k + j ] - ui * Qi[ j ] / u0;
   }
  Qi.swap( Qn );
  }

 // the data of the leaves: Q_{n1}^{-1} and - Q_{n1}^{-1} A_{n1}' / 2
 const int k = n - n1;
 f_Qil = Qi;
 f_Wl.assign( std::size_t( k ) * m , 0 );
 for( int j = 0 ; j < k ; ++j )
  for( int r = 0 ; r < m ; ++r ) {
   const double * a = f_A.data() + std::size_t( r ) * n + n1;
   double v = 0;
   for( int t = 0 ; t < k ; ++t )
    v += Qi[ std::size_t( j ) * k + t ] * a[ t ];
   f_Wl[ std::size_t( j ) * m + r ] = -v / 2;
   }

 // the workspace of FAST-QPA
 for( auto * w : { &w_xt , &w_gt , &w_d , &w_y , &w_r , &w_p , &w_Qp , &w_xb ,
      &w_gb ,
		   &w_xn , &w_gn , &w_x2 , &w_g2 , &w_pl } )
  w->assign( m , 0 );
 w_N.resize( m );

 return( true );

 }  // end( FASTQPAMILPSolver::preprocess )

/*--------------------------------------------------------------------------*/

double FASTQPAMILPSolver::max_eigenvalue( std::vector< double > M , int k )
{
 // cyclic Jacobi rotations until the off-diagonal part is negligible
 if( k == 0 )
  return( 0 );

 double fro = 0;
 for( auto v : M )
  fro += v * v;
 if( fro == 0 )
  return( 0 );

 for( int sweep = 0 ; sweep < 100 ; ++sweep ) {
  double off = 0;
  for( int p = 0 ; p < k ; ++p )
   for( int q = p + 1 ; q < k ; ++q )
    off += M[ std::size_t( p ) * k + q ] * M[ std::size_t( p ) * k + q ];
  if( off <= 1e-30 * fro )
   break;

  for( int p = 0 ; p < k ; ++p )
   for( int q = p + 1 ; q < k ; ++q ) {
    const double apq = M[ std::size_t( p ) * k + q ];
    if( apq == 0 )
     continue;
    const double app = M[ std::size_t( p ) * k + p ];
    const double aqq = M[ std::size_t( q ) * k + q ];
    const double theta = ( aqq - app ) / ( 2 * apq );
    const double t = ( theta >= 0 ? 1 : -1 ) /
		     ( std::abs( theta ) + std::sqrt( theta * theta + 1 ) );
    const double c = 1 / std::sqrt( t * t + 1 );
    const double s = t * c;
    for( int r = 0 ; r < k ; ++r ) {  // columns p and q
     const double mrp = M[ std::size_t( r ) * k + p ];
     const double mrq = M[ std::size_t( r ) * k + q ];
     M[ std::size_t( r ) * k + p ] = c * mrp - s * mrq;
     M[ std::size_t( r ) * k + q ] = s * mrp + c * mrq;
     }
    for( int r = 0 ; r < k ; ++r ) {  // rows p and q
     const double mpr = M[ std::size_t( p ) * k + r ];
     const double mqr = M[ std::size_t( q ) * k + r ];
     M[ std::size_t( p ) * k + r ] = c * mpr - s * mqr;
     M[ std::size_t( q ) * k + r ] = s * mpr + c * mqr;
     }
    }
  }

 double lmax = 0;
 for( int p = 0 ; p < k ; ++p )
  lmax = std::max( lmax , M[ std::size_t( p ) * k + p ] );

 // a margin for the rounding, an overestimate being harmless
 return( lmax * ( 1 + 1e-10 ) );
 }

/*--------------------------------------------------------------------------*/

double FASTQPAMILPSolver::eval_dual( int l , const double * ct ,
				     const double * x , double * g ) const
{
 // g = 2 Qt x + ct over the nonzeros of x, and q = x' ( g + ct ) / 2
 const int m = f_m;
 const double * Qt = f_Qt[ l ].data();
 std::copy( ct , ct + m , g );
 for( int j = 0 ; j < m ; ++j )
  if( x[ j ] != 0 ) {
   const double xj = 2 * x[ j ];
   const double * col = Qt + std::size_t( j ) * m;  // Qt is symmetric
   for( int i = 0 ; i < m ; ++i )
    g[ i ] += col[ i ] * xj;
   }

 double q = 0;
 for( int i = 0 ; i < m ; ++i )
  if( x[ i ] != 0 )
   q += x[ i ] * ( g[ i ] + ct[ i ] );
 return( q / 2 );
 }

/*--------------------------------------------------------------------------*/

int FASTQPAMILPSolver::fast_qpa( int l , const double * ct , double qcut ,
				 double tol , double * lam , char * act ,
				 bool use_act , double * g , double & qv ,
				 double lam_big )
{
 const int m = f_m;
 const double * Qt = f_Qt[ l ].data();
 const double eps = f_eps[ l ];
 const double eta = f_eta * f_lmax[ l ];

 double * xt = w_xt.data();
 double * gt = w_gt.data();
 double * d = w_d.data();
 double * y = w_y.data();
 double * r = w_r.data();
 double * p = w_p.data();
 double * Qp = w_Qp.data();
 double * xb = w_xb.data();
 double * gb = w_gb.data();
 double * pl = w_pl.data();
 int * N = w_N.data();

 // the active-set estimate of definition 2, x_i <= eps g_i
 auto estimate = [ & ]( const double * x , const double * gx ) {
  for( int i = 0 ; i < m ; ++i )
   act[ i ] = ! ( x[ i ] > eps * gx[ i ] );
  };

 qv = eval_dual( l , ct , lam , g );
 if( qv <= qcut )
  return( 1 );

 for( int it = 0 ; ; ++it ) {
  if( pg_norm( lam , g , m ) <= tol ) {
   estimate( lam , g );
   return( 0 );
   }

  if( it >= std::min( f_max_qp_iter , f_it_cap ) ) {
   estimate( lam , g );
   return( 2 );
   }

  ++f_qp_it;

  // step 3: the estimate, inherited from the father at the first iteration
  if( ! ( use_act && ( it == 0 ) ) )
   estimate( lam , g );

  // step 4: the active variables to 0, and the gradient there
  std::copy( lam , lam + m , xt );
  std::copy( g , g + m , gt );
  for( int j = 0 ; j < m ; ++j )
   if( act[ j ] && ( lam[ j ] != 0 ) ) {
    const double xj = 2 * lam[ j ];
    const double * col = Qt + std::size_t( j ) * m;
    for( int i = 0 ; i < m ; ++i )
     gt[ i ] -= col[ i ] * xj;
    xt[ j ] = 0;
    }
  double qt = 0;
  for( int i = 0 ; i < m ; ++i )
   if( xt[ i ] != 0 )
    qt += xt[ i ] * ( gt[ i ] + ct[ i ] );
  qt /= 2;

  if( qt <= qcut ) {
   std::copy( xt , xt + m , lam );
   std::copy( gt , gt + m , g );
   qv = qt;
   return( 1 );
   }

  int nN = 0;
  for( int i = 0 ; i < m ; ++i )
   if( ! act[ i ] )
    N[ nN++ ] = i;

  // step 6: algorithm 2, the conjugate gradient on min y' Qt_NN y + gt_N' y
  // giving the direction d_N (truncated at the first direction of small
  // curvature) and the candidate x_bar = x_tilde + y
  bool check = true;
  bool ray = false;
  int nl = 0;
  double rr = 0;
  for( int h = 0 ; h < nN ; ++h ) {
   r[ h ] = gt[ N[ h ] ];
   p[ h ] = -r[ h ];
   y[ h ] = 0;
   rr += r[ h ] * r[ h ];
   }
  const double rr0 = rr;

  for( ; ; ) {
   if( rr <= 1e-30 * std::max( rr0 , 1.0 ) )
    break;

   double pQp = 0;
   double pp = 0;
   for( int h = 0 ; h < nN ; ++h ) {
    const double * row = Qt + std::size_t( N[ h ] ) * m;
    double v = 0;
    for( int t = 0 ; t < nN ; ++t )
     v += row[ N[ t ] ] * p[ t ];
    Qp[ h ] = v;
    pQp += p[ h ] * v;
    pp += p[ h ] * p[ h ];
    }

   // the first direction of small curvature is kept as a candidate
   // direction of recession
   if( ( ! ray ) && ( pQp <= eta * pp ) ) {
    std::copy( p , p + nN , pl );
    ray = true;
    }

   if( ! ( pQp > 0 ) )
    break;

   if( check && ( pQp <= eta * pp ) ) {
    for( int h = 0 ; h < nN ; ++h )
     d[ h ] = ( nl == 0 ) ? -gt[ N[ h ] ] : y[ h ];
    check = false;
    }

   const double alpha = rr / ( 2 * pQp );
   double rrn = 0;
   for( int h = 0 ; h < nN ; ++h ) {
    y[ h ] += alpha * p[ h ];
    r[ h ] += 2 * alpha * Qp[ h ];
    rrn += r[ h ] * r[ h ];
    }
   ++nl;
   if( nl > 2 * nN + 5 )
    break;

   const double beta = rrn / rr;
   for( int h = 0 ; h < nN ; ++h )
    p[ h ] = -r[ h ] + beta * p[ h ];
   rr = rrn;
   }

  if( check )
   for( int h = 0 ; h < nN ; ++h )
    d[ h ] = ( nl == 0 ) ? -gt[ N[ h ] ] : y[ h ];

  double dd = 0;
  double gd = 0;
  for( int h = 0 ; h < nN ; ++h ) {
   dd += d[ h ] * d[ h ];
   gd += gt[ N[ h ] ] * d[ h ];
   }

  // step 8: the optimality of x_bar
  std::fill( xb , xb + m , 0.0 );
  bool xbfeas = true;
  for( int h = 0 ; h < nN ; ++h ) {
   const double v = xt[ N[ h ] ] + ( ( nl == 0 ) ? -gt[ N[ h ] ] : y[ h ] );
   xb[ N[ h ] ] = v;
   if( v < 0 )
    xbfeas = false;
   }
  if( xbfeas ) {
   const double qb = eval_dual( l , ct , xb , gb );
   if( ( qb <= qcut ) || ( pg_norm( xb , gb , m ) <= tol ) ) {
    std::copy( xb , xb + m , lam );
    std::copy( gb , gb + m , g );
    qv = qb;
    if( qb <= qcut )
     return( 1 );
    estimate( lam , g );
    return( 0 );
    }
   }

  // no descent direction on N: go on from x_tilde with a new estimate
  if( ( dd == 0 ) || ( ! ( gd < 0 ) ) ) {
   std::copy( xt , xt + m , lam );
   std::copy( gt , gt + m , g );
   qv = qt;
   use_act = false;
   continue;
   }

  // step 9: the projected Armijo line search along d
  double * xn = w_xn.data();
  double * gn = w_gn.data();
  double alpha = 1;
  double qn = 0;
  bool accepted = false;
  int j = 0;
  for( ; j < 60 ; ++j , alpha *= f_delta ) {
   std::fill( xn , xn + m , 0.0 );
   for( int h = 0 ; h < nN ; ++h )
    xn[ N[ h ] ] = std::max( 0.0 , xt[ N[ h ] ] + alpha * d[ h ] );
   qn = eval_dual( l , ct , xn , gn );
   if( qn <= qt + f_gamma * alpha * gd ) {
    accepted = true;
    break;
    }
   }

  if( ! accepted ) {
   std::copy( xt , xt + m , lam );
   std::copy( gt , gt + m , g );
   qv = qt;
   estimate( lam , g );
   return( 3 );
   }

  // the unit step accepted while the slope along d is still more than half
  // of that at x_tilde, as it is along a direction of recession of an
  // unbounded dual: the step is doubled as long as the value keeps
  // decreasing enough, so that the value of the dual reaches the cutoff in
  // a few iterations
  auto slope = [ & ]( const double * x , const double * gx ,
		      const double * dir ) {
   double sl = 0;
   for( int h = 0 ; h < nN ; ++h )
    if( ( x[ N[ h ] ] > 0 ) || ( dir[ h ] > 0 ) )
     sl += gx[ N[ h ] ] * dir[ h ];
   return( sl );
   };

  if( ( j == 0 ) && ( slope( xn , gn , d ) < gd / 2 ) ) {
   double * x2 = w_x2.data();
   double * g2 = w_g2.data();
   for( int t = 0 ; ( t < 100 ) && ( qn > qcut ) ; ++t ) {
    const double a2 = 2 * alpha;
    std::fill( x2 , x2 + m , 0.0 );
    for( int h = 0 ; h < nN ; ++h )
     x2[ N[ h ] ] = std::max( 0.0 , xt[ N[ h ] ] + a2 * d[ h ] );
    const double q2 = eval_dual( l , ct , x2 , g2 );
    if( ! ( ( q2 <= qt + f_gamma * a2 * gd ) && ( q2 < qn ) ) )
     break;
    alpha = a2;
    qn = q2;
    std::swap( w_xn , w_x2 );
    std::swap( w_gn , w_g2 );
    xn = w_xn.data();
    gn = w_gn.data();
    x2 = w_x2.data();
    g2 = w_g2.data();
    if( slope( xn , gn , d ) >= gd / 2 )
     break;
    }
   }

  // along the direction of small curvature, if the value keeps decreasing
  // along it from the new point, the step is doubled as long as it does:
  // along a direction of recession of an unbounded dual the value reaches
  // the cutoff, or the multipliers grow enough to be a certificate of
  // infeasibility, in a few iterations, rather than growing linearly
  if( ray && ( qn > qcut ) ) {
   double * x2 = w_x2.data();
   double * g2 = w_g2.data();
   double pn = 0;
   for( int h = 0 ; h < nN ; ++h )
    pn += pl[ h ] * pl[ h ];
   double beta = ( dd > 0 ) && ( pn > 0 ) ? std::sqrt( dd / pn ) : 1;
   for( int t = 0 ; ( t < 100 ) && ( qn > qcut ) ; ++t , beta *= 2 ) {
    const double s0 = slope( xn , gn , pl );
    if( ! ( s0 < 0 ) )
     break;
    std::copy( xn , xn + m , x2 );
    for( int h = 0 ; h < nN ; ++h )
     x2[ N[ h ] ] = std::max( 0.0 , xn[ N[ h ] ] + beta * pl[ h ] );
    const double q2 = eval_dual( l , ct , x2 , g2 );
    if( ! ( q2 <= qn + f_gamma * beta * s0 ) )
     break;
    qn = q2;
    std::swap( w_xn , w_x2 );
    std::swap( w_gn , w_g2 );
    xn = w_xn.data();
    gn = w_gn.data();
    x2 = w_x2.data();
    g2 = w_g2.data();
    }
   }

  std::copy( xn , xn + m , lam );
  std::copy( gn , gn + m , g );
  qv = qn;
  if( qv <= qcut )
   return( 1 );

  if( *std::max_element( lam , lam + m ) > lam_big ) {
   estimate( lam , g );
   return( 4 );
   }
  }
 }  // end( FASTQPAMILPSolver::fast_qpa )

/*--------------------------------------------------------------------------*/

bool FASTQPAMILPSolver::infeasible( int l , const double * ct ,
				    const double * xh ,
				    const double * lam ) const
{
 // u = lam / sum( lam ) proves the rows of the node infeasible if u' A_l =
 // 0 and u' b < 0, b being their right-hand sides after the fixings, i.e.,
 // u' ct + u' A_l xh with xh the minimum with no row
 const int n = f_n;
 const int m = f_m;
 const int k = n - l;

 double su = 0;
 for( int r = 0 ; r < m ; ++r )
  su += lam[ r ];
 if( ! ( su > 0 ) )
  return( false );

 double ub = 0;
 double amax = 0;
 std::vector< double > ua( k , 0 );
 for( int r = 0 ; r < m ; ++r ) {
  const double ur = lam[ r ] / su;
  ub += ur * ct[ r ];
  if( ur == 0 )
   continue;
  const double * a = f_A.data() + std::size_t( r ) * n + l;
  for( int j = 0 ; j < k ; ++j ) {
   ua[ j ] += ur * a[ j ];
   amax = std::max( amax , std::abs( a[ j ] ) );
   }
  }

 double res = 0;
 for( int j = 0 ; j < k ; ++j ) {
  res = std::max( res , std::abs( ua[ j ] ) );
  ub += ua[ j ] * xh[ j ];
  }

 return( ( res <= 1e-9 * std::max( amax , 1.0 ) ) && ( ub < -f_feas_tol ) );
 }

/*--------------------------------------------------------------------------*/

double FASTQPAMILPSolver::value_of( const std::vector< double > & x ) const
{
 const int n = f_n;
 double v = 0;
 for( int i = 0 ; i < n ; ++i ) {
  double qx = 0;
  const double * row = f_Q.data() + std::size_t( i ) * n;
  for( int j = 0 ; j < n ; ++j )
   qx += row[ j ] * x[ j ];
  v += x[ i ] * ( qx + f_c[ i ] );
  }
 return( v );
 }

/*--------------------------------------------------------------------------*/

double FASTQPAMILPSolver::prune_eps( double ub ) const
{
 const double abs_acc = is_finite( f_abs_acc ) ? f_abs_acc : 0;
 const double rel_acc = is_finite( f_rel_acc ) ? f_rel_acc : 0;
 return( std::max( abs_acc , rel_acc * std::max( 1.0 , std::abs( ub ) ) ) );
 }

/*--------------------------------------------------------------------------*/

bool FASTQPAMILPSolver::solve_leaf( const std::vector< double > & fix ,
				    std::vector< double > & lam , double tol ,
				    std::vector< double > & x )
{
 const int n = f_n;
 const int n1 = f_n1;
 const int m = f_m;
 const int k = n - n1;

 // the linear term of the free variables, c_j + 2 sum_i Q_ji r_i, and
 // the minimum with no row y = - Q_{n1}^{-1} c / 2
 std::vector< double > cb( k );
 for( int j = 0 ; j < k ; ++j ) {
  const double * row = f_Q.data() + std::size_t( n1 + j ) * n;
  double v = f_c[ n1 + j ];
  for( int i = 0 ; i < n1 ; ++i )
   v += 2 * row[ i ] * fix[ i ];
  cb[ j ] = v;
  }

 x.assign( n , 0 );
 for( int i = 0 ; i < n1 ; ++i )
  x[ i ] = fix[ i ];
 for( int j = 0 ; j < k ; ++j ) {
  double v = 0;
  for( int t = 0 ; t < k ; ++t )
   v += f_Qil[ std::size_t( j ) * k + t ] * cb[ t ];
  x[ n1 + j ] = -v / 2;
  }

 std::vector< double > ct( m );
 for( int r = 0 ; r < m ; ++r ) {
  const double * a = f_A.data() + std::size_t( r ) * n;
  double v = 0;
  for( int j = 0 ; j < n ; ++j )
   v += a[ j ] * x[ j ];
  ct[ r ] = f_b[ r ] - v;
  }

 std::vector< char > act( m , 0 );
 std::vector< double > g( m );
 double qv;
 lam.resize( m , 0 );
 const int res = fast_qpa( n1 , ct.data() , -Inf< double >() , tol ,
			   lam.data() ,
			   act.data() , false , g.data() , qv ,
			   Inf< double >() );
 if( ( res == 1 ) ||
     ( res == 3 && pg_norm( lam.data() , g.data() , m ) > f_feas_tol ) )
  return( false );

 // the gradient of the dual is the slack of the rows
 for( int r = 0 ; r < m ; ++r )
  if( g[ r ] < -f_feas_tol * f_bsc[ r ] )
   return( false );

 for( int j = 0 ; j < k ; ++j ) {
  const double * w = f_Wl.data() + std::size_t( j ) * m;
  double v = 0;
  for( int r = 0 ; r < m ; ++r )
   v += w[ r ] * lam[ r ];
  x[ n1 + j ] += v;
  }

 return( true );
 }

/*--------------------------------------------------------------------------*/

int FASTQPAMILPSolver::guts_of_compute( void )
{
 const auto t0 = std::chrono::steady_clock::now();
 auto elapsed = [ & ]() {
  return(
    std::chrono::duration< double >( std::chrono::steady_clock::now() - t0 )
      .count() );
  };

 f_has_sol = false;
 f_UB = Inf< double >();
 f_LB = -Inf< double >();
 f_nodes = f_qp_it = f_root_it = 0;
 f_x.clear();

 if( f_dirty ) {
  convert();
  f_dirty = false;
  f_prep = false;
  }

 if( ! f_prep ) {
  f_pd = preprocess();
  f_prep = true;
  }

 if( ! f_pd ) {
  if( f_log && ( log_verbosity > 0 ) )
   *f_log << "FASTQPAMILPSolver: the quadratic part of the Objective is "
	  << "not positive definite" << std::endl;
  return( f_status = kError );
  }

 if( f_bad_data ) {
  f_LB = Inf< double >();
  return( f_status = kInfeasible );
  }

 const int n = f_n;
 const int n1 = f_n1;
 const int m = f_m;

 // the cutoffs, in the internal form (minimization, no constant term)
 double upcut , lwcut;
 if( f_sign > 0 ) {
  upcut = f_up_cutoff - constant_value;
  lwcut = f_lw_cutoff - constant_value;
  }
 else {
  upcut = constant_value - f_lw_cutoff;
  lwcut = constant_value - f_up_cutoff;
  }
 if( std::isnan( upcut ) )
  upcut = Inf< double >();
 if( std::isnan( lwcut ) )
  lwcut = -Inf< double >();

 double UB = upcut;

 // a node is pruned when its lower bound reaches this
 auto cut = [ & ]() {
  return( is_finite( UB ) ? std::min( UB - prune_eps( UB ) , FQ_BIG )
			   : FQ_BIG );
  };

 // the state of the nodes along the current path, one per depth
 std::vector< std::vector< double > > xh( n1 + 1 );  // the minimum, no row
 std::vector< std::vector< double > > Axh( n1 + 1 ,
					   std::vector< double >( m ) );
 std::vector< std::vector< double > > lam( n1 + 1 ,
					   std::vector< double >( m ) );
 std::vector< std::vector< char > > act( n1 + 1 , std::vector< char >( m ) );
 for( int l = 0 ; l <= n1 ; ++l )
  xh[ l ].resize( n - l );
 std::vector< double > fu( n1 + 1 );  // the value of the minimum, no row
 std::vector< double > LB( n1 + 1 );  // the lower bound
 std::vector< double > v( n1 + 1 );   // the value of the branching variable
 std::vector< double > dn( n1 + 1 );  // the next child below v
 std::vector< double > up( n1 + 1 );  // the next child above v
 std::vector< char > dn_open( n1 + 1 );
 std::vector< char > up_open( n1 + 1 );
 std::vector< char > side( n1 + 1 );  // the side of the last child
 std::vector< double > fix( n1 );     // the values of the integer ones

 std::vector< double > ct( m );
 std::vector< double > g( m );
 std::vector< double > best_fix;
 std::vector< double > best_lam;
 std::vector< double > xbest;
 std::vector< double > xleaf( n );

 // the evaluation of the node at depth cl, whose data are set: returns its
 // lower bound, after having made it the incumbent if it is a leaf whose
 // solution is better than that
 bool accepted = false;
 bool postponed = false;
 long unsolved = 0;  // number of leaves neither pruned nor solved
 auto eval = [ & ]( int cl , bool use_act , bool capped ) -> double {
  ++f_nodes;
  accepted = postponed = false;
  for( int r = 0 ; r < m ; ++r )
   ct[ r ] = f_b[ r ] - Axh[ cl ][ r ];

  if( ( cl == n1 ) && ( n1 == n ) ) {  // a leaf with no free variable
   for( int r = 0 ; r < m ; ++r )
    if( ct[ r ] < -f_feas_tol * f_bsc[ r ] )
     return( Inf< double >() );
   if( fu[ cl ] < UB ) {
    xleaf.assign( fix.begin() , fix.end() );
    const double val = value_of( xleaf );
    if( val < UB ) {
     UB = val;
     xbest = xleaf;
     f_has_sol = accepted = true;
     }
    }
   return( fu[ cl ] );
   }

  const bool leaf = ( cl == n1 );
  const double cutv = cut();
  const double qcut = fu[ cl ] - ( f_early_prune ? cutv : FQ_BIG );
  const double tol = leaf ? std::min( f_opt_tol , f_feas_tol ) : f_opt_tol;
  // with no incumbent, a dual that does not converge quickly may be an
  // unbounded one, which no cutoff stops: its node is postponed
  capped = capped && ( ! is_finite( UB ) );
  f_it_cap = capped ? 3 * m + 20 : Inf< int >();
  double qv;
  int res;
  for( double lam_big = FQ_LAM_BIG ; ; lam_big *= 1e3 ) {
   res = fast_qpa( cl , ct.data() , qcut , tol , lam[ cl ].data() ,
		   act[ cl ].data() , use_act , g.data() , qv ,
		   lam_big < FQ_BIG ? lam_big : Inf< double >() );
   if( res != 4 )
    break;
   if( infeasible( cl , ct.data() , xh[ cl ].data() , lam[ cl ].data() ) ) {
    f_it_cap = Inf< int >();
    return( Inf< double >() );
    }
   if( capped )
    break;
   use_act = true;
   }
  f_it_cap = Inf< int >();
  if( ( m > 0 ) &&
      ( *std::max_element( lam[ cl ].begin() , lam[ cl ].end() ) >
	FQ_LAM_BIG ) &&
      infeasible( cl , ct.data() , xh[ cl ].data() , lam[ cl ].data() ) )
   return( Inf< double >() );
  const double lb = fu[ cl ] - qv;
  if( lb >= FQ_BIG )
   return( Inf< double >() );
  if( capped && ( ( res == 2 ) || ( res == 4 ) ) && ( lb < cutv ) ) {
   postponed = true;
   return( lb );
   }
  if( ( res == 1 ) || ( lb >= cutv ) || ( ! leaf ) )
   return( lb );

  // a leaf: the gradient of the dual is the slack of the rows; a leaf whose
  // dual is not solved and whose primal solution is not feasible is left
  // with its lower bound alone
  for( int r = 0 ; r < m ; ++r )
   if( g[ r ] < -f_feas_tol * f_bsc[ r ] ) {
    ++unsolved;
    return( lb );
    }

  // its value is fu + lam' Qt lam, i.e., fu + lam' ( g - ct ) / 2
  double val = fu[ cl ];
  for( int r = 0 ; r < m ; ++r )
   val += lam[ cl ][ r ] * ( g[ r ] - ct[ r ] ) / 2;
  if( ! ( val < UB ) )
   return( lb );

  const int k = n - n1;
  for( int i = 0 ; i < n1 ; ++i )
   xleaf[ i ] = fix[ i ];
  for( int j = 0 ; j < k ; ++j ) {
   const double * w = f_Wl.data() + std::size_t( j ) * m;
   double x = xh[ cl ][ j ];
   for( int r = 0 ; r < m ; ++r )
    x += w[ r ] * lam[ cl ][ r ];
   xleaf[ n1 + j ] = x;
   }
  val = value_of( xleaf );

  // a bound above the value of the solution is the sign of a dual that is
  // unbounded but for the rounding: the leaf is left with its bound alone
  if( lb > val + prune_eps( val ) ) {
   ++unsolved;
   return( lb );
   }

  if( val < UB ) {
   UB = val;
   xbest = xleaf;
   best_fix.assign( fix.begin() , fix.end() );
   best_lam = lam[ cl ];
   f_has_sol = accepted = true;
   }
  return( lb );
  };

 // the children of the node at depth l: the integer values of its
 // branching variable by increasing distance from v[ l ], within its window
 auto setup_children = [ & ]( int l ) {
  double vv = xh[ l ][ 0 ];
  for( int r = 0 ; r < m ; ++r )
   vv -= f_Az[ l ][ r ] * lam[ l ][ r ] / ( 2 * f_s[ l ] );
  v[ l ] = vv;
  const double fl = std::floor( vv );
  dn[ l ] = std::min( fl , f_ihi[ l ] );
  dn_open[ l ] = ( dn[ l ] >= f_ilo[ l ] );
  up[ l ] = std::max( fl + 1 , f_ilo[ l ] );
  up_open[ l ] = ( up[ l ] <= f_ihi[ l ] );
  };

 auto next_child = [ & ]( int l , double & r ) -> bool {
  bool down;
  if( dn_open[ l ] && up_open[ l ] )
   down = ( v[ l ] - dn[ l ] <= up[ l ] - v[ l ] );
  else if( dn_open[ l ] )
   down = true;
  else if( up_open[ l ] )
   down = false;
  else
   return( false );

  if( down ) {
   r = dn[ l ];
   dn[ l ] -= 1;
   dn_open[ l ] = ( dn[ l ] >= f_ilo[ l ] );
   side[ l ] = 0;
   }
  else {
   r = up[ l ];
   up[ l ] += 1;
   up_open[ l ] = ( up[ l ] <= f_ihi[ l ] );
   side[ l ] = 1;
   }
  return( true );
  };

 // the child of the node at depth l fixing its branching variable to r:
 // the minimum with no row moves along z^l, its value grows by alpha^2 s^l
 auto move_to_child = [ & ]( int l , double r ) {
  const int k = n - l;
  const double alpha = r - xh[ l ][ 0 ];
  const auto & z = f_z[ l ];
  auto & xc = xh[ l + 1 ];
  const auto & xp = xh[ l ];
  for( int j = 1 ; j < k ; ++j )
   xc[ j - 1 ] = xp[ j ] + alpha * z[ j ];
  for( int i = 0 ; i < m ; ++i )
   Axh[ l + 1 ][ i ] = Axh[ l ][ i ] + alpha * f_Az[ l ][ i ];
  fu[ l + 1 ] = fu[ l ] + alpha * alpha * f_s[ l ];
  fix[ l ] = r;
  };

 auto make_child = [ & ]( int l , double r ) {
  move_to_child( l , r );
  if( f_warm_start ) {
   lam[ l + 1 ] = lam[ l ];
   act[ l + 1 ] = act[ l ];
   }
  else {
   std::fill( lam[ l + 1 ].begin() , lam[ l + 1 ].end() , 0.0 );
   std::fill( act[ l + 1 ].begin() , act[ l + 1 ].end() , 0 );
   }
  };

 // the root
 std::copy( f_xroot.begin() , f_xroot.end() , xh[ 0 ].begin() );
 for( int r = 0 ; r < m ; ++r ) {
  const double * a = f_A.data() + std::size_t( r ) * n;
  double s = 0;
  for( int j = 0 ; j < n ; ++j )
   s += a[ j ] * f_xroot[ j ];
  Axh[ 0 ][ r ] = s;
  }
 fu[ 0 ] = f_froot;

 double LBmin = Inf< double >();  // the least bound of the pruned nodes
 int stop = kOK;
 int lvl = -1;

 // the sides of the nodes whose first child has been postponed: the node
 // at depth l, its fixings, multipliers, active set, bound and branching
 // value, the side and its next child
 struct Postponed
 {
  int l;
  std::vector< double > fix;
  std::vector< double > lam;
  std::vector< char > act;
  double lb;
  double v;
  bool down;
  double r;
  };
 std::vector< Postponed > later;
 std::size_t next_later = 0;

 auto close_side = [ & ]( int l ) {
  if( side[ l ] == 0 )
   dn_open[ l ] = false;
  else
   up_open[ l ] = false;
  };

 // the depth-first search of the subtree of the node at depth base, which
 // is lvl; with first_full, its first child is not postponed
 auto search = [ & ]( int base , bool first_full ) {
  for( long loop = 0 ; lvl >= base ; ++loop ) {
   if( f_has_sol && ( UB <= lwcut ) )
    return( int( kCutOff ) );
   if( f_nodes >= f_max_nodes )
    return( int( kStopIter ) );
   if( ( ( loop & 63 ) == 0 ) && ( elapsed() > f_max_time ) )
    return( int( kStopTime ) );

   // a new incumbent may prune a node of the path
   if( LB[ lvl ] >= cut() ) {
    LBmin = std::min( LBmin , LB[ lvl ] );
    --lvl;
    continue;
    }

   double r;
   if( ! next_child( lvl , r ) ) {
    --lvl;
    continue;
    }

   make_child( lvl , r );
   const int cl = lvl + 1;
   const double lbc = eval( cl , f_warm_start , ! first_full );
   first_full = false;

   // the side is postponed as a whole, its outer siblings being no less
   // likely to be infeasible
   if( postponed ) {
    later.push_back(
      { lvl , std::vector< double >( fix.begin() , fix.begin() + lvl ) ,
	lam[ lvl ] , act[ lvl ] , LB[ lvl ] , v[ lvl ] , side[ lvl ] == 0 ,
	r } );
    close_side( lvl );
    continue;
    }

   if( cl == n1 ) {
    LBmin = std::min( LBmin , lbc );
    if( accepted || ( lbc >= cut() ) )
     close_side( lvl );
    continue;
    }

   // by convexity, all the outer siblings of a pruned child are pruned
   if( lbc >= cut() ) {
    LBmin = std::min( LBmin , lbc );
    close_side( lvl );
    continue;
    }

   LB[ cl ] = lbc;
   setup_children( cl );
   lvl = cl;
   }
  return( int( kOK ) );
  };

 // the root
 const double lbroot = eval( 0 , false , false );
 f_root_it = f_qp_it;
 if( n1 == 0 )
  LBmin = lbroot;
 else if( lbroot >= cut() )
  LBmin = lbroot;
 else {
  LB[ 0 ] = lbroot;
  setup_children( 0 );
  lvl = 0;
  stop = search( 0 , false );
  }

 // the postponed sides, each from the node it belongs to, rebuilt along
 // its fixings: with an incumbent their nodes are mostly pruned at once
 while( ( stop == kOK ) && ( next_later < later.size() ) ) {
  const Postponed pp = later[ next_later++ ];
  for( int l = 0 ; l < pp.l ; ++l )
   move_to_child( l , pp.fix[ l ] );
  lam[ pp.l ] = pp.lam;
  act[ pp.l ] = pp.act;
  LB[ pp.l ] = pp.lb;
  v[ pp.l ] = pp.v;
  dn_open[ pp.l ] = up_open[ pp.l ] = false;
  if( pp.down ) {
   dn[ pp.l ] = pp.r;
   dn_open[ pp.l ] = true;
   }
  else {
   up[ pp.l ] = pp.r;
   up_open[ pp.l ] = true;
   }
  lvl = pp.l;
  stop = search( pp.l , true );
  }

 // what is still open: the path, and the postponed sides
 if( stop != kOK ) {
  for( int l = 0 ; l <= lvl ; ++l )
   LBmin = std::min( LBmin , LB[ l ] );
  for( std::size_t h = next_later ; h < later.size() ; ++h )
   LBmin = std::min( LBmin , later[ h ].lb );
  }
 else
  lvl = -1;

 // the leaf of the incumbent, solved again from scratch
 if( f_has_sol && f_post_solve && ( n1 < n ) ) {
  std::vector< double > x;
  auto lm = best_lam;
  if( solve_leaf( best_fix , lm , f_post_opt_tol , x ) ) {
   const double val = value_of( x );
   if( val <= UB + prune_eps( UB ) ) {
    UB = val;
    xbest.swap( x );
    }
   }
  }

 if( f_has_sol ) {
  f_UB = UB;
  f_LB = std::min( LBmin , UB );
  f_x.assign( n , 0 );
  for( int k = 0 ; k < n ; ++k )
   f_x[ f_perm[ k ] ] = xbest[ k ];
  }
 else
  f_LB = ( ( lvl < 0 ) && ( ! is_finite( upcut ) ) && ( unsolved == 0 ) )
	   ? Inf< double >()
	   : LBmin;

 // with the whole tree explored, the bound may still leave a gap if some
 // leaf has not been solved
 if( ( lvl < 0 ) && ( stop == kOK ) ) {
  if( f_has_sol )
   stop = ( f_LB >= UB - prune_eps( UB ) ) ? int( kOK ) : int( kLowPrecision );
  else if( unsolved > 0 )
   stop = kError;
  else
   stop = is_finite( upcut ) ? int( kCutOff ) : int( kInfeasible );
  }

 if( f_log && ( log_verbosity > 0 ) )
  *f_log << "FASTQPAMILPSolver: status " << stop << ", nodes " << f_nodes
	 << ", FAST-QPA iterations " << f_qp_it << " (root " << f_root_it
	 << "), rows " << m << ", integer variables " << n1 << " of " << n
	 << ", time " << elapsed() << std::endl;

 return( f_status = stop );

 }  // end( FASTQPAMILPSolver::guts_of_compute )

/*--------------------------------------------------------------------------*/
/*---------------------- METHODS FOR READING RESULTS -----------------------*/
/*--------------------------------------------------------------------------*/

Solver::OFValue FASTQPAMILPSolver::get_lb( void )
{
 if( f_status == kUnEval )
  return( -Inf< OFValue >() );

 // the internal value v is constant_value + v for a minimization and
 // constant_value - v for a maximization
 if( f_sign > 0 )
  return( constant_value + f_LB );
 return( constant_value - f_UB );
 }

/*--------------------------------------------------------------------------*/

Solver::OFValue FASTQPAMILPSolver::get_ub( void )
{
 if( f_status == kUnEval )
  return( Inf< OFValue >() );

 if( f_sign > 0 )
  return( constant_value + f_UB );
 return( constant_value - f_LB );
 }

/*--------------------------------------------------------------------------*/

bool FASTQPAMILPSolver::has_var_solution( void )
{
 return( f_has_sol );
 }

/*--------------------------------------------------------------------------*/

Solver::OFValue FASTQPAMILPSolver::get_var_value( void )
{
 if( ! f_has_sol )
  return( f_sign > 0 ? Inf< OFValue >() : -Inf< OFValue >() );

 return( constant_value + f_sign * f_UB );
 }

/*--------------------------------------------------------------------------*/

void FASTQPAMILPSolver::get_var_solution( Configuration * solc )
{
 if( ! f_has_sol )
  throw( std::logic_error( "FASTQPAMILPSolver::get_var_solution: "
			   "no solution available" ) );

 MILPSolver::write_var_solution( f_x );
 }

/*--------------------------------------------------------------------------*/
/*------------------ METHODS FOR HANDLING THE PARAMETERS -------------------*/
/*--------------------------------------------------------------------------*/

void FASTQPAMILPSolver::set_par( idx_type par , int value )
{
 switch( par ) {
  case( intMaxIter ):
   f_max_nodes = std::max( value , 0 );
   return;
  case( intLogVerb ):
   log_verbosity = value;
   break;
  case( intFQWarmStart ):
   f_warm_start = bool( value );
   return;
  case( intFQEarlyPrune ):
   f_early_prune = bool( value );
   return;
  case( intFQIntBoundsAsRows ):
   if( bool( value ) != f_int_bounds_rows ) {
    f_int_bounds_rows = bool( value );
    f_dirty = true;
    }
   return;
  case( intFQMaxQPIter ):
   f_max_qp_iter = std::max( value , 1 );
   return;
  case( intFQPostSolve ):
   f_post_solve = bool( value );
   return;
  }

 MILPSolver::set_par( par , value );
 }

/*--------------------------------------------------------------------------*/

void FASTQPAMILPSolver::set_par( idx_type par , double value )
{
 switch( par ) {
  case( dblMaxTime ):
   f_max_time = std::max( value , 0.0 );
   return;
  case( dblAbsAcc ):
   f_abs_acc = std::max( value , 0.0 );
   return;
  case( dblRelAcc ):
   f_rel_acc = std::max( value , 0.0 );
   return;
  case( dblUpCutOff ):
   f_up_cutoff = value;
   return;
  case( dblLwCutOff ):
   f_lw_cutoff = value;
   return;
  case( dblFAccSol ):
   if( value != f_feas_tol ) {
    f_feas_tol = std::max( value , 1e-12 );
    f_dirty = true;  // the empty rows are checked with it
    }
   return;
  case( dblFQOptTol ):
   f_opt_tol = std::max( value , 1e-14 );
   return;
  case( dblFQPostOptTol ):
   f_post_opt_tol = std::max( value , 1e-14 );
   return;
  case( dblFQEpsFactor ):
   value = std::clamp( value , 1e-12 , 1 - 1e-12 );
   if( value != f_eps_factor ) {
    f_eps_factor = value;
    f_prep = false;
    }
   return;
  case( dblFQDelta ):
   f_delta = std::clamp( value , 1e-3 , 1 - 1e-3 );
   return;
  case( dblFQGamma ):
   f_gamma = std::clamp( value , 1e-12 , 0.5 - 1e-12 );
   return;
  case( dblFQEta ):
   f_eta = std::max( value , 0.0 );
   return;
  case( dblFQPDTol ):
   value = std::max( value , 0.0 );
   if( value != f_pd_tol ) {
    f_pd_tol = value;
    f_prep = false;
    }
   return;
  }

 MILPSolver::set_par( par , value );
 }

/*--------------------------------------------------------------------------*/

int FASTQPAMILPSolver::get_dflt_int_par( idx_type par ) const
{
 switch( par ) {
  case( intFQWarmStart ):
   return( 1 );
  case( intFQEarlyPrune ):
   return( 1 );
  case( intFQIntBoundsAsRows ):
   return( 0 );
  case( intFQMaxQPIter ):
   return( 1000 );
  case( intFQPostSolve ):
   return( 1 );
  }

 return( MILPSolver::get_dflt_int_par( par ) );
 }

/*--------------------------------------------------------------------------*/

double FASTQPAMILPSolver::get_dflt_dbl_par( idx_type par ) const
{
 switch( par ) {
  case( dblAbsAcc ):
   return( 1e-6 );
  case( dblFAccSol ):
   return( 1e-6 );
  case( dblFQOptTol ):
   return( 1e-5 );
  case( dblFQPostOptTol ):
   return( 1e-9 );
  case( dblFQEpsFactor ):
   return( 0.5 );
  case( dblFQDelta ):
   return( 0.5 );
  case( dblFQGamma ):
   return( 1e-4 );
  case( dblFQEta ):
   return( 1e-9 );
  case( dblFQPDTol ):
   return( 1e-12 );
  }

 return( MILPSolver::get_dflt_dbl_par( par ) );
 }

/*--------------------------------------------------------------------------*/

int FASTQPAMILPSolver::get_int_par( idx_type par ) const
{
 switch( par ) {
  case( intMaxIter ):
   return( f_max_nodes >= Inf< int >() ? Inf< int >() : int( f_max_nodes ) );
  case( intLogVerb ):
   return( log_verbosity );
  case( intFQWarmStart ):
   return( f_warm_start );
  case( intFQEarlyPrune ):
   return( f_early_prune );
  case( intFQIntBoundsAsRows ):
   return( f_int_bounds_rows );
  case( intFQMaxQPIter ):
   return( f_max_qp_iter );
  case( intFQPostSolve ):
   return( f_post_solve );
  }

 return( MILPSolver::get_int_par( par ) );
 }

/*--------------------------------------------------------------------------*/

double FASTQPAMILPSolver::get_dbl_par( idx_type par ) const
{
 switch( par ) {
  case( dblMaxTime ):
   return( f_max_time );
  case( dblAbsAcc ):
   return( f_abs_acc );
  case( dblRelAcc ):
   return( f_rel_acc );
  case( dblUpCutOff ):
   return( f_up_cutoff );
  case( dblLwCutOff ):
   return( f_lw_cutoff );
  case( dblFAccSol ):
   return( f_feas_tol );
  case( dblFQOptTol ):
   return( f_opt_tol );
  case( dblFQPostOptTol ):
   return( f_post_opt_tol );
  case( dblFQEpsFactor ):
   return( f_eps_factor );
  case( dblFQDelta ):
   return( f_delta );
  case( dblFQGamma ):
   return( f_gamma );
  case( dblFQEta ):
   return( f_eta );
  case( dblFQPDTol ):
   return( f_pd_tol );
  }

 return( MILPSolver::get_dbl_par( par ) );
 }

/*--------------------------------------------------------------------------*/

Solver::idx_type
FASTQPAMILPSolver::int_par_str2idx( const std::string & name ) const
{
 auto it = std::find( FQ_int_pars_str.begin() , FQ_int_pars_str.end() , name );
 if( it != FQ_int_pars_str.end() )
  return( intLastAlgParMILP + ( it - FQ_int_pars_str.begin() ) );

 return( MILPSolver::int_par_str2idx( name ) );
 }

/*--------------------------------------------------------------------------*/

const std::string & FASTQPAMILPSolver::int_par_idx2str( idx_type idx ) const
{
 if( ( idx >= intLastAlgParMILP ) && ( idx < intLastAlgParFQ ) )
  return( FQ_int_pars_str[ idx - intLastAlgParMILP ] );

 return( MILPSolver::int_par_idx2str( idx ) );
 }

/*--------------------------------------------------------------------------*/

Solver::idx_type
FASTQPAMILPSolver::dbl_par_str2idx( const std::string & name ) const
{
 auto it = std::find( FQ_dbl_pars_str.begin() , FQ_dbl_pars_str.end() , name );
 if( it != FQ_dbl_pars_str.end() )
  return( dblLastAlgParMILP + ( it - FQ_dbl_pars_str.begin() ) );

 return( MILPSolver::dbl_par_str2idx( name ) );
 }

/*--------------------------------------------------------------------------*/

const std::string & FASTQPAMILPSolver::dbl_par_idx2str( idx_type idx ) const
{
 if( ( idx >= dblLastAlgParMILP ) && ( idx < dblLastAlgParFQ ) )
  return( FQ_dbl_pars_str[ idx - dblLastAlgParMILP ] );

 return( MILPSolver::dbl_par_idx2str( idx ) );
 }

/*--------------------------------------------------------------------------*/
/*--------------------- End File FASTQPAMILPSolver.cpp ---------------------*/
/*--------------------------------------------------------------------------*/
