/*--------------------------------------------------------------------------*/
/*---------------------- File QuadProgMILPSolver.cpp -----------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Implementation of the QuadProgMILPSolver class, which implements the
 * MILPSolver interface for strictly convex Quadratic Programs with linear
 * constraints using the dual active-set solver quadprog.
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

#include <Eigen/Dense>

#include "QuadProgMILPSolver.h"

/*--------------------------------------------------------------------------*/
/*------------------------------- QUADPROG ---------------------------------*/
/*--------------------------------------------------------------------------*/

// quadprog installs no header: this is the declaration of its solve.QP.c,
// which solves min { 1/2 x' G x - a' x : C' x = b (the first meq columns of
// C) , C' x >= b (the others) }, returning 0 at the optimum, 1 if the
// constraints are inconsistent and 2 if G is not positive definite; G is
// given as the inverse of its upper Cholesky factor if factorized != 0, and
// both G and av are overwritten

extern "C" int qpgen2_( double * G , double * av , int n ,
			double * xv , double * lagr , double * obj ,
			double * C , double * bv , int q , int meq ,
			int * iact , int * nact , int * iter ,
			double * work , int factorized );

/*--------------------------------------------------------------------------*/
/*-------------------------------- USING -----------------------------------*/
/*--------------------------------------------------------------------------*/

using namespace SMSpp_di_unipi_it;

/*--------------------------------------------------------------------------*/
/*-------------------------------- CONSTANTS -------------------------------*/
/*--------------------------------------------------------------------------*/

static const std::array< std::string , 1 > QuadProg_int_pars_str = {
 "intQuadProgReuseFactor" };

static const std::array< std::string , 1 > QuadProg_dbl_pars_str = {
 "dblQuadProgReg" };

/*--------------------------------------------------------------------------*/
/*------------------------- FACTORY MANAGEMENT -----------------------------*/
/*--------------------------------------------------------------------------*/

SMSpp_insert_in_factory_cpp_0( QuadProgMILPSolver );

/*--------------------------------------------------------------------------*/
/*--------------------- DERIVED METHODS OF BASE CLASS ----------------------*/
/*--------------------------------------------------------------------------*/

void QuadProgMILPSolver::clear_problem( unsigned int what )
{
 MILPSolver::clear_problem( what );
 d_factored = false;
 d_factH.clear();
 d_n = d_nc = d_meq = 0;
 d_status = kError;
 }

/*--------------------------------------------------------------------------*/

void QuadProgMILPSolver::load_problem( void )
{
 MILPSolver::load_problem();

 if( numquadrows > 0 )
  throw( std::invalid_argument( "QuadProgMILPSolver::load_problem: "
				"quadratic constraints are not supported" ) );

 if( ( int_vars > 0 ) && ( relax_int_vars == 0 ) )
  throw( std::invalid_argument( "QuadProgMILPSolver::load_problem: integer "
				"variables are only supported if relaxed "
				"(intRelaxIntVars)" ) );

 const int n = numcols;
 d_n = n;
 d_nc = numrows;
 d_sign = ( objsense == -1 ) ? -1 : 1;  // quadprog minimizes

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
 std::vector< double > A( std::size_t( d_nc ) * n , 0 );
 for( int j = 0 ; j < n ; ++j )
  for( int k = matbeg[ j ] ; k < matbeg[ j ] + matcnt[ j ] ; ++k )
   A[ std::size_t( matind[ k ] ) * n + j ] += matval[ k ];

 // the sides of the rows
 std::vector< double > lbA( d_nc );
 std::vector< double > ubA( d_nc );
 for( int i = 0 ; i < d_nc ; ++i )
  switch( sense[ i ] ) {
   case( 'L' ): lbA[ i ] = - Inf< double >();
                ubA[ i ] = rhs[ i ];
                break;
   case( 'G' ): lbA[ i ] = rhs[ i ];
                ubA[ i ] = Inf< double >();
                break;
   case( 'E' ): lbA[ i ] = ubA[ i ] = rhs[ i ];
                break;
   case( 'R' ): if( rngval[ i ] > 0 ) {
                 lbA[ i ] = rhs[ i ];
                 ubA[ i ] = rhs[ i ] + rngval[ i ];
                 }
                else {
                 lbA[ i ] = rhs[ i ] + rngval[ i ];
                 ubA[ i ] = rhs[ i ];
                 }
                break;
   default:
    throw( std::invalid_argument( "QuadProgMILPSolver::load_problem: row " +
				  std::to_string( i ) + " has sense " +
				  std::string( 1 , sense[ i ] ) ) );
   }

 // the constraints of quadprog, the equalities first: index j < n of d_y is
 // the bound of variable j, n + i is row i
 d_C.clear();
 d_b.clear();
 d_idx.clear();
 d_sgn.clear();

 auto add = [ & ]( int idx , double sgn , double side ) {
  const auto k = d_C.size();
  d_C.resize( k + n , 0 );
  if( idx < n )
   d_C[ k + idx ] = sgn;
  else
   for( int j = 0 ; j < n ; ++j )
    d_C[ k + j ] = sgn * A[ std::size_t( idx - n ) * n + j ];
  d_b.push_back( sgn * side );
  d_idx.push_back( idx );
  d_sgn.push_back( sgn );
  };

 for( int j = 0 ; j < n ; ++j )
  if( lb[ j ] == ub[ j ] )
   add( j , 1 , lb[ j ] );
 for( int i = 0 ; i < d_nc ; ++i )
  if( lbA[ i ] == ubA[ i ] )
   add( n + i , 1 , lbA[ i ] );
 d_meq = int( d_b.size() );

 for( int j = 0 ; j < n ; ++j )
  if( lb[ j ] != ub[ j ] ) {
   if( std::isfinite( lb[ j ] ) )
    add( j , 1 , lb[ j ] );
   if( std::isfinite( ub[ j ] ) )
    add( j , -1 , ub[ j ] );
   }
 for( int i = 0 ; i < d_nc ; ++i )
  if( lbA[ i ] != ubA[ i ] ) {
   if( std::isfinite( lbA[ i ] ) )
    add( n + i , 1 , lbA[ i ] );
   if( std::isfinite( ubA[ i ] ) )
    add( n + i , -1 , ubA[ i ] );
   }

 d_x.assign( n , 0 );
 d_y.assign( n + d_nc , 0 );
 d_status = kError;

 }  // end( QuadProgMILPSolver::load_problem )

/*--------------------------------------------------------------------------*/
/*-------------------- METHODS FOR SOLVING THE MODEL -----------------------*/
/*--------------------------------------------------------------------------*/

bool QuadProgMILPSolver::factorize( void )
{
 using RowMajor = Eigen::Matrix< double , Eigen::Dynamic , Eigen::Dynamic ,
				 Eigen::RowMajor >;
 Eigen::MatrixXd G = Eigen::Map< const RowMajor >( d_H.data() , d_n , d_n );
 if( d_reg > 0 )
  G.diagonal().array() += d_reg;

 d_factH = d_H;
 d_factReg = d_reg;
 d_factored = false;

 const Eigen::LLT< Eigen::MatrixXd > llt( G );
 if( llt.info() != Eigen::Success )
  return( false );

 // G = R' R with R upper triangular, and quadprog wants R^{-1}
 const Eigen::MatrixXd Rinv = llt.matrixU().solve(
				       Eigen::MatrixXd::Identity( d_n , d_n ) );
 d_Rinv.assign( Rinv.data() , Rinv.data() + std::size_t( d_n ) * d_n );
 return( d_factored = true );
 }

/*--------------------------------------------------------------------------*/

int QuadProgMILPSolver::guts_of_compute( void )
{
 d_fval = 0;

 // with no variable, the problem is feasible if and only if each
 // constraint admits 0
 if( d_n == 0 ) {
  d_status = kOK;
  for( int k = 0 ; k < int( d_b.size() ) ; ++k )
   if( ( d_b[ k ] > 0 ) || ( ( k < d_meq ) && ( d_b[ k ] < 0 ) ) )
    d_status = kInfeasible;
  d_y.assign( d_nc , 0 );
  return( sol_status = d_status );
  }

 // the factor of the previous solve is given if H has not changed since
 if( ! ( d_reuse && d_factored && ( d_factReg == d_reg ) &&
	 ( d_factH == d_H ) ) )
  if( ! factorize() )
   return( sol_status = d_status = kError );

 // quadprog overwrites the factor and the linear term, and wants at least
 // one constraint, which is 0' x >= -1 if there is none
 std::vector< double > G( d_Rinv );
 std::vector< double > a( d_n );
 for( int j = 0 ; j < d_n ; ++j )
  a[ j ] = - d_h[ j ];

 int q = int( d_b.size() );
 std::vector< double > C( d_C );
 std::vector< double > b( d_b );
 if( q == 0 ) {
  C.assign( d_n , 0 );
  b.assign( 1 , -1 );
  q = 1;
  }

 const int r = std::min( d_n , q );
 std::vector< double > work( 2 * std::size_t( d_n ) + 2 * std::size_t( q ) +
			     std::size_t( r ) * ( r + 5 ) / 2 );
 std::vector< double > lagr( q );
 std::vector< int > iact( q );
 int nact = 0;
 std::array< int , 2 > iter{};
 double obj = 0;

 const int ret = qpgen2_( G.data() , a.data() , d_n , d_x.data() ,
			  lagr.data() , &obj , C.data() , b.data() , q ,
			  d_meq , iact.data() , &nact , iter.data() ,
			  work.data() , 1 );

 if( ret == 1 )
  return( sol_status = d_status = kInfeasible );
 if( ret != 0 )
  return( sol_status = d_status = kError );

 d_status = kOK;

 // the bounds are constraints of quadprog, satisfied only up to the
 // accuracy of the dual method on the regularized H: x is projected on
 // them, so that a fixed variable gets exactly its value
 for( int j = 0 ; j < d_n ; ++j )
  d_x[ j ] = std::clamp( d_x[ j ] , double( lb[ j ] ) , double( ub[ j ] ) );

 // the value is computed on H as loaded, not regularized
 for( int i = 0 ; i < d_n ; ++i ) {
  double Hx = 0;
  for( int j = 0 ; j < d_n ; ++j )
   Hx += d_H[ std::size_t( i ) * d_n + j ] * d_x[ j ];
  d_fval += d_x[ i ] * ( d_h[ i ] + Hx / 2 );
  }

 // the multipliers: H x + h = sum_k lagr_k c_k, with c_k = sgn_k a_i for
 // a bound or a row a_i, hence y_i = sum of sgn_k lagr_k over the
 // constraints coming from i, > 0 on an active lower side and < 0 on an
 // active upper one
 d_y.assign( d_n + d_nc , 0 );
 for( int k = 0 ; k < int( d_b.size() ) ; ++k )
  d_y[ d_idx[ k ] ] += d_sgn[ k ] * lagr[ k ];

 return( sol_status = d_status );

 }  // end( QuadProgMILPSolver::guts_of_compute )

/*--------------------------------------------------------------------------*/
/*---------------------- METHODS FOR READING RESULTS -----------------------*/
/*--------------------------------------------------------------------------*/

Solver::OFValue QuadProgMILPSolver::get_lb( void )
{
 if( solved() )
  return( d_sign * d_fval + constant_value );

 if( sol_status == kInfeasible )
  return( objsense == -1 ? - Inf< OFValue >() : Inf< OFValue >() );

 return( - Inf< OFValue >() );
 }

/*--------------------------------------------------------------------------*/

Solver::OFValue QuadProgMILPSolver::get_ub( void )
{
 if( solved() )
  return( d_sign * d_fval + constant_value );

 if( sol_status == kInfeasible )
  return( objsense == -1 ? - Inf< OFValue >() : Inf< OFValue >() );

 return( Inf< OFValue >() );
 }

/*--------------------------------------------------------------------------*/

bool QuadProgMILPSolver::has_var_solution( void )
{
 return( solved() );
 }

/*--------------------------------------------------------------------------*/

Solver::OFValue QuadProgMILPSolver::get_var_value( void )
{
 return( objsense == -1 ? get_lb() : get_ub() );
 }

/*--------------------------------------------------------------------------*/

void QuadProgMILPSolver::get_var_solution( Configuration * solc )
{
 if( ! solved() )
  throw( std::logic_error( "QuadProgMILPSolver::get_var_solution: "
			   "no solution available" ) );

 MILPSolver::write_var_solution( d_x , solc );
 }

/*--------------------------------------------------------------------------*/

bool QuadProgMILPSolver::has_dual_solution( void )
{
 return( solved() );
 }

/*--------------------------------------------------------------------------*/

bool QuadProgMILPSolver::is_dual_feasible( void )
{
 return( solved() );
 }

/*--------------------------------------------------------------------------*/

void QuadProgMILPSolver::get_dual_solution( Configuration * solc )
{
 if( ! solved() )
  throw( std::logic_error( "QuadProgMILPSolver::get_dual_solution: "
			   "no solution available" ) );

 // the multipliers of the problem of the Block are y for a minimization and
 // - y for a maximization, whose Objective quadprog has been given negated
 std::vector< double > rc( d_n );
 for( int j = 0 ; j < d_n ; ++j )
  rc[ j ] = d_sign * d_y[ j ];

 std::vector< double > pi( d_nc );
 for( int i = 0 ; i < d_nc ; ++i )
  pi[ i ] = d_sign * d_y[ d_n + i ];

 MILPSolver::write_dual_solution( pi , rc , solc );
 }

/*--------------------------------------------------------------------------*/
/*------------------ METHODS FOR HANDLING THE PARAMETERS -------------------*/
/*--------------------------------------------------------------------------*/

void QuadProgMILPSolver::set_par( idx_type par , int value )
{
 if( par == intQuadProgReuseFactor ) {
  d_reuse = bool( value );
  return;
  }

 MILPSolver::set_par( par , value );
 }

/*--------------------------------------------------------------------------*/

void QuadProgMILPSolver::set_par( idx_type par , double value )
{
 if( par == dblQuadProgReg ) {
  d_reg = std::max( value , 0.0 );
  return;
  }

 MILPSolver::set_par( par , value );
 }

/*--------------------------------------------------------------------------*/

int QuadProgMILPSolver::get_dflt_int_par( idx_type par ) const
{
 if( par == intQuadProgReuseFactor )
  return( 1 );

 return( MILPSolver::get_dflt_int_par( par ) );
 }

/*--------------------------------------------------------------------------*/

double QuadProgMILPSolver::get_dflt_dbl_par( idx_type par ) const
{
 if( par == dblQuadProgReg )
  return( 0 );

 return( MILPSolver::get_dflt_dbl_par( par ) );
 }

/*--------------------------------------------------------------------------*/

int QuadProgMILPSolver::get_int_par( idx_type par ) const
{
 if( par == intQuadProgReuseFactor )
  return( d_reuse );

 return( MILPSolver::get_int_par( par ) );
 }

/*--------------------------------------------------------------------------*/

double QuadProgMILPSolver::get_dbl_par( idx_type par ) const
{
 if( par == dblQuadProgReg )
  return( d_reg );

 return( MILPSolver::get_dbl_par( par ) );
 }

/*--------------------------------------------------------------------------*/

Solver::idx_type QuadProgMILPSolver::int_par_str2idx(
					      const std::string & name ) const
{
 auto it = std::find( QuadProg_int_pars_str.begin() ,
		      QuadProg_int_pars_str.end() , name );
 if( it != QuadProg_int_pars_str.end() )
  return( intLastAlgParMILP + ( it - QuadProg_int_pars_str.begin() ) );

 return( MILPSolver::int_par_str2idx( name ) );
 }

/*--------------------------------------------------------------------------*/

const std::string & QuadProgMILPSolver::int_par_idx2str( idx_type idx ) const
{
 if( ( idx >= intLastAlgParMILP ) && ( idx < intLastAlgParQuadProg ) )
  return( QuadProg_int_pars_str[ idx - intLastAlgParMILP ] );

 return( MILPSolver::int_par_idx2str( idx ) );
 }

/*--------------------------------------------------------------------------*/

Solver::idx_type QuadProgMILPSolver::dbl_par_str2idx(
					      const std::string & name ) const
{
 auto it = std::find( QuadProg_dbl_pars_str.begin() ,
		      QuadProg_dbl_pars_str.end() , name );
 if( it != QuadProg_dbl_pars_str.end() )
  return( dblLastAlgParMILP + ( it - QuadProg_dbl_pars_str.begin() ) );

 return( MILPSolver::dbl_par_str2idx( name ) );
 }

/*--------------------------------------------------------------------------*/

const std::string & QuadProgMILPSolver::dbl_par_idx2str( idx_type idx ) const
{
 if( ( idx >= dblLastAlgParMILP ) && ( idx < dblLastAlgParQuadProg ) )
  return( QuadProg_dbl_pars_str[ idx - dblLastAlgParMILP ] );

 return( MILPSolver::dbl_par_idx2str( idx ) );
 }

/*--------------------------------------------------------------------------*/
/*-------------------- End File QuadProgMILPSolver.cpp ---------------------*/
/*--------------------------------------------------------------------------*/
