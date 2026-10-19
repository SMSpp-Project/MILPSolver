/*--------------------------------------------------------------------------*/
/*------------------------ File QPMADMILPSolver.cpp ------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Implementation of the QPMADMILPSolver class, which implements the
 * MILPSolver interface for strictly convex Quadratic Programs with linear
 * constraints using the dual active-set solver qpmad.
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

#include "QPMADMILPSolver.h"

#include <qpmad/solver.h>

/*--------------------------------------------------------------------------*/
/*-------------------------------- USING -----------------------------------*/
/*--------------------------------------------------------------------------*/

using namespace SMSpp_di_unipi_it;

/*--------------------------------------------------------------------------*/
/*-------------------------------- CONSTANTS -------------------------------*/
/*--------------------------------------------------------------------------*/

static const std::array< std::string , 1 > QPMAD_int_pars_str = {
 "intQPMADReuseFactor" };

static const std::array< std::string , 1 > QPMAD_dbl_pars_str = {
 "dblQPMADReg" };

/*--------------------------------------------------------------------------*/
/*----------------------------- QPMADData ----------------------------------*/
/*--------------------------------------------------------------------------*/

struct QPMADMILPSolver::QPMADData {
 qpmad::Solver solver;           ///< the solver of qpmad
 Eigen::MatrixXd Hf;             ///< H, then its Cholesky factor
 bool factored = false;          ///< true if Hf is the factor of d_factH
 };

/*--------------------------------------------------------------------------*/
/*------------------------- FACTORY MANAGEMENT -----------------------------*/
/*--------------------------------------------------------------------------*/

SMSpp_insert_in_factory_cpp_0( QPMADMILPSolver );

/*--------------------------------------------------------------------------*/
/*--------------------- CONSTRUCTOR AND DESTRUCTOR -------------------------*/
/*--------------------------------------------------------------------------*/

QPMADMILPSolver::QPMADMILPSolver( void ) : MILPSolver() ,
 d_qp( std::make_unique< QPMADData >() ) ,
 d_tol( qpmad::SolverParameters().tolerance_ ) {}

/*--------------------------------------------------------------------------*/

QPMADMILPSolver::~QPMADMILPSolver() = default;

/*--------------------------------------------------------------------------*/
/*--------------------- DERIVED METHODS OF BASE CLASS ----------------------*/
/*--------------------------------------------------------------------------*/

void QPMADMILPSolver::clear_problem( unsigned int what )
{
 MILPSolver::clear_problem( what );
 d_qp->factored = false;
 d_factH.clear();
 d_n = d_nc = 0;
 d_status = kError;
 }

/*--------------------------------------------------------------------------*/

void QPMADMILPSolver::load_problem( void )
{
 MILPSolver::load_problem();

 if( numquadrows > 0 )
  throw( std::invalid_argument( "QPMADMILPSolver::load_problem: "
				"quadratic constraints are not supported" ) );

 if( ( int_vars > 0 ) && ( relax_int_vars == 0 ) )
  throw( std::invalid_argument( "QPMADMILPSolver::load_problem: integer "
				"variables are only supported if relaxed "
				"(intRelaxIntVars)" ) );

 const int n = numcols;
 d_n = n;
 d_nc = numrows;
 d_sign = ( objsense == -1 ) ? -1 : 1;  // qpmad minimizes

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
    throw( std::invalid_argument( "QPMADMILPSolver::load_problem: row " +
				  std::to_string( i ) + " has sense " +
				  std::string( 1 , sense[ i ] ) ) );
   }

 d_x.assign( n , 0 );
 d_y.assign( n + d_nc , 0 );
 d_status = kError;

 }  // end( QPMADMILPSolver::load_problem )

/*--------------------------------------------------------------------------*/
/*-------------------- METHODS FOR SOLVING THE MODEL -----------------------*/
/*--------------------------------------------------------------------------*/

int QPMADMILPSolver::guts_of_compute( void )
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

 using RowMajor = Eigen::Matrix< double , Eigen::Dynamic , Eigen::Dynamic ,
				 Eigen::RowMajor >;
 const Eigen::Map< const Eigen::VectorXd > h( d_h.data() , d_n );
 const Eigen::Map< const Eigen::VectorXd > lb( d_lb.data() , d_n );
 const Eigen::Map< const Eigen::VectorXd > ub( d_ub.data() , d_n );
 const Eigen::Map< const RowMajor > A( d_A.data() , d_nc , d_n );
 const Eigen::Map< const Eigen::VectorXd > lbA( d_lbA.data() , d_nc );
 const Eigen::Map< const Eigen::VectorXd > ubA( d_ubA.data() , d_nc );

 qpmad::SolverParameters param;
 param.tolerance_ = d_tol;
 param.max_iter_ = ( d_max_iter < Inf< int >() ) ? d_max_iter : -1;

 // the factor of the previous solve is given if H has not changed since,
 // otherwise qpmad factorizes in place H, regularized if so asked
 if( d_reuse && d_qp->factored && ( d_factReg == d_reg ) &&
     ( d_factH == d_H ) )
  param.hessian_type_ = qpmad::SolverParameters::HESSIAN_CHOLESKY_FACTOR;
 else {
  d_qp->Hf = Eigen::Map< const RowMajor >( d_H.data() , d_n , d_n );
  if( d_reg > 0 )
   d_qp->Hf.diagonal().array() += d_reg;
  d_factH = d_H;
  d_factReg = d_reg;
  d_qp->factored = false;
  param.hessian_type_ = qpmad::SolverParameters::HESSIAN_LOWER_TRIANGULAR;
  }

 Eigen::VectorXd x( d_n );
 qpmad::Solver::ReturnStatus ret;
 try {
  ret = d_qp->solver.solve( x , d_qp->Hf , h , lb , ub , A , lbA , ubA ,
			    param );
  }
 catch( std::runtime_error & e ) {
  // qpmad throws if H has no Cholesky factorization, or if it finds the
  // constraints infeasible
  const std::string what = e.what();
  if( what.find( "Cholesky" ) != std::string::npos ) {
   d_qp->factored = false;
   d_status = kError;
   }
  else {
   d_qp->factored = true;
   if( ( what.find( "nfeasible" ) != std::string::npos ) ||
       ( what.find( "nconsistent" ) != std::string::npos ) )
    d_status = kInfeasible;
   else
    d_status = kError;
   }
  return( sol_status = d_status );
  }

 d_qp->factored = true;

 if( ret != qpmad::Solver::OK ) {
  if( ret == qpmad::Solver::MAXIMAL_NUMBER_OF_ITERATIONS )
   d_status = kStopIter;
  else
   d_status = kError;
  return( sol_status = d_status );
  }

 d_status = kOK;
 std::copy( x.data() , x.data() + d_n , d_x.begin() );

 // the value is computed here, qpmad giving none, and on H as loaded
 for( int i = 0 ; i < d_n ; ++i ) {
  double Hx = 0;
  for( int j = 0 ; j < d_n ; ++j )
   Hx += d_H[ std::size_t( i ) * d_n + j ] * d_x[ j ];
  d_fval += d_x[ i ] * ( d_h[ i ] + Hx / 2 );
  }

 compute_multipliers();

 return( sol_status = d_status );

 }  // end( QPMADMILPSolver::guts_of_compute )

/*--------------------------------------------------------------------------*/

void QPMADMILPSolver::compute_multipliers( void )
{
 // the active set: the equalities, which qpmad keeps active from the start
 // and whose multipliers it does not give, and the active inequalities
 std::vector< int > active;
 for( int j = 0 ; j < d_n ; ++j )
  if( std::abs( d_lb[ j ] - d_ub[ j ] ) <= d_tol )
   active.push_back( j );
 for( int i = 0 ; i < d_nc ; ++i )
  if( std::abs( d_lbA[ i ] - d_ubA[ i ] ) <= d_tol )
   active.push_back( d_n + i );

 Eigen::VectorXd dual;
 Eigen::Matrix< qpmad::MatrixIndex , Eigen::Dynamic , 1 > indices;
 Eigen::Matrix< bool , Eigen::Dynamic , 1 > is_lower;
 d_qp->solver.getInequalityDual( dual , indices , is_lower );
 for( Eigen::Index k = 0 ; k < indices.size() ; ++k )
  active.push_back( int( indices( k ) ) );

 std::sort( active.begin() , active.end() );
 active.erase( std::unique( active.begin() , active.end() ) , active.end() );

 // the multipliers y of the active constraints solve, in the least-squares
 // sense, (H + reg I) x + h = sum_k y_k a_k, with a_k the unit vector of a
 // bound or a row of A, which gives y > 0 on an active lower side and y < 0
 // on an active upper one
 Eigen::VectorXd g( d_n );
 for( int i = 0 ; i < d_n ; ++i ) {
  double Hx = d_reg * d_x[ i ];
  for( int j = 0 ; j < d_n ; ++j )
   Hx += d_H[ std::size_t( i ) * d_n + j ] * d_x[ j ];
  g( i ) = Hx + d_h[ i ];
  }

 d_y.assign( d_n + d_nc , 0 );
 if( active.empty() )
  return;

 Eigen::MatrixXd M = Eigen::MatrixXd::Zero( d_n , active.size() );
 for( std::size_t k = 0 ; k < active.size() ; ++k )
  if( active[ k ] < d_n )
   M( active[ k ] , k ) = 1;
  else
   for( int j = 0 ; j < d_n ; ++j )
    M( j , k ) = d_A[ std::size_t( active[ k ] - d_n ) * d_n + j ];

 const Eigen::VectorXd y = M.colPivHouseholderQr().solve( g );
 for( std::size_t k = 0 ; k < active.size() ; ++k )
  d_y[ active[ k ] ] = y( k );
 }

/*--------------------------------------------------------------------------*/
/*---------------------- METHODS FOR READING RESULTS -----------------------*/
/*--------------------------------------------------------------------------*/

Solver::OFValue QPMADMILPSolver::get_lb( void )
{
 if( solved() )
  return( d_sign * d_fval + constant_value );

 if( sol_status == kInfeasible )
  return( objsense == -1 ? - Inf< OFValue >() : Inf< OFValue >() );

 return( - Inf< OFValue >() );
 }

/*--------------------------------------------------------------------------*/

Solver::OFValue QPMADMILPSolver::get_ub( void )
{
 if( solved() )
  return( d_sign * d_fval + constant_value );

 if( sol_status == kInfeasible )
  return( objsense == -1 ? - Inf< OFValue >() : Inf< OFValue >() );

 return( Inf< OFValue >() );
 }

/*--------------------------------------------------------------------------*/

bool QPMADMILPSolver::has_var_solution( void )
{
 return( solved() );
 }

/*--------------------------------------------------------------------------*/

Solver::OFValue QPMADMILPSolver::get_var_value( void )
{
 return( objsense == -1 ? get_lb() : get_ub() );
 }

/*--------------------------------------------------------------------------*/

void QPMADMILPSolver::get_var_solution( Configuration * solc )
{
 if( ! solved() )
  throw( std::logic_error( "QPMADMILPSolver::get_var_solution: "
			   "no solution available" ) );

 MILPSolver::write_var_solution( d_x );
 }

/*--------------------------------------------------------------------------*/

bool QPMADMILPSolver::has_dual_solution( void )
{
 return( solved() );
 }

/*--------------------------------------------------------------------------*/

bool QPMADMILPSolver::is_dual_feasible( void )
{
 return( solved() );
 }

/*--------------------------------------------------------------------------*/

void QPMADMILPSolver::get_dual_solution( Configuration * solc )
{
 if( ! solved() )
  throw( std::logic_error( "QPMADMILPSolver::get_dual_solution: "
			   "no solution available" ) );

 // the multipliers of the problem of the Block are y for a minimization and
 // - y for a maximization, whose Objective qpmad has been given negated
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

void QPMADMILPSolver::set_par( idx_type par , int value )
{
 switch( par ) {
  case( intMaxIter ): d_max_iter = std::max( value , 0 ); return;
  case( intQPMADReuseFactor ): d_reuse = bool( value ); return;
  }

 MILPSolver::set_par( par , value );
 }

/*--------------------------------------------------------------------------*/

void QPMADMILPSolver::set_par( idx_type par , double value )
{
 switch( par ) {
  case( dblMaxTime ): d_max_time = std::max( value , 0.0 ); return;
  case( dblFAccSol ): d_tol = value; return;
  case( dblQPMADReg ): d_reg = std::max( value , 0.0 ); return;
  }

 MILPSolver::set_par( par , value );
 }

/*--------------------------------------------------------------------------*/

int QPMADMILPSolver::get_dflt_int_par( idx_type par ) const
{
 if( par == intQPMADReuseFactor )
  return( 1 );

 return( MILPSolver::get_dflt_int_par( par ) );
 }

/*--------------------------------------------------------------------------*/

double QPMADMILPSolver::get_dflt_dbl_par( idx_type par ) const
{
 switch( par ) {
  case( dblFAccSol ):  return( qpmad::SolverParameters().tolerance_ );
  case( dblQPMADReg ): return( 0 );
  }

 return( MILPSolver::get_dflt_dbl_par( par ) );
 }

/*--------------------------------------------------------------------------*/

int QPMADMILPSolver::get_int_par( idx_type par ) const
{
 switch( par ) {
  case( intMaxIter ):          return( d_max_iter );
  case( intQPMADReuseFactor ): return( d_reuse );
  }

 return( MILPSolver::get_int_par( par ) );
 }

/*--------------------------------------------------------------------------*/

double QPMADMILPSolver::get_dbl_par( idx_type par ) const
{
 switch( par ) {
  case( dblMaxTime ):  return( d_max_time );
  case( dblFAccSol ):  return( d_tol );
  case( dblQPMADReg ): return( d_reg );
  }

 return( MILPSolver::get_dbl_par( par ) );
 }

/*--------------------------------------------------------------------------*/

Solver::idx_type QPMADMILPSolver::int_par_str2idx(
					      const std::string & name ) const
{
 auto it = std::find( QPMAD_int_pars_str.begin() ,
		      QPMAD_int_pars_str.end() , name );
 if( it != QPMAD_int_pars_str.end() )
  return( intLastAlgParMILP + ( it - QPMAD_int_pars_str.begin() ) );

 return( MILPSolver::int_par_str2idx( name ) );
 }

/*--------------------------------------------------------------------------*/

const std::string & QPMADMILPSolver::int_par_idx2str( idx_type idx ) const
{
 if( ( idx >= intLastAlgParMILP ) && ( idx < intLastAlgParQPMAD ) )
  return( QPMAD_int_pars_str[ idx - intLastAlgParMILP ] );

 return( MILPSolver::int_par_idx2str( idx ) );
 }

/*--------------------------------------------------------------------------*/

Solver::idx_type QPMADMILPSolver::dbl_par_str2idx(
					      const std::string & name ) const
{
 auto it = std::find( QPMAD_dbl_pars_str.begin() ,
		      QPMAD_dbl_pars_str.end() , name );
 if( it != QPMAD_dbl_pars_str.end() )
  return( dblLastAlgParMILP + ( it - QPMAD_dbl_pars_str.begin() ) );

 return( MILPSolver::dbl_par_str2idx( name ) );
 }

/*--------------------------------------------------------------------------*/

const std::string & QPMADMILPSolver::dbl_par_idx2str( idx_type idx ) const
{
 if( ( idx >= dblLastAlgParMILP ) && ( idx < dblLastAlgParQPMAD ) )
  return( QPMAD_dbl_pars_str[ idx - dblLastAlgParMILP ] );

 return( MILPSolver::dbl_par_idx2str( idx ) );
 }

/*--------------------------------------------------------------------------*/
/*--------------------- End File QPMADMILPSolver.cpp -----------------------*/
/*--------------------------------------------------------------------------*/
