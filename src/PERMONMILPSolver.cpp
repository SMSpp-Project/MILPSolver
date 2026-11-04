/*--------------------------------------------------------------------------*/
/*------------------------ File PERMONMILPSolver.cpp -----------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Implementation of the PERMONMILPSolver class, which implements the
 * MILPSolver interface for convex Quadratic Programs with linear
 * constraints using the solvers of PermonQP.
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
#include <ostream>
#include <stdexcept>
#include <string>
#include <type_traits>

#include <Eigen/OrderingMethods>
#include <Eigen/SparseCore>
#include <Eigen/SparseQR>

#include <permonqps.h>

#include "PERMONMILPSolver.h"

/*--------------------------------------------------------------------------*/
/*-------------------------------- USING -----------------------------------*/
/*--------------------------------------------------------------------------*/

using namespace SMSpp_di_unipi_it;

static_assert( std::is_same_v< PetscScalar , double > ,
	       "PERMONMILPSolver needs PETSc built with real double scalars" );

/*--------------------------------------------------------------------------*/
/*-------------------------------- CONSTANTS -------------------------------*/
/*--------------------------------------------------------------------------*/

static const std::array< std::string , 1 > PERMON_int_pars_str = {
 "intPERMONWarmStart"
 };

static const std::array< std::string , 1 > PERMON_dbl_pars_str = {
 "dblPERMONRelTol"
 };

// the defaults of dblFAccSol and dblPERMONRelTol, tighter than the 1e-50
// and 1e-5 of PETSc, the latter leaving errors of that order relative to
// the norm of the linear term in the solution

static constexpr double PERMON_dflt_atol = 1e-9;
static constexpr double PERMON_dflt_rtol = 1e-9;

// the iteration limit of PermonQP when intMaxIter is not set, that of PETSc

static constexpr int PERMON_dflt_max_iter = 10000;

// the default of dblKKTTol: on a problem whose Hessian is singular along a
// direction of free variables, MPGP may drift along it to a point of norm
// 1e13 and stop there as converged, with a value that is only rounding
// noise; the check of the optimality conditions of MILPSolver turns such a
// solve into kError

static constexpr double PERMON_dflt_kkt_tol = 1e-6;

/*--------------------------------------------------------------------------*/
/*----------------------------- PERMONData ---------------------------------*/
/*--------------------------------------------------------------------------*/

// the objects of PETSc and PermonQP of one solve

struct PERMONMILPSolver::PERMONData
{

 Mat H = nullptr;    ///< Hessian, on the variables and the slacks
 Mat B = nullptr;    ///< [ A 0 ; C -I ]
 Vec rhs = nullptr;  ///< - g (PermonQP minimizes x' H x / 2 - rhs' x)
 Vec c = nullptr;    ///< [ b ; 0 ]
 Vec lb = nullptr;   ///< lower bounds of the variables and slacks
 Vec ub = nullptr;   ///< upper bounds of the variables and slacks
 Vec x = nullptr;    ///< the start, then the solution
 QP qp = nullptr;
 QPS qps = nullptr;

 // the data of the problem, kept between a load and the solve
 Eigen::SparseMatrix< double , Eigen::RowMajor , PetscInt > Hs , Bs;
 std::vector< double > vrhs , vc , vlb , vub;

 void free( void ) {
  QPSDestroy( &qps );
  QPDestroy( &qp );
  VecDestroy( &x );
  VecDestroy( &ub );
  VecDestroy( &lb );
  VecDestroy( &c );
  VecDestroy( &rhs );
  MatDestroy( &B );
  MatDestroy( &H );
  }

 ~PERMONData() { free(); }
 };

/*--------------------------------------------------------------------------*/
/*---------------------------- LOCAL FUNCTIONS -----------------------------*/
/*--------------------------------------------------------------------------*/

namespace
{

// the message of the first error of PETSc since the last solve began

std::string permon_error;

// the handler of the errors of PETSc: it keeps the message of the error
// where it arises and returns its code, rather than printing it and
// aborting

PetscErrorCode permon_error_handler( MPI_Comm , int line , const char * fun ,
				     const char * file , PetscErrorCode n ,
				     PetscErrorType p , const char * mess ,
				     void * )
{
 if( ( p == PETSC_ERROR_INITIAL ) && permon_error.empty() )
  permon_error = std::string( fun ) + "() at " + file + ":" +
		 std::to_string( line ) + ": " + ( mess ? mess : "" );
 return( n );
 }

// starts PETSc and PERMON the first time, unless the program has done it;
// the errors of PETSc then go to permon_error_handler(), PETSc does not
// catch the signals, and up to PERMON 3.26 (whose version is that of
// PETSc) the expansion step of MPGP (also the inner solver of SMALXE) is
// projcg unless the options of PETSc say otherwise: there, a proportioning
// step of the default one may leave the bounds, and MPGP may stop as
// converged at a point out of them

void permon_initialize( void )
{
 static bool done = false;
 if( done )
  return;

 if( PermonInitialize( nullptr , nullptr , nullptr , nullptr ) )
  throw( std::runtime_error( "PERMONMILPSolver: PERMON cannot be started" ) );
 PetscPushErrorHandler( permon_error_handler , nullptr );
 PetscPopSignalHandler();
#if PETSC_VERSION_LT( 3 , 27 , 0 )
 for( const char * opt :
      { "-qps_mpgp_expansion_type" , "-smalxe_qps_mpgp_expansion_type" } ) {
  PetscBool set = PETSC_FALSE;
  PetscOptionsHasName( nullptr , nullptr , opt , &set );
  if( ! set )
   PetscOptionsSetValue( nullptr , opt , "projcg" );
  }
#endif
 done = true;
 }

// a sequential matrix of PETSc with the entries of a row-major one

PetscErrorCode
make_mat( const Eigen::SparseMatrix< double , Eigen::RowMajor ,
	  PetscInt > & S ,
	  Mat * M )
{
 PetscFunctionBeginUser;
 PetscCall( MatCreate( PETSC_COMM_SELF , M ) );
 PetscCall( MatSetSizes( *M , S.rows() , S.cols() , S.rows() , S.cols() ) );
 PetscCall( MatSetType( *M , MATSEQAIJ ) );
 PetscCall( MatSeqAIJSetPreallocationCSR( *M , S.outerIndexPtr() ,
					  S.innerIndexPtr() , S.valuePtr() ) );
 PetscFunctionReturn( PETSC_SUCCESS );
 }

// a sequential vector of PETSc with the given entries

PetscErrorCode make_vec( const std::vector< double > & v , Vec * V )
{
 PetscFunctionBeginUser;
 PetscCall( VecCreateSeq( PETSC_COMM_SELF , PetscInt( v.size() ) , V ) );
 PetscScalar * a;
 PetscCall( VecGetArray( *V , &a ) );
 std::copy( v.begin() , v.end() , a );
 PetscCall( VecRestoreArray( *V , &a ) );
 PetscFunctionReturn( PETSC_SUCCESS );
 }

// copies a vector of PETSc

PetscErrorCode get_vec( Vec V , std::vector< double > & v )
{
 PetscFunctionBeginUser;
 PetscInt n;
 PetscCall( VecGetSize( V , &n ) );
 const PetscScalar * a;
 PetscCall( VecGetArrayRead( V , &a ) );
 v.assign( a , a + n );
 PetscCall( VecRestoreArrayRead( V , &a ) );
 PetscFunctionReturn( PETSC_SUCCESS );
 }

// the point of [ l , u ] closest to v

double clamp( double v , double l , double u )
{
 return( std::max( l , std::min( v , u ) ) );
 }

}  // end( namespace )

/*--------------------------------------------------------------------------*/
/*------------------------- FACTORY MANAGEMENT -----------------------------*/
/*--------------------------------------------------------------------------*/

SMSpp_insert_in_factory_cpp_0( PERMONMILPSolver );

/*--------------------------------------------------------------------------*/
/*--------------------- CONSTRUCTOR AND DESTRUCTOR -------------------------*/
/*--------------------------------------------------------------------------*/

PERMONMILPSolver::PERMONMILPSolver( void )
    : MILPSolver() , d_pd( std::make_unique< PERMONData >() ) ,
      d_atol( PERMON_dflt_atol ) , d_rtol( PERMON_dflt_rtol )
{
 MILPSolver::set_par( idx_type( dblKKTTol ) , PERMON_dflt_kkt_tol );
 }

/*--------------------------------------------------------------------------*/

PERMONMILPSolver::~PERMONMILPSolver() = default;

/*--------------------------------------------------------------------------*/
/*--------------------- DERIVED METHODS OF BASE CLASS ----------------------*/
/*--------------------------------------------------------------------------*/

void PERMONMILPSolver::clear_problem( unsigned int what )
{
 MILPSolver::clear_problem( what );
 d_pd->free();
 d_x.clear();
 d_y.clear();
 d_z.clear();
 d_n = d_neq = d_nin = d_nfix = 0;
 d_row_pos.clear();
 d_fix_row.clear();
 d_status = 0;
 }

/*--------------------------------------------------------------------------*/

void PERMONMILPSolver::load_problem( void )
{
 // the layout of the previous problem, whose solution is moved to the new
 // one in the next solve [see previous_column_index()]
 d_prev_row_pos = d_row_pos;
 d_prev_n = d_n;

 // the solution, which clear_problem() of MILPSolver::load_problem() would
 // throw away, is kept for the warm start of the next solve
 auto x = std::move( d_x );

 MILPSolver::load_problem();
 d_pd->free();

 d_x = std::move( x );

 if( numquadrows > 0 )
  throw( std::invalid_argument( "PERMONMILPSolver::load_problem: "
				"quadratic constraints are not supported" ) );

 if( ( int_vars > 0 ) && ( relax_int_vars == 0 ) )
  throw( std::invalid_argument( "PERMONMILPSolver::load_problem: integer "
				"variables are only supported if relaxed "
				"(intRelaxIntVars)" ) );

 const int n = numcols;
 d_n = n;
 d_sign = ( objsense == -1 ) ? -1 : 1;  // PermonQP minimizes

 // the sides of the rows: those whose two sides coincide are equalities,
 // all the others get a slack
 std::vector< double > lo( numrows ) , hi( numrows );
 d_row_pos.resize( numrows );
 d_neq = d_nin = 0;
 for( int i = 0 ; i < numrows ; ++i ) {
  switch( sense[ i ] ) {
   case( 'L' ):
    lo[ i ] = -Inf< double >();
    hi[ i ] = rhs[ i ];
    break;
   case( 'G' ):
    lo[ i ] = rhs[ i ];
    hi[ i ] = Inf< double >();
    break;
   case( 'E' ):
    lo[ i ] = hi[ i ] = rhs[ i ];
    break;
   case( 'R' ):
    if( rngval[ i ] > 0 ) {
     lo[ i ] = rhs[ i ];
     hi[ i ] = rhs[ i ] + rngval[ i ];
     }
    else {
     lo[ i ] = rhs[ i ] + rngval[ i ];
     hi[ i ] = rhs[ i ];
     }
    break;
   default:
    throw( std::invalid_argument( "PERMONMILPSolver::load_problem: row " +
				  std::to_string( i ) + " has sense " +
				  std::string( 1 , sense[ i ] ) ) );
   }
  if( lo[ i ] == hi[ i ] )
   d_row_pos[ i ] = d_neq++;
  else
   d_row_pos[ i ] = -( ++d_nin );
  }

 // the variables with equal bounds, which MPGP cannot handle in its box,
 // are fixed by an equality x_j = l_j and have no bounds
 d_fix_row.assign( n , -1 );
 d_nfix = 0;
 for( int j = 0 ; j < n ; ++j )
  if( lb[ j ] == ub[ j ] )
   d_fix_row[ j ] = d_neq + d_nfix++;

 const int nt = n + d_nin;       // the variables and the slacks
 const int ns = d_neq + d_nfix;  // the first row with a slack
 const int m = ns + d_nin;       // the equalities of PermonQP
 auto & pd = *d_pd;

 // the objective: c' x + sum_j q_j x_j^2 + sum_k ndq_k x_i x_j is
 // 1/2 x' H x + g' x with H_jj = 2 q_j and H_ij = H_ji = ndq_k, and
 // PermonQP is given rhs = - g; the slacks are not in the objective
 std::vector< Eigen::Triplet< double , PetscInt > > trp;
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
 pd.Hs.resize( nt , nt );
 pd.Hs.setFromTriplets( trp.begin() , trp.end() );
 pd.Hs.prune( 0.0 );
 pd.Hs.makeCompressed();

 pd.vrhs.assign( nt , 0 );
 for( int j = 0 ; j < n ; ++j )
  pd.vrhs[ j ] = -d_sign * objective[ j ];

 // the rows: [ A 0 ] x = b and [ C -I ] ( x , s ) = 0
 trp.clear();
 for( int j = 0 ; j < n ; ++j )
  for( int k = matbeg[ j ] ; k < matbeg[ j ] + matcnt[ j ] ; ++k ) {
   const int p = d_row_pos[ matind[ k ] ];
   trp.emplace_back( p >= 0 ? p : ns - p - 1 , j , matval[ k ] );
   }
 for( int j = 0 ; j < n ; ++j )
  if( d_fix_row[ j ] >= 0 )
   trp.emplace_back( d_fix_row[ j ] , j , 1.0 );
 for( int i = 0 ; i < d_nin ; ++i )
  trp.emplace_back( ns + i , n + i , -1.0 );
 pd.Bs.resize( m , nt );
 pd.Bs.setFromTriplets( trp.begin() , trp.end() );
 pd.Bs.makeCompressed();

 pd.vc.assign( m , 0 );
 pd.vlb.resize( nt );
 pd.vub.resize( nt );
 for( int j = 0 ; j < n ; ++j )
  if( d_fix_row[ j ] >= 0 ) {
   pd.vc[ d_fix_row[ j ] ] = lb[ j ];
   pd.vlb[ j ] = PETSC_NINFINITY;
   pd.vub[ j ] = PETSC_INFINITY;
   }
  else {
   pd.vlb[ j ] = lb[ j ] > -Inf< double >() ? lb[ j ] : PETSC_NINFINITY;
   pd.vub[ j ] = ub[ j ] < Inf< double >() ? ub[ j ] : PETSC_INFINITY;
   }
 for( int i = 0 ; i < numrows ; ++i )
  if( const int p = d_row_pos[ i ]; p >= 0 )
   pd.vc[ p ] = lo[ i ];
  else {
   pd.vlb[ n - p - 1 ] =
     lo[ i ] > -Inf< double >() ? lo[ i ] : PETSC_NINFINITY;
   pd.vub[ n - p - 1 ] = hi[ i ] < Inf< double >() ? hi[ i ] : PETSC_INFINITY;
   }

 d_status = 0;
 d_remap = true;

 }  // end( PERMONMILPSolver::load_problem )

/*--------------------------------------------------------------------------*/
/*-------------------- METHODS FOR SOLVING THE MODEL -----------------------*/
/*--------------------------------------------------------------------------*/

int PERMONMILPSolver::guts_of_compute( void )
{
 auto & pd = *d_pd;
 const int nt = d_n + d_nin;
 const int ns = d_neq + d_nfix;
 const int m = ns + d_nin;

 // a problem with no variables is solved here, PermonQP needing some
 if( nt == 0 ) {
  bool feas = true;
  for( int i = 0 ; i < m ; ++i )
   feas &= ( pd.vc[ i ] == 0 );
  d_status = feas ? int( KSP_CONVERGED_ATOL ) : int( KSP_DIVERGED_BREAKDOWN );
  d_x.clear();
  d_y.assign( m , 0 );
  d_z.clear();
  d_fval = 0;
  return( sol_status = feas ? kOK : kInfeasible );
  }

 // the solution of the previous problem is moved to the columns that are
 // still there, the new ones at 0
 if( d_remap && ( ! d_x.empty() ) && ( int( d_x.size() ) == d_prev_n ) ) {
  std::vector< double > x( d_n , 0 );
  for( int j = 0 ; j < d_n ; ++j )
   if( const int pj = previous_column_index( j ); pj >= 0 )
    x[ j ] = d_x[ pj ];
  d_x = std::move( x );
  }
 d_remap = false;

 // the start: the previous solution, or 0, in the bounds, and each slack
 // at the value of its row there
 std::vector< double > x0( nt , 0 );
 if( d_warm_start && ( int( d_x.size() ) == d_n ) )
  std::copy( d_x.begin() , d_x.end() , x0.begin() );
 for( int j = 0 ; j < d_n ; ++j )
  x0[ j ] = d_fix_row[ j ] >= 0 ? pd.vc[ d_fix_row[ j ] ]
				: clamp( x0[ j ] , pd.vlb[ j ] , pd.vub[ j ] );
 for( int i = 0 ; i < d_nin ; ++i ) {
  double r = 0;
  for( decltype( pd.Bs )::InnerIterator it( pd.Bs , ns + i ) ; it ; ++it )
   if( it.col() < d_n )
    r += it.value() * x0[ it.col() ];
  x0[ d_n + i ] = clamp( r , pd.vlb[ d_n + i ] , pd.vub[ d_n + i ] );
  }

 permon_initialize();

 const bool hasbox =
   std::any_of( pd.vlb.begin() , pd.vlb.end() ,
		[]( double v ) { return( v > PETSC_NINFINITY ); } ) ||
   std::any_of( pd.vub.begin() , pd.vub.end() ,
		[]( double v ) { return( v < PETSC_INFINITY ); } );

 // the objects of PETSc and PermonQP, made anew at each solve
 auto setup = [ & ]( void ) -> PetscErrorCode {
  PetscFunctionBeginUser;
  pd.free();
  PetscCall( make_mat( pd.Hs , &pd.H ) );
  PetscCall( make_vec( pd.vrhs , &pd.rhs ) );
  PetscCall( make_vec( x0 , &pd.x ) );
  PetscCall( QPCreate( PETSC_COMM_SELF , &pd.qp ) );
  PetscCall( QPSetOperator( pd.qp , pd.H ) );
  PetscCall( QPSetRhs( pd.qp , pd.rhs ) );
  PetscCall( QPSetInitialVector( pd.qp , pd.x ) );
  if( m > 0 ) {
   PetscCall( make_mat( pd.Bs , &pd.B ) );
   PetscCall( make_vec( pd.vc , &pd.c ) );
   PetscCall( QPSetEq( pd.qp , pd.B , pd.c ) );
   }
  if( hasbox ) {
   PetscCall( make_vec( pd.vlb , &pd.lb ) );
   PetscCall( make_vec( pd.vub , &pd.ub ) );
   PetscCall( QPSetBox( pd.qp , nullptr , pd.lb , pd.ub ) );
   }
  PetscCall( QPSetFromOptions( pd.qp ) );
  PetscCall( QPSCreate( PETSC_COMM_SELF , &pd.qps ) );
  PetscCall( QPSSetQP( pd.qps , pd.qp ) );
  PetscCall( QPSSetFromOptions( pd.qps ) );
  // the tolerances, given again after the set up, which gives those of the
  // method (the inner solver of SMALXE takes the first ones); no divergence
  // test, which stops a solve whose first residual is small (e.g., one
  // started from a nearly optimal point) as soon as a later one is some
  // orders of magnitude larger
  const PetscInt maxit =
    d_max_iter < Inf< int >() ? d_max_iter : PERMON_dflt_max_iter;
  PetscCall(
    QPSSetTolerances( pd.qps , d_rtol , d_atol , PETSC_MAX_REAL , maxit ) );
  PetscCall( QPSSetUp( pd.qps ) );
  PetscCall(
    QPSSetTolerances( pd.qps , d_rtol , d_atol , PETSC_MAX_REAL , maxit ) );
  PetscFunctionReturn( PETSC_SUCCESS );
  };

 // the solve, the solution and the multipliers
 auto solve = [ & ]( void ) -> PetscErrorCode {
  PetscFunctionBeginUser;
  PetscCall( QPSSolve( pd.qps ) );
  KSPConvergedReason reason;
  PetscCall( QPSGetConvergedReason( pd.qps , &reason ) );
  d_status = int( reason );
  PetscInt its;
  PetscCall( QPSGetIterationNumber( pd.qps , &its ) );
  elapsed_iter += long( its );
  if( d_status <= 0 )
   PetscFunctionReturn( PETSC_SUCCESS );

  Vec xs;
  PetscCall( QPGetSolutionVector( pd.qp , &xs ) );
  std::vector< double > xt;
  PetscCall( get_vec( xs , xt ) );
  d_x.assign( xt.begin() , xt.begin() + d_n );

  // the multipliers z = lub - llb of the bounds
  d_z.assign( nt , 0 );
  if( hasbox ) {
   PetscCall( QPComputeMissingBoxMultipliers( pd.qp ) );
   QPC qpc;
   PetscCall( QPGetQPC( pd.qp , &qpc ) );
   Vec llb = nullptr , lub = nullptr;
   if( qpc )
    PetscCall( QPCBoxGetMultipliers( qpc , &llb , &lub ) );
   std::vector< double > v;
   if( lub ) {
    PetscCall( get_vec( lub , v ) );
    for( int j = 0 ; j < nt ; ++j )
     d_z[ j ] += v[ j ];
    }
   if( llb ) {
    PetscCall( get_vec( llb , v ) );
    for( int j = 0 ; j < nt ; ++j )
     d_z[ j ] -= v[ j ];
    }
   }

  // PermonQP keeps the multipliers lambda of the rows only as B' lambda,
  // hence they are computed here from H x - rhs + B' lambda + z = 0: the
  // row C_i x - s_i = 0 of a slack has lambda_i = z of the slack, and the
  // equalities (those of the rows and of the fixed variables) the least
  // squares solution of A' lambda_A = rhs - H x - z - C' lambda_C on x
  d_y.assign( m , 0 );
  for( int i = 0 ; i < d_nin ; ++i )
   d_y[ ns + i ] = d_z[ d_n + i ];
  if( ns > 0 ) {
   Eigen::VectorXd v( d_n );
   for( int j = 0 ; j < d_n ; ++j )
    v[ j ] = pd.vrhs[ j ] - d_z[ j ];
   for( int i = 0 ; i < d_n ; ++i )
    for( decltype( pd.Hs )::InnerIterator it( pd.Hs , i ) ; it ; ++it )
     v[ i ] -= it.value() * xt[ it.col() ];
   std::vector< Eigen::Triplet< double > > trp;
   for( int i = 0 ; i < m ; ++i )
    for( decltype( pd.Bs )::InnerIterator it( pd.Bs , i ) ; it ; ++it )
     if( it.col() < d_n ) {
      if( i < ns )
       trp.emplace_back( int( it.col() ) , i , it.value() );
      else
       v[ it.col() ] -= it.value() * d_y[ i ];
      }
   Eigen::SparseMatrix< double > At( d_n , ns );
   At.setFromTriplets( trp.begin() , trp.end() );
   At.makeCompressed();
   Eigen::SparseQR< Eigen::SparseMatrix< double > ,
		    Eigen::COLAMDOrdering< int > >
     qr( At );
   if( qr.info() != Eigen::Success )
    PetscFunctionReturn( PETSC_ERR_LIB );
   const Eigen::VectorXd y = qr.solve( v );
   for( int i = 0 ; i < ns ; ++i )
    d_y[ i ] = y[ i ];
   }

  d_fval = 0;
  for( int j = 0 ; j < d_n ; ++j )
   d_fval -= pd.vrhs[ j ] * d_x[ j ];
  for( int i = 0 ; i < pd.Hs.outerSize() ; ++i )
   if( i < d_n )
    for( decltype( pd.Hs )::InnerIterator it( pd.Hs , i ) ; it ; ++it )
     if( it.col() < d_n )
      d_fval += d_x[ i ] * it.value() * d_x[ it.col() ] / 2;
  PetscFunctionReturn( PETSC_SUCCESS );
  };

 d_status = 0;
 permon_error.clear();
 if( setup() || solve() ) {
  if( f_log && ( log_verbosity > 0 ) )
   *f_log << "PERMONMILPSolver: error of PETSc in " << permon_error
	  << std::endl;
  d_status = 0;
  pd.free();
  return( sol_status = kError );
  }

 if( d_status > 0 )
  sol_status = kOK;
 else
  sol_status = ( d_status == int( KSP_DIVERGED_ITS ) ) ? kStopIter : kError;

 return( sol_status );

 }  // end( PERMONMILPSolver::guts_of_compute )

/*--------------------------------------------------------------------------*/

bool PERMONMILPSolver::solved( void ) const
{
 return( d_status > 0 );
 }

/*--------------------------------------------------------------------------*/
/*---------------------- METHODS FOR READING RESULTS -----------------------*/
/*--------------------------------------------------------------------------*/

Solver::OFValue PERMONMILPSolver::get_lb( void )
{
 if( solved() )
  return( d_sign * d_fval + constant_value );

 if( sol_status == kInfeasible )
  return( objsense == -1 ? -Inf< OFValue >() : Inf< OFValue >() );

 return( -Inf< OFValue >() );
 }

/*--------------------------------------------------------------------------*/

Solver::OFValue PERMONMILPSolver::get_ub( void )
{
 if( solved() )
  return( d_sign * d_fval + constant_value );

 if( sol_status == kInfeasible )
  return( objsense == -1 ? -Inf< OFValue >() : Inf< OFValue >() );

 return( Inf< OFValue >() );
 }

/*--------------------------------------------------------------------------*/

bool PERMONMILPSolver::has_var_solution( void )
{
 return( solved() );
 }

/*--------------------------------------------------------------------------*/

Solver::OFValue PERMONMILPSolver::get_var_value( void )
{
 return( objsense == -1 ? get_lb() : get_ub() );
 }

/*--------------------------------------------------------------------------*/

void PERMONMILPSolver::get_var_solution( Configuration * solc )
{
 if( ! solved() )
  throw( std::logic_error( "PERMONMILPSolver::get_var_solution: "
			   "no solution available" ) );

 MILPSolver::write_var_solution( d_x , solc );
 }

/*--------------------------------------------------------------------------*/

bool PERMONMILPSolver::has_dual_solution( void )
{
 return( solved() );
 }

/*--------------------------------------------------------------------------*/

bool PERMONMILPSolver::is_dual_feasible( void )
{
 return( solved() );
 }

/*--------------------------------------------------------------------------*/

void PERMONMILPSolver::get_dual_solution( Configuration * solc )
{
 if( ! solved() )
  throw( std::logic_error( "PERMONMILPSolver::get_dual_solution: "
			   "no solution available" ) );

 // H x + g + A' y + C' y_C = 0 for the problem PermonQP minimizes, the
 // multiplier y_C of the row C x - s = 0 of a slack being that of the
 // row of the Block (which equals the z of the bounds of the slack): the
 // multipliers of the problem of the Block are - y for a minimization and
 // y for a maximization, whose Objective PermonQP has been given negated
 std::vector< double > pi( d_row_pos.size() );
 for( std::size_t i = 0 ; i < d_row_pos.size() ; ++i ) {
  const int p = d_row_pos[ i ];
  pi[ i ] = -d_sign * d_y[ p >= 0 ? p : d_neq + d_nfix - p - 1 ];
  }

 // the multiplier of the bounds of a fixed variable is that of its row
 std::vector< double > rc( d_n );
 for( int j = 0 ; j < d_n ; ++j )
  rc[ j ] =
    -d_sign * ( d_fix_row[ j ] >= 0 ? d_y[ d_fix_row[ j ] ] : d_z[ j ] );

 MILPSolver::write_dual_solution( pi , rc , solc );
 }

/*--------------------------------------------------------------------------*/
/*------------------ METHODS FOR HANDLING THE PARAMETERS -------------------*/
/*--------------------------------------------------------------------------*/

void PERMONMILPSolver::set_par( idx_type par , int value )
{
 switch( par ) {
  case( intMaxIter ):
   d_max_iter = std::max( value , 1 );
   return;
  case( intPERMONWarmStart ):
   d_warm_start = bool( value );
   return;
  case( intLogVerb ):  // the base class does not keep it
   log_verbosity = value;
   break;
  }

 MILPSolver::set_par( par , value );
 }

/*--------------------------------------------------------------------------*/

void PERMONMILPSolver::set_par( idx_type par , double value )
{
 switch( par ) {
  case( dblFAccSol ):
   d_atol = value;
   return;
  case( dblPERMONRelTol ):
   d_rtol = value;
   return;
  }

 MILPSolver::set_par( par , value );
 }

/*--------------------------------------------------------------------------*/

int PERMONMILPSolver::get_dflt_int_par( idx_type par ) const
{
 if( par == intPERMONWarmStart )
  return( 1 );

 return( MILPSolver::get_dflt_int_par( par ) );
 }

/*--------------------------------------------------------------------------*/

double PERMONMILPSolver::get_dflt_dbl_par( idx_type par ) const
{
 switch( par ) {
  case( dblFAccSol ):
   return( PERMON_dflt_atol );
  case( dblPERMONRelTol ):
   return( PERMON_dflt_rtol );
  case( dblKKTTol ):
   return( PERMON_dflt_kkt_tol );
  }

 return( MILPSolver::get_dflt_dbl_par( par ) );
 }

/*--------------------------------------------------------------------------*/

int PERMONMILPSolver::get_int_par( idx_type par ) const
{
 switch( par ) {
  case( intMaxIter ):
   return( d_max_iter );
  case( intPERMONWarmStart ):
   return( d_warm_start );
  }

 return( MILPSolver::get_int_par( par ) );
 }

/*--------------------------------------------------------------------------*/

double PERMONMILPSolver::get_dbl_par( idx_type par ) const
{
 switch( par ) {
  case( dblFAccSol ):
   return( d_atol );
  case( dblPERMONRelTol ):
   return( d_rtol );
  }

 return( MILPSolver::get_dbl_par( par ) );
 }

/*--------------------------------------------------------------------------*/

Solver::idx_type
PERMONMILPSolver::int_par_str2idx( const std::string & name ) const
{
 auto it =
   std::find( PERMON_int_pars_str.begin() , PERMON_int_pars_str.end() , name );
 if( it != PERMON_int_pars_str.end() )
  return( intLastAlgParMILP + ( it - PERMON_int_pars_str.begin() ) );

 return( MILPSolver::int_par_str2idx( name ) );
 }

/*--------------------------------------------------------------------------*/

const std::string & PERMONMILPSolver::int_par_idx2str( idx_type idx ) const
{
 if( ( idx >= intLastAlgParMILP ) && ( idx < intLastAlgParPERMON ) )
  return( PERMON_int_pars_str[ idx - intLastAlgParMILP ] );

 return( MILPSolver::int_par_idx2str( idx ) );
 }

/*--------------------------------------------------------------------------*/

Solver::idx_type
PERMONMILPSolver::dbl_par_str2idx( const std::string & name ) const
{
 auto it =
   std::find( PERMON_dbl_pars_str.begin() , PERMON_dbl_pars_str.end() , name );
 if( it != PERMON_dbl_pars_str.end() )
  return( dblLastAlgParMILP + ( it - PERMON_dbl_pars_str.begin() ) );

 return( MILPSolver::dbl_par_str2idx( name ) );
 }

/*--------------------------------------------------------------------------*/

const std::string & PERMONMILPSolver::dbl_par_idx2str( idx_type idx ) const
{
 if( ( idx >= dblLastAlgParMILP ) && ( idx < dblLastAlgParPERMON ) )
  return( PERMON_dbl_pars_str[ idx - dblLastAlgParMILP ] );

 return( MILPSolver::dbl_par_idx2str( idx ) );
 }

/*--------------------------------------------------------------------------*/
/*--------------------- End File PERMONMILPSolver.cpp ----------------------*/
/*--------------------------------------------------------------------------*/
