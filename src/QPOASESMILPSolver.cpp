/*--------------------------------------------------------------------------*/
/*----------------------- File QPOASESMILPSolver.cpp -----------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Implementation of the QPOASESMILPSolver class, which implements the
 * MILPSolver interface for convex Quadratic Programs with linear
 * constraints using the online active-set solver qpOASES.
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
#include <stdexcept>

#include "QPOASESMILPSolver.h"

#include <qpOASES.hpp>

/*--------------------------------------------------------------------------*/
/*-------------------------------- USING -----------------------------------*/
/*--------------------------------------------------------------------------*/

using namespace SMSpp_di_unipi_it;

/*--------------------------------------------------------------------------*/
/*-------------------------------- CONSTANTS -------------------------------*/
/*--------------------------------------------------------------------------*/

static const std::array< std::string , 8 > QPOASES_int_pars_str = {
 "intQPOASESHotStart" , "intQPOASESPreset" , "intQPOASESHessianType" ,
 "intQPOASESEnableReg" , "intQPOASESRegSteps" , "intQPOASESFarBounds" ,
 "intQPOASESFlipBounds" , "intQPOASESRefineSteps" };

static const std::array< std::string , 5 > QPOASES_dbl_pars_str = {
 "dblQPOASESBoundTol" , "dblQPOASESEpsReg" , "dblQPOASESInitFarBounds" ,
 "dblQPOASESGrowFarBounds" , "dblQPOASESBoundRelax" };

/*--------------------------------------------------------------------------*/
/*-------------------------------- FUNCTIONS -------------------------------*/
/*--------------------------------------------------------------------------*/

// a side of SMS++ in the format of qpOASES, where an infinite one is INFTY

static double qpoases_side( double v )
{
 if( v >= qpOASES::INFTY )
  return( qpOASES::INFTY );
 if( v <= - qpOASES::INFTY )
  return( - qpOASES::INFTY );
 return( v );
 }

/*--------------------------------------------------------------------------*/

// the default Options of qpOASES

static const qpOASES::Options & default_options( void )
{
 static const qpOASES::Options options;
 return( options );
 }

/*--------------------------------------------------------------------------*/
/*------------------------- FACTORY MANAGEMENT -----------------------------*/
/*--------------------------------------------------------------------------*/

SMSpp_insert_in_factory_cpp_0( QPOASESMILPSolver );

/*--------------------------------------------------------------------------*/
/*--------------------- CONSTRUCTOR AND DESTRUCTOR -------------------------*/
/*--------------------------------------------------------------------------*/

QPOASESMILPSolver::QPOASESMILPSolver( void ) : MILPSolver() ,
 d_options( std::make_unique< qpOASES::Options >() )
{
 set_print_level();
 }

/*--------------------------------------------------------------------------*/

QPOASESMILPSolver::~QPOASESMILPSolver() = default;

/*--------------------------------------------------------------------------*/
/*--------------------- DERIVED METHODS OF BASE CLASS ----------------------*/
/*--------------------------------------------------------------------------*/

void QPOASESMILPSolver::clear_problem( unsigned int what )
{
 MILPSolver::clear_problem( what );
 d_qp.reset();
 d_qp_solved = false;
 d_n = d_nc = 0;
 d_status = kError;
 }

/*--------------------------------------------------------------------------*/

void QPOASESMILPSolver::load_problem( void )
{
 // the problem of qpOASES, which reads the copies of H and A of the
 // previous problem, is kept across the clear_problem() of MILPSolver, so
 // that the problem loaded again whole after a Modification is a hot start
 // from its solution
 auto qp = std::move( d_qp );
 const bool qp_solved = d_qp_solved;

 MILPSolver::load_problem();

 d_qp = std::move( qp );
 d_qp_solved = qp_solved;

 if( numquadrows > 0 )
  throw( std::invalid_argument( "QPOASESMILPSolver::load_problem: "
				"quadratic constraints are not supported" ) );

 if( ( int_vars > 0 ) && ( relax_int_vars == 0 ) )
  throw( std::invalid_argument( "QPOASESMILPSolver::load_problem: integer "
				"variables are only supported if relaxed "
				"(intRelaxIntVars)" ) );

 const int n = numcols;
 d_n = n;
 d_nc = numrows;
 d_sign = ( objsense == -1 ) ? -1 : 1;  // qpOASES minimizes

 // the objective: c' x + sum_j q_j x_j^2 + sum_k ndq_k x_i x_j is
 // 1/2 x' H x + g' x with H_jj = 2 q_j and H_ij = H_ji = ndq_k
 d_g.assign( n , 0 );
 for( int j = 0 ; j < n ; ++j )
  d_g[ j ] = d_sign * objective[ j ];

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
 d_lb.resize( n );
 d_ub.resize( n );
 for( int j = 0 ; j < n ; ++j ) {
  d_lb[ j ] = qpoases_side( lb[ j ] );
  d_ub[ j ] = qpoases_side( ub[ j ] );
  }

 d_lbA.resize( d_nc );
 d_ubA.resize( d_nc );
 for( int i = 0 ; i < d_nc ; ++i )
  switch( sense[ i ] ) {
   case( 'L' ): d_lbA[ i ] = - qpOASES::INFTY;
                d_ubA[ i ] = qpoases_side( rhs[ i ] );
                break;
   case( 'G' ): d_lbA[ i ] = qpoases_side( rhs[ i ] );
                d_ubA[ i ] = qpOASES::INFTY;
                break;
   case( 'E' ): d_lbA[ i ] = d_ubA[ i ] = rhs[ i ];
                break;
   case( 'R' ): if( rngval[ i ] > 0 ) {
                 d_lbA[ i ] = qpoases_side( rhs[ i ] );
                 d_ubA[ i ] = qpoases_side( rhs[ i ] + rngval[ i ] );
                 }
                else {
                 d_lbA[ i ] = qpoases_side( rhs[ i ] + rngval[ i ] );
                 d_ubA[ i ] = qpoases_side( rhs[ i ] );
                 }
                break;
   default:
    throw( std::invalid_argument( "QPOASESMILPSolver::load_problem: row " +
				  std::to_string( i ) + " has sense " +
				  std::string( 1 , sense[ i ] ) ) );
   }

 d_x.assign( n , 0 );
 d_y.assign( n + d_nc , 0 );
 d_status = kError;

 }  // end( QPOASESMILPSolver::load_problem )

/*--------------------------------------------------------------------------*/
/*-------------------- METHODS FOR SOLVING THE MODEL -----------------------*/
/*--------------------------------------------------------------------------*/

int QPOASESMILPSolver::guts_of_compute( void )
{
 d_fval = 0;

 // qpOASES needs at least one variable: with none, the problem is feasible
 // if and only if each row admits 0
 if( d_n == 0 ) {
  const double eps = d_options->terminationTolerance;
  d_status = kOK;
  for( int i = 0 ; i < d_nc ; ++i )
   if( ( d_lbA[ i ] > eps ) || ( d_ubA[ i ] < - eps ) )
    d_status = kInfeasible;
  d_y.assign( d_nc , 0 );
  return( sol_status = d_status );
  }

 // a hot start needs the problem of qpOASES solved, of the same size and
 // with the same Options; one that does not end with a solution, unless it
 // is stopped by the time limit, is repeated from scratch
 if( d_hot_start && d_qp && d_qp_solved && ( ! d_new_options ) &&
     ( d_qp->getNV() == d_n ) && ( d_qp->getNC() == d_nc ) ) {
  d_status = hot_start();
  if( ( d_status != kOK ) && ( d_status != kStopTime ) )
   d_status = solve_from_scratch();
  }
 else
  d_status = solve_from_scratch();

 d_qp_solved = ( d_status == kOK );

 if( d_qp_solved ) {
  d_qp->getPrimalSolution( d_x.data() );
  d_qp->getDualSolution( d_y.data() );

  // the value is computed here rather than taken from qpOASES, whose value
  // is that of the regularized problem if H has been regularized
  for( int i = 0 ; i < d_n ; ++i ) {
   double Hx = 0;
   for( int j = 0 ; j < d_n ; ++j )
    Hx += d_H[ std::size_t( i ) * d_n + j ] * d_x[ j ];
   d_fval += d_x[ i ] * ( d_g[ i ] + Hx / 2 );
   }
  }

 return( sol_status = d_status );

 }  // end( QPOASESMILPSolver::guts_of_compute )

/*--------------------------------------------------------------------------*/

int QPOASESMILPSolver::solve_from_scratch( void )
{
 qpOASES::HessianType ht = qpOASES::HST_UNKNOWN;
 if( d_hessian_type == 1 )
  ht = qpOASES::HST_POSDEF;
 else
  if( d_hessian_type == 2 )
   ht = qpOASES::HST_SEMIDEF;

 d_qp = std::make_unique< qpOASES::SQProblem >( d_n , d_nc , ht );
 d_qp->setOptions( *d_options );
 d_new_options = false;

 d_cur = 0;
 d_qpH[ 0 ] = d_oldH = d_H;
 d_qpA[ 0 ] = d_oldA = d_A;

 const int nwsr = max_nwsr();
 qpOASES::int_t done = nwsr;
 qpOASES::real_t time = d_max_time;
 const bool timed = d_max_time < Inf< double >();

 const auto ret = d_qp->init( d_qpH[ 0 ].data() , d_g.data() ,
			      d_nc ? d_qpA[ 0 ].data() : nullptr ,
			      d_lb.data() , d_ub.data() ,
			      d_nc ? d_lbA.data() : nullptr ,
			      d_nc ? d_ubA.data() : nullptr ,
			      done , timed ? & time : nullptr );
 elapsed_iter += done;

 return( decode( ret , nwsr , done ) );
 }

/*--------------------------------------------------------------------------*/

int QPOASESMILPSolver::hot_start( void )
{
 d_qp->setPrintLevel( d_options->printLevel );

 const int nwsr = max_nwsr();
 qpOASES::int_t done = nwsr;
 qpOASES::real_t time = d_max_time;
 const bool timed = d_max_time < Inf< double >();
 qpOASES::returnValue ret;

 if( ( d_H == d_oldH ) && ( d_A == d_oldA ) )
  // only the vectors have changed, if anything
  ret = d_qp->hotstart( d_g.data() , d_lb.data() , d_ub.data() ,
			d_nc ? d_lbA.data() : nullptr ,
			d_nc ? d_ubA.data() : nullptr ,
			done , timed ? & time : nullptr );
 else {
  // the new matrices go into the other copy, the current one being still
  // read by qpOASES while it moves to the new ones
  const int nxt = 1 - d_cur;
  d_qpH[ nxt ] = d_oldH = d_H;
  d_qpA[ nxt ] = d_oldA = d_A;
  ret = d_qp->hotstart( d_qpH[ nxt ].data() , d_g.data() ,
			d_nc ? d_qpA[ nxt ].data() : nullptr ,
			d_lb.data() , d_ub.data() ,
			d_nc ? d_lbA.data() : nullptr ,
			d_nc ? d_ubA.data() : nullptr ,
			done , timed ? & time : nullptr );
  d_cur = nxt;
  }
 elapsed_iter += done;

 return( decode( ret , nwsr , done ) );
 }

/*--------------------------------------------------------------------------*/

int QPOASESMILPSolver::decode( int ret , int nwsr , int done ) const
{
 if( ret == qpOASES::SUCCESSFUL_RETURN )
  return( kOK );

 if( ( ret == qpOASES::RET_QP_INFEASIBLE ) || d_qp->isInfeasible() )
  return( kInfeasible );

 if( ( ret == qpOASES::RET_QP_UNBOUNDED ) || d_qp->isUnbounded() )
  return( kUnbounded );

 switch( qpOASES::getSimpleStatus( qpOASES::returnValue( ret ) ) ) {
  case( 1 ):   // the working set recalculations, or the time, are over
   return( ( ( d_max_time < Inf< double >() ) && ( done < nwsr ) ) ?
	   kStopTime : kStopIter );
  case( -2 ):  return( kInfeasible );
  case( -3 ):  return( kUnbounded );
  default:     return( kError );
  }
 }

/*--------------------------------------------------------------------------*/

int QPOASESMILPSolver::max_nwsr( void ) const
{
 if( d_max_iter < Inf< int >() )
  return( d_max_iter );

 return( 100 * ( d_n + d_nc ) + 1000 );
 }

/*--------------------------------------------------------------------------*/

void QPOASESMILPSolver::set_print_level( void )
{
 if( log_verbosity <= 0 )
  d_options->printLevel = qpOASES::PL_NONE;
 else
  if( log_verbosity == 1 )
   d_options->printLevel = qpOASES::PL_LOW;
  else
   d_options->printLevel = ( log_verbosity == 2 ) ? qpOASES::PL_MEDIUM
                                                  : qpOASES::PL_HIGH;
 }

/*--------------------------------------------------------------------------*/

void QPOASESMILPSolver::apply_preset( void )
{
 switch( d_preset ) {
  case( 1 ): d_options->setToReliable(); break;
  case( 2 ): d_options->setToMPC();      break;
  default:   d_options->setToDefault();
  }

 for( const auto & [ par , value ] : d_int_set )
  set_option( par , value );
 for( const auto & [ par , value ] : d_dbl_set )
  set_option( par , value );

 set_print_level();
 d_new_options = true;
 }

/*--------------------------------------------------------------------------*/

bool QPOASESMILPSolver::set_option( idx_type par , int value )
{
 const auto bt = value ? qpOASES::BT_TRUE : qpOASES::BT_FALSE;

 switch( par ) {
  case( intQPOASESEnableReg ):  d_options->enableRegularisation = bt;
                                return( true );
  case( intQPOASESRegSteps ):   d_options->numRegularisationSteps = value;
                                return( true );
  case( intQPOASESFarBounds ):  d_options->enableFarBounds = bt;
                                return( true );
  case( intQPOASESFlipBounds ): d_options->enableFlippingBounds = bt;
                                return( true );
  case( intQPOASESRefineSteps ): d_options->numRefinementSteps = value;
                                 return( true );
  }

 return( false );
 }

/*--------------------------------------------------------------------------*/

bool QPOASESMILPSolver::set_option( idx_type par , double value )
{
 switch( par ) {
  case( dblFAccSol ):         d_options->terminationTolerance = value;
                              return( true );
  case( dblQPOASESBoundTol ): d_options->boundTolerance = value;
                              return( true );
  case( dblQPOASESEpsReg ):   d_options->epsRegularisation = value;
                              return( true );
  case( dblQPOASESInitFarBounds ): d_options->initialFarBounds = value;
                                   return( true );
  case( dblQPOASESGrowFarBounds ): d_options->growFarBounds = value;
                                   return( true );
  case( dblQPOASESBoundRelax ): d_options->boundRelaxation = value;
                                return( true );
  }

 return( false );
 }

/*--------------------------------------------------------------------------*/
/*---------------------- METHODS FOR READING RESULTS -----------------------*/
/*--------------------------------------------------------------------------*/

Solver::OFValue QPOASESMILPSolver::get_lb( void )
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

Solver::OFValue QPOASESMILPSolver::get_ub( void )
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

bool QPOASESMILPSolver::has_var_solution( void )
{
 return( solved() );
 }

/*--------------------------------------------------------------------------*/

Solver::OFValue QPOASESMILPSolver::get_var_value( void )
{
 return( objsense == -1 ? get_lb() : get_ub() );
 }

/*--------------------------------------------------------------------------*/

void QPOASESMILPSolver::get_var_solution( Configuration * solc )
{
 if( ! solved() )
  throw( std::logic_error( "QPOASESMILPSolver::get_var_solution: "
			   "no solution available" ) );

 MILPSolver::write_var_solution( d_x , solc );
 }

/*--------------------------------------------------------------------------*/

bool QPOASESMILPSolver::has_dual_solution( void )
{
 return( solved() );
 }

/*--------------------------------------------------------------------------*/

bool QPOASESMILPSolver::is_dual_feasible( void )
{
 return( solved() );
 }

/*--------------------------------------------------------------------------*/

void QPOASESMILPSolver::get_dual_solution( Configuration * solc )
{
 if( ! solved() )
  throw( std::logic_error( "QPOASESMILPSolver::get_dual_solution: "
			   "no solution available" ) );

 // qpOASES has H x + g = y_x + A' y_A, with y > 0 on an active lower side
 // and y < 0 on an active upper one, for the problem it minimizes: the
 // multipliers of the problem of the Block are y for a minimization and
 // - y for a maximization, whose Objective qpOASES has been given negated
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

void QPOASESMILPSolver::set_par( idx_type par , int value )
{
 switch( par ) {
  case( intMaxIter ): d_max_iter = std::max( value , 0 ); return;
  case( intLogVerb ):
   log_verbosity = value;
   set_print_level();
   break;  // the base class may want it too
  case( intQPOASESHotStart ): d_hot_start = bool( value ); return;
  case( intQPOASESPreset ):
   d_preset = value;
   apply_preset();
   return;
  case( intQPOASESHessianType ):
   d_hessian_type = value;
   d_new_options = true;
   return;
  default:
   if( set_option( par , value ) ) {
    d_int_set[ par ] = value;
    d_new_options = true;
    return;
    }
  }

 MILPSolver::set_par( par , value );
 }

/*--------------------------------------------------------------------------*/

void QPOASESMILPSolver::set_par( idx_type par , double value )
{
 if( par == dblMaxTime ) {
  d_max_time = std::max( value , 0.0 );
  return;
  }

 if( set_option( par , value ) ) {
  d_dbl_set[ par ] = value;
  d_new_options = true;
  return;
  }

 MILPSolver::set_par( par , value );
 }

/*--------------------------------------------------------------------------*/

int QPOASESMILPSolver::get_dflt_int_par( idx_type par ) const
{
 const auto & o = default_options();

 switch( par ) {
  case( intQPOASESHotStart ):    return( 1 );
  case( intQPOASESPreset ):      return( 0 );
  case( intQPOASESHessianType ): return( 0 );
  case( intQPOASESEnableReg ):   return( o.enableRegularisation );
  case( intQPOASESRegSteps ):    return( o.numRegularisationSteps );
  case( intQPOASESFarBounds ):   return( o.enableFarBounds );
  case( intQPOASESFlipBounds ):  return( o.enableFlippingBounds );
  case( intQPOASESRefineSteps ): return( o.numRefinementSteps );
  }

 return( MILPSolver::get_dflt_int_par( par ) );
 }

/*--------------------------------------------------------------------------*/

double QPOASESMILPSolver::get_dflt_dbl_par( idx_type par ) const
{
 const auto & o = default_options();

 switch( par ) {
  case( dblFAccSol ):              return( o.terminationTolerance );
  case( dblQPOASESBoundTol ):      return( o.boundTolerance );
  case( dblQPOASESEpsReg ):        return( o.epsRegularisation );
  case( dblQPOASESInitFarBounds ): return( o.initialFarBounds );
  case( dblQPOASESGrowFarBounds ): return( o.growFarBounds );
  case( dblQPOASESBoundRelax ):    return( o.boundRelaxation );
  }

 return( MILPSolver::get_dflt_dbl_par( par ) );
 }

/*--------------------------------------------------------------------------*/

int QPOASESMILPSolver::get_int_par( idx_type par ) const
{
 const auto & o = *d_options;

 switch( par ) {
  case( intMaxIter ):            return( d_max_iter );
  case( intLogVerb ):            return( log_verbosity );
  case( intQPOASESHotStart ):    return( d_hot_start );
  case( intQPOASESPreset ):      return( d_preset );
  case( intQPOASESHessianType ): return( d_hessian_type );
  case( intQPOASESEnableReg ):   return( o.enableRegularisation );
  case( intQPOASESRegSteps ):    return( o.numRegularisationSteps );
  case( intQPOASESFarBounds ):   return( o.enableFarBounds );
  case( intQPOASESFlipBounds ):  return( o.enableFlippingBounds );
  case( intQPOASESRefineSteps ): return( o.numRefinementSteps );
  }

 return( MILPSolver::get_int_par( par ) );
 }

/*--------------------------------------------------------------------------*/

double QPOASESMILPSolver::get_dbl_par( idx_type par ) const
{
 const auto & o = *d_options;

 switch( par ) {
  case( dblMaxTime ):              return( d_max_time );
  case( dblFAccSol ):              return( o.terminationTolerance );
  case( dblQPOASESBoundTol ):      return( o.boundTolerance );
  case( dblQPOASESEpsReg ):        return( o.epsRegularisation );
  case( dblQPOASESInitFarBounds ): return( o.initialFarBounds );
  case( dblQPOASESGrowFarBounds ): return( o.growFarBounds );
  case( dblQPOASESBoundRelax ):    return( o.boundRelaxation );
  }

 return( MILPSolver::get_dbl_par( par ) );
 }

/*--------------------------------------------------------------------------*/

Solver::idx_type QPOASESMILPSolver::int_par_str2idx(
					      const std::string & name ) const
{
 auto it = std::find( QPOASES_int_pars_str.begin() ,
		      QPOASES_int_pars_str.end() , name );
 if( it != QPOASES_int_pars_str.end() )
  return( intLastAlgParMILP + ( it - QPOASES_int_pars_str.begin() ) );

 return( MILPSolver::int_par_str2idx( name ) );
 }

/*--------------------------------------------------------------------------*/

const std::string & QPOASESMILPSolver::int_par_idx2str( idx_type idx ) const
{
 if( ( idx >= intLastAlgParMILP ) && ( idx < intLastAlgParQPOASES ) )
  return( QPOASES_int_pars_str[ idx - intLastAlgParMILP ] );

 return( MILPSolver::int_par_idx2str( idx ) );
 }

/*--------------------------------------------------------------------------*/

Solver::idx_type QPOASESMILPSolver::dbl_par_str2idx(
					      const std::string & name ) const
{
 auto it = std::find( QPOASES_dbl_pars_str.begin() ,
		      QPOASES_dbl_pars_str.end() , name );
 if( it != QPOASES_dbl_pars_str.end() )
  return( dblLastAlgParMILP + ( it - QPOASES_dbl_pars_str.begin() ) );

 return( MILPSolver::dbl_par_str2idx( name ) );
 }

/*--------------------------------------------------------------------------*/

const std::string & QPOASESMILPSolver::dbl_par_idx2str( idx_type idx ) const
{
 if( ( idx >= dblLastAlgParMILP ) && ( idx < dblLastAlgParQPOASES ) )
  return( QPOASES_dbl_pars_str[ idx - dblLastAlgParMILP ] );

 return( MILPSolver::dbl_par_idx2str( idx ) );
 }

/*--------------------------------------------------------------------------*/
/*-------------------- End File QPOASESMILPSolver.cpp ----------------------*/
/*--------------------------------------------------------------------------*/
