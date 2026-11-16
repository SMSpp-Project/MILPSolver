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

 // the fixed columns and the free rows are not given to qpOASES; the
 // positions in the problem of qpOASES of the columns and rows of the
 // previous problem are kept, to move to the new one its solution
 std::vector< int > ppos = std::move( d_dq.pos );
 std::vector< int > prpos = std::move( d_dq.rpos );
 const int qn = std::count_if( ppos.begin() , ppos.end() ,
			       []( int k ) { return( k >= 0 ); } );
 const int qm = std::count_if( prpos.begin() , prpos.end() ,
			       []( int k ) { return( k >= 0 ); } );
 d_sign = ( objsense == -1 ) ? -1 : 1;  // qpOASES minimizes
 make_dense_qp( d_dq , d_sign , qpOASES::INFTY );
 const int n = d_n = d_dq.n;
 d_nc = d_dq.m;

 // a solved problem with other columns or rows cannot be hot-started, but
 // its working set, solution and multipliers, moved to the columns and
 // rows that are still there [see previous_column_index()], are the guess
 // with which the next one is initialized, the new ones inactive and at 0;
 // an active side that is no longer finite is inactive
 d_guess = false;
 const auto [ pn , pm ] = previous_sizes();
 d_same = ( pn == numcols ) && ( pm == numrows ) &&
          ( ppos == d_dq.pos ) && ( prpos == d_dq.rpos );
 for( int j = 0 ; d_same && ( j < numcols ) ; ++j )
  d_same = ( previous_column_index( j ) == j );
 for( int i = 0 ; d_same && ( i < numrows ) ; ++i )
  d_same = ( previous_row_index( i ) == i );
 if( d_hot_start && d_qp && d_qp_solved && ( ! d_same ) &&
     ( d_qp->getNV() == qn ) && ( d_qp->getNC() == qm ) &&
     ( ppos.size() == std::size_t( pn ) ) &&
     ( prpos.size() == std::size_t( pm ) ) ) {
  qpOASES::Bounds ob;
  qpOASES::Constraints oc;
  d_qp->getBounds( ob );
  d_qp->getConstraints( oc );
  std::vector< double > ox( qn ) , oy( qn + qm );
  d_qp->getPrimalSolution( ox.data() );
  d_qp->getDualSolution( oy.data() );

  auto status = []( qpOASES::SubjectToStatus st , double l , double u ) {
   if( ( st == qpOASES::ST_LOWER ) && ( l > - qpOASES::INFTY ) )
    return( qpOASES::ST_LOWER );
   if( ( st == qpOASES::ST_UPPER ) && ( u < qpOASES::INFTY ) )
    return( qpOASES::ST_UPPER );
   return( qpOASES::ST_INACTIVE );
   };

  d_guess_x.assign( n , 0 );
  d_guess_y.assign( n + d_nc , 0 );
  d_guess_b = std::make_unique< qpOASES::Bounds >();
  d_guess_c = std::make_unique< qpOASES::Constraints >();
  d_guess_b->init( n );
  d_guess_c->init( d_nc );
  for( int j = 0 ; j < n ; ++j ) {
   const int pc = previous_column_index( d_dq.col[ j ] );
   const int pj = ( pc >= 0 ) ? ppos[ pc ] : -1;
   d_guess_b->setupBound( j , pj < 0 ? qpOASES::ST_INACTIVE :
                          status( ob.getStatus( pj ) , d_dq.lb[ j ] ,
                                  d_dq.ub[ j ] ) );
   if( pj >= 0 ) {
    d_guess_x[ j ] = ox[ pj ];
    d_guess_y[ j ] = oy[ pj ];
    }
   }
  for( int i = 0 ; i < d_nc ; ++i ) {
   const int pr = previous_row_index( d_dq.row[ i ] );
   const int pi = ( pr >= 0 ) ? prpos[ pr ] : -1;
   d_guess_c->setupConstraint( i , pi < 0 ? qpOASES::ST_INACTIVE :
			       status( oc.getStatus( pi ) , d_dq.rl[ i ] ,
				       d_dq.ru[ i ] ) );
   if( pi >= 0 )
    d_guess_y[ n + i ] = oy[ qn + pi ];
   }
  d_guess = true;
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
 d_fval = d_dq.fconst;

 // qpOASES needs at least one variable: with none, the problem is feasible
 // if and only if each row admits 0
 if( d_n == 0 ) {
  const double eps = d_options->terminationTolerance;
  d_status = kOK;
  for( int i = 0 ; i < d_nc ; ++i )
   if( ( d_dq.rl[ i ] > eps ) || ( d_dq.ru[ i ] < - eps ) )
    d_status = kInfeasible;
  d_y.assign( d_nc , 0 );
  return( sol_status = d_status );
  }

 // a hot start needs the problem of qpOASES solved, with the same columns
 // and rows in the same places and with the same Options, otherwise the
 // new problem is initialized from the guess made out of the previous one
 // by load_problem(), if any; one that does not end with a solution, unless
 // it is stopped by the time limit, is repeated from scratch
 if( d_hot_start && d_qp && d_qp_solved && ( ! d_new_options ) &&
     d_same && ( d_qp->getNV() == d_n ) && ( d_qp->getNC() == d_nc ) ) {
  d_status = hot_start();
  if( ( d_status != kOK ) && ( d_status != kStopTime ) )
   d_status = solve_from_scratch();
  }
 else {
  d_status = solve_from_scratch( d_guess && ( ! d_new_options ) );
  if( d_guess && ( d_status != kOK ) && ( d_status != kStopTime ) )
   d_status = solve_from_scratch( false );
  }
 d_guess = false;

 // with a singular H and no regularization qpOASES may cycle until the
 // limit on the working set recalculations, which is only a safeguard when
 // intMaxIter is not set (as on the dual of an SVM with more samples than
 // features): the solve is then repeated from scratch with the H
 // regularized, which the problem of qpOASES keeps for its hot starts
 if( ( d_status == kStopIter ) && ( d_max_iter == Inf< int >() ) &&
     ( d_hessian_type != 1 ) && ( ! d_options->enableRegularisation ) ) {
  d_options->enableRegularisation = qpOASES::BT_TRUE;
  d_status = solve_from_scratch();
  d_options->enableRegularisation = qpOASES::BT_FALSE;
  }

 d_qp_solved = ( d_status == kOK );

 if( d_qp_solved ) {
  d_qp->getPrimalSolution( d_x.data() );
  d_qp->getDualSolution( d_y.data() );

  // the value is computed here rather than taken from qpOASES, whose value
  // is that of the regularized problem if H has been regularized
  for( int i = 0 ; i < d_n ; ++i ) {
   double Hx = 0;
   for( int j = 0 ; j < d_n ; ++j )
    Hx += d_dq.H[ std::size_t( i ) * d_n + j ] * d_x[ j ];
   d_fval += d_x[ i ] * ( d_dq.f[ i ] + Hx / 2 );
   }
  }

 return( sol_status = d_status );

 }  // end( QPOASESMILPSolver::guts_of_compute )

/*--------------------------------------------------------------------------*/

int QPOASESMILPSolver::solve_from_scratch( bool guess )
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
 d_qpH[ 0 ] = d_oldH = d_dq.H;
 d_qpA[ 0 ] = d_oldA = d_dq.A;

 // a guess that is good takes few working set recalculations, and on a
 // singular H a bad one may take as many as from scratch, or more: a
 // guessed start is given up after half as many as variables and rows,
 // and the solve repeated from scratch
 const int nwsr = guess ? std::min( max_nwsr() ,
				    std::max( 10 , ( d_n + d_nc ) / 2 ) )
                        : max_nwsr();
 qpOASES::int_t done = nwsr;
 qpOASES::real_t time = d_max_time;
 const bool timed = d_max_time < Inf< double >();

 const auto ret = guess ?
  d_qp->init( d_qpH[ 0 ].data() , d_dq.f.data() ,
	      d_nc ? d_qpA[ 0 ].data() : nullptr ,
	      d_dq.lb.data() , d_dq.ub.data() ,
	      d_nc ? d_dq.rl.data() : nullptr ,
	      d_nc ? d_dq.ru.data() : nullptr ,
	      done , timed ? & time : nullptr ,
	      d_guess_x.data() , d_guess_y.data() ,
	      d_guess_b.get() , d_guess_c.get() ) :
  d_qp->init( d_qpH[ 0 ].data() , d_dq.f.data() ,
	      d_nc ? d_qpA[ 0 ].data() : nullptr ,
	      d_dq.lb.data() , d_dq.ub.data() ,
	      d_nc ? d_dq.rl.data() : nullptr ,
	      d_nc ? d_dq.ru.data() : nullptr ,
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

 if( ( d_dq.H == d_oldH ) && ( d_dq.A == d_oldA ) )
  // only the vectors have changed, if anything
  ret = d_qp->hotstart( d_dq.f.data() , d_dq.lb.data() , d_dq.ub.data() ,
			d_nc ? d_dq.rl.data() : nullptr ,
			d_nc ? d_dq.ru.data() : nullptr ,
			done , timed ? & time : nullptr );
 else {
  // the new matrices go into the other copy, the current one being still
  // read by qpOASES while it moves to the new ones
  const int nxt = 1 - d_cur;
  d_qpH[ nxt ] = d_oldH = d_dq.H;
  d_qpA[ nxt ] = d_oldA = d_dq.A;
  ret = d_qp->hotstart( d_qpH[ nxt ].data() , d_dq.f.data() ,
			d_nc ? d_qpA[ nxt ].data() : nullptr ,
			d_dq.lb.data() , d_dq.ub.data() ,
			d_nc ? d_dq.rl.data() : nullptr ,
			d_nc ? d_dq.ru.data() : nullptr ,
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

 // with H positive definite the problem cannot be unbounded: an unbounded
 // qpOASES only says that a projected H has turned out numerically singular
 const int unb = ( d_hessian_type == 1 ) ? kError : kUnbounded;

 if( ( ret == qpOASES::RET_QP_UNBOUNDED ) || d_qp->isUnbounded() )
  return( unb );

 switch( qpOASES::getSimpleStatus( qpOASES::returnValue( ret ) ) ) {
  case( 1 ):   // the working set recalculations, or the time, are over
   return( ( ( d_max_time < Inf< double >() ) && ( done < nwsr ) ) ?
	   kStopTime : kStopIter );
  case( -2 ):  return( kInfeasible );
  case( -3 ):  return( unb );
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

 MILPSolver::write_var_solution( dense_qp_x( d_dq , d_x.data() ) ,
				  solc );
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
 // - y for a maximization, whose Objective qpOASES has been given negated,
 // and 0 on a row that qpOASES has not been given
 std::vector< double > pi( numrows , 0 );
 for( int i = 0 ; i < d_nc ; ++i )
  pi[ d_dq.row[ i ] ] = d_sign * d_y[ d_n + i ];

 std::vector< double > rc( numcols , 0 );
 for( int j = 0 ; j < d_n ; ++j )
  rc[ d_dq.col[ j ] ] = d_sign * d_y[ j ];
 dense_qp_fixed_rc( d_dq , dense_qp_x( d_dq , d_x.data() ) , pi , rc );

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
