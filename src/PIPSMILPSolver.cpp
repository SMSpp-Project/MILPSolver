/*--------------------------------------------------------------------------*/
/*------------------------ File PIPSMILPSolver.cpp -------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Implementation of the PIPSMILPSolver class.
 *
 * \author Enrico Calandrini \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \author Antonio Frangioni \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \copyright &copy; by Enrico Calandrini, Antonio Frangioni
 */
/*--------------------------------------------------------------------------*/
/*---------------------------- IMPLEMENTATION ------------------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include <cstdlib>

#include <queue>

#include <LinearFunction.h>

#include <QuadFunction.h>

#include <mpi_debug.h>

#include "PIPSMILPSolver.h"

#include "PIPS_maps.h"

#ifdef MILPSOLVER_DEBUG
 #define DEBUG_LOG( stuff ) std::cout << "[PIPSMILPSolver DEBUG] " << stuff
#else
 #define DEBUG_LOG( stuff )
#endif

/*--------------------------------------------------------------------------*/
/*------------------------- NAMESPACE AND USING ----------------------------*/
/*--------------------------------------------------------------------------*/

using namespace SMSpp_di_unipi_it;

using pipsipmpp::DistributedInputTree;
using pipsipmpp::PIPSIPMppInterface;
using pipsipmpp::TerminationStatus;
using pipsipmpp::MPIErrorHandler;

/*--------------------------------------------------------------------------*/
/*-------------------------- FACTORY MANAGEMENT ----------------------------*/
/*--------------------------------------------------------------------------*/

SMSpp_insert_in_factory_cpp_0( PIPSMILPSolver );

/*--------------------------------------------------------------------------*/
/*--------------------- CONSTRUCTOR AND DESTRUCTOR -------------------------*/
/*--------------------------------------------------------------------------*/

PIPSMILPSolver::PIPSMILPSolver( void ) :
 MILPSolver() , pips_tree( nullptr ) , pips_interface( nullptr ) { }

/*--------------------------------------------------------------------------*/

PIPSMILPSolver::~PIPSMILPSolver()
{
 reset_pips_data();
}

/*--------------------------------------------------------------------------*/

void PIPSMILPSolver::reset_pips_data()
{
 delete pips_interface;
 pips_interface = nullptr;

 delete pips_tree;
 pips_tree = nullptr;

 n_nodes = 0;
 nodes_subtrees.clear();
 block_to_leaf.clear();

 n_var_node.clear();
 var_node.clear();
 n_eq_cons_node.clear();
 eq_cons_node.clear();
 n_ineq_cons_node.clear();
 ineq_cons_node.clear();

 n_link_eq_cons = 0;
 link_eq_cons.clear();
 n_link_ineq_cons = 0;
 link_ineq_cons.clear();

 var_to_node.clear();
}

/*--------------------------------------------------------------------------*/

void PIPSMILPSolver::clear_problem( unsigned int what )
{
 MILPSolver::clear_problem( what );
 reset_pips_data();
}

/*--------------------------------------------------------------------------*/
/*--------------------- DERIVED METHODS OF BASE CLASS ----------------------*/
/*--------------------------------------------------------------------------*/

void PIPSMILPSolver::set_Block( Block * block )
{
 if( block == f_Block )
  return;

 // MPI is initialised on demand: if the application has not done it already,
 // it is done here (funneled threading level) and finalised at program exit;
 // it is never finalised while any PIPSMILPSolver may still be constructed
 int finalized = 0;
 MPI_Finalized( & finalized );
 if( finalized )
  throw( std::runtime_error(
		   "PIPSMILPSolver::set_Block: MPI has been finalized" ) );

 int initialized = 0;
 MPI_Initialized( & initialized );
 if( ! initialized ) {
  int provided = 0;
  if( MPI_Init_thread( nullptr , nullptr , MPI_THREAD_FUNNELED ,
		       & provided ) != MPI_SUCCESS )
   throw( std::runtime_error(
		    "PIPSMILPSolver::set_Block: MPI_Init_thread failed" ) );

  MPIErrorHandler::init();  // install the MPI error handler

  std::atexit( []( void ) {
		int f = 0;
		MPI_Finalized( & f );
		if( ! f )
		 MPI_Finalize();
		} );
  }

 MILPSolver::set_Block( block );
}

/*--------------------------------------------------------------------------*/

void PIPSMILPSolver::load_problem( void )
{
 // clear the PIPS structures derived from any previously loaded Block; the
 // base class representation is kept alive by the PIPS callbacks, so it is
 // rebuilt from scratch as well
 reset_pips_data();

 // call MILPSolver to generate the entire LP problem
 MILPSolver::load_problem();

 if( numquadrows > 0 )
  throw( std::runtime_error( "PIPSMILPSolver::load_problem: PIPS-IPM++ "
			     "cannot solve quadratically constrained models"
			     ) );

 if( numnnzq > 0 )
  throw( std::runtime_error( "PIPSMILPSolver::load_problem: PIPS-IPM++ "
			     "cannot solve models with a quadratic objective"
			     ) );

 if( int_vars > 0 )
  throw( std::runtime_error( "PIPSMILPSolver::load_problem: only LP models "
	  "can be solved, but integer variables are present; set "
	  "intRelaxIntVars to solve the continuous relaxation" ) );

 /* Each node of the PIPS tree stores the elements (Variable and Constraint)
  * belonging to it; the callbacks then extract the portion of the base class
  * matrix corresponding to the rows and columns of the node and pass it to
  * PIPS. The Block is mapped to the PIPS tree as follows:
  * - f_Block is the root of the PIPS model;
  * - there is one leaf for each nested Block of f_Block, holding the whole
  *   subtree of Block rooted there. */

 // set the first node (root) to be the f_Block
 n_nodes = 1;
 Index n_blocks = 1;
 nodes_subtrees.push_back( { f_Block } );

 // collect all the Blocks representing the leaves of the PIPS tree
 for( auto leaf : f_Block->get_nested_Blocks() ) {
  std::vector< Block * > subtree;

  // recursively collect the whole subtree
  n_blocks += collect_subtree( leaf , subtree , n_nodes );
  subtree.shrink_to_fit();

  nodes_subtrees.push_back( std::move( subtree ) );
  n_nodes++;
  }
 nodes_subtrees.shrink_to_fit();

 if( n_nodes == 1 )
  throw( std::runtime_error( "PIPSMILPSolver::load_problem: the Block needs "
	  "at least one nested Block" ) );

 DEBUG_LOG( "Number of blocks      = " << n_blocks << std::endl );
 DEBUG_LOG( "Number of nodes      = " << n_nodes << std::endl );

 // PIPMILPSolver vectors allocation - - - - - - - - - - - -  - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 n_var_node.assign( n_nodes , 0 );
 n_eq_cons_node.assign( n_nodes , 0 );
 n_ineq_cons_node.assign( n_nodes , 0 );

 var_node.assign( n_nodes , {} );
 eq_cons_node.assign( n_nodes , {} );
 ineq_cons_node.assign( n_nodes , {} );

 var_to_node.reserve( numcols );

 // First scan: scan all the variables and assign them to the correct node
 DEBUG_LOG( "First Scan: storing variables per node" << std::endl );

 Index num_node = 0;
 for( auto node : nodes_subtrees ){
  for( auto qb : node ){
    for( const auto & i : qb->get_static_variables() ){
     // Call specific function to scan the new group of Variables
     scan_group( i , qb , num_node , un_any_type< ColVariable >() );
    }

    for( const auto & i : qb->get_dynamic_variables() ){
     // Call specific function to scan the new group of Variables
     auto push_var_to_node = [ this , num_node ]
      ( const ColVariable & c ) {
      n_var_node[ num_node ] += 1;
      var_node[ num_node ].push_back( &c );
      var_to_node.emplace( &c , num_node );
     };

     un_any_const_dynamic( i , push_var_to_node ,
                            un_any_type< ColVariable >() );
    }
  }
  num_node++;
 }

 // Second scan: scan all the constraints and assign them to the correct
 // group in the correct node
 DEBUG_LOG( "Second Scan: assigning constraints per group" << std::endl );

 num_node = 0;
 for( auto node : nodes_subtrees ){
  for( auto qb : node ){
    for( const auto & i : qb->get_static_constraints() ){
     // Call specific function to scan the new group of constraints
     scan_group( i , qb , num_node , un_any_type< FRowConstraint >() );
    }

    for( const auto & i : qb->get_dynamic_constraints() ){
     // Call specific function to scan the new group of constraints
     // Scanning a group of Constraints
     auto scan = [ this , num_node ]
      ( const FRowConstraint & c ) {
        scan_constraint( c , num_node );
     };

     un_any_const_dynamic( i , scan , un_any_type< FRowConstraint >() );
    }
  }
  num_node++;
 }

 // build the PIPS input tree- - - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // every node registers the same set of callbacks, which dispatch on the
 // node id; the quadratic objective callbacks are identically zero, since
 // only LP models are supported

 auto make_node = [ this ]( int id ) {
  return( std::make_unique< DistributedInputTree::DistributedInputNode >(
	   this , id ,
	   & no_var_in_node , & no_eq_cons_in_node , & no_link_eq_cons ,
	   & no_ineq_cons_in_node , & no_link_ineq_cons ,
	   & mat_all_zero , & nnz_all_zero , & obj_vars ,
	   & mat_eq_cons_vert , & nnz_eq_cons_vert ,
	   & mat_eq_cons_diag , & nnz_eq_cons_diag ,
	   & mat_link_eq_cons , & nnz_link_eq_cons ,
	   & rhs_eq_cons , & rhs_link_eq_cons ,
	   & mat_ineq_cons_vert , & nnz_ineq_cons_vert ,
	   & mat_ineq_cons_diag , & nnz_ineq_cons_diag ,
	   & mat_link_ineq_cons , & nnz_link_ineq_cons ,
	   & lhs_ineq_cons , & flag_lhs_ineq_cons ,
	   & rhs_ineq_cons , & flag_rhs_ineq_cons ,
	   & lhs_link_ineq_cons , & flag_lhs_link_ineq_cons ,
	   & rhs_link_ineq_cons , & flag_rhs_link_ineq_cons ,
	   & lb_vars , & flag_lb_vars , & ub_vars , & flag_ub_vars ,
	   nullptr , nullptr , false ) );
  };

 // the root is node 0, with one child per leaf
 auto * root = new DistributedInputTree( make_node( 0 ) );

 for( int id = 1 ; id < static_cast< int >( n_nodes ) ; ++id )
  root->add_child(
	 std::make_unique< DistributedInputTree >( make_node( id ) ) );

 pips_tree = root;
 pips_interface = new PIPSIPMppInterface( pips_tree , MPI_COMM_WORLD );
}

/*--------------------------------------------------------------------------*/

Solver::OFValue PIPSMILPSolver::get_lb( void )
{
 // PIPS-IPM++ only solves minimization problems, hence the objective value
 // is converted back to the original objective sense of the Block; a valid
 // bound is only available upon (possibly approximate) successful
 // termination, since the iterate at which an interior-point method stops
 // early need not be dual (or primal) feasible
 if( objsense == 1 ) {  // minimization: lower bound = dual bound
  switch( sol_status ) {
   case( kOK ):
    return( constant_value + pips_interface->getObjective() );
   case( kInfeasible ):
    return( Inf< OFValue >() );
   default:  // kUnbounded, early stops, errors: no valid lower bound
    return( -Inf< OFValue >() );
   }
  }

 // maximization: lower bound = value of the best feasible solution
 switch( sol_status ) {
  case( kOK ):
  case( kLowPrecision ):
   return( constant_value - pips_interface->getObjective() );
  case( kUnbounded ):
   return( Inf< OFValue >() );
  default:  // kInfeasible, early stops, errors: no feasible solution
   return( -Inf< OFValue >() );
  }
 }

/*--------------------------------------------------------------------------*/

Solver::OFValue PIPSMILPSolver::get_ub( void )
{
 // see the comment in get_lb()
 if( objsense == 1 ) {  // minimization: upper bound = value of the best
			// feasible solution
  switch( sol_status ) {
   case( kOK ):
   case( kLowPrecision ):
    return( constant_value + pips_interface->getObjective() );
   case( kUnbounded ):
    return( -Inf< OFValue >() );
   default:  // kInfeasible, early stops, errors: no feasible solution
    return( Inf< OFValue >() );
   }
  }

 // maximization: upper bound = dual bound
 switch( sol_status ) {
  case( kOK ):
   return( constant_value - pips_interface->getObjective() );
  case( kInfeasible ):
   return( -Inf< OFValue >() );
  default:  // kUnbounded, early stops, errors: no valid upper bound
   return( Inf< OFValue >() );
  }
 }

/*--------------------------------------------------------------------------*/

int PIPSMILPSolver::guts_of_compute( void )
{
 // locking, process_modifications() and the LP cut separation loop (when
 // intRelaxIntVars == 2) are all handled by MILPSolver::compute(); this
 // method is only responsible for the actual PIPS-IPM++ call
 sol_status = decode_pips_status( pips_interface->run() );

 return( sol_status );
}

/*--------------------------------------------------------------------------*/

int PIPSMILPSolver::decode_pips_status( TerminationStatus status )
{
 DEBUG_LOG( "pips_interface.run() returned " << static_cast< int >( status )
              << std::endl );

 /* The following are the symbols that may represent the status of
 * a PIPS solution as returned by pips_interface.run(). */

 switch( status ) {
  case( TerminationStatus::READ_ERROR ):
  case( TerminationStatus::DID_NOT_RUN ):
  case( TerminationStatus::NOT_FINISHED ):
  case( TerminationStatus::STOPPED_AFTER_PRESOLVE ):
   // some error happened
   return( kError );
  case( TerminationStatus::SLOW_CONVERGENCE ):
   // terminated with a usable, but not fully accurate, solution
   return( kLowPrecision );
  case( TerminationStatus::TIMELIMIT ):
   // time limit exceeded
   return( kStopTime );
  case( TerminationStatus::INFEASIBLE ):
   // problem is infeasible
   return( kInfeasible );
  case( TerminationStatus::UNBOUNDED ):
   // problem has been proven unbounded
   return( kUnbounded );
  case( TerminationStatus::MAX_ITS_EXCEEDED ):
   // iteration limit has been reached
   return( kStopIter );
  case( TerminationStatus::SUCCESSFUL_TERMINATION ):
   // solve terminated successfully
   return( kOK );
  default:;
  }

 throw( std::runtime_error(
	 "PIPSMILPSolver::decode_pips_status: unknown termination status " +
	 std::to_string( static_cast< int >( status ) ) ) );
}

/*--------------------------------------------------------------------------*/

void PIPSMILPSolver::get_var_solution( Configuration * solc )
{
 if( ! pips_interface )
  throw( std::runtime_error( "PIPSMILPSolver::get_var_solution: the PIPS "
			     "problem has not been loaded" ) );

 auto primalSolVec = pips_interface->gatherPrimalSolution();
 int rank = 0;
 MPI_Comm_rank( MPI_COMM_WORLD , &rank );
 if( rank != 0 )
  return;

 std::vector< double > x( numcols , 0.0 );
 std::size_t pips_pos = 0;

 for( Index id = 0 ; id < n_nodes ; ++id )
  for( auto * var : var_node[ id ] ) {
   if( pips_pos >= primalSolVec.size() )
    throw( std::runtime_error( "PIPSMILPSolver::get_var_solution: the PIPS "
			       "primal solution vector is too short" ) );

   auto col = index_of_variable( var );
   if( col < 0 || col >= numcols )
    throw( std::runtime_error( "PIPSMILPSolver::get_var_solution: a variable "
			       "in the PIPS primal solution has no column" ) );

   x[ col ] = primalSolVec[ pips_pos++ ];
  }

 if( pips_pos != primalSolVec.size() )
  throw( std::runtime_error( "PIPSMILPSolver::get_var_solution: the PIPS "
			     "primal solution vector has extra entries" ) );

 MILPSolver::write_var_solution( x );
}

/*--------------------------------------------------------------------------*/

bool PIPSMILPSolver::has_dual_solution( void )
{
 return( pips_interface &&
	 ( ( sol_status == kOK ) || ( sol_status == kLowPrecision ) ) );
}

/*--------------------------------------------------------------------------*/

bool PIPSMILPSolver::is_dual_feasible( void )
{
 return( pips_interface && ( sol_status == kOK ) );
}

/*--------------------------------------------------------------------------*/

void PIPSMILPSolver::get_dual_solution( Configuration * solc )
{
 if( ! pips_interface )
  throw( std::runtime_error( "PIPSMILPSolver::get_dual_solution: the PIPS "
			     "problem has not been loaded" ) );

 auto eq_dual = pips_interface->gatherDualSolutionEq();
 auto ineq_dual = pips_interface->gatherDualSolutionIneq();
 auto var_bound_dual = pips_interface->gatherDualSolutionVarBounds();

 int rank = 0;
 MPI_Comm_rank( MPI_COMM_WORLD , &rank );
 if( rank != 0 )
  return;

 std::vector< double > pi( numrows , 0.0 );
 std::vector< double > rc( numcols , 0.0 );
 // PIPS always receives a minimization model; rescale duals back to the
 // original SMS++ objective sense.
 const double dual_scale = objsense;

 auto scatter_rows = [ this , & pi , dual_scale ](
  const std::vector< double > & source , std::size_t & pips_pos ,
  const std::vector< const FRowConstraint * > & rows , const char * name ) {

  const auto global_rows = compute_cons_global_idxs( rows );

  for( auto row : global_rows ) {
   if( pips_pos >= source.size() )
    throw( std::runtime_error(
	    std::string( "PIPSMILPSolver::get_dual_solution: PIPS " ) + name +
	    " dual vector is too short" ) );

   if( row < 0 || row >= numrows )
    throw( std::runtime_error(
	    std::string( "PIPSMILPSolver::get_dual_solution: PIPS " ) + name +
	    " dual row has no MILPSolver row" ) );

   pi[ row ] = dual_scale * source[ pips_pos++ ];
  }
 };

 // PIPS gathers node-local rows first, then appends the global linking rows.
 std::size_t pips_pos = 0;
 for( Index id = 0 ; id < n_nodes ; ++id )
  scatter_rows( eq_dual , pips_pos , eq_cons_node[ id ] , "equality" );
 scatter_rows( eq_dual , pips_pos , link_eq_cons , "linking equality" );

 if( pips_pos != eq_dual.size() )
  throw( std::runtime_error( "PIPSMILPSolver::get_dual_solution: PIPS "
			     "equality dual vector has unexpected extra "
			     "entries" ) );

 pips_pos = 0;
 for( Index id = 0 ; id < n_nodes ; ++id )
  scatter_rows( ineq_dual , pips_pos , ineq_cons_node[ id ] , "inequality" );
 scatter_rows( ineq_dual , pips_pos , link_ineq_cons , "linking inequality" );

 if( pips_pos != ineq_dual.size() )
  throw( std::runtime_error( "PIPSMILPSolver::get_dual_solution: PIPS "
			     "inequality dual vector has unexpected extra "
			     "entries" ) );

 pips_pos = 0;
 for( Index id = 0 ; id < n_nodes ; ++id )
  for( auto * var : var_node[ id ] ) {
   if( pips_pos >= var_bound_dual.size() )
    throw( std::runtime_error( "PIPSMILPSolver::get_dual_solution: PIPS "
			       "reduced-cost vector is too short" ) );

   auto col = index_of_variable( var );
   if( col < 0 || col >= numcols )
    throw( std::runtime_error( "PIPSMILPSolver::get_dual_solution: PIPS "
			       "reduced-cost variable has no column" ) );

   rc[ col ] = dual_scale * var_bound_dual[ pips_pos++ ];
  }

 if( pips_pos != var_bound_dual.size() )
  throw( std::runtime_error( "PIPSMILPSolver::get_dual_solution: PIPS "
			     "reduced-cost vector has unexpected extra "
			     "entries" ) );

 MILPSolver::write_dual_solution( pi , rc );
}

/*--------------------------------------------------------------------------*/

void PIPSMILPSolver::var_modification( const VariableMod * )
{
 f_reset = true;
}

/*--------------------------------------------------------------------------*/

void PIPSMILPSolver::objective_modification( const ObjectiveMod * )
{
 f_reset = true;
}

/*--------------------------------------------------------------------------*/

void PIPSMILPSolver::const_modification( const ConstraintMod * )
{
 f_reset = true;
}

/*--------------------------------------------------------------------------*/

void PIPSMILPSolver::bound_modification( const OneVarConstraintMod * )
{
 f_reset = true;
}

/*--------------------------------------------------------------------------*/

void PIPSMILPSolver::objective_function_modification( const FunctionMod * )
{
 f_reset = true;
}

/*--------------------------------------------------------------------------*/

void PIPSMILPSolver::constraint_function_modification( const FunctionMod * )
{
 f_reset = true;
}

/*--------------------------------------------------------------------------*/

void PIPSMILPSolver::objective_fvars_modification( const FunctionModVars * )
{
 f_reset = true;
}

/*--------------------------------------------------------------------------*/

void PIPSMILPSolver::constraint_fvars_modification( const FunctionModVars * )
{
 f_reset = true;
}

/*--------------------------------------------------------------------------*/

void PIPSMILPSolver::dynamic_modification( const BlockModAD * )
{
 f_reset = true;
}

/*--------------------------------------------------------------------------*/

std::string PIPSMILPSolver::pips_int_par_map( idx_type par ) const
{
 switch( par ) {
  case( intMaxIter ): return( "IPM_MAX_ITER" );
  case( intLogVerb ): return( "SILENT" );
  }

 // PIPS parameters
 if( ( par >= intFirstPIPSPar ) && ( par < intLastAlgParPIPS ) ) {
  return( SMSpp_to_PIPS_int_pars[ par - intFirstPIPSPar ] );

  }

 return( "" );
 }

/*--------------------------------------------------------------------------*/

std::string PIPSMILPSolver::pips_dbl_par_map( idx_type par ) const
{
 switch( par ) {
  case( dblMaxTime ): return( "IPM_TIMELIMIT" );
  case( dblRelAcc ):  return( "OUTER_BICG_TOL" );
  }

 if( ( par >= dblFirstPIPSPar ) && ( par < dblLastAlgParPIPS ) )
  return( SMSpp_to_PIPS_dbl_pars[ par - dblFirstPIPSPar ] );

 return( "" );
 }

/*--------------------------------------------------------------------------*/
/*------------------- METHODS FOR HANDLING THE PARAMETERS ------------------*/
/*--------------------------------------------------------------------------*/

void PIPSMILPSolver::set_par( idx_type par , int value )
{
 // mirror intLogVerb into MILPSolver::log_verbosity (for the LP cut
 // separation loop logging) before letting PIPS consume it through the
 // mapping below
 if( par == intLogVerb )
  log_verbosity = value;

 std::string Pp = pips_int_par_map( par );
 if( Pp.size() > 0 ) {
  // the PIPS SILENT option has the opposite meaning of intLogVerb, hence
  // the value is inverted
  if( Pp == "SILENT" ) {
   pipsipmpp::options::set_parameter( Pp , value == 0 );
   return;
   }

  // the mapped option may be a bool one in PIPS, in which case setting it
  // as an int throws: retry setting it as a bool
  if( ( value == 0 ) || ( value == 1 ) ) {
   try {
    pipsipmpp::options::set_parameter( Pp , value );
    }
   catch( const std::runtime_error & ) {
    pipsipmpp::options::set_parameter( Pp , value != 0 );
    }
   return;
   }

  pipsipmpp::options::set_parameter( Pp , value );
  return;
  }

 MILPSolver::set_par( par , value );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

void PIPSMILPSolver::set_par( idx_type par , double value )
{
 std::string Pp = pips_dbl_par_map( par );
 if( Pp.size() > 0 ) {
  pipsipmpp::options::set_parameter( Pp , value );
  return;
  }

 MILPSolver::set_par( par , value );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

void PIPSMILPSolver::set_par( idx_type par , std::string && value )
{
 // PIPS parameters
 if( ( par >= strFirstPIPSPar ) && ( par < strLastAlgParPIPS ) ) {
  std::string pips_par = SMSpp_to_PIPS_str_pars[ par - strFirstPIPSPar ];
  pipsipmpp::options::set_parameter( pips_par , value );
  return;
  }

 // PIPS has no log-file option, but strLogFileName is still forwarded to
 // the base class so that the standard bookkeeping is performed; the
 // problem can be written to file via strOutputFile
 if( par == strOutputFile )
  pipsipmpp::options::set_parameter( "WRITE_ORIGINAL_PROBLEM_TO_LP" ,
				     value );

 MILPSolver::set_par( par , std::move( value ) );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

Solver::idx_type PIPSMILPSolver::get_num_int_par( void ) const {
 return( MILPSolver::get_num_int_par()
	 + intLastAlgParPIPS - intLastAlgParMILP );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

Solver::idx_type PIPSMILPSolver::get_num_dbl_par( void ) const {
 return( MILPSolver::get_num_dbl_par()
	 + dblLastAlgParPIPS - dblLastAlgParMILP );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

Solver::idx_type PIPSMILPSolver::get_num_str_par( void ) const {
 return( MILPSolver::get_num_str_par()
	 + strLastAlgParPIPS - strLastAlgParMILP );
 }

/*--------------------------------------------------------------------------*/

int PIPSMILPSolver::get_dflt_int_par( idx_type par ) const
{
 // note: for parameters mapped into PIPS options the current value of the
 // (process-global) option is returned, as PIPS does not expose defaults
 if( ( par >= intFirstPIPSPar ) && ( par < intLastAlgParPIPS ) )
  return( pipsipmpp::options::get_int_parameter(
			     SMSpp_to_PIPS_int_pars[ par - intFirstPIPSPar ]
			     ) );

 return( MILPSolver::get_dflt_int_par( par ) );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

double PIPSMILPSolver::get_dflt_dbl_par( idx_type par ) const
{
 // see the comment in get_dflt_int_par()
 if( ( par >= dblFirstPIPSPar ) && ( par < dblLastAlgParPIPS ) )
  return( pipsipmpp::options::get_double_parameter(
			     SMSpp_to_PIPS_dbl_pars[ par - dblFirstPIPSPar ]
			     ) );

 return( MILPSolver::get_dflt_dbl_par( par ) );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

const std::string & PIPSMILPSolver::get_dflt_str_par( idx_type par ) const
{
 // see the comment in get_dflt_int_par(); note that this is not thread safe
 static std::string value;

 if( ( par >= strFirstPIPSPar ) && ( par < strLastAlgParPIPS ) ) {
  value = pipsipmpp::options::get_string_parameter(
			     SMSpp_to_PIPS_str_pars[ par - strFirstPIPSPar ] );
  return( value );
  }

 return( MILPSolver::get_dflt_str_par( par ) );
 }

/*--------------------------------------------------------------------------*/

int PIPSMILPSolver::get_int_par( idx_type par ) const
{
 // intLogVerb is mirrored into the (inverted) PIPS SILENT option by
 // set_par(), hence the value kept in the base class is returned
 if( par == intLogVerb )
  return( log_verbosity );

 std::string Pp = pips_int_par_map( par );
 if( Pp.size() > 0 )
  return( pipsipmpp::options::get_int_parameter( Pp ) );

 return( MILPSolver::get_int_par( par ) );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

double PIPSMILPSolver::get_dbl_par( idx_type par ) const
{
 std::string Pp = pips_dbl_par_map( par );
 if( Pp.size() > 0 )
  return( pipsipmpp::options::get_double_parameter( Pp ) );

 return( MILPSolver::get_dbl_par( par ) );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

const std::string & PIPSMILPSolver::get_str_par( idx_type par ) const
{
 // note: static buffer, not thread safe
 static std::string value;

 if( ( par >= strFirstPIPSPar ) && ( par < strLastAlgParPIPS ) ) {
  value = pipsipmpp::options::get_string_parameter(
			     SMSpp_to_PIPS_str_pars[ par - strFirstPIPSPar ] );
  return( value );
  }

 return( MILPSolver::get_str_par( par ) );
 }

/*--------------------------------------------------------------------------*/

Solver::idx_type PIPSMILPSolver::int_par_str2idx(
					     const std::string & name ) const
{
 // check with MILPSolver first
 idx_type idx = MILPSolver::int_par_str2idx( name );
 if( idx < Inf< idx_type >() )
  return( idx );

 // PIPS parameters
 auto array_pos = std::find( SMSpp_to_PIPS_int_pars.begin() ,
			     SMSpp_to_PIPS_int_pars.end() , name );

 if( array_pos != SMSpp_to_PIPS_int_pars.end() ) {
  int pos = std::distance( SMSpp_to_PIPS_int_pars.begin() , array_pos );
  return( PIPS_to_SMSpp_int_pars[ pos ].second );
  }

 return( Inf< idx_type >() );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

const std::string & PIPSMILPSolver::int_par_idx2str( idx_type idx ) const
{
 // note: static buffer, not thread safe; the result must be used
 //       immediately after the call (prior to any other call to
 //       int_par_idx2str()), this may have to be improved upon
 static std::string par_name;

 if( ( idx >= intFirstPIPSPar ) && ( idx < intLastAlgParPIPS ) ) {
  par_name = SMSpp_to_PIPS_int_pars[ idx - intFirstPIPSPar ];
  return( par_name );
  }

 return( MILPSolver::int_par_idx2str( idx ) );
 }

/*--------------------------------------------------------------------------*/

Solver::idx_type PIPSMILPSolver::dbl_par_str2idx( const std::string & name )
 const
{
 // check with MILPSolver first
 idx_type idx = MILPSolver::dbl_par_str2idx( name );
 if( idx < Inf< idx_type >() )
  return( idx );

 // PIPS parameters
 auto array_pos = std::find( SMSpp_to_PIPS_dbl_pars.begin() ,
			     SMSpp_to_PIPS_dbl_pars.end() , name );

 if( array_pos != SMSpp_to_PIPS_dbl_pars.end() ) {
  int pos = std::distance( SMSpp_to_PIPS_dbl_pars.begin() , array_pos );
  return( PIPS_to_SMSpp_dbl_pars[ pos ].second );
  }

 return( Inf< idx_type >() );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

const std::string & PIPSMILPSolver::dbl_par_idx2str( idx_type idx ) const
{
 // note: static buffer, not thread safe; the result must be used
 //       immediately after the call (prior to any other call to
 //       dbl_par_idx2str()), this may have to be improved upon
 static std::string par_name;

 if( ( idx >= dblFirstPIPSPar ) && ( idx < dblLastAlgParPIPS ) ) {
  par_name = SMSpp_to_PIPS_dbl_pars[ idx - dblFirstPIPSPar ];
  return( par_name );
  }

 return( MILPSolver::dbl_par_idx2str( idx ) );
 }

/*--------------------------------------------------------------------------*/

Solver::idx_type PIPSMILPSolver::str_par_str2idx( const std::string & name )
 const
{
 // check with MILPSolver first
 idx_type idx = MILPSolver::str_par_str2idx( name );
 if( idx < Inf< idx_type >() )
  return( idx );

 // PIPS parameters
 auto array_pos = std::find( SMSpp_to_PIPS_str_pars.begin() ,
			     SMSpp_to_PIPS_str_pars.end() , name );

 if( array_pos != SMSpp_to_PIPS_str_pars.end() ) {
  int pos = std::distance( SMSpp_to_PIPS_str_pars.begin() , array_pos );
  return( PIPS_to_SMSpp_str_pars[ pos ].second );
  }

 return( Inf< idx_type >() );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

const std::string & PIPSMILPSolver::str_par_idx2str( idx_type idx ) const
{
 // note: static buffer, not thread safe; the result must be used
 //       immediately after the call (prior to any other call to
 //       str_par_idx2str()), this may have to be improved upon
 static std::string par_name;

 if( ( idx >= strFirstPIPSPar ) && ( idx < strLastAlgParPIPS ) ) {
  par_name = SMSpp_to_PIPS_str_pars[ idx - strFirstPIPSPar ];
  return( par_name );
  }

 return( MILPSolver::str_par_idx2str( idx ) );
 }

/*--------------------------------------------------------------------------*/
/*----------------------- PIPS CALLBACKS METHODS ---------------------------*/
/*--------------------------------------------------------------------------*/

PIPSMILPSolver::CSRMatrix PIPSMILPSolver::extract_block_matrix(
 const std::vector< const FRowConstraint * > & rows ,
 const std::vector< const ColVariable * > & cols
 ) const
{
 CSRMatrix mat;

 std::unordered_map< const ColVariable * , int > local_col;

 for( int j = 0 ; j < static_cast< int >( cols.size() ) ; ++j )
  local_col[ cols[ j ] ] = j;

 mat.krow.reserve( rows.size() + 1 );
 mat.krow.push_back( 0 );

 for( const auto * con : rows ) {
  const auto * f = con->get_function();

  if( const auto * lf = dynamic_cast< const LinearFunction * >( f ) ) {
   for( const auto & el : lf->get_v_var() ) {
    const auto * var =
     dynamic_cast< const ColVariable * >( std::get< 0 >( el ) );

    auto it = local_col.find( var );

    if( it != local_col.end() ) {
     mat.jcol.push_back( it->second );
     mat.val.push_back( std::get< 1 >( el ) );
     }
    }
   }

  mat.krow.push_back( mat.val.size() );
  }

 return( mat );
}

/*--------------------------------------------------------------------------*/

PIPSMILPSolver::CSRMatrix PIPSMILPSolver::extract_submatrix_to_csr(
 const std::vector< int > & selected_rows ,
 const int n_rows ,
 const std::vector< int > & selected_cols ,
 const int n_cols ) const
{
 CSRMatrix result;

 std::unordered_map< int , int > row_map;
 row_map.reserve( n_rows );

 for( int i = 0 ; i < n_rows ; ++i )
  row_map[ selected_rows[ i ] ] = i;

 std::vector< std::vector< std::pair< int , double > > > rows( n_rows );

 for( int local_col = 0 ; local_col < n_cols ; ++local_col ) {
  const int original_col = selected_cols[ local_col ];

  for( int k = matbeg[ original_col ] ;
       k < matbeg[ original_col ] + matcnt[ original_col ] ; ++k ) {
   const int original_row = matind[ k ];
   auto it = row_map.find( original_row );

   if( it != row_map.end() )
    rows[ it->second ].push_back( { local_col , matval[ k ] } );
   }
  }

 result.krow.reserve( n_rows + 1 );
 result.krow.push_back( 0 );

 for( int i = 0 ; i < n_rows ; ++i ) {
  for( const auto & el : rows[ i ] ) {
   result.jcol.push_back( el.first );
   result.val.push_back( el.second );
   }

  result.krow.push_back( result.val.size() );
  }

 return( result );
}

/*--------------------------------------------------------------------------*/

int PIPSMILPSolver::extract_matrix(
 int id , int * krowM , int * jcolM , double * M ,
 const std::vector< const FRowConstraint * > & node_cons ,
 const std::vector< const ColVariable * > & vars ) const
{
 if( ( id < 0 ) || ( id >= static_cast< int >( n_nodes ) ) )
  throw( std::runtime_error( "PIPSMILPSolver::extract_matrix: node ID "
			     "outside the expected range" ) );

 const auto n_cons = static_cast< int >( node_cons.size() );
 const auto n_vars = static_cast< int >( vars.size() );

 if( ( n_cons == 0 ) || ( n_vars == 0 ) ) {
  for( int i = 0 ; i <= n_cons ; ++i )
   krowM[ i ] = 0;

  return( 0 );
  }

 auto global_idxs_cons = compute_cons_global_idxs( node_cons );
 auto global_idxs_vars = compute_vars_global_idxs( vars );
 auto sub_matrix = extract_submatrix_to_csr( global_idxs_cons , n_cons ,
					     global_idxs_vars , n_vars );

 std::copy( sub_matrix.krow.begin() , sub_matrix.krow.end() , krowM );
 std::copy( sub_matrix.jcol.begin() , sub_matrix.jcol.end() , jcolM );
 std::copy( sub_matrix.val.begin() , sub_matrix.val.end() , M );

 return( 0 );
}

/*--------------------------------------------------------------------------*/

int PIPSMILPSolver::evaluate_nnz(
 int id , int * nnz ,
 const std::vector< const FRowConstraint * > & node_cons ,
 const std::vector< const ColVariable * > & vars ) const
{
 if( ( id < 0 ) || ( id >= static_cast< int >( n_nodes ) ) )
  throw( std::runtime_error( "PIPSMILPSolver::evaluate_nnz: node ID "
			     "outside the expected range" ) );

 const auto n_cons = static_cast< int >( node_cons.size() );
 const auto n_vars = static_cast< int >( vars.size() );

 if( ( n_cons == 0 ) || ( n_vars == 0 ) ) {
  *nnz = 0;
  return( 0 );
  }

 auto global_idxs_cons = compute_cons_global_idxs( node_cons );
 auto global_idxs_vars = compute_vars_global_idxs( vars );
 auto sub_matrix = extract_submatrix_to_csr( global_idxs_cons , n_cons ,
					     global_idxs_vars , n_vars );

 *nnz = sub_matrix.nnz();

 return( 0 );
}

/*--------------------------------------------------------------------------*/

std::vector< int > PIPSMILPSolver::compute_cons_global_idxs(
 const std::vector< const FRowConstraint * > & cons ) const
{
 std::vector< int > global_idxs( cons.size() , -1 );
 Index row = 0;

 for( auto con : cons ) {
  Index static_idx = index_of_static_constraint( con );

  if( static_idx < Inf< Index >() )
   global_idxs[ row ] = static_idx;
  else {
   Index dynamic_idx = index_of_dynamic_constraint( con );
   if( dynamic_idx < Inf< Index >() )
    global_idxs[ row ] = dynamic_idx;
   else
    throw( std::runtime_error(
	    "PIPSMILPSolver::compute_cons_global_idxs: the required "
	    "constraint has no available index" ) );
   }
  row++;
  }

 return( global_idxs );
}

/*--------------------------------------------------------------------------*/

std::vector< int > PIPSMILPSolver::compute_vars_global_idxs(
 const std::vector< const ColVariable * > & vars ) const
{
 std::vector< int > global_idxs( vars.size() , -1 );
 Index col = 0;

 for( auto var : vars ) {
  Index static_idx = index_of_static_variable( var );

  if( static_idx < Inf< Index >() )
   global_idxs[ col ] = static_idx;
  else {
   Index dynamic_idx = index_of_dynamic_variable( var );
   if( dynamic_idx < Inf< Index >() )
    global_idxs[ col ] = dynamic_idx;
   else
    throw( std::runtime_error(
	    "PIPSMILPSolver::compute_vars_global_idxs: the required "
	    "variable has no available index" ) );
   }
  col++;
  }

 return( global_idxs );
}

/*--------------------------------------------------------------------------*/

Block::Index PIPSMILPSolver::collect_subtree( Block * block ,
					      std::vector< Block * > & subtree ,
					      Index parent_leaf )
{
 // map the Block to its father leaf
 block_to_leaf[ block ] = parent_leaf;

 Index tree_size = 1;
 subtree.push_back( block );

 for( auto * child : block->get_nested_Blocks() )
  tree_size += collect_subtree( child , subtree , parent_leaf );

 return( tree_size );
}

/*--------------------------------------------------------------------------*/

template< typename T >
 void PIPSMILPSolver::scan_group( const boost::any & gr , Block * qb ,
				  Index num_node , un_any_type< T > )
{
 // search for the group type
 if( ( gr.type() == typeid( T * ) ) ||
     ( gr.type() == typeid( std::vector< T > * ) ) ||
     ( gr.type() == typeid( std::vector< std::vector< T > > * ) ) ) {
  // "simple" group
  scan_simple_group( gr , qb , num_node , un_any_type< T >() );
  }
 else {
  // "complex" group
  scan_multiarray_group( gr , qb , num_node , un_any_type< T >() );
  }
 }

/*--------------------------------------------------------------------------*/

template< typename T >
 void PIPSMILPSolver::scan_simple_group( const boost::any & gr , Block * qb ,
					 Index num_node , un_any_type< T > )
{
 if( typeid( T * ) == typeid( FRowConstraint * ) ) {
  // constraint group
  auto scan = [ this , num_node ]( const FRowConstraint & c ) {
   scan_constraint( c , num_node );
   };
  // scan all the constraints one at a time
  un_any_const_static( gr , scan , un_any_type< FRowConstraint >() );
  }
 else if( typeid( T * ) == typeid( ColVariable * ) ) {
  // variable group: simply add each variable to the corresponding node
  auto push_var_to_node = [ this , num_node ]( const ColVariable & c ) {
   n_var_node[ num_node ] += 1;
   var_node[ num_node ].push_back( & c );
   var_to_node.emplace( & c , num_node );
   };
  un_any_const_static( gr , push_var_to_node , un_any_type< ColVariable >() );
  }
 else
  throw( std::runtime_error( "PIPSMILPSolver::scan_simple_group: "
			     "unsupported group type" ) );
}

/*--------------------------------------------------------------------------*/

template< typename T >
 void PIPSMILPSolver::scan_multiarray_group( const boost::any & gr ,
	    Block * qb , Index num_node , un_any_type< T > )
{
 // get the multi_array dimensionality
 int ma_dim = get_multi_array_dim( gr , un_any_type< T >() ,
				   un_any_int< 2 >() );

 if( ma_dim == 2 ) {
  // 2D multi_array

  // get the multi_array type: see get_multi_array_type() in MILPSolver.h
  int type = get_multi_array_type( gr , un_any_type< T >() ,
				   un_any_int< 2 >() );

  // indices of the 2 dimensions
  int idx_0 = 0;
  int idx_1 = 0;

  if( type == 1 ) {
   // multi_arrays of type 1 (i.e. multi_array< std::vector< T * > >) do not
   // store elements in sequential cells: each std::vector is "unpacked" and
   // stored separately
   auto ma = get_multi_array1( gr , un_any_type< T >() ,
			       un_any_int< 2 >() );

   const auto dim_0 = static_cast< int >( ma->shape()[ 0 ] );
   const auto dim_1 = static_cast< int >( ma->shape()[ 1 ] );

   if( ( dim_0 == 0 ) || ( dim_1 == 0 ) )
    return;  // empty multi_array: nothing to scan

   if( typeid( T * ) == typeid( FRowConstraint * ) ) {
    // constraint group

    // scan the linearization of the array
    for( auto v = ma->data() ; idx_0 < dim_0 ; ++v ) {
     auto scan = [ this , num_node ]( const FRowConstraint & c ) {
      scan_constraint( c , num_node );
      };
     // scan all the constraints one at a time
     un_any_const_static( v , scan , un_any_type< FRowConstraint >() );

     // the linearization produced by ma->data() for the 2D multi_array
     // stores elements in row-major order: increment the column index
     // first, and when it exceeds the number of columns reset it and
     // increment the row index
     if( idx_1 < dim_1 - 1 )
      idx_1++;
     else {
      idx_1 = 0;
      idx_0++;
      }
     }
    }
   else if( typeid( T * ) == typeid( ColVariable * ) ) {
    // variable group

    // scan the linearization of the array
    for( auto v = ma->data() ; idx_0 < dim_0 ; ++v ) {
     // simply add each variable to the corresponding node
     auto push_var_to_node = [ this , num_node ]( const ColVariable & c ) {
      n_var_node[ num_node ] += 1;
      var_node[ num_node ].push_back( & c );
      var_to_node.emplace( & c , num_node );
      };
     un_any_const_static( v , push_var_to_node ,
			  un_any_type< ColVariable >() );

     // row-major order, see above
     if( idx_1 < dim_1 - 1 )
      idx_1++;
     else {
      idx_1 = 0;
      idx_0++;
      }
     }
    }
   else
    throw( std::runtime_error( "PIPSMILPSolver::scan_multiarray_group: "
			       "unsupported group type" ) );
   }
  else if( type == 0 ) {
   // multi_arrays of type 0 (i.e. multi_array< T >) store elements in
   // sequential cells, hence they can be stored as done for
   // std::vector< T > by only keeping track of the first element and
   // storing the number of non-empty cells in the structure
   auto ma = get_multi_array0( gr , un_any_type< T >() ,
			       un_any_int< 2 >() );

   const auto dim_0 = static_cast< int >( ma->shape()[ 0 ] );
   const auto dim_1 = static_cast< int >( ma->shape()[ 1 ] );

   if( ( dim_0 == 0 ) || ( dim_1 == 0 ) )
    return;  // empty multi_array: nothing to scan

   if( typeid( T * ) == typeid( FRowConstraint * ) ) {
    // constraint group

    // scan the linearization of the array
    for( auto v = ma->data() ; idx_0 < dim_0 ; ++v ) {
     auto scan = [ this , num_node ]( const FRowConstraint & c ) {
      scan_constraint( c , num_node );
      };
     // scan all the constraints one at a time
     un_any_const_static( v , scan , un_any_type< FRowConstraint >() );

     // row-major order, see above
     if( idx_1 < dim_1 - 1 )
      idx_1++;
     else {
      idx_1 = 0;
      idx_0++;
      }
     }
    }
   else if( typeid( T * ) == typeid( ColVariable * ) ) {
    // variable group

    // scan the linearization of the array
    for( auto v = ma->data() ; idx_0 < dim_0 ; ++v ) {
     // simply add each variable to the corresponding node
     auto push_var_to_node = [ this , num_node ]( const ColVariable & c ) {
      n_var_node[ num_node ] += 1;
      var_node[ num_node ].push_back( & c );
      var_to_node.emplace( & c , num_node );
      };
     un_any_const_static( v , push_var_to_node ,
			  un_any_type< ColVariable >() );

     // row-major order, see above
     if( idx_1 < dim_1 - 1 )
      idx_1++;
     else {
      idx_1 = 0;
      idx_0++;
      }
     }
    }
   else
    throw( std::runtime_error( "PIPSMILPSolver::scan_multiarray_group: "
			       "unsupported group type" ) );
   }
  else
   throw( std::runtime_error( "PIPSMILPSolver::scan_multiarray_group: "
			      "unsupported multi_array type" ) );
  }
 else if( ma_dim == 3 ) {
  // 3D multi_array

  // get the multi_array type: see get_multi_array_type() in MILPSolver.h
  int type = get_multi_array_type( gr , un_any_type< T >() ,
				   un_any_int< 3 >() );

  // indices of the 3 dimensions
  int idx_0 = 0;
  int idx_1 = 0;
  int idx_2 = 0;

  if( type == 1 ) {
   // multi_arrays of type 1, see above
   auto ma = get_multi_array1( gr , un_any_type< T >() ,
			       un_any_int< 3 >() );

   const auto dim_0 = static_cast< int >( ma->shape()[ 0 ] );
   const auto dim_1 = static_cast< int >( ma->shape()[ 1 ] );
   const auto dim_2 = static_cast< int >( ma->shape()[ 2 ] );

   if( ( dim_0 == 0 ) || ( dim_1 == 0 ) || ( dim_2 == 0 ) )
    return;  // empty multi_array: nothing to scan

   if( typeid( T * ) == typeid( FRowConstraint * ) ) {
    // constraint group

    // scan the linearization of the array
    for( auto v = ma->data() ; idx_0 < dim_0 ; ++v ) {
     auto scan = [ this , num_node ]( const FRowConstraint & c ) {
      scan_constraint( c , num_node );
      };
     // scan all the constraints one at a time
     un_any_const_static( v , scan , un_any_type< FRowConstraint >() );

     // the linearization produced by ma->data() for the 3D multi_array
     // stores elements in row-major order
     if( idx_2 < dim_2 - 1 )
      idx_2++;
     else if( idx_1 < dim_1 - 1 ) {
      idx_1++;
      idx_2 = 0;
      }
     else {
      idx_0++;
      idx_1 = 0;
      idx_2 = 0;
      }
     }
    }
   else if( typeid( T * ) == typeid( ColVariable * ) ) {
    // variable group

    // scan the linearization of the array
    for( auto v = ma->data() ; idx_0 < dim_0 ; ++v ) {
     // simply add each variable to the corresponding node
     auto push_var_to_node = [ this , num_node ]( const ColVariable & c ) {
      n_var_node[ num_node ] += 1;
      var_node[ num_node ].push_back( & c );
      var_to_node.emplace( & c , num_node );
      };
     un_any_const_static( v , push_var_to_node ,
			  un_any_type< ColVariable >() );

     // row-major order, see above
     if( idx_2 < dim_2 - 1 )
      idx_2++;
     else if( idx_1 < dim_1 - 1 ) {
      idx_1++;
      idx_2 = 0;
      }
     else {
      idx_0++;
      idx_1 = 0;
      idx_2 = 0;
      }
     }
    }
   else
    throw( std::runtime_error( "PIPSMILPSolver::scan_multiarray_group: "
			       "unsupported group type" ) );
   }
  else if( type == 0 ) {
   // multi_arrays of type 0, see above
   auto ma = get_multi_array0( gr , un_any_type< T >() ,
			       un_any_int< 3 >() );

   const auto dim_0 = static_cast< int >( ma->shape()[ 0 ] );
   const auto dim_1 = static_cast< int >( ma->shape()[ 1 ] );
   const auto dim_2 = static_cast< int >( ma->shape()[ 2 ] );

   if( ( dim_0 == 0 ) || ( dim_1 == 0 ) || ( dim_2 == 0 ) )
    return;  // empty multi_array: nothing to scan

   if( typeid( T * ) == typeid( FRowConstraint * ) ) {
    // constraint group

    // scan the linearization of the array
    for( auto v = ma->data() ; idx_0 < dim_0 ; ++v ) {
     auto scan = [ this , num_node ]( const FRowConstraint & c ) {
      scan_constraint( c , num_node );
      };
     // scan all the constraints one at a time
     un_any_const_static( v , scan , un_any_type< FRowConstraint >() );

     // row-major order, see above
     if( idx_2 < dim_2 - 1 )
      idx_2++;
     else if( idx_1 < dim_1 - 1 ) {
      idx_1++;
      idx_2 = 0;
      }
     else {
      idx_0++;
      idx_1 = 0;
      idx_2 = 0;
      }
     }
    }
   else if( typeid( T * ) == typeid( ColVariable * ) ) {
    // variable group

    // scan the linearization of the array
    for( auto v = ma->data() ; idx_0 < dim_0 ; ++v ) {
     // simply add each variable to the corresponding node
     auto push_var_to_node = [ this , num_node ]( const ColVariable & c ) {
      n_var_node[ num_node ] += 1;
      var_node[ num_node ].push_back( & c );
      var_to_node.emplace( & c , num_node );
      };
     un_any_const_static( v , push_var_to_node ,
			  un_any_type< ColVariable >() );

     // row-major order, see above
     if( idx_2 < dim_2 - 1 )
      idx_2++;
     else if( idx_1 < dim_1 - 1 ) {
      idx_1++;
      idx_2 = 0;
      }
     else {
      idx_0++;
      idx_1 = 0;
      idx_2 = 0;
      }
     }
    }
   else
    throw( std::runtime_error( "PIPSMILPSolver::scan_multiarray_group: "
			       "unsupported group type" ) );
   }
  else
   throw( std::runtime_error( "PIPSMILPSolver::scan_multiarray_group: "
			      "unsupported multi_array type" ) );
  }
 else
  throw( std::runtime_error(
	  "PIPSMILPSolver::scan_multiarray_group: unsupported multi_array "
	  "dimensionality " + std::to_string( ma_dim ) ) );
}

/*--------------------------------------------------------------------------*/

void PIPSMILPSolver::scan_constraint( const FRowConstraint & con ,
				      Index num_node )
{
 // check that the Constraint is not empty
 if( con.get_Block() == nullptr )
  throw( std::invalid_argument( "PIPSMILPSolver::scan_constraint: the "
				"provided constraint is empty" ) );

 // categorize the constraint into equality/inequality
 auto con_lhs = con.get_lhs();
 auto con_rhs = con.get_rhs();
 bool is_eq = ( con_lhs == con_rhs );

 // a row with lhs == rhs == +/- infinity is no equality: reject it
 if( is_eq && ( ( con_lhs == Inf< double >() ) ||
		( con_lhs == -Inf< double >() ) ) )
  throw( std::invalid_argument( "PIPSMILPSolver::scan_constraint: the "
				"provided constraint has lhs == rhs == "
				"+/- infinity" ) );

 // check whether the constraint involves only variables of the same node
 // (accepting also the root node) or also variables coming from other
 // leaves
 bool found_in_other_node = false;
 if( auto f = con.get_function() ) {
  if( auto lf = dynamic_cast< const LinearFunction * >( f ) ) {
   for( auto el : lf->get_v_var() ) {
    auto * v = dynamic_cast< ColVariable * >( std::get< 0 >( el ) );

    auto it = var_to_node.find( v );

    if( it == var_to_node.end() )
     throw( std::runtime_error( "PIPSMILPSolver::scan_constraint: variable "
				"not assigned to any PIPS node" ) );

    auto owner = static_cast< Index >( it->second );

    if( ( owner != num_node ) && ( owner != 0 ) ) {
     found_in_other_node = true;
     break;
     }
    }
   }
  else
   throw( std::invalid_argument( "PIPSMILPSolver::scan_constraint: "
				 "unexpected constraint type" ) );
  }

 // add the constraint to the correct group
 if( found_in_other_node && is_eq ) {
  // equality linking constraints
  n_link_eq_cons++;
  link_eq_cons.push_back( & con );
  }
 else if( found_in_other_node && ( ! is_eq ) ) {
  // inequality linking constraints
  n_link_ineq_cons++;
  link_ineq_cons.push_back( & con );
  }
 else if( ( ! found_in_other_node ) && is_eq ) {
  // node equality constraints
  n_eq_cons_node[ num_node ]++;
  eq_cons_node[ num_node ].push_back( & con );
  }
 else {
  // node inequality constraints
  n_ineq_cons_node[ num_node ]++;
  ineq_cons_node[ num_node ].push_back( & con );
  }
}

/*--------------------------------------------------------------------------*/

int PIPSMILPSolver::extract_rhs_vector( int id , double * vec , int len ,
	      const std::vector< const FRowConstraint * > & node_cons ,
	      const std::vector< double > & rhs ,
	      const std::vector< char > & sense ,
	      const std::vector< double > & ranges )
{
 if( ( id < 0 ) || ( id >= static_cast< int >( n_nodes ) ) )
  throw( std::runtime_error( "PIPSMILPSolver::extract_rhs_vector: node ID "
			     "outside the expected range" ) );

 if( node_cons.empty() ) {
  // PIPS may query nodes with no local entries: return a zero vector
  for( int i = 0 ; i < len ; ++i )
   vec[ i ] = 0.0;

  return( 0 );
  }

 // extract the global indices of the constraints
 std::vector< int > global_idxs_cons = compute_cons_global_idxs( node_cons );

 // extract the corresponding elements and store them
 std::vector< double > sub_rhs;
 sub_rhs.reserve( global_idxs_cons.size() );

 if( sense.empty() ) {
  // equality constraints: just select the corresponding rhs entries
  for( auto idx : global_idxs_cons ) {
   if( static_cast< std::size_t >( idx ) >= rhs.size() )
    throw( std::runtime_error( "PIPSMILPSolver::extract_rhs_vector: "
			       "constraint index outside the row data" ) );

   sub_rhs.push_back( rhs[ idx ] );
   }
  }
 else {
  // inequality constraints: distinguish among the several senses
  for( auto idx : global_idxs_cons ) {
   if( static_cast< std::size_t >( idx ) >= rhs.size() )
    throw( std::runtime_error( "PIPSMILPSolver::extract_rhs_vector: "
			       "constraint index outside the row data" ) );

   if( sense[ idx ] == 'L' )
    // <= RHS, store the element
    sub_rhs.push_back( rhs[ idx ] );
   else if( sense[ idx ] == 'G' )
    // >= LHS, so no elements to store
    sub_rhs.push_back( 0.0 );
   else if( sense[ idx ] == 'R' ) {
    // LHS <= f <= RHS, evaluate the sign of the range
    if( ranges[ idx ] > 0 )
     sub_rhs.push_back( rhs[ idx ] + ranges[ idx ] );
    else
     sub_rhs.push_back( rhs[ idx ] );
    }
   else
    throw( std::runtime_error( "PIPSMILPSolver::extract_rhs_vector: "
			       "unexpected inequality sense" ) );
   }
  }

 if( static_cast< int >( sub_rhs.size() ) != len )
  throw( std::runtime_error(
	  "PIPSMILPSolver::extract_rhs_vector: node " +
	  std::to_string( id ) + " built " +
	  std::to_string( sub_rhs.size() ) +
	  " values but PIPS requested length " + std::to_string( len ) ) );

 std::copy( sub_rhs.begin() , sub_rhs.end() , vec );

 return( 0 );
}

/*--------------------------------------------------------------------------*/

int PIPSMILPSolver::extract_lhs_vector( int id , double * vec , int len ,
	      const std::vector< const FRowConstraint * > & node_cons ,
	      const std::vector< double > & rhs ,
	      const std::vector< char > & sense ,
	      const std::vector< double > & ranges )
{
 if( ( id < 0 ) || ( id >= static_cast< int >( n_nodes ) ) )
  throw( std::runtime_error( "PIPSMILPSolver::extract_lhs_vector: node ID "
			     "outside the expected range" ) );

 if( node_cons.empty() ) {
  // zero-fill
  for( int i = 0 ; i < len ; ++i )
   vec[ i ] = 0.0;

  return( 0 );
  }

 // extract the global indices of the constraints
 std::vector< int > global_idxs_cons = compute_cons_global_idxs( node_cons );

 // extract the corresponding elements and store them
 std::vector< double > sub_rhs;
 sub_rhs.reserve( global_idxs_cons.size() );

 for( auto idx : global_idxs_cons ) {
  if( static_cast< std::size_t >( idx ) >= rhs.size() )
   throw( std::runtime_error( "PIPSMILPSolver::extract_lhs_vector: "
			      "constraint index outside the row data" ) );

  if( sense[ idx ] == 'G' )
   // >= LHS, store the element
   sub_rhs.push_back( rhs[ idx ] );
  else if( sense[ idx ] == 'L' )
   // <= RHS, so no elements to store
   sub_rhs.push_back( 0.0 );
  else if( sense[ idx ] == 'R' ) {
   // LHS <= f <= RHS, evaluate the sign of the range
   if( ranges[ idx ] > 0 )
    sub_rhs.push_back( rhs[ idx ] );
   else
    sub_rhs.push_back( rhs[ idx ] + ranges[ idx ] );
   }
  else
   throw( std::runtime_error( "PIPSMILPSolver::extract_lhs_vector: "
			      "unexpected inequality sense" ) );
  }

 if( static_cast< int >( sub_rhs.size() ) != len )
  throw( std::runtime_error(
	  "PIPSMILPSolver::extract_lhs_vector: node " +
	  std::to_string( id ) + " built " +
	  std::to_string( sub_rhs.size() ) +
	  " values but PIPS requested length " + std::to_string( len ) ) );

 std::copy( sub_rhs.begin() , sub_rhs.end() , vec );

 return( 0 );
}

/*--------------------------------------------------------------------------*/

int PIPSMILPSolver::extract_rhs_active_flag( int id , double * vec , int len ,
	      const std::vector< const FRowConstraint * > & node_cons ,
	      const std::vector< double > & rhs ,
	      const std::vector< char > & sense ,
	      const std::vector< double > & ranges )
{
 if( ( id < 0 ) || ( id >= static_cast< int >( n_nodes ) ) )
  throw( std::runtime_error( "PIPSMILPSolver::extract_rhs_active_flag: node "
			     "ID outside the expected range" ) );

 if( node_cons.empty() ) {
  // zero-fill
  for( int i = 0 ; i < len ; ++i )
   vec[ i ] = 0.0;

  return( 0 );
  }

 // extract the global indices of the constraints
 std::vector< int > global_idxs_cons = compute_cons_global_idxs( node_cons );

 // extract the corresponding elements and store them
 std::vector< double > sub_rhs;
 sub_rhs.reserve( global_idxs_cons.size() );

 if( sense.empty() ) {
  // equality constraints: just set the flag to true
  for( auto idx : global_idxs_cons ) {
   if( static_cast< std::size_t >( idx ) >= rhs.size() )
    throw( std::runtime_error( "PIPSMILPSolver::extract_rhs_active_flag: "
			       "constraint index outside the row data" ) );

   sub_rhs.push_back( 1.0 );
   }
  }
 else {
  // inequality constraints: distinguish among the several senses
  for( auto idx : global_idxs_cons ) {
   if( static_cast< std::size_t >( idx ) >= rhs.size() )
    throw( std::runtime_error( "PIPSMILPSolver::extract_rhs_active_flag: "
			       "constraint index outside the row data" ) );

   if( sense[ idx ] == 'G' )
    // >= LHS, set false
    sub_rhs.push_back( 0.0 );
   else
    // in all the other cases set true
    sub_rhs.push_back( 1.0 );
   }
  }

 if( static_cast< int >( sub_rhs.size() ) != len )
  throw( std::runtime_error(
	  "PIPSMILPSolver::extract_rhs_active_flag: node " +
	  std::to_string( id ) + " built " +
	  std::to_string( sub_rhs.size() ) +
	  " values but PIPS requested length " + std::to_string( len ) ) );

 std::copy( sub_rhs.begin() , sub_rhs.end() , vec );

 return( 0 );
}

/*--------------------------------------------------------------------------*/

int PIPSMILPSolver::extract_lhs_active_flag( int id , double * vec , int len ,
	      const std::vector< const FRowConstraint * > & node_cons ,
	      const std::vector< double > & rhs ,
	      const std::vector< char > & sense ,
	      const std::vector< double > & ranges )
{
 if( ( id < 0 ) || ( id >= static_cast< int >( n_nodes ) ) )
  throw( std::runtime_error( "PIPSMILPSolver::extract_lhs_active_flag: node "
			     "ID outside the expected range" ) );

 if( node_cons.empty() ) {
  // zero-fill
  for( int i = 0 ; i < len ; ++i )
   vec[ i ] = 0.0;

  return( 0 );
  }

 // extract the global indices of the constraints
 std::vector< int > global_idxs_cons = compute_cons_global_idxs( node_cons );

 // extract the corresponding elements and store them
 std::vector< double > sub_rhs;
 sub_rhs.reserve( global_idxs_cons.size() );

 for( auto idx : global_idxs_cons ) {
  if( static_cast< std::size_t >( idx ) >= rhs.size() )
   throw( std::runtime_error( "PIPSMILPSolver::extract_lhs_active_flag: "
			      "constraint index outside the row data" ) );

  if( sense[ idx ] == 'L' )
   // <= RHS, set false
   sub_rhs.push_back( 0.0 );
  else
   // in all the other cases set true
   sub_rhs.push_back( 1.0 );
  }

 if( static_cast< int >( sub_rhs.size() ) != len )
  throw( std::runtime_error(
	  "PIPSMILPSolver::extract_lhs_active_flag: node " +
	  std::to_string( id ) + " built " +
	  std::to_string( sub_rhs.size() ) +
	  " values but PIPS requested length " + std::to_string( len ) ) );

 std::copy( sub_rhs.begin() , sub_rhs.end() , vec );

 return( 0 );
}

/*--------------------------------------------------------------------------*/

int PIPSMILPSolver::extract_var_bounds( int id , double * vec , int len ,
	      const std::vector< const ColVariable * > & node_vars ,
	      const std::vector< double > & bounds )
{
 if( ( id < 0 ) || ( id >= static_cast< int >( n_nodes ) ) )
  throw( std::runtime_error( "PIPSMILPSolver::extract_var_bounds: node ID "
			     "outside the expected range" ) );

 if( node_vars.empty() ) {
  // zero-fill
  for( int i = 0 ; i < len ; ++i )
   vec[ i ] = 0.0;

  return( 0 );
  }

 // extract the global indices of the variables
 std::vector< int > global_idxs_vars = compute_vars_global_idxs( node_vars );

 // extract the corresponding elements and store them
 std::vector< double > sub_bounds;
 sub_bounds.reserve( global_idxs_vars.size() );

 for( auto idx : global_idxs_vars ) {
  if( static_cast< std::size_t >( idx ) >= bounds.size() )
   throw( std::runtime_error( "PIPSMILPSolver::extract_var_bounds: variable "
			      "index outside the column data" ) );

  if( ( bounds[ idx ] == -Inf< double >() ) ||
      ( bounds[ idx ] == Inf< double >() ) )
   sub_bounds.push_back( 0.0 );
  else
   sub_bounds.push_back( bounds[ idx ] );
  }

 if( static_cast< int >( sub_bounds.size() ) != len )
  throw( std::runtime_error(
	  "PIPSMILPSolver::extract_var_bounds: node " +
	  std::to_string( id ) + " built " +
	  std::to_string( sub_bounds.size() ) +
	  " values but PIPS requested length " + std::to_string( len ) ) );

 std::copy( sub_bounds.begin() , sub_bounds.end() , vec );

 return( 0 );
}

/*--------------------------------------------------------------------------*/

int PIPSMILPSolver::extract_flag_var_bounds( int id , double * vec , int len ,
	      const std::vector< const ColVariable * > & node_vars ,
	      const std::vector< double > & bounds )
{
 if( ( id < 0 ) || ( id >= static_cast< int >( n_nodes ) ) )
  throw( std::runtime_error( "PIPSMILPSolver::extract_flag_var_bounds: node "
			     "ID outside the expected range" ) );

 if( node_vars.empty() ) {
  // zero-fill
  for( int i = 0 ; i < len ; ++i )
   vec[ i ] = 0.0;

  return( 0 );
  }

 // extract the global indices of the variables
 std::vector< int > global_idxs_vars = compute_vars_global_idxs( node_vars );

 // extract the corresponding elements and store them
 std::vector< double > sub_bounds;
 sub_bounds.reserve( global_idxs_vars.size() );

 for( auto idx : global_idxs_vars ) {
  if( static_cast< std::size_t >( idx ) >= bounds.size() )
   throw( std::runtime_error( "PIPSMILPSolver::extract_flag_var_bounds: "
			      "variable index outside the column data" ) );

  if( ( bounds[ idx ] == -Inf< double >() ) ||
      ( bounds[ idx ] == Inf< double >() ) )
   sub_bounds.push_back( 0.0 );
  else
   sub_bounds.push_back( 1.0 );
  }

 if( static_cast< int >( sub_bounds.size() ) != len )
  throw( std::runtime_error(
	  "PIPSMILPSolver::extract_flag_var_bounds: node " +
	  std::to_string( id ) + " built " +
	  std::to_string( sub_bounds.size() ) +
	  " values but PIPS requested length " + std::to_string( len ) ) );

 std::copy( sub_bounds.begin() , sub_bounds.end() , vec );

 return( 0 );
}

/*--------------------------------------------------------------------------*/

int PIPSMILPSolver::extract_obj( int id , double * vec , int len ,
	      const std::vector< const ColVariable * > & node_vars ,
	      const std::vector< double > & obj_value , int objsense )
{
 if( ( id < 0 ) || ( id >= static_cast< int >( n_nodes ) ) )
  throw( std::runtime_error( "PIPSMILPSolver::extract_obj: node ID outside "
			     "the expected range" ) );

 if( node_vars.empty() ) {
  // zero-fill
  for( int i = 0 ; i < len ; ++i )
   vec[ i ] = 0.0;

  return( 0 );
  }

 // extract the global indices of the variables
 std::vector< int > global_idxs_vars = compute_vars_global_idxs( node_vars );

 // extract the corresponding elements and store them
 std::vector< double > sub_obj;
 sub_obj.reserve( global_idxs_vars.size() );

 for( auto idx : global_idxs_vars ) {
  if( static_cast< std::size_t >( idx ) >= obj_value.size() )
   throw( std::runtime_error( "PIPSMILPSolver::extract_obj: variable index "
			      "outside the objective data" ) );

  // PIPS always assumes to work with a minimization model, so the
  // MILPSolver objective values must be appropriately converted
  sub_obj.push_back( obj_value[ idx ] * objsense );
  }

 if( static_cast< int >( sub_obj.size() ) != len )
  throw( std::runtime_error(
	  "PIPSMILPSolver::extract_obj: node " +
	  std::to_string( id ) + " built " +
	  std::to_string( sub_obj.size() ) +
	  " values but PIPS requested length " + std::to_string( len ) ) );

 std::copy( sub_obj.begin() , sub_obj.end() , vec );

 return( 0 );
}

/*--------------------------------------------------------------------------*/

int PIPSMILPSolver::no_var_in_node( void * user_data , int id , int * nnz )
{
 // use the current class in the callback
 auto * solver = static_cast< PIPSMILPSolver * >( user_data );

 if( id >= static_cast< int >( solver->n_nodes ) )
  throw( std::runtime_error( "PIPSMILPSolver::no_var_in_node: tried to "
			     "query number of variables outside the current "
			     "node count" ) );

 *nnz = solver->n_var_node[ id ];

 return( 0 );
}

/*--------------------------------------------------------------------------*/

int PIPSMILPSolver::no_eq_cons_in_node( void * user_data , int id , int * nnz )
{
 // use the current class in the callback
 auto * solver = static_cast< PIPSMILPSolver * >( user_data );

 if( id >= static_cast< int >( solver->n_nodes ) )
  throw( std::runtime_error( "PIPSMILPSolver::no_eq_cons_in_node: tried to "
			     "query number of equality constraints outside "
			     "the current node count" ) );

 *nnz = solver->n_eq_cons_node[ id ];

 return( 0 );
}

/*--------------------------------------------------------------------------*/

int PIPSMILPSolver::no_ineq_cons_in_node( void * user_data , int id ,
					  int * nnz )
{
 // use the current class in the callback
 auto * solver = static_cast< PIPSMILPSolver * >( user_data );

 if( id >= static_cast< int >( solver->n_nodes ) )
  throw( std::runtime_error( "PIPSMILPSolver::no_ineq_cons_in_node: tried "
			     "to query number of inequality constraints "
			     "outside the current node count" ) );

 *nnz = solver->n_ineq_cons_node[ id ];

 return( 0 );
}

/*--------------------------------------------------------------------------*/

int PIPSMILPSolver::no_link_eq_cons( void * user_data , int id , int * nnz )
{
 // use the current class in the callback
 auto * solver = static_cast< PIPSMILPSolver * >( user_data );

 *nnz = solver->n_link_eq_cons;

 return( 0 );
}

/*--------------------------------------------------------------------------*/

int PIPSMILPSolver::no_link_ineq_cons( void * user_data , int id , int * nnz )
{
 // use the current class in the callback
 auto * solver = static_cast< PIPSMILPSolver * >( user_data );

 *nnz = solver->n_link_ineq_cons;

 return( 0 );
}

/*--------------------------------------------------------------------------*/

int PIPSMILPSolver::nnz_eq_cons_diag( void * user_data , int id , int * nnz )
{
 // use the current class in the callback
 auto * solver = static_cast< PIPSMILPSolver * >( user_data );

 solver->evaluate_nnz( id , nnz , solver->eq_cons_node[ id ] ,
		       solver->var_node[ id ] );

 return( 0 );
}

/*--------------------------------------------------------------------------*/

int PIPSMILPSolver::nnz_eq_cons_vert( void * user_data , int id , int * nnz )
{
 // extract the sub-matrix for a particular node, but considering the root
 // variables (i.e. T_1, T_2, ..., T_n)
 if( id == 0 ) {
  // T_0 can already be obtained through nnz_eq_cons_diag( node_id = 0 )
  *nnz = 0;
  return( 0 );
  }

 // use the current class in the callback
 auto * solver = static_cast< PIPSMILPSolver * >( user_data );

 solver->evaluate_nnz( id , nnz , solver->eq_cons_node[ id ] ,
		       solver->var_node[ 0 ] );

 return( 0 );
}

/*--------------------------------------------------------------------------*/

int PIPSMILPSolver::nnz_ineq_cons_diag( void * user_data , int id , int * nnz )
{
 // use the current class in the callback
 auto * solver = static_cast< PIPSMILPSolver * >( user_data );

 solver->evaluate_nnz( id , nnz , solver->ineq_cons_node[ id ] ,
		       solver->var_node[ id ] );

 return( 0 );
}

/*--------------------------------------------------------------------------*/

int PIPSMILPSolver::nnz_ineq_cons_vert( void * user_data , int id , int * nnz )
{
 // extract the sub-matrix for a particular node, but considering the root
 // variables (i.e. T_1, T_2, ..., T_n)
 if( id == 0 ) {
  // T_0 can already be obtained through nnz_ineq_cons_diag( node_id = 0 )
  *nnz = 0;
  return( 0 );
  }

 // use the current class in the callback
 auto * solver = static_cast< PIPSMILPSolver * >( user_data );

 solver->evaluate_nnz( id , nnz , solver->ineq_cons_node[ id ] ,
		       solver->var_node[ 0 ] );

 return( 0 );
}

/*--------------------------------------------------------------------------*/

int PIPSMILPSolver::nnz_link_eq_cons( void * user_data , int id , int * nnz )
{
 // use the current class in the callback
 auto * solver = static_cast< PIPSMILPSolver * >( user_data );

 solver->evaluate_nnz( id , nnz , solver->link_eq_cons ,
		       solver->var_node[ id ] );

 return( 0 );
}

/*--------------------------------------------------------------------------*/

int PIPSMILPSolver::nnz_link_ineq_cons( void * user_data , int id , int * nnz )
{
 // use the current class in the callback
 auto * solver = static_cast< PIPSMILPSolver * >( user_data );

 solver->evaluate_nnz( id , nnz , solver->link_ineq_cons ,
		       solver->var_node[ id ] );

 return( 0 );
}

/*--------------------------------------------------------------------------*/

int PIPSMILPSolver::nnz_all_zero( void * , int , int * nnz )
{
 *nnz = 0;

 return( 0 );
}

/*--------------------------------------------------------------------------*/

int PIPSMILPSolver::mat_eq_cons_diag( void * user_data , int id , int * krowM ,
				      int * jcolM , double * M )
{
 // use the current class in the callback
 auto * solver = static_cast< PIPSMILPSolver * >( user_data );

 solver->extract_matrix( id , krowM , jcolM , M ,
			 solver->eq_cons_node[ id ] , solver->var_node[ id ] );

 return( 0 );
}

/*--------------------------------------------------------------------------*/

int PIPSMILPSolver::mat_eq_cons_vert( void * user_data , int id , int * krowM ,
				      int * jcolM , double * M )
{
 // extract the sub-matrix for a particular node, but considering the root
 // variables (i.e. T_1, T_2, ..., T_n)
 if( id == 0 )
  // T_0 can already be obtained through mat_eq_cons_diag( node_id = 0 )
  return( 0 );

 // use the current class in the callback
 auto * solver = static_cast< PIPSMILPSolver * >( user_data );

 solver->extract_matrix( id , krowM , jcolM , M ,
			 solver->eq_cons_node[ id ] , solver->var_node[ 0 ] );

 return( 0 );
}

/*--------------------------------------------------------------------------*/

int PIPSMILPSolver::mat_ineq_cons_diag( void * user_data , int id ,
					int * krowM , int * jcolM , double * M )
{
 // use the current class in the callback
 auto * solver = static_cast< PIPSMILPSolver * >( user_data );

 solver->extract_matrix( id , krowM , jcolM , M ,
			 solver->ineq_cons_node[ id ] ,
			 solver->var_node[ id ] );

 return( 0 );
}

/*--------------------------------------------------------------------------*/

int PIPSMILPSolver::mat_ineq_cons_vert( void * user_data , int id ,
					int * krowM , int * jcolM , double * M )
{
 // extract the sub-matrix for a particular node, but considering the root
 // variables (i.e. T_1, T_2, ..., T_n)
 if( id == 0 )
  // T_0 can already be obtained through mat_ineq_cons_diag( node_id = 0 )
  return( 0 );

 // use the current class in the callback
 auto * solver = static_cast< PIPSMILPSolver * >( user_data );

 solver->extract_matrix( id , krowM , jcolM , M ,
			 solver->ineq_cons_node[ id ] ,
			 solver->var_node[ 0 ] );

 return( 0 );
}

/*--------------------------------------------------------------------------*/

int PIPSMILPSolver::mat_link_eq_cons( void * user_data , int id , int * krowM ,
				      int * jcolM , double * M )
{
 // use the current class in the callback
 auto * solver = static_cast< PIPSMILPSolver * >( user_data );

 solver->extract_matrix( id , krowM , jcolM , M ,
			 solver->link_eq_cons , solver->var_node[ id ] );

 return( 0 );
}

/*--------------------------------------------------------------------------*/

int PIPSMILPSolver::mat_link_ineq_cons( void * user_data , int id ,
					int * krowM , int * jcolM , double * M )
{
 // use the current class in the callback
 auto * solver = static_cast< PIPSMILPSolver * >( user_data );

 solver->extract_matrix( id , krowM , jcolM , M ,
			 solver->link_ineq_cons , solver->var_node[ id ] );

 return( 0 );
}

/*--------------------------------------------------------------------------*/

int PIPSMILPSolver::mat_all_zero( void * , int , int * , int * , double * )
{
 return( 0 );
}

/*--------------------------------------------------------------------------*/

int PIPSMILPSolver::obj_vars( void * user_data , int id , double * vec ,
			      int len )
{
 // use the current class in the callback
 auto * solver = static_cast< PIPSMILPSolver * >( user_data );

 solver->extract_obj( id , vec , len , solver->var_node[ id ] ,
		      solver->objective , solver->objsense );

 return( 0 );
}

/*--------------------------------------------------------------------------*/

int PIPSMILPSolver::rhs_eq_cons( void * user_data , int id , double * vec ,
				 int len )
{
 // use the current class in the callback
 auto * solver = static_cast< PIPSMILPSolver * >( user_data );

 solver->extract_rhs_vector( id , vec , len , solver->eq_cons_node[ id ] ,
			     solver->rhs , {} , {} );

 return( 0 );
}

/*--------------------------------------------------------------------------*/

int PIPSMILPSolver::rhs_ineq_cons( void * user_data , int id , double * vec ,
				   int len )
{
 // use the current class in the callback
 auto * solver = static_cast< PIPSMILPSolver * >( user_data );

 solver->extract_rhs_vector( id , vec , len , solver->ineq_cons_node[ id ] ,
			     solver->rhs , solver->sense , solver->rngval );

 return( 0 );
}

/*--------------------------------------------------------------------------*/

int PIPSMILPSolver::lhs_ineq_cons( void * user_data , int id , double * vec ,
				   int len )
{
 // use the current class in the callback
 auto * solver = static_cast< PIPSMILPSolver * >( user_data );

 solver->extract_lhs_vector( id , vec , len , solver->ineq_cons_node[ id ] ,
			     solver->rhs , solver->sense , solver->rngval );

 return( 0 );
}

/*--------------------------------------------------------------------------*/

int PIPSMILPSolver::rhs_link_eq_cons( void * user_data , int id ,
				      double * vec , int len )
{
 // use the current class in the callback
 auto * solver = static_cast< PIPSMILPSolver * >( user_data );

 solver->extract_rhs_vector( id , vec , len , solver->link_eq_cons ,
			     solver->rhs , {} , {} );

 return( 0 );
}

/*--------------------------------------------------------------------------*/

int PIPSMILPSolver::rhs_link_ineq_cons( void * user_data , int id ,
					double * vec , int len )
{
 // use the current class in the callback
 auto * solver = static_cast< PIPSMILPSolver * >( user_data );

 solver->extract_rhs_vector( id , vec , len , solver->link_ineq_cons ,
			     solver->rhs , solver->sense , solver->rngval );

 return( 0 );
}

/*--------------------------------------------------------------------------*/

int PIPSMILPSolver::lhs_link_ineq_cons( void * user_data , int id ,
					double * vec , int len )
{
 // use the current class in the callback
 auto * solver = static_cast< PIPSMILPSolver * >( user_data );

 solver->extract_lhs_vector( id , vec , len , solver->link_ineq_cons ,
			     solver->rhs , solver->sense , solver->rngval );

 return( 0 );
}

/*--------------------------------------------------------------------------*/

int PIPSMILPSolver::flag_rhs_ineq_cons( void * user_data , int id ,
					double * vec , int len )
{
 // use the current class in the callback
 auto * solver = static_cast< PIPSMILPSolver * >( user_data );

 solver->extract_rhs_active_flag( id , vec , len ,
				  solver->ineq_cons_node[ id ] ,
				  solver->rhs , solver->sense ,
				  solver->rngval );

 return( 0 );
}

/*--------------------------------------------------------------------------*/

int PIPSMILPSolver::flag_lhs_ineq_cons( void * user_data , int id ,
					double * vec , int len )
{
 // use the current class in the callback
 auto * solver = static_cast< PIPSMILPSolver * >( user_data );

 solver->extract_lhs_active_flag( id , vec , len ,
				  solver->ineq_cons_node[ id ] ,
				  solver->rhs , solver->sense ,
				  solver->rngval );

 return( 0 );
}

/*--------------------------------------------------------------------------*/

int PIPSMILPSolver::flag_rhs_link_ineq_cons( void * user_data , int id ,
					     double * vec , int len )
{
 // use the current class in the callback
 auto * solver = static_cast< PIPSMILPSolver * >( user_data );

 solver->extract_rhs_active_flag( id , vec , len , solver->link_ineq_cons ,
				  solver->rhs , solver->sense ,
				  solver->rngval );

 return( 0 );
}

/*--------------------------------------------------------------------------*/

int PIPSMILPSolver::flag_lhs_link_ineq_cons( void * user_data , int id ,
					     double * vec , int len )
{
 // use the current class in the callback
 auto * solver = static_cast< PIPSMILPSolver * >( user_data );

 solver->extract_lhs_active_flag( id , vec , len , solver->link_ineq_cons ,
				  solver->rhs , solver->sense ,
				  solver->rngval );

 return( 0 );
}

/*--------------------------------------------------------------------------*/

int PIPSMILPSolver::ub_vars( void * user_data , int id , double * vec ,
			     int len )
{
 // use the current class in the callback
 auto * solver = static_cast< PIPSMILPSolver * >( user_data );

 solver->extract_var_bounds( id , vec , len , solver->var_node[ id ] ,
			     solver->ub );

 return( 0 );
}

/*--------------------------------------------------------------------------*/

int PIPSMILPSolver::lb_vars( void * user_data , int id , double * vec ,
			     int len )
{
 // use the current class in the callback
 auto * solver = static_cast< PIPSMILPSolver * >( user_data );

 solver->extract_var_bounds( id , vec , len , solver->var_node[ id ] ,
			     solver->lb );

 return( 0 );
}

/*--------------------------------------------------------------------------*/

int PIPSMILPSolver::flag_ub_vars( void * user_data , int id , double * vec ,
				  int len )
{
 // use the current class in the callback
 auto * solver = static_cast< PIPSMILPSolver * >( user_data );

 solver->extract_flag_var_bounds( id , vec , len , solver->var_node[ id ] ,
				  solver->ub );

 return( 0 );
}

/*--------------------------------------------------------------------------*/

int PIPSMILPSolver::flag_lb_vars( void * user_data , int id , double * vec ,
				  int len )
{
 // use the current class in the callback
 auto * solver = static_cast< PIPSMILPSolver * >( user_data );

 solver->extract_flag_var_bounds( id , vec , len , solver->var_node[ id ] ,
				  solver->lb );

 return( 0 );
}

/*--------------------------------------------------------------------------*/
/*--------------------- End File PIPSMILPSolver.cpp ------------------------*/
