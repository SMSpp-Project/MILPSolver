/*--------------------------------------------------------------------------*/
/*------------------------- File PIPSMILPSolver.cpp -------------------------*/
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

#include <queue>

#include <LinearFunction.h>

#include <QuadFunction.h>

#include "PIPSMILPSolver.h"

#include "PIPS_maps.h"

#ifdef MILPSOLVER_DEBUG
 #define DEBUG_LOG( stuff ) std::cout << "[MILPSolver DEBUG] " << stuff
#else
 #define DEBUG_LOG( stuff )
#endif

/*--------------------------------------------------------------------------*/
/*------------------------- NAMESPACE AND USING ----------------------------*/
/*--------------------------------------------------------------------------*/

using namespace SMSpp_di_unipi_it;

/*
 * Anonymous namespace for file-local debug/printing utilities.
 * These helpers are used to inspect callback inputs and reconstructed
 * PIPS data structures while debugging the PIPSMILPSolver interface.
 */

#ifndef PIPS_CALLBACK_SANITY
 #define PIPS_CALLBACK_SANITY 0
#endif

namespace {

struct PIPSCallbackNullStream {
 template< class T >
 const PIPSCallbackNullStream & operator<<( const T & ) const { return *this; }

 const PIPSCallbackNullStream & operator<<( std::ostream & ( * )
        ( std::ostream & ) ) const
 {
  return *this;
 }
};

static const PIPSCallbackNullStream pips_callback_null_stream;

#if PIPS_CALLBACK_SANITY
 #define PIPS_CALLBACK_COUT std::cout
#else
 #define PIPS_CALLBACK_COUT pips_callback_null_stream
#endif

#if PIPS_CALLBACK_SANITY

void sanity_print_indices( const char * what , const std::vector< int > & idxs )
{
 PIPS_CALLBACK_COUT << "  " << what << ":";
 for( int idx : idxs )
  PIPS_CALLBACK_COUT << " " << idx;
 PIPS_CALLBACK_COUT << std::endl;
}

void sanity_check_id( const char * cb , int id , int n_nodes )
{
 PIPS_CALLBACK_COUT << "\n[PIPS CALLBACK] " << cb << " id=" << id << std::endl;
 if( id < 0 || id >= n_nodes )
  PIPS_CALLBACK_COUT << "  [SANITY ERROR] id outside [0," << n_nodes << ")" 
    << std::endl;
}

void sanity_check_count( const char * cb , int id , int value , 
                          int expected_min = 0 )
{
 PIPS_CALLBACK_COUT << "[PIPS CALLBACK] " << cb
           << " id=" << id
           << " returns " << value << std::endl;
 if( value < expected_min )
  PIPS_CALLBACK_COUT << "  [SANITY ERROR] negative/invalid count" 
    << std::endl;
}

template< class CSR >
void sanity_check_csr( const char * cb , int id , int nRows , int nCols ,
                       const CSR & mat ,
                       const std::vector< int > & globalRows ,
                       const std::vector< int > & globalCols )
{
 PIPS_CALLBACK_COUT << "\n[PIPS CALLBACK MATRIX] " << cb
           << " id=" << id
           << " nRows=" << nRows
           << " nCols=" << nCols
           << " nnz=" << mat.val.size() << std::endl;

 sanity_print_indices( "global rows", globalRows );
 sanity_print_indices( "global cols", globalCols );

 bool ok = true;

 if( static_cast< int >( mat.krow.size() ) != nRows + 1 ) {
  PIPS_CALLBACK_COUT << "  [SANITY ERROR] krow.size()=" << mat.krow.size()
            << " but expected " << nRows + 1 << std::endl;
  ok = false;
 }

 if( mat.jcol.size() != mat.val.size() ) {
  PIPS_CALLBACK_COUT << "  [SANITY ERROR] jcol.size()=" << mat.jcol.size()
            << " differs from val.size()=" << mat.val.size() << std::endl;
  ok = false;
 }

 if( ! mat.krow.empty() && mat.krow.front() != 0 ) {
  PIPS_CALLBACK_COUT << "  [SANITY ERROR] krow[0]=" << mat.krow.front()
            << " but expected 0" << std::endl;
  ok = false;
 }

 if( ! mat.krow.empty() &&
     mat.krow.back() != static_cast< int >( mat.val.size() ) ) {
  PIPS_CALLBACK_COUT << "  [SANITY ERROR] krow.back()=" << mat.krow.back()
            << " but nnz=" << mat.val.size() << std::endl;
  ok = false;
 }

 for( int i = 0 ; i + 1 < static_cast< int >( mat.krow.size() ) ; ++i ) {
  if( mat.krow[ i ] > mat.krow[ i + 1 ] ) {
   PIPS_CALLBACK_COUT << "  [SANITY ERROR] krow is decreasing at row " << i
             << ": " << mat.krow[ i ] << " > " << mat.krow[ i + 1 ]
             << std::endl;
   ok = false;
  }
 }

 std::vector< int > col_nnz( std::max( nCols , 0 ) , 0 );

 for( int i = 0 ; i < nRows && i + 1 < static_cast< int >( mat.krow.size() ) ; ++i ) {
  const int row_start = mat.krow[ i ];
  const int row_end = mat.krow[ i + 1 ];

  if( row_start < 0 || row_end < 0 ||
      row_start > static_cast< int >( mat.val.size() ) ||
      row_end > static_cast< int >( mat.val.size() ) ) {
   PIPS_CALLBACK_COUT << "  [SANITY ERROR] row " << i
             << " has invalid pointer range [" << row_start
             << "," << row_end << ")" << std::endl;
   ok = false;
   continue;
  }

  for( int k = row_start ; k < row_end ; ++k ) {
   if( mat.jcol[ k ] < 0 || mat.jcol[ k ] >= nCols ) {
    PIPS_CALLBACK_COUT << "  [SANITY ERROR] entry k=" << k
              << " has invalid local column " << mat.jcol[ k ]
              << " outside [0," << nCols << ")" << std::endl;
    ok = false;
   }
   else if( mat.jcol[ k ] < static_cast< int >( col_nnz.size() ) )
    ++col_nnz[ mat.jcol[ k ] ];

   if( ! std::isfinite( mat.val[ k ] ) ) {
    PIPS_CALLBACK_COUT << "  [SANITY ERROR] entry k=" << k
              << " is not finite: " << mat.val[ k ] << std::endl;
    ok = false;
   }
  }
 }

 for( int i = 0 ; i < nRows && i + 1 < static_cast< int >( mat.krow.size() ) ; ++i )
  PIPS_CALLBACK_COUT << "  row " << i << " nnz=" << mat.krow[ i + 1 ] - mat.krow[ i ]
            << std::endl;

 for( int j = 0 ; j < static_cast< int >( col_nnz.size() ) ; ++j )
  PIPS_CALLBACK_COUT << "  col " << j << " nnz=" << col_nnz[ j ] << std::endl;

 for( int i = 0 ; i < nRows && i + 1 < static_cast< int >( mat.krow.size() ) ; ++i ) {
  for( int k = mat.krow[ i ] ; k < mat.krow[ i + 1 ] ; ++k ) {
   PIPS_CALLBACK_COUT << "    localRow=" << i
             << " localCol=" << mat.jcol[ k ]
             << " val=" << mat.val[ k ] << std::endl;
  }
 }

 if( ok )
  PIPS_CALLBACK_COUT << "  [SANITY OK] CSR structure is consistent" << std::endl;
}

void sanity_check_vector( const char * cb , int id , const double * vec ,
                          int len , int expected_len )
{
 PIPS_CALLBACK_COUT << "\n[PIPS CALLBACK VECTOR] " << cb
           << " id=" << id
           << " len=" << len
           << " expected_len=" << expected_len << std::endl;

 if( len != expected_len )
  PIPS_CALLBACK_COUT << "  [SANITY WARNING] PIPS requested len=" << len
            << " but local expected length is " << expected_len << std::endl;

 for( int i = 0 ; i < len ; ++i ) {
  PIPS_CALLBACK_COUT << "  [" << i << "]=" << vec[ i ] << std::endl;
  if( std::isnan( vec[ i ] ) )
   PIPS_CALLBACK_COUT << "  [SANITY ERROR] NaN at position " << i << std::endl;
 }
}

#else

void sanity_check_id( const char * , int , int ) {}
void sanity_check_count( const char * , int , int , int = 0 ) {}
template< class CSR >
void sanity_check_csr( const char * , int , int , int ,
                       const CSR & ,
                       const std::vector< int > & ,
                       const std::vector< int > & ) {}
void sanity_check_vector( const char * , int , const double * , int , int ) {}

#endif

} // end anonymous namespace

/*--------------------------------------------------------------------------*/
/*-------------------------- FACTORY MANAGEMENT ----------------------------*/
/*--------------------------------------------------------------------------*/

SMSpp_insert_in_factory_cpp_0( PIPSMILPSolver );

/*--------------------------------------------------------------------------*/
/*--------------------- CONSTRUCTOR AND DESTRUCTOR -------------------------*/
/*--------------------------------------------------------------------------*/

PIPSMILPSolver::PIPSMILPSolver( void ) :
 MILPSolver() ,
 UpCutOff( Inf< double >() ) ,
 LwCutOff( -Inf< double >() ) ,

 /* Pips problem handling structures */
 pips_tree( nullptr ) ,
 pips_interface( nullptr ) { }

/*--------------------------------------------------------------------------*/

PIPSMILPSolver::~PIPSMILPSolver()
{
 for( auto el : v_ConfigDB )
  delete el;

 // Remove PIPS problem handling structures
 delete pips_interface;
 pips_interface = nullptr;

 delete pips_tree;
 pips_tree = nullptr;

 int finalized = 0;
 MPI_Finalized(&finalized);

 if (mpi_initialized_by_this_solver && !finalized) {
  MPI_Finalize();
 }

}

/*--------------------------------------------------------------------------*/
/*--------------------- DERIVED METHODS OF BASE CLASS ----------------------*/
/*--------------------------------------------------------------------------*/

void PIPSMILPSolver::set_Block( Block * block )
{
 if( block == f_Block )
  return;

 int initialized = 0;
 MPI_Initialized( &initialized );

 if( ! initialized ) {
  int err = MPI_Init( nullptr , nullptr );

  int initialized = 0;
  MPI_Initialized( &initialized );

  if( err != MPI_SUCCESS || ! initialized ) {
      throw std::runtime_error("MPI_Init failed");
  }
  MPIErrorHandler::init();  // Important: set error handlr
  mpi_initialized_by_this_solver = true;
 }

 MILPSolver::set_Block( block );

 UpCutOff = Inf< double >();
 LwCutOff = -Inf< double >();
}

/*--------------------------------------------------------------------------*/

void PIPSMILPSolver::load_problem( void )
{
 // Call MILPSolver to generate the entire LP problem
 MILPSolver::load_problem();

 if( numquadrows > 0 )
  throw( std::runtime_error( "PIPS cannot solve QCP models" ) );

 // Clear PIPS structures
 int status = 0;
  
 delete pips_interface;
 pips_interface = nullptr;

 delete pips_tree;
 pips_tree = nullptr;

 var_to_node.clear();

 /* The strategy in PIPSMILPSolver will simply be to correctly store the
  * elements belonging to each node. After that, by simply calling 
  * extractNode_sub_matrix, we will be able to extract the portion of
  * matrix for the set of constraints and variables that are involved in the
  * node, and pass it to PIPS. */

 // construct the set of involved Block - - - - - - - - - - - - - - - - - - -
 // At the moment, we will store the Blocks in the following way:
 // - The f_Block will be the root of the PIPS model;
 // - then, there will be one leaf for each sub-block of the f_Block. All the 
 //   sub-sub-...-blocks will be stored in the corresponding leaf.

 /*nodes_subtrees.clear();

 // Set the first node (root) to be the f_Block
 n_nodes = 1;
 Index n_blocks = 1; 
 nodes_subtrees.push_back( { f_Block } );

 // Collect all the Blocks representing the leaves of the PIPS tree
 for( auto leaf : f_Block->get_nested_Blocks() ){
  // Collect the subtree data
  std::vector< Block * > subtree;

  // recursively collect all the sub-tree;
  n_blocks += collect_subtree( leaf , subtree , n_nodes );
  subtree.shrink_to_fit();

  // Add to the global structure
  nodes_subtrees.push_back( std::move( subtree ) );
  n_nodes++;
 }
 nodes_subtrees.shrink_to_fit();*/

 // TEMP
 nodes_subtrees.clear();
 blockToleaf.clear();

 auto nested_blocks = f_Block->get_nested_Blocks();

 if( nested_blocks.empty() )
  throw( std::runtime_error( "PIPSMILPSolver requires at least one nested Block "
                            "to promote as PIPS root." ) );

 // Choose which child becomes the PIPS root.
 // For now: promote the last nested Block (DesignNetworkBlock).
 Block * promoted_root = nested_blocks.back();

 // Set the first node (root) to be the promoted block
 n_nodes = 1;
 Index n_blocks = 1; 
 nodes_subtrees.push_back( { promoted_root } );

 // All the sub_blocks of the promoted root become leaves
 for( auto leaf : promoted_root->get_nested_Blocks(); ) {
  std::vector< Block * > subtree;

  n_blocks += collect_subtree( leaf , subtree , n_nodes );

  subtree.shrink_to_fit();
  nodes_subtrees.push_back( std::move( subtree ) );

  n_nodes++;
 }

 // Other children of f_Block become PIPS leaves
 for( auto leaf : nested_blocks ) {
  if( leaf == promoted_root )
   continue;

  std::vector< Block * > subtree;

  n_blocks += collect_subtree( leaf , subtree , n_nodes );

  subtree.shrink_to_fit();
  nodes_subtrees.push_back( std::move( subtree ) );

  n_nodes++;
 }

 nodes_subtrees.shrink_to_fit();
 // END TEMP

 if( n_nodes == 1)
  throw( std::runtime_error( "PIPSMILPSolver requires at least two blocks "
    "in order to work." ) );

 DEBUG_LOG( "Number of blocks      = " << n_blocks << std::endl );
 DEBUG_LOG( "Number of nodes      = " << n_nodes << std::endl );

 // PIPMILPSolver vectors allocation - - - - - - - - - - - -  - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 n_varNode.resize( n_nodes, 0 );
 n_EqConsNode.resize( n_nodes, 0 );
 n_InEqConsNode.resize( n_nodes, 0 );

 varNode.resize( n_nodes );
 EqConsNode.resize( n_nodes );
 InEqConsNode.resize( n_nodes );

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
     auto push_var_toNode = [ this , num_node ]
      ( const ColVariable & c ) {
      n_varNode[ num_node ] += 1;
      varNode[ num_node ].push_back( &c );
      var_to_node.emplace( &c , num_node );
     };

     un_any_const_dynamic( i , push_var_toNode , 
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

 int rank;
 int size;
 MPI_Comm_rank(MPI_COMM_WORLD, &rank);
 MPI_Comm_size(MPI_COMM_WORLD, &size);

 // Now set all the callbacks - - - - - - - - - - - - - - - - - - - - - - - - 
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 // Number of elements
 FNNZ fNo_VarinNode = &No_VarinNode; // Number of variables
 FNNZ fNo_EqConsinNode = &No_EqConsinNode; 
              // Number of equality constraints per node
 FNNZ fNo_InEqConsinNode = &No_InEqConsinNode; 
              // Number of inequality constraints per node
 FNNZ fNo_LinkEqCons = &No_LinkEqCons; 
              // Number of global linking equality constraints
 FNNZ fNo_LinkInEqCons = &No_LinkInEqCons; 
              // Number of global linking inequality constraints

 // Nonzeros
 FNNZ fnnzQ = &nnzAllZero; 
              // TBD: number of quadratic nonzero terms in objective
 
 FNNZ fnnzEqConsDiag = &nnzEqConsDiag; 
              // Number of nonzeros in the diagonal matrices
 FNNZ fnnzEqConsVert = &nnzEqConsVert; 
              // Number of nonzeros in the vertical root matrices
 FNNZ fnnzInEqConsDiag = &nnzInEqConsDiag; 
              // Number of nonzeros in the equality diagonal matrices
 FNNZ fnnzInEqConsVert = &nnzInEqConsVert; 
              // Number of nonzeros in the inequality vertical root matrices
 FNNZ fnnzLinkEqCons = &nnzLinkEqCons; 
              // Number of nonzeros in the equality linking matrices
 FNNZ fnnzLinkInEqCons = &nnzLinkInEqCons;
              // Number of nonzeros in the inequality linking matrices
   
 // Vectors
 FVEC fRhsEqCons = &RhsEqCons; // rhs of equality constraints per node
 FVEC fRhsInEqCons = &RhsInEqCons; 
              // rhs of linking inequality constraints per node
 FVEC fLhsInEqCons = &LhsInEqCons; 
              // lhs of linking inequality constraints per node
 FVEC fRhsLinkEqCons = &RhsLinkEqCons; 
              // rhs of linking equality constraints
 FVEC fRhsLinkInEqCons = &RhsLinkInEqCons; 
              // rhs of linking inequality constraints
 FVEC fLhsLinkInEqCons = &LhsLinkInEqCons; 
              // lhs of linking inequality constraints

 FVEC fFlagRhsInEqCons = &FlagRhsInEqCons; 
              // active rhs of linking inequality constraints per node
 FVEC fFlagLhsInEqCons = &FlagLhsInEqCons; 
              // active lhs of linking inequality constraints per node
 FVEC fFlagRhsLinkInEqCons = &FlagRhsLinkInEqCons; 
              // active rhs of linking inequality constraints
 FVEC fFlagLhsLinkInEqCons = &FlagLhsLinkInEqCons; 
              // active lhs of linking inequality constraints

 FVEC fUBVars = &UBVars; // Upper bounds on the variables
 FVEC fLBVars = &LBVars; // Lower bounds on the variables
 FVEC fFlagUBVars = &FlagUBVars; // Upper bound flags on the variables
 FVEC fFlagLBVars = &FlagLBVars; // Lower bound flag on the variables

 FVEC fObjVars = &ObjVars; // Objective value

 // Matrices
 FMAT fMatEqConsDiag = &MatEqConsDiag; // Diagonal equality matrices
 FMAT fMatEqConsVert = &MatEqConsVert; // Vertical equality matrices
 FMAT fMatInEqConsDiag = &MatInEqConsDiag; // Diagonal inequality matrices
 FMAT fMatInEqConsVert = &MatInEqConsVert; // Vertical inequality matrices
 FMAT fMatLinkEqCons = &MatLinkEqCons; // Equality linking matrices
 FMAT fMatLinkInEqCons = &MatLinkInEqCons; // Inequality linking matrices

 FMAT fQ = &matAllZero; // TBD Quadratic objective matrix

 // Build the problem tree

 // Create the root node (depth = 0)
 std::unique_ptr<DistributedInputTree::DistributedInputNode> data_root = 
  std::make_unique<DistributedInputTree::DistributedInputNode>( this, 0, 
    fNo_VarinNode, fNo_EqConsinNode, fNo_LinkEqCons, fNo_InEqConsinNode, 
    fNo_LinkInEqCons, fQ, fnnzQ, fObjVars, fMatEqConsVert , fnnzEqConsVert , 
    fMatEqConsDiag, fnnzEqConsDiag, fMatLinkEqCons , fnnzLinkEqCons, 
    fRhsEqCons , fRhsLinkEqCons , fMatInEqConsVert , fnnzInEqConsVert , 
    fMatInEqConsDiag , fnnzInEqConsDiag , fMatLinkInEqCons , 
    fnnzLinkInEqCons , fLhsInEqCons , fFlagLhsInEqCons , fRhsInEqCons , 
    fFlagRhsInEqCons , fLhsLinkInEqCons , fFlagLhsLinkInEqCons , 
    fRhsLinkInEqCons ,fFlagRhsLinkInEqCons , fLBVars , fFlagLBVars , 
    fUBVars , fFlagUBVars , nullptr, nullptr, false );

 auto* root = new DistributedInputTree( std::move( data_root ) );

 for(int id = 1; id < n_nodes ; ++id ) {
  // Create the leaves (depth = id)
  std::unique_ptr<DistributedInputTree::DistributedInputNode> data_child = 
    std::make_unique<DistributedInputTree::DistributedInputNode>( this, id, 
      fNo_VarinNode, fNo_EqConsinNode, fNo_LinkEqCons, fNo_InEqConsinNode, 
      fNo_LinkInEqCons, fQ, fnnzQ, fObjVars, fMatEqConsVert , fnnzEqConsVert , 
      fMatEqConsDiag, fnnzEqConsDiag, fMatLinkEqCons , fnnzLinkEqCons, 
      fRhsEqCons , fRhsLinkEqCons , fMatInEqConsVert , fnnzInEqConsVert , 
      fMatInEqConsDiag , fnnzInEqConsDiag , fMatLinkInEqCons , 
      fnnzLinkInEqCons , fLhsInEqCons , fFlagLhsInEqCons , fRhsInEqCons , 
      fFlagRhsInEqCons , fLhsLinkInEqCons , fFlagLhsLinkInEqCons , 
      fRhsLinkInEqCons ,fFlagRhsLinkInEqCons , fLBVars , fFlagLBVars , 
      fUBVars , fFlagUBVars , nullptr, nullptr, false );

   // Add the child to the root
   root->add_child( 
    std::make_unique<DistributedInputTree>( std::move( data_child ) ) );
 }

 // Store the constructed elements in the SMS++ structure
 pips_tree = root;
 pips_interface = new PIPSIPMppInterface( pips_tree, MPI_COMM_WORLD );

 return;
}

/*--------------------------------------------------------------------------*/

Solver::OFValue PIPSMILPSolver::get_lb( void )
{
 OFValue lower_bound = constant_value;

 // PIPS only admits minimization problems, so we need to convert the result
 // based on the expected MILPSolver objective sense.
 if( objsense == 1 ){ // Minimization
  switch( sol_status ) {
   case( kUnbounded ):  lower_bound = -Inf< OFValue >(); break;
   case( kInfeasible ): lower_bound = Inf< OFValue >();  break;
   default:             lower_bound += pips_interface->getObjective();
    // TBD: Here we should retrieve a bound
  }
 }
 else{
  switch( sol_status ) {
   case( kUnbounded ):  lower_bound = Inf< OFValue >();  break;
   case( kInfeasible ): lower_bound = -Inf< OFValue >(); break;
   default:             lower_bound -= pips_interface->getObjective();
  }
 }

 return( lower_bound );
 }

/*--------------------------------------------------------------------------*/

Solver::OFValue PIPSMILPSolver::get_ub( void )
{
 OFValue upper_bound = constant_value;

 // PIPS only admits minimization problems, so we need to convert the result
 // based on the expected MILPSolver objective sense.
 if( objsense == 1 ){ // Minimization
  switch( sol_status ) {
   case( kUnbounded ):  upper_bound = -Inf< OFValue >(); break;
   case( kInfeasible ): upper_bound = Inf< OFValue >();  break;
   default:             upper_bound += pips_interface->getObjective();
  }
 }
 else{
  switch( sol_status ) {
   case( kUnbounded ):  upper_bound = Inf< OFValue >();  break;
   case( kInfeasible ): upper_bound = -Inf< OFValue >(); break;
   default:             upper_bound -= pips_interface->getObjective();
   // TBD: Here we should retrieve a bound
  }
 }

 return( upper_bound );
}

/*--------------------------------------------------------------------------*/

void PIPSMILPSolver::clear_problem( unsigned int what )
{
 MILPSolver::clear_problem( 0 );

 delete pips_interface;
 pips_interface = nullptr;

 delete pips_tree;
 pips_tree = nullptr;
}

/*--------------------------------------------------------------------------*/

int PIPSMILPSolver::guts_of_compute( void )
{
 // Note: locking, process_modifications() and the LP cut separation loop
 // (when intRelaxIntVars == 2) are all handled by MILPSolver::compute().
 // This method is only responsible for the actual PIPS-IPM++ call.

 // if required, write the problem to file- - - - - - - - - - - - - - - - - -
 // This should have already been done in set_par

 // the continuous case - - - - - - - - - - - - - - - - - - - - - - - - - - -
 sol_status = decode_pips_status( pips_interface->run() );

 Return_status:
 unlock();  // unlock the mutex
 return( sol_status );

}  

/*--------------------------------------------------------------------------*/

int PIPSMILPSolver::decode_pips_status( TerminationStatus status )
{
 DEBUG_LOG( "pips_interface.run() returned " << status << std::endl );

 /* The following are the symbols that may represent the status of
 * a PIPS solution as returned by pips_interface.run(). */

 switch( status ) {
  case( TerminationStatus::READ_ERROR ):
  case( TerminationStatus::UNKNOWN ) :
  case( TerminationStatus::DID_NOT_RUN ):
  case( TerminationStatus::NOT_FINISHED ):
  case( TerminationStatus::STOPPED_AFTER_PRESOLVE ):
   // Some error happened.
   return( kError );
  case( TerminationStatus::TIMELIMIT ):
   // Time limit exceeded
   return( kStopTime );
  case( TerminationStatus::INFEASIBLE ):
   // Problem is infeasible.
   return( kInfeasible );
  case( TerminationStatus::UNBOUNDED ):
   // Problem has been proven unbounded.
   return( kUnbounded );
  case( TerminationStatus::MAX_ITS_EXCEEDED ):
   // Iteration limit has been reached;
   return( kStopIter );
  case( TerminationStatus::SUCCESSFUL_TERMINATION ):
   // Compilation terminated succesfully
   return( kOK );
  default:;
  }

 throw( std::runtime_error( "pips_interface.run() returned unknown status." ) );
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
  // Firstly, if we are setting the output verbosity, we have to be careful 
  // that the SILENT parameter in PIPS has an opposite behaviour with respect
  // to the other solvers. So we must switch the value
  if( Pp == "SILENT" ){
    pipsipmpp_options::set_parameter( Pp , value == 0 );
    return;
  }

  // Now, we can set normally the other parameters but considering that the 
  // option could be a boolean one in PIPS.
  if( value == 0 || value == 1 ) {
   try {
    pipsipmpp_options::set_parameter( Pp , value );
    return;
   }
   catch( const std::runtime_error & ) {
    pipsipmpp_options::set_parameter( Pp , value != 0 );
    return;
   }
  }

  pipsipmpp_options::set_parameter( Pp , value );
  return;
  }

 MILPSolver::set_par( par, value );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

void PIPSMILPSolver::set_par( idx_type par , double value )
{
 // Solver parameters explicitly mapped in PIPS
 switch( par ) {
  case( dblUpCutOff ): UpCutOff = value; return;
  case( dblLwCutOff ): LwCutOff = value; return;
  }

 std::string Pp;
 Pp = pips_dbl_par_map( par );

 if( Pp.size() > 0 ) {
  pipsipmpp_options::set_parameter( Pp , value );
  return;
  }

 MILPSolver::set_par( par , value );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

void PIPSMILPSolver::set_par( idx_type par , std::string && value )
{
 // set the solver log to a specific file
 if( par == strLogFileName ) {
  // No currect option in PIPS
  return;
 }

 // PIPS parameters
 if( ( par >= strFirstPIPSPar ) && ( par < strLastAlgParPIPS ) ) {
  std::string pips_par = SMSpp_to_PIPS_str_pars[ par - strFirstPIPSPar ];
  pipsipmpp_options::set_parameter( pips_par , value );
  return;
  }

 if( par == strOutputFile )
  pipsipmpp_options::set_parameter( "WRITE_ORIGINAL_PROBLEM_TO_LP" , 
                                      value );

 MILPSolver::set_par( par, std::move( value ) );
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
 // intCutSepPar is now handled by MILPSolver base

 std::string Pp = pips_int_par_map( par );
 if( Pp.size() > 0 ) {
   int value = pipsipmpp_options::get_int_parameter( Pp );
   return( value );
  }

 return( MILPSolver::get_dflt_int_par( par ) );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

double PIPSMILPSolver::get_dflt_dbl_par( idx_type par ) const
{
 switch( par ) {
  case( dblUpCutOff ): return( Inf< double >() );
  case( dblLwCutOff ): return( -Inf< double >() );
  }

 std::string Pp = pips_dbl_par_map( par );
 if( Pp.size() > 0 ) {
   double value = pipsipmpp_options::get_double_parameter( Pp );
   return( value );
  }

 return( MILPSolver::get_dflt_dbl_par( par ) );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

const std::string & PIPSMILPSolver::get_dflt_str_par( idx_type par ) const
{
 static std::string value;

 if( ( par >= strFirstPIPSPar ) && ( par < strLastAlgParPIPS ) ) {
  std::string pips_par = SMSpp_to_PIPS_str_pars[ par - strFirstPIPSPar ];
  value.reserve( 512 );
  value = pipsipmpp_options::get_string_parameter( pips_par );

  return( value );
  }

 return( MILPSolver::get_str_par( par ) );
 }

/*--------------------------------------------------------------------------*/

int PIPSMILPSolver::get_int_par( idx_type par ) const
{
 // intCutSepPar is now handled by MILPSolver base

 std::string Pp = pips_int_par_map( par );
 if( Pp.size() > 0 ) {
   int value = pipsipmpp_options::get_int_parameter( Pp );
   return( value );
  }

 return( MILPSolver::get_int_par( par ) );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

double PIPSMILPSolver::get_dbl_par( idx_type par ) const
{
 switch( par ) {
  case( dblUpCutOff ): return( UpCutOff );
  case( dblLwCutOff ): return( LwCutOff );
  }

 std::string Pp = pips_dbl_par_map( par );
 if( Pp.size() > 0 ) {
   double value = pipsipmpp_options::get_double_parameter( Pp );
   return( value );
  }

 return( MILPSolver::get_dbl_par( par ) );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

const std::string & PIPSMILPSolver::get_str_par( idx_type par ) const
{
 static std::string value;

 if( ( par >= strFirstPIPSPar ) && ( par < strLastAlgParPIPS ) ) {
  std::string pips_par = SMSpp_to_PIPS_str_pars[ par - strFirstPIPSPar ];
  value.reserve( 512 );
  value = pipsipmpp_options::get_string_parameter( pips_par );

  return( value );
  }

 return( MILPSolver::get_str_par( par ) );
 }

/*--------------------------------------------------------------------------*/

Solver::idx_type PIPSMILPSolver::int_par_str2idx(
					     const std::string & name ) const
{
 // intCutSepPar is now handled by MILPSolver::int_par_str2idx()

 /* In PIPSMILPSolver::*_par_str2idx() methods we check with MILPSolver first */

 idx_type idx = MILPSolver::int_par_str2idx( name );
 if( idx < Inf< idx_type >() )
  return( idx );

 // PIPS parameters
 std::string pips_par = name;
 auto array_pos = std::find( SMSpp_to_PIPS_int_pars.begin() ,
                        SMSpp_to_PIPS_int_pars.end() ,
                        pips_par);

 if( array_pos != SMSpp_to_PIPS_int_pars.end() ) {
  int pos = std::distance( SMSpp_to_PIPS_int_pars.begin(), array_pos );
  auto idx_par = PIPS_to_SMSpp_int_pars[ pos ].second;
  return( idx_par );
  }

 return( Inf< idx_type >() );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

const std::string & PIPSMILPSolver::int_par_idx2str( idx_type idx ) const
{
 // intCutSepPar is now handled by MILPSolver::int_par_idx2str() (via the
 // fall-through at the end of this method)

 // note: this implementation is not thread safe, and it requires that the
 //       result is used immediately after the call (prior to any other call
 //       to int_par_idx2str()), this may have to be improved upon
 static std::string par_name;
 par_name.reserve( 512 );

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
 /* In PIPSMILPSolver::*_par_str2idx() methods we check with MILPSolver first */

 idx_type idx = MILPSolver::dbl_par_str2idx( name );
 if( idx < Inf< idx_type >() )
  return( idx );

 // PIPS parameters
 std::string pips_par = name;
 auto array_pos = std::find( SMSpp_to_PIPS_dbl_pars.begin() ,
                        SMSpp_to_PIPS_dbl_pars.end() ,
                        pips_par);

 if( array_pos != SMSpp_to_PIPS_dbl_pars.end() ) {
  int pos = std::distance( SMSpp_to_PIPS_dbl_pars.begin(), array_pos );
  auto idx_par = PIPS_to_SMSpp_dbl_pars[ pos ].second;
  return( idx_par );
  }

 return( Inf< idx_type >() );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

const std::string & PIPSMILPSolver::dbl_par_idx2str( idx_type idx ) const
{
 // note: this implementation is not thread safe, and it requires that the
 //       result is used immediately after the call (prior to any other call
 //       to int_par_idx2str()), this may have to be improved upon
 static std::string par_name;
 par_name.reserve( 512 );

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
 /* In PIPSMILPSolver::*_par_str2idx() methods we check with MILPSolver first */

 idx_type idx = MILPSolver::str_par_str2idx( name );
 if( idx < Inf< idx_type >() )
  return( idx );

 // PIPS parameters
 std::string pips_par = name;
 auto array_pos = std::find( SMSpp_to_PIPS_str_pars.begin() ,
                        SMSpp_to_PIPS_str_pars.end() ,
                        pips_par);

 if( array_pos != SMSpp_to_PIPS_str_pars.end() ) {
  int pos = std::distance( SMSpp_to_PIPS_str_pars.begin(), array_pos );
  auto idx_par = PIPS_to_SMSpp_str_pars[ pos ].second;
  return( idx_par );
  }

 return( Inf< idx_type >() );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

const std::string & PIPSMILPSolver::str_par_idx2str( idx_type idx ) const
{
 // note: this implementation is not thread safe, and it requires that the
 //       result is used immediately after the call (prior to any other call
 //       to int_par_idx2str()), this may have to be improved upon
 static std::string par_name;
 par_name.reserve( 512 );

 if( ( idx >= strFirstPIPSPar ) && ( idx < strLastAlgParPIPS ) ) {
  par_name = SMSpp_to_PIPS_str_pars[ idx - strFirstPIPSPar ];
  return( par_name );
  }

 return( MILPSolver::str_par_idx2str( idx ) );
 }

/*--------------------------------------------------------------------------*/
/*----------------------- PIPS CALLBACKS METHODS----------------------------*/
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

PIPSMILPSolver::CSRMatrix PIPSMILPSolver::extractSubmatrixToCRS(
 const std::vector< int > & selectedRows ,
 const int nRows ,
 const std::vector< int > & selectedCols ,
 const int nCols )
{
 CSRMatrix result;

 std::unordered_map< int , int > rowMap;
 rowMap.reserve( nRows );

 for( int i = 0 ; i < nRows ; ++i )
  rowMap[ selectedRows[ i ] ] = i;

 std::vector< std::vector< std::pair< int , double > > > rows( nRows );

 for( int localCol = 0 ; localCol < nCols ; ++localCol ) {
  const int originalCol = selectedCols[ localCol ];

  for( int k = matbeg[ originalCol ] ;
       k < matbeg[ originalCol ] + matcnt[ originalCol ] ;
       ++k ) {

   const int originalRow = matind[ k ];
   auto it = rowMap.find( originalRow );

   if( it != rowMap.end() )
    rows[ it->second ].push_back( { localCol , matval[ k ] } );
  }
 }

 result.krow.reserve( nRows + 1 );
 result.krow.push_back( 0 );

 for( int i = 0 ; i < nRows ; ++i ) {
  for( const auto & el : rows[ i ] ) {
   result.jcol.push_back( el.first );
   result.val.push_back( el.second );
  }

  result.krow.push_back( result.val.size() );
 }

 PIPS_CALLBACK_COUT << "\n[DEBUG extractSubmatrixToCRS]" << std::endl;
 PIPS_CALLBACK_COUT << "nRows = " << nRows
           << ", nCols = " << nCols
           << ", nnz = " << result.val.size()
           << std::endl;

 PIPS_CALLBACK_COUT << "selectedRows: ";
 for( int r : selectedRows )
  PIPS_CALLBACK_COUT << r << " ";
 PIPS_CALLBACK_COUT << std::endl;

 PIPS_CALLBACK_COUT << "selectedCols: ";
 for( int c : selectedCols )
  PIPS_CALLBACK_COUT << c << " ";
 PIPS_CALLBACK_COUT << std::endl;

 for( int i = 0 ; i < nRows ; ++i ) {
  for( int k = result.krow[ i ] ; k < result.krow[ i + 1 ] ; ++k ) {
   PIPS_CALLBACK_COUT << "  localRow = " << i
             << ", localCol = " << result.jcol[ k ]
             << ", val = " << result.val[ k ]
             << std::endl;
  }
 }

 return result;
}

/*--------------------------------------------------------------------------*/

int PIPSMILPSolver::ExtractMatrix( int id , int* krowM, int* jcolM, double* M , 
                              std::vector< const FRowConstraint * > node_cons ,
                              const int nCons , 
                              std::vector< const ColVariable * > vars ,
                              const int nVars ){ 
  if( id < n_nodes ){
   // Check if the node has own constraint or variables
   if( nCons == 0 || nVars == 0 ){
    for( int i = 0 ; i <= nCons ; ++i )
     krowM[ i ] = 0;

    return 0;
   }

   // Extract the global indices of constraints
   std::vector< int > global_idxs_cons = compute_cons_global_idxs( node_cons ,
                                                      nCons );

   // Extract the global indices of variables
   std::vector< int > global_idxs_vars = compute_vars_global_idxs( vars ,
                                                      nVars );
   
   // Extract the corresponding submatrix
   CSRMatrix sub_matrix = extractSubmatrixToCRS( global_idxs_cons ,
                                                      nCons , 
                                                      global_idxs_vars ,
                                                      nVars );
   sanity_check_csr( "ExtractMatrix", id, nCons, nVars, sub_matrix,
                     global_idxs_cons, global_idxs_vars );
   std::copy( sub_matrix.krow.begin(), sub_matrix.krow.end(), krowM );
   std::copy( sub_matrix.jcol.begin(), sub_matrix.jcol.end(), jcolM );
   std::copy( sub_matrix.val.begin(),  sub_matrix.val.end(),  M );
   return 0;
  }
  else
   throw( std::runtime_error( "Node ID outside the expected range" ) );
}

/*--------------------------------------------------------------------------*/

int PIPSMILPSolver::EvaluateNnz( int id , int* nnz , 
                              std::vector< const FRowConstraint * > node_cons ,
                              const int nCons , 
                              std::vector< const ColVariable * > vars ,
                              const int nVars ){ 
  if( id < n_nodes ){
   // Check if the node has own constraint or variables
   if( nCons == 0 || nVars == 0 ){
    *nnz = 0;
    return 0;
   }

   // Extract the global indices of constraints
   std::vector< int > global_idxs_cons = compute_cons_global_idxs( node_cons ,
                                                      nCons );

   // Extract the global indices of variables
   std::vector< int > global_idxs_vars = compute_vars_global_idxs( vars ,
                                                      nVars );
   
   // Extract the corresponding submatrix
   CSRMatrix sub_matrix = extractSubmatrixToCRS( global_idxs_cons ,
                                                      nCons , 
                                                      global_idxs_vars ,
                                                      nVars );
   sanity_check_csr( "EvaluateNnz", id, nCons, nVars, sub_matrix,
                     global_idxs_cons, global_idxs_vars );
   *nnz = sub_matrix.nnz();
   sanity_check_count( "EvaluateNnz", id, *nnz );
   return 0;
  }
  else
   throw( std::runtime_error( "Node ID outside the expected range" ) );
}

/*--------------------------------------------------------------------------*/

std::vector< int > PIPSMILPSolver::compute_cons_global_idxs( 
                                  std::vector< const FRowConstraint * > cons ,
                                  const int nCons ){
 // Initialize vector to -1
 std::vector< int > global_idxs( nCons , -1 );
 Index row = 0;

 for( auto con : cons ){
  Index static_idx = index_of_static_constraint( con );

  if( static_idx < Inf< int >() ){
    global_idxs[ row ] = static_idx;
  }
  else{
    Index dynamic_idx = index_of_dynamic_constraint( con );
    if( dynamic_idx < Inf< int >() )
     global_idxs[ row ] = dynamic_idx;
    else
     throw( std::runtime_error( "The required constraint has no available "
      "index" ) );
  }
  row++;
 }      
 return global_idxs;                            
}

/*--------------------------------------------------------------------------*/

std::vector< int > PIPSMILPSolver::compute_vars_global_idxs( 
                                  std::vector< const ColVariable * > vars ,
                                  const int nVars ){
 // Initialize vector to -1
 std::vector< int > global_idxs( nVars , -1 );
 Index col = 0;

 for( auto var : vars ){
  Index static_idx = index_of_static_variable( var );

  if( static_idx < Inf< int >() ){
    global_idxs[ col ] = static_idx;
  }
  else{
    Index dynamic_idx = index_of_dynamic_variable( var );
    if( dynamic_idx < Inf< int >() )
     global_idxs[ col ] = dynamic_idx;
    else
     throw( std::runtime_error( "The required variable has no available "
      "index" ) );
  }
  col++;
 }      
 return global_idxs;                    
}

/*--------------------------------------------------------------------------*/

int PIPSMILPSolver::collect_subtree(
 Block * block ,
 std::vector< Block * > & subtree ,
 Index parent_leaf
 )
{
 // Set the map from the block to the father leaf
 blockToleaf[ block ] = parent_leaf;

 Index tree_size = 1;
 subtree.push_back( block );

 for( auto * child : block->get_nested_Blocks() ){
  Index subtree_size = collect_subtree( child , subtree , parent_leaf );
  tree_size += subtree_size;
 }
 return tree_size;
}

/*--------------------------------------------------------------------------*/

template< typename T >
 void PIPSMILPSolver::scan_group( const boost::any & gr , Block * qb ,
                              Index num_node , un_any_type< T > )
{
 // Search for the group type
 if( gr.type() == typeid( T * ) ||
      gr.type() == typeid( std::vector< T > * ) ||
      gr.type() == typeid( std::vector< std::vector< T > > * ) ) {
  // "Simple" group
  scan_simple_group( gr , qb , num_node , un_any_type< T >() );
  }
 else {
  // "Complex" group
  scan_multiarray_group( gr , qb , num_node , un_any_type< T >() );
  }
 } 

/*--------------------------------------------------------------------------*/

template< typename T >
 void PIPSMILPSolver::scan_simple_group( const boost::any & gr , Block * qb , 
                                Index num_node , un_any_type< T > )
{
 if( typeid( T * ) == typeid( FRowConstraint * ) ) {
  // Scanning a group of Constraints
  auto scan = [ this , num_node ]
   ( const FRowConstraint & c ) {
    scan_constraint( c , num_node );
  };
  // Scan all the constraints one at a time
  un_any_const_static( gr , scan  , un_any_type< FRowConstraint >() );
 }
 else if( typeid( T * ) == typeid( ColVariable * ) ) {
  // Scanning a group of Variables
  // in this case we simply add the variable to the corresponding node
  auto push_var_toNode = [ this , num_node ]
    ( const ColVariable & c ) {
      n_varNode[ num_node ] += 1;
      varNode[ num_node ].push_back( &c );
      var_to_node.emplace( &c , num_node );
  };
  un_any_const_static( gr , push_var_toNode , un_any_type< ColVariable >() );
  }
 else
  throw( std::runtime_error( "Unsupported group type" ) );
} 

/*--------------------------------------------------------------------------*/

template< typename T >
 void PIPSMILPSolver::scan_multiarray_group( const boost::any & gr ,
            Block * qb , Index num_node , un_any_type< T > )
{
 // Get multi array number of dimension
 int ma_dim = get_multi_array_dim( gr , un_any_type< T >() , 
                                    un_any_int< 2 >() );

 if( ma_dim == 2 ) {
  // Use the 2D multi_array

  // Get type of multi array. See MILPSolver.h:1256 for further details.
  int type = get_multi_array_type( gr , un_any_type< T >() ,
                                  un_any_int< 2 >() );

  // Indices of the 2 dimensions
  int idx_0 = 0;
  int idx_1 = 0;

  if( type == 1 ) {
   // Multi arrays of type 1 (i.e. multi_array< std::vector < T * > >).
   // In this case elements are not stored in sequential cells. Thus, 
   // we have to "unpack" each std::vector and store them separately.
   auto ma = get_multi_array1( gr , un_any_type< T >() ,
                            un_any_int< 2 >() );

   if( typeid( T * ) == typeid( FRowConstraint * ) ) {
    // Constraint group

    // Scan the linearization of the array
    for( auto v = ma->data() ; idx_0 < ma->shape()[ 0 ] ; ++v ) {
     // Scanning a group of Constraints
     auto scan = [ this , num_node ]
      ( const FRowConstraint & c ) {
        scan_constraint( c , num_node );
     };
     // Scan all the constraints one at a time
     un_any_const_static( v , scan  , un_any_type< FRowConstraint >() );

     // The linearization produced by ma->data() for the 2D multi_array
     // stores elements in row-major order. Therefore, we should increment
     // the column index first, and when it exceeds the number of columns,
     // reset it and increment the row index.
     if( idx_1 < ma->shape()[ 1 ] - 1 )
      idx_1++;
     else {
      idx_1 = 0;
      idx_0++;
     }
    }
   }
   else if( typeid( T * ) == typeid( ColVariable * ) ) {
    // Variable group

    // Scan the linearization of the array
    for( auto v = ma->data() ; idx_0 < ma->shape()[ 0 ] ; ++v ) {
     // Scanning a group of Variables
     // in this case we simply add the variable to the corresponding node
     auto push_var_toNode = [ this , num_node ]
       ( const ColVariable & c ) {
        n_varNode[ num_node ] += 1;
        varNode[ num_node ].push_back( &c );
        var_to_node.emplace( &c , num_node );
     };
     un_any_const_static( v , push_var_toNode , 
      un_any_type< ColVariable >() );

     // The linearization produced by ma->data() for the 2D multi_array
     // stores elements in row-major order. Therefore, we should increment
     // the column index first, and when it exceeds the number of columns,
     // reset it and increment the row index.
     if( idx_1 < ma->shape()[ 1 ] - 1 )
      idx_1++;
     else {
      idx_1 = 0;
      idx_0++;
     }
    }
   }
   else
    throw( std::runtime_error( "Unsupported group type" ) );
  }
  else if( type == 0 ) {
   // Multi arrays of type 0 (i.e. multi_array< T >) store
   // elements in sequential cells. Thus, we can store them
   // as usually done for std::vector< T > by only keeping track
   // of the first element and storing the number of non empty
   // cells in the structure.
   auto ma = get_multi_array0( gr , un_any_type< T >() ,
                              un_any_int< 2 >() );

   if( typeid( T * ) == typeid( FRowConstraint * ) ) {
    // Constraint group

    // Scan the linearization of the array
    for( auto v = ma->data() ; idx_0 < ma->shape()[ 0 ] ; ++v ) {
     // Scanning a group of Constraints
     auto scan = [ this , num_node ]
      ( const FRowConstraint & c ) {
        scan_constraint( c , num_node );
     };
     // Scan all the constraints one at a time
     un_any_const_static( v , scan  , un_any_type< FRowConstraint >() );

     // The linearization produced by ma->data() for the 2D multi_array
     // stores elements in row-major order. Therefore, we should increment
     // the column index first, and when it exceeds the number of columns,
     // reset it and increment the row index.
     if( idx_1 < ma->shape()[ 1 ] - 1 )
      idx_1++;
     else {
      idx_1 = 0;
      idx_0++;
     }
    }
   }
   else if( typeid( T * ) == typeid( ColVariable * ) ) {
    // Variable group

    // Scan the linearization of the array
    for( auto v = ma->data() ; idx_0 < ma->shape()[ 0 ] ; ++v ) {
     // Scanning a group of Variables
     // in this case we simply add the variable to the corresponding node
     auto push_var_toNode = [ this , num_node ]
       ( const ColVariable & c ) {
        n_varNode[ num_node ] += 1;
        varNode[ num_node ].push_back( &c );
        var_to_node.emplace( &c , num_node );
     };
     un_any_const_static( v , push_var_toNode , 
      un_any_type< ColVariable >() );

     // The linearization produced by ma->data() for the 2D multi_array
     // stores elements in row-major order. Therefore, we should increment
     // the column index first, and when it exceeds the number of columns,
     // reset it and increment the row index.
     if( idx_1 < ma->shape()[ 1 ] - 1 )
      idx_1++;
     else {
      idx_1 = 0;
      idx_0++;
     }
    }
   }
   else
    throw( std::runtime_error( "Unsupported group type" ) );
  }
  else
   throw( std::runtime_error( "Unsupported multi-array type" ) );
 }
 else if( ma_dim == 3 ) {
  // Use the 3D multi_array

  // Get type of multi array. See MILPSolver.h:1256 for further details.
  int type = get_multi_array_type( gr , un_any_type< T >() ,
                un_any_int< 3 >() );

  // Indices of the 3 dimensions
  int idx_0 = 0;
  int idx_1 = 0;
  int idx_2 = 0;

  if( type == 1 ) {
   // Multi arrays of type 1 (i.e. multi_array< std::vector < T * > >).
   // In this case elements are not stored in sequential cells. Thus, 
   // we have to "unpack" each std::vector and store them separately.
   auto ma = get_multi_array1( gr , un_any_type< T >() ,
                            un_any_int< 3 >() );

   if( typeid( T * ) == typeid( FRowConstraint * ) ) {
    // Constraint group

    // Scan the linearization of the array
    for( auto v = ma->data() ; idx_0 < ma->shape()[ 0 ] ; ++v ) {
     // Scanning a group of Constraints
     auto scan = [ this , num_node ]
      ( const FRowConstraint & c ) {
        scan_constraint( c , num_node );
     };
     // Scan all the constraints one at a time
     un_any_const_static( v , scan  , un_any_type< FRowConstraint >() );

     // The linearization produced by ma->data() for the 3D multi_array
     // stores elements in row-major order.
     if( idx_2 < ma->shape()[ 2 ] - 1 )
      idx_2++;
     else if( idx_1 < ma->shape()[ 1 ] - 1 ) {
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
    // Variable group

    // Scan the linearization of the array
    for( auto v = ma->data() ; idx_0 < ma->shape()[ 0 ] ; ++v ) {
     // Scanning a group of Variables
     // in this case we simply add the variable to the corresponding node
     auto push_var_toNode = [ this , num_node ]
       ( const ColVariable & c ) {
        n_varNode[ num_node ] += 1;
        varNode[ num_node ].push_back( &c );
        var_to_node.emplace( &c , num_node );
     };
     un_any_const_static( v , push_var_toNode , 
      un_any_type< ColVariable >() );

     // The linearization produced by ma->data() for the 3D multi_array
     // stores elements in row-major order.
     if( idx_2 < ma->shape()[ 2 ] - 1 )
      idx_2++;
     else if( idx_1 < ma->shape()[ 1 ] - 1 ) {
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
    throw( std::runtime_error( "Unsupported group type" ) );
  }
  else if( type == 0 ) {
   // Multi arrays of type 0 (i.e., multi_array< T >) store
   // elements in sequential cells. Thus, we can store them
   // as usually done for std::vector< T > by only keeping track
   // of the first element and storing the number of non-empty
   // cells in the structure.
   auto ma = get_multi_array0( gr , un_any_type< T >() ,
                              un_any_int< 3 >() );

   if( typeid( T * ) == typeid( FRowConstraint * ) ) {
    // Constraint group

    // Scan the linearization of the array
    for( auto v = ma->data() ; idx_0 < ma->shape()[ 0 ] ; ++v ) {
     // Scanning a group of Constraints
     auto scan = [ this , num_node ]
      ( const FRowConstraint & c ) {
        scan_constraint( c , num_node );
     };
     // Scan all the constraints one at a time
     un_any_const_static( v , scan  , un_any_type< FRowConstraint >() );

     // The linearization produced by ma->data() for the 3D multi_array
     // stores elements in row-major order.
     if( idx_2 < ma->shape()[ 2 ] - 1 )
      idx_2++;
     else if( idx_1 < ma->shape()[ 1 ] - 1 ) {
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
    // Variable group

    // Scan the linearization of the array
    for( auto v = ma->data() ; idx_0 < ma->shape()[ 0 ] ; ++v ) {
     // Scanning a group of Variables
     // in this case we simply add the variable to the corresponding node
     auto push_var_toNode = [ this , num_node ]
       ( const ColVariable & c ) {
        n_varNode[ num_node ] += 1;
        varNode[ num_node ].push_back( &c );
        var_to_node.emplace( &c , num_node );
     };
     un_any_const_static( v , push_var_toNode , 
      un_any_type< ColVariable >() );

     // The linearization produced by ma->data() for the 3D multi_array
     // stores elements in row-major order.
     if( idx_2 < ma->shape()[ 2 ] - 1 )
      idx_2++;
     else if( idx_1 < ma->shape()[ 1 ] - 1 ) {
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
    throw( std::runtime_error( "Unsupported group type" ) );
  }
  else
   throw( std::runtime_error( "Unsupported multi-array type" ) );
 }
 else
    // Handle invalid or unsupported ma_dim 
    return; 

} 

/*--------------------------------------------------------------------------*/

void PIPSMILPSolver::scan_constraint( const FRowConstraint & con , 
                                      Index num_node )
{
 // Check if the Constraint is not empty
 if( con.get_Block() == nullptr )
  throw( std::invalid_argument( "The provided constraint is empty" ) );

 // First of all we can categorize the constraint into equality/inequality
 bool is_eq;
 auto con_lhs = con.get_lhs();
 auto con_rhs = con.get_rhs();

 if( con_lhs == con_rhs )
  is_eq = true;
 else
  is_eq = false;

 // Now we need to understand if the constraint involves only variables of the 
 // same node (accepting also root node) or also variables coming from 
 // multiple leaves

 // Initialize boolean variable to understand if exists any variable in the
 // constraint that belongs to another node (excluding root)
 bool found_inOtherNode = false;
 if( auto f = con.get_function() ) {
  if( auto lf = dynamic_cast< const LinearFunction * >( f ) ) {
   for( auto el : lf->get_v_var() ){
    auto * v = dynamic_cast< ColVariable * >( std::get< 0 >( el ) );

    auto it = var_to_node.find( v );

    if( it == var_to_node.end() )
      throw std::runtime_error( "Variable not assigned to any PIPS node" );

    int owner = it->second;

    if( owner != num_node && owner != 0 ) {
     found_inOtherNode = true;
     break;
    }
   }
  }
  else
      throw( std::invalid_argument( "Unexpected constraint type" ) );
 }
 
 // Now add the constraint to the correct group
 if( found_inOtherNode && is_eq ){
  // Add to the equality linking constraints
  n_LinkEqCons++;
  LinkEqCons.push_back( &con );
 }
 else if( found_inOtherNode && !is_eq ){
  // Add to the inequality linking constraints
  n_LinkInEqCons++;
  LinkInEqCons.push_back( &con );
 }
 else if( ! found_inOtherNode && is_eq ){
  // Add to the node equality constraints
  n_EqConsNode[ num_node ]++;
  EqConsNode[ num_node ].push_back( &con );
 }
 else{
  // Add to the node inequality constraints
  n_InEqConsNode[ num_node ]++;
  InEqConsNode[ num_node ].push_back( &con );
 }
}

/*--------------------------------------------------------------------------*/

int PIPSMILPSolver::ExtractRhsVector( int id , double* vec , int len ,
                              std::vector< const FRowConstraint * > node_cons ,
                              const int nCons ,
                              std::vector< double > rhs ,
                              std::vector< char > sense ,
                              std::vector< double > ranges ){ 
  if( id < n_nodes ){
   // Check if the node has own constraint
   if( nCons == 0 ){
    if( len > 0 ){
    // Strange, but we fill the expected vector with zeros
     for (int i = 0; i < len; i++)
      vec[i] = 0.0;
    }
    // Nothing to do
    return 0;
   }

   // Extract the global indices of constraints
   std::vector< int > global_idxs_cons = compute_cons_global_idxs( 
                                          node_cons , nCons );
   
   // Extract the corresponding elements and store them
   std::vector< double > sub_rhs;
   sub_rhs.reserve( global_idxs_cons.size() );

   if( sense.empty() ){
    // We are asking the rhs for equality constraints, just select them
    for( auto idx : global_idxs_cons )
     if( idx < rhs.size() )
      sub_rhs.push_back( rhs[ idx ] );
   }
   else{
    // Inequality constraints, we have to distinguish several cases
    for( int i = 0 ; i < global_idxs_cons.size() ; i++ ){
     int idx = global_idxs_cons[ i ]; // true constraint index
     if( idx < rhs.size() ){ 
      if( sense[ idx ] == 'L' )
       // <= RHS, store the element
       sub_rhs.push_back( rhs[ idx ] );
      else if( sense[ idx ] == 'G' )
       // >= LHS, so no elements to store
       sub_rhs.push_back( 0.0 );
      else if( sense[ idx ] == 'R' ){
       // LHS <= f <= RHS, we have to evaluate the sign of the range
       if( ranges[ idx ] > 0 )
        sub_rhs.push_back( rhs[ idx ] + ranges[ idx ] );
       else
        sub_rhs.push_back( rhs[ idx ] );
      }
      else
       throw( std::runtime_error( "Unexpected inequality sense" ) );
     }
    }
   }
   if( static_cast< int >( sub_rhs.size() ) != len )
    PIPS_CALLBACK_COUT << "[SANITY WARNING] ExtractRhsVector id=" << id
              << " built " << sub_rhs.size()
              << " values but PIPS requested len=" << len << std::endl;
   std::copy( sub_rhs.begin(), sub_rhs.end(), vec );
   sanity_check_vector( "ExtractRhsVector", id, vec, len, nCons );
   return 0;
  }
  else
   throw( std::runtime_error( "Node ID outside the expected range" ) );
}

/*--------------------------------------------------------------------------*/

int PIPSMILPSolver::ExtractLhsVector( int id , double* vec , int len ,
                              std::vector< const FRowConstraint * > node_cons ,
                              const int nCons ,
                              std::vector< double > rhs ,
                              std::vector< char > sense ,
                              std::vector< double > ranges ){ 
  if( id < n_nodes ){
   // Check if the node has own constraint
   if( nCons == 0 ){
    if( len > 0 ){
    // Strange, but we fill the expected vector with zeros
     for (int i = 0; i < len; i++)
      vec[i] = 0.0;
    }
    // Nothing to do
    return 0;
   }

   // Extract the global indices of constraints
   std::vector< int > global_idxs_cons = compute_cons_global_idxs( 
                                          node_cons , nCons );
   
   // Extract the corresponding elements and store them
   std::vector< double > sub_rhs;
   sub_rhs.reserve( global_idxs_cons.size() );

   for( int i = 0 ; i < global_idxs_cons.size() ; i++ ){
    int idx = global_idxs_cons[ i ]; // true constraint index
    if( idx < rhs.size() ){ 
     if( sense[ idx ] == 'G' )
      // >= LHS, store the element
      sub_rhs.push_back( rhs[ idx ] );
     else if( sense[ idx ] == 'L' )
      // <= RHS, so no elements to store
      sub_rhs.push_back( 0.0 );
     else if( sense[ idx ] == 'R' ){
      // LHS <= f <= RHS, we have to evaluate the sign of the range
      if( ranges[ idx ] > 0 )
       sub_rhs.push_back( rhs[ idx ] );
      else
       sub_rhs.push_back( rhs[ idx ] + ranges[ idx ] );
     }
     else
      throw( std::runtime_error( "Unexpected inequality sense" ) );
    }
   }
   if( static_cast< int >( sub_rhs.size() ) != len )
    PIPS_CALLBACK_COUT << "[SANITY WARNING] ExtractLhsVector id=" << id
              << " built " << sub_rhs.size()
              << " values but PIPS requested len=" << len << std::endl;
   std::copy( sub_rhs.begin(), sub_rhs.end(), vec );
   sanity_check_vector( "ExtractLhsVector", id, vec, len, nCons );
   return 0;
  }
  else
   throw( std::runtime_error( "Node ID outside the expected range" ) );
}

/*--------------------------------------------------------------------------*/

int PIPSMILPSolver::ExtractRhsActiveFlag( int id , double* vec , int len ,
                              std::vector< const FRowConstraint * > node_cons ,
                              const int nCons ,
                              std::vector< double > rhs ,
                              std::vector< char > sense ,
                              std::vector< double > ranges ){ 
  if( id < n_nodes ){
   // Check if the node has own constraint
   if( nCons == 0 ){
    if( len > 0 ){
    // Strange, but we fill the expected vector with zeros
     for (int i = 0; i < len; i++)
      vec[i] = 0.0;
    }
    // Nothing to do
    return 0;
   }

   // Extract the global indices of constraints
   std::vector< int > global_idxs_cons = compute_cons_global_idxs( 
                                          node_cons , nCons );
   
   // Extract the corresponding elements and store them
   std::vector< double > sub_rhs;
   sub_rhs.reserve( global_idxs_cons.size() );

   if( sense.empty() ){
    // We are asking the rhs for equality constraints, just set flag to true
    for( auto idx : global_idxs_cons )
     if( idx < rhs.size() )
      sub_rhs.push_back( 1.0 );
   }
   else{
    // Inequality constraints, we have to distinguish several cases
    for( int i = 0 ; i < global_idxs_cons.size() ; i++ ){
     int idx = global_idxs_cons[ i ]; // true constraint index
     if( idx < rhs.size() ){ 
      if( sense[ idx ] == 'G' )
       // >= LHS, set false
       sub_rhs.push_back( 0.0 );
      else 
       // In all the other cases set rue
       sub_rhs.push_back( 1.0 );
     }
    }
   }
   if( static_cast< int >( sub_rhs.size() ) != len )
    PIPS_CALLBACK_COUT << "[SANITY WARNING] ExtractRhsActiveFlag id=" << id
              << " built " << sub_rhs.size()
              << " values but PIPS requested len=" << len << std::endl;
   std::copy( sub_rhs.begin(), sub_rhs.end(), vec );
   sanity_check_vector( "ExtractRhsActiveFlag", id, vec, len, nCons );
   return 0;
  }
  else
   throw( std::runtime_error( "Node ID outside the expected range" ) );
}

/*--------------------------------------------------------------------------*/

int PIPSMILPSolver::ExtractLhsActiveFlag( int id , double* vec , int len ,
                              std::vector< const FRowConstraint * > node_cons ,
                              const int nCons ,
                              std::vector< double > rhs ,
                              std::vector< char > sense ,
                              std::vector< double > ranges ){ 
  if( id < n_nodes ){
   // Check if the node has own constraint
   if( nCons == 0 ){
    if( len > 0 ){
    // Strange, but we fill the expected vector with zeros
     for (int i = 0; i < len; i++)
      vec[i] = 0.0;
    }
    // Nothing to do
    return 0;
   }

   // Extract the global indices of constraints
   std::vector< int > global_idxs_cons = compute_cons_global_idxs( 
                                          node_cons , nCons );
   
   // Extract the corresponding elements and store them
   std::vector< double > sub_rhs;
   sub_rhs.reserve( global_idxs_cons.size() );

   for( int i = 0 ; i < global_idxs_cons.size() ; i++ ){
    int idx = global_idxs_cons[ i ]; // true constraint index
    if( idx < rhs.size() ){ 
     if( sense[ idx ] == 'L' )
      // <= RHS, set to false
      sub_rhs.push_back( 0.0 );
     else 
      // In all the other cases, set true
      sub_rhs.push_back( 1.0 );
    }
   }
   if( static_cast< int >( sub_rhs.size() ) != len )
    PIPS_CALLBACK_COUT << "[SANITY WARNING] ExtractLhsActiveFlag id=" << id
              << " built " << sub_rhs.size()
              << " values but PIPS requested len=" << len << std::endl;
   std::copy( sub_rhs.begin(), sub_rhs.end(), vec );
   sanity_check_vector( "ExtractLhsActiveFlag", id, vec, len, nCons );
   return 0;
  }
  else
   throw( std::runtime_error( "Node ID outside the expected range" ) );
}

/*--------------------------------------------------------------------------*/

int PIPSMILPSolver::ExtractVarBounds( int id , double* vec , int len ,
                              std::vector< const ColVariable * > node_vars ,
                              const int nVars ,
                              std::vector< double > bounds ){ 
  if( id < n_nodes ){
   // Check if the node has own constraint
   if( nVars == 0 ){
    if( len > 0 ){
    // Strange, but we fill the expected vector with zeros
     for (int i = 0; i < len; i++)
      vec[i] = 0.0;
    }
    // Nothing to do
    return 0;
   }

   // Extract the global indices of variables
   std::vector< int > global_idxs_vars = compute_vars_global_idxs( 
                                          node_vars , nVars );
   
   // Extract the corresponding elements and store them
   std::vector< double > sub_bounds;
   sub_bounds.reserve( global_idxs_vars.size() );
    
   for( auto idx : global_idxs_vars )
    if( idx < bounds.size() ){
     if( bounds[ idx ] == -Inf< double >() || 
          bounds[ idx ] == Inf< double >() )
      sub_bounds.push_back( 0.0 );
     else
      sub_bounds.push_back( bounds[ idx ] );
    }
   
   if( static_cast< int >( sub_bounds.size() ) != len )
    PIPS_CALLBACK_COUT << "[SANITY WARNING] ExtractVarBounds id=" << id
              << " built " << sub_bounds.size()
              << " values but PIPS requested len=" << len << std::endl;
   std::copy( sub_bounds.begin(), sub_bounds.end(), vec );
   sanity_check_vector( "ExtractVarBounds", id, vec, len, nVars );
   return 0;
  }
  else
   throw( std::runtime_error( "Node ID outside the expected range" ) );
}

/*--------------------------------------------------------------------------*/

int PIPSMILPSolver::ExtractFlagVarBounds( int id , double* vec , int len ,
                                std::vector< const ColVariable * > node_vars ,
                                const int nVars ,
                                std::vector< double > bounds ){ 
  if( id < n_nodes ){
   // Check if the node has own constraint
   if( nVars == 0 ){
    if( len > 0 ){
    // Strange, but we fill the expected vector with zeros
     for (int i = 0; i < len; i++)
      vec[i] = 0.0;
    }
    // Nothing to do
    return 0;
   }

   // Extract the global indices of variables
   std::vector< int > global_idxs_vars = compute_vars_global_idxs( 
                                          node_vars , nVars );
   
   // Extract the corresponding elements and store them
   std::vector< double > sub_bounds;
   sub_bounds.reserve( global_idxs_vars.size() );
    
   for( auto idx : global_idxs_vars )
    if( idx < bounds.size() ){
     if( bounds[ idx ] == -Inf< double >() || 
          bounds[ idx ] == Inf< double >() )
      sub_bounds.push_back( 0.0 );
     else
      sub_bounds.push_back( 1.0 );
    }
   
   if( static_cast< int >( sub_bounds.size() ) != len )
    PIPS_CALLBACK_COUT << "[SANITY WARNING] ExtractFlagVarBounds id=" << id
              << " built " << sub_bounds.size()
              << " values but PIPS requested len=" << len << std::endl;
   std::copy( sub_bounds.begin(), sub_bounds.end(), vec );
   sanity_check_vector( "ExtractFlagVarBounds", id, vec, len, nVars );
   return 0;
  }
  else
   throw( std::runtime_error( "Node ID outside the expected range" ) );
}

/*--------------------------------------------------------------------------*/

int PIPSMILPSolver::ExtractObj( int id , double* vec , int len ,
                                std::vector< const ColVariable * > node_vars ,
                                const int nVars ,
                                std::vector< double > obj_value ,
                                int objsense ){ 
  if( id < n_nodes ){
   // Check if the node has own constraint
   if( nVars == 0 ){
    if( len > 0 ){
    // Strange, but we fill the expected vector with zeros
     for (int i = 0; i < len; i++)
      vec[i] = 0.0;
    }
    // Nothing to do
    return 0;
   }

   // Extract the global indices of variables
   std::vector< int > global_idxs_vars = compute_vars_global_idxs( 
                                          node_vars , nVars );
   
   // Extract the corresponding elements and store them
   std::vector< double > sub_obj;
   sub_obj.reserve( global_idxs_vars.size() );
    
   for( auto idx : global_idxs_vars )
    if( idx < obj_value.size() )
     // PIPS always assume to work with a minimization model, so we must
     // approprtiately convert MILPSolver objective values
     sub_obj.push_back( obj_value[ idx ] * objsense );
   
   std::copy( sub_obj.begin(), sub_obj.end(), vec );
   sanity_check_vector( "ExtractObj", id, vec, len, nVars );
   return 0;
  }
  else
   throw( std::runtime_error( "Node ID outside the expected range" ) );
}

/*--------------------------------------------------------------------------*/

int PIPSMILPSolver::No_VarinNode( void * user_data, int id, int* nnz){
  // Use the current class in the callback
  auto * solver = static_cast< PIPSMILPSolver * >( user_data );
  sanity_check_id( "No_VarinNode", id, solver->n_nodes );
  if( id < solver->n_nodes ){
    *nnz = solver->n_varNode[ id ];
    sanity_check_count( "No_VarinNode", id, *nnz );
    return 0;
  }
  else
    throw( std::runtime_error( "Tried to query number of variables outside"
      " the current node count" ) );
}

/*--------------------------------------------------------------------------*/

int PIPSMILPSolver::No_EqConsinNode( void * user_data, int id, int* nnz){
  // Use the current class in the callback
  auto * solver = static_cast< PIPSMILPSolver * >( user_data );
  sanity_check_id( "No_EqConsinNode", id, solver->n_nodes );
  if( id < solver->n_nodes ){
    *nnz = solver->n_EqConsNode[ id ];
    sanity_check_count( "No_EqConsinNode", id, *nnz );
    return 0;
  }
  else
    throw( std::runtime_error( "Tried to query number of equality constraints"
      " outside the current node count" ) );
}

/*--------------------------------------------------------------------------*/

int PIPSMILPSolver::No_InEqConsinNode( void * user_data, int id, int* nnz){
  // Use the current class in the callback
  auto * solver = static_cast< PIPSMILPSolver * >( user_data );
  sanity_check_id( "No_InEqConsinNode", id, solver->n_nodes );
  if( id < solver->n_nodes ){
    *nnz = solver->n_InEqConsNode[ id ];
    sanity_check_count( "No_InEqConsinNode", id, *nnz );
    return 0;
  }
  else
    throw( std::runtime_error( "Tried to query number of equality constraints"
      " outside the current node count" ) );
}

/*--------------------------------------------------------------------------*/

int PIPSMILPSolver::No_LinkEqCons( void * user_data, int id, int* nnz){
  // Use the current class in the callback
  auto * solver = static_cast< PIPSMILPSolver * >( user_data );
  sanity_check_id( "No_LinkEqCons", id, solver->n_nodes );
  *nnz = solver->n_LinkEqCons;
  sanity_check_count( "No_LinkEqCons", id, *nnz );
  return 0;
}

/*--------------------------------------------------------------------------*/

int PIPSMILPSolver::No_LinkInEqCons( void * user_data, int id, int* nnz){
  // Use the current class in the callback
  auto * solver = static_cast< PIPSMILPSolver * >( user_data );
  sanity_check_id( "No_LinkInEqCons", id, solver->n_nodes );
  *nnz = solver->n_LinkInEqCons;
  sanity_check_count( "No_LinkInEqCons", id, *nnz );
  return 0;
}

/*--------------------------------------------------------------------------*/

int PIPSMILPSolver::nnzEqConsDiag( void * user_data, int id , int* nnz ){
  // Use the current class in the callback
  auto * solver = static_cast< PIPSMILPSolver * >( user_data );
  solver->EvaluateNnz( id , nnz , solver->EqConsNode[ id ] , 
    solver->n_EqConsNode[ id ] , solver->varNode[ id ] , 
    solver->n_varNode[ id ] );
  sanity_check_count( "nnzEqConsDiag", id, *nnz );
  return 0;
}

/*--------------------------------------------------------------------------*/

int PIPSMILPSolver::nnzEqConsVert( void * user_data, int id , int* nnz ){
  // In this case we are extracting the sub-matrix for a particular Node but
  // considering the root variables (i.e. T_1, T_2, ..., T_n)
  if( id == 0 ){
    // T_0 can already be obtained through nnzEqConsDiag( node_id = 0 )
    *nnz = 0;
    return 0;
  }
  else{
   // Use the current class in the callback
   auto * solver = static_cast< PIPSMILPSolver * >( user_data );
   solver->EvaluateNnz( id , nnz , solver->EqConsNode[ id ] , 
    solver->n_EqConsNode[ id ] , solver->varNode[ 0 ] , 
    solver->n_varNode[ 0 ] );
  sanity_check_count( "nnzEqConsVert", id, *nnz );
   return 0;
  }
}

/*--------------------------------------------------------------------------*/

int PIPSMILPSolver::nnzInEqConsDiag( void * user_data, int id , int* nnz ){
  // Use the current class in the callback
  auto * solver = static_cast< PIPSMILPSolver * >( user_data );
  solver->EvaluateNnz( id , nnz , solver->InEqConsNode[ id ] , 
    solver->n_InEqConsNode[ id ] , solver->varNode[ id ] , 
    solver->n_varNode[ id ] );
  sanity_check_count( "nnzInEqConsDiag", id, *nnz );
  return 0;
}

/*--------------------------------------------------------------------------*/

int PIPSMILPSolver::nnzInEqConsVert( void * user_data, int id , int* nnz ){
  // In this case we are extracting the sub-matrix for a particular Node but
  // considering the root variables (i.e. T_1, T_2, ..., T_n)
  if( id == 0 ){
    // T_0 can already be obtained through nnzEqConsDiag( node_id = 0 )
    *nnz = 0;
    return 0;
  }
  else{
   // Use the current class in the callback
   auto * solver = static_cast< PIPSMILPSolver * >( user_data );
   solver->EvaluateNnz( id , nnz , solver->InEqConsNode[ id ] , 
    solver->n_InEqConsNode[ id ] , solver->varNode[ 0 ] , 
    solver->n_varNode[ 0 ] );
  sanity_check_count( "nnzInEqConsVert", id, *nnz );
   return 0;
  }
}

/*--------------------------------------------------------------------------*/

int PIPSMILPSolver::nnzLinkEqCons( void * user_data, int id , int* nnz ){
  // Use the current class in the callback
  auto * solver = static_cast< PIPSMILPSolver * >( user_data );
  solver->EvaluateNnz( id , nnz , solver->LinkEqCons , 
    solver->n_LinkEqCons , solver->varNode[ id ] , 
    solver->n_varNode[ id ] );
  sanity_check_count( "nnzLinkEqCons", id, *nnz );
  return 0;
}

/*--------------------------------------------------------------------------*/

int PIPSMILPSolver::nnzLinkInEqCons( void * user_data, int id , int* nnz ){
  // Use the current class in the callback
  auto * solver = static_cast< PIPSMILPSolver * >( user_data );
  solver->EvaluateNnz( id , nnz , solver->LinkInEqCons , 
    solver->n_LinkInEqCons , solver->varNode[ id ] , 
    solver->n_varNode[ id ] );
  sanity_check_count( "nnzLinkInEqCons", id, *nnz );
  return 0;
}

/*--------------------------------------------------------------------------*/

int PIPSMILPSolver::nnzAllZero(void*, int, int* nnz) {
   *nnz = 0;
   return 0;
}

/*--------------------------------------------------------------------------*/

int PIPSMILPSolver::MatEqConsDiag( void * user_data, int id , int* krowM, 
                                    int* jcolM, double* M ){

  // Use the current class in the callback
  auto * solver = static_cast< PIPSMILPSolver * >( user_data );

  PIPS_CALLBACK_COUT << "\n[PIPS CALLBACK ENTER] MatEqConsDiag id=" << id 
    << std::endl;
  solver->ExtractMatrix( id , krowM, jcolM , M , solver->EqConsNode[ id ] , 
                solver->n_EqConsNode[ id ] , solver->varNode[ id ] , 
                solver->n_varNode[ id ] );
  return 0;
}

/*--------------------------------------------------------------------------*/

int PIPSMILPSolver::MatEqConsVert( void * user_data, int id , int* krowM, 
                                    int* jcolM, double* M ){
  // In this case we are extracting the sub-matrix for a particular Node but
  // considering the root variables (i.e. T_1, T_2, ..., T_n)
  if( id == 0 ){
    // T_0 can already be obtained through nnzEqConsDiag( node_id = 0 )
    // Nothing to do
    return 0;
  }
  else{
   // Use the current class in the callback
   auto * solver = static_cast< PIPSMILPSolver * >( user_data );
   PIPS_CALLBACK_COUT << "\n[PIPS CALLBACK ENTER] MatEqConsVert id=" << id 
    << std::endl;
  solver->ExtractMatrix( id , krowM, jcolM , M , solver->EqConsNode[ id ] , 
                solver->n_EqConsNode[ id ] , solver->varNode[ 0 ] , 
                solver->n_varNode[ 0 ] );
   return 0;
  }
}

/*--------------------------------------------------------------------------*/

int PIPSMILPSolver::MatInEqConsDiag( void * user_data, int id , int* krowM, 
                                    int* jcolM, double* M ){
  // Use the current class in the callback
  auto * solver = static_cast< PIPSMILPSolver * >( user_data );

  PIPS_CALLBACK_COUT << "\n[PIPS CALLBACK ENTER] MatInEqConsDiag id=" << id 
    << std::endl;
  solver->ExtractMatrix( id , krowM, jcolM , M , solver->InEqConsNode[ id ] , 
                solver->n_InEqConsNode[ id ] , solver->varNode[ id ] , 
                solver->n_varNode[ id ] );
  return 0;
}

/*--------------------------------------------------------------------------*/

int PIPSMILPSolver::MatInEqConsVert( void * user_data, int id , int* krowM, 
                                    int* jcolM, double* M ){
  // In this case we are extracting the sub-matrix for a particular Node but
  // considering the root variables (i.e. T_1, T_2, ..., T_n)
  if( id == 0 ){
    // T_0 can already be obtained through nnzEqConsDiag( node_id = 0 )
    // Nothing to do
    return 0;
  }
  else{
   // Use the current class in the callback
   auto * solver = static_cast< PIPSMILPSolver * >( user_data );

   PIPS_CALLBACK_COUT << "\n[PIPS CALLBACK ENTER] MatInEqConsVert id=" << id 
    << std::endl;
  solver->ExtractMatrix( id , krowM, jcolM , M , solver->InEqConsNode[ id ] , 
                solver->n_InEqConsNode[ id ] , solver->varNode[ 0 ] , 
                solver->n_varNode[ 0 ] );
   return 0;
  }
}

/*--------------------------------------------------------------------------*/

int PIPSMILPSolver::MatLinkEqCons( void * user_data, int id , int* krowM, 
                                    int* jcolM, double* M ){
  // Use the current class in the callback
  auto * solver = static_cast< PIPSMILPSolver * >( user_data );

  PIPS_CALLBACK_COUT << "\n[PIPS CALLBACK ENTER] MatLinkEqCons id=" << id 
    << std::endl;
  solver->ExtractMatrix( id , krowM, jcolM , M , solver->LinkEqCons , 
                solver->n_LinkEqCons , solver->varNode[ id ] , 
                solver->n_varNode[ id ] );
  return 0;
}

/*--------------------------------------------------------------------------*/

int PIPSMILPSolver::MatLinkInEqCons( void * user_data, int id , int* krowM, 
                                    int* jcolM, double* M ){
  // Use the current class in the callback
  auto * solver = static_cast< PIPSMILPSolver * >( user_data );

  PIPS_CALLBACK_COUT << "\n[PIPS CALLBACK ENTER] MatLinkInEqCons id=" << id 
    << std::endl;
  solver->ExtractMatrix( id , krowM, jcolM , M , solver->LinkInEqCons , 
                solver->n_LinkInEqCons , solver->varNode[ id ] , 
                solver->n_varNode[ id ] );
  return 0;
}

/*--------------------------------------------------------------------------*/

int PIPSMILPSolver::matAllZero(void*, int, int*, int*, double*) {
   return 0;
}

/*--------------------------------------------------------------------------*/

int PIPSMILPSolver::ObjVars( void * user_data, int id , double* vec,
                                int len ){

  PIPS_CALLBACK_COUT << "\n[PIPS CALLBACK ENTER] ObjVars id=" << id << " len=" 
    << len << std::endl;

  // Use the current class in the callback
  auto * solver = static_cast< PIPSMILPSolver * >( user_data );

  solver->ExtractObj( id , vec, len , solver->varNode[ id ] , 
                solver->n_varNode[ id ] , solver->objective ,
                solver->objsense );
  return 0;
}

/*--------------------------------------------------------------------------*/

int PIPSMILPSolver::RhsEqCons( void * user_data, int id , double* vec,
                                int len ){
  PIPS_CALLBACK_COUT << "\n[PIPS CALLBACK ENTER] RhsEqCons id=" << id 
    << " len=" << len << std::endl;

  // Use the current class in the callback
  auto * solver = static_cast< PIPSMILPSolver * >( user_data );

  solver->ExtractRhsVector( id , vec, len , solver->EqConsNode[ id ] , 
                solver->n_EqConsNode[ id ] , solver-> rhs ,
                {} , {} );
  return 0;
}

/*--------------------------------------------------------------------------*/

int PIPSMILPSolver::RhsInEqCons( void * user_data, int id , double* vec,
                                int len ){
  PIPS_CALLBACK_COUT << "\n[PIPS CALLBACK ENTER] RhsInEqCons id=" << id 
    << " len=" << len << std::endl;

  // Use the current class in the callback
  auto * solver = static_cast< PIPSMILPSolver * >( user_data );

  solver->ExtractRhsVector( id , vec, len , solver->InEqConsNode[ id ] , 
                solver->n_InEqConsNode[ id ] , solver-> rhs ,
                solver->sense , solver->rngval );
  return 0;
}

/*--------------------------------------------------------------------------*/

int PIPSMILPSolver::LhsInEqCons( void * user_data, int id , double* vec,
                                int len ){
  PIPS_CALLBACK_COUT << "\n[PIPS CALLBACK ENTER] LhsInEqCons id=" << id 
    << " len=" << len << std::endl;

  // Use the current class in the callback
  auto * solver = static_cast< PIPSMILPSolver * >( user_data );

  solver->ExtractLhsVector( id , vec, len , solver->InEqConsNode[ id ] , 
                solver->n_InEqConsNode[ id ] , solver-> rhs ,
                solver->sense , solver->rngval );
  return 0;
}

/*--------------------------------------------------------------------------*/

int PIPSMILPSolver::RhsLinkEqCons( void * user_data, int id , double* vec,
                                int len ){
  PIPS_CALLBACK_COUT << "\n[PIPS CALLBACK ENTER] RhsLinkEqCons id=" << id 
    << " len=" << len << std::endl;

  // Use the current class in the callback
  auto * solver = static_cast< PIPSMILPSolver * >( user_data );

  solver->ExtractRhsVector( id , vec, len , solver->LinkEqCons , 
                solver->n_LinkEqCons , solver-> rhs , {} , {} );
  return 0;
}

/*--------------------------------------------------------------------------*/

int PIPSMILPSolver::RhsLinkInEqCons( void * user_data, int id , double* vec,
                                int len ){
  PIPS_CALLBACK_COUT << "\n[PIPS CALLBACK ENTER] RhsLinkInEqCons id=" << id 
    << " len=" << len << std::endl;

  // Use the current class in the callback
  auto * solver = static_cast< PIPSMILPSolver * >( user_data );

  solver->ExtractRhsVector( id , vec, len , solver->LinkInEqCons , 
                solver->n_LinkInEqCons , solver->rhs ,
                solver->sense , solver->rngval );
  return 0;
}

/*--------------------------------------------------------------------------*/

int PIPSMILPSolver::LhsLinkInEqCons( void * user_data, int id , double* vec,
                                int len ){
  PIPS_CALLBACK_COUT << "\n[PIPS CALLBACK ENTER] LhsLinkInEqCons id=" << id 
    << " len=" << len << std::endl;

  // Use the current class in the callback
  auto * solver = static_cast< PIPSMILPSolver * >( user_data );

  solver->ExtractLhsVector( id , vec, len , solver->LinkInEqCons , 
                solver->n_LinkInEqCons , solver-> rhs ,
                solver->sense , solver->rngval );
  return 0;
}

/*--------------------------------------------------------------------------*/

int PIPSMILPSolver::FlagRhsInEqCons( void * user_data, int id , double* vec,
                                int len ){
  PIPS_CALLBACK_COUT << "\n[PIPS CALLBACK ENTER] FlagRhsInEqCons id=" << id 
    << " len=" << len << std::endl;

  // Use the current class in the callback
  auto * solver = static_cast< PIPSMILPSolver * >( user_data );

  solver->ExtractRhsActiveFlag( id , vec, len , solver->InEqConsNode[ id ] , 
                solver->n_InEqConsNode[ id ] , solver-> rhs ,
                solver->sense , solver->rngval );
  return 0;
}

/*--------------------------------------------------------------------------*/

int PIPSMILPSolver::FlagLhsInEqCons( void * user_data, int id , double* vec,
                                int len ){
  PIPS_CALLBACK_COUT << "\n[PIPS CALLBACK ENTER] FlagLhsInEqCons id=" << id 
    << " len=" << len << std::endl;

  // Use the current class in the callback
  auto * solver = static_cast< PIPSMILPSolver * >( user_data );

  solver->ExtractLhsActiveFlag( id , vec, len , solver->InEqConsNode[ id ] , 
                solver->n_InEqConsNode[ id ] , solver-> rhs ,
                solver->sense , solver->rngval );
  return 0;
}

/*--------------------------------------------------------------------------*/

int PIPSMILPSolver::FlagRhsLinkInEqCons( void * user_data, int id , double* vec,
                                int len ){
  PIPS_CALLBACK_COUT << "\n[PIPS CALLBACK ENTER] FlagRhsLinkInEqCons id=" 
    << id << " len=" << len << std::endl;

  // Use the current class in the callback
  auto * solver = static_cast< PIPSMILPSolver * >( user_data );

  solver->ExtractRhsActiveFlag( id , vec, len , solver->LinkInEqCons , 
                solver->n_LinkInEqCons , solver-> rhs ,
                solver->sense , solver->rngval );
  return 0;
}

/*--------------------------------------------------------------------------*/

int PIPSMILPSolver::FlagLhsLinkInEqCons( void * user_data, int id , double* vec,
                                int len ){
  PIPS_CALLBACK_COUT << "\n[PIPS CALLBACK ENTER] FlagLhsLinkInEqCons id=" 
    << id << " len=" << len << std::endl;

  // Use the current class in the callback
  auto * solver = static_cast< PIPSMILPSolver * >( user_data );

  solver->ExtractLhsActiveFlag( id , vec, len , solver->LinkInEqCons , 
                solver->n_LinkInEqCons , solver-> rhs ,
                solver->sense , solver->rngval );
  return 0;
}

/*--------------------------------------------------------------------------*/

int PIPSMILPSolver::UBVars( void * user_data, int id , double* vec,
                                int len ){
  PIPS_CALLBACK_COUT << "\n[PIPS CALLBACK ENTER] UBVars id=" << id << " len=" 
    << len << std::endl;

  // Use the current class in the callback
  auto * solver = static_cast< PIPSMILPSolver * >( user_data );

  solver->ExtractVarBounds( id , vec, len , solver->varNode[ id ] , 
                solver->n_varNode[ id ] , solver->ub );
  return 0;
}

/*--------------------------------------------------------------------------*/

int PIPSMILPSolver::LBVars( void * user_data, int id , double* vec,
                                int len ){
  PIPS_CALLBACK_COUT << "\n[PIPS CALLBACK ENTER] LBVars id=" << id << " len=" 
    << len << std::endl;

  // Use the current class in the callback
  auto * solver = static_cast< PIPSMILPSolver * >( user_data );

  solver->ExtractVarBounds( id , vec, len , solver->varNode[ id ] , 
                solver->n_varNode[ id ] , solver->lb );
  return 0;
}

/*--------------------------------------------------------------------------*/

int PIPSMILPSolver::FlagUBVars( void * user_data, int id , double* vec,
                                int len ){
  PIPS_CALLBACK_COUT << "\n[PIPS CALLBACK ENTER] FlagUBVars id=" << id 
    << " len=" << len << std::endl;

  // Use the current class in the callback
  auto * solver = static_cast< PIPSMILPSolver * >( user_data );

  solver->ExtractFlagVarBounds( id , vec, len , solver->varNode[ id ] , 
                solver->n_varNode[ id ] , solver->ub );
  return 0;
}

/*--------------------------------------------------------------------------*/

int PIPSMILPSolver::FlagLBVars( void * user_data, int id , double* vec,
                                int len ){
  PIPS_CALLBACK_COUT << "\n[PIPS CALLBACK ENTER] FlagLBVars id=" << id 
    << " len=" << len << std::endl;

  // Use the current class in the callback
  auto * solver = static_cast< PIPSMILPSolver * >( user_data );

  solver->ExtractFlagVarBounds( id , vec, len , solver->varNode[ id ] , 
                solver->n_varNode[ id ] , solver->lb );
  return 0;
}

/*--------------------------------------------------------------------------*/
/*--------------------- End File PIPSMILPSolver.cpp ------------------------*/
/*--------------------------------------------------------------------------*/