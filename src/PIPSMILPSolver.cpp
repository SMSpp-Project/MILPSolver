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

#ifdef MILPSOLVER_DEBUG
 #define DEBUG_LOG( stuff ) std::cout << "[MILPSolver DEBUG] " << stuff
#else
 #define DEBUG_LOG( stuff )
#endif

/*--------------------------------------------------------------------------*/
/*------------------------- NAMESPACE AND USING ----------------------------*/
/*--------------------------------------------------------------------------*/

using namespace SMSpp_di_unipi_it;

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

void PIPSMILPSolver::clear_problem( unsigned int what )
{
 MILPSolver::clear_problem( 0 );

 delete pips_interface;
 pips_interface = nullptr;

 delete pips_tree;
 pips_tree = nullptr;
}

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

int PIPSMILPSolver::No_VarinNode( void * user_data, int id, int* nnz){
  // Use the current class in the callback
  auto * solver = static_cast< PIPSMILPSolver * >( user_data );
  if( id < solver->n_nodes ){
    *nnz = solver->n_varNode[ id ];
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
  if( id < solver->n_nodes ){
    *nnz = solver->n_EqConsNode[ id ];
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
  if( id < solver->n_nodes ){
    *nnz = solver->n_InEqConsNode[ id ];
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
  *nnz = solver->n_LinkEqCons;
  return 0;
}

/*--------------------------------------------------------------------------*/

int PIPSMILPSolver::No_LinkInEqCons( void * user_data, int id, int* nnz){
  // Use the current class in the callback
  auto * solver = static_cast< PIPSMILPSolver * >( user_data );
  *nnz = solver->n_LinkInEqCons;
  return 0;
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
   *nnz = sub_matrix.nnz();
   return 0;
  }
  else
   throw( std::runtime_error( "Node ID outside the expected range" ) );
}

/*--------------------------------------------------------------------------*/

int PIPSMILPSolver::nnzEqConsDiag( void * user_data, int id , int* nnz ){
  // Use the current class in the callback
  auto * solver = static_cast< PIPSMILPSolver * >( user_data );
  solver->EvaluateNnz( id , nnz , solver->EqConsNode[ id ] , 
    solver->n_EqConsNode[ id ] , solver->varNode[ id ] , 
    solver->n_varNode[ id ] );
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
  return 0;
}

/*--------------------------------------------------------------------------*/

int PIPSMILPSolver::nnzLinkInEqCons( void * user_data, int id , int* nnz ){
  // Use the current class in the callback
  auto * solver = static_cast< PIPSMILPSolver * >( user_data );
  solver->EvaluateNnz( id , nnz , solver->LinkInEqCons , 
    solver->n_LinkInEqCons , solver->varNode[ id ] , 
    solver->n_varNode[ id ] );
  return 0;
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
       sub_rhs.push_back( Inf< double >() );
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
   std::copy( sub_rhs.begin(), sub_rhs.end(), vec );
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
      sub_rhs.push_back( -Inf< double >() );
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
   std::copy( sub_rhs.begin(), sub_rhs.end(), vec );
   return 0;
  }
  else
   throw( std::runtime_error( "Node ID outside the expected range" ) );
}

/*--------------------------------------------------------------------------*/

int PIPSMILPSolver::RhsEqCons( void * user_data, int id , double* vec,
                                int len ){

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

  // Use the current class in the callback
  auto * solver = static_cast< PIPSMILPSolver * >( user_data );

  solver->ExtractRhsVector( id , vec, len , solver->LinkEqCons , 
                solver->n_LinkEqCons , solver-> rhs , {} , {} );
  return 0;
}

/*--------------------------------------------------------------------------*/

int PIPSMILPSolver::RhsLinkInEqCons( void * user_data, int id , double* vec,
                                int len ){

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

  // Use the current class in the callback
  auto * solver = static_cast< PIPSMILPSolver * >( user_data );

  solver->ExtractLhsVector( id , vec, len , solver->LinkInEqCons , 
                solver->n_LinkInEqCons , solver-> rhs ,
                solver->sense , solver->rngval );
  return 0;
}

/*--------------------------------------------------------------------------*/

int PIPSMILPSolver::ExtractRhsActiveFlag( int id , double* vec , int len ,
                                std::vector< const FRowConstraint * > node_cons ,
                                const int nCons ,
                                std::vector< double > rhs ,
                                std::vector< char > sense ,
                                std::vector< double > ranges ){ 
  ExtractRhsVector( id , vec , len , node_cons , nCons ,
                    rhs , sense , ranges );

  // Now vec stores the true rhs, convert them in flag
  for( int i = 0 ; i < len ; ++i ) 
   if( vec[ i ] == Inf< double >() )
    vec[ i ] = 0.0;
   else
    vec[ i ] = 1.0;

  return 0;
}

/*--------------------------------------------------------------------------*/

int PIPSMILPSolver::ExtractLhsActiveFlag( int id , double* vec , int len ,
                                std::vector< const FRowConstraint * > node_cons ,
                                const int nCons ,
                                std::vector< double > rhs ,
                                std::vector< char > sense ,
                                std::vector< double > ranges ){ 
  ExtractLhsVector( id , vec , len , node_cons , nCons ,
                    rhs , sense , ranges );

  // Now vec stores the true rhs, convert them in flag
  for( int i = 0 ; i < len ; ++i ) 
   if( vec[ i ] == -Inf< double >() )
    vec[ i ] = 0.0;
   else
    vec[ i ] = 1.0;

  return 0;
}

/*--------------------------------------------------------------------------*/

int PIPSMILPSolver::FlagRhsInEqCons( void * user_data, int id , double* vec,
                                int len ){

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

  // Use the current class in the callback
  auto * solver = static_cast< PIPSMILPSolver * >( user_data );

  solver->ExtractLhsActiveFlag( id , vec, len , solver->LinkInEqCons , 
                solver->n_LinkInEqCons , solver-> rhs ,
                solver->sense , solver->rngval );
  return 0;
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
    
   // We are asking the rhs for equality constraints, just select them
   for( auto idx : global_idxs_vars )
    if( idx < bounds.size() )
     sub_bounds.push_back( bounds[ idx ] );
   
   std::copy( sub_bounds.begin(), sub_bounds.end(), vec );
   return 0;
  }
  else
   throw( std::runtime_error( "Node ID outside the expected range" ) );
}

/*--------------------------------------------------------------------------*/

int PIPSMILPSolver::UBVars( void * user_data, int id , double* vec,
                                int len ){

  // Use the current class in the callback
  auto * solver = static_cast< PIPSMILPSolver * >( user_data );

  solver->ExtractVarBounds( id , vec, len , solver->varNode[ id ] , 
                solver->n_varNode[ id ] , solver->ub );
  return 0;
}

/*--------------------------------------------------------------------------*/

int PIPSMILPSolver::LBVars( void * user_data, int id , double* vec,
                                int len ){

  // Use the current class in the callback
  auto * solver = static_cast< PIPSMILPSolver * >( user_data );

  solver->ExtractVarBounds( id , vec, len , solver->varNode[ id ] , 
                solver->n_varNode[ id ] , solver->lb );
  return 0;
}

/*--------------------------------------------------------------------------*/

int PIPSMILPSolver::ObjVars( void * user_data, int id , double* vec,
                                int len ){

  // Use the current class in the callback
  auto * solver = static_cast< PIPSMILPSolver * >( user_data );

  solver->ExtractVarBounds( id , vec, len , solver->varNode[ id ] , 
                solver->n_varNode[ id ] , solver->objective );
  return 0;
}

/*--------------------------------------------------------------------------*/

int PIPSMILPSolver::ExtractFlagVarBounds( int id , double* vec , int len ,
                                std::vector< const ColVariable * > node_vars ,
                                const int nVars ,
                                std::vector< double > bounds ){ 
  ExtractVarBounds( id , vec , len , node_vars , nVars , bounds );

  // Now vec stores the true bounds, convert them in flag
  for( int i = 0 ; i < len ; ++i ) 
   if( vec[ i ] == -Inf< double >() || vec[ i ] == Inf< double >() )
    vec[ i ] = 0.0;
   else
    vec[ i ] = 1.0;

  return 0;
}

/*--------------------------------------------------------------------------*/

int PIPSMILPSolver::FlagUBVars( void * user_data, int id , double* vec,
                                int len ){

  // Use the current class in the callback
  auto * solver = static_cast< PIPSMILPSolver * >( user_data );

  solver->ExtractFlagVarBounds( id , vec, len , solver->varNode[ id ] , 
                solver->n_varNode[ id ] , solver->ub );
  return 0;
}

/*--------------------------------------------------------------------------*/

int PIPSMILPSolver::FlagLBVars( void * user_data, int id , double* vec,
                                int len ){

  // Use the current class in the callback
  auto * solver = static_cast< PIPSMILPSolver * >( user_data );

  solver->ExtractFlagVarBounds( id , vec, len , solver->varNode[ id ] , 
                solver->n_varNode[ id ] , solver->lb );
  return 0;
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
   std::copy( sub_matrix.krow.begin(), sub_matrix.krow.end(), krowM );
   std::copy( sub_matrix.jcol.begin(), sub_matrix.jcol.end(), jcolM );
   std::copy( sub_matrix.val.begin(),  sub_matrix.val.end(),  M );
   return 0;
  }
  else
   throw( std::runtime_error( "Node ID outside the expected range" ) );
}

/*--------------------------------------------------------------------------*/

int PIPSMILPSolver::MatEqConsDiag( void * user_data, int id , int* krowM, 
                                    int* jcolM, double* M ){

  // Use the current class in the callback
  auto * solver = static_cast< PIPSMILPSolver * >( user_data );

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

  solver->ExtractMatrix( id , krowM, jcolM , M , solver->LinkInEqCons , 
                solver->n_LinkInEqCons , solver->varNode[ id ] , 
                solver->n_varNode[ id ] );
  return 0;
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

 return result;
}

/*--------------------------------------------------------------------------*/

int PIPSMILPSolver::matAllZero(void*, int, int*, int*, double*) {
   return 0;
}

/*--------------------------------------------------------------------------*/

int PIPSMILPSolver::nnzAllZero(void*, int, int* nnz) {
   *nnz = 0;
   return 0;
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

 nodes_subtrees.clear();

 // Set the first node (root) to be the f_Block
 n_nodes = 1;
 Index n_blocks = 1; 
 nodes_subtrees.push_back( { f_Block } );

 // Collect all the Blocks representing the leaves of the PIPS tree
 for( auto leaf : f_Block->get_nested_Blocks() ){
  // Collect the subtree data
  std::vector< Block * > subtree;

  // Add the leaf
  subtree.push_back( leaf );

  // recursively collect all the sub-tree;
  n_blocks += collect_subtree( leaf , subtree , n_nodes );
  subtree.shrink_to_fit();

  // Add to the global structure
  nodes_subtrees.push_back( std::move( subtree ) );
  n_nodes++;
 }
 nodes_subtrees.shrink_to_fit();

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

 // First scan: scan all the variables and assign them to the correct node
 DEBUG_LOG( "First Scan: storing variables per node" << std::endl );

 Index num_node = 0;
 for( auto node : nodes_subtrees ){
  for( auto qb : node ){
    for( const auto & i : qb->get_static_variables() ){
     // Call specific function to scan the new group of Variables
     scan_group( i , qb , num_node , true , un_any_type< ColVariable >() );
    }

    for( const auto & i : qb->get_dynamic_variables() ){
     // Call specific function to scan the new group of Variables
     scan_group( i , qb , num_node , false , un_any_type< ColVariable >() );
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
     scan_group( i , qb , num_node , true , un_any_type< FRowConstraint >() );
    }

    for( const auto & i : qb->get_dynamic_constraints() ){
     // Call specific function to scan the new group of constraints
     scan_group( i , qb , num_node , false , un_any_type< FRowConstraint >() );
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
 FNNZ fNo_EqConsinNode = &No_EqConsinNode; // Number of equality constraints per node
 FNNZ fNo_InEqConsinNode = &No_InEqConsinNode; // Number of inequality constraints per node
 FNNZ fNo_LinkEqCons = &No_LinkEqCons; // Number of global linking equality constraints
 FNNZ fNo_LinkInEqCons = &No_LinkInEqCons; // Number of global linking inequality constraints

 // Nonzeros
 FNNZ fnnzQ = &nnzAllZero; // TBD: number of quadratic nonzero terms in objective
 
 FNNZ fnnzEqConsDiag = &nnzEqConsDiag; // Number of nonzeros in the diagonal matrices
 FNNZ fnnzEqConsVert = &nnzEqConsVert; // Number of nonzeros in the vertical root matrices

 FNNZ fnnzInEqConsDiag = &nnzInEqConsDiag; // Number of nonzeros in the equality diagonal matrices
 FNNZ fnnzInEqConsVert = &nnzInEqConsVert; // Number of nonzeros in the inequality vertical root matrices

 FNNZ fnnzLinkEqCons = &nnzLinkEqCons; // Number of nonzeros in the equality linking matrices
 FNNZ fnnzLinkInEqCons = &nnzLinkInEqCons; // Number of nonzeros in the inequality linking matrices
   
 // Vectors (TBD)
 FVEC fRhsEqCons = &RhsEqCons; // rhs of equality constraints per node
 FVEC fRhsInEqCons = &RhsInEqCons; // rhs of linking inequality constraints per node
 FVEC fLhsInEqCons = &LhsInEqCons; // lhs of linking inequality constraints per node

 FVEC fRhsLinkEqCons = &RhsLinkEqCons; // rhs of linking equality constraints
 FVEC fRhsLinkInEqCons = &RhsLinkInEqCons; // rhs of linking inequality constraints
 FVEC fLhsLinkInEqCons = &LhsLinkInEqCons; // lhs of linking inequality constraints

 FVEC fFlagRhsInEqCons = &FlagRhsInEqCons; // active rhs of linking inequality constraints per node
 FVEC fFlagLhsInEqCons = &FlagLhsInEqCons; // active lhs of linking inequality constraints per node

 FVEC fFlagRhsLinkInEqCons = &FlagRhsLinkInEqCons; // active rhs of linking inequality constraints
 FVEC fFlagLhsLinkInEqCons = &FlagLhsLinkInEqCons; // active lhs of linking inequality constraints

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
 std::unique_ptr<DistributedInputTree::DistributedInputNode> data_root = 
  std::make_unique<DistributedInputTree::DistributedInputNode>( this, 0, fNo_VarinNode, 
    fNo_EqConsinNode, fNo_LinkEqCons, fNo_InEqConsinNode, fNo_LinkInEqCons, fQ, fnnzQ, fObjVars, 
    fMatEqConsDiag, fnnzEqConsDiag, fMatEqConsVert , fnnzEqConsVert , fMatLinkEqCons , fnnzLinkEqCons, 
    fRhsEqCons , fRhsLinkEqCons , fMatInEqConsDiag , fnnzInEqConsDiag , fMatInEqConsVert , 
    fnnzInEqConsVert , fMatLinkInEqCons , fnnzLinkInEqCons , fLhsInEqCons , fFlagLhsInEqCons ,
    fRhsInEqCons , fFlagRhsInEqCons , fLhsLinkInEqCons , fFlagLhsLinkInEqCons , fRhsLinkInEqCons ,
    fFlagRhsLinkInEqCons , fLBVars , fFlagLBVars , fUBVars , fFlagUBVars , nullptr, nullptr, false );

 auto* root = new DistributedInputTree( std::move( data_root ) );

 for(int id = 1; id < n_nodes ; ++id ) {
  // Build the problem tree
  std::unique_ptr<DistributedInputTree::DistributedInputNode> data_child = 
    std::make_unique<DistributedInputTree::DistributedInputNode>( this, id, fNo_VarinNode, 
      fNo_EqConsinNode, fNo_LinkEqCons, fNo_InEqConsinNode, fNo_LinkInEqCons, fQ, fnnzQ, fObjVars, 
      fMatEqConsDiag, fnnzEqConsDiag, fMatEqConsVert , fnnzEqConsVert , fMatLinkEqCons , fnnzLinkEqCons, 
      fRhsEqCons , fRhsLinkEqCons , fMatInEqConsDiag , fnnzInEqConsDiag , fMatInEqConsVert , 
      fnnzInEqConsVert , fMatLinkInEqCons , fnnzLinkInEqCons , fLhsInEqCons , fFlagLhsInEqCons ,
      fRhsInEqCons , fFlagRhsInEqCons , fLhsLinkInEqCons , fFlagLhsLinkInEqCons , fRhsLinkInEqCons ,
      fFlagRhsLinkInEqCons , fLBVars , fFlagLBVars , fUBVars , fFlagUBVars , nullptr, nullptr, false );

   root->add_child( std::make_unique<DistributedInputTree>( std::move( data_child ) ) );
 }

 if (rank == 0)
  std::cout << "Using a total of " << size << " MPI processes.\n";

/* use BiCGStab for outer solve */
//pipsipmpp_options::set_parameter("PRESOLVE", false);
//pipsipmpp_options::set_parameter("SCALER", "geometricmean");

 pips_tree = root;
 pips_interface = new PIPSIPMppInterface( pips_tree, MPI_COMM_WORLD );

 return;
}

/*--------------------------------------------------------------------------*/

template< typename T >
 void PIPSMILPSolver::scan_group( const boost::any & gr , Block * qb ,
                              Index num_node , bool is_static , 
                              un_any_type< T > )
{
 // Search for the group type
 if( gr.type() == typeid( T * ) ||
      gr.type() == typeid( std::vector< T > * ) ||
      gr.type() == typeid( std::vector< std::vector< T > > * ) ) {
  // "Simple" group
  scan_simple_group( gr , qb , num_node , is_static, un_any_type< T >() );
  }
 else {
  // "Complex" group
  scan_multiarray_group( gr , qb , num_node , is_static , un_any_type< T >() );
  }
 } // end( PIPSMILPSolver::scan_group )

/*--------------------------------------------------------------------------*/

template< typename T >
 void PIPSMILPSolver::scan_simple_group( const boost::any & gr , Block * qb , 
                                Index num_node , bool is_static , 
                                un_any_type< T > )
{
 if( typeid( T * ) == typeid( FRowConstraint * ) ) {
  // Scanning a group of Constraints
  auto scan = [ this , num_node ]
   ( const FRowConstraint & c ) {
    scan_constraint( c , num_node );
  };
  // Scan all the constraints one at a time
  if( is_static )
   un_any_const_static( gr , scan  , un_any_type< FRowConstraint >() );
  else
   un_any_const_dynamic( gr , scan , un_any_type< FRowConstraint >() );
 }
 else if( typeid( T * ) == typeid( ColVariable * ) ) {
  // Scanning a group of Variables
  // in this case we simply add the variable to the corresponding node
  auto push_var_toNode = [ this , num_node ]
    ( const ColVariable & c ) {
      n_varNode[ num_node ] += 1;
      varNode[ num_node ].push_back( &c );
  };
  if( is_static )
   un_any_const_static( gr , push_var_toNode , un_any_type< ColVariable >() );
  else
   un_any_const_dynamic( gr , push_var_toNode , un_any_type< ColVariable >() );
  }
 else
  throw( std::runtime_error( "Unsupported group type" ) );
} // end( PIPSMILPSolver::scan_simple_group )

/*--------------------------------------------------------------------------*/

template< typename T >
 void PIPSMILPSolver::scan_multiarray_group( const boost::any & gr ,
            Block * qb , Index num_node , bool is_static , un_any_type< T > )
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
     if( is_static )
      un_any_const_static( gr , scan  , un_any_type< FRowConstraint >() );
     else
      un_any_const_dynamic( gr , scan , un_any_type< FRowConstraint >() );

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
     };
     if( is_static )
      un_any_const_static( gr , push_var_toNode , un_any_type< ColVariable >() );
     else
      un_any_const_dynamic( gr , push_var_toNode , un_any_type< ColVariable >() );

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
     if( is_static )
      un_any_const_static( gr , scan  , un_any_type< FRowConstraint >() );
     else
      un_any_const_dynamic( gr , scan , un_any_type< FRowConstraint >() );

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
     };
     if( is_static )
      un_any_const_static( gr , push_var_toNode , un_any_type< ColVariable >() );
     else
      un_any_const_dynamic( gr , push_var_toNode , un_any_type< ColVariable >() );

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
     if( is_static )
      un_any_const_static( gr , scan  , un_any_type< FRowConstraint >() );
     else
      un_any_const_dynamic( gr , scan , un_any_type< FRowConstraint >() );

     // The linearization produced by ma->data() for the 3D multi_array
     // stores elements in row-major order.
     if( idx_2 < ma->shape()[ 1 ] - 1 )
      idx_2++; // Move third counter
     else if( idx_1 < ma->shape()[ 1 ] - 1 ) {
      idx_1++; // Move second counter
      idx_2 = 0; // Reset third counter
     }
     else {
      idx_0++; // Move first counter
      idx_1 = 0; // Reset second counter
      idx_2 = 0; // Reset third counter
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
     };
     if( is_static )
      un_any_const_static( gr , push_var_toNode , un_any_type< ColVariable >() );
     else
      un_any_const_dynamic( gr , push_var_toNode , un_any_type< ColVariable >() );

     // The linearization produced by ma->data() for the 3D multi_array
     // stores elements in row-major order.
     if( idx_2 < ma->shape()[ 1 ] - 1 )
      idx_2++; // Move third counter
     else if( idx_1 < ma->shape()[ 1 ] - 1 ) {
      idx_1++; // Move second counter
      idx_2 = 0; // Reset third counter
     }
     else {
      idx_0++; // Move first counter
      idx_1 = 0; // Reset second counter
      idx_2 = 0; // Reset third counter
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
     if( is_static )
      un_any_const_static( gr , scan  , un_any_type< FRowConstraint >() );
     else
      un_any_const_dynamic( gr , scan , un_any_type< FRowConstraint >() );

     // The linearization produced by ma->data() for the 3D multi_array
     // stores elements in row-major order.
     if( idx_2 < ma->shape()[ 1 ] - 1 )
      idx_2++; // Move third counter
     else if( idx_1 < ma->shape()[ 1 ] - 1 ) {
      idx_1++; // Move second counter
      idx_2 = 0; // Reset third counter
     }
     else {
      idx_0++; // Move first counter
      idx_1 = 0; // Reset second counter
      idx_2 = 0; // Reset third counter
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
     };
     if( is_static )
      un_any_const_static( gr , push_var_toNode , un_any_type< ColVariable >() );
     else
      un_any_const_dynamic( gr , push_var_toNode , un_any_type< ColVariable >() );

     // The linearization produced by ma->data() for the 3D multi_array
     // stores elements in row-major order.
     if( idx_2 < ma->shape()[ 1 ] - 1 )
      idx_2++; // Move third counter
     else if( idx_1 < ma->shape()[ 1 ] - 1 ) {
      idx_1++; // Move second counter
      idx_2 = 0; // Reset third counter
     }
     else {
      idx_0++; // Move first counter
      idx_1 = 0; // Reset second counter
      idx_2 = 0; // Reset third counter
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

} // end( PIPSMILPSolver::scan_multiarray_group )

/*--------------------------------------------------------------------------*/

void PIPSMILPSolver::scan_constraint( const FRowConstraint & con , Index num_node )
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

    // Understand if v belongs to the same node
    bool found_inNode = 
      std::find( varNode[ num_node ].begin() , varNode[ num_node ].end() , v ) 
        != varNode[ num_node ].end();

    // Understand if v belongs to the root
    bool found_inRoot = 
      std::find( varNode[ 0 ].begin() , varNode[ 0 ].end() , v ) 
        != varNode[ 0 ].end();

    if( ! found_inNode && ! found_inRoot ) {
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

int PIPSMILPSolver::compute( bool changedvars )
{
 lock();  // lock the mutex: this is done again inside MILPSolver::compute,
          // but that's OK since the mutex is recursive

 // process Modification: this is driven by MILPSolver- - - - - - - - - - - -
 if( MILPSolver::compute( changedvars ) != kOK )
  throw( std::runtime_error( "an error occurred in MILPSolver::compute()" ) );

 // if required, write the problem to file- - - - - - - - - - - - - - - - - -
 // Not possible in PIPS

 // the continuous case - - - - - - - - - - - - - - - - - - - - - - - - - - -
 sol_status = decode_pips_status( pips_interface->run() );

 Return_status:
 unlock();  // unlock the mutex
 return( sol_status );

}  // end( PIPSMILPSolver::compute )

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