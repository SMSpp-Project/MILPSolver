/*--------------------------------------------------------------------------*/
/*-------------------------- File MILPSolver.cpp ---------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Implementation of the MILPSolver class.
 *
 * \author Antonio Frangioni \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \author Niccolo' Iardella \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * Copyright &copy; by Antonio Frangioni, Niccolo' Iardella
 */
/*--------------------------------------------------------------------------*/
/*---------------------------- IMPLEMENTATION ------------------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*------------------------------- MACROS -----------------------------------*/
/*--------------------------------------------------------------------------*/
/* If the macro MILPSOLVER_DEBUG is externally defined, then some costly
 * checks on the data astructures of MILPSolver are performed and debug
 * information printed. Also, the method check_status() is defined and
 * used to check the whole set of data structures. */

#ifdef MILPSOLVER_DEBUG
 #define DEBUG_LOG( stuff ) std::cout << "[MILPSolver DEBUG] " << stuff
#else
 #define DEBUG_LOG( stuff )
#endif

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include <queue>

#include <LinearFunction.h>

#include <DQuadFunction.h>

#include "MILPSolver.h"

/*--------------------------------------------------------------------------*/
/*------------------------- NAMESPACE AND USING ----------------------------*/
/*--------------------------------------------------------------------------*/

using namespace SMSpp_di_unipi_it;

/*--------------------------------------------------------------------------*/
/*-------------------------- FACTORY MANAGEMENT ----------------------------*/
/*--------------------------------------------------------------------------*/

SMSpp_insert_in_factory_cpp_0( MILPSolver );

/*--------------------------------------------------------------------------*/
/*----------------------------- FUNCTIONS ----------------------------------*/
/*--------------------------------------------------------------------------*/

template< typename T >
std::string log_vector( const std::vector< T > & v, int limit = 10 );

template<>
std::string log_vector( const std::vector< char > & v, int limit );

/*--------------------------------------------------------------------------*/
/*-------------------------------- SET_BLOCK -------------------------------*/
/*--------------------------------------------------------------------------*/

void MILPSolver::set_Block( Block * block )
{
 if( f_Block == block )  // registering to the same Block
  return;                // cowardly and silently return

 Solver::set_Block( block );

 if( block ) {
  bool owned = block->is_owned_by( f_id );
  if( ( ! owned ) && ( ! block->lock( f_id ) ) )
   throw( std::runtime_error( "Unable to lock the Block" ) );

  // Generate abstract representation
  block->generate_abstract_variables( nullptr );
  block->generate_abstract_constraints( nullptr );
  block->generate_objective( nullptr );

  if( ! owned )
   block->unlock( f_id );

  load_problem();
  }
 }

/*--------------------------------------------------------------------------*/
/*------------------------------- CLEAR/LOAD -------------------------------*/
/*--------------------------------------------------------------------------*/

void MILPSolver::clear_problem( unsigned int what )
{
 if( what & 1u ) {
  matbeg.clear();
  matcnt.clear();
  matind.clear();
  matval.clear();
  xctype.clear();

  for( auto & i: colname )
   delete[] i;
  for( auto & i: rowname )
   delete[] i;
  colname.clear();
  rowname.clear();
 }

 if( what & 2u ) {
  objective.clear();
  q_objective.clear();
 }

 if( what & 4u ) {
  sense.clear();
  rhs.clear();
  rngval.clear();
 }

 if( what & 8u ) {
  lb.clear();
  ub.clear();
 }
}

/*--------------------------------------------------------------------------*/

void MILPSolver::load_problem( void )
{
 numrows = 0;
 numcols = 0;
 static_vars = 0;
 static_cons = 0;
 nzelements = 0;

 std::queue< Block * > Q;

 // Locking the Block
 bool owned = f_Block->is_owned_by( f_id );
 if( ( ! owned ) && ( ! f_Block->read_lock() ) )
  throw( std::runtime_error( "Unable to lock the Block" ) );

 int num_block = 0;       // Counter for the blocks
 int row = 0;             // Counter for the rows
 int col = 0;             // Counter for the columns
 int static_con_grps = 0; // Counter for static constraint groups
 int static_var_grps = 0; // Counter for static variable groups

 // Count variables and constraints
 // --------------------------------------------------------------------------

 Q.push( f_Block );
 while( ! Q.empty() ) {
  Block * q_Block = Q.front();
  Q.pop();
  DEBUG_LOG( "Processing Block " << num_block << " ["
	                         << q_Block << "]" << std::endl );
  DEBUG_LOG( *q_Block << std::endl );
 
  for( auto * i : q_Block->get_nested_Blocks() )
   Q.push( i );

  for( const auto & i : q_Block->get_static_constraints() ) {
   // Singles
   if( un_any_thing_0( FRowConstraint , i ,
                       {
                        ++numrows;
                        ++static_cons;
                        ++static_con_grps;
                        }
		       ) )
    continue;

   // Vectors
   if( un_any_thing_1( FRowConstraint , i ,
                       {
                        numrows += var.size();
                        static_cons += var.size();
                        ++static_con_grps;
                        }
		       ) )
    continue;

   // Multiarrays
   if( un_any_thing_K( FRowConstraint, i,
                       {
                        numrows += var.num_elements();
                        static_cons += var.num_elements();
                        ++static_con_grps;
                        }
		       ) )
    continue;
   }

  for( const auto & i : q_Block->get_dynamic_constraints() ) {
   // Single lists
   if( un_any_thing_0( std::list< FRowConstraint > , i ,
                       { numrows += var.size(); } ) )
    continue;

   // Vectors of lists
   if( un_any_thing_1( std::list< FRowConstraint > , i ,
                       {
                        for( auto & el: var )
                         numrows += el.size();
                        }
		       ) )
    continue;

   // Multiarrays of lists
   if( un_any_thing_K( std::list< FRowConstraint > , i ,
                       {
                        auto it = var.data();
                        for( auto i = var.num_elements() ; i-- ; ++it )
                         numrows += it->size();
		        }
		       ) )
    continue;
   }

  for( const auto & i : q_Block->get_static_variables() ) {
   // Singles
   if( un_any_thing_0( ColVariable , i ,
                       {
                        ++numcols;
                        ++static_vars;
                        ++static_var_grps;
                        }
		       ) )
    continue;

   // Vectors
   if( un_any_thing_1( ColVariable , i ,
                       {
                        numcols += var.size();
                        static_vars += var.size();
                        ++static_var_grps;
                        }
		       ) )
    continue;

   // Multiarrays
   if( un_any_thing_K( ColVariable , i ,
                       {
                        numcols += var.num_elements();
                        static_vars += var.num_elements();
                        ++static_var_grps;
                        }
		       ) )
    continue;
   }

  for( const auto & i : q_Block->get_dynamic_variables() ) {
   // Single lists
   if( un_any_thing_0( std::list< ColVariable > , i ,
                       { numcols += var.size(); } ) )
    continue;

   // Vectors of lists
   if( un_any_thing_1( std::list< ColVariable > , i ,
                       {
                        for( auto & el: var )
                         numcols += el.size();
                        }
		       ) )
    continue;

   // Multiarrays of lists
   if( un_any_thing_K( std::list< ColVariable > , i ,
                       {
                        auto it = var.data();
                        for( auto i = var.num_elements() ; i-- ; ++it )
                         numcols += it->size();
                        }
		       ) )
    continue;
   }

  auto counter = [ this ]( ColVariable & var ) {
   for( auto * i : var.active_stuff() ) {
    auto * row = dynamic_cast< FRowConstraint * >( i );
    if( ! row )
     continue;

    bool is_mine = false;
    auto b = row->get_Block();
    do {
     if( b == f_Block ) {
      is_mine = true;
      break;
      }
     else
      b = b->get_f_Block();
     } while( b );

    if( is_mine )
     ++nzelements;
    }
   };

  for( const auto & i : q_Block->get_static_variables() )
   un_any_const_static( i , counter , un_any_type< ColVariable >() );

  for( const auto & i : q_Block->get_dynamic_variables() )
   un_any_const_dynamic( i , counter , un_any_type< ColVariable >() );

  ++num_block;
  }

 DEBUG_LOG( "Number of blocks      = " << num_block << std::endl );
 DEBUG_LOG( "numrows (constraints) = " << numrows << " (S:" << static_cons
	    << "/D:" << numrows - static_cons << ")" << std::endl );
 DEBUG_LOG( "numcols (variables)   = " << numcols << " (S:" << static_vars
	    << "/D:" << numcols - static_vars << ")" << std::endl );
 DEBUG_LOG( "nzelements            = " << nzelements << std::endl );

 // LP vector allocation
 // --------------------------------------------------------------------------

 // The +1 is needed by generic interface
 matbeg.resize( numcols + 1, 0 );
 matbeg[ numcols ] = nzelements;

 matcnt.resize( numcols, 0 );
 matind.resize( nzelements, 0 );
 matval.resize( nzelements, 0 );
 rhs.resize( numrows, 0 );
 rngval.resize( numrows, 0 );
 sense.resize( numrows, 0 );
 objective.resize( numcols, 0 );
 q_objective.resize( numcols, 0 );
 lb.resize( numcols, 0 );
 ub.resize( numcols, 0 );
 xctype.resize( numcols, 0 );
 colname.resize( numcols , nullptr );
 rowname.resize( numrows , nullptr );

 svar_to_idx.clear();
 idx_to_svar.clear();
 scon_to_idx.clear();
 idx_to_scon.clear();
 dvar_to_idx.clear();
 idx_to_dvar.clear();
 dcon_to_idx.clear();
 idx_to_dcon.clear();

 svar_to_idx.reserve( static_var_grps );
 idx_to_svar.reserve( static_var_grps );
 dvar_to_idx.reserve( numcols - static_vars );
 idx_to_dvar.reserve( numcols - static_vars );
 scon_to_idx.reserve( static_con_grps );
 idx_to_scon.reserve( static_con_grps );
 dcon_to_idx.reserve( numrows - static_cons );
 idx_to_dcon.reserve( numrows - static_cons );

 // Scan the static constraints
 // --------------------------------------------------------------------------

 num_block = 0;
 Q.push( f_Block );
 while( ! Q.empty() ) {
  Block * q_Block = Q.front();
  Q.pop();

  for( auto i : q_Block->get_nested_Blocks() )
   Q.push( i );

  int set = 0; // Counter for the constraint groups

  for( const auto & i : q_Block->get_static_constraints() ) {
   int elements = 0; // Counter for group elements
   int start = row;

   auto scan = [ this, &elements, &row ]( FRowConstraint & c ) {
    scan_constraint( c, elements, row );
    };
   un_any_const_static( i , scan , un_any_type< FRowConstraint >() );

   // Write names
   auto base = q_Block->get_s_const_name()[ set ];
   int end = row - start;
   for( int n = 0 ; n < end ; ++n ) {
    std::string name;
    if( base.empty() )
     name = "cs_" + std::to_string( num_block )
          + "_" + std::to_string( set ) + "_" + std::to_string( n );
    else
     name = base + "_" + std::to_string( num_block )
          + "_" + std::to_string( n );

    rowname[ start + n ] = strcpy( new char[ name.length() + 1 ] ,
				   name.c_str() );
    }
   set++;
   if( elements )
    std::get< 2 >( scon_to_idx.back() ) = elements;
   }
  num_block++;
  }

 std::sort( scon_to_idx.begin() , scon_to_idx.end() );

 // Scan the dynamic constraints
 // --------------------------------------------------------------------------

 num_block = 0;
 Q.push( f_Block );
 while( ! Q.empty() ) {
  Block * q_Block = Q.front();
  Q.pop();

  for( auto i : q_Block->get_nested_Blocks() )
   Q.push( i );

  int set = 0; // Counter for the constraint groups

  for( const auto & i : q_Block->get_dynamic_constraints() ) {
   int elements = -1;
   int start = row;

   auto scan = [ this , &elements , &row ]( FRowConstraint & c ) {
    scan_constraint( c , elements , row );
    };
   un_any_const_dynamic( i , scan , un_any_type< FRowConstraint >() );

   // Write names
   auto base = q_Block->get_d_const_name()[ set ];
   int end = row - start;
   for( int n = 0 ; n < end ; ++n ) {
    std::string name;
    if( base.empty() )
     name = "cd_" + std::to_string( num_block )
          + "_" + std::to_string( set ) + "_" + std::to_string( n );
    else
     name = base + "_" + std::to_string( num_block )
          + "_" + std::to_string( n );

    rowname[ start + n ] = strcpy( new char[ name.length() + 1 ] ,
				   name.c_str() );
    }
   set++;
   }
  num_block++;
  }

 std::sort( dcon_to_idx.begin() , dcon_to_idx.end() );

 // Scan the static variables
 // --------------------------------------------------------------------------

 num_block = 0;
 Q.push( f_Block );

 while( ! Q.empty() ) {
  Block * q_Block = Q.front();
  Q.pop();

  for( auto * i : q_Block->get_nested_Blocks() )
   Q.push( i );

  int set = 0;   // Counter for the variable groups

  for( const auto & i : q_Block->get_static_variables() ) {
   int elements = 0; // Counter for group elements
   int start = col;

   auto scan = [ this , &elements , &col ]( ColVariable & v ) {
    scan_variable( v , elements , col );
    };
   un_any_const_static( i , scan , un_any_type< ColVariable >() );

   // Write names
   auto base = q_Block->get_s_var_name()[ set ];
   int end = col - start;
   for( int n = 0 ; n < end ; ++n ) {
    std::string name;
    if( base.empty() )
     name = "xs_" + std::to_string( num_block )
          + "_" + std::to_string( set ) + "_" + std::to_string( n );
    else
     name = base + "_" + std::to_string( num_block )
          + "_" + std::to_string( n );

    colname[ start + n ] = strcpy( new char[ name.length() + 1 ] ,
				   name.c_str() );
    }
   set++;
   if( elements )
    std::get< 2 >( svar_to_idx.back() ) = elements;
   }
  num_block++;
  }

 std::sort( svar_to_idx.begin() , svar_to_idx.end() );

 // Scan the dynamic variables
 // --------------------------------------------------------------------------

 num_block = 0;
 Q.push( f_Block );
 while( ! Q.empty() ) {
  Block * q_Block = Q.front();
  Q.pop();

  for( auto * i : q_Block->get_nested_Blocks() )
   Q.push( i );

  int set = 0;   // Counter for the variable groups

  for( const auto & i : q_Block->get_dynamic_variables() ) {
   int elements = -1;
   int start = col;

   auto scan = [ this , &elements , &col ]( ColVariable & v ) {
    scan_variable( v, elements, col );
    };
   un_any_const_dynamic( i , scan , un_any_type< ColVariable >() );

   // Write names
   auto base = q_Block->get_d_var_name()[ set ];
   int end = col - start;
   for( int n = 0 ; n < end ; ++n ) {
    std::string name;
    if( base.empty() )
     name = "xv_" + std::to_string( num_block )
            + "_" + std::to_string( set ) + "_" + std::to_string( n );
    else
     name = base + "_" + std::to_string( num_block )
          + "_" + std::to_string( n );

    colname[ start + n ] = strcpy( new char[ name.length() + 1 ] ,
				   name.c_str() );
    }
   set++;
  }
  num_block++;
  }

 std::sort( dvar_to_idx.begin() , dvar_to_idx.end() );

 // Scan the objective
 // --------------------------------------------------------------------------

 switch( f_Block->get_objective_sense() ) {
  case( Objective::eMax ): objsense = -1; break;
  case( Objective::eMin ): objsense = 1;  break;
  default:                 objsense = 0;
  }

 Q.push( f_Block );
 while( ! Q.empty() ) {
  Block * q_Block = Q.front();
  Q.pop();

  for( auto * i : q_Block->get_nested_Blocks() )
   Q.push( i );

  auto * obj = dynamic_cast< FRealObjective * >( q_Block->get_objective() );
  if( obj )
   scan_objective( obj );
  }

 DEBUG_LOG( "objective   = " << log_vector( objective ) << std::endl );
 DEBUG_LOG( "q_objective = " << log_vector( q_objective ) << std::endl );
 DEBUG_LOG( "rhs         = " << log_vector( rhs ) << std::endl );
 DEBUG_LOG( "rngval      = " << log_vector( rngval ) << std::endl );
 DEBUG_LOG( "sense       = " << log_vector( sense ) << std::endl );
 DEBUG_LOG( "matbeg      = " << log_vector( matbeg ) << std::endl );
 DEBUG_LOG( "matcnt      = " << log_vector( matcnt ) << std::endl );
 DEBUG_LOG( "matind      = " << log_vector( matind ) << std::endl );
 DEBUG_LOG( "matval      = " << log_vector( matval ) << std::endl );
 DEBUG_LOG( "lb          = " << log_vector( lb ) << std::endl );
 DEBUG_LOG( "ub          = " << log_vector( ub ) << std::endl );
 DEBUG_LOG( "xctype      = " << log_vector( xctype ) << std::endl );

 // Unlock the Block
 if( ! owned )
  f_Block->read_unlock();
 }

/*--------------------------------------------------------------------------*/

double MILPSolver::get_problem_lb( const ColVariable & var )
{
 double b = var.get_lb();

 for( auto * i : var.active_stuff() )
  if( auto box = dynamic_cast< OneVarConstraint * >( i ) )
   b = std::max( b , box->get_lhs() );

 return( b );
 }

/*--------------------------------------------------------------------------*/

double MILPSolver::get_problem_ub( const ColVariable & var )
{
 double b = var.get_ub();

 for( auto * i : var.active_stuff() )
  if( auto box = dynamic_cast< OneVarConstraint * >( i ) )
   b = std::min( b , box->get_rhs() );

 return( b );
 }

/*--------------------------------------------------------------------------*/

std::vector< FRowConstraint * >
MILPSolver::get_active_constraints( const ColVariable & var )
{
 std::vector< FRowConstraint * > active_constraints;
 for( auto * i : var.active_stuff() )
  if( auto * row = dynamic_cast<FRowConstraint *>( i ) ) {
   bool is_mine = false;
   auto b = row->get_Block();
   do {
    if( b == f_Block ) {
     is_mine = true;
     break;
     }
    else
     b = b->get_f_Block();

    } while( b );

   if( is_mine )
    active_constraints.push_back( row );
   }
 
 return( active_constraints );
 }

/*--------------------------------------------------------------------------*/

std::vector< OneVarConstraint * >
MILPSolver::get_active_bounds( const ColVariable & var )
{
 std::vector< OneVarConstraint * > active_bounds;
 for( auto * i : var.active_stuff() )
  if( auto row = dynamic_cast< OneVarConstraint * >( i ) )
   active_bounds.push_back( row );

 return( active_bounds );
 }

/*--------------------------------------------------------------------------*/
/*-------------------- METHODS FOR PROBLEM DESCRIPTION ---------------------*/
/*--------------------------------------------------------------------------*/

int MILPSolver::index_of_variable( const ColVariable * var )
{
 auto i = index_of_static_variable( var );
 return( i < Inf< int >() ? i : index_of_dynamic_variable( var ) );
 }

/*--------------------------------------------------------------------------*/

int MILPSolver::index_of_static_variable( const ColVariable * var )
{
 if( svar_to_idx.empty() )
  return( Inf< int >() );

 #ifdef MILPSOLVER_DEBUG
  assert( std::is_sorted( svar_to_idx.begin(), svar_to_idx.end() ) );
 #endif
 auto it = upper_bound( svar_to_idx.begin(), svar_to_idx.end(),
                        std::make_tuple( var, 0, 0 ),
                        [ & ]( auto & p1, auto & p2 ) {
                         return std::get< 0 >( p1 ) < std::get< 0 >( p2 );
                        } );

 // Now it refers to the first (group of) element(s) greater than var
 if( it == svar_to_idx.begin() )
  return( Inf< int >() );

 --it;

 // First element of the variable group
 auto first = const_cast< const ColVariable * >( std::get< 0 >( *it ) );
 int distance = std::distance( first , var );

 if( ( distance >= 0 ) && ( distance < std::get< 2 >( *it ) ) )
  // The element belongs to this group
  return( std::get< 1 >( *it ) + distance );

 // The element doesn't exist
 return( Inf< int >() );
 }

/*--------------------------------------------------------------------------*/

int MILPSolver::index_of_dynamic_variable( const ColVariable * var )
{
 #ifdef MILPSOLVER_DEBUG
  assert( std::is_sorted( dvar_to_idx.begin(), dvar_to_idx.end() ) );
 #endif
 auto it = lower_bound( dvar_to_idx.begin(), dvar_to_idx.end(),
                        std::make_pair( var, 0 ),
                        [ & ]( auto & p1, auto & p2 ) {
                         return p1.first < p2.first;
                        } );

 if( it != dvar_to_idx.end() && it->first == var )
  return( it->second );

 return( Inf< int >() );
 }

/*--------------------------------------------------------------------------*/

int MILPSolver::index_of_constraint( const FRowConstraint * con )
{
 auto i = index_of_static_constraint( con );
 return( i < Inf< int >() ? i : index_of_dynamic_constraint( con ) );
 }

/*--------------------------------------------------------------------------*/

int MILPSolver::index_of_static_constraint( const FRowConstraint * con )
{
 if( scon_to_idx.empty() )
  return( Inf< int >() );

 #ifdef MILPSOLVER_DEBUG
  assert( std::is_sorted( scon_to_idx.begin(), scon_to_idx.end() ) );
 #endif
 auto it = upper_bound( scon_to_idx.begin(), scon_to_idx.end(),
                        std::make_tuple( con, 0, 0 ),
                        [ & ]( auto & p1, auto & p2 ) {
                         return std::get< 0 >( p1 ) < std::get< 0 >( p2 );
                        } );

 // Now it refers to the first (group of) element(s) greater than var
 if( it == scon_to_idx.begin() )
  return( Inf< int >() );

 --it;

 // First element of the constraint group
 auto first = const_cast<const FRowConstraint *>(std::get< 0 >( *it ));
 int distance = std::distance( first , con );

 if( ( distance >= 0 ) && ( distance < std::get< 2 >( *it ) ) )
  // The element belongs to this group
  return( std::get< 1 >( *it ) + distance );

 // The element doesn't exist
 return( Inf< int >() );
 }

/*--------------------------------------------------------------------------*/

int MILPSolver::index_of_dynamic_constraint( const FRowConstraint * con )
{
 #ifdef MILPSOLVER_DEBUG
  assert( std::is_sorted( dcon_to_idx.begin(), dcon_to_idx.end() ) );
 #endif
 auto it = lower_bound( dcon_to_idx.begin(), dcon_to_idx.end(),
                        std::make_pair( con, 0 ),
                        [ & ]( auto & p1, auto & p2 ) {
                         return p1.first < p2.first;
                        } );

 if( it != dcon_to_idx.end() && it->first == con )
  return( it->second );

 return( Inf< int >() );
 }

/*--------------------------------------------------------------------------*/

ColVariable * MILPSolver::variable_with_index( int i )
{
 if( i < static_vars )
  return( static_variable_with_index( i ) );

 return( dynamic_variable_with_index( i ) );
 }

/*--------------------------------------------------------------------------*/

ColVariable * MILPSolver::static_variable_with_index( int i )
{
 if( idx_to_svar.empty() )
  return( nullptr );

 #ifdef MILPSOLVER_DEBUG
  assert( std::is_sorted( idx_to_svar.begin(), idx_to_svar.end() ) );
 #endif
 auto it = upper_bound( idx_to_svar.begin(), idx_to_svar.end(),
                        std::make_pair( i, nullptr ),
                        [ & ]( auto & p1, auto & p2 ) {
                         return p1.first < p2.first;
                        } );

 // Now it refers to the first (group of) element(s) greater than i
 if( it == idx_to_svar.begin() )
  return( nullptr );

 --it;

 int distance = i - it->first;
 return( it->second + distance );
 }

/*--------------------------------------------------------------------------*/

ColVariable * MILPSolver::dynamic_variable_with_index( int i )
{
 if( ( static_vars < i ) || ( i < numcols ) )
  return( idx_to_dvar[ i - static_vars ] );

 return( nullptr );
 }

/*--------------------------------------------------------------------------*/

FRowConstraint * MILPSolver::constraint_with_index( int i )
{
 if( i < static_cons )
  return( static_constraint_with_index( i ) );

 return( dynamic_constraint_with_index( i ) );
 }

/*--------------------------------------------------------------------------*/

FRowConstraint * MILPSolver::static_constraint_with_index( int i )
{
 if( idx_to_scon.empty() )
  return( nullptr );

 #ifdef MILPSOLVER_DEBUG
  assert( std::is_sorted( idx_to_scon.begin(), idx_to_scon.end() ) );
 #endif
 auto it = upper_bound( idx_to_scon.begin(), idx_to_scon.end(),
                        std::make_pair( i, nullptr ),
                        [ & ]( auto & p1, auto & p2 ) {
                         return p1.first < p2.first;
                        } );

 // Now it refers to the first (group of) element(s) greater than i
 if( it == idx_to_scon.begin() )
  return( nullptr );

 --it;

 int distance = i - it->first;
 return( it->second + distance );
 }

/*--------------------------------------------------------------------------*/

FRowConstraint * MILPSolver::dynamic_constraint_with_index( int i )
{
 if( ( static_cons < i ) || ( i < numrows ) )
  return( idx_to_dcon[ i - static_cons ] );

 return( nullptr );
 }

/*--------------------------------------------------------------------------*/
/*------------- AUXILIARY METHODS FOR POPULATING THE PROBLEM  --------------*/
/*--------------------------------------------------------------------------*/

void MILPSolver::scan_variable( ColVariable & var , int & n , int & col )
{
 if( n < 0 ) {
  // DEBUG_LOG( "MILPSolver::scan_variable(): D#" << col << " " << var );
  dvar_to_idx.emplace_back( &var, col );
  idx_to_dvar.emplace_back( &var );
  }
 else {
  // DEBUG_LOG( "MILPSolver::scan_variable(): S#"
  //             << n << "/" << col << " " << var );

  if( n == 0 ) {  // The tuple's third field will be filled later
   svar_to_idx.emplace_back( &var, col, 0 );
   idx_to_svar.emplace_back( col, &var );
   }
  ++n;
  }

 if( var.is_fixed() ) {
  lb[ col ] = std::max( get_problem_lb( var ) , var.get_value() );
  ub[ col ] = std::min( get_problem_ub( var ) , var.get_value() );
  }
 else {
  lb[ col ] = get_problem_lb( var );
  ub[ col ] = get_problem_ub( var );
  }

 if( var.is_integer() && ( ! relax_int_vars ) ) {
  ++int_vars;
  if( var.is_unitary() && var.is_positive() )
   xctype[ col ] = 'B'; // Binary
  else
   xctype[ col ] = 'I'; // Integer
  }
 else
  xctype[ col ] = 'C';  // Continuous

 /* From CPLEX documentation about matval, matbeg, matcnt, and matind:
  *
  * CPLEX needs to know only the nonzero coefficients.
  * These are grouped by column in the array matval.
  * The nonzero elements of every column must be stored in sequential
  * locations in this array with matbeg[ j ] containing the index of the
  * beginning of column j and matcnt[ j ] containing the number of entries
  * in column j. The components of matbeg must be in ascending order.
  * For each k, matind[ k ] specifies the row number of the
  * corresponding coefficient, matval[ k ]. */

 auto active_constraints = get_active_constraints( var );
 int nz_elements = active_constraints.size();
 matcnt[ col ] = nz_elements;

 if( col == 0 )
  matbeg[ col ] = 0;
 else
  matbeg[ col ] = matbeg[ col - 1 ] + matcnt[ col - 1 ];

 for( int j = 0; j < nz_elements; ++j ) {
  auto * con = active_constraints[ j ];
  const auto * f = dynamic_cast<const LinearFunction *> ( con->get_function()
							  );
  // TODO: Is dynamic_cast necessary?
  if( ! f )  // TODO: Throw exception?
   continue;

  auto it = std::find_if( f->get_v_var().begin() , f->get_v_var().end() ,
                          [ & ]( LinearFunction::coeff_pair pair ) {
                           return( pair.first == &var );
                           } );

  if( it != f->get_v_var().end() ) {
   matval[ matbeg[ col ] + j ] = it->second;
   matind[ matbeg[ col ] + j ] = index_of_constraint( con );
   }
  else
   // This should never happen since we are looping on the active contraints
   throw( std::invalid_argument(
	"This ColVariable is not active in the examined FRowConstraint" ) );
  }
 ++col;
 }

/*--------------------------------------------------------------------------*/

void MILPSolver::scan_constraint( FRowConstraint & con , int & n ,
				  int & row )
{
 if( n < 0 ) {
  // DEBUG_LOG( "MILPSolver::scan_constraint(): D#" << row << " " << con );
  dcon_to_idx.emplace_back( &con, row );
  idx_to_dcon.emplace_back( &con );
  }
 else {
  // DEBUG_LOG( "MILPSolver::scan_constraint(): S#"
  //             << n << "/" << row << " " << con );

  if( n == 0 ) {
   // The tuple's third field will be filled later
   scon_to_idx.emplace_back( &con, row, 0 );
   idx_to_scon.emplace_back( row, &con );
   }
  ++n;
  }

 const auto * lf = dynamic_cast< const LinearFunction * >(
						     con.get_function() );
 if( ! lf )
  throw( std::invalid_argument( "The Constraint is not linear" ) );

 /* We need to define the sense of the constraints as requested by CPLEX.
  *
  * In SMS++ FRowConstraints are defined as:
  * LHS <= ( some function from Variables to reals ) <= RHS
  *
  * CPLEX uses rhs and rngval arrays as documented in CPXcopylp reference. */

 auto con_lhs = con.get_lhs();
 auto con_rhs = con.get_rhs();

 if( con.is_relaxed() ) {
  // A relaxed constraint becomes: function >= -Inf
  sense[ row ] = 'G';
  rhs[ row ] = -Inf< double >();
  }
 else
  if( con_lhs == con_rhs ) {
   // LHS <= function <= RHS, with LHS = RHS
   // becomes: function = RHS
   sense[ row ] = 'E';
   rhs[ row ] = con_rhs;
   }
  else
   if( con_lhs == -Inf< double >() ) {
    // -inf <= function <= RHS
    // becomes: function <= RHS
    sense[ row ] = 'L';
    rhs[ row ] = con_rhs;
    }
   else if( con_rhs == Inf< double >() ) {
    // LHS <= function <= inf
    // becomes function >= LHS
    sense[ row ] = 'G';
    rhs[ row ] = con_lhs;
    }
   else {
    // LHS <= function <= RHS
    // becomes: LHS <= function <= LHS + (range),
    // with range = RHS - LHS
    sense[ row ] = 'R';
    rhs[ row ] = con_lhs;
    rngval[ row ] = con_rhs - con_lhs;
    }
 ++row;
 }

/*--------------------------------------------------------------------------*/

void MILPSolver::scan_objective( const FRealObjective * obj )
{
 // DEBUG_LOG( "MILPSolver::scan_objective() " << *obj );

 if( const auto * lf = dynamic_cast< const LinearFunction * >(
						 obj->get_function() ) ) {
  for( auto el : lf->get_v_var() )
   objective[ index_of_variable( el.first ) ] = el.second;

  return;
  }

 if( const auto * qf = dynamic_cast< const DQuadFunction * >(
						 obj->get_function() ) ) {
  for( auto el : qf->get_v_var() ) {
   // DQuadFunction::get_v_var() returns std::tuples of 3 elements
   auto k = index_of_variable( std::get< 0 >( el ) );
   objective[ k ] = std::get< 1 >( el );
   q_objective[ k ] = std::get< 2 >( el );
   }

  return;
  }

 throw( std::invalid_argument( "Unknown type of Objective Function" ) );
 }

/*--------------------------------------------------------------------------*/
/*----------------------------- MODIFICATIONS ------------------------------*/
/*--------------------------------------------------------------------------*/

int MILPSolver::compute( bool changedvars )
{
 MILPSolver::process_modifications();
 return( kOK );
 }

/*--------------------------------------------------------------------------*/

void MILPSolver::process_modifications( void )
{
 /* This function processes one modification after another, without any
  * attempt of optimization, moreover you have to be CAREFUL to write all
  * the cases in order from the most specialized to the more generic,
  * e.g., OneVarConstraintMod before RowConstraintMod before ConstraintMod,
  * otherwise the generic case will intercept the more specialized Mods. */

 // A function like this is needed to be called
 // recursively with GroupModifications
 std::function< void( sp_Mod ) > f;

 f = [ this, &f ]( const sp_Mod & mod ) {
  // DEBUG_LOG( *mod );

  if( const auto gm = std::dynamic_pointer_cast< GroupModification >( mod )
      ) {
   // DEBUG_LOG( "GroupModification containing:" << std::endl );
   for( const auto & submod : gm->sub_Modifications() )
    f( submod );

   return;
   }

  if( const auto vm = std::dynamic_pointer_cast< VariableMod >( mod ) ) {
   var_modification( vm.get() );
   return;
   }

  if( const auto om = std::dynamic_pointer_cast< ObjectiveMod >( mod ) ) {
   objective_modification( om.get() );
   return;
   }

  if( const auto bm = std::dynamic_pointer_cast< OneVarConstraintMod >(
								 mod ) ) {
   bound_modification( bm.get() );
   return;
   }

  if( const auto tmod = std::dynamic_pointer_cast< RowConstraintMod >( mod
								       ) ) {
   const_modification( tmod.get() );
   return;
   }

  if( const auto cm = std::dynamic_pointer_cast< ConstraintMod >( mod ) ) {
   const_modification( cm.get() );
   return;
   }

  if( const auto fm = std::dynamic_pointer_cast< FunctionMod >( mod ) ) {
   if( is_of( fm->function() ) )
    objective_function_modification( fm.get() );
   else
    constraint_function_modification( fm.get() );

   return;
   }

  if( const auto fvm = std::dynamic_pointer_cast< FunctionModVars >( mod )
      ) {
   if( is_of( fvm->function() ) )
    objective_fvars_modification( fvm.get() );
   else
    constraint_fvars_modification( fvm.get() );

   return;
   }

  if( const auto dm = std::dynamic_pointer_cast< BlockModAD >( mod ) ) {
   dynamic_modification( dm.get() );
   return;
   }

  if( std::dynamic_pointer_cast< NBModification >( mod ) )
   load_problem();
  };

 // Process the Modifications
 for( auto mod = front(); mod; mod = front() ) {
  f( mod );
  pop_front();

  if( std::dynamic_pointer_cast< NBModification >( mod ) ) {
   // An NBModification has just been handled.
   // All the remaining Modifications must be ignored.
   while( front() )
    pop_front();
   break;
   }
  }

 #ifdef MILPSOLVER_DEBUG
  check_status();
 #endif
 }

/*--------------------------------------------------------------------------*/

bool MILPSolver::is_of( Function * f )
{
 return( dynamic_cast< Objective * >( f->get_Observer() ) );
 }

/*--------------------------------------------------------------------------*/

void MILPSolver::var_modification( VariableMod * mod )
{
 auto * var = dynamic_cast< const ColVariable * >( mod->variable() );
 // TODO: Is dynamic_cast necessary?
 if( ! var ) // TODO: Throw exception?
  return;

 // Update the number of integer variables
 if( var->is_integer( mod->old_state() ) !=
     var->is_integer( mod->new_state() ) ) {
  if( var->is_integer( mod->new_state() ) )
   ++int_vars;
  else    
   --int_vars;
  }

 int idx = index_of_variable( var );

 // Update bounds
 if( ( ! lb.empty() ) && ( ! ub.empty() ) ) {
  lb[ idx ] = get_problem_lb( *var );
  ub[ idx ] = get_problem_ub( *var );
  }
 else
  throw( std::logic_error( "Bound representation is empty" ) );

 // Update variable type
 if( ! xctype.empty() ) {
  if( var->is_integer() && ( ! relax_int_vars ) ) {
   if( var->is_unitary() && var->is_positive() )
    xctype[ idx ] = 'B';
   else
    xctype[ idx ] = 'I';
   }
  else
   xctype[ idx ] = 'C';
  }
 else
  throw( std::logic_error( "Variable type representation is empty" ) );
 }

/*--------------------------------------------------------------------------*/

void MILPSolver::objective_modification( ObjectiveMod * mod )
{
 switch( mod->type() ) {
  case ObjectiveMod::eSetMin: objsense = 1; break;
  case ObjectiveMod::eSetMax: objsense = -1; break;
  default:
   throw( std::invalid_argument( "Invalid type of ObjectiveMod" ) );
  }
 }

/*--------------------------------------------------------------------------*/

void MILPSolver::const_modification( ConstraintMod * mod )
{
 // This modification has nothing to do if these are empty
 if( sense.empty() || rhs.empty() || rngval.empty() )
  throw( std::logic_error( "Constraint range representation is empty" ) );

 auto * con = dynamic_cast<FRowConstraint *>(mod->constraint());
 // TODO: Is dynamic_cast necessary?
 if( ! con )  // TODO: Throw exception?
  return;

 int idx = index_of_constraint( con );

 RowConstraint::RHSValue con_lhs = NAN;
 RowConstraint::RHSValue con_rhs = NAN;

 switch( mod->type() ) {
  case ConstraintMod::eRelaxConst:
   sense[ idx ] = 'G';
   rhs[ idx ] = -Inf< double >();
   rngval[ idx ] = 0;
   break;

  case ConstraintMod::eEnforceConst:
  case RowConstraintMod::eChgLHS:
  case RowConstraintMod::eChgRHS:
  case RowConstraintMod::eChgBTS:

   con_lhs = con->get_lhs();
   con_rhs = con->get_rhs();

   if( con_lhs == con_rhs ) {
    sense[ idx ] = 'E';
    rhs[ idx ] = con_rhs;
    rngval[ idx ] = 0;
    }
   else
    if( con_lhs == -Inf< double >() ) {
     sense[ idx ] = 'L';
     rhs[ idx ] = con_rhs;
     rngval[ idx ] = 0;
     }
    else
     if( con_rhs == Inf< double >() ) {
      sense[ idx ] = 'G';
      rhs[ idx ] = con_lhs;
      rngval[ idx ] = 0;
      }
     else {
      sense[ idx ] = 'R';
      rhs[ idx ] = con_lhs;
      rngval[ idx ] = con_rhs - con_lhs;
      }

   break;
  default:
   throw( std::invalid_argument( "Invalid type of ConstraintMod" ) );
  }
 }

/*--------------------------------------------------------------------------*/

void MILPSolver::bound_modification( OneVarConstraintMod * mod )
{
 // This modification has nothing to do if these are empty
 if( lb.empty() || ub.empty() )
  throw( std::logic_error( "Bound representation is empty" ) );

 auto con = static_cast< OneVarConstraint * >( mod->constraint() );
 auto var = static_cast< ColVariable * >( con->get_active_var( 0 ) );
 int idx = index_of_variable( var );

 switch( mod->type() ) {
  case RowConstraintMod::eChgLHS:
   lb[ idx ] = get_problem_lb( *var );
   break;
  case RowConstraintMod::eChgRHS:
   ub[ idx ] = get_problem_ub( *var );
   break;
  case RowConstraintMod::eChgBTS:
   lb[ idx ] = get_problem_lb( *var );
   ub[ idx ] = get_problem_ub( *var );
   break;
  default:
   throw( std::invalid_argument( "Invalid type of OneVarConstraintMod" ) );
  }
 }

/*--------------------------------------------------------------------------*/

void MILPSolver::objective_function_modification( FunctionMod * mod )
{
 // This modification has nothing to do if these are empty
 if( objective.empty() )
  throw( std::logic_error( "Objective representation is empty" ) );

 auto * f = mod->function();

 // C05FunctionModLin
 // --------------------------------------------------------------------------
 if( const auto * modl = dynamic_cast< C05FunctionModLin * >( mod ) ) {
  if( const auto * lf = dynamic_cast<const LinearFunction *>( f ) ) {
   for( Block::Index i = 0 ; i < modl->vars().size() ; ++i ) {
    auto var = static_cast<const ColVariable *>(modl->vars()[ i ]);
    objective[ index_of_variable( var ) ] += modl->delta()[ i ];
    }
   return;
   }

  // if( const auto * qf = dynamic_cast<const DQuadFunction *> (f) ) {
  //
  //  // This may happen if we change from LP to QP
  //  if( q_objective.empty() ) {
  //   q_objective.resize( numcols );
  //  }
  //
  //  for( auto i : sbst->subset() ) {
  //   auto var = static_cast<const ColVariable *>(qf->get_active_var( i ));
  //   auto idx = index_of_variable( var );
  //   objective[ idx ] = qf->get_linear_coefficient( i );
  //   q_objective[ idx ] = qf->get_quadratic_coefficient( i );
  //  }
  //
  //  return;
  // }

  // This should never happen
  throw( std::invalid_argument( "Unknown type of Objective Function" ) );
  }

 // Fallback method - Update all costs
 // --------------------------------------------------------------------------
 if( const auto * lf = dynamic_cast< const LinearFunction * >( f ) ) {
  // Linear objective function
  for( auto el : lf->get_v_var() )
   objective[ index_of_variable( el.first ) ] = el.second;

  return;
  }

 if( const auto * qf = dynamic_cast< const DQuadFunction * >( f ) ) {
  // Quadratic objective function

  // This may happen if we change from LP to QP
  if( q_objective.empty() )
   q_objective.resize( numcols );

  for( auto el : qf->get_v_var() ) {
   int idx = index_of_variable( std::get< 0 >( el ) );
   objective[ idx ] = std::get< 1 >( el );
   q_objective[ idx ] = std::get< 2 >( el );
   }

  return;
  }

 // This should never happen
 throw( std::invalid_argument( "Unsupported type of Objective function" ) );
 }

/*--------------------------------------------------------------------------*/

void MILPSolver::constraint_function_modification( FunctionMod * mod )
{
 if( matval.empty() )
  throw std::logic_error( "Constraint representation is empty" );

 auto * f = mod->function();

 const auto * lf = dynamic_cast< const LinearFunction * >( f );
 if( ! lf )
  return;

 auto * con = dynamic_cast< FRowConstraint * >( lf->get_Observer() );
 // TODO: Is dynamic_cast necessary?
 if( ! con )  // TODO: Throw exception?
  return;

 auto row = index_of_constraint( con );
 // TODO: update constraint matrix
 // C05FunctionModLin
 // --------------------------------------------------------------------------
 // if( const auto * modl = dynamic_cast<C05FunctionModLin *>(mod) ) {
 //
 //  for( int i = 0; i < modl->vars().size(); ++i ) {
 //   auto var = static_cast<const ColVariable *>(modl->vars()[ i ]);
 //   auto col = index_of_variable( var );
 //
 //   auto it = lower_bound( matind.begin() + matbeg[ col ],
 //                          matind.begin() + matbeg[ col + 1 ],
 //                          row );
 //   if( *it == row ) {
 //    matval[ std::distance( matind.begin(), it ) ] += modl->delta()[ i ];
 //   }
 //  }
 // }

 // Fallback method - Reload all coefficients
 // --------------------------------------------------------------------------
 // for( auto var : lf->get_v_var() ) {
 //  auto col = index_of_variable( var.first );
 //
 //  auto it = lower_bound( matind.begin() + matbeg[ col ],
 //                         matind.begin() + matbeg[ col + 1 ],
 //                         row );
 //  if( *it == row ) {
 //   matval[ std::distance( matind.begin(), it ) ] = var.second;
 //  }
 // }
 throw( std::logic_error(
		  "constraint_function_modification not implememted yet" ) );
 }

/*--------------------------------------------------------------------------*/

void MILPSolver::objective_fvars_modification( FunctionModVars * mod )
{
 // This modification has nothing to do if these are empty
 if( objective.empty() )
  throw( std::logic_error( "Objective representation is empty" ) );

 auto * f = mod->function();

 // Check the modification type
 if( dynamic_cast<C05FunctionModVarsAddd *>( mod ) == nullptr &&
     dynamic_cast<C05FunctionModVarsRngd *>( mod ) == nullptr &&
     dynamic_cast<C05FunctionModVarsSbst *>( mod ) == nullptr )
  throw( std::invalid_argument(
			 "This type of FunctionModVars is not handled" ) );

 if( const auto * lf = dynamic_cast< const LinearFunction * >( f ) ) {
  // Linear objective function

  for( auto * it1 : mod->vars() )
   for( auto it2: lf->get_v_var() )
    if( it1 == it2.first ) {
     if( mod->added() )
      objective[ index_of_variable( it2.first ) ] = it2.second;
     else
      objective[ index_of_variable( it2.first ) ] = 0;

     break;
     }

  return;
  }

 if( const auto * qf = dynamic_cast< const DQuadFunction * >( f ) ) {
  // Quadratic objective function

  for( auto * it1 : mod->vars() )
   for( auto it2: qf->get_v_var() )
    if( it1 == std::get< 0 >( it2 ) ) {
     int idx = index_of_variable( std::get< 0 >( it2 ) );
     if( mod->added() ) {
      objective[ idx ] = std::get< 1 >( it2 );
      q_objective[ idx ] = std::get< 2 >( it2 );
      }
     else {
      objective[ idx ] = 0;
      q_objective[ idx ] = 0;
      }
     break;
     }

  return;
  }

 // This should never happen
 throw( std::invalid_argument( "Unsupported type of Objective function" ) );
 }

/*--------------------------------------------------------------------------*/

void MILPSolver::constraint_fvars_modification( FunctionModVars * mod )
{
 if( matval.empty() )
  throw( std::logic_error( "Constraint representation is empty" ) );

 auto * f = mod->function();
 const auto * lf = dynamic_cast< const LinearFunction * >( f );
 if( ! lf )
  return;

 // Check the modification type
 if( dynamic_cast<C05FunctionModVarsAddd *>( mod ) == nullptr &&
     dynamic_cast<C05FunctionModVarsRngd *>( mod ) == nullptr &&
     dynamic_cast<C05FunctionModVarsSbst *>( mod ) == nullptr )
  throw( std::invalid_argument(
			 "This type of FunctionModVars is not handled" ) );

 // Get indices and coefficients
 auto * con = dynamic_cast< FRowConstraint * >( lf->get_Observer() );
 // TODO: Is dynamic_cast necessary?
 if( ! con )  // TODO: Throw exception?
  return;

 auto row = index_of_constraint( con );
 // TODO: Update coefficient matrix
 throw( std::logic_error(
		 "constraint_fvars_modification not implemented yet" ) );
 }

/*--------------------------------------------------------------------------*/

void MILPSolver::dynamic_modification( BlockModAD * mod )
{
 if( auto tmod = dynamic_cast< BlockModAdd< FRowConstraint > * >( mod ) ) {
  for( auto * i : tmod->added() )
   add_dynamic_constraint( i );

  return;
  }

 if( auto tmod = dynamic_cast< BlockModRmv< FRowConstraint > * >( mod ) ) {
  for( const auto & i : tmod->removed() )
   remove_dynamic_constraint( & i );

  return;
  }

 if( auto tmod = dynamic_cast< BlockModAdd< ColVariable > * >( mod ) ) {
  for( auto * i : tmod->added() )
   add_dynamic_variable( i );

  return;
  }

 if( auto tmod = dynamic_cast< BlockModRmv< ColVariable > * >( mod ) ) {
  for( const auto & i : tmod->removed() )
   remove_dynamic_variable( &i );

  return;
  }

 // FIXME: It should be a generic OneVarConstraint!
 if( auto tmod = dynamic_cast< BlockModAdd< LB0Constraint > * >( mod ) ) {
  for( auto * i : tmod->added() )
   add_dynamic_bound( i );

  return;
  }

 // FIXME: It should be a generic OneVarConstraint!
 if( auto tmod = dynamic_cast< BlockModRmv< LB0Constraint > * >( mod ) ) {
  for( const auto & i : tmod->removed() )
   remove_dynamic_bound( &i );

  return;
  }

 throw( std::invalid_argument( "Unknown type of BlockModAD" ) );
 }

/*--------------------------------------------------------------------------*/

void MILPSolver::add_dynamic_constraint( FRowConstraint * con )
{
 auto f = dynamic_cast< const LinearFunction * >( con->get_function() );
 if( ! f )
  throw( std::invalid_argument( "The Constraint is not linear" ) );

 // Update the dictionaries
 auto it = lower_bound( dcon_to_idx.begin() , dcon_to_idx.end() , con ,
                        []( con_int pair , FRowConstraint * c ) {
                         return pair.first < c;
                         } );
 dcon_to_idx.insert( it , { con , numrows } );
 idx_to_dcon.emplace_back( con );

 // Update the counter
 ++numrows;

 // Update the vectors
 if( sense.empty() || rhs.empty() || rngval.empty() )
  throw( std::logic_error( "Constraint range representation is empty" ) );

 auto con_lhs = con->get_lhs();
 auto con_rhs = con->get_rhs();

 if( con_lhs == con_rhs ) {
  sense.emplace_back( 'E' );
  rhs.emplace_back( con_rhs );
  rngval.emplace_back( 0 );
  }
 else
  if( con_lhs == -Inf< double >() ) {
   sense.emplace_back( 'L' );
   rhs.emplace_back( con_rhs );
   rngval.emplace_back( 0 );
   }
  else
   if( con_rhs == Inf< double >() ) {
    sense.emplace_back( 'G' );
    rhs.emplace_back( con_lhs );
    rngval.emplace_back( 0 );
    }
   else {
    sense.emplace_back( 'R' );
    rhs.emplace_back( con_lhs );
    rngval.emplace_back( con_rhs - con_lhs );
    }


 if( matval.empty() )
  throw( std::logic_error( "Constraint representation is empty" ) );

 throw( std::logic_error(
		"MILPSolver::add_dynamic_constraint not implememted yet" ) );
 }

/*--------------------------------------------------------------------------*/

void MILPSolver::add_dynamic_variable( ColVariable * var )
{
 // Update the dictionaries
 auto it = lower_bound( dvar_to_idx.begin() , dvar_to_idx.end() , var ,
                        []( var_int pair , ColVariable * v ) {
                         return pair.first < v;
                         } );

 dvar_to_idx.insert( it, { var, numcols } );
 idx_to_dvar.emplace_back( var );

 // Update the counters
 ++numcols;
 if( var->is_integer() )
  ++int_vars;

 // Update the bound vectors
 if( lb.empty() || ub.empty() )
  throw( std::logic_error( "Bound representation is empty" ) );

 lb.emplace_back( get_problem_lb( *var ) );
 ub.emplace_back( get_problem_ub( *var ) );

 // Update the variable type vector
 if( xctype.empty() )
  throw( std::logic_error( "Variable type representation is empty" ) );

 if( var->is_integer() && ( ! relax_int_vars ) ) {
  if( var->is_unitary() && var->is_positive() )
   xctype.emplace_back( 'B' );
  else
   xctype.emplace_back( 'I' );
  }
 else
  xctype.emplace_back( 'C' );

 // Update the objective vectors
 if( objective.empty() )
  throw( std::logic_error( "Objective representation is empty" ) );

 objective.push_back( 0 );
 if( ! q_objective.empty() )
  q_objective.push_back( 0 );

 if( matval.empty() )
  throw( std::logic_error( "Constraint representation is empty" ) );

 // TODO: update constraint matrix
 // matbeg.push_back( matbeg.back() );
 // matcnt.push_back( 0 );
 throw( std::logic_error(
		 "MILPSolver::add_dynamic_variable not implememted yet" ) );
 }

/*--------------------------------------------------------------------------*/

void MILPSolver::add_dynamic_bound( OneVarConstraint * con )
{
 // This modification has nothing to do if these are empty
 if( lb.empty() || ub.empty() )
  throw( std::logic_error( "Bound representation is empty" ) );

 auto var = static_cast< ColVariable * >( con->get_active_var( 0 ) );
 int idx = index_of_variable( var );
 lb[ idx ] = get_problem_lb( *var );
 ub[ idx ] = get_problem_ub( *var );
 }

/*--------------------------------------------------------------------------*/

void MILPSolver::remove_dynamic_constraint( const FRowConstraint * con )
{
 // TODO: Implement remove_dynamic_with_index(i)?

 // Update the dictionaries
 int index = 0;
 auto it1 = lower_bound( dcon_to_idx.begin() , dcon_to_idx.end() ,
                         std::make_pair( con , 0 ) ,
                         []( auto & p1 , auto & p2 ) {
                          return p1.first < p2.first;
                          } );

 if( ( it1 != dcon_to_idx.end() ) && ( it1->first == con ) ) {
  index = it1->second;
  dcon_to_idx.erase( it1 );
  idx_to_dcon.erase( idx_to_dcon.begin() + index - static_cons );
  }
 else
  throw( std::runtime_error( "Dynamic constraint not found" ) );

 for( auto & it: dcon_to_idx )
  if( it.second > index )
   it.second--;

 // Update the counter
 --numrows;

 // Update the vectors
 if( sense.empty() || rhs.empty() || rngval.empty() )
  throw( std::logic_error( "Constraint range representation is empty" ) );

 sense.erase( sense.begin() + index );
 rhs.erase( rhs.begin() + index );
 rngval.erase( rngval.begin() + index );


 if( matval.empty() )
  throw( std::logic_error( "Constraint representation is empty" ) );

 // TODO: update constraint matrix
 throw( std::logic_error(
	    "MILPSolver::remove_dynamic_constraint not implememted yet" ) );
 }

/*--------------------------------------------------------------------------*/

void MILPSolver::remove_dynamic_variable( const ColVariable * var )
{
 // TODO: Implement remove_dynamic_with_index(i)?

 // Update the dictionaries
 int index = 0;
 auto it1 = lower_bound( dvar_to_idx.begin() , dvar_to_idx.end() ,
                         std::make_pair( var , 0 ),
                         []( auto & p1 , auto & p2 ) {
                          return p1.first < p2.first;
                          } );

 if( ( it1 != dvar_to_idx.end() ) && ( it1->first == var ) ) {
  index = it1->second;
  dvar_to_idx.erase( it1 );
  idx_to_dvar.erase( idx_to_dvar.begin() + index - static_vars );
  }
 else
  throw( std::runtime_error( "Dynamic variable not found" ) );

 for( auto & it : dvar_to_idx )
  if( it.second > index )
   it.second--;

 // Update the counters
 --numcols;
 if( var->is_integer() )
  --int_vars;

 // Update the bound vectors
 if( lb.empty() || ub.empty() )
  throw( std::logic_error( "Bound representation is empty" ) );

 lb.erase( lb.begin() + index );
 ub.erase( ub.begin() + index );

 // Update the variable type vector
 if( xctype.empty() )
  throw( std::logic_error( "Variable type representation is empty" ) );

 xctype.erase( xctype.begin() + index );

 // Update the objective vectors
 if( objective.empty() )
  throw( std::logic_error( "Objective representation is empty" ) );

 objective.erase( objective.begin() + index );
 if( ! q_objective.empty() )
  q_objective.erase( q_objective.begin() + index );

 if( matval.empty() )
  throw( std::logic_error( "Constraint representation is empty" ) );

 // TODO: update constraint matrix
 throw( std::logic_error(
	       "MILPSolver::remove_dynamic_variable not implememted yet" ) );
 }

/*--------------------------------------------------------------------------*/

void MILPSolver::remove_dynamic_bound( const OneVarConstraint * con )
{
 // This modification has nothing to do if these are empty
 if( lb.empty() || ub.empty() )
  throw( std::logic_error( "Bound representation is empty" ) );

 auto var = static_cast< ColVariable * >( con->get_active_var( 0 ) );
 int idx = index_of_variable( var );
 lb[ idx ] = get_problem_lb( * var );
 ub[ idx ] = get_problem_ub( * var );
 }

/*--------------------------------------------------------------------------*/
/*------------------------------- PARAMETERS -------------------------------*/
/*--------------------------------------------------------------------------*/

void MILPSolver::set_par( idx_type par , int value )
{
 if( par == intUseCustomNames ) {
  use_custom_names = bool( value );
  return;
  }
 if( par == intRelaxIntVars ) {
  relax_int_vars = bool( value );
  return;
  }

 CDASolver::set_par( par, value );
 }

/*----------------------------------------------------------------------------

void MILPSolver::set_par( idx_type par , double value )
{
 CDASolver::set_par( par, value );
 }

----------------------------------------------------------------------------*/

void MILPSolver::set_par( idx_type par , std::string && value )
{
 if( par == strProblemName ) {
  prob_name = std::move( value );
  return;
  }
 if( par == strOutputFile ) {
  output_file = std::move( value );
  return;
  }

 CDASolver::set_par( par, std::move( value ) );
 }

/*--------------------------------------------------------------------------*/

ThinComputeInterface::idx_type MILPSolver::get_num_int_par( void ) const {
 return( CDASolver::get_num_int_par() + intLastAlgParMILP - intLastParCDAS );
 }

ThinComputeInterface::idx_type MILPSolver::get_num_dbl_par( void ) const {
 return( CDASolver::get_num_dbl_par() + dblLastAlgParMILP - dblLastParCDAS );
 }

ThinComputeInterface::idx_type MILPSolver::get_num_str_par( void ) const {
 return( CDASolver::get_num_str_par() + strLastAlgParMILP - strLastParCDAS );
 }

/*--------------------------------------------------------------------------*/

int MILPSolver::get_dflt_int_par( idx_type par ) const
{
 if( par == intUseCustomNames )
  return( 1 );

 if( par == intRelaxIntVars )
  return( 0 );

 return( CDASolver::get_dflt_int_par( par ) );
 }

/*----------------------------------------------------------------------------

double MILPSolver::get_dflt_dbl_par( idx_type par ) const
{
 return( CDASolver::get_dflt_dbl_par( par ) );
 }

----------------------------------------------------------------------------*/

const std::string & MILPSolver::get_dflt_str_par( idx_type par ) const
{
 static const std::vector< std::string > vals = { "MILPSolver_prob", "" };
 if( par == strProblemName )
  return( vals[ 0 ] );

 if( par == strOutputFile )
  return( vals[ 1 ] );

 return( CDASolver::get_dflt_str_par( par ) );
 }

/*--------------------------------------------------------------------------*/

int MILPSolver::get_int_par( idx_type par ) const
{
 if( par == intUseCustomNames )
  return( use_custom_names );

 if( par == intRelaxIntVars )
  return( relax_int_vars );

 return( CDASolver::get_int_par( par ) );
 }

/*----------------------------------------------------------------------------

double MILPSolver::get_dbl_par( idx_type par ) const
{
 return( CDASolver::get_dbl_par( par ) );
 }

----------------------------------------------------------------------------*/

const std::string & MILPSolver::get_str_par( idx_type par ) const
{
 if( par == strProblemName )
  return( prob_name );

 if( par == strOutputFile )
  return( output_file );

 return( CDASolver::get_str_par( par ) );
 }

/*--------------------------------------------------------------------------*/

ThinComputeInterface::idx_type
MILPSolver::int_par_str2idx( const std::string & name ) const
{
 if( name == "intUseCustomNames" )
  return( intUseCustomNames );

 if( name == "intRelaxIntVars" )
  return( intRelaxIntVars );

 return( CDASolver::int_par_str2idx( name ) );
 }

/*--------------------------------------------------------------------------*/

const std::string & MILPSolver::int_par_idx2str( idx_type idx ) const
{
 static const std::vector< std::string > pars = { "intUseCustomNames",
                                                  "intRelaxIntVars" };
 if( idx == intUseCustomNames )
  return( pars[ 0 ] );

 if( idx == intRelaxIntVars )
  return( pars[ 1 ] );

 return( CDASolver::int_par_idx2str( idx ) );
 }

/*----------------------------------------------------------------------------

ThinComputeInterface::idx_type
MILPSolver::dbl_par_str2idx( const std::string & name ) const {
 return( CDASolver::dbl_par_str2idx( name ) );
 }

------------------------------------------------------------------------------

const std::string &
MILPSolver::dbl_par_idx2str( idx_type idx ) const {
 return( CDASolver::dbl_par_idx2str( idx ) );
 }

----------------------------------------------------------------------------*/

ThinComputeInterface::idx_type
MILPSolver::str_par_str2idx( const std::string & name ) const
{
 if( name == "strProblemName" )
  return( strProblemName );

 if( name == "strOutputFile" )
  return( strOutputFile );

 return( CDASolver::str_par_str2idx( name ) );
 }

/*--------------------------------------------------------------------------*/

const std::string & MILPSolver::str_par_idx2str( idx_type idx ) const
{
 static const std::vector< std::string > pars = { "strProblemName" ,
                                                  "strOutputFile" };
 if( idx == strProblemName )
  return( pars[ 0 ] );

 if( idx == strOutputFile )
  return( pars[ 1 ] );

 return( CDASolver::str_par_idx2str( idx ) );
 }

/*--------------------------------------------------------------------------*/
/*---------------------------------- DEBUG ---------------------------------*/
/*--------------------------------------------------------------------------*/

template< typename T >
std::string log_vector( const std::vector< T > & v, int limit ) {
 std::string temp_log = "[";
 for( auto i : v ) {
  temp_log += " " + std::to_string( i );
  if( --limit == 0 ) {
   temp_log += " ... ";
   break;
  }
 }
 temp_log += "]";
 return temp_log;
}

template<>
std::string log_vector( const std::vector< char > & v, int limit ) {
 std::string temp_log = "[";
 for( auto i : v ) {
  temp_log += " ";
  temp_log += i;
  if( --limit == 0 ) {
   temp_log += " ... ";
   break;
  }
 }
 temp_log += "]";
 return temp_log;
}

/*--------------------------------------------------------------------------*/

#ifdef MILPSOLVER_DEBUG

void MILPSolver::check_status( void )
{
 int v = 0;
 int c = 0;
 int sv = 0;
 int sc = 0;
 int svg = 0;
 int scg = 0;
 int dv = 0;
 int dc = 0;

 // ------------------ Count everything -------------------
 std::queue< Block * > Q;

 // Locking the Block
 bool owned = f_Block->is_owned_by( f_id );
 if( ( ! owned ) && ( ! f_Block->read_lock() ) )
  throw( std::runtime_error( "Unable to lock the Block" ) );

 Q.push( f_Block );
 while( ! Q.empty() ) {
  Block * q_Block = Q.front();
  Q.pop();

  for( auto * i : q_Block->get_nested_Blocks() ) {
   Q.push( i );
  }

  for( const auto & i : q_Block->get_static_constraints() ) {
   // Singles
   if( un_any_thing_0( FRowConstraint , i ,
                       {
                        ++scg;
                        ++c;
                        ++sc;
                        }
		       ) )
    continue;

   // Vectors
   if( un_any_thing_1( FRowConstraint , i ,
                       {
                        ++scg;
                        c += var.size();
                        sc += var.size();
                        }
		       ) )
    continue;

   // Multiarrays
   if( un_any_thing_K( FRowConstraint, i,
                       {
                        ++scg;
                        c += var.num_elements();
                        sc += var.num_elements();
                        }
		       ) )
    continue;
   }

  for( const auto & i : q_Block->get_dynamic_constraints() ) {
   // Single lists
   if( un_any_thing_0( std::list< FRowConstraint > , i ,
                       { c += var.size();  } ) )
    continue;

   // Vectors of lists
   if( un_any_thing_1( std::list< FRowConstraint >, i,
                       {
                        for( auto & el: var )
                         c += el.size();
                        }
		       ) )
    continue;

   // Multiarrays of lists
   if( un_any_thing_K( std::list< FRowConstraint > , i ,
                       {
                        auto it = var.data();
                        for( auto i = var.num_elements() ; i-- ; ++it )
                         c += it->size();
                        }
		       ) )
    continue;
   }

  dc = c - sc;

  for( const auto & i : q_Block->get_static_variables() ) {
   // Singles
   if( un_any_thing_0( ColVariable , i ,
                       {
                        ++svg;
                        ++v;
                        ++sv;
                        }
		       ) )
    continue;

   // Vectors
   if( un_any_thing_1( ColVariable , i ,
                       {
                        ++svg;
                        v += var.size();
                        sv += var.size();
                        }
		       ) )
    continue;

   // Multiarrays
   if( un_any_thing_K( ColVariable , i ,
                       {
                        ++svg;
                        v += var.num_elements();
                        sv += var.num_elements();
                        }
		       ) )
    continue;
   }

  for( const auto & i : q_Block->get_dynamic_variables() ) {
   // Single lists
   if( un_any_thing_0( std::list< ColVariable > , i ,
                       { v += var.size(); } ) )
    continue;

   // Vectors of lists
   if( un_any_thing_1( std::list< ColVariable > , i ,
                       {
                        for( auto & el: var )
                         v += el.size();
                        }
		       ) )
    continue;

   // Multiarrays of lists
   if( un_any_thing_K( std::list< ColVariable > , i ,
                       {
                        auto it = var.data();
                        for( auto i = var.num_elements() ; i-- ; ++it )
                         v += it->size();
                        }
		       ) )
    continue;
   }

  dv = v - sv;
  }

 // Unlock the Block
 if( ! owned )
  f_Block->read_unlock();

 if( numcols != v )
  DEBUG_LOG( "numcols is " << numcols << ", it should be " << v
	     << std::endl );

 if( numrows != c )
  DEBUG_LOG( "numrows is " << numrows << ", it should be " << c
	     << std::endl );

 if( static_vars != sv )
  DEBUG_LOG( "static_vars is " << static_vars
                               << ", it should be " << sv << std::endl );

 if( static_cons != sc )
  DEBUG_LOG( "static_cons is " << static_cons
                               << ", it should be " << sc << std::endl );

 // ------------------ Svar dictionaries ------------------

 for( auto & i: idx_to_svar ) {
  auto j = std::find_if( svar_to_idx.begin(), svar_to_idx.end(),
                         [ & ]( auto & pair ) {
                          return std::get< 1 >( pair ) == i.first &&
                                 std::get< 0 >( pair ) == i.second;
                         } );
  if( j == svar_to_idx.end() ) {
   DEBUG_LOG( "Element [" << i.first << ", " << i.second
                          << "] of idx_to_svar was not found in svar_to_idx"
                          << std::endl );
  }
 }

 for( auto & i: svar_to_idx ) {
  auto j = std::find_if( idx_to_svar.begin(), idx_to_svar.end(),
                         [ & ]( auto & pair ) {
                          return std::get< 0 >( i ) == pair.second &&
                                 std::get< 1 >( i ) == pair.first;
                         } );
  if( j == idx_to_svar.end() ) {
   DEBUG_LOG( ", " << std::get< 1 >( i ) <<
                   "] of svar_to_idx was not found in idx_to_svar"
                   << std::endl );
  }
 }

 if( idx_to_svar.size() != svg ) {
  DEBUG_LOG( "Size of idx_to_svar is " << idx_to_svar.size()
                                       << ", it should be " << sv
                                       << std::endl );
 }

 if( svar_to_idx.size() != svg ) {
  DEBUG_LOG( "Size of svar_to_idx is " << svar_to_idx.size()
                                       << ", it should be " << svg
                                       << std::endl );
 }

 // ------------------ Dvar dictionaries ------------------

 for( auto & i: idx_to_dvar ) {
  auto j = std::find_if( dvar_to_idx.begin(), dvar_to_idx.end(),
                         [ & ]( auto & pair ) { return pair.first == i; } );
  if( j == dvar_to_idx.end() ) {
   DEBUG_LOG( "Element [" << i
                          << "] of idx_to_dvar was not found in dvar_to_idx"
                          << std::endl );
  }
 }

 for( auto & i: dvar_to_idx ) {
  auto j = std::find_if( idx_to_dvar.begin(), idx_to_dvar.end(),
                         [ & ]( auto & var ) { return i.first == var; } );
  if( j == idx_to_dvar.end() ) {
   DEBUG_LOG( "Element [" << i.first << ", " << i.second
                          << "] of dvar_to_idx was not found in idx_to_dvar"
                          << std::endl );
  }
 }

 if( idx_to_dvar.size() != dv ) {
  DEBUG_LOG( "Size of idx_to_dvar is " << idx_to_dvar.size()
                                       << ", it should be " << dv
                                       << std::endl );
 }

 if( dvar_to_idx.size() != dv ) {
  DEBUG_LOG( "Size of dvar_to_idx is " << dvar_to_idx.size()
                                       << ", it should be " << dv
                                       << std::endl );
 }

 // ------------------ Scon dictionaries ------------------

 for( auto & i: idx_to_scon ) {
  auto j = std::find_if( scon_to_idx.begin(), scon_to_idx.end(),
                         [ & ]( auto & pair ) {
                          return std::get< 1 >( pair ) == i.first &&
                                 std::get< 0 >( pair ) == i.second;
                         } );
  if( j == scon_to_idx.end() ) {
   DEBUG_LOG( "Element [" << i.first << ", " << i.second
                          << "] of idx_to_scon was not found in scon_to_idx"
                          << std::endl );
  }
 }

 for( auto & i: scon_to_idx ) {
  auto j = std::find_if( idx_to_scon.begin(), idx_to_scon.end(),
                         [ & ]( auto & pair ) {
                          return std::get< 0 >( i ) == pair.second &&
                                 std::get< 1 >( i ) == pair.first;
                         } );
  if( j == idx_to_scon.end() ) {
   DEBUG_LOG( ", " << std::get< 1 >( i )
                   << "] of scon_to_idx was not found in idx_to_scon"
                   << std::endl );
  }
 }

 if( idx_to_scon.size() != scg ) {
  DEBUG_LOG( "Size of idx_to_scon is " << idx_to_scon.size()
                                       << ", it should be " << scg
                                       << std::endl );
 }

 if( scon_to_idx.size() != scg ) {
  DEBUG_LOG( "Size of scon_to_idx is " << scon_to_idx.size()
                                       << ", it should be " << scg
                                       << std::endl );
 }

 // ------------------ Dcon dictionaries ------------------

 for( auto & i: idx_to_dcon ) {
  auto j = std::find_if( dcon_to_idx.begin(), dcon_to_idx.end(),
                         [ & ]( auto & pair ) { return pair.first == i; } );
  if( j == dcon_to_idx.end() ) {
   DEBUG_LOG( "Element [" << i
                          << "] of idx_to_dcon was not found in dcon_to_idx"
                          << std::endl );
  }
 }
 for( auto & i: dcon_to_idx ) {
  auto j = std::find_if( idx_to_dcon.begin(), idx_to_dcon.end(),
                         [ & ]( auto & con ) { return i.first == con; } );
  if( j == idx_to_dcon.end() ) {
   DEBUG_LOG( "Element [" << i.first << ", " << i.second
                          << "] of dcon_to_idx was not found in idx_to_dcon"
                          << std::endl );
  }
 }

 if( idx_to_dcon.size() != dc ) {
  DEBUG_LOG( "Size of idx_to_dcon is " << idx_to_dcon.size()
                                       << ", it should be " << dc
                                       << std::endl );
 }

 if( dcon_to_idx.size() != dc ) {
  DEBUG_LOG( "Size of dcon_to_idx is " << dcon_to_idx.size()
                                       << ", it should be " << dc
                                       << std::endl );
 }


 // ------------------ Check int_vars ------------------

 int actual_int_vars = 0;

}

#endif // MILPSOLVER_DEBUG

/*--------------------------------------------------------------------------*/
/*----------------------- End File MILPSolver.cpp --------------------------*/
/*--------------------------------------------------------------------------*/ 
