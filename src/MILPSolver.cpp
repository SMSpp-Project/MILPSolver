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
 * \copyright &copy; by Antonio Frangioni, Niccolo' Iardella
 */
/*--------------------------------------------------------------------------*/
/*---------------------------- IMPLEMENTATION ------------------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*------------------------------- MACROS -----------------------------------*/
/*--------------------------------------------------------------------------*/
/* If the macro MILPSolver_DEBUG is externally defined, then some costly
 * checks on the data structures of MILPSolver are performed and debug
 * information printed. Also, the method check_status() is defined and
 * used to check the whole set of data structures. */

#ifdef MILPSolver_DEBUG
 #include <queue>
 #define DEBUG_LOG( stuff ) std::cout << "[MILPSolver DEBUG] " << stuff
#else
 #define DEBUG_LOG( stuff )
#endif

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include <map>

#include <set>

#include <LinearFunction.h>
#include <PolyhedralFunction.h>
#include <C05Function.h>

#include <iostream>

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
   throw( std::runtime_error(
              "MILPSolver::set_Block: unable to lock the Block" ) );

  // generate abstract representation
  block->generate_abstract_variables();
  block->generate_abstract_constraints();
  block->generate_objective();

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
  q_part.clear();

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
  ndq_objective.clear();
  ndq_rowind.clear();
  ndq_colind.clear();
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
 // Clean any left-up structures
 clear_problem( 15 );
 
 numrows = 0;
 numcols = 0;
 static_vars = 0;
 static_cons = 0;
 static_quadcons = 0;
 constant_value = 0;
 Index nzelements = 0;
 Index nst_linrow = 0;
 Index nst_quadrow = 0;
 Index ndy_linrow = 0;
 Index ndy_quadrow = 0;

 int_vars = 0;
 numquadrows = 0;
 numnnzq = 0;
 objsense = 0;

 // locking the Block
 bool owned = f_Block->is_owned_by( f_id );
 if( ( ! owned ) && ( ! f_Block->read_lock() ) )
  throw( std::runtime_error(
             "MILPSolver::load_problem: unable to lock the Block" ) );

 // construct the set of involved Block in BFS order, skipping any
 // sub-Block marked as ignored. The exclusion list is owned by the
 // Solver base (Solver::f_excluded) and stored as the *minimal*
 // user-supplied set (no eager descendant enumeration): the literal
 // .count() test at every immediate child is therefore sufficient,
 // because once a child is excluded its subtree is never enqueued and
 // the BFS naturally prunes (no need to walk the father chain here)
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 const auto & ignored = get_excluded_blocks();

 v_BFS.clear();
 if( ignored.count( f_Block ) == 0 )
  v_BFS.push_back( f_Block );
 for( Index i = 0 ; i < v_BFS.size() ; ++i )
  for( auto el : v_BFS[ i ]->get_nested_Blocks() )
   if( ignored.count( el ) == 0 )
    v_BFS.push_back( el );

 v_BFS.shrink_to_fit();

 // count variables and constraints - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // note that we only count FRowConstraint, but we also allow
 // OneVarConstraint which are used to set bounds but do not produce any row
 // in the coefficient matrix. for this reason we don't throw exception if
 // the constraint is not a FRowConstraint (we should check if it is a
 // OneVarConstraint and throw exception if not, but this is pesky because
 // we'd need to check all derived classes to OneVarConstraint so we avoid
 // it). conversely we only allow ColVariable, so we immediately throw
 // exception if any Variable, be it static or dynamic, is not a ColVariable

 Index num_block = 0;        // counter for the blocks
 Index row = 0;              // counter for the rows
 Index col = 0;              // counter for the columns
 Index static_con_grps = 0;  // counter for static constraint groups
 Index static_var_grps = 0;  // counter for static variable groups

 for( auto qb : v_BFS ) {
  DEBUG_LOG( "Processing Block " << num_block << " [" << qb << "]"
	     << std::endl );
  DEBUG_LOG( *qb << std::endl );

  // Static constraints
  for( const auto & i : qb->get_static_constraints() ) {
   auto count = un_any_thing_count_static( FRowConstraint , i );
   if( count != Inf< std::size_t >() ) {
    numrows += count;
    static_cons += count;
    ++static_con_grps;
    continue;
    }

   // if it's not FRowConstraint, accept any known OneVarConstraint silently
   // note that the second argument of the macro is empty
   if( un_any_thing_OneVarConstraint_static( i , ) )
    continue;

   throw( std::invalid_argument(
              "MILPSolver::load_problem: static constraint is neither "
              "FRowConstraint nor OneVarConstraint" ) );
   }

  // Dynamic constraints
  for( const auto & i : qb->get_dynamic_constraints() ) {
   auto count = un_any_thing_count_dynamic( FRowConstraint , i );
   if( count != Inf< std::size_t >() ) {
    numrows += count;
    continue;
    }

   // if it's not FRowConstraint, accept any known OneVarConstraint silently
   // note that the second argument of the macro is empty
   if( un_any_thing_OneVarConstraint_dynamic( i , ) )
    continue;

   throw( std::invalid_argument(
              "MILPSolver::load_problem: dynamic constraint is neither "
              "FRowConstraint nor OneVarConstraint" ) );
   }

  auto counter_static_lin_quad_row = [ this , & nst_linrow, & nst_quadrow ]
                              ( FRowConstraint & cons ) {
   if( dynamic_cast< LinearFunction * >( cons.get_function() ) )
    ++nst_linrow;
   else
    if( dynamic_cast< DQuadFunction * >( cons.get_function() ) )
     ++nst_quadrow;
    else
     if( dynamic_cast< QuadFunction * >( cons.get_function() ) )
      ++nst_quadrow;
   };

  auto counter_dynamic_lin_quad_row = [ this , & ndy_linrow, & ndy_quadrow ]
                              ( FRowConstraint & cons ) {
   if( dynamic_cast< LinearFunction * >( cons.get_function() ) )
    ++ndy_linrow;
   else
    if( dynamic_cast< DQuadFunction * >( cons.get_function() ) )
     ++ndy_quadrow;
    else
     if( dynamic_cast< QuadFunction * >( cons.get_function() ) )
      ++ndy_quadrow;
   };

  for( const auto & i : qb->get_static_constraints() )
   un_any_const_static( i , counter_static_lin_quad_row ,
			un_any_type< FRowConstraint >() );

  for( const auto & i : qb->get_dynamic_constraints() )
   un_any_const_dynamic( i , counter_dynamic_lin_quad_row , 
			 un_any_type< FRowConstraint >() );

  // Fill number of quadratic rows
  numquadrows = nst_quadrow + ndy_quadrow;
  static_quadcons = nst_quadrow;

  for( const auto & i : qb->get_static_variables() ) {
   auto count = un_any_thing_count_static( ColVariable , i );
   if( count == Inf< std::size_t >() )
    throw( std::invalid_argument(
               "MILPSolver::load_problem: not a ColVariable" ) );
   numcols += count;
   static_vars += count;
   ++static_var_grps;
   }

  for( const auto & i : qb->get_dynamic_variables() ) {
   auto count = un_any_thing_count_dynamic( ColVariable , i );
   if( count == Inf< std::size_t >() )
    throw( std::invalid_argument(
               "MILPSolver::load_problem: not a ColVariable" ) );
   numcols += count;
   }

  auto counter = [ this , & nzelements ]( ColVariable & var ) {
   // count each constraint the variable is active in ONCE: a variable may
   // appear more than once in the same constraint (so active_stuff() lists
   // it more than once), but it contributes a single matrix entry (with the
   // coefficients summed); see the matching coalescing in scan_variable()
   std::set< const FRowConstraint * > seen;
   for( auto * i : var.active_stuff() )
    if( auto * row = dynamic_cast< FRowConstraint * >( i ) )
      if( is_mine( row->get_Block() ) )
       if( seen.insert( row ).second )
        ++nzelements;
   };

  for( const auto & i : qb->get_static_variables() )
   un_any_const_static( i , counter , un_any_type< ColVariable >() );

  for( const auto & i : qb->get_dynamic_variables() )
   un_any_const_dynamic( i , counter , un_any_type< ColVariable >() );

  ++num_block;
  }

 DEBUG_LOG( "Number of blocks      = " << num_block << std::endl );
 DEBUG_LOG( "numrows (constraints) = " << numrows << " (S:" << static_cons
	    << "/D:" << numrows - static_cons << ")" << std::endl );
 DEBUG_LOG( "numcols (variables)   = " << numcols << " (S:" << static_vars
	    << "/D:" << numcols - static_vars << ")" << std::endl );
 DEBUG_LOG( "nzelements            = " << nzelements << std::endl );

 // MILP vectors allocation - - - - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 // The +1 is needed by generic interface

 // If no quadratic constraints are in the model, then the coefficient matrix
 // is grouped by columns. Otherwise, it is grouped by rows.
 if( numquadrows == 0 ) {
  matbeg.resize( numcols + 1, 0 );
  matbeg[ numcols ] = nzelements;
  matcnt.resize( numcols, 0 );
  }
 else {
  matbeg.resize( numrows + 1, 0 );
  matbeg[ numrows ] = nzelements;
  matcnt.resize( numrows , 0 );
  q_part.resize( numrows );
  }

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

 // Acccount also for link between variables and bound
 if( single_bound ) {
  svar_to_bound.clear();
  dvar_to_bound.clear();
  svar_to_bound.reserve( static_vars );
  dvar_to_bound.reserve( numcols - static_vars );
  }

 /* Now we have to check wheter there are any quadratic constraints. 
  * If this is the case, then matbeg, matcnt, ... will represent the 
  * coefficient matrix by rows. Thus, we will firstly scan the variables
  * and after the constraints. Otherwise, if no quadratic costraints 
  * are in the model, the coefficient matrix will be grouped by column,
  * and so we will firstly scan the constraint and after the variables. */

 if( numquadrows == 0 ) {
 // scan the static constraints - - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 num_block = 0;
 // quad_row = nst_linrow + ndy_linrow;
 // quadratic rows starts after linear ones
 for( auto qb : v_BFS ) {
  Index set = 0;  // counter for the constraint groups

  for( const auto & i : qb->get_static_constraints() ) {
    // Call specific function to scan the new group of Constraints
    scan_group( i , qb , num_block, set, row, un_any_type< FRowConstraint >() );

    set++;
   }
  num_block++;
  }
 std::sort( scon_to_idx.begin() , scon_to_idx.end() );

 // scan the dynamic constraints- - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 num_block = 0;
 for( auto qb : v_BFS ) {
  int set = 0; // Counter for the constraint groups

  for( const auto & i : qb->get_dynamic_constraints() ) {
   Index start = row;

   auto scan = [ this , & row ]
    ( const FRowConstraint & c ) { scan_dynamic_constraint( c , row ); };
   un_any_const_dynamic( i , scan , un_any_type< FRowConstraint >() );

   //  write names
   auto base = qb->get_d_const_name()[ set ];
   Index end = row - start;
   for( Index n = 0 ; n < end ; ++n ) {
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
   }
  num_block++;
  }

 std::sort( dcon_to_idx.begin() , dcon_to_idx.end() );

 } // end if( numquadrows == 0 )

 // scan the static variables - - - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 num_block = 0;
 for( auto qb : v_BFS ) {
  Index set = 0;   // counter for the variable groups

  for( const auto & i : qb->get_static_variables() ) {
   // Call specific function to scan the new group of Constraints
   scan_group( i , qb , num_block, set, col , un_any_type< ColVariable >() );

   set++;
   }
  num_block++;
  }

 std::sort( svar_to_idx.begin() , svar_to_idx.end() );

 // scan the dynamic variables- - - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 num_block = 0;
 for( auto qb : v_BFS ) {
  Index set = 0;   // Counter for the variable groups

  for( const auto & i : qb->get_dynamic_variables() ) {
   Index start = col;
   auto scan = [ this , & col ]( const ColVariable & v ) {
    scan_dynamic_variable( v ,  col );
    };
   un_any_const_dynamic( i , scan , un_any_type< ColVariable >() );

   // write names
   auto base = qb->get_d_var_name()[ set ];
   Index end = col - start;
   for( Index n = 0 ; n < end ; ++n ) {
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

   // If the option single_bound is true, we have to check that maximum one
   // OneVarConstraint is associated with a single variable.
   // Moreover, the vector linking the variable with the associated bound,
   // needs to be filled.
   if( single_bound ) {
    auto scan_bound = [ this ]( const ColVariable & v ) {
      scan_dynamic_variable_bound( v );
    };
    
    un_any_const_dynamic( i , scan_bound , un_any_type< ColVariable >() );
    }
  }
  num_block++;
  }

 std::sort( dvar_to_idx.begin() , dvar_to_idx.end() );

 if( numquadrows != 0 ) {
  // scan the static constraints - - - - - - - - - - - - - - - - - - - - - - -
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  num_block = 0;
  //quad_row = nst_linrow + ndy_linrow; // quadratic rows starts after linear ones
  for( auto qb : v_BFS ) {
   Index set = 0;  // counter for the constraint groups

   for( const auto & i : qb->get_static_constraints() ) {
    // Call specific function to scan the new group of Constraints
    scan_group( i , qb , num_block, set, row, un_any_type< FRowConstraint >() );

    set++;
   }
  num_block++;
  }
 std::sort( scon_to_idx.begin() , scon_to_idx.end() );

 // scan the dynamic constraints- - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 num_block = 0;
 for( auto qb : v_BFS ) {
  int set = 0; // Counter for the constraint groups

  for( const auto & i : qb->get_dynamic_constraints() ) {
   Index start = row;

   auto scan = [ this , & row ]
      ( const FRowConstraint & c ) {
        scan_dynamic_constraint( c , row );
    };
   un_any_const_dynamic( i , scan , un_any_type< FRowConstraint >() );

   //  write names
   auto base = qb->get_d_const_name()[ set ];
   Index end = row - start;
   for( Index n = 0 ; n < end ; ++n ) {
    std::string name;
    if( base.empty() )
     name = "cs_" + std::to_string( num_block )
          + "_" + std::to_string( set ) + "_" + std::to_string( n );
    else
     name = base + "_" + std::to_string( num_block )
          + "_" + std::to_string( n );

    rowname[ start + n ] = strcpy( new char[ name.length() + 1 ] , name.c_str() );
    }
   set++;
   }
  num_block++;
  }

 std::sort( dcon_to_idx.begin() , dcon_to_idx.end() );
 } // end if( numquadrows != 0)

 // scan the objective- - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 objsense = 0;
 for( auto qb : v_BFS ) {
  switch( qb->get_objective_sense() ) {
   case( Objective::eMax ):
    if( objsense == 1 )
     throw( std::invalid_argument(
                "MILPSolver::load_problem: mixed max/min "
                "Objective not supported" ) );
    objsense = -1; break;
   case( Objective::eMin ):
    if( objsense == -1 )
     throw( std::invalid_argument(
                "MILPSolver::load_problem: mixed max/min "
                "Objective not supported" ) );
     objsense = 1;
   }

  if( auto * obj = dynamic_cast< FRealObjective * >( qb->get_objective() ) )
   scan_objective( obj );
  }

 // unlock the Block- - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 if( ! owned )
  f_Block->read_unlock();

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

 if( ! objsense )  // not defined anywhere in the Block
  objsense = 1;    // take one pick (minimization)

 }  // end( MILPSolver::load_problem )

/*--------------------------------------------------------------------------*/

template< typename T >
 void MILPSolver::scan_group( const boost::any & gr , Block * qb ,
                              Index num_block , Index set ,
                              Index & counter , un_any_type< T > )
{
 // Search for the group type
 if( gr.type() == typeid( T * ) ||
      gr.type() == typeid( std::vector< T > * ) ||
      gr.type() == typeid( std::vector< std::vector< T > > * ) ) {
  // "Simple" group
  scan_st_group( gr , qb , num_block , set , counter , un_any_type< T >() );
  }
 else {
  // "Complex" group
  scan_multiarray_st_group( gr , qb , num_block , set , counter , 
    un_any_type< T >() );
  }
 } // end( MILPSolver::scan_group )

/*--------------------------------------------------------------------------*/

template< typename T >
 void MILPSolver::scan_st_group( const boost::any & gr , Block * qb , 
                                Index num_block , Index set , Index & counter , 
                                un_any_type< T > )
{
 Index elements = 0;  // counter for group elements
 Index start = counter;

 if( typeid( T * ) == typeid( FRowConstraint * ) ) {
  // Scanning a group of Constraints
  auto scan = [ this , & elements , & counter ]
   ( const FRowConstraint & c ) {
    scan_static_constraint( c , elements , counter );
  };
  un_any_const_static( gr , scan , un_any_type< FRowConstraint >() );

  //  write names
  auto base = qb->get_s_const_name()[ set ];
  Index end = counter - start;
  for( Index n = 0 ; n < end ; ++n ) {
   std::string name;
   if( base.empty() )
    name = "cs_" + std::to_string( num_block )
      + "_" + std::to_string( set ) + "_" + std::to_string( n );
   else
    name = base + "_" + std::to_string( num_block )
      + "_" + std::to_string( n );

    rowname[ start + n ] = strcpy( new char[ name.length() + 1 ] , name.c_str() );
  }

  // Add number of elements of the group
  if( elements )
    std::get< 2 >( scon_to_idx.back() ) = elements;
 }
 else if( typeid( T * ) == typeid( ColVariable * ) ) {
  // Scanning a group of Variables
  auto scan = [ this , & elements , & counter ]
   ( const ColVariable & c ) {
    scan_static_variable( c , elements , counter );
  };
  un_any_const_static( gr , scan , un_any_type< ColVariable >() );

  //  write names
  auto base = qb->get_s_var_name()[ set ];
  Index end = counter - start;
  for( Index n = 0 ; n < end ; ++n ) {
   std::string name;
   if( base.empty() )
    name = "xs_" + std::to_string( num_block )
     + "_" + std::to_string( set ) + "_" + std::to_string( n );
   else
    name = base + "_" + std::to_string( num_block )
     + "_" + std::to_string( n );

   colname[ start + n ] = strcpy( new char[ name.length() + 1 ] , name.c_str() );
  }

  // Add number of elements of the group
  if( elements )
   std::get< 2 >( svar_to_idx.back() ) = elements;

  // If the option single_bound is true, we have to check that maximum one
  // OneVarConstraint is associated with a single variable.
  // Moreover, the vector linking the variable with the associated bound,
  // needs to be filled.
  if( single_bound == true ) {
   auto scan_bound = [ this ]( const ColVariable & v ) {
    scan_static_variable_bound( v );
   };
    
   un_any_const_static( gr , scan_bound , un_any_type< ColVariable >() );
  }
 }
 else
  throw( std::runtime_error(
             "MILPSolver::scan_st_group: unsupported group type" ) );
} // end( MILPSolver::scan_st_group )

/*--------------------------------------------------------------------------*/

template< typename T >
 void MILPSolver::scan_multiarray_st_group( const boost::any & gr ,
            Block * qb , Index num_block , Index set , 
            Index & counter , un_any_type< T > )
{
 // Get multi array number of dimension
 int ma_dim = get_multi_array_dim( gr , un_any_type< T >() , 
                                    un_any_int< 2 >() );

 // Name of the group
 auto base = (typeid(T*) == typeid(FRowConstraint*))
                ? qb->get_s_const_name()[set]
                : qb->get_s_var_name()[set];

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
     Index elements = 0;  // counter for group elements
     Index start = counter;

     auto scan = [ this , & elements , & counter ]
      ( const FRowConstraint & c ) {
        scan_static_constraint( c , elements , counter );
     };
     // Scan a single group
     un_any_const_static( v , scan , un_any_type< FRowConstraint >() );

     //  write names
     Index end = counter - start;
     for( Index n = 0 ; n < end ; ++n ) {
      std::string name;
      if( base.empty() )
       name = "cs_" + std::to_string( num_block )
        + "_" + std::to_string( set ) + "_" 
        + std::to_string( idx_0 ) + "_" 
        + std::to_string( idx_1 ) + "_" + std::to_string( n );
      else
       name = base + "_" + std::to_string( num_block ) + "_"
        + std::to_string( idx_0 ) + "_" 
        + std::to_string( idx_1 ) + "_" + std::to_string( n );

     rowname[ start + n ] = strcpy( new char[ name.length() + 1 ] , name.c_str() );
     }

     // Add number of elements of the group
     if( elements )
      std::get< 2 >( scon_to_idx.back() ) = elements;

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
     Index start = counter;
     Index elements = 0;  // counter for group elements

     auto scan = [ this , & elements , & counter ]
      ( const ColVariable & c ) {
       scan_static_variable( c , elements , counter );
     };
     un_any_const_static( v , scan , un_any_type< ColVariable >() );

     //  write names
     Index end = counter - start;
     for( Index n = 0 ; n < end ; ++n ) {
      std::string name;
      if( base.empty() )
       name = "xs_" + std::to_string( num_block )
        + "_" + std::to_string( set ) + "_" 
        + std::to_string( idx_0 ) + "_" 
        + std::to_string( idx_1 ) + "_" + std::to_string( n );
      else
       name = base + "_" + std::to_string( num_block ) + "_"
        + std::to_string( idx_0 ) + "_" 
        + std::to_string( idx_1 ) + "_" + std::to_string( n );

      colname[ start + n ] = strcpy( new char[ name.length() + 1 ] , name.c_str() );
     }
      
     // Add number of elements of the group
     if( elements )
      std::get< 2 >( svar_to_idx.back() ) = elements;

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

     // If the option single_bound is true, we have to check that maximum one
     // OneVarConstraint is associated with a single variable.
     // Moreover, the vector linking the variable with the associated bound,
     // needs to be filled.
     if( single_bound == true ) {
      auto scan_bound = [ this ]( const ColVariable & c ) {
       scan_static_variable_bound( c );
      };
    
     un_any_const_static( v , scan_bound , un_any_type< ColVariable >() );
     }
    }
   }
   else
    throw( std::runtime_error(
               "MILPSolver::scan_multiarray_st_group: "
               "unsupported group type" ) );
  }
  else if( type == 0 ) {
   // Multi arrays of type 0 (i.e. multi_array< T >) store
   // elements in sequential cells. Thus, we can store them
   // as usually done for std::vector< T > by only keeping track
   // of the first element and storing the number of non empty
   // cells in the structure.
   auto ma = get_multi_array0( gr , un_any_type< T >() ,
                              un_any_int< 2 >() );
   Index elements = 0;  // counter for group elements

   if( typeid( T * ) == typeid( FRowConstraint * ) ) {
    // Constraint group

    // Scan the linearization of the array
    for( auto v = ma->data() ; idx_0 < ma->shape()[ 0 ] ; ++v ) {
     auto scan = [ this , & elements , & counter ]
      ( const FRowConstraint & c ) {
       scan_static_constraint( c , elements , counter );
     };
     // Scan a single group
     un_any_const_static( v , scan , un_any_type< FRowConstraint >() );

     //  write names
     std::string name;
     if( base.empty() )
      name = "cs_" + std::to_string( num_block )
       + "_" + std::to_string( set ) + "_" 
       + std::to_string( idx_0 ) + "_" 
       + std::to_string( idx_1 );
     else
      name = base + "_" + std::to_string( num_block ) + "_"
       + std::to_string( idx_0 ) + "_" 
       + std::to_string( idx_1 );

     rowname[ counter - 1 ] = strcpy( new char[ name.length() + 1 ] , name.c_str() );

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
    
    // Add number of elements of the group
    if( elements )
     std::get< 2 >( scon_to_idx.back() ) = elements;
   }
   else if( typeid( T * ) == typeid( ColVariable * ) ) {
    // Variable group

    // Scan the linearization of the array
    for( auto v = ma->data() ; idx_0 < ma->shape()[ 0 ] ; ++v ) {
     auto scan = [ this , & elements , & counter ]
      ( const ColVariable & c ) {
       scan_static_variable( c , elements , counter );
     };
     un_any_const_static( v , scan , un_any_type< ColVariable >() );

     //  write names
     std::string name;
     if( base.empty() )
      name = "xs_" + std::to_string( num_block )
       + "_" + std::to_string( set ) + "_" 
       + std::to_string( idx_0 ) + "_" 
       + std::to_string( idx_1 );
     else
      name = base + "_" + std::to_string( num_block ) + "_"
       + std::to_string( idx_0 ) + "_" 
       + std::to_string( idx_1 );

     colname[ counter - 1 ] = strcpy( new char[ name.length() + 1 ] , name.c_str() );

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

     // If the option single_bound is true, we have to check that maximum one
     // OneVarConstraint is associated with a single variable.
     // Moreover, the vector linking the variable with the associated bound,
     // needs to be filled.
     if( single_bound == true ) {
      auto scan_bound = [ this ]( const ColVariable & c ) {
       scan_static_variable_bound( c );
      };
      un_any_const_static( v , scan_bound , un_any_type< ColVariable >() );
     }
    }

    // Add number of elements of the group
    if( elements )
     std::get< 2 >( svar_to_idx.back() ) = elements;
   }
   else
    throw( std::runtime_error(
               "MILPSolver::scan_multiarray_st_group: "
               "unsupported group type" ) );
  }
  else
   throw( std::runtime_error(
              "MILPSolver::scan_multiarray_st_group: "
              "unsupported multi-array type" ) );
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
     Index start = counter;
     Index elements = 0;  // counter for group elements

     auto scan = [ this , & elements , & counter ]
      ( const FRowConstraint & c ) {
       scan_static_constraint( c , elements , counter );
     };
     // Scan a single group
     un_any_const_static( v , scan , un_any_type< FRowConstraint >() );

     //  write names
     Index end = counter - start;
     for( Index n = 0 ; n < end ; ++n ) {
      std::string name;
      if( base.empty() )
       name = "cs_" + std::to_string( num_block )
        + "_" + std::to_string( set ) + "_" 
        + std::to_string( idx_0 ) + "_" 
        + std::to_string( idx_1 ) + "_" 
        + std::to_string( idx_2 ) + "_" + std::to_string( n );
      else
       name = base + "_" + std::to_string( num_block ) + "_"
        + std::to_string( idx_0 ) + "_" 
        + std::to_string( idx_1 ) + "_" 
        + std::to_string( idx_2 ) + "_" + std::to_string( n );

      rowname[ start + n ] = strcpy( new char[ name.length() + 1 ] , name.c_str() );
     }

     // Add number of elements of the group
     if( elements )
      std::get< 2 >( scon_to_idx.back() ) = elements;

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
     Index start = counter;
     Index elements = 0;  // counter for group elements

     auto scan = [ this , & elements , & counter ]
      ( const ColVariable & c ) {
       scan_static_variable( c , elements , counter );
     };
     un_any_const_static( v , scan , un_any_type< ColVariable >() );

     //  write names
     Index end = counter - start;
     for( Index n = 0 ; n < end ; ++n ) {
      std::string name;
      if( base.empty() )
       name = "xs_" + std::to_string( num_block )
        + "_" + std::to_string( set ) + "_" 
        + std::to_string( idx_0 ) + "_" 
        + std::to_string( idx_1 ) + "_" 
        + std::to_string( idx_2 ) + "_" + std::to_string( n );
      else
       name = base + "_" + std::to_string( num_block ) + "_"
        + std::to_string( idx_0 ) + "_" 
        + std::to_string( idx_1 ) + "_" 
        + std::to_string( idx_2 ) + "_" + std::to_string( n );

      colname[ start + n ] = strcpy( new char[ name.length() + 1 ] , name.c_str() );
     }

     // Add number of elements of the group
     if( elements )
      std::get< 2 >( svar_to_idx.back() ) = elements;

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

     // If the option single_bound is true, we have to check that maximum one
     // OneVarConstraint is associated with a single variable.
     // Moreover, the vector linking the variable with the associated bound,
     // needs to be filled.
     if( single_bound == true ) {
      auto scan_bound = [ this ]( const ColVariable & c ) {
       scan_static_variable_bound( c );
      };
        
      un_any_const_static( v , scan_bound , un_any_type< ColVariable >() );
     }
    }
   }
   else
    throw( std::runtime_error(
               "MILPSolver::scan_multiarray_st_group: "
               "unsupported group type" ) );
  }
  else if( type == 0 ) {
   // Multi arrays of type 0 (i.e., multi_array< T >) store
   // elements in sequential cells. Thus, we can store them
   // as usually done for std::vector< T > by only keeping track
   // of the first element and storing the number of non-empty
   // cells in the structure.
   auto ma = get_multi_array0( gr , un_any_type< T >() ,
                              un_any_int< 3 >() );

   Index elements = 0;  // counter for group elements

   if( typeid( T * ) == typeid( FRowConstraint * ) ) {
    // Constraint group

    // Scan the linearization of the array
    for( auto v = ma->data() ; idx_0 < ma->shape()[ 0 ] ; ++v ) {
     auto scan = [ this , & elements , & counter ]
      ( const FRowConstraint & c ) {
       scan_static_constraint( c , elements , counter );
      };
     // Scan a single group
     un_any_const_static( v , scan , un_any_type< FRowConstraint >() );

     //  write names
     std::string name;
     if( base.empty() )
      name = "cs_" + std::to_string( num_block )
       + "_" + std::to_string( set ) + "_" 
       + std::to_string( idx_0 ) + "_" 
       + std::to_string( idx_1 ) + "_" 
       + std::to_string( idx_2 );
     else
      name = base + "_" + std::to_string( num_block ) + "_"
       + std::to_string( idx_0 ) + "_" 
       + std::to_string( idx_1 ) + "_" 
       + std::to_string( idx_2 );

     rowname[ counter - 1 ] = strcpy( new char[ name.length() + 1 ] , name.c_str() );

     // The linearization produced by ma->data() for the 3D multi_array
     // stores elements in row-major order.
     if( idx_2 < ma->shape()[ 2 ] - 1 )
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
    
    // Add number of elements of the group
    if( elements )
     std::get< 2 >( scon_to_idx.back() ) = elements;
   }
   else if( typeid( T * ) == typeid( ColVariable * ) ) {
    // Variable group

    // Scan the linearization of the array
    for( auto v = ma->data() ; idx_0 < ma->shape()[ 0 ] ; ++v ) {
     auto scan = [ this , & elements , & counter ]
      ( const ColVariable & c ) {
       scan_static_variable( c , elements , counter );
     };
     un_any_const_static( v , scan , un_any_type< ColVariable >() );

     //  write names
     std::string name;
     if( base.empty() )
      name = "xs_" + std::to_string( num_block )
       + "_" + std::to_string( set ) + "_" 
       + std::to_string( idx_0 ) + "_" 
       + std::to_string( idx_1 ) + "_" 
       + std::to_string( idx_2 );
     else
      name = base + "_" + std::to_string( num_block ) + "_"
       + std::to_string( idx_0 ) + "_" 
       + std::to_string( idx_1 ) + "_" 
       + std::to_string( idx_2 );

     colname[ counter - 1 ] = strcpy( new char[ name.length() + 1 ] , name.c_str() );

     // The linearization produced by ma->data() for the 3D multi_array
     // stores elements in row-major order.
     if( idx_2 < ma->shape()[ 2 ] - 1 )
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

     // If the option single_bound is true, we have to check that maximum one
     // OneVarConstraint is associated with a single variable.
     // Moreover, the vector linking the variable with the associated bound,
     // needs to be filled.
     if( single_bound == true ) {
      auto scan_bound = [ this ]( const ColVariable & c ) {
       scan_static_variable_bound( c );
      };
    
      un_any_const_static( v , scan_bound , un_any_type< ColVariable >() );
     }
    }
    
    // Add number of elements of the group
    if( elements )
     std::get< 2 >( svar_to_idx.back() ) = elements;
   }
   else
    throw( std::runtime_error(
               "MILPSolver::scan_multiarray_st_group: "
               "unsupported group type" ) );
  }
  else
   throw( std::runtime_error(
              "MILPSolver::scan_multiarray_st_group: "
              "unsupported multi-array type" ) );
 }
 else
    // Handle invalid or unsupported ma_dim 
    return; 

} // end( MILPSolver::scan_multiarray_st_group )

/*--------------------------------------------------------------------------*/

double MILPSolver::get_problem_lb( const ColVariable & var ) const
{
 double b = var.get_lb();

 for( auto * i : var.active_stuff() )
  if( auto box = dynamic_cast< OneVarConstraint * >( i ) )
   b = std::max( b , box->get_lhs() );

 return( b );
 }

/*--------------------------------------------------------------------------*/

double MILPSolver::get_problem_ub( const ColVariable & var ) const
{
 double b = var.get_ub();

 for( auto * i : var.active_stuff() )
  if( auto box = dynamic_cast< OneVarConstraint * >( i ) )
   b = std::min( b , box->get_rhs() );

 return( b );
 }

/*--------------------------------------------------------------------------*/

std::array< double , 2 > MILPSolver::get_problem_bounds(
					      const ColVariable & var ) const
{
 std::array< double , 2 > ret;
 ret[ 0 ] = var.get_lb();
 ret[ 1 ] = var.get_ub();

 for( auto * i : var.active_stuff() )
  if( auto box = dynamic_cast< OneVarConstraint * >( i ) ) {
   ret[ 0 ] = std::max( ret[ 0 ] , box->get_lhs() );
   ret[ 1 ] = std::min( ret[ 1 ] , box->get_rhs() );
   }

 return( ret );
 }

/*--------------------------------------------------------------------------*/

std::vector< FRowConstraint * > MILPSolver::get_active_constraints(
					      const ColVariable & var ) const
{
 std::vector< FRowConstraint * > active_constraints;
 for( auto * i : var.active_stuff() )
  if( auto * row = dynamic_cast< FRowConstraint * >( i ) )
    if( is_mine( row->get_Block() ) )
     active_constraints.push_back( row );
 
 return( active_constraints );
 }

/*--------------------------------------------------------------------------*/

std::vector< OneVarConstraint * > MILPSolver::get_active_bounds(
					      const ColVariable & var , bool first_scan ) const
{
 std::vector< OneVarConstraint * > active_bounds;
 
 /* If the option intSingleBound is activated, then each constraint can have 
 *  just a single OneVarConstraint associated, which is stored in the vectors
 *  svar_to_bound and dvar_to_bound. 
 *  We also check if we are calling this function from the load_problem: in
 *  this case (i.e. first_scan = true) we still have to fill the dictionaries.
 */
 if( single_bound && ( ! first_scan ) ) {
  int idx = index_of_variable( &var ); // get variable index

  if( idx < static_vars ) { // the variable is static
   if( svar_to_bound[ idx ] != nullptr )
    active_bounds.push_back( 
      const_cast< OneVarConstraint * >( svar_to_bound[ idx ] ) );
   }
  else { // the variable is dynamic
   if( dvar_to_bound[ idx - static_vars ] != nullptr )
    active_bounds.push_back( 
     const_cast< OneVarConstraint * >( dvar_to_bound[ idx - static_vars ] ) );
   }
  }
 else { // The option is not activated, scan all active stuff
  for( auto * i : var.active_stuff() )
   if( auto row = dynamic_cast< OneVarConstraint * >( i ) )
    active_bounds.push_back( row );
  }

 return( active_bounds );
 }

/*--------------------------------------------------------------------------*/
/*-------------------- METHODS FOR PROBLEM DESCRIPTION ---------------------*/
/*--------------------------------------------------------------------------*/

int MILPSolver::index_of_variable( const ColVariable * var ) const
{
 auto i = index_of_static_variable( var );
 return( i < Inf< int >() ? i : index_of_dynamic_variable( var ) );
 }

/*--------------------------------------------------------------------------*/

int MILPSolver::index_of_static_variable( const ColVariable * var ) const
{
 if( svar_to_idx.empty() )
  return( Inf< int >() );

 #ifdef MILPSolver_DEBUG
  assert( std::is_sorted( svar_to_idx.begin() , svar_to_idx.end() ) );
 #endif
 auto it = upper_bound( svar_to_idx.begin() , svar_to_idx.end(),
                        std::make_tuple( var , 0 , 0 ) ,
                        [ & ]( auto & p1 , auto & p2 ) {
                         return( std::get< 0 >( p1 ) < std::get< 0 >( p2 ) );
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

int MILPSolver::index_of_dynamic_variable( const ColVariable * var ) const
{
 #ifdef MILPSolver_DEBUG
  assert( std::is_sorted( dvar_to_idx.begin() , dvar_to_idx.end() ) );
 #endif
 auto it = lower_bound( dvar_to_idx.begin() , dvar_to_idx.end(),
                        std::make_pair( var , 0 ) ,
                        [ & ]( auto & p1 , auto & p2 ) {
                         return( p1.first < p2.first );
                        } );

 if( ( it != dvar_to_idx.end() ) && ( it->first == var ) )
  return( it->second );

 return( Inf< int >() );
 }

/*--------------------------------------------------------------------------*/

int MILPSolver::index_of_constraint( const FRowConstraint * con ) const
{
 auto i = index_of_static_constraint( con );
 return( i < Inf< int >() ? i : index_of_dynamic_constraint( con ) );
 }

/*--------------------------------------------------------------------------*/

int MILPSolver::index_of_static_constraint( const FRowConstraint * con ) const
{
 if( scon_to_idx.empty() )
  return( Inf< int >() );

 #ifdef MILPSolver_DEBUG
  assert( std::is_sorted( scon_to_idx.begin(), scon_to_idx.end() ) );
 #endif
 auto it = upper_bound( scon_to_idx.begin() , scon_to_idx.end() ,
                        std::make_tuple( con , 0 , 0 , false ) ,
                        [ & ]( auto & p1 , auto & p2 ) {
                         return( std::get< 0 >( p1 ) < std::get< 0 >( p2 ) );
                        } );

 // Now it refers to the first (group of) element(s) greater than var
 if( it == scon_to_idx.begin() )
  return( Inf< int >() );

 --it;

 // First element of the constraint group
 auto first = const_cast< const FRowConstraint * >( std::get< 0 >( *it ) );
 int distance = std::distance( first , con );

 if( ( distance >= 0 ) && ( distance < std::get< 2 >( *it ) ) )
  // The element belongs to this group
  return( std::get< 1 >( *it ) + distance );

 // The element doesn't exist
 return( Inf< int >() );
 }

/*--------------------------------------------------------------------------*/

int MILPSolver::index_of_dynamic_constraint( const FRowConstraint * con )
 const {
 #ifdef MILPSolver_DEBUG
  assert( std::is_sorted( dcon_to_idx.begin() , dcon_to_idx.end() ) );
 #endif
 auto it = lower_bound( dcon_to_idx.begin() , dcon_to_idx.end() ,
                        std::make_pair( con , 0 ),
                        [ & ]( auto & p1 , auto & p2 ) {
                         return( p1.first < p2.first );
                        } );

 if( ( it != dcon_to_idx.end() ) && ( it->first == con ) )
  return( it->second );

 return( Inf< int >() );
 }

/*--------------------------------------------------------------------------*/

const ColVariable * MILPSolver::variable_with_index( int i ) const
{
 if( i < static_vars )
  return( static_variable_with_index( i ) );

 return( dynamic_variable_with_index( i ) );
 }

/*--------------------------------------------------------------------------*/

const ColVariable * MILPSolver::static_variable_with_index( int i ) const
{
 if( idx_to_svar.empty() )
  return( nullptr );

 #ifdef MILPSolver_DEBUG
  assert( std::is_sorted( idx_to_svar.begin(), idx_to_svar.end() ) );
 #endif
 auto it = upper_bound( idx_to_svar.begin(), idx_to_svar.end(),
                        std::make_pair( i, nullptr ),
                        [ & ]( auto & p1, auto & p2 ) {
                         return( p1.first < p2.first );
                        } );

 // Now it refers to the first (group of) element(s) greater than i
 if( it == idx_to_svar.begin() )
  return( nullptr );

 --it;

 int distance = i - it->first;
 return( it->second + distance );
 }

/*--------------------------------------------------------------------------*/

const ColVariable * MILPSolver::dynamic_variable_with_index( int i ) const
{
 if( ( static_vars < i ) || ( i < numcols ) )
  return( idx_to_dvar[ i - static_vars ] );

 return( nullptr );
 }

/*--------------------------------------------------------------------------*/

const FRowConstraint * MILPSolver::constraint_with_index( int i ) const
{
 if( i < static_cons )
  return( static_constraint_with_index( i ) );

 return( dynamic_constraint_with_index( i ) );
 }

/*--------------------------------------------------------------------------*/

const FRowConstraint * MILPSolver::static_constraint_with_index( int i )
 const {
 if( idx_to_scon.empty() )
  return( nullptr );

 #ifdef MILPSolver_DEBUG
  assert( std::is_sorted( idx_to_scon.begin() , idx_to_scon.end() ) );
 #endif
 auto it = upper_bound( idx_to_scon.begin() , idx_to_scon.end() ,
                        std::make_pair( i , nullptr ) ,
                        [ & ]( auto & p1 , auto & p2 ) {
                         return( p1.first < p2.first );
                        } );

 // Now it refers to the first (group of) element(s) greater than i
 if( it == idx_to_scon.begin() )
  return( nullptr );

 --it;

 int distance = i - it->first;
 return( it->second + distance );
 }

/*--------------------------------------------------------------------------*/

const FRowConstraint * MILPSolver::dynamic_constraint_with_index( int i )
 const {
 if( ( static_cons < i ) || ( i < numrows ) )
  return( idx_to_dcon[ i - static_cons ] );

 return( nullptr );
 }

/*--------------------------------------------------------------------------*/
/*------------- AUXILIARY METHODS FOR POPULATING THE PROBLEM  --------------*/
/*--------------------------------------------------------------------------*/

void MILPSolver::scan_static_variable( const ColVariable & var , Index & n ,
				       Index & col )
{
 if( ! n++ ) {  // the tuple's third field will be filled later
  svar_to_idx.emplace_back( & var , col , 0 );
  idx_to_svar.emplace_back( col , & var );
  }

 scan_variable( var , col );
 }

/*--------------------------------------------------------------------------*/

void MILPSolver::scan_dynamic_variable( const ColVariable & var ,
					Index & col )
{
 dvar_to_idx.emplace_back( & var , col );
 idx_to_dvar.emplace_back( & var );
 scan_variable( var , col );
 }

/*--------------------------------------------------------------------------*/

void MILPSolver::scan_variable( const ColVariable & var , Index & col )
{
 // Check if the Variable is not empty
 if( var.get_Block() == nullptr )
  throw( std::invalid_argument(
             "MILPSolver::scan_variable: the provided variable is empty" ) );

 auto bd = MILPSolver::get_problem_bounds( var );
 if( var.is_fixed() ) {
  lb[ col ] = std::max( bd[ 0 ] , var.get_value() );
  ub[ col ] = std::min( bd[ 1 ] , var.get_value() );
  }
 else {
  lb[ col ] = bd[ 0 ];
  ub[ col ] = bd[ 1 ];
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

 // Get linear active constraints
 auto active_constraints = get_active_constraints( var );

 /* We have to check wheter we have quadratic constraints in the model or not.
  * If the model is LP, then matbeg, matcnt, ... store the matrix coefficients
  * grouped by column. */
 if( numquadrows == 0 ) {
  // coalesce the column: one matrix entry per DISTINCT constraint the
  // variable is active in, with coefficient = the SUM of all of the
  // variable's coefficients in that constraint. A variable may legitimately
  // appear more than once in a single LinearFunction (e.g. in some AC
  // network formulations): emitting one (row, col) entry per occurrence is
  // both wrong (the coefficients must add up, not repeat) and rejected by
  // some backends (HiGHS >= 1.14 errors on duplicate row indices). The
  // std::map keeps the entries sorted by row index, as the solvers expect;
  // it stays consistent with the distinct count of the nzelements counter.
  std::map< int , double > col_entries;
  for( auto * con : active_constraints ) {
   const int row = index_of_constraint( con );
   if( col_entries.count( row ) )
    continue;                  // this constraint has already been summed
   auto * f = static_cast< const LinearFunction * >( con->get_function() );
   double coeff = 0;
   bool found = false;
   for( const auto & pr : f->get_v_var() )
    if( pr.first == & var ) {
     coeff += pr.second;
     found = true;
     }
   if( ! found )
    // this should never happen since we loop on the active constraints
    throw( std::invalid_argument(
               "MILPSolver::scan_variable: "
               "this ColVariable is not active in the examined "
               "FRowConstraint" ) );
   col_entries[ row ] = coeff;
   }

  matcnt[ col ] = col_entries.size();
  matbeg[ col ] = ( col == 0 ) ? 0 : matbeg[ col - 1 ] + matcnt[ col - 1 ];

  Index j = 0;
  for( const auto & [ row , coeff ] : col_entries ) {
   matind[ matbeg[ col ] + j ] = row;
   matval[ matbeg[ col ] + j ] = coeff;
   ++j;
   }
  }
 ++col;
 }

/*--------------------------------------------------------------------------*/

void MILPSolver::scan_static_constraint( const FRowConstraint & con ,
					 Index & n , Index & row )
{
 if( ! n++ ) {  // the tuple's third field will be filled later
  scon_to_idx.emplace_back( & con , row , 0 );
  idx_to_scon.emplace_back( row , & con );
  }

 scan_constraint( con , row );
 }

/*--------------------------------------------------------------------------*/

void MILPSolver::scan_dynamic_constraint( const FRowConstraint & con ,
					  Index & row )
{
 dcon_to_idx.emplace_back( & con , row );
 idx_to_dcon.emplace_back( & con );

 scan_constraint( con , row );
 }

/*--------------------------------------------------------------------------*/

void MILPSolver::scan_constraint( const FRowConstraint & con , Index & row  )
{
 // Check if the Constraint is not empty
 if( con.get_Block() == nullptr )
  throw( std::invalid_argument(
             "MILPSolver::scan_constraint: "
             "the provided constraint is empty" ) );

 /* We have to check wheter we have quadratic constraints in the model or not.
  * If the model is QP, then matbeg, matcnt, ... store the matrix coefficients
  * grouped by rows. */
 if( numquadrows != 0 ) {
  if( auto f = con.get_function() ) {

   /* If the DEBUG is activated, we can check wether a variable appears 
   * multiple times in a single constraint, which will clearly produce
   * an error in later stage of the process. */
   #ifdef MILPSolver_DEBUG
    // Initialize vector to store active variables indices
    std::vector< int > idxs_av;
   #endif
    
    if( row == 0 )
      matbeg[ row ] = 0;
    else
      matbeg[ row ] = matbeg[ row - 1 ] + matcnt[ row - 1 ];

    if( auto lf = dynamic_cast< const LinearFunction * >( f ) ) {
      matcnt[ row ] = lf->get_v_var().size();
      int j = 0;

      for( auto el : lf->get_v_var() ) {
        auto * v = dynamic_cast< ColVariable * >( std::get< 0 >( el ) );
        auto idx_v = index_of_variable( v ); 

        matval[ matbeg[ row ] + j ] = std::get< 1 >( el );;
        matind[ matbeg[ row ] + j ] = idx_v;

        #ifdef MILPSolver_DEBUG
          // Add the variable index to the safety vector
          idxs_av.push_back( idx_v );
        #endif

        j++;
      }
    }
    else if( auto qf = dynamic_cast< const QuadFunction * >( f ) ) {
      int nnz = 0;

      /* The quadratic part of the constraint will be represented as an
       * Eigen::SparseMatrix, as already done in QuadFunction. However, we
       * need to translate the local indices stored in a specific QuadFunction
       * into the global one of the model. */

      // Retrieve sparse quadratic matrix
      auto local_qmatrix = qf->get_matrix();
      std::map< mat_indices , float > global_qmatrix;
      //Qmat global_qmatrix;
      //global_qmatrix.reserve(5);
      //std::vector<Eigen::Triplet<Coefficient>> vv_nd( local_qmatrix.nonZeros() );

      // Create map from local indices to global ones
      std::vector< int > map_local_to_global( qf->get_num_active_var() );
      int num_var = 0;

      for( auto el : qf->get_v_var() ) {
        // Fill linear part of the constraint
        auto * v = dynamic_cast< ColVariable * >( std::get< 0 >( el ) );
        auto idx_v = index_of_variable( v );
        map_local_to_global[ num_var ] = idx_v;

        // If the linear coefficient is nonzero
        if( std::get< 1 >( el ) != 0 ) {
          matval[ matbeg[ row ] + nnz ] = std::get< 1 >( el );
          matind[ matbeg[ row ] + nnz ] = idx_v;
          nnz++;

          #ifdef MILPSolver_DEBUG
            // Add the variable index to the safety vector
            idxs_av.push_back( idx_v );
          #endif
        }

        // Note: the quadratic matrix does not contain the diagonal elements.
        // Thus, if there are nonzeros, we have to insert them.
        double q_coeff = std::get< 2 >( el );
        if( q_coeff != 0 ) {
          global_qmatrix[ { idx_v , idx_v } ] = q_coeff;
          //Eigen::Triplet< Coefficient > term( idx_v , idx_v , q_coeff );
          //vv_nd.push_back( term ); 
        }
        ++num_var;
      }
      matcnt[ row ] = nnz;

      // Now update local quadratic matrix into global one
      int k_term = 0;
      for (int k=0; k < local_qmatrix.outerSize(); ++k) {
        for( Qmat::InnerIterator it( local_qmatrix , k ) ; it ; ++it ) {
          int glob_idx1 = map_local_to_global[ it.row() ];
          int glob_idx2 = map_local_to_global[ it.col() ];
          global_qmatrix[ { glob_idx1 , glob_idx2 } ] = it.value();
          //Eigen::Triplet< Coefficient > term( glob_idx1 , glob_idx2 , it.value() );
          //vv_nd[ k_term ] = term;
          ++k_term;
        }
      }

      //global_qmatrix.setFromTriplets( vv_nd.begin(), vv_nd.end() );
      //global_qmatrix.makeCompressed();

      q_part[ row ] = global_qmatrix;
    }
    else if( auto dqf = dynamic_cast< const DQuadFunction * >( f ) ) {
      int nnz = 0;

      /* The quadratic part of the constraint will be represented as an
       * Eigen::SparseMatrix, as already done in QuadFunction. However, we
       * need to translate the local indices stored in a specific DQuadFunction
       * into the global one of the model. */

      // In this case we can simply insert the non-zero diagonal element
      std::map< mat_indices , float > global_qmatrix;
      //Qmat global_qmatrix;
      //std::vector<Eigen::Triplet<Coefficient>> vv_nd;

      for( auto el : dqf->get_v_var() ) {
        // Fill linear part of the constraint
        auto * v = dynamic_cast< ColVariable * >( std::get< 0 >( el ) );
        auto idx_v = index_of_variable( v );

        // If the linear coefficient is nonzero
        if( std::get< 1 >( el ) != 0 ) {
          matval[ matbeg[ row ] + nnz ] = std::get< 1 >( el );
          matind[ matbeg[ row ] + nnz ] = idx_v;
          nnz++;

          #ifdef MILPSolver_DEBUG
            // Add the variable index to the safety vector
            idxs_av.push_back( idx_v );
          #endif
        }

        // Check if the diagonal quadratic coefficient is nonzero
        double q_coeff = std::get< 2 >( el );
        if( q_coeff != 0 ) {
          global_qmatrix[ { idx_v , idx_v } ] = q_coeff;
          //Eigen::Triplet< Coefficient > term( idx_v , idx_v , q_coeff );
          //vv_nd.push_back( term ); 
        }
      }
      matcnt[ row ] = nnz;

      //global_qmatrix.setFromTriplets( vv_nd.begin(), vv_nd.end() );
      //global_qmatrix.makeCompressed();

      q_part[ row ] = global_qmatrix;
    }
    else
      throw( std::invalid_argument(
                  "MILPSolver::scan_constraint: "
                  "unexpected constraint type" ) );

    #ifdef MILPSolver_DEBUG
      // Now perform the actual sanity check of not having repeated variables
      // in the constraint

      // Sort the vector with respect to the variable indices
      std::sort( idxs_av.begin(), idxs_av.end() );

      // Now check there are no repeated indices
      if( idxs_av.size() > 0 ){
        for( int j = 0 ; j < idxs_av.size() - 1 ; ++j ){
          if( idxs_av[ j ] == idxs_av[ j + 1 ] ){
            // Error message
            std::string msg = std::string("MILPSolver Error [")
              + __func__ + "]: Variable with index " + 
              std::to_string( idxs_av[ j ] ) + 
              + " repeated multiple times in constraint with index " +
              std::to_string( row ) + ".\n"; 

            // Print warning message in MILPSolver DEBUG
            DEBUG_LOG( msg.c_str() );
          }
        }
      }
    #endif
   }
  }

 /* We need to define the sense of the constraint.
  * In SMS++ FRowConstraints are defined as:
  * LHS <= ( some function from Variables to reals ) <= RHS
  * Here we rather use the rhs + rngval approach */

 auto con_lhs = con.get_lhs();
 auto con_rhs = con.get_rhs();

 if( con.is_relaxed() ) {  // a relaxed constraint becomes: function >= -Inf
  sense[ row ] = 'L';
  rhs[ row ] = Inf< double >();
  }
 else
  if( con_lhs == con_rhs ) {
   // LHS <= function <= RHS, with LHS = RHS becomes: function = RHS
   sense[ row ] = 'E';
   rhs[ row ] = con_rhs;
   }
  else
   if( con_lhs == -Inf< double >() ) {
    // -inf <= function <= RHS becomes: function <= RHS
    sense[ row ] = 'L';
    rhs[ row ] = con_rhs;
    }
   else if( con_rhs == Inf< double >() ) {
    // LHS <= function <= inf becomes function >= LHS
    sense[ row ] = 'G';
    rhs[ row ] = con_lhs;
    }
   else {
    // LHS <= function <= RHS becomes: LHS <= function <= LHS + (range),
    // with range = RHS - LHS
    sense[ row ] = 'R';
    rhs[ row ] = con_lhs;
    rngval[ row ] = con_rhs - con_lhs;
    }
  ++row;
 }

/*--------------------------------------------------------------------------*/

void MILPSolver::scan_static_variable_bound( const ColVariable & var )
{
 // Get vector of bounds
 auto active_bounds = get_active_bounds( var , true );

 // This function should be called only if the parameter single_bound is
 // set to true. Thus, we have to check that maximum a single OneVarConstraint
 // is contained in the vector.
 if( active_bounds.size() > 1 )
   throw( std::logic_error(
              "MILPSolver::scan_static_variable_bound: only a single "
              "OneVarConstraint can be associated to a variable when "
              "the option intSingleBound is set to 1" ) );

 // Otherwise, if a single OneVarConstraint exists, we have to fill the
 // dictionary.
 if( active_bounds.size() == 1 )
  svar_to_bound.push_back( active_bounds[ 0 ] );
 else
  svar_to_bound.push_back( nullptr );
 }

/*--------------------------------------------------------------------------*/

void MILPSolver::scan_dynamic_variable_bound( const ColVariable & var )
{
 // Get vector of bounds
 auto active_bounds = get_active_bounds( var );

 // This function should be called only if the parameter single_bound is
 // set to true. Thus, we have to check that maximum a single OneVarConstraint
 // is contained in the vector.
 if( active_bounds.size() > 1 )
   throw( std::logic_error(
              "MILPSolver::scan_dynamic_variable_bound: only a single "
              "OneVarConstraint can be associated to a variable when "
              "the option intSingleBound is set to 1" ) );

 // Otherwise, if a single OneVarConstraint exists, we have to fill the 
 // dictionary.
 if( active_bounds.size() == 1 )
  dvar_to_bound.push_back( active_bounds[ 0 ] );
 else
  dvar_to_bound.push_back( nullptr );
 }

/*--------------------------------------------------------------------------*/

void MILPSolver::scan_objective( const FRealObjective * obj )
{
 // DEBUG_LOG( "MILPSolver::scan_objective() " << *obj );

 constant_value += obj->get_constant_term();

 // important note: a ColVariable may appear in multiple FRealObjective, this
 // meaning that all the corresponding terms will be present in the final
 // objective summed together (the Objective of a Block is the sum of its own
 // Objective plus the Objective of all its sub-Block, recursively); this is
 // why the "+=" below, each time a new linear or quadratic coefficient is
 // found it has to be *added* to ones already there (initialised to 0, hence
 // if the ColVariable appears only once this is equivalent to the fact that
 // the final coefficients are the ones found there)

 if( auto * lf = dynamic_cast< const LinearFunction * >(
						 obj->get_function() ) ) {
  /* Linear objective function */
  for( auto el : lf->get_v_var() )
   objective[ index_of_variable( el.first ) ] += el.second;

  return;
  }

 if( auto * qf = dynamic_cast< const QuadFunction * >(
						 obj->get_function() ) ) {
  /* Non separable Quadratic objective function */

  // Firstly, fill the linear and diagonal part of the
  // quadratic function.
  for( auto el : qf->get_v_var() ) {
   auto k = index_of_variable( std::get< 0 >( el ) );
   objective[ k ] += std::get< 1 >( el );
   q_objective[ k ] += std::get< 2 >( el );
   }

  // Now get sparse matrix related to the off diagonal terms
  v_off_diag_term nd_terms; // vector of tuples in the form (Idx,Idx,value)
  qf->get_v_nd_var( nd_terms );
  numnnzq = nd_terms.size();

  // Quadratic vectors allocation - - - - - - - - - -- - - - - - - - - - - -
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ndq_objective.resize( numnnzq, 0 );
  ndq_rowind.resize( numnnzq, 0 );
  ndq_colind.resize( numnnzq, 0 );

  for( int i = 0 ; i < numnnzq ; ++i  ) {
   auto el_nd = nd_terms[ i ];

   // Start from local idx and retrieve global index
   auto * v1 = qf->get_active_var( std::get< 0 >( el_nd ) );
   auto * v2 = qf->get_active_var( std::get< 1 >( el_nd ) );

   int global_idx1 = index_of_variable( dynamic_cast< ColVariable * >( v1 ) );
   int global_idx2 = index_of_variable( dynamic_cast< ColVariable * >( v2 ) );

   ndq_rowind[ i ] = global_idx1;
   ndq_colind[ i ] = global_idx2;
   ndq_objective[ i ] = std::get< 2 >( el_nd );
   }

  return;
  }

 if( auto * qf = dynamic_cast< const DQuadFunction * >(
						 obj->get_function() ) ) {
  /* Separable Quadratic objective function */
  for( auto el : qf->get_v_var() ) {
   auto k = index_of_variable( std::get< 0 >( el ) );
   objective[ k ] += std::get< 1 >( el );
   q_objective[ k ] += std::get< 2 >( el );
   }

  return;
  }

 throw( std::invalid_argument(
            "MILPSolver::scan_objective: "
            "unknown type of Objective Function" ) );
 }

/*--------------------------------------------------------------------------*/
/*----------------------------- MODIFICATIONS ------------------------------*/
/*--------------------------------------------------------------------------*/

int MILPSolver::compute( bool changedvars )
{
 lock();  // lock the Solver mutex

 // whatever certificate the previous solve left behind is stale
 f_dual_direction_value = std::numeric_limits< OFValue >::quiet_NaN();

 // separation may happen during this compute() either via the derived
 // Solver's cut callback (CutSepPar > 0) or via the explicit LP cut
 // separation loop driven here when intRelaxIntVars == 2; in both cases
 // we need write access to the Block (the cut callback can write the LP
 // solution into the Variable; the explicit loop does the same via
 // get_var_solution()), so a full lock() is acquired. Otherwise the
 // Block is only read_lock()-ed for the duration of process_modifications
 const bool may_separate = ( CutSepPar > 0 ) || ( relax_int_vars == 2 );

 bool owned = f_Block->is_owned_by( f_id );
 if( ! owned ) {
  if( may_separate ) {
   if( ! f_Block->lock( f_id ) )
    throw( std::runtime_error( "MILPSolver::compute: unable to lock Block" ) );
   }
  else
   if( ! f_Block->read_lock() )
    throw( std::runtime_error( "MILPSolver::compute: unable to read_lock Block"
                               ) );
  }

 MILPSolver::process_modifications();

 if( ( ! owned ) && ( ! may_separate ) )
  f_Block->read_unlock();  // no separation: release the read lock early

 // optional logging: gated by both a non-null log stream and the
 // SMS++-semantic intLogVerb >= 2 (captured at set_par() time as
 // log_verbosity, see comment in MILPSolver.h)
 const int log_verb = ( f_log ) ? log_verbosity : 0;

 // dispatch to the derived class for the actual back-end solve - - - - - - -
 int sts = guts_of_compute();

 // intRelaxIntVars == 2: explicit LP cut separation loop - - - - - - - - - -
 if( relax_int_vars == 2 ) {
  if( log_verb >= 2 )
   *f_log << "MILPSolver::compute: starting LP cut separation loop "
	  << "(intRelaxIntVars == 2, max " << max_cut_passes << " passes)"
	  << std::endl;

  // save and restore precision so we don't perturb the caller's stream
  auto savprec = ( log_verb >= 2 ) ? f_log->precision() : std::streamsize( 0 );

  for( int pass = 1 ; ( sts == kOK ) && ( pass < max_cut_passes ) ; ++pass ) {
   // write the LP solution into the Block Variable (so the separator can
   // see the right point); get_var_solution() with a nullptr Configuration
   // is the existing public API to do this
   get_var_solution( nullptr );

   // call the Block's separator; new dynamic Constraint reach this Solver
   // as BlockModAdd< FRowConstraint > queued in v_mod
   auto nM = v_mod.size();
   f_Block->generate_dynamic_constraints();
   auto new_mods = v_mod.size() - nM;

   if( log_verb >= 2 ) {
    f_log->precision( 10 );
    *f_log << "  pass " << pass << ": guts_of_compute -> sts = " << sts
	   << ", value = " << get_var_value()
	   << "; generate_dynamic_constraints -> " << new_mods
	   << " new Modification" << std::endl;
    }

   if( new_mods == 0 )
    break;                       // no new constraint: converged

   // process the new Modification (each will call the derived class's
   // add_dynamic_constraint() and update the back-end model)
   MILPSolver::process_modifications();

   // re-solve with the augmented model
   sts = guts_of_compute();
   }

  if( log_verb >= 2 ) {
   f_log->precision( 10 );
   *f_log << "MILPSolver::compute: LP cut separation loop ended, sts = "
	  << sts << ", value = " << get_var_value() << std::endl;
   f_log->precision( savprec );
   }
  }

 if( ( ! owned ) && may_separate )
  f_Block->unlock( f_id );       // release the write lock at the end

 unlock();  // unlock the Solver mutex
 return( sts );
 }

/*--------------------------------------------------------------------------*/

void MILPSolver::process_modifications( void )
{
 /* This function processes one modification after another, without any
  * attempt of optimization. Moreover, you have to be CAREFUL to write all
  * the cases in order from the most specialized to the more generic,
  * e.g., OneVarConstraintMod before RowConstraintMod before ConstraintMod,
  * otherwise the generic cases will intercept the more specialized ones. */

 for( ; ; )                 // process all the Modification loop
  if( auto mod = pop() ) {  // get next Modification, if any
   auto pmod = mod.get();   // down to regular Modification *
   if( dynamic_cast< const NBModification * >( pmod ) ) {
    f_reset = false;
    load_problem();         // an NBModification: reload everything
    mod_clear();            // all the remaining Modification must be ignored
    break;                  // all done
    }

   guts_of_process_modifications( pmod );  // process the Modification

   if( f_reset ) {
    f_reset = false;
    load_problem();         // reload everything
    mod_clear();            // all the remaining Modification must be ignored
    break;                  // all done
    }
   }
  else                      // no more Modification to process
   break;                   // all done

 #ifdef MILPSolver_DEBUG
  check_status();
 #endif

 }  // end( MILPSolver::process_modifications )

/*--------------------------------------------------------------------------*/

void MILPSolver::guts_of_process_modifications( const p_Mod mod )
{
 // silently drop Modification originating from a sub-Block the user has
 // asked to ignore via Solver::set_excluded_blocks(): is_excluded() walks
 // get_f_Block() chain of mod->get_Block(), so a Modification born deep
 // inside the excluded sub-tree is filtered out at any depth
 if( is_excluded( const_cast< Block * >( mod->get_Block() ) ) )
  return;

 // for a GroupModification, re-dispatch itself to all the sub-Modification
 if( auto gm = dynamic_cast< const GroupModification * >( mod ) ) {
  // DEBUG_LOG( "GroupModification containing:" << std::endl );
  for( const auto & submod : gm->sub_Modifications() )
   guts_of_process_modifications( submod.get() );
  return;
  }

 // for all other Modification, dispatch the appropriate virtual method
 if( auto vm = dynamic_cast< const VariableMod * >( mod ) ) {
  var_modification( vm );
  return;
  }

 if( auto om = dynamic_cast< const ObjectiveMod * >( mod ) ) {
  objective_modification( om );
  return;
  }

 if( auto bm = dynamic_cast< const OneVarConstraintMod * >( mod ) ) {
  bound_modification( bm );
  return;
  }

 if( auto tmod = dynamic_cast< const RowConstraintMod * >( mod ) ) {
  const_modification( tmod );
  return;
  }

 if( auto cm = dynamic_cast< const ConstraintMod * >( mod ) ) {
  const_modification( cm );
  return;
  }

 if( auto fm = dynamic_cast< const FunctionMod * >( mod ) ) {
  if( is_of( fm->function() ) )
   objective_function_modification( fm );
  else
   constraint_function_modification( fm );
  return;
  }

 if( auto fvm = dynamic_cast< const FunctionModVars * >( mod ) ) {
  if( is_of( fvm->function() ) )
   objective_fvars_modification( fvm );
  else
   constraint_fvars_modification( fvm );
  return;
  }

 if( auto dm = dynamic_cast< const BlockModAD * >( mod ) )
  dynamic_modification( dm );

 // PolyhedralFunctionMod (Addd / Rngd / Sbst) is intentionally not handled
 // here: the PolyhedralFunctionBlock owning the PolyhedralFunction reacts
 // to it inside its own add_Modification() and re-emits the equivalent
 // BlockModAdd<ColVariable> / BlockModRmv<ColVariable> /
 // BlockModAdd<FRowConstraint> + FunctionModVars{Addd,Rngd,Sbst} on the
 // affected obj_lf / normcns_lf / coupling[j] LinearFunctions, all of
 // which are picked up by the dispatch above. The PFB-level mod itself
 // is therefore safely ignored at the solver level
 // any other Modification is ignored

 }  // end( MILPSolver::guts_of_process_modifications )

/*--------------------------------------------------------------------------*/

bool MILPSolver::is_of( Function * f )
{
 return( dynamic_cast< Objective * >( f->get_Observer() ) );
 }

/*--------------------------------------------------------------------------*/

void MILPSolver::var_modification( const VariableMod * mod )
{
 auto var = static_cast< const ColVariable * >( mod->variable() );

 // update the number of integer variables
 if( ColVariable::is_integer( mod->old_state() ) !=
     ColVariable::is_integer( mod->new_state() ) ) {
  if( ColVariable::is_integer( mod->new_state() ) )
   ++int_vars;
  else    
   --int_vars;
  }

 if( lb.empty() && xctype.empty() )
  return;
 
 int idx = index_of_variable( var );

 // update bounds (if any)
 if( ! lb.empty() ) {
  auto bd = MILPSolver::get_problem_bounds( *var );
  lb[ idx ] = bd[ 0 ];
  ub[ idx ] = bd[ 1 ];
  }

 // update variable type (if any)
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
 } // end( MILPSolver::var_modification )

/*--------------------------------------------------------------------------*/

void MILPSolver::objective_modification( const ObjectiveMod * mod )
{
 switch( mod->type() ) {
  case( ObjectiveMod::eSetMin ): objsense = 1; break;
  case( ObjectiveMod::eSetMax ): objsense = -1; break;
  default:
   throw( std::invalid_argument(
              "MILPSolver::objective_modification: "
              "invalid type of ObjectiveMod" ) );
  }
 }

/*--------------------------------------------------------------------------*/

void MILPSolver::const_modification( const ConstraintMod * mod )
{
 // this Modification has nothing to do if these are empty
 if( rhs.empty() )
  return;

 auto * con = dynamic_cast< FRowConstraint * >( mod->constraint() );
 if( ! con )  // non-FRowConstraint Mods are silently ignored
  return;

 int idx = index_of_constraint( con );

 RowConstraint::RHSValue con_lhs = NAN;
 RowConstraint::RHSValue con_rhs = NAN;

 switch( mod->type() ) {
  case( ConstraintMod::eRelaxConst ):
   sense[ idx ] = 'G';
   rhs[ idx ] = -Inf< double >();
   rngval[ idx ] = 0;
   break;

  case( ConstraintMod::eEnforceConst ):
  case( RowConstraintMod::eChgLHS ):
  case( RowConstraintMod::eChgRHS ):
  case( RowConstraintMod::eChgBTS ):

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
   throw( std::invalid_argument(
              "MILPSolver::const_modification: "
              "invalid type of ConstraintMod" ) );
  }
 }  // end( MILPSolver::const_modification )

/*--------------------------------------------------------------------------*/

void MILPSolver::bound_modification( const OneVarConstraintMod * mod )
{
 // this modification has nothing to do if these are empty
 if( lb.empty() )
  return;

 auto con = static_cast< OneVarConstraint * >( mod->constraint() );
 auto var = static_cast< ColVariable * >( con->get_active_var( 0 ) );
 int idx = index_of_variable( var );

 switch( mod->type() ) {
  case( RowConstraintMod::eChgLHS ):
   lb[ idx ] = get_problem_lb( *var );
   break;
  case( RowConstraintMod::eChgRHS ):
   ub[ idx ] = get_problem_ub( *var );
   break;
  case( RowConstraintMod::eChgBTS ): {
   auto bd = MILPSolver::get_problem_bounds( *var );
   lb[ idx ] = bd[ 0 ];
   ub[ idx ] = bd[ 1 ];
   break;
   }
  default:
   throw( std::invalid_argument(
              "MILPSolver::bound_modification: "
              "invalid type of OneVarConstraintMod" ) );
  }
 }

/*--------------------------------------------------------------------------*/

void MILPSolver::objective_function_modification( const FunctionMod * mod )
{
 // this modification has nothing to do if these are empty
 if( objective.empty() )
  return;

 auto * f = mod->function();

 // important note: a ColVariable may appear in multiple FRealObjective, this
 // meaning that all the corresponding terms will be present in the final
 // objective summed together (the Objective of a Block is the sum of its own
 // Objective plus the Objective of all its sub-Block, recursively); this has
 // to be properly taken into account

 // C05FunctionModLin - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // we exploit the delta() vector of C05FunctionModLin, giving the difference
 // between the new and the old value of the linear coefficient, to update
 // the objective[] values without having to recompute them: since they are
 // (potentially) a sum of terms, recomputing them would require fetching
 // back all the terms, while the delta() can just be applied to the sum

 if( auto * modl = dynamic_cast< const C05FunctionModLin * >( mod ) ) {
  if( auto * lf = dynamic_cast< const LinearFunction * >( f ) ) {
   for( Block::Index i = 0 ; i < modl->vars().size() ; ++i ) {
    auto var = static_cast< const ColVariable * >( modl->vars()[ i ] );
    objective[ index_of_variable( var ) ] += modl->delta()[ i ];
    }
   return;
   }

  // if( const auto * qf = dynamic_cast< const DQuadFunction * > (f) ) {
  //
  //  // This may happen if we change from LP to QP
  //  if( q_objective.empty() ) {
  //   q_objective.resize( numcols );
  //  }
  //
  //  for( auto i : sbst->subset() ) {
  //   auto var = static_cast< const ColVariable * >( qf->get_active_var( i ) );
  //   auto idx = index_of_variable( var );
  //   objective[ idx ] = qf->get_linear_coefficient( i );
  //   q_objective[ idx ] = qf->get_quadratic_coefficient( i );
  //  }
  //
  //  return;
  // }

  // this should never happen
  throw( std::invalid_argument(
             "MILPSolver::objective_function_modification: "
             "unknown type of Objective Function" ) );
  }

 // fallback method - update all costs
 if( const auto * lf = dynamic_cast< const LinearFunction * >( f ) ) {
  // Linear objective function
  for( auto el : lf->get_v_var() )
   objective[ index_of_variable( el.first ) ] = el.second;

  return;
  }

 if( const auto * qf = dynamic_cast< const DQuadFunction * >( f ) ) {
  // quadratic objective function

  // this may happen if we change from LP to QP
  if( q_objective.empty() )
   q_objective.resize( numcols );

  for( auto el : qf->get_v_var() ) {
   int idx = index_of_variable( std::get< 0 >( el ) );
   objective[ idx ] = std::get< 1 >( el );
   q_objective[ idx ] = std::get< 2 >( el );
   }

  return;
  }

 // this should never happen
 throw( std::invalid_argument(
            "MILPSolver::objective_function_modification: "
            "unsupported type of Objective Function" ) );

 }  // end( MILPSolver::objective_function_modification )

/*--------------------------------------------------------------------------*/

void MILPSolver::constraint_function_modification( const FunctionMod * mod )
{
 // this modification has nothing to do if these are empty
 if( matval.empty() )
  return;

 auto * f = mod->function();

 auto * lf = dynamic_cast< const LinearFunction * >( f );
 if( ! lf )
  return;

 auto * con = dynamic_cast< const FRowConstraint * >( lf->get_Observer() );
 if( ! con )  // non-FRowConstraint Mods are silently ignored
  return;

 auto row = index_of_constraint( con );
 // Note: the per-backend overrides update the constraint matrix in place;
 // the base method only validates the lookup. The draft below illustrates
 // the shape of the implementation each backend specialises.
 // C05FunctionModLin
 // --------------------------------------------------------------------------
 // if( const auto * modl = dynamic_cast< C05FunctionModLin * >( mod ) ) {
 //
 //  for( int i = 0; i < modl->vars().size(); ++i ) {
 //   auto var = static_cast< const ColVariable * >( modl->vars()[ i ] );
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
            "MILPSolver::constraint_function_modification: "
            "not implemented yet" ) );

 }  // end( MILPSolver::constraint_function_modification )

/*--------------------------------------------------------------------------*/

void MILPSolver::objective_fvars_modification( const FunctionModVars * mod )
{
 // this modification has nothing to do if these are empty
 if( objective.empty() )
  return;

 auto * f = mod->function();

 // Check the modification type
 if( ( ! dynamic_cast< const C05FunctionModVarsAddd * >( mod ) ) &&
     ( ! dynamic_cast< const C05FunctionModVarsRngd * >( mod ) ) &&
     ( ! dynamic_cast< const C05FunctionModVarsSbst * >( mod ) ) )
  throw( std::invalid_argument(
             "MILPSolver::objective_fvars_modification: "
             "this type of FunctionModVars is not handled" ) );

 if( auto * lf = dynamic_cast< const LinearFunction * >( f ) ) {
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

 if( auto * qf = dynamic_cast< const DQuadFunction * >( f ) ) {
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
 throw( std::invalid_argument(
            "MILPSolver::objective_fvars_modification: "
            "unsupported type of Objective Function" ) );

 }  // end( MILPSolver::objective_fvars_modification )

/*--------------------------------------------------------------------------*/

void MILPSolver::constraint_fvars_modification( const FunctionModVars * mod )
{
 // this modification has nothing to do if these are empty
 if( matval.empty() )
  return;

 auto * lf = dynamic_cast< const LinearFunction * >( mod->function() );
 if( ! lf )
  return;

 // Check the modification type
 if( ( ! dynamic_cast< const C05FunctionModVarsAddd * >( mod ) ) &&
     ( ! dynamic_cast< const C05FunctionModVarsRngd * >( mod ) ) &&
     ( ! dynamic_cast< const C05FunctionModVarsSbst * >( mod ) ) )
  throw( std::invalid_argument(
             "MILPSolver::constraint_fvars_modification: "
             "this type of FunctionModVars is not handled" ) );

 // Get indices and coefficients
 auto * con = dynamic_cast< FRowConstraint * >( lf->get_Observer() );
 if( ! con )  // non-FRowConstraint Mods are silently ignored
  return;

 auto row = index_of_constraint( con );
 // Add/Remove of single coefficients within an existing row is not
 // wired in the base class; the per-backend overrides are expected to
 // handle it. Reaching here means a backend forgot to override.
 throw( std::logic_error(
            "MILPSolver::constraint_fvars_modification: "
            "not implemented in the base class" ) );

 }  // end( MILPSolver::constraint_fvars_modification )

/*--------------------------------------------------------------------------*/

void MILPSolver::dynamic_modification( const BlockModAD * mod )
{
 if( auto tmod = dynamic_cast< const BlockModAdd< FRowConstraint > * >( mod
									) ) {
  const auto & added = tmod->added();
  if( added.size() <= 1 ) {
   // single-row fast path: no need to materialise a vector copy when
   // there is at most one constraint to insert, and the derived
   // :MILPSolver's add_dynamic_constraints() default would just loop
   // over a length-1 vector to the same effect
   for( auto * i : added )
    add_dynamic_constraint( i );
   }
  else {
   // batch path: hand the whole list to add_dynamic_constraints() so
   // back-ends that expose a multi-row insertion API (CPLEX
   // CPXaddrows / Gurobi GRBaddconstrs / HiGHS Highs_addRows) can push
   // every row in a single call instead of paying the per-row overhead
   // of repeatedly touching their internal model. The typical trigger
   // is a PolyhedralFunctionModAddd retraduced by the owning
   // PolyhedralFunctionBlock into a single batched
   // add_dynamic_constraints(f_const, newc, ...)
   std::vector< const FRowConstraint * > batch;
   batch.reserve( added.size() );
   for( auto * i : added )
    batch.push_back( i );
   add_dynamic_constraints( batch );
   }

  return;
  }

 if( auto tmod = dynamic_cast< const BlockModRmv< FRowConstraint > * >( mod
									) ) {
  for( const auto & i : tmod->removed() )
   remove_dynamic_constraint( & i );

  return;
  }

 if( auto tmod = dynamic_cast< const BlockModAdd< ColVariable > * >( mod
								     ) ) {
  for( auto * i : tmod->added() )
   add_dynamic_variable( i );

  return;
  }

 if( auto tmod = dynamic_cast< const BlockModRmv< ColVariable > * >( mod
								    ) ) {
  for( const auto & i : tmod->removed() )
   remove_dynamic_variable( &i );

  return;
  }

 // all that has remained to do is to deal with OneVarConstraint
 // we deal with this using the base class BlockModAD, which is *not*
 // template, to be able to deal with all derived classes from
 // OneVarConstraint at once
 if( auto tmod = dynamic_cast< const BlockModAD * >( mod ) ) {
  if( mod->is_variable() )
   throw( std::invalid_argument(
              "MILPSolver::dynamic_modification: "
              "adding non-ColVariable not supported" ) );

  std::vector< Constraint * > mdcns;
  mod->get_elements( mdcns );
  if( mdcns.empty() )  // this should not happen, but ...
   return;             // there you go

  if( ! dynamic_cast< OneVarConstraint * >( mdcns.front() ) )
   throw( std::invalid_argument(
              "MILPSolver::dynamic_modification: "
              "adding unsupported Constraint type" ) );
  if( mod->is_added() )
   for( auto cnst : mdcns )
    add_dynamic_bound( static_cast< OneVarConstraint * >( cnst ) );
  else
   for( auto cnst : mdcns )
    remove_dynamic_bound( static_cast< OneVarConstraint * >( cnst ) );

  return;
  }

 throw( std::invalid_argument(
            "MILPSolver::dynamic_modification: "
            "unknown type of BlockModAD" ) );

 }  // end( MILPSolver::dynamic_modification )

/*--------------------------------------------------------------------------*/

void MILPSolver::add_dynamic_constraint( const FRowConstraint * con )
{
 // update the dictionaries
 auto it = lower_bound( dcon_to_idx.begin() , dcon_to_idx.end() , con ,
                        []( auto & pair , const FRowConstraint * c ) {
                         return( pair.first < c );
                         } );
 dcon_to_idx.insert( it , { con , numrows } );
 idx_to_dcon.emplace_back( con );

 // update the counter
 ++numrows;

 // update the vectors (if any)
 if( ! rhs.empty() ) {
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
  }

 // The CSC matrix (matbeg/matcnt/matind/matval) is the snapshot produced
 // by load_problem(); we do NOT mutate it here. Inserting a row into a
 // column-major CSC means scanning the LinearFunction of `con` and pushing
 // an entry into every involved column's slice, which is O(nnz_total) and
 // not worth the bookkeeping: the snapshot is only used to push the model
 // to the back-end at load_problem() time. After this point each derived
 // :MILPSolver (CPX/GRB/SCIP/HiGHS) holds its own internal model and
 // applies the row directly via its native API (CPXaddrows / GRBaddconstr
 // / SCIPcreateConsLinear / HighsAddRow) in its own override of
 // add_dynamic_constraint(). Therefore leaving the CSC stale is safe as
 // long as load_problem() is not invoked again on this Solver instance

 }  // end( MILPSolver::add_dynamic_constraint )

/*--------------------------------------------------------------------------*/

void MILPSolver::add_dynamic_variable( const ColVariable * var )
{
 // update the dictionaries
 auto it = lower_bound( dvar_to_idx.begin() , dvar_to_idx.end() , var ,
                        []( auto & pair , const ColVariable * v ) {
                         return( pair.first < v );
                         } );

 dvar_to_idx.insert( it , { var , numcols } );
 idx_to_dvar.emplace_back( var );

 // update the counters
 ++numcols;
 if( var->is_integer() )
  ++int_vars;

 // update the LB/UB vectors (if any)
 if( ! lb.empty() ) {
  auto bd = MILPSolver::get_problem_bounds( *var );
  lb.emplace_back( bd[ 0 ] );
  ub.emplace_back( bd[ 1 ] );
  }

 // update the objective vectors
 if( ! objective.empty() ) {
  objective.push_back( 0 );
  if( ! q_objective.empty() )
   q_objective.push_back( 0 );
  }

 // update the variable type vector
 if( ! xctype.empty() ) {
  if( var->is_integer() && ( ! relax_int_vars ) ) {
   if( var->is_unitary() && var->is_positive() )
    xctype.emplace_back( 'B' );
   else
    xctype.emplace_back( 'I' );
   }
  else
   xctype.emplace_back( 'C' );
  }

 /* If the check on SingleBound is active, scan the
 *  OneVarConstraint associated with the variable */
 if( single_bound )
  scan_dynamic_variable_bound( *var );

 // The CSC matrix (matbeg/matcnt/matind/matval) is left untouched: see the
 // comment in add_dynamic_constraint() for the rationale. The derived
 // :MILPSolver (CPX/GRB/SCIP/HiGHS) plumbs the new column into its own
 // back-end model via its override (CPXaddcols / GRBaddvar /
 // SCIPcreateVar / HighsAddCol). Leaving the CSC stale here is safe as
 // long as load_problem() is not invoked again on this Solver instance

 }  // end( MILPSolver::add_dynamic_variable )

/*--------------------------------------------------------------------------*/

void MILPSolver::add_dynamic_bound( const OneVarConstraint * con )
{
 // this modification has nothing to do if these are empty
 if( lb.empty() )
  return;

 auto var = static_cast< ColVariable * >( con->get_active_var( 0 ) );
 if( ! var )
  throw( std::logic_error(
             "MILPSolver::add_dynamic_bound: "
             "added a bound on no Variable" ) );

 int idx = index_of_variable( var );
 if( idx == Inf< int >() )
  throw( std::logic_error(
             "MILPSolver::add_dynamic_bound: "
             "added a bound on unknown Variable" ) );

 auto bd = MILPSolver::get_problem_bounds( *var );
 lb[ idx ] = bd[ 0 ];
 ub[ idx ] = bd[ 1 ];

 /* If the check on SingleBound is active, it is important to check that 
 *  no other OneVarConstraint are already associated with the variable.
 *  If this is the case, then add the new bound to the dictionary. */
 if( single_bound ) {
  if( idx < static_vars ) { // the variable is static
   if( svar_to_bound[ idx ] != nullptr ) // There was already a bound set
     throw( std::logic_error(
                "MILPSolver::add_dynamic_bound: only a single "
                "OneVarConstraint can be associated to a variable when "
                "the option intSingleBound is set to 1" ) );

   else // No OneVarConstraint was previously associated with the variable.
    svar_to_bound[ idx ] = con; // Update the dictionary
   }
  else { // the variable is dynamic
   if( dvar_to_bound[ idx - static_vars ] != nullptr ) // There was already a bound set
     throw( std::logic_error(
                "MILPSolver::add_dynamic_bound: only a single "
                "OneVarConstraint can be associated to a variable when "
                "the option intSingleBound is set to 1" ) );

   else // No OneVarConstraint was previously associated with the variable.
    dvar_to_bound[ idx - static_vars ] = con; // Update the dictionary
   }
  }
 }

/*--------------------------------------------------------------------------*/

void MILPSolver::remove_dynamic_constraint( const FRowConstraint * con )
{

 // update the dictionaries
 int index = 0;
 auto it1 = lower_bound( dcon_to_idx.begin() , dcon_to_idx.end() ,
                         std::make_pair( con , 0 ) ,
                         []( auto & p1 , auto & p2 ) {
                          return( p1.first < p2.first );
                          } );

 if( ( it1 != dcon_to_idx.end() ) && ( it1->first == con ) ) {
  index = it1->second;
  dcon_to_idx.erase( it1 );
  idx_to_dcon.erase( idx_to_dcon.begin() + ( index - static_cons ) );
  }
 else
  throw( std::runtime_error(
             "MILPSolver::remove_dynamic_constraint: "
             "dynamic constraint not found" ) );

 for( auto & it : dcon_to_idx )
  if( it.second > index )
   it.second--;

 // update the counter
 --numrows;

 // update the vectors, if any
 if( ! rhs.empty() ) {
  sense.erase( sense.begin() + index );
  rhs.erase( rhs.begin() + index );
  rngval.erase( rngval.begin() + index );
  }

 // update the matrix, if any
 if( matval.empty() )
  return;

 throw( std::logic_error(
            "MILPSolver::remove_dynamic_constraint: "
            "not fully implemented yet" ) );

 }  // end( MILPSolver::remove_dynamic_constraint )

/*--------------------------------------------------------------------------*/

void MILPSolver::remove_dynamic_variable( const ColVariable * var )
{

 // update the dictionaries
 int index = 0;
 auto it1 = lower_bound( dvar_to_idx.begin() , dvar_to_idx.end() ,
                         std::make_pair( var , 0 ),
                         []( auto & p1 , auto & p2 ) {
                          return( p1.first < p2.first );
                          } );

 if( ( it1 != dvar_to_idx.end() ) && ( it1->first == var ) ) {
  index = it1->second;
  dvar_to_idx.erase( it1 );
  idx_to_dvar.erase( idx_to_dvar.begin() + ( index - static_vars ) );
  if( single_bound )
    dvar_to_bound.erase( dvar_to_bound.begin() + ( index - static_vars ) );
  }
 else
  throw( std::runtime_error(
             "MILPSolver::remove_dynamic_variable: "
             "dynamic variable not found" ) );

 for( auto & it : dvar_to_idx )
  if( it.second > index )
   it.second--;

 // update the counters
 --numcols;
 if( var->is_integer() )
  --int_vars;

 // update the bound vectors, if any
 if( ! lb.empty() ) {
  lb.erase( lb.begin() + index );
  ub.erase( ub.begin() + index );
  }

 // update the objective vectors, if any
 if( ! objective.empty() ) {
  objective.erase( objective.begin() + index );
  if( ! q_objective.empty() )
   q_objective.erase( q_objective.begin() + index );
  }

 // update the matrix, if any
 if( matval.empty() )
  return;

 /* this would be easy, but the rest is not
 xctype.erase( xctype.begin() + index ); */

 throw( std::logic_error(
            "MILPSolver::remove_dynamic_variable: "
            "not fully implemented yet" ) );

 }  // end( MILPSolver::remove_dynamic_variable )

/*--------------------------------------------------------------------------*/

void MILPSolver::remove_dynamic_bound( const OneVarConstraint * con )
{
 // this Modification has nothing to do if these are empty
 if( lb.empty() )
  return;

 // note: this only works because remove_dynamic_constraint[ s ]() do *not*
 //       clear the removed OneVarConstraint, and therefore we can easily
 //       reconstruct which ColVariable it was about
 auto var = static_cast< const ColVariable * >( con->get_active_var( 0 ) );
 if( ! var )  // this should never happen
  return;     // but in case, there is nothing to do

 int idx = index_of_variable( var );
 if( idx == Inf< int >() )  // the ColVariable has been removed
  return;                   // is strange, but there is nothing to do
 
 auto bd = get_problem_bounds( * var );
 lb[ idx ] = bd[ 0 ];
 ub[ idx ] = bd[ 1 ];

 /* If the check on SingleBound is active, we have to remove the pointer to 
 *  the OneVarConstraint from the svar_to_bound or svar_to_bound dictionaries. 
 */
 if( single_bound ) {
  if( idx < static_vars ) // the variable is static
   svar_to_bound[ idx ] = nullptr;
  else // the variable is dynamic
   dvar_to_bound[ idx - static_vars ] = nullptr;
  }
 }

/*--------------------------------------------------------------------------*/
/*------------------------------- PARAMETERS -------------------------------*/
/*--------------------------------------------------------------------------*/

void MILPSolver::set_par( idx_type par , int value )
{
 if( par == intThrowReducedCostException ) {
  throw_reduced_cost_exception = bool( value );
  return;
  }
 if( par == intUseCustomNames ) {
  use_custom_names = bool( value );
  return;
  }
 if( par == intRelaxIntVars ) {
  relax_int_vars = value;
  return;
  }
 if( par == intSingleBound ) {
  single_bound = bool( value );
  return;
 }
 if( par == intConsModification ) {
  cons_modification = bool( value );
  return;
 }
 if( par == intCutSepPar ) {
  CutSepPar = (unsigned char)( value );
  return;
  }
 if( par == intMaxCutPasses ) {
  max_cut_passes = value;
  return;
  }

 CDASolver::set_par( par, value );
 }

/*--------------------------------------------------------------------------*/

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
 if( par == strWarmStartVariables ) {
  warmstart_variables = std::move( value );
  return;
  }
 if( par == strWarmStartSolution ) {
  warmstart_solution = std::move( value );
  return;
  }

 CDASolver::set_par( par, std::move( value ) );
 }

/*--------------------------------------------------------------------------*/

void MILPSolver::set_par( idx_type par , std::vector< std::string > && value )
{
 // MILPSolver itself defines no vstr_par; the path-based ignore list
 // previously here (vstrMILPIgnSBlks) is now a typed
 // std::unordered_set<Block*> installed via Solver::set_excluded_blocks()
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

ThinComputeInterface::idx_type MILPSolver::get_num_vstr_par( void ) const {
 return( CDASolver::get_num_vstr_par() + vstrLastAlgParMILP - vstrLastParCDAS );
 }

/*--------------------------------------------------------------------------*/

int MILPSolver::get_dflt_int_par( idx_type par ) const
{
 if( par == intThrowReducedCostException )
  return( 1 );

 if( par == intUseCustomNames )
  return( 1 );

 if( par == intRelaxIntVars )
  return( 0 );

 if( par == intSingleBound )
  return( 0 );

 if( par == intConsModification )
  return( 1 );

 if( par == intCutSepPar )
  return( 0 );

 if( par == intMaxCutPasses )
  return( 1000 );

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

 if( par == strWarmStartVariables )
  return( vals[ 1 ] );

 if( par == strWarmStartSolution )
  return( vals[ 1 ] );

 return( CDASolver::get_dflt_str_par( par ) );
 }

/*--------------------------------------------------------------------------*/

const std::vector< std::string > &
MILPSolver::get_dflt_vstr_par( idx_type par ) const
{
 // MILPSolver itself defines no vstr_par (see set_par(vstr) above)
 return( CDASolver::get_dflt_vstr_par( par ) );
 }

/*--------------------------------------------------------------------------*/

int MILPSolver::get_int_par( idx_type par ) const
{
 if( par == intThrowReducedCostException )
  return( throw_reduced_cost_exception );

 if( par == intUseCustomNames )
  return( use_custom_names );

 if( par == intRelaxIntVars )
  return( relax_int_vars );

 if( par == intSingleBound )
  return( single_bound );

 if( par == intConsModification )
  return( cons_modification );

 if( par == intCutSepPar )
  return( CutSepPar );

 if( par == intMaxCutPasses )
  return( max_cut_passes );

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

 if( par == strWarmStartVariables )
  return( warmstart_variables );

 if( par == strWarmStartSolution )
  return( warmstart_solution );

 return( CDASolver::get_str_par( par ) );
 }

/*--------------------------------------------------------------------------*/

const std::vector< std::string > &
MILPSolver::get_vstr_par( idx_type par ) const
{
 // MILPSolver itself defines no vstr_par (see set_par(vstr) above)
 return( CDASolver::get_vstr_par( par ) );
 }

/*--------------------------------------------------------------------------*/

Solver::idx_type MILPSolver::int_par_str2idx( const std::string & name ) const
{
 if( name == "intThrowReducedCostException" )
  return( intThrowReducedCostException );

 if( name == "intUseCustomNames" )
  return( intUseCustomNames );

 if( name == "intRelaxIntVars" )
  return( intRelaxIntVars );

 if( name == "intSingleBound" )
  return( intSingleBound );

 if( name == "intConsModification" )
  return( intConsModification );

 if( name == "intCutSepPar" )
  return( intCutSepPar );

 if( name == "intMaxCutPasses" )
  return( intMaxCutPasses );

 return( CDASolver::int_par_str2idx( name ) );
 }

/*--------------------------------------------------------------------------*/

const std::string & MILPSolver::int_par_idx2str( idx_type idx ) const
{
 static const std::vector< std::string > pars = { "intThrowReducedCostException",
                                                  "intUseCustomNames",
                                                  "intRelaxIntVars" ,
                                                  "intSingleBound" ,
                                                  "intConsModification" ,
                                                  "intCutSepPar" ,
                                                  "intMaxCutPasses" };
 if( idx == intThrowReducedCostException )
  return( pars[ 0 ] );

 if( idx == intUseCustomNames )
  return( pars[ 1 ] );

 if( idx == intRelaxIntVars )
  return( pars[ 2 ] );

 if( idx == intSingleBound )
  return( pars[ 3 ] );

 if( idx == intConsModification )
  return( pars[ 4 ] );

 if( idx == intCutSepPar )
  return( pars[ 5 ] );

 if( idx == intMaxCutPasses )
  return( pars[ 6 ] );

 return( CDASolver::int_par_idx2str( idx ) );
 }

/*----------------------------------------------------------------------------

Solver::idx_type MILPSolver::dbl_par_str2idx( const std::string & name ) const
{
 return( CDASolver::dbl_par_str2idx( name ) );
 }

------------------------------------------------------------------------------

const std::string & MILPSolver::dbl_par_idx2str( idx_type idx ) const {
 return( CDASolver::dbl_par_idx2str( idx ) );
 }

----------------------------------------------------------------------------*/

Solver::idx_type MILPSolver::str_par_str2idx( const std::string & name ) const
{
 if( name == "strProblemName" )
  return( strProblemName );

 if( name == "strOutputFile" )
  return( strOutputFile );

 if( name == "strWarmStartVariables" )
  return( strWarmStartVariables );

 if( name == "strWarmStartSolution" )
  return( strWarmStartSolution );

 return( CDASolver::str_par_str2idx( name ) );
 }

/*--------------------------------------------------------------------------*/

const std::string & MILPSolver::str_par_idx2str( idx_type idx ) const
{
 static const std::vector< std::string > pars = { "strProblemName" ,
                                                  "strOutputFile" ,
                                                  "strWarmStartVariables" ,
                                                  "strWarmStartSolution" };
 if( idx == strProblemName )
  return( pars[ 0 ] );

 if( idx == strOutputFile )
  return( pars[ 1 ] );

 if( idx == strWarmStartVariables )
  return( pars[ 2 ] );

 if( idx == strWarmStartSolution )
  return( pars[ 3 ] );

 return( CDASolver::str_par_idx2str( idx ) );
 }

/*--------------------------------------------------------------------------*/

Solver::idx_type
MILPSolver::vstr_par_str2idx( const std::string & name ) const
{
 // MILPSolver itself defines no vstr_par (see set_par(vstr) above)
 return( CDASolver::vstr_par_str2idx( name ) );
 }

/*--------------------------------------------------------------------------*/

const std::string & MILPSolver::vstr_par_idx2str( idx_type idx ) const
{
 // MILPSolver itself defines no vstr_par (see set_par(vstr) above)
 return( CDASolver::vstr_par_idx2str( idx ) );
 }

/*--------------------------------------------------------------------------*/
/*---------------------------------- DEBUG ---------------------------------*/
/*--------------------------------------------------------------------------*/

template< typename T >
std::string log_vector( const std::vector< T > & v , int limit ) {
 std::string temp_log = "[";
 for( auto i : v ) {
  temp_log += " " + std::to_string( i );
  if( --limit == 0 ) {
   temp_log += " ... ";
   break;
  }
 }
 temp_log += "]";
 return( temp_log );
 }

template<>
std::string log_vector( const std::vector< char > & v , int limit ) {
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
 return( temp_log );
 }

/*--------------------------------------------------------------------------*/

#ifdef MILPSolver_DEBUG

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
  throw( std::runtime_error(
             "MILPSolver::check_status: unable to lock the Block" ) );

 Q.push( f_Block );
 while( ! Q.empty() ) {
  Block * q_Block = Q.front();
  Q.pop();

  for( auto * i : q_Block->get_nested_Blocks() )
   Q.push( i );

  for( const auto & i : q_Block->get_static_constraints() ) {
   auto count = un_any_thing_count_static( FRowConstraint , i );
   if( count != Inf< std::size_t >() ) {
    if( count > 0 )
      ++scg;
    
    c += count;
    sc += count;
    continue;
   }

   // if it's not FRowConstraint, accept any known OneVarConstraint silently
   if( un_any_thing_OneVarConstraint_static( i , [](){}() ) )
    continue;

   throw( std::invalid_argument(
              "MILPSolver::check_status: static Constraint is neither "
              "FRowConstraint nor OneVarConstraint" ) );
  }

  for( const auto & i : q_Block->get_dynamic_constraints() ) {
   auto count = un_any_thing_count_dynamic( FRowConstraint , i );
   if( count != Inf< std::size_t >() ) {
    c += count;
    continue;
   }

   // if it's not FRowConstraint, accept any known OneVarConstraint silently
   if( un_any_thing_OneVarConstraint_dynamic( i , [](){}() ) )
    continue;

   throw( std::invalid_argument(
              "MILPSolver::check_status: dynamic Constraint is neither "
              "FRowConstraint nor OneVarConstraint" ) );
  }

  dc = c - sc;

  for( const auto & i : q_Block->get_static_variables() ) {
   auto count = un_any_thing_count_static( ColVariable , i );
   if( count == Inf< std::size_t >() )
    throw( std::invalid_argument(
               "MILPSolver::check_status: not a ColVariable" ) );

   if( count > 0 )
    ++svg;

    v += count;
   sv += count;
  }

  for( const auto & i : q_Block->get_dynamic_variables() ) {
   auto count = un_any_thing_count_dynamic( ColVariable , i );
   if( count == Inf< std::size_t >() )
    throw( std::invalid_argument(
               "MILPSolver::check_status: not a ColVariable" ) );
   v += count;
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
                          return( ( std::get< 1 >( pair ) == i.first ) &&
                                  ( std::get< 0 >( pair ) == i.second ) );
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
                          return( ( std::get< 0 >( i ) == pair.second ) &&
                                  ( std::get< 1 >( i ) == pair.first ) );
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
                         [ & ]( auto & pair ) { return( pair.first == i ); } );
  if( j == dvar_to_idx.end() ) {
   DEBUG_LOG( "Element [" << i
                          << "] of idx_to_dvar was not found in dvar_to_idx"
                          << std::endl );
  }
 }

 for( auto & i: dvar_to_idx ) {
  auto j = std::find_if( idx_to_dvar.begin(), idx_to_dvar.end(),
                         [ & ]( auto & var ) { return( i.first == var ); } );
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
                          return( ( std::get< 1 >( pair ) == i.first ) &&
                                  ( std::get< 0 >( pair ) == i.second ) );
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
                          return( ( std::get< 0 >( i ) == pair.second ) &&
                                  ( std::get< 1 >( i ) == pair.first ) );
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
                         [ & ]( auto & pair ) { return( pair.first == i ); } );
  if( j == dcon_to_idx.end() ) {
   DEBUG_LOG( "Element [" << i
                          << "] of idx_to_dcon was not found in dcon_to_idx"
                          << std::endl );
  }
 }
 for( auto & i: dcon_to_idx ) {
  auto j = std::find_if( idx_to_dcon.begin(), idx_to_dcon.end(),
                         [ & ]( auto & con ) { return( i.first == con ); } );
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

#endif // MILPSolver_DEBUG

/*--------------------------------------------------------------------------*/

void MILPSolver::write_var_solution( const std::vector< double > & x )
{
 if( x.empty() )  // actually no solution passed
  return;         // silently return

 if( x.size() < get_numcols() )
  throw( std::invalid_argument(
             "MILPSolver::write_var_solution: x too short" ) );

 int col = 0;

 auto set = [ & x , & col ]( ColVariable & v ) {
  v.set_value( x[ col++ ] );
  };

 for( auto qb : v_BFS ) {
  for( const auto & vi : qb->get_static_variables() )
   un_any_const_static( vi , set , un_any_type< ColVariable >() );
  }

 // Dynamic columns are appended to the solver in modification-arrival order,
 // which is generally different from the block/BFS order above (for example,
 // when bundle cuts are added to different PolyhedralFunctionBlock-s over
 // time). idx_to_dvar is maintained in the actual solver-column order by
 // add_dynamic_variable() / remove_dynamic_variable(), so it is the only
 // reliable map for writing the dynamic part of the solution back.
 for( std::size_t i = 0 ; i < idx_to_dvar.size() ; ++i )
  const_cast< ColVariable * >( idx_to_dvar[ i ] )->set_value(
                                                x[ static_vars + i ] );
 }  // end( MILPSolver::write_var_solution )

/*--------------------------------------------------------------------------*/

void MILPSolver::write_dual_solution( const std::vector< double > & pi ,
				      const std::vector< double > & rc )
{
 // handle dual variables, if any- - - - - - - - - - - - - - - - - - - - - -
 if( ! pi.empty() ) {  // there actually is a dual solution

  if( pi.size() < get_numrows() )
   throw( std::invalid_argument(
              "MILPSolver::write_dual_solution: pi too short" ) );

  // NOTE: this only supports pi written for linear constraints!
  // TODO: extend pi to quadratic constraints
  int row = 0;

  auto set = [ & pi , & row ]( FRowConstraint & c ) {
    if( dynamic_cast< LinearFunction * >( c.get_function() ) )
      c.set_dual( - pi[ row++ ] );
    else
      // Skip quadratic rows
      row++;
   };

  /* The dynamic rows are appended to the solver in modification-arrival
   * order, which is generally different from the block/BFS order the static
   * ones are written in: a running counter would therefore scramble the
   * duals across the dynamic groups, exactly as it did for the dynamic
   * columns [see write_var_solution()]. The row of each of them is asked
   * for instead. */

  auto set_dynamic = [ this , & pi ]( FRowConstraint & c ) {
    if( ! dynamic_cast< LinearFunction * >( c.get_function() ) )
     return;   // skip quadratic rows

    const int row = index_of_constraint( & c );
    if( row < int( get_numrows() ) )
     c.set_dual( - pi[ row ] );
   };

  for( auto qb : v_BFS ) {
   for( const auto & ci : qb->get_static_constraints() )
    un_any_const_static( ci , set , un_any_type< FRowConstraint >() );

   for( const auto & ci : qb->get_dynamic_constraints() )
    un_any_const_dynamic( ci , set , un_any_type< FRowConstraint >() );
   }
  }  // end( ! p.empty() )

 // handle reduced costs, if any - - - - - - - - - - - - - - - - - - - - - -
 if( rc.empty() )  // actually no reduced costs passed
  return;         // all done

 if( rc.size() < get_numcols() )
  throw( std::invalid_argument(
             "MILPSolver::write_dual_solution: rc too short" ) );

 int col = 0;

 auto set_bound = [ this , & rc , & col ]( ColVariable & var ) {
  auto active_bounds = get_active_bounds( var );

  // bounds that will have the dual value set.
  OneVarConstraint * lhs_con = nullptr;
  OneVarConstraint * rhs_con = nullptr;

  auto var_lb = var.get_lb();
  auto var_ub = var.get_ub();

  const auto var_is_fixed = var.is_fixed();
  if( var_is_fixed ) {
   /* The Variable is fixed. There should be at least one OneVarConstraint
    * (for this Variable) whose lower and upper bounds are equal to the value
    * of this Variable. If such OneVarConstraint exists, the reduced cost of
    * this Variable will be dual of that OneVarConstraint. If there is no such
    * OneVarConstraint, the reduced cost of this variable will be lost. */
   var_lb = var.get_value();
   var_ub = var.get_value(); //REMOVE

   for( auto b : active_bounds ) {
    b->set_dual( 0 );
    if( ( b->get_lhs() == var_lb ) && ( b->get_rhs() == var_lb ) ) {
     lhs_con = b;
     rhs_con = b;
     }
    }

   assert( lhs_con == rhs_con );
   }
  else {  // a non-fixed Variable
   for( auto b: active_bounds ) {
    b->set_dual( 0 );
 
    if( b->get_lhs() >= var_lb ) {
     var_lb = b->get_lhs();
     lhs_con = b;
     }

    if( b->get_rhs() <= var_ub ) {
     var_ub = b->get_rhs();
     rhs_con = b;
     }
    }
   }

  if( lhs_con && ( rc[ col ] >= 0 ) )
   lhs_con->set_dual( - rc[ col ] );
  else
   if( rhs_con && ( rc[ col ] <= 0 ) )
    rhs_con->set_dual( - rc[ col ] );
   else
    // interior-point solutions carry sign noise up to the dual feasibility
    // tolerance on reduced costs: a tiny value of the "wrong" sign for the
    // only existing bound means 0, which the duals of all the active bounds
    // have just been set to, so only complain on a significant mismatch
    if( ( lhs_con || rhs_con ) && ( std::abs( rc[ col ] ) > 1e-6 ) )
     throw( std::logic_error(
	       "MILPSolver::write_dual_solution: invalid dual value" ) );

  if( throw_reduced_cost_exception ) {
    if( var_is_fixed && ( ! lhs_con ) && ( var_lb != 0 ) ) {
     // the Variable is fixed, but it has no associated OneVarConstraint
     // with both bounds equal to the value of the Variable

     throw( std::logic_error(
      "MILPSolver::write_dual_solution: variable with index " +
      std::to_string( col ) + " is fixed to " +
      std::to_string( var.get_value() ) + ", but it has no OneVarConstraint" +
      "with both bounds equal to the value of this variable." ) );
     }
    else
     if( ( ! var_is_fixed ) && ( ! lhs_con ) && ( ! rhs_con ) ) {
      // the Variable is not fixed, and it has no associated OneVarConstraint
      // an exception is thrown if it has a finite nonzero bound
      if( ( ( var_lb != 0 ) && ( std::abs( var_lb ) < Inf< double >() ) ) ||
	  ( ( var_ub != 0 ) && ( std::abs( var_ub ) < Inf< double >() ) ) )
       throw( std::logic_error(
                "MILPSolver::get_dual_solution: variable with index " +
		 std::to_string( col ) + " has no OneVarConstraint." ) );
      }
    }

   col += 1;  // update variable counter
   };

 for( auto qb : v_BFS ) {
  // write all static variables first
  for( const auto & vi : qb->get_static_variables() )
   un_any_const_static( vi , set_bound , un_any_type< ColVariable >() );
  // write all dynamic variables
  for( const auto & vi : qb->get_dynamic_variables() )
   un_any_const_dynamic( vi , set_bound , un_any_type< ColVariable >() );
  }
 }  // end( MILPSolver::write_dual_solution )

/*--------------------------------------------------------------------------*/
/*----------------------- End File MILPSolver.cpp --------------------------*/
/*--------------------------------------------------------------------------*/ 
