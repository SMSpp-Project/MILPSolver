/*--------------------------------------------------------------------------*/
/*-------------------------- File MILPSolver.cpp ---------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Implementation of the MILPSolver class.
 *
 * \author Antonio Frangioni \n
 *         Operations Research Group \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \author Niccolò Iardella \n
 *         Operations Research Group \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \author Kostas Tavlaridis-Gyparakis \n
 *         Operations Research Group \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * Copyright &copy by Antonio Frangioni, Kostas Tavlaridis-Gyparakis, Niccolò Iardella
 */
/*--------------------------------------------------------------------------*/
/*---------------------------- IMPLEMENTATION ------------------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include <functional>
#include <queue>

#include <boost/log/core.hpp>
#include <boost/log/trivial.hpp>
#include <boost/log/expressions.hpp>

#include <Block.h>
#include <OneVarConstraint.h>
#include <LinearFunction.h>
#include <DQuadFunction.h>

#include "MILPSolver.h"

/*--------------------------------------------------------------------------*/
/*------------------------- NAMESPACE AND USING ----------------------------*/
/*--------------------------------------------------------------------------*/

using namespace SMSpp_di_unipi_it;

SMSpp_insert_in_factory_cpp_0( MILPSolver );

/*--------------------------------------------------------------------------*/
/*--------------------- CONSTRUCTOR AND DESTRUCTOR -------------------------*/
/*--------------------------------------------------------------------------*/

MILPSolver::MILPSolver() : CDASolver() {
#ifdef MILPSLVR_DEBUG
 boost::log::core::get()->set_filter(
  boost::log::trivial::severity >= boost::log::trivial::debug
 );
#else
 boost::log::core::get()->set_filter(
  boost::log::trivial::severity >= boost::log::trivial::info
 );
#endif
}

MILPSolver::~MILPSolver() {
 for( auto & i: colname )
  delete i;
 for( auto & i: rowname )
  delete i;
}

/*--------------------------------------------------------------------------*/
/*------------------------------- GETTERS ----------------------------------*/
/*--------------------------------------------------------------------------*/

int MILPSolver::get_numcols() const {
 return numcols;
}

int MILPSolver::get_numrows() const {
 return numrows;
}

int MILPSolver::get_nzelements() const {
 return nzelements;
}

int MILPSolver::get_objsense() const {
 return objsense;
}

const std::vector< double > & MILPSolver::get_objective() const {
 return objective;
}

const std::vector< double > & MILPSolver::get_q_objective() const {
 return q_objective;
}

const std::vector< double > & MILPSolver::get_rhs() const {
 return rhs;
}

const std::vector< double > & MILPSolver::get_rngval() const {
 return rngval;
}

const std::vector< char > & MILPSolver::get_sense() const {
 return sense;
}

const std::vector< int > & MILPSolver::get_matbeg() const {
 return matbeg;
}

const std::vector< int > & MILPSolver::get_matcnt() const {
 return matcnt;
}

const std::vector< int > & MILPSolver::get_matind() const {
 return matind;
}

const std::vector< double > & MILPSolver::get_matval() const {
 return matval;
}

const std::vector< double > & MILPSolver::get_var_lb() const {
 return lb;
}

const std::vector< double > & MILPSolver::get_var_ub() const {
 return ub;
}

const std::vector< char > & MILPSolver::get_xctype() const {
 return xctype;
}

const std::vector< char* > & MILPSolver::get_rowname() const {
 return rowname;
}

const std::vector< char* > & MILPSolver::get_colname() const {
 return colname;
}

int MILPSolver::get_nodes() const {
 return nodes;
}

int MILPSolver::get_num_integer_vars() const {
 return int_vars;
}

/*--------------------------------------------------------------------------*/
/*-------------------------------- SET_BLOCK -------------------------------*/
/*--------------------------------------------------------------------------*/

void MILPSolver::set_Block( Block * block ) {

 Solver::set_Block( block );

 if( block) {

  bool owned = block->is_owned_by( f_id );
  if( !owned && !block->lock( f_id ) ) {
   throw std::runtime_error( "Unable to lock the Block" );
  }

  // Generate abstract representation
  block->generate_abstract_variables( nullptr );
  block->generate_abstract_constraints( nullptr );
  block->generate_objective( nullptr );

  if( !owned ) {
   block->unlock( f_id );
  }

  clear_problem();
  load_problem();
 }
}

/*--------------------------------------------------------------------------*/

void MILPSolver::clear_problem() {
 numrows = 0;
 numcols = 0;
 nzelements = 0;
 objsense = 0;

 matbeg.clear();
 matcnt.clear();
 matind.clear();
 matval.clear();
 rhs.clear();
 rngval.clear();
 sense.clear();
 objective.clear();
 q_objective.clear();
 lb.clear();
 ub.clear();
 xctype.clear();

 for( auto & i: colname )
  delete i;
 for( auto & i: rowname )
  delete i;
 colname.clear();
 rowname.clear();

 svar_to_idx.clear();
 idx_to_svar.clear();
 scon_to_idx.clear();
 idx_to_scon.clear();
 dvar_to_idx.clear();
 idx_to_dvar.clear();
 dcon_to_idx.clear();
 idx_to_dcon.clear();

 for( auto i: active_constraints )
  i.clear();
 active_constraints.clear();
 for( auto i: active_bounds )
  i.clear();
 active_bounds.clear();
}

/*--------------------------------------------------------------------------*/

void MILPSolver::load_problem() {
 /*
  * Passing all the data of the Block to the LP.
  *
  * Following a Breadth First Search we proceed with scanning the received
  * Block and all of each corresponding children if any, in order to populate
  * the LP data.
  */
 std::queue< Block * > Q;

 // Locking the Block
 bool owned = f_Block->is_owned_by( f_id );
 if( !owned && !f_Block->read_lock() ) {
  throw std::runtime_error( "Unable to lock the Block" );
 }

 // First loop on the queue to count variables and constraints
 Q.push( f_Block );

 int var = 0;       // Variable counter for ncount_nzelements()
 int num_block = 0; // Counter for the blocks

 while( !Q.empty() ) {
  Block * q_Block = Q.front();
  Q.pop();
  BOOST_LOG_TRIVIAL( debug ) << "Processing Block " << num_block << " ["
                             << q_Block << "]";
  BOOST_LOG_TRIVIAL( trace ) << *q_Block;

  for( auto * i : q_Block->get_nested_Blocks() ) {
   Q.push( i );
  }

  BOOST_LOG_TRIVIAL( trace )
   << "MILPSolver::set_Block() counting static constraints";
  for( const auto & i : q_Block->get_static_constraints() ) {
   auto f1 = std::bind( &MILPSolver::count_constraints,
                        this,
                        std::placeholders::_1,
                        std::ref( numrows ) );
   un_any_const_static( i, f1, un_any_type< FRowConstraint >() );
  }
  static_cons = numrows;

  BOOST_LOG_TRIVIAL( trace )
   << "MILPSolver::set_Block() counting dynamic constraints";
  for( const auto & i : q_Block->get_dynamic_constraints() ) {
   auto f1 = std::bind( &MILPSolver::count_constraints,
                        this,
                        std::placeholders::_1,
                        std::ref( numrows ) );
   un_any_const_dynamic( i, f1, un_any_type< FRowConstraint >() );
  }

  BOOST_LOG_TRIVIAL( trace )
   << "MILPSolver::set_Block() counting static variables";
  for( const auto & i : q_Block->get_static_variables() ) {
   auto f1 = std::bind( &MILPSolver::count_variables,
                        this,
                        std::placeholders::_1,
                        std::ref( numcols ) );
   un_any_const_static( i, f1, un_any_type< ColVariable >() );
  }
  static_vars = numcols;

  BOOST_LOG_TRIVIAL( trace )
   << "MILPSolver::set_Block() counting dynamic variables";
  for( const auto & i : q_Block->get_dynamic_variables() ) {
   auto f1 = std::bind( &MILPSolver::count_variables,
                        this,
                        std::placeholders::_1,
                        std::ref( numcols ) );
   un_any_const_dynamic( i, f1, un_any_type< ColVariable >() );
  }

  BOOST_LOG_TRIVIAL( trace ) << "MILPSolver::set_Block() nonzero elements";
  for( const auto & i : q_Block->get_static_variables() ) {
   auto f1 = std::bind( &MILPSolver::count_nzelements,
                        this,
                        std::placeholders::_1,
                        std::ref( nzelements ),
                        std::ref( var ) );
   un_any_const_static( i, f1, un_any_type< ColVariable >() );
  }

  for( const auto & i : q_Block->get_dynamic_variables() ) {
   auto f1 = std::bind( &MILPSolver::count_nzelements,
                        this,
                        std::placeholders::_1,
                        std::ref( nzelements ),
                        std::ref( var ) );
   un_any_const_dynamic( i, f1, un_any_type< ColVariable >() );
  }
  ++num_block;
 } // End of while loop on Block queue

 BOOST_LOG_TRIVIAL( debug ) << "numrows (constraints) = " << numrows;
 BOOST_LOG_TRIVIAL( debug ) << "numcols (variables)   = " << numcols;
 BOOST_LOG_TRIVIAL( debug ) << "nzelements            = " << nzelements;
 BOOST_LOG_TRIVIAL( debug ) << "static constraints    = " << static_cons;
 BOOST_LOG_TRIVIAL( debug ) << "static variables      = " << static_vars;

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
 colname.resize( numcols, nullptr );
 rowname.resize( numrows, nullptr );

 // Second loop to scan the constraints
 Q.push( f_Block );
 num_block = 0; // Counter for the blocks
 int row = 0;   // Counter for the rows

 while( !Q.empty() ) {
  Block * q_Block = Q.front();
  Q.pop();

  for( auto i : q_Block->get_nested_Blocks() ) {
   Q.push( i );
  }

  /*
   * Scanning and passing all the data of the examined block
   * to the corresponding CPLEX data. We do this by first scanning
   * the static part of the problem and then the dynamic part.
   * This is because we want to have an order in the columns
   * and rows of the CPLEX coeff matrix where the static part is
   * being followed by the dynamic one
   */

  int set = 0; // Counter for the constraint groups

  for( const auto & i : q_Block->get_static_constraints() ) {
   int elements = 0; // Counter for group elements
   int start = row;

   auto f1 = std::bind( &MILPSolver::scan_constraint,
                        this,
                        std::placeholders::_1,
                        std::ref( elements ),
                        std::ref( row ) );
   un_any_const_static( i, f1, un_any_type< FRowConstraint >() );

   // Write names
   auto base = q_Block->get_s_const_name()[ set ];
   int end = row - start;
   for( int n = 0; n < end; ++n ) {
    std::string name;
    if( base.empty() ) {
     name = "cs_" + std::to_string( num_block )
            + "_" + std::to_string( set )
            + "_" + std::to_string( n );
    } else {
     name = base
            + "_" + std::to_string( num_block )
            + "_" + std::to_string( n );
    }
    rowname[ start + n ] = strcpy( new char[name.length() + 1], name.c_str() );
   }
   set++;
   if( elements ) {
    std::get< 2 >( scon_to_idx.back() ) = elements;
   }
  }

  set = 0;
  for( const auto & i : q_Block->get_dynamic_constraints() ) {
   int start = row;

   auto f1 = std::bind( &MILPSolver::scan_constraint,
                        this,
                        std::placeholders::_1,
                        -1,
                        std::ref( row ) );
   un_any_const_dynamic( i, f1, un_any_type< FRowConstraint >() );

   // Write names
   auto base = q_Block->get_d_const_name()[ set ];
   int end = row - start;
   for( int n = 0; n < end; ++n ) {
    std::string name;
    if( base.empty() ) {
     name = "cd_" + std::to_string( num_block )
            + "_" + std::to_string( set )
            + "_" + std::to_string( n );
    } else {
     name = base
            + "_" + std::to_string( num_block )
            + "_" + std::to_string( n );
    }
    rowname[ start + n ] = strcpy( new char[name.length() + 1], name.c_str() );
   }
   set++;
  }
  num_block++;
 } // End of while loop on Block queue

 // Keep the dictionaries sorted
 std::sort( scon_to_idx.begin(), scon_to_idx.end() );
 std::sort( dcon_to_idx.begin(), dcon_to_idx.end() );

 // These are already sorted at this point
 // std::sort( v_int_s_const.begin(), v_int_s_const.end() );
 // std::sort( v_int_d_const.begin(), v_int_d_const.end() );

 // Third loop to scan the variables
 Q.push( f_Block );
 num_block = 0;
 int col = 0; // Counter for the columns

 while( !Q.empty() ) {
  Block * q_Block = Q.front();
  Q.pop();

  for( auto *i : q_Block->get_nested_Blocks() ) {
   Q.push( i );
  }

  int set = 0;   // Counter for the variable groups

  for( const auto & i : q_Block->get_static_variables() ) {
   int elements = 0; // Counter for group elements
   int start = col;

   auto f1 = std::bind( &MILPSolver::scan_variable,
                        this,
                        std::placeholders::_1,
                        std::ref( elements ),
                        std::ref( col ) );
   un_any_const_static( i, f1, un_any_type< ColVariable >() );

   // Write names
   auto base = q_Block->get_s_var_name()[ set ];
   int end = col - start;
   for( int n = 0; n < end; ++n ) {
    std::string name;
    if( base.empty() ) {
     name = "xs_" + std::to_string( num_block )
            + "_" + std::to_string( set )
            + "_" + std::to_string( n );
    } else {
     name = base
            + "_" + std::to_string( num_block )
            + "_" + std::to_string( n );
    }
    colname[ start + n ] = strcpy( new char[name.length() + 1], name.c_str() );
   }
   set++;
   if( elements ) {
    std::get< 2 >( svar_to_idx.back() ) = elements;
   }
  }

  set = 0;
  for( const auto & i : q_Block->get_dynamic_variables() ) {
   int start = col;

   auto f1 = std::bind( &MILPSolver::scan_variable,
                        this,
                        std::placeholders::_1,
                        -1,
                        std::ref( col ) );
   un_any_const_dynamic( i, f1, un_any_type< ColVariable >() );

   // Write names
   auto base = q_Block->get_d_var_name()[ set ];
   int end = col - start;
   for( int n = 0; n < end; ++n ) {
    std::string name;
    if( base.empty() ) {
     name = "xv_" + std::to_string( num_block )
            + "_" + std::to_string( set )
            + "_" + std::to_string( n );
    } else {
     name = base
            + "_" + std::to_string( num_block )
            + "_" + std::to_string( n );
    }
    colname[ start + n ] = strcpy( new char[name.length() + 1], name.c_str() );
   }
   set++;
  }
  num_block++;
 } // End of while loop on Block queue

 // Keep the dictionaries sorted
 std::sort( svar_to_idx.begin(), svar_to_idx.end() );
 std::sort( dvar_to_idx.begin(), dvar_to_idx.end() );

 // These are already sorted at this point
 // std::sort( v_int_s_var.begin(), v_int_s_var.end() );
 // std::sort( v_int_d_var.begin(), v_int_d_var.end() );

 switch( f_Block->get_objective_sense() ) {
  case ( Objective::eMax ):
   objsense = -1;
   break;
  case ( Objective::eMin ):
   objsense = 1;
   break;
  default:
   objsense = 0;
   break;
 }

 // Fourth loop to scan the objective(s?)
 Q.push( f_Block );
 while( !Q.empty() ) {
  Block * q_Block = Q.front();
  Q.pop();

  for( auto *i : q_Block->get_nested_Blocks() ) {
   Q.push( i );
  }

  auto * p_obj = dynamic_cast< FRealObjective * >( q_Block->get_objective() );
  if( p_obj ) {
   scan_objective( p_obj );
  }
 } // End of while loop on Block queue

 BOOST_LOG_TRIVIAL(debug) << "objective   = " << log_vector(objective);
 BOOST_LOG_TRIVIAL(debug) << "q_objective = " << log_vector(q_objective );
 BOOST_LOG_TRIVIAL(debug) << "rhs         = " << log_vector(rhs );
 BOOST_LOG_TRIVIAL(debug) << "rngval      = " << log_vector(rngval );
 BOOST_LOG_TRIVIAL(debug) << "sense       = " << log_vector(sense );
 BOOST_LOG_TRIVIAL(debug) << "matbeg      = " << log_vector(matbeg );
 BOOST_LOG_TRIVIAL(debug) << "matcnt      = " << log_vector(matcnt );
 BOOST_LOG_TRIVIAL(debug) << "matind      = " << log_vector(matind );
 BOOST_LOG_TRIVIAL(debug) << "matval      = " << log_vector(matval );
 BOOST_LOG_TRIVIAL(debug) << "lb          = " << log_vector(lb );
 BOOST_LOG_TRIVIAL(debug) << "ub          = " << log_vector(ub );
 BOOST_LOG_TRIVIAL(debug) << "xctype      = " << log_vector(xctype );

 // Unlock the Block
 if( !owned ) {
  f_Block->read_unlock();
 }
}

/*--------------------------------------------------------------------------*/
/*-------------------- METHODS FOR PROBLEM DESCRIPTION ---------------------*/
/*--------------------------------------------------------------------------*/

int MILPSolver::index_of_variable( const ColVariable * p_var ) {
 auto i = index_of_static_variable( p_var );
 if( i < Inf< int >() ) {
  return i;
 } else {
  return index_of_dynamic_variable( p_var );
 }
}

/*--------------------------------------------------------------------------*/

int MILPSolver::index_of_static_variable( const ColVariable * p_var ) {

 if( svar_to_idx.empty() ) {
  return Inf< int >();
 }

 assert( std::is_sorted( svar_to_idx.begin(), svar_to_idx.end() ) );
 auto it = upper_bound( svar_to_idx.begin(), svar_to_idx.end(),
                        std::make_tuple( p_var, 0, 0 ),
                        [ & ]( auto & p1, auto & p2 ) {
                         return std::get<0>(p1) < std::get<0>(p2);
                        } );

 // Now it refers to the first (group of) element(s) greater than p_var
 if( it == svar_to_idx.begin() ) {
  return Inf< int >();
 }
 --it;

 // First element of the variable group
 auto first = const_cast<const ColVariable *>(std::get< 0 >( *it ));
 int distance = std::distance( first, p_var );

 if( distance < std::get< 2 >( *it ) ) {
  // The element belongs to this group
  return std::get< 1 >( *it ) + distance;
 } else {
  // The element doesn't exist
  return Inf< int >();
 }
}

/*--------------------------------------------------------------------------*/

int MILPSolver::index_of_dynamic_variable( const ColVariable * p_var ) {

 assert( std::is_sorted( dvar_to_idx.begin(), dvar_to_idx.end() ) );
 auto it = lower_bound( dvar_to_idx.begin(), dvar_to_idx.end(),
                        std::make_pair( p_var, 0 ),
                        [ & ]( auto & p1, auto & p2 ) {
                         return p1.first < p2.first;
                        } );

 if( it != dvar_to_idx.end() && it->first == p_var ) {
  return it->second;
 } else {
  return Inf< int >();
 }
}

/*--------------------------------------------------------------------------*/

int MILPSolver::index_of_constraint( const FRowConstraint * p_const ) {
 auto i = index_of_static_constraint( p_const );
 if( i < Inf< int >() ) {
  return i;
 } else {
  return index_of_dynamic_constraint( p_const );
 }
}

/*--------------------------------------------------------------------------*/

int MILPSolver::index_of_static_constraint( const FRowConstraint * p_const ) {

 if( scon_to_idx.empty() ) {
  return Inf< int >();
 }

 assert( std::is_sorted( scon_to_idx.begin(), scon_to_idx.end() ) );
 auto it = upper_bound( scon_to_idx.begin(), scon_to_idx.end(),
                        std::make_tuple( p_const, 0, 0 ),
                        [ & ]( auto & p1, auto & p2 ) {
                         return std::get<0>(p1) < std::get<0>(p2);
                        } );

 // Now it refers to the first (group of) element(s) greater than p_var
 if( it == scon_to_idx.begin() ) {
  return Inf< int >();
 }
 --it;

 // First element of the constraint group
 auto first = const_cast<const FRowConstraint *>(std::get< 0 >( *it ));
 int distance = std::distance( first, p_const );

 if( distance < std::get< 2 >( *it ) ) {
  // The element belongs to this group
  return std::get< 1 >( *it ) + distance;
 } else {
  // The element doesn't exist
  return Inf< int >();
 }
}

/*--------------------------------------------------------------------------*/

int MILPSolver::index_of_dynamic_constraint( const FRowConstraint * p_const ) {

 assert( std::is_sorted( dcon_to_idx.begin(), dcon_to_idx.end() ) );
 auto it = lower_bound( dcon_to_idx.begin(), dcon_to_idx.end(),
                        std::make_pair( p_const, 0 ),
                        [ & ]( auto & p1, auto & p2 ) {
                         return p1.first < p2.first;
                        } );

 if( it != dcon_to_idx.end() && it->first == p_const ) {
  return it->second;
 } else {
  return Inf< int >();
 }
}

/*--------------------------------------------------------------------------*/

ColVariable * MILPSolver::variable_with_index( int i ) {
 if( i < static_vars ) {
  return static_variable_with_index( i );
 } else {
  return dynamic_variable_with_index( i );
 }
}

/*--------------------------------------------------------------------------*/

ColVariable * MILPSolver::static_variable_with_index( int i ) {

 if( idx_to_svar.empty() ) {
  return nullptr;
 }

 assert( std::is_sorted( idx_to_svar.begin(), idx_to_svar.end() ) );
 auto it = upper_bound( idx_to_svar.begin(), idx_to_svar.end(),
                        std::make_pair( i, nullptr ),
                        [ & ]( auto & p1, auto & p2 ) {
                         return p1.first < p2.first;
                        } );

 // Now it refers to the first (group of) element(s) greater than i
 if( it == idx_to_svar.begin() ) {
  return nullptr;
 }
 --it;

 int distance = i - it->first;
 return it->second + distance;
}

/*--------------------------------------------------------------------------*/

ColVariable * MILPSolver::dynamic_variable_with_index( int i ) {

 assert( std::is_sorted( idx_to_dvar.begin(), idx_to_dvar.end() ) );
 auto it = lower_bound( idx_to_dvar.begin(), idx_to_dvar.end(),
                        std::make_pair( i, nullptr ),
                        [ & ]( auto & p1, auto & p2 ) {
                         return p1.first < p2.first;
                        } );

 if( it != idx_to_dvar.end() && it->first == i ) {
  return it->second;
 } else {
  return nullptr;
 }
}

/*--------------------------------------------------------------------------*/

FRowConstraint * MILPSolver::constraint_with_index( int i ) {
 if( i < static_cons ) {
  return static_constraint_with_index( i );
 } else {
  return dynamic_constraint_with_index( i );
 }
}

/*--------------------------------------------------------------------------*/

FRowConstraint * MILPSolver::static_constraint_with_index( int i ) {

 if( idx_to_scon.empty() ) {
  return nullptr;
 }

 assert( std::is_sorted( idx_to_scon.begin(), idx_to_scon.end() ) );
 auto it = upper_bound( idx_to_scon.begin(), idx_to_scon.end(),
                        std::make_pair( i, nullptr ),
                        [ & ]( auto & p1, auto & p2 ) {
                         return p1.first < p2.first;
                        } );

 // Now it refers to the first (group of) element(s) greater than i
 if( it == idx_to_scon.begin() ) {
  return nullptr;
 }
 --it;

 int distance = i - it->first;
 return it->second + distance;
}

/*--------------------------------------------------------------------------*/

FRowConstraint * MILPSolver::dynamic_constraint_with_index( int i ) {

 assert( std::is_sorted( idx_to_dcon.begin(), idx_to_dcon.end() ) );
 auto it = lower_bound( idx_to_dcon.begin(), idx_to_dcon.end(),
                        std::make_pair( i, nullptr ),
                        [ & ]( auto & p1, auto & p2 ) {
                         return p1.first < p2.first;
                        } );

 if( it != idx_to_dcon.end() && it->first == i ) {
  return it->second;
 } else {
  return nullptr;
 }
}

/*--------------------------------------------------------------------------*/
/*------------- AUXILIARY METHODS FOR POPULATING THE PROBLEM  --------------*/
/*--------------------------------------------------------------------------*/

void MILPSolver::count_constraints( FRowConstraint & constraint, int & n_rows ) {
 BOOST_LOG_TRIVIAL(trace) << "MILPSolver::count_constraints(): row " << n_rows << " " << constraint;
 const auto *fun = dynamic_cast<const LinearFunction *>(constraint.get_function());
 if( fun != nullptr ) {
  ++n_rows;
 } else {
  throw ( std::invalid_argument( "The Constraint is not linear" ) );
 }
}

/*--------------------------------------------------------------------------*/

void MILPSolver::count_variables( ColVariable & variable, int & n_cols ) {
 BOOST_LOG_TRIVIAL(trace) << "MILPSolver::count_variables(): col " << n_cols << " " << variable;
 ++n_cols;
}

/*--------------------------------------------------------------------------*/

void MILPSolver::count_nzelements( ColVariable & variable,
                                   int & nz_elements,
                                   int & var ) {
 BOOST_LOG_TRIVIAL( trace ) << "MILPSolver::count_nzelements(): nz/cnt "
                            << nz_elements << "/" << var << " " << variable;

 /*
  * Since counting non-zero elements requires checking if each active thing
  * is a FRowConstraint, we populate active_constraints and
  * active_bounds here so we don't have to loop over active stuff
  * once again later.
  */

 if( active_constraints.size() < numcols ) {
  active_constraints.resize( static_cast<unsigned long>(numcols) );
  active_bounds.resize( static_cast<unsigned long>(numcols) );
 }

 for( auto * i : variable.active_stuff() ) {
  auto * row = dynamic_cast<FRowConstraint *>(i);
  if( row != nullptr ) {
   active_constraints[ var ].push_back( row );
   ++nz_elements;
  }
  auto * box = dynamic_cast<OneVarConstraint *>(i);
  if( box != nullptr ) {
   active_bounds[ var ].push_back( box );
  }
  // auto *obj = dynamic_cast<Objective *>(i);
  // if( obj != nullptr ) {
  // }
 }
 ++var;
}

/*--------------------------------------------------------------------------*/

void MILPSolver::scan_variable( ColVariable & var, int & n, int & col ) {
 if( n < 0 ) {
  BOOST_LOG_TRIVIAL( trace ) << "MILPSolver::scan_variable(): D#"
                             << col << " " << var;
  dvar_to_idx.emplace_back( &var, col );
  idx_to_dvar.emplace_back( col, &var );
 } else {
  BOOST_LOG_TRIVIAL( trace ) << "MILPSolver::scan_variable(): S#"
                             << n << "/" << col << " " << var;

  if( n == 0 ) {
   // The tuple's third field will be filled later
   svar_to_idx.emplace_back( &var, col, 0 );
   idx_to_svar.emplace_back( col, &var );
  }
  ++n;
 }

 if( var.is_fixed() ) {
  lb[ col ] = var.get_value();
  ub[ col ] = var.get_value();
 } else {

  lb[ col ] = var.get_lb();
  ub[ col ] = var.get_ub();

  int num_bounds = static_cast<int>(active_bounds[ col ].size());
  for( int j = 0; j < num_bounds; ++j ) {
   auto * bound = active_bounds[ col ][ j ];
   if( lb[ col ] < bound->get_lhs() ) {
    lb[ col ] = bound->get_lhs();
   }
   if( ub[ col ] > bound->get_rhs() ) {
    ub[ col ] = bound->get_rhs();
   }
  }
 }

 if( var.is_integer() ) {
  ++int_vars;
  if( var.is_unitary() && var.is_positive() ) {
   xctype[ col ] = 'B'; // Binary
  } else {
   xctype[ col ] = 'I'; // Integer
  }
 } else {
  xctype[ col ] = 'C';  // Continuous
 }

 /*
  * From CPLEX documentation about matval, matbeg, matcnt, and matind:
  *
  * CPLEX needs to know only the nonzero coefficients.
  * These are grouped by column in the array matval.
  * The nonzero elements of every column must be stored in sequential locations
  * in this array with matbeg[j] containing the index of the beginning of
  * column j and matcnt[j] containing the number of entries in column j.
  * The components of matbeg must be in ascending order.
  * For each k, matind[k] specifies the row number of the
  * corresponding coefficient, matval[k].
  */

 int nz_elements = static_cast<int>(active_constraints[ col ].size());
 matcnt[ col ] = nz_elements;

 if( col == 0 ) {
  matbeg[ col ] = 0;
 } else {
  matbeg[ col ] = matbeg[ col - 1 ] + matcnt[ col - 1 ];
 }

 for( int j = 0; j < nz_elements; ++j ) {

  auto * p_const = dynamic_cast<FRowConstraint *> (active_constraints[ col ][ j ]);
  const auto * p_fun = dynamic_cast<const LinearFunction *> (p_const
   ->get_function());

  auto it = std::find_if( p_fun->get_v_var().begin(),
                          p_fun->get_v_var().end(),
                          [ & ]( LinearFunction::coeff_pair pair ) {
                           return pair.first == &var;
                          } );

  if( it != p_fun->get_v_var().end() ) {
   matval[ matbeg[ col ] + j ] = it->second;
   matind[ matbeg[ col ] + j ] = index_of_constraint( p_const );
  } else {
   // This should never happen because we are looping on the active contraints
   throw ( std::invalid_argument( "This ColVariable is not active in the examined FRowConstraint" ) );
  }
 }
 ++col;
}

/*--------------------------------------------------------------------------*/

void MILPSolver::scan_constraint( FRowConstraint & con, int & n, int & row ) {
 if( n < 0 ) {
  BOOST_LOG_TRIVIAL( trace ) << "MILPSolver::scan_constraint(): D#"
                             << row << " " << con;
  dcon_to_idx.emplace_back( &con, row );
  idx_to_dcon.emplace_back( row, &con );
 } else {
  BOOST_LOG_TRIVIAL( trace ) << "MILPSolver::scan_constraint(): S#"
                             << n << "/" << row << " " << con;

  if( n == 0 ) {
   // The tuple's third field will be filled later
   scon_to_idx.emplace_back( &con, row, 0 );
   idx_to_scon.emplace_back( row, &con );
  }
  ++n;
 }

 const auto * lf = dynamic_cast<const LinearFunction *>(con.get_function());
 if( lf == nullptr ) {
  throw ( std::invalid_argument( "The Constraint is not linear" ) );
 }

/*
 * We need to define the sense of the constraints as requested by CPLEX.
 *
 * In SMS++ FRowConstraints are defined as:
 * LHS <= ( some function from Variables to reals ) <= RHS
 *
 * CPLEX uses rhs and rngval arrays as documented in CPXcopylp reference.
 */

 auto const_lhs = con.get_lhs();
 auto const_rhs = con.get_rhs();

 if( const_lhs == const_rhs ) {
  // LHS <= function <= RHS, with LHS = RHS
  // becomes:
  // function = RHS
  sense[ row ] = 'E';
  rhs[ row ] = const_rhs;

 } else if( const_lhs == -Inf< double >() ) {
  // -inf <= function <= RHS
  // becomes:
  // function <= RHS
  sense[ row ] = 'L';
  rhs[ row ] = const_rhs;

 } else if( const_rhs == Inf< double >() ) {
  // LHS <= function <= inf
  // becomes:
  // function >= LHS
  sense[ row ] = 'G';
  rhs[ row ] = const_lhs;

 } else {
  // LHS <= function <= RHS
  // becomes:
  // LHS <= function <= LHS + (range),
  // with range = RHS - LHS
  sense[ row ] = 'R';
  rhs[ row ] = const_lhs;
  rngval[ row ] = const_rhs - const_lhs;
 }
 ++row;
}

/*--------------------------------------------------------------------------*/

void MILPSolver::scan_objective( const FRealObjective * obj ) {
 BOOST_LOG_TRIVIAL( trace ) << "MILPSolver::scan_objective() " << *obj;

 const auto *lin_fun = dynamic_cast<const LinearFunction *> (obj->get_function());
 int k = 0;

 if( lin_fun != nullptr ) {
  for( auto el : lin_fun->get_v_var() ) {
   k = index_of_variable( el.first );
   objective[ k ] = el.second;
  }
 } else {
  const auto *dquad_fun = dynamic_cast<const DQuadFunction *> (obj->get_function());
  if( dquad_fun != nullptr ) {
   for( auto el : dquad_fun->get_v_var() ) {
    // DQuadFunction::get_v_var() returns std::tuples of 3 elements
    k = index_of_variable( std::get< 0 >( el ) );
    objective[ k ] = std::get< 1 >( el );
    q_objective[ k ] = std::get< 2 >( el );
   }
  } else {
   throw ( std::invalid_argument( "Unknown type of Objective Function" ) );
  }
 }
}

/*--------------------------------------------------------------------------*/

void MILPSolver::process_modifications() {
 /*
  * This function processes one modification after another, without
  * any attempt of optimization, moreover you have to be CAREFUL to write
  * all the cases in order from the most specialized to the more generic,
  * e.g., OneVarConstraintMod before RowConstraintMod before ConstraintMod,
  * otherwise the generic case will intercept the more specialized Mods.
  */
 for( auto mod = front() ; mod ; mod = front() ) {
 // while( !v_mod.empty() ) {
  // auto mod = v_mod.front();

  // A function like this is needed to be called recursively with GroupModifications
  std::function< void( sp_Mod ) > f;
  f = [ this, &f ]( const sp_Mod& mod ) {

   {
    const auto tmod = std::dynamic_pointer_cast< GroupModification >( mod );
    if( tmod ) {
     BOOST_LOG_TRIVIAL( trace ) << "GroupModification containing:";
     for( const auto & submod : tmod->sub_Modifications() ) {
      f( submod );
     }
     return;
    }
   }
   {
    const auto tmod = std::dynamic_pointer_cast< VariableMod >( mod );
    if( tmod ) {
     BOOST_LOG_TRIVIAL( trace ) << *mod;
     var_modification( tmod.get() );
     return;
    }
   }
   {
    const auto tmod = std::dynamic_pointer_cast< ObjectiveMod >( mod );
    if( tmod ) {
     BOOST_LOG_TRIVIAL( trace ) << *mod;
     of_modification( tmod.get() );
     return;
    }
   }
   {
    const auto tmod = std::dynamic_pointer_cast< OneVarConstraintMod >( mod );
    if( tmod ) {
     BOOST_LOG_TRIVIAL( trace ) << *mod;
     bound_modification( tmod.get() );
     return;
    }
   }
   {
    const auto tmod = std::dynamic_pointer_cast< RowConstraintMod >( mod );
    if( tmod ) {
     BOOST_LOG_TRIVIAL( trace ) << *mod;
     const_modification( tmod.get() );
     return;
    }
   }
   {
    const auto tmod = std::dynamic_pointer_cast< ConstraintMod >( mod );
    if( tmod ) {
     BOOST_LOG_TRIVIAL( trace ) << *mod;
     const_modification( tmod.get() );
     return;
    }
   }
   {
    const auto tmod = std::dynamic_pointer_cast< FunctionMod >( mod );
    if( tmod ) {
     BOOST_LOG_TRIVIAL( trace ) << *mod;
     function_modification( tmod.get() );
     return;
    }
   }
   {
    const auto tmod = std::dynamic_pointer_cast< FunctionModVars >( mod );
    if( tmod ) {
     BOOST_LOG_TRIVIAL( trace ) << *mod;
     function_vars_modification( tmod.get() );
     return;
    }
   }
   {
    const auto tmod = std::dynamic_pointer_cast< BlockModAD >( mod );
    if( tmod ) {
     BOOST_LOG_TRIVIAL( trace ) << *mod;
     dynamic_modification( tmod.get() );
     return;
    }
   }
   {
    const auto tmod = std::dynamic_pointer_cast< NBModification >( mod );
    if( tmod ) {
     BOOST_LOG_TRIVIAL( trace ) << "\033[1;33m" << *mod << "\033[0m";
     clear_problem();
     load_problem();
    }
   }
  };

  f( mod );
  // v_mod.pop_front();
  pop_front();
 }
}

/*--------------------------------------------------------------------------*/

void MILPSolver::var_modification( VariableMod * mod ){
 const auto var = dynamic_cast< const ColVariable * >( mod->variable() );
 // int idx = index_of_variable( var );

 // Update the number of integer variables
 if( var->is_integer( mod->old_state() ) ) {
  --int_vars;
 }

 if( var->is_integer() ) {
  ++int_vars;
 }

 // TODO: Update MILPSolver bounds
 // TODO: Update MILPSolver ctypes?
}

/*--------------------------------------------------------------------------*/

void MILPSolver::of_modification( ObjectiveMod * mod ){
 // TODO: Update MILPSolver objsense
}

/*--------------------------------------------------------------------------*/

void MILPSolver::const_modification( ConstraintMod * mod ){
 auto * p_const = dynamic_cast<FRowConstraint *>(mod->constraint());

 // TODO: Update MILPSolver rhs
 // TODO: Update MILPSolver sense
 // TODO: Update MILPSolver rngval
}

/*--------------------------------------------------------------------------*/

void MILPSolver::bound_modification( OneVarConstraintMod * mod ){
 auto * p_const = dynamic_cast<OneVarConstraint *>(mod->constraint());
 auto * p_var = dynamic_cast<ColVariable *>(p_const->get_active_var( 0 ));

 // TODO: Update MILPSolver bounds
}

/*--------------------------------------------------------------------------*/

void MILPSolver::function_modification( FunctionMod * mod ){
 // TODO: update objective [and q_objective], or constraint coefficient
}

/*--------------------------------------------------------------------------*/

void MILPSolver::function_vars_modification( FunctionModVars * mod ){
 // TODO: update objective [and q_objective], or constraint coefficient
}

/*--------------------------------------------------------------------------*/

void MILPSolver::dynamic_modification( BlockModAD * mod ) {
 auto * addcon_mod = dynamic_cast<BlockModAdd< FRowConstraint > *>(mod);
 if( addcon_mod ) {
  for( auto * i : addcon_mod->added() ) {
   add_dynamic_constraint( i );
  }
  return;
 }

 auto * rmvcon_mod = dynamic_cast<BlockModRmv< FRowConstraint > *>(mod);
 if( rmvcon_mod ) {
  for( const auto & i : rmvcon_mod->removed() ) {
   remove_dynamic_constraint( &i );
  }
  return;
 }

 auto * addvar_mod = dynamic_cast<BlockModAdd< ColVariable > *>(mod);
 if( addvar_mod ) {
  for( auto * i : addvar_mod->added() ) {
   add_dynamic_variable( i );
  }
  return;
 }

 auto * rmvvar_mod = dynamic_cast<BlockModRmv< ColVariable > *>(mod);
 if( rmvvar_mod ) {
  for( const auto & i : rmvvar_mod->removed() ) {
   remove_dynamic_variable( &i );
  }
  return;
 }

 auto * addbnd_mod = dynamic_cast<BlockModAdd< LB0Constraint > *>(mod);
 if( addbnd_mod ) {
  for( auto * i : addbnd_mod->added() ) {
   add_dynamic_bound( i );
  }
  return;
 }

 auto * rmvbnd_mod = dynamic_cast<BlockModRmv< LB0Constraint > *>(mod);
 if( rmvbnd_mod ) {
  for( const auto & i : rmvbnd_mod->removed() ) {
   remove_dynamic_bound( &i );
  }
  return;
 }

 throw std::invalid_argument( "Unknown type of BlockAD" );
}

/*--------------------------------------------------------------------------*/

void MILPSolver::add_dynamic_constraint( FRowConstraint * p_const ) {

 const auto * p_fun =
  dynamic_cast<const LinearFunction *>(p_const->get_function());
 if( p_fun == nullptr ) {
  throw std::invalid_argument( "The Constraint is not linear" );
 }

 // Set the constraint as active for all its variables
 int nzcnt = p_const->get_num_active_var();
 for( int i = 0; i < nzcnt; ++i ) {
  auto * p_var = dynamic_cast<ColVariable *>(p_fun->get_active_var( i ));
  active_constraints[ index_of_variable( p_var ) ].push_back( p_const );
 }

 // Update the dictionaries
 auto it = lower_bound( dcon_to_idx.begin(),
                        dcon_to_idx.end(),
                        p_const,
                        [ & ]( const_int pair, FRowConstraint * c ) {
                         return pair.first < c;
                        } );
 dcon_to_idx.insert( it, { p_const, numrows } );
 idx_to_dcon.emplace_back( numrows, p_const );
 ++numrows;
}

/*--------------------------------------------------------------------------*/

void MILPSolver::add_dynamic_variable( ColVariable * p_var ) {
 // Get the constraints and bounds of the new variable
 std::vector< FRowConstraint * > var_constraints;
 std::vector< OneVarConstraint * > var_bounds;


 for( auto * stuff : p_var->active_stuff() ) {
  auto * constraint = dynamic_cast<FRowConstraint *>(stuff);
  if( constraint != nullptr ) {
   var_constraints.push_back( constraint );
  }
  auto * bound = dynamic_cast<OneVarConstraint *>(stuff);
  if( bound != nullptr ) {
   var_bounds.push_back( bound );
  }
 }

 active_constraints.emplace_back( var_constraints );
 active_bounds.emplace_back( var_bounds );

 // Update the number of integer vars
 if ( p_var->is_integer() ) {
  ++int_vars;
 }

 // Update the dictionaries
 auto it = lower_bound( dvar_to_idx.begin(),
                        dvar_to_idx.end(),
                        p_var,
                        [ & ]( var_int pair, ColVariable * v ) {
                         return pair.first < v;
                        } );
 dvar_to_idx.insert( it, { p_var, numcols } );
 idx_to_dvar.emplace_back( numcols, p_var );
 ++numcols;
}

/*--------------------------------------------------------------------------*/

void MILPSolver::add_dynamic_bound( OneVarConstraint * p_bound ) {
 auto * p_var = dynamic_cast<ColVariable *>(p_bound->get_active_var( 0 ));
 auto active_bnds = active_bounds[ index_of_variable( p_var ) ];

 // Look if the bound is already there (say, added with the Variable)
 auto it = std::find( active_bnds.begin(), active_bnds.end(), p_bound );
 if( it < active_bnds.end() ) {
  return;
 }

 // Add the bound
 active_bnds.emplace_back( p_bound );
}

/*--------------------------------------------------------------------------*/

void MILPSolver::remove_dynamic_constraint( const FRowConstraint * p_const ){

 // Remove the constraint from the dictionaries
 int index = 0;
 assert( std::is_sorted( dcon_to_idx.begin(), dcon_to_idx.end() ) );
 assert( std::is_sorted( idx_to_dcon.begin(), idx_to_dcon.end() ) );

 auto it1 = lower_bound( dcon_to_idx.begin(), dcon_to_idx.end(),
                         std::make_pair( p_const, 0 ),
                         [ & ]( auto & p1, auto & p2 ) {
                          return p1.first < p2.first;
                         } );

 if( it1 != dcon_to_idx.end() && it1->first == p_const ) {
  index = it1->second;
  dcon_to_idx.erase( it1 );
 } else {
  return; // TODO Not sure if it's ok
 }

 auto it2 = lower_bound( idx_to_dcon.begin(), idx_to_dcon.end(),
                         std::make_pair( index, 0 ),
                         [ & ]( auto & p1, auto & p2 ) {
                          return p1.first < p2.first;
                         } );

 if( it2 != idx_to_dcon.end() && it2->second == p_const ) {
  idx_to_dcon.erase( it2 );
 }

 // Update the other indices
 for( auto & it: dcon_to_idx ) {
  if( it.second > index ) {
   it.second--;
  }
 }
 for( auto & it: idx_to_dcon ) {
  if( it.first > index ) {
   it.first--;
  }
 }

 // Remove the constraint from the active contraints
 for( auto & constraints: active_constraints ) {
  auto constraint = find( constraints.begin(), constraints.end(), p_const );
  if( constraint != constraints.end() ) {
   constraints.erase( constraint );
  }
 }

 --numrows;
}

/*--------------------------------------------------------------------------*/

void MILPSolver::remove_dynamic_variable( const ColVariable * p_var ){

 // Remove the constraint from the dictionaries
 int index = 0;
 assert( std::is_sorted( dvar_to_idx.begin(), dvar_to_idx.end() ) );
 assert( std::is_sorted( idx_to_dvar.begin(), idx_to_dvar.end() ) );

 auto it1 = lower_bound( dvar_to_idx.begin(), dvar_to_idx.end(),
                         std::make_pair( p_var, 0 ),
                         [ & ]( auto & p1, auto & p2 ) {
                          return p1.first < p2.first;
                         } );

 if( it1 != dvar_to_idx.end() && it1->first == p_var ) {
  index = it1->second;
  dvar_to_idx.erase( it1 );
 } else {
  return; // TODO Not sure if it's ok
 }

 auto it2 = lower_bound( idx_to_dvar.begin(), idx_to_dvar.end(),
                         std::make_pair( index, 0 ),
                         [ & ]( auto & p1, auto & p2 ) {
                          return p1.first < p2.first;
                         } );

 if( it2 != idx_to_dvar.end() && it2->second == p_var ) {
  idx_to_dvar.erase( it2 );
 }

 // Update the other indices
 for( auto & it: dvar_to_idx ) {
  if( it.second > index ) {
   it.second--;
  }
 }
 for( auto & it: idx_to_dvar ) {
  if( it.first > index ) {
   it.first--;
  }
 }

 // Remove the variable's active contraints
 active_constraints.erase( active_constraints.begin() + index );
 active_bounds.erase( active_bounds.begin() + index );

 --numcols;

 // Update the number of integer vars
 if( p_var->is_integer() ) {
  --int_vars;
 }
}

/*--------------------------------------------------------------------------*/

void MILPSolver::remove_dynamic_bound( const OneVarConstraint * p_bound ) {
 // Remove the bound from the active bounds
 for( auto & bounds: active_bounds ) {
  auto bound = find( bounds.begin(), bounds.end(), p_bound );
  if( bound != bounds.end() ) {
   bounds.erase( bound );
  }
 }
}

/*--------------------------------------------------------------------------*/

template< typename T >
std::string MILPSolver::log_vector( std::vector< T > v ) {
 std::string temp_log = "[";
 for( auto i : v ) {
  temp_log += " " + std::to_string(i);
 }
 temp_log += "]";
 return temp_log;
}

/*--------------------------------------------------------------------------*/
/*----------------------- End File MILPSolver.cpp --------------------------*/
/*--------------------------------------------------------------------------*/ 
