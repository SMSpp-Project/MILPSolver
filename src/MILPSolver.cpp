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
 * Copyright &copy; by Antonio Frangioni, Niccolò Iardella
 */
/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include <functional>
#include <queue>

#include <Block.h>
#include <OneVarConstraint.h>
#include <LinearFunction.h>
#include <DQuadFunction.h>

#include "MILPSolver.h"

/*--------------------------------------------------------------------------*/
/*------------------------- NAMESPACE AND USING ----------------------------*/
/*--------------------------------------------------------------------------*/

using namespace SMSpp_di_unipi_it;

SMSpp_insert_in_factory_cpp_0( MILPSolver )

/*--------------------------------------------------------------------------*/
/*--------------------- CONSTRUCTOR AND DESTRUCTOR -------------------------*/
/*--------------------------------------------------------------------------*/

MILPSolver::MILPSolver() : CDASolver() {}

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

const std::vector< char * > & MILPSolver::get_rowname() const {
 return rowname;
}

const std::vector< char * > & MILPSolver::get_colname() const {
 return colname;
}

int MILPSolver::get_num_integer_vars() const {
 return int_vars;
}

/*--------------------------------------------------------------------------*/
/*-------------------------------- SET_BLOCK -------------------------------*/
/*--------------------------------------------------------------------------*/

void MILPSolver::set_Block( Block * block ) {

 Solver::set_Block( block );

 if( block ) {

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

  load_problem();
 }
}

/*--------------------------------------------------------------------------*/

void MILPSolver::clear_problem( unsigned int what ) {

 if( what & 1u ) {
  matbeg.clear();
  matcnt.clear();
  matind.clear();
  matval.clear();
  xctype.clear();

  for( auto & i: colname )
   delete i;
  for( auto & i: rowname )
   delete i;
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

void MILPSolver::load_problem() {

 numrows = 0;
 numcols = 0;
 nzelements = 0;

 std::queue< Block * > Q;

 // Locking the Block
 bool owned = f_Block->is_owned_by( f_id );
 if( !owned && !f_Block->read_lock() ) {
  throw std::runtime_error( "Unable to lock the Block" );
 }

 int num_block = 0; // Counter for the blocks
 int row = 0;       // Counter for the rows
 int col = 0;       // Counter for the columns

 // Count variables and constraints
 // --------------------------------------------------------------------------

 Q.push( f_Block );
 while( !Q.empty() ) {
  Block * q_Block = Q.front();
  Q.pop();
  DEBUG_LOG( "Processing Block " << num_block << " ["
                                 << q_Block << "]" << std::endl );
  DEBUG_LOG( *q_Block << std::endl );

  for( auto * i : q_Block->get_nested_Blocks() ) {
   Q.push( i );
  }

  DEBUG_LOG( "MILPSolver::set_Block() counting static constraints"
              << std::endl );
  for( const auto & i : q_Block->get_static_constraints() ) {
   // Singles
   if( un_any_thing_0( FRowConstraint, i,
                       {
                        ++numrows;
                        ++static_cons;
                       }
   ) ) {
    continue;
   }
   // Vectors
   if( un_any_thing_1( FRowConstraint, i,
                       {
                        numrows += var.size();
                        static_cons += var.size();
                       }
   ) ) {
    continue;
   }
   // Multiarrays
   if( un_any_thing_K( FRowConstraint, i,
                       {
                        numrows += var.num_elements();
                        static_cons += var.num_elements();
                       }
   ) ) {
    continue;
   }
  }

  DEBUG_LOG( "MILPSolver::set_Block() counting dynamic constraints"
              << std::endl );
  for( const auto & i : q_Block->get_dynamic_constraints() ) {
   // Single lists
   if( un_any_thing_0( std::list< FRowConstraint >, i,
                       {
                        numrows += var.size();
                       }
   ) ) {
    continue;
   }
   // Vectors of lists
   if( un_any_thing_1( std::list< FRowConstraint >, i,
                       {
                        for( auto & el: var ) {
                         numrows += el.size();
                        }
                       }
   ) ) {
    continue;
   }
   // Multiarrays of lists
   if( un_any_thing_K( std::list< FRowConstraint >, i,
                       {
                        auto it = var.data();
                        for( auto i = var.num_elements(); i--; ++it ) {
                         numrows += it->size();
                        }
                       }
   ) ) {
    continue;
   }
  }

  DEBUG_LOG( "MILPSolver::set_Block() counting static variables" << std::endl );
  for( const auto & i : q_Block->get_static_variables() ) {
   // Singles
   if( un_any_thing_0( ColVariable, i,
                       {
                        ++numcols;
                        ++static_vars;
                       }
   ) ) {
    continue;
   }
   // Vectors
   if( un_any_thing_1( ColVariable, i,
                       {
                        numcols += var.size();
                        static_vars += var.size();
                       }
   ) ) {
    continue;
   }
   // Multiarrays
   if( un_any_thing_K( ColVariable, i,
                       {
                        numcols += var.num_elements();
                        static_vars += var.num_elements();
                       }
   ) ) {
    continue;
   }
  }

  DEBUG_LOG( "MILPSolver::set_Block() counting dynamic variables"
              << std::endl );
  for( const auto & i : q_Block->get_dynamic_variables() ) {
   // Single lists
   if( un_any_thing_0( std::list< ColVariable >, i,
                       {
                        numcols += var.size();
                       }
   ) ) {
    continue;
   }
   // Vectors of lists
   if( un_any_thing_1( std::list< ColVariable >, i,
                       {
                        for( auto & el: var ) {
                         numcols += el.size();
                        }
                       }
   ) ) {
    continue;
   }
   // Multiarrays of lists
   if( un_any_thing_K( std::list< ColVariable >, i,
                       {
                        auto it = var.data();
                        for( auto i = var.num_elements(); i--; ++it ) {
                         numcols += it->size();
                        }
                       }
   ) ) {
    continue;
   }
  }

  DEBUG_LOG( "MILPSolver::set_Block() nonzero elements" << std::endl );
  auto counter = [ this ]( ColVariable & var ) {
   for( auto * i : var.active_stuff() ) {
    auto * row = dynamic_cast<FRowConstraint *>(i);
    if( row != nullptr ) {
     ++nzelements;
    }
   }
  };

  for( const auto & i : q_Block->get_static_variables() ) {
   un_any_const_static( i, counter, un_any_type< ColVariable >() );
  }

  for( const auto & i : q_Block->get_dynamic_variables() ) {
   un_any_const_dynamic( i, counter, un_any_type< ColVariable >() );
  }
  ++num_block;
 }

 DEBUG_LOG( "Number of blocks      = " << num_block << std::endl );
 DEBUG_LOG( "numrows (constraints) = " << numrows
                                       << " (S:" << static_cons
                                       << "/D:" << numrows - static_cons << ")"
                                       << std::endl );
 DEBUG_LOG( "numcols (variables)   = " << numcols
                                       << " (S:" << static_vars
                                       << "/D:" << numcols - static_vars << ")"
                                       << std::endl );
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
 colname.resize( numcols, nullptr );
 rowname.resize( numrows, nullptr );

 svar_to_idx.clear();
 idx_to_svar.clear();
 scon_to_idx.clear();
 idx_to_scon.clear();
 dvar_to_idx.clear();
 idx_to_dvar.clear();
 dcon_to_idx.clear();
 idx_to_dcon.clear();

 svar_to_idx.reserve( numcols );
 idx_to_svar.reserve( numcols );
 dvar_to_idx.reserve( numcols );
 idx_to_dvar.reserve( numcols );
 scon_to_idx.reserve( numrows );
 idx_to_scon.reserve( numrows );
 dcon_to_idx.reserve( numrows );
 idx_to_dcon.reserve( numrows );

 // Scan the static constraints
 // --------------------------------------------------------------------------

 num_block = 0;
 Q.push( f_Block );
 while( !Q.empty() ) {
  Block * q_Block = Q.front();
  Q.pop();

  for( auto i : q_Block->get_nested_Blocks() ) {
   Q.push( i );
  }

  int set = 0; // Counter for the constraint groups

  for( const auto & i : q_Block->get_static_constraints() ) {
   int elements = 0; // Counter for group elements
   int start = row;

   auto scan = [ this, &elements, &row ]( FRowConstraint & c ) {
    scan_constraint( c, elements, row );
   };
   un_any_const_static( i, scan, un_any_type< FRowConstraint >() );

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
  num_block++;
 }

 std::sort( scon_to_idx.begin(), scon_to_idx.end() );

 // Scan the dynamic constraints
 // --------------------------------------------------------------------------

 num_block = 0;
 Q.push( f_Block );
 while( !Q.empty() ) {
  Block * q_Block = Q.front();
  Q.pop();

  for( auto i : q_Block->get_nested_Blocks() ) {
   Q.push( i );
  }

  int set = 0; // Counter for the constraint groups

  for( const auto & i : q_Block->get_dynamic_constraints() ) {
   int elements = -1;
   int start = row;

   auto scan = [ this, &elements, &row ]( FRowConstraint & c ) {
    scan_constraint( c, elements, row );
   };
   un_any_const_dynamic( i, scan, un_any_type< FRowConstraint >() );

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
 }

 std::sort( dcon_to_idx.begin(), dcon_to_idx.end() );

 // Scan the static variables
 // --------------------------------------------------------------------------

 num_block = 0;
 Q.push( f_Block );

 while( !Q.empty() ) {
  Block * q_Block = Q.front();
  Q.pop();

  for( auto * i : q_Block->get_nested_Blocks() ) {
   Q.push( i );
  }

  int set = 0;   // Counter for the variable groups

  for( const auto & i : q_Block->get_static_variables() ) {
   int elements = 0; // Counter for group elements
   int start = col;

   auto scan = [ this, &elements, &col ]( ColVariable & v ) {
    scan_variable( v, elements, col );
   };
   un_any_const_static( i, scan, un_any_type< ColVariable >() );

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
  num_block++;
 }

 std::sort( svar_to_idx.begin(), svar_to_idx.end() );

// Scan the dynamic variables
 // --------------------------------------------------------------------------

 num_block = 0;
 Q.push( f_Block );
 while( !Q.empty() ) {
  Block * q_Block = Q.front();
  Q.pop();

  for( auto * i : q_Block->get_nested_Blocks() ) {
   Q.push( i );
  }

  int set = 0;   // Counter for the variable groups

  for( const auto & i : q_Block->get_dynamic_variables() ) {
   int elements = -1;
   int start = col;

   auto scan = [ this, &elements, &col ]( ColVariable & v ) {
    scan_variable( v, elements, col );
   };
   un_any_const_dynamic( i, scan, un_any_type< ColVariable >() );

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
 }

 std::sort( dvar_to_idx.begin(), dvar_to_idx.end() );

 // Scan the objective
 // --------------------------------------------------------------------------

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

 Q.push( f_Block );
 while( !Q.empty() ) {
  Block * q_Block = Q.front();
  Q.pop();

  for( auto * i : q_Block->get_nested_Blocks() ) {
   Q.push( i );
  }

  auto * p_obj = dynamic_cast< FRealObjective * >( q_Block->get_objective() );
  if( p_obj ) {
   scan_objective( p_obj );
  }
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
 if( !owned ) {
  f_Block->read_unlock();
 }
}

/*--------------------------------------------------------------------------*/

double MILPSolver::get_problem_lb( const ColVariable & var ) {
 double b = var.get_lb();

 for( auto * i : var.active_stuff() ) {
  auto * box = dynamic_cast<OneVarConstraint *>(i);
  if( box != nullptr ) {
   b = b < box->get_lhs() ? box->get_lhs() : b;
  }
 }

 return b;
}

/*--------------------------------------------------------------------------*/

double MILPSolver::get_problem_ub( const ColVariable & var ) {
 double b = var.get_ub();

 for( auto * i : var.active_stuff() ) {
  auto * box = dynamic_cast<OneVarConstraint *>(i);
  if( box != nullptr ) {
   b = b > box->get_rhs() ? box->get_rhs() : b;
  }
 }

 return b;
}

/*--------------------------------------------------------------------------*/

std::vector< FRowConstraint * >
MILPSolver::get_active_constraints( const ColVariable & var ) {
 std::vector< FRowConstraint * > active_constraints;
 for( auto * i : var.active_stuff() ) {
  auto * row = dynamic_cast<FRowConstraint *>(i);
  if( row != nullptr ) {
   active_constraints.push_back( row );
  }
 }
 return active_constraints;
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
                         return std::get< 0 >( p1 ) < std::get< 0 >( p2 );
                        } );

 // Now it refers to the first (group of) element(s) greater than p_var
 if( it == svar_to_idx.begin() ) {
  return Inf< int >();
 }
 --it;

 // First element of the variable group
 auto first = const_cast<const ColVariable *>(std::get< 0 >( *it ));
 int distance = std::distance( first, p_var );

 if( distance >= 0 && distance < std::get< 2 >( *it ) ) {
  // The element belongs to this group
  return std::get< 1 >( *it ) + distance;
 }
 // The element doesn't exist
 return Inf< int >();
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
                         return std::get< 0 >( p1 ) < std::get< 0 >( p2 );
                        } );

 // Now it refers to the first (group of) element(s) greater than p_var
 if( it == scon_to_idx.begin() ) {
  return Inf< int >();
 }
 --it;

 // First element of the constraint group
 auto first = const_cast<const FRowConstraint *>(std::get< 0 >( *it ));
 int distance = std::distance( first, p_const );

 if( distance >= 0 && distance < std::get< 2 >( *it ) ) {
  // The element belongs to this group
  return std::get< 1 >( *it ) + distance;
 }

 // The element doesn't exist
 return Inf< int >();
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

void MILPSolver::scan_variable( ColVariable & var, int & n, int & col ) {
 if( n < 0 ) {
  DEBUG_LOG( "MILPSolver::scan_variable(): D#" << col << " " << var );
  dvar_to_idx.emplace_back( &var, col );
  idx_to_dvar.emplace_back( col, &var );
 } else {
  DEBUG_LOG( "MILPSolver::scan_variable(): S#"
              << n << "/" << col << " " << var );

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
  lb[ col ] = get_problem_lb( var );
  ub[ col ] = get_problem_ub( var );
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

 auto active_constraints = get_active_constraints( var );
 int nz_elements = active_constraints.size();
 matcnt[ col ] = nz_elements;

 if( col == 0 ) {
  matbeg[ col ] = 0;
 } else {
  matbeg[ col ] = matbeg[ col - 1 ] + matcnt[ col - 1 ];
 }

 for( int j = 0; j < nz_elements; ++j ) {
  auto * p_const = active_constraints[ j ];
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
  DEBUG_LOG( "MILPSolver::scan_constraint(): D#" << row << " " << con );
  dcon_to_idx.emplace_back( &con, row );
  idx_to_dcon.emplace_back( row, &con );
 } else {
  DEBUG_LOG( "MILPSolver::scan_constraint(): S#"
              << n << "/" << row << " " << con );

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
 DEBUG_LOG( "MILPSolver::scan_objective() " << *obj );

 const auto * lf = dynamic_cast<const LinearFunction *> (obj->get_function());
 int k = 0;

 if( lf != nullptr ) {
  for( auto el : lf->get_v_var() ) {
   k = index_of_variable( el.first );
   objective[ k ] = el.second;
  }
 } else {
  const auto * qf = dynamic_cast<const DQuadFunction *> (obj->get_function());
  if( qf != nullptr ) {
   for( auto el : qf->get_v_var() ) {
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
 for( auto mod = front(); mod; mod = front() ) {

  // A function like this is needed to be called
  // recursively with GroupModifications
  std::function< void( sp_Mod ) > f;

  f = [ this, &f ]( const sp_Mod & mod ) {
   DEBUG_LOG( *mod );

   const auto gm = std::dynamic_pointer_cast< GroupModification >( mod );
   if( gm ) {
    DEBUG_LOG( "GroupModification containing:" << std::endl );
    for( const auto & submod : gm->sub_Modifications() ) {
     f( submod );
    }
    return;
   }

   const auto vm = std::dynamic_pointer_cast< VariableMod >( mod );
   if( vm ) {
    var_modification( vm.get() );
    return;
   }

   const auto om = std::dynamic_pointer_cast< ObjectiveMod >( mod );
   if( om ) {
    of_modification( om.get() );
    return;
   }


   const auto bm = std::dynamic_pointer_cast< OneVarConstraintMod >( mod );
   if( bm ) {
    bound_modification( bm.get() );
    return;
   }


   const auto tmod = std::dynamic_pointer_cast< RowConstraintMod >( mod );
   if( tmod ) {
    const_modification( tmod.get() );
    return;
   }

   const auto cm = std::dynamic_pointer_cast< ConstraintMod >( mod );
   if( cm ) {
    const_modification( cm.get() );
    return;
   }

   const auto fm = std::dynamic_pointer_cast< FunctionMod >( mod );
   if( fm ) {
    function_modification( fm.get() );
    return;
   }

   const auto fvm = std::dynamic_pointer_cast< FunctionModVars >( mod );
   if( fvm ) {
    function_vars_modification( fvm.get() );
    return;
   }

   const auto dm = std::dynamic_pointer_cast< BlockModAD >( mod );
   if( dm ) {
    dynamic_modification( dm.get() );
    return;
   }

   const auto nm = std::dynamic_pointer_cast< NBModification >( mod );
   if( nm ) {
    load_problem();
   }
  };

  f( mod );
  pop_front();
 }

#ifdef MILPSOLVER_DEBUG
 check_status();
#endif
}

/*--------------------------------------------------------------------------*/

void MILPSolver::var_modification( VariableMod * mod ) {
 const auto var = dynamic_cast< const ColVariable * >( mod->variable() );
 int idx = index_of_variable( var );

 if( !lb.empty() ) {
  lb[ idx ] = get_problem_lb( *var );
 }
 if( !ub.empty() ) {
  ub[ idx ] = get_problem_ub( *var );
 }

 if( !xctype.empty() ) {
  if( var->is_integer() ) {
   if( var->is_unitary() && var->is_positive() ) {
    xctype[ idx ] = 'B';
   } else {
    xctype[ idx ] = 'I';
   }
  } else {
   xctype[ idx ] = 'C';
  }
 }

 // Update the number of integer variables
 if( var->is_integer( mod->old_state() ) ) {
  --int_vars;
 }

 if( var->is_integer() ) {
  ++int_vars;
 }
}

/*--------------------------------------------------------------------------*/

void MILPSolver::of_modification( ObjectiveMod * mod ) {
 switch( mod->type() ) {

  case ObjectiveMod::eSetMin:
   objsense = 1;
   break;
  case ObjectiveMod::eSetMax:
   objsense = -1;
   break;
  default:
   throw std::invalid_argument( "Invalid type of ObjectiveMod" );
 }
}

/*--------------------------------------------------------------------------*/

void MILPSolver::const_modification( ConstraintMod * mod ) {
 auto * p_const = dynamic_cast<FRowConstraint *>(mod->constraint());
 int idx = index_of_constraint( p_const );

 RowConstraint::RHSValue const_lhs = NAN;
 RowConstraint::RHSValue const_rhs = NAN;

 switch( mod->type() ) {
  case ConstraintMod::eRelaxConst:
   if( !sense.empty() ) {
    sense[ idx ] = 'G';
   }
   if( !rhs.empty() ) {
    rhs[ idx ] = -Inf< double >();
   }
   if( !rngval.empty() ) {
    rngval[ idx ] = 0;
   }
   break;
  case ConstraintMod::eEnforceConst:
  case RowConstraintMod::eChgLHS:
  case RowConstraintMod::eChgRHS:
  case RowConstraintMod::eChgBTS:

   const_lhs = p_const->get_lhs();
   const_rhs = p_const->get_rhs();

   if( const_lhs == const_rhs ) {
    if( !sense.empty() ) {
     sense[ idx ] = 'E';
    }
    if( !rhs.empty() ) {
     rhs[ idx ] = const_rhs;
    }
    if( !rngval.empty() ) {
     rngval[ idx ] = 0;
    }

   } else if( const_lhs == -Inf< double >() ) {
    if( !sense.empty() ) {
     sense[ idx ] = 'L';
    }
    if( !rhs.empty() ) {
     rhs[ idx ] = const_rhs;
    }
    if( !rngval.empty() ) {
     rngval[ idx ] = 0;
    }

   } else if( const_rhs == Inf< double >() ) {
    if( !sense.empty() ) {
     sense[ idx ] = 'G';
    }
    if( !rhs.empty() ) {
     rhs[ idx ] = const_lhs;
    }
    if( !rngval.empty() ) {
     rngval[ idx ] = 0;
    }

   } else {
    if( !sense.empty() ) {
     sense[ idx ] = 'R';
    }
    if( !rhs.empty() ) {
     rhs[ idx ] = -Inf< double >();
    }
    if( !rngval.empty() ) {
     rngval[ idx ] = const_rhs - const_lhs;
    }
   }

   break;
  default:
   throw std::invalid_argument( "Invalid type of ConstraintMod" );
 }
}

/*--------------------------------------------------------------------------*/

void MILPSolver::bound_modification( OneVarConstraintMod * mod ) {
 auto * p_const = dynamic_cast<OneVarConstraint *>(mod->constraint());
 auto * p_var = dynamic_cast<ColVariable *>(p_const->get_active_var( 0 ));
 int idx = index_of_variable( p_var );

 switch( mod->type() ) {
  case RowConstraintMod::eChgLHS:
   if( !lb.empty() ) {
    lb[ idx ] = get_problem_lb( *p_var );
   }
   break;
  case RowConstraintMod::eChgRHS:
   if( !ub.empty() ) {
    ub[ idx ] = get_problem_ub( *p_var );
   }
   break;
  case RowConstraintMod::eChgBTS:
   if( !lb.empty() ) {
    lb[ idx ] = get_problem_lb( *p_var );
   }
   if( !ub.empty() ) {
    ub[ idx ] = get_problem_ub( *p_var );
   }
   break;
  default:
   throw std::invalid_argument( "Invalid type of OneVarConstraintMod" );
 }
}

/*--------------------------------------------------------------------------*/

void MILPSolver::function_modification( FunctionMod * mod ) {
 // TODO: Change only involved variables, see function_vars_modification()

 auto * mod_f = mod->function();
 bool changing_of = false;
 const auto * lf = dynamic_cast<const LinearFunction *> (mod_f);
 const auto * qf = dynamic_cast<const DQuadFunction *> (mod_f);

 // Check if OF or a Constraint is involved
 // --------------------------------------------------------------------------

 std::queue< Block * > Q;

 // Locking the Block
 bool owned = f_Block->is_owned_by( f_id );
 if( !owned && !f_Block->read_lock() ) {
  throw std::runtime_error( "Unable to lock the Block" );
 }

 Q.push( f_Block );
 while( !Q.empty() ) {
  Block * q_Block = Q.front();
  Q.pop();

  for( auto * i : q_Block->get_nested_Blocks() ) {
   Q.push( i );
  }

  auto * p_obj = dynamic_cast< FRealObjective * >( q_Block->get_objective() );
  if( p_obj != nullptr ) {
   auto * of = p_obj->get_function();
   if( of == mod_f ) {
    changing_of = true;
    break;
   }
  }
 }

 // Unlock the Block
 if( !owned ) {
  f_Block->read_unlock();
 }

 // Change the coefficients
 // --------------------------------------------------------------------------

 if( changing_of ) {
  // Changing the coefficients of the objective function

  if( lf != nullptr && !objective.empty() ) {
   // Linear objective function
   objective.resize( lf->get_num_active_var() );
   for( auto el : lf->get_v_var() ) {
    objective[ index_of_variable( el.first ) ] = el.second;
   }

  } else if( qf != nullptr && !q_objective.empty() ) {
   // Quadratic objective function
   objective.resize( qf->get_num_active_var() );
   q_objective.resize( qf->get_num_active_var() );
   for( auto el : qf->get_v_var() ) {
    objective[ index_of_variable( std::get< 0 >( el ) ) ] =
     std::get< 1 >( el );
    q_objective[ index_of_variable( std::get< 0 >( el ) ) ] =
     std::get< 2 >( el );
   }

  } else {
   // This should never happen
   throw std::invalid_argument( "Unknown type of Objective Function" );
  }
 } else {
  // Changing coefficients of a constraint
  if( lf != nullptr ) {
   // TODO: update constraint matrix
  }
 }
}

/*--------------------------------------------------------------------------*/

void MILPSolver::function_vars_modification( FunctionModVars * mod ) {

 auto * mod_f = mod->function();
 bool changing_of = false;
 const auto * lf = dynamic_cast<const LinearFunction *> (mod_f);
 const auto * qf = dynamic_cast<const DQuadFunction *> (mod_f);

 // Check if OF or a Constraint is involved
 // --------------------------------------------------------------------------

 std::queue< Block * > Q;

 // Locking the Block
 bool owned = f_Block->is_owned_by( f_id );
 if( !owned && !f_Block->read_lock() ) {
  throw std::runtime_error( "Unable to lock the Block" );
 }

 Q.push( f_Block );
 while( !Q.empty() ) {
  Block * q_Block = Q.front();
  Q.pop();

  for( auto * i : q_Block->get_nested_Blocks() ) {
   Q.push( i );
  }

  auto * p_obj = dynamic_cast< FRealObjective * >( q_Block->get_objective() );
  if( p_obj != nullptr ) {
   auto * of = p_obj->get_function();
   if( of == mod_f ) {
    changing_of = true;
    break;
   }
  }
 }

 // Unlock the Block
 if( !owned ) {
  f_Block->read_unlock();
 }

 // Modify the coefficients
 // --------------------------------------------------------------------------

 // Check the modification type
 auto * add = dynamic_cast<C05FunctionModVarsAddd *>( mod );
 auto * rmvr = dynamic_cast<C05FunctionModVarsRngd *>( mod );
 auto * rmvs = dynamic_cast<C05FunctionModVarsSbst *>( mod );

 if( add == nullptr && rmvr == nullptr && rmvs == nullptr ) {
  throw std::invalid_argument( "This type of FunctionModVars is not handled" );
 }

 if( changing_of ) {
  // Changing the coefficients of the objective function

  if( lf != nullptr && !objective.empty() ) {
   // Linear objective function
   objective.resize( mod->vars().size() );

   for( auto * it1 : mod->vars() ) {
    for( auto it2: lf->get_v_var() ) {
     if( it1 == it2.first ) {
      if( add ) {
       objective[ index_of_variable( it2.first ) ] = it2.second;
      } else {
       objective[ index_of_variable( it2.first ) ] = 0;
      }
      break;
     }
    }
   }

  } else if( qf != nullptr && !q_objective.empty() ) {
   // Quadratic objective function
   objective.resize( mod->vars().size() );
   q_objective.resize( mod->vars().size() );

   for( auto * it1 : mod->vars() ) {
    for( auto it2: qf->get_v_var() ) {
     if( it1 == std::get< 0 >( it2 ) ) {
      if( add ) {
       objective[ index_of_variable( std::get< 0 >( it2 ) ) ] =
        std::get< 1 >( it2 );
       q_objective[ index_of_variable( std::get< 0 >( it2 ) ) ] =
        std::get< 2 >( it2 );
      } else {
       objective[ index_of_variable( std::get< 0 >( it2 ) ) ] = 0;
       q_objective[ index_of_variable( std::get< 0 >( it2 ) ) ] = 0;
      }
      break;
     }
    }
   }
  } else {
   // This should never happen
   throw std::invalid_argument( "Unknown type of Objective Function" );
  }
 } else {
  // Changing coefficients of a constraint
  if( lf != nullptr ) {
   // TODO: update constraint matrix
  }
 }
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

 // Update the dictionaries
 auto it = lower_bound( dcon_to_idx.begin(), dcon_to_idx.end(), p_const,
                        [ & ]( const_int pair, FRowConstraint * c ) {
                         return pair.first < c;
                        } );
 dcon_to_idx.insert( it, { p_const, numrows } );
 idx_to_dcon.emplace_back( numrows, p_const );

 auto const_lhs = p_const->get_lhs();
 auto const_rhs = p_const->get_rhs();

 if( const_lhs == const_rhs ) {
  if( !sense.empty() ) {
   sense.emplace_back( 'E' );
  }
  if( !rhs.empty() ) {
   rhs.emplace_back( const_rhs );
  }
  if( !rngval.empty() ) {
   rngval.emplace_back( 0 );
  }

 } else if( const_lhs == -Inf< double >() ) {
  if( !sense.empty() ) {
   sense.emplace_back( 'L' );
  }
  if( !rhs.empty() ) {
   rhs.emplace_back( const_rhs );
  }
  if( !rngval.empty() ) {
   rngval.emplace_back( 0 );
  }

 } else if( const_rhs == Inf< double >() ) {
  if( !sense.empty() ) {
   sense.emplace_back( 'G' );
  }
  if( !rhs.empty() ) {
   rhs.emplace_back( const_lhs );
  }
  if( !rngval.empty() ) {
   rngval.emplace_back( 0 );
  }

 } else {
  if( !sense.empty() ) {
   sense.emplace_back( 'R' );
  }
  if( !rhs.empty() ) {
   rhs.emplace_back( const_lhs );
  }
  if( !rngval.empty() ) {
   rngval.emplace_back( const_rhs - const_lhs );
  }
 }

 // TODO: update constraint matrix
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

 // Update the number of integer vars
 if( p_var->is_integer() ) {
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

 int idx = numcols;
 if( !lb.empty() ) {
  lb.emplace_back( get_problem_lb( *p_var ) );
 }
 if( !ub.empty() ) {
  ub.emplace_back( get_problem_ub( *p_var ) );
 }

 if( !xctype.empty() ) {
  if( p_var->is_integer() ) {
   if( p_var->is_unitary() && p_var->is_positive() ) {
    xctype.emplace_back( 'B' );
   } else {
    xctype.emplace_back( 'I' );
   }
  } else {
   xctype.emplace_back( 'C' );
  }
 }

 // TODO: update constraint matrix
 ++numcols;
}

/*--------------------------------------------------------------------------*/

void MILPSolver::add_dynamic_bound( OneVarConstraint * p_bound ) {

 auto * p_var = dynamic_cast<ColVariable *>(p_bound->get_active_var( 0 ));
 if( p_var != nullptr ) {
  int idx = index_of_variable( p_var );
  if( !lb.empty() ) {
   lb[ idx ] = get_problem_lb( *p_var );
  }
  if( !ub.empty() ) {
   ub[ idx ] = get_problem_ub( *p_var );
  }
 }
}

/*--------------------------------------------------------------------------*/

void MILPSolver::remove_dynamic_constraint( const FRowConstraint * p_const ) {
 // TODO: Implement remove_dynamic_with_index(i)

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
  throw std::runtime_error( "Dynamic constraint not found" );
 }

 for( auto & it: dcon_to_idx ) {
  if( it.second > index ) {
   it.second--;
  }
 }

 auto it2 = lower_bound( idx_to_dcon.begin(), idx_to_dcon.end(),
                         std::make_pair( index, 0 ),
                         [ & ]( auto & p1, auto & p2 ) {
                          return p1.first < p2.first;
                         } );

 if( it2 != idx_to_dcon.end() && it2->second == p_const ) {
  for( auto it = it2; it != idx_to_dcon.end(); ++it ) {
   it->first--;
  }
  idx_to_dcon.erase( it2 );
 } else {
  throw std::runtime_error( "Dynamic constraint not found" );
 }

 if( !sense.empty() ) {
  sense.erase( sense.begin() + index );
 }
 if( !rhs.empty() ) {
  rhs.erase( rhs.begin() + index );
 }
 if( !rngval.empty() ) {
  rngval.erase( rngval.begin() + index );
 }

 // TODO: update constraint matrix
 --numrows;
}

/*--------------------------------------------------------------------------*/

void MILPSolver::remove_dynamic_variable( const ColVariable * p_var ) {
 // TODO: Implement remove_dynamic_with_index(i)

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
  throw std::runtime_error( "Dynamic variable not found" );
 }

 for( auto & it: dvar_to_idx ) {
  if( it.second > index ) {
   it.second--;
  }
 }

 auto it2 = lower_bound( idx_to_dvar.begin(), idx_to_dvar.end(),
                         std::make_pair( index, 0 ),
                         [ & ]( auto & p1, auto & p2 ) {
                          return p1.first < p2.first;
                         } );

 if( it2 != idx_to_dvar.end() && it2->second == p_var ) {
  for( auto it = it2; it != idx_to_dvar.end(); ++it ) {
   it->first--;
  }
  idx_to_dvar.erase( it2 );
 } else {
  throw std::runtime_error( "Dynamic constraint not found" );
 }

 if( !lb.empty() ) {
  lb.erase( lb.begin() + index );
 }
 if( !ub.empty() ) {
  ub.erase( ub.begin() + index );
 }
 if( !xctype.empty() ) {
  xctype.erase( xctype.begin() + index );
 }

 // TODO: update constraint matrix

 --numcols;

 // Update the number of integer vars
 if( p_var->is_integer() ) {
  --int_vars;
 }
}

/*--------------------------------------------------------------------------*/

void MILPSolver::remove_dynamic_bound( const OneVarConstraint * p_bound ) {

 auto * p_var = dynamic_cast<ColVariable *>(p_bound->get_active_var( 0 ));
 if( p_var != nullptr ) {
  int idx = index_of_variable( p_var );
  if( !lb.empty() ) {
   lb[ idx ] = get_problem_lb( *p_var );
  }
  if( !ub.empty() ) {
   ub[ idx ] = get_problem_ub( *p_var );
  }
 }
}

/*--------------------------------------------------------------------------*/

template< typename T >
std::string MILPSolver::log_vector( const std::vector< T > & v, int limit ) {
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
std::string MILPSolver::log_vector( const std::vector< char > & v, int limit ) {
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

void MILPSolver::check_status() {
 int v = 0;
 int c = 0;
 int sv = 0;
 int sc = 0;
 int svg = 0;
 int scg = 0;
 int dv = 0;
 int dc = 0;

 DEBUG_LOG( "Checking MILPSolver dictionaries" << std::endl );

 // ------------------ Count everything -------------------
 std::queue< Block * > Q;

 // Locking the Block
 bool owned = f_Block->is_owned_by( f_id );
 if( !owned && !f_Block->read_lock() ) {
  throw std::runtime_error( "Unable to lock the Block" );
 }

 Q.push( f_Block );
 while( !Q.empty() ) {
  Block * q_Block = Q.front();
  Q.pop();

  for( auto * i : q_Block->get_nested_Blocks() ) {
   Q.push( i );
  }

  for( const auto & i : q_Block->get_static_constraints() ) {
   // Singles
   if( un_any_thing_0( FRowConstraint, i,
                       {
                        ++scg;
                        ++c;
                        ++sc;
                       }
   ) ) {
    continue;
   }
   // Vectors
   if( un_any_thing_1( FRowConstraint, i,
                       {
                        ++scg;
                        c += var.size();
                        sc += var.size();
                       }
   ) ) {
    continue;
   }
   // Multiarrays
   if( un_any_thing_K( FRowConstraint, i,
                       {
                        ++scg;
                        c += var.num_elements();
                        sc += var.num_elements();
                       }
   ) ) {
    continue;
   }
  }

  for( const auto & i : q_Block->get_dynamic_constraints() ) {
   // Single lists
   if( un_any_thing_0( std::list< FRowConstraint >, i,
                       {
                        c += var.size();
                       }
   ) ) {
    continue;
   }
   // Vectors of lists
   if( un_any_thing_1( std::list< FRowConstraint >, i,
                       {
                        for( auto & el: var ) {
                         c += el.size();
                        }
                       }
   ) ) {
    continue;
   }
   // Multiarrays of lists
   if( un_any_thing_K( std::list< FRowConstraint >, i,
                       {
                        auto it = var.data();
                        for( auto i = var.num_elements(); i--; ++it ) {
                         c += it->size();
                        }
                       }
   ) ) {
    continue;
   }
  }
  dc = c - sc;

  for( const auto & i : q_Block->get_static_variables() ) {
   // Singles
   if( un_any_thing_0( ColVariable, i,
                       {
                        ++svg;
                        ++v;
                        ++sv;
                       }
   ) ) {
    continue;
   }
   // Vectors
   if( un_any_thing_1( ColVariable, i,
                       {
                        ++svg;
                        v += var.size();
                        sv += var.size();
                       }
   ) ) {
    continue;
   }
   // Multiarrays
   if( un_any_thing_K( ColVariable, i,
                       {
                        ++svg;
                        v += var.num_elements();
                        sv += var.num_elements();
                       }
   ) ) {
    continue;
   }
  }

  for( const auto & i : q_Block->get_dynamic_variables() ) {
   // Single lists
   if( un_any_thing_0( std::list< ColVariable >, i,
                       {
                        v += var.size();
                       }
   ) ) {
    continue;
   }
   // Vectors of lists
   if( un_any_thing_1( std::list< ColVariable >, i,
                       {
                        for( auto & el: var ) {
                         v += el.size();
                        }
                       }
   ) ) {
    continue;
   }
   // Multiarrays of lists
   if( un_any_thing_K( std::list< ColVariable >, i,
                       {
                        auto it = var.data();
                        for( auto i = var.num_elements(); i--; ++it ) {
                         v += it->size();
                        }
                       }
   ) ) {
    continue;
   }
  }
  dv = v - sv;
 }

 // Unlock the Block
 if( !owned ) {
  f_Block->read_unlock();
 }

 if( numcols != v ) {
  DEBUG_LOG( "numcols is " << numcols << ", it should be " << v << std::endl );
 }
 if( numrows != c ) {
  DEBUG_LOG( "numrows is " << numrows << ", it should be " << c << std::endl );
 }
 if( static_vars != sv ) {
  DEBUG_LOG( "static_vars is " << static_vars
                               << ", it should be " << sv << std::endl );
 }
 if( static_cons != sc ) {
  DEBUG_LOG( "static_cons is " << static_cons
                               << ", it should be " << sc << std::endl );
 }

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
                         [ & ]( auto & pair ) {
                          return pair.second == i.first &&
                                 pair.first == i.second;
                         } );
  if( j == dvar_to_idx.end() ) {
   DEBUG_LOG( "Element [" << i.first << ", " << i.second
                          << "] of idx_to_dvar was not found in dvar_to_idx"
                          << std::endl );
  }
 }

 for( auto & i: dvar_to_idx ) {
  auto j = std::find_if( idx_to_dvar.begin(), idx_to_dvar.end(),
                         [ & ]( auto & pair ) {
                          return i.first == pair.second &&
                                 i.second == pair.first;
                         } );
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
                         [ & ]( auto & pair ) {
                          return pair.second == i.first &&
                                 pair.first == i.second;
                         } );
  if( j == dcon_to_idx.end() ) {
   DEBUG_LOG( "Element [" << i.first << ", " << i.second
                          << "] of idx_to_dcon was not found in dcon_to_idx"
                          << std::endl );
  }
 }
 for( auto & i: dcon_to_idx ) {
  auto j = std::find_if( idx_to_dcon.begin(), idx_to_dcon.end(),
                         [ & ]( auto & pair ) {
                          return i.first == pair.second &&
                                 i.second == pair.first;
                         } );
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
