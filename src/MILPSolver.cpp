/*--------------------------------------------------------------------------*/
/*-------------------------- File MILPSolver.cpp ---------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Implementation of the MILPSolver class.
 *
 * \version 0.90
 *
 * \date 14 - 06 - 2019
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

#include <cstdlib>
#include <functional>
#include <queue>

#include <Block.h>
#include <OneVarConstraint.h>
#include <LinearFunction.h>
#include <DQuadFunction.h>

#include "MILPSolver.h"

// TODO: Remove this and all the printouts when done
#ifdef MILPSLVR_DEBUG
#define LOG( stuff ) std::cout << stuff
#define LOG_VEC( stuff ) std::cout << "[";      \
                       for (auto i : stuff)    \
                        std::cout << " " << i; \
                       std::cout << "]\n"
#else
#define LOG(stuff)
#define LOG_VEC(stuff)
#endif

/*--------------------------------------------------------------------------*/
/*------------------------- NAMESPACE AND USING ----------------------------*/
/*--------------------------------------------------------------------------*/

using namespace SMSpp_di_unipi_it;

SMSpp_insert_in_factory_cpp_0( MILPSolver );

/*--------------------------------------------------------------------------*/
/*--------------------- CONSTRUCTOR AND DESTRUCTOR -------------------------*/
/*--------------------------------------------------------------------------*/

MILPSolver::MILPSolver() : CDASolver() {}

MILPSolver::~MILPSolver() = default;

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

const std::vector< double > & MILPSolver::get_lb() const {
 return lb;
}

const std::vector< double > & MILPSolver::get_ub() const {
 return ub;
}

const std::vector< char > & MILPSolver::get_xctype() const {
 return xctype;
}

int MILPSolver::get_nodes() const {
 return nodes;
}

/*--------------------------------------------------------------------------*/
/*-------------------------------- SET_BLOCK -------------------------------*/
/*--------------------------------------------------------------------------*/

void MILPSolver::set_Block( Block * block ) {
 if( block == f_Block ) {
  return;
 }
 Solver::set_Block( block );
 if( block) {
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

 v_s_var_int.clear();
 v_int_s_var.clear();
 v_s_const_int.clear();
 v_int_s_const.clear();
 v_d_var_int.clear();
 v_int_d_var.clear();
 v_d_const_int.clear();
 v_int_d_const.clear();

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

 // First loop on the queue to count variables and constraints
 Q.push( f_Block );

 while( !Q.empty() ) {
  Block * q_Block = Q.front();
  Q.pop();

  for( auto i : q_Block->get_nested_Blocks() ) {
   Q.push( i );
  }

  LOG( "[DEBUG] ========= MILPSolver::set_Block() counting static constraints\n" );
  for( const auto & i : q_Block->get_static_constraints() ) {
   auto f1 = std::bind( &MILPSolver::count_constraints,
                        this,
                        std::placeholders::_1,
                        std::ref( numrows ) );
   un_any_const_static( i, f1, un_any_type< FRowConstraint >() );
  }

  LOG( "[DEBUG] ========= MILPSolver::set_Block() counting dynamic constraints\n" );
  for( const auto & i : q_Block->get_dynamic_constraints() ) {
   auto f1 = std::bind( &MILPSolver::count_constraints,
                        this,
                        std::placeholders::_1,
                        std::ref( numrows ) );
   un_any_const_dynamic( i, f1, un_any_type< FRowConstraint >() );
  }

  LOG( "[DEBUG] ========= MILPSolver::set_Block() counting static variables\n" );
  for( const auto & i : q_Block->get_static_variables() ) {
   auto f1 = std::bind( &MILPSolver::count_variables,
                        this,
                        std::placeholders::_1,
                        std::ref( numcols ) );
   un_any_const_static( i, f1, un_any_type< ColVariable >() );
  }

  LOG( "[DEBUG] ========= MILPSolver::set_Block() counting dynamic variables\n" );
  for( const auto & i : q_Block->get_dynamic_variables() ) {
   auto f1 = std::bind( &MILPSolver::count_variables,
                        this,
                        std::placeholders::_1,
                        std::ref( numcols ) );
   un_any_const_dynamic( i, f1, un_any_type< ColVariable >() );
  }

  LOG( "[DEBUG] ========= MILPSolver::set_Block() nonzero elements\n" );
  int cnt = 0;
  for( const auto & i : q_Block->get_static_variables() ) {
   auto f1 = std::bind( &MILPSolver::count_nzelements,
                        this,
                        std::placeholders::_1,
                        std::ref( nzelements ),
                        std::ref( cnt ) );
   un_any_const_static( i, f1, un_any_type< ColVariable >() );
  }

  for( const auto & i : q_Block->get_dynamic_variables() ) {
   auto f1 = std::bind( &MILPSolver::count_nzelements,
                        this,
                        std::placeholders::_1,
                        std::ref( nzelements ),
                        std::ref( cnt ) );
   un_any_const_dynamic( i, f1, un_any_type< ColVariable >() );
  }
 } // End of while loop on Block queue

 LOG( "[DEBUG] ========= MILPSolver::set_Block() after counting\n" );
 LOG( "constraints/numrows = " << numrows << std::endl );
 LOG( "variables/numcols =   " << numcols << std::endl );
 LOG( "nzelements =          " << nzelements << std::endl );

 // The +1 is needed by generic interface
 matbeg.resize( numcols + 1 );
 matbeg[ numcols ] = nzelements;

 matcnt.resize( numcols );
 matind.resize( nzelements );
 matval.resize( nzelements );
 rhs.resize( numrows );
 rngval.resize( numrows );
 sense.resize( numrows );
 objective.resize( numcols );
 q_objective.resize( numcols );
 std::fill( objective.begin(), objective.end(), 0 );
 std::fill( objective.begin(), objective.end(), 0 );
 lb.resize( numcols );
 ub.resize( numcols );
 xctype.resize( numcols );

 // Second loop to scan the constraints
 Q.push( f_Block );

 int col = 0;

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

  for( const auto & i : q_Block->get_static_constraints() ) {
   // Variable used to locate the first element of each type of constraint
   int first = 0;
   auto f1 = std::bind( &MILPSolver::scan_static_constraint,
                        this,
                        std::placeholders::_1,
                        std::ref( first ),
                        std::ref( col ) );
   un_any_const_static( i, f1, un_any_type< FRowConstraint >() );
  }

  for( const auto & i : q_Block->get_dynamic_constraints() ) {
   auto f1 = std::bind( &MILPSolver::scan_dynamic_constraint,
                        this,
                        std::placeholders::_1,
                        std::ref( col ) );
   un_any_const_dynamic( i, f1, un_any_type< FRowConstraint >() );
  }
 } // End of while loop on Block queue

 std::sort( v_s_const_int.begin(), v_s_const_int.end() );
 std::sort( v_int_s_const.begin(), v_int_s_const.end() );
 std::sort( v_d_const_int.begin(), v_d_const_int.end() );
 std::sort( v_int_d_const.begin(), v_int_d_const.end() );

 LOG( "[DEBUG] ========= MILPSolver::set_Block() after constraint scan\n" );
 // Third loop to scan the variables
 Q.push( f_Block );

 int var = 0;
 while( !Q.empty() ) {
  Block * q_Block = Q.front();
  Q.pop();

  for( auto i : q_Block->get_nested_Blocks() ) {
   Q.push( i );
  }

  for( const auto & i : q_Block->get_static_variables() ) {
   int first = 0;
   auto f1 = std::bind( &MILPSolver::scan_static_variable,
                        this,
                        std::placeholders::_1,
                        std::ref( first ),
                        std::ref( var ) );
   un_any_const_static( i, f1, un_any_type< ColVariable >() );
  }

  for( const auto & i : q_Block->get_dynamic_variables() ) {
   auto f1 = std::bind( &MILPSolver::scan_dynamic_variable,
                        this,
                        std::placeholders::_1,
                        std::ref( var ) );
   un_any_const_dynamic( i, f1, un_any_type< ColVariable >() );
  }
 } // End of while loop on Block queue

 std::sort( v_s_var_int.begin(), v_s_var_int.end() );
 std::sort( v_int_s_var.begin(), v_int_s_var.end() );
 std::sort( v_d_var_int.begin(), v_d_var_int.end() );
 std::sort( v_int_d_var.begin(), v_int_d_var.end() );

 const FRealObjective * p_obj;
 p_obj = dynamic_cast< FRealObjective * >( f_Block->get_objective() );
 switch( p_obj->get_sense() ) {
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

  for( auto i : q_Block->get_nested_Blocks() ) {
   Q.push( i );
  }

  p_obj = dynamic_cast< FRealObjective * >( q_Block->get_objective() );
  scan_objective( p_obj );
 } // End of while loop on Block queue

 LOG( "[DEBUG] ========= MILPSolver::set_Block() after objective scan\n" );
 LOG( "[DEBUG] objective   = " );
 LOG_VEC( objective );
 LOG( "[DEBUG] q_objective = " );
 LOG_VEC( q_objective );
 LOG( "[DEBUG] rhs         = " );
 LOG_VEC( rhs );
 LOG( "[DEBUG] rngval      = " );
 LOG_VEC( rngval );
 LOG( "[DEBUG] sense       = " );
 LOG_VEC( sense );
 LOG( "[DEBUG] matbeg      = " );
 LOG_VEC( matbeg );
 LOG( "[DEBUG] matcnt      = " );
 LOG_VEC( matcnt );
 LOG( "[DEBUG] matind      = " );
 LOG_VEC( matind );
 LOG( "[DEBUG] matval      = " );
 LOG_VEC( matval );
 LOG( "[DEBUG] lb          = " );
 LOG_VEC( lb );
 LOG( "[DEBUG] ub          = " );
 LOG_VEC( ub );
 LOG( "[DEBUG] xctype      = " );
 LOG_VEC( xctype );
}

/*--------------------------------------------------------------------------*/
/*-------------------- METHODS FOR PROBLEM DESCRIPTION ---------------------*/
/*--------------------------------------------------------------------------*/

int MILPSolver::index_of_variable( ColVariable * p_var ) {
 try {
  return index_of_dynamic_variable( p_var );
 } catch( ... ) {
  return index_of_static_variable( p_var );
 }
}

int MILPSolver::index_of_static_variable( ColVariable * p_var ) {

 int i;
 auto it = lower_bound( v_s_var_int.begin(),
                        v_s_var_int.end(),
                        p_var,
                        [ & ]( var_int pair, ColVariable * pvar ) {
                         return pair.first < pvar;
                        } );

 if( it == v_s_var_int.end() ) {
  it = ( v_s_var_int.rbegin() + 1 ).base();
 } else if( it != v_s_var_int.begin() && it->first > p_var ) {
  --it;
 }

 if( it != v_s_var_int.end() ) {
  i = static_cast<int>(it->second + std::distance( it->first, p_var ));
 } else {
  throw ( std::invalid_argument( "Variable not found" ) );
 }
 return i;
}

int MILPSolver::index_of_dynamic_variable( ColVariable * p_var ) {
 auto it = find_if( v_d_var_int.begin(),
                    v_d_var_int.end(),
                    [ & ]( var_int pair ) {
                     return pair.first == p_var;
                    } );
 if( it != v_d_var_int.end() ) {
  return it->second;
 } else {
  throw ( std::invalid_argument( "Variable not found" ) );
 }
}

/*--------------------------------------------------------------------------*/
int MILPSolver::index_of_constraint( FRowConstraint * p_const ) {
 try {
  return index_of_dynamic_constraint( p_const );
 } catch( ... ) {
  return index_of_static_constraint( p_const );
 }
}

int MILPSolver::index_of_static_constraint( FRowConstraint * p_const ) {

 int i = 0;

 auto it = lower_bound( v_s_const_int.begin(),
                        v_s_const_int.end(),
                        p_const,
                        [ & ]( const_int pair, FRowConstraint * pconst ) {
                         return pair.first < pconst;
                        } );

 if( it == v_s_const_int.end() ) {
  it = ( v_s_const_int.rbegin() + 1 ).base();
 } else if( it != v_s_const_int.begin() && it->first > p_const ) {
  --it;
 }

 if( it != v_s_const_int.end() ) {
  i = static_cast<int>(it->second + std::distance( it->first, p_const ));
 } else {
  throw ( std::invalid_argument( "Constraint not found" ) );
 }

 return i;
}

int MILPSolver::index_of_dynamic_constraint( FRowConstraint * p_const ) {
 auto it = find_if( v_d_const_int.begin(),
                    v_d_const_int.end(),
                    [ & ]( const_int pair ) {
                     return pair.first == p_const;
                    } );
 if( it != v_d_const_int.end() ) {
  return it->second;
 } else {
  throw ( std::invalid_argument( "Constraint not found" ) );
 }
}

/*--------------------------------------------------------------------------*/

ColVariable * MILPSolver::static_variable_with_index( int i ) {

 ColVariable * p_var;
 auto it = lower_bound( v_int_s_var.begin(),
                        v_int_s_var.end(),
                        i,
                        [ & ]( int_var pair, int i ) {
                         return pair.first < i;
                        } );

 if( it == v_int_s_var.end() ) {
  it = ( v_int_s_var.rbegin() + 1 ).base();
 } else if( it != v_int_s_var.begin() && it->first > i ) {
  --it;
 }

 if( it != v_int_s_var.end() ) {
  p_var = it->second;
 } else {
  throw ( std::invalid_argument( "Index not found" ) );
 }

 return p_var;
}

/*--------------------------------------------------------------------------*/

FRowConstraint * MILPSolver::static_constraint_with_index( int i ) {

 FRowConstraint * p_const;
 auto it = lower_bound( v_int_s_const.begin(),
                        v_int_s_const.end(),
                        i,
                        [ & ]( int_const pair, int i ) {
                         return pair.first < i;
                        } );

 if( it == v_int_s_const.end() ) {
  it = ( v_int_s_const.rbegin() + 1 ).base();
 } else if( it != v_int_s_const.begin() && it->first > i ) {
  --it;
 }

 if( it != v_int_s_const.end() ) {
  p_const = it->second;
 } else {
  throw ( std::invalid_argument( "Index not found" ) );
 }

 return p_const;
}

/*--------------------------------------------------------------------------*/
/*------------- AUXILIARY METHODS FOR POPULATING THE PROBLEM  --------------*/
/*--------------------------------------------------------------------------*/

void MILPSolver::count_constraints( FRowConstraint & constraint, int & n_rows ) {
 LOG( "[DEBUG] ========= MILPSolver::count_constraints()  " << n_rows << " " << constraint );
 auto fun = dynamic_cast<const LinearFunction *>(constraint.get_function());
 if( fun != nullptr ) {
  ++n_rows;
 } else {
  throw ( std::invalid_argument( "The Constraint is not linear" ) );
 }
}

/*--------------------------------------------------------------------------*/

void MILPSolver::count_variables( ColVariable & variable, int & n_cols ) {
 LOG( "[DEBUG] ========= MILPSolver::count_variables()    " << n_cols << " " << variable );
 ++n_cols;
}

/*--------------------------------------------------------------------------*/

void MILPSolver::count_nzelements( ColVariable & variable,
                                   int & nz_elements,
                                   int & cnt ) {
 LOG( "[DEBUG] ========= MILPSolver::count_nzelements()   " << cnt << " " << variable );
 LOG( "[DEBUG] The active stuff is:\n" );

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

 for( auto i : variable.active_stuff() ) {
  auto row = dynamic_cast<FRowConstraint *>(i);
  if( row != nullptr ) {
   LOG( "[DEBUG] " << *row );
   active_constraints[ cnt ].push_back( row );
   ++nz_elements;
  }
  auto box = dynamic_cast<OneVarConstraint *>(i);
  if( box != nullptr ) {
   LOG( "[DEBUG] " << *box );
   active_bounds[ cnt ].push_back( box );
  }
  auto obj = dynamic_cast<Objective *>(i);
  if( obj != nullptr ) {
   LOG( "[DEBUG] " << *obj );
  }
 }
 ++cnt;
}

/*--------------------------------------------------------------------------*/

void MILPSolver::scan_static_variable( ColVariable & var, int & first, int & i ) {
 LOG( "[DEBUG] ========= MILPSolver::scan_static_variable() " );
 LOG( i << " " << var );


 if( first == 0 ) {
  v_s_var_int.emplace_back( &var, i );
  v_int_s_var.emplace_back( i, &var );
 }

 lb[ i ] = var.get_lb();
 ub[ i ] = var.get_ub();

 int num_bounds = static_cast<int>(active_bounds[ i ].size());
 for( int j = 0; j < num_bounds; ++j ) {
  auto bound = active_bounds[ i ][ j ];
  lb[ i ] = lb[ i ] > bound->get_lhs() ? lb[ i ] : bound->get_lhs();
  ub[ i ] = ub[ i ] < bound->get_rhs() ? ub[ i ] : bound->get_rhs();
 }

 if( var.is_integer() ) {
  if( var.is_unitary() && var.is_positive() ) {
   xctype[ i ] = 'B'; // Binary
  } else {
   xctype[ i ] = 'I'; // Integer
  }
 } else {
  xctype[ i ] = 'C';  // Continuous
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

 int nz_elements = static_cast<int>(active_constraints[ i ].size());
 matcnt[ i ] = nz_elements;

 if( i == 0 ) {
  matbeg[ i ] = 0;
 } else {
  matbeg[ i ] = matbeg[ i - 1 ] + matcnt[ i - 1 ];
 }

 for( int j = 0; j < nz_elements; ++j ) {

  auto p_const = dynamic_cast<FRowConstraint *> (active_constraints[ i ][ j ]);
  auto p_fun = dynamic_cast<const LinearFunction *> (p_const->get_function());

  auto it = std::find_if( p_fun->get_v_var().begin(),
                          p_fun->get_v_var().end(),
                          [ & ]( LinearFunction::coeff_pair pair ) {
                           return pair.first == &var;
                          } );

  if( it != p_fun->get_v_var().end() ) {
   matval[ matbeg[ i ] + j ] = it->second;
   matind[ matbeg[ i ] + j ] = index_of_constraint( p_const );
  } else {
   // This should never happen because we are looping on the active contraints
   throw ( std::invalid_argument( "This ColVariable is not active in the examined FRowConstraint" ) );
  }
 }
 ++first;
 ++i;
}

/*--------------------------------------------------------------------------*/

void MILPSolver::scan_dynamic_variable( ColVariable & var, int & i ) {
 LOG( "[DEBUG] ========= MILPSolver::scan_dynamic_variable() " );
 LOG( i << " " << var );

 v_d_var_int.emplace_back( &var, i );
 v_int_d_var.emplace_back( i, &var );

 lb[ i ] = var.get_lb();
 ub[ i ] = var.get_ub();

 int num_bounds = static_cast<int>(active_bounds[ i ].size());
 for( int j = 0; j < num_bounds; ++j ) {
  auto bound = active_bounds[ i ][ j ];
  lb[ i ] = lb[ i ] > bound->get_lhs() ? lb[ i ] : bound->get_lhs();
  ub[ i ] = ub[ i ] < bound->get_rhs() ? ub[ i ] : bound->get_rhs();
 }

 if( var.is_integer() ) {
  if( var.is_unitary() && var.is_positive() ) {
   xctype[ i ] = 'B'; // Binary
  } else {
   xctype[ i ] = 'I'; // Integer
  }
 } else {
  xctype[ i ] = 'C';  // Continuous
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

 int nz_elements = static_cast<int>(active_constraints[ i ].size());
 matcnt[ i ] = nz_elements;

 if( i == 0 ) {
  matbeg[ i ] = 0;
 } else {
  matbeg[ i ] = matbeg[ i - 1 ] + matcnt[ i - 1 ];
 }

 for( int j = 0; j < nz_elements; ++j ) {

  auto p_const = dynamic_cast<FRowConstraint *> (active_constraints[ i ][ j ]);
  auto p_fun = dynamic_cast<const LinearFunction *> (p_const->get_function());

  auto it = std::find_if( p_fun->get_v_var().begin(),
                          p_fun->get_v_var().end(),
                          [ & ]( LinearFunction::coeff_pair pair ) {
                           return pair.first == &var;
                          } );

  if( it != p_fun->get_v_var().end() ) {
   matval[ matbeg[ i ] + j ] = it->second;
   matind[ matbeg[ i ] + j ] = index_of_constraint( p_const );
  } else {
   // This should never happen because we are looping on the active contraints
   throw ( std::invalid_argument( "This ColVariable is not active in the examined FRowConstraint" ) );
  }
 }
 ++i;
}

/*--------------------------------------------------------------------------*/

void MILPSolver::scan_static_constraint( FRowConstraint & p_const, int & first, int & i ) {
 LOG( "[DEBUG] ========= MILPSolver::scan_static_constraint()  " << p_const );

 auto lin_fun = dynamic_cast<const LinearFunction *>(p_const.get_function());
 if( lin_fun == nullptr ) {
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

 auto const_lhs = p_const.get_lhs();
 auto const_rhs = p_const.get_rhs();

 if( const_lhs == const_rhs ) {
  // LHS <= function <= RHS, with LHS = RHS
  // becomes:
  // function = RHS
  sense[ i ] = 'E';
  rhs[ i ] = const_rhs;

 } else if( const_lhs == -Inf< double >() ) {
  // -inf <= function <= RHS
  // becomes:
  // function <= RHS
  sense[ i ] = 'L';
  rhs[ i ] = const_rhs;

 } else if( const_rhs == Inf< double >() ) {
  // LHS <= function <= inf
  // becomes:
  // function >= LHS
  sense[ i ] = 'G';
  rhs[ i ] = const_lhs;

 } else {
  // LHS <= function <= RHS
  // becomes:
  // LHS <= function <= LHS + (range),
  // with range = RHS - LHS
  sense[ i ] = 'R';
  rhs[ i ] = const_lhs;
  rngval[ i ] = const_rhs - const_lhs;
 }

 if( first == 0 ) {
  v_s_const_int.emplace_back( &p_const, i );
  v_int_s_const.emplace_back( i, &p_const );
 }

 ++first;
 ++i;
}

/*--------------------------------------------------------------------------*/

void MILPSolver::scan_dynamic_constraint( FRowConstraint & p_const, int & i ) {
 LOG( "[DEBUG] ========= MILPSolver::scan_dynamic_constraint() " << p_const );

 auto lin_fun = dynamic_cast<const LinearFunction *>(p_const.get_function());
 if( lin_fun == nullptr ) {
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

 auto const_lhs = p_const.get_lhs();
 auto const_rhs = p_const.get_rhs();

 if( const_lhs == const_rhs ) {
  // LHS <= function <= RHS, with LHS = RHS
  // becomes:
  // function = RHS
  sense[ i ] = 'E';
  rhs[ i ] = const_rhs;

 } else if( const_lhs == -Inf< double >() ) {
  // -inf <= function <= RHS
  // becomes:
  // function <= RHS
  sense[ i ] = 'L';
  rhs[ i ] = const_rhs;

 } else if( const_rhs == Inf< double >() ) {
  // LHS <= function <= inf
  // becomes:
  // function >= LHS
  sense[ i ] = 'G';
  rhs[ i ] = const_lhs;

 } else {
  // LHS <= function <= RHS
  // becomes:
  // LHS <= function <= LHS + (range),
  // with range = RHS - LHS
  sense[ i ] = 'R';
  rhs[ i ] = const_lhs;
  rngval[ i ] = const_rhs - const_lhs;
 }

 v_d_const_int.emplace_back( &p_const, i );
 v_int_d_const.emplace_back( i, &p_const );
 ++i;
}

/*--------------------------------------------------------------------------*/

void MILPSolver::scan_objective( const FRealObjective * obj ) {
 LOG( "[DEBUG] ========= MILPSolver::scan_objective() " << *obj );

 auto lin_fun = dynamic_cast<const LinearFunction *> (obj->get_function());
 int k;

 if( lin_fun != nullptr ) {
  for( auto el : lin_fun->get_v_var() ) {
   k = index_of_variable( el.first );
   objective[ k ] = el.second;
  }
 } else {
  auto dquad_fun = dynamic_cast<const DQuadFunction *> (obj->get_function());
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

 while( !v_mod.empty() ) {
  auto mod = v_mod.front();

  // A function like this is needed to be called recursively with GroupModifications
  std::function< void( sp_Mod ) > f;
  f = [ this, &f ]( sp_Mod mod ) {

   {
    const auto tmod = std::dynamic_pointer_cast< GroupModification >( mod );
    if( tmod ) {
     LOG("GroupModification containing: " << std::endl);
     for( const auto & submod : tmod->v_sub_Modifications ) {
      f( submod );
     }
     return;
    }
   }
   {
    const auto tmod = std::dynamic_pointer_cast< VariableMod >( mod );
    if( tmod ) {
     LOG(*mod);
     var_modification( tmod.get() );
     return;
    }
   }
   {
    const auto tmod = std::dynamic_pointer_cast< ObjectiveMod >( mod );
    if( tmod ) {
     LOG(*mod);
     of_modification( tmod.get() );
     return;
    }
   }
   {
    const auto tmod = std::dynamic_pointer_cast< OneVarConstraintMod >( mod );
    if( tmod ) {
     LOG(*mod);
     bound_modification( tmod.get() );
     return;
    }
   }
   {
    const auto tmod = std::dynamic_pointer_cast< RowConstraintMod >( mod );
    if( tmod ) {
     LOG(*mod);
     const_modification( tmod.get() );
     return;
    }
   }
   {
    const auto tmod = std::dynamic_pointer_cast< ConstraintMod >( mod );
    if( tmod ) {
     LOG(*mod);
     const_modification( tmod.get() );
     return;
    }
   }
   {
    const auto tmod = std::dynamic_pointer_cast< FunctionMod >( mod );
    if( tmod ) {
     LOG(*mod);
     function_modification( tmod.get() );
     return;
    }
   }
   {
    const auto tmod = std::dynamic_pointer_cast< BlockModAD >( mod );
    if( tmod ) {
     LOG(*mod);
     dynamic_modification( tmod.get() );
     return;
    }
   }
   {
    const auto tmod = std::dynamic_pointer_cast< NBModification >( mod );
    if( tmod ) {
     LOG("\033[1;33m" << *mod << "\033[0m");
     clear_problem();
     load_problem();
    }
   }
  };

  f( mod );
  v_mod.pop_front();
 }
}

/*--------------------------------------------------------------------------*/
/*----------------------- End File MILPSolver.cpp --------------------------*/
/*--------------------------------------------------------------------------*/ 
