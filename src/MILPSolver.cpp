/*--------------------------------------------------------------------------*/
/*-------------------------- File MILPSolver.cpp ---------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Implementation of the MILPSolver class.
 *
 * \version 0.10
 *
 * \date 22 - 12 - 2016
 *
 * \author Antonio Frangioni \n
 *         Operations Research Group \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \author Kostas Tavlaridis-Gyparakis \n
 *         Operations Research Group \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * Copyright &copy by Antonio Frangioni, Kostas Tavlaridis-Gyparakis
 */
/*--------------------------------------------------------------------------*/
/*---------------------------- IMPLEMENTATION ------------------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include <queue>
#include <functional>
#include <chrono>
#include <cstdlib>
#include <OneVarConstraint.h>

#include "MILPSolver.h"
#include "Block.h"
#include "LinearFunction.h"
#include "DQuadFunction.h"


/*--------------------------------------------------------------------------*/
/*------------------------- NAMESPACE AND USING ----------------------------*/
/*--------------------------------------------------------------------------*/

using namespace SMSpp_di_unipi_it;

// extern int t_pc;


int mycallback(CPXCENVptr env,
               void* cbdata,
               int wherefrom,
               void* cbhandle,
               int* useraction_p) {

 return static_cast<MILPSolver*>(cbhandle)->Callback(env, cbdata, wherefrom, useraction_p);
}

SMSpp_insert_in_factory_cpp_0( MILPSolver );

/*--------------------------------------------------------------------------*/
/*--------------------------------- METHODS --------------------------------*/
/*--------------------------------------------------------------------------*/

void MILPSolver::count_constraints(FRowConstraint& constraint, int& numrows) {
 std::cout << "[DEBUG] ========= MILPSolver::count_constraints()" << std::endl;
 std::cout << "[DEBUG] numrows = " << numrows << std::endl;
 std::cout << "[DEBUG] " << constraint;

 // Counter is incremented only if constraint is described by a linear function
 auto fun = dynamic_cast<const LinearFunction*>(constraint.get_function());
 if (fun != nullptr) {
  ++numrows;
 } else {
  throw (std::invalid_argument("The Constraint is not linear"));
 }
}

/*--------------------------------------------------------------------------*/

void MILPSolver::count_variables(ColVariable& variable, int& numcols, int& nzelements) {
 std::cout << "[DEBUG] ========= MILPSolver::count_variables()" << std::endl;
 std::cout << "[DEBUG] numcols    = " << numrows << std::endl;
 std::cout << "[DEBUG] nzelements = " << numrows << std::endl;
 std::cout << "[DEBUG] " << variable;
 std::cout << "[DEBUG] The active stuff is:" << std::endl;

 ++numcols;


 // TODO resize these somewhere else
 if (active_row_constraints.size() < numcols) {
  // The two vectors are resized together so they should have the same size
  active_row_constraints.resize(static_cast<unsigned long>(numcols));
  active_box_constraints.resize(static_cast<unsigned long>(numcols));
 }

 // Each ColVariable has a vector of active things
 // so we check if said things are row constraints or box constraints
 for (auto i : variable.active_stuff()) {
  auto row = dynamic_cast<FRowConstraint*>(i);
  if (row != nullptr) {
   std::cout << "[DEBUG] " << *row;
   active_row_constraints[numcols - 1].push_back(row);
   ++nzelements;
  }
  auto box = dynamic_cast<OneVarConstraint*>(i);
  if (box != nullptr) {
   std::cout << "[DEBUG] " << *box;
   active_box_constraints[numcols - 1].push_back(box);
  }
  auto obj = dynamic_cast<Objective*>(i);
  if (obj != nullptr) {
   std::cout << "[DEBUG] " << *obj;
  }
 }
}

/*--------------------------------------------------------------------------*/

int MILPSolver::index_of_constraint(FRowConstraint* p_const) {

 auto p_fun = dynamic_cast<const LinearFunction*> (p_const->get_function());
 if (p_fun == nullptr) {
  throw (std::invalid_argument("The Constraint is not linear"));
 }

 int i = 0;

 auto it = lower_bound(v_s_const_int.begin(),
                       v_s_const_int.end(),
                       p_const,
                       [&](const_int pair, FRowConstraint* pconst) {
                        return pair.first < pconst;
                       });

 if (it == v_s_const_int.end()) {
  it = (v_s_const_int.rbegin() + 1).base();
 } else if (it != v_s_const_int.begin() && it->first > p_const) {
  --it;
 }

 if (it != v_s_const_int.end()) {
  i = static_cast<int>(it->second + std::distance(it->first, p_const));
 } else {
  throw (std::invalid_argument("Current Constraint is not defined in the CPLEX Coeff Matrix"));
 }

 return i;
}

/*--------------------------------------------------------------------------*/

FRowConstraint* MILPSolver::constraint_with_index(int i) {

 FRowConstraint* p_const;
 auto it = lower_bound(v_int_s_const.begin(),
                       v_int_s_const.end(),
                       i,
                       [&](int_const pair, int i) {
                        return pair.first < i;
                       });

 if (it == v_int_s_const.end()) {
  it = (v_int_s_const.rbegin() + 1).base();
 } else if (it != v_int_s_const.begin() && it->first > i) {
  --it;
 }

 if (it != v_int_s_const.end()) {
  p_const = it->second;
 } else {
  throw (std::invalid_argument("Current index is not defined in the CPLEX Coeff Matrix"));
 }

 return p_const;
}

/*--------------------------------------------------------------------------*/

int MILPSolver::index_of_variable(ColVariable* p_var) {

 int i;
 auto it = lower_bound(v_s_var_int.begin(),
                       v_s_var_int.end(),
                       p_var,
                       [&](var_int pair, ColVariable* pvar) {
                        return pair.first < pvar;
                       });

 if (it == v_s_var_int.end()) {
  it = (v_s_var_int.rbegin() + 1).base();
 } else if (it != v_s_var_int.begin() && it->first > p_var) {
  --it;
 }

 if (it != v_s_var_int.end()) {
  i = static_cast<int>(it->second + std::distance(it->first, p_var));
 } else {
  throw (std::invalid_argument("Current Variable is not defined in the CPLEX Coeff Matrix"));
 }

 return i;
}

/*--------------------------------------------------------------------------*/

ColVariable* MILPSolver::variable_with_index(int i) {

 ColVariable* p_var;
 auto it = lower_bound(v_int_s_var.begin(),
                       v_int_s_var.end(),
                       i,
                       [&](int_var pair, int i) {
                        return pair.first < i;
                       });

 if (it == v_int_s_var.end()) {
  it = (v_int_s_var.rbegin() + 1).base();
 } else if (it != v_int_s_var.begin() && it->first > i) {
  --it;
 }

 if (it != v_int_s_var.end()) {
  p_var = it->second;
 } else {
  throw (std::invalid_argument("Current index is not defined in the CPLEX Coeff Matrix"));
 }

 return p_var;
}

/*--------------------------------------------------------------------------*/

void MILPSolver::scan_static_constraints(FRowConstraint& lconst,
                                         char* sense,
                                         double* rhs,
                                         int& first,
                                         int& i) {
 // We have to check that FRowConstraint has a Linear Function
 // If not, the problem is not solvable with MILP
 auto lin_fun = dynamic_cast<const LinearFunction*>(lconst.get_function());
 if (lin_fun == nullptr) {
  throw (std::invalid_argument("The Constraint is not linear"));
 }
 std::cout << "[DEBUG] ========= MILPSolver::scan_static_constraints()" << std::endl;
 // Remember: i is the index of the constraint/row from the CPLEX matrix,
 // only FRowConstraints with LinearFunctions are counted.

/*
 * In order to define the type of the operation of the constraint we need to
 * check and compare the lhs and the rhs of the Constraint due to the fact
 * that is built in the following form:
 *
 * LHS <= ( some function from Variables to reals ) <= RHS
 */

 auto lconst_lhs = lconst.get_lhs();
 auto lconst_rhs = lconst.get_rhs();

 if (lconst_lhs == lconst_rhs) {
  // equality
  sense[i] = 'E';
  // setting the rhs
  rhs[i] = lconst_rhs;
 } else if (lconst_lhs == -Inf<double>()) {
  // inequality (smaller/equal)
  sense[i] = 'L';
  //setting the rhs, which is equal to lhs
  rhs[i] = lconst_rhs;
 } else if (lconst_rhs == Inf<double>()) {
  // inequality (greate/equal)
  sense[i] = 'G';
  // setting the rhs, which is equal to rhs
  rhs[i] = lconst_lhs;
 }

 std::cout << "[DEBUG] i          = " << i << std::endl;
 std::cout << "[DEBUG] lconst     = " << lconst;
 std::cout << "[DEBUG] sense[i]   = " << i << std::endl;
 std::cout << "[DEBUG] lconst_lhs = " << lconst_lhs << std::endl;
 std::cout << "[DEBUG] lconst_rhs = " << lconst_rhs << std::endl;


 // Check if we are in dynamic or static part of the problem
 if (first == 0) {
  // emplace_back() spares a copy operation
  // v_s_const_int.push_back( std::make_pair( &lconst, i ));
  // v_int_s_const.push_back( std::make_pair( i, &lconst ));
  v_s_const_int.emplace_back(&lconst, i);
  v_int_s_const.emplace_back(i, &lconst);
 }

 // Increasing the counter, since we are in static part
 first++;
 i++;
}

/*--------------------------------------------------------------------------*/

void MILPSolver::scan_dynamic_constraints(FRowConstraint& lconst,
                                          char* sense,
                                          double* rhs,
                                          int& i) {

 // We have to check that FRowConstraint has a Linear Function
 // If not, the problem is not solvable with MILP
 auto lin_fun = dynamic_cast<const LinearFunction*>(lconst.get_function());
 if (lin_fun == nullptr) {
  throw (std::invalid_argument("The Constraint is not linear"));
 }

 std::cout << "[DEBUG] ========= MILPSolver::scan_dynamic_constraints()" << std::endl;
 // Remember: i is the index of the constraint/row from the CPLEX matrix,
 // only FRowConstraints with LinearFunctions are counted.

/*
 * In order to define the type of the operation of the constraint we need to
 * check and compare the lhs and the rhs of the Constraint due to the fact
 * that is built in the following form:
 *
 * LHS <= ( some function from Variables to reals ) <= RHS
 */

 auto lconst_lhs = lconst.get_lhs();
 auto lconst_rhs = lconst.get_rhs();

 if (lconst_lhs == lconst_rhs) {
  // equality
  sense[i] = 'E';
  // setting the rhs
  rhs[i] = lconst_rhs;
 } else if (lconst_lhs == -Inf<double>()) {
  // inequality (smaller/equal)
  sense[i] = 'L';
  // setting the rhs, which is equal to lhs
  rhs[i] = lconst_rhs;
 } else if (lconst_rhs == Inf<double>()) {
  // inequality (greate/equal)
  sense[i] = 'G';
  //setting the rhs, which is equal to rhs
  rhs[i] = lconst_lhs;
 }

 std::cout << "[DEBUG] i          = " << i << std::endl;
 std::cout << "[DEBUG] lconst     = " << lconst;
 std::cout << "[DEBUG] sense[i]   = " << i << std::endl;
 std::cout << "[DEBUG] lconst_lhs = " << lconst_lhs << std::endl;
 std::cout << "[DEBUG] lconst_rhs = " << lconst_rhs << std::endl;

 // emplace_back() spares a copy operation
 // v_d_const_int.push_back( std::make_pair( &lconst, i ));
 // v_int_d_const.push_back( std::make_pair( i, &lconst ));
 v_d_const_int.emplace_back(&lconst, i);
 v_int_d_const.emplace_back(i, &lconst);
 i++;
}

/*--------------------------------------------------------------------------*/

void MILPSolver::scan_static_variables(ColVariable& var,
                                       int* matbeg,
                                       int* matcnt,
                                       int* matind,
                                       double* matval,
                                       // double *objective,
                                       // double *q_objective,
                                       double* lb,
                                       double* ub,
                                       char* xctype,
                                       int& first,
                                       int& i) {

 std::cout << "[DEBUG] ========= MILPSolver::scan_static_variables()" << std::endl;

 if (first == 0) {
  // emplace_back() spares a copy operation
  // v_s_var_int.push_back( std::make_pair( &var, i ));
  // v_int_s_var.push_back( std::make_pair( i, &var ));
  v_s_var_int.emplace_back(&var, i);
  v_int_s_var.emplace_back(i, &var);
 }

 // Setting the lower and upper bounds of the variable
 // First check the bounds from the ColVariable class
 if (var.get_lb() == -Inf<double>())
  lb[i] = -CPX_INFBOUND;
 else
  lb[i] = var.get_lb();

 if (var.get_ub() == Inf<double>())
  ub[i] = CPX_INFBOUND;
 else
  ub[i] = var.get_ub();
 // Then scan the box constraints
 int bounds = static_cast<int>(active_box_constraints[i].size());
 for (int j = 0; j < bounds; ++j) {
  auto box = active_box_constraints[i][j];
  lb[i] = lb[i] > box->get_lhs()? lb[i] : box->get_lhs();
  ub[i] = ub[i] < box->get_rhs()? ub[i] : box->get_rhs();
 }

 // Setting the type of the variable
 if (var.is_integer()) {
  if (var.is_unitary() && var.is_positive()) {
   xctype[i] = 'B'; // Binary
  } else {
   xctype[i] = 'I'; // Integer
  }
 } else {
  xctype[i] = 'C';  // Continuous
 }

 std::cout << "[DEBUG] i         = " << i << std::endl;
 std::cout << "[DEBUG] var       = " << var;
 std::cout << "[DEBUG] lb[i]     = " << lb[i] << std::endl;
 std::cout << "[DEBUG] ub[i]     = " << ub[i] << std::endl;
 std::cout << "[DEBUG] xctype[i] = " << xctype[i] << std::endl;

 // Remember: i is the index of the variable/column from the CPLEX matrix,
 // all the ColVariables are counted, even if they don't have active FRowConstraints

 // Setting the number of non-zero elements for this Variable in the CPLEX coeff matrix
 int nzelements = static_cast<int>(active_row_constraints[i].size());
 matcnt[i] = nzelements;

 // Setting the index of index of the beginning of column in the CPLEX coefficient matrix
 if (i == 0) {
  matbeg[i] = 0;
 } else {
  matbeg[i] = matbeg[i - 1] + matcnt[i - 1];
 }

 // We need to find where is the index in matval
 // and matind for the corresponding variable
 // FIXME I suspect that indexed[] and matbeg[] do the same thing
 if (i == 0) {
  indexed[i] = 0;
 } else if (i > 0) {
  indexed[i] = indexed[i - 1] + matcnt[i - 1];
 }

 // Setting the coefficients and the indexes of the corresponging rows in the CPLEX coeff matrix
 for (int j = 0; j < nzelements; ++j) {

  // Retreive the function of each active constraint for the variable
  auto p_const = dynamic_cast<FRowConstraint*> ( active_row_constraints[i][j]);
  auto p_fun = dynamic_cast<const LinearFunction*> (p_const->get_function());

  // Retrieve the coefficient of the variable for that constraint
  auto it = std::find_if(p_fun->get_v_var().begin(),
                         p_fun->get_v_var().end(),
                         [&](LinearFunction::coeff_pair pair) {
                          return pair.first == &var;
                         });

  if (it != p_fun->get_v_var().end()) {
   matval[indexed[i] + j] = it->second;
  } else {
   throw (std::invalid_argument("Variable is not active in the examined Constraint"));
  }

  /* Passing the index of the constraint to the corresponding vector matind,
   * note that here we need to check if we are in the dynamic or static part
   * of the problem, since we will have to use the vector of pairs in order to
   * locate the index of the examined constraint in the CPLEX coeff matrix.
   */

  // Locating the index of the constraint in CPLEX
  // TODO check that this still works after changing index_of_constraint()
  matind[indexed[i] + j] = index_of_constraint(p_const);
 } // for j

 first++; //increasing the counter, since we are in static part
 i++;
}

/*--------------------------------------------------------------------------*/

void MILPSolver::scan_dynamic_variables( ColVariable &var,
                                         int *matbeg,
                                         int *matcnt,
                                         int *matind,
                                         double *matval,
                                         double *lb,
                                         double *ub,
                                         char *xctype,
                                         int &i ){

 std::cout << "[DEBUG] ========= MILPSolver::scan_dynamic_variables()" << std::endl;

 // emplace_back() spares a copy operation
 // v_d_var_int.push_back( std::make_pair( &lvar, i ));
 // v_int_d_var.push_back( std::make_pair( i, &lvar ));
 v_d_var_int.emplace_back( &var, i );
 v_int_s_var.emplace_back( i, &var );

 // Setting the lower and upper bounds of the variable
 // First check the bounds from the ColVariable class
 if (var.get_lb() == -Inf<double>())
  lb[i] = -CPX_INFBOUND;
 else
  lb[i] = var.get_lb();

 if (var.get_ub() == Inf<double>())
  ub[i] = CPX_INFBOUND;
 else
  ub[i] = var.get_ub();
 // Then scan the box constraints
 int bounds = static_cast<int>(active_box_constraints[i].size());
 for (int j = 0; j < bounds; ++j) {
  auto box = active_box_constraints[i][j];
  lb[i] = lb[i] > box->get_lhs()? lb[i] : box->get_lhs();
  ub[i] = ub[i] < box->get_rhs()? ub[i] : box->get_rhs();
 }

 // Setting the type of the variable
 if (var.is_integer()) {
  if (var.is_unitary() && var.is_positive()) {
   xctype[i] = 'B'; // Binary
  } else {
   xctype[i] = 'I'; // Integer
  }
 } else {
  xctype[i] = 'C';  // Continuous
 }

 std::cout << "[DEBUG] i         = " << i << std::endl;
 std::cout << "[DEBUG] var       = " << var;
 std::cout << "[DEBUG] lb[i]     = " << lb[i] << std::endl;
 std::cout << "[DEBUG] ub[i]     = " << ub[i] << std::endl;
 std::cout << "[DEBUG] xctype[i] = " << xctype[i] << std::endl;

 // Remember: i is the index of the variable/column from the CPLEX matrix,
 // all the ColVariables are counted, even if they don't have active FRowConstraints

 // Setting the index of variable to the CPLEX coefficient matrix
 matbeg[i] = i;

 // Setting the number of non-zero elements for this Variable in the CPLEX coeff matrix
 int nzelements = static_cast<int>(active_row_constraints[i].size());
 matcnt[i] = nzelements;

 // We need to find where is the index in matval
 // and matind for the corresponding variable
 int index = 0;
 if( i > 0 ) {
  for( int j = 0; j < i; ++j ) {
   index = index + matcnt[j];
  }
 }

 // Setting the coefficients and the indexes of the corresponging rows in the CPLEX coeff matrix
 for (int j = 0; j < nzelements; ++j) {

  // Retreive the function of each active constraint for the variable
  auto p_const = dynamic_cast<FRowConstraint*> ( active_row_constraints[i][j]);
  auto p_fun = dynamic_cast<const LinearFunction*> (p_const->get_function());

  // Retrieve the coefficient of the variable for that constraint
  auto it = std::find_if(p_fun->get_v_var().begin(),
                         p_fun->get_v_var().end(),
                         [&](LinearFunction::coeff_pair pair) {
                          return pair.first == &var;
                         });

  if (it != p_fun->get_v_var().end()) {
   matval[index + j] = (it->second);
  } else {
   throw (std::invalid_argument("Variable is not active in the examined Constraint"));
  }

  /*
   * Passing the index of the constraint to the corresponding vector matind,
   * note that here we need to check if we are in the dynamic or static part
   * of the problem, since we will have to use the vector of pairs in order to
   * locate the index of the examined constraint in the CPLEX coeff matrix.
   */

  // Locating the index of the constraint in CPLEX
  auto it1 = find_if(v_d_const_int.begin(),
                     v_d_const_int.end(),
                     [&](const_int pair) {
                      return pair.first == p_const;
                     });

  if (it1 != v_d_const_int.end()) {
   matind[index + j] = static_cast<int>(it->second);
  } else {
   throw (std::invalid_argument("Current Constraint is not defined in the CPLEX Coeff Matrix"));
  }
 } // for j


 // TODO I do not understand why checking the objective function in a scan_variables()
 // TODO Apparently this doesn't work but the scan_objective() does, so I have to understand why this code was here in the first place

 // Setting the coefficient of the variable for the objective function
 // Checking the type of the objective function
 // auto p_obj = boost::any_cast<const FRealObjective*>(var.get_Block()->get_objective());
 // auto lin_fun = dynamic_cast<const LinearFunction*> (p_obj->get_function());
 //
 // if (lin_fun != nullptr) {
 //  LinearFunction::Index var_index = lin_fun->is_active(&var);
 //  if (var_index != std::numeric_limits<LinearFunction::Index>::infinity()) {
 //   objective[i] = lin_fun->get_coefficient(var_index);
 //  } else {
 //   objective[i] = 0;
 //  }
 // } else {
 //  auto dquad_fun = dynamic_cast<const DQuadFunction*> (p_obj->get_function());
 //  if (dquad_fun != nullptr) {
 //   DQuadFunction::Index var_index = dquad_fun->is_active(&var);
 //   if (var_index != std::numeric_limits<DQuadFunction::Index>::infinity()) {
 //    objective[i] = dquad_fun->get_linear_coefficient(var_index);
 //    q_objective[i] = dquad_fun->get_quadratic_coefficient(var_index);
 //   } else {
 //    objective[i] = 0;
 //    q_objective[i] = 0;
 //   }
 //  } else {
 //   // Throw exception
 //   throw (std::invalid_argument("Unknown type of Objective Function"));
 //  }
 // }
 ++i;
}

void MILPSolver::scan_objective(const FRealObjective *obj, double *objective, double *q_objective ) {
 std::cout << "[DEBUG] ========= MILPSolver::scan_objective()" << std::endl;
 std::cout << "[DEBUG] " << *obj;

 auto lin_fun = dynamic_cast<const LinearFunction *> (obj->get_function());
 int k;

 if(lin_fun != nullptr) {
  for (auto el : lin_fun->get_v_var()) {
   k = index_of_variable(el.first);
   objective[k] = el.second;
  }
 } else {
  auto dquad_fun = dynamic_cast<const DQuadFunction *> (obj->get_function());
  if( dquad_fun != nullptr ) {
   for( auto el : dquad_fun->get_v_var()) {
    // DQuadFunction::get_v_var() returns std::tuples of 3 elements
    k = index_of_variable(std::get<0>(el));
    objective[k] = std::get<1>( el );
    q_objective[k] = std::get<2>( el );
   }
  } else {
   // Throw exception
   throw ( std::invalid_argument( "Unknown type of Objective Function" ));
  }
 }
}

/*--------------------------------------------------------------------------*/

void MILPSolver::set_Block( Block *block ) {

 // If it was attached to some other Block, unregister and initialize
 if (f_Block) {
  f_Block->unregister_Solver(this);
  CPXfreeprob(env, &milp);
  CPXcloseCPLEX(&env);
  milp = nullptr;
  env = nullptr;
 }

 f_Block = block;
 // f_Block->register_Solver( this );

 // Initialisation the CPLEX environment and problem
 int status;
 env = CPXopenCPLEX(&status);
 milp = CPXcreateprob(env, &status, "MILPCPX");

 /*
  * Passing all the data of the Block to the CPLEX Problem
  *
  * Following a Breadth First Search we proceed with scanning the received
  * Block and all of each corresponding children if any, in order to populate
  * the CPLEX data needed to define the corresponding MILP. This is done by
  * the use of the following two CPLEX commands:
  *
  * 1)  int CPXcopylp(CPXCENVptr env, CPXLPptr lp, int numcols, int numrows,
  *     int objsense, double const * objective, double const * rhs, char const *
  *     sense, int const * matbeg, int const * matcnt, int const * matind,
  *     double const * matval, double const * lb, double const * ub, double
  *     const * rngval)
  * 2)  int CPXcopyctype(CPXCENVptr env, CPXLPptr lp, char const * xctype)
  *
  *
  * The following vectors are used in order to retreive from the Block all the
  * information needed in order to construct and pass to CPLEX the corresponding
  * coeff matrix and all the rest of the information needed
  *
  * Where:
  * - objective: refers to a vector with the total coefficients of the objective
  *   function
  * - rhs: refers to a vector with the rhs of all the constraints
  * - sense: refers to a vector with the type of all the constraints
  * - matbeg: refers to a vector with the indices of all the variables
  * - matcnt: refers to a vector with the non-zero elements to correspond to each
  *   variable, keeping track of the indicing from matbeg
  * - matval: refers to a vector with the coefficients that correspond to each
  *   variable, keeping track of the indicing from matbeg
  * - matind: refers to a vector that asociates the coefficients that correspond
  *   to each variable with the corresponding row that they refer to
  * - lb, ub: refers to two array with the lower and upper bounds respectively
  *   of each variable
  * - xctype: refers to a vector that describes the type of each variable
  * - first: is an integer used to denote if the examined constraint/variable is
  *   static or dynamic and in case is static if its the first element of the
  *   examined type of variable.
  * - obj_fun: pointer to the objective function of the Block
  * - numcols, numrows, nzelements: total number of variables, constraints and
  *   non-zero elements of the Block respectively
  */

 numrows = 0;
 numcols = 0;
 nzelements = 0;

 // Queue of blocks to be processed
 std::queue<Block*> Q;

 // LOOP ON THE QUEUE TO COUNT VARIABLES AND CONSTRAINTS
 // Adding the father block to the queue
 Q.push(f_Block);

 while (!Q.empty()) {
  Block* q_Block = Q.front();
  Q.pop();

  // Insert all the children of this block in the queue
  for (auto i : q_Block->get_nested_Blocks()) {
   Q.push(i);
  }

  //We count all the static linear constraints/rows
  for (const auto& i : q_Block->get_static_constraints()) {
   auto f1 = std::bind(&MILPSolver::count_constraints,
                       this,
                       std::placeholders::_1,
                       std::ref(numrows));
   un_any_const_static(i, f1, un_any_type<FRowConstraint>());
  }

  // We count all the dynamic linear constraints/rows
  for (const auto& i : q_Block->get_dynamic_constraints()) {
   auto f1 = std::bind(&MILPSolver::count_constraints,
                       this,
                       std::placeholders::_1,
                       std::ref(numrows));
   un_any_const_dynamic(i, f1, un_any_type<FRowConstraint>());
  }

  // We count all the static Variables/columns
  for (const auto& i : q_Block->get_static_variables()) {
   auto f1 = std::bind(&MILPSolver::count_variables,
                       this,
                       std::placeholders::_1,
                       std::ref(numcols),
                       std::ref(nzelements));
   un_any_const_static(i, f1, un_any_type<ColVariable>());
  }

  // We count all the dynamic Variables/columns
  for (const auto& i : q_Block->get_dynamic_variables()) {
   auto f1 = std::bind(&MILPSolver::count_variables,
                       this,
                       std::placeholders::_1,
                       std::ref(numcols),
                       std::ref(nzelements));
   un_any_const_dynamic(i, f1, un_any_type<ColVariable>());
  }
 } // End of while loop on Block queue

 std::cout << "[DEBUG] ========= MILPSolver::set_Block() after counting loop" << std::endl;
 std::cout << "constraints/numrows = " << numrows << std::endl;
 std::cout << "variables/numcols =   " << numcols << std::endl;
 std::cout << "nzelements =          " << nzelements << std::endl;

 int* matbeg = new int[numcols];
 int* matcnt = new int[numcols];
 int* matind = new int[nzelements];
 auto* matval = new double[nzelements];
 indexed.resize(static_cast<unsigned long>(numcols));
 auto* rhs = new double[numrows];
 char* sense = new char[numrows];
 auto* objective = new double[numcols];
 auto* q_objective = new double[numcols];
 for (int i = 0; i < numcols; i++) {
  objective[i] = 0;
  q_objective[i] = 0;
 }
 auto* lb = new double[numcols];
 auto* ub = new double[numcols];
 char* xctype = new char[numcols];

/*
 * Proceeding with Breadth First Scan/Search of the Block and its children
 */

 // LOOP ON THE QUEUE TO SCAN CONSTRAINTS
 Q.push(f_Block);

 // Counters used to locate the index in the arrays of constraints and
 // variables that need to be populated with the Block data
 int var = 0;
 int col = 0;

 while (!Q.empty()) {
  Block* q_Block = Q.front();
  Q.pop();

  // Insert all the children of this block in the queue
  for (auto i : q_Block->get_nested_Blocks()) {
   Q.push(i);
  }

  /*
   * Scanning and passing all the data of the examined block
   * to the corresponding CPLEX data. We do this by first scanning
   * the static part of the problem and then the dynamic part.
   * This is because we want to have an order in the columns
   * and rows of the CPLEX coeff matrix where the static part is
   * being followed by the dynamic one
   */

  // We scan all the static constraints
  for (const auto& i : q_Block->get_static_constraints()) {
   // Variable used to locate the first element of each type of constraint
   int first = 0;
   auto f1 = std::bind(&MILPSolver::scan_static_constraints,
                       this,
                       std::placeholders::_1,
                       std::ref(sense),
                       std::ref(rhs),
                       std::ref(first),
                       std::ref(col));
   un_any_const_static(i, f1, un_any_type<FRowConstraint>());
  }

  //We scan all the dynamic constraints
  for (const auto& i : q_Block->get_dynamic_constraints()) {
   auto f1 = std::bind(&MILPSolver::scan_dynamic_constraints,
                       this,
                       std::placeholders::_1,
                       std::ref(sense),
                       std::ref(rhs),
                       std::ref(col));
   un_any_const_dynamic(i, f1, un_any_type<FRowConstraint>());
  }
 } // End of while loop on Block queue

 std::cout << "[DEBUG] ========= MILPSolver::set_Block() after constraints scan" << std::endl;
 std::sort(v_s_const_int.begin(), v_s_const_int.end());
 std::sort(v_int_s_const.begin(), v_int_s_const.end());
 std::sort(v_d_const_int.begin(), v_d_const_int.end());
 std::sort(v_int_d_const.begin(), v_int_d_const.end());


 // LOOP ON THE QUEUE TO SCAN VARIABLES
 Q.push(f_Block);
 while (!Q.empty()) {
  Block* q_Block = Q.front();
  Q.pop();

  // Insert all the children of this block in the queue
  for (auto i : q_Block->get_nested_Blocks()) {
   Q.push(i);
  }

  //We scan all the static Variables
  for (const auto& i : q_Block->get_static_variables()) {
   // Variable used to locate the first element of each type of constraint
   int first = 0;
   auto f1 = std::bind(&MILPSolver::scan_static_variables,
                       this,
                       std::placeholders::_1,
                       std::ref(matbeg),
                       std::ref(matcnt),
                       std::ref(matind),
                       std::ref(matval),
                       // std::ref(objective),
                       // std::ref(q_objective),
                       std::ref(lb),
                       std::ref(ub),
                       std::ref(xctype),
                       std::ref(first),
                       std::ref(var));
   un_any_const_static(i, f1, un_any_type<ColVariable>());
  }

  //We scan all the dynamic Variables
  for (const auto& i : q_Block->get_dynamic_variables()) {
   auto f1 = std::bind(&MILPSolver::scan_dynamic_variables,
                       this,
                       std::placeholders::_1,
                       std::ref(matbeg),
                       std::ref(matcnt),
                       std::ref(matind),
                       std::ref(matval),
                       std::ref(lb),
                       std::ref(ub),
                       std::ref(xctype),
                       std::ref(var));
   un_any_const_dynamic(i, f1, un_any_type<ColVariable>());
  }
 } // End of while loop on Block queue

 std::cout << "[DEBUG] ========= MILPSolver::set_Block() after variables scan" << std::endl;
 /* Order all the different vectors of pairs in ascending order based
on the adress of the constraints/variables or the index of CPLEX
rows/columns respectively. */
 std::sort(v_s_var_int.begin(), v_s_var_int.end());
 std::sort(v_int_s_var.begin(), v_int_s_var.end());
 std::sort(v_d_var_int.begin(), v_d_var_int.end());
 std::sort(v_int_d_var.begin(), v_int_d_var.end());

 // Checking the type of the objective function
 try {
  auto tmp_obj = boost::any_cast<FRealObjective*>(f_Block->get_objective());
  switch (tmp_obj->get_sense()) {
   case (Objective::eMax):
    obj_type = -1;
    break;
   case (Objective::eMin):
    obj_type = 1;
    break;
   default:
    obj_type = 0;
    break;
  }
 } catch (boost::bad_any_cast&) {
  std::cout << "Bad any_cast" << std::endl;
  // TODO Do something if objective is not a FRealObjective
 }

 std::cout << "[DEBUG] ========= MILPSolver::set_Block() before objective scan" << std::endl;

 std::cout << "[DEBUG] objective   = [";
 for (int i = 0; i < numcols; ++i) {
  std::cout << " " << objective[i];
 }
 std::cout << "]" << std::endl;

 std::cout << "[DEBUG] q_objective = [";
 for (int i = 0; i < numcols; ++i) {
  std::cout << " " << q_objective[i];
 }
 std::cout << "]" << std::endl;

 // LOOP ON THE QUEUE TO SCAN OBJECTIVES
 Q.push(f_Block);
 int sel = 0;
 const FRealObjective* p_obj = nullptr;
 while (!Q.empty()) {
  Block* q_Block = Q.front();
  Q.pop();

  for (auto i : q_Block->get_nested_Blocks()) {
   Q.push(i);
  }

  // TODO I think that this is already done in scan_dynamic_variables()
  try {
   p_obj = boost::any_cast<FRealObjective*>(q_Block->get_objective());
  } catch (boost::bad_any_cast&) {
   std::cout << "Bad any_cast" << std::endl;
   // TODO Do something if objective is not a FRealObjective
  }

  scan_objective(p_obj, objective, q_objective);
  sel++;
 } // End of while loop on Block queue

/*
 * Due to the fact that CPLEX receives as arguments arrays and our data
 * has been stored with the usage of vectors (due to the fact that they
 * can be dynamically resized) we need to create temporary arrays to which
 * we will pass the vectors by pointers.
 */

 std::cout << "[DEBUG] ========= MILPSolver::set_Block() after objective scan" << std::endl;

 std::cout << "[DEBUG] objective   = [";
 for (int i = 0; i < numcols; ++i) {
  std::cout << " " << objective[i];
 }
 std::cout << "]" << std::endl;

 std::cout << "[DEBUG] q_objective = [";
 for (int i = 0; i < numcols; ++i) {
  std::cout << " " << q_objective[i];
 }
 std::cout << "]" << std::endl;

 std::cout << "[DEBUG] rhs         = [";
 for (int i = 0; i < numrows; ++i) {
  std::cout << " " << rhs[i];
 }
 std::cout << "]" << std::endl;

 std::cout << "[DEBUG] sense       = [";
 for (int i = 0; i < numrows; ++i) {
  std::cout << " " << sense[i];
 }
 std::cout << "]" << std::endl;

 std::cout << "[DEBUG] matbeg      = [";
 for (int i = 0; i < numcols; ++i) {
  std::cout << " " << matbeg[i];
 }
 std::cout << "]" << std::endl;

 std::cout << "[DEBUG] matcnt      = [";
 for (int i = 0; i < numcols; ++i) {
  std::cout << " " << matcnt[i];
 }
 std::cout << "]" << std::endl;

 std::cout << "[DEBUG] matind      = [";
 for (int i = 0; i < nzelements; ++i) {
  std::cout << " " << matind[i];
 }
 std::cout << "]" << std::endl;

 std::cout << "[DEBUG] matval      = [";
 for (int i = 0; i < nzelements; ++i) {
  std::cout << " " << matval[i];
 }
 std::cout << "]" << std::endl;

 std::cout << "[DEBUG] lb          = [";
 for (int i = 0; i < numcols; ++i) {
  std::cout << " " << lb[i];
 }
 std::cout << "]" << std::endl;

 std::cout << "[DEBUG] ub          = [";
 for (int i = 0; i < numcols; ++i) {
  std::cout << " " << ub[i];
 }
 std::cout << "]" << std::endl;

 std::cout << "[DEBUG] xctype      = [";
 for (int i = 0; i < numcols; ++i) {
  std::cout << " " << xctype[i];
 }
 std::cout << "]" << std::endl;

 // Creation of the problem
 CPXcopylp(env, milp, numcols, numrows, obj_type, objective, rhs, sense,
           matbeg, matcnt, matind, matval, lb, ub, nullptr);

 // Copy the objective function coeff matrix
 auto p_dquad_fun = dynamic_cast<const DQuadFunction*> (p_obj->get_function());
 if (p_dquad_fun != nullptr) {
  CPXcopyqpsep(env, milp, q_objective);
 }
 CPXcopyctype(env, milp, xctype);

 // TODO
 // Adding cuts
 // cuts = t_pc;
 // if (cuts == 1) {
 //  // Variables are referred according to original model in callback function
 //  CPXsetintparam(env, CPX_PARAM_MIPCBREDLP, CPX_OFF);
 //  CPXsetintparam(env, CPX_PARAM_PRELINEAR, CPX_OFF);
 //
 //  status = CPXsetusercutcallbackfunc(env, mycallback, this);
 //  status = CPXsetlazyconstraintcallbackfunc(env, mycallback, this);
 // }

 delete[]matbeg;
 delete[]matcnt;
 delete[]matind;
 delete[]matval;

 delete[]rhs;
 delete[]sense;
 delete[]objective;
 delete[]q_objective;
 delete[]lb;
 delete[]ub;
 delete[]xctype;
}


/*--------------------------------------------------------------------------*/

Solver::OFValue MILPSolver::get_lb() {

 OFValue lower_bound = 0;

 switch( obj_type ) {
  case 1: // Minimization problem
   switch( sol_status ) {
    case kUnbounded:  lower_bound = -Inf<OFValue>(); break;
    case kInfeasible: lower_bound = Inf<OFValue>(); break;
    default:          CPXgetbestobjval( env, milp, &lower_bound ); break;
   }
   break;
  case -1: // Maximization problem
   switch( sol_status ) {
    case kUnbounded:  lower_bound = Inf<OFValue>(); break;
    case kInfeasible: lower_bound = -Inf<OFValue>(); break;
    default:          CPXgetbestobjval( env, milp, &lower_bound ); break;
   }
   break;
  default:
   throw ( std::runtime_error( "Objective type not yet defined" ));
   break;
 }

 return lower_bound;
}

 /*--------------------------------------------------------------------------*/

Solver::OFValue MILPSolver::get_ub(){

 OFValue upper_bound = 0;
 switch( obj_type ) {
  case 1: // Minimization problem
  switch( sol_status ) {
   case kUnbounded:  upper_bound = -Inf<OFValue>(); break;
   case kInfeasible: upper_bound = Inf<OFValue>(); break;
   default:          CPXgetobjval( env, milp, &upper_bound ); break;
  }
  break;
  case -1: // Maximization problem
  switch( sol_status ) {
   case kUnbounded:  upper_bound = Inf<OFValue>(); break;
   case kInfeasible: upper_bound = -Inf<OFValue>(); break;
   default:          CPXgetobjval( env, milp, &upper_bound ); break;
  }
  break;
  default:
   throw ( std::runtime_error( "Objective type not yet defined" ));
   break;
 }

 return upper_bound;
}

/*--------------------------------------------------------------------------*/
// TODO all the three set_par() have to be revised
// I'm not sure it is ideal to split them according the type of the value
void MILPSolver::set_par(const int par, const int value) {
 switch (par) {
  case intLogVerb:
   f_log_verb = value;
   CPXsetintparam(env, CPX_PARAM_SCRIND, value);
   break; //1
  case intMaxSol:
   f_max_sol = value;
   CPXsetintparam(env, CPXPARAM_MIP_Limits_Solutions, value);
   break;
  case intMaxIter:
   f_max_iter = value;
   CPXsetlongparam(env, CPXPARAM_MIP_Limits_Nodes, value);
  default:
   CPXsetintparam(env, par - intLastAlgPar, value);
   break;
 }
}

/*--------------------------------------------------------------------------*/

// TODO: Use ThinComputeInterface setters
void MILPSolver::set_par( const int par, const double value ) {
 switch (par) {
  case dblMaxTime:
   f_max_time = value;
   CPXsetdblparam(env, CPXPARAM_TimeLimit, value);
   break; //10
  case dblRelAcc:
   f_rel_acc = value;
   // CPXsetdblparam(env, ,value );
   break;
  case dblAbsAcc:
   f_abs_acc = value;
   // CPXsetdblparam(env, ,value );
   break;
  case dblUpCutOff:
   f_up_cutoff = value;
   CPXsetdblparam(env, CPXPARAM_MIP_Tolerances_UpperCutoff, value);
   break;
  case dblLwCutOff:
   f_lw_cutoff = value;
   CPXsetdblparam(env, CPXPARAM_MIP_Tolerances_LowerCutoff, value);
   break;
  case dblRAccSol:
   f_r_acc_sol = value;
   CPXsetdblparam(env, CPXPARAM_MIP_Pool_RelGap, value);
   break;
  case dblAAccSol:
   f_a_acc_sol = value;
   CPXsetdblparam(env, CPXPARAM_MIP_Pool_AbsGap, value);
   break;
  case dblFAccSol:
   f_f_acc_sol = value;
   // CPXsetdblparam(env, ,value );
   break;
  default:
   CPXsetdblparam(env, par - dblLastAlgPar, value);
   break;
 }
 }

// void MILPSolver::set_par( const int par , const long value )
// {
//  switch( par ) {
//   case( kMaxIter ): f_max_iter = value; CPXsetlongparam(env, CPXPARAM_MIP_Limits_Nodes, value); break; // 0 or max
//   default:                              CPXsetlongparam(env, par-kLastAlgPar,value ); break;
//   }
//  }

/*--------------------------------------------------------------------------*/

int MILPSolver::solve() {

 // Write the problem on file
	// TODO Add this feature as configurable
 // CPXwriteprob(env, milp, "mipex1.lp", nullptr);

 // Solve the problem
 CPXmipopt(env, milp);

 // Get the number of nodes used to solve the problem
 nodes = CPXgetnodecnt(env, milp);

 // Write the solution on the Block
 // TODO pass configuration object
 get_var_solution(nullptr);

 // Get the status
 int status = CPXgetstat(env, milp);

 // Obtaining appropriate sol_status
 switch (status) {
  case 101:
  case 102:
  case 104:
   sol_status = kOK;
   break;
  case 103:
   sol_status = kInfeasible;
   break;
  case 105:
  case 106:
   sol_status = kStopIter;
   break;
  case 107:
  case 108:
   sol_status = kStopTime;
   break;
  case 109:
  case 110:
   sol_status = kError;
   break;
  case 118:
   sol_status = kUnbounded;
   break;
  default:
   sol_status = status;
   break;
 }
 return sol_status;
}

int MILPSolver::Callback(CPXCENVptr env, void* cbdata, int wherefrom, int* useraction_p) {

 // Pass to Block the solution of current node
 double* tmpx = new double[numcols];

 *useraction_p = CPX_CALLBACK_DEFAULT;

 // Obtain solution from CPLEX
 CPXgetcallbacknodex(env, cbdata, wherefrom, tmpx, 0, numcols - 1);

  // Retrieve Solution for all the Variables with respect to the order
 int col = 0;

 std::queue<Block*> Q; // Creating the queue
 Q.push(f_Block);      // Passing the root, i.e. the father block

 while (!Q.empty())    // Iterating for the father block and all children
 {
  Block* q_Block = Q.front(); // Block to be examined
  Q.pop(); //take out from the queue the examined block

  for (auto i : q_Block->get_nested_Blocks()) {
   Q.push(i);
  }

  // Passing values for all static Variables
  for (const auto& i : q_Block->get_static_variables()) {
   auto f1 = std::bind(&MILPSolver::set_var_value,
                       this,
                       std::placeholders::_1,
                       std::ref(tmpx),
                       std::ref(col));
   un_any_const_static(i, f1, un_any_type<ColVariable>());
  }

  // We scan all the dynamic Variables
  for (const auto& i : q_Block->get_dynamic_variables()) {
   auto f1 = std::bind(&MILPSolver::set_var_value,
                       this,
                       std::placeholders::_1,
                       std::ref(tmpx),
                       std::ref(col));
   un_any_const_dynamic(i, f1, un_any_type<ColVariable>());
  }

  // we call this method where is expected to be written the mechanism
  // for constructing new constraints
  q_Block->generate_dynamic_constraints();
 } // while loop

 // Checking if new constraints have been creating
 // by looking at the list of modifications

 // TODO: Port modifications in the new version
 /*
 if (!v_mod.empty()) {

  int mynzcnt = 0;
  double myrhs = 0;
  char sense = ('L');

  std::cout << std::endl << "New Modification received to Solver its size is: " << v_mod.size();
  std::cout << std::endl << "//";
  for (auto it = v_mod.begin(); it != v_mod.end(); it++) {

   shared_ptr <BlockModificationAD> dmod = dynamic_pointer_cast<BlockModificationAD>(*it);

   if (dmod) {

    //if(dmod->f_type == BlockModificationAD::eAddConst){ // we need to add all attached new constraints to CPLEX
    std::cout << std::endl << "Beefore stuff ";
    std::cout << std::endl << "//";
    auto v_cuts = boost::any_cast<std::list<LinearConstraint>*>(dmod->whc_list);
    std::cout << std::endl << "No. of Constraints to be added: " << v_cuts->size();
    for (auto ct = v_cuts->begin(); ct != v_cuts->end(); ct++) {
     std::cout << std::endl << "suppa";
     //pass the new cuts to CPLEX
     mynzcnt = ct->get_num_active_var();
     if (ct->get_lhs() == ct->get_rhs()) { //equality
      sense = ('E');
      myrhs = (ct->get_rhs());//setting the rhs
     } else if (ct->get_lhs() == -Inf<double>()) { // inequality (smaller/equal)
      sense = ('L');
      myrhs = (ct->get_rhs());//setting the rhs, which is equal to lhs
     } else if (ct->get_rhs() == Inf<double>()) {  // inequality (greate/equal)
      sense = ('G');
      myrhs = (ct->get_lhs());//setting the rhs, which is equal to rhs
     }

     int* mycutind = new int[mynzcnt];
     double* mycutval = new double[mynzcnt];

     int i = 0;
     for (auto it = ct->get_v_var()->begin(); it != ct->get_v_var()->end(); ++it) {
      mycutind[i] = index_of_variable(it->first);
      mycutval[i] = it->second;
      i++;
     }

     CPXcutcallbackadd(env, cbdata, wherefrom, mynzcnt, myrhs, sense, mycutind, mycutval, CPX_USECUT_PURGE);

     delete[]mycutind;
     delete[]mycutval;


    } // for (auto ct:v_cuts)

   }// if dmod check
   else {
    std::cout << "PROBLEEEEEMOOOOOOOOOO";
   }

  }//for-loops of modifications


  int j = v_mod.size();
  auto it = v_mod.begin();
  for (int i = 0; i < j && it != v_mod.end(); i++) {
   it = v_mod.erase(it);
  }
  std::cout << std::endl << "check clear = " << v_mod.size();
 }//if existing modifications
 */



 return( 0 );

 }


/*--------------------------------------------------------------------------*/

void MILPSolver::set_var_value(ColVariable& lvar, double* tmpx, int& i) {
 std::cout << "[DEBUG] ========= MILPSolver::set_var_value()" << std::endl;
 std::cout << "[DEBUG] i     = " << i << std::endl;
 std::cout << "[DEBUG] value = " << tmpx[i] << std::endl;
 lvar.set_value(tmpx[i]);
 i++;
}

/*--------------------------------------------------------------------------*/

void MILPSolver::get_var_solution(Configuration* solc) {
 std::cout << "[DEBUG] ========= MILPSolver::get_var_solution()" << std::endl;
 // TODO Support Configuration

 auto* tmpx = new double[numcols];

 // Obtain solution from CPLEX
 CPXgetmipx(env, milp, tmpx, 0, numcols - 1);

 // Retrieve Solution for all the Variables with respect to the order
 int col = 0;

 std::queue<Block*> Q;
 Q.push(f_Block);

 while (!Q.empty()) {
  Block* q_Block = Q.front();
  Q.pop();

  for (auto i : q_Block->get_nested_Blocks()) {
   Q.push(i);
  }

  // Passing values for all static Variables
  for (const auto& i : q_Block->get_static_variables()) {
   auto f1 = std::bind(&MILPSolver::set_var_value,
                       this,
                       std::placeholders::_1,
                       std::ref(tmpx),
                       std::ref(col));
   un_any_const_static(i, f1, un_any_type<ColVariable>());
  }

  // We scan all the dynamic Variables
  for (const auto& i : q_Block->get_dynamic_variables()) {
   auto f1 = std::bind(&MILPSolver::set_var_value,
                       this,
                       std::placeholders::_1,
                       std::ref(tmpx),
                       std::ref(col));
   un_any_const_dynamic(i, f1, un_any_type<ColVariable>());
  }

 } // While queue is not empty

 // The following code evaluates the objective function
 try {
  auto p_obj = boost::any_cast<FRealObjective*>(f_Block->get_objective());
  p_obj->compute();
 } catch (boost::bad_any_cast&) {
  std::cout << "Bad any_cast" << std::endl;
  // TODO Do something if objective is not a FRealObjective
 }

 // TODO The following code was used to write the value of the objective function
 // double objval;
 // CPXgetobjval(env, milp, &objval);
 //
 // try {
 //  auto p_obj = boost::any_cast<FRealObjective*>(f_Block->get_objective());
 //  auto p_lin_fun = dynamic_cast<const LinearFunction*> (p_obj->get_function());
 //
 //  if (p_lin_fun != nullptr) {
 //   // TODO
 //
 //  } else {
 //   auto p_dquad_fun = dynamic_cast<const DQuadFunction*> (p_obj->get_function());
 //   if (p_dquad_fun != nullptr) {
 //    // TODO
 //   } else {
 //    throw (std::invalid_argument("Unknown type of Objective Function"));
 //   }
 //  }
 // } catch (boost::bad_any_cast&) {
 //  std::cout << "Bad any_cast" << std::endl;
 //  // TODO Do something if objective is not a FRealObjective
 // }

// if( f_Block->get_objective_function().type() ==  typeid(LinearObjectiveFunction *) ){
//
// 	LinearObjectiveFunction * obj_fun = boost::any_cast< LinearObjectiveFunction * >(f_Block->get_objective_function());
// 	obj_fun->set_value(objval);
//
// }// if linear objective function
// else if( f_Block->get_objective_function().type() ==  typeid(DQuadObjectiveFunction *) ){
//
// 	DQuadObjectiveFunction * obj_fun = boost::any_cast< DQuadObjectiveFunction * >(f_Block->get_objective_function());
//     obj_fun->set_value(objval);
// }
  delete []tmpx;

}


void MILPSolver::add_modifications(sp_Mod& mod) {

 auto var_mod = dynamic_cast< VariableMod* >(mod.get());
 auto obj_mod = dynamic_cast< ObjectiveMod* >(mod.get());
 auto const_mod = dynamic_cast< ConstraintMod* >(mod.get());
 auto dyn_mod = dynamic_cast< BlockModAD* >(mod.get());

 // Checking if Modification refers to Obj_Function, Constraint or Variable
 if (var_mod) {
  // var_modification(var_mod);
 } else if (obj_mod) {
  // of_modification(obj_mod);
 } else if (const_mod) {
  // const_modification(const_mod);
 } else if (dyn_mod) {
  dynamic_modification(dyn_mod);
 } else {
  throw (std::invalid_argument("Unknown type of Modification"));
 }
}



void MILPSolver::dynamic_modification(BlockModAD* mod) {

 // TODO The casts in this method will probably need a lot of debugging
 switch (mod->f_type) {


  case (BlockModAD::eAddConst): {
   // Add one or more dynamic constraints

   if (mod->mod_list.type() == typeid(std::vector<FRowConstraint*>)) {
    auto v_const = boost::any_cast<std::vector<FRowConstraint*> >(mod->mod_list);
    for (auto& i : v_const) {
     add_dynamic_constraint(i);
    }

   } else if (mod->mod_list.type() == typeid(FRowConstraint*)) {
    auto p_const = boost::any_cast<FRowConstraint*>(mod->mod_list);
    add_dynamic_constraint(p_const);

   } else {
    throw (std::invalid_argument("Received unexpected type of Constraint to be added"));
   }
   break;
  }

/*
  case (BlockModAD::eAddVar): {
   // Add one or more dynamic variables

   if (mod->mod_list.type() == typeid(std::vector<ColVariable*>)) {

    auto v_vars = boost::any_cast<std::vector<ColVariable*> >(mod->mod_list);
    for (auto& v_var : v_vars) {
     add_dynamic_variable(v_var);
    }

   } else if (mod->mod_list.type() == typeid(ColVariable*)) {

    auto p_var = boost::any_cast<ColVariable*>(mod->mod_list);
    add_dynamic_variable(p_var);

   } else {
    throw (std::invalid_argument("Received unexpected type of Variable to be added"));
   }
   break;
  }


  case (BlockModAD::eDelConst): {
   // Remove one or more dynamic constraints

   if (mod->mod_list.type() == typeid(std::list<FRowConstraint*>)) {

    auto l_const = boost::any_cast<std::list<FRowConstraint*> >(mod->mod_list);
    for (auto& it : l_const) {
     remove_dynamic_constraint(it);
    }

   } else if (mod->mod_list.type() == typeid(FRowConstraint*)) {

    auto* p_const = boost::any_cast<FRowConstraint*>(mod->mod_list);
    remove_dynamic_constraint(p_const);

   } else {
    throw (std::invalid_argument("Received unexpected type of Constraint to be removed"));
   }
   break;
  }

  case (BlockModAD::eDelVar): {
   // Remove one or more dynamic variables

   if (mod->mod_list.type() == typeid(std::list<ColVariable*>)) {

    auto l_var = boost::any_cast<std::list<ColVariable*> >(mod->mod_list);
    for (auto& it : l_var) {
     remove_dynamic_variable(it);
    }

   } else if (mod->mod_list.type() == typeid(ColVariable*)) {

    auto* p_const = boost::any_cast<ColVariable*>(mod->mod_list);
    remove_dynamic_variable(p_const);

   } else {
    throw (std::invalid_argument("Received unexpected type of Constraint to be removed"));
   }
   break;
  }
  */

  default:
   break;
 }
}


void MILPSolver::add_dynamic_constraint(FRowConstraint* r_const) {

 // Retrieve the function
 // TODO Handle exceptions and other classes of Functions
 auto p_fun = dynamic_cast<const LinearFunction*>(r_const->get_function());

 // In order to add a new dynamic constraint to CPLEX we use the
 // routine addrows that needs the following information:

 int nz = r_const->get_num_active_var();

 int matbeg[2] = {0, nz};
 int matind[nz];
 double matval[nz];
 double rhs[1];
 char sense[1];

 // Loop through all the active variables for the constraint and
 // retrieve the corresponding CPLEX indexes
 int i = 0;
 for (auto it = r_const->begin(); it != r_const->end(); ++it) {

  // Retrieve the variable
  // TODO Handle exceptions
  auto p_var = dynamic_cast<ColVariable*>(&*it);

  // Search among dynamic variables
  auto it1 = find_if(v_int_d_var.begin(),
                     v_int_d_var.end(),
                     [&](MILPSolver::int_var pair) {
                      return pair.second == p_var;
                     });

  if (it1 != v_int_d_var.end()) {
   // Found
   matind[i] = (it1->first);
  } else {
   // Search among static variables
   matind[i] = index_of_variable(p_var);
  }

  matval[i] = p_fun->get_coefficient(i);
  i++;
 }

 auto const_lhs = r_const->get_lhs();
 auto const_rhs = r_const->get_rhs();
 if (const_lhs == const_lhs) {
  // equality
  sense[0] = ('E');
  rhs[0] = const_rhs;
 } else if (const_lhs == -Inf<double>()) {
  // inequality (smaller/equal)
  sense[0] = ('L');
  rhs[0] = const_rhs;
 } else if (const_rhs == Inf<double>()) {
  // inequality (greate/equal)
  sense[0] = ('G');
  rhs[0] = const_lhs;
 }

 // Update the vectors of pairs
 v_d_const_int.emplace_back(r_const, v_d_const_int.back().second + 1);
 v_int_d_const.emplace_back(v_int_d_const.back().first + 1, r_const);

 // Update CPLEX problem
 CPXaddrows(env, milp, 0, 1, nz, rhs, sense, matbeg, matind, matval, nullptr, nullptr);
}

/*
void MILPSolver::add_dynamic_variable(ColVariable* r_var) {

 // In order to add a new dynamic variable to CPLEX we use the
 // routine addrows that needs the following information:

 int nz = r_var->get_num_active_const();
 int matbeg[2] = {0, nz};
 int matind[nz];
 double matval[nz];

 int i = 0;
 for (auto it = r_var->active_constraints().begin(); it != r_var->active_constraints().end(); ++it) {

  // locating the index in CPLEX searching first in the dynamic
  // and then in the static part
  LinearConstraint* p_const = dynamic_cast<LinearConstraint*>(*it);
  auto it1 = find_if(v_int_d_const.begin(), v_int_d_const.end(), [&](MILPSolver::int_const pair) {
   return pair.second == p_const;
  });
  if (it1 != v_int_d_const.end())
   matind[i] = (it1->first);
  else
   matind[i] = ind_const(p_const);

  //find the coefficient of the constraint that refers to the examined variable
  auto it2 = find_if(p_const->get_v_var()->begin(), p_const->get_v_var()->end(), [&](LinearConstraint::coeff_pair pair) {
   return pair.first == r_var;
  });

  matval[i] = it2->second;
  i++;
 }


//setting the bounds
 double lb[1];
 double ub[1];

 lb[0] = r_var->get_lb();
 ub[0] = r_var->get_ub();

 v_d_var_int.push_back(std::make_pair(r_var, v_d_var_int.back().second + 1));
 v_int_d_var.push_back(std::make_pair(v_int_d_var.back().first + 1, r_var));

 CPXaddcols(env, milp, 1, nz, NULL, matbeg, matind, matval, lb, ub, NULL);

}

void MILPSolver::remove_dynamic_constraint(LinearConstraint * r_const)
{

//all we need is to find the index of the constraint in CPLEX matrix
int i=0;
auto it = find_if (v_int_d_const.begin(),v_int_d_const.end(), [&](MILPSolver::int_const pair) {
                    return pair.second == r_const;
                    });
          if (it != v_int_d_const.end())
                i = it->first;
          else
        throw( std::invalid_argument( "Constraint is not inside the CPLEX matrix" ) );

CPXdelrows (env, milp, i, i+1);

//remove constraint from the corresponding vector of pairs
v_int_d_const.erase(it);
auto it2 = find_if (v_d_const_int.begin(),v_d_const_int.end(), [&](MILPSolver::const_int pair) {
                    return pair.first == r_const;
                    });

 if (it2 != v_d_const_int.end())
                v_d_const_int.erase(it2);
          else
        throw( std::invalid_argument( "Constraint is not inside the CPLEX matrix" ) );

}


void MILPSolver::remove_dynamic_variable(ColVariable * r_var)
{

//all we need is to find the index of the constraint in CPLEX matrix
int i=0;
auto it = find_if (v_int_d_var.begin(),v_int_d_var.end(), [&](MILPSolver::int_var pair) {
                    return pair.second == r_var;
                    });
          if (it != v_int_d_var.end())
                i = it->first;
          else
        throw( std::invalid_argument( "Constraint is not inside the CPLEX matrix" ) );

CPXdelcols (env, milp, i, i+1);

//remove constraint from the corresponding vector of pairs
v_int_d_var.erase(it);
auto it2 = find_if (v_d_var_int.begin(),v_d_var_int.end(), [&](MILPSolver::var_int pair) {
                    return pair.first == r_var;
                    });

 if (it2 != v_d_var_int.end())
                v_d_var_int.erase(it2);
          else
        throw( std::invalid_argument( "Constraint is not inside the CPLEX matrix" ) );

}


void MILPSolver::var_modification( VariableModification* mod ) //TODO
{


        ColVariable* colvar = dynamic_cast< ColVariable* >( mod->f_variable );

		switch( mod->f_type ) {
		case( VariableModification::eFixVar ):   {
			//in case the Variable is set to be fixed we set the
                        //  lower and upper bound equal to the Variable value

			int cnt = 1; //total number of bounds to be changed
			int * indices = new int [cnt]; //array with the CPLEX index of variable
			indices[0] = ind_var(colvar);
			char * lu = new char [cnt] ; //array with the type of bound for each element
                        lu[0] = 'B' ;
			double * bd = new double [cnt]; //array with the value of bounds
			bd[0] = colvar->get_value();
			CPXchgbds (env, milp, cnt, indices, lu, bd);
            delete []indices; delete []lu; delete []bd;
            break;
        }//case
		case( VariableModification::eUnFixVar ):   {
			//in case the Variable is set to be unfixed we set the
                        //lower and upper bound equal to their original value

			int cnt = 2; //total number of bounds to be changed
			int * indices = new int [cnt]; //array with the CPLEX index of variable
			indices[0] = ind_var(colvar);
			indices[1] = ind_var(colvar);
			char * lu = new char [cnt]; //array with the type of bound for each element
                        lu[0] = 'L' ;
                        lu[1] = 'U' ;
			double * bd = new double [cnt]; //array with the value of bounds
			bd[0] = colvar->get_lb();
			bd[1] = colvar->get_ub();
			CPXchgbds (env, milp, cnt, indices, lu, bd);

            delete []indices; delete []lu; delete []bd;
                        break;
        }//case
		case( ColVariableModification::eChgLB ):   {
                        //changing the LB
			int cnt = 1; //total number of bounds to be changed
			int * indices = new int [cnt]; //array with the CPLEX index of variable
			indices[0] = ind_var(colvar);
			char * lu = new char [cnt]; //array with the type of bound for each element
                        lu[0] = 'L' ;
			double * bd = new double [cnt]; //array with the value of bounds
			bd[0] = colvar->get_lb();
			CPXchgbds (env, milp, cnt, indices, lu, bd);

            delete []indices; delete []lu; delete []bd;
                        break;
        }//case
		case( ColVariableModification::eChgUB ):   {
                        //changing the UB

			int cnt = 1; //total number of bounds to be changed
			int * indices = new int [cnt]; //array with the CPLEX index of variable
			indices[0] = ind_var(colvar);
			char * lu = new char [cnt]; //array with the type of bound for each element
                        lu[0] = 'U' ;
			double * bd = new double [cnt]; //array with the value of bounds
			bd[0] = colvar->get_ub();
			CPXchgbds (env, milp, cnt, indices, lu, bd);

            delete []indices; delete []lu; delete []bd;
                        break;
        }//case
		case( ColVariableModification::eChgType ):   {
                        //in case the type of Variable is set to be changed, we
                        //need to upgrade equivalently the lower and upper bound
                        //as well

			int cnt = 1; //total number of variables that change type
			int * indices = new int [cnt]; //array with the CPLEX index of variable
			indices[0] = ind_var(colvar);
   			char * ctype = new char [cnt]; //array with the ne2 type of the variable
                        switch( colvar->get_type() ) {
			  case( ColVariable::integer ):      ctype[0] = ('I'); break;
			  case( ColVariable::binary ):       ctype[0] = ('B'); break;
			  case( ColVariable::continuous ):   ctype[0] = ('C');
			  }
			CPXchgctype (env, milp, cnt, indices, ctype);

			int cnt2 = 2; //total number of bounds to be changed
			int * indices2 = new int [cnt2]; //array with the CPLEX index of variable
			indices2[0] = ind_var(colvar);
			indices2[1] = ind_var(colvar);
			char * lu = new char [cnt]; //array with the type of bound for each element
                        lu[0] = 'L' ;
                        lu[1] = 'U' ;
			double * bd = new double [cnt2]; //array with the value of bounds
			bd[0] = colvar->get_lb();
			bd[1] = colvar->get_ub();
			CPXchgbds (env, milp, cnt2, indices2, lu, bd);

            delete []indices;  delete []indices2; delete []lu; delete []bd; delete []ctype;

        }//case

			}//switch



}


void MILPSolver::of_modification( ObjFunModification* mod ) //TODO
{


		switch( mod->f_type ) {
			case( ObjFunModification::eSetMin ):  {
			//setting the objective function to minimize
			CPXchgobjsen (env, milp, 1);
                        break;
        }//case
		case( ObjFunModification::eSetMax ):   {
			//setting the objective function to maximize
			CPXchgobjsen (env, milp, -1);
                        break;
        }//case
		case( LinearOFModification::eAddVar ):  {
                        //adding linear coefficients
            LinearOFModification* l_mod = dynamic_cast< LinearOFModification* >( mod );
			int cnt = l_mod->v_variables->size();
			int * indices = new int [cnt]; //array with the CPLEX index of variables
			double * values = new double [cnt]; //array with the new values of variables

            int i=0;
			for (auto it = l_mod->v_variables->begin(); it< l_mod->v_variables->end(); ++it){
			indices[i]= index_of_variable (it->first) ;
			values[i]= it->second;
            i++;
			}
			CPXchgobj(env, milp, cnt, indices, values);
            delete []indices;  delete []values;
                        break;
        }//case
		case( LinearOFModification::eRemoveVar ):   {
                        //deleting linear coefficients
            LinearOFModification* l_mod = dynamic_cast< LinearOFModification* >( mod );
			int cnt = l_mod->v_variables->size();
            int * indices = new int [cnt]; //array with the CPLEX index of variables
			double * values = new double [cnt]; //array with the new values of variables

            int i=0;
			for (auto it = l_mod->v_variables->begin(); it< l_mod->v_variables->end(); ++it){
			indices[i]= ind_var (it->first) ;
			values[i]= 0;
            i++;
			}
			CPXchgobj(env, milp, cnt, indices, values);
            delete []indices;  delete []values;
                        break;
        }//case
		case( LinearOFModification::eModifyCoeff ):   {
                       //changin linear coefficients
            LinearOFModification* l_mod = dynamic_cast< LinearOFModification* >( mod );
			int cnt = l_mod->v_variables->size();
			int * indices = new int [cnt]; //array with the CPLEX index of variables
			double * values = new double [cnt]; //array with the new values of variables

            int i=0;
			for (auto it = l_mod->v_variables->begin(); it< l_mod->v_variables->end(); ++it){
			indices[i]= ind_var (it->first) ;
			values[i]= 0;
            i++;
			}
			CPXchgobj(env, milp, cnt, indices, values);
            delete []indices;  delete []values;
			break;
        }//case

		case( DQOFModification::eqAddVar ):   {
                        //adding quad coefficients
            DQOFModification* q_mod = dynamic_cast< DQOFModification* >( mod );
            for (auto it = q_mod->v_variables->begin(); it< q_mod->v_variables->end(); ++it)
            	CPXchgqpcoef(env, milp,ind_var (it->first),ind_var (it->first), it->second ) ;
            break;
        }//case
		case( DQOFModification::eqRemoveVar ): {
                        //deleting quad coefficients
            DQOFModification* q_mod = dynamic_cast< DQOFModification* >( mod );
			for (auto it = q_mod->v_variables->begin(); it< q_mod->v_variables->end(); ++it)
            	CPXchgqpcoef(env, milp,ind_var (it->first),ind_var (it->first), 0 ) ;
            break;
        }//case
		case( DQOFModification::eqModifyCoeff ):  {
                       //changing quad coefficients
            DQOFModification* q_mod = dynamic_cast< DQOFModification* >( mod );
			for (auto it = q_mod->v_variables->begin(); it< q_mod->v_variables->end(); ++it)
            	CPXchgqpcoef(env, milp,ind_var (it->first),ind_var (it->first), it->second ) ;
            break;
		}//case
} //switch


}


void MILPSolver::const_modification( ConstraintModification* mod )
{

		switch( mod->f_type ) {

		case (ConstraintModification::eRelaxConst):{
            LinearConstraint* r_const = dynamic_cast< LinearConstraint* > ( mod->f_constraint );
            // In order to relax the constraint all we do is transform it to an inequality with
            // rhs equal to infinity

            int cnt = 1; //total number of bounds to be changed
			int * indices = new int [cnt];  //array with the CPLEX index of constraint
			double * values = new double [cnt]; //array with the new rhs of constraint
			char * sense = new char [cnt]; //array with the new sense of constraint
            sense[0] = ('G');
            values[0] =  -Inf<double>();
   			indices[0] = ind_const(r_const);
            CPXchgrhs (env, milp, cnt, indices, values);
			CPXchgsense (env, milp, cnt, indices, sense);
			delete []indices;  delete []values; delete []sense;
			break;
            }//case

		case (ConstraintModification::eEnforceConst):{
            LinearConstraint* r_const = dynamic_cast< LinearConstraint* > ( mod->f_constraint );
            // In order to enforce a relaxed constraint all we need to do is reverse the process
            //  of relaxing it, by changing the sense and the rhs back to the original form of the
            //  constraint

             int cnt = 1; //total number of bounds to be changed
			int * indices = new int [cnt];  //array with the CPLEX index of constraint
			double * values = new double [cnt]; //array with the new rhs of constraint
			char * sense = new char [cnt]; //array with the new sense of constraint

    		if( r_const->get_lhs() == r_const->get_rhs() ){ //equality
				sense[0] = ('E');
				values[0] = (r_const->get_rhs());//setting the rhs
			}
			else if(r_const->get_lhs() == -Inf<double>() ){ // inequality (smaller/equal)
				sense[0] = ('L');
				values[0] = (r_const->get_rhs());//setting the rhs, which is equal to lhs
			}
			else if( r_const->get_rhs() == Inf<double>() ){  // inequality (greate/equal)
				sense[0] = ('G');
				values[0] = (r_const->get_lhs());//setting the rhs, which is equal to rhs
			}

			indices[0] = ind_const(r_const);
            CPXchgrhs (env, milp, cnt, indices, values);
			CPXchgsense (env, milp, cnt, indices, sense);
			delete []indices;  delete []values; delete []sense;
			break;

            }//case

		case (RowConstraintModification::eChgLHS): {
            RowConstraintModification* l_mod = dynamic_cast< RowConstraintModification* >( mod );
            LinearConstraint* r_const = dynamic_cast< LinearConstraint* > ( l_mod->f_constraint );
			int cnt = 1; //total number of bounds to be changed
			int * indices = new int [cnt];  //array with the CPLEX index of constraint
			double * values = new double [cnt]; //array with the new rhs of constraint
			char * sense = new char [cnt]; //array with the new sense of constraint

			if( r_const->get_lhs() == r_const->get_rhs() ){ //equality
				sense[0] = ('E');
				values[0] = (r_const->get_rhs());//setting the rhs
			}
			else if(r_const->get_lhs() == -Inf<double>() ){ // inequality (smaller/equal)
				sense[0] = ('L');
				values[0] = (r_const->get_rhs());//setting the rhs, which is equal to lhs
			}
			else if( r_const->get_rhs() == Inf<double>() ){  // inequality (greate/equal)
				sense[0] = ('G');
				values[0] = (r_const->get_lhs());//setting the rhs, which is equal to rhs
			}

			indices[0] = ind_const(r_const);
            CPXchgrhs (env, milp, cnt, indices, values);
			CPXchgsense (env, milp, cnt, indices, sense);
			delete []indices;  delete []values; delete []sense;
    		break;
        } //case
		case (RowConstraintModification::eChgRHS): {

            RowConstraintModification* l_mod = dynamic_cast< RowConstraintModification* >( mod );
            LinearConstraint* r_const = dynamic_cast< LinearConstraint* > ( l_mod->f_constraint );

			int cnt = 1; //total number of bounds to be changed
			int * indices = new int [cnt];  //array with the CPLEX index of constraint
			double * values = new double [cnt]; //array with the new rhs of constraint
			char * sense = new char [cnt]; //array with the new sense of constraint

			if(  r_const->get_lhs() == r_const->get_rhs() ){ //equality
				sense[0] = ('E');
				values[0] = (r_const->get_rhs());//setting the rhs
			}
			else if( r_const->get_lhs() == -Inf<double>() ){ // inequality (smaller/equal)
				sense[0] = ('L');
				values[0] = (r_const->get_rhs());//setting the rhs, which is equal to lhs
			}
			else if( r_const->get_rhs() == Inf<double>() ){  // inequality (greate/equal)
				sense[0] = ('G');
				values[0] = (r_const->get_lhs());//setting the rhs, which is equal to rhs
			}

			indices[0] = ind_const(r_const);
            CPXchgrhs (env, milp, cnt, indices, values);
			CPXchgsense (env, milp, cnt, indices, sense);
            delete []indices;  delete []values; delete []sense;
		    break;
        } //case

		case( LinearConstraintModification::eAddVar ):   {
                        //adding linear coefficients

            LinearConstraintModification* l_mod = dynamic_cast< LinearConstraintModification* >( mod );
            LinearConstraint* r_const = dynamic_cast< LinearConstraint* > ( l_mod->f_constraint );

           for (auto it = l_mod->v_variables->begin(); it< l_mod->v_variables->end(); ++it)
				 CPXchgcoef (env, milp, ind_const (r_const), ind_var (it->first), it->second);
		   break;
        } //case

		case( LinearConstraintModification::eRemoveVar ):   {

            LinearConstraintModification* l_mod = dynamic_cast< LinearConstraintModification* >( mod );
            LinearConstraint* r_const = dynamic_cast< LinearConstraint* > ( l_mod->f_constraint );
                        //deleting linear coefficients
			for (auto it = l_mod->v_variables->begin(); it< l_mod->v_variables->end(); ++it)
				 CPXchgcoef (env, milp, ind_const (r_const), ind_var (it->first), 0);
		   break;
        } //case

		case( LinearConstraintModification::eModifyCoeff ):   {

            LinearConstraintModification* l_mod = dynamic_cast< LinearConstraintModification* >( mod );
            LinearConstraint* r_const = dynamic_cast< LinearConstraint* > ( l_mod->f_constraint );
                       //changin linear coefficients
			for (auto it = l_mod->v_variables->begin(); it< l_mod->v_variables->end(); ++it)
				 CPXchgcoef (env, milp, index_of_constraint (r_const), ind_var (it->first), it->second);
		   break;
        } //case

}//switch


}
*/
/*--------------------------------------------------------------------------*/
/*----------------------- End File MILPSolver.cpp --------------------------*/
/*--------------------------------------------------------------------------*/ 
