/*--------------------------------------------------------------------------*/
/*-------------------------- File MILPSolver.cpp ---------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Implementation of the MILPSolver class.
 *
 * \version 0.20
 *
 * \date 26 - 03 - 2019
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

#define DEBUG_COUT 0 // TODO: Remove this and all the printouts when done

/*--------------------------------------------------------------------------*/
/*------------------------- NAMESPACE AND USING ----------------------------*/
/*--------------------------------------------------------------------------*/

using namespace SMSpp_di_unipi_it;

// extern int t_pc;

// TODO: Implement this functionality
int mycallback(CPXCENVptr env,
               void* cbdata,
               int wherefrom,
               void* cbhandle,
               int* useraction_p) {
 return static_cast<MILPSolver*>(cbhandle)->Callback(env, cbdata, wherefrom, useraction_p);
}

SMSpp_insert_in_factory_cpp_0(MILPSolver);

/*--------------------------------------------------------------------------*/
/*--------------------------------- METHODS --------------------------------*/
/*--------------------------------------------------------------------------*/

void MILPSolver::count_constraints(FRowConstraint& constraint, int& n_rows) {
#if DEBUG_COUT
 std::cout << "[DEBUG] ========= MILPSolver::count_constraints()" << std::endl;
 std::cout << "[DEBUG] " << constraint;
#endif

 // Counter is incremented only if constraint is described by a linear function
 auto fun = dynamic_cast<const LinearFunction*>(constraint.get_function());
 if (fun != nullptr) {
  ++n_rows;
 } else {
  throw (std::invalid_argument("The Constraint is not linear"));
 }
}

/*--------------------------------------------------------------------------*/

void MILPSolver::count_variables(ColVariable& variable, int& n_cols) {
#if DEBUG_COUT
 std::cout << "[DEBUG] ========= MILPSolver::count_variables()" << std::endl;
 std::cout << "[DEBUG] " << variable;
#endif
 ++n_cols;
}

/*--------------------------------------------------------------------------*/

void MILPSolver::count_nzelements(ColVariable& variable, int& nz_elements) {
#if DEBUG_COUT
 std::cout << "[DEBUG] ========= MILPSolver::count_nzelements()" << std::endl;
 std::cout << "[DEBUG] " << variable;
 std::cout << "[DEBUG] The active stuff is:" << std::endl;
#endif

 static int counter = 0;

 if (active_row_constraints.size() < numcols) {
  // The two vectors are resized together so they should have the same size
  active_row_constraints.resize(static_cast<unsigned long>(numcols));
  active_box_constraints.resize(static_cast<unsigned long>(numcols));
 }

 // Each ColVariable has a vector of active things
 // so we check if such things are row constraints or box constraints
 for (auto i : variable.active_stuff()) {
  auto row = dynamic_cast<FRowConstraint*>(i);
  if (row != nullptr) {
#if DEBUG_COUT
   std::cout << "[DEBUG] " << *row;
#endif
   active_row_constraints[counter].push_back(row);
   ++nz_elements;
  }
  auto box = dynamic_cast<OneVarConstraint*>(i);
  if (box != nullptr) {
#if DEBUG_COUT
   std::cout << "[DEBUG] " << *box;
#endif
   active_box_constraints[counter].push_back(box);
  }
#if DEBUG_COUT
  auto obj = dynamic_cast<Objective*>(i);
  if (obj != nullptr) {
   std::cout << "[DEBUG] " << *obj;
  }
#endif
 }
 ++counter;
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

#if DEBUG_COUT
 std::cout << "[DEBUG] ========= MILPSolver::scan_static_constraints()" << std::endl;
#endif
 auto lin_fun = dynamic_cast<const LinearFunction*>(lconst.get_function());
 if (lin_fun == nullptr) {
  throw (std::invalid_argument("The Constraint is not linear"));
 }

/*
 * In order to define the type of the operation of the constraint we need to
 * check and compare the lhs and the rhs of the Constraint due to the fact
 * that is built in the following form:
 *
 * LHS <= ( some function from Variables to reals ) <= RHS
 *
 * Remember: i is the index of the constraint/row from the CPLEX matrix,
 * and only FRowConstraints with LinearFunctions are counted.
 */

 auto lconst_lhs = lconst.get_lhs();
 auto lconst_rhs = lconst.get_rhs();

 if (lconst_lhs == lconst_rhs) {
  sense[i] = 'E'; // equality
  rhs[i] = lconst_rhs;
 } else if (lconst_lhs == -Inf<double>()) {
  sense[i] = 'L'; // smaller/equal
  rhs[i] = lconst_rhs;
 } else if (lconst_rhs == Inf<double>()) {
  sense[i] = 'G'; // greater/equal
  rhs[i] = lconst_lhs;
 }

 if (first == 0) {
  // We are in the static part
  v_s_const_int.emplace_back(&lconst, i);
  v_int_s_const.emplace_back(i, &lconst);
 }

 ++first;
 ++i;
}

/*--------------------------------------------------------------------------*/

void MILPSolver::scan_dynamic_constraints(FRowConstraint& lconst,
                                          char* sense,
                                          double* rhs,
                                          int& i) {

#if DEBUG_COUT
 std::cout << "[DEBUG] ========= MILPSolver::scan_dynamic_constraints()" << std::endl;
#endif
 auto lin_fun = dynamic_cast<const LinearFunction*>(lconst.get_function());
 if (lin_fun == nullptr) {
  throw (std::invalid_argument("The Constraint is not linear"));
 }

/*
 * In order to define the type of the operation of the constraint we need to
 * check and compare the lhs and the rhs of the Constraint due to the fact
 * that is built in the following form:
 *
 * LHS <= ( some function from Variables to reals ) <= RHS
 *
 * Remember: i is the index of the constraint/row from the CPLEX matrix,
 * and only FRowConstraints with LinearFunctions are counted.
 */

 auto lconst_lhs = lconst.get_lhs();
 auto lconst_rhs = lconst.get_rhs();

 if (lconst_lhs == lconst_rhs) {
  sense[i] = 'E'; // equality
  rhs[i] = lconst_rhs;
 } else if (lconst_lhs == -Inf<double>()) {
  sense[i] = 'L'; // smaller/equal
  rhs[i] = lconst_rhs;
 } else if (lconst_rhs == Inf<double>()) {
  sense[i] = 'G'; // greater/equal
  rhs[i] = lconst_lhs;
 }

 v_d_const_int.emplace_back(&lconst, i);
 v_int_d_const.emplace_back(i, &lconst);
 ++i;
}

/*--------------------------------------------------------------------------*/

void MILPSolver::scan_static_variables(ColVariable& var,
                                       int* matbeg,
                                       int* matcnt,
                                       int* matind,
                                       double* matval,
                                       double* lb,
                                       double* ub,
                                       char* xctype,
                                       int& first,
                                       int& i) {

#if DEBUG_COUT
 std::cout << "[DEBUG] ========= MILPSolver::scan_static_variables()" << std::endl;
#endif

 if (first == 0) {
  v_s_var_int.emplace_back(&var, i);
  v_int_s_var.emplace_back(i, &var);
 }

 // Setting the lower and upper bounds of the variable
 // First check the bounds from the ColVariable class
 // FIXME: I am not sure if ColVariable bounds are used here
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
  lb[i] = lb[i] > box->get_lhs() ? lb[i] : box->get_lhs();
  ub[i] = ub[i] < box->get_rhs() ? ub[i] : box->get_rhs();
 }

 if (var.is_integer()) {
  if (var.is_unitary() && var.is_positive()) {
   xctype[i] = 'B'; // Binary
  } else {
   xctype[i] = 'I'; // Integer
  }
 } else {
  xctype[i] = 'C';  // Continuous
 }

 // Remember: i is the index of the variable/column from the CPLEX matrix,
 // all the ColVariables are counted, even if they don't have active FRowConstraints

 // Number of non-zero elements for this Variable in the CPLEX coeff matrix
 int nz_elements = static_cast<int>(active_row_constraints[i].size());
 matcnt[i] = nz_elements;

 // Index of the beginning of column in the CPLEX coefficient matrix
 if (i == 0) {
  matbeg[i] = 0;
 } else {
  matbeg[i] = matbeg[i - 1] + matcnt[i - 1];
 }

 // We need to find where is the index in matval
 // and matind for the corresponding variable
 // FIXME: I suspect that indexed[] and matbeg[] do the same thing
 if (i == 0) {
  indexed[i] = 0;
 } else if (i > 0) {
  indexed[i] = indexed[i - 1] + matcnt[i - 1];
 }

 // Setting the coefficients and the indexes of the corresponging rows in the CPLEX coeff matrix
 for (int j = 0; j < nz_elements; ++j) {

  // Retrieve the function of each active constraint for the variable
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
  matind[indexed[i] + j] = index_of_constraint(p_const);
 }

 ++first;
 ++i;
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

#if DEBUG_COUT
 std::cout << "[DEBUG] ========= MILPSolver::scan_dynamic_variables()" << std::endl;
#endif

 v_d_var_int.emplace_back( &var, i );
 v_int_s_var.emplace_back( i, &var );

 // Setting the lower and upper bounds of the variable
 // First check the bounds from the ColVariable class
 // FIXME: I am not sure if ColVariable bounds are used here
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

 if (var.is_integer()) {
  if (var.is_unitary() && var.is_positive()) {
   xctype[i] = 'B'; // Binary
  } else {
   xctype[i] = 'I'; // Integer
  }
 } else {
  xctype[i] = 'C';  // Continuous
 }


 // Remember: i is the index of the variable/column from the CPLEX matrix,
 // all the ColVariables are counted, even if they don't have active FRowConstraints

 // Index of variable to the CPLEX coefficient matrix
 matbeg[i] = i;

 // Nnumber of non-zero elements for this Variable in the CPLEX coeff matrix
 int nz_elements = static_cast<int>(active_row_constraints[i].size());
 matcnt[i] = nz_elements;

 // We need to find where is the index in matval
 // and matind for the corresponding variable
 int index = 0;
 if( i > 0 ) {
  for( int j = 0; j < i; ++j ) {
   index = index + matcnt[j];
  }
 }

 // Setting the coefficients and the indexes of the corresponging rows in the CPLEX coeff matrix
 for (int j = 0; j < nz_elements; ++j) {

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

 // TODO: Apparently this doesn't work but the scan_objective() (see) does, so I have to figure out why this code existed in the old version

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
 //   throw (std::invalid_argument("Unknown type of Objective Function"));
 //  }
 // }
 ++i;
}

void MILPSolver::scan_objective(const FRealObjective *obj, double *objective, double *q_objective ) {
#if DEBUG_COUT
 std::cout << "[DEBUG] ========= MILPSolver::scan_objective()" << std::endl;
 std::cout << "[DEBUG] " << *obj;
#endif

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
    objective[k] = std::get<1>(el);
    q_objective[k] = std::get<2>(el);
   }
  } else {
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
                       std::ref(numcols));
   un_any_const_static(i, f1, un_any_type<ColVariable>());
  }

  // We count all the dynamic Variables/columns
  for (const auto& i : q_Block->get_dynamic_variables()) {
   auto f1 = std::bind(&MILPSolver::count_variables,
                       this,
                       std::placeholders::_1,
                       std::ref(numcols));
   un_any_const_dynamic(i, f1, un_any_type<ColVariable>());
  }

  // We count all the static active contraints/non-zero elements
  for (const auto& i : q_Block->get_static_variables()) {
   auto f1 = std::bind(&MILPSolver::count_nzelements,
                       this,
                       std::placeholders::_1,
                       std::ref(nzelements));
   un_any_const_static(i, f1, un_any_type<ColVariable>());
  }

  // We count all the dynamic active contraints/non-zero elements
  for (const auto& i : q_Block->get_dynamic_variables()) {
   auto f1 = std::bind(&MILPSolver::count_nzelements,
                       this,
                       std::placeholders::_1,
                       std::ref(nzelements));
   un_any_const_dynamic(i, f1, un_any_type<ColVariable>());
  }
 } // End of while loop on Block queue

#if DEBUG_COUT
 std::cout << "[DEBUG] ========= MILPSolver::set_Block() after counting loop" << std::endl;
 std::cout << "constraints/numrows = " << numrows << std::endl;
 std::cout << "variables/numcols =   " << numcols << std::endl;
 std::cout << "nzelements =          " << nzelements << std::endl;
#endif

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

 std::sort(v_s_var_int.begin(), v_s_var_int.end());
 std::sort(v_int_s_var.begin(), v_int_s_var.end());
 std::sort(v_d_var_int.begin(), v_d_var_int.end());
 std::sort(v_int_d_var.begin(), v_int_d_var.end());

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
  throw (std::invalid_argument("Objective is not a FRealObjective"));
 }

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

  // TODO: I think that this is already done in scan_dynamic_variables()
  try {
   p_obj = boost::any_cast<FRealObjective*>(q_Block->get_objective());
  } catch (boost::bad_any_cast&) {
   throw (std::invalid_argument("Objective is not a FRealObjective"));
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

#if DEBUG_COUT
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
#endif

 CPXcopylp(env, milp, numcols, numrows, obj_type, objective, rhs, sense,
           matbeg, matcnt, matind, matval, lb, ub, nullptr);

 // Copy the objective function coeff matrix
 auto p_dquad_fun = dynamic_cast<const DQuadFunction*> (p_obj->get_function());
 if (p_dquad_fun != nullptr) {
  CPXcopyqpsep(env, milp, q_objective);
 }
 CPXcopyctype(env, milp, xctype);

 // TODO: Implement cuts
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

// TODO: All the three set_par() have to be revised
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

 process_modifications();

	// TODO: Add this feature as configurable
 CPXwriteprob(env, milp, "problem.lp", nullptr);

 CPXmipopt(env, milp);

 nodes = CPXgetnodecnt(env, milp);

 // TODO: Pass configuration object
 get_var_solution(nullptr);

 int status = CPXgetstat(env, milp);

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
#if DEBUG_COUT
 std::cout << "[DEBUG] ========= MILPSolver::set_var_value()" << std::endl;
 std::cout << "[DEBUG] i     = " << i << std::endl;
 std::cout << "[DEBUG] value = " << tmpx[i] << std::endl;
#endif
 lvar.set_value(tmpx[i]);
 i++;
}

/*--------------------------------------------------------------------------*/

void MILPSolver::get_var_solution(Configuration* solc) {
 // TODO: Support Configuration
#if DEBUG_COUT
 std::cout << "[DEBUG] ========= MILPSolver::get_var_solution()" << std::endl;
#endif

 auto* tmpx = new double[numcols];

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

 }

 // The following code evaluates the objective function
 try {
  auto p_obj = boost::any_cast<FRealObjective*>(f_Block->get_objective());
  p_obj->compute();
 } catch (boost::bad_any_cast&) {
  throw (std::invalid_argument("Objective is not a FRealObjective"));
 }

 // FIXME: The following code was used to write the value of the objective function, figure out if we can toss it
 // double objval;
 // CPXgetobjval(env, milp, &objval);
 //
 // try {
 //  auto p_obj = boost::any_cast<FRealObjective*>(f_Block->get_objective());
 //  auto p_lin_fun = dynamic_cast<const LinearFunction*> (p_obj->get_function());
 //
 //  if (p_lin_fun != nullptr) {
 //   //
 //
 //  } else {
 //   auto p_dquad_fun = dynamic_cast<const DQuadFunction*> (p_obj->get_function());
 //   if (p_dquad_fun != nullptr) {
 //    //
 //   } else {
 //    throw (std::invalid_argument("Unknown type of Objective Function"));
 //   }
 //  }
 // } catch (boost::bad_any_cast&) {
 //  std::cout << "Bad any_cast" << std::endl;
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

/*--------------------------------------------------------------------------*/
void MILPSolver::process_modifications() {
 /*
  * TODO: Write a better demux function
  * This function processes one modification after another, without
  * any attempt of optimization, moreover you have to be CAREFUL to write
  * all the cases in order from the most specialized to the more generic,
  * e.g., OneVarConstraintMod before RowConstraintMod before ConstraintMod,
  * otherwise the generic case will intercept the more specialized Mods.
  */

 while (!v_mod.empty()) {
  auto mod = v_mod.front();

  // A function like this is needed to be called recursively with GroupModifications
  std::function<void(sp_Mod)> f;
  f = [this, &f](sp_Mod mod) {
   {
    const auto tmod = std::dynamic_pointer_cast<GroupModification>(mod);
    if (tmod) {
     for (const auto& submod : tmod->v_sub_Modifications) {
      f(submod);
     }
     return;
    }
   }
   {
    const auto tmod = std::dynamic_pointer_cast<VariableMod>(mod);
    if (tmod) {
     var_modification(tmod.get());
     return;
    }
   }
   {
    const auto tmod = std::dynamic_pointer_cast<ObjectiveMod>(mod);
    if (tmod) {
     of_modification(tmod.get());
     return;
    }
   }
   {
    const auto tmod = std::dynamic_pointer_cast<OneVarConstraintMod>(mod);
    if (tmod) {
     bound_modification(tmod.get());
     return;
    }
   }
   {
    const auto tmod = std::dynamic_pointer_cast<RowConstraintMod>(mod);
    if (tmod) {
     const_modification(tmod.get());
     return;
    }
   }
   {
    const auto tmod = std::dynamic_pointer_cast<ConstraintMod>(mod);
    if (tmod) {
     const_modification(tmod.get());
     return;
    }
   }
   {
    const auto tmod = std::dynamic_pointer_cast<FunctionMod>(mod);
    if (tmod) {
     function_modification(tmod.get());
     return;
    }
   }
   {
    const auto tmod = std::dynamic_pointer_cast<BlockModAD>(mod);
    if (tmod) {
     dynamic_modification(tmod.get());
     return;
    }
   }
  };

  f(mod);
  v_mod.pop_front();
 }
}

/*--------------------------------------------------------------------------*/
void MILPSolver::var_modification(VariableMod* mod) {

 /*
  * VariableMod class does not include any modification types, so we "refresh"
  * all the CPLEX information for the variable. In particular we:
  *
  *  - reset the ctype (binary, integer, continuous)
  *  - if the variable is fixed, we fix it in CPLEX by doing LHS = RHS
  *  - if the variable is not fixed, we reset LHS and RHS values
  */

 auto* var = dynamic_cast<ColVariable*>(mod->f_variable);

 int indices[2];
 indices[0] = index_of_variable(var);
 indices[1] = indices[0];

 char ctype[1];

 if (var->is_integer()) {
  if (var->is_unitary() && var->is_positive()) {
   ctype[0] = 'B'; // Binary
  } else {
   ctype[0] = 'I'; // Integer
  }
 } else {
  ctype[0] = 'C';  // Continuous
 }
 CPXchgctype(env, milp, 1, indices, ctype);

 char* lu;
 double* bd;

 if (var->is_fixed()) {
  lu = new char[1];
  bd = new double[1];
  lu[0] = 'B';
  bd[0] = var->get_value();
  CPXchgbds(env, milp, 1, indices, lu, bd);
 } else {
  lu = new char[2];
  bd = new double[2];

  int bounds = static_cast<int>(active_box_constraints[indices[0]].size());
  lu[0] = 'L';
  lu[1] = 'U';
  bd[0] = -CPX_INFBOUND;
  bd[1] = CPX_INFBOUND;
  for (int i = 0; i < bounds; ++i) {
   auto box = active_box_constraints[indices[0]][i];
   bd[0] = bd[0] > box->get_lhs() ? bd[0] : box->get_lhs();
   bd[1] = bd[1] < box->get_rhs() ? bd[1] : box->get_rhs();
  }
  CPXchgbds(env, milp, 2, indices, lu, bd);
 }

 delete[] lu;
 delete[] bd;
}

/*--------------------------------------------------------------------------*/
void MILPSolver::of_modification(ObjectiveMod* mod) {

 /*
  * ObjectiveMod class does not include any modification types except
  * for eSetMin and eSetMax.
  * To change OF coefficents, a FunctionMod must be used.
  */

 switch (mod->f_type) {

  case (ObjectiveMod::eSetMin): {
   CPXchgobjsen(env, milp, 1);
   break;
  }

  case (ObjectiveMod::eSetMax): {
   //setting the objective function to maximize
   CPXchgobjsen(env, milp, -1);
   break;
  }

  default: {
   throw (std::invalid_argument("Invalid type of ObjectiveMod"));
   break;
  }
 }
}

/*--------------------------------------------------------------------------*/
void MILPSolver::const_modification(ConstraintMod* mod) {

 /*
  * To change the coefficents, a FunctionMod must be used.
  */

 auto p_const = dynamic_cast<FRowConstraint*>(mod->f_constraint);
 auto lf = dynamic_cast<const LinearFunction*>(p_const->get_function());

 int num_bounds;
 int* indices;
 double* values;
 char* sense;

 switch (mod->f_type) {

  case (ConstraintMod::eRelaxConst): {
   // In order to relax the constraint all we do is transform it
   // into an inequality with RHS equal to infinity

   num_bounds = 1;
   indices = new int[num_bounds];
   values = new double[num_bounds];
   sense = new char[num_bounds];

   sense[0] = ('G');
   values[0] = -Inf<double>();
   indices[0] = index_of_constraint(p_const);
   CPXchgrhs(env, milp, num_bounds, indices, values);
   CPXchgsense(env, milp, num_bounds, indices, sense);

   break;
  }

  case (ConstraintMod::eEnforceConst): {
   // In order to enforce a relaxed constraint all we need to do is
   // reverse the process of relaxing it, by changing the sense and
   // the rhs back to the original form of the constraint

   num_bounds = 1;
   indices = new int[num_bounds];
   values = new double[num_bounds];
   sense = new char[num_bounds];

   auto lhs = p_const->get_lhs();
   auto rhs = p_const->get_rhs();

   if (lhs == rhs) {
    sense[0] = 'E';
    values[0] = rhs;
   } else if (lhs == -Inf<double>()) {
    sense[0] = 'L';
    values[0] = rhs;
   } else if (rhs == Inf<double>()) {
    sense[0] = 'G';
    values[0] = lhs;
   }

   indices[0] = index_of_constraint(p_const);
   CPXchgrhs(env, milp, num_bounds, indices, values);
   CPXchgsense(env, milp, num_bounds, indices, sense);

   break;

  }

  default: {
   // TODO: Specifically handle RHS/LHS/BTS modification

   num_bounds = 1;
   indices = new int[num_bounds];
   values = new double[num_bounds];
   sense = new char[num_bounds];

   auto lhs = p_const->get_lhs();
   auto rhs = p_const->get_rhs();

   if (lhs == rhs) {
    sense[0] = 'E';
    values[0] = rhs;
   } else if (lhs == -Inf<double>()) {
    sense[0] = 'L';
    values[0] = rhs;
   } else if (rhs == Inf<double>()) {
    sense[0] = 'G';
    values[0] = lhs;
   }

   indices[0] = index_of_constraint(p_const);
   CPXchgrhs(env, milp, num_bounds, indices, values);
   CPXchgsense(env, milp, num_bounds, indices, sense);
   break;
  }
 }

 delete[]indices;
 delete[]values;
 delete[]sense;
}

/*--------------------------------------------------------------------------*/
void MILPSolver::bound_modification(OneVarConstraintMod* mod) {

 /*
  * The same ColVariable can have more active OneVarConstraints,
  * so each time we modify one of them we have to check if LHS and RHS
  * of the Variable change.
  */

 auto p_const = dynamic_cast<OneVarConstraint*>(mod->f_constraint);
 auto p_var = dynamic_cast<ColVariable*>(p_const->get_active_var(0));

 int indices[2];
 indices[0] = index_of_variable(p_var);
 indices[1] = indices[0];
 int num_bounds = static_cast<int>(active_box_constraints[indices[0]].size());

 char* lu = nullptr;
 double* bd = nullptr;

 switch (mod->f_type) {

  case (RowConstraintMod::eChgLHS):
   lu = new char[1];
   bd = new double[1];

   lu[0] = 'L';
   bd[0] = -CPX_INFBOUND;

   for (int i = 0; i < num_bounds; ++i) {
    auto box = active_box_constraints[indices[0]][i];
    bd[0] = bd[0] > box->get_lhs() ? bd[0] : box->get_lhs();
   }

   CPXchgbds(env, milp, 1, indices, lu, bd);
   break;

  case (RowConstraintMod::eChgRHS):
   lu = new char[1];
   bd = new double[1];

   lu[0] = 'U';
   bd[0] = CPX_INFBOUND;

   for (int i = 0; i < num_bounds; ++i) {
    auto box = active_box_constraints[indices[0]][i];
    bd[0] = bd[0] < box->get_rhs() ? bd[0] : box->get_rhs();
   }

   CPXchgbds(env, milp, 1, indices, lu, bd);
   break;

  case (RowConstraintMod::eChgBTS):
   lu = new char[2];
   bd = new double[2];
   lu[0] = 'L';
   lu[1] = 'U';
   bd[0] = -CPX_INFBOUND;
   bd[1] = CPX_INFBOUND;

   for (int i = 0; i < num_bounds; ++i) {
    auto box = active_box_constraints[indices[0]][i];
    bd[0] = bd[0] > box->get_lhs() ? bd[0] : box->get_lhs();
    bd[1] = bd[1] < box->get_rhs() ? bd[1] : box->get_rhs();
   }

   CPXchgbds(env, milp, 2, indices, lu, bd);
   break;

  default:
   break;
 }
 delete[] lu;
 delete[] bd;
}

/*--------------------------------------------------------------------------*/
void MILPSolver::function_modification(FunctionMod* mod) {

 /*
  * This function is used when changing coefficents for OFs or constraints.
  */
 auto mod_f = mod->f_function;

 // TODO: Handle exceptions here
 auto p_obj = boost::any_cast<FRealObjective*>(f_Block->get_objective());
 auto of = p_obj->get_function();

 if (of == mod_f) {
  auto lf = dynamic_cast<const LinearFunction*> (mod_f);
  auto qf = dynamic_cast<const DQuadFunction*> (mod_f);

  int num_vars;
  int* indices;
  double* values;

  if (lf != nullptr) {
   num_vars = static_cast<int>(lf->get_v_var().size());
   indices = new int[num_vars];
   values = new double[num_vars];

   int i = 0;
   for (auto el : lf->get_v_var()) {
    indices[i] = index_of_variable(el.first);
    values[i] = el.second;
    ++i;
   }
   CPXchgobj(env, milp, num_vars, indices, values);

  } else if (qf != nullptr) {
   num_vars = static_cast<int>(qf->get_v_var().size());
   indices = new int[num_vars];
   values = new double[num_vars];

   int i = 0;
   for (auto el : qf->get_v_var()) {
    // Linear coefficients can be changed all at once with CPXchgobj
    indices[i] = index_of_variable(std::get<0>(el));
    values[i] = std::get<1>(el);

    // Quadratic coefficients can be changed one at a time
    CPXchgqpcoef(env, milp, indices[i], indices[i], std::get<2>(el));
    ++i;
   }

   CPXchgobj(env, milp, num_vars, indices, values);

  } else {
   throw (std::invalid_argument("Unknown type of Objective Function"));
  }

  delete[] indices;
  delete[] values;

 } else {
  // Assume Function is a constraint, so it can be only linear
  auto lf = dynamic_cast<const LinearFunction*> (mod_f);

  if (lf != nullptr) {
   auto p_const = (FRowConstraint*)lf->get_Observer();

   int indices[1];
   indices[0] = index_of_constraint(p_const);

   for (auto el : lf->get_v_var()) {
    CPXchgcoef(env, milp, indices[0], index_of_variable(el.first), el.second);
   }
  } else {
   throw (std::invalid_argument("Unknown type of Function"));
  }
 }
}

/*--------------------------------------------------------------------------*/
void MILPSolver::dynamic_modification(BlockModAD* mod) {

 // TODO: The casts in this method will probably need a lot of debugging
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

  default:
   break;
 }
}

/*--------------------------------------------------------------------------*/
void MILPSolver::add_dynamic_constraint(FRowConstraint* r_const) {
 // TODO: We have to do the same thing with OneVarConstraints

 // Retrieve the function
 // TODO: Handle exceptions and other classes of Functions
 auto p_fun = dynamic_cast<const LinearFunction*>(r_const->get_function());
 if (p_fun == nullptr) {
  throw (std::invalid_argument("The Constraint is not linear"));
 }

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
  // TODO: Handle exceptions
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
 std::sort(v_d_const_int.begin(), v_d_const_int.end());
 std::sort(v_int_d_const.begin(), v_int_d_const.end());

 // We also need to add the constraint to active_row_constraints
 active_row_constraints[i].push_back(r_const);

 // Update CPLEX problem
 CPXaddrows(env, milp, 0, 1, nz, rhs, sense, matbeg, matind, matval, nullptr, nullptr);
}

/*--------------------------------------------------------------------------*/
void MILPSolver::add_dynamic_variable(ColVariable* r_var) {

 std::vector<FRowConstraint *>   var_row_constraints;
 std::vector<OneVarConstraint *> var_box_constraints;
 int nz = 0;

 // We count the nz elements,
 // and build the vectors of active constraints and bounds
 for (auto i : r_var->active_stuff()) {
  auto row = dynamic_cast<FRowConstraint*>(i);
  if (row != nullptr) {
   var_row_constraints.push_back(row);
   ++nz;
  }
  auto box = dynamic_cast<OneVarConstraint*>(i);
  if (box != nullptr) {
   var_box_constraints.push_back(box);
  }
 }

 // In order to add a new dynamic variable to CPLEX we use the
 // routine addrows that needs the following information

 int matbeg[2] = {0, nz};
 int matind[nz];
 double matval[nz];

 int i = 0;
 for (auto* p_const : var_row_constraints) {

  // Locate the index in CPLEX,
  // searching first in the dynamic and then in the static part
  auto p_fun = dynamic_cast<const LinearFunction*>(p_const->get_function());
  if (p_fun == nullptr) {
   throw (std::invalid_argument("The Constraint is not linear"));
  }

  auto it1 = find_if(v_int_d_const.begin(),
                     v_int_d_const.end(),
                     [&](MILPSolver::int_const pair) {
                      return pair.second == p_const;
                     });
  if (it1 != v_int_d_const.end()) {
   matind[i] = (it1->first);
  } else {
   matind[i] = index_of_constraint(p_const);
  }

  matval[i] = p_fun->get_coefficient(i);
  ++i;
 }


 // Setting the lower and upper bounds of the variable
 double lb[1];
 double ub[1];

 // First check the bounds from the ColVariable class
 if (r_var->get_lb() == -Inf<double>()) {
  lb[0] = -CPX_INFBOUND;
 } else {
  lb[0] = r_var->get_lb();
 }

 if (r_var->get_ub() == Inf<double>()) {
  ub[0] = CPX_INFBOUND;
 } else {
  ub[0] = r_var->get_ub();
 }

 // Then scan the box constraints
 int bounds = static_cast<int>(var_box_constraints.size());
 for (int j = 0; j < bounds; ++j) {
  auto box = var_box_constraints[j];
  lb[0] = lb[0] > box->get_lhs() ? lb[0] : box->get_lhs();
  ub[0] = ub[0] < box->get_rhs() ? ub[0] : box->get_rhs();
 }

 // Update all the vectors
 active_row_constraints.emplace_back(var_row_constraints);
 active_box_constraints.emplace_back(var_box_constraints);
 v_d_var_int.emplace_back(r_var, v_d_var_int.back().second + 1);
 v_int_d_var.emplace_back(v_int_d_var.back().first + 1, r_var);

 // Update the CPLEX problem
 CPXaddcols(env, milp, 1, nz, nullptr, matbeg, matind, matval, lb, ub, nullptr);
}

/*--------------------------------------------------------------------------*/
void MILPSolver::remove_dynamic_constraint(FRowConstraint* r_const) {

// Find the index of the constraint in CPLEX matrix
 int i = 0;
 auto it1 = find_if(v_int_d_const.begin(),
                   v_int_d_const.end(),
                   [&](MILPSolver::int_const pair) {
                    return pair.second == r_const;
                   });
 if (it1 != v_int_d_const.end()) {
  i = it1->first;
 } else {
  throw (std::invalid_argument("Constraint is not inside the CPLEX matrix"));
 }

 // Now i is the index of the constraint in the CPLEX matrix
 CPXdelrows(env, milp, i, i + 1);

// Remove constraint from the vectors of pairs
 v_int_d_const.erase(it1);
 auto it2 = find_if(v_d_const_int.begin(),
                    v_d_const_int.end(),
                    [&](MILPSolver::const_int pair) {
                     return pair.first == r_const;
                    });
 if (it2 != v_d_const_int.end()) {
  v_d_const_int.erase(it2);
 } else {
  throw (std::invalid_argument("Constraint is not inside the CPLEX matrix"));
 }

 // We have to update the indexes greater than i in the vectors of pairs
 // TODO: This is quite ugly but it has to be done, think of something smarter
 for (auto it: v_d_const_int) {
  if (it.second > i) {
   it.second--;
  }
 }
 for (auto it: v_int_d_const) {
  if (it.first > i) {
   it.first--;
  }
 }
 for (auto it: v_s_const_int) {
  if (it.second > i) {
   it.second--;
  }
 }
 for (auto it: v_int_s_const) {
  if (it.first > i) {
   it.first--;
  }
 }

 // We have to remove the constraint from the active_row_contraints[]
 // TODO: Look if we can avoiding use active_row_contraints[] altogether
 for (auto& constraints: active_row_constraints) {
  auto constraint = find(constraints.begin(), constraints.end(), r_const);
  if (constraint != constraints.end()) {
   constraints.erase(constraint);
  }
 }
}

/*--------------------------------------------------------------------------*/
void MILPSolver::remove_dynamic_variable(ColVariable* r_var) {

// Find the index of the variable in CPLEX matrix
 int i = 0;
 auto it = find_if(v_int_d_var.begin(),
                   v_int_d_var.end(),
                   [&](MILPSolver::int_var pair) {
                    return pair.second == r_var;
                   });
 if (it != v_int_d_var.end()) {
  i = it->first;
 } else {
  throw (std::invalid_argument("Constraint is not inside the CPLEX matrix"));
 }

 CPXdelcols(env, milp, i, i + 1);

// Remove variable from the vectors of pairs
 v_int_d_var.erase(it);
 auto it2 = find_if(v_d_var_int.begin(),
                    v_d_var_int.end(),
                    [&](MILPSolver::var_int pair) {
                     return pair.first == r_var;
                    });

 if (it2 != v_d_var_int.end()) {
  v_d_var_int.erase(it2);
 } else {
  throw (std::invalid_argument("Constraint is not inside the CPLEX matrix"));
 }

 // Also, remove from the active constraints vectors
 // TODO: Verify that these vectors are kept in the right order
 active_row_constraints.erase(active_row_constraints.begin() + i);
 active_box_constraints.erase(active_box_constraints.begin() + i);
}

/*--------------------------------------------------------------------------*/
/*----------------------- End File MILPSolver.cpp --------------------------*/
/*--------------------------------------------------------------------------*/ 
