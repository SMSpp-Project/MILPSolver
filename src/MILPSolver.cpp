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

#define DEBUG_COUT 1 // TODO: Remove this and all the printouts when done

/*--------------------------------------------------------------------------*/
/*------------------------- NAMESPACE AND USING ----------------------------*/
/*--------------------------------------------------------------------------*/

using namespace SMSpp_di_unipi_it;
SMSpp_insert_in_factory_cpp_0(MILPSolver);

/*--------------------------------------------------------------------------*/
/*--------------------- CONSTRUCTOR AND DESTRUCTOR -------------------------*/
/*--------------------------------------------------------------------------*/

MILPSolver::MILPSolver() : Solver() {}

MILPSolver::~MILPSolver() {}

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

const std::vector<double>& MILPSolver::get_objective() const {
 return objective;
}

const std::vector<double>& MILPSolver::get_q_objective() const {
 return q_objective;
}

const std::vector<double>& MILPSolver::get_rhs() const {
 return rhs;
}

const std::vector<double>& MILPSolver::get_rngval() const {
 return rngval;
}

const std::vector<char>& MILPSolver::get_sense() const {
 return sense;
}

const std::vector<int>& MILPSolver::get_matbeg() const {
 return matbeg;
}

const std::vector<int>& MILPSolver::get_matcnt() const {
 return matcnt;
}

const std::vector<int>& MILPSolver::get_matind() const {
 return matind;
}

const std::vector<double>& MILPSolver::get_matval() const {
 return matval;
}

const std::vector<double>& MILPSolver::get_lb() const {
 return lb;
}

const std::vector<double>& MILPSolver::get_ub() const {
 return ub;
}

const std::vector<char>& MILPSolver::get_xctype() const {
 return xctype;
}

/*--------------------------------------------------------------------------*/
/*-------------------------------- SET_BLOCK -------------------------------*/
/*--------------------------------------------------------------------------*/

void MILPSolver::set_Block(Block* block) {

 if (f_Block) {
  f_Block->unregister_Solver(this);
  clear_matrices();
 }

 f_Block = block;

 /*
  * Passing all the data of the Block to the LP.
  *
  * Following a Breadth First Search we proceed with scanning the received
  * Block and all of each corresponding children if any, in order to populate
  * the LP data.
  */

 numrows = 0;
 numcols = 0;
 nzelements = 0;

 std::queue<Block*> Q;

 // First loop on the queue to count variables and constraints
 Q.push(f_Block);

 while (!Q.empty()) {
  Block* q_Block = Q.front();
  Q.pop();

  for (auto i : q_Block->get_nested_Blocks()) {
   Q.push(i);
  }

#if DEBUG_COUT
  std::cout << "[DEBUG] ========= MILPSolver::set_Block() counting static constraints" << std::endl;
#endif
  for (const auto& i : q_Block->get_static_constraints()) {
   auto f1 = std::bind(&MILPSolver::count_constraints,
                       this,
                       std::placeholders::_1,
                       std::ref(numrows));
   un_any_const_static(i, f1, un_any_type<FRowConstraint>());
  }

#if DEBUG_COUT
  std::cout << "[DEBUG] ========= MILPSolver::set_Block() counting dynamic constraints" << std::endl;
#endif
  for (const auto& i : q_Block->get_dynamic_constraints()) {
   auto f1 = std::bind(&MILPSolver::count_constraints,
                       this,
                       std::placeholders::_1,
                       std::ref(numrows));
   un_any_const_dynamic(i, f1, un_any_type<FRowConstraint>());
  }

#if DEBUG_COUT
  std::cout << "[DEBUG] ========= MILPSolver::set_Block() counting static variables" << std::endl;
#endif
  for (const auto& i : q_Block->get_static_variables()) {
   auto f1 = std::bind(&MILPSolver::count_variables,
                       this,
                       std::placeholders::_1,
                       std::ref(numcols));
   un_any_const_static(i, f1, un_any_type<ColVariable>());
  }

#if DEBUG_COUT
  std::cout << "[DEBUG] ========= MILPSolver::set_Block() counting dynamic variables" << std::endl;
#endif
  for (const auto& i : q_Block->get_dynamic_variables()) {
   auto f1 = std::bind(&MILPSolver::count_variables,
                       this,
                       std::placeholders::_1,
                       std::ref(numcols));
   un_any_const_dynamic(i, f1, un_any_type<ColVariable>());
  }
#if DEBUG_COUT
  std::cout << "[DEBUG] ========= MILPSolver::set_Block() nonzero elements" << std::endl;
#endif
  int cnt = 0;
  for (const auto& i : q_Block->get_static_variables()) {
   auto f1 = std::bind(&MILPSolver::count_nzelements,
                       this,
                       std::placeholders::_1,
                       std::ref(nzelements),
                       std::ref(cnt));
   un_any_const_static(i, f1, un_any_type<ColVariable>());
  }

  for (const auto& i : q_Block->get_dynamic_variables()) {
   auto f1 = std::bind(&MILPSolver::count_nzelements,
                       this,
                       std::placeholders::_1,
                       std::ref(nzelements),
                       std::ref(cnt));
   un_any_const_dynamic(i, f1, un_any_type<ColVariable>());
  }
 } // End of while loop on Block queue

#if DEBUG_COUT
 std::cout << "[DEBUG] ========= MILPSolver::set_Block() after counting" << std::endl;
 std::cout << "constraints/numrows = " << numrows << std::endl;
 std::cout << "variables/numcols =   " << numcols << std::endl;
 std::cout << "nzelements =          " << nzelements << std::endl;
#endif

 // The +1 is needed by generic interface
 matbeg.resize(numcols + 1);
 matbeg[numcols] = nzelements;

 matcnt.resize(numcols);
 matind.resize(nzelements);
 matval.resize(nzelements);
 // indexed.resize(static_cast<unsigned long>(numcols));
 rhs.resize(numrows);
 rngval.resize(numrows);
 sense.resize(numrows);
 objective.resize(numcols);
 q_objective.resize(numcols);
 std::fill(objective.begin(), objective.end(), 0);
 std::fill(objective.begin(), objective.end(), 0);
 lb.resize(numcols);
 ub.resize(numcols);
 xctype.resize(numcols);

 // Second loop to scan the constraints
 Q.push(f_Block);

 int col = 0;

 while (!Q.empty()) {
  Block* q_Block = Q.front();
  Q.pop();

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

  for (const auto& i : q_Block->get_static_constraints()) {
   // Variable used to locate the first element of each type of constraint
   int first = 0;
   auto f1 = std::bind(&MILPSolver::scan_static_constraint,
                       this,
                       std::placeholders::_1,
                       std::ref(first),
                       std::ref(col));
   un_any_const_static(i, f1, un_any_type<FRowConstraint>());
  }

  for (const auto& i : q_Block->get_dynamic_constraints()) {
   auto f1 = std::bind(&MILPSolver::scan_dynamic_constraint,
                       this,
                       std::placeholders::_1,
                       std::ref(col));
   un_any_const_dynamic(i, f1, un_any_type<FRowConstraint>());
  }
 } // End of while loop on Block queue

 std::sort(v_s_const_int.begin(), v_s_const_int.end());
 std::sort(v_int_s_const.begin(), v_int_s_const.end());
 std::sort(v_d_const_int.begin(), v_d_const_int.end());
 std::sort(v_int_d_const.begin(), v_int_d_const.end());

#if DEBUG_COUT
 std::cout << "[DEBUG] ========= MILPSolver::set_Block() after constraint scan" << std::endl;
#endif

 // Third loop to scan the variables
 Q.push(f_Block);

 int var = 0;
 while (!Q.empty()) {
  Block* q_Block = Q.front();
  Q.pop();

  for (auto i : q_Block->get_nested_Blocks()) {
   Q.push(i);
  }

  for (const auto& i : q_Block->get_static_variables()) {
   int first = 0;
   auto f1 = std::bind(&MILPSolver::scan_static_variable,
                       this,
                       std::placeholders::_1,
                       std::ref(first),
                       std::ref(var));
   un_any_const_static(i, f1, un_any_type<ColVariable>());
  }

  for (const auto& i : q_Block->get_dynamic_variables()) {
   auto f1 = std::bind(&MILPSolver::scan_dynamic_variable,
                       this,
                       std::placeholders::_1,
                       std::ref(var));
   un_any_const_dynamic(i, f1, un_any_type<ColVariable>());
  }
 } // End of while loop on Block queue

 std::sort(v_s_var_int.begin(), v_s_var_int.end());
 std::sort(v_int_s_var.begin(), v_int_s_var.end());
 std::sort(v_d_var_int.begin(), v_d_var_int.end());
 std::sort(v_int_d_var.begin(), v_int_d_var.end());

 const FRealObjective* p_obj;
 try {
  p_obj = boost::any_cast<FRealObjective*>(f_Block->get_objective());
  switch (p_obj->get_sense()) {
   case (Objective::eMax):
    objsense = -1;
    break;
   case (Objective::eMin):
    objsense = 1;
    break;
   default:
    objsense = 0;
    break;
  }
 } catch (boost::bad_any_cast&) {
  throw (std::invalid_argument("Objective is not a FRealObjective"));
 }

 // Fourth loop to scan the objective(s?)
 // FIXME: I am not sure that this scan should be done in a loop for the subs
 Q.push(f_Block);
 while (!Q.empty()) {
  Block* q_Block = Q.front();
  Q.pop();

  for (auto i : q_Block->get_nested_Blocks()) {
   Q.push(i);
  }

  try {
   p_obj = boost::any_cast<FRealObjective*>(q_Block->get_objective());
  } catch (boost::bad_any_cast&) {
   throw (std::invalid_argument("Objective is not a FRealObjective"));
  }

  scan_objective(p_obj);
 } // End of while loop on Block queue

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

 std::cout << "[DEBUG] rngval      = [";
 for (int i = 0; i < numrows; ++i) {
  std::cout << " " << rngval[i];
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
}
/*--------------------------------------------------------------------------*/
/*-------------------- METHODS FOR PROBLEM DESCRIPTION ---------------------*/
/*--------------------------------------------------------------------------*/

int MILPSolver::index_of_variable(ColVariable* p_var) {
 try {
  return index_of_dynamic_variable(p_var);
 } catch (...) {
  return index_of_static_variable(p_var);
 }
}

int MILPSolver::index_of_static_variable(ColVariable* p_var) {

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
  throw (std::invalid_argument("Variable not found"));
 }
 return i;
}

int MILPSolver::index_of_dynamic_variable(ColVariable* p_var) {
 auto it = find_if(v_d_var_int.begin(),
                   v_d_var_int.end(),
                   [&](var_int pair) {
                    return pair.first == p_var;
                   });
 if (it != v_d_var_int.end()) {
  return it->second;
 } else {
  throw (std::invalid_argument("Variable not found"));
 }
}

/*--------------------------------------------------------------------------*/
int MILPSolver::index_of_constraint(FRowConstraint* p_const) {
 try {
  return index_of_dynamic_constraint(p_const);
 } catch (...) {
  return index_of_static_constraint(p_const);
 }
}

int MILPSolver::index_of_static_constraint(FRowConstraint* p_const) {

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
  throw (std::invalid_argument("Constraint not found"));
 }

 return i;
}

int MILPSolver::index_of_dynamic_constraint(FRowConstraint* p_const) {
 auto it = find_if(v_d_const_int.begin(),
                   v_d_const_int.end(),
                   [&](const_int pair) {
                    return pair.first == p_const;
                   });
 if (it != v_d_const_int.end()) {
  return it->second;
 } else {
  throw (std::invalid_argument("Constraint not found"));
 }
}

/*--------------------------------------------------------------------------*/

ColVariable* MILPSolver::static_variable_with_index(int i) {

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
  throw (std::invalid_argument("Index not found"));
 }

 return p_var;
}

/*--------------------------------------------------------------------------*/

FRowConstraint* MILPSolver::static_constraint_with_index(int i) {

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
  throw (std::invalid_argument("Index not found"));
 }

 return p_const;
}

/*--------------------------------------------------------------------------*/
/*------------- AUXILIARY METHODS FOR POPULATING THE PROBLEM  --------------*/
/*--------------------------------------------------------------------------*/

void MILPSolver::count_constraints(FRowConstraint& constraint, int& n_rows) {
#if DEBUG_COUT
 std::cout << "[DEBUG] ========= MILPSolver::count_constraints()  " << n_rows << " " << constraint;
#endif

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
 std::cout << "[DEBUG] ========= MILPSolver::count_variables()    " << n_cols << " " << variable;
#endif
 ++n_cols;
}

/*--------------------------------------------------------------------------*/

void MILPSolver::count_nzelements(ColVariable& variable,
                                  int& nz_elements,
                                  int& cnt) {
#if DEBUG_COUT
 std::cout << "[DEBUG] ========= MILPSolver::count_nzelements()   " << cnt << " " << variable;
 std::cout << "[DEBUG] The active stuff is:" << std::endl;
#endif

 /*
  * Since counting non-zero elements requires checking if each active thing
  * is a FRowConstraint, we populate active_constraints and
  * active_bounds here so we don't have to loop over active stuff
  * once again later.
  */

 if (active_constraints.size() < numcols) {
  active_constraints.resize(static_cast<unsigned long>(numcols));
  active_bounds.resize(static_cast<unsigned long>(numcols));
 }

 for (auto i : variable.active_stuff()) {
  auto row = dynamic_cast<FRowConstraint*>(i);
  if (row != nullptr) {
#if DEBUG_COUT
   std::cout << "[DEBUG] " << *row;
#endif
   active_constraints[cnt].push_back(row);
   ++nz_elements;
  }
  auto box = dynamic_cast<OneVarConstraint*>(i);
  if (box != nullptr) {
#if DEBUG_COUT
   std::cout << "[DEBUG] " << *box;
#endif
   active_bounds[cnt].push_back(box);
  }
#if DEBUG_COUT
  auto obj = dynamic_cast<Objective*>(i);
  if (obj != nullptr) {
   std::cout << "[DEBUG] " << *obj;
  }
#endif
 }
 ++cnt;
}

/*--------------------------------------------------------------------------*/

void MILPSolver::scan_static_variable(ColVariable& var, int& first, int& i) {

#if DEBUG_COUT
 std::cout << "[DEBUG] ========= MILPSolver::scan_static_variable() ";
 std::cout << i << " " << var;
#endif

 if (first == 0) {
  v_s_var_int.emplace_back(&var, i);
  v_int_s_var.emplace_back(i, &var);
 }

 lb[i] = var.get_lb();
 ub[i] = var.get_ub();

 int num_bounds = static_cast<int>(active_bounds[i].size());
 for (int j = 0; j < num_bounds; ++j) {
  auto bound = active_bounds[i][j];
  lb[i] = lb[i] > bound->get_lhs() ? lb[i] : bound->get_lhs();
  ub[i] = ub[i] < bound->get_rhs() ? ub[i] : bound->get_rhs();
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

 int nz_elements = static_cast<int>(active_constraints[i].size());
 matcnt[i] = nz_elements;

 if (i == 0) {
  matbeg[i] = 0;
 } else {
  matbeg[i] = matbeg[i - 1] + matcnt[i - 1];
 }

 for (int j = 0; j < nz_elements; ++j) {

  auto p_const = dynamic_cast<FRowConstraint*> (active_constraints[i][j]);
  auto p_fun = dynamic_cast<const LinearFunction*> (p_const->get_function());

  auto it = std::find_if(p_fun->get_v_var().begin(),
                         p_fun->get_v_var().end(),
                         [&](LinearFunction::coeff_pair pair) {
                          return pair.first == &var;
                         });

  if (it != p_fun->get_v_var().end()) {
   matval[matbeg[i] + j] = it->second;
   matind[matbeg[i] + j] = index_of_constraint(p_const);
  } else {
   // This should never happen because we are looping on the active contraints
   throw (std::invalid_argument("This ColVariable is not active in the examined FRowConstraint"));
  }
 }
 ++first;
 ++i;
}

/*--------------------------------------------------------------------------*/

void MILPSolver::scan_dynamic_variable(ColVariable& var, int& i){

#if DEBUG_COUT
 std::cout << "[DEBUG] ========= MILPSolver::scan_dynamic_variable() ";
 std::cout << i << " " << var;
#endif

 v_d_var_int.emplace_back(&var, i);
 v_int_d_var.emplace_back(i, &var);

 lb[i] = var.get_lb();
 ub[i] = var.get_ub();

 int num_bounds = static_cast<int>(active_bounds[i].size());
 for (int j = 0; j < num_bounds; ++j) {
  auto bound = active_bounds[i][j];
  lb[i] = lb[i] > bound->get_lhs()? lb[i] : bound->get_lhs();
  ub[i] = ub[i] < bound->get_rhs()? ub[i] : bound->get_rhs();
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

 int nz_elements = static_cast<int>(active_constraints[i].size());
 matcnt[i] = nz_elements;

 if (i == 0) {
  matbeg[i] = 0;
 } else {
  matbeg[i] = matbeg[i - 1] + matcnt[i - 1];
 }

 for (int j = 0; j < nz_elements; ++j) {

  auto p_const = dynamic_cast<FRowConstraint*> (active_constraints[i][j]);
  auto p_fun = dynamic_cast<const LinearFunction*> (p_const->get_function());

  auto it = std::find_if(p_fun->get_v_var().begin(),
                         p_fun->get_v_var().end(),
                         [&](LinearFunction::coeff_pair pair) {
                          return pair.first == &var;
                         });

  if (it != p_fun->get_v_var().end()) {
   matval[matbeg[i] + j] = it->second;
   matind[matbeg[i] + j] = index_of_constraint(p_const);
  } else {
   // This should never happen because we are looping on the active contraints
   throw (std::invalid_argument("This ColVariable is not active in the examined FRowConstraint"));
  }
 }
 ++i;
}

/*--------------------------------------------------------------------------*/

void MILPSolver::scan_static_constraint(FRowConstraint& p_const, int& first, int& i) {

#if DEBUG_COUT
 std::cout << "[DEBUG] ========= MILPSolver::scan_static_constraint()  " << p_const;
#endif
 auto lin_fun = dynamic_cast<const LinearFunction*>(p_const.get_function());
 if (lin_fun == nullptr) {
  throw (std::invalid_argument("The Constraint is not linear"));
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

 if (const_lhs == const_rhs) {
  // LHS <= function <= RHS, with LHS = RHS
  // becomes:
  // function = RHS
  sense[i] = 'E';
  rhs[i] = const_rhs;

 } else if (const_lhs == -Inf<double>()) {
  // -inf <= function <= RHS
  // becomes:
  // function <= RHS
  sense[i] = 'L';
  rhs[i] = const_rhs;

 } else if (const_rhs == Inf<double>()) {
  // LHS <= function <= inf
  // becomes:
  // function >= LHS
  sense[i] = 'G';
  rhs[i] = const_lhs;

 } else {
  // LHS <= function <= RHS
  // becomes:
  // LHS <= function <= LHS + (range),
  // with range = RHS - LHS
  sense[i] = 'R';
  rhs[i] = const_lhs;
  rngval[i] = const_rhs - const_lhs;
 }

 if (first == 0) {
  v_s_const_int.emplace_back(&p_const, i);
  v_int_s_const.emplace_back(i, &p_const);
 }

 ++first;
 ++i;
}

/*--------------------------------------------------------------------------*/

void MILPSolver::scan_dynamic_constraint(FRowConstraint& p_const, int& i) {

#if DEBUG_COUT
 std::cout << "[DEBUG] ========= MILPSolver::scan_dynamic_constraint() " << p_const;
#endif
 auto lin_fun = dynamic_cast<const LinearFunction*>(p_const.get_function());
 if (lin_fun == nullptr) {
  throw (std::invalid_argument("The Constraint is not linear"));
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

 if (const_lhs == const_rhs) {
  // LHS <= function <= RHS, with LHS = RHS
  // becomes:
  // function = RHS
  sense[i] = 'E';
  rhs[i] = const_rhs;

 } else if (const_lhs == -Inf<double>()) {
  // -inf <= function <= RHS
  // becomes:
  // function <= RHS
  sense[i] = 'L';
  rhs[i] = const_rhs;

 } else if (const_rhs == Inf<double>()) {
  // LHS <= function <= inf
  // becomes:
  // function >= LHS
  sense[i] = 'G';
  rhs[i] = const_lhs;

 } else {
  // LHS <= function <= RHS
  // becomes:
  // LHS <= function <= LHS + (range),
  // with range = RHS - LHS
  sense[i] = 'R';
  rhs[i] = const_lhs;
  rngval[i] = const_rhs - const_lhs;
 }

 v_d_const_int.emplace_back(&p_const, i);
 v_int_d_const.emplace_back(i, &p_const);
 ++i;
}

/*--------------------------------------------------------------------------*/

void MILPSolver::scan_objective(const FRealObjective* obj) {
#if DEBUG_COUT
 std::cout << "[DEBUG] ========= MILPSolver::scan_objective() " << *obj;
#endif

 auto lin_fun = dynamic_cast<const LinearFunction*> (obj->get_function());
 int k;

 if (lin_fun != nullptr) {
  for (auto el : lin_fun->get_v_var()) {
   k = index_of_variable(el.first);
   objective[k] = el.second;
  }
 } else {
  auto dquad_fun = dynamic_cast<const DQuadFunction*> (obj->get_function());
  if (dquad_fun != nullptr) {
   for (auto el : dquad_fun->get_v_var()) {
    // DQuadFunction::get_v_var() returns std::tuples of 3 elements
    k = index_of_variable(std::get<0>(el));
    objective[k] = std::get<1>(el);
    q_objective[k] = std::get<2>(el);
   }
  } else {
   throw (std::invalid_argument("Unknown type of Objective Function"));
  }
 }
}

/*--------------------------------------------------------------------------*/
void MILPSolver::process_modifications() {
 /* TODO: Write a better demux function
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

   std::cout << *mod;
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

void MILPSolver::clear_matrices() {
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
}

/*--------------------------------------------------------------------------*/

// TODO: Support cuts callback
// extern int t_pc;
// int mycallback(CPXCENVptr env,
//                void* cbdata,
//                int wherefrom,
//                void* cbhandle,
//                int* useraction_p) {
//  return static_cast<MILPSolver*>(cbhandle)->Callback(env, cbdata, wherefrom, useraction_p);
// }

// I'm not sure it is ideal to split them according the type of the value
// void MILPSolver::set_par(const int par, const int value) {
//  switch (par) {
//   case intLogVerb:
//    f_log_verb = value;
//    CPXsetintparam(env, CPX_PARAM_SCRIND, value);
//    break; //1
//   case intMaxSol:
//    f_max_sol = value;
//    CPXsetintparam(env, CPXPARAM_MIP_Limits_Solutions, value);
//    break;
//   case intMaxIter:
//    f_max_iter = value;
//    CPXsetlongparam(env, CPXPARAM_MIP_Limits_Nodes, value);
//   default:
//    CPXsetintparam(env, par - intLastAlgPar, value);
//    break;
//  }
// }

/*--------------------------------------------------------------------------*/

// void MILPSolver::set_par( const int par, const double value ) {
//  switch (par) {
//   case dblMaxTime:
//    f_max_time = value;
//    CPXsetdblparam(env, CPXPARAM_TimeLimit, value);
//    break; //10
//   case dblRelAcc:
//    f_rel_acc = value;
//    // CPXsetdblparam(env, ,value );
//    break;
//   case dblAbsAcc:
//    f_abs_acc = value;
//    // CPXsetdblparam(env, ,value );
//    break;
//   case dblUpCutOff:
//    f_up_cutoff = value;
//    CPXsetdblparam(env, CPXPARAM_MIP_Tolerances_UpperCutoff, value);
//    break;
//   case dblLwCutOff:
//    f_lw_cutoff = value;
//    CPXsetdblparam(env, CPXPARAM_MIP_Tolerances_LowerCutoff, value);
//    break;
//   case dblRAccSol:
//    f_r_acc_sol = value;
//    CPXsetdblparam(env, CPXPARAM_MIP_Pool_RelGap, value);
//    break;
//   case dblAAccSol:
//    f_a_acc_sol = value;
//    CPXsetdblparam(env, CPXPARAM_MIP_Pool_AbsGap, value);
//    break;
//   case dblFAccSol:
//    f_f_acc_sol = value;
//    // CPXsetdblparam(env, ,value );
//    break;
//   default:
//    CPXsetdblparam(env, par - dblLastAlgPar, value);
//    break;
//  }
//  }

/*--------------------------------------------------------------------------*/

// void MILPSolver::set_par( const int par , const long value )
// {
//  switch( par ) {
//   case( kMaxIter ): f_max_iter = value; CPXsetlongparam(env, CPXPARAM_MIP_Limits_Nodes, value); break; // 0 or max
//   default:                              CPXsetlongparam(env, par-kLastAlgPar,value ); break;
//   }
//  }

/*--------------------------------------------------------------------------*/

// int MILPSolver::Callback(CPXCENVptr env, void* cbdata, int wherefrom, int* useraction_p) {
//
//  // Pass to Block the solution of current node
//  double* tmpx = new double[numcols];
//
//  *useraction_p = CPX_CALLBACK_DEFAULT;
//
//  // Obtain solution from CPLEX
//  CPXgetcallbacknodex(env, cbdata, wherefrom, tmpx, 0, numcols - 1);
//
//   // Retrieve Solution for all the Variables with respect to the order
//  int col = 0;
//
//  std::queue<Block*> Q; // Creating the queue
//  Q.push(f_Block);      // Passing the root, i.e. the father block
//
//  while (!Q.empty())    // Iterating for the father block and all children
//  {
//   Block* q_Block = Q.front(); // Block to be examined
//   Q.pop(); //take out from the queue the examined block
//
//   for (auto i : q_Block->get_nested_Blocks()) {
//    Q.push(i);
//   }
//
//   // Passing values for all static Variables
//   for (const auto& i : q_Block->get_static_variables()) {
//    auto f1 = std::bind(&MILPSolver::set_var_value,
//                        this,
//                        std::placeholders::_1,
//                        std::ref(tmpx),
//                        std::ref(col));
//    un_any_const_static(i, f1, un_any_type<ColVariable>());
//   }
//
//   // We scan all the dynamic Variables
//   for (const auto& i : q_Block->get_dynamic_variables()) {
//    auto f1 = std::bind(&MILPSolver::set_var_value,
//                        this,
//                        std::placeholders::_1,
//                        std::ref(tmpx),
//                        std::ref(col));
//    un_any_const_dynamic(i, f1, un_any_type<ColVariable>());
//   }
//
//   // we call this method where is expected to be written the mechanism
//   // for constructing new constraints
//   q_Block->generate_dynamic_constraints();
//  } // while loop
//
//  // Checking if new constraints have been creating
//  // by looking at the list of modifications

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
      mycutind[i] = index_of_static_variable(it->first);
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
 //
 // return( 0 );
 // }

/*--------------------------------------------------------------------------*/
/*----------------------- End File MILPSolver.cpp --------------------------*/
/*--------------------------------------------------------------------------*/ 
