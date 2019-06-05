/*--------------------------------------------------------------------------*/
/*------------------------- File CPXMILPSolver.cpp -------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Implementation of the MILPSolver class.
 *
 * \version 0.20
 *
 * \date 31 - 05 - 2019
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
#include <MILPSolver.h>
#include <OneVarConstraint.h>
#include <LinearFunction.h>
#include <DQuadFunction.h>

#include "CPXMILPSolver.h"

/*--------------------------------------------------------------------------*/
/*------------------------- NAMESPACE AND USING ----------------------------*/
/*--------------------------------------------------------------------------*/

using namespace SMSpp_di_unipi_it;
SMSpp_insert_in_factory_cpp_0(CPXMILPSolver);

/*--------------------------------------------------------------------------*/
/*--------------------- CONSTRUCTOR AND DESTRUCTOR -------------------------*/
/*--------------------------------------------------------------------------*/

CPXMILPSolver::CPXMILPSolver() : MILPSolver() {
 env = nullptr;
 milp = nullptr;
}

CPXMILPSolver::~CPXMILPSolver() {
 if (env) {
  CPXfreeprob(env, &milp);
  CPXcloseCPLEX(&env);
 }
}

/*--------------------------------------------------------------------------*/
/*--------------------- DERIVED METHODS OF BASE CLASS ----------------------*/
/*--------------------------------------------------------------------------*/

int CPXMILPSolver::compute(bool changedvars) {

 int status;
 env = CPXopenCPLEX(&status);
 milp = CPXcreateprob(env, &status, "MILPCPX");

 for (int i = 0; i < numcols; ++i) {
  if (lb[i] == -Inf<double>()) {
   lb[i] = -CPX_INFBOUND;
  }
  if (ub[i] == Inf<double>()) {
   lb[i] = CPX_INFBOUND;
  }
 }

 CPXcopylp(env, milp,
           numcols,
           numrows,
           objsense,
           objective.data(),
           rhs.data(),
           sense.data(),
           matbeg.data(),
           matcnt.data(),
           matind.data(),
           matval.data(),
           lb.data(),
           ub.data(),
           rngval.data());

 const FRealObjective* p_obj;
 p_obj = boost::any_cast<FRealObjective*>(f_Block->get_objective());
 auto p_dquad_fun = dynamic_cast<const DQuadFunction*> (p_obj->get_function());
 if (p_dquad_fun != nullptr) {
  CPXcopyqpsep(env, milp, q_objective.data());
 }
 CPXcopyctype(env, milp, xctype.data());

 process_modifications();

 // TODO: Add this feature as configurable
 CPXwriteprob(env, milp, "problem.lp", nullptr);

 CPXmipopt(env, milp);

 // FIXME: nodes value is not read anywhere
 nodes = CPXgetnodecnt(env, milp);

 // TODO: Support configuration object
 get_var_solution(nullptr);

 status = CPXgetstat(env, milp);

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

/*--------------------------------------------------------------------------*/

Solver::OFValue CPXMILPSolver::get_lb() {

 OFValue lower_bound = 0;

 switch (objsense) {
  case 1: // Minimization problem
   switch (sol_status) {
    case kUnbounded:
     lower_bound = -Inf<OFValue>();
     break;
    case kInfeasible:
     lower_bound = Inf<OFValue>();
     break;
    default:
     CPXgetbestobjval(env, milp, &lower_bound);
     break;
   }
   break;
  case -1: // Maximization problem
   switch (sol_status) {
    case kUnbounded:
     lower_bound = Inf<OFValue>();
     break;
    case kInfeasible:
     lower_bound = -Inf<OFValue>();
     break;
    default:
     CPXgetbestobjval(env, milp, &lower_bound);
     break;
   }
   break;
  default:
   throw (std::runtime_error("Objective type not yet defined"));
   break;
 }

 return lower_bound;
}

/*--------------------------------------------------------------------------*/

Solver::OFValue CPXMILPSolver::get_ub() {

 OFValue upper_bound = 0;

 switch (objsense) {
  case 1: // Minimization problem
   switch (sol_status) {
    case kUnbounded:
     upper_bound = -Inf<OFValue>();
     break;
    case kInfeasible:
     upper_bound = Inf<OFValue>();
     break;
    default:
     CPXgetobjval(env, milp, &upper_bound);
     break;
   }
   break;
  case -1: // Maximization problem
   switch (sol_status) {
    case kUnbounded:
     upper_bound = Inf<OFValue>();
     break;
    case kInfeasible:
     upper_bound = -Inf<OFValue>();
     break;
    default:
     CPXgetobjval(env, milp, &upper_bound);
     break;
   }
   break;
  default:
   throw (std::runtime_error("Objective type not yet defined"));
   break;
 }

 return upper_bound;
}

/*--------------------------------------------------------------------------*/
/*-------------------- PROTECTED FIELDS OF THE CLASS -----------------------*/
/*--------------------------------------------------------------------------*/

void CPXMILPSolver::get_var_solution(Configuration* solc) {

#if DEBUG_COUT
 std::cout << "[DEBUG] ========= MILPSolver::get_var_solution()" << std::endl;
#endif

 auto* tmpx = new double[numcols];

 CPXgetmipx(env, milp, tmpx, 0, numcols - 1);

 int col = 0;

 std::queue<Block*> Q;
 Q.push(f_Block);

 while (!Q.empty()) {
  Block* q_Block = Q.front();
  Q.pop();

  for (auto i : q_Block->get_nested_Blocks()) {
   Q.push(i);
  }

  for (const auto& i : q_Block->get_static_variables()) {
   auto f1 = std::bind(&CPXMILPSolver::set_var_value,
                       this,
                       std::placeholders::_1,
                       std::ref(tmpx),
                       std::ref(col));
   un_any_const_static(i, f1, un_any_type<ColVariable>());
  }

  for (const auto& i : q_Block->get_dynamic_variables()) {
   auto f1 = std::bind(&CPXMILPSolver::set_var_value,
                       this,
                       std::placeholders::_1,
                       std::ref(tmpx),
                       std::ref(col));
   un_any_const_dynamic(i, f1, un_any_type<ColVariable>());
  }

 }

 // After the Objective is computed (evaluated), the solution can be retrieved
 // directly from there.
 try {
  auto p_obj = boost::any_cast<FRealObjective*>(f_Block->get_objective());
  p_obj->compute();
 } catch (boost::bad_any_cast&) {
  throw (std::invalid_argument("Objective is not a FRealObjective"));
 }
  delete []tmpx;
}

/*--------------------------------------------------------------------------*/

void CPXMILPSolver::var_modification(VariableMod* mod) {

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

  int bounds = static_cast<int>(active_bounds[indices[0]].size());
  lu[0] = 'L';
  lu[1] = 'U';
  bd[0] = -CPX_INFBOUND;
  bd[1] = CPX_INFBOUND;
  for (int i = 0; i < bounds; ++i) {
   auto box = active_bounds[indices[0]][i];
   bd[0] = bd[0] > box->get_lhs() ? bd[0] : box->get_lhs();
   bd[1] = bd[1] < box->get_rhs() ? bd[1] : box->get_rhs();
  }
  CPXchgbds(env, milp, 2, indices, lu, bd);
 }

 delete[] lu;
 delete[] bd;
}

/*--------------------------------------------------------------------------*/

void CPXMILPSolver::of_modification(ObjectiveMod* mod) {

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

void CPXMILPSolver::const_modification(ConstraintMod* mod) {

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

void CPXMILPSolver::bound_modification(OneVarConstraintMod* mod) {

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
 int num_bounds = static_cast<int>(active_bounds[indices[0]].size());

 char* lu = nullptr;
 double* bd = nullptr;

 switch (mod->f_type) {

  case (RowConstraintMod::eChgLHS):
   lu = new char[1];
   bd = new double[1];

   lu[0] = 'L';
   bd[0] = -CPX_INFBOUND;

   for (int i = 0; i < num_bounds; ++i) {
    auto box = active_bounds[indices[0]][i];
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
    auto box = active_bounds[indices[0]][i];
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
    auto box = active_bounds[indices[0]][i];
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

void CPXMILPSolver::function_modification(FunctionMod* mod) {

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

void CPXMILPSolver::dynamic_modification(BlockModAD* mod) {

 switch (mod->f_type) {

  case (BlockModAD::eAddConst): {
   // TODO: Handle OneVarConstraints
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
   // TODO: Handle OneVarConstraints
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
   throw (std::invalid_argument("Unknown type of BlockAD"));
   break;
 }
}

/*--------------------------------------------------------------------------*/

void CPXMILPSolver::add_dynamic_constraint(FRowConstraint* p_const) {

 auto p_fun = dynamic_cast<const LinearFunction*>(p_const->get_function());
 if (p_fun == nullptr) {
  throw (std::invalid_argument("The Constraint is not linear"));
 }

 int nzcnt = p_const->get_num_active_var();
 std::array<int, 2> rmatbeg = {0, nzcnt};
 std::vector<int> rmatind(nzcnt);
 std::vector<double> rmatval(nzcnt);

 std::array<double, 1> rhs{};
 std::array<double, 1> rngval{};
 std::array<int, 1> indices{};
 std::array<char, 1> sense{};

 int i = 0;

 // Get the coefficients to fill the matrix
 for (auto it = p_const->begin(); it != p_const->end(); ++it) {

  auto p_var = dynamic_cast<ColVariable*>(&*it);
  auto it1 = find_if(v_int_d_var.begin(),
                     v_int_d_var.end(),
                     [&](MILPSolver::int_var pair) {
                      return pair.second == p_var;
                     });
  if (it1 != v_int_d_var.end()) {
   rmatind[i] = (it1->first);
  } else {
   rmatind[i] = index_of_variable(p_var);
  }

  rmatval[i] = p_fun->get_coefficient(p_var);
  active_constraints[rmatind[i]].push_back(p_const);
  ++i;
 }

 auto const_lhs = p_const->get_lhs();
 auto const_rhs = p_const->get_rhs();

 if (const_lhs == const_rhs) {
  sense[0] = 'E';
  rhs[0] = const_rhs;

 } else if (const_lhs == -Inf<double>()) {
  sense[0] = 'L';
  rhs[0] = const_rhs;

 } else if (const_rhs == Inf<double>()) {
  sense[0] = 'G';
  rhs[0] = const_lhs;

 } else {
  sense[0] = 'R';
  rhs[0] = const_lhs;
  rngval[0] = const_rhs - const_lhs;
  indices[0] = v_d_const_int.back().second + 1;
 }

 v_d_const_int.emplace_back(p_const, v_d_const_int.back().second + 1);
 v_int_d_const.emplace_back(v_int_d_const.back().first + 1, p_const);
 std::sort(v_d_const_int.begin(), v_d_const_int.end());
 std::sort(v_int_d_const.begin(), v_int_d_const.end());

 CPXaddrows(env, milp, 0, 1, nzcnt, rhs.data(),
                                    sense.data(),
                                    rmatbeg.data(),
                                    rmatind.data(),
                                    rmatval.data(),
                                    nullptr,
                                    nullptr);
 if (sense[0] == 'R') {
  CPXchgrngval (env, milp, 1, indices.data(), rngval.data());
 }
}

/*--------------------------------------------------------------------------*/

void CPXMILPSolver::add_dynamic_variable(ColVariable* p_var) {

 std::vector<FRowConstraint*> var_constraints;
 std::vector<OneVarConstraint*> var_bounds;

 int nzcnt = 0;
 for (auto stuff : p_var->active_stuff()) {
  auto constraint = dynamic_cast<FRowConstraint*>(stuff);
  if (constraint != nullptr) {
   var_constraints.push_back(constraint);
   ++nzcnt;
  }
  auto bound = dynamic_cast<OneVarConstraint*>(stuff);
  if (bound != nullptr) {
   var_bounds.push_back(bound);
  }
 }

 std::array<int, 2> cmatbeg = {0, nzcnt};
 std::vector<int> cmatind(nzcnt);
 std::vector<double> cmatval(nzcnt);

 int i = 0;

 // We need the coefficients for this variable in each constraint
 for (auto* p_const : var_constraints) {

  auto p_fun = dynamic_cast<const LinearFunction*>(p_const->get_function());
  if (p_fun == nullptr) {
   throw (std::invalid_argument("The Constraint is not linear"));
  }

  auto it = find_if(v_int_d_const.begin(),
                    v_int_d_const.end(),
                    [&](MILPSolver::int_const pair) {
                     return pair.second == p_const;
                    });
  if (it != v_int_d_const.end()) {
   cmatind[i] = it->first;
  } else {
   cmatind[i] = index_of_constraint(p_const);
  }

  cmatval[i] = p_fun->get_coefficient(p_var);
  ++i;
 }

 std::array<double, 1> lb{};
 std::array<double, 1> ub{};

 lb[0] = p_var->get_lb() == -Inf<double>() ? -CPX_INFBOUND : p_var->get_lb();
 ub[0] = p_var->get_ub() == Inf<double>() ? CPX_INFBOUND : p_var->get_ub();

 for (auto bnd : var_bounds) {
  lb[0] = lb[0] > bnd->get_lhs() ? lb[0] : bnd->get_lhs();
  ub[0] = ub[0] < bnd->get_rhs() ? ub[0] : bnd->get_rhs();
 }

 active_constraints.emplace_back(var_constraints);
 active_bounds.emplace_back(var_bounds);

 v_d_var_int.emplace_back(p_var, v_d_var_int.back().second + 1);
 v_int_d_var.emplace_back(v_int_d_var.back().first + 1, p_var);
 std::sort(v_d_var_int.begin(), v_d_var_int.end());
 std::sort(v_int_d_var.begin(), v_int_d_var.end());

 CPXaddcols(env, milp, 1, nzcnt, nullptr, cmatbeg.data(),
            cmatind.data(), cmatval.data(), lb.data(), ub.data(), nullptr);
}

/*--------------------------------------------------------------------------*/

void CPXMILPSolver::remove_dynamic_constraint(FRowConstraint* p_const) {

 int index = 0;
 auto it1 = find_if(v_int_d_const.begin(),
                    v_int_d_const.end(),
                    [&](MILPSolver::int_const pair) {
                     return pair.second == p_const;
                    });
 if (it1 != v_int_d_const.end()) {
  index = it1->first;
  v_int_d_const.erase(it1);
 } else {
  throw (std::invalid_argument("Cannot find the Constraint"));
 }

 auto it2 = find_if(v_d_const_int.begin(),
                    v_d_const_int.end(),
                    [&](MILPSolver::const_int pair) {
                     return pair.first == p_const;
                    });
 if (it2 != v_d_const_int.end()) {
  v_d_const_int.erase(it2);
 } else {
  throw (std::invalid_argument("Cannot find the Constraint"));
 }

 CPXdelrows(env, milp, index, index + 1);

 for (auto it: v_d_const_int) {
  if (it.second > index) {
   it.second--;
  }
 }
 for (auto it: v_int_d_const) {
  if (it.first > index) {
   it.first--;
  }
 }

 for (auto& constraints: active_constraints) {
  auto constraint = find(constraints.begin(), constraints.end(), p_const);
  if (constraint != constraints.end()) {
   constraints.erase(constraint);
  }
 }
}

/*--------------------------------------------------------------------------*/

void CPXMILPSolver::remove_dynamic_variable(ColVariable* p_var) {

 int index = 0;
 auto it1 = find_if(v_int_d_var.begin(),
                   v_int_d_var.end(),
                   [&](MILPSolver::int_var pair) {
                    return pair.second == p_var;
                   });
 if (it1 != v_int_d_var.end()) {
  index = it1->first;
  v_int_d_var.erase(it1);
 } else {
  throw (std::invalid_argument("Cannot find the Variable"));
 }

 auto it2 = find_if(v_d_var_int.begin(),
                    v_d_var_int.end(),
                    [&](MILPSolver::var_int pair) {
                     return pair.first == p_var;
                    });

 if (it2 != v_d_var_int.end()) {
  v_d_var_int.erase(it2);
 } else {
  throw (std::invalid_argument("Cannot find the Variable"));
 }

 CPXdelcols(env, milp, index, index + 1);

 for (auto it: v_d_var_int) {
  if (it.second > index) {
   it.second--;
  }
 }
 for (auto it: v_int_d_var) {
  if (it.first > index) {
   it.first--;
  }
 }

 active_constraints.erase(active_constraints.begin() + index);
 active_bounds.erase(active_bounds.begin() + index);
}

/*--------------------------------------------------------------------------*/
/*--------------------- PRIVATE FIELDS OF THE CLASS ------------------------*/
/*--------------------------------------------------------------------------*/

void CPXMILPSolver::set_var_value(ColVariable& lvar, double* tmpx, int& i) {
#if DEBUG_COUT
 std::cout << "[DEBUG] ========= MILPSolver::set_var_value()" << std::endl;
 std::cout << "[DEBUG] i     = " << i << std::endl;
 std::cout << "[DEBUG] value = " << tmpx[i] << std::endl;
#endif
 lvar.set_value(tmpx[i]);
 i++;
}

/*--------------------------------------------------------------------------*/
/*--------------------- End File CPXMILPSolver.cpp -------------------------*/
/*--------------------------------------------------------------------------*/
