//
// Created by Niccolò Iardella on 2019-03-07.
//

#include <iostream>
#include <fstream>

#include "SimpleMILPBlock.h"
#include "MILPSolver.h"
#include "CPXMILPSolver.h"
#include "LinearFunction.h"
#include "Configuration.h"

using namespace std;
using namespace SMSpp_di_unipi_it;

int main(int argc, char** argv) {

 if (argc != 2) {
  std::cerr << "Usage: " << argv[0] << " MILP_file_name" << std::endl;
  return (1);
 }

 std::ifstream file(argv[1]);
 if (!file.is_open()) {
  std::cerr << "Error: cannot open file " << argv[1] << std::endl;
  return (1);
 }

 auto block = Block::new_Block("SimpleMILPBlock");
 file >> *block;
 std::cout << *block;

 // Solver* solver = Solver::new_Solver("CPXMILPSolver");
 Solver* solver = new CPXMILPSolver();

 ComputeConfig conf;
 std::pair<std::string, std::string> problem_name = {"strProblemName", "testCPX"};
 std::pair<std::string, std::string> output_file = {"strOutputFile", "output.lp"};
 conf.str_pars.emplace_back(problem_name);
 conf.str_pars.emplace_back(output_file);

 block->register_Solver(solver);

 // First solve
 int status = solver->compute();
 if (solver->has_var_solution()) {
  solver->get_var_solution();
 }

 auto smilpblock = dynamic_cast<SimpleMILPBlock*>(block);
 auto obj = dynamic_cast<FRealObjective*>(smilpblock->get_objective());
 auto obj_f = obj->get_function();
 obj_f->compute();

 std::cout << "Status = " << status << std::endl;
 auto vars = smilpblock->get_x();
 for (auto &i : vars) {
  std::cout << "Variable value =  " << i.get_value() << std::endl;
 }
 std::cout << "Function value =  " << obj_f->get_value() << std::endl;

 // Testing FunctionMod on Objective
 auto lf0 = dynamic_cast<LinearFunction*>(obj_f);
 lf0->modify_coefficient(0, -8);
 status = solver->compute();

 // Testing ObjectiveMod, maximize
 obj->set_sense(Objective::eMax);
 status = solver->compute();

 // Testing ObjectiveMod, minimize
 obj->set_sense(Objective::eMin);
 status = solver->compute();

 // Testing RowConstraintMod, set RHS
 auto constraints = boost::any_cast<std::vector<FRowConstraint> *>(smilpblock->get_static_constraints()[0]);
 (*constraints)[0].set_rhs(20);
 status = solver->compute();

 // Testing FunctionMod on FRowConstraint
 auto frow_f = (*constraints)[0].get_function();
 auto lf1 = dynamic_cast<LinearFunction*>(frow_f);
 lf1->modify_coefficient(0, 60);
 status = solver->compute();

 // Testing VariableMod, fix
 auto x = boost::any_cast<std::vector<ColVariable> *>(smilpblock->get_static_variables()[0]);
 (*x)[0].is_fixed(true);
 status = solver->compute();

 // Testing VariableMod, unfix
 (*x)[0].is_fixed(false);
 status = solver->compute();

 // Testing RowConstraintMod, relax
 (*constraints)[0].relax(true);
 status = solver->compute();

 // Testing RowConstraintMod, enforce
 (*constraints)[0].relax(false);
 status = solver->compute();

 // Testing OneVarConstraint, enforce
 auto bounds = boost::any_cast<std::vector<BoxConstraint> *>(smilpblock->get_static_constraints()[1]);
 (*bounds)[0].set_lhs(1);
 (*bounds)[0].set_rhs(5);
 (*bounds)[1].set_lhs(2);
 (*bounds)[1].set_rhs(6);
 (*bounds)[1].set_both(0);
 status = solver->compute();

 std::cout << "Quitting" << std::endl;
 delete block;

 return 0;
}

