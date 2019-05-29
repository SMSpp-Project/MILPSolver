//
// Created by Niccolò Iardella on 2019-03-07.
//

#include <iostream>
#include <fstream>

#include "SimpleMILPBlock.h"
#include "MILPSolver.h"
#include "LinearFunction.h"

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

 block->register_Solver(Solver::new_Solver("MILPSolver"));

 // First solve
 auto solver = (block->get_registered_solvers()).front();
 int status = solver->compute();

 auto smilpblock = dynamic_cast<SimpleMILPBlock*>(block);
 auto obj = boost::any_cast<FRealObjective*>(smilpblock->get_objective());
 auto obj_f = obj->get_function();

 std::cout << "Status = " << status << std::endl;
 auto vars = smilpblock->get_x();
 for (auto &i : vars) {
  std::cout << "Variable value =  " << i.get_value() << std::endl;
 }
 std::cout << "Function value =  " << obj_f->get_value() << std::endl;

 // // Testing FunctionMod on Objective
 // LinearFunction::v_coeff nc0 = {-8};
 // auto lf0 = dynamic_cast<LinearFunction*>(obj_f);
 // lf0->modify_coefficients(nc0.begin(), 0, 1);
 // status = solver->compute();
 //
 // // Testing ObjectiveMod, maximize
 // obj->set_sense(Objective::eMax);
 // status = solver->compute();
 //
 // // Testing ObjectiveMod, minimize
 // obj->set_sense(Objective::eMin);
 // status = solver->compute();
 //
 // // Testing RowConstraintMod, set RHS
 // auto constraints = boost::any_cast<std::vector<FRowConstraint> *>(smilpblock->get_static_constraints()[0]);
 // (*constraints)[0].set_rhs(20);
 // status = solver->compute();
 //
 // // Testing FunctionMod on FRowConstraint
 // auto frow_f = (*constraints)[0].get_function();
 // auto lf1 = dynamic_cast<LinearFunction*>(frow_f);
 // LinearFunction::v_coeff nc1 = {60};
 // lf1->modify_coefficients(nc1.begin(), 0, 1);
 // status = solver->compute();
 //
 // // Testing VariableMod, fix
 // auto x = boost::any_cast<std::vector<ColVariable> *>(smilpblock->get_static_variables()[0]);
 // (*x)[0].is_fixed(true);
 // status = solver->compute();
 //
 // // Testing VariableMod, unfix
 // (*x)[0].is_fixed(false);
 // status = solver->compute();
 //
 // // Testing RowConstraintMod, relax
 // (*constraints)[0].relax(true);
 // status = solver->compute();
 //
 // // Testing RowConstraintMod, enforce
 // (*constraints)[0].relax(false);
 // status = solver->compute();
 //
 // // Testing OneVarConstraint, enforce
 // auto bounds = boost::any_cast<std::vector<BoxConstraint> *>(smilpblock->get_static_constraints()[1]);
 // (*bounds)[0].set_lhs(1);
 // (*bounds)[0].set_rhs(5);
 // (*bounds)[1].set_lhs(2);
 // (*bounds)[1].set_rhs(6);
 // (*bounds)[1].set_both(0);
 // status = solver->compute();

 // Testing dynamic modifications, add constraint
 auto x = boost::any_cast<std::vector<ColVariable> *>(smilpblock->get_static_variables()[0]);

 auto x1 = &x[0][0];
 auto x2 = &x[0][1];

 LinearFunction::v_coeff_pair v_cp = {{x1, 4}, {x2, 4}};
 auto p_lf = new LinearFunction(std::move(v_cp));


 std::list<FRowConstraint> new_constraints(1);
 new_constraints.front().set_function(p_lf);



 smilpblock->add_dynamic_constraint(new_constraints);

 auto cmod = std::make_shared<BlockModAD>( BlockModAD::eAddConst );
 cmod->mod_list = &(new_constraints);
 smilpblock->add_Modification(cmod);
 status = solver->compute();


 std::cout << "Quitting" << std::endl;
 delete block;

 return 0;
}

