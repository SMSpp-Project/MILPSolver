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

 // Load problem in the block
 auto block = Block::new_Block("SimpleMILPBlock");
 file >> *block;
 std::cout << *block;

 // Register MILP Solver
 block->register_Solver(Solver::new_Solver("MILPSolver"));

 // First solve
 auto solver = (block->get_registered_solvers()).front();
 int status = solver->compute();

 // Print results
 std::cout << "[DEBUG] Status = " << status << std::endl;
 std::cout << "SOLUTION" << std::endl;
 auto smilpblock = dynamic_cast<SimpleMILPBlock*>(block);
 auto vars = smilpblock->get_x();
 for (auto &i : vars) {
  std::cout << "Variable value =  " << i.get_value() << std::endl;
 }
 auto obj = boost::any_cast<FRealObjective*>(smilpblock->get_objective());
 auto obj_f = obj->get_function();
 std::cout << "Function value =  " << obj_f->get_value() << std::endl;

 // Modification test: Change OF
 LinearFunction::v_coeff nc = {-8};
 auto lf = dynamic_cast<LinearFunction*>(obj_f);
 lf->modify_coefficients(nc.begin(), 0, 1);
 status = solver->compute();
 std::cout << "[DEBUG] Status = " << status << std::endl;
 std::cout << "SOLUTION" << std::endl;
 for (auto &i : vars) {
  std::cout << "Variable value =  " << i.get_value() << std::endl;
 }
 std::cout << "Function value =  " << obj_f->get_value() << std::endl;

 // Modification test: Change a constraint
 auto constraints = boost::any_cast<std::vector<FRowConstraint> *>((smilpblock->get_static_constraints())[0]);
 (*constraints)[0].set_rhs(20);
 status = solver->compute();
 std::cout << "[DEBUG] Status = " << status << std::endl;
 std::cout << "SOLUTION" << std::endl;
 for (auto &i : vars) {
  std::cout << "Variable value =  " << i.get_value() << std::endl;
 }
 std::cout << "Function value =  " << obj_f->get_value() << std::endl;




 std::cout << "Quitting" << std::endl;
 delete block;

 return 0;
}

