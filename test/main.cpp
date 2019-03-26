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

 auto solver = (block->get_registered_solvers()).front();
 int status = solver->compute();

 std::cout << "[DEBUG] Status = " << status << std::endl;


 std::cout << "SOLUTION" << std::endl;
 auto smilpblock = dynamic_cast<SimpleMILPBlock*>(block);
 auto vars = smilpblock->get_x();

 for (auto &i : vars) {
  std::cout << "Variable value =  " << i.get_value() << std::endl;
 }

 auto obj = boost::any_cast<FRealObjective*>(smilpblock->get_objective());
 auto fun = obj->get_function();
 std::cout << "Function value =  " << fun->get_value() << std::endl;


 std::cout << "Quitting" << std::endl;
 delete block;

 return 0;
}

