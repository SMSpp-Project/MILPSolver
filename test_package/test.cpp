#include "CPXMILPSolver.h"

using namespace SMSpp_di_unipi_it;

int main(int argc, char** argv) {

 Solver* solver = new CPXMILPSolver();
 delete solver;

 return( 0 );
}

