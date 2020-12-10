/** @file
 * SMS++ MILP solver.
 *
 * A tool that loads an SimpleMILPBlock from a .milp file,
 * optionally configures it with a BlockConfig and a BlockSolverConfig,
 * and solves it with all the loaded Solvers.
 *
 * Optionally, it writes back the Block, the BlockConfig and the
 * BlockSolverConfig on a SMS++ nc4 problem file.
 *
 * \author Niccolo' Iardella \n
 *         Operations Research Group \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * Copyright &copy; by Niccolo' Iardella
 */

#include <iostream>
#include <iomanip>

#include <CPXMILPSolver.h>
#include <BlockSolverConfig.h>

#include "solver_utils.h"
// TODO: Remove SimpleMILPBlock
#include "../test/SimpleMILPBlock.h"

using namespace SMSpp_di_unipi_it;

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/


int main(int argc, char** argv) {

 // Manage options and help, see common_utils.h
 docopt_desc = "SMS++ MILP solver.\n";
 exe = get_filename( argv[ 0 ] );
 process_args( argc, argv );

 std::ifstream file(filename);
 if (!file.is_open()) {
  std::cerr << exe << ": cannot open file " << filename << std::endl;
  exit (1);
 }

 // Read block
 auto block = new SimpleMILPBlock();
 file >> *block;
 std::cout << *block;

 // Configure block
 BlockConfig * b_config;
 if( !bconf_file.empty() ) {
  b_config = get_blockconfig( bconf_file );
  if( b_config == nullptr ) {
   std::cerr << exe << ": Block configuration not valid" << std::endl;
   exit( 1 );
  }
  b_config->apply( block );
 }

 // Configure solver
 BlockSolverConfig * s_config;
 if( !sconf_file.empty() ) {
  s_config = get_blocksolverconfig( sconf_file );
  if( s_config == nullptr ) {
   std::cerr << exe << ": Block configuration not valid" << std::endl;
   exit( 1 );
  }
 } else {
  std::cout << "Using a default Solver configuration" << std::endl;
  s_config = default_configure_solver( solvVerbose );
 }
 s_config->apply( block );

 // Write nc4 problem
 if( writeprob ) {
  write_nc4problem( block, b_config, s_config );
 }

 // Solve
 std::cout.setf( std::ios::scientific, std::ios::floatfield );
 std::cout << std::setprecision( 8 );
 solve_all( block );

 // Print solution
 for( auto solver : block->get_registered_solvers() ) {
  auto status = solver->compute();
  if( solver->compute() == Solver::kOK ) {
   solver->get_var_solution();
   std::string s = "[";
   for( auto & i : block->get_x() ) {
    auto x = i.get_value();
    s += " " + std::to_string( x );
   }
   s += " ]";
   std::cout << "Solution = " << s << std::endl;
  }
  if( static_cast<CPXMILPSolver *>(solver)->has_dual_solution() ) {
   static_cast<CPXMILPSolver *>(solver)->get_dual_solution();
   std::string s = "[";

   auto set = [ &s ]( FRowConstraint & c ) {
    auto pi = c.get_dual();
    s += " " + std::to_string( pi );
   };

   for( const auto & i : block->get_static_constraints() ) {
    un_any_const_static( i, set, un_any_type< FRowConstraint >() );
   }

   for( const auto & i : block->get_dynamic_constraints() ) {
    un_any_const_dynamic( i, set, un_any_type< FRowConstraint >() );
   }
   s += " ]";
   std::cout << "Dual solution = " << s << std::endl;
  }
 }
 return 0;
}
