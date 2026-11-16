/*--------------------------------------------------------------------------*/
/*------------------------- File test_fail_feas.cpp ------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Checks the status and the bounds that a :MILPSolver gives when its
 * solver ends a MIP with a solution and no accuracy promise attached to it:
 * the solution is there, and so is its value, hence the status has to be
 * kLowPrecision, as that of a continuous problem whose solution is not
 * proved optimal, and not kError, which says that there is nothing to read;
 * the bound of the search is given where the solver still vouches for it.
 *
 * Both problems are the second-order cone formulation of a
 * delay-constrained routing problem, with binary variables and quadratic
 * rows, written in the LP format by the solver itself:
 *
 * - fail-feas.lp, 110 binary variables and 150 quadratic rows, on which
 *   CPLEX 22.1.1 on one thread (MILPCfg-fail-feas.txt) stops its search
 *   for an error at the root node (CPXMIP_FAIL_FEAS), after the
 *   linearization of the cones, with an incumbent of value 15283.82,
 *   which is the optimum: CPXMILPSolver has to answer kLowPrecision with
 *   that solution (or kOK, with a CPLEX that meets no error);
 *
 * - grb-suboptimal.lp, on which Gurobi 13.0 on 4 threads with its
 *   automatic numerical care (MILPCfg-grb-suboptimal.txt) ends the search
 *   as sub-optimal (GRB_SUBOPTIMAL): GRBMILPSolver has to answer
 *   kLowPrecision (or kOK) with a feasible solution and keep the bound of
 *   the search, which Gurobi still gives there, i.e., get_lb() has to be
 *   finite and not above get_ub(). Only a model stopped on numerical
 *   troubles (GRB_NUMERIC) has its bound dropped, which no small problem
 *   reproduces from a file.
 *
 * A case whose :MILPSolver is not in the build is skipped; if no case is
 * left the test exits with 77, the status ctest reads as "skipped" [see
 * test_common.h].
 *
 * \author Donato Meoli \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 */

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include <cmath>
#include <fstream>
#include <iostream>
#include <memory>
#include <string>

#include "AbstractBlock.h"

#include "CDASolver.h"

#include "BlockSolverConfig.h"

/*--------------------------------------------------------------------------*/

using namespace SMSpp_di_unipi_it;

/*--------------------------------------------------------------------------*/
/*----------------------------- CONSTANTS ----------------------------------*/
/*--------------------------------------------------------------------------*/

/// a case: the :MILPSolver, the problem, the ComputeConfig and the optimal
/// value, NaN if the check is on the bounds rather than on the value

struct Case {
 std::string solver;
 std::string problem;
 std::string config;
 double optimum;
 };

static const Case Cases[] = {
 { "CPXMILPSolver" , "fail-feas.lp" , "MILPCfg-fail-feas.txt" ,
   15283.8234 } ,
 { "GRBMILPSolver" , "grb-suboptimal.lp" , "MILPCfg-grb-suboptimal.txt" ,
   std::nan( "" ) } };

/// the relative gap of the checks, that of the configurations with leeway
static constexpr double RelGap = 1e-5;

/*--------------------------------------------------------------------------*/
/*------------------------------ FUNCTIONS ---------------------------------*/
/*--------------------------------------------------------------------------*/

/// solves the case and checks what the :MILPSolver says; true if it passes

static bool check( const Case & c )
{
 std::cout << c.solver << " on " << c.problem << ": ";

 std::ifstream file( c.problem );
 if( ! file ) {
  std::cout << "cannot open the problem -> KO" << std::endl;
  return( false );
  }

 std::unique_ptr< Configuration > cfg(
			       Configuration::deserialize( c.config ) );
 auto cc = dynamic_cast< ComputeConfig * >( cfg.get() );
 if( ! cc ) {
  std::cout << c.config << " is not a ComputeConfig -> KO" << std::endl;
  return( false );
  }

 auto block = std::make_unique< AbstractBlock >();
 block->load( file , 'L' );

 auto solver = Solver::new_Solver( c.solver );
 solver->set_ComputeConfig( cc );
 block->register_Solver( solver );

 const auto status = solver->compute();
 bool passed = true;

 std::cout << "status " << status;
 if( status == Solver::kError ) {
  std::cout << " -> KO, kError with a solution at hand";
  passed = false;
  }
 else if( ( status == Solver::kOK ) ||
	  ( status == Solver::kLowPrecision ) ) {
  auto cda = dynamic_cast< CDASolver * >( solver );
  const double value = cda->get_var_value();
  const double lb = solver->get_lb();
  const double ub = solver->get_ub();
  std::cout << ( status == Solver::kOK ? " (kOK" : " (kLowPrecision" )
	    << "), value " << value << " in [ " << lb << " , " << ub << " ]";
  const double tol = RelGap * std::max( 1.0 , std::abs( value ) );
  if( ! ( cda->has_var_solution() && cda->is_var_feasible() ) ) {
   std::cout << " -> KO, no feasible solution";
   passed = false;
   }
  else if( ! std::isnan( c.optimum ) ) {
   if( std::abs( value - c.optimum ) > RelGap * std::abs( c.optimum ) ) {
    std::cout << " -> KO, the optimum is " << c.optimum;
    passed = false;
    }
   }
  else if( ! std::isfinite( lb ) ) {
   std::cout << " -> KO, the bound of the search is lost";
   passed = false;
   }
  else if( ( lb > ub + tol ) || ( lb > value + tol ) ) {
   std::cout << " -> KO, the bound is above the solution";
   passed = false;
   }
  if( passed )
   std::cout << " -> OK";
  }
 else {
  std::cout << " -> KO, the problem has an optimal solution";
  passed = false;
  }
 std::cout << std::endl;

 block->unregister_Solver( solver );
 delete solver;

 return( passed );
 }

/*--------------------------------------------------------------------------*/
/*-------------------------------- MAIN ------------------------------------*/
/*--------------------------------------------------------------------------*/

int main( void )
{
 bool passed = true;
 bool any = false;

 for( const auto & c : Cases ) {
  if( ! Solver::has_Solver( c.solver ) ) {
   std::cout << c.solver << " is not in this build, skipped" << std::endl;
   continue;
   }
  any = true;
  passed &= check( c );
  }

 if( ! any ) {
  std::cout << "no :MILPSolver of the cases in this build, nothing to check"
	    << std::endl;
  return( 77 );
  }

 if( passed )
  std::cout << "All tests passed!!" << std::endl;
 else
  std::cout << "Shit happened!!" << std::endl;

 return( passed ? 0 : 1 );
 }

/*--------------------------------------------------------------------------*/
/*------------------------- File test_fail_feas.cpp ------------------------*/
/*--------------------------------------------------------------------------*/
