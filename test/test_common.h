/*--------------------------------------------------------------------------*/
/*------------------------- File test_common.h -----------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * What the testers of this module share: a configuration of this directory
 * names the four `:MILPSolver` backends, whichever of them a build has, and
 * keep_available_Solvers() takes out of it the ones that are not there before
 * it is applied, so that a tester runs with what the machine offers rather
 * than with the one whose name happens to be first in the file.
 *
 * A tester needing more backends than the build has cannot compare anything:
 * it has nothing to do rather than something to fail, and exits with 77, the
 * status ctest reads as "skipped" [see SKIP_RETURN_CODE in CMakeLists.txt].
 *
 * \author Donato Meoli \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \copyright &copy; by Donato Meoli
 */
/*--------------------------------------------------------------------------*/
/*----------------------------- DEFINITIONS --------------------------------*/
/*--------------------------------------------------------------------------*/

#ifndef __MILPSolver_test_common
 #define __MILPSolver_test_common

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "BlockSolverConfig.h"
#include "Solver.h"

#include <cstdlib>
#include <iostream>
#include <string>

/*--------------------------------------------------------------------------*/
/*----------------------------- FUNCTIONS ----------------------------------*/
/*--------------------------------------------------------------------------*/
/// keeps of \p bsc the first \p wanted Solver this build has
/** Takes out of \p bsc the Solver whose name is not in the factory, keeps the
 * first \p wanted of those that are left, each with the ComputeConfig that
 * goes with it, and returns how many there are. Since the ComputeConfig are
 * positional, the entry of a Solver that goes has to go with it, which is why
 * remove_ComputeConfig() is what takes a name out.
 *
 * If fewer than \p wanted are left the run is given up with the status 77,
 * \p fn being said in the message so that the log tells which configuration
 * asked for what the machine does not have. */

inline SMSpp_di_unipi_it::Block::Index keep_available_Solvers(
			       SMSpp_di_unipi_it::BlockSolverConfig * bsc ,
			       const std::string & fn ,
			       SMSpp_di_unipi_it::Block::Index wanted = 1 )
{
 using namespace SMSpp_di_unipi_it;

 if( ! bsc ) {
  std::cerr << "Error: " << fn << " is not a BlockSolverConfig" << std::endl;
  exit( 1 );
  }

 const auto asked = bsc->num_ComputeConfig();

 for( Block::Index i = asked ; i-- ; )
  if( ! Solver::has_Solver( bsc->get_SolverName( i ) ) )
   bsc->remove_ComputeConfig( i );

 // what is left beyond the ones that are wanted is not attached: a Solver
 // that nobody computes would only cost the time of building its matrix
 while( bsc->num_ComputeConfig() > wanted )
  bsc->remove_ComputeConfig( bsc->num_ComputeConfig() - 1 );

 const auto left = bsc->num_ComputeConfig();

 if( left < wanted ) {
  std::cerr << "this build has " << left << " of the " << asked
	    << " :MILPSolver that " << fn << " names, and this test needs "
	    << wanted << ": nothing to do" << std::endl;
  exit( 77 );
  }

 for( Block::Index i = 0 ; i < left ; ++i )
  std::cout << "using " << bsc->get_SolverName( i ) << std::endl;

 return( left );
 }

/*--------------------------------------------------------------------------*/

#endif  /* __MILPSolver_test_common */

/*--------------------------------------------------------------------------*/
/*--------------------- End File test_common.h -----------------------------*/
/*--------------------------------------------------------------------------*/
