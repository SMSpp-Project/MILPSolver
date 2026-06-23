/*--------------------------------------------------------------------------*/
/*------------------------- File test_pips.cpp -----------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Small standalone driver for exercising PIPSMILPSolver on one nc4 instance.
 *
 * The program:
 * - deserializes a Block from an nc4 file;
 * - applies a BlockSolverConfig, plain or classname-based meta configuration;
 * - runs the first Solver registered on the root Block;
 * - prints status, time, iterations, and objective information;
 * - asks the Solver to write the primal solution back into the Block variables;
 * - dumps all ColVariable values in the Block tree to a text file;
 * - optionally dumps all FRowConstraint dual values in the Block tree to a text file.
 */
/*--------------------------------------------------------------------------*/

#include <chrono>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <list>
#include <map>
#include <memory>
#include <stdexcept>
#include <string>

#include <Block.h>
#include <BlockInspection.h>
#include <BlockSolverConfig.h>
#include <CDASolver.h>
#include <ColVariable.h>
#include <FRowConstraint.h>
#include <Configuration.h>
#include <Solver.h>

using namespace SMSpp_di_unipi_it;

/*--------------------------------------------------------------------------*/

namespace {

using Clock = std::chrono::steady_clock;

/*--------------------------------------------------------------------------*/

std::string status_to_string( int status )
{
 switch( status ) {
  case Solver::kOK:             return( "kOK" );
  case Solver::kLowPrecision:   return( "kLowPrecision" );
  case Solver::kInfeasible:     return( "kInfeasible" );
  case Solver::kUnbounded:      return( "kUnbounded" );
  case Solver::kStopTime:       return( "kStopTime" );
  case Solver::kStopIter:       return( "kStopIter" );
  case Solver::kError:          return( "kError" );
  case Solver::kUnEval:         return( "kUnEval" );
  default:                      return( "status_" + std::to_string( status ) );
  }
}

/*--------------------------------------------------------------------------*/

bool has_primal_solution( int status )
{
 return( ( ( status >= Solver::kOK ) && ( status < Solver::kError ) &&
           ( status != Solver::kUnbounded ) &&
           ( status != Solver::kInfeasible ) ) ||
         ( status == Solver::kLowPrecision ) );
}

/*--------------------------------------------------------------------------*/

void apply_solver_config( Block * block , Configuration * config ,
                          const std::string & filename )
{
 if( auto * meta =
     dynamic_cast< SimpleConfiguration< std::map< std::string ,
                                                  Configuration * > > * >(
      config ) ) {
  std::list< Block * > bfs;
  bfs.push_back( block );

  for( auto it = bfs.begin() ; it != bfs.end() ; ++it )
   for( auto nested : ( *it )->get_nested_Blocks() )
    bfs.push_back( nested );

  auto & map = meta->f_value;
  for( auto b : bfs ) {
   const auto cfg_it = map.find( b->classname() );
   if( cfg_it == map.end() )
    continue;

   auto * bsc = dynamic_cast< BlockSolverConfig * >( cfg_it->second );
   if( ! bsc )
    throw( std::runtime_error(
     "meta configuration entry for " + cfg_it->first +
     " in " + filename + " is not a BlockSolverConfig" ) );

   bsc->apply( b );
   }

  for( auto & entry : map )
   entry.second->clear();

  return;
  }

 auto * bsc = dynamic_cast< BlockSolverConfig * >( config );
 if( ! bsc )
  throw( std::runtime_error(
   filename + " does not contain a valid [meta]BlockSolverConfig" ) );

 bsc->apply( block );
 bsc->clear();
}

/*--------------------------------------------------------------------------*/

void dump_block_variables( const Block * block , const std::string & path ,
                           std::ofstream & out , std::size_t & count )
{
 const auto & static_vars = block->get_static_variables();
 for( Block::Index group = 0 ; group < static_vars.size() ; ++group ) {
  const auto size =
   inspection::get_element_size< ColVariable >( block , true , group );

  for( Block::Index index = 0 ; index < size ; ++index ) {
   auto * var =
    inspection::get_element< ColVariable >( block , true , group , index );
   if( ! var )
    continue;

   out << path << "\tstatic\t" << group << "\t" << index << "\t"
       << static_cast< const void * >( var ) << "\t"
       << std::setprecision( 17 ) << var->get_value() << "\n";
   ++count;
   }
  }

 const auto & dynamic_vars = block->get_dynamic_variables();
 for( Block::Index group = 0 ; group < dynamic_vars.size() ; ++group ) {
  const auto size =
   inspection::get_element_size< ColVariable >( block , false , group );

  for( Block::Index index = 0 ; index < size ; ++index ) {
   auto * var =
    inspection::get_element< ColVariable >( block , false , group , index );
   if( ! var )
    continue;

   out << path << "\tdynamic\t" << group << "\t" << index << "\t"
       << static_cast< const void * >( var ) << "\t"
       << std::setprecision( 17 ) << var->get_value() << "\n";
   ++count;
   }
  }

 const auto & nested = block->get_nested_Blocks();
 for( Block::Index index = 0 ; index < nested.size() ; ++index ) {
  if( ! nested[ index ] )
   continue;

  dump_block_variables( nested[ index ] ,
                        path + "/" + std::to_string( index ) + ":" +
                         nested[ index ]->classname() ,
                        out , count );
  }
}

/*--------------------------------------------------------------------------*/

std::size_t dump_primal_solution( const Block * block ,
                                  const std::string & solution_file )
{
 std::ofstream out( solution_file );
 if( ! out )
  throw( std::runtime_error( "cannot open solution file " + solution_file ) );

 out << "# block_path\tkind\tgroup\tindex\taddress\tvalue\n";

 std::size_t count = 0;
 dump_block_variables( block , "root:" + block->classname() , out , count );
 return( count );
}

/*--------------------------------------------------------------------------*/

void dump_block_constraints( const Block * block , const std::string & path ,
                             std::ofstream & out , std::size_t & count )
{
 const auto & static_cons = block->get_static_constraints();
 for( Block::Index group = 0 ; group < static_cons.size() ; ++group ) {
  const auto size =
   inspection::get_element_size< FRowConstraint >( block , true , group );
  if( size == Inf< Block::Index >() )
   continue;

  for( Block::Index index = 0 ; index < size ; ++index ) {
   auto * con =
    inspection::get_element< FRowConstraint >( block , true , group , index );
   if( ! con )
    continue;

   out << path << "\tstatic\t" << group << "\t" << index << "\t"
       << static_cast< const void * >( con ) << "\t"
       << std::setprecision( 17 ) << con->get_dual() << "\n";
   ++count;
   }
  }

 const auto & dynamic_cons = block->get_dynamic_constraints();
 for( Block::Index group = 0 ; group < dynamic_cons.size() ; ++group ) {
  const auto size =
   inspection::get_element_size< FRowConstraint >( block , false , group );
  if( size == Inf< Block::Index >() )
   continue;

  for( Block::Index index = 0 ; index < size ; ++index ) {
   auto * con =
    inspection::get_element< FRowConstraint >( block , false , group , index );
   if( ! con )
    continue;

   out << path << "\tdynamic\t" << group << "\t" << index << "\t"
       << static_cast< const void * >( con ) << "\t"
       << std::setprecision( 17 ) << con->get_dual() << "\n";
   ++count;
   }
  }

 const auto & nested = block->get_nested_Blocks();
 for( Block::Index index = 0 ; index < nested.size() ; ++index ) {
  if( ! nested[ index ] )
   continue;

  dump_block_constraints( nested[ index ] ,
                          path + "/" + std::to_string( index ) + ":" +
                           nested[ index ]->classname() ,
                          out , count );
  }
}

/*--------------------------------------------------------------------------*/

std::size_t dump_dual_solution( const Block * block ,
                                const std::string & dual_solution_file )
{
 std::ofstream out( dual_solution_file );
 if( ! out )
  throw( std::runtime_error( "cannot open dual solution file " +
                             dual_solution_file ) );

 out << "# block_path\tkind\tgroup\tindex\taddress\tdual\n";

 std::size_t count = 0;
 dump_block_constraints( block , "root:" + block->classname() , out , count );
 return( count );
}

/*--------------------------------------------------------------------------*/

void usage( const char * exe )
{
 std::cerr
  << "usage: " << exe
  << " [instance.nc4] [BlockSolverConfig.txt] [primal.txt] [dual.txt]\n\n"
  << "defaults:\n"
  << "  instance.nc4           TSSB_EC_CO_Test.nc4\n"
  << "  BlockSolverConfig.txt  BSCfg1.txt\n"
  << "  primal.txt             primal_solution.txt\n"
  << "  dual.txt               dual_solution.txt\n";
}

/*--------------------------------------------------------------------------*/

}  // namespace

/*--------------------------------------------------------------------------*/

int main( int argc , char ** argv )
{
 try {
  if( argc > 1 && ( std::string( argv[ 1 ] ) == "-h" ||
                    std::string( argv[ 1 ] ) == "--help" ) ) {
   usage( argv[ 0 ] );
   return( 0 );
   }

  const std::string instance_file = argc > 1 ? argv[ 1 ] : "EC_CO_Test.nc4";
  const std::string config_file = argc > 2 ? argv[ 2 ] : "BSCfg1.txt";
  const std::string solution_file =
   argc > 3 ? argv[ 3 ] : "primal_solution.txt";
  const std::string dual_solution_file =
   argc > 4 ? argv[ 4 ] : "dual_solution.txt";

  std::cout << "Instance:      " << instance_file << "\n"
            << "Solver config: " << config_file << "\n"
            << "Primal file:   " << solution_file << "\n"
            << "Dual file:     " << dual_solution_file << "\n";

  std::unique_ptr< Block > block( Block::deserialize( instance_file ) );
  if( ! block )
   throw( std::runtime_error( "Block::deserialize failed for " +
                              instance_file ) );

  std::unique_ptr< Configuration > config(
   Configuration::deserialize( config_file ) );
  if( ! config )
   throw( std::runtime_error( "Configuration::deserialize failed for " +
                              config_file ) );

  apply_solver_config( block.get() , config.get() , config_file );

  const auto & solvers = block->get_registered_solvers();
  if( solvers.empty() )
   throw( std::runtime_error( "no Solver registered on the root Block" ) );

  Solver * solver = solvers.front();
  solver->set_log( &std::cout );

  std::cout << "Block class:    " << block->classname() << "\n"
            << "Solver class:   " << solver->classname() << "\n"
            << "Solving...\n";

  const auto start = Clock::now();
  const int status = solver->compute( false );
  const auto stop = Clock::now();
  const std::chrono::duration< double > elapsed = stop - start;

  std::cout << std::setprecision( 10 )
            << "Status:        " << status << " ("
            << status_to_string( status ) << ")\n"
            << "Time_s:        " << elapsed.count() << "\n"
            << "Iterations:    " << solver->get_elapsed_iterations() << "\n"
            << "Lower_bound:   " << solver->get_lb() << "\n"
            << "Upper_bound:   " << solver->get_ub() << "\n";

  if( has_primal_solution( status ) ) {
   std::cout << "Objective:     " << solver->get_var_value() << "\n";
   solver->get_var_solution();
   const auto n_values = dump_primal_solution( block.get() , solution_file );
   std::cout << "Wrote " << n_values << " primal variable values to "
             << solution_file << "\n";

   auto * cda_solver = dynamic_cast< CDASolver * >( solver );
   if( cda_solver && cda_solver->has_dual_solution() ) {
    cda_solver->get_dual_solution();
    const auto n_duals =
     dump_dual_solution( block.get() , dual_solution_file );
    std::cout << "Wrote " << n_duals << " row dual values to "
              << dual_solution_file << "\n";
    }
   else {
    std::cout << "No dual solution available; dual solution file not written.\n";
    }
   }
  else {
   std::cout << "No primal solution available; solution files not written.\n";
   }

  block->unregister_Solvers( true );
  return( has_primal_solution( status ) ? 0 : 1 );
  }
 catch( const std::exception & e ) {
  std::cerr << "error: " << e.what() << std::endl;
  return( 1 );
  }
}

/*--------------------------------------------------------------------------*/
/*------------------------- End File test_pips.cpp -------------------------*/
