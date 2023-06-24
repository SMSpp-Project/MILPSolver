/** @file
 * Some common utilities for SMS++ tools.
 *
 * \author Niccolo' Iardella \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \copyright &copy; by Niccolo' Iardella
 */

#ifndef __COMMON_UTILS
#define __COMMON_UTILS

#include <getopt.h>

#include <Block.h>
#include <RBlockConfig.h>

using namespace SMSpp_di_unipi_it;

/**
 * @name Global variables used by every tool
 * @{
 */

std::string filename{};    ///< Input filename
std::string bconf_file{};  ///< BlockConfig filename
std::string sconf_file{};  ///< BlockSolverConfig filename
bool solvVerbose = false;  ///< If the solver should be verbose
bool writeprob = false;    ///< If the problem should be written back
std::string exe{};         ///< Name of the executable file
std::string docopt_desc{}; ///< Tool description
/// @}

/*--------------------------------------------------------------------------*/

/// Gets the name of the executable from its full path
std::string get_filename( const std::string & fullpath ) {
 std::size_t found = fullpath.find_last_of( "/\\" );
 return( fullpath.substr( found + 1 ) );
}

/*--------------------------------------------------------------------------*/

/// Prints the tool description and usage
void docopt() {
 // http://docopt.org
 std::cout << docopt_desc << std::endl;
 std::cout << "Usage:\n"
           << "  " << exe << " [options] <file>\n"
           << "  " << exe << " -h | --help\n"
           << std::endl
           << "Options:\n"
           << "  -B, --blockcfg <file>    Block configuration.\n"
           << "  -S, --solvercfg <file>   Solver configuration.\n"
           << "  -n, --nc4problem <file>  Write nc4 problem on file.\n"
           << "  -v, --verbose            Make the solver verbose.\n"
           << "  -h, --help               Print this help.\n";
}

/*--------------------------------------------------------------------------*/

/// Processes the command line arguments
void process_args( int argc, char ** argv ) {

 if( argc < 2 ) {
  std::cout << exe << ": no input file\n"
            << "Try " << exe << "' --help' for more information.\n";
  exit( 1 );
 }

 const char * const short_opts = "B:S:nvh";
 const option long_opts[] = {
  { "blockcfg",   required_argument, nullptr, 'B' },
  { "solvercfg",  required_argument, nullptr, 'S' },
  { "nc4problem", no_argument,       nullptr, 'n' },
  { "verbose",    no_argument,       nullptr, 'v' },
  { "help",       no_argument,       nullptr, 'h' },
  { nullptr,      no_argument,       nullptr, 0 }
 };

 // Options
 while( true ) {
  const auto opt = getopt_long( argc, argv, short_opts, long_opts, nullptr );

  if( -1 == opt ) {
   break;
  }
  switch( opt ) {
   case 'B':
    bconf_file = std::string( optarg );
    break;
   case 'S':
    sconf_file = std::string( optarg );
    break;
   case 'n':
    writeprob = true;
    break;
   case 'v':
    solvVerbose = true;
    break;
   case 'h':
    docopt();
    exit( 0 );
   case '?':
   default:
    std::cout << "Try " << exe << "' --help' for more information.\n";
    exit( 1 );
  }
 }

 // Last argument
 if( optind < argc ) {
  filename = std::string( argv[ optind ] );
 } else {
  std::cout << exe << ": no input file\n"
            << "Try " << exe << "' --help' for more information.\n";
  exit( 1 );
 }
}

/*--------------------------------------------------------------------------*/

/// Returns a default Solver configuration
BlockSolverConfig * default_configure_solver( int verbose ) {
 auto s_config = new BlockSolverConfig;
 auto c_config = new ComputeConfig;

 if( verbose )
  c_config->set_par( "intLogVerb" , 1 );

 s_config->add_ComputeConfig( "CPXMILPSolver" , c_config );
 return( s_config );
}

/*--------------------------------------------------------------------------*/

/// Prints the status in a human-readable form
void print_status( int status ) {
 std::cout << "Status = " << status << " (";

 switch( status ) {
  case Solver::kOK:
   std::cout << "Success)" << std::endl;
   break;
  case Solver::kError:
   std::cout << "Error)" << std::endl;
   break;
  case Solver::kInfeasible:
   std::cout << "Infeasible)" << std::endl;
   break;
  case Solver::kUnbounded:
   std::cout << "Unbounded)" << std::endl;
   break;
  case Solver::kStopTime:
   std::cout << "Stopped for time limit)" << std::endl;
   break;
  case Solver::kStopIter:
   std::cout << "Stopped for iteration limit)" << std::endl;
   break;
  default:;
 }
}

/*--------------------------------------------------------------------------*/

/// Solves the problem with all available solvers
void solve_all( Block * block ) {
 for( auto solver : block->get_registered_solvers() ) {
  std::cout << "Solver: " << solver->classname() << std::endl;
  auto status = solver->compute();
  auto ub = solver->get_ub();
  auto lb = solver->get_lb();
  print_status( status );
  std::cout << "Upper bound = " << ub << std::endl;
  std::cout << "Lower bound = " << lb << std::endl;
 }
}

/*--------------------------------------------------------------------------*/

/// Gets a BlockConfig from a BlockConfig file
BlockConfig * get_blockconfig( const std::string & conf_file ) {
 BlockConfig * b_config = nullptr;
 std::ifstream bcf;

 bcf.open( conf_file, std::ifstream::in );
 if( ! bcf.is_open() ) {
  return( nullptr );
 }

 std::string name;
 bcf >> eatcomments >> name;
 b_config = dynamic_cast< BlockConfig * > ( Configuration::new_Configuration( name ) );

 if( ! b_config ) {
  return( nullptr );
 }

 try {
  bcf >> *b_config;
 } catch( const std::exception & e ) {
  return( nullptr );
 }

 return( b_config );
}

/*--------------------------------------------------------------------------*/

/// Gets a BlockSolverConfig from a BlockSolverConfig file
BlockSolverConfig *
get_blocksolverconfig( const std::string & conf_file ) {
 BlockSolverConfig * s_config = nullptr;
 std::ifstream scf;

 scf.open( conf_file, std::ifstream::in );
 if( ! scf.is_open() ) {
  return( nullptr );
 }

 std::string name;
 scf >> eatcomments >> name;
 s_config = dynamic_cast< BlockSolverConfig * > ( Configuration::new_Configuration( name ) );

 if( ! s_config ) {
  return( nullptr );
 }

 try {
  scf >> *s_config;
 } catch( const std::exception & e ) {
  return( nullptr );
 }

 return( s_config );
}

/*--------------------------------------------------------------------------*/

/// Writes a new nc4 problem using the block and its configurations
void write_nc4problem( Block * block,
                       BlockConfig * b_config,
                       BlockSolverConfig * s_config ) {

 std::size_t found = filename.find_last_of( '.' );
 std::string nc4_file = filename.substr( 0, found ) + "_problem.nc4";

 netCDF::NcFile outfile;
 try {
  outfile.open( nc4_file, netCDF::NcFile::replace );
 } catch( netCDF::exceptions::NcException & e ) {
  std::cerr << exe << ": cannot open nc4 file " << nc4_file << std::endl;
  exit( 1 );
 }
 outfile.putAtt( "SMS++_file_type", netCDF::NcInt(), eProbFile );

 block->Block::serialize( outfile, eProbFile );
 netCDF::NcGroup g = outfile.getGroup( "Prob_0" );

 auto new_bc = g.addGroup( "BlockConfig" );
 if( b_config ) {
  b_config->serialize( new_bc );
 }

 auto new_bsc = g.addGroup( "BlockSolver" );
 if( s_config ) {
  s_config->serialize( new_bsc );
 }
 outfile.close();
}

/*--------------------------------------------------------------------------*/


#endif //__COMMON_UTILS
