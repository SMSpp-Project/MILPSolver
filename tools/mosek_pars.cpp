/*--------------------------------------------------------------------------*/
/*--------------------------- File mosek_pars.cpp --------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Small tool for parsing the macros and the standard std::maps that enable
 * the support for plain MOSEK parameters into MOSEKMILPSolver.
 *
 * The tool generates two files, MOSEK<MOSEK_VERSION>_defs.h and
 * MOSEK<MOSEK_VERSION>_maps.h in the specified path, MOSEK_VERSION being
 * the major, minor and revision numbers of MOSEK one after the other.
 * The path should be the include directory of the MILPSolver source tree.
 * The parameters are listed by MSK_getparamname() on a task, which needs
 * no license.
 *
 * \author Donato Meoli \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \copyright &copy; by Donato Meoli
 */
/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include <iostream>
#include <fstream>
#include <string>
#include <cstring>
#include <map>
#include <vector>
#include <getopt.h>

#include <mosek.h>

/*--------------------------------------------------------------------------*/

bool verbose = false;          ///< If the tool should be verbose
std::string path{};           ///< Path for output files
std::string exe{};            ///< Name of the executable file
std::string docopt_desc{};    ///< Tool description

/*--------------------------------------------------------------------------*/

/// Gets the name of the executable from its full path
std::string get_filename( const std::string & fullpath )
{
 std::size_t found = fullpath.find_last_of( "/\\" );
 return( fullpath.substr( found + 1 ) );
}

/*--------------------------------------------------------------------------*/

/// Prints the tool description and usage
void docopt()
{
 // http://docopt.org
 std::cout << docopt_desc << std::endl;
 std::cout << "Usage:\n"
           << "  " << exe << " [-v] <path>\n"
           << "  " << exe << " -h | --help\n"
           << std::endl
           << "Options:\n"
           << "  -v, --verbose  Make the tool verbose.\n"
           << "  -h, --help     Print this help.\n";
}

/*--------------------------------------------------------------------------*/

/// Processes the command line arguments
void process_args( int argc , char ** argv )
{

 const char * const short_opts = "vh";
 const option long_opts[] = {
  { "verbose" , no_argument , nullptr , 'v' } ,
  { "help" ,    no_argument , nullptr , 'h' } ,
  { nullptr ,   no_argument , nullptr , 0 }
 };

 // Options
 while( true ) {
  const auto opt = getopt_long( argc , argv , short_opts , long_opts ,
                                nullptr );

  if( -1 == opt ) {
   break;
  }
  switch( opt ) {
   case( 'v' ):
    verbose = true;
    break;
   case( 'h' ):
    docopt();
    exit( 0 );
   case( '?' ):
   default:
    std::cout << "Try " << exe << "' --help' for more information.\n";
    exit( 1 );
  }
 }

 // Last argument
 if( optind < argc ) {
  path = std::string( argv[ optind ] );
 }
}

/*--------------------------------------------------------------------------*/

/// Custom terminate function to print the exception message

void smspp_terminate( void ) {
 std::cerr << "Uncaught exception in executing SMS++:\n";
 try {
  std::rethrow_exception( std::current_exception() );
 }
 catch( const std::exception & e ) {
  std::cerr << "\tException type: " << typeid( e ).name() << "\n";
  std::cerr << "\tException message: " << e.what() << "\n";
 } catch( ... ) {
  std::cerr << "\tUnknown exception" << std::endl;
 }
 std::abort(); // or exit(1)
}

/*--------------------------------------------------------------------------*/

int main( int argc , char ** argv )
{
 // override the default terminate handler to print the exception message
 std::set_terminate( smspp_terminate );

 // Manage options and help
 path = "../include";
 docopt_desc = "MOSEK parameter map generator.\n";
 exe = get_filename( argv[ 0 ] );
 process_args( argc , argv );

 std::string MOSEK_VERSION =
  std::to_string( MSK_VERSION_MAJOR ) + std::to_string( MSK_VERSION_MINOR ) +
  std::to_string( MSK_VERSION_REVISION );

 auto defs_path = path + "/MOSEK" + MOSEK_VERSION + "_defs.h";
 auto maps_path = path + "/MOSEK" + MOSEK_VERSION + "_maps.h";

 if( verbose )
  std::cout << "MOSEK_VERSION is " << MOSEK_VERSION << std::endl;

 MSKenv_t env = nullptr;
 MSKtask_t task = nullptr;
 if( ( MSK_makeenv( & env , nullptr ) != MSK_RES_OK ) ||
     ( MSK_maketask( env , 0 , 0 , & task ) != MSK_RES_OK ) ) {
  std::cerr << "Unable to create a MOSEK environment and task" << std::endl;
  return( 1 );
  }

 // the parameters of each type, in the order of their enum
 auto list = [ & ]( MSKparametertypee type , int begin , int end ) {
  std::vector< std::string > names;
  char name[ MSK_MAX_STR_LEN ];
  for( int i = begin ; i < end ; ++i )
   if( MSK_getparamname( task , type , i , name ) == MSK_RES_OK )
    names.emplace_back( name );
  return( names );
  };

 const auto int_parameters = list( MSK_PAR_INT_TYPE , MSK_IPAR_BEGIN ,
                                   MSK_IPAR_END );
 const auto dbl_parameters = list( MSK_PAR_DOU_TYPE , MSK_DPAR_BEGIN ,
                                   MSK_DPAR_END );
 const auto str_parameters = list( MSK_PAR_STR_TYPE , MSK_SPAR_BEGIN ,
                                   MSK_SPAR_END );

 MSK_deletetask( & task );
 MSK_deleteenv( & env );

 // Generate defs file
 std::ofstream defs_file( defs_path );

 defs_file << "/* FILE GENERATED AUTOMATICALLY, DO NOT EDIT */" << std::endl
           << std::endl
           << "#ifndef __MOSEK" << MOSEK_VERSION << "_DEFS" << std::endl
           << "#define __MOSEK" << MOSEK_VERSION << "_DEFS" << std::endl
           << std::endl
           << "#define MOSEK_NUM_INT_PARS " << int_parameters.size()
           << std::endl
           << "#define MOSEK_NUM_DBL_PARS " << dbl_parameters.size()
           << std::endl
           << "#define MOSEK_NUM_STR_PARS " << str_parameters.size()
           << std::endl << std::endl
           << "#endif //__MOSEK" << MOSEK_VERSION << "_DEFS" << std::endl;

 defs_file.close();
 std::cout << "Defs file written on " << defs_path << std::endl;

 // Generate maps file
 std::ofstream maps_file( maps_path );

 maps_file << "/* FILE GENERATED AUTOMATICALLY, DO NOT EDIT */" << std::endl
           << std::endl
           << "#include <mosek.h>" << std::endl
           << "#include \"MOSEKMILPSolver.h\"" << std::endl
           << std::endl
           << "using namespace SMSpp_di_unipi_it;" << std::endl
           << std::endl;

 // SMSpp_to_MOSEK_***_pars and MOSEK_to_SMSpp_***_pars maps
 auto write = [ & ]( const std::vector< std::string > & names ,
                     const std::string & type , const std::string & num ,
                     const std::string & first ) {
  maps_file << "const std::array< std::string, " << num << " >"
            << std::endl
            << " MOSEKMILPSolver::SMSpp_to_MOSEK_" << type << "_pars{"
            << std::endl;
  for( const auto & n : names )
   maps_file << " \"" << n << "\"," << std::endl;
  maps_file << "};" << std::endl << std::endl;

  maps_file << "const std::array< std::pair< std::string, int >, " << num
            << " >" << std::endl
            << " MOSEKMILPSolver::MOSEK_to_SMSpp_" << type << "_pars{"
            << std::endl << " {" << std::endl;
  for( std::size_t i = 0 ; i < names.size() ; ++i )
   maps_file << "  { \"" << names[ i ] << "\", " << first << " + " << i
             << " }," << std::endl;
  maps_file << " }" << std::endl << "};" << std::endl << std::endl;
  };

 write( int_parameters , "int" , "MOSEK_NUM_INT_PARS" , "intFirstMOSEKPar" );
 write( dbl_parameters , "dbl" , "MOSEK_NUM_DBL_PARS" , "dblFirstMOSEKPar" );
 write( str_parameters , "str" , "MOSEK_NUM_STR_PARS" , "strFirstMOSEKPar" );

 maps_file.close();
 std::cout << "Maps file written on " << maps_path << std::endl;

 return( 0 );
}

/*--------------------------------------------------------------------------*/
/*------------------------ End File mosek_pars.cpp -------------------------*/
/*--------------------------------------------------------------------------*/
