/*--------------------------------------------------------------------------*/
/*---------------------------- File cpx_pars.cpp ---------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Small tool for parsing the macros and the standard std::maps that enable
 * the support for plain CPLEX parameters into CPXMILPSolver.
 *
 * The tool generates two files, CPX<CPX_VERSION>_defs.h and
 * CPX<CPX_VERSION>_maps.h in the specified path.
 * The path should be the include directory of the MILPSolver source tree.
 *
 * \author Niccolo' Iardella \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \copyright &copy; by Niccolo' Iardella
 */
/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include <iostream>
#include <fstream>
#include <string>
#include <cstring>
#include <map>
#include <getopt.h>
// #include <filesystem>

#include <ilcplex/cplex.h>

/*--------------------------------------------------------------------------*/

bool verbose = false;         ///< If the tool should be verbose
// std::filesystem::path path{};
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

int main( int argc , char ** argv )
{

 // Manage options and help
 path = "../include";
 // path = std::filesystem::current_path();
 docopt_desc = "CPLEX parameter map generator.\n";
 exe = get_filename( argv[ 0 ] );
 process_args( argc , argv );

 // if (! std::filesystem::exists(path)) {
 //  std::filesystem::create_directory(path);
 // }

 std::string defs_filename = "CPX" + std::to_string( CPX_VERSION ) + "_defs.h";
 std::string maps_filename = "CPX" + std::to_string( CPX_VERSION ) + "_maps.h";

 // auto defs_path = path / defs_filename;
 // auto maps_path = path / maps_filename;

 auto defs_path = path + "/" + defs_filename;
 auto maps_path = path + "/" + maps_filename;

 std::ofstream defs_file;
 std::ofstream maps_file;

 CPXENVptr env;
 int status;
 char name[CPX_STR_PARAM_MAX];

 std::map< int , std::string > int_parameters;
 std::map< int , std::string > dbl_parameters;
 std::map< int , std::string > str_parameters;

 int int_counter = 0;
 int dbl_counter = 0;
 int str_counter = 0;

 env = CPXopenCPLEX( &status );
 if( verbose ) {
  std::cout << "CPX_VERSION is " << CPX_VERSION << std::endl;
 }

 for( int i = CPX_PARAM_ALL_MIN ; i <= CPX_PARAM_ALL_MAX ; ++i ) {

#if CPX_VERSION < 12090000
  status = CPXgetparamname( env, i, name );
#else
  status = CPXgetparamhiername( env , i , name );
#endif

  if( status == CPXERR_BAD_PARAM_NUM ) {
   continue;
  }

  if( strlen( name ) == 0 ) {
   continue;
  }

  if( status == 0 ) {
   int type;
   CPXgetparamtype( env , i , &type );

   switch( type ) {
    case( CPX_PARAMTYPE_INT ):
    case( CPX_PARAMTYPE_LONG ):
     int_parameters.insert( { int_counter++ , std::string( name ) } );
     break;

    case( CPX_PARAMTYPE_DOUBLE ):
     // Remove unsupported internal parameters
     if( strcmp( name , "CPXPARAM_Internal_cfilemul" ) == 0 ||
         strcmp( name , "CPXPARAM_Internal_rfilemul" ) == 0 ||
         strcmp( name , "CPXPARAM_Internal_singtol" ) == 0 ||
         strcmp( name , "CPX_PARAM_CFILEMUL" ) == 0 ||
         strcmp( name , "CPX_PARAM_RFILEMUL" ) == 0 ||
         strcmp( name , "CPX_PARAM_SINGTOL" ) == 0 ) {
      break;
     }
     dbl_parameters.insert( { dbl_counter++ , std::string( name ) } );
     break;

    case( CPX_PARAMTYPE_STRING ):
     str_parameters.insert( { str_counter++ , std::string( name ) } );
     break;

    default:
     std::cerr << "Unknown type from CPXgetparamtype()" << std::endl;
     return( 1 );
   }

  } else {
   std::cerr << "Unknown error in CPXgetparamhiername()" << std::endl;
   return( 1 );
  }
 }

 // Generate defs file
 defs_file.open( defs_path );

 defs_file << "/* FILE GENERATED AUTOMATICALLY, DO NOT EDIT */" << std::endl
           << std::endl
           << "#ifndef __CPX" << std::to_string( CPX_VERSION ) << "_DEFS"
           << std::endl
           << "#define __CPX" << std::to_string( CPX_VERSION ) << "_DEFS"
           << std::endl << std::endl
           << "#define CPX_NUM_INT_PARS " << int_counter << std::endl
           << "#define CPX_NUM_DBL_PARS " << dbl_counter << std::endl
           << "#define CPX_NUM_STR_PARS " << str_counter << std::endl
           << std::endl
           << "#endif //__CPX" << std::to_string( CPX_VERSION ) << "_DEFS"
           << std::endl;

 defs_file.close();
 std::cout << "Defs file written on " << defs_path << std::endl;

 // Generate maps file
 maps_file.open( maps_path );
 maps_file << "/* FILE GENERATED AUTOMATICALLY, DO NOT EDIT */" << std::endl
           << std::endl
           << "#include <ilcplex/cplex.h>" << std::endl
           << "#include \"CPXMILPSolver.h\"" << std::endl
           << std::endl
           << "using namespace SMSpp_di_unipi_it;" << std::endl
           << std::endl;

 // SMSpp_to_CPLEX_***_pars maps
 maps_file
  << "const std::array< int, CPX_NUM_INT_PARS >"
  << " CPXMILPSolver::SMSpp_to_CPLEX_int_pars{"
  << std::endl;
 for( const auto & i : int_parameters ) {
  maps_file << " " << i.second << "," << std::endl;
 }
 maps_file << "};" << std::endl;
 maps_file << std::endl;

 maps_file
  << "const std::array< int, CPX_NUM_DBL_PARS >"
  << " CPXMILPSolver::SMSpp_to_CPLEX_dbl_pars{"
  << std::endl;
 for( const auto & i : dbl_parameters ) {
  maps_file << " " << i.second << "," << std::endl;
 }
 maps_file << "};" << std::endl;
 maps_file << std::endl;

 maps_file
  << "const std::array< int, CPX_NUM_STR_PARS >"
  << " CPXMILPSolver::SMSpp_to_CPLEX_str_pars{"
  << std::endl;
 for( const auto & i : str_parameters ) {
  maps_file << " " << i.second << "," << std::endl;
 }
 maps_file << "};" << std::endl;
 maps_file << std::endl;

 // Reverse CPLEX_to_SMSpp_***_pars maps
 maps_file
  << "const std::array< std::pair< int, int >, CPX_NUM_INT_PARS >" << std::endl
  << " CPXMILPSolver::CPLEX_to_SMSpp_int_pars{" << std::endl
  << " {" << std::endl;
 for( const auto & i : int_parameters ) {
  maps_file << "  { " << i.second << ", intFirstCPLEXPar + " << i.first << " },"
            << std::endl;
 }
 maps_file
  << " }" << std::endl
  << "};" << std::endl
  << std::endl;

 maps_file
  << "const std::array< std::pair< int, int >, CPX_NUM_DBL_PARS >" << std::endl
  << " CPXMILPSolver::CPLEX_to_SMSpp_dbl_pars{" << std::endl
  << " {" << std::endl;
 for( const auto & i : dbl_parameters ) {
  maps_file << "  { " << i.second << ", dblFirstCPLEXPar + " << i.first << " },"
            << std::endl;
 }
 maps_file
  << " }" << std::endl
  << "};" << std::endl
  << std::endl;

 maps_file
  << "const std::array< std::pair< int, int >, CPX_NUM_STR_PARS >" << std::endl
  << " CPXMILPSolver::CPLEX_to_SMSpp_str_pars{" << std::endl
  << " {" << std::endl;
 for( const auto & i : str_parameters ) {
  maps_file << "  { " << i.second << ", strFirstCPLEXPar + " << i.first << " },"
            << std::endl;
 }
 maps_file
  << " }" << std::endl
  << "};" << std::endl;
 maps_file.close();
 std::cout << "Maps file written on " << maps_path << std::endl;

 return( 0 );
}
