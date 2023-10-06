/*--------------------------------------------------------------------------*/
/*---------------------------- File highs_pars.cpp ---------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Small tool for parsing the macros and the standard std::maps that enable
 * the support for plain HiGHS options into HiGHSMILPSolver.
 *
 * The tool generates two files, HiGHS<HiGHS_VERSION>_defs.h and
 * HiGHS<HiGHS_VERSION>_maps.h in the specified path.
 * The path should be the include directory of the MILPSolver source tree.
 *
 * \author Niccolo' Iardella \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \author Enrico Calandrini \n
 *         Dipartimento di Matematica \n
 *         Universita' di Pisa \n
 *
 * Copyright &copy; Niccolo' Iardella
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

#include <interfaces/highs_c_api.h>

/*--------------------------------------------------------------------------*/

bool verbose = true;          ///< If the tool should be verbose
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
 docopt_desc = "HiGHS options map generator.\n";
 exe = get_filename( argv[ 0 ] );
 process_args( argc , argv );

 // if( ! std::filesystem::exists( path ) ) {
 //  std::filesystem::create_directory( path );
 // }

 std::string HiGHS_VERSION =
  std::to_string( HIGHS_VERSION_MAJOR ) + std::to_string( HIGHS_VERSION_MINOR ) +
  std::to_string( HIGHS_VERSION_PATCH );

 std::string defs_filename = "HiGHS" + HiGHS_VERSION + "_defs.h";
 std::string maps_filename = "HiGHS" + HiGHS_VERSION + "_maps.h";

 // auto defs_path = path / defs_filename;
 // auto maps_path = path / maps_filename;

 auto defs_path = path + "/" + defs_filename;
 auto maps_path = path + "/" + maps_filename;

 std::ofstream defs_file;
 std::ofstream maps_file;

 // Create a Highs instance
 void* highsenv = Highs_create();
 int status;
 char * name;

 std::map< int , std::string > int_options;
 std::map< int , std::string > dbl_options;
 std::map< int , std::string > str_options;

 int int_counter = 0;
 int dbl_counter = 0;
 int str_counter = 0;

 if( verbose ) {
  std::cout << "HiGHS_VERSION is " << HiGHS_VERSION << std::endl;
 }

 for( int i = 0 ; Highs_getOptionName( highsenv , i , &name ) == kHighsStatusOk ;
             ++i ) {

  if( strlen( name ) == 0 ) {
   break;
  }

  int type;
  Highs_getOptionType( highsenv , name , &type);

  switch( type ) {
   case( kHighsOptionTypeInt ):
   case( kHighsOptionTypeBool ):
    int_options.insert( { int_counter++ , std::string( name ) } );
    break;

   case( kHighsOptionTypeDouble ):
    dbl_options.insert( { dbl_counter++ , std::string( name ) } );
    break;

   case( kHighsOptionTypeString ):
    str_options.insert( { str_counter++ , std::string( name ) } );
    break;

   default:
    std::cerr << "Unknown type from Highs_getOptionName()" << std::endl;
    return( 1 );

  }
 }

 // Generate defs file
 defs_file.open( defs_path );

 defs_file << "/* FILE GENERATED AUTOMATICALLY, DO NOT EDIT */" << std::endl
           << std::endl
           << "#ifndef __HiGHS" << HiGHS_VERSION << "_DEFS"
           << std::endl
           << "#define __HiGHS" << HiGHS_VERSION << "_DEFS"
           << std::endl << std::endl
           << "#define HiGHS_NUM_INT_PARS " << int_counter << std::endl
           << "#define HiGHS_NUM_DBL_PARS " << dbl_counter << std::endl
           << "#define HiGHS_NUM_STR_PARS " << str_counter << std::endl
           << std::endl
           << "#endif //__HiGHS" << HiGHS_VERSION << "_DEFS"
           << std::endl;

 defs_file.close();
 std::cout << "Defs file written on " << defs_path << std::endl;


// Generate maps file
 maps_file.open( maps_path );
 maps_file << "/* FILE GENERATED AUTOMATICALLY, DO NOT EDIT */" << std::endl
           << std::endl
           << "#include <interfaces/highs_c_api.h>" << std::endl
           << "#include \"HiGHSMILPSolver.h\"" << std::endl
           << std::endl
           << "using namespace SMSpp_di_unipi_it;" << std::endl
           << std::endl;

 // SMSpp_to_HiGHS_***_pars maps
 maps_file
  << "const std::array< std::string, HiGHS_NUM_INT_PARS >"
  << " HiGHSMILPSolver::SMSpp_to_HiGHS_int_pars{"
  << std::endl;
 for( const auto & i : int_options ) {
  maps_file << " \"" << i.second << "\"," << std::endl;
 }
 maps_file << "};" << std::endl;
 maps_file << std::endl;

 maps_file
  << "const std::array< std::string, HiGHS_NUM_DBL_PARS >"
  << " HiGHSMILPSolver::SMSpp_to_HiGHS_dbl_pars{"
  << std::endl;
 for( const auto & i : dbl_options ) {
  maps_file << " \"" << i.second << "\"," << std::endl;
 }
 maps_file << "};" << std::endl;
 maps_file << std::endl;

 maps_file
  << "const std::array< std::string, HiGHS_NUM_STR_PARS >"
  << " HiGHSMILPSolver::SMSpp_to_HiGHS_str_pars{"
  << std::endl;
 for( const auto & i : str_options ) {
  maps_file << " \"" << i.second << "\"," << std::endl;
 }
 maps_file << "};" << std::endl;
 maps_file << std::endl;

 // Reverse HiGHS_to_SMSpp_***_pars maps
 maps_file
  << "const std::array< std::pair< std::string, int >, HiGHS_NUM_INT_PARS >"
  << std::endl
  << " HiGHSMILPSolver::HiGHS_to_SMSpp_int_pars{" << std::endl
  << " {" << std::endl;
 for( const auto & i : int_options ) {
  maps_file << "  { \"" << i.second << "\", intFirstHiGHSPar + " << i.first
            << " },"
            << std::endl;
 }
 maps_file
  << " }" << std::endl
  << "};" << std::endl
  << std::endl;

 maps_file
  << "const std::array< std::pair< std::string, int >, HiGHS_NUM_DBL_PARS >"
  << std::endl
  << " HiGHSMILPSolver::HiGHS_to_SMSpp_dbl_pars{" << std::endl
  << " {" << std::endl;
 for( const auto & i : dbl_options ) {
  maps_file << "  { \"" << i.second << "\", dblFirstHiGHSPar + " << i.first
            << " },"
            << std::endl;
 }
 maps_file
  << " }" << std::endl
  << "};" << std::endl
  << std::endl;

 maps_file
  << "const std::array< std::pair< std::string, int >, HiGHS_NUM_STR_PARS >"
  << std::endl
  << " HiGHSMILPSolver::HiGHS_to_SMSpp_str_pars{" << std::endl
  << " {" << std::endl;
 for( const auto & i : str_options ) {
  maps_file << "  { \"" << i.second << "\", strFirstHiGHSPar + " << i.first
            << " },"
            << std::endl;
 }
 maps_file
  << " }" << std::endl
  << "};" << std::endl;

 maps_file.close();
 std::cout << "Maps file written on " << maps_path << std::endl;

 return( 0 );
}