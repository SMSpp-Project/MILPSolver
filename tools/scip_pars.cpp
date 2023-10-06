/*--------------------------------------------------------------------------*/
/*---------------------------- File cpx_pars.cpp ---------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Small tool for parsing the macros and the standard std::maps that enable
 * the support for plain SCIP parameters into SCIPMILPSolver.
 *
 * The tool generates two files, SCIP<SCIP_VERSION>_defs.h and
 * SCIP<SCIP_VERSION>_maps.h in the specified path.
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

#include <scip/scip.h>

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
 docopt_desc = "SCIP parameter map generator.\n";
 exe = get_filename( argv[ 0 ] );
 process_args( argc , argv );

 // if (! std::filesystem::exists(path)) {
 //  std::filesystem::create_directory(path);
 // }

 std::string defs_filename =
  "SCIP" + std::to_string( SCIP_VERSION ) + "_defs.h";
 std::string maps_filename =
  "SCIP" + std::to_string( SCIP_VERSION ) + "_maps.h";

 // auto defs_path = path / defs_filename;
 // auto maps_path = path / maps_filename;

 auto defs_path = path + "/" + defs_filename;
 auto maps_path = path + "/" + maps_filename;

 std::ofstream defs_file;
 std::ofstream maps_file;

 SCIP * scip{};
 const char * name;

 std::map< int , std::string > int_parameters;
 std::map< int , std::string > dbl_parameters;
 std::map< int , std::string > str_parameters;

 int int_counter = 0;
 int dbl_counter = 0;
 int str_counter = 0;

 SCIP_CALL_ABORT( SCIPcreate( &scip ) );
 if( verbose ) {
  std::cout << "SCIP_VERSION is " << SCIP_VERSION << std::endl;
 }

 SCIP_PARAM ** params = SCIPgetParams( scip );
 int nparams = SCIPgetNParams( scip );

 for( int i = 0 ; i < nparams ; ++i ) {

  SCIP_PARAM * param = params[ i ];
  name = SCIPparamGetName( param );
  auto type = SCIPparamGetType( param );

  switch( type ) {
   case( SCIP_PARAMTYPE_BOOL ):
   case( SCIP_PARAMTYPE_INT ):
   case( SCIP_PARAMTYPE_LONGINT ):
    int_parameters.insert( { int_counter++ , std::string( name ) } );
    break;

   case( SCIP_PARAMTYPE_REAL ):
    dbl_parameters.insert( { dbl_counter++ , std::string( name ) } );
    break;

   case( SCIP_PARAMTYPE_CHAR ):
   case( SCIP_PARAMTYPE_STRING ):
    str_parameters.insert( { str_counter++ , std::string( name ) } );
    break;

   default:
    std::cerr << "Unknown type from SCIPparamGetType()" << std::endl;
    return( 1 );
  }
 }

 // Generate defs file
 defs_file.open( defs_path );

 defs_file << "/* FILE GENERATED AUTOMATICALLY, DO NOT EDIT */" << std::endl
           << std::endl
           << "#ifndef __SCIP" << std::to_string( SCIP_VERSION ) << "_DEFS"
           << std::endl
           << "#define __SCIP" << std::to_string( SCIP_VERSION ) << "_DEFS"
           << std::endl << std::endl
           << "#define SCIP_NUM_INT_PARS " << int_counter << std::endl
           << "#define SCIP_NUM_DBL_PARS " << dbl_counter << std::endl
           << "#define SCIP_NUM_STR_PARS " << str_counter << std::endl
           << std::endl
           << "#endif //__SCIP" << std::to_string( SCIP_VERSION ) << "_DEFS"
           << std::endl;

 defs_file.close();
 std::cout << "Defs file written on " << defs_path << std::endl;

 // Generate maps file
 maps_file.open( maps_path );
 maps_file << "/* FILE GENERATED AUTOMATICALLY, DO NOT EDIT */" << std::endl
           << std::endl
           << "#include <string>" << std::endl
           << "#include <scip/scip.h>" << std::endl
           << "#include \"SCIPMILPSolver.h\"" << std::endl
           << std::endl
           << "using namespace SMSpp_di_unipi_it;" << std::endl
           << std::endl;

 // SMSpp_to_SCIP_***_pars maps
 maps_file
  << "const std::array< std::string, SCIP_NUM_INT_PARS >"
  << " SCIPMILPSolver::SMSpp_to_SCIP_int_pars{"
  << std::endl;
 for( const auto & i : int_parameters ) {
  maps_file << " \"" << i.second << "\"," << std::endl;
 }
 maps_file << "};" << std::endl;
 maps_file << std::endl;

 maps_file
  << "const std::array< std::string, SCIP_NUM_DBL_PARS >"
  << " SCIPMILPSolver::SMSpp_to_SCIP_dbl_pars{"
  << std::endl;
 for( const auto & i : dbl_parameters ) {
  maps_file << " \"" << i.second << "\"," << std::endl;
 }
 maps_file << "};" << std::endl;
 maps_file << std::endl;

 maps_file
  << "const std::array< std::string, SCIP_NUM_STR_PARS >"
  << " SCIPMILPSolver::SMSpp_to_SCIP_str_pars{"
  << std::endl;
 for( const auto & i : str_parameters ) {
  maps_file << " \"" << i.second << "\"," << std::endl;
 }
 maps_file << "};" << std::endl;
 maps_file << std::endl;

 // Reverse SCIP_to_SMSpp_***_pars maps
 maps_file
  << "const std::array< std::pair< std::string, int >, SCIP_NUM_INT_PARS >"
  << std::endl
  << " SCIPMILPSolver::SCIP_to_SMSpp_int_pars{" << std::endl
  << " {" << std::endl;
 for( const auto & i : int_parameters ) {
  maps_file << "  { \"" << i.second << "\", intFirstSCIPPar + " << i.first
            << " },"
            << std::endl;
 }
 maps_file
  << " }" << std::endl
  << "};" << std::endl
  << std::endl;

 maps_file
  << "const std::array< std::pair< std::string, int >, SCIP_NUM_DBL_PARS >"
  << std::endl
  << " SCIPMILPSolver::SCIP_to_SMSpp_dbl_pars{" << std::endl
  << " {" << std::endl;
 for( const auto & i : dbl_parameters ) {
  maps_file << "  { \"" << i.second << "\", dblFirstSCIPPar + " << i.first
            << " },"
            << std::endl;
 }
 maps_file
  << " }" << std::endl
  << "};" << std::endl
  << std::endl;

 maps_file
  << "const std::array< std::pair< std::string, int >, SCIP_NUM_STR_PARS >"
  << std::endl
  << " SCIPMILPSolver::SCIP_to_SMSpp_str_pars{" << std::endl
  << " {" << std::endl;
 for( const auto & i : str_parameters ) {
  maps_file << "  { \"" << i.second << "\", strFirstSCIPPar + " << i.first
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
