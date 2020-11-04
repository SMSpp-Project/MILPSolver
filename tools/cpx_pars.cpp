/*--------------------------------------------------------------------------*/
/*---------------------------- File cpx_pars.cpp ---------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Small tool for parsing the macros and the standard std::maps that enable
 * the support for plain CPLEX parameters into CPXMILPSolver.
 *
 * The tool generates two files, CPX<CPX_VERSION>_defs.h and
 * CPX<CPX_VERSION>_maps.h. These are automatically placed in the include
 * directory of the MILPSolver project (assumed to be ../include).
 *
 * \author Niccolo' Iardella \n
 *         Operations Research Group \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * Copyright &copy; Niccolo' Iardella
 */
/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include <iostream>
#include <fstream>
#include <cstring>
#include <ilcplex/cplex.h>
#include <map>

int main( int argc, char ** argv ) {

 std::string defs_file_name(
  "../include/CPX" + std::to_string( CPX_VERSION ) + "_defs.h" );
 std::string maps_file_name(
  "../include/CPX" + std::to_string( CPX_VERSION ) + "_maps.h" );

 std::ofstream defs_file;
 std::ofstream maps_file;

 CPXENVptr env;
 int status;
 char name[CPX_STR_PARAM_MAX];

 std::map< int, std::string > int_parameters;
 std::map< int, std::string > dbl_parameters;
 std::map< int, std::string > str_parameters;

 int int_counter = 0;
 int dbl_counter = 0;
 int str_counter = 0;
 int total = 0;

 env = CPXopenCPLEX( &status );
 std::cout << "CPX_VERSION is " << CPX_VERSION << std::endl;

 for( int i = CPX_PARAM_ALL_MIN; i <= CPX_PARAM_ALL_MAX; ++i ) {

#if CPX_VERSION < 12090000
  status = CPXgetparamname(env, i, name);
#else
  status = CPXgetparamhiername( env, i, name );
#endif

  if( status == CPXERR_BAD_PARAM_NUM ) {
   continue;
  }

  if( strlen( name ) == 0 ) {
   continue;
  }

  if( status == 0 ) {
   int type;
   CPXgetparamtype( env, i, &type );

   switch( type ) {
    case CPX_PARAMTYPE_INT:
    case CPX_PARAMTYPE_LONG:
     int_parameters.insert( { int_counter++, std::string( name ) } );
     break;

    case CPX_PARAMTYPE_DOUBLE:
     // Remove unsupported internal parameters
     if( strcmp( name, "CPXPARAM_Internal_cfilemul" ) == 0 ||
         strcmp( name, "CPXPARAM_Internal_rfilemul" ) == 0 ||
         strcmp( name, "CPXPARAM_Internal_singtol" ) == 0 ) {
      break;
     }
     dbl_parameters.insert( { dbl_counter++, std::string( name ) } );
     break;

    case CPX_PARAMTYPE_STRING:
     str_parameters.insert( { str_counter++, std::string( name ) } );
     break;

    default:
     std::cerr << "Unknown type from CPXgetparamtype()" << std::endl;
     return 1;
   }

  } else {
   std::cerr << "Unknown error in CPXgetparamhiername()" << std::endl;
   return 1;
  }
 }

 // Generate defs file
 defs_file.open( defs_file_name );

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
 std::cout << "Defs file written on " << defs_file_name << std::endl;

 // Generate maps file
 maps_file.open( maps_file_name );
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
 for( const auto & i: int_parameters ) {
  maps_file << " " << i.second << "," << std::endl;
 }
 maps_file << "};" << std::endl;
 maps_file << std::endl;

 maps_file
  << "const std::array< int, CPX_NUM_DBL_PARS >"
  << " CPXMILPSolver::SMSpp_to_CPLEX_dbl_pars{"
  << std::endl;
 for( const auto & i: dbl_parameters ) {
  maps_file << " " << i.second << "," << std::endl;
 }
 maps_file << "};" << std::endl;
 maps_file << std::endl;

 maps_file
  << "const std::array< int, CPX_NUM_STR_PARS >"
  << " CPXMILPSolver::SMSpp_to_CPLEX_str_pars{"
  << std::endl;
 for( const auto & i: str_parameters ) {
  maps_file << " " << i.second << "," << std::endl;
 }
 maps_file << "};" << std::endl;
 maps_file << std::endl;

 // Reverse CPLEX_to_SMSpp_***_pars maps
 maps_file
  << "const std::array< std::pair< int, int >, CPX_NUM_INT_PARS >" << std::endl
  << " CPXMILPSolver::CPLEX_to_SMSpp_int_pars{" << std::endl
  << " {" << std::endl;
 for( const auto & i: int_parameters ) {
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
 for( const auto & i: dbl_parameters ) {
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
 for( const auto & i: str_parameters ) {
  maps_file << "  { " << i.second << ", strFirstCPLEXPar + " << i.first << " },"
            << std::endl;
 }
 maps_file
  << " }" << std::endl
  << "};" << std::endl;
 maps_file.close();
 std::cout << "Maps file written on " << maps_file_name << std::endl;

 return 0;
}
