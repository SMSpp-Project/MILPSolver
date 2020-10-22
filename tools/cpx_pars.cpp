/*--------------------------------------------------------------------------*/
/*---------------------------- File cpx_pars.cpp ---------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Small tool for parsing the standard std::maps that add the support for
 * plain CPLEX parameters to CPXMILPSolver.
 *
 * The output of the tool is some C++ code that must be copied and pasted
 * in CPXMILPSolver_pars.cpp and CPXMILPSolver.h.
 * The output depends on CPLEX version.
 *
 * In theory, users shouldn't need to do this, as we developers plan to
 * update those files when needed (that is, at each new CPLEX version.
 * This tool is provided just in case.
 *
 * \author Niccolò Iardella \n
 *         Operations Research Group \n
 *         Dipartimento di Informatica \n
 *         Università di Pisa \n
 *
 * Copyright &copy; Niccolò Iardella
 */
/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include <iostream>
#include <cstring>
#include <ilcplex/cplex.h>
#include <map>

int main() {

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
 std::cout
  << "-------------------- COPY THE FOLLOWING IN THE HEADER --------------------"
  << std::endl
  << "#if CPX_VERSION == " << CPX_VERSION << std::endl
  << "#define CPX_NUM_INT_PARS " << int_counter << std::endl
  << "#define CPX_NUM_DBL_PARS " << dbl_counter << std::endl
  << "#define CPX_NUM_STR_PARS " << str_counter << std::endl
  << "#endif" << std::endl
  << "--------------------------------------------------------------------------"
  << std::endl
  << std::endl
  << "--------------------- COPY THE FOLLOWING IN THE BODY ---------------------"
  << std::endl
  << "#if CPX_VERSION == " << CPX_VERSION << std::endl;

 // SMSpp_to_CPLEX_***_pars maps
 std::cout
  << "const std::map< int, int > CPXMILPSolver::SMSpp_to_CPLEX_int_pars{"
  << std::endl;
 for( const auto & i: int_parameters ) {
  std::cout << "{ intFirstCPLEXPar + " << i.first << ", " << i.second << " },"
            << std::endl;
 }
 std::cout << "};" << std::endl;
 std::cout << std::endl;

 std::cout
  << "const std::map< int, int > CPXMILPSolver::SMSpp_to_CPLEX_dbl_pars{"
  << std::endl;
 for( const auto & i: dbl_parameters ) {
  std::cout << "{ dblFirstCPLEXPar + " << i.first << ", " << i.second << " },"
            << std::endl;
 }
 std::cout << "};" << std::endl;
 std::cout << std::endl;

 std::cout
  << "const std::map< int, int > CPXMILPSolver::SMSpp_to_CPLEX_str_pars{"
  << std::endl;
 for( const auto & i: str_parameters ) {
  std::cout << "{ strFirstCPLEXPar + " << i.first << ", " << i.second << " },"
            << std::endl;
 }
 std::cout << "};" << std::endl;
 std::cout << std::endl;

 // Reverse CPLEX_to_SMSpp_***_pars maps
 std::cout
  << "const std::map< int, int > CPXMILPSolver::CPLEX_to_SMSpp_int_pars{"
  << std::endl;
 for( const auto & i: int_parameters ) {
  std::cout << "{ " << i.second << ", intFirstCPLEXPar + " << i.first << " },"
            << std::endl;
 }
 std::cout << "};" << std::endl;
 std::cout << std::endl;

 std::cout
  << "const std::map< int, int > CPXMILPSolver::CPLEX_to_SMSpp_dbl_pars{"
  << std::endl;
 for( const auto & i: dbl_parameters ) {
  std::cout << "{ " << i.second << ", dblFirstCPLEXPar + " << i.first << " },"
            << std::endl;
 }
 std::cout << "};" << std::endl;
 std::cout << std::endl;

 std::cout
  << "const std::map< int, int > CPXMILPSolver::CPLEX_to_SMSpp_str_pars{"
  << std::endl;
 for( const auto & i: str_parameters ) {
  std::cout << "{ " << i.second << ", strFirstCPLEXPar + " << i.first << " },"
            << std::endl;
 }
 std::cout << "};" << std::endl;
 std::cout
  << "--------------------------------------------------------------------------"
  << std::endl
  << "#endif" << std::endl;
 return 0;
}
