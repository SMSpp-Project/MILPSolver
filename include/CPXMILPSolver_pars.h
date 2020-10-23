/*--------------------------------------------------------------------------*/
/*----------------------- File CPXMILPSolver_pars.h ------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the CPLEX parameters.
 *
 * This file contains some macros that enable CPXMILPSolver support to
 * plain CPLEX parameters. For each new CPLEX version, the proper defines
 * must be added, using the cpx_pars tool for convenience.
 *
 * \author Antonio Frangioni \n
 *         Operations Research Group \n
 *         Dipartimento di Informatica \n
 *         Università di Pisa \n
 *
 * \author Niccolò Iardella \n
 *         Operations Research Group \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \copyright &copy; Antonio Frangioni, Niccolò Iardella
 */

/*--------------------------------------------------------------------------*/

#ifndef __CPXMILPSOLVER_PARS
#define __CPXMILPSOLVER_PARS

#if CPX_VERSION == 12080000
#define CPX_NUM_INT_PARS 135
#define CPX_NUM_DBL_PARS 44
#define CPX_NUM_STR_PARS 5

#elif CPX_VERSION == 12090000
#define CPX_NUM_INT_PARS 137
#define CPX_NUM_DBL_PARS 44
#define CPX_NUM_STR_PARS 5

#elif CPX_VERSION >= 12100000
#define CPX_NUM_INT_PARS 137
#define CPX_NUM_DBL_PARS 43
#define CPX_NUM_STR_PARS 5
#endif

#endif //__CPXMILPSOLVER_PARS
