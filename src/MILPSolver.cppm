/*--------------------------------------------------------------------------*/
/*-------------------------- File MILPSolver.cppm --------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Interface of the C++20 module smspp.milpsolver, which exports every name
 * that the headers of MILPSolver declare at namespace scope, together with
 * the modules it depends on; a translation unit can then write
 *
 *     import smspp.milpsolver;
 *
 * in place of the #include of the headers. The macros (e.g., those of the
 * factory in SMSFactory.h) are not exported by a module, so a file that uses
 * them still includes the header that defines them.
 *
 * This file is written by SMS++/cmake/smspp_module_exports.py, which reads
 * the names from the AST of the headers: it is not to be edited by hand.
 *
 * The headers *_defs.h, *_maps.h are not in this module: they need a library
 * of their own, or only a source includes them.
 *
 * \author Enrico Calandrini \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \author Antonio Frangioni \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \author Niccolo' Iardella \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \author Donato Meoli \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \copyright &copy; by Enrico Calandrini, Antonio Frangioni,
 *                   Niccolo' Iardella,
 *                   Donato Meoli
 */
/*--------------------------------------------------------------------------*/
/*------------------------- GLOBAL MODULE FRAGMENT -------------------------*/
/*--------------------------------------------------------------------------*/

module;

#if defined( MILPSolver_MODULE_CPLEX )
 #include "CPXMILPSolver.h"
#endif
#if defined( MILPSolver_MODULE_Clarabel )
 #include "ClarabelMILPSolver.h"
#endif
#if defined( MILPSolver_MODULE_DAQP )
 #include "DAQPMILPSolver.h"
#endif
#if defined( MILPSolver_MODULE_GUROBI )
 #include "GRBMILPSolver.h"
#endif
#if defined( MILPSolver_MODULE_HiGHS )
 #include "HiGHSMILPSolver.h"
#endif
#include "MILPSolver.h"
#if defined( MILPSolver_MODULE_MOSEK )
 #include "MOSEKMILPSolver.h"
#endif
#if defined( MILPSolver_MODULE_OSQP )
 #include "OSQPMILPSolver.h"
#endif
#if defined( MILPSolver_MODULE_PIPS )
 #include "PIPSMILPSolver.h"
#endif
#if defined( MILPSolver_MODULE_PIQP )
 #include "PIQPMILPSolver.h"
#endif
#if defined( MILPSolver_MODULE_ProxQP )
 #include "ProxQPMILPSolver.h"
#endif
#if defined( MILPSolver_MODULE_qpmad )
 #include "QPMADMILPSolver.h"
#endif
#if defined( MILPSolver_MODULE_qpOASES )
 #include "QPOASESMILPSolver.h"
#endif
#if defined( MILPSolver_MODULE_quadprog )
 #include "QuadProgMILPSolver.h"
#endif
#if defined( MILPSolver_MODULE_SCIP )
 #include "SCIPMILPSolver.h"
#endif

/*--------------------------------------------------------------------------*/
/*----------------------------- MODULE PURVIEW -----------------------------*/
/*--------------------------------------------------------------------------*/

export module smspp.milpsolver;

export import smspp;

export namespace SMSpp_di_unipi_it {

 // CPXMILPSolver.h
#if defined( MILPSolver_MODULE_CPLEX )
 using SMSpp_di_unipi_it::CPXMILPSolver;
#endif

 // ClarabelMILPSolver.h
#if defined( MILPSolver_MODULE_Clarabel )
 using SMSpp_di_unipi_it::ClarabelMILPSolver;
#endif

 // DAQPMILPSolver.h
#if defined( MILPSolver_MODULE_DAQP )
 using SMSpp_di_unipi_it::DAQPMILPSolver;
#endif

 // GRBMILPSolver.h
#if defined( MILPSolver_MODULE_GUROBI )
 using SMSpp_di_unipi_it::GRBMILPSolver;
#endif

 // HiGHSMILPSolver.h
#if defined( MILPSolver_MODULE_HiGHS )
 using SMSpp_di_unipi_it::HiGHSMILPSolver;
#endif

 // MILPSolver.h
 using SMSpp_di_unipi_it::MILPSolver;

 // MOSEKMILPSolver.h
#if defined( MILPSolver_MODULE_MOSEK )
 using SMSpp_di_unipi_it::MOSEKMILPSolver;
#endif

 // OSQPMILPSolver.h
#if defined( MILPSolver_MODULE_OSQP )
 using SMSpp_di_unipi_it::OSQPMILPSolver;
#endif

 // PIPSMILPSolver.h
#if defined( MILPSolver_MODULE_PIPS )
 using SMSpp_di_unipi_it::PIPSMILPSolver;
#endif

 // PIQPMILPSolver.h
#if defined( MILPSolver_MODULE_PIQP )
 using SMSpp_di_unipi_it::PIQPMILPSolver;
#endif

 // ProxQPMILPSolver.h
#if defined( MILPSolver_MODULE_ProxQP )
 using SMSpp_di_unipi_it::ProxQPMILPSolver;
#endif

 // QPMADMILPSolver.h
#if defined( MILPSolver_MODULE_qpmad )
 using SMSpp_di_unipi_it::QPMADMILPSolver;
#endif

 // QPOASESMILPSolver.h
#if defined( MILPSolver_MODULE_qpOASES )
 using SMSpp_di_unipi_it::QPOASESMILPSolver;
#endif

 // QuadProgMILPSolver.h
#if defined( MILPSolver_MODULE_quadprog )
 using SMSpp_di_unipi_it::QuadProgMILPSolver;
#endif

 // SCIPMILPSolver.h
#if defined( MILPSolver_MODULE_SCIP )
 using SMSpp_di_unipi_it::SCIPMILPSolver;
 using SMSpp_di_unipi_it::SCIPMILPSolver_CBK;
#endif
}

export namespace proxsuite {
 namespace proxqp {
  namespace dense {
  }
 }
}

export namespace proxsuite {
 namespace proxqp {
  namespace sparse {
  }
 }
}

export namespace qpOASES {
}

/*--------------------------------------------------------------------------*/
/*------------------------ End File MILPSolver.cppm ------------------------*/
/*--------------------------------------------------------------------------*/
