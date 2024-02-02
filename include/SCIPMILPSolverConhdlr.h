/*--------------------------------------------------------------------------*/
/*--------------------- File SCIPMILPSolver_Conhdlr.h ----------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the SCIPMILPSolver_Conhdlr class.
 *
 * SCIPMILPSolver_Conhdlr is a tool developed for SCIPMILPSolver, in order 
 * to generate valid inequalities or even facets of the polyhedron described 
 * by a single constraint or a subset of the constraints of a single 
 * constraint class. It is essentially used to for user cuts / lazy
 * constraint separation.
 *
 * \author Enrico Calandrini \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \copyright &copy; by Enrico Calandrini
 * 
 */
/*--------------------------------------------------------------------------*/
/*----------------------------- DEFINITIONS --------------------------------*/
/*--------------------------------------------------------------------------*/

#ifndef __SCIPMILPSOLVER_CONHDLR_H
#define __SCIPMILPSOLVER_CONHDLR_H

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "scip/scip.h"
#include "objscip/objscip.h"
#include "SCIPMILPSolver.h"

/*--------------------------------------------------------------------------*/
/*----------------------------- NAMESPACE ----------------------------------*/
/*--------------------------------------------------------------------------*/

/// namespace for the Structured Modeling System++ (SMS++)
namespace SMSpp_di_unipi_it
{
/*--------------------------------------------------------------------------*/
/*------------------ CLASS SCIPMILPSolver_Conhdlr --------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------- GENERAL NOTES --------------------------------*/
/*--------------------------------------------------------------------------*/
/// Class for add lazy constraints / user cut to a SCIP instance handled by 
/// SCIPMILPSolver.

/** The SCIPMILPSolver_Conhdlr class derives from scip::ObjConshdlr and it
  * creates specific user cut/ or lazy constraint to be added within a
  * SMS++ model handled by SCIPMILPSolver.  */

class SCIPMILPSolver_Conhdlr : public scip::ObjConshdlr
{
/*--------------------------------------------------------------------------*/
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

 public:

/*--------------------------------------------------------------------------*/
/*--------------------- CONSTRUCTOR AND DESTRUCTOR -------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Constructor and Destructor
 *  @{ */

 SCIPMILPSolver_Conhdlr( SCIP* scip, /**< SCIP data structure */
    SMSpp_di_unipi_it::SCIPMILPSolver* scipmilpsolver, /**< "parent" 
                                    * SCIPMILPSolver from which the Constraint 
                                    * handler has been called. */
    unsigned char SeparationPar    /**< separation decider parameter */
    );

 ~SCIPMILPSolver_Conhdlr() override;

/** @} ---------------------------------------------------------------------*/
/*--------------------- FUNDAMENTAL CALLBACK METHODS -----------------------*/
/*--------------------------------------------------------------------------*/

/** constraint enforcing method of constraint handler for LP solutions
*
*  The method is called at the end of the node processing loop for a node 
*  where the LP was solved. The LP solution has to be checked for 
*  feasibility.
*
*  In this function we add new lazy constraints (lc) in the model.
*  NOTE: we are sure that a new lc can be added, because first a SCIP_CHECK
*  method has declared that separation is possible. 
*
*  Possible return values for *result:
*  - SCIP_CUTOFF     : the node is infeasible in the variable's bounds and 
*                      can be cut off
*  - SCIP_SEPARATED  : a cutting plane was generated
*  - SCIP_INFEASIBLE : at least one constraint is infeasible, but it was not 
*                      resolved
*  - SCIP_FEASIBLE   : all constraints of the handler are feasible
*/
 virtual SCIP_DECL_CONSENFOLP(scip_enfolp);

/** separation method of constraint handler for LP solution
*
*  Separates all constraints of the constraint handler. The method is called in 
*  the LP solution loop, which means that a valid LP solution exists.
*
*  In this function we add new user cut (uc) in the model.
*  NOTE: we are sure that a new uc can be added, because first a SCIP_CHECK
*  method has declared that separation is possible. 
*
*  possible return values for *result (if more than one applies, the first in 
*  the list should be used):
*  - SCIP_CUTOFF     : the node is infeasible in the variable's bounds and 
*                      can be cut off
*  - SCIP_SEPARATED  : a cutting plane was generated
*  - SCIP_INFEASIBLE : at least one constraint is infeasible, but it was not 
*                      resolved
*  - SCIP_FEASIBLE   : all constraints of the handler are feasible
*/
   virtual SCIP_DECL_CONSSEPALP(scip_sepalp);

/** constraint enforcing method of constraint handler for pseudo solutions
*
*  The method is called at the end of the node processing loop for a node 
*  where the LP was not solved. The pseudo solution has to be checked for 
*  feasibility. If possible, an infeasibility should be resolved by
*  branching, reducing a variable's domain to exclude the solution or adding 
*  an additional constraint. Separation is not possible, since the LP is 
*  not processed at the current node. All LP informations like
*  LP solution, slack values, or reduced costs are invalid and must not 
*  be accessed.
*
*  NOTE: At the moment in the solving loop of the algorithm SMS++ can't
*  separate pseudo-solution. For this reason this function is not yet
*  implemented
*
* Possible return values for *result:
*  - SCIP_DIDNOTRUN  : the enforcement was skipped 
*/
   virtual SCIP_DECL_CONSENFOPS(scip_enfops);

/** feasibility check method of constraint handler for primal solutions
*
*  The given solution has to be checked for feasibility.
*
*  In this method we check if possible new user cut / lazy constraint
*  can be added to the model.
* 
*  Possible return values for *result:
*  - SCIP_INFEASIBLE : at least one constraint of the handler is infeasible
*                      (i.e. a new cut can be added) 
*  - SCIP_FEASIBLE   : all constraints of the handler are feasible
*/
   virtual SCIP_DECL_CONSCHECK(scip_check);

/** variable rounding lock method of constraint handler
*
*  This method is called, after a constraint is added or removed from 
*  the transformed problem. It should update the rounding locks of all 
*  associated variables with calls to SCIPaddVarLocksType(),
*  depending on the way, the variable is involved in the constraint:
*  - If the constraint may get violated by decreasing the value of a 
*    variable, it should call SCIPaddVarLocksType(scip, var, 
*    SCIP_LOCKTYPE_MODEL, nlockspos, nlocksneg), saying that rounding 
*    down is potentially rendering the (positive) constraint infeasible 
*    and rounding up is potentially rendering the negation of the constraint 
*    infeasible.
*  - If the constraint may get violated by increasing the value of a variable, 
*    it should call SCIPaddVarLocksType(scip, var, SCIP_LOCKTYPE_MODEL, 
*    nlocksneg, nlockspos), saying that rounding up is potentially rendering 
*    the constraint's negation infeasible and rounding up is potentially 
*    rendering the constraint itself infeasible.
*  - If the constraint may get violated by changing the variable in any direction,
*    it should call SCIPaddVarLocksType(scip, var, SCIP_LOCKTYPE_MODEL, 
*    nlockspos + nlocksneg, nlockspos + nlocksneg).
*
*    NOTE: when creating the caallback method, we don't know which variables 
*    may get involved in future cuts or lazy constraints and in which direction.
*    Thus, the last method will be always called on all the problem variables,
*    in order to avoid SCIP from fixing them. 
*/
   virtual SCIP_DECL_CONSLOCK(scip_lock);

/** @} ---------------------------------------------------------------------*/
/*--------------------- ADDITIONAL CALLBACK METHODS ------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Additional Callback Methods
 *  @{ */

/** transforms constraint data into data belonging to the transformed problem */
   virtual SCIP_DECL_CONSTRANS(scip_trans);

/** frees specific constraint data */
   virtual SCIP_DECL_CONSDELETE(scip_delete);

 protected:

/*--------------------------------------------------------------------------*/
/*-------------------- PROTECTED METHODS OF THE CLASS ----------------------*/
/*--------------------------------------------------------------------------*/

 /** bitwise-encoded parameter for deciding if and when separation of user
  * cuts and lazy constraints is performed */
 unsigned char CutSepPar;

 /* parent *milpsolver from which the separator* is called */
 SMSpp_di_unipi_it::SCIPMILPSolver* parent_scipmilpsolver;

/*--------------------------------------------------------------------------*/

 };  // end( class SCIPMILPSolver_Conhdlr )

/** @} ---------------------------------------------------------------------*/
/*----------------------- METHODS FOR LAZY CONSTRAINT ----------------------*/
/*--------------------------------------------------------------------------*/
/** @name Methods for add a new constraint which will eventually be 
 *  enforced or separated producing new cuts or lazy constraints.
 *  @{ */

/** creates and captures a constraint used which will be used as a separator */
 SCIP_RETCODE SCIPcreateSCIPMILPSolver_cb(
   SCIP*        scip,               /**< SCIP data structure */
   SCIP_CONS**  cons,               /**< pointer to hold the created 
                                         constraint */
   const char*  name,               /**< name of constraint */
   std::vector< SCIP_VAR * > vars,  /**< SCIP vars */
   SCIP_Bool    initial,            /**< should the LP relaxation of 
                                         constraint be in the initial LP? */
   SCIP_Bool    separate,           /**< should the constraint be 
                                         separated during LP processing? */
   SCIP_Bool    enforce,            /**< should the constraint be enforced 
                                         during node processing? */
   SCIP_Bool    check,              /**< should the constraint be checked 
                                         for feasibility? */
   SCIP_Bool    propagate,          /**< should the constraint be propagated 
                                         during node processing? */
   SCIP_Bool    local,              /**< is constraint only valid locally? */
   SCIP_Bool    modifiable,         /**< is constraint modifiable (subject 
                                         to column generation)? */
   SCIP_Bool    dynamic,            /**< is constraint dynamic? */
   SCIP_Bool    removable           /**< should the constraint be removed 
                                         from the LP due to aging or cleanup? */
   );

/** creates and captures a a constraint which will be used as a separator
 *  with all its constraint flags set to their default values */
SCIP_RETCODE SCIPcreateSCIPMILPSolver_basiccb(
   SCIP*        scip,               /**< SCIP data structure */
   SCIP_CONS**  cons,               /**< pointer to hold the created constraint */
   const char*  name,               /**< name of constraint */
   std::vector< SCIP_VAR * > vars   /**< SCIP vars */
   );

/** @} ---------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/

}  // end( namespace SMSpp_di_unipi_it )

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

#endif  /* SCIPMILPSolver_Conhdlr.h included */

/*--------------------------------------------------------------------------*/
/*------------------ End File SCIPMILPSolver_Conhdlr.h ---------------------*/
/*--------------------------------------------------------------------------*/
