/*-------------------------------------------------------------------------------*/
/*------------------------- File SCIPMILPSolver_Conhdlr.cpp ---------------------*/
/*-------------------------------------------------------------------------------*/
/** @file
 * Implementation of the SCIPMILPSolver_Conhdlr class.
 *
 * \author Enrico Calandrini \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \copyright &copy; by Enrico Calandrini
 */
/*-------------------------------------------------------------------------------*/
/*---------------------------- IMPLEMENTATION -----------------------------------*/
/*-------------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ---------------------------------------*/
/*-------------------------------------------------------------------------------*/

#include <assert.h>

#include "SCIPMILPSolverConhdlr.h"
#include <scip/cons_linear.h>

/*-------------------------------------------------------------------------------*/
/*---------------------------- DEFINITIONS --------------------------------------*/
/*-------------------------------------------------------------------------------*/

/* fundamental constraint handler properties */
#define CONSHDLR_NAME          "SCIPMILPSolver_Conhdlr"
#define CONSHDLR_DESC          "Constraint handler for SCIPMILPSolver. This \
particular class will provide callback method inside a SCIPMILPSolver object, \
allowing it to generate user cuts or lazy constraints. "
#define CONSHDLR_ENFOPRIORITY  1.0  // priority of the constraint handler for
                                    // constraint enforcing
#define CONSHDLR_CHECKPRIORITY -1   // priority of the constraint handler for 
                                    // checking feasibility
#define CONSHDLR_EAGERFREQ     -1   // frequency for using all instead of only 
                                    // the useful constraints in separation,
                                    //propagation and enforcement
#define CONSHDLR_NEEDSCONS     TRUE // should the constraint handler be skipped, 
                                    // if no constraints are available? 
#define CONSHDLR_SEPAPRIORITY  1.0  // priority of the constraint handler for 
                                     // separation
#define CONSHDLR_SEPAFREQ      1.0   // frequency for separating cuts; zero 
                                     // means to separate only in the root node 
#define CONSHDLR_DELAYSEPA     FALSE // should separation method be delayed, if 
                                     // other separators found cuts?
#define CONSHDLR_PROPFREQ      -1    // frequency for propagating domains; 
                                     // zero means only preprocessing propagation
#define CONSHDLR_DELAYPROP     FALSE // should propagation method be delayed, 
                                     // if other propagators found reductions?
#define CONSHDLR_PROP_TIMING   SCIP_PROPTIMING_BEFORELP 
#define CONSHDLR_PRESOLTIMING  SCIP_PRESOLTIMING_MEDIUM 
#define CONSHDLR_MAXPREROUNDS  -1    // maximal number of presolving rounds the
                                     // constraint handler participates in 
                                     // (-1: no limit) */

/*-------------------------------------------------------------------------------*/
/*---------------------------- NAMESPACE AND USING ------------------------------*/
/*-------------------------------------------------------------------------------*/

using namespace SMSpp_di_unipi_it;

/*-------------------------------------------------------------------------------*/
/*------------------------------ DATA STRUCTURES --------------------------------*/
/*-------------------------------------------------------------------------------*/

/** constraint data for lazy constraints */
struct SCIP_ConsData
{
   SCIP_Bool new_cut; // Parameter used to tell if the constraint handler found 
                      // with the CHECK function a possible cut
};

/*-------------------------------------------------------------------------------*/
/*------------------------- CONSTRUCTOR AND DESTRUCTOR --------------------------*/
/*-------------------------------------------------------------------------------*/

SCIPMILPSolver_Conhdlr::SCIPMILPSolver_Conhdlr( SCIP* scip,
    SMSpp_di_unipi_it::SCIPMILPSolver* scipmilpsolver,
     unsigned char SeparationPar
    ) : ObjConshdlr( scip , CONSHDLR_NAME , CONSHDLR_DESC, CONSHDLR_SEPAPRIORITY,
         CONSHDLR_ENFOPRIORITY, CONSHDLR_CHECKPRIORITY, CONSHDLR_SEPAFREQ, 
         CONSHDLR_PROPFREQ, CONSHDLR_EAGERFREQ , CONSHDLR_MAXPREROUNDS,
         CONSHDLR_DELAYSEPA, CONSHDLR_DELAYPROP, CONSHDLR_NEEDSCONS,
         CONSHDLR_PROP_TIMING, CONSHDLR_PRESOLTIMING )
   {
      parent_scipmilpsolver = scipmilpsolver; // set parent scipmilpsolver

      CutSepPar = SeparationPar; // set parameter for deciding if and when 
                                 // separation should be performed
   }

 SCIPMILPSolver_Conhdlr::~SCIPMILPSolver_Conhdlr(){
 }

/*-------------------------------------------------------------------------------*/
/*--------------------------------- LOCAL METHODS -------------------------------*/
/*-------------------------------------------------------------------------------*/

/** local method used to perform the separation within the constraint handler.
*
*  The method is called from the verify_separation or scipmilpsolver_separation 
*  function, in order to create the structures to add a lazy constraint/user cut. 
*/
 static SCIP_RETCODE perform_separation( 
   SCIP*                 scip,                         /* SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,                     /* the constraint handler 
                                                         itself */
   SMSpp_di_unipi_it::SCIPMILPSolver* scipmilpsolver,  /* the milpsolver model from
                                                        * which the conhdlr has been 
                                                        * called */
   SCIP_CONS**           conss,                        /* array of constraints 
                                                        * to process */
   SCIP_SOL*             sol,                          /* primal solution that should 
                                                        * be separated */
   SCIP_Bool             enforce,                      /* whether we are in 
                                                        *enforcing */
   SCIP_RESULT*          result,                        /* pointer to store the result 
                                                        * of the separation call */
	std::vector< int >   & rmatbeg ,
	std::vector< int >   & rmatind ,
	std::vector< double > & rmatval ,
	std::vector< double > & rhs ,
   std::vector< double > & lhs 
   )
 {
 std::vector< SCIP_VAR * > scip_vars;
 int nvars;

 // get variables data from the scipmilpsolver 
 /** NOTE: This is essential, because using SCIPMILPSolver we are adding auxiliary 
 * variables which don't have to be considered in the computation of the best 
 * solution. */
 scip_vars =  scipmilpsolver->get_SCIP_var();
 nvars = scip_vars.size();

 // this is a critical section where different SCIP threads may compete
 // for access to the Block: ensure mutual exclusion
 scipmilpsolver->set_f_cb_mutex();

 // get the feasible solution
 std::vector< double > x( nvars );
 for( int i = 0 ; i < nvars ; ++i ){
   SCIP_VAR* var = scip_vars[ i ];
   SCIP_Real best_sol = SCIPgetSolVal( scip , sol, var );
   x[i] = best_sol;
 }

 // write it in the Variable of the Block
 scipmilpsolver->get_var_solution( x );

 // now perform the lazy constraint separation with the right Configuration

 if( enforce == FALSE ){
  // Adding a user cut
  int depth;
  depth = SCIPgetSubscipDepth( scip ); // find the depth of the current node

  scipmilpsolver->perform_separation( 
   scipmilpsolver->get_cfg( depth ? 1 : 0 ) ,
   rmatbeg , rmatind , rmatval , rhs , lhs );
  }
 else{
  // Adding a lazy constraint
  scipmilpsolver->perform_separation( 
   scipmilpsolver->get_cfg( 2 ) ,
   rmatbeg , rmatind , rmatval , rhs , lhs );
  }

 // critical section ends here, release the mutex
 scipmilpsolver->unset_f_cb_mutex();

 return SCIP_OKAY;
}

/** local method used to check if a user cut/lazy constraint exists
*
*  The method is called from the scip_check function, in order to state if a 
*  lazy constraint/user cut is available and can be added. 
*/
 static SCIP_RETCODE verify_separation( 
   SCIP*                 scip,                         /* SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,                     /* the constraint handler 
                                                         itself */
   SMSpp_di_unipi_it::SCIPMILPSolver* scipmilpsolver,  /* the milpsolver model from
                                                        * which the conhdlr has been 
                                                        * called */
   SCIP_CONS**           conss,                        /* array of constraints 
                                                        * to process */
   SCIP_SOL*             sol,                          /* primal solution that should 
                                                        * be separated */
   SCIP_Bool             enforce,                      /* whether we are in 
                                                        *enforcing */
   SCIP_RESULT*          result                        /* pointer to store the result 
                                                        * of the separation call */
   )
 {
 assert(result != NULL);

 // get all required structures
 SCIP_CONSDATA* consdata;

 consdata = SCIPconsGetData(conss[0]);
 assert(consdata != NULL);

 std::vector< SCIP_VAR * > scip_vars;
 int nvars;

 // get variables data from the scipmilpsolver 
 /** NOTE: This is essential, because using SCIPMILPSolver we are adding auxiliary 
 * variables which don't have to be considered in the computation of the best 
 * solution. */
 scip_vars =  scipmilpsolver->get_SCIP_var();
 nvars = scip_vars.size();

 // now perform the lazy constraint separation with the right Configuration
 std::vector< int > rmatbeg;
 std::vector< int > rmatind;
 std::vector< double > rmatval;
 std::vector< double > rhs;
 std::vector< double > lhs;

 perform_separation( scip , conshdlr , scipmilpsolver , conss , sol ,
                     enforce , result , rmatbeg , rmatind , rmatval ,
                     rhs , lhs );

 // Understand if any user cut/ lazy constraint are available
 if( ! rmatbeg.empty() ){
    // at least acutting plane has been found
    *result = SCIP_INFEASIBLE;
    consdata->new_cut = TRUE; // inform constraint data that 
                              // a new cut has been found
  }
 else{
    // the solution is already feasible
    *result = SCIP_FEASIBLE;
    consdata->new_cut = TRUE; // inform constraint data that 
                              // a no cut is available
 }

 return SCIP_OKAY;
}


/** local method used to add a lazy constraint or a user cut
*
*  The method is called from a ENFOLP/SEPALP, respectively to add a new
*  lazy constraint or a user cut to the model. The type of constraint is
*  specified with the 'enforce' parameter.
*/
 static
 SCIP_RETCODE scipmilpsolver_separation( 
   SCIP*                 scip,                         /* SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,                     /* the constraint handler 
                                                         itself */
   SMSpp_di_unipi_it::SCIPMILPSolver* scipmilpsolver,  /* the milpsolver model from
                                                        * which the conhdlr has been 
                                                        * called */
   SCIP_CONS**           conss,                        /* array of constraints 
                                                        * to process */
   SCIP_SOL*             sol,                          /* primal solution that should 
                                                        * be separated */
   SCIP_Bool             enforce,                      /* whether we are in 
                                                        *enforcing */
   SCIP_RESULT*          result                        /* pointer to store the result 
                                                        * of the separation call */
   )
 {
 assert(result != NULL);

 // get all required structures
 SCIP_CONSDATA* consdata;

 consdata = SCIPconsGetData(conss[0]);
 assert(consdata != NULL);

 /* if a new cut is available, the constraint data must be already informed */
 if( !consdata->new_cut )
  // strange, but nothing to do
  return SCIP_OKAY;

 std::vector< SCIP_VAR * > scip_vars;
 int nvars;

 // get variables data from the scipmilpsolver 
 /** NOTE: This is essential, because using SCIPMILPSolver we are adding auxiliary 
 * variables which don't have to be considered in the computation of the best 
 * solution. */
 //scip_vars = consdata->vars;
 scip_vars =  scipmilpsolver->get_SCIP_var();
 nvars = scip_vars.size();

 // now perform the lazy constraint separation with the right Configuration
 std::vector< int > rmatbeg;
 std::vector< int > rmatind;
 std::vector< double > rmatval;
 std::vector< double > rhs;
 std::vector< double > lhs;

 perform_separation( scip , conshdlr , scipmilpsolver , conss , sol ,
                     enforce , result , rmatbeg , rmatind , rmatval ,
                     rhs , lhs );

 // if any lazy constraint/user cut was generated, add them
 if( ! rmatbeg.empty() )
   for( int c = 0 ; c < rhs.size() ; ++c ){
      SCIP_ROW* row;
      SCIP_CALL( SCIPcreateEmptyRowConshdlr( scip, &row, conshdlr, 
                     "scipmilpsolver_cut", lhs[ c ] , rhs[ c ], 
                     FALSE, FALSE, TRUE) );

      SCIP_CALL( SCIPcacheRowExtensions(scip, row) );

      int nnz; // number of nonzero coefficients in the actual lazy constraint
      int beg_idx = rmatbeg[ c ]; // idx from where new coefficients begin 
         
      // retrieve number of nonzero
      if( c < rhs.size() - 1)
         nnz = rmatbeg[ c + 1 ] - rmatbeg[ c ]; 
      else
         nnz = rmatind.size() - rmatbeg[ c ];

      for( int counter = 0 ; counter < nnz ; ++counter ){
         int var_idx = rmatind[ beg_idx + counter ];

         SCIP_CALL( SCIPaddVarToRow( scip , row , scip_vars[ var_idx ] , 
                        rmatval[ beg_idx + counter ] ) );
         }

      SCIP_CALL( SCIPflushRowExtensions(scip, row) );
      //SCIP_CALL( SCIPprintRow(scip, row, NULL));

      // Add violated cut. If we are enforcing, then this is enough to add 
      // the cut. Otherwise (we are separating), we check whether the
      // cut is efficacious.
      if( enforce || SCIPisCutEfficacious( scip , sol , row ) ){
         SCIP_Bool infeasible;
         SCIP_CALL( SCIPaddRow( scip , row , FALSE , &infeasible) );
         if ( infeasible )
            *result = SCIP_CUTOFF;
         else
            *result = SCIP_SEPARATED;
         }
      
      SCIP_CALL( SCIPreleaseRow(scip, &row) );
      }

   return SCIP_OKAY;
}

/*-------------------------------------------------------------------------------*/
/*------------------- CALLBACK METHODS OF CONSTRAINT HANDLER --------------------*/
/*-------------------------------------------------------------------------------*/

/*-------------------------------------------------------------------------------*/

/** constraint enforcing method of constraint handler for LP solutions */
SCIP_DECL_CONSENFOLP(SCIPMILPSolver_Conhdlr::scip_enfolp)
{  /*lint --e{715}*/

 assert( result != NULL );

 if( ! ( CutSepPar & 4 ) )  // but we don't do lazy constraint separation
   return SCIP_OKAY;         // nothing to do

 SCIP_CALL( scipmilpsolver_separation ( scip, conshdlr, parent_scipmilpsolver,
                conss , NULL , TRUE, result ) );

 return SCIP_OKAY;

}

/*-------------------------------------------------------------------------------*/

/** separation method of constraint handler for LP solution */
SCIP_DECL_CONSSEPALP(SCIPMILPSolver_Conhdlr::scip_sepalp)
{
 assert( result != NULL );

 if( ! ( CutSepPar & 3 ) )  // but we don't do user cut separation
   return SCIP_OKAY;        // nothing to do

 int depth;
 depth = SCIPgetSubscipDepth( scip ); // find the depth of the current node

 // if we are at a depth for which separation is not enabled
 if( ( ( ! depth ) && ( ! ( CutSepPar & 1 ) ) ) ||
   ( depth && ( ! ( CutSepPar & 2 ) ) ) )
     return SCIP_OKAY;     // nothing to do

 SCIP_CALL( scipmilpsolver_separation ( scip, conshdlr , parent_scipmilpsolver,
                conss , NULL , FALSE, result) );

 return SCIP_OKAY;
}

/*-------------------------------------------------------------------------------*/

/** constraint enforcing method of constraint handler for pseudo solutions */
SCIP_DECL_CONSENFOPS(SCIPMILPSolver_Conhdlr::scip_enfops)
{  /*lint --e{715}*/
   *result = SCIP_DIDNOTRUN;
   return SCIP_OKAY;
}

/** feasibility check method of constraint handler for primal solutions */
SCIP_DECL_CONSCHECK(SCIPMILPSolver_Conhdlr::scip_check)
{  /*lint --e{715}*/

   // Retrieve actual SCIP stage
   SCIP_STAGE current_stage = SCIPgetStage( scip );

   int sep = 0; // parameter to decide if separation is enabled

   sep = verify_separation( scip, conshdlr, parent_scipmilpsolver,
                conss , sol , TRUE, result);
   
   if( sep == 0 ) // No separation is available
    *result = SCIP_FEASIBLE;
   else // we can add some user cut/lazy constraint
    *result = SCIP_INFEASIBLE;

   return SCIP_OKAY;
}

/** variable rounding lock method of constraint handler */
SCIP_DECL_CONSLOCK(SCIPMILPSolver_Conhdlr::scip_lock)
{  /*lint --e{715}*/

 std::vector< SCIP_VAR * > scip_vars;
 int nvars;

 // get variables data from the scipmilpsolver 
 /** NOTE: This is essential, because using SCIPMILPSolver we are adding auxiliary 
 * variables which don't have to be considered in the computation of the best 
 * solution. */
 scip_vars =  parent_scipmilpsolver->get_SCIP_var();
 nvars = scip_vars.size();

 for( int i = 0; i < nvars; i++){
      SCIP_CALL( SCIPaddVarLocksType(scip, scip_vars[i], locktype, nlockspos + nlocksneg, nlockspos + nlocksneg) );
   }

 return SCIP_OKAY;
}

/** transforms constraint data into data belonging to the transformed problem */
SCIP_DECL_CONSTRANS(SCIPMILPSolver_Conhdlr::scip_trans){
   SCIP_CONSDATA* sourcedata;
   SCIP_CONSDATA* targetdata = NULL;

   sourcedata = SCIPconsGetData(sourcecons);
   assert( sourcedata != NULL );

   SCIP_CALL( SCIPallocBlockMemory(scip, &targetdata) );
   targetdata->new_cut = FALSE;

   /* create target constraint */
   SCIP_CALL( SCIPcreateCons(scip, targetcons, SCIPconsGetName(sourcecons), conshdlr, targetdata,
         SCIPconsIsInitial(sourcecons), SCIPconsIsSeparated(sourcecons), SCIPconsIsEnforced(sourcecons),
         SCIPconsIsChecked(sourcecons), SCIPconsIsPropagated(sourcecons),  SCIPconsIsLocal(sourcecons),
         SCIPconsIsModifiable(sourcecons), SCIPconsIsDynamic(sourcecons), SCIPconsIsRemovable(sourcecons),
         SCIPconsIsStickingAtNode(sourcecons)) );

   return SCIP_OKAY;
}

/** frees specific constraint data */
SCIP_DECL_CONSDELETE(SCIPMILPSolver_Conhdlr::scip_delete){  /*lint --e{715}*/
   
   assert(consdata != NULL);
   SCIPfreeBlockMemory(scip, consdata);

   return SCIP_OKAY;
}


/** creates and captures a lazy constraint with all its constraint 
 *  flags set to their default values */
SCIP_RETCODE SMSpp_di_unipi_it::SCIPcreateSCIPMILPSolver_basiccb(
   SCIP*        scip,               /**< SCIP data structure */
   SCIP_CONS**  cons,               /**< pointer to hold the created constraint */
   const char*  name,               /**< name of constraint */
   std::vector< SCIP_VAR * > vars   /**< SCIP vars */
   ){

   SCIP_CALL( SCIPcreateSCIPMILPSolver_cb(scip, cons, name, vars ,
         FALSE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, TRUE) );

   return SCIP_OKAY;
}

/** creates and captures a constraint used which will be used as a separator */
 SCIP_RETCODE SMSpp_di_unipi_it::SCIPcreateSCIPMILPSolver_cb(
   SCIP*        scip,               /**< SCIP data structure */
   SCIP_CONS**  cons,              /**< pointer to hold the created 
                                        constraint */
   const char*  name,               /**< name of constraint */
   std::vector< SCIP_VAR * > vars,   /**< SCIP vars */
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
   ){
   
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSDATA* consdata;
   int nvars;

   /* find the subtour constraint handler */
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   if( conshdlr == NULL )
   {
      SCIPerrorMessage("scipmilpsolver constraint handler not found\n");
      return SCIP_PLUGINNOTFOUND;
   }

   /* create constraint data */
   SCIP_CALL( SCIPallocBlockMemory(scip, &consdata) ); /*lint !e530*/

   consdata->new_cut = FALSE;

   /* create constraint */
   SCIP_CALL( SCIPcreateCons(scip, cons, name, conshdlr, consdata, initial, 
         separate, enforce, check, propagate, local, modifiable, dynamic, 
         removable, FALSE) );

   return SCIP_OKAY;
}

/*--------------------------------------------------------------------------*/
/*----------------- End File SCIPMILPSolverConhdlr.cpp ---------------------*/
/*--------------------------------------------------------------------------*/