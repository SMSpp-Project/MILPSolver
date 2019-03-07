/*--------------------------------------------------------------------------*/
/*---------------------------- File MILPSolver.h ---------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for MILP Solver that defines a General Purpose Solver that is
 * able to tackle a MILP Problem expressed by a Block via the use of CLPEX
 *
 *
 * \version 0.10
 *
 * \date 22 - 12 - 2016
 *
 * \author Antonio Frangioni \n
 *         Operations Research Group \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \author Kostas Tavlaridis-Gyparakis \n
 *         Operations Research Group \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * Copyright &copy by Antonio Frangioni, Kostas Tavlaridis-Gyparakis
 */
/*--------------------------------------------------------------------------*/
/*----------------------------- DEFINITIONS --------------------------------*/
/*--------------------------------------------------------------------------*/

#ifndef __MILPSolver
#define __MILPSolver
/* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

// #include "SMSTypedefs.h"
// #include "ColVariable.h"
//#include "Observer.h"
#include "FRowConstraint.h"
#include "FRealObjective.h"
//#include "DQuadObjectiveFunction.h"
//#include "Block.h"
// #include "Solver.h"
//!!#include <ilcplex/cplex.h>
//!!#include <ilcplex/cplexcheck.h>
#include <cplex.h>
  //#include <cplexcheck.h>

/*--------------------------------------------------------------------------*/
/*----------------------------- NAMESPACE ----------------------------------*/
/*--------------------------------------------------------------------------*/

/// namespace for the Structured Modeling System++ (SMS++)
namespace SMSpp_di_unipi_it {

class Block;  ///< forward definition of class Block

/*--------------------------------------------------------------------------*/
/*------------------------- CLASS MILPSolver -------------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------- GENERAL NOTES --------------------------------*/
/*--------------------------------------------------------------------------*/

/// derived class for solving MILP Problems via CPLEX
/** The MILPSolver class derives from Solver [see Solver.h] and extends the 
 * interface of the base class to be able to efficiently handle and solve
 * MILP Problems. The main funcionalities of the MILPSolver can be summarised 
 * as follows. The MILPSolver is capable of receiving any kind of Block, ta-
 * king as assumption that this exceeds a MILP formulation ofcourse, and is
 * responsible to pass this Block to CPLEX and then retrieve back the corres-
 * ponding solution. On top of that MILPSolver is capable of upgrading the 
 * CPLEX Solver for any possible modifications that may have occured to the
 * examined Block.
 *
 * More specifically in order for the above to happen the first fundamendal
 * thing that one has to take care is the connectivity between the Block and
 * the CPLEX coefficient Matrix. The proper connection between the two is be-
 * ing constructed via the method of set_Block(). In this method a Breadth 
 * First Search is conducted, that scans the received Block and all of each 
 * corresponding children if any, in order to populate the CPLEX data needed 
 * to define the corresponding MILP. 
 *
 * Following the Block structure the CPLEX coefficient matrix is being devided
 * in two parts:
 * - The static, i.e. the part of the rows (constraints) and columns (variab-
 *   les) that is not going to be deleted for all the optimisation problem,
 *   and is being processed first.
 * - The dynamic, i.e. the part of the rows (constraints) and columns (variab-
 *   les) that expected to be potentially stop existing throughout the optimi-
 *   sation problem. And in this part is also expected for new rows and/or co-
 *   lumns to be added. This part is being processed second.
 * 
 * An other fundamendal element that is used in order to establish the conne-
 * ctivity between the Block and the CPLEX coeff. matrix are the vectors of
 * pairs that keep track between the variables and constraints of the Blocks
 * and their corresponding indexes in the CPLEX coeff. matrix. 
 *
 * For the static part this achieved with the introduction of the following
 * elements:
 *
 * -v_s_var_int,v_s_const_int : vectors of pair that store the adress of the first
 *  element of each different type of variavle and constraint respectively and
 *  the corresponding index in CPLEX coeff matrix and are stored in an ascending
 *  order based on the adresses.
 *
 * -v_int_s_var,v_int_s_const : vectors of pair that store the adress of the first
 *  element of each different type of variavle and constraint respectively and
 *  the corresponding index in CPLEX coeff matrix and are stored in an ascending
 *  order based on the indexes of CPLEX.
 *
 * Via the usage of these vectors of pair we can at any time locate the position
 * of each individual static constraint and variable with (probably) the minimum
 * calculation cost, via the usage of the corresponding methods.
 *
 * For the dynamic part this achieved with the introduction of the following
 * elements:
 *
 * -v_d_var_int,v_d_const_int : vectors of pair that store the adresses of all 
 *  elements of each different type of variavle and constraint respectively and
 *  the corresponding index in CPLEX coeff matrix and are stored in an ascending
 *  order based on the adresses.
 *
 * -v_int_d_var,v_int_d_const : vectors of pair that store the adresses of all
 *  elements of each different type of variavle and constraint respectively and
 *  the corresponding index in CPLEX coeff matrix and are stored in an ascending
 *  order based on the indexes of CPLEX.
 *
 * Via the usage of these vectors of pair we can at any time locate the position
 * of each individual static constraint and variable. Note that here the size of
 * the vectors is bigger since the dynamic variables and constraints are stored 
 * in lists and can at any time appear and dissapear from the problem and as a
 * result we can not treat them in the same fashion as the static variables and
 * constraint.
 *
 * After this point the CPLEX coeff matrix is assumed to have all the information
 * needed from the Block in order to tackle down the optimisation MILP problem
 * that has received. This is done via the method solve() that is followed by
 * the method get_var_solution(), which is responsible to write back the obtained
 * solution to the Block structure.
 *
 * The next thing that should be mentioned is the handling of the modulations in
 * the MILP Solver. Initially each new modification is being processed in the me-
 * thod of add_modification(), that is responsible to proparly direct the modifi-
 * cation to the corresponding method based on its nature, where one can have:
 * 
 * -Variable Modification: Is refered to the static part of the problem and is be-
 *  ing processed via the use of the method var_modification(). Where the method
 *  is responsible for communicating the variable's modification to the CPLEX coeff
 *  matrix and proceeding with the corresponding needed changes on it. 
 *
 * -Constraint Modification: Is refered to the static part of the problem and is
 *  being processed via the use of the method const_modification(). Where the met-
 *  hod is responsible for communicating the constraint's modification to the CPLEX 
 *  coeff matrix and proceeding with the corresponding needed changes on it.
 *
 *
 * -Objective Function Modification: Is refered to the static part of the problem 
 *  and is being processed via the use of the method of_modification(). Where the
 *  method is responsible for communicating the o.f. modification to the CPLEX 
 *  coeff matrix and proceeding with the corresponding needed changes on it.
 *
 * -Dynamic Modification: Is refered as the title says to the dynamic part of the
 *  problem and more specifically in the addition or removal of a dynamic Constra-
 *  int or Variable. These modifications are handled via the method of 
 *  dynamic_modification() that is responsible of receiving the modification and
 *  for all the corresponding Constraints/Variables taking care to add/remove them
 *  from the CPLEX coeff matrix and upgrade equivelantely the corresponding vectors
 *  of pairs.
 *
 */

class MILPSolver : public Solver {

 public:

/*--------------------------------------------------------------------------*/
/*---------------------- PUBLIC TYPES OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

 typedef std::pair<ColVariable *, int> var_int;

 typedef std::pair<int, ColVariable *> int_var;

 typedef std::pair<FRowConstraint *, int> const_int;

 typedef std::pair<int, FRowConstraint *> int_const;

/*--------------------------------------------------------------------------*/
/*--------------------- CONSTRUCTOR AND DESTRUCTOR -------------------------*/
/*--------------------------------------------------------------------------*/

/** @name Constructor and Destructor
 *  @{ */

 MILPSolver() : Solver() {
  env = nullptr;
  milp = nullptr;
 }

 ~MILPSolver() override {
  if( env ) {
   CPXfreeprob( env, &milp );
   CPXcloseCPLEX( &env );
  }
 }

/*@}------------------------------------------------------------------------*/
/*--------------------- PUBLIC METHODS OF THE CLASS ------------------------*/
/*--------------------------------------------------------------------------*/

/*@}------------------------------------------------------------------------*/
/*--------------------- DERIVED METHODS OF BASE CLASS ----------------------*/
/*--------------------------------------------------------------------------*/

/** @name Public Methods derived of the Base Class
 *  @{ */
 void set_Block( Block *block ) override;
 ///< method for setting the Block and building the corresponding CPLEX problem

 int solve();
 ///< method for solving the problem with CPLEX

 void get_var_solution( Configuration *solc) override;
 ///< method for writing the solution in Block

 bool new_var_solution() override { return ( false ); }

/*@}------------------------------------------------------------------------*/
/*----------------------- METHODS FOR READING DATA  ------------------------*/
/*--------------------------------------------------------------------------*/

/** @name Public Methods for reading the data of the Block
 *  @{ */

 void scan_static_constraints( FRowConstraint &lconst,
                               char *sense,
                               double *rhs,
                               int &first, int &i );
 /***< method used to scan all the data of the each individual static constraint and
    stores them in the appropriate vectors in order to pass them to CPLEX
 */

 void scan_dynamic_constraints( FRowConstraint &lconst,
                                char *sense,
                                double *rhs,
                                int &i );
 /***< method used to scan all the data of the each individual dynamic constraint and
    stores them in the appropriate vectors in order to pass them to CPLEX
 */

 // TODO Fix this
 void scan_static_variables( ColVariable &lvar,
                             int *matbeg,
                             int *matcnt,
                             int *matind,
                             double *matval,
                             // double *objective,
                             // double *q_objective,
                             double *lb,
                             double *ub,
                             char *xctype,
                             int &first,
                             int &i );
 /***< method used to scan all the data of the each individual static variable and
    stores them in the appropriate vectors in order to pass them to CPLEX
 */

 void scan_dynamic_variables( ColVariable &lvar,
                              int *matbeg,
                              int *matcnt,
                              int *matind,
                              double *matval,
                              double *objective,
                              double *q_objective,
                              double *lb,
                              double *ub,
                              char *xctype,
                              int &i );
 /***< method used to scan all the data of the each individual dynamic variable and
    stores them in the appropriate vectors in order to pass them to CPLEX
 */

 void scan_objective( const FRealObjective *obj, double *objective, double *q_objective );

/*@}------------------------------------------------------------------------*/
/*-------------------- METHODS FOR MODIFYING THE PROBLEM -------------------*/
/*--------------------------------------------------------------------------*/

/** @name Public Methods for modifying the constructed CPLEX Problem
 *  @{ */


 // void add_modifications( sp_Mod &mod );
 // ///< method for adding and handling a Modification of the Problem
 //
 // void var_modification( VariableMod *mod );
 // ///< method for adding and handling a Variable Modification
 //
 // void of_modification( ObjectiveMod *mod );
 // ///< method for adding and handling an Objective Function Modification
 //
 // void const_modification( ConstraintMod *mod );
 // ///< method for adding and handling a Constraint Modification
 //
 // void dynamic_modification( BlockModAD *mod );
 // ///< method for handling a dynamic Modification
 //
 // void add_dynamic_constraint( FRowConstraint *r_const );
 // ///< method for adding a single new dynamic constraint to CPLEX
 //
 // void add_dynamic_variable( ColVariable *r_var );
 // ///< method for adding a single new dynamic variable to CPLEX
 //
 // void remove_dynamic_constraint( FRowConstraint *r_const );
 // ///< method for deleting a single dynamic constraint to CPLEX
 //
 // void remove_dynamic_variable( ColVariable *r_var );
 // ///< method for removing a single dynamic variable to CPLEX

/*@}------------------------------------------------------------------------*/
/*-------------------------- SUPPLEMENTARY METHODS  ------------------------*/
/*--------------------------------------------------------------------------*/

/** @name Supplementary Methods that are used in order to assist the above methods
 *  @{ */

 void count_constraints( FRowConstraint &lconst, int &count );
///< method used to count the total number of rows

 void count_variables( ColVariable &lvar, int &count, int &count2 );
 ///< method used to count the total number of columns and non-zero elements

 void set_var_value( ColVariable &lvar, double *tmpx, int &i );
///< method used to pass the solution to the variables

 int ind_var( ColVariable *p_var );
///< method for returning the CPLEX coeff-matrix index of the examined Variable

 int ind_const( FRowConstraint *p_const );
///< method for returning the CPLEX coeff-matrix index of the examined Constraint

 ColVariable *var_ind( int i );
///< method for returning the Variable of the examined CPLEX coeff-matrix index

 FRowConstraint *const_ind( int i );
///< method for returning the Constraint of the examined CPLEX coeff-matrix index

 OFValue get_lb() override;

 OFValue get_ub() override;

 void set_par( int par, int value );

 void set_par( int par, double value );

 // void set_par( int par, long value );


 int sol_status;
 int nodes;
 int obj_type;
 std::vector<int> indexed;

 int cuts; //parameter for adding cuts

 int Callback( CPXCENVptr env, void *cbdata, int wherefrom, int *useraction_p );

/*@}------------------------------------------------------------------------*/
/*-------------------- PROTECTED FIELDS OF THE CLASS -----------------------*/
/*--------------------------------------------------------------------------*/
/** @name Protected Fields of the MILPSolver
 *  @{ */

protected: 

/*** The following vectors are used in order to keep track between the static vari-
  ables and constraints of the Block and the CPLEX coeff matrix. Where we have:

  -v_s_var_int,v_s_const_int : vectors of pair that store the adress of the first
   element of each different type of variavle and constraint respectively and
   the corresponding index in CPLEX coeff matrix and are stored in an ascending
   order based on the adresses.

  -v_int_s_var,v_int_s_const : vectors of pair that store the adress of the first
   element of each different type of variavle and constraint respectively and
   the corresponding index in CPLEX coeff matrix and are stored in an ascending
   order based on the indexes of CPLEX.

  Via the usage of these vectors of pair we can at any time locate the position
  of each individual static constraint and variable with (probably) the minimum
  calculation cost, via the usage of the corresponding methods.

*/

 std::vector<var_int> v_s_var_int;
 std::vector<int_var> v_int_s_var;

 std::vector<const_int> v_s_const_int;
 std::vector<int_const> v_int_s_const;

/*** The following vectors are used in order to keep track between the dynamic 
  variables and constraints of the Block and the CPLEX coeff matrix. Where we have:

  -v_d_var_int,v_d_const_int : vectors of pair that store the adresses of all 
   elements of each different type of variavle and constraint respectively and
   the corresponding index in CPLEX coeff matrix and are stored in an ascending
   order based on the adresses.

  -v_int_d_var,v_int_d_const : vectors of pair that store the adresses of all
   elements of each different type of variavle and constraint respectively and
   the corresponding index in CPLEX coeff matrix and are stored in an ascending
   order based on the indexes of CPLEX.

  Via the usage of these vectors of pair we can at any time locate the position
  of each individual static constraint and variable. Note that here the size of
  the vectors is bigger since the dynamic variables and constraints are stored 
  in lists and can at any time appear and dissapear from the problem and as a
  result we can not treat them in the same fashion as the static variables and
  constraint.
*/

 std::vector<var_int> v_d_var_int;
 std::vector<int_var> v_int_d_var;

 std::vector<const_int> v_d_const_int;
 std::vector<int_const> v_int_d_const;

/***Two vectors that store the initial index in CPLEX and the initial adress of 
each different type of variable keeping a row-wise order with respect with the
variable types, in the sense that it will first scan the variables found in the
block and then will proceed with the variables in each different sublock. 
Knowing the exact order gives us the advantage of being able to read all the
parts of the solution in the fastest possible way. */

 std::vector<int> v_ind;
 std::vector<ColVariable *> v_pt;

 // TODO Check if already available typedef
 std::vector<std::vector<FRowConstraint *> > p_active_constraints;

 // TODO The following fields are for supporting set_var() I have to figure out what to use in new version of SMS++
 double f_max_time;    ///< maximum time for each call to solve()
 int f_max_iter;       ///< maximum iterations in each call to solve()
 int f_log_verb;       ///< "verbosity" of the log
 double f_rel_acc;     ///< relative objective function accuracy
 double f_abs_acc;    ///< absolute objective function accuracy
 double f_up_cutoff;  ///< upper cutoff
 double f_lw_cutoff;  ///< lower cutoff
 int f_max_sol;        ///< max number of solutions for each call to solve()
 double f_r_acc_sol;  ///< max relative error of a solution
 double f_a_acc_sol;  ///< max absolute error of a solution
 double f_f_acc_sol;   ///< max relative constraint violation of a solution

/*** Following integer variables are used in order to store information needed to
    to populate the CPLEX Matrix
*/

 int numrows;
 int numcols;
 int nzelements;

/// Following two elements are used to initialize the CPLEX environement
 CPXENVptr env;
 CPXLPptr milp;

/*@}*/

};   // end( class MILPSolver )

}; // end( namespace SMSpp_di_unipi_it )

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

#endif  /* MILPSolver.h included */

