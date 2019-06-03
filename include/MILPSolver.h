/*--------------------------------------------------------------------------*/
/*---------------------------- File MILPSolver.h ---------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the MILP Solver class, a base class for any other
 * General Purpose Solvers that are able to tackle a MILP problem.
 * It's main purpose is to describe the LP through a series of vectors that
 * can be used by other solvers like CPLEX.
 *
 * \version 0.20
 *
 * \date 31 - 05 - 2019
 *
 * \author Antonio Frangioni \n
 *         Operations Research Group \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \author Niccolò Iardella \n
 *         Operations Research Group \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \author Kostas Tavlaridis-Gyparakis \n
 *         Operations Research Group \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * Copyright &copy by Antonio Frangioni, Kostas Tavlaridis-Gyparakis, Niccolò Iardella
 */
/*--------------------------------------------------------------------------*/
/*----------------------------- DEFINITIONS --------------------------------*/
/*--------------------------------------------------------------------------*/

#ifndef __MILPSolver
#define __MILPSolver

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include <SMSTypedefs.h>
#include <Observer.h>
#include <Block.h>
#include <Solver.h>
#include <ColVariable.h>
#include <FRealObjective.h>
#include <FRowConstraint.h>
#include <OneVarConstraint.h>

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

/// derived class for describing MILP Problems
/**
 * The MILPSolver class derives from Solver [see Solver.h] and extends the
 * interface of the base class to be able to efficiently handle MILP Problems.
 * This class alone does not solve problems, but it serves as base class to
 * other MILP Solver like CPXMILPSolver. Nonetheless it can be used by itself
 * wherever a description of the LP in matricial form is needed.
 *
 * The MILPSolver can be registered to any kind of Block, assuming that
 * it contains a MILP formulation, and it generates a collection of vectors
 * that describes the LP and can be used by an external framework, e.g., CPLEX.
 * Moreover, MILPSolver provides an interface that derived classes can
 * implement to process modifications that may have occured to the Block.
 *
 * The first thing that one has to take care is the correspondance between
 * the constraints and variables of the Block and the constraint matrix.
 * The proper connection between the two is built via the method of set_Block().
 * In this method a Breadth First Search is conducted, that scans the Block and
 * all its children, if any, populating the vectors needed to define the LP.
 *
 * Accordingly to the Block structure we build the constraint matrix
 * in two steps:
 * First, the static part, i.e. the part of the rows (constraints) and
 * columns (variables) that is not going to be deleted for all the lifecycle
 * of the optimization problem.
 * Second, the dynamic part, i.e. the part of the rows (constraints) and
 * columns (variables) that can potentially change or stop existing.
 *
 * After set_Block() has been called the constraint matrix is assumed to have
 * all the information needed from the Block to solve the problem.
 * The information can be retrieved by a library of getters.
 * The methods compute(), get_var_solution() and new_var_solution(), derived
 * from the Solver class, do nothing and must be impelmented by derived classes.
 */
class MILPSolver : public Solver {

/*--------------------------------------------------------------------------*/
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

 public:

/*--------------------------------------------------------------------------*/
/*--------------------- CONSTRUCTOR AND DESTRUCTOR -------------------------*/
/*--------------------------------------------------------------------------*/

/** @name Constructor and Destructor
 *  @{ */

 MILPSolver();

 ~MILPSolver() override;

/*@}------------------------------------------------------------------------*/
/*--------------------- PUBLIC METHODS OF THE CLASS ------------------------*/
/*--------------------------------------------------------------------------*/

/** @name Public Methods
 *  @{ */

 /// Getter for numcols
 int get_numcols() const;

 /// Getter for numrows
 int get_numrows() const;

 /// Getter for nzelements
 int get_nzelements() const;

 /// Getter for objsense
 int get_objsense() const;

 /// Getter for objective
 const std::vector<double>& get_objective() const;

 /// Getter for q_objective
 const std::vector<double>& get_q_objective() const;

 /// Getter for rhs
 const std::vector<double>& get_rhs() const;

 /// Getter for sense
 const std::vector<char>& get_sense() const;

 /// Getter for matbeg
 const std::vector<int>& get_matbeg() const;

 /// Getter for matcnt
 const std::vector<int>& get_matcnt() const;

 /// Getter for matind
 const std::vector<int>& get_matind() const;

 /// Getter for matval
 const std::vector<double>& get_matval() const;

 /// Getter for lb
 const std::vector<double>& get_lb() const;

 /// Getter for ub
 const std::vector<double>& get_ub() const;

 /// Getter for xctype
 const std::vector<char>& get_xctype() const;

/*@}------------------------------------------------------------------------*/
/*--------------------- DERIVED METHODS OF BASE CLASS ----------------------*/
/*--------------------------------------------------------------------------*/

/** @name Public Methods derived of the Base Class
 *  @{ */

 int compute(bool changedvars) override { return 0; }

 void set_Block(Block* block) override;

 void get_var_solution(Configuration* solc) override {}

 bool new_var_solution() override { return (false); }

/*@}------------------------------------------------------------------------*/
/*-------------------- PROTECTED FIELDS OF THE CLASS -----------------------*/
/*--------------------------------------------------------------------------*/

/** @name Protected Fields of the MILPSolver
 *  @{ */

protected:

/*--------------------------------------------------------------------------*/
/*---------------- VARIABLE AND CONSTRAINT TRACKING VECTORS ----------------*/
/*--------------------------------------------------------------------------*/

 typedef std::pair<ColVariable*, int> var_int;
 typedef std::pair<int, ColVariable*> int_var;
 typedef std::pair<FRowConstraint*, int> const_int;
 typedef std::pair<int, FRowConstraint*> int_const;

 /**
  * The following vectors are used in order to keep track between the
  * variables and constraints of the Block and the constraint matrix.
  *
  *  - v_s_var_int, v_s_const_int : vectors of pairs that store the address
  *    of the first element of each different type of static variable and
  *    constraint, respectively, and the corresponding index in constraint
  *    matrix. Vectors are kept sorted in ascending order by address.
  *
  *  - v_d_var_int, v_d_const_int : vectors of pairs that store the addresses
  *    of all the dynamic variables and constraints respectively and the
  *    corresponding index in constraint matrix. and are stored in an ascending
  *    Vectors are kept sorted in ascending order by address.
  *
  *  - v_int_s_var, v_int_s_const : vectors of pairs that store the indices of
  *    columns and rows, respectively, of the constraint matrix and the
  *    address of the corresponding (group of) static variables and constraints.
  *    Vectors are kept sorted in ascending order by index.
  *
  *  - v_int_d_var, v_int_d_const : vectors of pairs that store the indices of
  *    columns and rows, respectively, of the constraint matrix and the
  *    address of the corresponding dynamic variables and constraints.
  *    Vectors are kept sorted in ascending order by index.
  *
  * Using these vectors of pair we can at any time locate the index
  * of each constraint and variable within the constraint matrix, and viceversa.
  * Note that the Block stores static variables and constraints grouped, so
  * the vectors that store addresses of static elements (v_s_var_int,
  * v_int_s_var, v_s_const_int, v_int_s_const) contain the address of the first
  * element of each group.
  * The methods that use these vectors need to be aware of that.
  */

 std::vector<var_int> v_s_var_int; ///< From static variable to index
 std::vector<int_var> v_int_s_var; ///< From index to static variable

 std::vector<const_int> v_s_const_int; ///< From static constraint to index
 std::vector<int_const> v_int_s_const; ///< From index to static constraint

 std::vector<var_int> v_d_var_int; ///< From dynamic variable to index
 std::vector<int_var> v_int_d_var; ///< From index to dynamic variable

 std::vector<const_int> v_d_const_int; ///< From dynamic constraint to index
 std::vector<int_const> v_int_d_const; ///< From index to dynamic constraint

 /**
  * The following two vectors exist because when building and updating
  * the problem we treat differently actual constraints and bounds of
  * the variables so it is useful to keep track of their addresses.
  */

 std::vector<std::vector<FRowConstraint*> > active_row_constraints;
 std::vector<std::vector<OneVarConstraint*> > active_box_constraints;

/*--------------------------------------------------------------------------*/
/*--------------------- FIELDS FOR PROBLEM DESCRIPTION ---------------------*/
/*--------------------------------------------------------------------------*/

 /**
  * The following variables and vectors are used to describe the LP.
  * The format is the one required by CPLEX Callable Library, but we decided
  * to maintain them in a class that does not depend on CPLEX as other solvers
  * may require the same description.
  */

 /**
  * An integer that specifies the number of columns in the constraint matrix,
  * or equivalently, the number of variables in the problem object.
  */
 int numcols{};

 /**
  * An integer that specifies the number of rows in the constraint matrix,
  * not including the objective function or bounds on the variables.
  */
 int numrows{};

 /**
  * An integer that specifies the number of nonzero coefficients in the
  * constraint matrix, used to properly size matind and matval.
  */
 int nzelements{};
 
 /**
 * An integer that specifies whether the problem is a minimization or
 * maximization problem.
 */
 int objsense{};

 /**
  * An array of length at least numcols containing the objective function
  * coefficients.
  */
 std::vector<double> objective;

 /**
  * An array of length numcolrs containing the quadratic coefficients of
  * the separable quadratic objective.
  */
 std::vector<double> q_objective;

 /**
 * An array of length at least numrows containing the righthand side value
 * for each constraint in the constraint matrix.
 */
 std::vector<double> rhs;

 /**
  * An array of length at least numrows containing the sense of each constraint
  * in the constraint matrix.
  */
 std::vector<char> sense;

 /**
 * The following vectors define the constraint matrix.
 *
 * CPLEX needs to know only the nonzero coefficients. These are grouped by
 * column in the array matval. The nonzero elements of every column must be
 * stored in sequential locations in this array with matbeg[j] containing the
 * index of the beginning of column j and matcnt[j] containing the number of
 * entries in column j. The components of matbeg must be in ascending order.
 * For each k, matind[k] specifies the row number of the corresponding
 * coefficient, matval[k].
 */
 std::vector<int> matbeg;
 std::vector<int> matcnt;
 std::vector<int> matind;
 std::vector<double> matval;
 
 /**
  * An array of length at least numcols containing the lower bound on each
  * of the variables.
  */
 std::vector<double> lb;

 /**
  * An array of length at least numcols containing the upper bound on each
  * of the variables.
  */
 std::vector<double> ub;

 /**
  * An array of length numcols containing the type of each column in
  * the constraint matrix. Possible values:
  * CPX_CONTINUOUS 'C' continuous variable
  * CPX_BINARY     'B' binary variable
  * CPX_INTEGER    'I' general integer variable
  * CPX_SEMICONT   'S' semi-continuous variable
  * CPX_SEMIINT    'N' semi-integer variable
  */
 std::vector<char> xctype;

 int sol_status{};     ///< Solution status (OK, Infeasible, Unbounded, ...)
 int nodes{};          ///< Number of nodes used to solve the problem
 

 /* TODO: Support CPLEX parameters
  * The following fields are for used in set_var() functions
  * I have to figure out what to use in new version of SMS++
  */
 // double f_max_time{};   ///< maximum time for each call to solve()
 // int f_max_iter{};      ///< maximum iterations in each call to solve()
 // int f_log_verb{};      ///< "verbosity" of the log
 // double f_rel_acc{};    ///< relative objective function accuracy
 // double f_abs_acc{};    ///< absolute objective function accuracy
 // double f_up_cutoff{};  ///< upper cutoff
 // double f_lw_cutoff{};  ///< lower cutoff
 // int f_max_sol{};       ///< max number of solutions for each call to solve()
 // double f_r_acc_sol{};  ///< max relative error of a solution
 // double f_a_acc_sol{};  ///< max absolute error of a solution
 // double f_f_acc_sol{};  ///< max relative constraint violation of a solution

/*--------------------------------------------------------------------------*/

 /**
  * The following methods use the tracking vectors to get the indices of the
  * Variables/Constraints from the pointers and viceversa.
  *
  * FIXME: The following methods only look in the static part
  * I think that it's because index search is rarely done for the dynamic part.
  * For the moment it works well and it's not worth changing it,
  * variable_with_index() and constraint_with_index() are not even used.
  */

 /**
  * It returns the constraint matrix column index of a given variable.
  *
  * @param p_var a pointer to a ColVariable
  * @return the corresponding constraint matrix column index
  */
 int index_of_variable(ColVariable* p_var);

 /**
  * It returns the constraint matrix row index of the given constraint.
  *
  * @param p_const a pointer to a FRowConstraint
  * @return the corresponding constraint matrix row index
  */
 int index_of_constraint(FRowConstraint* p_const);

 /**
  * It returns the variable corresponding to the given constraint matrix index.
  * @param i a constraint matrix column index
  * @return a pointer to the corresponding ColVariable
  */
 ColVariable* variable_with_index(int i);

 /**
 * It returns the constraint corresponding to the given constraint matrix index.
 * @param i a constraint matrix row index
 * @return a pointer to the corresponding FRowConstraint
 */
 FRowConstraint* constraint_with_index(int i);

/*@}------------------------------------------------------------------------*/
/*--------------------- PRIVATE FIELDS OF THE CLASS ------------------------*/
/*--------------------------------------------------------------------------*/

 private:

/*--------------------------------------------------------------------------*/
/*------------- AUXILIARY METHODS FOR POPULATING THE PROBLEM  --------------*/
/*--------------------------------------------------------------------------*/

/**
 * @name Private auxiliary methods for reading the data from the Block
 *  @{ */

 /**
  * All of these auxiliary methods, except scan_objective(), take integer
  * counters as input parameters. That's because they are meant to be used
  * by un_any_const_static() and un_any_const_dynamic() template functions
  * on boost::any containers, and the counters keep track of the elements
  * inside the containers.
  */

 /**
  * Accumulating procedure that helps getting the total number of rows.
  * If the constraint is linear, n_rows is incremented.
  *
  * @param constraint FRowConstraint to be counted
  * @param n_rows accumulator for rows
  */
 void count_constraints(FRowConstraint& constraint, int& n_rows);

 /**
  * Accumulating procedure that helps getting the total number of columns.
  *
  * @param variable ColVariable to be counted
  * @param n_cols accumulator for columns
  */
 void count_variables(ColVariable& variable, int& n_cols);

 /**
  * Accumulating procedure that helps getting the total number of nonzero
  * elements (that is, the coefficients of the matrix).
  * Since it requires getting all the active stuff for each ColVariable and
  * checking if they are constraints or bounds or else, we use this procedure
  * to populate the active_row_constraints and active_box_constraints.
  *
  * @param variable ColVariable to be checked for active constraints
  * @param nz_elements accumulator for nonzero elements
  * @param cnt counter that keeps track of variable index
  */
 void count_nzelements(ColVariable& variable, int& nz_elements, int& cnt);

 /**
  * It scans a static ColVariable or a group of static ColVariables and fills
  * the vectors of the LP accordingly.
  *
  * @param var a reference to a [vector of] static ColVariable[s]
  * @param first an counter that should be 0 when var is the first
  *              element of a vector of ColVariables.
  * @param i a counter for variables/columns shared w/ scan_dynamic_variable()
  */
 void scan_static_variable(ColVariable& var, int& first, int& i);

 /**
  * It scans a dynamic ColVariable and fills the vectors of the LP accordingly.
  *
  * @param var a reference to a dynamic ColVariable
  * @param i a counter for variables/columns shared with scan_static_variable()
  */
 void scan_dynamic_variable(ColVariable& lvar, int& i);

 /**
  * It scans a static FRowConstraint or a group of static FRowConstraints and
  * fills the vectors of the LP accordingly.
  * @param lconst a reference to a [vector of] static FRowConstraint[s]
  * @param first an counter that should be 0 when lconst is the first
  *              element of a vector of FRowConstraints.
  * @param i a counter for constraints/rows shared w/ scan_dynamic_constraint()
  */
 void scan_static_constraint(FRowConstraint& lconst, int& first, int& i);

 /**
  * It scans a dyn FRowConstraint and fills the vectors of the LP accordingly.
  *
  * @param lconst a reference to a dynamic FRowConstraint
  * @param i a counter for constraints/rows shared w/ scan_static_constraint()
  */
 void scan_dynamic_constraint(FRowConstraint& lconst, int& i);

 /**
  * It scans a FRealObjective and fills the vectors of the LP accordingly.
  * @param obj a FRealObjective
  */
 void scan_objective(const FRealObjective* obj);

 /**
  * It clears all the LP vectors
  */
 void clear_matrices();

/*@}------------------------------------------------------------------------*/

 // TODO: Support CPLEX parameters
 // void set_par(int par, int value);
 // void set_par(int par, double value);
 // void set_par(int par, long value);

 // std::vector<int> indexed;

 // TODO: Support cuts callback
 // int cuts{}; //parameter for adding cuts
 // int Callback(CPXCENVptr env, void* cbdata, int wherefrom, int* useraction_p);

 SMSpp_insert_in_factory_h;
};   // end( class MILPSolver )

}; // end( namespace SMSpp_di_unipi_it )

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

#endif

