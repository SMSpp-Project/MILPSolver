/*--------------------------------------------------------------------------*/
/*---------------------------- File MILPSolver.h ---------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the MILPSolver class.
 *
 * The class is a base class for any other general purpose solvers
 * that are able to tackle a MILP problem. Its purpose is to describe the LP
 * problem with a series of vectors that can be used by other solvers like
 * CPLEX.
 *
 * \author Antonio Frangioni \n
 *         Operations Research Group \n
 *         Dipartimento di Informatica \n
 *         Università di Pisa \n
 *
 * \author Niccolò Iardella \n
 *         Operations Research Group \n
 *         Dipartimento di Informatica \n
 *         Università di Pisa \n
 *
 * \author Kostas Tavlaridis-Gyparakis \n
 *         Operations Research Group \n
 *         Dipartimento di Informatica \n
 *         Università di Pisa \n
 *
 * \copyright &copy; Antonio Frangioni, Kostas Tavlaridis-Gyparakis, Niccolò Iardella
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
#include <CDASolver.h>
#include <ColVariable.h>
#include <FRealObjective.h>
#include <FRowConstraint.h>
#include <OneVarConstraint.h>

/*--------------------------------------------------------------------------*/
/*----------------------------- NAMESPACE ----------------------------------*/
/*--------------------------------------------------------------------------*/

/// namespace for the Structured Modeling System++ (SMS++)
namespace SMSpp_di_unipi_it {

class Block;

/*--------------------------------------------------------------------------*/
/*------------------------- CLASS MILPSolver -------------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------- GENERAL NOTES --------------------------------*/
/*--------------------------------------------------------------------------*/

/// Class for describing MILP problems.
/**
 * The MILPSolver class derives from Solver and extends the
 * interface of the base class to be able to efficiently handle MILP problems.
 * This class alone does not solve problems, but it serves as base class to
 * other solvers like CPXMILPSolver. Nonetheless it can be used by itself
 * wherever a description of the LP problem in matricial form is needed.
 *
 * The MILPSolver can be registered to any kind of Block (assuming that
 * it contains a MILP formulation) and it generates a collection of vectors
 * that describes the LP problem and can be used by an external framework,
 * e.g., CPLEX. Moreover, MILPSolver provides an interface that derived classes
 * can implement to process modifications that may have occured to the Block.
 *
 * The main thing that this class has to take care is the correspondance between
 * the Constraints and Variables of the Block and the constraint matrix.
 * The corrispondance is built via set_Block(), that conducts a Breadth
 * First Search is conducted, scanning the Block and all its children, if any,
 * populating the vectors needed to define the LP problem.
 *
 * Following the Block structure we build the constraint matrix
 * in two steps:
 * 1. We build the static part, i.e. the part of the rows (Constraints) and
 *    columns (Variables) that is not going to be deleted for all the lifecycle
 *    of the problem.
 * 2. Then, the dynamic part, i.e. the part of the rows (Constraints) and
 *    columns (Variables) that can potentially change or stop existing.
 *
 * After set_Block() has been called, the constraint matrix is assumed to have
 * all the information needed from the Block to solve the problem.
 * The information can be retrieved by a library of getters.
 * The methods compute(), get_var_solution() and new_var_solution(), derived
 * from the Solver class, do nothing and should be implemented by derived
 * classes.
 *
 * The class defines also an interface that the derived classes should implement
 * to support modifications. The method process_modifications() is already
 * implemented and it is the one that dispatches the modifications to the other
 * methods accordingly.
 */
class MILPSolver : public CDASolver {

/*--------------------------------------------------------------------------*/
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

 public:

/*--------------------------------------------------------------------------*/
/*---------------------------- PUBLIC TYPES --------------------------------*/
/*--------------------------------------------------------------------------*/

 /** Types of integer parameters.
  * Public enum describing the different types of algorithmic parameters
  * of "int" type that the MILPSolver might have. At the moment it has none,
  * but since CPXMILPSolver has them, we add this as compatibility.
  */
 enum int_par_type_MILP {
  intLastAlgParMILP = intLastAlgPar
 };

 /** Types of string parameters.
  * Public enum describing the different types of algorithmic parameters
  * of "string" type that the MILPSolver might have. At the moment it has none,
  * but since CPXMILPSolver has them, we add this as compatibility.
  */
 enum str_par_type_MILP {
  strLastAlgParMILP = strLastAlgPar
 };

 /** Types of string parameters.
  * Public enum describing the different types of algorithmic parameters
  * of "double" type that the MILPSolver might have. At the moment it has none,
  * but since CPXMILPSolver has them, we add this as compatibility.
 */
 enum dbl_par_type_MILP {
  dblLastAlgParMILP = strLastAlgPar
 };

/*--------------------------------------------------------------------------*/
/*--------------------- CONSTRUCTOR AND DESTRUCTOR -------------------------*/
/*--------------------------------------------------------------------------*/

/**
 * @name Constructor and destructor
 * @{
 */

 MILPSolver();

 ~MILPSolver() override;
 /// @}

/*--------------------------------------------------------------------------*/
/*--------------------- PUBLIC METHODS OF THE CLASS ------------------------*/
/*--------------------------------------------------------------------------*/

/** @name Getters for the vectors of the LP problem
 *
 * @{
 */

 [[nodiscard]] int get_numcols() const;

 [[nodiscard]] int get_numrows() const;

 [[nodiscard]] int get_nzelements() const;

 [[nodiscard]] int get_objsense() const;

 [[nodiscard]] const std::vector< double > & get_objective() const;

 [[nodiscard]] const std::vector< double > & get_q_objective() const;

 [[nodiscard]] const std::vector< double > & get_rhs() const;

 [[nodiscard]] const std::vector< double > & get_rngval() const;

 [[nodiscard]] const std::vector< char > & get_sense() const;

 [[nodiscard]] const std::vector< int > & get_matbeg() const;

 [[nodiscard]] const std::vector< int > & get_matcnt() const;

 [[nodiscard]] const std::vector< int > & get_matind() const;

 [[nodiscard]] const std::vector< double > & get_matval() const;

 /// returns the lower bounds on the Variables
 [[nodiscard]] const std::vector< double > & get_var_lb() const;

 /// returns the upper bounds on the Variables
 [[nodiscard]] const std::vector< double > & get_var_ub() const;

 [[nodiscard]] const std::vector< char > & get_xctype() const;

 [[nodiscard]] const std::vector< char * > & get_rowname() const;

 [[nodiscard]] const std::vector< char * > & get_colname() const;

 [[nodiscard]] int get_nodes() const;

 [[nodiscard]] int get_num_integer_vars() const;

 /// @}

 /** @name Methods that use the correspondance vectors
  *
  * The following methods use the tracking vectors to get the indices of the
  * Variables/Constraints from the pointers and viceversa.
  *
  * We provide separate methods for looking into static, dynamic or both parts
  * of the problem, so we can reduce searching time when possible
  *
  * @{
  */

 /**
  * It returns the matrix column index of a given variable.
  *
  * @param p_var a pointer to a ColVariable
  * @return the corresponding matrix column index
  * @throws std::invalid_argument if no indices are associated to that variable
  */
 int index_of_variable( ColVariable * p_var );

 /**
  * It returns the matrix column index of a given static variable.
  *
  * @param p_var a pointer to a ColVariable
  * @return the corresponding matrix column index
  * @throws std::invalid_argument if no indices are associated to the variable
  */
 int index_of_static_variable( ColVariable * p_var );

 /**
  * It returns the matrix column index of a given dynamic variable.
  *
  * @param p_var a pointer to a ColVariable
  * @return the corresponding matrix column index
  * @throws std::invalid_argument if no indices are associated to the variable
  */
 int index_of_dynamic_variable( ColVariable * p_var );

 /**
  * It returns the matrix row index of the given constraint.
  *
  * @param p_const a pointer to a FRowConstraint
  * @return the corresponding matrix row index
  */
 int index_of_constraint( FRowConstraint * p_const );

 /**
  * It returns the matrix row index of the given static constraint.
  *
  * @param p_const a pointer to a FRowConstraint
  * @return the corresponding matrix row index
  * @throws std::invalid_argument if no indices are associated to the constraint
  */
 int index_of_static_constraint( FRowConstraint * p_const );

 /**
  * It returns the matrix row index of the given dynamic constraint.
  *
  * @param p_const a pointer to a FRowConstraint
  * @return the corresponding matrix row index
  * @throws std::invalid_argument if no indices are associated to the constraint
  */
 int index_of_dynamic_constraint( FRowConstraint * p_const );

 /**
  * It returns the variable corresponding to a variable matrix column index.
  *
  * @param i a constraint matrix column index
  * @return a pointer to the corresponding ColVariable
  * @throws std::invalid_argument if the index doesn't correspond to a variable
  */
 ColVariable * variable_with_index( int i );

 /**
  * It returns the static variable corresponding to the given variable
  * matrix column index.
  *
  * @param i a constraint matrix column index
  * @return a pointer to the corresponding static ColVariable
  * @throws std::invalid_argument if the index doesn't correspond to
  *                               a static variable
  */
 ColVariable * static_variable_with_index( int i );

 /**
  * It returns the dynamic variable corresponding to the given variable
  * matrix column index.
  *
  * @param i a constraint matrix column index
  * @return a pointer to the corresponding dynamic ColVariable
  * @throws std::invalid_argument if the index doesn't correspond to
  *                               a dynamic variable.
  */
 ColVariable * dynamic_variable_with_index( int i );

 /**
  * It returns the constraint corresponding to a constraint matrix row index.
  *
  * @param i a constraint matrix row index
  * @return a pointer to the corresponding FRowConstraint
  * @throws std::invalid_argument if the index doesn't correspond to
  *                               a constraint
  */
 FRowConstraint * constraint_with_index( int i );

 /**
  * It returns the static constraint corresponding to the given constraint
  * matrix row index.
  *
  * @param i a constraint matrix row index
  * @return a pointer to the corresponding FRowConstraint
  * @throws std::invalid_argument if the index doesn't correspond to
  *                               a static constraint
  */
 FRowConstraint * static_constraint_with_index( int i );

 /**
 * It returns the dynamic constraint corresponding to the given constraint
 * matrix row index.
 *
 * @param i a constraint matrix row index
 * @return a pointer to the corresponding FRowConstraint
 * @throws std::invalid_argument if the index doesn't correspond to
 *                               a dynamic constraint
 */
 FRowConstraint * dynamic_constraint_with_index( int i );
 /// @}

 virtual void write_lp(const std::string & filename) {}

/*--------------------------------------------------------------------------*/
/*--------------------- DERIVED METHODS OF BASE CLASS ----------------------*/
/*--------------------------------------------------------------------------*/

/**
 * @name Public Methods derived from base classes
 * @{
 */

 /// It sets the Block that the Solver has to solve and build the LP vectors.
 void set_Block( Block * block ) override;

 /// It does nothing as there is nothing to do.
 int compute( bool changedvars ) override { return 0; }

 void get_var_solution( Configuration * solc ) override {}

 void get_dual_solution( Configuration * solc ) override {}
 /// @}

/*--------------------------------------------------------------------------*/
/*-------------------- PROTECTED FIELDS OF THE CLASS -----------------------*/
/*--------------------------------------------------------------------------*/

 protected:

/*--------------------------------------------------------------------------*/
/*---------------- VARIABLE AND CONSTRAINT TRACKING VECTORS ----------------*/
/*--------------------------------------------------------------------------*/

 typedef std::pair< ColVariable *, int > var_int;
 typedef std::pair< int, ColVariable * > int_var;
 typedef std::pair< FRowConstraint *, int > const_int;
 typedef std::pair< int, FRowConstraint * > int_const;

 /** @name Correspondance vectors
  *
  * The following vectors are used in order to keep track between the
  * Variables and Constraints of the Block and the constraint matrix.
  *
  *  - #v_s_var_int, #v_s_const_int : vectors of pairs that store the address
  *    of the first element of each different type of static variable and
  *    constraint, respectively, and the corresponding index in constraint
  *    matrix. Vectors are kept sorted in ascending order by address.
  *
  *  - #v_d_var_int, #v_d_const_int : vectors of pairs that store the addresses
  *    of all the dynamic variables and constraints respectively and the
  *    corresponding index in constraint matrix. and are stored in an ascending
  *    Vectors are kept sorted in ascending order by address.
  *
  *  - #v_int_s_var, #v_int_s_const : vectors of pairs that store the indices of
  *    columns and rows, respectively, of the constraint matrix and the
  *    address of the corresponding (group of) static variables and constraints.
  *    Vectors are kept sorted in ascending order by index.
  *
  *  - #v_int_d_var, #v_int_d_const : vectors of pairs that store the indices of
  *    columns and rows, respectively, of the constraint matrix and the
  *    address of the corresponding dynamic variables and constraints.
  *    Vectors are kept sorted in ascending order by index.
  *
  * Using these vectors of pair we can at any time locate the index
  * of each constraint and variable within the constraint matrix, and viceversa.
  * Note that the Block stores static variables and constraints grouped, so
  * the vectors that store addresses of static elements (#v_s_var_int,
  * #v_int_s_var, #v_s_const_int, #v_int_s_const) contain the address of the first
  * element of each group.
  * The methods that use these vectors need to be aware of that.
  *
  * @{
  */
 std::vector< var_int > v_s_var_int; ///< From static variable to index
 std::vector< int_var > v_int_s_var; ///< From index to static variable

 std::vector< const_int > v_s_const_int; ///< From static constraint to index
 std::vector< int_const > v_int_s_const; ///< From index to static constraint

 std::vector< var_int > v_d_var_int; ///< From dynamic variable to index
 std::vector< int_var > v_int_d_var; ///< From index to dynamic variable

 std::vector< const_int > v_d_const_int; ///< From dynamic constraint to index
 std::vector< int_const > v_int_d_const; ///< From index to dynamic constraint
 /// @}

 /**
  * @name Active stuff vectors
  *
  * The following two vectors exist because when building and updating
  * the problem we need to treat differently actual constraints and bounds of
  * the variables so it is useful to keep track of their addresses.
  *
  * @{
  */

 /// Active constraints for each Variable
 std::vector< std::vector< FRowConstraint * > > active_constraints;

 /// Active bounds for each Variable
 std::vector< std::vector< OneVarConstraint * > > active_bounds;
 /// @}

/*--------------------------------------------------------------------------*/
/*--------------------- FIELDS FOR PROBLEM DESCRIPTION ---------------------*/
/*--------------------------------------------------------------------------*/

 /** @name LP problem description
  *
  * The following variables and vectors are used to describe the LP problem.
  * The format is the one required by CPLEX Callable Library, but we decided
  * to maintain them in a class that does not depend on CPLEX as other solvers
  * may require the same description.
  *
  * #matbeg, #matcnt, #matind and #matval define the constraint matrix.
  * CPLEX needs to know only the nonzero coefficients. These are grouped by
  * column in the array matval. The nonzero elements of every column must be
  * stored in sequential locations in this array with matbeg[j] containing the
  * index of the beginning of column j and matcnt[j] containing the number of
  * entries in column j. The components of matbeg must be in ascending order.
  * For each k, matind[k] specifies the row number of the corresponding
  * coefficient, matval[k].
  *
  * @{
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
 std::vector< double > objective;

 /**
  * An array of length numcols containing the quadratic coefficients of
  * the separable quadratic objective.
  */
 std::vector< double > q_objective;

 /**
  * An array of length at least numrows containing the righthand side value
  * for each constraint in the constraint matrix.
  */
 std::vector< double > rhs;

 /**
  * An array of length at least numrows containing the range value of each
  * ranged constraint. Ranged rows are those designated by 'R' in the sense
  * array. If the row is not ranged, the rngval array entry is ignored.
  * If rngval[i] > 0, then row i activity is in [rhs[i],rhs[i]+rngval[i]],
  * and if rngval[i] <= 0,then row i activity is in [rhs[i]+rngval[i],rhs[i]].
  */
 std::vector< double > rngval;

 /**
  * An array of length at least numrows containing the sense of each constraint
  * in the constraint matrix.
  */
 std::vector< char > sense;

 std::vector< int > matbeg;     ///< Beginnings of constraint matrix columns
 std::vector< int > matcnt;     ///< Sizes of constraint matrix columns
 std::vector< int > matind;     ///< Indices of rows for each coefficient
 std::vector< double > matval;  ///< All nonzero coefficients

 /**
  * An array of length at least numcols containing the lower bound on each
  * of the variables.
  */
 std::vector< double > lb;

 /**
  * An array of length at least numcols containing the upper bound on each
  * of the variables.
  */
 std::vector< double > ub;

 /**
  * An array of length numcols containing the type of each column in
  * the constraint matrix. Possible values:
  *
  * | Symbol           | Value | Type of variable |
  * | ---------------- | :---: | ---------------- |
  * | `CPX_CONTINUOUS` |  'C'  | continuous       |
  * | `CPX_BINARY`     |  'B'  | binary           |
  * | `CPX_INTEGER`    |  'I'  | general integer  |
  * | `CPX_SEMICONT`   |  'S'  | semi-continuous  |
  * | `CPX_SEMIINT`    |  'N'  | semi-integer     |
  */
 std::vector< char > xctype;

 /**
  * An array of length at least numcols containing pointers to character
  * strings containing the names of the variables.
  */
 std::vector< char* > colname;

 /**
  * An array of length at least numrows containing pointers to character
  * strings containing the names of the constraints.
  */
 std::vector< char* > rowname;

 int sol_status{};      ///< Solution status (OK, Infeasible, Unbounded, ...)
 int nodes{};           ///< Number of nodes used to solve the problem
 int int_vars = 0;      ///< Number of integer variables
 /// @}

 /** @name Clear and load the problem
  *
  * The following two methods include the main logic of the class.
  * They are called in the set_Block() method to build the problem when the
  * Solver is [re]registered to a Block, but also when a NBModification is
  * processed.
  * Note: A derived class can override these methods but should call the base
  * versions if it wants to use the LP vectors (see CPXMILPSolver).
  *
  * @{
  */

 virtual void clear_problem(); ///< It clears all the LP vectors

 virtual void load_problem();  ///< It loads all the LP vectors
 /// @}

/*--------------------------------------------------------------------------*/
/*----------------- INTERFACE FOR SUPPORTING MODIFICATIONS ---------------- */
/*--------------------------------------------------------------------------*/

/**
 * @name Methods for modifying the constructed problem
 *
 * This is the API for supporting Modifications.
 * The method process_modifications() is the one that checks the Modification
 * queue for pending modifications and calls the proper method.
 * A derived class that wants to support modifications should implement the
 * virtual methods of this group.
 *
 * @{
 */

 /// It processes all the pending modifications
 void process_modifications();

 /// It handles a variable modification
 virtual void var_modification( VariableMod * mod );

 /// It handles an objective modification
 virtual void of_modification( ObjectiveMod * mod );

 /// It handles a constraint modification
 virtual void const_modification( ConstraintMod * mod );

 /// It handles a bound modification
 virtual void bound_modification( OneVarConstraintMod * mod );

 /// It handles a function modification
 virtual void function_modification( FunctionMod * mod );

 /// It handles a function vars modification
 virtual void function_vars_modification( FunctionModVars * mod );

 /// It handles a dynamic modification
 virtual void dynamic_modification( BlockModAD * mod );

 /// It adds a single new dynamic constraint
 virtual void add_dynamic_constraint( FRowConstraint * p_const );

 /// It adds a single new dynamic variable
 virtual void add_dynamic_variable( ColVariable * p_var );

 /// It adds a single new dynamic bound
 virtual void add_dynamic_bound( OneVarConstraint * p_bound );

 /// It removes a single dynamic constraint
 virtual void remove_dynamic_constraint( const FRowConstraint * p_const );

 /// It removes a single dynamic variable
 virtual void remove_dynamic_variable( const ColVariable * p_var );

 /// It removes a single dynamic bound
 virtual void remove_dynamic_bound( const OneVarConstraint * p_bound );
 /// @}

/*--------------------------------------------------------------------------*/
/*--------------------- PRIVATE FIELDS OF THE CLASS ------------------------*/
/*--------------------------------------------------------------------------*/

 private:

/*--------------------------------------------------------------------------*/
/*------------- AUXILIARY METHODS FOR POPULATING THE PROBLEM  --------------*/
/*--------------------------------------------------------------------------*/

/**
 * @name Private auxiliary methods
 *
 * These methods are used in load_problem() to read data from the Block.
 * All of them, except scan_objective(), take integer counters as input
 * parameters. That's because they are meant to be used by
 * un_any_const_static() and un_any_const_dynamic() template functions
 * on boost::any containers, and the counters keep track of the elements
 * inside the containers.
 *
 * @{
 */

 /**
  * Accumulating procedure that helps getting the total number of rows.
  * If the constraint is linear, n_rows is incremented.
  *
  * @param constraint FRowConstraint to be counted
  * @param n_rows accumulator for rows
  */
 void count_constraints( FRowConstraint & constraint, int & n_rows );

 /**
  * Accumulating procedure that helps getting the total number of columns.
  *
  * @param variable ColVariable to be counted
  * @param n_cols accumulator for columns
  */
 void count_variables( ColVariable & variable, int & n_cols );

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
 void count_nzelements( ColVariable & variable, int & nz_elements, int & cnt );

 /**
  * It scans a static ColVariable or a group of static ColVariables and fills
  * the vectors of the LP accordingly.
  *
  * @param var a reference to a [vector of] static ColVariable[s]
  * @param first an counter that should be 0 when var is the first
  *              element of a vector of ColVariables.
  * @param i a counter for variables/columns shared w/ scan_dynamic_variable()
  */
 void scan_static_variable( ColVariable & var, int & first, int & i );

 /**
  * It scans a dynamic ColVariable and fills the vectors of the LP accordingly.
  *
  * @param var a reference to a dynamic ColVariable
  * @param i a counter for variables/columns shared with scan_static_variable()
  */
 void scan_dynamic_variable( ColVariable & lvar, int & i );

 /**
  * It scans a static FRowConstraint or a group of static FRowConstraints and
  * fills the vectors of the LP accordingly.
  * @param p_const a reference to a [vector of] static FRowConstraint[s]
  * @param first an counter that should be 0 when lconst is the first
  *              element of a vector of FRowConstraints.
  * @param i a counter for constraints/rows shared w/ scan_dynamic_constraint()
  */
 void scan_static_constraint( FRowConstraint & p_const, int & first, int & i );

 /**
  * It scans a dyn FRowConstraint and fills the vectors of the LP accordingly.
  *
  * @param p_const a reference to a dynamic FRowConstraint
  * @param i a counter for constraints/rows shared w/ scan_static_constraint()
  */
 void scan_dynamic_constraint( FRowConstraint & p_const, int & i );

 /**
  * It scans a FRealObjective and fills the vectors of the LP accordingly.
  * @param obj a FRealObjective
  */
 void scan_objective( const FRealObjective * obj );
 /// @}

 /// Returns a loggable representation of a vector
 template<typename T> std::string log_vector(std::vector<T> v);

 SMSpp_insert_in_factory_h;
};   // end( class MILPSolver )

} // end( namespace SMSpp_di_unipi_it )

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

#endif

