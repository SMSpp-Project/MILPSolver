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
 * \copyright &copy; Antonio Frangioni, Niccolò Iardella
 */
/*--------------------------------------------------------------------------*/
/*----------------------------- DEFINITIONS --------------------------------*/
/*--------------------------------------------------------------------------*/

#ifndef __MILPSolver
#define __MILPSolver

#ifdef MILPSOLVER_DEBUG
#define DEBUG_LOG( stuff ) std::cout << "[MILPSolver DEBUG] " << stuff
#else
#define DEBUG_LOG( stuff )
#endif

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

 /// Types of integer parameters
 enum int_par_type_MILP {
  /// Use custom names for rows/columns
  intUseCustomNames = intLastParCDAS,
  /// First allowed new int parameter for derived classes
  intLastAlgParMILP
 };

 /// Types of double parameters
 enum dbl_par_type_MILP {
  /// First allowed new double parameter for derived classes
  dblLastAlgParMILP = dblLastParCDAS
 };

 /// Types of string parameters
 enum str_par_type_MILP {
  /// Problem name
  strProblemName = strLastParCDAS,
  /// Output file
  strOutputFile,
  /// First allowed new string parameter for derived classes
  strLastAlgParMILP
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

 /**
  * It clears the LP vectors.
  *
  * The input parameter is a bitwise value that allows to specify which
  * vectors should be cleared. From the LSB to the MSB:
  *
  * - 1 clears the costraint matrix, xctype the column/row names
  * - 2 clears the OF related vectors (objective and q_objective)
  * - 4 clears rhs, rngval and sense
  * - 8 clears lb and ub.
  *
  * To clear everything, use what = 15.
  */
 virtual void clear_problem( unsigned int what );

 /// It loads the problem from the Block into the LP vectors
 virtual void load_problem();

 /// @}

/** @name Getters for the vectors of the LP problem.
 *
 * The following methods return the data that define the LP problem as
 * specified in the CPLEX CPXcopylp(), CPXcopyctype(),  CPXcopyqpsep() and
 * CPXcopylpwnames() documentation.
 *
 * @{
 */

 /// Returns the number of variables/columns
 [[nodiscard]] int get_numcols() const;

 /// Returns the number of constraints/rows
 [[nodiscard]] int get_numrows() const;

 /// Returns the number of non-zero elements
 [[nodiscard]] int get_nzelements() const;

 /// Returns the sense of the objective function, see CPXchgobjsen()
 [[nodiscard]] int get_objsense() const;

 /// Returns the linear cofficients of the objective function
 [[nodiscard]] const std::vector< double > & get_objective() const;

 /// Returns the quadratic cofficients of the objective function
 [[nodiscard]] const std::vector< double > & get_q_objective() const;

 /// Returns the RHS values of the constraints
 [[nodiscard]] const std::vector< double > & get_rhs() const;

 /// Returns the range values of the ranged constraints
 [[nodiscard]] const std::vector< double > & get_rngval() const;

 /// Returns the sense of the constraints, see  CPXchgsense()
 [[nodiscard]] const std::vector< char > & get_sense() const;

 /// Returns matbeg, one of the arrays that define the constraint matrix
 [[nodiscard]] const std::vector< int > & get_matbeg() const;

 /// Returns matcnt, one of the arrays that define the constraint matrix
 [[nodiscard]] const std::vector< int > & get_matcnt() const;

 /// Returns matind, one of the arrays that define the constraint matrix
 [[nodiscard]] const std::vector< int > & get_matind() const;

 /// Returns matval, one of the arrays that define the constraint matrix
 [[nodiscard]] const std::vector< double > & get_matval() const;

 /// Returns the lower bounds on the variables
 [[nodiscard]] const std::vector< double > & get_var_lb() const;

 /// Returns the upper bounds on the variables
 [[nodiscard]] const std::vector< double > & get_var_ub() const;

 /// Returns the types of the variables, see CPXcopyctype()
 [[nodiscard]] const std::vector< char > & get_xctype() const;

 /// Returns the names of the constraints/rows
 [[nodiscard]] const std::vector< char * > & get_rowname() const;

 /// Returns the names of the variables/columns
 [[nodiscard]] const std::vector< char * > & get_colname() const;

 /// @}

 /** @name Methods that use the dictionaries
  *
  * The following methods use the dictionaries to get the indices of the
  * Variables/Constraints from the pointers and viceversa.
  *
  * We provide separate methods for looking into static, dynamic or both parts
  * of the problem, so we can reduce searching time when possible.
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
 int index_of_variable( const ColVariable * p_var );

 /**
  * It returns the matrix column index of a given static variable.
  *
  * @param p_var a pointer to a ColVariable
  * @return the corresponding matrix column index
  * @throws std::invalid_argument if no indices are associated to the variable
  */
 int index_of_static_variable( const ColVariable * p_var );

 /**
  * It returns the matrix column index of a given dynamic variable.
  *
  * @param p_var a pointer to a ColVariable
  * @return the corresponding matrix column index
  * @throws std::invalid_argument if no indices are associated to the variable
  */
 int index_of_dynamic_variable( const ColVariable * p_var );

 /**
  * It returns the matrix row index of the given constraint.
  *
  * @param p_const a pointer to a FRowConstraint
  * @return the corresponding matrix row index
  */
 int index_of_constraint( const FRowConstraint * p_const );

 /**
  * It returns the matrix row index of the given static constraint.
  *
  * @param p_const a pointer to a FRowConstraint
  * @return the corresponding matrix row index
  * @throws std::invalid_argument if no indices are associated to the constraint
  */
 int index_of_static_constraint( const FRowConstraint * p_const );

 /**
  * It returns the matrix row index of the given dynamic constraint.
  *
  * @param p_const a pointer to a FRowConstraint
  * @return the corresponding matrix row index
  * @throws std::invalid_argument if no indices are associated to the constraint
  */
 int index_of_dynamic_constraint( const FRowConstraint * p_const );

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

 /// Writes the LP on the specified file
 virtual void write_lp( const std::string & filename ) {}

 /// Returns the number of nodes used to solve a MIP
 [[nodiscard]] virtual int get_nodes() const { return 0; }

 /// Returns the number of integer variables
 [[nodiscard]] int get_num_integer_vars() const;

#ifdef MILPSOLVER_DEBUG

 /// Check the dictionaries for inconsistencies
 virtual void check_status();

#endif

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

 /// It does nothing as there is nothing to do.
 void get_var_solution( Configuration * solc ) override {}

 /// It does nothing as there is nothing to do.
 void get_dual_solution( Configuration * solc ) override {}
 /// @}

 /*--------------------------------------------------------------------------*/
/*------------------- METHODS FOR HANDLING THE PARAMETERS ------------------*/
/*--------------------------------------------------------------------------*/
/**
 * @name Methods for handling parameters
 * @{
 */

 /// Sets an integer parameter with the given value
 void set_par( idx_type par, int value ) override;

 /// Sets a double parameter with the given value
 void set_par( idx_type par, double value ) override;

 /// Sets a string parameter with the given value
 void set_par( idx_type par, const std::string & value ) override;

 /// Gets the number of integer parameters
 [[nodiscard]] idx_type get_num_int_par() const override;

 /// Gets the number of double parameters
 [[nodiscard]] idx_type get_num_dbl_par() const override;

 /// Gets the number of string parameters
 [[nodiscard]] idx_type get_num_str_par() const override;

 /// Gets the default value of the specified integer parameter
 [[nodiscard]] int get_dflt_int_par( idx_type par ) const override;

 /// Gets the default value of the specified double parameter
 [[nodiscard]] double get_dflt_dbl_par( idx_type par ) const override;

 /// Gets the default value of the specified string parameter
 [[nodiscard]] const std::string &
 get_dflt_str_par( idx_type par ) const override;

 /// Gets the value of the specified integer parameter
 [[nodiscard]] int get_int_par( idx_type par ) const override;

 /// Gets the value of the specified double parameter
 [[nodiscard]] double get_dbl_par( idx_type par ) const override;

 /// Gets the value of the specified string parameter
 [[nodiscard]] const std::string & get_str_par( idx_type par ) const override;

 /// Returns the index of the int parameter with the specified name
 [[nodiscard]] idx_type
 int_par_str2idx( const std::string & name ) const override;

 /// Returns the name of the int parameter with the specified index
 [[nodiscard]] const std::string &
 int_par_idx2str( idx_type idx ) const override;

 /// Returns the index of the double parameter with the specified name
 [[nodiscard]] idx_type
 dbl_par_str2idx( const std::string & name ) const override;

 /// Returns the name of the double parameter with the specified index
 [[nodiscard]] const std::string &
 dbl_par_idx2str( idx_type idx ) const override;

 /// Returns the index of the string parameter with the specified name
 [[nodiscard]] idx_type
 str_par_str2idx( const std::string & name ) const override;

 /// Returns the name of the string parameter with the specified index
 [[nodiscard]] const std::string &
 str_par_idx2str( idx_type idx ) const override;
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
 typedef std::tuple< ColVariable *, int, int > var_int_int;
 typedef std::tuple< FRowConstraint *, int, int > con_int_int;

 /** @name Variable and Constraint dictionaries
  *
  * The following vectors are used in order to keep track between the
  * Variables and Constraints of the Block and the constraint matrix.
  *
  *  - svar_to_idx, scon_to_idx : vectors of tuples that store 1) the address
  *    of the first element of each group of static variables and
  *    constraints, respectively, 2) the corresponding index in constraint
  *    matrix (column or row), and 3) the number of elements in the group.
  *    Vectors are kept sorted in ascending order by address.
  *
  *  - dvar_to_idx, dcon_to_idx : vectors of pairs that store the addresses
  *    of all the dynamic variables and constraints respectively and the
  *    corresponding index in constraint matrix. and are stored in an ascending
  *    Vectors are kept sorted in ascending order by address.
  *
  *  - idx_to_svar, idx_to_scon : vectors of pairs that store the indices of
  *    columns and rows, respectively, of the constraint matrix and the
  *    address of the corresponding (group of) static variables and constraints.
  *    Vectors are kept sorted in ascending order by index.
  *
  *  - idx_to_dvar, idx_to_dcon : vectors of pairs that store the indices of
  *    columns and rows, respectively, of the constraint matrix and the
  *    address of the corresponding dynamic variables and constraints.
  *    Vectors are kept sorted in ascending order by index.
  *
  * Using these vectors of pair we can at any time locate the index
  * of each constraint and variable within the constraint matrix, and viceversa.
  * Note that the Block stores static variables and constraints grouped, so
  * the vectors that store addresses of static elements contain the address
  * of the first element of each group.
  * The methods that use these vectors need to be aware of that.
  *
  * @{
  */
 std::vector< var_int_int > svar_to_idx; ///< From static variable to index
 std::vector< int_var > idx_to_svar;     ///< From index to static variable

 std::vector< con_int_int > scon_to_idx; ///< From static constraint to index
 std::vector< int_const > idx_to_scon;   ///< From index to static constraint

 std::vector< var_int > dvar_to_idx;     ///< From dynamic variable to index
 std::vector< int_var > idx_to_dvar;     ///< From index to dynamic variable

 std::vector< const_int > dcon_to_idx;   ///< From dynamic constraint to index
 std::vector< int_const > idx_to_dcon;   ///< From index to dynamic constraint
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

 std::string prob_name;    ///< Problem name
 std::string output_file;  ///< Output file

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


 /// Use Variable/Constraint custom names
 bool use_custom_names = true;

 /**
  * An array of length at least numcols containing pointers to character
  * strings containing the names of the variables.
  */
 std::vector< char * > colname;

 /**
  * An array of length at least numrows containing pointers to character
  * strings containing the names of the constraints.
  */
 std::vector< char * > rowname;

 int sol_status = kUnEval; ///< Solution status (OK, Infeasible, Unbounded, ...)
 int int_vars{};           ///< Number of integer variables
 int static_vars{};        ///< Number of static variables
 int static_cons{};        ///< Number of static constraints
 /// @}

 /** @name Get variable bounds for the problem
  *
  * The following two methods retrieve the upper and lower bound for the
  * given variable considering both the Variable bounds and all the active
  * OneVarConstraints active for that Variable.
  *
  * @{
  */

 /// Gets the LB fot the given variable in the problem
 virtual double get_problem_lb( const ColVariable & var );

 /// Gets the UB fot the given variable in the problem
 virtual double get_problem_ub( const ColVariable & var );
 /// @}

 /// Gets the active constraints for the specified variable
 // TODO: This should be temporary
 static std::vector< FRowConstraint * >
 get_active_constraints( const ColVariable & var );

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
  * It scans a ColVariable and fills the vectors of the LP accordingly.
  *
  * @param var a reference to a ColVariable
  * @param n   an counter that should be 0 when var is the first
  *              element of a vector of static ColVariables,
  *              and -1 when it's a dynamic ColVariable
  * @param col a counter for variables/columns
  */
 void scan_variable( ColVariable & var, int & n, int & col );

 /**
  * It scans a FRowConstraint and fills the vectors of the LP accordingly.
  * @param con a reference to a FRowConstraint
  * @param n   an counter that should be 0 when lconst is the first
  *              element of a vector of static FRowConstraints,
  *              and -1 when it's a dynamic FRowConstraint
  * @param row a counter for constraints/rows
  */
 void scan_constraint( FRowConstraint & con, int & n, int & row );


 /**
  * It scans a FRealObjective and fills the vectors of the LP accordingly.
  * @param obj a FRealObjective
  */
 void scan_objective( const FRealObjective * obj );
 /// @}

 SMSpp_insert_in_factory_h;
};

}

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

#endif

