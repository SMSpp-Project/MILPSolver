/*--------------------------------------------------------------------------*/
/*---------------------------- File MILPSolver.h ---------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the MILPSolver class, which implements parts of the
 * Solver concept for MILP solvers. The class does not directly provide
 * solving capabilities, but rather a first layer for analysing a "MILP
 * Block" and constructing/maintaining a standard "sparse matrix as a vector
 * of doubles + two vectors of int, plus vectors for costs, bounds and
 * lhs/rhs of constraints" representation of the problem. This is thought to
 * be used as an input for derived classes that use it to provide actual
 * solving capabilities.
 *
 * However, MILPSolver is also engineered to be able to be used in isolation
 * (that is, it is not an abstract class) for just keeping up the standard
 * "sparse matrix as a vector of doubles + two vectors of int, plus vectors
 * for costs, bounds and lhs/rhs of constraints" representation of the
 * problem. 
 *
 * \author Antonio Frangioni \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \author Niccolo' Iardella \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \copyright &copy; by Antonio Frangioni, Niccolo' Iardella
 */
/*--------------------------------------------------------------------------*/
/*----------------------------- DEFINITIONS --------------------------------*/
/*--------------------------------------------------------------------------*/

#ifndef __MILPSolver
 #define __MILPSolver /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

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
namespace SMSpp_di_unipi_it
{
/*--------------------------------------------------------------------------*/
/*------------------------- CLASS MILPSolver -------------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------- GENERAL NOTES --------------------------------*/
/*--------------------------------------------------------------------------*/
/// base class for solving MILP problems.
/** The MILPSolver class derives from Solver and extends the interface of
 * the base class to be able to efficiently handle MILP problems. This class
 * alone does not solve problems, but it serves as base class to other
 * solvers, to be implemented in derived classes. Nonetheless it can be used
 * by itself wherever a description of the MILP problem in matricial form is
 * needed.
 *
 * In fact, since LP is a special case of MILP and therefore MILP solvers
 * necessarily need to have LP solving capabilities,  MILPSolver derives
 * from CDASolver and implements the interface for obtaining dual solutions
 * as well, although this can of course only used when no integer variables
 * are present.
 *
 * The MILPSolver can be registered to any kind of Block (assuming that it
 * contains a MILP formulation) and it generates a collection of vectors
 * that describes the MILP problem in the usual form "sparse matrix as a
 * vector of doubles + two vectors of int, plus vectors for costs, bounds
 * and lhs/rhs of constraints". This makes it easy to construct derived
 * classes that interface with standard solvers.
 *
 * The main thing that this class has to take care is the correspondence
 * between the Constraints and Variables of the Block and the constraint
 * matrix. The correspondence is built via set_Block(), that conducts a
 * Breadth First Search, scanning the Block and all its children, if any,
 * populating the vectors needed to define the MILP problem.
 *
 * Following the Block structure we build the constraint matrix in two
 * steps:
 *
 * 1. We build the static part, i.e., the part of the rows (Constraints)
 *    and columns (Variables) that is not going to be deleted for all the
 *    lifecycle of the problem.
 *
 * 2. Then, the dynamic part, i.e., the part of the rows (Constraints) and
 *    columns (Variables) that can potentially change or stop existing.
 *
 * After set_Block() has been called, the constraint matrix is assumed to
 * have all the information needed from the Block to solve the problem.
 * The information can be retrieved by a library of getters. The methods
 * compute(), get_var_solution() and new_var_solution(), derived from the
 * Solver class, do nothing and should be implemented by derived classes.
 *
 * The class defines also an interface that the derived classes should
 * implement to support Modification. The method process_modifications() is
 * already implemented and it is the one that dispatches the Modifications
 * to the other methods accordingly. Provided that the derived classes
 * properly implement the virtual methods for handling the individual kinds
 * of Modification a(n admittedly, rather rough on that it processes each
 * Modification individually without any attempt at optimizing the sequence
 * and reducing the operations on the MILP problem) complete solution for
 * handling the Modification is obtained.
 *
 * Note that it is generally assumed that the "sparse matrix as a vector of
 * doubles + two vectors of int, plus vectors for costs, bounds and lhs/rhs
 * of constraints" representation of the MILP problem is only useful when
 * initialising the actual MILP solver in derived classes, and can be
 * deleted right after that: any change in the formulation should directly
 * happen in the solver internal representation. However, MILPSolver still
 * has a role at least in keeping updated the dictionaries that map the
 * [Col]Variable and [FRow]Constraint of the Block into the columns and
 * rows of the coefficient matrix and back. Furthermore, MILPSolver is also
 * engineered to be able to be used in isolation (that is, it is not an
 * abstract class) for just keeping up the standard matrix-based
 * representation of the MILP. For this reason, it is possible to freely
 * choose which parts of the matrix-based representation will be discarded
 * after that the problem has been loaded in the underlying solver (if any)
 * and which ones will rather be kept, and properly updated. */

class MILPSolver : public CDASolver
{
/*--------------------------------------------------------------------------*/
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

 public:

/*--------------------------------------------------------------------------*/
/*---------------------------- PUBLIC TYPES --------------------------------*/
/*--------------------------------------------------------------------------*/

 /// enum for int parameters
 enum int_par_type_MILP {
  intUseCustomNames = intLastParCDAS , ///< use custom names for rows/columns
  /// Relax [M]ILP by removing integrality constraints for integer variables
  intRelaxIntVars ,
  intLastAlgParMILP  ///< 1st allowed new int parameter for derived classes
  };

 /// enum for double parameters
 enum dbl_par_type_MILP {
  /// First allowed new double parameter for derived classes
  dblLastAlgParMILP = dblLastParCDAS
  };

 /// enum for string parameters
 enum str_par_type_MILP {
  strProblemName = strLastParCDAS ,  ///< problem name
  strOutputFile ,                    ///< output filename
  strLastAlgParMILP  ///< 1st allowed new string parameter for derived classes
  };

 /// enum for vector-of-int parameters
 enum vint_par_type_MILP {
  ///< first allowed new vector-of-int parameter for derived classes
  vintLastAlgParMILP = vintLastParCDAS
  };

 /// enum for vector-of-double parameters
 enum vdbl_par_type_MILP {
  /// first allowed new vector-of-double parameter for derived classes
  vdblLastAlgParMILP = vdblLastParCDAS
  };

 /// enum for vector-of-string parameters
 enum vstr_par_type_MILP {
  /// first allowed new vector-of-double parameter for derived classes
  vstrLastAlgParMILP = vstrLastParCDAS
  };

 using Index = Block::Index;  // "import" Index from Block

/*--------------------------------------------------------------------------*/
/*--------------------- CONSTRUCTOR AND DESTRUCTOR -------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Constructor and destructor
 * @{ */

 MILPSolver( void ) : CDASolver() {}

/*--------------------------------------------------------------------------*/

 ~MILPSolver() override {
  for( auto & i : colname )
   delete[] i;
  for( auto & i : rowname )
   delete[] i;
  }

/** @} ---------------------------------------------------------------------*/
/*--------------------- PUBLIC METHODS OF THE CLASS ------------------------*/
/*--------------------------------------------------------------------------*/
 /** @name Clear and load the problem
  *
  * The following two methods include the main logic of the class. They are
  * called in the set_Block() method to build the problem when the MILPSolver
  * is [re]registered to a Block, but also when a NBModification is processed.
  * Note: A derived class can override these methods but must call the base
  * versions if it wants to use the matrix-based representation (which is
  * the only rationale for deriving from MILPSolver in the first place).
  *
  * Note that it is generally assumed that the "sparse matrix as a vector of
  * doubles + two vectors of int, plus vectors for costs, bounds and lhs/rhs
  * of constraints" representation of the MILP problem is only useful when
  * initialising the actual MILP solver in derived classes, and can be
  * deleted right after that: any change in the formulation should directly
  * happen in the solver internal representation. However, MILPSolver still
  * has a role at least in keeping updated the dictionaries that map the
  * [Col]Variable and [FRow]Constraint of the Block into the columns and
  * rows of the coefficient matrix and back. Furthermore, MILPSolver is also
  * engineered to be able to be used in isolation (that is, it is not an
  * abstract class) for just keeping up the standard matrix-based
  * representation of the MILP. For this reason, it is possible to freely
  * choose which parts of the matrix-based representation will be discarded
  * after that the problem has been loaded in the underlying solver (if any)
  * and which ones will rather be kept, and properly updated. This is
  * controlled by the parameter of the clear_problem() method.
  *
  * IMPORTANT NOTE: retaining the information related to the Objective, the
  *                 LHS/RHS of the Constraint and the LB/UB of the Variable
  *                 and changing it when Modification are received is
  *                 properly implemented, but doing the same for the matrix
  *                 is currently *not*. Hence, exceptions will be thrown if
  *                 the matrix is kept and a Modification changing it occurs.
  * @{ */

 /// clears the description of the MILP
 /** Clears the matrix-based representation of the MILP.
  *
  * The input parameter is a bitwise value that allows to specify which
  * vectors should be cleared. From the LSB to the MSB:
  *
  * - 1 clears the constraint matrix, xctype the column/row names
  * - 2 clears the OF related vectors (objective and q_objective)
  * - 4 clears rhs, rngval and sense
  * - 8 clears lb and ub.
  *
  * To clear everything, use what = 15.
  *
  * Note: this method is provided so the user can clear the stuff that
  * she is sure IT WILL NOT BE CHANGED. If some vectors are cleared and
  * a method tries to change them, it will throw exception! */

 virtual void clear_problem( unsigned int what );

/*--------------------------------------------------------------------------*/
 /// loads the problem from the Block into the MILP vectors

 virtual void load_problem( void );

/** @} ---------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Getters for the vectors of the MILP problem.
 *
 * The following methods return the data that define the MILP problem.
 * @{ */

 /// returns the number of variables/columns
 [[nodiscard]] int get_numcols( void ) const { return( numcols ); }

 /// returns the number of constraints/rows
 [[nodiscard]] int get_numrows( void ) const { return( numrows ); }

 /// returns the number of non-zero elements
 [[nodiscard]] int get_nzelements( void ) const { return( matval.size() ); }

 /// returns the sense of the objective function, see CPXchgobjsen()
 [[nodiscard]] int get_objsense( void ) const { return( objsense ); }

 /// returns the linear coefficients of the objective function
 [[nodiscard]] const std::vector< double > & get_objective( void ) const {
  return( objective );
  }

 /// returns the quadratic coefficients of the objective function
 [[nodiscard]] const std::vector< double > & get_q_objective( void ) const {
  return( q_objective );
  }

 /// returns the RHS values of the constraints
 [[nodiscard]] const std::vector< double > & get_rhs( void ) const {
  return( rhs );
  }

 /// returns the range values of the ranged constraints
 [[nodiscard]] const std::vector< double > & get_rngval( void ) const {
  return( rngval );
  }

 /// returns the sense of the constraints, see  CPXchgsense()
 [[nodiscard]] const std::vector< char > & get_sense( void ) const {
  return( sense );
  }

 /// returns matbeg, one of the arrays that define the constraint matrix
 [[nodiscard]] const std::vector< int > & get_matbeg( void ) const {
  return( matbeg );
  }

 /// returns matcnt, one of the arrays that define the constraint matrix
 [[nodiscard]] const std::vector< int > & get_matcnt( void ) const {
  return( matcnt );
  }

 /// returns matind, one of the arrays that define the constraint matrix
 [[nodiscard]] const std::vector< int > & get_matind( void ) const {
  return( matind );
  }

 /// returns matval, one of the arrays that define the constraint matrix
 [[nodiscard]] const std::vector< double > & get_matval( void ) const {
  return( matval );
  }

 /// returns the lower bounds on the variables
 [[nodiscard]] const std::vector< double > & get_var_lb( void ) const {
  return( lb );
  }

 /// returns the upper bounds on the variables
 [[nodiscard]] const std::vector< double > & get_var_ub( void ) const {
  return( ub );
  }

 /// returns the types of the variables
 [[nodiscard]] const std::vector< char > & get_xctype( void ) const {
  return( xctype );
  }

 /// returns the names of the constraints/rows
 [[nodiscard]] const std::vector< char * > & get_rowname( void ) const {
  return( rowname );
  }

 /// returns the names of the variables/columns
 [[nodiscard]] const std::vector< char * > & get_colname( void ) const {
  return( colname );
  }

/** @} ---------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
 /** @name Methods that use the dictionaries
  *
  * The following methods use the dictionaries to get the indices of the
  * Variables/Constraints from the pointers and viceversa.
  *
  * We provide separate methods for looking into static, dynamic or both parts
  * of the problem, so we can reduce searching time when possible.
  * @{ */

 /// returns the matrix column index of a given variable
 /** Returns the matrix column index of a given variable.
  *
  * @param var a pointer to a ColVariable
  * @return the corresponding matrix column index, or Inf< int >() if the
  *         ColVariable was not found */

 int index_of_variable( const ColVariable * var ) const;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// returns the matrix column index of a given static variable
 /** Returns the matrix column index of a given static variable.
  *
  * @param var a pointer to a ColVariable
  * @return the corresponding matrix column index, or Inf< int >() if the
  *         ColVariable was not found among the static ones */

 int index_of_static_variable( const ColVariable * var ) const;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /** Returns the matrix column index of a given dynamic variable.
  *
  * @param var a pointer to a ColVariable
  * @return the corresponding matrix column index, or Inf< int >() if the
  *         ColVariable was not found among the dynamic ones */

 int index_of_dynamic_variable( const ColVariable * var ) const;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// returns the matrix row index of the given constraint
 /** Returns the matrix row index of the given constraint.
  *
  * @param con a pointer to a FRowConstraint
  * @return the corresponding matrix row index, or Inf< int >() if the
  *         FRowConstraint was not found */

 int index_of_constraint( const FRowConstraint * con ) const;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// returns the matrix row index of the given static constraint
 /** Returns the matrix row index of the given static constraint.
  *
  * @param con a pointer to a FRowConstraint
  * @return the corresponding matrix row index, or Inf< int >() if the
  *         FRowConstraint was not found among the static ones */

 int index_of_static_constraint( const FRowConstraint * con ) const;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// returns the matrix row index of the given dynamic constraint
 /** Returns the matrix row index of the given dynamic constraint.
  *
  * @param con a pointer to a FRowConstraint
  * @return the corresponding matrix row index, or Inf< int >() if the
  *         FRowConstraint was not found among the dynamic ones */

 int index_of_dynamic_constraint( const FRowConstraint * con ) const;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// returns the variable corresponding to the given column index
 /** Returns the variable corresponding to a variable matrix column index.
  *
  * @param i a constraint matrix column index
  * @return a pointer to the corresponding ColVariable, or nullptr if \p i
  *         is an invalid variable index */

 const ColVariable * variable_with_index( int i ) const;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// returns the static variable corresponding to the given column index
 /** Returns the static variable corresponding to the given variable
  * matrix column index.
  *
  * @param i a constraint matrix column index
  * @return a pointer to the corresponding ColVariable, or nullptr if \p i
  *         is an invalid static variable index */

 const ColVariable * static_variable_with_index( int i ) const;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// returns the dynamic variable corresponding to the given column index
 /** Returns the dynamic variable corresponding to the given variable
  * matrix column index.
  *
  * @param i a constraint matrix column index
  * @return a pointer to the corresponding ColVariable, or nullptr if \p i
  *         is an invalid dynamic variable index */

 const ColVariable * dynamic_variable_with_index( int i ) const;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// returns the constraint corresponding to the given matrix row index
 /** Returns the constraint corresponding to the given matrix row index.
  *
  * @param i a constraint matrix row index
  * @return a pointer to the corresponding FRowConstraint, or nullptr if \p i
  *         is an invalid constraint index */

 const FRowConstraint * constraint_with_index( int i ) const;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// returns the static constraint corresponding to the given row index
 /** Returns the static constraint corresponding to the given constraint
  * matrix row index.
  *
  * @param i a constraint matrix row index
  * @return a pointer to the corresponding FRowConstraint, or nullptr if \p i
  *         is an invalid static constraint index */

 const FRowConstraint * static_constraint_with_index( int i ) const;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// returns the dynamic constraint corresponding to the given row index
 /** Returns the dynamic constraint corresponding to the given constraint
  * matrix row index.
  *
  * @param i a constraint matrix row index
  * @return a pointer to the corresponding FRowConstraint, or nullptr if \p i
  *         is an invalid dynamic constraint index */

 const FRowConstraint * dynamic_constraint_with_index( int i ) const;

/** @} ---------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

 /// Writes the LP on the specified file
 virtual void write_lp( const std::string & filename ) {}

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// Returns the number of nodes used to solve a MIP
 [[nodiscard]] virtual int get_nodes( void ) const { return( 0 ); }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// Returns the number of integer variables
 [[nodiscard]] int get_num_integer_vars( void ) const { return( int_vars ); }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 #ifdef MILPSOLVER_DEBUG
  /// Check the dictionaries for inconsistencies
  virtual void check_status( void );
 #endif

/*--------------------------------------------------------------------------*/
/*--------------------- DERIVED METHODS OF BASE CLASS ----------------------*/
/*--------------------------------------------------------------------------*/
/** @name Public Methods derived from base classes
 *  @{ */

 /// sets the Block that the Solver has to solve and build the MILP vectors
 void set_Block( Block * block ) override;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// does nothing as there is nothing to do
 int compute( bool changedvars = true ) override;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// does nothing as there is nothing to do
 void get_var_solution( Configuration * solc ) override {}

 /// does nothing as there is nothing to do
 void get_dual_solution( Configuration * solc ) override {}

/** @} --------------------------------------------------------------------*/
/*------------------- METHODS FOR HANDLING THE PARAMETERS ------------------*/
/*--------------------------------------------------------------------------*/
/** @name Methods for handling parameters
 * @{ */

 /// sets an integer parameter with the given value
 void set_par( idx_type par , int value ) override;

 // sets a double parameter with the given value
 // void set_par( idx_type par , double value ) override;

 /// sets a string parameter with the given value
 void set_par( idx_type par , std::string && value ) override;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// gets the number of integer parameters
 [[nodiscard]] idx_type get_num_int_par( void ) const override;

 /// gets the number of double parameters
 [[nodiscard]] idx_type get_num_dbl_par( void ) const override;

 /// gets the number of string parameters
 [[nodiscard]] idx_type get_num_str_par( void ) const override;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// gets the default value of the specified integer parameter
 [[nodiscard]] int get_dflt_int_par( idx_type par ) const override;

 // gets the default value of the specified double parameter
 // [[nodiscard]] double get_dflt_dbl_par( idx_type par ) const override;

 /// returns the default value of the specified string parameter
 [[nodiscard]] const std::string & get_dflt_str_par( idx_type par )
  const override;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// returns the value of the specified integer parameter
 [[nodiscard]] int get_int_par( idx_type par ) const override;

 // returns the value of the specified double parameter
 // [[nodiscard]] double get_dbl_par( idx_type par ) const override;

 /// returns the value of the specified string parameter
 [[nodiscard]] const std::string & get_str_par( idx_type par ) const override;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// returns the index of the int parameter with the specified name
 [[nodiscard]] idx_type int_par_str2idx( const std::string & name )
  const override;

 /// returns the name of the int parameter with the specified index
 [[nodiscard]] const std::string & int_par_idx2str( idx_type idx )
  const override;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 // returns the index of the double parameter with the specified name
 // [[nodiscard]] idx_type
 // dbl_par_str2idx( const std::string & name ) const override;

 // returns the name of the double parameter with the specified index
 // [[nodiscard]] const std::string &
 // dbl_par_idx2str( idx_type idx ) const override;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// returns the index of the string parameter with the specified name
 [[nodiscard]] idx_type str_par_str2idx( const std::string & name )
  const override;

 /// Returns the name of the string parameter with the specified index
 [[nodiscard]] const std::string & str_par_idx2str( idx_type idx )
  const override;

/** @} ---------------------------------------------------------------------*/
/*-------------------- PROTECTED FIELDS OF THE CLASS -----------------------*/
/*--------------------------------------------------------------------------*/

 protected:

/*--------------------------------------------------------------------------*/
/*---------------- VARIABLE AND CONSTRAINT TRACKING VECTORS ----------------*/
/*--------------------------------------------------------------------------*/

 using var_int = std::pair< const ColVariable * , int >;
 using int_var = std::pair< int , const ColVariable * >;
 using con_int = std::pair< const FRowConstraint * , int >;
 using int_con = std::pair< int , const FRowConstraint * >;
 using var_int_int = std::tuple< const ColVariable * , int , int >;
 using con_int_int = std::tuple< const FRowConstraint * , int , int >;

 /** @name Variable and Constraint dictionaries
  *
  * The following vectors are used in order to keep track between the
  * Variables and Constraints of the Block and the constraint matrix.
  *
  *  - svar_to_idx, scon_to_idx : vectors of tuples that store 1) the
  *    address of the first element of each group of static variables and
  *    constraints, respectively, 2) the corresponding index in constraint
  *    matrix (column or row), and 3) the number of elements in the group.
  *    The vectors are kept sorted in ascending order by address.
  *
  *  - dvar_to_idx, dcon_to_idx : vectors of pairs that store the addresses
  *    of all the dynamic variables and constraints respectively and the
  *    corresponding index in constraint matrix. and are stored in an
  *    ascending. The vectors are kept sorted in ascending order by address.
  *
  *  - idx_to_svar, idx_to_scon : vectors of pairs that store the indices
  *    of columns and rows, respectively, of the constraint matrix and the
  *    address of the corresponding (group of) static variables and
  *    constraints. The vectors are kept sorted in ascending order by index.
  *
  *  - idx_to_dvar, idx_to_dcon : vectors that store the addresses of dynamic
  *    variables and constraints, respectively. The element at idx_to_d*[ i ]
  *    has index ( i + static_*s ), that is, it is the column/row
  *    ( i + static_*s ) of the constraint matrix.
  *
  * Using these vectors of pair we can at any time locate the index of each
  * constraint and variable within the constraint matrix, and viceversa.
  * Note that the Block stores static variables and constraints grouped, so
  * the vectors that store addresses of static elements contain the address
  * of the first element of each group. The methods that use these vectors
  * need to be aware of that.
  * @{ */

 std::vector< var_int_int > svar_to_idx; ///< from static variable to index
 std::vector< int_var > idx_to_svar;     ///< from index to static variable

 std::vector< con_int_int > scon_to_idx; ///< from static constraint to index
 std::vector< int_con> idx_to_scon;      ///< from index to static constraint

 std::vector< var_int > dvar_to_idx;     ///< from dynamic variable to index
 std::vector< const ColVariable * > idx_to_dvar;
                                         ///< from index to dynamic variable

 std::vector< con_int > dcon_to_idx;     ///< from dynamic constraint to index
 std::vector< const FRowConstraint * > idx_to_dcon;
                                         ///< from index to dynamic constraint

/** @} ---------------------------------------------------------------------*/
/*--------------------- FIELDS FOR PROBLEM DESCRIPTION ---------------------*/
/*--------------------------------------------------------------------------*/
 /** @name MILP problem description
  *
  * The following fields are used to describe the MILP problem in the
  * standard "sparse matrix as a vector of doubles + two vectors of int,
  * plus vectors for costs, bounds and lhs/rhs of constraints" format.
  *
  * #matbeg, #matcnt, #matind and #matval define the (sparse) constraint
  * matrix by its nonzero coefficients. These are grouped by column in the
  * array matval. The nonzero elements of every column must be stored in
  * sequential locations in this array with matbeg[ j ] containing the index
  * of the beginning of column j and matcnt[ j ] containing the number of
  * entries in column j. The components of matbeg must be in ascending
  * order. For each k, matind[ k ] specifies the row number of the
  * corresponding coefficient, matval[ k ].
  * @{  */

 /** Pointers to the currently registered Block and all its descendants
  * arranged along a BFS order. */
 std::vector< Block * > v_BFS;

 std::string prob_name;    ///< problem name
 std::string output_file;  ///< output file

 /** An integer that specifies the number of columns in the constraint matrix,
  * or equivalently, the number of variables in the problem object. */
 int numcols{};

 /** An integer that specifies the number of rows in the constraint matrix,
  * not including the objective function or bounds on the variables. */
 int numrows{};

 /** An integer that specifies whether the problem is a minimization or
  * maximization problem. */
 int objsense{};

 /// A double that specify the summation of the constant terms of all blocks.
 OFValue constant_value{};

 /** An array of length at least numcols containing the objective function
  * coefficients. */
 std::vector< double > objective;

 /** An array of length numcols containing the quadratic coefficients of
  * the separable quadratic objective. */
 std::vector< double > q_objective;

 /** An array of length at least numrows containing the righthand side value
  * for each constraint in the constraint matrix. */
 std::vector< double > rhs;

 /** An array of length at least numrows containing the range value of each
  * ranged constraint. Ranged rows are those designated by 'R' in the sense
  * array. If the row is not ranged, the rngval array entry is ignored.
  * If rngval[ i ] > 0, then row i activity is in
  * [ rhs[ i ] , rhs[ i ] + rngval[ i ] ], while if rngval[ i ] <= 0 then
  * row i activity is in [ rhs[ i ] + rngval[ i ] , rhs[ i ] ] */
 std::vector< double > rngval;

 /** An array of length at least numrows containing the sense of each constraint
  * in the constraint matrix. */
 std::vector< char > sense;

 std::vector< int > matbeg;     ///< Beginnings of constraint matrix columns
 std::vector< int > matcnt;     ///< Sizes of constraint matrix columns
 std::vector< int > matind;     ///< Indices of rows for each coefficient
 std::vector< double > matval;  ///< All nonzero coefficients

 /** An array of length at least numcols containing the lower bound on each
  * of the variables. */
 std::vector< double > lb;

 /** An array of length at least numcols containing the upper bound on each
  * of the variables. */
 std::vector< double > ub;

 /** An array of length numcols containing the type of each column in
  * the constraint matrix. Possible values:
  *
  * | Value | Type of variable |
  * | :---: | ---------------- |
  * |  'C'  | continuous       |
  * |  'B'  | binary           |
  * |  'I'  | general integer  |
  * |  'S'  | semi-continuous  |
  * |  'N'  | semi-integer     | */
 std::vector< char > xctype;

 /// if true, use Variable/Constraint custom names
 bool use_custom_names = true;

 /// if true, relax [M]ILP by removing integrality constraints
 bool relax_int_vars = false;

 /** An array of length at least numcols containing pointers to character
  * strings containing the names of the variables. */
 std::vector< char * > colname;

 /** An array of length at least numrows containing pointers to character
  * strings containing the names of the constraints. */
 std::vector< char * > rowname;

 int sol_status = kUnEval; ///< Solution status (OK, Infeasible, Unbounded, ...)
 int int_vars{};           ///< Number of integer variables
 int static_vars{};        ///< Number of static variables
 int static_cons{};        ///< Number of static constraints

/** @} ---------------------------------------------------------------------*/
/*--------------------------- PROTECTED METHODS ----------------------------*/
/*--------------------------------------------------------------------------*/
 /** @name Get variable bounds for the problem
  *
  * The following two methods retrieve the upper and lower bound for the
  * given variable considering both the Variable bounds and all the active
  * OneVarConstraints active for that Variable.
  * @{ */

 /// gets the LB for the given variable in the problem
 virtual double get_problem_lb( const ColVariable & var ) const;

 /// gets the UB for the given variable in the problem
 virtual double get_problem_ub( const ColVariable & var ) const;

 /// gets both bounds for the given variable in the problem
 virtual std::array< double , 2 > get_problem_bounds(
					     const ColVariable & var ) const;

/** @} ---------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
 /// returns true if b is "mine" (f_Block or one of its descendants)

 bool is_mine( Block * b ) const {
  while( b ) {
   if( b == f_Block )
    return( true );
   b = b->get_f_Block();
   }
  return( false );
  }

/*--------------------------------------------------------------------------*/
 /// gets the active constraints for the specified variable
 // TODO: This should be temporary
 std::vector< FRowConstraint * > get_active_constraints(
					     const ColVariable & var ) const;

 /// gets the active bounds for the specified variable
 // TODO: This should be temporary
 std::vector< OneVarConstraint * > get_active_bounds(
					     const ColVariable & var ) const;

/*--------------------------------------------------------------------------*/
/*----------------- INTERFACE FOR SUPPORTING MODIFICATIONS ---------------- */
/*--------------------------------------------------------------------------*/
/** @name Methods for modifying the constructed problem
 *
 * This is the API for supporting Modifications.
 * The method process_modifications() is the one that checks the Modification
 * queue for pending modifications and calls the proper method.
 * A derived class that wants to support modifications should implement the
 * virtual methods of this group.
 * @{ */

 /// processes all the pending modifications
 void process_modifications( void );

 /// process all not-GroupModification (bulk of the work)
 void guts_of_process_modifications( const p_Mod mod );

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// Checks if the given function is an objective function
 /** This method is meant to be used by process_modifications() when a
  * FunctionMod is catched, in order to discriminate between a modification
  * of the objective function and one of a constraint.
  * The method checks against the block's and all the sub-blocks' OFs.
  *
  * @param f a function
  * @return true if the function belongs to the objective, false otherwise
  */
 static bool is_of( Function * f );

 /// handles a variable modification
 virtual void var_modification( const VariableMod * mod );

 /// handles an objective modification
 virtual void objective_modification( const ObjectiveMod * mod );

 /// handles a constraint modification
 virtual void const_modification( const ConstraintMod * mod );

 /// handles a bound modification
 virtual void bound_modification( const OneVarConstraintMod * mod );

 /// handles a function modification applied to the objective
 virtual void objective_function_modification( const FunctionMod * mod );

 /// handles a function modification applied to a constraint
 virtual void constraint_function_modification( const FunctionMod * mod );

 /// handles a function vars modification to the objective
 virtual void objective_fvars_modification( const FunctionModVars * mod );

 /// handles a function vars modification to a constraint
 virtual void constraint_fvars_modification( const FunctionModVars * mod );

 /// handles a dynamic modification
 virtual void dynamic_modification( const BlockModAD * mod );

 /// adds a single new dynamic constraint
 /** Notice that empty constraints, i.e., constraints with null function, are by
  * definition equals to zero, so as in some cases it might be useful to
  * handle them, if FRowConstraint::get_function() returns nullptr, then the
  * constraint will be added and the respective row will be generated since,
  * formally speaking, an empty constraint is linear since the identical
  * function zero is.
  *
  * @param con a reference to a FRowConstraint
  */
 virtual void add_dynamic_constraint( const FRowConstraint * con );

 /// adds a single new dynamic variable
 virtual void add_dynamic_variable( const ColVariable * var );

 /// adds a single new dynamic bound
 virtual void add_dynamic_bound( const OneVarConstraint * con );

 /// removes a single dynamic constraint
 virtual void remove_dynamic_constraint( const FRowConstraint * con );

 /// removes a single dynamic variable
 virtual void remove_dynamic_variable( const ColVariable * var );

 /// removes a single dynamic bound
 virtual void remove_dynamic_bound( const OneVarConstraint * con );

/** @} ---------------------------------------------------------------------*/
/*--------------------- PRIVATE FIELDS OF THE CLASS ------------------------*/
/*--------------------------------------------------------------------------*/

 private:

/*--------------------------------------------------------------------------*/
/*------------- AUXILIARY METHODS FOR POPULATING THE PROBLEM  --------------*/
/*--------------------------------------------------------------------------*/
/** @name Private auxiliary methods
 *
 * These methods are used in load_problem() to read data from the Block.
 * All of them, except scan_objective(), take integer counters as input
 * parameters. That's because they are meant to be used by
 * un_any_const_static() and un_any_const_dynamic() template functions
 * on boost::any containers, and the counters keep track of the elements
 * inside the containers.
 * @{ */

 /** Scans a static ColVariable and fills the dictionaries accordingly
  *
  * @param var a reference to a ColVariable
  * @param n   an counter that should be 0 when var is the first
  *            element of a vector of static ColVariables
  * @param col a counter for variables/columns */

 void scan_static_variable( const ColVariable & var , Index & n ,
			    Index & col );

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /** Scans a dynamic ColVariable and fills the dictionaries accordingly
  *
  * @param var a reference to a ColVariable
  * @param col a counter for variables/columns */

 void scan_dynamic_variable( const ColVariable & var , Index & col );

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// common part of scan_static_variable() and scan_dynamic_variable()

 void scan_variable( const ColVariable & var , Index & col );

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /** Scans a static FRowConstraint and fills the dictionaries accordingly.
  *
  * Notice that empty constraints, i.e., constraints with null function, are by
  * definition equals to zero, so as in some cases it might be useful to
  * handle them, if FRowConstraint::get_function() returns nullptr, then the
  * constraint will be considered since, formally speaking, an empty
  * constraint is linear since the identical function zero is.
  *
  * @param con a reference to a FRowConstraint
  * @param n   an counter that should be 0 when lconst is the first
  *            element of a vector of static FRowConstraints
  * @param row a counter for constraints/rows */

 void scan_static_constraint( const FRowConstraint & con , Index & n ,
			      Index & col );

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /** Scans a dynamic FRowConstraint and fills the dictionaries accordingly.
  *
  * Notice that empty constraints, i.e., constraints with null function, are by
  * definition equals to zero, so as in some cases it might be useful to
  * handle them, if FRowConstraint::get_function() returns nullptr, then the
  * constraint will be considered since, formally speaking, an empty
  * constraint is linear since the identical function zero is.
  *
  * @param con a reference to a FRowConstraint
  * @param row a counter for constraints/rows */

 void scan_dynamic_constraint( const FRowConstraint & con , Index & row );

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// common part of scan_static_constraint() and scan_dynamic_constraint()

 void scan_constraint( const FRowConstraint & con , Index & row );

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /** Scans a FRealObjective and fills the vectors of the LP accordingly.
  * Moreover, since both the CPLEX and SCIP C API does not support the concept
  * of "constant term", all of them, for each Block of the problem, are
  * accumulated in the homonymous variable to provide the updated OF value.
  *
  * @param obj a FRealObjective */

 void scan_objective( const FRealObjective * obj );

/** @} ---------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

 SMSpp_insert_in_factory_h;

/*--------------------------------------------------------------------------*/
 
 };  // end( class( MILPSolver ) )

/*--------------------------------------------------------------------------*/

}  // end( namespace SMSpp_di_unipi_it )

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

#endif  /* MILPSolver.h included */

/*--------------------------------------------------------------------------*/
/*------------------------- End File MILPSolver.h --------------------------*/
/*--------------------------------------------------------------------------*/
