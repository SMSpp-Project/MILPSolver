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
 * \copyright &copy; by Enrico Calandrini, Antonio Frangioni,
 *                   Niccolo' Iardella
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

#include <QuadFunction.h>

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
  /// throws exception if there is inconsistency when storing a reduced cost
  intThrowReducedCostException = intLastParCDAS ,
  intUseCustomNames , ///< use custom names for rows/columns
  /// Relax [M]ILP by removing integrality constraints for integer variables
  /**< Controls how the integrality of integer Variable is treated. Three
   * semantics are supported:
   *
   * - intRelaxIntVars == 0 [default]: the problem is solved as a MILP, with
   *   integer Variable kept integer.
   *
   * - intRelaxIntVars == 1: the problem is solved as the LP relaxation,
   *   i.e., integer Variable are treated as continuous. The Solver uses
   *   the MILP engine on a problem whose integrality has been relaxed
   *   (the cut callback, if any, still fires as for a MIP). This is
   *   appropriate when the user wants the LP relaxation but the Solver's
   *   MIP machinery (in particular the user-cut callback at the root
   *   node) is required.
   *
   * - intRelaxIntVars == 2: the problem is solved as a pure LP with an
   *   explicit user-cut separation loop driven by MILPSolver. After each
   *   LP solve, the Block's generate_dynamic_constraints() is called and,
   *   if any new dynamic Constraint is added, the LP is re-solved; the
   *   loop terminates when no new constraint is generated or
   *   intMaxCutPasses rounds have been performed. This mode does *not*
   *   install the MIP cut callback (which would not fire on a pure LP
   *   solve) and is therefore the only one suitable for getting the LP
   *   relaxation value of a Block with non-trivial dynamic constraints
   *   when the underlying Solver does not natively support user-cut
   *   separation on continuous problems. */
  intRelaxIntVars ,
  intSingleBound , // Force that at maximum one OneVarConstraint can be
                   // associated to a single variable
  intConsModification , // Enable/Disable constraint modifications
  /// parameter for deciding if/when cut separation is done
  /**< If > 0, the Solver may perform user-cut and/or lazy-constraint
   * separation during compute(), calling generate_dynamic_constraints()
   * on the attached Block. The specific bit-coded semantics (which
   * sub-events fire the separation, when in the search tree, ...) is
   * defined by the derived class; the base class only inspects whether
   * the value is non-zero, which has two effects:
   *
   * - when non-zero, compute() holds the Block under a write lock() for
   *   its full duration (since the cut callback may modify the Block);
   *   when zero, compute() only acquires a read_lock() for the duration
   *   of process_modifications() and then releases it.
   *
   * - when non-zero, the loop driven by intRelaxIntVars == 2 (see) is
   *   executed (only if intRelaxIntVars == 2; otherwise the derived
   *   class governs the callback installation as appropriate). */
  intCutSepPar ,
  /// maximum number of LP cut separation passes (when intRelaxIntVars == 2)
  /**< Caps the number of "solve LP - separate user cuts - re-solve"
   * iterations performed by the base compute() when intRelaxIntVars == 2.
   * Defaults to 1000. Increase if the separator is known to require very
   * many passes to converge; reduce to limit the time spent in the
   * separation loop when fast (possibly weaker) bounds are acceptable. */
  intMaxCutPasses ,
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
  strWarmStartSolution ,             ///< warm start solution filename
  strWarmStartVariables ,            ///< warm start variables filename
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

 using v_off_diag_term = QuadFunction::v_off_diag_term; 
                              // "import" v_off_diag_term from QuadFunction
 using Qmat = QuadFunction::Qmat; // "import" Qmat from QuadFunction
 using Coefficient = QuadFunction::Coefficient;

/*--------------------------------------------------------------------------*/
/*--------------------- CONSTRUCTOR AND DESTRUCTOR -------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Constructor and destructor
 * @{ */

 MILPSolver( void ) : CDASolver() , throw_reduced_cost_exception( 0 ) {}

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

/*--------------------------------------------------------------------------*/
 /// methods to scan a single group of Constraints or Variables

 /** Scans a "simple" static group of FRowConstraint or ColVariable.
 *
 * This function is called from load_problem() whenever a new group of
 * FRowConstraint or ColVariable is encountered. The group should be contained
 * in an "easy" structure (e.g. T *, std::vector< T* >).
 *
 * @param gr        Reference to the group being scanned.
 * @param qb        The block from which the group originated.
 * @param num_block Sequential number of the block in MILPSolver.
 * @param set       Index of the constraint set within qb.
 * @param counter   Current number of elements of type T that have been
 *                  scanned.
 * @param T         The element type of the group. Should be either 
 *                  FRowConstraint or ColVariable, depending on the group 
 *                  being scanned.
 */
 
 template< typename T >
 void scan_st_group( const boost::any & gr , Block * qb , Index num_block ,
                    Index set , Index & counter , un_any_type< T > );

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 /** Scans a "complex" static group of FRowConstraint or ColVariable.
 *
 * This function is called from load_problem() whenever a new group of
 * FRowConstraint or ColVariable is encountered. The group should be contained
 * in a more "complex" structure (e.g. multi_array< T* >, 
 * multi_array< std::vector < T* > >).
 *
 * @param gr        Reference to the group being scanned.
 * @param qb        The block from which the group originated.
 * @param num_block Sequential number of the block in MILPSolver.
 * @param set       Index of the constraint set within qb.
 * @param counter   Current number of elements of type T that have been
 *                  scanned.
 * @param T         The element type of the group. Should be either 
 *                  FRowConstraint or ColVariable, depending on the group 
 *                  being scanned.
 */
 
 template< typename T >
 void scan_multiarray_st_group( const boost::any & gr , Block * qb ,
            Index num_block , Index set , Index & counter , un_any_type< T > );

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// common part of scan_st_group() and scan_multiarray_st_group()

 template< typename T >
 void scan_group( const boost::any & gr , Block * qb , Index num_block ,
                  Index set , Index & row , un_any_type< T > );

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

  /// returns the number of quadratic constraints/rows
 [[nodiscard]] int get_numquadrows( void ) const { return( numquadrows ); }

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
 /** @name Informative methods about different aspect of the problem
 *
 * The following methods are used to provide useful information about the
 * internal status of the problem (e.g. number of nodes explored ).
 * These methods are not properly implemented in the base class MILPSolver
 * as they should be overwritten in derived classes.
 * @{ */

 /// Returns the number of nodes explored so far.
 [[nodiscard]] virtual int get_explored_nodes( void ) const {
  throw( std::runtime_error( "Function get_explored_nodes is not supported "
    "by the current *MILPSolver" ) );
  
  return( 0 );
  }

 /// Returns the estimated number of nodes to explore.
 [[nodiscard]] virtual long get_left_nodes( void ) const {
  throw( std::runtime_error( "Function get_left_nodes is not supported "
    "by the current *MILPSolver" ) );
  
  return( 0 );
  }

 /// Returns a true value if a feasible solution is known, 
 //  false otherwise.
 [[nodiscard]] virtual bool has_feasible_sol( void ) {
  throw( std::runtime_error( "Function has_feasible_sol is not supported "
    "by the current *MILPSolver" ) );
  
  return( 0 );
  }

 /// Returns elapsed solver runtime (in second).
 [[nodiscard]] virtual double get_runtime( void ) const {
  throw( std::runtime_error( "Function get_runtime is not supported "
    "by the current *MILPSolver" ) );
  
  return( 0 );
  }

 /// Returns a unique identifier for the node currently being explored  
 //  in the branch-and-bound algorithm for a MIP problem.  
 //  
 /// NOTE: This method should only be called during the callback process  
 //  and in specific situations (e.g., when a new incumbent solution is found,  
 //  and you need to identify the node from which it originates).  
 [[nodiscard]] virtual long get_id_node( void ) const {
  throw( std::runtime_error( "Function get_id_node is not supported "
    "by the current *MILPSolver" ) );
  
  return( 0 );
  }

 /** 
 * Adds multiple MIP starts to a MIP problem. This function allows the solver 
 * to receive multiple sets of starting values by providing vectors of variable 
 * indices and corresponding values for each start.
 * 
 * NOTE: Partial solutions are allowed. In such cases, the solver will attempt 
 * to infer values for the unspecified variables.
 */
 virtual void add_mip_starts( 
  std::vector< std::vector< int > > varidxs ,
  std::vector< std::vector< double > > varvalues ) {

    throw( std::runtime_error( "Function add_mip_starts is not supported "
      "by the current *MILPSolver" ) );

    return;
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

/** 
 * Adds a new warm start to the solver. See Solver.h for further details.
 * 
 * NOTE: Partial solutions are allowed. In such cases, the solver will attempt 
 * to infer values for the unspecified variables.  */
 //void add_warm_start( const Solution * sol ,
  //                    const std::vector< AbstractPath > varpaths ) override;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

/** 
 * Adds a set of warm start to the solver. See Solver.h for further details.
 * 
 * NOTE: Partial solutions are allowed. In such cases, the solver will attempt 
 * to infer values for the unspecified variables.  */
//void add_warm_start( const std::vector< Solution * > sols ,
//                     const std::vector< AbstractPath > varpaths ) override;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 #ifdef MILPSolver_DEBUG
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
 /// compute() entry point shared by all the *MILPSolver derived classes
 /** The base compute() takes care of:
  *
  * - acquiring the appropriate locks (Solver mutex; on the Block, write
  *   lock() if cut separation may happen during the solve, otherwise
  *   read_lock() only for the duration of process_modifications());
  *
  * - dispatching pending Modification (via process_modifications());
  *
  * - dispatching to the derived class via guts_of_compute() to perform
  *   the actual back-end solve;
  *
  * - if intRelaxIntVars == 2, wrapping guts_of_compute() in an LP
  *   cut-separation loop: after each solve the LP solution is written
  *   into the Block Variable (via get_var_solution()),
  *   generate_dynamic_constraints() is called on the Block, and the
  *   resulting new Modification (if any) is processed before the next
  *   guts_of_compute() iteration. The loop terminates when no new
  *   Modification is produced or intMaxCutPasses rounds have been done.
  *
  * Derived classes should *not* override compute(); they should instead
  * override guts_of_compute(). */

 int compute( bool changedvars = true ) override;

/*--------------------------------------------------------------------------*/
 /// does nothing as there is nothing to do
 void get_var_solution( Configuration * solc ) override {}

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// takes a get_numcols()-vector of doubles and writes it as var solution
 /** The get_numcols()-vector \p x is supposed to encode a var solution in the
  * natural format, i.e., x[ i ] is the value of the ColVariable corresponding
  * to the i-th column in the coefficient matrix as constructed by
  * MILPSolver; this method writes it in the Block. It can obviously be used
  * to implement get_var_solution(). */
 
 void write_var_solution( const std::vector< double > & x );

/*--------------------------------------------------------------------------*/
 /// does nothing as there is nothing to do
 void get_dual_solution( Configuration * solc ) override {}

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// takes two vectors of doubles and writes them as dual solution
 /** The get_numrows()-vector \p pi and the get_numcols()-vector \p rc are
  * supposed to encode dual a solution in the natural format, i.e., pi[ j ]
  * is the value of the dual variable of to the FRowConstraint corresponding
  * to the j-th row in the coefficient matrix as constructed by
  * MILPSolver, while rc[ i ] is the value of the reduced cost of the
  * ColVariable corresponding to the i-th column in the coefficient matrix
  * as constructed by MILPSolver; this method writes them in the Block if
  * they are provided, i.e., both pi and rc may be empty(), in which case
  * they are ignored.
  *
  * Note that the reduced cost part is nontrivial, in that there is not
  * really such a thing as the reduced cost of a ColVariable; this has to
  * be implemented as the dual variable of a BoxConstraint using that
  * ColVariable. Currently, the method just looks if there is one, and if
  * not it ignores the thing; maybe we'll add a general mechanism to throw
  * exception instead.
  *
  * This method can obviously be used to implement get_dual_solution(). */
 
 void write_dual_solution( const std::vector< double > & pi ,
			   const std::vector< double > & rc );

/** @} ---------------------------------------------------------------------*/
/*------------------- METHODS FOR HANDLING THE PARAMETERS ------------------*/
/*--------------------------------------------------------------------------*/
/** @name Methods for handling parameters
 * @{ */

 /// sets an integer parameter with the given value
  /** Set the "int" parameters specific of MILPSolver, together with the
  * parameters of MILPSolver that various solver actually "listens to":
  *
  * - intThrowReducedCostException [0]: it indicates whether an exception must
  *                                     be thrown if there is an inconsistency
  *   when a reduced cost is being stored during a call to get_dual_solution()
  *   or get_dual_direction(). The reduced cost of a Variable is stored in at
  *   most one OneVarConstraint on that Variable. It may happen that a
  *   Variable has no OneVarConstraint, in which case its reduced cost will
  *   not be stored and will be lost. Usually, the reduced cost of a Variable
  *   is of interest if the Variable has a finite nonzero lower or upper
  *   bound. In this case, if a OneVariableConstraint for that Variable is not
  *   found, an exception is thrown. More specifically, there are two cases in
  *   which an exception is thrown:
  *
  *   1) The Variable is fixed to a finite nonzero value and there is no
  *      OneVarConstraint on that Variable whose lower and upper bounds are
  *      both equal to the value of that Variable.
  *
  *   2) The Variable is not fixed, it has a finite nonzero lower or upper
  *      bound and there is no OneVarConstraint on that Variable whose lower
  *      or upper bound match the bounds of the Variable. */
 void set_par( idx_type par , int value ) override;

 /// sets a double parameter with the given value
 // although it actually does nothing, it has to be there since, due to an
 // excess of caution, *MILPSolver::set_par( double ) call it, and if it's
 // not there then set_par( int ) ends up being called with duuble -> int
 // conversion that can go awry
 void set_par( idx_type par , double value ) override {
  CDASolver::set_par( par, value );
  }

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
/*------------------------- PROTECTED METHODS ------------------------------*/
/*--------------------------------------------------------------------------*/

 /// derived-class hook performing the actual back-end solve
 /** Virtual method called by MILPSolver::compute() after the locks have
  * been acquired and any pending Modification have been processed. The
  * derived class implements only the part that is specific to its
  * underlying solver: calling the LP/MIP optimisation function, decoding
  * the return code, and storing it into f_status or the analogous member.
  *
  * The default base-class implementation returns kOK and does nothing,
  * which makes MILPSolver usable as a stand-alone class for the sole
  * purpose of keeping the matrix-based representation in sync with the
  * Block (no actual optimisation is performed). */

 virtual int guts_of_compute( void ) { return( kOK ); }

/*--------------------------------------------------------------------------*/
/*---------------- VARIABLE AND CONSTRAINT TRACKING VECTORS ----------------*/
/*--------------------------------------------------------------------------*/

 using var_int = std::pair< const ColVariable * , int >;
 using int_var = std::pair< int , const ColVariable * >;
 using con_int = std::pair< const FRowConstraint * , int >;
 using int_con = std::pair< int , const FRowConstraint * >;
 using var_int_int = std::tuple< const ColVariable * , int , int >;
 using con_int_int = std::tuple< const FRowConstraint * , int , int >;
 using c_v_coeff_pair = DQuadFunction::c_v_coeff_pair;
 using mat_indices = std::pair< int , int >;

 /** @name Variable and Constraint dictionaries
  *
  * The following vectors are used in order to keep track between the
  * Variables and Constraints of the Block and the constraint matrix.
  *
  *  - svar_to_idx, scon_to_idx, scon_to_idx : vectors of tuples 
  *    that store 1) the address of the first element of each group of static 
  *    variables and constraints, respectively, 2) the 
  *    corresponding index in constraint matrix (column or row) and 3) 
  *    the number of elements in the group.
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
  *  - svar_to_bound : vector of single OneVarConstraint * associated to 
  *    each static variable. The order of the variables used to sort the 
  *    array is the one given by Block::get_static_variables.
  *    NOTE: this is used only if the option intSingleBound is set to 1.
  *
  *  - dvar_to_bound : vector of single OneVarConstraint * associated to 
  *    each dynamic variable. The order of the variables used to sort the 
  *    array is the one given by Block::get_static_variables.
  *    NOTE: this is used only if the option intSingleBound is set to 1.
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
 std::vector< int_con> idx_to_scon;      ///< from index to linear static constraint

 std::vector< var_int > dvar_to_idx;     ///< from dynamic variable to index
 std::vector< const ColVariable * > idx_to_dvar;
                                         ///< from index to dynamic variable

 std::vector< con_int > dcon_to_idx;     ///< from dynamic constraint to index
 std::vector< const FRowConstraint * > idx_to_dcon;     
                                        ///< from index to dynamic constraint

 std::vector< const OneVarConstraint * > svar_to_bound; 
                                         ///< from static variable to bound
 std::vector< const OneVarConstraint * > dvar_to_bound; 
                                         /// from dynamic variable to bound

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
  * matrix associated to linear constraints by its nonzero coefficients. 
  * These are grouped by column in the array matval. The nonzero elements 
  * of every column must be stored in sequential locations in this array 
  * with matbeg[ j ] containing the index of the beginning of column j 
  * and matcnt[ j ] containing the number of entries in column j. 
  * The components of matbeg must be in ascending order. For each k, 
  * matind[ k ] specifies the row number of the corresponding coefficient, 
  * matval[ k ].
  *
  * NOTE: the above mentioned structures are grouped by column when no 
  * quadratic constraint is in the initial load of the model. Otherwise, 
  * the representation is switched to the rows.
  *
  * For a quadratic constraint i, we separtely store all the nonzeros linear
  * terms in the above mentioned matbeg, matcnt, ... structures (grouped by
  * rows). The quadratic part of the constraint is instead represented with a
  * map< mat_indices,float > stored in the i-th position of the vector q_part.
  * 
  * The same procedure is applied for the objective function, with the linear
  * coefficients stored in objective, the diagonal coefficients of the quadratic
  * matrix stored in q_objective and the off-diagonal ones stored using three 
  * vectors ndq_rowind, ndq_colind, ndq_objective.
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

 /** An integer that specifies the number of quadratic rows in the set of
  * constraints, not including the objective function or bounds on the 
  * variables. */
 int numquadrows{};

 /** An integer that specifies the number of nonzero coefficients outside
  * the diagonal of the quadratic objective matrix. */
 int numnnzq{};

 /** An integer that specifies whether the problem is a minimization or
  * maximization problem. */
 int objsense{};

 /// A double that specify the summation of the constant terms of all blocks.
 OFValue constant_value{};

 /** An array of length at least numcols containing the linear objective function
  * coefficients. */
 std::vector< double > objective;

 /** An array of length numcols containing the quadratic coefficients along
  * the diagonal of the quadratic objective matrix. */
 std::vector< double > q_objective;

 /** An array of length numnnzq containing the quadratic coefficients outside
  * the diagonal of the quadratic objective matrix. */
 std::vector< double > ndq_objective;

 std::vector< int > ndq_rowind; ///< Indices of rows for each quadratic coefficient
 std::vector< int > ndq_colind; ///< Indices of columns for each quadratic coefficient

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

 // Linear constraints container
 std::vector< int > matbeg;     ///< Beginnings of constraint matrix columns
 std::vector< int > matcnt;     ///< Sizes of constraint matrix columns
 std::vector< int > matind;     ///< Indices of rows for each coefficient
 std::vector< double > matval;  ///< All nonzero coefficients

 /* Quadratic constraints container */ 

 // Linear coefficients (already stored in matbeg,...)

 // Quadratic coefficients
 std::vector< std::map< mat_indices , float > > q_part;

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

 /// how to handle the integrality of integer Variable
 /**< 0 = solve as MILP, 1 = MILP engine on the relaxed problem, 2 = LP
  * with explicit user-cut separation loop; see the comments to the
  * intRelaxIntVars parameter for details. */
 int relax_int_vars = 0;

 /// parameter for deciding if/when cut separation is done; see intCutSepPar
 unsigned char CutSepPar = 0;

 /// maximum number of LP cut separation passes (see intMaxCutPasses)
 int max_cut_passes = 1000;

 /* if true, no more than one OneVarConstraint can be associated to a
 *  single variable. 
 *  Moreover, the vectors svar_to_bound and dvar_to_bound are activated
 *  to guarantee a direct link between variables and bound. */
 bool single_bound = false;

  /* if true, modification on constraints are enabled. This parameter can
  *  be useful when dealing with quadratic constraint, where Modification 
  * from some Solver (e.g. CPLEX) are not allowed. */
 bool cons_modification = true;

 /** This variable indicates whether an exception must be thrown if there is
 * an inconsistency when a reduced cost is being stored during a call to
 * get_dual_solution() or get_dual_direction(). */
 bool throw_reduced_cost_exception;

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
 int static_quadcons{};    ///< Number of static quadratic constraints

 /** Warm start structures used by the solver.
 *
 * These fields store the initial solution data that can be provided to
 * the solver. As described in Solver.h, each warm start is represented
 * as a tuple: (Solution*, std::vector<AbstractPath>).
 *
 * Since multiple warm starts may be supplied, two vectors are maintained:
 * - v_warmstart_sol: stores the Solution* instances
 * - v_warmstart_vars: stores the AbstractPath vectors
 *
 * @note A single std::vector<AbstractPath> can be shared across multiple
 * Solution * instances. To track this association, a std::vector<int> 
 * — matching the length of v_warmstart_sol— is used to record which :AbstractPath
 * vector is referenced by each Solution *.
 *
 * Currently, warm starts can be provided in two ways:
 *
 * 1. **Programmatically**, by storing the necessary data structures
 *    directly in the Solver instance using the methods `add_warm_start()`
 *    or `add_warm_starts()`.
 *
 * 2. **Via input files**, using two files of type `eWarmStartFile` and
 *    `eSolutionFile`. The corresponding filenames can be set using the
 *    string parameters:
 *      - "strWarmStartVariables" — path to the file defining set of variables
 *      for which we would like to provide a warm start.
 *      - "strWarmStartSolution"  — path to the file containing the
 *      warm start solution (i.e. initial values for the specified set of
 *      variables).
 *
 *    If a filename follows the format "filename[idx]", only the indexed 
 *    structure within the file will be used. 
 *    (Behavior for unspecified indices or other formats is TBD.)
 */

 std::string warmstart_variables; // warm start variables filename
 std::string warmstart_solution;  // warm start solution filename

 //std::vector< std::vector< AbstractPath >> v_warmstart_vars;
                                  // warm start variables
 //std::vector< Solution * > v_warmstart_sol;               
                                  // warm start solutions
 
 //std::vector< int > v_warmstart_sol2vars; 
                                  // corresponding variables for each solution

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
					     const ColVariable & var ,
               bool first_scan = false ) const;

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
  * @param n a counter that should be 0 when row is the first
  *            element of a vector of linear static FRowConstraints
  * @param row a counter for constraints/rows */

 void scan_static_constraint( const FRowConstraint & con , Index & n,
			      Index & row );

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
 /** Scans a static ColVariable to found the associated bound and fills the 
  *  dictionaries accordingly
  *
  * @param var a reference to a ColVariable */

 void scan_static_variable_bound( const ColVariable & var );

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /** Scans a dynamic ColVariable to found the associated bound and fills the 
  *  dictionaries accordingly
  *
  * @param var a reference to a ColVariable */

 void scan_dynamic_variable_bound( const ColVariable & var );

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /** Scans a FRealObjective and fills the vectors of the LP accordingly.
  * Moreover, since both the CPLEX and SCIP C API does not support the concept
  * of "constant term", all of them, for each Block of the problem, are
  * accumulated in the homonymous variable to provide the updated OF value.
  *
  * @param obj a FRealObjective */

 void scan_objective( const FRealObjective * obj );

/*--------------------------------------------------------------------------*/
/*--------------- AUXILIARY METHODS FOR MULTI-ARRAY GROUP  -----------------*/
/*--------------------------------------------------------------------------*/
/** @name Multi-array methods
 *
 * These methods are used in load_problem() to read data from complex
 * multi_array<> structures.
 * Each method is templated with:
 *  1) T - the type of elements in the group, expected to be either
 *     ColVariable or FRowConstraint.
 *  2) K - the number of dimensions of the multi_array.
 *
 * NOTE: Currently, only 2D or 3D arrays are supported.
 * @{ */

 /** Scans a multi_array structure and returns its number of dimensions.
  *
  * This method attempts to cast a boost::any element to a boost::multi_array.
  * It should be used as a recursive method, as it will try to cast an 
  * increasing number of dimensions until the cast succeeds.
  * If the cast succeeds, it returns the number of dimensions of the array.
  * A default maximum of K = 9 dimensions is used when attempting the cast.
  * 
  * @param any the reference to the multi_array
  * @param T the basic type of the multi_array
  * @param K the number of dimensions of the multi array 
 */

 template< typename T , unsigned short K >
  int get_multi_array_dim( const boost::any & any ,
                           un_any_type< T > , un_any_int< K > );

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 template< typename T >
  int get_multi_array_dim( const boost::any & ,
                          un_any_type< T > , un_any_int< 9 > ) {
  return( -1 );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 /** Scans a multi_array structure and returns its type.
  *
  * This method attempts to cast a boost::any element to a boost::multi_array
  * with fixed number of dimensions K.
  * If the cast succeeds, it returns the type of the array.
  * In SMS++ currently two different types of multi_array are available:
  *
  * - boost::multi_array< T > -> type 0
  * - boost::multi_array< std::vector < T > > -> type 1
  * 
  * @param any the reference to the multi_array
  * @param T the basic type of the multi_array
  * @param K the number of dimensions of the multi array 
 */
 template< typename T , unsigned short K >
  int get_multi_array_type( 
                         const boost::any & any ,
                         un_any_type< T > , 
                         un_any_int< K > );

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 /** These methods attempt to cast a multi_array with specific attributes.
 * If the cast is successful, they return a pointer to the resulting 
 * structure. 
 * 
 * @param any the reference to the multi_array
 * @param T the basic type of the multi_array
 * @param K the number of dimensions of the multi array */

 // Cast to a 2D multi_array of type 0
 template< typename T >
 boost::multi_array< T , 2 > * get_multi_array0( 
                          const boost::any & ,
                          un_any_type< T > , un_any_int< 2 > );

 // Cast to a 2D multi_array of type 1
 template< typename T >
 boost::multi_array< std::vector< T >, 2 > * get_multi_array1( 
                          const boost::any & ,
                          un_any_type< T > , un_any_int< 2 > );

 // Cast to a 3D multi_array of type 0
 template< typename T >
 boost::multi_array< T , 3 > * get_multi_array0( 
                          const boost::any & ,
                          un_any_type< T > , un_any_int< 3 > );

 // Cast to a 3D multi_array of type 1
 template< typename T >
 boost::multi_array< std::vector< T >, 3 > * get_multi_array1( 
                          const boost::any & ,
                          un_any_type< T > , un_any_int< 3 > );

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
