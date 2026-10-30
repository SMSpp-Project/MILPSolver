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
 * \author Donato Meoli \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \copyright &copy; by Enrico Calandrini, Antonio Frangioni,
 *                   Niccolo' Iardella,
 *                   Donato Meoli
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

#include <LinearFunction.h>

#include <OneVarConstraint.h>

#include <QuadFunction.h>

#include <algorithm>

#include <array>

#include <atomic>

#include <cmath>

#include <cstddef>

#include <functional>

#include <limits>

#include <map>

#include <typeinfo>

#include <unordered_map>

#include <unordered_set>

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
 * The derived classes implement get_var_solution(), get_dual_solution() and
 * the two methods of the directions with write_var_solution() and
 * write_dual_solution(), to which they hand the Configuration they are
 * given: a SimpleConfiguration< std::vector< int > > listing sub-Blocks of
 * the Block, and -1 for the Block itself, has only that part of the solution
 * written, any other Configuration all of it [see write_var_solution()].
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
   * - intRelaxIntVars == 1: the problem is solved as its continuous
   *   relaxation: the integer Variable are given to the Solver as
   *   continuous ones, so that the model has no integer column and it is
   *   solved as an LP (or QP), as a problem without integer Variable
   *   would be. No cut callback is installed and no cut is separated,
   *   whatever intCutSepPar is: for the relaxation of a Block with
   *   dynamic Constraint, see intRelaxIntVars == 2.
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
   *   class governs the callback installation as appropriate).
   *
   * The dynamic Constraint that the Block adds while the callback of the
   * derived class separates are given to the back-end as user cuts or
   * lazy constraints of that solve, and do not become rows of the model
   * as long as the following solves separate lazy constraints: they are
   * kept in a pool that the callback gives back when they are violated
   * [see separate_in_callback()]. */
  intCutSepPar ,
  /// maximum number of LP cut separation passes (when intRelaxIntVars == 2)
  /**< Caps the number of "solve LP - separate user cuts - re-solve"
   * iterations performed by the base compute() when intRelaxIntVars == 2.
   * Defaults to 1000. Increase if the separator is known to require very
   * many passes to converge; reduce to limit the time spent in the
   * separation loop when fast (possibly weaker) bounds are acceptable. */
  intMaxCutPasses ,
  /// the multipliers of the columns of a dual direction are homogeneous
  /**< What get_dual_direction() writes as the multiplier of a column is, by
   * default, the reduced cost \f$ c - A' y \f$: that is what the multiplier
   * of a dual *point* is, and it is what whoever reads the duals of an
   * optimum needs. A dual *direction*, though, is a ray of the homogeneous
   * system, and the objective has no part in it: its column multipliers are
   * \f$ - A' y \f$. The two differ on every column whose objective
   * coefficient is not zero, which on a real model is most of them, and only
   * the homogeneous ones describe the certificate.
   *
   * - intHomogeneousDirection == 0 [default]: get_dual_direction() writes
   *   \f$ c - A' y \f$, as it has always done.
   *
   * - intHomogeneousDirection == 1: it writes \f$ - A' y \f$.
   *
   * get_dual_solution() is not affected either way. The derived classes that
   * compute the multipliers themselves out of the duals honour the parameter
   * [CPXMILPSolver and GRBMILPSolver]; SCIPMILPSolver asks its back-end for
   * the Farkas coefficients, which are homogeneous to begin with, so there is
   * nothing for it to do; HiGHSMILPSolver, whose back-end only hands out the
   * reduced costs of the iterate it stopped at, refuses the value 1 rather
   * than take it and go on writing multipliers that carry the Objective. */
  intHomogeneousDirection ,
  intLastAlgParMILP  ///< 1st allowed new int parameter for derived classes
  };

 /// enum for double parameters
 enum dbl_par_type_MILP {
  /// tolerance for rounding objective coefficients after a removal delta
  /**< When a Variable is removed from a LinearFunction objective and the
   * Modification carries its coefficient delta, the derived Solver adds
   * that delta to the current total objective coefficient. A finite result
   * is rounded to zero if its absolute value is at most
   *
   *   dblObjRemovalTol * max( 1 , abs( old coefficient ) , abs( delta ) ).
   *
   * The default is 1e-12. The value must be finite and nonnegative; zero
   * disables rounding of nonzero residuals. The tolerance applies only to
   * removals carrying a delta, not to additions or coefficient changes.
   * Removals without a delta retain their existing zeroing behaviour.
   * CPXMILPSolver, GRBMILPSolver, SCIPMILPSolver and HiGHSMILPSolver use it.
   * Changing it affects subsequent removals without reloading the model. */
  dblObjRemovalTol = dblLastParCDAS ,
  /// tolerance of the check of the optimality conditions, 0 = no check
  /**< If positive, whenever the solve returns kOK compute() writes the
   * solution in the Block [see get_var_solution() and get_dual_solution()]
   * and checks it against the Block itself, with the convention of the duals
   * of write_dual_solution(). For the problem written as
   * \f$ \min \{ f(x) : l \leq A x \leq u , \underline{x} \leq x \leq
   * \overline{x} \} \f$ (a maximization is checked as the minimization of
   * \f$ - f \f$), with \f$ \pi = - y \f$ the multipliers of the rows,
   * \f$ y \f$ being their get_dual(), and
   * \f$ z = \nabla f(x) - A^\top \pi \f$ the ones of the bounds, three
   * residuals are measured:
   *
   * - primal: the largest violation of a row or of a bound (and of the
   *   integrality of an integer Variable in a MIP), each divided by
   *   \f$ \max\{ 1 , | \text{side} | \} \f$;
   *
   * - dual: the largest part of a \f$ \pi_i \f$ or \f$ z_j \f$ whose sign
   *   points to an infinite side or bound, divided by
   *   \f$ \max\{ 1 , \| \nabla f(x) \|_\infty \} \f$;
   *
   * - gap: \f$ \sum_i | \pi_i ( A_i x - b_i ) | + \sum_j | z_j ( x_j -
   *   \beta_j ) | \f$, with \f$ b_i \f$ and \f$ \beta_j \f$ the side
   *   and the bound the sign of the multiplier points to, divided by
   *   \f$ \max\{ 1 , | f(x) | \} \f$; for a convex QP this is the
   *   difference between the primal and the dual value.
   *
   * The dual residual and the gap are only measured on a continuous problem
   * (or with intRelaxIntVars != 0) with linear rows, an Objective whose
   * Function is a C05Function and a Solver that has a dual solution; for a
   * MIP, or when the back-end gives no duals, only the primal residual is.
   * If one of them is larger than the value, compute() returns kError, says
   * which on the log, and the solution can still be read; the residuals of
   * the last check are given by get_kkt_residuals(). Since the solution is
   * written in the Block, compute() holds the Block under a write lock()
   * when the check is on. The default is 0. */
  dblKKTTol ,
  /// First allowed new double parameter for derived classes
  dblLastAlgParMILP
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
  /// first allowed new vector-of-string parameter for derived classes
  /// (MILPSolver itself defines none: the path-based sub-Block ignore
  /// list previously here as vstrMILPIgnSBlks has been promoted to the
  /// Solver base interface via Solver::set_excluded_blocks(), which
  /// takes a typed std::unordered_set<Block*> instead of slash-separated
  /// paths; see Solver.h)
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

 /** Scans a static group of FRowConstraint or ColVariable.
 *
 * This function is called from load_problem() for each group of static
 * FRowConstraint or ColVariable of a Block, whatever the shape of the group:
 * the elements are visited in storage order and each run of contiguous ones
 * is recorded, so that the index of an element can be found back from its
 * address. A group of any other type, e.g., of :OneVarConstraint, is left
 * alone.
 *
 * @param group     The group being scanned.
 * @param qb        The block from which the group originated.
 * @param num_block Sequential number of the block in MILPSolver.
 * @param set       Index of the group within qb.
 * @param counter   Current number of elements of type T that have been
 *                  scanned.
 * @tparam T        The element type of the group. Should be either
 *                  FRowConstraint or ColVariable, depending on the group
 *                  being scanned.
 */

 template< typename T >
 void scan_group( const BaseGroup & group , Block * qb , Index num_block ,
                  Index set , Index & counter );

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// scans a dynamic group of FRowConstraint or ColVariable
 /** The dynamic counterpart of scan_group(): the elements of a dynamic group
  * come and go one by one, so each of them gets its own entry of the
  * dictionaries. */

 template< typename T >
 void scan_dynamic_group( const BaseGroup & group , Block * qb ,
                          Index num_block , Index set , Index & counter );

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

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// returns the index of a column in the problem loaded before this one
 /** Each load_problem() but the first one finds, for each column of the
  * problem it loads, the one of the same ColVariable in the problem loaded
  * before it, so that a back-end that loads the problem again whole at a
  * Modification [see reload_on_modification()] can start its solve from
  * the previous one even if columns have been added or removed. The
  * ColVariable are told apart by their address, hence a Variable created
  * where a deleted one was is taken for it: what is carried over is only
  * the start of a solve, which this cannot make wrong.
  *
  * @param j a column index of the problem loaded last
  * @return its index in the problem loaded before, or -1 if it was not
  *         there (or there was no such problem) */

 [[nodiscard]] int previous_column_index( int j ) const {
  return( ( j >= 0 ) && ( std::size_t( j ) < v_prev_col.size() ) ?
          v_prev_col[ j ] : -1 );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// returns the index of a row in the problem loaded before this one
 /** The same as previous_column_index() for the rows, told apart by the
  * address of their FRowConstraint.
  *
  * @param i a row index of the problem loaded last
  * @return its index in the problem loaded before, or -1 if it was not
  *         there (or there was no such problem) */

 [[nodiscard]] int previous_row_index( int i ) const {
  return( ( i >= 0 ) && ( std::size_t( i ) < v_prev_row.size() ) ?
          v_prev_row[ i ] : -1 );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// the number of columns and rows of the problem loaded before this one
 /** The sizes of the problem the indices of previous_column_index() and
  * previous_row_index() refer to, 0 and 0 if there was none. */

 [[nodiscard]] std::pair< int , int > previous_sizes( void ) const {
  return( std::make_pair( prev_numcols , prev_numrows ) );
  }

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
  * A MIP solve that the callback of a derived class interrupts because
  * the incumbent or the bound of the search has reached dblUpCutOff or
  * dblLwCutOff [see Solver.h] returns kCutOff, which is in [ kOK , kError )
  * as kStopTime is, rather than the error status of an interrupted run.
  * Nothing is then proved optimal: as after kStopTime, get_lb() and get_ub()
  * give the value of the incumbent, if any, on its side and the bound of the
  * search on the other, either of them infinite if the solve has none, and
  * get_var_solution() gives the incumbent. The cutoffs are not checked on
  * a continuous problem, nor with intRelaxIntVars == 2.
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
  * to implement get_var_solution().
  *
  * \p solc selects the part of the solution that is written, as the
  * Configuration of get_var_solution() [see Solver.h] does. If it is a
  * SimpleConfiguration< std::vector< int > >, only the ColVariable of the
  * Block that its f_value lists are written: an index k >= 0 stands for the
  * k-th sub-Block of the Block, i.e., get_nested_Block( k ), together with
  * all its sub-Blocks, recursively, while -1 stands for the ColVariable of
  * the Block itself, without those of its sub-Blocks. Any other
  * Configuration, nullptr included, asks for the whole solution. An index
  * that is neither -1 nor that of a sub-Block throws std::invalid_argument.
  * The columns not asked for are skipped a group at a time, so writing a
  * part costs in proportion to the part rather than to the whole Block. */

 void write_var_solution( const std::vector< double > & x ,
                          Configuration * solc = nullptr );

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
  * This method can obviously be used to implement get_dual_solution().
  *
  * \p solc selects the part of the solution that is written as in
  * write_var_solution(): the dual values of the FRowConstraint and the
  * reduced costs of the ColVariable of the Block it lists, a row or column
  * belonging to the Block its Constraint or Variable belongs to. */

 void write_dual_solution( const std::vector< double > & pi ,
			   const std::vector< double > & rc ,
			   Configuration * solc = nullptr );

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// true if the solver handed over the value of the dual direction
 /** The value is the one the solver computes while producing the certificate,
  * and it is only there if the solver both computes it and gives it out: the
  * derived classes that do record it during get_dual_direction() [see
  * f_dual_direction_value]. */

 [[nodiscard]] bool has_dual_direction_value( void ) override {
  return( ! std::isnan( f_dual_direction_value ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// the value of the dual direction produced by the last get_dual_direction()
 /** What the back-end reports for its own certificate, which is CPLEX's
  * dualfarkas proof and Gurobi's FarkasProof; the other two back-ends do not
  * have it [see has_dual_direction_value()].
  *
  * THIS IS NOT THE CONSTANT OF THE CUT THE CERTIFICATE GIVES, and the two are
  * different numbers on any but the smallest models. What gives that constant
  * is walking the Block and accumulating \f$ - \pi_i b_i \f$ over rows and
  * bounds, each taken on the side its multiplier points to (the multipliers
  * reach the Constraint negated, see write_dual_solution()), and that is what
  * a consumer has to do. On the small infeasible instance of the MILPSolver
  * test the two happen to agree, CPLEX and Gurobi reporting 0.6 against
  * multipliers summing to -0.6; on a real one they do not, measured 106000
  * against a certificate worth 6000, the difference being nine of the twelve
  * columns the certificate uses aggregated on the opposite bound. */

 [[nodiscard]] OFValue get_dual_direction_value( void ) override {
  return( f_dual_direction_value );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// the iterations of the back-end in the last compute()
 /** The iterations that the back-end has made in the last compute(), all
  * its solves summed [see elapsed_iter], what an iteration is being the
  * back-end's own (e.g., a change of the working set for an active-set
  * method); 0 for a back-end that does not count them. */

 [[nodiscard]] long get_elapsed_iterations( void ) const override {
  return( elapsed_iter );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// the residuals of the last check of the optimality conditions
 /** The primal residual, the dual residual and the gap of the solution
 * returned by the last compute() [see dblKKTTol], NaN for those that were
 * not measured, all of them if the check is off or the solve did not
 * return kOK. */

 [[nodiscard]] const std::array< double , 3 > & get_kkt_residuals( void )
  const { return( kkt_residuals ); }

/** @} ---------------------------------------------------------------------*/
/*------------------- METHODS FOR HANDLING THE PARAMETERS ------------------*/
/*--------------------------------------------------------------------------*/
/** @name Methods for handling parameters
 * @{ */

 /// sets an integer parameter with the given value
  /** Set the "int" parameters specific of MILPSolver, together with the
  * parameters of MILPSolver that various solver actually "listens to":
  *
  * - intThrowReducedCostException [1]: it indicates whether an exception must
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
  *      or upper bound match the bounds of the Variable.
  *
  *   The default 1 is what get_dflt_int_par() returns, i.e., the value that
  *   a ComputeConfig not giving the parameter sets; a Solver to which no
  *   ComputeConfig has been applied starts with 0. */
 void set_par( idx_type par , int value ) override;

 /// sets a double parameter with the given value
 void set_par( idx_type par , double value ) override;

 /// sets a string parameter with the given value
 void set_par( idx_type par , std::string && value ) override;

 /// sets a vector-of-string parameter with the given value
 void set_par( idx_type par , std::vector< std::string > && value ) override;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// gets the number of integer parameters
 [[nodiscard]] idx_type get_num_int_par( void ) const override;

 /// gets the number of double parameters
 [[nodiscard]] idx_type get_num_dbl_par( void ) const override;

 /// gets the number of string parameters
 [[nodiscard]] idx_type get_num_str_par( void ) const override;

 /// gets the number of vector-of-string parameters
 [[nodiscard]] idx_type get_num_vstr_par( void ) const override;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// gets the default value of the specified integer parameter
 [[nodiscard]] int get_dflt_int_par( idx_type par ) const override;

 /// gets the default value of the specified double parameter
 [[nodiscard]] double get_dflt_dbl_par( idx_type par ) const override;

 /// returns the default value of the specified string parameter
 [[nodiscard]] const std::string & get_dflt_str_par( idx_type par )
  const override;

 /// returns the default value of the specified vector-of-string parameter
 [[nodiscard]] const std::vector< std::string > & get_dflt_vstr_par(
  idx_type par ) const override;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// returns the value of the specified integer parameter
 [[nodiscard]] int get_int_par( idx_type par ) const override;

 /// returns the value of the specified double parameter
 [[nodiscard]] double get_dbl_par( idx_type par ) const override;

 /// returns the value of the specified string parameter
 [[nodiscard]] const std::string & get_str_par( idx_type par ) const override;

 /// returns the value of the specified vector-of-string parameter
 [[nodiscard]] const std::vector< std::string > & get_vstr_par(
  idx_type par ) const override;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// returns the index of the int parameter with the specified name
 [[nodiscard]] idx_type int_par_str2idx( const std::string & name )
  const override;

 /// returns the name of the int parameter with the specified index
 [[nodiscard]] const std::string & int_par_idx2str( idx_type idx )
  const override;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// returns the index of the double parameter with the specified name
 [[nodiscard]] idx_type
 dbl_par_str2idx( const std::string & name ) const override;

 /// returns the name of the double parameter with the specified index
 [[nodiscard]] const std::string &
 dbl_par_idx2str( idx_type idx ) const override;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// returns the index of the string parameter with the specified name
 [[nodiscard]] idx_type str_par_str2idx( const std::string & name )
  const override;

 /// Returns the name of the string parameter with the specified index
 [[nodiscard]] const std::string & str_par_idx2str( idx_type idx )
  const override;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// returns the index of the vstr parameter with the specified name
 [[nodiscard]] idx_type vstr_par_str2idx( const std::string & name )
  const override;

 /// returns the name of the vstr parameter with the specified index
 [[nodiscard]] const std::string & vstr_par_idx2str( idx_type idx )
  const override;

/** @} ---------------------------------------------------------------------*/
/*-------------------- PROTECTED FIELDS OF THE CLASS -----------------------*/
/*--------------------------------------------------------------------------*/

 protected:

/*--------------------------------------------------------------------------*/
/*------------------------- PROTECTED METHODS ------------------------------*/
/*--------------------------------------------------------------------------*/

 /// the Block of #v_BFS whose solution solc asks for
 /** Returns false if \p solc asks for the whole solution, otherwise true and
  * the Block of #v_BFS whose ColVariable and FRowConstraint have to be
  * written in \p sel [see write_var_solution()]; \p who is the name of the
  * method in the exception thrown on an invalid index. */

 bool selected_Blocks( Configuration * solc , const char * who ,
                       std::unordered_set< const Block * > & sel ) const;

/*--------------------------------------------------------------------------*/
 /// the weight of the Objective whose Function is \p f
 /** Returns the factor by which the coefficients of \p f, the Function of
  * the FRealObjective of a Block in the tree of the Block of this Solver,
  * enter the objective of the problem, i.e., the product of the weights of
  * the Objective of that Block and of its ancestors up to, and excluding,
  * the Block of this Solver [see Block::get_objective_weight()], 1 for a
  * Function that is not that of the Objective of such a Block; the
  * back-ends multiply by it whatever they take out of a Modification of
  * \p f. */

 double objective_weight( const Function * f ) const;

/*--------------------------------------------------------------------------*/
 /// the weight of the coefficient of \p var in \p con
 /** Returns the factor by which the coefficient of the Variable \p var in
  * the Constraint \p con, both in the tree of the Block of this Solver,
  * enters the problem, i.e., the product of the weights of the Variable
  * [see Block::get_variable_weight( const Block * )] from the Block of
  * \p var up to the closest Block that also contains the Block of \p con,
  * the latter excluded; it is 1, with no computation, when no Block of the
  * tree has a weight of the Variable other than 1. The back-ends multiply
  * by it whatever coefficient of a row they take out of a Modification. */

 double variable_weight( const Variable * var , const Constraint * con ) const;

/*--------------------------------------------------------------------------*/
 /// MILPSolver handles a change of the weights of a sub-Block
 /** A BlockWeightMod coming from a sub-Block [see Block::set_weight()] has
  * the problem loaded again, with the new weights, as a NBModification. */

 [[nodiscard]] bool handles_weight_change( void ) const override {
  return( true );
  }

/*--------------------------------------------------------------------------*/
 /// MILPSolver applies the weights of the Variable [see variable_weight()]

 [[nodiscard]] bool handles_variable_weight( void ) const override {
  return( true );
  }

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
 /// scatter a LinearFunction into a CSR row of (index, value) buffers
 /** Iterate over the (variable, coefficient) pairs of \p lf and push,
  * into the back of \p out_ind / \p out_val, the (back-end column
  * index, coefficient) of every term whose variable is currently
  * registered with the underlying model. The mapping from
  * ColVariable * to back-end column index is provided by the callable
  * \p idx_fn (typically the derived solver's `*_index_of_variable`
  * member, wrapped in a lambda). Terms pointing at variables that are
  * not registered (\p idx_fn returns Inf<int>()) are silently skipped:
  * the back-end is then free to assume the row touches only declared
  * columns, which mirrors the per-row handling that the legacy
  * implementations of add_dynamic_constraint() were performing inline.
  * \p out_ind and \p out_val are not cleared, so the helper can be
  * called repeatedly to concatenate the rows of a batch into a single
  * pair of buffers (in conjunction with a separate rmatbeg vector).
  *
  * The coefficients of a row are multiplied by their weight [see
  * variable_weight()], hence the helper is a member. It is inline so that
  * the resulting code is identical to the hand-rolled loop the derived
  * solvers were carrying around, and a template so the
  * LinearFunction-iteration boilerplate is stated once and the per-back-end
  * idx_fn stays a zero-overhead abstraction. */

 template< typename IdxFn >
 inline void scatter_lf_to_csr( const LinearFunction * lf , IdxFn && idx_fn ,
                                std::vector< int > & out_ind ,
                                std::vector< double > & out_val ) const {
  const auto con = dynamic_cast< const Constraint * >( lf->get_Observer() );
  for( auto & el : lf->get_v_var() )
   if( auto idx = idx_fn( el.first ) ; idx < Inf< int >() ) {
    out_ind.push_back( idx );
    out_val.push_back( con ? el.second * variable_weight( el.first , con )
                           : el.second );
    }
  }

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
  *  - svar_to_idx, scon_to_idx : vectors of tuples that store 1) the
  *    address of the first element of each run of static variables and
  *    constraints, respectively, 2) the corresponding index in constraint
  *    matrix (column or row) and 3) the number of elements in the run. A run
  *    is a set of elements of a group that are contiguous in memory: a group
  *    that is a single array is one run, a group made of many arrays (e.g.,
  *    a std::vector of std::vector) has one run per array.
  *    The vectors are kept sorted in ascending order by address.
  *
  *  - dvar_to_idx, dcon_to_idx : vectors of pairs that store the addresses
  *    of all the dynamic variables and constraints respectively and the
  *    corresponding index in constraint matrix. and are stored in an
  *    ascending. The vectors are kept sorted in ascending order by address.
  *
  *  - idx_to_svar, idx_to_scon : vectors of pairs that store the indices
  *    of columns and rows, respectively, of the constraint matrix and the
  *    address of the corresponding (run of) static variables and
  *    constraints. The vectors are kept sorted in ascending order by index.
  *
  *  - idx_to_dvar, idx_to_dcon : vectors that store the addresses of dynamic
  *    variables and constraints, respectively. The element at idx_to_d*[ i ]
  *    has index ( i + static_*s ), that is, it is the column/row
  *    ( i + static_*s ) of the constraint matrix.
  * 
  *  - svar_to_bound : a OneVarConstraint * associated to each static
  *    variable (nullptr if it has none), in the order of the columns, so that
  *    the reduced costs are written without looking for the bounds among the
  *    active stuff of each variable [see get_column_bounds()]. If the option
  *    intSingleBound is set to 1, each variable has at most one.
  *
  *  - dvar_to_bound : the same for the dynamic variables, in the order of
  *    idx_to_dvar.
  *
  *  - m_more_bounds : the other OneVarConstraint * of the few variables that
  *    have more than one, by variable rather than by index, so that it stays
  *    valid when a dynamic variable is removed and the indices shift. With
  *    one pointer per column the memory of a problem with many columns is
  *    the same as with intSingleBound.
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

 std::vector< OneVarConstraint * > svar_to_bound;
                                         ///< from static variable to bound
 std::vector< OneVarConstraint * > dvar_to_bound;
                                         ///< from dynamic variable to bound

 std::vector< int > v_prev_col;  ///< previous index of each column, or -1
 std::vector< int > v_prev_row;  ///< previous index of each row, or -1
 int prev_numcols = 0;           ///< columns of the previous problem
 int prev_numrows = 0;           ///< rows of the previous problem

 /// the bounds of a variable beyond the one in svar_to_bound/dvar_to_bound
 std::unordered_map< const ColVariable * ,
                     std::vector< OneVarConstraint * > > m_more_bounds;

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
  * map< mat_indices , double > stored in the i-th position of the vector q_part.
  * 
  * The same procedure is applied for the objective function, with the linear
  * coefficients stored in objective, the diagonal coefficients of the quadratic
  * matrix stored in q_objective and the off-diagonal ones stored using three 
  * vectors ndq_rowind, ndq_colind, ndq_objective.
  * @{  */

 /** Pointers to the currently registered Block and all its descendants
  * arranged along a BFS order. Sub-Blocks listed in
  * Solver::get_excluded_blocks() are excluded from #v_BFS during
  * load_problem() (along with their entire subtree, by virtue of the
  * BFS naturally pruning at the excluded root): they are invisible to
  * every method that traverses the attached Block through #v_BFS. */
 std::vector< Block * > v_BFS;

 /** The same Block pointers as #v_BFS, stored for constant-time membership
  * queries while constructing and updating the MILP matrix. */
 std::unordered_set< Block * > BFS_set;

 /** The weight with which the Objective of each Block of #v_BFS enters the
  * problem [see objective_weight()], for those whose weight is not 1: it is
  * computed going down the tree in load_problem(), each Block having that
  * of its father times its own, and the Block of this Solver 1. */
 std::unordered_map< const Block * , double > f_obj_weight;

 /// true if a Block of #v_BFS below the Block of this Solver has a weight
 /// of the Variable other than 1 [see variable_weight()]
 bool f_var_weighted = false;

 /// the weight of the Variable of a Block in a Block [see
 /// variable_weight()], by pair of Block, as it is asked
 mutable std::map< std::pair< const Block * , const Block * > , double >
  f_var_weight_cache;

 /// buffers of the scan of the columns, kept for the whole load
 /** Each column is scanned with these rather than with containers of its
  * own [see scan_variable()]: the rows it is active in, ( row , position )
  * pairs sorted by row, the ( row , coefficient ) entries of the column and
  * the bounds it has. */
 std::vector< FRowConstraint * > f_scan_rows;
 std::vector< std::pair< int , Index > > f_scan_pos;
 std::vector< std::pair< int , double > > f_scan_sorted;
 std::vector< OneVarConstraint * > f_scan_bounds;

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
 std::vector< std::map< mat_indices , double > > q_part;

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

 /// objective coefficient removal tolerance (see dblObjRemovalTol)
 double obj_removal_tol = 1e-12;

 /// tolerance of the check of the optimality conditions (see dblKKTTol)
 double kkt_tol = 0;

 /// the residuals of the last check [see get_kkt_residuals()]
 std::array< double , 3 > kkt_residuals = {
  std::numeric_limits< double >::quiet_NaN() ,
  std::numeric_limits< double >::quiet_NaN() ,
  std::numeric_limits< double >::quiet_NaN() };

 /// the rows of the pool [see separate_in_callback()]
 std::vector< const FRowConstraint * > v_cut_pool;

 /// the set of Variable each row of v_cut_pool is written on
 std::vector< Index > v_cut_support;

 /// the position of each row of the pool in v_cut_pool
 std::unordered_map< const FRowConstraint * , Index > m_cut_pool;

 /// the index of each set of Variable some row of the pool is written on
 std::map< std::vector< const Variable * > , Index > m_cut_supports;

 /// SMS++-semantic value of intLogVerb captured at set_par() time
 /**< Mirrors the value most recently set via set_par( intLogVerb , v ).
  * The derived classes each map intLogVerb to their backend log
  * parameter (e.g., LogToConsole for Gurobi, CPXPARAM_ScreenOutput for
  * CPLEX), which may clamp or otherwise alter the value; this member
  * preserves the original integer the caller asked for, so that base
  * methods (in particular compute() with its intRelaxIntVars == 2 loop)
  * can take it as the verbosity intended for SMS++-side logging. */
 int log_verbosity = 0;

 /* if true, no more than one OneVarConstraint can be associated to a
 *  single variable. */
 bool single_bound = false;

 /* if true, the model should be cleared and re-loaded. This parameter
  * can be useful when some Modifications are issued but they are not
  * handled by the specialized *MILPSolver. In this case, the model
  * should be reconstructed from scratch. */
 bool f_reset = false;

  /* if true, modification on constraints are enabled. This parameter can
  *  be useful when dealing with quadratic constraint, where Modification 
  * from some Solver (e.g. CPLEX) are not allowed. */
 bool cons_modification = true;

 /** This variable indicates whether an exception must be thrown if there is
 * an inconsistency when a reduced cost is being stored during a call to
 * get_dual_solution() or get_dual_direction(). */
 bool throw_reduced_cost_exception;

 /** true if get_dual_direction() writes the homogeneous multipliers of the
  * columns, - A' y, instead of the reduced costs c - A' y
  * [see intHomogeneousDirection]. */
 bool homogeneous_direction = false;

 /** The value the dual objective takes along the unbounded dual direction
  * the last call to get_dual_direction() wrote in the Block, in the sign
  * convention of the dual solution this Solver produces. A derived class
  * whose back-end computes it records it there; NaN means that no direction
  * has been produced since the last compute(), or that the back-end does not
  * hand the value over, which is what has_dual_direction_value() reports. */
 OFValue f_dual_direction_value = std::numeric_limits< OFValue >::quiet_NaN();

 /** An array of length at least numcols containing pointers to character
  * strings containing the names of the variables. */
 std::vector< char * > colname;

 /** An array of length at least numrows containing pointers to character
  * strings containing the names of the constraints. */
 std::vector< char * > rowname;

 int sol_status = kUnEval; ///< Solution status (OK, Infeasible, Unbounded, ...)

 /// the iterations of the back-end in the current compute()
 /** Cleared by compute() before the first solve, and increased by each
  * guts_of_compute() of a derived class that counts the iterations of its
  * back-end [see get_elapsed_iterations()]. */
 long elapsed_iter = 0;

 /// whether the callback of the cutoffs has interrupted the current solve
 /** Set by the callback of a derived class (possibly from several threads
  * of the back-end) when the incumbent or the bound of the search reaches
  * dblUpCutOff or dblLwCutOff and it asks the back-end to stop; compute()
  * clears it before the solve and reads it after, giving kCutOff rather
  * than the error status of an interrupted run. */
 std::atomic< bool > f_cutoff_stop = false;

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
 *    structure within the file will be used. Without an explicit index
 *    the loader picks the first structure in the file.
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
 /// checks the optimality conditions of the solution [see dblKKTTol]
 /** Writes the solution in the Block, measures the residuals in
  * #kkt_residuals and returns the largest of them. */

 double check_kkt( void );

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

 virtual std::array< double , 2 > get_problem_bounds(
					     const ColVariable & var ) const;

 /// restricts the bounds \p bd of the column of \p var to its value if fixed
 static void fix_bounds( const ColVariable & var ,
			 std::array< double , 2 > & bd ) {
  if( var.is_fixed() ) {
   bd[ 0 ] = std::max( bd[ 0 ] , double( var.get_value() ) );
   bd[ 1 ] = std::min( bd[ 1 ] , double( var.get_value() ) );
   }
  }

/** @} ---------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
 /// returns true if b belongs to the Block subtree loaded in this Solver

 bool is_mine( Block * b ) const {
  return( b && ( BFS_set.count( b ) > 0 ) );
  }

/*--------------------------------------------------------------------------*/
 /// gets the active constraints for the specified variable
 /** The constraints are written in \p active_constraints, which is cleared
  * first, so that a caller scanning many variables keeps one buffer. */
 void get_active_constraints( const ColVariable & var ,
               std::vector< FRowConstraint * > & active_constraints ) const;

 /// the dynamic_cast of \p i to T *, made once per dynamic type
 /** Gives what dynamic_cast< T * >( \p i ) gives. The active stuff of a
  * ColVariable is taken to an FRowConstraint or to a OneVarConstraint for
  * each ColVariable that is loaded, and a dynamic_cast in a hierarchy with
  * multiple inheritance spread over shared libraries compares the names of
  * the types, which a large problem pays tens of millions of times; what
  * the cast gives only depends on the dynamic type of the object, so it is
  * made the first time a type is seen, and the offset from the active
  * stuff to the object of the cast (the same for all the objects of a
  * type, which have the same layout) is kept in a table of the thread,
  * keyed by the std::type_info of the type. */
 template< class T >
 static T * cast_active( ThinVarDepInterface * i ) {
  // PTRDIFF_MIN if the cast fails
  static thread_local std::unordered_map< const std::type_info * ,
					  std::ptrdiff_t > cache;
  auto [ it , fresh ] = cache.try_emplace( & typeid( * i ) , PTRDIFF_MIN );
  if( fresh )
   if( auto * t = dynamic_cast< T * >( i ) )
    it->second = reinterpret_cast< char * >( t ) -
                 reinterpret_cast< char * >( i );
  if( it->second == PTRDIFF_MIN )
   return( nullptr );
  return( reinterpret_cast< T * >( reinterpret_cast< char * >( i ) +
				   it->second ) );
  }

 /// gets the active bounds for the specified variable
 /** As get_active_constraints(), the bounds go in \p active_bounds. They are
  * found among the active stuff of the variable, which is what
  * svar_to_bound and dvar_to_bound keep track of. */
 void get_active_bounds( const ColVariable & var ,
                       std::vector< OneVarConstraint * > & active_bounds ) const;

 /// gets the bounds of a column, from svar_to_bound and dvar_to_bound
 /** As get_active_bounds(), but in constant time: the bounds of the
  * column \p col, whose ColVariable is \p var, go in \p bounds. */
 void get_column_bounds( int col , const ColVariable & var ,
                         std::vector< OneVarConstraint * > & bounds ) const;

 /// records the bounds of a column found by get_active_bounds()
 /** Writes the first of \p bounds in \p slot and the others in
  * m_more_bounds, for the ColVariable \p var. */
 void set_column_bounds( OneVarConstraint * & slot , const ColVariable & var ,
                         const std::vector< OneVarConstraint * > & bounds );

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

 /// whether a pending Modification makes the problem be loaded again whole
 /** If this returns true, process_modifications() does not change the
  * problem of the back-end Modification by Modification: as soon as one is
  * pending, the whole problem is loaded again from the Block, as for an
  * NBModification. This is for a back-end that is not reliable on a problem
  * changed in place with the algorithm it has been asked to use. */

 virtual bool reload_on_modification( void ) { return( false ); }

 /// process all not-GroupModification (bulk of the work)
 void guts_of_process_modifications( const p_Mod mod );

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// executes a GroupModification whole, if the back-end knows how to
 /** A GroupModification is a bunch of Modification that the :Block has
  * declared to belong together, and some of these bunches are one single
  * operation of the underlying solver rather than many: the obvious case is
  * a whole column, which the "abstract" representation describes as one new
  * Variable plus one coefficient change for each row the Variable appears
  * in [see VariableGroupMod in Modification.h]. This method is asked whether
  * \p gmod is one of those before the group is taken apart: if it answers
  * true the group has been dealt with and nothing else is done with it, if
  * it answers false the sub-Modification are dispatched one by one exactly
  * as they are when they arrive alone.
  *
  * The implementation here recognises the shapes that the core describes and
  * hands each of them to a method of its own [see add_columns(),
  * remove_columns(), change_coefficients() and change_bounds()], which is
  * what a :MILPSolver overrides to turn the operation into one call of its
  * own API. Since those answer false unless
  * overridden, a back-end that does nothing keeps exactly the behaviour it
  * has today.
  *
  * @param gmod the GroupModification to be executed whole
  * @return true if the group has been executed, false to take it apart */

 virtual bool process_group_modification( const GroupModification * gmod );

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// adds the whole column of each of the given Variable in one operation
 /** Called when a VariableGroupMod says that \p vars have been added to the
  * Block, with all the coefficient changes that the addition has caused
  * inside the same group [see process_group_modification()]. The Variable
  * are in the Block already, hence each of them knows the Function it is
  * active in and each Function knows its coefficient: the column is
  * therefore available without reading the sub-Modification at all.
  *
  * The implementation here does nothing and answers false, which is how a
  * back-end says that the column has to be built one Modification at a time
  * as usual.
  *
  * @param vars the Variable whose columns are being added
  * @return true if the columns have been added */

 virtual bool add_columns( const std::vector< Variable * > & vars ,
                           const GroupModification * gmod ) {
  return( false );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// reads the entries of a group that declares a column
 /** Reads what a group of the column shape holds, without touching the
  * model, so that a group which turns out not to be a plain column can
  * still be handed back to the one-by-one path: the entries of the new
  * columns come back as the three parallel vectors \p econs (the row an
  * entry belongs to), \p evars (which of \p vars the entry belongs to) and
  * \p evals (its value), and the cost of each of \p vars, summed over the
  * Objective it appears in, comes back in \p cost.
  *
  * Returns false, leaving the four untouched, if the group is not a plain
  * addition of whole columns: something other than a linear Function, a
  * removal inside an addition, or an entry of a Variable that was there
  * already, whose cost would have to be read back from the model.
  *
  * The entries are replayed from the group rather than read from the Block,
  * which by the time a Modification is processed has moved on: this is what
  * the handlers do one at a time [see constraint_fvars_modification() and
  * objective_fvars_modification()]. */

 bool read_column_group( const std::vector< Variable * > & vars ,
                         const GroupModification * gmod ,
                         std::vector< const FRowConstraint * > & econs ,
                         std::vector< Index > & evars ,
                         std::vector< double > & evals ,
                         std::vector< double > & cost );

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// removes the whole column of each of the given Variable in one operation
 /** The counterpart of add_columns(): \p vars are being removed from the
  * Block, and with them every coefficient they have in any row, which is
  * what deleting the column does by itself. Note that the Variable are
  * still in the data structures of the :MILPSolver when this is called, and
  * they are still those of the model, so their indices are the right ones.
  *
  * The implementation here does nothing and answers false, which is how a
  * back-end says that the column has to be emptied one Modification at a
  * time as usual.
  *
  * @param vars the Variable whose columns are being removed
  * @return true if the columns have been removed */

 virtual bool remove_columns( const std::vector< Variable * > & vars ,
                              const GroupModification * gmod ) {
  return( false );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// changes the linear coefficients of \p mods in one operation
 /** Called when a group contains nothing but changes of the linear
  * coefficients of some Function, which is what a :Block issues when one
  * datum of it enters many rows, or many entries of the same row: each of
  * them is a Modification of its own, and the model has to be told all of
  * them, but it can be told once [see process_group_modification()].
  *
  * The implementation here does nothing and answers false, which is how a
  * back-end says that the changes have to be made one at a time as usual.
  *
  * @param mods the changes, in the order they were issued
  * @return true if the changes have been made */

 virtual bool change_coefficients(
                     const std::vector< const FunctionMod * > & mods ) {
  return( false );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// changes the bounds of \p mods in one operation
 /** The counterpart of change_coefficients() for the bounds: a group of
  * OneVarConstraintMod, which is what a :Block issues when it fixes or frees
  * a whole set of Variable at once.
  *
  * The implementation here does nothing and answers false, which is how a
  * back-end says that the bounds have to be written one at a time as usual.
  *
  * @param mods the changes, in the order they were issued
  * @return true if the bounds have been written */

 virtual bool change_bounds(
              const std::vector< const OneVarConstraintMod * > & mods ) {
  return( false );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// changes the sides of the rows of \p mods in one operation
 /** The counterpart of change_bounds() for the rows: a group of
  * RowConstraintMod, which is what a :Block issues when one datum of it is
  * the right-hand side of a whole vector of Constraint, one per time
  * instant.
  *
  * The implementation here does nothing and answers false, which is how a
  * back-end says that the sides have to be written one at a time as usual.
  *
  * @param mods the changes, in the order they were issued
  * @return true if the sides have been written */

 /// changes the state of a whole set of Variable in one operation
 /** Executes a batch of VariableMod, i.e., of changes of the state of a
  * Variable, all of which the Solver has collected out of one
  * GroupModification [see process_group_modification()]; returns true if
  * it has done so, false if the batch has to be executed one Modification
  * at a time, which is what the implementation in the base class does.
  *
  * What a back-end can do in one operation here is the fixing and the
  * unfixing, which is a change of the bounds of the columns; a batch
  * carrying a change of integrality is refused, that one being a call per
  * column anyway. A column can be named several times in the same batch,
  * e.g., a Variable fixed and unfixed again, and what it has to end up with
  * is the state of the last of them, exactly as when the same changes are
  * executed one at a time.
  *
  * A batch of these and a batch of changes of the bounds [see
  * change_bounds()] are kept apart and executed in the order they were
  * issued, both being about the same attribute of the model. */

 virtual bool change_variables(
                        const std::vector< const VariableMod * > & mods ) {
  return( false );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 virtual bool change_sides(
              const std::vector< const RowConstraintMod * > & mods ) {
  return( false );
  }

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

 /// true if \p mod is, or has inside, a Modification of an ignored Objective
 /** Returns true if \p mod is an ObjectiveMod, a FunctionMod or a
  * FunctionModVars of the Objective of a Block given to
  * Solver::set_excluded_objectives(), or a GroupModification having one of
  * them inside at any depth. */
 [[nodiscard]] bool touches_excluded_objective( const Modification * mod ) const;

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

 /// handles a change of the weights of a Block [see BlockWeightMod]
 /** The Objective of the problem, and its constant, are built again when
  * the weight of the Objective has changed [see set_linear_objective()], and
  * the coefficients of the Variable of the Block and of those below it in
  * the rows outside of it when the weight of the Variable has [see
  * set_matrix_coefficients()]; f_reset is set, so that the problem is loaded
  * again whole, if a Block below the changed one has a quadratic Objective,
  * if a row whose coefficients change is not linear, or if the back-end
  * does not change them in place. */
 void weight_modification( const BlockWeightMod * mod );

 /// gives the back-end the linear coefficients \p obj of the Objective
 /** \p obj has one coefficient per column, in the order of the columns,
  * the weights comprised. The method of the base class sets f_reset, so
  * that the problem is loaded again: a back-end that can change the
  * coefficients of its Objective in place redefines it. */
 virtual void set_linear_objective( const std::vector< double > & obj ) {
  f_reset = true;
  }

 /// gives the back-end new coefficients of the matrix
 /** The coefficient of the column cols[ i ] in the linear row rows[ i ]
  * becomes vals[ i ], the weights comprised. The method of the base class
  * sets f_reset, so that the problem is loaded again: a back-end that can
  * change the coefficients in place redefines it. */
 virtual void set_matrix_coefficients(
	    const std::vector< const FRowConstraint * > & rows ,
	    const std::vector< int > & cols ,
	    const std::vector< double > & vals ) {
  f_reset = true;
  }

 /// computes the weights of the Block of v_BFS
 /** Fills f_obj_weight and f_var_weighted out of the weights of the Block
  * of v_BFS [see objective_weight() and variable_weight()]. */
 void refresh_weights( void );

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// calls \p f on every addition of FRowConstraint inside \p mod
 /** The rows that a Block generates may well arrive inside a
  * GroupModification, which a scan of the first level alone does not see,
  * and the rows would then silently disappear: this walks the group,
  * however deep it is, and calls \p f on each BlockModAdd< FRowConstraint >
  * in the order they were issued [see perform_separation()]. */

 static void for_each_row_addition( const Modification * mod ,
  const std::function< void( const BlockModAdd< FRowConstraint > * ) > & f );

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// the rows that a callback gives the back-end at the current point
 /** To be called from within the callback of the back-end, with the Block
  * lock()-ed and the point (be it a relaxed or a feasible one) written in
  * its Variable. It calls generate_dynamic_constraints( cfg ) on the Block
  * and writes in \p rows, which is supposed to be empty, the rows that the
  * callback has to give the back-end at this point: (i) all the dynamic
  * FRowConstraint that the Block has added in the call, and (ii) for each
  * set of Variable on which none of (i) is written, the most violated row
  * of the pool written on that set, if any is violated.
  *
  * The pool is the set of the rows that the Block has added during the
  * callbacks of this Solver: they are dynamic Constraint of the Block, and
  * this Solver knows them, but they are not rows of the model. A row that
  * a callback gives a back-end as a user cut or a lazy constraint belongs
  * to that solve, and a back-end is free to discard it and to show the
  * same point again later (another thread, a cut not deemed useful): a
  * Block that does not give again a Constraint it already has is then
  * served by (ii). When the BlockModAdd of a row of the pool is processed,
  * the row is not added to the model; when its BlockModRmv is, the row
  * leaves the pool [see dynamic_modification()].
  *
  * Hence, with lazy constraint separation the model keeps the rows of the
  * formulation alone, and those separated are given back by the callback
  * whenever they are violated, solve after solve. If the next solve does
  * not separate lazy constraints [see cut_pool_enforced()], before it the
  * pool becomes rows of the model, as the BlockModAdd would have made it;
  * so does load_problem(), which reads all the dynamic Constraint of the
  * Block. */

 void separate_in_callback( Configuration * cfg ,
                            std::vector< const FRowConstraint * > & rows );

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// whether the callback of the next solve gives back the pool
 /** True if the next solve will run the callback that separates lazy
  * constraints, which is what makes it safe to keep the pool out of the
  * model: the default says so for the MIP path, i.e., with integer columns
  * and intRelaxIntVars == 0, when intCutSepPar has bit 2 set; a back-end
  * that cannot give lazy constraints from its callback answers false. */

 virtual bool cut_pool_enforced( void ) const {
  return( ( CutSepPar & 4 ) && ( int_vars > 0 ) && ( ! relax_int_vars ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// adds all the rows of the pool to the model, and empties it

 void flush_cut_pool( void );

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// empties the pool without adding its rows to the model

 void clear_cut_pool( void );

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// removes \p con from the pool, returning false if it is not there

 bool erase_from_cut_pool( const FRowConstraint * con );

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

 /// batch-adds a sequence of new dynamic constraints
 /** Default implementation simply loops over \p cons calling
  * add_dynamic_constraint() once per row, which keeps every derived
  * :MILPSolver functionally correct. Derived classes whose back-end
  * exposes a multi-row insertion entry point (e.g. CPLEX
  * CPXaddrows( ..., nrows, ... ), Gurobi GRBaddconstrs, HiGHS
  * Highs_addRows) override this method to push the whole batch in a
  * single API call, which avoids the per-row overhead of repeatedly
  * touching the back-end's internal data structures when a
  * BlockModAdd<FRowConstraint> arrives with many rows at once (the
  * typical case for a PolyhedralFunctionModAddd retraduced by the
  * owning PolyhedralFunctionBlock into a single batched
  * add_dynamic_constraints(f_const, newc, ...)). The contract is the
  * same as add_dynamic_constraint() applied row-by-row: the order in
  * \p cons is preserved, every row gets the next free slot in numrows,
  * and the dictionaries (dcon_to_idx / idx_to_dcon) are updated
  * accordingly. */

 virtual void add_dynamic_constraints(
                       const std::vector< const FRowConstraint * > & cons ) {
  for( auto * con : cons )
   add_dynamic_constraint( con );
  }

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
 * parameters. That's because they are meant to be used while walking the
 * groups of the Block, and the counters keep track of the elements inside
 * them.
 * @{ */

 /** Scans a dynamic ColVariable and fills the dictionaries accordingly
  *
  * @param var a reference to a ColVariable
  * @param col a counter for variables/columns */

 void scan_dynamic_variable( const ColVariable & var , Index & col );

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// fills the dictionaries and the matrix with one ColVariable

 void scan_variable( const ColVariable & var , Index & col );

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
 /// fills the dictionaries and the matrix with one FRowConstraint

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
  * All the coefficients and the constant term are multiplied by \p w, the
  * weight of the Objective [see objective_weight()].
  *
  * @param obj a FRealObjective
  * @param w the weight of its Block */

 void scan_objective( const FRealObjective * obj , double w );

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
