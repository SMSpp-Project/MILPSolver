/*--------------------------------------------------------------------------*/
/*-------------------------- File OSQPMILPSolver.h -------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the OSQPMILPSolver class, implementing the MILPSolver
 * interface for convex Quadratic Programs with linear constraints using the
 * operator splitting solver OSQP.
 *
 * \author Donato Meoli \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \author Antonio Frangioni \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \copyright &copy; by Donato Meoli, Antonio Frangioni
 */
/*--------------------------------------------------------------------------*/
/*----------------------------- DEFINITIONS --------------------------------*/
/*--------------------------------------------------------------------------*/

#ifndef __OSQPMILPSOLVER_H
 #define __OSQPMILPSOLVER_H
                      /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include <string>
#include <vector>

#include <osqp.h>

#include "MILPSolver.h"

/*--------------------------------------------------------------------------*/
/*--------------------------- NAMESPACE ------------------------------------*/
/*--------------------------------------------------------------------------*/

/// namespace for the Structured Modeling System++ (SMS++)
namespace SMSpp_di_unipi_it
{

/*--------------------------------------------------------------------------*/
/*------------------------ CLASS OSQPMILPSolver ----------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------- GENERAL NOTES --------------------------------*/
/*--------------------------------------------------------------------------*/
/// class for solving convex QP Blocks via the solver OSQP
/** The OSQPMILPSolver implements the MILPSolver interface using OSQP, the
 * open-source solver for convex Quadratic Programs based on the alternating
 * direction method of multipliers (ADMM)
 *
 * \f[
 *   \min \{ \tfrac{1}{2} x^\top P x + q^\top x \,:\,
 *           l \leq A x \leq u \} \,,
 * \f]
 *
 * where \f$ P \f$ only needs to be positive semidefinite (a zero \f$ P \f$,
 * i.e., that of a Linear Program, is given to OSQP as such). It is
 * registered in the Solver factory, hence it can be selected by name
 * ("OSQPMILPSolver") in a BlockSolverConfig.
 *
 * OSQP has no bounds on the variables: the rows of \f$ A \f$ are those of
 * the problem with at least one finite side (as OSQP sees it, i.e., below
 * OSQP_INFTY in absolute value), followed by a unit row for each variable
 * with at least one finite bound, whose multiplier is the reduced cost of
 * the variable. The equalities are the rows with \f$ l_i = u_i \f$.
 *
 * Whenever a Modification is received the problem is loaded again whole
 * from the Block [see MILPSolver::reload_on_modification()]; if the sizes
 * and the sparsity patterns of \f$ P \f$ and \f$ A \f$ have not changed,
 * the new data are given to the solver of OSQP already set up, which is
 * set up again otherwise. After the first solve, each one starts from the
 * primal and dual solution of the previous one, as long as the size of the
 * problem has not changed (see intOSQPWarmStart). After the ADMM has
 * converged, the solution is polished (see intOSQPPolishing), i.e., the
 * problem restricted to the active constraints of the ADMM solution is
 * solved, which gives the solution and the multipliers to high accuracy
 * if the active set has been guessed right, and leaves the ADMM ones
 * otherwise.
 *
 * Unlike active-set and interior point solvers, OSQP detects the dual
 * infeasibility of a problem also when \f$ P \f$ is singular, e.g., an
 * unbounded Linear Program, hence such a problem ends with kUnbounded.
 *
 * The following restrictions apply:
 *
 * - quadratic constraints are not supported, and integer variables are only
 *   accepted if relaxed via intRelaxIntVars;
 *
 * - OSQP reaches the tolerances of the ADMM in many iterations of low
 *   accuracy; if intMaxIter or dblMaxTime stop it before, the iterate is
 *   given as an approximate primal and dual solution if it is within ten
 *   times the tolerances (SOLVED_INACCURATE in OSQP), but compute()
 *   returns kStopIter or kStopTime and both get_lb() and get_ub() are
 *   infinite, the solution not being certified;
 *
 * - dblMaxTime is only handled if OSQP has been built with its profiling
 *   (OSQP_ENABLE_PROFILING, the default), and includes the time of the
 *   set up of OSQP after a change of the sparsity patterns;
 *
 * - without polishing, or when it fails, the solution and the multipliers
 *   are accurate to the tolerances of the ADMM only (dblFAccSol and
 *   dblOSQPEpsRel), and the value of the solution is computed from it;
 *
 * - no unbounded ray or infeasibility certificate is provided, i.e., the
 *   has_var_direction() / has_dual_direction() interface of CDASolver is
 *   not implemented. */

class OSQPMILPSolver : public MILPSolver
{
/*--------------------------------------------------------------------------*/
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

 public:

/*--------------------------------------------------------------------------*/
/*---------------------------- PUBLIC TYPES --------------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Public types of OSQPMILPSolver
 * @{ */

 /// public enum for the integer algorithmic parameters
 /** Public enum describing the different algorithmic parameters of int type
  * that OSQPMILPSolver has in addition to those of MILPSolver. The value
  * intLastAlgParOSQP is provided so that the list can be easily further
  * extended by derived classes. */
 enum int_par_type_OSQP {
  /// whether a solve starts from the solution of the previous one
  /**< If 1 [default], each solve after the first starts from the primal
   * and dual solution of the previous one, as long as the number of
   * variables and of rows of OSQP has not changed; if 0, each solve
   * starts from zero. */
  intOSQPWarmStart = intLastAlgParMILP ,
  /// whether the solution of the ADMM is polished (default 1)
  intOSQPPolishing ,
  /// number of iterative refinement steps of the polishing (default 3)
  intOSQPPolishRefineIter ,
  /// number of iterations of the scaling of the data, 0 for none
  /**< The default is that of OSQP, 10; a change has the solver of OSQP set
   * up again at the next solve. */
  intOSQPScaling ,
  /// the rule of the adaptation of rho
  /**< 0 none, 1 every some iterations (the default of OSQP), 2 on the time
   * of the set up, 3 [default] on the decrease of the error of the
   * optimality conditions, with which the ADMM converges on degenerate
   * Linear Programs where with 1 rho oscillates; a change has the solver
   * of OSQP set up again at the next solve. */
  intOSQPAdaptiveRho ,
  /// number of iterations between two checks of termination, 0 for none
  /**< The default is that of OSQP, 25; with 0, the ADMM only stops at
   * intMaxIter or dblMaxTime. */
  intOSQPCheckTermination ,
  /// first allowed new int parameter for derived classes
  intLastAlgParOSQP
  };

 /// public enum for the double algorithmic parameters
 /** Public enum describing the different algorithmic parameters of double
  * type that OSQPMILPSolver has in addition to those of MILPSolver. The
  * value dblLastAlgParOSQP is provided so that the list can be easily
  * further extended by derived classes. Note that the absolute tolerance
  * on the primal and dual residuals (eps_abs of OSQP) is the standard
  * dblFAccSol. */
 enum dbl_par_type_OSQP {
  /// relative tolerance on the primal and dual residuals (eps_rel)
  dblOSQPEpsRel = dblLastAlgParMILP ,
  /// tolerance of the certificate of primal infeasibility (eps_prim_inf)
  dblOSQPEpsPrimInf ,
  /// tolerance of the certificate of dual infeasibility (eps_dual_inf)
  dblOSQPEpsDualInf ,
  /// initial step size of the ADMM (rho)
  /**< A change has the solver of OSQP set up again at the next solve. */
  dblOSQPRho ,
  /// regularization of the primal variables of the ADMM (sigma)
  /**< A change has the solver of OSQP set up again at the next solve. */
  dblOSQPSigma ,
  /// relaxation parameter of the ADMM, in ( 0 , 2 ) (alpha)
  dblOSQPAlpha ,
  /// regularization of the linear system of the polishing (delta)
  dblOSQPDelta ,
  /// first allowed new double parameter for derived classes
  dblLastAlgParOSQP
  };

/** @} ---------------------------------------------------------------------*/
/*------------------- CONSTRUCTING AND DESTRUCTING -------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Constructing and destructing OSQPMILPSolver
 * @{ */

 /// void constructor
 OSQPMILPSolver( void );

/*--------------------------------------------------------------------------*/
 /// destructor: releases the solver of OSQP

 ~OSQPMILPSolver() override;

/** @} ---------------------------------------------------------------------*/
/*--------------------- DERIVED METHODS OF BASE CLASS ----------------------*/
/*--------------------------------------------------------------------------*/
/** @name Methods of MILPSolver that are redefined
 * @{ */

 /// clears the problem, together with the solver of OSQP

 void clear_problem( unsigned int what ) override;

/*--------------------------------------------------------------------------*/
 /// loads the problem from the Block into the data of OSQP

 void load_problem( void ) override;

/** @} ---------------------------------------------------------------------*/
/*-------------------- METHODS FOR READING RESULTS -------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Accessing the found solutions (if any)
 * @{ */

 /// returns a lower bound on the optimal objective function value

 OFValue get_lb( void ) override;

 /// returns an upper bound on the optimal objective function value

 OFValue get_ub( void ) override;

 /// returns true if a primal solution is available

 bool has_var_solution( void ) override;

 /// returns the objective function value of the primal solution

 OFValue get_var_value( void ) override;

 /// writes the primal solution in the Variable of the Block

 void get_var_solution( Configuration * solc = nullptr ) override;

 /// returns true if a dual solution is available

 bool has_dual_solution( void ) override;

 /// returns true if the dual solution is feasible

 bool is_dual_feasible( void ) override;

 /// writes the dual solution in the Constraint of the Block

 void get_dual_solution( Configuration * solc = nullptr ) override;

/** @} ---------------------------------------------------------------------*/
/*------------------ METHODS FOR HANDLING THE PARAMETERS -------------------*/
/*--------------------------------------------------------------------------*/
/** @name Handling the parameters of OSQPMILPSolver
 *
 * Besides those of MILPSolver and those of OSQP above, the standard
 * intMaxIter (the iteration limit of OSQP, none by default rather than the
 * 4000 of OSQP, which do not suffice to many Linear Programs), dblMaxTime
 * (its time limit) and dblFAccSol (its absolute tolerance on the primal
 * and dual residuals) are handled.
 * @{ */

 void set_par( idx_type par , int value ) override;

 void set_par( idx_type par , double value ) override;

 [[nodiscard]] idx_type get_num_int_par( void ) const override {
  return( MILPSolver::get_num_int_par()
	  + intLastAlgParOSQP - intLastAlgParMILP );
  }

 [[nodiscard]] idx_type get_num_dbl_par( void ) const override {
  return( MILPSolver::get_num_dbl_par()
	  + dblLastAlgParOSQP - dblLastAlgParMILP );
  }

 [[nodiscard]] int get_dflt_int_par( idx_type par ) const override;

 [[nodiscard]] double get_dflt_dbl_par( idx_type par ) const override;

 [[nodiscard]] int get_int_par( idx_type par ) const override;

 [[nodiscard]] double get_dbl_par( idx_type par ) const override;

 [[nodiscard]] idx_type int_par_str2idx( const std::string & name )
  const override;

 [[nodiscard]] const std::string & int_par_idx2str( idx_type idx )
  const override;

 [[nodiscard]] idx_type dbl_par_str2idx( const std::string & name )
  const override;

 [[nodiscard]] const std::string & dbl_par_idx2str( idx_type idx )
  const override;

/** @} ---------------------------------------------------------------------*/
/*--------------------- PROTECTED PART OF THE CLASS ------------------------*/
/*--------------------------------------------------------------------------*/

 protected:

/*--------------------------------------------------------------------------*/
/*-------------------------- PROTECTED METHODS -----------------------------*/
/*--------------------------------------------------------------------------*/

 /// sets up OSQP on the loaded problem, if needed, and solves it

 int guts_of_compute( void ) override;

 /// OSQP is not given the changes of a problem: it is always loaded again

 bool reload_on_modification( void ) override { return( true ); }

/*--------------------------------------------------------------------------*/
/*---------------------------- PRIVATE PART --------------------------------*/
/*--------------------------------------------------------------------------*/

 private:

 /// translates the status of the last solve into that of compute()

 [[nodiscard]] int decode_status( void ) const;

 /// releases the solver of OSQP

 void free_solver( void );

 /// true if the last solve has given an optimal solution

 [[nodiscard]] bool solved( void ) const;

 /// true if the last solve has given a solution, possibly approximate

 [[nodiscard]] bool has_solution( void ) const;

/*--------------------------------------------------------------------------*/
/*--------------------------- PRIVATE FIELDS -------------------------------*/
/*--------------------------------------------------------------------------*/

 // the problem in the format of OSQP: P (its upper triangle) and A in
 // compressed sparse columns, with l <= A x <= u

 OSQPInt d_n = 0;                ///< number of variables
 OSQPInt d_m = 0;                ///< number of rows of A
 double d_sign = 1;              ///< -1 if a maximization was negated

 std::vector< OSQPInt > d_Pp;    ///< column starts of P
 std::vector< OSQPInt > d_Pi;    ///< row indices of P
 std::vector< OSQPFloat > d_Px;  ///< values of P
 std::vector< OSQPFloat > d_q;   ///< linear term
 std::vector< OSQPInt > d_Ap;    ///< column starts of A
 std::vector< OSQPInt > d_Ai;    ///< row indices of A
 std::vector< OSQPFloat > d_Ax;  ///< values of A
 std::vector< OSQPFloat > d_l;   ///< lower sides of the rows of A
 std::vector< OSQPFloat > d_u;   ///< upper sides of the rows of A

 /// position of each row of the problem in A, or -1 if it is left out
 std::vector< OSQPInt > d_row_pos;
 /// position in A of the bounds of each variable, or -1
 std::vector< OSQPInt > d_box_pos;

 OSQPSolver * d_solver = nullptr;  ///< the solver of OSQP, if any
 bool d_loaded = false;            ///< true if the problem has changed
 bool d_setup = false;             ///< true if a set up is needed
 bool d_reloading = false;         ///< true within load_problem()

 // the data the solver of OSQP has been set up with, for the update
 std::vector< OSQPInt > d_sPp;   ///< column starts of P
 std::vector< OSQPInt > d_sPi;   ///< row indices of P
 std::vector< OSQPFloat > d_sPx; ///< values of P
 std::vector< OSQPInt > d_sAp;   ///< column starts of A
 std::vector< OSQPInt > d_sAi;   ///< row indices of A
 std::vector< OSQPFloat > d_sAx; ///< values of A

 std::vector< OSQPFloat > d_x;   ///< primal solution
 std::vector< OSQPFloat > d_y;   ///< multipliers of the rows of A
 double d_fval = 0;              ///< value of the last solution
 int d_status = -1;              ///< status of the last solve
 bool d_iter_limit = false;      ///< true if it has hit intMaxIter

 OSQPSettings d_settings;            ///< the parameters of OSQP
 int d_max_iter = Inf< int >();      ///< intMaxIter
 double d_max_time = Inf< double >(); ///< dblMaxTime
 bool d_warm_start = true;           ///< intOSQPWarmStart

/*--------------------------------------------------------------------------*/

 SMSpp_insert_in_factory_h;

/*--------------------------------------------------------------------------*/

 };  // end( class( OSQPMILPSolver ) )

/*--------------------------------------------------------------------------*/

}  // end( namespace SMSpp_di_unipi_it )

/*--------------------------------------------------------------------------*/

#endif  /* OSQPMILPSolver.h included */

/*--------------------------------------------------------------------------*/
/*---------------------- End File OSQPMILPSolver.h -------------------------*/
/*--------------------------------------------------------------------------*/
