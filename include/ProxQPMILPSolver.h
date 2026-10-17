/*--------------------------------------------------------------------------*/
/*------------------------- File ProxQPMILPSolver.h ------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the ProxQPMILPSolver class, implementing the MILPSolver
 * interface for convex Quadratic Programs with linear constraints using the
 * proximal augmented Lagrangian solver ProxQP of the ProxSuite library.
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

#ifndef __PROXQPMILPSOLVER_H
 #define __PROXQPMILPSOLVER_H
                      /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include <memory>
#include <string>
#include <vector>

#include <Eigen/Core>
#include <Eigen/SparseCore>

#include "MILPSolver.h"

/*--------------------------------------------------------------------------*/
/*----------------------- FORWARD DECLARATIONS -----------------------------*/
/*--------------------------------------------------------------------------*/

// the solvers of ProxQP, whose headers are only included by the .cpp

namespace proxsuite::proxqp::dense { template< typename T > struct QP; }
namespace proxsuite::proxqp::sparse {
 template< typename T , typename I > struct QP; }

/*--------------------------------------------------------------------------*/
/*--------------------------- NAMESPACE ------------------------------------*/
/*--------------------------------------------------------------------------*/

/// namespace for the Structured Modeling System++ (SMS++)
namespace SMSpp_di_unipi_it
{

/*--------------------------------------------------------------------------*/
/*----------------------- CLASS ProxQPMILPSolver ---------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------- GENERAL NOTES --------------------------------*/
/*--------------------------------------------------------------------------*/
/// class for solving convex QP Blocks via the solver ProxQP
/** The ProxQPMILPSolver implements the MILPSolver interface using ProxQP,
 * the open-source proximal augmented Lagrangian solver of the ProxSuite
 * library for convex Quadratic Programs
 *
 * \f[
 *   \min \{ \tfrac{1}{2} x^\top H x + g^\top x \,:\,
 *           A x = b \,,\, l \leq C x \leq u \,,\,
 *           l_x \leq x \leq u_x \} \,,
 * \f]
 *
 * where \f$ H \f$ only needs to be positive semidefinite (a zero \f$ H \f$,
 * i.e., that of a Linear Program, is given to ProxQP as such). It is
 * registered in the Solver factory, hence it can be selected by name
 * ("ProxQPMILPSolver") in a BlockSolverConfig.
 *
 * The rows of the problem whose two sides coincide are the equalities
 * \f$ A x = b \f$, all the others the inequalities \f$ l \leq C x \leq u \f$.
 * Two backends of ProxQP can be used (see intProxQPBackend): the dense one
 * [default], which keeps the bounds on the variables as box constraints of
 * its own and is meant for the small and medium problems solved many times
 * with small changes, and the sparse one, which has no box constraints, so
 * that each variable with a finite bound is a further row of \f$ C \f$.
 * After the first solve, each one starts from the primal and dual solution
 * of the previous one, as long as the size of the problem has not changed
 * (see intProxQPWarmStart).
 *
 * The following restrictions apply:
 *
 * - quadratic constraints are not supported, and integer variables are only
 *   accepted if relaxed via intRelaxIntVars;
 *
 * - ProxQP is not given the changes of the loaded problem: whenever a
 *   Modification is received the problem is loaded again whole from the
 *   Block [see MILPSolver::reload_on_modification()], and if its size has
 *   not changed the solution of the previous solve is the start of the
 *   next one;
 *
 * - ProxQP has no time limit, hence dblMaxTime is not handled and
 *   compute() never returns kStopTime;
 *
 * - the infeasibility criteria of ProxQP do not always succeed: an
 *   unbounded problem, and an infeasible one with the dense backend, may
 *   end with kStopIter at the iteration limit rather than with kUnbounded
 *   or kInfeasible; conversely, too loose thresholds (see
 *   dblProxQPEpsPrimalInf) have a feasible problem declared infeasible;
 *
 * - the solution and the multipliers are accurate to the tolerances of
 *   ProxQP (dblFAccSol and dblProxQPEpsRel), e.g., a multiplier of 1 may
 *   come out as 0.999999999, and the value of the solution is computed
 *   from it;
 *
 * - no unbounded ray or infeasibility certificate is provided, i.e., the
 *   has_var_direction() / has_dual_direction() interface of CDASolver is
 *   not implemented. */

class ProxQPMILPSolver : public MILPSolver
{
/*--------------------------------------------------------------------------*/
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

 public:

/*--------------------------------------------------------------------------*/
/*---------------------------- PUBLIC TYPES --------------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Public types of ProxQPMILPSolver
 * @{ */

 /// public enum for the integer algorithmic parameters
 /** Public enum describing the different algorithmic parameters of int type
  * that ProxQPMILPSolver has in addition to those of MILPSolver. The value
  * intLastAlgParProxQP is provided so that the list can be easily further
  * extended by derived classes. */
 enum int_par_type_ProxQP {
  /// whether a solve starts from the solution of the previous one
  /**< If 1 [default], each solve after the first starts from the primal
   * and dual solution of the previous one, as long as the number of
   * variables, of equalities and of inequalities has not changed; if 0,
   * each solve starts from the initial guess of intProxQPInitialGuess. */
  intProxQPWarmStart = intLastAlgParMILP ,
  /// the start of a solve that is not warm-started
  /**< If 0, the primal and dual solutions start from zero; if 1 [default],
   * the primal solution starts from that of the problem with the
   * equalities alone (the initial_guess EQUALITY_CONSTRAINED_INITIAL_GUESS
   * of ProxQP). */
  intProxQPInitialGuess ,
  /// the backend of ProxQP: 0 [default] dense, 1 sparse
  intProxQPBackend ,
  /// the factorization of the dense backend
  /**< 0 [default] chosen by ProxQP on the size of the problem, 1 that of
   * the whole regularized KKT system (PrimalDualLDLT), 2 that of the
   * primal Hessian with the constraints added to it (PrimalLDLT); ignored
   * by the sparse backend. */
  intProxQPDenseBackend ,
  /// whether the duality gap is also checked by the stopping criterion
  intProxQPCheckDualityGap ,
  /// first allowed new int parameter for derived classes
  intLastAlgParProxQP
  };

 /// public enum for the double algorithmic parameters
 /** Public enum describing the different algorithmic parameters of double
  * type that ProxQPMILPSolver has in addition to those of MILPSolver. The
  * value dblLastAlgParProxQP is provided so that the list can be easily
  * further extended by derived classes. Note that the absolute tolerance
  * on the primal and dual residuals (eps_abs of ProxQP) is the standard
  * dblFAccSol. */
 enum dbl_par_type_ProxQP {
  /// relative tolerance on the primal and dual residuals (eps_rel)
  dblProxQPEpsRel = dblLastAlgParMILP ,
  /// threshold under which the problem is declared primal infeasible
  /**< The default 1e-8 is tighter than the 1e-4 of ProxQP, with which
   * feasible problems whose costs are large (e.g., 1000) are declared
   * infeasible. */
  dblProxQPEpsPrimalInf ,
  /// threshold under which the problem is declared dual infeasible
  /**< The default is 1e-8, as for dblProxQPEpsPrimalInf. */
  dblProxQPEpsDualInf ,
  /// proximal step size on the primal variables (rho)
  dblProxQPRho ,
  /// initial penalty of the equality constraints (mu_eq)
  dblProxQPMuEq ,
  /// initial penalty of the inequality constraints (mu_in)
  dblProxQPMuIn ,
  /// first allowed new double parameter for derived classes
  dblLastAlgParProxQP
  };

/** @} ---------------------------------------------------------------------*/
/*------------------- CONSTRUCTING AND DESTRUCTING -------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Constructing and destructing ProxQPMILPSolver
 * @{ */

 /// void constructor
 ProxQPMILPSolver( void );

/*--------------------------------------------------------------------------*/
 /// destructor: releases the solver of ProxQP

 ~ProxQPMILPSolver() override;

/** @} ---------------------------------------------------------------------*/
/*--------------------- DERIVED METHODS OF BASE CLASS ----------------------*/
/*--------------------------------------------------------------------------*/
/** @name Methods of MILPSolver that are redefined
 * @{ */

 /// clears the problem, together with the solver of ProxQP

 void clear_problem( unsigned int what ) override;

/*--------------------------------------------------------------------------*/
 /// loads the problem from the Block into the data of ProxQP

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
/** @name Handling the parameters of ProxQPMILPSolver
 *
 * Besides those of MILPSolver and those of ProxQP above, the standard
 * intMaxIter (the iteration limit of ProxQP) and dblFAccSol (its absolute
 * tolerance on the primal and dual residuals) are handled.
 * @{ */

 void set_par( idx_type par , int value ) override;

 void set_par( idx_type par , double value ) override;

 [[nodiscard]] idx_type get_num_int_par( void ) const override {
  return( MILPSolver::get_num_int_par()
	  + intLastAlgParProxQP - intLastAlgParMILP );
  }

 [[nodiscard]] idx_type get_num_dbl_par( void ) const override {
  return( MILPSolver::get_num_dbl_par()
	  + dblLastAlgParProxQP - dblLastAlgParMILP );
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

 /// sets up ProxQP on the loaded problem, if needed, and solves it

 int guts_of_compute( void ) override;

 /// ProxQP is not given the changes of a problem: it is always loaded again

 bool reload_on_modification( void ) override { return( true ); }

/*--------------------------------------------------------------------------*/
/*---------------------------- PRIVATE PART --------------------------------*/
/*--------------------------------------------------------------------------*/

 private:

 /// translates a status of ProxQP into a return status of compute()

 static int decode_status( int status );

 /// releases the solver of ProxQP

 void free_solver( void );

 /// true if the last solve has given a primal solution

 [[nodiscard]] bool solved( void ) const;

/*--------------------------------------------------------------------------*/
/*--------------------------- PRIVATE FIELDS -------------------------------*/
/*--------------------------------------------------------------------------*/

 using SpMat = Eigen::SparseMatrix< double , Eigen::ColMajor , int >;
 using DQP = proxsuite::proxqp::dense::QP< double >;
 using SQP = proxsuite::proxqp::sparse::QP< double , int >;

 // the problem in the format of ProxQP

 int d_n = 0;                    ///< number of variables
 int d_neq = 0;                  ///< number of equalities
 int d_nin = 0;                  ///< number of inequalities
 double d_sign = 1;              ///< -1 if a maximization was negated

 SpMat d_H;                      ///< Hessian
 Eigen::VectorXd d_g;            ///< linear term
 SpMat d_A;                      ///< equality rows
 Eigen::VectorXd d_b;            ///< right-hand side of the equalities
 SpMat d_C;                      ///< inequality rows
 Eigen::VectorXd d_l;            ///< lower sides of the inequalities
 Eigen::VectorXd d_u;            ///< upper sides of the inequalities
 Eigen::VectorXd d_lbox;         ///< lower bounds (dense backend)
 Eigen::VectorXd d_ubox;         ///< upper bounds (dense backend)

 /// position of each row: k >= 0 in A, or - k - 1 in C
 std::vector< int > d_row_pos;
 /// position in C of the bounds of each variable (sparse backend), or -1
 std::vector< int > d_box_pos;

 std::unique_ptr< DQP > d_dqp;   ///< the dense solver, if any
 std::unique_ptr< SQP > d_sqp;   ///< the sparse solver, if any
 bool d_loaded = false;          ///< true if the problem has changed
 int d_backend = 0;              ///< backend of the loaded problem

 Eigen::VectorXd d_x;            ///< primal solution
 Eigen::VectorXd d_y;            ///< multipliers of the equalities
 Eigen::VectorXd d_z;            ///< of the inequalities, then of the box
 int d_ws_backend = -1;          ///< backend that gave d_x, d_y and d_z
 double d_fval = 0;              ///< value of the last solution
 int d_status = -1;              ///< the status of the last solve

 bool d_warm_start = true;           ///< intProxQPWarmStart
 int d_initial_guess = 1;            ///< intProxQPInitialGuess
 int d_par_backend = 0;              ///< intProxQPBackend
 int d_dense_backend = 0;            ///< intProxQPDenseBackend
 bool d_check_gap = false;           ///< intProxQPCheckDualityGap
 int d_max_iter = Inf< int >();      ///< intMaxIter
 double d_eps_abs;                   ///< dblFAccSol
 double d_eps_rel;                   ///< dblProxQPEpsRel
 double d_eps_primal_inf;            ///< dblProxQPEpsPrimalInf
 double d_eps_dual_inf;              ///< dblProxQPEpsDualInf
 double d_rho;                       ///< dblProxQPRho
 double d_mu_eq;                     ///< dblProxQPMuEq
 double d_mu_in;                     ///< dblProxQPMuIn

/*--------------------------------------------------------------------------*/

 SMSpp_insert_in_factory_h;

/*--------------------------------------------------------------------------*/

 };  // end( class( ProxQPMILPSolver ) )

/*--------------------------------------------------------------------------*/

}  // end( namespace SMSpp_di_unipi_it )

/*--------------------------------------------------------------------------*/

#endif  /* ProxQPMILPSolver.h included */

/*--------------------------------------------------------------------------*/
/*--------------------- End File ProxQPMILPSolver.h ------------------------*/
/*--------------------------------------------------------------------------*/
