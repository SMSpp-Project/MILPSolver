/*--------------------------------------------------------------------------*/
/*-------------------------- File PIQPMILPSolver.h -------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the PIQPMILPSolver class, implementing the MILPSolver
 * interface for convex Quadratic Programs with linear constraints using the
 * proximal interior point solver PIQP.
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

#ifndef __PIQPMILPSOLVER_H
 #define __PIQPMILPSOLVER_H
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
/*--------------------------- NAMESPACE ------------------------------------*/
/*--------------------------------------------------------------------------*/

/// namespace for the Structured Modeling System++ (SMS++)
namespace SMSpp_di_unipi_it
{

/*--------------------------------------------------------------------------*/
/*------------------------ CLASS PIQPMILPSolver ----------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------- GENERAL NOTES --------------------------------*/
/*--------------------------------------------------------------------------*/
/// class for solving convex QP Blocks via the solver PIQP
/** The PIQPMILPSolver implements the MILPSolver interface using PIQP, the
 * open-source proximal interior point solver for convex Quadratic Programs
 *
 * \f[
 *   \min \{ \tfrac{1}{2} x^\top P x + c^\top x \,:\,
 *           A x = b \,,\, h_l \leq G x \leq h_u \,,\,
 *           x_l \leq x \leq x_u \} \,,
 * \f]
 *
 * where \f$ P \f$ only needs to be positive semidefinite (a zero \f$ P \f$,
 * i.e., that of a Linear Program, is given to PIQP as such). It is
 * registered in the Solver factory, hence it can be selected by name
 * ("PIQPMILPSolver") in a BlockSolverConfig.
 *
 * The rows of the problem whose two sides coincide are the equalities
 * \f$ A x = b \f$, those with no finite side are left out, all the others
 * are the inequalities \f$ h_l \leq G x \leq h_u \f$, and the bounds on the
 * variables are the box constraints \f$ x_l \leq x \leq x_u \f$ of PIQP.
 * Two backends of PIQP can be used (see intPIQPBackend): the dense one
 * [default], meant for the small and medium problems, and the sparse one,
 * for the large and sparse ones.
 *
 * The following restrictions apply:
 *
 * - quadratic constraints are not supported, and integer variables are only
 *   accepted if relaxed via intRelaxIntVars;
 *
 * - PIQP is not given the changes of the loaded problem: whenever a
 *   Modification is received the problem is loaded again whole from the
 *   Block [see MILPSolver::reload_on_modification()];
 *
 * - PIQP has no warm start, hence each solve starts from the initial point
 *   of the interior point method, whatever the previous one;
 *
 * - PIQP has no time limit, hence dblMaxTime is not handled and compute()
 *   never returns kStopTime;
 *
 * - the infeasibility criterion of PIQP does not always succeed: an
 *   infeasible problem may end with kStopIter at the iteration limit
 *   rather than with kInfeasible, e.g., when the conflict is between rows
 *   and bounds; a lower dblPIQPInfeasThreshold catches some of these, but
 *   has unbounded problems declared infeasible;
 *
 * - the solution and the multipliers are accurate to the tolerances of
 *   PIQP (dblFAccSol and dblPIQPEpsRel), e.g., a multiplier of 1 may come
 *   out as 0.999999999, and the value of the solution is computed from it;
 *
 * - no unbounded ray or infeasibility certificate is provided, i.e., the
 *   has_var_direction() / has_dual_direction() interface of CDASolver is
 *   not implemented. */

class PIQPMILPSolver : public MILPSolver
{
/*--------------------------------------------------------------------------*/
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

 public:

/*--------------------------------------------------------------------------*/
/*---------------------------- PUBLIC TYPES --------------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Public types of PIQPMILPSolver
 * @{ */

 /// public enum for the integer algorithmic parameters
 /** Public enum describing the different algorithmic parameters of int type
  * that PIQPMILPSolver has in addition to those of MILPSolver. The value
  * intLastAlgParPIQP is provided so that the list can be easily further
  * extended by derived classes. */
 enum int_par_type_PIQP {
  /// the backend of PIQP: 0 [default] dense, 1 sparse
  intPIQPBackend = intLastAlgParMILP ,
  /// whether the duality gap is also checked by the stopping criterion
  intPIQPCheckDualityGap ,
  /// maximum number of iterations of the Ruiz equilibration of the data
  intPIQPPrecondIter ,
  /// whether the Ruiz equilibration also scales the cost
  intPIQPPrecondScaleCost ,
  /// whether the iterative refinement is run at every iteration
  /**< If 0 [default], the iterative refinement of the solution of the KKT
   * system is only run when its factorization fails; if 1, it is run at
   * every iteration. */
  intPIQPIterRefinement ,
  /// first allowed new int parameter for derived classes
  intLastAlgParPIQP
  };

 /// public enum for the double algorithmic parameters
 /** Public enum describing the different algorithmic parameters of double
  * type that PIQPMILPSolver has in addition to those of MILPSolver. The
  * value dblLastAlgParPIQP is provided so that the list can be easily
  * further extended by derived classes. Note that the absolute tolerance
  * on the primal and dual residuals (eps_abs of PIQP) is the standard
  * dblFAccSol. */
 enum dbl_par_type_PIQP {
  /// relative tolerance on the primal and dual residuals (eps_rel)
  dblPIQPEpsRel = dblLastAlgParMILP ,
  /// absolute tolerance on the duality gap (eps_duality_gap_abs)
  dblPIQPEpsGapAbs ,
  /// relative tolerance on the duality gap (eps_duality_gap_rel)
  dblPIQPEpsGapRel ,
  /// threshold of the infeasibility detection (infeasibility_threshold)
  dblPIQPInfeasThreshold ,
  /// initial proximal penalty on the primal variables (rho_init)
  dblPIQPRhoInit ,
  /// initial augmented Lagrangian penalty (delta_init)
  dblPIQPDeltaInit ,
  /// first allowed new double parameter for derived classes
  dblLastAlgParPIQP
  };

/** @} ---------------------------------------------------------------------*/
/*------------------- CONSTRUCTING AND DESTRUCTING -------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Constructing and destructing PIQPMILPSolver
 * @{ */

 /// void constructor
 PIQPMILPSolver( void );

/*--------------------------------------------------------------------------*/
 /// destructor: releases the solver of PIQP

 ~PIQPMILPSolver() override;

/** @} ---------------------------------------------------------------------*/
/*--------------------- DERIVED METHODS OF BASE CLASS ----------------------*/
/*--------------------------------------------------------------------------*/
/** @name Methods of MILPSolver that are redefined
 * @{ */

 /// clears the problem, together with the solver of PIQP

 void clear_problem( unsigned int what ) override;

/*--------------------------------------------------------------------------*/
 /// loads the problem from the Block into the data of PIQP

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
/** @name Handling the parameters of PIQPMILPSolver
 *
 * Besides those of MILPSolver and those of PIQP above, the standard
 * intMaxIter (the iteration limit of PIQP) and dblFAccSol (its absolute
 * tolerance on the primal and dual residuals) are handled.
 * @{ */

 void set_par( idx_type par , int value ) override;

 void set_par( idx_type par , double value ) override;

 [[nodiscard]] idx_type get_num_int_par( void ) const override {
  return( MILPSolver::get_num_int_par()
	  + intLastAlgParPIQP - intLastAlgParMILP );
  }

 [[nodiscard]] idx_type get_num_dbl_par( void ) const override {
  return( MILPSolver::get_num_dbl_par()
	  + dblLastAlgParPIQP - dblLastAlgParMILP );
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

 /// sets up PIQP on the loaded problem, if needed, and solves it

 int guts_of_compute( void ) override;

 /// PIQP is not given the changes of a problem: it is always loaded again

 bool reload_on_modification( void ) override { return( true ); }

/*--------------------------------------------------------------------------*/
/*---------------------------- PRIVATE PART --------------------------------*/
/*--------------------------------------------------------------------------*/

 private:

 /// translates a status of PIQP into a return status of compute()

 static int decode_status( int status );

 /// releases the solver of PIQP

 void free_solver( void );

 /// true if the last solve has given a primal solution

 [[nodiscard]] bool solved( void ) const;

/*--------------------------------------------------------------------------*/
/*--------------------------- PRIVATE FIELDS -------------------------------*/
/*--------------------------------------------------------------------------*/

 using SpMat = Eigen::SparseMatrix< double , Eigen::ColMajor , int >;

 /// the solvers of PIQP, whose headers are only included by the .cpp
 struct Solvers;

 // the problem in the format of PIQP

 int d_n = 0;                    ///< number of variables
 int d_neq = 0;                  ///< number of equalities
 int d_nin = 0;                  ///< number of inequalities
 double d_sign = 1;              ///< -1 if a maximization was negated

 SpMat d_P;                      ///< Hessian
 Eigen::VectorXd d_c;            ///< linear term
 SpMat d_A;                      ///< equality rows
 Eigen::VectorXd d_b;            ///< right-hand side of the equalities
 SpMat d_G;                      ///< inequality rows
 Eigen::VectorXd d_hl;           ///< lower sides of the inequalities
 Eigen::VectorXd d_hu;           ///< upper sides of the inequalities
 Eigen::VectorXd d_xl;           ///< lower bounds
 Eigen::VectorXd d_xu;           ///< upper bounds

 /// position of each row: k >= 0 in A, - k - 1 in G, or the minimum int
 /// if the row has no finite side and is left out of PIQP
 std::vector< int > d_row_pos;

 std::unique_ptr< Solvers > d_slv;  ///< the solver of PIQP, if any
 bool d_loaded = false;          ///< true if the problem has changed
 int d_backend = 0;              ///< backend of the loaded problem

 Eigen::VectorXd d_x;            ///< primal solution
 Eigen::VectorXd d_y;            ///< multipliers of the equalities
 Eigen::VectorXd d_z;            ///< of the inequalities, upper - lower
 Eigen::VectorXd d_zb;           ///< of the bounds, upper - lower
 double d_fval = 0;              ///< value of the last solution
 int d_status = 0;               ///< the status of the last solve, 0 none

 int d_par_backend = 0;              ///< intPIQPBackend
 bool d_check_gap;                   ///< intPIQPCheckDualityGap
 int d_precond_iter;                 ///< intPIQPPrecondIter
 bool d_precond_cost;                ///< intPIQPPrecondScaleCost
 bool d_iter_ref;                    ///< intPIQPIterRefinement
 int d_max_iter = Inf< int >();      ///< intMaxIter
 double d_eps_abs;                   ///< dblFAccSol
 double d_eps_rel;                   ///< dblPIQPEpsRel
 double d_eps_gap_abs;               ///< dblPIQPEpsGapAbs
 double d_eps_gap_rel;               ///< dblPIQPEpsGapRel
 double d_infeas_thr;                ///< dblPIQPInfeasThreshold
 double d_rho_init;                  ///< dblPIQPRhoInit
 double d_delta_init;                ///< dblPIQPDeltaInit

/*--------------------------------------------------------------------------*/

 SMSpp_insert_in_factory_h;

/*--------------------------------------------------------------------------*/

 };  // end( class( PIQPMILPSolver ) )

/*--------------------------------------------------------------------------*/

}  // end( namespace SMSpp_di_unipi_it )

/*--------------------------------------------------------------------------*/

#endif  /* PIQPMILPSolver.h included */

/*--------------------------------------------------------------------------*/
/*---------------------- End File PIQPMILPSolver.h -------------------------*/
/*--------------------------------------------------------------------------*/
