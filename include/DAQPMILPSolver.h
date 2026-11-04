/*--------------------------------------------------------------------------*/
/*-------------------------- File DAQPMILPSolver.h -------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the DAQPMILPSolver class, implementing the MILPSolver
 * interface for convex Quadratic Programs with linear constraints using the
 * dual active-set solver DAQP.
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

#ifndef __DAQPMILPSOLVER_H
 #define __DAQPMILPSOLVER_H
                      /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include <string>
#include <vector>

#include "MILPSolver.h"

#include <daqp/api.h>

/*--------------------------------------------------------------------------*/
/*--------------------------- NAMESPACE ------------------------------------*/
/*--------------------------------------------------------------------------*/

/// namespace for the Structured Modeling System++ (SMS++)
namespace SMSpp_di_unipi_it
{

/*--------------------------------------------------------------------------*/
/*------------------------ CLASS DAQPMILPSolver ----------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------- GENERAL NOTES --------------------------------*/
/*--------------------------------------------------------------------------*/
/// class for solving convex QP Blocks via the dual active-set solver DAQP
/** The DAQPMILPSolver implements the MILPSolver interface using DAQP, an
 * open-source dual active-set solver for convex Quadratic Programs
 *
 * \f[
 *   \min \{ \tfrac{1}{2} x^\top H x + f^\top x \,:\,
 *           l \leq x \leq u \,,\, b_l \leq A x \leq b_u \} \,,
 * \f]
 *
 * where \f$ H \f$ is dense and only needs to be positive semidefinite (a
 * singular or zero \f$ H \f$, e.g., that of a Linear Program, is dealt with
 * by the proximal-point iterations of DAQP). It is registered in the Solver
 * factory, hence it can be selected by name ("DAQPMILPSolver") in a
 * BlockSolverConfig.
 *
 * The data are kept dense: DAQP is meant for the small and medium problems,
 * typically solved many times with small changes, where an active-set
 * method is at its best, since the working set of the previous solution is
 * a good start for the next one (see intDAQPWarmStart); the iterations of
 * DAQP in the last compute() are given by get_elapsed_iterations().
 *
 * The following restrictions apply:
 *
 * - quadratic constraints are not supported, and integer variables are only
 *   accepted if relaxed via intRelaxIntVars;
 *
 * - DAQP does not support changing the loaded problem in place: whenever a
 *   Modification is received the problem is loaded again whole from the
 *   Block [see MILPSolver::reload_on_modification()], and the working set
 *   of the previous solution, on the variables and rows that are still
 *   there, is the start of the next one, but for the constraints that
 *   either problem fixes (a variable with equal bounds, an equality row),
 *   which are active in the new one if and only if it fixes them;
 *
 * - a problem whose H is singular, e.g., a Linear Program, is solved by
 *   proximal-point iterations, whose solution and multipliers are accurate
 *   to about the regularization (see dblDAQPEpsProx), e.g., a multiplier
 *   of 1 may come out as 0.999999; the value of the solution is computed
 *   from it rather than taken from DAQP;
 *
 * - no unbounded ray or infeasibility certificate is provided, i.e., the
 *   has_var_direction() / has_dual_direction() interface of CDASolver is
 *   not implemented. */

class DAQPMILPSolver : public MILPSolver
{
/*--------------------------------------------------------------------------*/
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

 public:

/*--------------------------------------------------------------------------*/
/*---------------------------- PUBLIC TYPES --------------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Public types of DAQPMILPSolver
 * @{ */

 /// public enum for the integer algorithmic parameters
 /** Public enum describing the different algorithmic parameters of int type
  * that DAQPMILPSolver has in addition to those of MILPSolver. The value
  * intLastAlgParDAQP is provided so that the list can be easily further
  * extended by derived classes. */
 enum int_par_type_DAQP {
  /// whether a solve starts from the working set of the previous one
  /**< If 1 [default], each solve after the first starts from the working
   * set, i.e., the set of the constraints active at the solution, of the
   * previous one, on the variables and rows that are still there [see
   * MILPSolver::previous_column_index()], the new ones inactive; if 0,
   * each solve starts from the empty working set. */
  intDAQPWarmStart = intLastAlgParMILP ,
  /// how the equality constraints are eliminated before the solve
  /**< The value of eq_reduction in DAQPSettings: -1 never, 0 [default]
   * only in a solve from scratch, 1 also in a warm-started one; ignored by
   * the versions of DAQP (up to 0.9) that have no elimination. */
  intDAQPEqReduction ,
  /// number of non-improving iterations after which DAQP declares cycling
  intDAQPCycleTol ,
  /// first allowed new int parameter for derived classes
  intLastAlgParDAQP
  };

 /// public enum for the double algorithmic parameters
 /** Public enum describing the different algorithmic parameters of double
  * type that DAQPMILPSolver has in addition to those of MILPSolver. The
  * value dblLastAlgParDAQP is provided so that the list can be easily
  * further extended by derived classes. Note that the tolerance on the
  * primal feasibility is the standard dblFAccSol. */
 enum dbl_par_type_DAQP {
  /// tolerance on the dual feasibility, i.e., on the sign of multipliers
  dblDAQPDualTol = dblLastAlgParMILP ,
  /// value below which a number is taken as zero
  dblDAQPZeroTol ,
  /// tolerance on the pivots of the updates of the factorization
  dblDAQPPivotTol ,
  /// regularization of the proximal-point iterations for a singular H
  /**< The value of eps_prox in DAQPSettings: if negative [default], the
   * proximal-point iterations are used, with regularization equal to the
   * absolute value, only if H is found to be singular; if positive, they
   * are always used; if 0, never. */
  dblDAQPEpsProx ,
  /// tolerance of the proximal-point iterations (< 0 for automatic)
  dblDAQPEtaProx ,
  /// first allowed new double parameter for derived classes
  dblLastAlgParDAQP
  };

/** @} ---------------------------------------------------------------------*/
/*------------------- CONSTRUCTING AND DESTRUCTING -------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Constructing and destructing DAQPMILPSolver
 * @{ */

 /// void constructor
 DAQPMILPSolver( void );

/*--------------------------------------------------------------------------*/
 /// destructor: releases the DAQP workspace

 ~DAQPMILPSolver() override;

/** @} ---------------------------------------------------------------------*/
/*--------------------- DERIVED METHODS OF BASE CLASS ----------------------*/
/*--------------------------------------------------------------------------*/
/** @name Methods of MILPSolver that are redefined
 * @{ */

 /// clears the problem, together with the DAQP workspace

 void clear_problem( unsigned int what ) override;

/*--------------------------------------------------------------------------*/
 /// loads the problem from the Block into the dense data of DAQP

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
/** @name Handling the parameters of DAQPMILPSolver
 *
 * Besides those of MILPSolver and those of DAQP above, the standard
 * intMaxIter (the iteration limit of DAQP), dblMaxTime (its time limit, in
 * seconds) and dblFAccSol (its tolerance on the primal feasibility) are
 * handled.
 * @{ */

 void set_par( idx_type par , int value ) override;

 void set_par( idx_type par , double value ) override;

 [[nodiscard]] idx_type get_num_int_par( void ) const override {
  return( MILPSolver::get_num_int_par()
	  + intLastAlgParDAQP - intLastAlgParMILP );
  }

 [[nodiscard]] idx_type get_num_dbl_par( void ) const override {
  return( MILPSolver::get_num_dbl_par()
	  + dblLastAlgParDAQP - dblLastAlgParMILP );
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

 /// sets up DAQP on the loaded problem, if needed, and solves it

 int guts_of_compute( void ) override;

 /// DAQP does not change a problem in place: it is always loaded again

 bool reload_on_modification( void ) override { return( true ); }

/*--------------------------------------------------------------------------*/
/*---------------------------- PRIVATE PART --------------------------------*/
/*--------------------------------------------------------------------------*/

 private:

 /// translates an exit flag of DAQP into a return status of compute()

 static int decode_exitflag( int exitflag );

 /// frees the DAQP workspace, keeping its working set for the next start

 void free_workspace( void );

 /// true if the last solve has given a primal solution

 [[nodiscard]] bool solved( void ) const;

/*--------------------------------------------------------------------------*/
/*--------------------------- PRIVATE FIELDS -------------------------------*/
/*--------------------------------------------------------------------------*/

 // the problem in the format of DAQP, kept as long as the workspace is,
 // which points to it

 int d_n = 0;                    ///< number of variables
 int d_m = 0;                    ///< number of constraints, bounds included
 double d_sign = 1;              ///< -1 if a maximization was negated

 std::vector< double > d_H;      ///< dense Hessian, row-major
 std::vector< double > d_f;      ///< linear term
 std::vector< double > d_A;      ///< dense rows, row-major
 std::vector< double > d_bu;     ///< upper sides: bounds, then rows
 std::vector< double > d_bl;     ///< lower sides: bounds, then rows
 std::vector< int > d_sense;     ///< sense of each constraint
 std::vector< int > d_fixed;     ///< sense the problem fixes (equalities)

 DAQPProblem d_qp{};             ///< the problem DAQP is given
 DAQPSettings d_settings{};      ///< the settings DAQP is given
 DAQPWorkspace d_work{};         ///< the workspace of DAQP
 bool d_has_work = false;        ///< true if d_work is set up
 bool d_loaded = false;          ///< true if the problem has changed

 std::vector< int > d_ws;        ///< working set of the previous solve
 int d_ws_n = -1;                ///< number of variables when d_ws was taken

 std::vector< double > d_x;      ///< primal solution
 std::vector< double > d_lam;    ///< multipliers, bounds then rows
 DAQPResult d_res{};             ///< the result of the last solve
 double d_fval = 0;              ///< value of the last solution
 int d_exitflag = 0;             ///< the exit flag of the last solve

 bool d_warm_start = true;       ///< intDAQPWarmStart
 int d_max_iter = Inf< int >();       ///< intMaxIter
 double d_max_time = Inf< double >(); ///< dblMaxTime

/*--------------------------------------------------------------------------*/

 SMSpp_insert_in_factory_h;

/*--------------------------------------------------------------------------*/

 };  // end( class( DAQPMILPSolver ) )

/*--------------------------------------------------------------------------*/

}  // end( namespace SMSpp_di_unipi_it )

/*--------------------------------------------------------------------------*/

#endif  /* DAQPMILPSolver.h included */

/*--------------------------------------------------------------------------*/
/*---------------------- End File DAQPMILPSolver.h -------------------------*/
/*--------------------------------------------------------------------------*/
