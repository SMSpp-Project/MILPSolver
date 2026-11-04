/*--------------------------------------------------------------------------*/
/*------------------------- File PERMONMILPSolver.h ------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the PERMONMILPSolver class, implementing the MILPSolver
 * interface for convex Quadratic Programs with linear constraints using the
 * solvers of PermonQP, the QP module of the PERMON library built on PETSc.
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

#ifndef __PERMONMILPSOLVER_H
 #define __PERMONMILPSOLVER_H
/* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include <memory>
#include <string>
#include <vector>

#include "MILPSolver.h"

/*--------------------------------------------------------------------------*/
/*--------------------------- NAMESPACE ------------------------------------*/
/*--------------------------------------------------------------------------*/

/// namespace for the Structured Modeling System++ (SMS++)
namespace SMSpp_di_unipi_it
{

/*--------------------------------------------------------------------------*/
/*----------------------- CLASS PERMONMILPSolver ---------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------- GENERAL NOTES --------------------------------*/
/*--------------------------------------------------------------------------*/
/// class for solving convex QP Blocks via the solvers of PermonQP
/** The PERMONMILPSolver implements the MILPSolver interface using PermonQP,
 * the open-source library of solvers for large sparse convex Quadratic
 * Programs of PERMON, which is built on PETSc. It is registered in the
 * Solver factory, hence it can be selected by name ("PERMONMILPSolver") in
 * a BlockSolverConfig.
 *
 * The solvers of PermonQP handle a quadratic objective, linear equalities
 * and bounds on the variables: the MPGP method (modified proportioning with
 * reduced gradient projections) the bounds alone, and the SMALXE method (a
 * semi-monotonic augmented Lagrangian, whose inner problems are solved by
 * MPGP) the equalities as well. The problem is therefore given to PermonQP
 * in the form
 *
 * \f[
 *   \min \{ \tfrac{1}{2} x^\top H x + g^\top x \,:\,
 *           A x = b \,,\, C x - s = 0 \,,\,
 *           l_x \leq x \leq u_x \,,\, l \leq s \leq u \} \,,
 * \f]
 *
 * where \f$ A x = b \f$ are the rows whose two sides coincide and each other
 * row \f$ l_i \leq C_i x \leq u_i \f$ has a slack \f$ s_i \f$ with the sides
 * of the row as bounds, while a variable whose two bounds coincide is fixed
 * by a further equality (MPGP cannot handle it in the bounds); PermonQP then
 * chooses SMALXE if there are rows and MPGP otherwise. The method and all
 * the settings of PETSc and PERMON that are not parameters of this class can
 * be given in the options database of PETSc (e.g., in the environment
 * variable PETSC_OPTIONS), with the prefix "qps_" of PermonQP.
 *
 * The following restrictions apply:
 *
 * - quadratic constraints are not supported, and integer variables are only
 *   accepted if relaxed via intRelaxIntVars;
 *
 * - the methods of PermonQP are meant for strictly convex problems, or at
 *   least for problems whose Hessian is positive definite on the null space
 *   of the equalities (as the dual of the training of a Support Vector
 *   Machine); a Linear Program, or a QP whose Hessian is singular on a
 *   direction along which the problem is unbounded, may not converge;
 *
 * - PermonQP has no infeasibility or unboundedness detection: such a
 *   problem ends with kStopIter at the iteration limit, or with kError,
 *   and never with kInfeasible or kUnbounded;
 *
 * - up to PERMON 3.26, with its default expansion step MPGP may stop as
 *   converged at a point out of the bounds, hence there the expansion step
 *   is projcg unless the options of PETSc say otherwise
 *   (-qps_mpgp_expansion_type, and -smalxe_qps_mpgp_expansion_type for the
 *   inner solver of SMALXE); also,
 *   on a problem whose Hessian is singular along a direction of free
 *   variables MPGP may drift along it and stop there as converged, with a
 *   value that is only rounding noise, hence dblKKTTol is 1e-6 by default
 *   (rather than the 0 of MILPSolver), so that such a solve is checked and
 *   ends with kError;
 *
 * - PermonQP is not given the changes of the loaded problem: whenever a
 *   Modification is received the problem is loaded again whole from the
 *   Block [see MILPSolver::reload_on_modification()], and the solution of
 *   the previous solve, on the variables that are still there, is the start
 *   of the next one (see intPERMONWarmStart);
 *
 * - PermonQP has no time limit, hence dblMaxTime is not handled and
 *   compute() never returns kStopTime;
 *
 * - the solution and the multipliers are accurate to the tolerances of
 *   PermonQP (dblFAccSol and dblPERMONRelTol), and the value of the
 *   solution is computed from it; PermonQP gives the multipliers of the
 *   bounds, while those of the rows are computed from the stationarity
 *   condition (exactly for a row with a slack, by least squares for the
 *   equalities);
 *
 * - PETSc is initialized (sequential, on PETSC_COMM_SELF) the first time a
 *   problem is solved, unless the program has already done it, and it is
 *   never finalized. */

class PERMONMILPSolver : public MILPSolver
{
/*--------------------------------------------------------------------------*/
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

 public:
/*--------------------------------------------------------------------------*/
/*---------------------------- PUBLIC TYPES --------------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Public types of PERMONMILPSolver
 * @{ */

 /// public enum for the integer algorithmic parameters
 /** Public enum describing the different algorithmic parameters of int type
  * that PERMONMILPSolver has in addition to those of MILPSolver. The value
  * intLastAlgParPERMON is provided so that the list can be easily further
  * extended by derived classes. */
 enum int_par_type_PERMON {
  /// whether a solve starts from the solution of the previous one
  /**< If 1 [default], each solve after the first starts from the primal
   * solution of the previous one, on the variables that are still there
   * [see MILPSolver::previous_column_index()], the new ones at 0 (or at
   * their bound closest to 0), and the slack of each row at the value of
   * the row there; if 0, each solve starts from the point of the variables
   * at 0 (or at their bound closest to 0). */
  intPERMONWarmStart = intLastAlgParMILP ,
  /// first allowed new int parameter for derived classes
  intLastAlgParPERMON
  };

 /// public enum for the double algorithmic parameters
 /** Public enum describing the different algorithmic parameters of double
  * type that PERMONMILPSolver has in addition to those of MILPSolver. The
  * value dblLastAlgParPERMON is provided so that the list can be easily
  * further extended by derived classes. Note that the absolute tolerance
  * of PermonQP on the norm of the projected gradient (atol of PETSc) is
  * the standard dblFAccSol. */
 enum dbl_par_type_PERMON {
  /// relative tolerance of PermonQP (rtol of PETSc)
  /**< PermonQP stops when the norm of the projected gradient (and, with
   * SMALXE, that of the residual of the equalities) is below the larger of
   * dblFAccSol and this tolerance times the norm of the linear term. */
  dblPERMONRelTol = dblLastAlgParMILP ,
  /// first allowed new double parameter for derived classes
  dblLastAlgParPERMON
  };

/** @} ---------------------------------------------------------------------*/
/*------------------- CONSTRUCTING AND DESTRUCTING -------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Constructing and destructing PERMONMILPSolver
 * @{ */

 /// void constructor
 PERMONMILPSolver( void );

/*--------------------------------------------------------------------------*/
 /// destructor: releases the objects of PermonQP

 ~PERMONMILPSolver() override;

/** @} ---------------------------------------------------------------------*/
/*--------------------- DERIVED METHODS OF BASE CLASS ----------------------*/
/*--------------------------------------------------------------------------*/
/** @name Methods of MILPSolver that are redefined
 * @{ */

 /// clears the problem, together with the objects of PermonQP

 void clear_problem( unsigned int what ) override;

/*--------------------------------------------------------------------------*/
 /// loads the problem from the Block into the data of PermonQP

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
/** @name Handling the parameters of PERMONMILPSolver
 *
 * Besides those of MILPSolver and those of PermonQP above, the standard
 * intMaxIter (the iteration limit of PermonQP, i.e., the outer iterations
 * of SMALXE or those of MPGP) and dblFAccSol (its absolute tolerance) are
 * handled.
 * @{ */

 void set_par( idx_type par , int value ) override;

 void set_par( idx_type par , double value ) override;

 [[nodiscard]] idx_type get_num_int_par( void ) const override {
  return( MILPSolver::get_num_int_par() + intLastAlgParPERMON -
	   intLastAlgParMILP );
  }

 [[nodiscard]] idx_type get_num_dbl_par( void ) const override {
  return( MILPSolver::get_num_dbl_par() + dblLastAlgParPERMON -
	   dblLastAlgParMILP );
  }

 [[nodiscard]] int get_dflt_int_par( idx_type par ) const override;

 [[nodiscard]] double get_dflt_dbl_par( idx_type par ) const override;

 [[nodiscard]] int get_int_par( idx_type par ) const override;

 [[nodiscard]] double get_dbl_par( idx_type par ) const override;

 [[nodiscard]] idx_type
 int_par_str2idx( const std::string & name ) const override;

 [[nodiscard]] const std::string &
 int_par_idx2str( idx_type idx ) const override;

 [[nodiscard]] idx_type
 dbl_par_str2idx( const std::string & name ) const override;

 [[nodiscard]] const std::string &
 dbl_par_idx2str( idx_type idx ) const override;

/** @} ---------------------------------------------------------------------*/
/*--------------------- PROTECTED PART OF THE CLASS ------------------------*/
/*--------------------------------------------------------------------------*/

 protected:
/*--------------------------------------------------------------------------*/
/*-------------------------- PROTECTED METHODS -----------------------------*/
/*--------------------------------------------------------------------------*/

 /// sets up PermonQP on the loaded problem and solves it

 int guts_of_compute( void ) override;

 /// PermonQP is not given the changes of a problem: it is always loaded again

 bool reload_on_modification( void ) override { return( true ); }

/*--------------------------------------------------------------------------*/
/*---------------------------- PRIVATE PART --------------------------------*/
/*--------------------------------------------------------------------------*/

 private:
 /// true if the last solve has given a primal solution

 [[nodiscard]] bool solved( void ) const;

/*--------------------------------------------------------------------------*/
/*--------------------------- PRIVATE FIELDS -------------------------------*/
/*--------------------------------------------------------------------------*/

 /// the objects of PETSc and PermonQP, defined in the .cpp, so that the
 /// headers of PETSc do not reach whoever includes this one
 struct PERMONData;
 std::unique_ptr< PERMONData > d_pd;

 int d_n = 0;     ///< number of variables
 int d_neq = 0;   ///< number of rows with equal sides
 int d_nin = 0;   ///< number of the other rows (slacks)
 int d_nfix = 0;  ///< number of fixed variables

 /// the equality of PermonQP that fixes each variable, or -1
 std::vector< int > d_fix_row;
 double d_sign = 1;  ///< -1 if a maximization was negated

 /// position of each row: k >= 0 among the equalities, or - k - 1 among
 /// the rows with a slack, which follow those that fix the variables
 std::vector< int > d_row_pos;
 std::vector< int > d_prev_row_pos;  ///< d_row_pos of the previous problem
 int d_prev_n = 0;                   ///< d_n of the previous problem
 bool d_remap = false;               ///< true if d_x is to be moved

 std::vector< double > d_x;  ///< primal solution (without the slacks)
 std::vector< double > d_y;  ///< multipliers of the rows
 std::vector< double > d_z;  ///< multipliers of the bounds
 double d_fval = 0;          ///< value of the last solution
 int d_status = 0;           ///< reason of the last stop (KSP), or 0

 bool d_warm_start = true;       ///< intPERMONWarmStart
 int d_max_iter = Inf< int >();  ///< intMaxIter
 double d_atol;                  ///< dblFAccSol
 double d_rtol;                  ///< dblPERMONRelTol

/*--------------------------------------------------------------------------*/

 SMSpp_insert_in_factory_h;

/*--------------------------------------------------------------------------*/

 };  // end( class( PERMONMILPSolver ) )

/*--------------------------------------------------------------------------*/

}  // end( namespace SMSpp_di_unipi_it )

/*--------------------------------------------------------------------------*/

#endif /* PERMONMILPSolver.h included */

/*--------------------------------------------------------------------------*/
/*--------------------- End File PERMONMILPSolver.h ------------------------*/
/*--------------------------------------------------------------------------*/
