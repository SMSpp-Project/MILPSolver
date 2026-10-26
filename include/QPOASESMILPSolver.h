/*--------------------------------------------------------------------------*/
/*------------------------ File QPOASESMILPSolver.h ------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the QPOASESMILPSolver class, implementing the MILPSolver
 * interface for convex Quadratic Programs with linear constraints using the
 * online active-set solver qpOASES.
 *
 * \author Donato Meoli \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \copyright &copy; by Donato Meoli
 */
/*--------------------------------------------------------------------------*/
/*----------------------------- DEFINITIONS --------------------------------*/
/*--------------------------------------------------------------------------*/

#ifndef __QPOASESMILPSOLVER_H
 #define __QPOASESMILPSOLVER_H
                      /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include <map>
#include <memory>
#include <string>
#include <vector>

#include "MILPSolver.h"

/*--------------------------------------------------------------------------*/
/*----------------------- FORWARD DECLARATIONS -----------------------------*/
/*--------------------------------------------------------------------------*/

// the classes of qpOASES are only declared here, so that the macros of its
// headers do not reach whoever includes this one

namespace qpOASES
{
 class Options;
 class SQProblem;
 }

/*--------------------------------------------------------------------------*/
/*--------------------------- NAMESPACE ------------------------------------*/
/*--------------------------------------------------------------------------*/

/// namespace for the Structured Modeling System++ (SMS++)
namespace SMSpp_di_unipi_it
{

/*--------------------------------------------------------------------------*/
/*---------------------- CLASS QPOASESMILPSolver ---------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------- GENERAL NOTES --------------------------------*/
/*--------------------------------------------------------------------------*/
/// class for solving convex QP Blocks via the active-set solver qpOASES
/** The QPOASESMILPSolver implements the MILPSolver interface using
 * qpOASES, an open-source implementation of the online active-set strategy
 * for convex Quadratic Programs
 *
 * \f[
 *   \min \{ \tfrac{1}{2} x^\top H x + g^\top x \,:\,
 *           l \leq x \leq u \,,\, l_A \leq A x \leq u_A \} \,,
 * \f]
 *
 * where \f$ H \f$ is dense and only needs to be positive semidefinite: a
 * direction of zero curvature found along the homotopy is followed until
 * a constraint blocks it, and if none does the problem is declared
 * unbounded, while a zero \f$ H \f$, i.e., a Linear Program, is solved
 * with the regularization of qpOASES (see intQPOASESEnableReg). It is
 * registered in the Solver factory, hence it can be selected by name
 * ("QPOASESMILPSolver") in a BlockSolverConfig.
 *
 * The data are kept dense: qpOASES is meant for the small and medium
 * problems, typically solved many times with small changes, where its
 * homotopy from the solution of the previous problem to that of the new
 * one is at its best (see intQPOASESHotStart); the iterations of qpOASES in
 * the last compute(), i.e., its working set recalculations, those of a hot
 * start repeated from scratch included, are given by
 * get_elapsed_iterations().
 *
 * The following restrictions apply:
 *
 * - quadratic constraints are not supported, and integer variables are only
 *   accepted if relaxed via intRelaxIntVars;
 *
 * - qpOASES does not support adding or removing variables or constraints:
 *   whenever a Modification is received the problem is loaded again whole
 *   from the Block [see MILPSolver::reload_on_modification()], and if its
 *   size has not changed the next solve is a hot start from the solution of
 *   the previous one, otherwise it starts from scratch; a hot start that
 *   does not end with an optimal solution, or stopped by the time limit, is
 *   repeated from scratch, so that a declaration of infeasibility or of
 *   unboundedness only comes from a solve from scratch;
 *
 * - the value of the solution is computed from it rather than taken from
 *   qpOASES, whose value is that of the regularized problem whenever the
 *   regularization is used;
 *
 * - no unbounded ray or infeasibility certificate is provided, i.e., the
 *   has_var_direction() / has_dual_direction() interface of CDASolver is
 *   not implemented. */

class QPOASESMILPSolver : public MILPSolver
{
/*--------------------------------------------------------------------------*/
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

 public:

/*--------------------------------------------------------------------------*/
/*---------------------------- PUBLIC TYPES --------------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Public types of QPOASESMILPSolver
 * @{ */

 /// public enum for the integer algorithmic parameters
 /** Public enum describing the different algorithmic parameters of int type
  * that QPOASESMILPSolver has in addition to those of MILPSolver. The value
  * intLastAlgParQPOASES is provided so that the list can be easily further
  * extended by derived classes. Except intQPOASESHotStart and
  * intQPOASESHessianType, they are fields of the Options of qpOASES. */
 enum int_par_type_QPOASES {
  /// whether a solve is a hot start from the solution of the previous one
  /**< If 1 [default], each solve after the first is a hot start of qpOASES
   * from the solution of the previous one, as long as that has been found
   * and the number of variables and constraints has not changed; if 0,
   * each solve starts from scratch. */
  intQPOASESHotStart = intLastAlgParMILP ,
  /// the set of Options of qpOASES the other parameters change
  /**< 0 [default] for the default Options of qpOASES, 1 for those it calls
   * reliable, 2 for those it calls MPC (faster, but with a singular H they
   * may give a wrong solution, or a finite one of an unbounded problem).
   * Setting it brings every Option back to the set, and then applies again
   * the other parameters of qpOASES that have been explicitly set. */
  intQPOASESPreset ,
  /// what qpOASES is told about H in a solve from scratch
  /**< 0 [default] if qpOASES finds it out itself (zero, identity, or else
   * taken as positive definite until a zero curvature shows up), 1 if H is
   * positive definite, 2 if it is positive semidefinite, which has it
   * regularized in a solve from scratch if intQPOASESEnableReg is 1. A hot
   * start with a new H always has qpOASES find it out. */
  intQPOASESHessianType ,
  /// whether a semidefinite or zero H is regularized (enableRegularisation)
  /**< If 1, H is added dblQPOASESEpsReg times its norm on the diagonal
   * whenever it is known to be semidefinite or zero; a zero H has it 1 in
   * any case. Note that the regularized problem of an unbounded one has a
   * finite solution, far from the origin. */
  intQPOASESEnableReg ,
  /// number of refinements of the solution of a regularized problem
  intQPOASESRegSteps ,
  /// whether the far bounds of the homotopy are used (enableFarBounds)
  intQPOASESFarBounds ,
  /// whether the bounds are flipped when H is singular
  /**< The value of enableFlippingBounds, which keeps the projected H
   * nonsingular by moving a variable from one of its bounds to the other
   * rather than freeing it. */
  intQPOASESFlipBounds ,
  /// number of iterative refinements of each solve (numRefinementSteps)
  intQPOASESRefineSteps ,
  /// first allowed new int parameter for derived classes
  intLastAlgParQPOASES
  };

 /// public enum for the double algorithmic parameters
 /** Public enum describing the different algorithmic parameters of double
  * type that QPOASESMILPSolver has in addition to those of MILPSolver. The
  * value dblLastAlgParQPOASES is provided so that the list can be easily
  * further extended by derived classes. Note that the tolerance of the
  * optimality of the solution (terminationTolerance) is the standard
  * dblFAccSol. */
 enum dbl_par_type_QPOASES {
  /// two sides closer than this are an equality (boundTolerance)
  dblQPOASESBoundTol = dblLastAlgParMILP ,
  /// relative regularization of H (epsRegularisation)
  dblQPOASESEpsReg ,
  /// initial value of the far bounds of the homotopy (initialFarBounds)
  dblQPOASESInitFarBounds ,
  /// factor by which the far bounds grow when reached (growFarBounds)
  dblQPOASESGrowFarBounds ,
  /// initial relaxation of the bounds of the homotopy (boundRelaxation)
  dblQPOASESBoundRelax ,
  /// first allowed new double parameter for derived classes
  dblLastAlgParQPOASES
  };

/** @} ---------------------------------------------------------------------*/
/*------------------- CONSTRUCTING AND DESTRUCTING -------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Constructing and destructing QPOASESMILPSolver
 * @{ */

 /// void constructor
 QPOASESMILPSolver( void );

/*--------------------------------------------------------------------------*/
 /// destructor: releases the problem of qpOASES

 ~QPOASESMILPSolver() override;

/** @} ---------------------------------------------------------------------*/
/*--------------------- DERIVED METHODS OF BASE CLASS ----------------------*/
/*--------------------------------------------------------------------------*/
/** @name Methods of MILPSolver that are redefined
 * @{ */

 /// clears the problem, together with that of qpOASES

 void clear_problem( unsigned int what ) override;

/*--------------------------------------------------------------------------*/
 /// loads the problem from the Block into the dense data of qpOASES

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
/** @name Handling the parameters of QPOASESMILPSolver
 *
 * Besides those of MILPSolver and those of qpOASES above, the standard
 * intMaxIter (the limit on the working set recalculations of qpOASES, if
 * not set 100 times the number of variables and constraints plus 1000),
 * dblMaxTime (its time limit, in seconds), dblFAccSol (its termination
 * tolerance) and intLogVerb (its printLevel: none if 0, low if 1, medium
 * if 2, high if more) are handled.
 * @{ */

 void set_par( idx_type par , int value ) override;

 void set_par( idx_type par , double value ) override;

 [[nodiscard]] idx_type get_num_int_par( void ) const override {
  return( MILPSolver::get_num_int_par()
	  + intLastAlgParQPOASES - intLastAlgParMILP );
  }

 [[nodiscard]] idx_type get_num_dbl_par( void ) const override {
  return( MILPSolver::get_num_dbl_par()
	  + dblLastAlgParQPOASES - dblLastAlgParMILP );
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

 /// solves the loaded problem, by a hot start if possible

 int guts_of_compute( void ) override;

 /// qpOASES cannot change the size of a problem: it is always loaded again

 bool reload_on_modification( void ) override { return( true ); }

/*--------------------------------------------------------------------------*/
/*---------------------------- PRIVATE PART --------------------------------*/
/*--------------------------------------------------------------------------*/

 private:

 /// solves the problem from scratch, in a new problem of qpOASES

 int solve_from_scratch( void );

 /// solves the problem by a hot start of the current problem of qpOASES

 int hot_start( void );

 /// translates the outcome of a solve into a return status of compute()

 int decode( int ret , int nwsr , int done ) const;

 /// brings the Options back to the preset, then applies those set

 void apply_preset( void );

 /// sets an Option of qpOASES from an int parameter, if it is one

 bool set_option( idx_type par , int value );

 /// sets an Option of qpOASES from a double parameter, if it is one

 bool set_option( idx_type par , double value );

 /// sets the printLevel of the Options from intLogVerb

 void set_print_level( void );

 /// the limit on the working set recalculations of a solve

 [[nodiscard]] int max_nwsr( void ) const;

 /// true if the last solve has given a primal solution

 [[nodiscard]] bool solved( void ) const { return( d_status == kOK ); }

/*--------------------------------------------------------------------------*/
/*--------------------------- PRIVATE FIELDS -------------------------------*/
/*--------------------------------------------------------------------------*/

 int d_n = 0;                    ///< number of variables
 int d_nc = 0;                   ///< number of rows
 double d_sign = 1;              ///< -1 if a maximization was negated

 // the problem as loaded from the Block

 std::vector< double > d_H;      ///< dense Hessian, row-major
 std::vector< double > d_g;      ///< linear term
 std::vector< double > d_A;      ///< dense rows, row-major
 std::vector< double > d_lb;     ///< lower bounds on the variables
 std::vector< double > d_ub;     ///< upper bounds on the variables
 std::vector< double > d_lbA;    ///< lower sides of the rows
 std::vector< double > d_ubA;    ///< upper sides of the rows

 // qpOASES works on the H and A it is given rather than on copies, and
 // regularizes H in place, hence it is given copies of them; a hot start
 // with new matrices still reads the old ones, hence there are two copies,
 // d_qpH[ d_cur ] and d_qpA[ d_cur ] being those of the current problem

 std::vector< double > d_qpH[ 2 ];  ///< copies of H for qpOASES
 std::vector< double > d_qpA[ 2 ];  ///< copies of A for qpOASES
 int d_cur = 0;                     ///< the copy qpOASES is working on
 std::vector< double > d_oldH;  ///< H of the current problem, as loaded
 std::vector< double > d_oldA;  ///< A of the current problem, as loaded

 std::unique_ptr< qpOASES::SQProblem > d_qp;   ///< the problem of qpOASES
 std::unique_ptr< qpOASES::Options > d_options;  ///< its Options
 bool d_qp_solved = false;       ///< true if d_qp has been solved
 bool d_new_options = false;     ///< true if d_options changed after it

 std::map< idx_type , int > d_int_set;     ///< int Options explicitly set
 std::map< idx_type , double > d_dbl_set;  ///< double Options explicitly set

 std::vector< double > d_x;      ///< primal solution
 std::vector< double > d_y;      ///< multipliers, bounds then rows
 double d_fval = 0;              ///< value of the last solution
 int d_status = kError;          ///< status of the last solve

 bool d_hot_start = true;        ///< intQPOASESHotStart
 int d_preset = 0;               ///< intQPOASESPreset
 int d_hessian_type = 0;         ///< intQPOASESHessianType
 int d_max_iter = Inf< int >();       ///< intMaxIter
 double d_max_time = Inf< double >(); ///< dblMaxTime

/*--------------------------------------------------------------------------*/

 SMSpp_insert_in_factory_h;

/*--------------------------------------------------------------------------*/

 };  // end( class( QPOASESMILPSolver ) )

/*--------------------------------------------------------------------------*/

}  // end( namespace SMSpp_di_unipi_it )

/*--------------------------------------------------------------------------*/

#endif  /* QPOASESMILPSolver.h included */

/*--------------------------------------------------------------------------*/
/*-------------------- End File QPOASESMILPSolver.h ------------------------*/
/*--------------------------------------------------------------------------*/
