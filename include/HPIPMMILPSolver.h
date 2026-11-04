/*--------------------------------------------------------------------------*/
/*------------------------- File HPIPMMILPSolver.h -------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the HPIPMMILPSolver class, implementing the MILPSolver
 * interface for convex Quadratic Programs with linear constraints using the
 * dense interior-point solver of HPIPM.
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

#ifndef __HPIPMMILPSOLVER_H
 #define __HPIPMMILPSOLVER_H
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
/*----------------------- CLASS HPIPMMILPSolver ----------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------- GENERAL NOTES --------------------------------*/
/*--------------------------------------------------------------------------*/
/// class for solving convex QP Blocks via the dense QP solver of HPIPM
/** The HPIPMMILPSolver implements the MILPSolver interface using the dense
 * QP solver of HPIPM, an open-source interior-point method built on the
 * linear algebra of BLASFEO for the small and dense problems of Model
 * Predictive Control, for convex Quadratic Programs
 *
 * \f[
 *   \min \{ \tfrac{1}{2} x^\top H x + h^\top x \,:\,
 *           l \leq x \leq u \,,\, l_A \leq A x \leq u_A \} \,,
 * \f]
 *
 * where \f$ H \f$ is dense and positive semidefinite. It is registered in
 * the Solver factory, hence it can be selected by name ("HPIPMMILPSolver")
 * in a BlockSolverConfig.
 *
 * The problem is given to HPIPM in its own form: the rows with equal sides
 * and the variables with equal bounds are its equality constraints, the
 * bounds of the other variables its box constraints and the other rows its
 * general constraints, an infinite side being masked out.
 *
 * The following restrictions apply:
 *
 * - quadratic constraints are not supported, and integer variables are only
 *   accepted if relaxed via intRelaxIntVars;
 *
 * - HPIPM has no infeasibility or unboundedness detection: an infeasible or
 *   unbounded problem ends with the iteration limit (kStopIter) or with a
 *   too short step (kError), except for inconsistent equality constraints,
 *   which give kInfeasible;
 *
 * - whenever a Modification is received the problem is loaded again whole
 *   from the Block [see MILPSolver::reload_on_modification()]; if
 *   intHPIPMWarmStart is set, the solve starts from the previous
 *   solution, on the variables that are still there, and one that does not
 *   end so is repeated from the initial point of HPIPM;
 *
 * - HPIPM has no time limit, hence dblMaxTime is not handled;
 *
 * - the value of the solution is computed from it, on \f$ H \f$ as loaded;
 *
 * - no unbounded ray or infeasibility certificate is provided, i.e., the
 *   has_var_direction() / has_dual_direction() interface of CDASolver is
 *   not implemented. */

class HPIPMMILPSolver : public MILPSolver
{
/*--------------------------------------------------------------------------*/
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

 public:

/*--------------------------------------------------------------------------*/
/*---------------------------- PUBLIC TYPES --------------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Public types of HPIPMMILPSolver
 * @{ */

 /// public enum for the integer algorithmic parameters
 /** Public enum describing the different algorithmic parameters of int type
  * that HPIPMMILPSolver has in addition to those of MILPSolver. The value
  * intLastAlgParHPIPM is provided so that the list can be easily further
  * extended by derived classes. */
 enum int_par_type_HPIPM {
  /// the mode of HPIPM, which sets its other settings
  /**< 0 is SPEED_ABS, 1 SPEED, 2 BALANCE and 3 [default] ROBUST, the modes
   * of HPIPM from the fastest to the most reliable one; each of them sets
   * the iteration limit (15, 15, 30 and 100) when intMaxIter is not set. */
  intHPIPMMode = intLastAlgParMILP ,
  /// whether a solve starts from the previous solution
  /**< If 1, a solve after the first starts from the primal solution of the
   * previous one, on the variables that are still there [see
   * MILPSolver::previous_column_index()], the new ones at 0; if 0
   * [default], each solve starts from the initial point of HPIPM. */
  intHPIPMWarmStart ,
  /// first allowed new int parameter for derived classes
  intLastAlgParHPIPM
  };

/** @} ---------------------------------------------------------------------*/
/*------------------- CONSTRUCTING AND DESTRUCTING -------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Constructing and destructing HPIPMMILPSolver
 * @{ */

 /// void constructor
 HPIPMMILPSolver( void );

/*--------------------------------------------------------------------------*/
 /// destructor: releases the memory of HPIPM

 ~HPIPMMILPSolver() override;

/** @} ---------------------------------------------------------------------*/
/*--------------------- DERIVED METHODS OF BASE CLASS ----------------------*/
/*--------------------------------------------------------------------------*/
/** @name Methods of MILPSolver that are redefined
 * @{ */

 /// clears the problem, together with the memory of HPIPM

 void clear_problem( unsigned int what ) override;

/*--------------------------------------------------------------------------*/
 /// loads the problem from the Block into dense data

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
/** @name Handling the parameters of HPIPMMILPSolver
 *
 * Besides those of MILPSolver and those of HPIPM above, the standard
 * intMaxIter (the limit on the iterations of HPIPM; that of intHPIPMMode if
 * not set) and dblFAccSol (the tolerance of HPIPM on the residuals of
 * stationarity, equality, inequality and complementarity) are handled.
 * @{ */

 void set_par( idx_type par , int value ) override;

 void set_par( idx_type par , double value ) override;

 [[nodiscard]] idx_type get_num_int_par( void ) const override {
  return( MILPSolver::get_num_int_par()
	  + intLastAlgParHPIPM - intLastAlgParMILP );
  }

 [[nodiscard]] int get_dflt_int_par( idx_type par ) const override;

 [[nodiscard]] double get_dflt_dbl_par( idx_type par ) const override;

 [[nodiscard]] int get_int_par( idx_type par ) const override;

 [[nodiscard]] double get_dbl_par( idx_type par ) const override;

 [[nodiscard]] idx_type int_par_str2idx( const std::string & name )
  const override;

 [[nodiscard]] const std::string & int_par_idx2str( idx_type idx )
  const override;

/** @} ---------------------------------------------------------------------*/
/*--------------------- PROTECTED PART OF THE CLASS ------------------------*/
/*--------------------------------------------------------------------------*/

 protected:

/*--------------------------------------------------------------------------*/
/*-------------------------- PROTECTED METHODS -----------------------------*/
/*--------------------------------------------------------------------------*/

 /// solves the loaded problem

 int guts_of_compute( void ) override;

 /// HPIPM has no problem to change: it is always loaded again

 bool reload_on_modification( void ) override { return( true ); }

/*--------------------------------------------------------------------------*/
/*---------------------------- PRIVATE PART --------------------------------*/
/*--------------------------------------------------------------------------*/

 private:

 /// true if the last solve has given a primal solution

 [[nodiscard]] bool solved( void ) const { return( d_status == kOK ); }

/*--------------------------------------------------------------------------*/
/*--------------------------- PRIVATE FIELDS -------------------------------*/
/*--------------------------------------------------------------------------*/

 int d_n = 0;                    ///< number of variables
 int d_nc = 0;                   ///< number of rows
 double d_sign = 1;              ///< -1 if a maximization was negated

 // the problem as loaded from the Block

 std::vector< double > d_H;      ///< dense Hessian
 std::vector< double > d_h;      ///< linear term
 std::vector< double > d_A;      ///< dense rows, row-major
 std::vector< double > d_lb;     ///< lower bounds on the variables
 std::vector< double > d_ub;     ///< upper bounds on the variables
 std::vector< double > d_lbA;    ///< lower sides of the rows
 std::vector< double > d_ubA;    ///< upper sides of the rows

 /// the data of HPIPM (its structures and their memory), defined in the
 /// .cpp, so that the headers of HPIPM and BLASFEO do not reach whoever
 /// includes this
 struct HPIPMData;
 std::unique_ptr< HPIPMData > d_hp;

 std::vector< double > d_x;      ///< primal solution
 std::vector< double > d_y;      ///< multipliers, bounds then rows
 double d_fval = 0;              ///< value of the last solution
 int d_status = kError;          ///< status of the last solve

 int d_mode = 3;                 ///< intHPIPMMode
 bool d_warm = false;            ///< intHPIPMWarmStart
 std::vector< double > d_guess_x;  ///< start of a warm-started solve
 bool d_has_guess = false;       ///< true if d_guess_x is a solution
 double d_tol;                   ///< dblFAccSol
 int d_max_iter = Inf< int >();       ///< intMaxIter
 double d_max_time = Inf< double >(); ///< dblMaxTime, not used

/*--------------------------------------------------------------------------*/

 SMSpp_insert_in_factory_h;

/*--------------------------------------------------------------------------*/

 };  // end( class( HPIPMMILPSolver ) )

/*--------------------------------------------------------------------------*/

}  // end( namespace SMSpp_di_unipi_it )

/*--------------------------------------------------------------------------*/

#endif  /* HPIPMMILPSolver.h included */

/*--------------------------------------------------------------------------*/
/*--------------------- End File HPIPMMILPSolver.h -------------------------*/
/*--------------------------------------------------------------------------*/
