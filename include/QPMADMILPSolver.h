/*--------------------------------------------------------------------------*/
/*------------------------- File QPMADMILPSolver.h -------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the QPMADMILPSolver class, implementing the MILPSolver
 * interface for strictly convex Quadratic Programs with linear constraints
 * using the dual active-set solver qpmad.
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

#ifndef __QPMADMILPSOLVER_H
 #define __QPMADMILPSOLVER_H
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
/*----------------------- CLASS QPMADMILPSolver ----------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------- GENERAL NOTES --------------------------------*/
/*--------------------------------------------------------------------------*/
/// class for solving strictly convex QP Blocks via the solver qpmad
/** The QPMADMILPSolver implements the MILPSolver interface using qpmad, an
 * open-source implementation of the dual active-set method of Goldfarb and
 * Idnani for strictly convex Quadratic Programs
 *
 * \f[
 *   \min \{ \tfrac{1}{2} x^\top H x + h^\top x \,:\,
 *           l \leq x \leq u \,,\, l_A \leq A x \leq u_A \} \,,
 * \f]
 *
 * where \f$ H \f$ is dense and has to be positive definite, since the method
 * starts from the unconstrained minimum \f$ - H^{-1} h \f$, which it computes
 * from the Cholesky factorization of \f$ H \f$. It is registered in the
 * Solver factory, hence it can be selected by name ("QPMADMILPSolver") in a
 * BlockSolverConfig.
 *
 * The data are kept dense: the method is meant for the small and dense
 * problems, e.g., those of Model Predictive Control, where it is among the
 * fastest ones.
 *
 * The following restrictions apply:
 *
 * - quadratic constraints are not supported, and integer variables are only
 *   accepted if relaxed via intRelaxIntVars;
 *
 * - \f$ H \f$ has to be positive definite: if its Cholesky factorization
 *   fails, which is the case of a singular \f$ H \f$ (e.g., a Linear
 *   Program), the solve ends with kError, unless \f$ H \f$ is regularized
 *   (see dblQPMADReg), in which case the regularized problem is solved; as a
 *   consequence, an unbounded problem is never told as such;
 *
 * - qpmad has no warm start of the active set: whenever a Modification is
 *   received the problem is loaded again whole from the Block [see
 *   MILPSolver::reload_on_modification()], and each solve starts from the
 *   unconstrained minimum, but if \f$ H \f$ has not changed its Cholesky
 *   factor is not computed again (see intQPMADReuseFactor);
 *
 * - qpmad has no time limit, hence dblMaxTime is not handled;
 *
 * - qpmad gives the multipliers of the active inequality constraints only,
 *   hence all the multipliers are computed from the active set it ends with
 *   as the least-squares solution of the stationarity condition;
 *
 * - the value of the solution is computed from it, qpmad giving none;
 *
 * - no unbounded ray or infeasibility certificate is provided, i.e., the
 *   has_var_direction() / has_dual_direction() interface of CDASolver is
 *   not implemented. */

class QPMADMILPSolver : public MILPSolver
{
/*--------------------------------------------------------------------------*/
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

 public:

/*--------------------------------------------------------------------------*/
/*---------------------------- PUBLIC TYPES --------------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Public types of QPMADMILPSolver
 * @{ */

 /// public enum for the integer algorithmic parameters
 /** Public enum describing the different algorithmic parameters of int type
  * that QPMADMILPSolver has in addition to those of MILPSolver. The value
  * intLastAlgParQPMAD is provided so that the list can be easily further
  * extended by derived classes. */
 enum int_par_type_QPMAD {
  /// whether the Cholesky factor of H is kept from a solve to the next
  /**< If 1 [default], the Cholesky factor of H computed by a solve is given
   * to the next one, as long as H (and the regularization) has not changed;
   * if 0, each solve factorizes H again. */
  intQPMADReuseFactor = intLastAlgParMILP ,
  /// first allowed new int parameter for derived classes
  intLastAlgParQPMAD
  };

 /// public enum for the double algorithmic parameters
 /** Public enum describing the different algorithmic parameters of double
  * type that QPMADMILPSolver has in addition to those of MILPSolver. The
  * value dblLastAlgParQPMAD is provided so that the list can be easily
  * further extended by derived classes. Note that the tolerance of qpmad,
  * on the violation of the constraints and on the linear dependence of the
  * active ones, is the standard dblFAccSol. */
 enum dbl_par_type_QPMAD {
  /// the regularization of H, added to its diagonal
  /**< If positive, this value is added to each diagonal entry of H before
   * its factorization, so that a positive semidefinite H becomes positive
   * definite; note that the regularized problem of an unbounded one has a
   * finite solution, far from the origin. 0 [default] means no
   * regularization. */
  dblQPMADReg = dblLastAlgParMILP ,
  /// first allowed new double parameter for derived classes
  dblLastAlgParQPMAD
  };

/** @} ---------------------------------------------------------------------*/
/*------------------- CONSTRUCTING AND DESTRUCTING -------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Constructing and destructing QPMADMILPSolver
 * @{ */

 /// void constructor
 QPMADMILPSolver( void );

/*--------------------------------------------------------------------------*/
 /// destructor: releases the solver of qpmad

 ~QPMADMILPSolver() override;

/** @} ---------------------------------------------------------------------*/
/*--------------------- DERIVED METHODS OF BASE CLASS ----------------------*/
/*--------------------------------------------------------------------------*/
/** @name Methods of MILPSolver that are redefined
 * @{ */

 /// clears the problem, together with the factor of H

 void clear_problem( unsigned int what ) override;

/*--------------------------------------------------------------------------*/
 /// loads the problem from the Block into the dense data of qpmad

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
/** @name Handling the parameters of QPMADMILPSolver
 *
 * Besides those of MILPSolver and those of qpmad above, the standard
 * intMaxIter (the limit on the iterations of qpmad, i.e., on the changes of
 * its active set; none if not set) and dblFAccSol (its tolerance) are
 * handled.
 * @{ */

 void set_par( idx_type par , int value ) override;

 void set_par( idx_type par , double value ) override;

 [[nodiscard]] idx_type get_num_int_par( void ) const override {
  return( MILPSolver::get_num_int_par()
	  + intLastAlgParQPMAD - intLastAlgParMILP );
  }

 [[nodiscard]] idx_type get_num_dbl_par( void ) const override {
  return( MILPSolver::get_num_dbl_par()
	  + dblLastAlgParQPMAD - dblLastAlgParMILP );
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

 /// solves the loaded problem

 int guts_of_compute( void ) override;

 /// qpmad has no problem to change: it is always loaded again

 bool reload_on_modification( void ) override { return( true ); }

/*--------------------------------------------------------------------------*/
/*---------------------------- PRIVATE PART --------------------------------*/
/*--------------------------------------------------------------------------*/

 private:

 /// computes all the multipliers from the active set of the last solve

 void compute_multipliers( void );

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
 std::vector< double > d_h;      ///< linear term
 std::vector< double > d_A;      ///< dense rows, row-major
 std::vector< double > d_lb;     ///< lower bounds on the variables
 std::vector< double > d_ub;     ///< upper bounds on the variables
 std::vector< double > d_lbA;    ///< lower sides of the rows
 std::vector< double > d_ubA;    ///< upper sides of the rows

 /// the data of qpmad (its solver and the factor of H), defined in the
 /// .cpp, so that the headers of qpmad do not reach whoever includes this
 struct QPMADData;
 std::unique_ptr< QPMADData > d_qp;

 std::vector< double > d_factH;  ///< H whose factor is in d_qp, as loaded
 double d_factReg = 0;           ///< the regularization of that factor

 std::vector< double > d_x;      ///< primal solution
 std::vector< double > d_y;      ///< multipliers, bounds then rows
 double d_fval = 0;              ///< value of the last solution
 int d_status = kError;          ///< status of the last solve

 bool d_reuse = true;            ///< intQPMADReuseFactor
 double d_reg = 0;               ///< dblQPMADReg
 double d_tol;                   ///< dblFAccSol
 int d_max_iter = Inf< int >();       ///< intMaxIter
 double d_max_time = Inf< double >(); ///< dblMaxTime, not used

/*--------------------------------------------------------------------------*/

 SMSpp_insert_in_factory_h;

/*--------------------------------------------------------------------------*/

 };  // end( class( QPMADMILPSolver ) )

/*--------------------------------------------------------------------------*/

}  // end( namespace SMSpp_di_unipi_it )

/*--------------------------------------------------------------------------*/

#endif  /* QPMADMILPSolver.h included */

/*--------------------------------------------------------------------------*/
/*--------------------- End File QPMADMILPSolver.h -------------------------*/
/*--------------------------------------------------------------------------*/
