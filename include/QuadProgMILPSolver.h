/*--------------------------------------------------------------------------*/
/*----------------------- File QuadProgMILPSolver.h ------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the QuadProgMILPSolver class, implementing the MILPSolver
 * interface for strictly convex Quadratic Programs with linear constraints
 * using the dual active-set solver quadprog.
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

#ifndef __QUADPROGMILPSOLVER_H
 #define __QUADPROGMILPSOLVER_H
                      /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

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
/*--------------------- CLASS QuadProgMILPSolver ---------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------- GENERAL NOTES --------------------------------*/
/*--------------------------------------------------------------------------*/
/// class for solving strictly convex QP Blocks via the solver quadprog
/** The QuadProgMILPSolver implements the MILPSolver interface using
 * quadprog, the implementation of the dual active-set method of Goldfarb and
 * Idnani for strictly convex Quadratic Programs
 *
 * \f[
 *   \min \{ \tfrac{1}{2} x^\top H x + h^\top x \,:\,
 *           l \leq x \leq u \,,\, l_A \leq A x \leq u_A \} \,,
 * \f]
 *
 * where \f$ H \f$ is dense and has to be positive definite, since the method
 * starts from the unconstrained minimum \f$ - H^{-1} h \f$. It is registered
 * in the Solver factory, hence it can be selected by name
 * ("QuadProgMILPSolver") in a BlockSolverConfig.
 *
 * quadprog only knows constraints of the form \f$ c^\top x = b \f$ and
 * \f$ c^\top x \geq b \f$, given as the dense columns of a matrix, hence each
 * finite bound of a variable and each finite side of a row is one column of
 * it, and an equality is one column only. The data are kept dense: the
 * method is meant for the small and dense problems, e.g., those of Model
 * Predictive Control, where it is among the fastest ones.
 *
 * Note that quadprog is distributed under the GNU General Public License
 * (version 2 or later), hence a library of MILPSolver built with it is
 * covered by that license too.
 *
 * The following restrictions apply:
 *
 * - quadratic constraints are not supported, and integer variables are only
 *   accepted if relaxed via intRelaxIntVars;
 *
 * - \f$ H \f$ has to be positive definite: if its Cholesky factorization
 *   fails, which is the case of a singular \f$ H \f$ (e.g., a Linear
 *   Program), the solve ends with kError, unless \f$ H \f$ is regularized
 *   (see dblQuadProgReg), in which case the regularized problem is solved;
 *   as a consequence, an unbounded problem is never told as such;
 *
 * - quadprog has no warm start of the active set: whenever a Modification is
 *   received the problem is loaded again whole from the Block [see
 *   MILPSolver::reload_on_modification()], and each solve starts from the
 *   unconstrained minimum, but if \f$ H \f$ has not changed the inverse of
 *   its Cholesky factor is not computed again (see intQuadProgReuseFactor);
 *
 * - quadprog has no limit on the iterations or on the time and no tolerance,
 *   hence intMaxIter, dblMaxTime and dblFAccSol are not handled;
 *
 * - linearly dependent equality constraints (e.g., the same row twice) may
 *   have the problem declared infeasible;
 *
 * - the value of the solution is computed from it, on the \f$ H \f$ of the
 *   Block, i.e., not regularized;
 *
 * - no unbounded ray or infeasibility certificate is provided, i.e., the
 *   has_var_direction() / has_dual_direction() interface of CDASolver is
 *   not implemented. */

class QuadProgMILPSolver : public MILPSolver
{
/*--------------------------------------------------------------------------*/
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

 public:

/*--------------------------------------------------------------------------*/
/*---------------------------- PUBLIC TYPES --------------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Public types of QuadProgMILPSolver
 * @{ */

 /// public enum for the integer algorithmic parameters
 /** Public enum describing the different algorithmic parameters of int type
  * that QuadProgMILPSolver has in addition to those of MILPSolver. The value
  * intLastAlgParQuadProg is provided so that the list can be easily further
  * extended by derived classes. */
 enum int_par_type_QuadProg {
  /// whether the inverse Cholesky factor of H is kept between solves
  /**< If 1 [default], the inverse of the Cholesky factor of H computed by a
   * solve is given to the next one, as long as H (and the regularization)
   * has not changed; if 0, each solve factorizes H again. */
  intQuadProgReuseFactor = intLastAlgParMILP ,
  /// first allowed new int parameter for derived classes
  intLastAlgParQuadProg
  };

 /// public enum for the double algorithmic parameters
 /** Public enum describing the different algorithmic parameters of double
  * type that QuadProgMILPSolver has in addition to those of MILPSolver. The
  * value dblLastAlgParQuadProg is provided so that the list can be easily
  * further extended by derived classes. */
 enum dbl_par_type_QuadProg {
  /// the regularization of H, added to its diagonal
  /**< If positive, this value is added to each diagonal entry of H before
   * its factorization, so that a positive semidefinite H becomes positive
   * definite; note that the regularized problem of an unbounded one has a
   * finite solution, far from the origin. 0 [default] means no
   * regularization. */
  dblQuadProgReg = dblLastAlgParMILP ,
  /// first allowed new double parameter for derived classes
  dblLastAlgParQuadProg
  };

/** @} ---------------------------------------------------------------------*/
/*------------------- CONSTRUCTING AND DESTRUCTING -------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Constructing and destructing QuadProgMILPSolver
 * @{ */

 /// void constructor
 QuadProgMILPSolver( void ) : MILPSolver() {}

/*--------------------------------------------------------------------------*/
 /// destructor

 ~QuadProgMILPSolver() override = default;

/** @} ---------------------------------------------------------------------*/
/*--------------------- DERIVED METHODS OF BASE CLASS ----------------------*/
/*--------------------------------------------------------------------------*/
/** @name Methods of MILPSolver that are redefined
 * @{ */

 /// clears the problem, together with the factor of H

 void clear_problem( unsigned int what ) override;

/*--------------------------------------------------------------------------*/
 /// loads the problem from the Block into the dense data of quadprog

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
/** @name Handling the parameters of QuadProgMILPSolver
 *
 * Besides those of MILPSolver, only those of quadprog above are handled.
 * @{ */

 void set_par( idx_type par , int value ) override;

 void set_par( idx_type par , double value ) override;

 [[nodiscard]] idx_type get_num_int_par( void ) const override {
  return( MILPSolver::get_num_int_par()
	  + intLastAlgParQuadProg - intLastAlgParMILP );
  }

 [[nodiscard]] idx_type get_num_dbl_par( void ) const override {
  return( MILPSolver::get_num_dbl_par()
	  + dblLastAlgParQuadProg - dblLastAlgParMILP );
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

 /// quadprog has no problem to change: it is always loaded again

 bool reload_on_modification( void ) override { return( true ); }

/*--------------------------------------------------------------------------*/
/*---------------------------- PRIVATE PART --------------------------------*/
/*--------------------------------------------------------------------------*/

 private:

 /// computes in d_Rinv the inverse Cholesky factor of H, regularized

 bool factorize( void );

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

 // the constraints c_k' x = b_k (the first d_meq ones) and c_k' x >= b_k of
 // quadprog: c_k is column k of d_C (d_n entries), d_idx[ k ] the index in
 // d_y of the bound or row it comes from, and d_sgn[ k ] is -1 if c_k is the
 // opposite of that row (an upper side), +1 otherwise

 std::vector< double > d_C;      ///< columns of the constraints
 std::vector< double > d_b;      ///< right-hand sides of the constraints
 std::vector< int > d_idx;       ///< bound or row of each constraint
 std::vector< double > d_sgn;    ///< +1 or -1 for each constraint
 int d_meq = 0;                  ///< number of equality constraints

 std::vector< double > d_Rinv;   ///< inverse Cholesky factor, column-major
 std::vector< double > d_factH;  ///< H whose factor is in d_Rinv
 double d_factReg = 0;           ///< the regularization of that factor
 bool d_factored = false;        ///< true if d_Rinv is the factor of d_factH

 std::vector< double > d_x;      ///< primal solution
 std::vector< double > d_y;      ///< multipliers, bounds then rows
 double d_fval = 0;              ///< value of the last solution
 int d_status = kError;          ///< status of the last solve

 bool d_reuse = true;            ///< intQuadProgReuseFactor
 double d_reg = 0;               ///< dblQuadProgReg

/*--------------------------------------------------------------------------*/

 SMSpp_insert_in_factory_h;

/*--------------------------------------------------------------------------*/

 };  // end( class( QuadProgMILPSolver ) )

/*--------------------------------------------------------------------------*/

}  // end( namespace SMSpp_di_unipi_it )

/*--------------------------------------------------------------------------*/

#endif  /* QuadProgMILPSolver.h included */

/*--------------------------------------------------------------------------*/
/*-------------------- End File QuadProgMILPSolver.h -----------------------*/
/*--------------------------------------------------------------------------*/
