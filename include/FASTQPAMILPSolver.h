/*--------------------------------------------------------------------------*/
/*------------------------ File FASTQPAMILPSolver.h ------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the FASTQPAMILPSolver class, implementing the MILPSolver
 * interface for strictly convex Mixed-Integer Quadratic Programs with linear
 * constraints by a depth-first branch-and-bound that fixes the integer
 * variables in a predetermined order and bounds each node by the dual of its
 * continuous relaxation, solved by a feasible active-set method.
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

#ifndef __FASTQPAMILPSOLVER_H
 #define __FASTQPAMILPSOLVER_H
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
/*---------------------- CLASS FASTQPAMILPSolver ---------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------- GENERAL NOTES --------------------------------*/
/*--------------------------------------------------------------------------*/
/// branch-and-bound for strictly convex MIQPs with a dual active-set method
/** The FASTQPAMILPSolver implements the MILPSolver interface for the
 * problems
 *
 * \f[
 *   \min \{ x^\top Q x + c^\top x \,:\, A x \leq b \,,\,
 *           x_i \in \mathbb{Z} \,,\, i \in I \}
 * \f]
 *
 * where \f$ Q \f$ is positive definite on all the variables, by the
 * branch-and-bound of "A Feasible Active Set Method with Reoptimization for
 * Convex Quadratic Mixed-Integer Programming", SIAM J. Optim. 26(3), 2016,
 * which extends that of "An effective branch-and-bound algorithm for convex
 * quadratic integer programming", Math. Program. 135, 2012, to linear
 * constraints. It is registered in the Solver factory, hence it can be
 * selected by name ("FASTQPAMILPSolver") in a BlockSolverConfig.
 *
 * The integer variables are fixed one at a time, always in the same order
 * (that of the columns, the integer ones first), so that all the data of a
 * node that do not depend on the values of the fixings are computed once
 * per depth before the search: the inverse of the Hessian of the free
 * variables (by Schur complements of the previous one), the direction
 * along which the continuous minimum moves when the next variable is fixed,
 * the Hessian of the dual of the relaxation and its largest eigenvalue.
 * A node fixes the next variable to every integer value in turn, by
 * increasing distance from its value in the optimum of the relaxation, and
 * since the value of the relaxation is convex in the fixed value, the first
 * child of each side whose bound reaches the incumbent closes that side.
 * The search is depth-first, and the only upper bounds are those of the
 * leaves, where all the integer variables are fixed.
 *
 * The relaxation of a node, \f$ \min \{ x^\top Q_l x + \bar{c}^\top x \,:\,
 * A_l x \leq \bar{b} \} \f$ over the free variables, is bounded by its
 * dual \f$ \max \{ f(\hat{x}) - \lambda^\top \tilde{Q}_l \lambda -
 * \tilde{c}^\top \lambda \,:\, \lambda \geq 0 \} \f$, with \f$ \tilde{Q}_l
 * = \frac{1}{4} A_l Q_l^{-1} A_l^\top \f$, \f$ \hat{x} \f$ the minimum of
 * the objective with no constraint and \f$ \tilde{c} = \bar{b} - A_l
 * \hat{x} \f$: it has as many variables as rows, whatever the depth, and
 * only their signs as constraints. It is solved by the feasible active-set
 * method FAST-QPA of the paper (an estimate of the active set by multiplier
 * functions, a truncated conjugate gradient on the free multipliers and a
 * projected Armijo line search), starting from the optimal multipliers of
 * the father and from the active set estimated there. Every iterate is
 * dual feasible, hence a lower bound, and the node is pruned as soon as one
 * of them reaches the incumbent. The gradient of the dual at the optimum is
 * the slack of the rows at the primal optimum \f$ \hat{x} - \frac{1}{2}
 * Q_l^{-1} A_l^\top \lambda \f$, which is how the solution of a leaf is
 * recovered; at the end the leaf of the incumbent of a problem with
 * continuous variables is solved again, from scratch and to a tighter
 * tolerance.
 *
 * The problem of the Block is put in the above form as follows: a
 * maximization is negated, a row with sense 'G' is negated, an equality or
 * a ranged row gives two rows, and a finite bound of a continuous variable
 * gives a row. The bounds of the integer variables only restrict the values
 * the children fix them to, unless intFQIntBoundsAsRows says otherwise.
 *
 * The following restrictions apply:
 *
 * - the quadratic part of the Objective must be positive definite on all
 *   the variables, continuous ones included (a variable appearing only
 *   linearly is not allowed): compute() returns kError otherwise, after a
 *   failed Cholesky factorization, and quadratic constraints are not
 *   supported (load_problem() throws);
 *
 * - with unbounded integer variables and no incumbent, the enumeration of
 *   a problem whose relaxation is feasible but which has no integer
 *   solution need not end: a finite cutoff [see dblUpCutOff] or finite
 *   bounds make it finite, the node and time limits stop it in any case;
 *
 * - an infeasible node has an unbounded dual, which only a cutoff stops: as
 *   long as there is none, a dual whose multipliers grow past 1e6 is checked
 *   for a certificate of the infeasibility of the rows (u >= 0 with u' A_l
 *   = 0 and u' b < 0), a dual whose value reaches 1e20 is taken as
 *   unbounded, and a side of a node whose first child has a dual that does
 *   not converge within 3 m + 20 iterations is postponed, to be explored
 *   from its node once the rest of the tree is (with an incumbent by then,
 *   in all likelihood);
 *
 * - no dual solution is provided, not even if there is no integer
 *   variable, and neither is any unbounded ray or infeasibility
 *   certificate;
 *
 * - whenever a Modification is received the problem is loaded again whole
 *   from the Block [see MILPSolver::reload_on_modification()], and the
 *   preprocessing is done again at the next compute(). */

class FASTQPAMILPSolver : public MILPSolver
{
/*--------------------------------------------------------------------------*/
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

 public:
/*--------------------------------------------------------------------------*/
/*---------------------------- PUBLIC TYPES --------------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Public types of FASTQPAMILPSolver
 * @{ */

 /// public enum for the integer algorithmic parameters
 /** Public enum describing the different algorithmic parameters of int type
  * that FASTQPAMILPSolver has in addition to those of MILPSolver. The value
  * intLastAlgParFQ is provided so that the list can be easily further
  * extended by derived classes. */
 enum int_par_type_FQ {
  /// whether a child starts from the optimal multipliers of its father
  /**< If 1 [default], the dual of a child starts from the optimal
   * multipliers of its father and from the active set estimated there (as
   * in section 3.3 of the paper); if 0, from the zero multipliers and the
   * active set estimated at them. */
  intFQWarmStart = intLastAlgParMILP ,
  /// whether a node is pruned as soon as a dual iterate reaches the cutoff
  /**< If 1 [default], the dual of a node stops, and the node is pruned, as
   * soon as an iterate gives a lower bound reaching the incumbent (as in
   * section 3.2 of the paper); if 0, the dual is solved to optimality and
   * then compared with the incumbent. */
  intFQEarlyPrune ,
  /// whether the finite bounds of the integer variables are rows
  /**< If 0 [default], the finite bounds of an integer variable only
   * restrict the values the children of a node fix it to, so that the
   * relaxations do not see them; if 1, they are also rows of A, which makes
   * the relaxations stronger and the duals larger. */
  intFQIntBoundsAsRows ,
  /// maximum number of iterations of FAST-QPA at a node
  /**< The default is 1000. A node that reaches it is not pruned unless its
   * lower bound says so, and the multipliers it has reached are used as
   * they are; the leaf of the incumbent is solved again anyway. */
  intFQMaxQPIter ,
  /// whether the leaf of the incumbent is solved again at the end
  /**< If 1 [default], at the end the continuous variables of the incumbent
   * are computed again by solving the dual of its leaf from scratch to the
   * tolerance dblFQPostOptTol, as in section 3.5 of the paper; ignored if
   * all the variables are integer. */
  intFQPostSolve ,
  /// first allowed new int parameter for derived classes
  intLastAlgParFQ
  };

 /// public enum for the double algorithmic parameters
 /** Public enum describing the different algorithmic parameters of double
  * type that FASTQPAMILPSolver has in addition to those of MILPSolver. The
  * value dblLastAlgParFQ is provided so that the list can be easily further
  * extended by derived classes. */
 enum dbl_par_type_FQ {
  /// tolerance of the optimality of the dual of a node
  /**< The dual of a node is optimal when the 2-norm of its projected
   * gradient is at most this; the default 1e-5 is the one of section 3.2 of
   * the paper. At a leaf the tolerance is also at most dblFAccSol, the
   * gradient of the dual being the slack of the rows. */
  dblFQOptTol = dblLastAlgParMILP ,
  /// tolerance of the dual of the leaf solved again at the end
  /**< The default is 1e-9 [see intFQPostSolve]. */
  dblFQPostOptTol ,
  /// the parameter of the active-set estimate, as a fraction of its limit
  /**< The estimate is \f$ A(\lambda) = \{ i \,:\, \lambda_i \leq
   * \varepsilon g_i(\lambda) \} \f$, g being the gradient of the dual,
   * with \f$ 0 < \varepsilon < 1 / ( 2 \lambda_{\max}( \tilde{Q}_l ) ) \f$
   * (assumption 1 of the paper); \f$ \varepsilon \f$ is this, which must be
   * in ( 0 , 1 ), times that limit. The paper gives no value; the default
   * is 0.5. */
  dblFQEpsFactor ,
  /// the reduction factor of the step of the line search, in ( 0 , 1 )
  /**< The \f$ \delta \f$ of algorithm 1 of the paper, which gives no
   * value; the default is 0.5. */
  dblFQDelta ,
  /// the sufficient decrease of the line search, in ( 0 , 1 / 2 )
  /**< The \f$ \gamma \f$ of algorithm 1 of the paper, which gives no
   * value; the default is 1e-4. */
  dblFQGamma ,
  /// the threshold of the curvature test of the conjugate gradient
  /**< A direction p with \f$ p^\top \tilde{Q}_l p \leq \eta \| p \|^2 \f$
   * truncates the direction of the conjugate gradient (algorithm 2 of the
   * paper); \f$ \eta \f$ is this times \f$ \lambda_{\max}( \tilde{Q}_l )
   * \f$. The paper gives no value; the default is 1e-9. */
  dblFQEta ,
  /// the relative tolerance of the Cholesky factorization of Q
  /**< Q is taken as not positive definite if a pivot of its Cholesky
   * factorization is at most this times the largest of its diagonal
   * elements. The default is 1e-12. */
  dblFQPDTol ,
  /// first allowed new double parameter for derived classes
  dblLastAlgParFQ
  };

/** @} ---------------------------------------------------------------------*/
/*------------------- CONSTRUCTING AND DESTRUCTING -------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Constructing and destructing FASTQPAMILPSolver
 * @{ */

 /// void constructor
 FASTQPAMILPSolver( void );

/*--------------------------------------------------------------------------*/
 /// destructor

 ~FASTQPAMILPSolver() override;

/** @} ---------------------------------------------------------------------*/
/*--------------------- DERIVED METHODS OF BASE CLASS ----------------------*/
/*--------------------------------------------------------------------------*/
/** @name Methods of MILPSolver that are redefined
 * @{ */

 /// clears the problem, together with the data of the search

 void clear_problem( unsigned int what ) override;

/*--------------------------------------------------------------------------*/
 /// loads the problem from the Block, to be preprocessed at the next solve

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

 /// no dual solution is ever available

 bool has_dual_solution( void ) override { return( false ); }

 /// returns the number of nodes of the last solve, the root included

 [[nodiscard]] int get_explored_nodes( void ) const override {
  return( int( f_nodes ) );
  }

 /// returns the number of nodes of the last solve, the root included

 [[nodiscard]] int get_nodes( void ) const override {
  return( int( f_nodes ) );
  }

 /// returns the number of iterations of FAST-QPA at the root

 [[nodiscard]] long get_root_qp_iterations( void ) const {
  return( f_root_it );
  }

 /// returns the number of iterations of FAST-QPA at all the nodes

 [[nodiscard]] long get_qp_iterations( void ) const { return( f_qp_it ); }

/** @} ---------------------------------------------------------------------*/
/*------------------ METHODS FOR HANDLING THE PARAMETERS -------------------*/
/*--------------------------------------------------------------------------*/
/** @name Handling the parameters of FASTQPAMILPSolver
 *
 * Besides those of MILPSolver and those above, the following standard
 * parameters are handled:
 *
 * - intMaxIter: the maximum number of nodes, the root included;
 *
 * - dblMaxTime: the time limit, in seconds;
 *
 * - dblAbsAcc and dblRelAcc: a node is pruned if its lower bound is at
 *   least the incumbent minus the larger of dblAbsAcc and dblRelAcc times
 *   the absolute value of the incumbent (at least 1); the default of
 *   dblAbsAcc is 1e-6, the absolute optimality tolerance of section 4 of
 *   the paper;
 *
 * - dblUpCutOff and dblLwCutOff: the former, for a minimization, is the
 *   initial incumbent, and the latter stops the search as soon as a
 *   solution at least as good is found; the roles are swapped for a
 *   maximization;
 *
 * - dblFAccSol: the maximum violation of a row, relative to the absolute
 *   value of its right-hand side (at least 1), of the solution of a leaf;
 *   the default is 1e-6, the relative feasibility tolerance of section 4
 *   of the paper;
 *
 * - intLogVerb: if positive, a summary of the solve is written in the log.
 * @{ */

 void set_par( idx_type par , int value ) override;

 void set_par( idx_type par , double value ) override;

 [[nodiscard]] idx_type get_num_int_par( void ) const override {
  return( MILPSolver::get_num_int_par() + intLastAlgParFQ -
	   intLastAlgParMILP );
  }

 [[nodiscard]] idx_type get_num_dbl_par( void ) const override {
  return( MILPSolver::get_num_dbl_par() + dblLastAlgParFQ -
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

 /// preprocesses the loaded problem, if needed, and runs the search

 int guts_of_compute( void ) override;

 /// the preprocessing depends on all the data: the problem is loaded again

 bool reload_on_modification( void ) override { return( true ); }

/*--------------------------------------------------------------------------*/
/*---------------------------- PRIVATE PART --------------------------------*/
/*--------------------------------------------------------------------------*/

 private:
/*--------------------------------------------------------------------------*/
/*--------------------------- PRIVATE METHODS ------------------------------*/
/*--------------------------------------------------------------------------*/

 /// puts the problem of MILPSolver in the form min x'Qx + c'x , Ax <= b

 void convert( void );

 /// computes the data of every depth; false if Q is not positive definite

 bool preprocess( void );

 /// solves the dual of a node at depth \p l by FAST-QPA
 /** Solves \f$ \min \{ q( \lambda ) = \lambda^\top \tilde{Q}_l \lambda +
  * \tilde{c}^\top \lambda \,:\, \lambda \geq 0 \} \f$ starting from \p lam
  * and, if \p use_act, from the active set \p act, stopping as soon as \f$
  * q( \lambda ) \leq \f$ \p qcut. On return \p lam, \p act, \p g (the
  * gradient) and \p qv (the value) are those of the last iterate. Returns
  * 0 if optimal within \p tol, 1 if stopped by \p qcut, 2 if stopped by the
  * iteration limit, 3 if the line search has failed and 4 if a multiplier
  * has gone past \p lam_big. */

 int fast_qpa( int l , const double * ct , double qcut , double tol ,
	       double * lam ,
	       char * act , bool use_act , double * g , double & qv ,
	       double lam_big );

 /// true if \p lam proves the relaxation of a node at depth \p l empty
 /** The node has right-hand sides \p ct of its dual and minimum with no row
  * \p xh; \p lam is taken as a certificate of the infeasibility of the
  * rows, i.e., u = lam / sum( lam ) has u' A_l = 0 and u' b < 0 within the
  * tolerances, b being the right-hand sides of the node. */

 bool infeasible( int l , const double * ct , const double * xh ,
		  const double * lam ) const;

 /// the gradient \p g and the value of the dual at \p x, given \p ct

 double eval_dual( int l , const double * ct , const double * x ,
		   double * g ) const;

 /// the value x'Qx + c'x of a point in the order of the internal columns

 double value_of( const std::vector< double > & x ) const;

 /// the tolerance of pruning against the incumbent \p ub

 double prune_eps( double ub ) const;

 /// the solution of the leaf of the integer values \p fix, from scratch
 /** Writes in \p x the solution of the leaf whose integer variables have
  * the values in \p fix, the dual starting from \p lam and solved to the
  * tolerance \p tol; returns false if the dual is not solved or its primal
  * solution is not feasible within dblFAccSol. */

 bool solve_leaf( const std::vector< double > & fix ,
		  std::vector< double > & lam , double tol ,
		  std::vector< double > & x );

 /// the largest eigenvalue of the symmetric \p k x \p k matrix \p M

 static double max_eigenvalue( std::vector< double > M , int k );

/*--------------------------------------------------------------------------*/
/*--------------------------- PRIVATE FIELDS -------------------------------*/
/*--------------------------------------------------------------------------*/

 // the parameters

 bool f_warm_start = true;               ///< intFQWarmStart
 bool f_early_prune = true;              ///< intFQEarlyPrune
 bool f_int_bounds_rows = false;         ///< intFQIntBoundsAsRows
 int f_max_qp_iter = 1000;               ///< intFQMaxQPIter
 int f_it_cap = Inf< int >();            ///< a tighter limit, if any
 bool f_post_solve = true;               ///< intFQPostSolve
 double f_opt_tol = 1e-5;                ///< dblFQOptTol
 double f_post_opt_tol = 1e-9;           ///< dblFQPostOptTol
 double f_eps_factor = 0.5;              ///< dblFQEpsFactor
 double f_delta = 0.5;                   ///< dblFQDelta
 double f_gamma = 1e-4;                  ///< dblFQGamma
 double f_eta = 1e-9;                    ///< dblFQEta
 double f_pd_tol = 1e-12;                ///< dblFQPDTol
 long f_max_nodes = Inf< int >();        ///< intMaxIter
 double f_max_time = Inf< double >();    ///< dblMaxTime
 double f_abs_acc = 1e-6;                ///< dblAbsAcc
 double f_rel_acc = 1e-6;                ///< dblRelAcc
 double f_up_cutoff = Inf< double >();   ///< dblUpCutOff
 double f_lw_cutoff = -Inf< double >();  ///< dblLwCutOff
 double f_feas_tol = 1e-6;               ///< dblFAccSol

 // the problem in the internal form: the integer columns first, minimize

 bool f_dirty = true;      ///< true if convert() is to be called
 bool f_prep = false;      ///< true if preprocess() has been called
 bool f_pd = true;         ///< false if Q is not positive definite
 bool f_bad_data = false;  ///< true if the data are infeasible

 int f_n = 0;        ///< number of variables
 int f_n1 = 0;       ///< number of integer variables
 int f_m = 0;        ///< number of rows
 double f_sign = 1;  ///< -1 if a maximization was negated

 std::vector< int > f_perm;    ///< the column of each internal variable
 std::vector< double > f_Q;    ///< the Hessian / 2, dense, row-major
 std::vector< double > f_c;    ///< the linear term
 std::vector< double > f_A;    ///< the rows, dense, row-major
 std::vector< double > f_b;    ///< the right-hand sides
 std::vector< double > f_bsc;  ///< max( 1 , | b | ), to scale violations
 std::vector< double > f_ilo;  ///< lower bounds of the integer variables
 std::vector< double > f_ihi;  ///< upper bounds of the integer variables

 // the data of each depth l: those of l < f_n1 are the ones of the
 // branching, the last of f_Qt , f_lmax and f_eps is the one of the leaves

 std::vector< std::vector< double > > f_z;   ///< direction of the minimum
 std::vector< std::vector< double > > f_Az;  ///< A_l z^l
 std::vector< double > f_s;                  ///< z^l' Q_l z^l
 std::vector< std::vector< double > > f_Qt;  ///< Hessian of the dual
 std::vector< double > f_lmax;  ///< largest eigenvalue of f_Qt[ l ]
 std::vector< double > f_eps;   ///< parameter of the active-set estimate

 std::vector< double > f_xroot;  ///< minimum of the objective with no row
 double f_froot = 0;             ///< its value
 std::vector< double > f_Qil;    ///< Q_l^{-1} of the leaves
 std::vector< double > f_Wl;     ///< - Q_l^{-1} A_l' / 2 of the leaves

 // the workspace of FAST-QPA

 std::vector< double > w_xt , w_gt , w_d , w_y , w_r , w_p , w_Qp , w_xb ,
   w_gb , w_xn ,
   w_gn , w_x2 , w_g2 , w_pl;
 std::vector< int > w_N;

 // the results of the last solve

 int f_status = kUnEval;          ///< what the last compute() returned
 bool f_has_sol = false;          ///< true if an incumbent is there
 double f_UB = Inf< double >();   ///< the internal value of the incumbent
 double f_LB = -Inf< double >();  ///< the internal lower bound
 std::vector< double > f_x;       ///< the incumbent, in the original order
 long f_nodes = 0;                ///< number of nodes
 long f_qp_it = 0;                ///< number of iterations of FAST-QPA
 long f_root_it = 0;              ///< number of them at the root

/*--------------------------------------------------------------------------*/

 SMSpp_insert_in_factory_h;

/*--------------------------------------------------------------------------*/

 };  // end( class( FASTQPAMILPSolver ) )

/*--------------------------------------------------------------------------*/

}  // end( namespace SMSpp_di_unipi_it )

/*--------------------------------------------------------------------------*/

#endif /* FASTQPAMILPSolver.h included */

/*--------------------------------------------------------------------------*/
/*---------------------- End File FASTQPAMILPSolver.h ----------------------*/
/*--------------------------------------------------------------------------*/
