/*--------------------------------------------------------------------------*/
/*------------------------ File ClarabelMILPSolver.h -----------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the ClarabelMILPSolver class, implementing the MILPSolver
 * interface for convex Quadratic Programs with linear constraints using the
 * interior-point solver Clarabel.
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

#ifndef __CLARABELMILPSOLVER_H
 #define __CLARABELMILPSOLVER_H
                      /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include <cstdint>
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
/*---------------------- CLASS ClarabelMILPSolver --------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------- GENERAL NOTES --------------------------------*/
/*--------------------------------------------------------------------------*/
/// class for solving convex QP Blocks via the interior-point solver Clarabel
/** The ClarabelMILPSolver implements the MILPSolver interface using
 * Clarabel, an open-source interior-point solver for convex conic programs
 * based on a homogeneous embedding, here only used for the convex Quadratic
 * Programs
 *
 * \f[
 *   \min \{ \tfrac{1}{2} x^\top P x + q^\top x \,:\,
 *           A x + s = b \,,\, s \in \mathcal{K} \} \,,
 * \f]
 *
 * where \f$ P \f$ is sparse and only needs to be positive semidefinite (a
 * singular or zero \f$ P \f$, e.g., that of a Linear Program, is allowed),
 * and \f$ \mathcal{K} \f$ is the product of the zero cone, which holds the
 * equality constraints and the fixed variables, and of the nonnegative
 * orthant, which holds each finite side of the other constraints and of
 * the bounds as a row of its own (an infinite side giving no row). It is
 * registered in the Solver factory, hence it can be selected by name
 * ("ClarabelMILPSolver") in a BlockSolverConfig.
 *
 * Because of the homogeneous embedding, an unbounded problem is detected as
 * such (Clarabel says it is dual infeasible) rather than being run until
 * the iteration limit, and a singular P needs no regularization of the
 * problem, which makes Clarabel a candidate for the master problems of
 * bundle methods, whose Hessian is often singular and which are unbounded
 * in the first iterations if vertical cuts are present.
 *
 * The following restrictions apply:
 *
 * - quadratic constraints are not supported, and integer variables are only
 *   accepted if relaxed via intRelaxIntVars;
 *
 * - Clarabel has no warm start, and changing the loaded problem in place
 *   only covers the values of a fixed structure: whenever a Modification
 *   is received the problem is loaded again whole from the Block [see
 *   MILPSolver::reload_on_modification()], and each solve starts from
 *   scratch;
 *
 * - the value of the solution is computed from it rather than taken from
 *   Clarabel, and both get_lb() and get_ub() return it;
 *
 * - no unbounded ray or infeasibility certificate is provided, i.e., the
 *   has_var_direction() / has_dual_direction() interface of CDASolver is
 *   not implemented. */

class ClarabelMILPSolver : public MILPSolver
{
/*--------------------------------------------------------------------------*/
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

 public:

/*--------------------------------------------------------------------------*/
/*---------------------------- PUBLIC TYPES --------------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Public types of ClarabelMILPSolver
 * @{ */

 /// public enum for the double algorithmic parameters
 /** Public enum describing the different algorithmic parameters of double
  * type that ClarabelMILPSolver has in addition to those of MILPSolver.
  * The value dblLastAlgParClarabel is provided so that the list can be
  * easily further extended by derived classes. Note that the tolerances on
  * the duality gap and on the residuals are the standard dblAbsAcc,
  * dblRelAcc and dblFAccSol. */
 enum dbl_par_type_Clarabel {
  /// absolute tolerance of the certificates of infeasibility
  dblClarabelInfeasAbsTol = dblLastAlgParMILP ,
  /// relative tolerance of the certificates of infeasibility
  dblClarabelInfeasRelTol ,
  /// first allowed new double parameter for derived classes
  dblLastAlgParClarabel
  };

/** @} ---------------------------------------------------------------------*/
/*------------------- CONSTRUCTING AND DESTRUCTING -------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Constructing and destructing ClarabelMILPSolver
 * @{ */

 /// void constructor
 ClarabelMILPSolver( void );

/*--------------------------------------------------------------------------*/
 /// destructor

 ~ClarabelMILPSolver() override = default;

/** @} ---------------------------------------------------------------------*/
/*--------------------- DERIVED METHODS OF BASE CLASS ----------------------*/
/*--------------------------------------------------------------------------*/
/** @name Methods of MILPSolver that are redefined
 * @{ */

 /// clears the problem, together with the data of Clarabel

 void clear_problem( unsigned int what ) override;

/*--------------------------------------------------------------------------*/
 /// loads the problem from the Block into the sparse data of Clarabel

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
/** @name Handling the parameters of ClarabelMILPSolver
 *
 * Besides those of MILPSolver and those of Clarabel above, the standard
 * intMaxIter (the iteration limit of Clarabel), intLogVerb (any positive
 * value has Clarabel print its log), dblMaxTime (its time limit, in
 * seconds), dblAbsAcc and dblRelAcc (its absolute and relative tolerances
 * on the duality gap, tol_gap_abs and tol_gap_rel) and dblFAccSol (its
 * relative tolerance on the primal and dual residuals, tol_feas) are
 * handled, their default being those of Clarabel.
 * @{ */

 void set_par( idx_type par , int value ) override;

 void set_par( idx_type par , double value ) override;

 [[nodiscard]] idx_type get_num_dbl_par( void ) const override {
  return( MILPSolver::get_num_dbl_par()
	  + dblLastAlgParClarabel - dblLastAlgParMILP );
  }

 [[nodiscard]] int get_dflt_int_par( idx_type par ) const override;

 [[nodiscard]] double get_dflt_dbl_par( idx_type par ) const override;

 [[nodiscard]] int get_int_par( idx_type par ) const override;

 [[nodiscard]] double get_dbl_par( idx_type par ) const override;

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

 /// sets up Clarabel on the loaded problem and solves it

 int guts_of_compute( void ) override;

 /// Clarabel does not change a problem in place: it is always loaded again

 bool reload_on_modification( void ) override { return( true ); }

/*--------------------------------------------------------------------------*/
/*---------------------------- PRIVATE PART --------------------------------*/
/*--------------------------------------------------------------------------*/

 private:

 /// translates a status of Clarabel into a return status of compute()

 static int decode_status( int status );

 /// true if the last solve has given a primal solution

 [[nodiscard]] bool solved( void ) const;

/*--------------------------------------------------------------------------*/
/*--------------------------- PRIVATE FIELDS -------------------------------*/
/*--------------------------------------------------------------------------*/

 // the problem in the format of Clarabel: P and A are in compressed sparse
 // column format, P only with its upper triangle, and the rows of A are
 // the d_meq ones of the zero cone followed by those of the nonnegative
 // orthant

 using csc_idx = std::uintptr_t;  ///< the type of the indices of Clarabel

 int d_n = 0;                    ///< number of variables
 int d_m = 0;                    ///< number of rows of A
 int d_meq = 0;                  ///< number of rows in the zero cone
 double d_sign = 1;              ///< -1 if a maximization was negated

 std::vector< csc_idx > d_Pbeg;  ///< column starts of P
 std::vector< csc_idx > d_Pind;  ///< row indices of P
 std::vector< double > d_Pval;   ///< values of P
 std::vector< double > d_q;      ///< linear term

 std::vector< csc_idx > d_Abeg;  ///< column starts of A
 std::vector< csc_idx > d_Aind;  ///< row indices of A
 std::vector< double > d_Aval;   ///< values of A
 std::vector< double > d_b;      ///< right-hand side of A x + s = b

 // the rows of A of each constraint: the bounds of the variables first,
 // then the rows of MILPSolver; an equality has its row in d_up, and -1
 // says that there is none

 std::vector< int > d_up;        ///< row of the upper side
 std::vector< int > d_lo;        ///< row of the lower side

 bool d_loaded = false;          ///< true if the problem has been loaded

 std::vector< double > d_x;      ///< primal solution
 std::vector< double > d_z;      ///< multipliers of the rows of A
 double d_fval = 0;              ///< value of the last solution
 int d_status = 0;               ///< the status of Clarabel of the last solve

 // the parameters of Clarabel

 int d_max_iter;                 ///< intMaxIter
 double d_max_time = Inf< double >(); ///< dblMaxTime
 double d_tol_gap_abs;           ///< dblAbsAcc
 double d_tol_gap_rel;           ///< dblRelAcc
 double d_tol_feas;              ///< dblFAccSol
 double d_tol_infeas_abs;        ///< dblClarabelInfeasAbsTol
 double d_tol_infeas_rel;        ///< dblClarabelInfeasRelTol

/*--------------------------------------------------------------------------*/

 SMSpp_insert_in_factory_h;

/*--------------------------------------------------------------------------*/

 };  // end( class( ClarabelMILPSolver ) )

/*--------------------------------------------------------------------------*/

}  // end( namespace SMSpp_di_unipi_it )

/*--------------------------------------------------------------------------*/

#endif  /* ClarabelMILPSolver.h included */

/*--------------------------------------------------------------------------*/
/*-------------------- End File ClarabelMILPSolver.h -----------------------*/
/*--------------------------------------------------------------------------*/
