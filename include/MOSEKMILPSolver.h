/*--------------------------------------------------------------------------*/
/*------------------------- File MOSEKMILPSolver.h -------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the MOSEKMILPSolver class, which implements the
 * MILPSolver interface for Mixed-Integer Linear and convex Quadratic
 * Programs using the Optimizer API for C of MOSEK. Only supports MOSEK
 * versions >= 11.0.1.
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

#ifndef __MOSEKMILPSOLVER_H
 #define __MOSEKMILPSOLVER_H
                      /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include <array>

#include <mutex>

#include <mosek.h>

#include "MILPSolver.h"

// Include the proper MOSEK parameter mapping
#include <boost/preprocessor/cat.hpp>
#include <boost/preprocessor/stringize.hpp>
#include BOOST_PP_STRINGIZE( BOOST_PP_CAT( BOOST_PP_CAT( BOOST_PP_CAT( BOOST_PP_CAT( MOSEK , MSK_VERSION_MAJOR ) , MSK_VERSION_MINOR ) , MSK_VERSION_REVISION ) , _defs.h ) )

/*--------------------------------------------------------------------------*/
/*----------------------------- NAMESPACE ----------------------------------*/
/*--------------------------------------------------------------------------*/

/// namespace for the Structured Modeling System++ (SMS++)
namespace SMSpp_di_unipi_it {

 class LinearFunction;  // forward declaration of LinearFunction

/*--------------------------------------------------------------------------*/
/*----------------------- CLASS MOSEKMILPSolver ----------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------- GENERAL NOTES --------------------------------*/
/*--------------------------------------------------------------------------*/
/// class for solving MILP and convex MIQP problems via MOSEK
/** The MOSEKMILPSolver class derives from MILPSolver and extends the base
 * class to solve, using the Optimizer API for C of MOSEK, the problems with
 * linear constraints and a linear or convex quadratic objective, with or
 * without integer Variable, i.e., LP, QP, MILP and MIQP.
 *
 * The MOSEKMILPSolver can be registered to any kind of Block (assuming that
 * it contains such a formulation): it uses the base class to build the
 * matricial representation of the problem, which it gives to a task of
 * MOSEK, and then it implements the interface that MILPSolver provides for
 * processing the Modification on the task itself. The bounds of the
 * columns, the sides of the rows, the linear costs, the coefficients of the
 * matrix, the integrality of the columns, the addition and the removal of
 * rows and columns and the terms of the quadratic objective are all changed
 * in place, with the call of MOSEK that does it; a task changed in place
 * keeps its last basic solution, which the simplex algorithms of MOSEK use
 * to start the next solve from (MOSEK solves an LP with its interior point
 * algorithm unless MSK_IPAR_OPTIMIZER says otherwise, and that one does not
 * start from a basis).
 *
 * MOSEK has no range sense to emulate: each row and each column has a
 * bound key and two sides, hence a ranged row, a free row and an empty
 * range are given to it as they are. A quadratic objective
 * \f$ c^\top x + \frac{1}{2} x^\top Q x \f$ is given to MOSEK as the lower
 * triangle of \f$ Q \f$, which MOSEK only takes as a whole: the
 * Modification change a copy of it, which is given to MOSEK before the next
 * solve. A non-convex \f$ Q \f$ makes the solve fail, MOSEK being a convex
 * solver.
 *
 * Constraints with a quadratic Function are not supported, and neither are
 * semi-continuous or semi-integer columns: load_problem() throws
 * std::invalid_argument if the Block has any. MOSEK has no callback through
 * which user cuts or lazy constraints can be given to the branch-and-bound,
 * hence a MIP solve with intCutSepPar asking for separation throws
 * std::invalid_argument as well; with intRelaxIntVars == 2 the separation
 * loop of MILPSolver::compute() works as for any other :MILPSolver, since
 * it only needs LP solves.
 *
 * A solve that ends with an infeasible or an unbounded problem leaves, if
 * the problem is continuous, the certificate MOSEK computes along with the
 * status, which get_dual_direction() and get_var_direction() write in the
 * Block: no second solve is needed to have it.
 *
 * The interior point of MOSEK may stop on slow progress (MSK_RES_TRM_STALL)
 * when what it has is already an answer, typically the certificate of an
 * infeasible or an unbounded problem, which is then taken as such; a stall
 * with no certificate and no solution is kError. On the master problems of
 * a bundle that are unbounded, a stall of the second kind is avoided by
 * relaxing the tolerance of the certificates of a QP
 * (MSK_DPAR_INTPNT_QO_TOL_INFEAS, e.g., 1e-8).
 *
 * Besides the parameters of MILPSolver, the user can include in the
 * configuration all the parameters of MOSEK, with their symbolic names (say,
 * MSK_IPAR_OPTIMIZER or MSK_DPAR_INTPNT_TOL_REL_GAP) as the names of the
 * int, double and string parameters of SMS++: they are listed by a tool
 * of MILPSolver/tools out of the version of MOSEK being compiled against
 * [see set_par( idx_type , int )]. The parameters of Solver that MOSEK
 * has a counterpart of are mapped onto it:
 *
 * - intMaxIter    = MSK_IPAR_MIO_MAX_NUM_BRANCHES
 * - intMaxSol     = MSK_IPAR_MIO_MAX_NUM_SOLUTIONS
 * - intMaxThread  = MSK_IPAR_NUM_THREADS
 * - intLogVerb    = MSK_IPAR_LOG, the log being printed on std::cout only
 *                   if it is > 0, or on the file of strLogFileName
 * - dblMaxTime    = MSK_DPAR_OPTIMIZER_MAX_TIME
 * - dblRelAcc     = MSK_DPAR_MIO_TOL_REL_GAP
 * - dblAbsAcc     = MSK_DPAR_MIO_TOL_ABS_GAP
 * - dblFAccSol    = the primal feasibility tolerances of all the algorithms
 *                   of MOSEK, MSK_DPAR_INTPNT_TOL_PFEAS,
 *                   MSK_DPAR_INTPNT_QO_TOL_PFEAS, MSK_DPAR_BASIS_TOL_X and
 *                   MSK_DPAR_MIO_TOL_FEAS, the first giving the value read
 *
 * while dblUpCutOff and dblLwCutOff stop a MIP solve from the callback of
 * MOSEK as soon as the bound or the incumbent of the search reaches them
 * [see MILPSolver::compute()].
 *
 * The environment of MOSEK, which holds the license once a solve has
 * checked it out, is one for the whole process, shared by all the
 * MOSEKMILPSolver, each of which has a task of its own. */

class MOSEKMILPSolver : public MILPSolver {

/*--------------------------------------------------------------------------*/
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

 public:

/*--------------------------------------------------------------------------*/
/*---------------------------- PUBLIC TYPES --------------------------------*/
/*--------------------------------------------------------------------------*/
/** @name public types of MOSEKMILPSolver
 * @{ */

 /// enum for integer parameters
 enum int_par_type_MOSEKS {
  intFirstMOSEKPar = intLastAlgParMILP ,  ///< first MOSEK int parameter
  /// first allowed new int parameter for derived classes
  intLastAlgParMOSEKS = intFirstMOSEKPar + MOSEK_NUM_INT_PARS
  };

 /// enum for double parameters
 enum dbl_par_type_MOSEKS {
  /// first MOSEK double parameter
  dblFirstMOSEKPar = dblLastAlgParMILP ,
  /// first allowed new double parameter for derived classes
  dblLastAlgParMOSEKS = dblFirstMOSEKPar + MOSEK_NUM_DBL_PARS
  };

 /// enum for string parameters
 enum str_par_type_MOSEKS {
  /// first MOSEK string parameter
  strFirstMOSEKPar = strLastAlgParMILP ,
  /// first allowed new string parameter for derived classes
  strLastAlgParMOSEKS = strFirstMOSEKPar + MOSEK_NUM_STR_PARS
  };

/*--------------------------------------------------------------------------*/
 // "importing" a few types from Block

 using Subset = Block::Subset;

 using c_Subset = Block::c_Subset;

 using Range = Block::Range;

/** @} ---------------------------------------------------------------------*/
/*--------------------- CONSTRUCTOR AND DESTRUCTOR -------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Constructor and Destructor
 * @{ */

 MOSEKMILPSolver( void );

 ~MOSEKMILPSolver() override;

/** @} ---------------------------------------------------------------------*/
/*--------------------- DERIVED METHODS OF BASE CLASS ----------------------*/
/*--------------------------------------------------------------------------*/
/** @name Public Methods derived from base classes
 * @{ */

 /// returns a valid lower bound on the optimal objective function value
 OFValue get_lb( void ) override;

 /// returns a valid upper bound on the optimal objective function value
 OFValue get_ub( void ) override;

 /// returns the value of the current solution, if any
 OFValue get_var_value( void ) override;

 /// tells whether a solution is available
 bool has_var_solution( void ) override;

 /// writes the current solution in the Block
 void get_var_solution( Configuration * solc = nullptr ) override;

 /// tells whether an unbounded direction is available
 /** True if the last solve has found a continuous problem unbounded, and
  * MOSEK has given the primal ray certifying it. */
 bool has_var_direction( void ) override;

 /// writes the current unbounded direction in the Block
 void get_var_direction( Configuration * dirc = nullptr ) override;

 /// tells whether a dual solution is available
 /** True if the problem is continuous and the solution of MOSEK has a dual
  * part that MOSEK says feasible. */
 bool has_dual_solution( void ) override;

 /// tells whether the current dual solution is feasible
 bool is_dual_feasible( void ) override;

 /// writes the current dual solution in the Block
 /** The dual value of a row is the y of MOSEK and the reduced cost of a
  * column is its \f$ s_l^x - s_u^x \f$, which is \f$ c - A^\top y \f$ (plus
  * \f$ Q x \f$ on a QP) for both a minimization and a maximization, i.e.,
  * the convention of the other :MILPSolver [see
  * MILPSolver::write_dual_solution()]. */
 void get_dual_solution( Configuration * solc = nullptr ) override;

 /// tells whether a dual unbounded direction is available
 /** True if the last solve has found a continuous problem infeasible, and
  * MOSEK has given the certificate of it. */
 bool has_dual_direction( void ) override;

 /// writes the current dual unbounded direction in the Block
 /** The multipliers of the rows are the y of the certificate of MOSEK,
  * which has \f$ y_i \geq 0 \f$ on a row whose lower side is used and
  * \f$ y_i \leq 0 \f$ on one whose upper side is, as the Farkas
  * certificates of the other :MILPSolver; those of the columns are
  * \f$ c - A^\top y \f$, or \f$ - A^\top y \f$ with
  * intHomogeneousDirection [see MILPSolver]. The value of the
  * certificate, \f$ (l^c)^\top s_l^c - (u^c)^\top s_u^c + (l^x)^\top s_l^x
  * - (u^x)^\top s_u^x \f$ over the finite sides, is recorded as well [see
  * MILPSolver::get_dual_direction_value()]. */
 void get_dual_direction( Configuration * dirc = nullptr ) override;

 /// writes the problem on the specified file
 /** The format is the one MOSEK associates with the extension of the file
  * (.lp, .mps, .opf, .ptf, .task, ...); a name without an extension is
  * written as an LP. A format that cannot hold a free or a ranged row, as
  * the LP one, gets the model without the free rows, which bind nothing,
  * and with each ranged row split in its two sides. */
 void write_lp( const std::string & filename ) override;

 /// returns the number of nodes used to solve a MIP
 [[nodiscard]] int get_nodes( void ) const override;

 /// clears the problem, both in MILPSolver and in the task of MOSEK
 void clear_problem( unsigned int what ) override;

 /// loads the problem into the task of MOSEK
 void load_problem( void ) override;

 /// returns the number of nodes used to solve a MIP
 [[nodiscard]] int get_explored_nodes( void ) const override;

 /// returns the number of nodes of the search not yet explored
 [[nodiscard]] long get_left_nodes( void ) const override;

 /// returns true if a feasible solution is known
 [[nodiscard]] bool has_feasible_sol( void ) override;

 /// returns the time (in seconds) taken by the last solve
 [[nodiscard]] double get_runtime( void ) const override;

 /// gives MOSEK a starting solution of a MIP
 /** MOSEK takes one starting solution: only the first one is used, and a
  * warning is printed if there are more. A partial solution is completed
  * with the values of the last integer solution, or with 0 if there is
  * none, and MOSEK is asked to construct a solution out of the values of
  * the integer columns (MSK_IPAR_MIO_CONSTRUCT_SOL). */
 void add_mip_starts( std::vector< std::vector< int > > varidxs ,
                      std::vector< std::vector< double > > varvalues )
  override;

 #ifdef MILPSolver_DEBUG
  /// check the dictionaries for inconsistencies
  void check_status( void ) override;
 #endif

/** @} ---------------------------------------------------------------------*/
/*------------------- METHODS FOR HANDLING THE PARAMETERS ------------------*/
/*--------------------------------------------------------------------------*/
/** @name Methods for handling parameters
 * @{ */

 /// sets an integer parameter with the given value
 /** Besides those of MILPSolver and those of Solver that are mapped onto a
  * parameter of MOSEK [see the class], every integer parameter of MOSEK
  * can be set with its symbolic name, e.g., MSK_IPAR_OPTIMIZER with the
  * value of the MSKoptimizertypee it has to take in the version of MOSEK
  * of the build, since the values of an enumeration may change from a
  * version to the other (MSK_OPTIMIZER_INTPNT is 4 in MOSEK 10 and 11, 6
  * in MOSEK 12). A value MOSEK refuses throws std::invalid_argument.
  *
  * intCutSepPar can be set, but a MIP solve asking for the separation of
  * user cuts or lazy constraints throws [see the class]. */
 void set_par( idx_type par , int value ) override;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// sets a double parameter with the given value
 void set_par( idx_type par , double value ) override;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// sets a string parameter with the given value
 void set_par( idx_type par , std::string && value ) override;

/*--------------------------------------------------------------------------*/
 /// returns the number of integer parameters
 [[nodiscard]] idx_type get_num_int_par( void ) const override;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// returns the number of double parameters
 [[nodiscard]] idx_type get_num_dbl_par( void ) const override;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// returns the number of string parameters
 [[nodiscard]] idx_type get_num_str_par( void ) const override;

/*--------------------------------------------------------------------------*/
 /// returns the default value of the specified integer parameter
 [[nodiscard]] int get_dflt_int_par( idx_type par ) const override;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// returns the default value of the specified double parameter
 [[nodiscard]] double get_dflt_dbl_par( idx_type par ) const override;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /** Returns the default value of the specified string parameter
  * @note
  * Due to a limit in the implementation, the string referenced by
  * the return value is *overwritten* each time the method is called with
  * par as a MOSEK parameter. */
 [[nodiscard]] const std::string & get_dflt_str_par( idx_type par )
  const override;

/*--------------------------------------------------------------------------*/
 /// returns the value of the specified integer parameter
 [[nodiscard]] int get_int_par( idx_type par ) const override;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// returns the value of the specified double parameter
 [[nodiscard]] double get_dbl_par( idx_type par ) const override;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /** returns the value of the specified string parameter
  * @note
  * Due to a limit in the implementation, the string referenced by
  * the return value is *overwritten* each time the method is called with
  * par as a MOSEK parameter. */
 [[nodiscard]] const std::string & get_str_par( idx_type par ) const override;

/*--------------------------------------------------------------------------*/
 /// returns the index of the int parameter with the specified name
 [[nodiscard]] idx_type int_par_str2idx( const std::string & name )
  const override;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// returns the name of the int parameter with the specified index
 [[nodiscard]] const std::string & int_par_idx2str( idx_type idx )
  const override;

/*--------------------------------------------------------------------------*/
 /// returns the index of the double parameter with the specified name
 [[nodiscard]] idx_type dbl_par_str2idx( const std::string & name )
  const override;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// returns the name of the double parameter with the specified index
 [[nodiscard]] const std::string & dbl_par_idx2str( idx_type idx )
  const override;

/*--------------------------------------------------------------------------*/
 /// returns the index of the string parameter with the specified name
 [[nodiscard]] idx_type str_par_str2idx( const std::string & name )
  const override;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// returns the name of the string parameter with the specified index
 [[nodiscard]] const std::string & str_par_idx2str( idx_type idx )
  const override;

/*--------------------------------------------------------------------------*/

 double up_cut_off( void ) const { return( UpCutOff ); }

 double lw_cut_off( void ) const { return( LwCutOff ); }

/*--------------------------------------------------------------------------*/
 /// the callback of MOSEK, implemented as a method of the class
 /** Called by MOSEK during a MIP solve, only when a cutoff is set: it stops
  * the solve as soon as the bound of the search reaches dblUpCutOff (when
  * minimizing; dblLwCutOff when maximizing) or the incumbent reaches
  * dblLwCutOff (dblUpCutOff), returning nonzero. */
 int callback( MSKcallbackcodee caller , const MSKrealt * douinf ,
               const MSKint32t * intinf );

/** @} ---------------------------------------------------------------------*/
/*--------------------- PROTECTED PART OF THE CLASS ------------------------*/
/*--------------------------------------------------------------------------*/

 protected:

/*--------------------------------------------------------------------------*/
/*-------------------- PROTECTED METHODS OF THE CLASS ----------------------*/
/*--------------------------------------------------------------------------*/

 /// MOSEK back-end solve, called by MILPSolver::compute()
 int guts_of_compute( void ) override;

/*--------------------------------------------------------------------------*/
/** @name Methods for modifying the problem in the task of MOSEK
 *
 *  These methods implement the "Modification interface" of MILPSolver, so
 *  that all Modification are applied to the task of MOSEK.
 *  @{ */

 /// handles a Variable Modification
 void var_modification( const VariableMod * mod ) override;

 /// handles an Objective Modification
 void objective_modification( const ObjectiveMod * mod ) override;

 /// handles a Constraint Modification
 void const_modification( const ConstraintMod * mod ) override;

 /// handles a bound (OneVarConstraint) Modification
 void bound_modification( const OneVarConstraintMod * mod ) override;

 /// handles a Function Modification applied to the Objective
 void objective_function_modification( const FunctionMod * mod ) override;

 /// handles a Function Modification applied to a Constraint
 void constraint_function_modification( const FunctionMod * mod ) override;

 /// handles a Function Variable Modification applied to the Objective
 void objective_fvars_modification( const FunctionModVars * mod )
  override;

 /// handles a Function Variable Modification applied to a Constraint
 void constraint_fvars_modification( const FunctionModVars * mod )
  override;

 /// adds the whole column of each of the given Variable in one operation
 /** One MSK_appendvars() for all the columns, one MSK_putaijlist() for
  * their entries and one MSK_putclist() for their costs, the entries being
  * read from the sub-Modification of the group [see
  * MILPSolver::add_columns()]. */
 bool add_columns( const std::vector< Variable * > & vars ,
                   const GroupModification * gmod ) override;

 /// removes the whole column of each of the given Variable
 bool remove_columns( const std::vector< Variable * > & vars ,
                      const GroupModification * gmod ) override;

 /// changes the linear coefficients of a whole group in one operation
 /** One MSK_putaijlist() for the entries of the rows and one MSK_putclist()
  * for the costs, the old value of each cost being read once for the whole
  * group [see MILPSolver::change_coefficients()]. */
 bool change_coefficients(
               const std::vector< const FunctionMod * > & mods ) override;

 /// writes the bounds of a whole group in one MSK_putvarboundlist()
 bool change_bounds(
        const std::vector< const OneVarConstraintMod * > & mods ) override;

 /// one MSK_putvarboundlist() for a whole set of fixings
 /** A batch carrying a change of integrality is refused [see
  * MILPSolver::change_variables()]. */
 bool change_variables(
             const std::vector< const VariableMod * > & mods ) override;

 /// writes the sides of a whole group of rows in one MSK_putconboundlist()
 bool change_sides(
        const std::vector< const RowConstraintMod * > & mods ) override;

 /// adds a single new dynamic FRowConstraint
 void add_dynamic_constraint( const FRowConstraint * con ) override;

 /// adds a sequence of new dynamic FRowConstraint in one operation
 void add_dynamic_constraints(
          const std::vector< const FRowConstraint * > & cons ) override;

 /// adds a single new dynamic bound (OneVarConstraint)
 void add_dynamic_bound( const OneVarConstraint * con ) override;

 /// adds a single new dynamic ColVariable
 void add_dynamic_variable( const ColVariable * var ) override;

 /// removes a single dynamic FRowConstraint
 void remove_dynamic_constraint( const FRowConstraint * con ) override;

 /// removes a single dynamic ColVariable
 void remove_dynamic_variable( const ColVariable * var ) override;

 /// removes a single dynamic bound (OneVarConstraint)
 void remove_dynamic_bound( const OneVarConstraint * con ) override;

/** @} ---------------------------------------------------------------------*/
 /// the name of the MOSEK parameter a Solver integer parameter is mapped on
 /** Returns the empty string if there is none. */
 std::string msk_int_par_map( idx_type par ) const;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// the names of the MOSEK parameters a Solver double parameter is mapped on
 /** All of them are set to the value of the parameter, the first one gives
  * the value read; the vector is empty if there is none. */
 std::vector< std::string > msk_dbl_par_map( idx_type par ) const;

/*--------------------------------------------------------------------------*/
/*-------------------- PROTECTED FIELDS OF THE CLASS -----------------------*/
/*--------------------------------------------------------------------------*/

 MSKenv_t env;    ///< MOSEK environment, shared by all the MOSEKMILPSolver
 MSKtask_t task;  ///< the task of MOSEK holding the problem

 /// how many MOSEKMILPSolver are using the shared environment
 static unsigned int f_env_users;

 /// the shared environment itself, and the lock that guards its life
 static MSKenv_t f_shared_env;
 static std::mutex f_env_mutex;

 /// the solution of the task that the last solve has produced
 /** The integer one for a MIP, the basic one for an LP solved by a simplex
  * algorithm or by the interior point followed by the basis
  * identification, the interior point one otherwise [see
  * guts_of_compute()]. */
 MSKsoltypee f_soltype = MSK_SOL_ITR;

 bool f_callback_set = false;  ///< true if the callback has been set

 /// the bound and the incumbent of the search, as seen by the callback
 double f_cb_bound;
 double f_cb_incumbent;

 /// the name of the file the log goes to, if any (strLogFileName)
 std::string f_log_file;

 /// the Hessian of the objective, one element per column of the model
 /** MOSEK only takes the quadratic part of the objective as a whole, so
  * this is the copy that every change is applied to before it is given
  * again [see pass_hessian()]. The element j holds the pairs ( i , value )
  * of column j with i >= j, ordered by i, i.e., the lower triangle with
  * the diagonal first; the elements follow the columns of the model as
  * they are added and removed, and so do the indices in the pairs. */
 std::vector< std::vector< std::pair< int , double > > > q_hessian;

 /// the number of nonzeros in q_hessian
 int q_hessian_nnz = 0;

 /// true if q_hessian has changed since it was last given to MOSEK
 bool f_hessian_changed = false;

 /** @name Handling of MOSEK parameters
  *
  * The following maps keep the relationship between the parameter system
  * of SMS++ and the parameters of MOSEK, which can therefore be used as if
  * they were SMS++ parameters with the same (symbolic) names, for example
  * in configuration files. They are generated by MILPSolver/tools/mosek_pars
  * for the version of MOSEK being compiled against.
  * @{ */

 const static std::array< std::string , MOSEK_NUM_INT_PARS >
  SMSpp_to_MOSEK_int_pars;
 const static std::array< std::string , MOSEK_NUM_DBL_PARS >
  SMSpp_to_MOSEK_dbl_pars;
 const static std::array< std::string , MOSEK_NUM_STR_PARS >
  SMSpp_to_MOSEK_str_pars;

 const static std::array< std::pair< std::string , int > ,
                          MOSEK_NUM_INT_PARS > MOSEK_to_SMSpp_int_pars;
 const static std::array< std::pair< std::string , int > ,
                          MOSEK_NUM_DBL_PARS > MOSEK_to_SMSpp_dbl_pars;
 const static std::array< std::pair< std::string , int > ,
                          MOSEK_NUM_STR_PARS > MOSEK_to_SMSpp_str_pars;

/** @} ---------------------------------------------------------------------*/

 double UpCutOff;  ///< externally set upper cutoff to terminate
 double LwCutOff;  ///< externally set lower cutoff to terminate

/*--------------------------------------------------------------------------*/
/*---------------------- PRIVATE PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

 private:

 /// throws std::runtime_error if r is an error of MOSEK, saying where
 void check( MSKrescodee r , const char * where ) const;

 /// the SMS++ status of the last solve, out of what MOSEK has returned
 int decode_status( MSKrescodee r , MSKrescodee trm );

 /// the solution status of f_soltype, MSK_SOL_STA_UNKNOWN if undefined
 MSKsolstae solution_status( void ) const;

 /// true if the solution f_soltype has a primal part that can be used
 bool has_primal( void ) const;

 /// the value of the solution f_soltype, constant_value comprised
 OFValue primal_value( void ) const;

 /// the bound of the search after a stopped MIP, none if there is none
 OFValue stopped_dual_bound( OFValue none ) const;

 /// the bound key of MOSEK for the two sides lo and up
 static MSKboundkeye bound_key( double lo , double up );

 /// the two sides of a row of the Block, as MOSEK takes them
 static std::array< double , 2 > row_sides( const FRowConstraint * con );

 /// sets the bounds of the column j of the task
 void put_var_bound( int j , const std::array< double , 2 > & bd );

 /// sets the bounds of a set of columns of the task in one call
 void put_var_bounds( const std::vector< int > & idx ,
                      const std::vector< std::array< double , 2 > > & bd );

 /// sets the sides of a set of rows of the task in one call
 void put_con_bounds( const std::vector< int > & idx ,
                      const std::vector< std::array< double , 2 > > & sd );

 /// the type of the column of var in the task
 MSKvariabletypee var_type( const ColVariable * var ) const;

 /// links the log of the task to std::cout or to f_log_file, or unlinks it
 void link_log( void );

 /// builds q_hessian out of the quadratic objective of MILPSolver
 void generate_qobj_hessian( void );

 /// adds delta to the entry ( i , j ) of the Hessian, removing it if zero
 void add_to_hessian( int i , int j , double delta );

 /// sets the entry ( i , j ) of the Hessian to value, zero meaning removed
 void set_hessian( int i , int j , double value );

 /// the column c of the model is gone: so is it in the Hessian
 void remove_hessian_column( int c );

 /// gives q_hessian to MOSEK if it has changed
 void pass_hessian( void );

 /// changes the costs of the columns idx, adding val to them or setting them
 /** Changes the costs of the columns in idx in one call, reading the old ones
  * in one call as well if add is true. The indices need not be ordered, and
  * a column named twice gets the sum of its values if add is true and the
  * last one otherwise. If add and removal are both true, a finite result
  * is rounded to zero when its magnitude is at most
  * dblObjRemovalTol * max( 1 , abs( old cost ) , abs( summed delta ) ). */
 void change_obj_costs( std::vector< int > & idx ,
                        std::vector< double > & val ,
                        bool add , bool removal = false );

/*--------------------------------------------------------------------------*/

 SMSpp_insert_in_factory_h;

/*--------------------------------------------------------------------------*/

 };  // end( class( MOSEKMILPSolver ) )

/*--------------------------------------------------------------------------*/

}  // end( namespace SMSpp_di_unipi_it )

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

#endif  /* MOSEKMILPSolver.h included */

/*--------------------------------------------------------------------------*/
/*----------------------- End File MOSEKMILPSolver.h -----------------------*/
/*--------------------------------------------------------------------------*/
