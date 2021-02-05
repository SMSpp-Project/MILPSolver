/*--------------------------------------------------------------------------*/
/*---------------------- File SimpleMILPBlock.h ----------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the *concrete* class SimpleMILPBlock, which implements the
 * Block concept [see Block.h] for a very simple "flat" problem having only
 * a fixed number of linear constraints on a fixed number of variables,
 * everything only represented in "abstract form.
 *
 * \version 0.21
 *
 * \date 31 - 01 - 2019
 *
 * \author Antonio Frangioni \n
 *         Operations Research Group \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * Copyright &copy; by Antonio Frangioni
 */
/*--------------------------------------------------------------------------*/
/*----------------------------- DEFINITIONS --------------------------------*/
/*--------------------------------------------------------------------------*/

#ifndef __SimpleMILPBlock
 #define __SimpleMILPBlock
                      /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "Block.h"
#include "LinearFunction.h"
#include "FRowConstraint.h"
#include "FRealObjective.h"
#include "OneVarConstraint.h"

/*--------------------------------------------------------------------------*/
/*--------------------------- NAMESPACE ------------------------------------*/
/*--------------------------------------------------------------------------*/

/// namespace for the Structured Modeling System++ (SMS++)
namespace SMSpp_di_unipi_it
{
/*--------------------------------------------------------------------------*/
/*-------------------- SimpleMILPBlock-RELATED TYPES -----------------------*/
/*--------------------------------------------------------------------------*/
/** @defgroup SimpleMILPBlock_TYPES SimpleMILPBlock-related types
 *  @{ */

/** @} end( group( SimpleMILPBlock_TYPES ) ) */ 
/*--------------------------------------------------------------------------*/
/*------------------------------- CLASSES ----------------------------------*/
/*--------------------------------------------------------------------------*/
/** @defgroup SimpleMILPBlock_CLASSES Classes in SimpleMILPBlock.h
 *  @{ */

/*--------------------------------------------------------------------------*/
/*---------------------- CLASS SimpleMILPBlock -----------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------- GENERAL NOTES --------------------------------*/
/*--------------------------------------------------------------------------*/
/// implementation of a simple MILP Block concept
/** The SimpleMILPBlock class implements the Block concept [see Block.h] for
 * a very simple "flat" problem having only a fixed number of linear
 * constraints on a fixed number of variables, everything only represented
 * in "abstract form.
 *
 * In particular, the problem represented by a SimpleMILPBlock is a 
 * Mixed-Integer Linear Program of the generic form
 * \f[
 *  \min \{ \bar{c} + cx : \underline{b} \leq Ax \leq \bar{b} ,
 *                         l \leq x \leq u \}
 * \f]
 * where \f$ c, l, u \f$ are real n-vectors, \f$ \underline{b}, \bar{b} \f$
 * are real m-vectors, A is a real m-times-n matrix, and \f$ \bar{c} \f$ is
 * a real constant. Actually, entries of \f$ \underline{b}, l \f$ can be
 * \f$ -\infty \f$, and entries of \f$ \bar{b}, u \f$ can be \f$ +\infty \f$.
 * The coefficient matrix A is usually very sparse, which justifies the
 * choice of a standard sparse (row-wise) format to represent it. */

class SimpleMILPBlock : public Block {

/*--------------------------------------------------------------------------*/
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

public:

/*--------------------------------------------------------------------------*/
/*---------------------------- PUBLIC TYPES --------------------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*--------------------- PUBLIC METHODS OF THE CLASS ------------------------*/
/*--------------------------------------------------------------------------*/
/*---------------------------- CONSTRUCTOR ---------------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Constructor and Destructor
 *  @{ */

 /// constructor of SimpleMILPBlock
 /** Constructor of SimpleMILPBlock. It accepts a pointer to the father
  * Block, which can be of any type. */

 SimpleMILPBlock( Block *father = nullptr ) : Block( father ) {}

/*--------------------------------------------------------------------------*/
 /// destructor of SimpleMILPBlock
 /** Destructor of SimpleMILPBlock. Since all Variable, Constraint and
  * Objective live in fields of the class, it explicitly clear() all
  * Constraint and Objective to ensure that they do not make any reference to
  * no-longer-existing Variable while they are destroyed. */

 virtual ~SimpleMILPBlock() {
  for( int i = Bx.size() ; i-- ; )
   Bx[ i ].clear();

  for( int i = A.size() ; i-- ; )
   A[ i ].clear();

  f.clear();
  }

/**@} ----------------------------------------------------------------------*/
/*-------------------------------- ACCESSORS -------------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Accessors
 *  @{ */

 std::vector<ColVariable> & get_x( void ) { return( x ); }

/**@} ----------------------------------------------------------------------*/
/*-------------------- PROTECTED PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

 protected:

/*--------------------------------------------------------------------------*/
/*-------------------------- PROTECTED METHODS -----------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Protected methods for inserting and extracting
    @{ */

 virtual void print( std::ostream &output ) const override;
 ///< print the SimpleMILPBlock on an ostream with the given verbosity

/*--------------------------------------------------------------------------*/

 virtual void load( std::istream &input ) override;

 ///< load the SimpleMILPBlock out of an istream
 /**< Load the SimpleMILPBlock out of an istream. The format is:
  *
  * number of variables
  *
  * number of constraints
  *
  * for_each( variable )
  *   lower bound   upper bound
  *
  * number of nonzeroes in the objective function
  * for_each( nonzero )
  *   index of variable   real coefficient
  *
  * for_each( constraint )
  *   lower bound   upper bound
  *   number of nonzeroes in the constraint
  *   for_each( nonzero )
  *     index of variable , real coefficient
  *
  * Infinity values are "INF", "inf", (actually, anything beginning with "I"
  * or "i" is taken as +INF), "-INF", and "-inf" (...).
  *
  * The values can be separated by any whitespace (space, tab, newline). If a
  * '#' is found, all the rest of the line is treated as a comment and
  * discarded. */

/*--------------------------------------------------------------------------*/
 /// extends Block::serialize( netCDF::NcGroup )
 /** Extends Block::serialize( netCDF::NcGroup ) to the specific format of a
  * SimpleMILPBlock. Besides what is managed by the serialize() method of
  * the base Block class, the group will contain the following:
  *
  * - the dimension "num_vars" containing the number of Variables (columns)
  *   of the SimpleMILPBlock;
  *
  * - the dimension "num_constraints" containing the number of Constraints
  *   (rows) of the SimpleMILPBlock;
  *
  * - the dimension "nzcount" containing the total number of nonzero
  *   coefficients in the Constraints of the SimpleMILPBlock;
  *
  * - the variable "rmatval", of type double and indexed over the dimension
  *   "nzcount", that contains the nonzeros in the coefficient matrix of the
  *   MILP, arranged row-wise (see rmatbeg);
  *
  * - the variable "rmatind", of type UInt64 and indexed over the dimension
  *   "nzcount", that contains the column indices of the nonzeros in the
  *   coefficient matrix of the MILP, arranged row-wise; that is, the i-th
  *   entry of rmatind indicates to which variable (column) the (allegedly,
  *   nonzero) coefficient to be found in the i-th entry of rmatval refers;
  *
  * - the variable "rmatbeg", of type UInt64 and indexed over the dimension
  *   "num_constraints", that contains the indices into the rmatval/rmatind
  *   vectors where each constraint (row of the coefficient matrix of the
  *   MILP) begins: all the nonzeros corresponding to constraint h = 0, ...,
  *   num_constraints - 1 are to be found in rmatval[ rmatbeg[ h ] ], ...,
  *   rmatval[ rmatbeg[ h + 1 ] - 1 ], with their column (variable) indices
  *   to be found in the corresponding entries of rmatind[];
  *
  * - the variable "x_lb", of type double and indexed over the dimension
  *   "num_vars", whose i-th entry is the (possibly, - infinite) lower bound
  *   on the feasible value of the i-th variable;
  *
  * - the variable "x_ub", of type double and indexed over the dimension
  *   "num_vars", whose i-th entry is the (possibly, + infinite) upper bound
  *   on the feasible value of the i-th variable;
  *
  * - the variable "row_lb", of type double and indexed over the dimension
  *   ""num_constraints", whose i-th entry is the (possibly, - infinite)
  *   lower bound on the feasible value that the linear expression of the
  *   i-th constraint (as specified by the coefficients in rmatval, rmatind,
  *   and rmatbeg) can take;
  *
  * - the variable "row_ub", of type double and indexed over the dimension
  *   ""num_constraints", whose i-th entry is the (possibly, + infinite)
  *   upper bound on the feasible value that the linear expression of the
  *   i-th constraint (as specified by the coefficients in rmatval, rmatind,
  *   and rmatbeg) can take;
  *
  * - the variable "obj", of type double and indexed over the dimension
  *   "num_vars", whose i-th entry is thecoefficient of the i-th variable in
  *   the linear objective function of the MILP;
  *
  * - the single variable "obj_offset", of type double (and, since it is a
  *   scalar, not indexed over any dimension) giving the fixed offset in
  *   the linear objective function of the MILP.
  *
  * All these are mandatory. */

 virtual void serialize( netCDF::NcGroup & file ) const override;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// extends Block::deserialize( netCDF::NcGroup )
 /** Extends Block::deserialize( netCDF::NcGroup ) to the specific format of
  * a SimpleMILPBlock. See SimpleMILPBlock::serialize( netCDF::NcGroup )
  * for details of the format of the read netCDF group.
  *
  * The method dispatches the protected guts_of_deserialize(), and then the
  * method of the base class, in this order. This is done so that the
  * necessary NBModification is issued last, while allowing derived classes
  * to call the (equivalent to) SimpleMILPBlock::deserialize() at any point
  * in their deserialize() is they so need; see comments to
  * Block::deserialize() for a complete discussion. */

 void deserialize( const netCDF::NcGroup & group ) override
 {
  guts_of_deserialize( group );
  Block::deserialize( group );
  }

/**@} ----------------------------------------------------------------------*/
/*--------------------------- PROTECTED FIELDS  ----------------------------*/
/*--------------------------------------------------------------------------*/

 protected:

/*--------------------------------------------------------------------------*/

 void guts_of_deserialize( const netCDF::NcGroup & group );
 
/*--------------------------------------------------------------------------*/
/*-------------------------- PROTECTED METHODS  ----------------------------*/
/*--------------------------------------------------------------------------*/

 std::vector<ColVariable> x;       /// the variables
 std::vector<FRowConstraint> A;    /// the constraints
 std::vector<BoxConstraint> Bx;    /// the box constraints

 FRealObjective f;                 /// the objective function

/*--------------------------------------------------------------------------*/
/*--------------------- PRIVATE PART OF THE CLASS --------------------------*/
/*--------------------------------------------------------------------------*/

 private:

/*--------------------------------------------------------------------------*/
/*-------------------------- PRIVATE METHODS -------------------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*---------------------------- PRIVATE FIELDS ------------------------------*/
/*--------------------------------------------------------------------------*/

 SMSpp_insert_in_factory_h;        // insert it in the Block factory

/*--------------------------------------------------------------------------*/

};  // end( class( SimpleMILPBlock ) )

/*--------------------------------------------------------------------------*/

/*@}  end( group( SimpleMILPBlock_CLASSES ) ) ------------------------------*/
/*--------------------------------------------------------------------------*/

 }  // end( namespace SMSpp_di_unipi_it )

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

#endif  /* SimpleMILPBlock.h included */

/*--------------------------------------------------------------------------*/
/*--------------------- End File SimpleMILPBlock.h -------------------------*/
/*--------------------------------------------------------------------------*/
