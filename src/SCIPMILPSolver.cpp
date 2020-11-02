/*--------------------------------------------------------------------------*/
/*------------------------- File SCIPMILPSolver.cpp -------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Implementation of the SCIPMILPSolver class.
 *
 * \author Antonio Frangioni \n
 *         Operations Research Group \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \author Niccolò Iardella \n
 *         Operations Research Group \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \copyright &copy by Antonio Frangioni, Niccolò Iardella
 */
/*--------------------------------------------------------------------------*/
/*---------------------------- IMPLEMENTATION ------------------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include <functional>
#include <queue>
#include <cstdio>
#include <cmath>

#include <Block.h>
#include <MILPSolver.h>
#include <OneVarConstraint.h>
#include <LinearFunction.h>
#include <DQuadFunction.h>

#include "SCIPMILPSolver.h"
#include <scip/scipdefplugins.h>
#include <scip/cons_linear.h>

// Logging
#include <boost/log/core.hpp>
#include <boost/log/trivial.hpp>
#include <boost/log/expressions.hpp>

/*--------------------------------------------------------------------------*/
/*------------------------- NAMESPACE AND USING ----------------------------*/
/*--------------------------------------------------------------------------*/

using namespace SMSpp_di_unipi_it;

SMSpp_insert_in_factory_cpp_0( SCIPMILPSolver );

/*--------------------------------------------------------------------------*/
/*--------------------- CONSTRUCTOR AND DESTRUCTOR -------------------------*/
/*--------------------------------------------------------------------------*/

SCIPMILPSolver::SCIPMILPSolver() : MILPSolver() {
 SCIP_CALL_ABORT( SCIPcreate( &scip ) );
 SCIP_CALL_ABORT( SCIPincludeDefaultPlugins( scip ) );
}

SCIPMILPSolver::~SCIPMILPSolver() {
 SCIPfree( &scip );
}

/*--------------------------------------------------------------------------*/
/*--------------------- DERIVED METHODS OF BASE CLASS ----------------------*/
/*--------------------------------------------------------------------------*/

void SCIPMILPSolver::set_Block( Block * block ) {
 if( block == f_Block ) {
  return;
 }
 MILPSolver::set_Block( block );
}

/*--------------------------------------------------------------------------*/

void SCIPMILPSolver::clear_problem() {
 MILPSolver::clear_problem();

 SCIP_CALL_ABORT( SCIPfreeProb( scip ) );
}

/*--------------------------------------------------------------------------*/

void SCIPMILPSolver::load_problem() {
 MILPSolver::load_problem();

 SCIP_CALL_ABORT( SCIPcreateProbBasic( scip, prob_name.c_str() ) );
 SCIP_CALL_ABORT( SCIPsetObjsense( scip, ( SCIP_OBJSENSE ) objsense ) );
 vars.resize( numcols );

 bool is_mip = std::any_of( xctype.begin(),
                            xctype.end(),
                            []( char c ) { return c != 'C'; } );

 for( int i = 0; i < numcols; ++i ) {
  char * name = use_custom_names ? colname[ i ] : nullptr;
  SCIP_Real collb = ( lb[ i ] == -Inf< double >() ) ?
                    -SCIPinfinity( scip ) : lb[ i ];
  SCIP_Real colub = ( ub[ i ] == Inf< double >() ) ?
                    SCIPinfinity( scip ) : ub[ i ];

  SCIP_VARTYPE vartype = SCIP_VARTYPE_CONTINUOUS;

  if( is_mip ) {
   switch( xctype[ i ] ) {
    default:
    case 'C':
     break;
    case 'B':
     collb = 0;
     colub = 1;
     vartype = SCIP_VARTYPE_BINARY;
     break;
    case 'I':
     vartype = SCIP_VARTYPE_INTEGER;
     break;
    case 'S':
    case 'N':
     SCIPABORT();
   }
  }
  SCIP_VAR * var = nullptr;
  SCIP_CALL_ABORT( SCIPcreateVarBasic( scip, &var, name, collb, colub,
                                       objective[ i ],
                                       SCIP_VARTYPE_CONTINUOUS ) );
  SCIP_CALL_ABORT( SCIPaddVar( scip, var ) );
  vars[ i ] = var;
  SCIP_CALL_ABORT( SCIPreleaseVar( scip, &var ) );
 }

 conss.resize( numrows );
 for( int i = 0; i < numrows; ++i ) {
  SCIP_Real conslhs = NAN;
  SCIP_Real consrhs = NAN;

  char * name = use_custom_names ? rowname[ i ] : nullptr;

  switch( sense[ i ] ) {
   case 'L':
    conslhs = -SCIPinfinity( scip );
    consrhs = rhs[ i ];
    break;
   case 'E':
    conslhs = rhs[ i ];
    consrhs = rhs[ i ];
    break;
   case 'G':
    conslhs = rhs[ i ];
    consrhs = SCIPinfinity( scip );
    break;
   case 'R':
    conslhs = rhs[ i ];
    consrhs = rhs[ i ];
    if( rngval[ i ] > 0 )
     consrhs += rngval[ i ];
    else
     conslhs += rngval[ i ];
  }

  SCIP_CONS * cons = nullptr;
  SCIP_CALL_ABORT( SCIPcreateConsBasicLinear( scip, &cons, name, 0,
                                              nullptr, nullptr,
                                              conslhs, consrhs ) );
  SCIP_CALL_ABORT( SCIPaddCons( scip, cons ) );
  conss[ i ] = cons;
  SCIP_CALL_ABORT( SCIPreleaseCons( scip, &cons ) );
 }

 for( int c = 0; c < numcols; ++c ) {
  for( int i = matbeg[ c ]; i < matbeg[ c + 1 ]; ++i )
   SCIP_CALL_ABORT( SCIPaddCoefLinear( scip,
                                       conss[ matind[ i ] ],
                                       vars[ c ],
                                       matval[ i ] ) );
 }

 bool is_qp = std::any_of( q_objective.begin(),
                           q_objective.end(),
                           []( double d ) { return d != 0; } );
 if( is_qp ) {
  SCIPABORT();
 }
}

/*--------------------------------------------------------------------------*/

double SCIPMILPSolver::get_problem_lb( const ColVariable & var ) {
 double b = MILPSolver::get_problem_lb( var );
 if( b == -Inf< double >() ) {
  b = -SCIPinfinity( scip );
 }
 return b;
}

/*--------------------------------------------------------------------------*/

double SCIPMILPSolver::get_problem_ub( const ColVariable & var ) {
 double b = MILPSolver::get_problem_ub( var );
 if( b == Inf< double >() ) {
  b = SCIPinfinity( scip );;
 }
 return b;
}

/*--------------------------------------------------------------------------*/

int SCIPMILPSolver::compute( bool changedvars ) {
 process_modifications();

 if( !output_file.empty() ) {
  SCIP_CALL_ABORT( SCIPwriteOrigProblem( scip, output_file
   .c_str(), "mps", FALSE ) );
 }

 SCIP_CALL_ABORT( SCIPsolve( scip ) );

 SCIP_STATUS status = SCIPgetStatus( scip );
 switch( status ) {
  case SCIP_STATUS_OPTIMAL:
  case SCIP_STATUS_GAPLIMIT:
  case SCIP_STATUS_SOLLIMIT:
   sol_status = kOK;
   break;
  case SCIP_STATUS_INFEASIBLE:
   sol_status = kInfeasible;
   break;
  case SCIP_STATUS_NODELIMIT:
   sol_status = kStopIter;
   break;
  case SCIP_STATUS_TIMELIMIT:
   sol_status = kStopTime;
   break;
  default:
   sol_status = kError;
   break;
  case SCIP_STATUS_INFORUNBD:
  case SCIP_STATUS_UNBOUNDED:
   sol_status = kUnbounded;
 }
 nodes = static_cast<int>(SCIPgetNTotalNodes( scip ));

 return sol_status;
}

/*--------------------------------------------------------------------------*/

Solver::OFValue SCIPMILPSolver::get_lb() {

 OFValue lower_bound = 0;

 switch( objsense ) {
  case 1: // Minimization problem
   switch( sol_status ) {
    case kUnbounded:
     lower_bound = -Inf< OFValue >();
     break;
    case kInfeasible:
     lower_bound = Inf< OFValue >();
     break;
    default:
     lower_bound = SCIPgetDualbound( scip );
     break;
   }
   break;
  case -1: // Maximization problem
   switch( sol_status ) {
    case kUnbounded:
     lower_bound = Inf< OFValue >();
     break;
    case kInfeasible:
     lower_bound = -Inf< OFValue >();
     break;
    default:
     lower_bound = SCIPgetDualbound( scip );
     break;
   }
   break;
  default:
   throw std::runtime_error( "Objective type not yet defined" );
   break;
 }

 return lower_bound;
}

/*--------------------------------------------------------------------------*/

Solver::OFValue SCIPMILPSolver::get_ub() {

 OFValue upper_bound = 0;

 switch( objsense ) {
  case 1: // Minimization problem
   switch( sol_status ) {
    case kUnbounded:
     upper_bound = -Inf< OFValue >();
     break;
    case kInfeasible:
     upper_bound = Inf< OFValue >();
     break;
    default:
     upper_bound = SCIPgetPrimalbound( scip );
     break;
   }
   break;
  case -1: // Maximization problem
   switch( sol_status ) {
    case kUnbounded:
     upper_bound = Inf< OFValue >();
     break;
    case kInfeasible:
     upper_bound = -Inf< OFValue >();
     break;
    default:
     upper_bound = SCIPgetPrimalbound( scip );
     break;
   }
   break;
  default:
   throw std::runtime_error( "Objective type not yet defined" );
   break;
 }

 return upper_bound;
}

/*--------------------------------------------------------------------------*/

bool SCIPMILPSolver::has_var_solution() {
 switch( sol_status ) {
  case ( kOK ):
  case ( kInfeasible ):
   return ( true );
  default:
   return ( false );
 }
}

/*--------------------------------------------------------------------------*/

bool SCIPMILPSolver::is_var_feasible() {
 switch( sol_status ) {
  case ( kInfeasible ):
   return ( false );
  default:
   return ( true );
 }
}

/*--------------------------------------------------------------------------*/

Solver::OFValue SCIPMILPSolver::get_var_value() {
 switch( objsense ) {
  case 1: // Minimization problem
   return get_ub();
   break;
  case -1: // Maximization problem
   return get_lb();
   break;
  default:
   throw std::runtime_error( "Objective type not yet defined" );
   break;
 }
}

/*--------------------------------------------------------------------------*/

void SCIPMILPSolver::get_var_solution( Configuration * solc ) {

 SCIP_SOL * sol = SCIPgetBestSol( scip );

 if( sol == nullptr ) {
  return;
 }

 auto * x = new double[numcols];
 SCIP_RETCODE status = SCIPgetSolVals( scip, sol, numcols, vars.data(), x );
 if( status != SCIP_OKAY ) {
  delete[] x;
  throw std::runtime_error( "Unable to get the solution values" );
 }

 int col = 0;

 std::queue< Block * > Q;

 bool owned = f_Block->is_owned_by( f_id );
 if( !owned && !f_Block->lock( f_id ) ) {
  throw std::runtime_error( "Unable to lock the Block" );
 }

 Q.push( f_Block );

 while( !Q.empty() ) {
  Block * q_Block = Q.front();
  Q.pop();

  for( auto * i : q_Block->get_nested_Blocks() ) {
   Q.push( i );
  }

  for( const auto & i : q_Block->get_static_variables() ) {
   auto f1 = std::bind( &SCIPMILPSolver::set_var_value,
                        this,
                        std::placeholders::_1,
                        std::ref( x ),
                        std::ref( col ) );
   un_any_const_static( i, f1, un_any_type< ColVariable >() );
  }

  for( const auto & i : q_Block->get_dynamic_variables() ) {
   auto f1 = std::bind( &SCIPMILPSolver::set_var_value,
                        this,
                        std::placeholders::_1,
                        std::ref( x ),
                        std::ref( col ) );
   un_any_const_dynamic( i, f1, un_any_type< ColVariable >() );
  }
 }

 if( !owned ) {
  f_Block->unlock( f_id );
 }

 delete[]x;
}

/*--------------------------------------------------------------------------*/

bool SCIPMILPSolver::has_dual_solution() {
 return false;
}

/*--------------------------------------------------------------------------*/

void SCIPMILPSolver::get_dual_solution( Configuration * solc ) {
 SCIPABORT();
}

/*--------------------------------------------------------------------------*/

void SCIPMILPSolver::write_lp( const std::string & filename ) {
 SCIP_CALL_ABORT( SCIPwriteOrigProblem( scip,
                                        filename.c_str(),
                                        "mps", FALSE ) );
}

/*--------------------------------------------------------------------------*/
/*-------------------- PROTECTED FIELDS OF THE CLASS -----------------------*/
/*--------------------------------------------------------------------------*/

void SCIPMILPSolver::var_modification( VariableMod * mod ) {
 MILPSolver::var_modification( mod );

 if( SCIPisTransformed( scip ) ) {
  SCIP_CALL_ABORT( SCIPfreeTransform( scip ) );
 }

 auto * var = dynamic_cast<ColVariable *>(mod->variable());
 int idx = index_of_variable( var );
 SCIP_VARTYPE oldtype = SCIPvarGetType( vars[ idx ] );

 SCIP_Bool infeas = FALSE;

 if( var->is_fixed() ) {
  // Fix the variable
  SCIP_Bool fixed = 0;
  SCIP_CALL_ABORT( SCIPfixVar( scip, vars[ idx ],
                               var->get_value(), &infeas, &fixed ) );
  assert( fixed );

 } else if( var->is_integer() ) {
  // Unfix or refresh the variable
  if( var->is_unitary() && var->is_positive() ) {
   if( oldtype != SCIP_VARTYPE_BINARY ) {
    SCIP_CALL_ABORT( SCIPchgVarType( scip, vars[ idx ],
                                     SCIP_VARTYPE_BINARY, &infeas ) );
   }
  } else if( oldtype != SCIP_VARTYPE_INTEGER ) {
   SCIP_CALL_ABORT( SCIPchgVarType( scip, vars[ idx ],
                                    SCIP_VARTYPE_INTEGER, &infeas ) );
  }
 } else if( oldtype != SCIP_VARTYPE_CONTINUOUS ) {
  SCIP_CALL_ABORT( SCIPchgVarType( scip, vars[ idx ],
                                   SCIP_VARTYPE_CONTINUOUS, &infeas ) );
 }
 assert( !infeas );

 SCIP_Real lb = var->get_lb() == -Inf< double >() ?
                -SCIPinfinity( scip ) : var->get_lb();
 SCIP_Real ub = var->get_ub() == Inf< double >() ?
                SCIPinfinity( scip ) : var->get_ub();

 if( SCIPvarGetLbOriginal( vars[ idx ] ) != lb )
  SCIP_CALL_ABORT( SCIPchgVarLb( scip, vars[ idx ], lb ) );

 if( SCIPvarGetUbOriginal( vars[ idx ] ) != ub )
  SCIP_CALL_ABORT( SCIPchgVarUb( scip, vars[ idx ], ub ) );
}

/*--------------------------------------------------------------------------*/

void SCIPMILPSolver::of_modification( ObjectiveMod * mod ) {
 MILPSolver::of_modification( mod );

 if( SCIPisTransformed( scip ) )
  SCIP_CALL_ABORT( SCIPfreeTransform( scip ) );
 /*
  * ObjectiveMod class does not include any modification types except
  * for eSetMin and eSetMax.
  * To change OF coefficents, a FunctionMod must be used.
  */

 switch( mod->type() ) {

  case ObjectiveMod::eSetMin:
   SCIP_CALL_ABORT( SCIPsetObjsense( scip, SCIP_OBJSENSE_MINIMIZE ) );
   break;

  case ObjectiveMod::eSetMax:
   SCIP_CALL_ABORT( SCIPsetObjsense( scip, SCIP_OBJSENSE_MAXIMIZE ) );
   break;

  default:
   throw std::invalid_argument( "Invalid type of ObjectiveMod" );
 }
}

/*--------------------------------------------------------------------------*/

void SCIPMILPSolver::const_modification( ConstraintMod * mod ) {
 MILPSolver::const_modification( mod );

 /*
  * To change the coefficents, a FunctionMod must be used.
  */

 if( SCIPisTransformed( scip ) )
  SCIP_CALL_ABORT( SCIPfreeTransform( scip ) );

 auto * p_const = dynamic_cast<FRowConstraint *>(mod->constraint());

 RowConstraint::RHSValue const_lhs = NAN;
 RowConstraint::RHSValue const_rhs = NAN;

 SCIP_CONS * cons = conss[ index_of_constraint( p_const ) ];

 switch( mod->type() ) {

  case ConstraintMod::eRelaxConst:
   // In order to relax the constraint all we do is transform it
   // into an inequality with RHS equal to infinity
   SCIP_CALL_ABORT( SCIPchgLhsLinear( scip, cons, -SCIPinfinity( scip ) ) );
   SCIP_CALL_ABORT( SCIPchgRhsLinear( scip, cons, SCIPinfinity( scip ) ) );
   break;

  case ConstraintMod::eEnforceConst:
  case RowConstraintMod::eChgLHS:
  case RowConstraintMod::eChgRHS:
  case RowConstraintMod::eChgBTS:
   // In order to enforce a relaxed constraint all we need to do is
   // reverse the process of relaxing it, by changing the sense and
   // the rhs back to the original form of the constraint
   // Moreover, for the way the LP vectors are built, handling
   // LHS/RHS/BTS cases separately is not worth it.

   const_lhs =
    p_const->get_lhs() == -Inf< double >() ?
    -SCIPinfinity( scip ) : p_const->get_lhs();
   const_rhs = p_const->get_rhs() == Inf< double >() ?
               SCIPinfinity( scip ) : p_const->get_rhs();

   SCIP_CALL_ABORT( SCIPchgLhsLinear( scip, cons, const_lhs ) );
   SCIP_CALL_ABORT( SCIPchgRhsLinear( scip, cons, const_rhs ) );

   break;

  default:
   throw std::invalid_argument( "Invalid type of ObjectiveMod" );
 }
}

/*--------------------------------------------------------------------------*/

void SCIPMILPSolver::bound_modification( OneVarConstraintMod * mod ) {
 MILPSolver::bound_modification( mod );

 /*
  * The same ColVariable can have more active OneVarConstraints,
  * so each time we modify one of them we have to check if LHS and RHS
  * of the Variable change.
  */

 if( SCIPisTransformed( scip ) )
  SCIP_CALL_ABORT( SCIPfreeTransform( scip ) );

 auto * p_const = dynamic_cast<OneVarConstraint *>(mod->constraint());
 auto * p_var = dynamic_cast<ColVariable *>(p_const->get_active_var( 0 ));

 int vidx = index_of_variable( p_var );
 std::vector< int > indices( 2, index_of_variable( p_var ) );
 std::vector< char > lu;
 std::vector< double > bd;

 SCIP_Real lb;
 SCIP_Real ub;
 SCIP_VAR * var = vars[ index_of_variable( p_var ) ];

 switch( mod->type() ) {

  case RowConstraintMod::eChgLHS:
   lb = get_problem_lb( *p_var );
   SCIP_CALL_ABORT( SCIPchgVarLb( scip, var, lb ) );
   break;

  case RowConstraintMod::eChgRHS:
   ub = get_problem_ub( *p_var );
   SCIP_CALL_ABORT( SCIPchgVarUb( scip, var, ub ) );
   break;

  case RowConstraintMod::eChgBTS:
   lb = get_problem_lb( *p_var );
   ub = get_problem_ub( *p_var );
   SCIP_CALL_ABORT( SCIPchgVarLb( scip, var, lb ) );
   SCIP_CALL_ABORT( SCIPchgVarUb( scip, var, ub ) );
   break;

  default:
   throw std::invalid_argument( "Invalid type of OneVarConstraintMod" );
 }
}

/*--------------------------------------------------------------------------*/

void SCIPMILPSolver::function_modification( FunctionMod * mod ) {
 MILPSolver::function_modification( mod );

 if( SCIPisTransformed( scip ) )
  SCIP_CALL_ABORT( SCIPfreeTransform( scip ) );

 /*
  * This function is used when changing coefficents for OFs or constraints.
  */

 // Check if OF or a Constraint is involved
 auto * mod_f = mod->function();
 bool changing_of = false;
 const auto * lf = dynamic_cast<const LinearFunction *> (mod_f);
 const auto * qf = dynamic_cast<const DQuadFunction *> (mod_f);

 auto * p_obj = dynamic_cast< FRealObjective * >( f_Block->get_objective() );
 if( p_obj != nullptr ) {
  auto * of = p_obj->get_function();
  if( of == mod_f ) {
   changing_of = true;
  }
 }

 if( changing_of ) {
  // Changing objective function

  if( lf != nullptr ) {
   // Linear objective function

   for( auto el : lf->get_v_var() ) {
    SCIP_VAR * var = vars[ index_of_variable( el.first ) ];
    SCIP_CALL_ABORT( SCIPchgVarObj( scip, var, el.second ) );
   }
  } else if( qf != nullptr ) {
   // Quadratic objective function
   SCIPABORT();

  } else {
   // This should never happen
   throw std::invalid_argument( "Unknown type of Objective Function" );
  }

 } else {
// Changing coefficients of a constraint
  if( lf != nullptr ) {
   auto * p_const = ( FRowConstraint * ) lf->get_Observer();

   SCIP_CONS * cons = conss[ index_of_constraint( p_const ) ];

   for( auto el : lf->get_v_var() ) {
    SCIP_VAR * var = vars[ index_of_variable( el.first ) ];
    SCIP_CALL_ABORT( SCIPchgCoefLinear( scip, cons, var, el.second ) );
   }
  }
 }
}

/*--------------------------------------------------------------------------*/

void SCIPMILPSolver::function_vars_modification( FunctionModVars * mod ) {
 MILPSolver::function_vars_modification( mod );

 /*
 * This function is used when adding coefficents to OFs or constraints.
 */

 // Check if OF or a Constraint is involved
 auto * mod_f = mod->function();
 bool changing_of = false;
 const auto * lf = dynamic_cast<const LinearFunction *> (mod_f);
 const auto * qf = dynamic_cast<const DQuadFunction *> (mod_f);

 auto * p_obj = dynamic_cast< FRealObjective * >( f_Block->get_objective() );
 if( p_obj != nullptr ) {
  auto * of = p_obj->get_function();
  if( of == mod_f ) {
   changing_of = true;
  }
 }

 // C05FunctionModVarsAddd (Add coefficients)
 // --------------------------------------------------------
 auto * add = dynamic_cast<C05FunctionModVarsAddd *>( mod );
 if( add ) {

 }


 // C05FunctionModVarsRngd (Remove coefficients in a range)
 // ---------------------------------------------------------
 auto * rmvr = dynamic_cast<C05FunctionModVarsRngd *>( mod );
 if( rmvr ) {

 }


 // C05FunctionModVarsRngd (Remove coefficients in a subset)
 // ---------------------------------------------------------
 auto * rmvs = dynamic_cast<C05FunctionModVarsSbst *>( mod );
 if( rmvs ) {

 }


 throw std::invalid_argument( "This type of FunctionModVars is not handled" );
}

/*--------------------------------------------------------------------------*/

void SCIPMILPSolver::dynamic_modification( BlockModAD * mod ) {
 MILPSolver::dynamic_modification( mod );
}

/*--------------------------------------------------------------------------*/

void SCIPMILPSolver::add_dynamic_constraint( FRowConstraint * p_const ) {
 MILPSolver::add_dynamic_constraint( p_const );

 const auto * p_fun =
  dynamic_cast<const LinearFunction *>(p_const->get_function());
 if( p_fun == nullptr ) {
  throw std::invalid_argument( "The Constraint is not linear" );
 }

 if( SCIPisTransformed( scip ) )
  SCIP_CALL_ABORT( SCIPfreeTransform( scip ) );

 int nzcnt = p_const->get_num_active_var();

 SCIP_CONS * cons = nullptr;

 SCIP_Real const_lhs =
  p_const->get_lhs() == -Inf< double >() ?
  -SCIPinfinity( scip ) : p_const->get_lhs();
 SCIP_Real const_rhs =
  p_const->get_rhs() == Inf< double >() ?
  SCIPinfinity( scip ) : p_const->get_rhs();

 char name[32];
 std::snprintf(name, sizeof(name), "%p", (void*)p_const);
 SCIP_CALL_ABORT( SCIPcreateConsBasicLinear( scip, &cons, name, 0,
                                             nullptr, nullptr,
                                             const_lhs, const_rhs ) );

 // Get the coefficients to fill the matrix
 for( int i = 0; i < nzcnt; ++i ) {
  auto * p_var = dynamic_cast<ColVariable *>(p_fun->get_active_var( i ));
  SCIP_VAR * var = vars[ index_of_variable( p_var ) ];
  SCIP_Real coef = p_fun->get_coefficient( i );
  SCIP_CALL_ABORT( SCIPaddCoefLinear( scip, cons, var, coef ) );
 }

 SCIP_CALL_ABORT( SCIPaddCons( scip, cons ) );
 assert( conss.size() == dcon_to_idx.back().second + 1 );
 conss.push_back( cons );
 SCIP_CALL_ABORT( SCIPreleaseCons( scip, &cons ) );
}

/*--------------------------------------------------------------------------*/

void SCIPMILPSolver::add_dynamic_variable( ColVariable * p_var ) {
 MILPSolver::add_dynamic_variable( p_var );

 if( SCIPisTransformed( scip ) )
  SCIP_CALL_ABORT( SCIPfreeTransform( scip ) );

 SCIP_Real lb = NAN;
 SCIP_Real ub = NAN;
 SCIP_VARTYPE vartype = SCIP_VARTYPE_BINARY;

 lb =
  p_var->get_lb() == -Inf< double >() ?
  -SCIPinfinity( scip ) : p_var->get_lb();
 ub = p_var->get_ub() == Inf< double >() ?
      SCIPinfinity( scip ) : p_var->get_ub();

 // Variable type

 if( p_var->is_integer() ) {
  // ++mip;
  if( p_var->is_unitary() && p_var->is_positive() ) {
   vartype = SCIP_VARTYPE_BINARY;
  } else {
   vartype = SCIP_VARTYPE_INTEGER;
  }
 } else {
  vartype = SCIP_VARTYPE_CONTINUOUS;
 }

 SCIP_VAR * var = nullptr;

 SCIP_CALL_ABORT( SCIPcreateVarBasic( scip, &var, nullptr,
                                      lb, ub, 0.0, vartype ) );
 vars.push_back( var );
 SCIP_CALL_ABORT( SCIPaddVar( scip, var ) );

 int i = 0;

 // We need the coefficients for this variable in each constraint
 for( auto * p_const : active_constraints.back() ) {

  const auto * p_fun =
   dynamic_cast<const LinearFunction *>(p_const->get_function());
  if( p_fun == nullptr ) {
   throw ( std::invalid_argument( "The Constraint is not linear" ) );
  }
  SCIP_CONS * cons = conss[ index_of_constraint( p_const ) ];
  SCIP_Real coeff = p_fun->get_coefficient( i );
  SCIP_CALL_ABORT( SCIPaddCoefLinear( scip, cons, var, coeff ) );
  ++i;
 }

 SCIP_CALL_ABORT( SCIPreleaseVar( scip, &var ) );
}

/*--------------------------------------------------------------------------*/

void
SCIPMILPSolver::add_dynamic_bound( OneVarConstraint * p_bound ) {
 MILPSolver::add_dynamic_bound( p_bound );
 // TODO
}

/*--------------------------------------------------------------------------*/

void
SCIPMILPSolver::remove_dynamic_constraint( const FRowConstraint * p_const ) {
 if( SCIPisTransformed( scip ) )
  SCIP_CALL_ABORT( SCIPfreeTransform( scip ) );

 int index = index_of_dynamic_constraint( p_const );

 auto it = conss.begin() + index;
 SCIP_CONS * cons = *it;
 conss.erase( it );
 SCIP_CALL_ABORT( SCIPdelCons( scip, cons ) );

 MILPSolver::remove_dynamic_constraint( p_const );
}

/*--------------------------------------------------------------------------*/

void SCIPMILPSolver::remove_dynamic_variable( const ColVariable * p_var ) {
 if( SCIPisTransformed( scip ) )
  SCIP_CALL_ABORT( SCIPfreeTransform( scip ) );

 int index = index_of_dynamic_variable( p_var );

 auto it = vars.begin() + index;
 SCIP_VAR * var = *it;
 vars.erase( it );

 SCIP_Bool deleted = 0;
 SCIP_CALL_ABORT( SCIPdelVar( scip, var, &deleted ) );
 assert( deleted );

 MILPSolver::remove_dynamic_variable( p_var );
}

/*--------------------------------------------------------------------------*/

void
SCIPMILPSolver::remove_dynamic_bound( const OneVarConstraint * p_bound ) {
 MILPSolver::remove_dynamic_bound( p_bound );
 // TODO
}

/*--------------------------------------------------------------------------*/
/*------------------- METHODS FOR HANDLING THE PARAMETERS ------------------*/
/*--------------------------------------------------------------------------*/

void SCIPMILPSolver::set_par( const idx_type par, const int value ) {
 switch( par ) {
  // case intMaxIter:
  //  SCIP_CALL_ABORT( SCIPsetLongintParam( scip, "limits/nodes", value ) );
  //  break;
  // case intMaxSol:
  //  SCIP_CALL_ABORT( SCIPsetIntParam( scip, "limits/solutions", value ) );
  //  break;
  // case intLogVerb:
  //  SCIP_CALL_ABORT( SCIPsetIntParam( scip, "display/verblevel", value ) );
  //  break;
  case intUseCustomNames:
   use_custom_names = bool( value );
   break;
  default:
   MILPSolver::set_par( par, value );
 }
}

/*--------------------------------------------------------------------------*/

void SCIPMILPSolver::set_par( idx_type par, const double value ) {
 switch( par ) {
  // case dblMaxTime:
  //  SCIP_CALL_ABORT( SCIPsetRealParam( scip, "limits/time", value ) );
  //  break;
  // case dblRelAcc: // TODO
  //  break;
  // case dblAbsAcc: // TODO
  //  break;
  // case dblUpCutOff:
  //  if( objsense == 1 )
  //   SCIP_CALL_ABORT( SCIPsetObjlimit( scip, value ) );
  //  break;
  // case dblLwCutOff:
  //  if( objsense == -1 )
  //   SCIP_CALL_ABORT( SCIPsetObjlimit( scip, value ) );
  //  break;
  // case dblRAccSol:
  //  SCIP_CALL_ABORT( SCIPsetRealParam( scip, "limits/gap", value ) );
  //  break;
  // case dblAAccSol:
  //  SCIP_CALL_ABORT( SCIPsetRealParam( scip, "limits/absgap", value ) );
  //  break;
  // case dblFAccSol:
  //  SCIP_CALL_ABORT( SCIPsetRealParam( scip, "numerics/feastol", value ) );
  //  break;
  default:
   MILPSolver::set_par( par, value );
 }
}

/*--------------------------------------------------------------------------*/

void SCIPMILPSolver::set_par( idx_type par, const std::string & value ) {
 switch( par ) {
  case strProblemName:
   prob_name = value;
   break;
  case strOutputFile:
   output_file = value;
   break;
  default:
   MILPSolver::set_par( par, value );
 }
}

/*--------------------------------------------------------------------------*/

ThinComputeInterface::idx_type SCIPMILPSolver::get_num_int_par() const {
 return MILPSolver::get_num_int_par() + intLastAlgParSCPS - intLastAlgParMILP;
}

ThinComputeInterface::idx_type SCIPMILPSolver::get_num_str_par() const {
 return MILPSolver::get_num_str_par() + strLastAlgParSCPS - strLastAlgParMILP;
}

ThinComputeInterface::idx_type SCIPMILPSolver::get_num_dbl_par() const {
 return MILPSolver::get_num_dbl_par() + dblLastAlgParSCPS - dblLastAlgParMILP;
}

/*--------------------------------------------------------------------------*/

int SCIPMILPSolver::get_int_par( idx_type par ) const {
 int value;
 SCIP_Longint long_value;

 switch( par ) {
  // case intMaxIter:
  //  SCIP_CALL_ABORT( SCIPgetLongintParam( scip, "limits/nodes", &long_value ) );
  //  return ( int ) long_value;
  // case intMaxSol:
  //  SCIP_CALL_ABORT( SCIPgetIntParam( scip, "limits/solutions", &value ) );
  //  return value;
  // case intLogVerb:
  //  SCIP_CALL_ABORT( SCIPgetIntParam( scip, "display/verblevel", &value ) );
  //  return value;
  case intUseCustomNames:
   return use_custom_names;
  default:
   return MILPSolver::get_int_par( par );
 }
}

/*--------------------------------------------------------------------------*/

double SCIPMILPSolver::get_dbl_par( idx_type par ) const {
 double value;

 switch( par ) {
  // case dblMaxTime:
  //  SCIP_CALL_ABORT( SCIPgetRealParam( scip, "limits/time", &value ) );
  //  return value;
  // case dblRelAcc:   // TODO
  //  return 1e-6;
  // case dblAbsAcc:   // TODO
  //  return Inf< OFValue >();
  // case dblUpCutOff: // TODO
  //  if( objsense == 1 )
  //   return SCIPgetObjlimit( scip );
  //  return Inf< OFValue >();
  // case dblLwCutOff: // TODO
  //  if( objsense == -1 )
  //   return SCIPgetObjlimit( scip );
  //  return -Inf< OFValue >();
  // case dblRAccSol:
  //  SCIP_CALL_ABORT( SCIPgetRealParam( scip, "limits/gap", &value ) );
  //  return value;
  // case dblAAccSol:
  //  SCIP_CALL_ABORT( SCIPgetRealParam( scip, "limits/absgap", &value ) );
  //  return value;
  // case dblFAccSol:
  //  SCIP_CALL_ABORT( SCIPgetRealParam( scip, "numerics/feastol", &value ) );
  //  return value;
  default:
   return MILPSolver::get_dbl_par( par );
 }
}

/*--------------------------------------------------------------------------*/

const std::string &
SCIPMILPSolver::get_str_par( const idx_type par ) const {
 switch( par ) {
  case strProblemName:
   return prob_name;
  case strOutputFile:
   return output_file;
  default:
   return MILPSolver::get_str_par( par );
 }
}

/*--------------------------------------------------------------------------*/

int SCIPMILPSolver::get_dflt_int_par( const idx_type par ) const {
 int value;
 SCIP_Longint long_value;

 switch( par ) {
  // case intMaxIter:
  //  return ( int ) long_value;
  // case intMaxSol:
  //  return value;
  // case intLogVerb:
  //  return value;
  case intUseCustomNames:
   return 1;
  default:
   return MILPSolver::get_dflt_int_par( par );
 }
}

/*--------------------------------------------------------------------------*/

double SCIPMILPSolver::get_dflt_dbl_par( const idx_type par ) const {
 double value;

 switch( par ) {
  // case dblMaxTime:
  //  return value;
  // case dblRelAcc:   // TODO
  //  return 1e-6;
  // case dblAbsAcc:   // TODO
  //  return Inf< OFValue >();
  // case dblUpCutOff: // TODO
  //  if( objsense == 1 )
  //   return SCIPgetObjlimit( scip );
  //  return Inf< OFValue >();
  // case dblLwCutOff: // TODO
  //  if( objsense == -1 )
  //   return SCIPgetObjlimit( scip );
  //  return -Inf< OFValue >();
  // case dblRAccSol:
  //  return value;
  // case dblAAccSol:
  //  return value;
  // case dblFAccSol:
  //  return value;
  default:
   return MILPSolver::get_dflt_dbl_par( par );
 }
}

/*--------------------------------------------------------------------------*/

const std::string &
SCIPMILPSolver::get_dflt_str_par( const idx_type par ) const {

 if( par == strProblemName ) {
  return std::move( std::string( "SCIPMILPSolver_prob" ) );
 }

 if( par == strOutputFile ) {
  return std::move( std::string( "output.lp" ) );
 }

 return MILPSolver::get_dflt_str_par( par );
}

/*--------------------------------------------------------------------------*/

ThinComputeInterface::idx_type
SCIPMILPSolver::int_par_str2idx( const std::string & name ) const {
 if( name == "intUseCustomNames" )
  return intUseCustomNames;
 return MILPSolver::int_par_str2idx( name );
}

/*--------------------------------------------------------------------------*/

const std::string &
SCIPMILPSolver::int_par_idx2str( const idx_type idx ) const {
 if( idx == intUseCustomNames ) {
  return std::move( std::string( "intUseCustomNames" ) );
 }

 return MILPSolver::int_par_idx2str( idx );
}

/*--------------------------------------------------------------------------*/

ThinComputeInterface::idx_type
SCIPMILPSolver::dbl_par_str2idx( const std::string & name ) const {
 return MILPSolver::dbl_par_str2idx( name );
}

/*--------------------------------------------------------------------------*/

const std::string &
SCIPMILPSolver::dbl_par_idx2str( const idx_type idx ) const {
 return MILPSolver::dbl_par_idx2str( idx );
}

/*--------------------------------------------------------------------------*/


ThinComputeInterface::idx_type
SCIPMILPSolver::str_par_str2idx( const std::string & name ) const {
 if( name == "strProblemName" )
  return strProblemName;
 if( name == "strOutputFile" )
  return strOutputFile;
 return MILPSolver::str_par_str2idx( name );
}

/*--------------------------------------------------------------------------*/

const std::string &
SCIPMILPSolver::str_par_idx2str( const idx_type idx ) const {
 static const std::vector< std::string > pars = { "strProblemName",
                                                  "strOutputFile" };
 switch( idx ) {
  case strProblemName:
   return pars[ 0 ];
  case strOutputFile:
   return pars[ 1 ];
  default:
   return MILPSolver::str_par_idx2str( idx );
 }
}

/*--------------------------------------------------------------------------*/
/*--------------------- PRIVATE FIELDS OF THE CLASS ------------------------*/
/*--------------------------------------------------------------------------*/

void SCIPMILPSolver::set_var_value( ColVariable & lvar, double * x, int & i ) {

 lvar.set_value( x[ i++ ] );
}

/*--------------------------------------------------------------------------*/

void SCIPMILPSolver::set_dual_value( FRowConstraint & lconst,
                                     double * pi,
                                     int & i ) {

 lconst.set_dual( pi[ i++ ] );
}

/*--------------------------------------------------------------------------*/

// void SCIPMILPSolver::fix_integer_vars() {
//  int probtype = CPXgetprobtype( env, milp );
//  switch( probtype ) {
//   case CPXPROB_MILP:
//    probtype = CPXPROB_FIXEDMILP;
//    break;
//   case CPXPROB_MIQP:
//    probtype = CPXPROB_FIXEDMIQP;
//    break;
//   default:
//    throw std::runtime_error( "Wrong problem type from CPXgetprobtype()" );
//  }
//  int status = CPXchgprobtype( env, milp, probtype );
//  if( status ) {
//   throw std::runtime_error( "Unable to change problem type with CPXchgprobtype()" );
//  }
//
//  status = CPXprimopt( env, milp );
//  if( status ) {
//   throw std::runtime_error( "An error occurred in CPXprimopt()" );
//  }
// }
/*--------------------------------------------------------------------------*/
/*--------------------- End File SCIPMILPSolver.cpp -------------------------*/
/*--------------------------------------------------------------------------*/
