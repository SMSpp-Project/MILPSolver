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

// Include the proper SCIP parameter mapping
#include <boost/preprocessor/cat.hpp>
#include <boost/preprocessor/stringize.hpp>
#include BOOST_PP_STRINGIZE( BOOST_PP_CAT( BOOST_PP_CAT( SCIP, SCIP_VERSION ), _maps.h ) )

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

void SCIPMILPSolver::clear_problem( unsigned int what ) {
 MILPSolver::clear_problem( what );

 SCIP_CALL_ABORT( SCIPfreeProb( scip ) );
}

/*--------------------------------------------------------------------------*/

void SCIPMILPSolver::load_problem() {
 MILPSolver::load_problem();

 SCIP_CALL_ABORT( SCIPcreateProbBasic( scip, prob_name.c_str() ) );

 // Set objective sense
 if( objsense == 1 ) {
  SCIP_CALL_ABORT( SCIPsetObjsense( scip, SCIP_OBJSENSE_MINIMIZE ) );
 } else if( objsense == -1 ) {
  SCIP_CALL_ABORT( SCIPsetObjsense( scip, SCIP_OBJSENSE_MAXIMIZE ) );
 } else {
  SCIPABORT();
 }

 // Add variables
 vars.resize( numcols );
 for( int i = 0; i < numcols; ++i ) {

  SCIP_Real collb = ( lb[ i ] == -Inf< double >() ) ?
                    -SCIPinfinity( scip ) : lb[ i ];
  SCIP_Real colub = ( ub[ i ] == Inf< double >() ) ?
                    SCIPinfinity( scip ) : ub[ i ];

  SCIP_VARTYPE vartype;
  switch( xctype[ i ] ) {
   case 'C':
    vartype = SCIP_VARTYPE_CONTINUOUS;
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
   default:
    SCIPABORT();
  }


  SCIP_VAR * var = nullptr;
  char * name = use_custom_names ? colname[ i ] : nullptr;
  SCIP_CALL_ABORT( SCIPcreateVarBasic( scip, &var, name, collb, colub,
                                       objective[ i ],
                                       vartype ) );
  SCIP_CALL_ABORT( SCIPaddVar( scip, var ) );
  vars[ i ] = var;
  SCIP_CALL_ABORT( SCIPreleaseVar( scip, &var ) );
 }

 // Add constraints
 cons.resize( numrows );
 for( int i = 0; i < numrows; ++i ) {
  SCIP_Real con_lhs = NAN;
  SCIP_Real con_rhs = NAN;

  switch( sense[ i ] ) {
   case 'L':
    con_lhs = -SCIPinfinity( scip );
    con_rhs = rhs[ i ];
    break;
   case 'E':
    con_lhs = rhs[ i ];
    con_rhs = rhs[ i ];
    break;
   case 'G':
    con_lhs = rhs[ i ];
    con_rhs = SCIPinfinity( scip );
    break;
   case 'R':
    con_lhs = rhs[ i ];
    con_rhs = rhs[ i ];
    // TODO Check this
    if( rngval[ i ] > 0 )
     con_rhs += rngval[ i ];
    else
     con_lhs += rngval[ i ];
  }

  SCIP_CONS * con = nullptr;
  char * name = use_custom_names ? rowname[ i ] : nullptr;
  SCIP_CALL_ABORT( SCIPcreateConsBasicLinear( scip, &con, name, 0,
                                              nullptr, nullptr,
                                              con_lhs, con_rhs ) );
  SCIP_CALL_ABORT( SCIPaddCons( scip, con ) );
  cons[ i ] = con;
  SCIP_CALL_ABORT( SCIPreleaseCons( scip, &con ) );
 }

 // Add constraint coefficients
 for( int c = 0; c < numcols; ++c ) {
  for( int i = matbeg[ c ]; i < matbeg[ c + 1 ]; ++i ) {
   SCIP_CALL_ABORT( SCIPaddCoefLinear( scip,
                                       cons[ matind[ i ] ],
                                       vars[ c ],
                                       matval[ i ] ) );
  }
 }

 bool is_qp = std::any_of( q_objective.begin(),
                           q_objective.end(),
                           []( double d ) { return d != 0; } );
 if( is_qp ) {
  // For each col add create a new aux_var z.
  // In the objective z goes with the q_objective value.
  // Then add a new quad aux_constraint z - xˆ2 >= 0

  aux_vars.resize( numcols );
  aux_cons.resize( numcols );

  for( int i = 0; i < numcols; ++i ) {

   // Add auxiliary variables
   SCIP_Real z_lb = ( lb[ i ] == -Inf< double >() ) ?
                    -SCIPinfinity( scip ) : lb[ i ];
   SCIP_Real z_ub = ( ub[ i ] == Inf< double >() ) ?
                    SCIPinfinity( scip ) : ub[ i ];

   SCIP_VARTYPE z_type;
   switch( xctype[ i ] ) {
    case 'C':
     z_type = SCIP_VARTYPE_CONTINUOUS;
     break;
    case 'B':
     z_lb = 0;
     z_ub = 1;
     z_type = SCIP_VARTYPE_BINARY;
     break;
    case 'I':
     z_type = SCIP_VARTYPE_INTEGER;
     break;
    case 'S':
    case 'N':
    default:
     SCIPABORT();
   }

   SCIP_VAR * z = nullptr;
   SCIP_CALL_ABORT( SCIPcreateVarBasic( scip, &z, nullptr, z_lb, z_ub,
                                        q_objective[ i ], z_type ) );
   SCIP_CALL_ABORT( SCIPaddVar( scip, z ) );
   aux_vars[ i ] = z;
   SCIP_CALL_ABORT( SCIPreleaseVar( scip, &z ) );


   // Add auxiliary constraints
   SCIP_Real con_lhs = NAN;
   SCIP_Real con_rhs = NAN;

   if( SCIPgetObjsense( scip ) == SCIP_OBJSENSE_MINIMIZE ) {
    con_lhs = 0;
    con_rhs = SCIPinfinity( scip );
   } else if( SCIPgetObjsense( scip ) == SCIP_OBJSENSE_MAXIMIZE ) {
    con_lhs = -SCIPinfinity( scip );
    con_rhs = 0;
   } else {
    SCIPABORT();
   }

   SCIP_CONS * con = nullptr;
   SCIP_CALL_ABORT( SCIPcreateConsBasicQuadratic( scip, &con, nullptr, 0,
                                                  nullptr, nullptr, 0,
                                                  nullptr, nullptr, nullptr,
                                                  con_lhs, con_rhs ) );
   SCIP_CALL_ABORT( SCIPaddCons( scip, con ) );
   aux_cons[ i ] = con;
   SCIP_CALL_ABORT( SCIPreleaseCons( scip, &con ) );


   // Add constraint coefficients
   SCIP_CALL_ABORT( SCIPaddCoefLinear( scip,
                                       aux_cons[ i ],
                                       aux_vars[ i ],
                                       1 ) );
   SCIP_CALL_ABORT( SCIPaddQuadVarQuadratic( scip,
                                             aux_cons[ i ],
                                             aux_vars[ i ],
                                             0,
                                             -1 ) );
  }
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
  b = SCIPinfinity( scip );
 }
 return b;
}

/*--------------------------------------------------------------------------*/

int SCIPMILPSolver::compute( bool changedvars ) {
 if( MILPSolver::compute( changedvars ) != kOK ) {
  // This should never happen
  throw std::runtime_error( "An error occurred in MILPSolver::compute()" );
 }

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

 return sol_status;
}

/*--------------------------------------------------------------------------*/

Solver::OFValue SCIPMILPSolver::get_lb() {

 OFValue lower_bound = 0;

 switch( SCIPgetObjsense( scip ) ) {
  case SCIP_OBJSENSE_MINIMIZE:
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
  case SCIP_OBJSENSE_MAXIMIZE:
   switch( sol_status ) {
    case kUnbounded:
     lower_bound = Inf< OFValue >();
     break;
    case kInfeasible:
     lower_bound = -Inf< OFValue >();
     break;
    default:
     lower_bound = SCIPgetPrimalbound( scip );
     break;
   }
   break;
  default:
   throw std::runtime_error( "Objective type not yet defined" );
 }

 return lower_bound;
}

/*--------------------------------------------------------------------------*/

Solver::OFValue SCIPMILPSolver::get_ub() {

 OFValue upper_bound = 0;

 switch( SCIPgetObjsense( scip ) ) {
  case SCIP_OBJSENSE_MINIMIZE:
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
  case SCIP_OBJSENSE_MAXIMIZE:
   switch( sol_status ) {
    case kUnbounded:
     upper_bound = Inf< OFValue >();
     break;
    case kInfeasible:
     upper_bound = -Inf< OFValue >();
     break;
    default:
     upper_bound = SCIPgetDualbound( scip );
     break;
   }
   break;
  default:
   throw std::runtime_error( "Objective type not yet defined" );
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
 switch( SCIPgetObjsense( scip ) ) {
  case SCIP_OBJSENSE_MINIMIZE:
   return get_ub();
  case SCIP_OBJSENSE_MAXIMIZE:
   return get_lb();
  default:
   throw std::runtime_error( "Objective type not yet defined" );
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

  auto set = [ &x, &col ]( ColVariable & v ) {
   v.set_value( x[ col++ ] );
  };

  for( const auto & i : q_Block->get_static_variables() ) {
   un_any_const_static( i, set, un_any_type< ColVariable >() );
  }

  for( const auto & i : q_Block->get_dynamic_variables() ) {
   un_any_const_dynamic( i, set, un_any_type< ColVariable >() );
  }
 }

 if( !owned ) {
  f_Block->unlock( f_id );
 }

 delete[]x;
}

/*--------------------------------------------------------------------------*/

bool SCIPMILPSolver::has_dual_solution() {
 // TODO
 return false;
}

/*--------------------------------------------------------------------------*/

bool SCIPMILPSolver::is_dual_feasible() {
 // TODO
 return false;
}

/*--------------------------------------------------------------------------*/

void SCIPMILPSolver::get_dual_solution( Configuration * solc ) {
 // TODO
 SCIPABORT();
}

/*--------------------------------------------------------------------------*/

bool SCIPMILPSolver::has_dual_direction() {
 // TODO
 return false;
}

/*--------------------------------------------------------------------------*/

void SCIPMILPSolver::get_dual_direction( Configuration * dirc ) {
 // TODO
 SCIPABORT();
}

/*--------------------------------------------------------------------------*/

void SCIPMILPSolver::write_lp( const std::string & filename ) {
 SCIP_CALL_ABORT( SCIPwriteOrigProblem( scip,
                                        filename.c_str(),
                                        "mps", FALSE ) );
}

int SCIPMILPSolver::get_nodes() const {
 return static_cast<int>(SCIPgetNTotalNodes( scip ));
}

/*--------------------------------------------------------------------------*/
/*-------------------- PROTECTED FIELDS OF THE CLASS -----------------------*/
/*--------------------------------------------------------------------------*/

void SCIPMILPSolver::var_modification( VariableMod * mod ) {
 try {
  MILPSolver::var_modification( mod );
 } catch( std::logic_error & e ) {}

 if( SCIPisTransformed( scip ) ) {
  SCIP_CALL_ABORT( SCIPfreeTransform( scip ) );
 }

 auto * var = dynamic_cast< ColVariable * >( mod->variable() );
 // TODO: Is dynamic_cast necessary?
 if( var == nullptr ) {
  // TODO: Throw exception?
  return;
 }

 int idx = index_of_variable( var );

 // Read old variable type
 SCIP_VARTYPE oldtype = SCIPvarGetType( vars[ idx ] );

 SCIP_Bool infeas = FALSE;

 if( var->is_fixed() ) {
  // Fix the variable
  SCIP_Bool fixed = 0;
  SCIP_CALL_ABORT(
   SCIPfixVar( scip, vars[ idx ], var->get_value(), &infeas, &fixed )
  );
  assert( fixed );

 } else if( var->is_integer() ) {
  // Unfix or refresh the variable
  // TODO: Check this, it seems unnecessary
  if( var->is_unitary() && var->is_positive() ) {
   if( oldtype != SCIP_VARTYPE_BINARY ) {
    SCIP_CALL_ABORT(
     SCIPchgVarType( scip, vars[ idx ], SCIP_VARTYPE_BINARY, &infeas )
    );
   }
  } else if( oldtype != SCIP_VARTYPE_INTEGER ) {
   SCIP_CALL_ABORT(
    SCIPchgVarType( scip, vars[ idx ], SCIP_VARTYPE_INTEGER, &infeas )
   );
  }
 } else if( oldtype != SCIP_VARTYPE_CONTINUOUS ) {
  SCIP_CALL_ABORT(
   SCIPchgVarType( scip, vars[ idx ], SCIP_VARTYPE_CONTINUOUS, &infeas )
  );
 }
 assert( !infeas );

 SCIP_Real lb = get_problem_lb( *var );
 SCIP_Real ub = get_problem_ub( *var );

 if( SCIPvarGetLbOriginal( vars[ idx ] ) != lb ) {
  SCIP_CALL_ABORT( SCIPchgVarLb( scip, vars[ idx ], lb ) );
 }

 if( SCIPvarGetUbOriginal( vars[ idx ] ) != ub ) {
  SCIP_CALL_ABORT( SCIPchgVarUb( scip, vars[ idx ], ub ) );
 }
}

/*--------------------------------------------------------------------------*/

void SCIPMILPSolver::objective_modification( ObjectiveMod * mod ) {
 try {
  MILPSolver::objective_modification( mod );
 } catch( std::logic_error & e ) {}

 if( SCIPisTransformed( scip ) ) {
  SCIP_CALL_ABORT( SCIPfreeTransform( scip ) );
 }
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
 try {
  MILPSolver::const_modification( mod );
 } catch( std::logic_error & e ) {}

 if( SCIPisTransformed( scip ) ) {
  SCIP_CALL_ABORT( SCIPfreeTransform( scip ) );
 }

 auto * con = dynamic_cast<FRowConstraint *>(mod->constraint());
 // TODO: Is dynamic_cast necessary?
 if( con == nullptr ) {
  // TODO: Throw exception?
  return;
 }

 RowConstraint::RHSValue con_lhs = NAN;
 RowConstraint::RHSValue con_rhs = NAN;

 SCIP_CONS * scip_con = cons[ index_of_constraint( con ) ];

 switch( mod->type() ) {

  case ConstraintMod::eRelaxConst:
   // In order to relax the constraint all we do is transform it
   // into an inequality with RHS equal to infinity
   SCIP_CALL_ABORT( SCIPchgLhsLinear( scip, scip_con, -SCIPinfinity( scip ) ) );
   SCIP_CALL_ABORT( SCIPchgRhsLinear( scip, scip_con, SCIPinfinity( scip ) ) );
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

   con_lhs = con->get_lhs() == -Inf< double >() ?
             -SCIPinfinity( scip ) : con->get_lhs();
   con_rhs = con->get_rhs() == Inf< double >() ?
             SCIPinfinity( scip ) : con->get_rhs();

   SCIP_CALL_ABORT( SCIPchgLhsLinear( scip, scip_con, con_lhs ) );
   SCIP_CALL_ABORT( SCIPchgRhsLinear( scip, scip_con, con_rhs ) );

   break;

  default:
   throw std::invalid_argument( "Invalid type of ObjectiveMod" );
 }
}

/*--------------------------------------------------------------------------*/

void SCIPMILPSolver::bound_modification( OneVarConstraintMod * mod ) {
 try {
  MILPSolver::bound_modification( mod );
 } catch( std::logic_error & e ) {}

 /*
  * The same ColVariable can have more active OneVarConstraints,
  * so each time we modify one of them we have to check if LHS and RHS
  * of the Variable change.
  */

 if( SCIPisTransformed( scip ) ) {
  SCIP_CALL_ABORT( SCIPfreeTransform( scip ) );
 }

 auto * con = static_cast<OneVarConstraint *>(mod->constraint());
 auto * var = static_cast<ColVariable *>(con->get_active_var( 0 ));

 SCIP_Real lb;
 SCIP_Real ub;
 SCIP_VAR * scip_var = vars[ index_of_variable( var ) ];

 switch( mod->type() ) {

  case RowConstraintMod::eChgLHS:
   lb = get_problem_lb( *var );
   SCIP_CALL_ABORT( SCIPchgVarLb( scip, scip_var, lb ) );
   break;

  case RowConstraintMod::eChgRHS:
   ub = get_problem_ub( *var );
   SCIP_CALL_ABORT( SCIPchgVarUb( scip, scip_var, ub ) );
   break;

  case RowConstraintMod::eChgBTS:
   lb = get_problem_lb( *var );
   ub = get_problem_ub( *var );
   SCIP_CALL_ABORT( SCIPchgVarLb( scip, scip_var, lb ) );
   SCIP_CALL_ABORT( SCIPchgVarUb( scip, scip_var, ub ) );
   break;

  default:
   throw std::invalid_argument( "Invalid type of OneVarConstraintMod" );
 }
}

/*--------------------------------------------------------------------------*/

void SCIPMILPSolver::objective_function_modification( FunctionMod * mod ) {
 try {
  MILPSolver::objective_function_modification( mod );
 } catch( std::logic_error & e ) {}

 if( SCIPisTransformed( scip ) ) {
  SCIP_CALL_ABORT( SCIPfreeTransform( scip ) );
 }

 auto * f = mod->function();

 // C05FunctionModLin
 // --------------------------------------------------------------------------

 // Fallback method - Update all costs
 // --------------------------------------------------------------------------
 if( const auto * lf = dynamic_cast<const LinearFunction *> (f) ) {
  // Linear objective function
  for( auto el : lf->get_v_var() ) {
   SCIP_VAR * var = vars[ index_of_variable( el.first ) ];
   SCIP_CALL_ABORT( SCIPchgVarObj( scip, var, el.second ) );
  }
  return;
 }

 if( const auto * qf = dynamic_cast<const DQuadFunction *> (f) ) {
  // Quadratic objective function
  // TODO
  SCIPABORT();
  return;
 }

 // This should never happen
 throw std::invalid_argument( "Unknown type of Objective Function" );
}

/*--------------------------------------------------------------------------*/

// TODO: Change only involved variables
void SCIPMILPSolver::constraint_function_modification( FunctionMod * mod ) {
 try {
  MILPSolver::constraint_function_modification( mod );
 } catch( std::logic_error & e ) {}

 auto * f = mod->function();
 const auto * lf = dynamic_cast<const LinearFunction *> (f);
 if( lf == nullptr ) {
  return;
 }

 auto * con = dynamic_cast<FRowConstraint *>(lf->get_Observer());
 // TODO: Is dynamic_cast necessary?
 if( con == nullptr ) {
  // TODO: Throw exception?
  return;
 }
 SCIP_CONS * scip_con = cons[ index_of_constraint( con ) ];

 if( SCIPisTransformed( scip ) ) {
  SCIP_CALL_ABORT( SCIPfreeTransform( scip ) );
 }

 // C05FunctionModLin
 // --------------------------------------------------------------------------

 // Fallback method - Reload all coefficients
 // --------------------------------------------------------------------------
 for( auto el : lf->get_v_var() ) {
  SCIP_VAR * var = vars[ index_of_variable( el.first ) ];
  SCIP_CALL_ABORT( SCIPchgCoefLinear( scip, scip_con, var, el.second ) );
 }
}

/*--------------------------------------------------------------------------*/

void SCIPMILPSolver::objective_fvars_modification( FunctionModVars * mod ) {
 try {
  MILPSolver::objective_fvars_modification( mod );
 } catch( std::logic_error & e ) {}

 auto * f = mod->function();

 // Check the modification type
 if( dynamic_cast<C05FunctionModVarsAddd *>( mod ) == nullptr &&
     dynamic_cast<C05FunctionModVarsRngd *>( mod ) == nullptr &&
     dynamic_cast<C05FunctionModVarsSbst *>( mod ) == nullptr ) {
  throw std::invalid_argument( "This type of FunctionModVars is not handled" );
 }

 if( SCIPisTransformed( scip ) ) {
  SCIP_CALL_ABORT( SCIPfreeTransform( scip ) );
 }

 if( const auto * lf = dynamic_cast<const LinearFunction *> (f) ) {
  // Linear objective function

  for( auto * it1 : mod->vars() ) {
   for( auto it2: lf->get_v_var() ) {
    if( it1 == it2.first ) {
     SCIP_VAR * scip_var = vars[ index_of_variable( it2.first ) ];
     if( mod->added() ) {
      SCIP_CALL_ABORT( SCIPchgVarObj( scip, scip_var, it2.second ) );
     } else {
      SCIP_CALL_ABORT( SCIPchgVarObj( scip, scip_var, 0 ) );
     }
     break;
    }
   }
  }

  return;
 }

 if( const auto * qf = dynamic_cast<const DQuadFunction *> (f) ) {
  // TODO
  SCIPABORT();
  return;
 }

 // This should never happen
 throw std::invalid_argument( "Unknown type of Objective Function" );
}

/*--------------------------------------------------------------------------*/

void SCIPMILPSolver::constraint_fvars_modification( FunctionModVars * mod ) {
 try {
  MILPSolver::constraint_fvars_modification( mod );
 } catch( std::logic_error & e ) {}

 auto * f = mod->function();
 const auto * lf = dynamic_cast<const LinearFunction *> (f);
 if( lf == nullptr ) {
  return;
 }

 // Check the modification type
 if( dynamic_cast<C05FunctionModVarsAddd *>( mod ) == nullptr &&
     dynamic_cast<C05FunctionModVarsRngd *>( mod ) == nullptr &&
     dynamic_cast<C05FunctionModVarsSbst *>( mod ) == nullptr ) {
  throw std::invalid_argument( "This type of FunctionModVars is not handled" );
 }

 auto * con = dynamic_cast<FRowConstraint *>(lf->get_Observer());
 // TODO: Is dynamic_cast necessary?
 if( con == nullptr ) {
  // TODO: Throw exception?
  return;
 }
 SCIP_CONS * scip_con = cons[ index_of_constraint( con ) ];

 if( SCIPisTransformed( scip ) ) {
  SCIP_CALL_ABORT( SCIPfreeTransform( scip ) );
 }

 for( auto * it1 : mod->vars() ) {
  for( auto it2: lf->get_v_var() ) {
   if( it1 == it2.first ) {
    SCIP_VAR * scip_var = vars[ index_of_variable( it2.first ) ];
    if( mod->added() ) {
     SCIP_CALL_ABORT( SCIPchgCoefLinear( scip, scip_con, scip_var, it2
      .second ) );
    } else {
     SCIP_CALL_ABORT( SCIPchgCoefLinear( scip, scip_con, scip_var, 0 ) );
    }
    break;
   }
  }
 }
}

/*--------------------------------------------------------------------------*/

void SCIPMILPSolver::dynamic_modification( BlockModAD * mod ) {
 try {
  MILPSolver::dynamic_modification( mod );
 } catch( std::logic_error & e ) {}
}

/*--------------------------------------------------------------------------*/

void SCIPMILPSolver::add_dynamic_constraint( FRowConstraint * con ) {
 try {
  MILPSolver::add_dynamic_constraint( con );
 } catch( std::logic_error & e ) {}

 const auto * f = dynamic_cast<const LinearFunction *>(con->get_function());
 if( f == nullptr ) {
  throw std::invalid_argument( "The Constraint is not linear" );
 }

 if( SCIPisTransformed( scip ) ) {
  SCIP_CALL_ABORT( SCIPfreeTransform( scip ) );
 }

 SCIP_CONS * scip_con = nullptr;

 SCIP_Real con_lhs = con->get_lhs() == -Inf< double >() ?
                     -SCIPinfinity( scip ) : con->get_lhs();
 SCIP_Real con_rhs = con->get_rhs() == Inf< double >() ?
                     SCIPinfinity( scip ) : con->get_rhs();

 char name[32];
 std::snprintf( name, sizeof( name ), "%p", ( void * ) con );
 SCIP_CALL_ABORT( SCIPcreateConsBasicLinear( scip, &scip_con, name, 0,
                                             nullptr, nullptr,
                                             con_lhs, con_rhs ) );

 // Get the coefficients to fill the matrix
 for( Block::Index i = 0; i < con->get_num_active_var(); ++i ) {
  auto * var = dynamic_cast<ColVariable *>( f->get_active_var( i ) );
  // TODO: Is dynamic_cast necessary?
  if( var == nullptr ) {
   // TODO: Throw exception?
   continue;
  }
  SCIP_VAR * scip_var = vars[ index_of_variable( var ) ];
  SCIP_Real coef = f->get_coefficient( i );
  SCIP_CALL_ABORT( SCIPaddCoefLinear( scip, scip_con, scip_var, coef ) );
 }

 SCIP_CALL_ABORT( SCIPaddCons( scip, scip_con ) );
 cons.push_back( scip_con );
 SCIP_CALL_ABORT( SCIPreleaseCons( scip, &scip_con ) );
}

/*--------------------------------------------------------------------------*/

void SCIPMILPSolver::add_dynamic_variable( ColVariable * var ) {
 try {
  MILPSolver::add_dynamic_variable( var );
 } catch( std::logic_error & e ) {}

 if( SCIPisTransformed( scip ) ) {
  SCIP_CALL_ABORT( SCIPfreeTransform( scip ) );
 }

 SCIP_Real lb = get_problem_lb( *var );
 SCIP_Real ub = get_problem_ub( *var );
 SCIP_VARTYPE vartype = SCIP_VARTYPE_BINARY;

 // Variable type

 if( var->is_integer() ) {
  if( var->is_unitary() && var->is_positive() ) {
   vartype = SCIP_VARTYPE_BINARY;
  } else {
   vartype = SCIP_VARTYPE_INTEGER;
  }
 } else {
  vartype = SCIP_VARTYPE_CONTINUOUS;
 }

 SCIP_VAR * scip_var = nullptr;

 SCIP_CALL_ABORT( SCIPcreateVarBasic( scip, &scip_var, nullptr,
                                      lb, ub, 0.0, vartype ) );
 vars.push_back( scip_var );
 SCIP_CALL_ABORT( SCIPaddVar( scip, scip_var ) );

 int i = 0;

 // We need the coefficients for this variable in each constraint
 auto active_constraints = get_active_constraints( *var );
 for( auto * con : active_constraints ) {

  const auto * f = dynamic_cast<const LinearFunction *>(con->get_function());
  if( f == nullptr ) {
   throw ( std::invalid_argument( "The Constraint is not linear" ) );
  }
  SCIP_CONS * scip_con = cons[ index_of_constraint( con ) ];
  SCIP_Real coeff = f->get_coefficient( i );
  SCIP_CALL_ABORT( SCIPaddCoefLinear( scip, scip_con, scip_var, coeff ) );
  ++i;
 }

 SCIP_CALL_ABORT( SCIPreleaseVar( scip, &scip_var ) );
}

/*--------------------------------------------------------------------------*/

void
SCIPMILPSolver::add_dynamic_bound( OneVarConstraint * con ) {
 try {
  MILPSolver::add_dynamic_bound( con );
 } catch( std::logic_error & e ) {}

 if( SCIPisTransformed( scip ) ) {
  SCIP_CALL_ABORT( SCIPfreeTransform( scip ) );
 }

 auto * var = dynamic_cast<ColVariable *>(con->get_active_var( 0 ));
 // TODO: Is dynamic_cast necessary?
 if( var == nullptr ) {
  // TODO: Throw exception?
  return;
 }

 SCIP_VAR * scip_var = vars[ index_of_variable( var ) ];

 SCIP_Real lb = get_problem_lb( *var );
 SCIP_Real ub = get_problem_ub( *var );
 SCIP_CALL_ABORT( SCIPchgVarLb( scip, scip_var, lb ) );
 SCIP_CALL_ABORT( SCIPchgVarUb( scip, scip_var, ub ) );
}

/*--------------------------------------------------------------------------*/

void
SCIPMILPSolver::remove_dynamic_constraint( const FRowConstraint * con ) {
 if( SCIPisTransformed( scip ) ) {
  SCIP_CALL_ABORT( SCIPfreeTransform( scip ) );
 }

 int index = index_of_dynamic_constraint( con );

 auto it = cons.begin() + index;
 SCIP_CONS * scip_con = *it;
 cons.erase( it );
 SCIP_CALL_ABORT( SCIPdelCons( scip, scip_con ) );

 try {
  MILPSolver::remove_dynamic_constraint( con );
 } catch( std::logic_error & e ) {}
}

/*--------------------------------------------------------------------------*/

void SCIPMILPSolver::remove_dynamic_variable( const ColVariable * var ) {
 if( SCIPisTransformed( scip ) ) {
  SCIP_CALL_ABORT( SCIPfreeTransform( scip ) );
 }

 int index = index_of_dynamic_variable( var );

 auto it = vars.begin() + index;
 SCIP_VAR * scip_var = *it;
 vars.erase( it );

 SCIP_Bool deleted = 0;
 SCIP_CALL_ABORT( SCIPdelVar( scip, scip_var, &deleted ) );
 assert( deleted );

 try {
  MILPSolver::remove_dynamic_variable( var );
 } catch( std::logic_error & e ) {}
}

/*--------------------------------------------------------------------------*/

void
SCIPMILPSolver::remove_dynamic_bound( const OneVarConstraint * con ) {
 try {
  MILPSolver::remove_dynamic_bound( con );
 } catch( std::logic_error & e ) {}

 if( SCIPisTransformed( scip ) ) {
  SCIP_CALL_ABORT( SCIPfreeTransform( scip ) );
 }

 auto * var = dynamic_cast<ColVariable *>(con->get_active_var( 0 ));
 // TODO: Is dynamic_cast necessary?
 if( var == nullptr ) {
  // TODO: Throw exception?
  return;
 }

 SCIP_VAR * scip_var = vars[ index_of_variable( var ) ];

 SCIP_Real lb = get_problem_lb( *var );
 SCIP_Real ub = get_problem_ub( *var );
 SCIP_CALL_ABORT( SCIPchgVarLb( scip, scip_var, lb ) );
 SCIP_CALL_ABORT( SCIPchgVarUb( scip, scip_var, ub ) );
}

/*--------------------------------------------------------------------------*/
/*------------------- METHODS FOR HANDLING THE PARAMETERS ------------------*/
/*--------------------------------------------------------------------------*/

void SCIPMILPSolver::set_par( const idx_type par, const int value ) {

 // Solver parameters explicitly mapped in SCIP
 switch( par ) {
  case intMaxIter:
   SCIP_CALL_ABORT( SCIPsetLongintParam( scip, "limits/nodes", value ) );
   return;
  case intMaxSol:
   SCIP_CALL_ABORT( SCIPsetIntParam( scip, "limits/solutions", value ) );
   return;
  case intLogVerb:
   SCIP_CALL_ABORT( SCIPsetIntParam( scip, "display/verblevel", value ) );
   return;
  default:;
 }

 // SCIP parameters
 if( par >= intFirstSCIPPar && par < intLastAlgParSCPS ) {
  const std::string & scip_par =
   SMSpp_to_SCIP_int_pars[ par - intFirstSCIPPar ];

  // Bool, int and long SCIP parameters are handled as SMS++ int parameters
  SCIP_PARAM * param = SCIPgetParam( scip, scip_par.c_str() );
  SCIP_PARAMTYPE type = SCIPparamGetType( param );

  if( type == SCIP_PARAMTYPE_BOOL ) {
   SCIP_CALL_ABORT( SCIPsetBoolParam( scip, scip_par.c_str(), value ) );
  } else if( type == SCIP_PARAMTYPE_INT ) {
   SCIP_CALL_ABORT( SCIPsetIntParam( scip, scip_par.c_str(), value ) );
  } else if( type == SCIP_PARAMTYPE_LONGINT ) {
   SCIP_CALL_ABORT( SCIPsetLongintParam( scip, scip_par.c_str(), value ) );
  }
  return;
 }

 MILPSolver::set_par( par, value );
}

/*--------------------------------------------------------------------------*/

void SCIPMILPSolver::set_par( idx_type par, const double value ) {

 // Solver parameters explicitly mapped in SCIP
 switch( par ) {
  case dblMaxTime:
   SCIP_CALL_ABORT( SCIPsetRealParam( scip, "limits/time", value ) );
   return;
   // case dblRelAcc: // TODO
   //  return;
   // case dblAbsAcc: // TODO
   //  return;
   // case dblUpCutOff:
   //  if( SCIPgetObjsense( scip ) == SCIP_OBJSENSE_MINIMIZE ) {
   //   SCIP_CALL_ABORT( SCIPsetObjlimit( scip, value ) );
   //  }
   //  return;
   // case dblLwCutOff:
   //  if( SCIPgetObjsense( scip ) == SCIP_OBJSENSE_MAXIMIZE ) {
   //   SCIP_CALL_ABORT( SCIPsetObjlimit( scip, value ) );
   //  }
   //  return;
  case dblRAccSol:
   SCIP_CALL_ABORT( SCIPsetRealParam( scip, "limits/gap", value ) );
   return;
  case dblAAccSol:
   SCIP_CALL_ABORT( SCIPsetRealParam( scip, "limits/absgap", value ) );
   return;
  case dblFAccSol:
   SCIP_CALL_ABORT( SCIPsetRealParam( scip, "numerics/feastol", value ) );
   return;
  default:;
 }

 // SCIP parameters
 if( par >= dblFirstSCIPPar && par < dblLastAlgParSCPS ) {
  const std::string & scip_par =
   SMSpp_to_SCIP_dbl_pars[ par - dblFirstSCIPPar ];
  SCIP_CALL_ABORT( SCIPsetRealParam( scip, scip_par.c_str(), value ) );
  return;
 }

 MILPSolver::set_par( par, value );
}

/*--------------------------------------------------------------------------*/

void SCIPMILPSolver::set_par( idx_type par, const std::string & value ) {

 // SCIP parameters
 if( par >= strFirstSCIPPar && par < strLastAlgParSCPS ) {
  const std::string & scip_par =
   SMSpp_to_SCIP_str_pars[ par - strFirstSCIPPar ];

  // Char and string SCIP parameters are handled as SMS++ int parameters
  SCIP_PARAM * param = SCIPgetParam( scip, scip_par.c_str() );
  SCIP_PARAMTYPE type = SCIPparamGetType( param );

  if( type == SCIP_PARAMTYPE_CHAR ) {
   SCIP_CALL_ABORT( SCIPsetCharParam( scip, scip_par.c_str(), value[ 0 ] ) );
  } else if( type == SCIP_PARAMTYPE_STRING ) {
   SCIP_CALL_ABORT( SCIPsetStringParam( scip, scip_par.c_str(),
                                        value.c_str() ) );
  }
  return;
 }

 MILPSolver::set_par( par, value );
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
 SCIP_Bool bool_val;
 SCIP_Longint long_val;

 // Solver parameters explicitly mapped in SCIP
 switch( par ) {
  case intMaxIter:
   SCIP_CALL_ABORT( SCIPgetLongintParam( scip, "limits/nodes", &long_val ) );
   return ( int ) long_val;
  case intMaxSol:
   SCIP_CALL_ABORT( SCIPgetIntParam( scip, "limits/solutions", &value ) );
   return value;
  case intLogVerb:
   SCIP_CALL_ABORT( SCIPgetIntParam( scip, "display/verblevel", &value ) );
   return value;
  default:;
 }

 // SCIP parameters
 if( par >= intFirstSCIPPar && par < intLastAlgParSCPS ) {
  const std::string & scip_par =
   SMSpp_to_SCIP_int_pars[ par - intFirstSCIPPar ];

  // Bool, int and long SCIP parameters are handled as SMS++ int parameters
  SCIP_PARAM * param = SCIPgetParam( scip, scip_par.c_str() );
  SCIP_PARAMTYPE type = SCIPparamGetType( param );

  switch( type ) {
   case SCIP_PARAMTYPE_BOOL:
    SCIP_CALL_ABORT( SCIPgetBoolParam( scip, scip_par.c_str(), &bool_val ) );
    return ( int ) bool_val;
   case SCIP_PARAMTYPE_INT:
    SCIP_CALL_ABORT( SCIPgetIntParam( scip, scip_par.c_str(), &value ) );
    return value;
   case SCIP_PARAMTYPE_LONGINT:
    SCIP_CALL_ABORT( SCIPgetLongintParam( scip, scip_par.c_str(), &long_val ) );
    return ( int ) long_val;
   default:;
  }
 }

 return MILPSolver::get_int_par( par );
}

/*--------------------------------------------------------------------------*/

double SCIPMILPSolver::get_dbl_par( idx_type par ) const {
 double value;

 // Solver parameters explicitly mapped in SCIP
 switch( par ) {
  case dblMaxTime:
   SCIP_CALL_ABORT( SCIPgetRealParam( scip, "limits/time", &value ) );
   return value;
   // case dblRelAcc:   // TODO
   //  return 1e-6;
   // case dblAbsAcc:   // TODO
   //  return Inf< OFValue >();
   // case dblUpCutOff:
   //  if( SCIPgetObjsense( scip ) == SCIP_OBJSENSE_MINIMIZE ) {
   //   return SCIPgetObjlimit( scip );
   //  }
   //  return Inf< OFValue >();
   // case dblLwCutOff:
   //  if( SCIPgetObjsense( scip ) == SCIP_OBJSENSE_MAXIMIZE ) {
   //   return SCIPgetObjlimit( scip );
   //  }
   //  return -Inf< OFValue >();
  case dblRAccSol:
   SCIP_CALL_ABORT( SCIPgetRealParam( scip, "limits/gap", &value ) );
   return value;
  case dblAAccSol:
   SCIP_CALL_ABORT( SCIPgetRealParam( scip, "limits/absgap", &value ) );
   return value;
  case dblFAccSol:
   SCIP_CALL_ABORT( SCIPgetRealParam( scip, "numerics/feastol", &value ) );
   return value;
  default:;
 }

 // SCIP parameters
 if( par >= dblFirstSCIPPar && par < dblLastAlgParSCPS ) {
  const std::string & scip_par =
   SMSpp_to_SCIP_dbl_pars[ par - dblFirstSCIPPar ];
  SCIP_CALL_ABORT( SCIPgetRealParam( scip, scip_par.c_str(), &value ) );
  return value;
 }

 return MILPSolver::get_dbl_par( par );
}

/*--------------------------------------------------------------------------*/

const std::string &
SCIPMILPSolver::get_str_par( const idx_type par ) const {
 // TODO: check if all these are necessary
 char char_val;
 char * str_val;
 std::string return_str;

 // SCIP parameters
 if( par >= strFirstSCIPPar && par < strLastAlgParSCPS ) {
  const std::string & scip_par =
   SMSpp_to_SCIP_str_pars[ par - strFirstSCIPPar ];

  // Char and string SCIP parameters are handled as SMS++ string parameters
  SCIP_PARAM * param = SCIPgetParam( scip, scip_par.c_str() );
  SCIP_PARAMTYPE type = SCIPparamGetType( param );

  switch( type ) {
   case SCIP_PARAMTYPE_CHAR:
    SCIP_CALL_ABORT( SCIPgetCharParam( scip, scip_par.c_str(), &char_val ) );
    return std::move( std::to_string( char_val ) );
   case SCIP_PARAMTYPE_STRING:
    SCIP_CALL_ABORT( SCIPgetStringParam( scip, scip_par.c_str(), &str_val ) );
    return_str = str_val;
    delete[] str_val;
    return std::move( return_str );
   default:;
  }
 }

 return MILPSolver::get_str_par( par );
}

/*--------------------------------------------------------------------------*/

int SCIPMILPSolver::get_dflt_int_par( const idx_type par ) const {
 int value;
 SCIP_Longint long_val;
 SCIP_PARAM * param;

 // Solver parameters explicitly mapped in SCIP
 switch( par ) {
  case intMaxIter:
   param = SCIPgetParam( scip, "limits/nodes" );
   return ( int ) SCIPparamGetLongintDefault( param );
  case intMaxSol:
   param = SCIPgetParam( scip, "limits/solutions" );
   return SCIPparamGetIntDefault( param );
  case intLogVerb:
   param = SCIPgetParam( scip, "display/verblevel" );
   return SCIPparamGetIntDefault( param );
  default:;
 }

 // SCIP parameters
 if( par >= intFirstSCIPPar && par < intLastAlgParSCPS ) {
  const std::string & scip_par =
   SMSpp_to_SCIP_int_pars[ par - intFirstSCIPPar ];

  // Bool, int and long SCIP parameters are handled as SMS++ int parameters
  param = SCIPgetParam( scip, scip_par.c_str() );
  SCIP_PARAMTYPE type = SCIPparamGetType( param );

  switch( type ) {
   case SCIP_PARAMTYPE_BOOL:
    return ( int ) SCIPparamGetBoolDefault( param );
   case SCIP_PARAMTYPE_INT:
    return SCIPparamGetIntDefault( param );
   case SCIP_PARAMTYPE_LONGINT:
    return ( int ) SCIPparamGetLongintDefault( param );
   default:;
  }
 }

 return MILPSolver::get_dflt_int_par( par );
}

/*--------------------------------------------------------------------------*/

double SCIPMILPSolver::get_dflt_dbl_par( const idx_type par ) const {
 double value;
 SCIP_PARAM * param;

 switch( par ) {
  case dblMaxTime:
   param = SCIPgetParam( scip, "limits/time" );
   return SCIPparamGetRealDefault( param );
   // case dblRelAcc:   // TODO
   //  return 1e-6;
   // case dblAbsAcc:   // TODO
   //  return Inf< OFValue >();
   // case dblUpCutOff:
   //  if( SCIPgetObjsense( scip ) == SCIP_OBJSENSE_MINIMIZE ) {
   //   return -Inf< OFValue >();
   //  }
   //  return Inf< OFValue >();
   // case dblLwCutOff:
   //  if( SCIPgetObjsense( scip ) == SCIP_OBJSENSE_MAXIMIZE ) {
   //   return Inf< OFValue >();
   //  }
   //  return -Inf< OFValue >();
  case dblRAccSol:
   param = SCIPgetParam( scip, "limits/gap" );
   return SCIPparamGetRealDefault( param );
  case dblAAccSol:
   param = SCIPgetParam( scip, "limits/absgap" );
   return SCIPparamGetRealDefault( param );
  case dblFAccSol:
   param = SCIPgetParam( scip, "numerics/feastol" );
   return SCIPparamGetRealDefault( param );
  default:;
 }

 // SCIP parameters
 if( par >= dblFirstSCIPPar && par < dblLastAlgParSCPS ) {
  const std::string & scip_par =
   SMSpp_to_SCIP_dbl_pars[ par - dblFirstSCIPPar ];

  param = SCIPgetParam( scip, scip_par.c_str() );
  SCIP_PARAMTYPE type = SCIPparamGetType( param );
  return SCIPparamGetRealDefault( param );
 }

 return MILPSolver::get_dflt_dbl_par( par );
}

/*--------------------------------------------------------------------------*/

const std::string &
SCIPMILPSolver::get_dflt_str_par( const idx_type par ) const {
 // TODO: check if all these are necessary
 std::string char_val( 1, '\0' );
 char * str_val;
 std::string return_str;

 // SCIP parameters
 if( par >= strFirstSCIPPar && par < strLastAlgParSCPS ) {
  const std::string & scip_par =
   SMSpp_to_SCIP_str_pars[ par - strFirstSCIPPar ];

  // Char and string SCIP parameters are handled as SMS++ string parameters
  SCIP_PARAM * param = SCIPgetParam( scip, scip_par.c_str() );
  SCIP_PARAMTYPE type = SCIPparamGetType( param );

  switch( type ) {
   case SCIP_PARAMTYPE_CHAR:
    char_val[ 0 ] = SCIPparamGetCharDefault( param );
    return std::move( char_val );
   case SCIP_PARAMTYPE_STRING:
    str_val = SCIPparamGetStringDefault( param );
    return_str = str_val;
    delete[] str_val;
    return std::move( return_str );
   default:;
  }
 }

 return MILPSolver::get_dflt_str_par( par );
}

/*--------------------------------------------------------------------------*/

ThinComputeInterface::idx_type
SCIPMILPSolver::int_par_str2idx( const std::string & name ) const {

 // SCIP parameters
 auto it = lower_bound( SCIP_to_SMSpp_int_pars.begin(),
                        SCIP_to_SMSpp_int_pars.end(),
                        std::make_pair( name, 0 ) );

 if( it != SCIP_to_SMSpp_int_pars.end() && it->first == name ) {
  return it->second;
 }

 return MILPSolver::int_par_str2idx( name );
}

/*--------------------------------------------------------------------------*/

const std::string &
SCIPMILPSolver::int_par_idx2str( const idx_type idx ) const {

 // SCIP parameters
 if( idx >= intFirstSCIPPar && idx < intLastAlgParSCPS ) {
  return SMSpp_to_SCIP_int_pars[ idx - intFirstSCIPPar ];
 }

 return MILPSolver::int_par_idx2str( idx );
}

/*--------------------------------------------------------------------------*/

ThinComputeInterface::idx_type
SCIPMILPSolver::dbl_par_str2idx( const std::string & name ) const {

 // SCIP parameters
 auto it = lower_bound( SCIP_to_SMSpp_dbl_pars.begin(),
                        SCIP_to_SMSpp_dbl_pars.end(),
                        std::make_pair( name, 0 ) );

 if( it != SCIP_to_SMSpp_dbl_pars.end() && it->first == name ) {
  return it->second;
 }

 return MILPSolver::dbl_par_str2idx( name );
}

/*--------------------------------------------------------------------------*/

const std::string &
SCIPMILPSolver::dbl_par_idx2str( const idx_type idx ) const {

 // SCIP parameters
 if( idx >= dblFirstSCIPPar && idx < dblLastAlgParSCPS ) {
  return SMSpp_to_SCIP_dbl_pars[ idx - dblFirstSCIPPar ];
 }

 return MILPSolver::dbl_par_idx2str( idx );
}

/*--------------------------------------------------------------------------*/


ThinComputeInterface::idx_type
SCIPMILPSolver::str_par_str2idx( const std::string & name ) const {

 // SCIP parameters
 auto it = lower_bound( SCIP_to_SMSpp_str_pars.begin(),
                        SCIP_to_SMSpp_str_pars.end(),
                        std::make_pair( name, 0 ) );

 if( it != SCIP_to_SMSpp_str_pars.end() && it->first == name ) {
  return it->second;
 }

 return MILPSolver::str_par_str2idx( name );
}

/*--------------------------------------------------------------------------*/

const std::string &
SCIPMILPSolver::str_par_idx2str( const idx_type idx ) const {

 // SCIP parameters
 if( idx >= strFirstSCIPPar && idx < strLastAlgParSCPS ) {
  return SMSpp_to_SCIP_str_pars[ idx - strFirstSCIPPar ];
 }

 return MILPSolver::str_par_idx2str( idx );
}

/*--------------------------------------------------------------------------*/
/*--------------------- End File SCIPMILPSolver.cpp -------------------------*/
/*--------------------------------------------------------------------------*/
