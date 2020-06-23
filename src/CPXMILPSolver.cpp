/*--------------------------------------------------------------------------*/
/*------------------------- File CPXMILPSolver.cpp -------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Implementation of the CPXMILPSolver class.
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
 * \author Kostas Tavlaridis-Gyparakis \n
 *         Operations Research Group \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * Copyright &copy by Antonio Frangioni, Kostas Tavlaridis-Gyparakis, Niccolò Iardella
 */
/*--------------------------------------------------------------------------*/
/*---------------------------- IMPLEMENTATION ------------------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include <functional>
#include <queue>

#include <Block.h>
#include <MILPSolver.h>
#include <OneVarConstraint.h>
#include <LinearFunction.h>
#include <DQuadFunction.h>

#include "CPXMILPSolver.h"

/*--------------------------------------------------------------------------*/
/*------------------------- NAMESPACE AND USING ----------------------------*/
/*--------------------------------------------------------------------------*/

using namespace SMSpp_di_unipi_it;

SMSpp_insert_in_factory_cpp_0( CPXMILPSolver );

/*--------------------------------------------------------------------------*/
/*--------------------- CONSTRUCTOR AND DESTRUCTOR -------------------------*/
/*--------------------------------------------------------------------------*/

CPXMILPSolver::CPXMILPSolver() : MILPSolver() {
 int status = 0;
 env = CPXopenCPLEX( &status );
 if( env == nullptr ) {
  throw std::runtime_error( "CPXopenCPLEX returned with status " +
                            std::to_string( status ) );
 }
 lp = nullptr;
}

CPXMILPSolver::~CPXMILPSolver() {
 if( lp ) {
  CPXfreeprob( env, &lp );
 }
 CPXcloseCPLEX( &env );
}

/*--------------------------------------------------------------------------*/
/*--------------------- DERIVED METHODS OF BASE CLASS ----------------------*/
/*--------------------------------------------------------------------------*/

void CPXMILPSolver::set_Block( Block * block ) {
 if( block == f_Block ) {
  return;
 }
 MILPSolver::set_Block( block );
}

/*--------------------------------------------------------------------------*/

void CPXMILPSolver::clear_problem() {
 MILPSolver::clear_problem();

 if( lp ) {
  CPXfreeprob( env, &lp );
 }
}

/*--------------------------------------------------------------------------*/

void CPXMILPSolver::load_problem() {
 MILPSolver::load_problem();

 int status = 0;
 lp = CPXcreateprob( env, &status, prob_name.c_str() );

 for( int i = 0; i < numcols; ++i ) {
  if( lb[ i ] == -Inf< double >() ) {
   lb[ i ] = -CPX_INFBOUND;
  }
  if( ub[ i ] == Inf< double >() ) {
   ub[ i ] = CPX_INFBOUND;
  }
 }

 if( use_custom_names ) {
  CPXcopylpwnames( env, lp,
                   numcols,
                   numrows,
                   objsense,
                   objective.data(),
                   rhs.data(),
                   sense.data(),
                   matbeg.data(),
                   matcnt.data(),
                   matind.data(),
                   matval.data(),
                   lb.data(),
                   ub.data(),
                   rngval.data(),
                   colname.data(),
                   rowname.data() );
 } else {
  CPXcopylp( env, lp,
             numcols,
             numrows,
             objsense,
             objective.data(),
             rhs.data(),
             sense.data(),
             matbeg.data(),
             matcnt.data(),
             matind.data(),
             matval.data(),
             lb.data(),
             ub.data(),
             rngval.data() );
 }

 if( qp ) {
  // CPLEX evaluates the corresponding objective with a factor
  // of 0.5 in front of the quadratic objective term.
  std::vector< double > double_q_obj = q_objective;
  for( auto & i: double_q_obj ) {
   i = i * 2;
  }
  CPXcopyqpsep( env, lp, double_q_obj.data() );
 }

 if( mip ) {
  // Adding ctype information automatically changes the problem type
  // from continuous to mixed integer
  CPXcopyctype( env, lp, xctype.data() );
 }
}

/*--------------------------------------------------------------------------*/

int CPXMILPSolver::compute( bool changedvars ) {

 int status = 0;
 process_modifications();
 // std::stringstream lpname;
 // lpname << std::setfill( '0' ) << std::setw( 2 ) << print_debug++
 //        << "_computeaftermods.lp";
 // output_file = lpname.str();

 if( !output_file.empty() )
  CPXwriteprob( env, lp, output_file.c_str(), "LP" );

#if MILPSLVR_DEBUG
 int probtype = CPXgetprobtype( env, lp );
 LOG( "[DEBUG] ========= Problem type:" );
 switch( probtype ) {
  case CPXPROB_LP :
   LOG( "CPXPROB_LP" );
   break;
  case CPXPROB_MILP :
   LOG( "CPXPROB_MILP" );
   break;
  case CPXPROB_FIXEDMILP :
   LOG( "CPXPROB_FIXEDMILP" );
   break;
  case CPXPROB_QP :
   LOG( "CPXPROB_QP" );
   break;
  case CPXPROB_MIQP :
   LOG( "CPXPROB_MIQP" );
   break;
  case CPXPROB_FIXEDMIQP :
   LOG( "CPXPROB_FIXEDMIQP" );
   break;
  case CPXPROB_QCP :
   LOG( "CPXPROB_QCP" );
   break;
  case CPXPROB_MIQCP :
   LOG( "CPXPROB_MIQCP" );
   break;
  default:
   throw std::runtime_error( "Undefined CPLEX problem type" );
 }
 LOG( std::endl );
#endif

 if( mip ) {
  status = CPXmipopt( env, lp );
  if( status ) {
   throw std::runtime_error( "CPXmipopt() encountered an error" );
  }

  status = CPXgetstat( env, lp );
  switch( status ) {
   case CPXMIP_OPTIMAL:
   case CPXMIP_OPTIMAL_TOL:
   case CPXMIP_SOL_LIM:
    sol_status = kOK;
    break;
   case CPXMIP_INFEASIBLE:
    sol_status = kInfeasible;
    break;
   case CPXMIP_NODE_LIM_FEAS:
   case CPXMIP_NODE_LIM_INFEAS:
    sol_status = kStopIter;
    break;
   case CPXMIP_TIME_LIM_FEAS:
   case CPXMIP_TIME_LIM_INFEAS:
    sol_status = kStopTime;
    break;
   case CPXMIP_FAIL_FEAS:
   case CPXMIP_FAIL_INFEAS:
    sol_status = kError;
    break;
   case CPXMIP_UNBOUNDED:
    sol_status = kUnbounded;
    break;
   default:
    sol_status = status;
    break;
  }
  nodes = CPXgetnodecnt( env, lp );

 } else {
  status = CPXlpopt( env, lp );
  if( status ) {
   throw std::runtime_error( "CPXlpopt() encountered an error" );
  }

  status = CPXgetstat( env, lp );
  switch( status ) {
   case CPX_STAT_OPTIMAL :
    sol_status = kOK;
    break;
   case CPX_STAT_INFEASIBLE :
    sol_status = kInfeasible;
    break;
   case CPX_STAT_ABORT_IT_LIM :
    sol_status = kStopIter;
    break;
   case CPX_STAT_ABORT_TIME_LIM :
    sol_status = kStopTime;
    break;
   case CPX_STAT_UNBOUNDED :
   case CPX_STAT_INForUNBD :
    sol_status = kUnbounded;
    break;
   default:
    sol_status = status;
    break;
  }
 }

 return sol_status;
}

/*--------------------------------------------------------------------------*/

Solver::OFValue CPXMILPSolver::get_lb() {

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
     CPXgetbestobjval( env, lp, &lower_bound );
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
     CPXgetbestobjval( env, lp, &lower_bound );
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

Solver::OFValue CPXMILPSolver::get_ub() {

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
     CPXgetobjval( env, lp, &upper_bound );
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
     CPXgetobjval( env, lp, &upper_bound );
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

bool CPXMILPSolver::has_var_solution() {
 switch( sol_status ) {
  case ( kOK ):
  case ( kInfeasible ):
   return ( true );
  default:
   return ( false );
 }
}

/*--------------------------------------------------------------------------*/

bool CPXMILPSolver::is_var_feasible() {
 switch( sol_status ) {
  case ( kInfeasible ):
   return ( false );
  default:
   return ( true );
 }
}

/*--------------------------------------------------------------------------*/

Solver::OFValue CPXMILPSolver::get_var_value() {
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

void CPXMILPSolver::get_var_solution( Configuration * solc ) {

#if MILPSLVR_DEBUG
 std::cout << "[DEBUG] ========= MILPSolver::get_var_solution()" << std::endl;
#endif

 auto * x = new double[numcols];
 int status = CPXgetx( env, lp, x, 0, numcols - 1 );
 if( status ) {
  delete[] x;
  throw std::runtime_error( "Unable to get the solution values with CPXgetx()" );
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
   auto f1 = std::bind( &CPXMILPSolver::set_var_value,
                        this,
                        std::placeholders::_1,
                        std::ref( x ),
                        std::ref( col ) );
   un_any_const_static( i, f1, un_any_type< ColVariable >() );
  }

  for( const auto & i : q_Block->get_dynamic_variables() ) {
   auto f1 = std::bind( &CPXMILPSolver::set_var_value,
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

 // After the Objective is computed (evaluated), the solution can be retrieved
 // directly from there.
 // auto p_obj = dynamic_cast< FRealObjective * >( f_Block->get_objective() );
 // if( p_obj ) {
 //  p_obj->compute();
 // }

 delete[]x;
}

/*--------------------------------------------------------------------------*/

bool CPXMILPSolver::has_dual_solution() {
 // FIXME: Use CPXsolninfo()
 auto * pi = new double[numrows];
 int status = CPXgetpi( env, lp, pi, 0, numrows - 1 );
 delete[]pi;
 return status == 0;
}

/*--------------------------------------------------------------------------*/

void CPXMILPSolver::get_dual_solution( Configuration * solc ) {

#if MILPSLVR_DEBUG
 std::cout << "[DEBUG] ========= MILPSolver::get_dual_solution()" << std::endl;
#endif

 int method = CPXgetmethod( env, lp );
 int status = CPXgetstat( env, lp );
 double * res = nullptr;
 // std::vector<double> res;

 if( method == CPX_ALG_PRIMAL ) {
  // Primal simplex optimizer is used
  if( status == CPX_STAT_UNBOUNDED ) {
   // The model is primal unbounded/dual infeasible

   res = new double[numcols];
   status = CPXgetray( env, lp, res );
  }
  if( status == CPX_STAT_INFEASIBLE ) {
   // The model is primal infeasible/dual unbounded
   res = new double[numrows];
   status = CPXgetpi( env, lp, res, 0, numrows - 1 );
  }
 }

 if( method == CPX_ALG_DUAL ) {
  // Dual simplex optimizer is used
  if( status == CPX_STAT_INFEASIBLE ) {
   // The model is dual infeasible/primal unbounded
   res = new double[numcols];
   status = CPXgetray( env, lp, res );
  }
  if( status == CPX_STAT_UNBOUNDED ) {
   // The model is dual unbounded/primal infeasible
   res = new double[numrows];
   status = CPXdualfarkas( env, lp, res, nullptr );
  }
 }

 // // TODO ?
 // if( status || res.size() != numrows ) {
 //  return;
 // }

 int row = 0;

 std::queue< Block * > Q;
 Q.push( f_Block );

 while( !Q.empty() ) {
  Block * q_Block = Q.front();
  Q.pop();

  for( auto * i : q_Block->get_nested_Blocks() ) {
   Q.push( i );
  }

  for( const auto & i : q_Block->get_static_constraints() ) {
   auto f1 = std::bind( &CPXMILPSolver::set_dual_value,
                        this,
                        std::placeholders::_1,
                        std::ref( res ),
                        std::ref( row ) );
   un_any_const_static( i, f1, un_any_type< FRowConstraint >() );
  }

  for( const auto & i : q_Block->get_dynamic_constraints() ) {
   auto f1 = std::bind( &CPXMILPSolver::set_dual_value,
                        this,
                        std::placeholders::_1,
                        std::ref( res ),
                        std::ref( row ) );
   un_any_const_dynamic( i, f1, un_any_type< FRowConstraint >() );
  }
 }
 delete[] res;
}

/*--------------------------------------------------------------------------*/

void CPXMILPSolver::write_lp( const std::string & filename ) {
 CPXwriteprob( env, lp, filename.c_str(), "LP" );
}

/*--------------------------------------------------------------------------*/
/*-------------------- PROTECTED FIELDS OF THE CLASS -----------------------*/
/*--------------------------------------------------------------------------*/

void CPXMILPSolver::var_modification( VariableMod * mod ) {

 /*
  * VariableMod class does not include any modification types, so we "refresh"
  * all the CPLEX information for the variable. In particular we:
  *
  *  - reset the ctype (binary, integer, continuous)
  *  - if the variable is fixed, we fix it in CPLEX by doing LHS = RHS
  *  - if the variable is not fixed, we reset LHS and RHS values
  *
  * The CPLEX problem type must also be updated if we remove the last integer
  * variable of a MIP or add a integer variable to a LP/QP.
  *
  */

 auto * var = dynamic_cast<ColVariable *>(mod->variable());

 int idx = index_of_variable( var );
 std::vector< int > indices( 2, idx );
 std::array< char, 1 > ctype{};
 std::vector< char > lu;
 std::vector< double > bd;

 int status = 0;
 int probtype = CPXgetprobtype( env, lp );

 // Read old type

 status = CPXgetctype( env, lp, ctype.data(), idx, idx );
 if( status == 0 ) {
  // Problem is MIP
  if( ctype[ 0 ] == 'B' || ctype[ 0 ] == 'I' ) {
   --mip;
  }
 } else if( status == CPXERR_NOT_MIP ) {
  // Problem is LP/QP, nothing to do
 } else {
  throw std::runtime_error( "CPXgetctype() returned " +
                            std::to_string( status ) );
 }

 // Read new type

 if( var->is_integer() ) {
  ++mip;
  if( var->is_unitary() && var->is_positive() ) {
   ctype[ 0 ] = 'B'; // Binary
  } else {
   ctype[ 0 ] = 'I'; // Integer
  }
 } else {
  ctype[ 0 ] = 'C';  // Continuous
 }

 // Update problem type

 if( mip == 0 ) {
  // The last integer variable was removed, or the problem stays continuous
  // This call removes all ctype values
  switch( probtype ) {
   case CPXPROB_LP :
   case CPXPROB_QP :
    break;
   case CPXPROB_MILP :
   case CPXPROB_FIXEDMILP :
    CPXchgprobtype( env, lp, CPXPROB_LP );
    break;
   case CPXPROB_MIQP :
   case CPXPROB_FIXEDMIQP :
    CPXchgprobtype( env, lp, CPXPROB_QP );
    break;
   default:
    throw std::runtime_error( "Wrong CPLEX problem type" );
  }

 } else if( mip == 1 ) {
  // The first integer variable was added
  // All ctype values must be [re]added to the problem
  switch( probtype ) {
   case CPXPROB_LP :
    CPXchgprobtype( env, lp, CPXPROB_MILP );
    break;
   case CPXPROB_QP :
    CPXchgprobtype( env, lp, CPXPROB_MIQP );
    break;
   default:
    throw std::runtime_error( "Wrong CPLEX problem type" );
  }

  std::vector< char > new_ctype( CPXgetnumcols( env, lp ), 'C' );
  new_ctype[ idx ] = ctype[ 0 ];
  CPXcopyctype( env, lp, new_ctype.data() );

 } else {
  // The problem stays a MIP, update only the one variable
  CPXchgctype( env, lp, 1, indices.data(), ctype.data() );
 }

 if( var->is_fixed() ) {
  lu.resize( 1 );
  bd.resize( 1 );
  lu[ 0 ] = 'B';
  bd[ 0 ] = var->get_value();
  CPXchgbds( env, lp, 1, indices.data(), lu.data(), bd.data() );

 } else {
  lu.resize( 2 );
  bd.resize( 2 );
  lu[ 0 ] = 'L';
  lu[ 1 ] = 'U';
  bd[ 0 ] = -CPX_INFBOUND;
  bd[ 1 ] = CPX_INFBOUND;
  for( auto * bnd : active_bounds[ indices[ 0 ] ] ) {
   bd[ 0 ] = bd[ 0 ] > bnd->get_lhs() ? bd[ 0 ] : bnd->get_lhs();
   bd[ 1 ] = bd[ 1 ] < bnd->get_rhs() ? bd[ 1 ] : bnd->get_rhs();
  }
  CPXchgbds( env, lp, 2, indices.data(), lu.data(), bd.data() );
 }
}

/*--------------------------------------------------------------------------*/

void CPXMILPSolver::of_modification( ObjectiveMod * mod ) {

 /*
  * ObjectiveMod class does not include any modification types except
  * for eSetMin and eSetMax.
  * To change OF coefficents, a FunctionMod must be used.
  */

 switch( mod->type() ) {

  case ObjectiveMod::eSetMin:
   CPXchgobjsen( env, lp, 1 );
   break;

  case ObjectiveMod::eSetMax:
   CPXchgobjsen( env, lp, -1 );
   break;

  default:
   throw std::invalid_argument( "Invalid type of ObjectiveMod" );
 }
}

/*--------------------------------------------------------------------------*/

void CPXMILPSolver::const_modification( ConstraintMod * mod ) {

 /*
  * To change the coefficents, a FunctionMod must be used.
  */

 auto * p_const = dynamic_cast<FRowConstraint *>(mod->constraint());

 const int cnt = 1;
 std::array< int, cnt > indices{};
 std::array< double, cnt > values{};
 std::array< char, cnt > sense{};
 std::array< double, cnt > rngval{};

 RowConstraint::RHSValue const_lhs = NAN;
 RowConstraint::RHSValue const_rhs = NAN;

 switch( mod->type() ) {

  case ConstraintMod::eRelaxConst:
   // In order to relax the constraint all we do is transform it
   // into an inequality with RHS equal to infinity

   sense[ 0 ] = ( 'G' );
   values[ 0 ] = -Inf< double >();
   indices[ 0 ] = index_of_constraint( p_const );
   CPXchgrhs( env, lp, cnt, indices.data(), values.data() );
   CPXchgsense( env, lp, cnt, indices.data(), sense.data() );
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

   const_lhs = p_const->get_lhs();
   const_rhs = p_const->get_rhs();

   if( const_lhs == const_rhs ) {
    sense[ 0 ] = 'E';
    values[ 0 ] = const_rhs;
   } else if( const_lhs == -Inf< double >() ) {
    sense[ 0 ] = 'L';
    values[ 0 ] = const_rhs;
   } else if( const_rhs == Inf< double >() ) {
    sense[ 0 ] = 'G';
    values[ 0 ] = const_lhs;
   } else {
    sense[ 0 ] = 'R';
    values[ 0 ] = const_rhs;
    rngval[ 0 ] = const_rhs - const_lhs;
   }

   indices[ 0 ] = index_of_constraint( p_const );
   CPXchgrhs( env, lp, cnt, indices.data(), values.data() );
   CPXchgsense( env, lp, cnt, indices.data(), sense.data() );
   if( sense[ 0 ] == 'R' ) {
    CPXchgrngval( env, lp, cnt, indices.data(), rngval.data() );
   }
   break;

  default:
   throw std::invalid_argument( "Invalid type of ObjectiveMod" );
 }
}

/*--------------------------------------------------------------------------*/

void CPXMILPSolver::bound_modification( OneVarConstraintMod * mod ) {

 /*
  * The same ColVariable can have more active OneVarConstraints,
  * so each time we modify one of them we have to check if LHS and RHS
  * of the Variable change.
  */

 auto * p_const = dynamic_cast<OneVarConstraint *>(mod->constraint());
 auto * p_var = dynamic_cast<ColVariable *>(p_const->get_active_var( 0 ));

 std::vector< int > indices( 2, index_of_variable( p_var ) );
 std::vector< char > lu;
 std::vector< double > bd;

 switch( mod->type() ) {

  case RowConstraintMod::eChgLHS:
   lu.resize( 1 );
   bd.resize( 1 );
   lu[ 0 ] = 'L';
   bd[ 0 ] = -CPX_INFBOUND;

   for( auto * bnd : active_bounds[ indices[ 0 ] ] ) {
    bd[ 0 ] = bd[ 0 ] > bnd->get_lhs() ? bd[ 0 ] : bnd->get_lhs();
   }

   CPXchgbds( env, lp, 1, indices.data(), lu.data(), bd.data() );
   break;

  case RowConstraintMod::eChgRHS:
   lu.resize( 1 );
   bd.resize( 1 );
   lu[ 0 ] = 'U';
   bd[ 0 ] = CPX_INFBOUND;

   for( auto * bnd : active_bounds[ indices[ 0 ] ] ) {
    bd[ 0 ] = bd[ 0 ] < bnd->get_rhs() ? bd[ 0 ] : bnd->get_rhs();
   }

   CPXchgbds( env, lp, 1, indices.data(), lu.data(), bd.data() );
   break;

  case RowConstraintMod::eChgBTS:
   lu.resize( 2 );
   bd.resize( 2 );
   lu[ 0 ] = 'L';
   lu[ 1 ] = 'U';
   bd[ 0 ] = -CPX_INFBOUND;
   bd[ 1 ] = CPX_INFBOUND;

   for( auto * bnd : active_bounds[ indices[ 0 ] ] ) {
    bd[ 0 ] = bd[ 0 ] > bnd->get_lhs() ? bd[ 0 ] : bnd->get_lhs();
    bd[ 1 ] = bd[ 1 ] < bnd->get_rhs() ? bd[ 1 ] : bnd->get_rhs();
   }

   CPXchgbds( env, lp, 2, indices.data(), lu.data(), bd.data() );
   break;

  default:
   throw std::invalid_argument( "Invalid type of OneVarConstraintMod" );
 }
}

/*--------------------------------------------------------------------------*/

void CPXMILPSolver::function_modification( FunctionMod * mod ) {
 // TODO: Change only involved variables, see function_vars_modification() below
 /*
  * This function is used when changing coefficents for OFs or constraints.
  */
 auto * mod_f = mod->function();
 auto * p_obj = dynamic_cast< FRealObjective * >( f_Block->get_objective() );
 auto * of = p_obj->get_function();

 if( of == mod_f ) {
  // Changing objective function

  const auto * lf = dynamic_cast<const LinearFunction *> (mod_f);
  const auto * qf = dynamic_cast<const DQuadFunction *> (mod_f);

  int num_vars = 0;
  std::vector< int > indices;
  std::vector< double > values;
  int probtype = CPXgetprobtype( env, lp );

  if( lf != nullptr ) {
   // Linear objective function

   num_vars = static_cast<int>(lf->get_v_var().size());
   indices.resize( num_vars );
   values.resize( num_vars );

   // Update problem type
   switch( probtype ) {
    case CPXPROB_LP :
    case CPXPROB_MILP :
    case CPXPROB_FIXEDMILP :
     break;
    case CPXPROB_QP :
     CPXchgprobtype( env, lp, CPXPROB_LP );
     break;
    case CPXPROB_MIQP :
     CPXchgprobtype( env, lp, CPXPROB_MILP );
     break;
    case CPXPROB_FIXEDMIQP :
     CPXchgprobtype( env, lp, CPXPROB_FIXEDMILP );
     break;
    default:
     throw std::runtime_error( "Wrong CPLEX problem type" );
   }

   // Update objective coefficients
   int i = 0;
   for( auto el : lf->get_v_var() ) {
    indices[ i ] = index_of_variable( el.first );
    values[ i ] = el.second;
    ++i;
   }
   CPXchgobj( env, lp, num_vars, indices.data(), values.data() );

  } else if( qf != nullptr ) {
   // Quadratic objective function

   num_vars = static_cast<int>(qf->get_v_var().size());
   indices.resize( num_vars );
   values.resize( num_vars );

   // Update problem type
   switch( probtype ) {
    case CPXPROB_LP :
     CPXchgprobtype( env, lp, CPXPROB_QP );
     break;
    case CPXPROB_MILP :
     CPXchgprobtype( env, lp, CPXPROB_MIQP );
     break;
    case CPXPROB_FIXEDMILP :
     CPXchgprobtype( env, lp, CPXPROB_FIXEDMIQP );
     break;
    case CPXPROB_QP :
    case CPXPROB_MIQP :
    case CPXPROB_FIXEDMIQP :
     break;
    default:
     throw std::runtime_error( "Wrong CPLEX problem type" );
   }

   int i = 0;
   for( auto el : qf->get_v_var() ) {
    // Linear coefficients can be changed all at once with CPXchgobj
    indices[ i ] = index_of_variable( std::get< 0 >( el ) );
    values[ i ] = std::get< 1 >( el );

    // Quadratic coefficients can be changed one at a time
    CPXchgqpcoef( env, lp, indices[ i ], indices[ i ], std::get< 2 >( el ) );
    ++i;
   }

   CPXchgobj( env, lp, num_vars, indices.data(), values.data() );

  } else {
   throw std::invalid_argument( "Unknown type of Objective Function" );
  }

 } else {
  // Changing a constraint function, so it can be only linear
  const auto * lf = dynamic_cast<const LinearFunction *> (mod_f);

  if( lf != nullptr ) {
   auto * p_const = dynamic_cast<FRowConstraint *>(lf->get_Observer());

   int indices[1];
   indices[ 0 ] = index_of_constraint( p_const );

   for( auto el : lf->get_v_var() ) {
    CPXchgcoef( env,
                lp,
                indices[ 0 ],
                index_of_variable( el.first ),
                el.second );
   }
  } else {
   throw std::invalid_argument( "Unknown type of Function" );
  }
 }
}

/*--------------------------------------------------------------------------*/

void CPXMILPSolver::function_vars_modification( FunctionModVars * mod ) {

 /*
  * This function is used when adding coefficents to OFs or constraints.
  */
 auto * mod_f = mod->function();
 auto * p_obj = dynamic_cast< FRealObjective * >( f_Block->get_objective() );
 auto * of = p_obj->get_function();

 auto * add = dynamic_cast<C05FunctionModVarsAddd *>( mod );
 if( add ) {

  int num_vars = static_cast<int>(add->vars().size());
  std::vector< int > indices( num_vars );
  std::vector< double > values( num_vars );

  if( of == mod_f ) {
   // Adding the coefficients to the objective function

   const auto * lf = dynamic_cast<const LinearFunction *> (mod_f);
   const auto * qf = dynamic_cast<const DQuadFunction *> (mod_f);

   if( lf != nullptr ) {
    // Linear objective function

    // Get indices and coefficients
    int i = 0;
    for( auto * it1 : add->vars() ) {
     for( auto it2: lf->get_v_var() ) {
      if( it1 == it2.first ) {
       indices[ i ] = index_of_variable( it2.first );
       values[ i ] = it2.second;
       break;
      }
     }
     ++i;
    }

    // Update the coefficients
    CPXchgobj( env, lp, num_vars, indices.data(), values.data() );

   } else if( qf != nullptr ) {
    // Quadratic objective function
    // TODO
   } else {
    throw std::invalid_argument( "Unknown type of Objective Function" );
   }

  } else {
   // Adding coefficients to a Constraint
   // TODO: Now it works because they are already added with a var
  }
  return;
 } // add

 auto * rmv = dynamic_cast<C05FunctionModVarsRngd *>( mod );
 if( rmv ) {

  int num_vars = static_cast<int>(rmv->vars().size());
  std::vector< int > indices( num_vars );
  std::vector< double > values( num_vars );

  if( of == mod_f ) {
   // Removing coefficients from the objective function

   const auto * lf = dynamic_cast<const LinearFunction *> (mod_f);
   const auto * qf = dynamic_cast<const DQuadFunction *> (mod_f);

   if( lf != nullptr ) {
    // TODO: The following doesn't find the vars in lf->get_v_var()?
    //       It works because var is removed with remove_dynamic_variable
    // // Linear objective function
    //
    // // Get indices and coefficients (all zeroes)
    // int i = 0;
    // for( auto * it1 : rmv->vars() ) {
    //  for( auto it2: lf->get_v_var() ) {
    //   if( it1 == it2.first ) {
    //    indices[ i ] = index_of_variable( it2.first );
    //    values[ i ] = 0;
    //    break;
    //   }
    //  }
    //  ++i;
    // }
    //
    // // Update the coefficients (all zeroes)
    // CPXchgobj( env, lp, num_vars, indices.data(), values.data() );

   } else if( qf != nullptr ) {
    // TODO
   } else {
    throw std::invalid_argument( "Unknown type of Objective Function" );
   }
  } else {
   // Removing coefficients from a constraint

   // const auto * lf = dynamic_cast<const LinearFunction *> (mod_f);
   // auto * p_const = dynamic_cast<FRowConstraint *>(lf->get_Observer());
   // std::vector< int > j( num_vars, index_of_constraint( p_const ) );
   //
   // // Get indices and coefficients (all zeroes)
   // int i = 0;
   // bool found = false;
   // for( auto * it1 : rmv->vars() ) {
   //  for( auto it2: lf->get_v_var() ) {
   //   if( it1 == it2.first ) {
   //    indices[ i ] = index_of_variable( it2.first );
   //    values[ i ] = 0;
   //    found = true;
   //    break;
   //   }
   //  }
   //  ++i;
   // }
   //
   // // Update the coefficients (all zeroes)
   // if (found)
   // CPXchgcoeflist( env, lp, num_vars, indices.data(), j.data(), values.data() );
   //
   // std::stringstream lpname;
   // lpname << std::setfill( '0' ) << std::setw( 3 ) << print_debug << "_function_vars_modification_RMV-CS.lp";
   // write_lp(lpname.str());
   // print_debug++;
  }
  return;
 } // rmv

 throw std::invalid_argument( "This type of FunctionModVars is not handled" );
}

/*--------------------------------------------------------------------------*/

void CPXMILPSolver::dynamic_modification( BlockModAD * mod ) {

 auto * addcon_mod = dynamic_cast<BlockModAdd< FRowConstraint > *>(mod);
 if( addcon_mod ) {
  for( auto * i : addcon_mod->added() ) {
   add_dynamic_constraint( i );
  }
  return;
 }

 auto * rmvcon_mod = dynamic_cast<BlockModRmv< FRowConstraint > *>(mod);
 if( rmvcon_mod ) {
  for( const auto & i : rmvcon_mod->removed() ) {
   remove_dynamic_constraint( &i );
  }
  return;
 }

 auto * addvar_mod = dynamic_cast<BlockModAdd< ColVariable > *>(mod);
 if( addvar_mod ) {
  for( auto * i : addvar_mod->added() ) {
   add_dynamic_variable( i );
  }
  return;
 }

 auto * rmvvar_mod = dynamic_cast<BlockModRmv< ColVariable > *>(mod);
 if( rmvvar_mod ) {
  for( const auto & i : rmvvar_mod->removed() ) {
   remove_dynamic_variable( &i );
  }
  return;
 }

 auto * addbnd_mod = dynamic_cast<BlockModAdd< LB0Constraint > *>(mod);
 if( addbnd_mod ) {
  for( auto * i : addbnd_mod->added() ) {
   add_dynamic_bound( i );
  }
  return;
 }

 auto * rmvbnd_mod = dynamic_cast<BlockModRmv< LB0Constraint > *>(mod);
 if( rmvbnd_mod ) {
  for( const auto & i : rmvbnd_mod->removed() ) {
   remove_dynamic_bound( &i );
  }
  return;
 }

 throw std::invalid_argument( "Unknown type of BlockAD" );
}

/*--------------------------------------------------------------------------*/

void CPXMILPSolver::add_dynamic_constraint( FRowConstraint * p_const ) {

 const auto * p_fun =
  dynamic_cast<const LinearFunction *>(p_const->get_function());
 if( p_fun == nullptr ) {
  throw std::invalid_argument( "The Constraint is not linear" );
 }

 int nzcnt = p_const->get_num_active_var();
 std::array< int, 2 > rmatbeg = { 0, nzcnt };
 std::vector< int > rmatind( nzcnt );
 std::vector< double > rmatval( nzcnt );

 std::array< double, 1 > rhs{};
 std::array< double, 1 > rngval{};
 std::array< int, 1 > indices{};
 std::array< char, 1 > sense{};

 int i = 0;

 // Get the coefficients to fill the matrix
 for( int it = 0; it < nzcnt; ++it ) {
  auto * p_var = dynamic_cast<ColVariable *>(p_fun->get_active_var( it ));
  rmatind[ i ] = index_of_variable( p_var );
  rmatval[ i ] = p_fun->get_coefficient( it );
  active_constraints[ rmatind[ i ] ].push_back( p_const );
  ++i;
 }

 auto const_lhs = p_const->get_lhs();
 auto const_rhs = p_const->get_rhs();

 int new_index = numrows++;

 if( const_lhs == const_rhs ) {
  sense[ 0 ] = 'E';
  rhs[ 0 ] = const_rhs;

 } else if( const_lhs == -Inf< double >() ) {
  sense[ 0 ] = 'L';
  rhs[ 0 ] = const_rhs;

 } else if( const_rhs == Inf< double >() ) {
  sense[ 0 ] = 'G';
  rhs[ 0 ] = const_lhs;

 } else {
  sense[ 0 ] = 'R';
  rhs[ 0 ] = const_lhs;
  rngval[ 0 ] = const_rhs - const_lhs;
  indices[ 0 ] = new_index;
 }

 v_d_const_int.emplace_back( p_const, new_index );
 v_int_d_const.emplace_back( new_index, p_const );
 std::sort( v_d_const_int.begin(), v_d_const_int.end() );
 std::sort( v_int_d_const.begin(), v_int_d_const.end() );

 CPXaddrows( env, lp, 0, 1, nzcnt, rhs.data(),
             sense.data(),
             rmatbeg.data(),
             rmatind.data(),
             rmatval.data(),
             nullptr,
             nullptr );
 if( sense[ 0 ] == 'R' ) {
  CPXchgrngval( env, lp, 1, indices.data(), rngval.data() );
 }
}

/*--------------------------------------------------------------------------*/

void CPXMILPSolver::add_dynamic_variable( ColVariable * p_var ) {

 std::vector< FRowConstraint * > var_constraints;
 std::vector< OneVarConstraint * > var_bounds;

 int nzcnt = 0;
 for( auto * stuff : p_var->active_stuff() ) {
  auto * constraint = dynamic_cast<FRowConstraint *>(stuff);
  if( constraint != nullptr ) {
   var_constraints.push_back( constraint );
   ++nzcnt;
  }
  auto * bound = dynamic_cast<OneVarConstraint *>(stuff);
  if( bound != nullptr ) {
   var_bounds.push_back( bound );
  }
 }

 std::array< int, 2 > cmatbeg = { 0, nzcnt };
 std::vector< int > cmatind( nzcnt );
 std::vector< double > cmatval( nzcnt );

 int i = 0;

 // We need the coefficients for this variable in each constraint
 for( auto * p_const : var_constraints ) {

  const auto * p_fun =
   dynamic_cast<const LinearFunction *>(p_const->get_function());
  if( p_fun == nullptr ) {
   throw std::invalid_argument( "The Constraint is not linear" );
  }
  cmatind[ i ] = index_of_constraint( p_const );
  cmatval[ i ] = p_fun->get_coefficient( i );
  ++i;
 }

 std::array< double, 1 > lb{};
 std::array< double, 1 > ub{};
 std::array< char, 1 > ctype{};

 lb[ 0 ] =
  p_var->get_lb() == -Inf< double >() ? -CPX_INFBOUND : p_var->get_lb();
 ub[ 0 ] = p_var->get_ub() == Inf< double >() ? CPX_INFBOUND : p_var->get_ub();

 for( auto * bnd : var_bounds ) {
  lb[ 0 ] = lb[ 0 ] > bnd->get_lhs() ? lb[ 0 ] : bnd->get_lhs();
  ub[ 0 ] = ub[ 0 ] < bnd->get_rhs() ? ub[ 0 ] : bnd->get_rhs();
 }

 active_constraints.emplace_back( var_constraints );
 active_bounds.emplace_back( var_bounds );

 v_d_var_int.emplace_back( p_var, numcols );
 v_int_d_var.emplace_back( numcols, p_var );
 std::sort( v_d_var_int.begin(), v_d_var_int.end() );
 std::sort( v_int_d_var.begin(), v_int_d_var.end() );
 ++numcols;

 CPXaddcols( env, lp, 1, nzcnt, nullptr, cmatbeg.data(),
             cmatind.data(), cmatval.data(), lb.data(), ub.data(), nullptr );

 // Variable type

 if( p_var->is_integer() ) {
  ++mip;
  if( p_var->is_unitary() && p_var->is_positive() ) {
   ctype[ 0 ] = 'B'; // Binary
  } else {
   ctype[ 0 ] = 'I'; // Integer
  }
 } else {
  ctype[ 0 ] = 'C';  // Continuous
 }

 // Update problem type, if necessary
 switch( mip ) {
  case 0:
   // The problem wasn't and still isn't a MIP
   break;
  case 1:
   // The first integer variable was added
   // All ctype values must be [re]added to the problem
   switch( CPXgetprobtype( env, lp ) ) {
    case CPXPROB_LP :
     CPXchgprobtype( env, lp, CPXPROB_MILP );
     break;
    case CPXPROB_QP :
     CPXchgprobtype( env, lp, CPXPROB_MIQP );
     break;
    default:
     throw std::runtime_error( "Wrong CPLEX problem type" );
   }
   break;
  default:
   // The problem stays a MIP, update only the one variable
   std::array< int, 1 > indices = { numcols - 1 };
   CPXchgctype( env, lp, 1, indices.data(), ctype.data() );
 }
}

/*--------------------------------------------------------------------------*/

void
CPXMILPSolver::add_dynamic_bound( OneVarConstraint * p_bound ) {

 auto * p_var = dynamic_cast<ColVariable *>(p_bound->get_active_var( 0 ));
 auto active_bnds = active_bounds[ index_of_variable( p_var ) ];

 // Look if the bound is already there (say, added with the Variable)
 auto it = std::find( active_bnds.begin(), active_bnds.end(), p_bound );
 if( it < active_bnds.end() ) {
  return;
 }

 // Add the bound
 active_bnds.emplace_back( p_bound );

 // Change CPLEX

 std::vector< int > indices( 2, index_of_variable( p_var ) );
 std::vector< char > lu( 2 );
 std::vector< double > bd( 2 );

 lu.resize( 2 );
 bd.resize( 2 );
 lu[ 0 ] = 'L';
 lu[ 1 ] = 'U';
 bd[ 0 ] = -CPX_INFBOUND;
 bd[ 1 ] = CPX_INFBOUND;

 for( auto * bnd : active_bnds ) {
  bd[ 0 ] = bd[ 0 ] > bnd->get_lhs() ? bd[ 0 ] : bnd->get_lhs();
  bd[ 1 ] = bd[ 1 ] < bnd->get_rhs() ? bd[ 1 ] : bnd->get_rhs();
 }

 CPXchgbds( env, lp, 2, indices.data(), lu.data(), bd.data() );
}

/*--------------------------------------------------------------------------*/

void
CPXMILPSolver::remove_dynamic_constraint( const FRowConstraint * p_const ) {

 int index = 0;
 auto it1 = find_if( v_int_d_const.begin(),
                     v_int_d_const.end(),
                     [ & ]( MILPSolver::int_const pair ) {
                      return pair.second == p_const;
                     } );
 if( it1 != v_int_d_const.end() ) {
  index = it1->first;
  v_int_d_const.erase( it1 );
 } else {
  throw std::invalid_argument( "Cannot find the Constraint" );
 }

 auto it2 = find_if( v_d_const_int.begin(),
                     v_d_const_int.end(),
                     [ & ]( MILPSolver::const_int pair ) {
                      return pair.first == p_const;
                     } );
 if( it2 != v_d_const_int.end() ) {
  v_d_const_int.erase( it2 );
 } else {
  throw std::invalid_argument( "Cannot find the Constraint" );
 }

 CPXdelrows( env, lp, index, index );
 numrows--;

 for( auto & it: v_d_const_int ) {
  if( it.second > index ) {
   it.second--;
  }
 }
 for( auto & it: v_int_d_const ) {
  if( it.first > index ) {
   it.first--;
  }
 }

 for( auto & constraints: active_constraints ) {
  auto constraint = find( constraints.begin(), constraints.end(), p_const );
  if( constraint != constraints.end() ) {
   constraints.erase( constraint );
  }
 }
}

/*--------------------------------------------------------------------------*/

void
CPXMILPSolver::remove_dynamic_variable( const ColVariable * p_var ) {

 int index = 0;
 auto it1 = find_if( v_int_d_var.begin(),
                     v_int_d_var.end(),
                     [ & ]( MILPSolver::int_var pair ) {
                      return pair.second == p_var;
                     } );
 if( it1 != v_int_d_var.end() ) {
  index = it1->first;
  v_int_d_var.erase( it1 );
 } else {
  throw std::invalid_argument( "Cannot find the Variable" );
 }

 auto it2 = find_if( v_d_var_int.begin(),
                     v_d_var_int.end(),
                     [ & ]( MILPSolver::var_int pair ) {
                      return pair.first == p_var;
                     } );

 if( it2 != v_d_var_int.end() ) {
  v_d_var_int.erase( it2 );
 } else {
  throw std::invalid_argument( "Cannot find the Variable" );
 }

 CPXdelcols( env, lp, index, index );
 numcols--;

 for( auto & it: v_d_var_int ) {
  if( it.second > index ) {
   it.second--;
  }
 }
 for( auto & it: v_int_d_var ) {
  if( it.first > index ) {
   it.first--;
  }
 }

 active_constraints.erase( active_constraints.begin() + index );
 active_bounds.erase( active_bounds.begin() + index );
}

/*--------------------------------------------------------------------------*/

void
CPXMILPSolver::remove_dynamic_bound( const OneVarConstraint * p_bound ) {
 // TODO
}

/*--------------------------------------------------------------------------*/
/*------------------- METHODS FOR HANDLING THE PARAMETERS ------------------*/
/*--------------------------------------------------------------------------*/
void CPXMILPSolver::set_par( const ThinComputeInterface::idx_type par,
                             const int value ) {
 switch( par ) {
  case intMaxIter:
   CPXsetlongparam( env, CPXPARAM_MIP_Limits_Nodes, value );
   break;
  case intMaxSol:
   CPXsetintparam( env, CPXPARAM_MIP_Limits_Solutions, value );
   break;
  case intLogVerb:
   CPXsetintparam( env, CPX_PARAM_SCRIND, value );
   break;
  case intUseCustomNames:
   use_custom_names = value; // use_custom_names is bool!
   break;
  default:
   // We assume that the symbolic constant is defined in CPLEX instead of SMS++
   CPXsetintparam( env, par, value );
 }
}

void CPXMILPSolver::set_par( ThinComputeInterface::idx_type par,
                             const double value ) {
 switch( par ) {
  case dblMaxTime:
   CPXsetdblparam( env, CPXPARAM_TimeLimit, value );
   break;
  case dblRelAcc:
   // TODO
   break;
  case dblAbsAcc:
   // TODO
   break;
  case dblUpCutOff:
   CPXsetdblparam( env, CPXPARAM_MIP_Tolerances_UpperCutoff, value );
   break;
  case dblLwCutOff:
   CPXsetdblparam( env, CPXPARAM_MIP_Tolerances_LowerCutoff, value );
   break;
  case dblRAccSol:
   CPXsetdblparam( env, CPX_PARAM_EPGAP, value );
   break;
  case dblAAccSol:
   CPXsetdblparam( env, CPX_PARAM_EPAGAP, value );
   break;
  case dblFAccSol:
   CPXsetdblparam( env, CPXPARAM_Simplex_Tolerances_Feasibility, value );
   break;
  default:
   // We assume that the symbolic constant is defined in CPLEX instead of SMS++
   CPXsetdblparam( env, par, value );
 }
}

void CPXMILPSolver::set_par( ThinComputeInterface::idx_type par,
                             const std::string & value ) {
 switch( par ) {
  case strProblemName:
   prob_name = value;
   break;
  case strOutputFile:
   output_file = value;
   break;
  default:
   // We assume that the symbolic constant is defined in CPLEX instead of SMS++
   CPXsetstrparam( env, par, value.c_str() );
 }
}

ThinComputeInterface::idx_type CPXMILPSolver::get_num_int_par() const {
 return MILPSolver::get_num_int_par() + intLastAlgParCPXS - intLastAlgParMILP;
}

ThinComputeInterface::idx_type CPXMILPSolver::get_num_str_par() const {
 return MILPSolver::get_num_str_par() + strLastAlgParCPXS - strLastAlgParMILP;
}

int CPXMILPSolver::get_int_par( idx_type par ) const {
 switch( par ) {
  case intUseCustomNames:
   return use_custom_names;
  default:
   return MILPSolver::get_int_par( par );
 }
}

const std::string &
CPXMILPSolver::get_str_par( const ThinComputeInterface::idx_type par ) const {
 switch( par ) {
  case strProblemName:
   return prob_name;
  case strOutputFile:
   return output_file;
  default:
   return MILPSolver::get_str_par( par );
 }
}

ThinComputeInterface::idx_type
CPXMILPSolver::int_par_str2idx( const std::string & name ) const {
 if( name == "intUseCustomNames" )
  return ( intUseCustomNames );
 return ( MILPSolver::str_par_str2idx( name ) );
}

const std::string &
CPXMILPSolver::int_par_idx2str( const ThinComputeInterface::idx_type idx ) const {
 // It is convoluted for extendability
 static const std::vector< std::string > pars = { "intUseCustomNames" };
 switch( idx ) {
  case intUseCustomNames:
   return pars[ 0 ];
  default:
   return MILPSolver::str_par_idx2str( idx );
 }
}

ThinComputeInterface::idx_type
CPXMILPSolver::str_par_str2idx( const std::string & name ) const {
 if( name == "strProblemName" )
  return ( strProblemName );
 if( name == "strOutputFile" )
  return ( strOutputFile );
 return ( MILPSolver::str_par_str2idx( name ) );
}

const std::string &
CPXMILPSolver::str_par_idx2str( const ThinComputeInterface::idx_type idx ) const {
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

void CPXMILPSolver::set_var_value( ColVariable & lvar, double * x, int & i ) {
#if MILPSLVR_DEBUG
 LOG( "[DEBUG] ========= MILPSolver::set_var_value():"
       << " index = " << std::setw( 4 ) << i
       << ", value = " << x[ i ] << std::endl );
#endif
 lvar.set_value( x[ i++ ] );
}

void CPXMILPSolver::set_dual_value( FRowConstraint & lconst,
                                    double * pi,
                                    int & i ) {
#if MILPSLVR_DEBUG
 LOG( "[DEBUG] ========= MILPSolver::set_dual_value():"
       << " index = " << std::setw( 4 ) << i
       << ", value = " << pi[ i ] << std::endl );
#endif
 lconst.set_dual( pi[ i++ ] );
}

// void CPXMILPSolver::fix_integer_vars() {
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
/*--------------------- End File CPXMILPSolver.cpp -------------------------*/
/*--------------------------------------------------------------------------*/
