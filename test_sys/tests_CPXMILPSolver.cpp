/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include <gtest/gtest.h>
#include "SimpleMILPBlock.h"
#include "CPXMILPSolver.h"

/*--------------------------------------------------------------------------*/
/*-------------------------------- USING -----------------------------------*/
/*--------------------------------------------------------------------------*/

using namespace SMSpp_di_unipi_it;

/*--------------------------------------------------------------------------*/
/*------------------------------- GLOBALS ----------------------------------*/
/*--------------------------------------------------------------------------*/

typedef std::string TestFile;
typedef std::vector< double > OptVars;
typedef double OptSolution;
typedef std::tuple< TestFile, OptVars, OptSolution > TestParameter;

/*--------------------------------------------------------------------------*/
/*------------------------- PARAMETRIZED FIXTURE ---------------------------*/
/*--------------------------------------------------------------------------*/
class MILPSolverTest :
 public ::testing::TestWithParam< TestParameter > {
 protected:

 MILPSolverTest() = default;

 ~MILPSolverTest() override = default;

 void SetUp() override {
  block = dynamic_cast<SimpleMILPBlock *>(Block::new_Block( "SimpleMILPBlock" ));
  EXPECT_TRUE( block != nullptr );

  std::ifstream file( std::get< 0 >( GetParam() ) );
  ASSERT_TRUE( file.is_open() );
  file >> *block;
 }

 void TearDown() override {
  delete block;
 }

 SimpleMILPBlock * block{};
};

/*--------------------------------------------------------------------------*/
/*------------------------ PARAMETRIZED TEST CASES -------------------------*/
/*--------------------------------------------------------------------------*/
TEST_P( MILPSolverTest, SimpleSolve ) {

 // block->register_Solver( Solver::new_Solver( "CPXMILPSolver" ) );
 // auto solver = ( block->get_registered_solvers() ).front();

 // Solve
 Solver * solver = new CPXMILPSolver();
 block->register_Solver( solver );
 int status = solver->compute();
 ASSERT_EQ( status, Solver::kOK );

 // Check the variable values
 solver->get_var_solution();
 const auto & opt_vars = std::get< 1 >( GetParam() );
 for( int i = 0; i < block->get_x().size(); ++i ) {
  auto x = block->get_x()[ i ].get_value();
  ASSERT_EQ( x, opt_vars[ i ] );
 }

 // Check the objective function value
 auto obj = dynamic_cast<FRealObjective *>(block->get_objective());
 obj->get_function()->compute();
 auto of = obj->get_function()->get_value();
 ASSERT_EQ( of, std::get< 2 >( GetParam() ) );
}

/*--------------------------------------------------------------------------*/
/*------------------------- TEST CASE INSTANCES ----------------------------*/
/*--------------------------------------------------------------------------*/

INSTANTIATE_TEST_CASE_P( CPXMILPSolverTests,
                         MILPSolverTest,
                         ::testing::Values(
                          TestParameter( TestFile( "data/test.milp" ), OptVars{ 4, 1 }, OptSolution( -22 ) )
                         ) );

/*--------------------------------------------------------------------------*/
/*---------------------------------- MAIN ----------------------------------*/
/*--------------------------------------------------------------------------*/

int main( int argc, char ** argv ) {
 ::testing::InitGoogleTest( &argc, argv );
 return RUN_ALL_TESTS();
}
