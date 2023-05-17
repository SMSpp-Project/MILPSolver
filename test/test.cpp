/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include <gtest/gtest.h>
#include <AbstractBlock.h>
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
/*----------------------- PARAMETERIZED TEST FIXTURE -----------------------*/
/*--------------------------------------------------------------------------*/

class MILPSolverTest :
 public ::testing::TestWithParam< TestParameter > {
 protected:

 MILPSolverTest() = default;

 ~MILPSolverTest() override = default;

 void SetUp() override {
  const auto filename = std::get< 0 >( GetParam() );
  std::ifstream istream( filename );
  if( ! istream.is_open() )
   throw( std::runtime_error( "Failed to open file " + filename ) );

  block = new AbstractBlock();
  EXPECT_TRUE( block != nullptr );
  block->load( istream );
 }

 void TearDown() override {
  delete block;
 }

 AbstractBlock * block{};
};

/*--------------------------------------------------------------------------*/
/*--------------------------- PARAMETERIZED TESTS --------------------------*/
/*--------------------------------------------------------------------------*/

TEST_P( MILPSolverTest, SimpleSolve ) {
 // Solve
 Solver * solver = new CPXMILPSolver();
 block->register_Solver( solver );
 int status = solver->compute();
 ASSERT_EQ( status, Solver::kOK );

 // Check the variable values
 solver->get_var_solution();
 const auto & opt_vars = std::get< 1 >( GetParam() );
 for( int i = 0;
      i < block->get_static_variable_v< ColVariable >( 0 )->size(); ++i ) {
  auto x = ( *block->get_static_variable_v< ColVariable >( 0 ) )[ i ]
   .get_value();
  ASSERT_NEAR( x, opt_vars[ i ], 1e-6 );
 }

 // Check the objective function value
 auto obj = dynamic_cast< FRealObjective * >( block->get_objective() );
 obj->get_function()->compute();
 auto of = obj->get_function()->get_value();
 ASSERT_NEAR( of, std::get< 2 >( GetParam() ), 1e-6 );
}

/*--------------------------------------------------------------------------*/
/*------------------------- TEST SUITE INSTANCES ---------------------------*/
/*--------------------------------------------------------------------------*/

INSTANTIATE_TEST_SUITE_P( CPXMILPSolverTests,
                          MILPSolverTest,
                          ::testing::Values(
                           TestParameter( TestFile( "test0.mps" ),
                                          OptVars{ 4, 1 },
                                          OptSolution( -22 ) ),
                           TestParameter( TestFile( "test1.mps" ),
                                          OptVars{ 0.333333, 0, 0.333333, 2 },
                                          OptSolution( 3.33333333e+00 ) ) ) );

/*--------------------------------------------------------------------------*/
/*---------------------------------- MAIN ----------------------------------*/
/*--------------------------------------------------------------------------*/

int main( int argc, char ** argv ) {
 ::testing::InitGoogleTest( &argc, argv );
 return( RUN_ALL_TESTS() );
}
