/*--------------------------------------------------------------------------*/
/*-------------------------- File test.cpp ---------------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Main for testing MILPSolver via MCFBlock/MCFSolver. It creates one
 * MCFBlock, possibly clones it, then attaches it one MCFSolver and one
 * MILPSolver (or attaches them separately to each copy). The costs /
 * capacities / deficits of the MCFBlock are then repeatedly changed, and
 * arcs are opened / closed; if two MCFBlocks are created, Modifications are
 * mapped between the two. The two Solver are then used to solve the (two
 * allegedly identical) MCFBlock(s), and the results are compared.
 *
 * \version 0.10
 *
 * \date 07 - 03 - 2019
 *
 * \author Antonio Frangioni \n
 *         Operations Research Group \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * Copyright &copy by Antonio Frangioni
 */
/*--------------------------------------------------------------------------*/
/*------------------------------ DEFINES -----------------------------------*/
/*--------------------------------------------------------------------------*/
/* If any of the following macros is defined, then the corresponding
 * :MCFClass solver is included and the corresponding version of
 * MCFSolver<> can be tested.
 *
 * - HAVE_CSCL2      for the CS2 class
 *
 * - HAVE_CPLEX      for the MCFCplex class
 *
 * - HAVE_MFSMX      for the MCFSimplex class
 *
 * - HAVE_MFZIB      for the MCFZIB class
 *
 * - HAVE_RELAX      for the RelaxIV class
 *
 * - HAVE_SPTRE      for the MCFCplex class
 *
 * - HAVE_CPLEX      for the SPTree class; note that SPTree cannot solve
 *                   most MCF instances, except those with SPT structure
 *
 * Thus, the choice of the specific :MCFClass solver can be done in the
 * makefile with a simple -DHAVE_* argument to the compiler.
 */

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include <fstream>
#include <sstream>
#include <iomanip>

#ifdef HAVE_CSCL2
 #include "CS2.h"
 #define MCFC CS2
#endif

#ifdef HAVE_CPLEX
 #include "MCFCplex.h"
 #define MCFC MCFCplex
#endif

#ifdef HAVE_MFSMX
 #include "MCFSimplex.h"
 #define MCFC MCFSimplex
#endif

#ifdef HAVE_MFZIB
 #include "MCFZIB.h"
 #define MCFC MCFZIB
#endif

#ifdef HAVE_RELAX
 #include "RelaxIV.h"
 #define MCFC RelaxIV
#endif

#ifdef HAVE_SPTRE
 #include "SPTree.h"
 #define MCFC SPTree
#endif

#include "MCFBlock.h"
#include "MCFSolver.h"
#include "MILPSolver.h"
#include "FakeSolver.h"

/*--------------------------------------------------------------------------*/
/*-------------------------------- MACROS ----------------------------------*/
/*--------------------------------------------------------------------------*/

#define arg_2_str( s ) #s

#define solver_name( snm ) "MCFSolver<" arg_2_str( snm ) ">"

/*--------------------------------------------------------------------------*/
/*-------------------------------- USING -----------------------------------*/
/*--------------------------------------------------------------------------*/

#if( OPT_USE_NAMESPACES )
 using namespace MCFClass_di_unipi_it;
#else
 using namespace std;
#endif

using namespace SMSpp_di_unipi_it;

/*--------------------------------------------------------------------------*/
/*------------------------------- GLOBALS ----------------------------------*/
/*--------------------------------------------------------------------------*/

int mode = 0;                  // what is modified, what is solved

MCFBlock * oMCFB = nullptr;    // original MCFBlock
MCFBlock * dMCFB = nullptr;    // "derived" (i.e., R3B) MCFBlock

vector<MCFBlock::Index> open;  // names of open (and closed) arcs
int opened;                    // number of opened arcs

/*--------------------------------------------------------------------------*/
/*------------------------------ FUNCTIONS ---------------------------------*/
/*--------------------------------------------------------------------------*/

template<class T>
static inline void Str2Sthg( const char* const str , T &sthg )
{
 istringstream( str ) >> sthg;
 }

/*--------------------------------------------------------------------------*/

static inline double rndfctr( void )
{
 // return a random number between 0.5 and 2, with 50% probability of being
 // < 1
 double fctr = drand48() - 0.5;
 return( fctr < 0 ? - fctr : fctr * 4 );
 }

/*--------------------------------------------------------------------------*/

static inline void load( char * fn )
{
 ifstream iFile( fn );
 if( ! iFile ) {
  cerr << "Can't open input file " << fn << endl;
  exit( 1 );
  }

 try {
  // load the "original" MCFBlock
  iFile >> *oMCFB;

  if( mode & 4 ) {
   // if so instructed, generate abstract representation for oMCFB
   oMCFB->generate_abstract_constraints();
   oMCFB->generate_objective();
   }

  if( ( mode & 8 ) && dMCFB ) {
   // if so instructed, generate abstract representation for dMCFB (if any)
   dMCFB->generate_abstract_constraints();
   dMCFB->generate_objective();
   }
  }
 catch( exception &e ) {
  cerr << "MCFClass: " << e.what() << endl;
  exit( 1 );
  }
 catch(...) {
  cerr << "Error: unknown exception thrown" << endl;
  exit( 1 );
  }

 iFile.close();

 // open[ 0 .. opened - 1 ] = names of open arcs
 // open[ opened ... m ] = names of closed arcs
 open.resize( opened = oMCFB->get_NArcs() );  // all arcs are open
 for( int i = 0 ; i < opened ; i++ )
  open[ i ] = i;

 }  // end( load )

/*--------------------------------------------------------------------------*/

static inline bool SolveMCF( void ) 
{
 try {
  // when the modified MCFBlock is not the solved one, pass the - - - - - - -
  // Modification from one to the other
  // note that also "abstract" Modification (if any) are issued with the
  // "eNoBlck" setting to avoid that they generate a "physical" Modification;
  // this is clearly useless, as the "physical" Modification corresponding to
  // the "abstract" one must already be in the Modification queue of the
  // FakeSolver

  if( dMCFB ) {
   FakeSolver * fs = dynamic_cast<FakeSolver *>(
				  (dMCFB->get_registered_solvers()).back() );
   assert( fs );
   Lst_sp_Mod & modlist = fs->get_Modification_list();

   for( auto mod : modlist ) {
    //!! std::cout << *mod << std::endl;
    oMCFB->map_forward_Modification( dMCFB , mod );
    }

   modlist.clear();  // clear the processed Modification
   }

  // solve the MCFBlock(s) - - - - - - - - - - - - - - - - - - - - - - - - - -
  Solver * slvr1 = (oMCFB->get_registered_solvers()).front();
  Solver * slvr2 = dMCFB ? (dMCFB->get_registered_solvers()).front()
                         : (oMCFB->get_registered_solvers()).back();

  int rtrn1 = slvr1->compute( false );
  int rtrn2 = slvr2->compute( false );

  if( rtrn1 == rtrn2 ) {
   if( ( rtrn1 >= Solver::kOK ) && ( rtrn1 < Solver::kError ) ) {
    auto fo1 = slvr1->get_ub();
    auto fo2 = slvr2->get_ub();
    if( abs( fo1 - fo2 )
	<= 1e-9 *  max( double( 1 ) , abs( max( fo1 , fo2 ) ) ) ) {
     cout << "OK(f)" << endl;
     return( true );
     }
    }

   if( rtrn1 == Solver::kInfeasible ) {
    cout << "OK(e)" << endl;
    return( false );
    }

   if( rtrn1 == Solver::kUnbounded ) {
    cout << "OK(u)" << endl;
    return( false );
    }
   }

  if( ( mode & 3 ) == 1 )  // MILP attached to original
   cout << "MILPSolver = ";
  else
   cout << "MCFSolver = ";

  if( ( rtrn1 >= Solver::kOK ) && ( rtrn1 < Solver::kError ) ) {
   cout << slvr1->get_ub() << endl;
   return( false );
   }
  
  switch( rtrn1 ) {
   case( Solver::kInfeasible ):   cout << "        +INF";
                                  break;
   case( Solver::kUnbounded ):    cout << "        -INF";
                                  break;
   default:                       cout << "      Error!";
   }

  if( ( mode & 3 ) == 1 )  // MILP attached to original
   cout << "MCFSolver = ";
  else
   cout << "MILPSolver = ";

  if( ( rtrn2 >= Solver::kOK ) && ( rtrn2 < Solver::kError ) ) {
   cout << slvr2->get_ub() << endl;
   return( false );
   }
  
  switch( rtrn2 ) {
   case( Solver::kInfeasible ):   cout << "        +INF";
                                  break;
   case( Solver::kUnbounded ):    cout << "        -INF";
                                  break;
   default:                       cout << "      Error!";
   }

  cout << endl;

  return( false );
  }
 catch( exception &e ) {
  cerr << e.what() << endl;
  exit( 1 );
  }
 catch(...) {
  cerr << "Error: unknown exception thrown" << endl;
  exit( 1 );
  }
 }

/*--------------------------------------------------------------------------*/

int main( int argc , char **argv )
{
 // reading command line parameters - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 long int seed = 1;
 int wchg = 31;
 double p_change = 0.4;
 MCFBlock::Index n_change = 10;
 MCFBlock::Index n_repeat = 40;

 switch( argc ) {
  case( 8 ): Str2Sthg( argv[ 7 ] , p_change );
  case( 7 ): Str2Sthg( argv[ 6 ] , n_change );
  case( 6 ): Str2Sthg( argv[ 5 ] , n_repeat );
  case( 5 ): Str2Sthg( argv[ 4 ] , wchg );
  case( 4 ): Str2Sthg( argv[ 3 ] , mode );
  case( 3 ): Str2Sthg( argv[ 2 ] , seed );
  case( 2 ): break;
 default: cerr << "Usage: " << argv[ 0 ] <<
	   " <dmx file> [seed mode wchg #rounds #chng %chng optns]"
		<< endl <<
	   "       mode: 0 = only one, 1 = orig <- MILP, 2 = copy <- MILP"
		<< endl <<
	   "             +4 = abstract orig, +8 = abstract copy,"
		<< endl <<
 	   "             +16 = also change abstract representation"
		<< endl <<
           "       wchg: what to change, coded bit-wise "
		<< endl <<
           "             0 = cost, 1 = cap, 2 = dfct, 3 = o.arc, 4 = c.arc"
	        << endl;
	   return( 1 );
  }

 if( mode & 16 )  // if the abstract representations are changed
  mode |= 12;     // ensure they exist in the first place

 if( ( mode & 3 ) == 1 )  // MILP attached to original
  mode |= 4;              // abstract representation need be there

 if( ( mode & 3 ) == 2 )  // MILP attached to copy
  mode |= 8;              // abstract representation need be there

 // construction and loading of the objects - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 // construct the "original" MCFBlock - - - - - - - - - - - - - - - - - - - -

 oMCFB = dynamic_cast<MCFBlock *>( Block::new_Block( "MCFBlock" ) );
 assert( oMCFB );

 // load the instance - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 load( argv[ 1 ] );

 // if so instructed, construct the R3 MCFBlock = copy- - - - - - - - - - - -

 if( mode & 3 ) {
  dMCFB = dynamic_cast<MCFBlock *>( oMCFB->get_R3_Block() );
  assert( dMCFB );           // excess of caution (we know it is)

  if( mode & 4 ) {
   // if so instructed, generate abstract representation for oMCFB
   oMCFB->generate_abstract_constraints();
   oMCFB->generate_objective();
   }

  if( mode & 8 ) {
   // if so instructed, also generate abstract representation for dMCFB
   dMCFB->generate_abstract_constraints();
   dMCFB->generate_objective();
   }

  if( ( mode & 3 ) == 1 ) {  // original <- MILPSolver, copy <- MCFSolver
   oMCFB->register_Solver( Solver::new_Solver( "MILPSolver" ) );
   dMCFB->register_Solver( Solver::new_Solver( solver_name( MCFC ) ) );
   }
  else {                     // copy <- MILPSolver, original <- MCFSolver
   dMCFB->register_Solver( Solver::new_Solver( "MILPSolver" ) );
   oMCFB->register_Solver( Solver::new_Solver( solver_name( MCFC ) ) );
   }

  // attach a FakeSolver to the original one to syphoon off Modification
  oMCFB->register_Solver( Solver::new_Solver( "FakeSolver" ) );
  }
 else {                      // just use one MCFBlock
  dMCFB = nullptr;
  oMCFB->register_Solver( Solver::new_Solver( "MILPSolver" ) );
  oMCFB->register_Solver( Solver::new_Solver( solver_name( MCFC ) ) );
  }

 // compute min/max cost & max deficit- - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 MCFBlock::c_Index n = oMCFB->get_NNodes();
 MCFBlock::c_Index m = oMCFB->get_NArcs();

 cout << ", n = " << n << ", m = " << m << endl;
 if( n_change > m )
  n_change = m;

 MCFBlock::CNumber c_max = - SMSpp_di_unipi_it::Inf<MCFBlock::CNumber>();
                                                        // max cost
 MCFBlock::CNumber c_min = - c_max;                     // min cost
 MCFBlock::FNumber u_max = 0;                           // max capacity

 for( MCFBlock::Index i = 0 ; i < m ; i++ ) {
  MCFBlock::c_CNumber ci = oMCFB->get_C( i );
  if( ci < c_min )
   c_min = ci;

  if( ci > c_max )
   c_max = ci;

  MCFBlock::c_FNumber ui = oMCFB->get_U( i );
  if( ui > u_max )
   u_max = ui;
  }

 bool nzdfct = false;

 for( MCFBlock::Index i = 0 ; i < n ; )
  if( oMCFB->get_B( i++ ) > 0 ) {
   nzdfct = true;
   break;
   }

 // first solver call - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 cout << "First call: ";
 cout.setf( ios::scientific, ios::floatfield );
 cout << setprecision( 6 );

 SolveMCF();
 
 // main loop - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // now, for n_repeat times:
 // - up tp n_change costs are changed, then the two problems are re-solved;
 // - up to n_change capacities are changed, then the two problems are
 //   re-solved, then the original capacities are restored;
 // - if the problem is not a circulation problem, 2 deficits are modified
 //   (adding and subtracting the same number), then the two problems are
 //   re-solved, then the original deficits are restored;
 // - up to n_change arcs are closed, then the two problems are re-solved;
 //   the same arcs arcs are re-opened, then the two problems are re-solved

 srand48( seed );  // seed the pseudo-random number generator

 MCFBlock::Vec_CNumber newcsts( n_change );
 MCFBlock::Vec_FNumber newcaps( max( n_change , n ) );

 while( n_repeat-- ) {

  cout << "Changing: ";

  // change costs - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  if( ( wchg & 1 ) && ( drand48() <= p_change ) ) {
   MCFBlock::Index tochange = max( double( 1 ) , drand48() * n_change );
   cout << tochange << " cost";

   if( tochange == 1 ) {
    MCFBlock::CNumber newcst = c_min +
                          MCFBlock::CNumber( drand48() * ( c_max - c_min ) );

    MCFBlock::Index arc = MCFBlock::Index( drand48() * ( m - 1 ) );

    if( ( mode & 16 ) && ( drand48() < 0.5 ) ) {
     // change via abstract representation
     cout << "(a)";
     auto obj = boost::any_cast<FRealObjective *>( oMCFB->get_objective() );
     assert( obj );
     auto lf = dynamic_cast<LinearFunction *>( obj->get_function() );
     assert( lf );
     LinearFunction::v_coeff nc = { newcst };
     lf->modify_coefficients( nc.begin() , arc , arc + 1 );
     }
    else  // change via call to chg_* method
     oMCFB->chg_cost( newcst , arc );

    cout << " - ";
    }
   else {
    for( MCFBlock::Index i = 0 ; i < tochange ; i++ )
     newcsts[ i ] = c_min +
                          MCFBlock::CNumber( drand48() * ( c_max - c_min ) );

    // in 50% of the cases do a ranged change, in the others a sparse change
    if( drand48() <= 0.5 ) {
     MCFBlock::Index strt = drand48() * ( m - tochange );
     MCFBlock::Index stp = strt + tochange;

     if( ( mode & 16 ) && ( drand48() < 0.5 ) ) {
      // change via abstract representation
      cout << "s(r,a) - ";
      auto obj = boost::any_cast<FRealObjective *>( oMCFB->get_objective() );
      assert( obj );
      auto lf = dynamic_cast<LinearFunction *>( obj->get_function() );
      assert( lf );
      lf->modify_coefficients( newcsts.begin() , strt , stp );
      }
     else {  // change via call to chg_* method
      oMCFB->chg_costs( newcsts.begin() , strt , stp );
      cout << "s(r) - ";
      }
     }
    else {
     MCFBlock::Vec_Index nms( m );
     for( MCFBlock::Index i = 0 ; i < m ; i++ )
      nms[ i ] = i;

     for( MCFBlock::Index i = 0 ; i < tochange ; i++ )
      swap( nms[ i ] , nms[ i + drand48() * ( m - i ) ] );

     sort( nms.begin() , nms.begin() + tochange );
     nms.resize( tochange );

     if( ( mode & 16 ) && ( drand48() < 0.5 ) ) {
      // change via abstract representation
      cout << "s(s,a) - ";
      auto obj = boost::any_cast<FRealObjective *>( oMCFB->get_objective() );
      assert( obj );
      auto lf = dynamic_cast<LinearFunction *>( obj->get_function() );
      assert( lf );
      lf->modify_coefficients( newcsts.begin() , nms );
      }
     else {  // change via call to chg_* method
      oMCFB->chg_costs( newcsts.begin() , std::move( nms ) , true );
      cout << "s(s) - ";
      }
     }
    }
   }  // end( if( change costs ) )

  // change capacities- - - - - - - - - - - - - - - - - - - - - - - - - - - -

  if( ( wchg & 2 ) && ( drand48() <= p_change ) ) {
   MCFBlock::Index tochange = max( double( 1 ) , drand48() * n_change );
   cout << tochange << " capacit";

   if( tochange == 1 ) {
    MCFBlock::Index arc = MCFBlock::Index( drand48() * ( m - 1 ) );
    MCFBlock::CNumber newcap = oMCFB->get_U( arc ) * rndfctr();

    if( ( mode & 16 ) && ( drand48() < 0.5 ) ) {
     // change via abstract representation
     cout << "y(a) - ";
     auto bnd = boost::any_cast<std::vector<LB0Constraint> *>(
				    (oMCFB->get_static_constraints())[ 1 ] );
     assert( bnd );
     (*bnd)[ arc ].set_rhs( newcap );
     }
    else {  // change via call to chg_* method
     oMCFB->chg_ucap( newcap , arc );
     cout << "y - ";
     }
    }
   else
    // in 50% of the cases do a ranged change, in the others a sparse change
    if( drand48() <= 0.5 ) {
     MCFBlock::Index strt = drand48() * ( m - tochange );
     MCFBlock::Index stp = strt + tochange;
     for( MCFBlock::Index i = 0 ; i < tochange ; ++i )
      newcaps[ i ] = oMCFB->get_U( i + strt ) * rndfctr();

     if( ( mode & 16 ) && ( drand48() < 0.5 ) ) {
      // change via abstract representation
      cout << "ies(a,r) - ";
      auto bnd = boost::any_cast<std::vector<LB0Constraint> *>(
				    (oMCFB->get_static_constraints())[ 1 ] );
      assert( bnd );
      for( MCFBlock::Index i = 0 ; i < tochange ; ++i )
       (*bnd)[ i + strt ].set_rhs( newcaps[ i ] );
      }
     else {  // change via call to chg_* method
      oMCFB->chg_ucaps( newcaps.begin() , strt , stp );
      cout << "ies(r) - ";
      }
     }
    else {
     MCFBlock::Vec_Index nms( m );
     for( MCFBlock::Index i = 0 ; i < m ; i++ )
      nms[ i ] = i;

     for( MCFBlock::Index i = 0 ; i < tochange ; i++ ) {
      swap( nms[ i ] , nms[ i + drand48() * ( m - i ) ] );
      newcaps[ i ] = oMCFB->get_U( nms[ i ]  ) * rndfctr();
      }

     sort( nms.begin() , nms.begin() + tochange );
     nms.resize( tochange );

     if( ( mode & 16 ) && ( drand48() < 0.5 ) ) {
      // change via abstract representation
      cout << "ies(a,s) - ";
      auto bnd = boost::any_cast<std::vector<LB0Constraint> *>(
				    (oMCFB->get_static_constraints())[ 1 ] );
      assert( bnd );
      for( MCFBlock::Index i = 0 ; i < tochange ; ++i )
       (*bnd)[ nms[ i ] ].set_rhs( newcaps[ i ] );
      }
     else {  // change via call to chg_* method
      oMCFB->chg_ucaps( newcaps.begin() , std::move( nms ) , true );
      cout << "ies(s) - ";
      }
     }
   }  // end( if( change capacities ) )

  // change deficits- - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  if( ( wchg & 4 ) && ( drand48() <= p_change ) ) {
   cout << "2 deficits";

   MCFBlock::Index posn;
   MCFBlock::Index negn;
   MCFBlock::FNumber posd;
   MCFBlock::FNumber negd;

   if( nzdfct ) {  // if there are nonzero deficits
    newcaps = oMCFB->get_B();

    do
     posn = MCFBlock::Index( drand48() * n );  // select node with positive
    while( newcaps[ posn ] <= 0 );             // deficit (one must exist)
    posd = newcaps[ posn ];

    do
     negn = MCFBlock::Index( drand48() * n );  // select node with negative
    while( newcaps[ negn ] >= 0 );             // deficit (one must exist)
    negd = newcaps[ negn ];
    }
   else {
    posn = MCFBlock::Index( drand48() * n );   // just select at random
    negn = MCFBlock::Index( drand48() * n );
    posd = negd = 0;
    }

   MCFBlock::FNumber Dlt = u_max / 5;
   if( drand48() <= 0.5 ) {  // in 50% of cases up, in 50% of cases down
    posd += Dlt;
    negd -= Dlt;
    }
   else {
    Dlt = min( Dlt , max( max( posd , - negd ) / 2 , double( 1 ) ) );
    posd -= Dlt;
    negd += Dlt;
    }

   if( ( mode & 16 ) && ( drand48() < 0.5 ) ) {
    // change via abstract representation
    cout << "(a)";
    auto flw = boost::any_cast<std::vector<FRowConstraint> *>(
				    (oMCFB->get_static_constraints())[ 0 ] );
    assert( flw );
    (*flw)[ posn ].set_both( posd );
    (*flw)[ negn ].set_both( negd );
    }
   else {  // change via call to chg_* method
    oMCFB->chg_dfct( posd , posn );
    oMCFB->chg_dfct( negd , negn );
    }

   cout << " - ";

   }  // end( change deficits )

  // closing arcs- - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  if( ( wchg & 8 ) && ( drand48() <= p_change ) ) {
   // at most half of the open ones
   MCFBlock::Index tochange = min( MCFBlock::Index( opened / 2 ) ,
				   MCFBlock::Index( drand48() * n_change ) );
   if( tochange ) {
    cout << tochange << " close";

    MCFBlock::Vec_Index nms( tochange );
    for( MCFBlock::Index i = 0 ; i < tochange ; ++i ) {
     MCFBlock::Index pos = drand48() * opened;
     if( pos >= opened )
      pos = opened - 1;
     MCFBlock::Index arc = open[ pos ];
     open[ pos ] = open[ --opened ];
     open[ opened ] = arc;
     nms[ i ] = arc;
     }

    if( ( mode & 16 ) && ( drand48() < 0.5 ) ) {
     // change via abstract representation
     cout << "(a)";
     auto x = boost::any_cast<std::vector<ColVariable> *>(
				      (oMCFB->get_static_variables())[ 0 ] );
     assert( x );
     for( MCFBlock::Index i = 0 ; i < tochange ; ++i ) {
      (*x)[ nms[ i ] ].set_value( 0 );
      (*x)[ nms[ i ] ].is_fixed( true );
      }
     }
    else  // change via call to chg_* method
     oMCFB->close_arcs( std::move( nms ) );

    cout << " - ";
    }
   }

  // re-opening arcs - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  if( ( wchg & 16 ) && ( drand48() <= p_change ) ) {
   // at most half of the closed ones
   MCFBlock::Index tochange = min( MCFBlock::Index( ( m - opened ) / 2 ) ,
				   MCFBlock::Index( drand48() * n_change ) );
   if( tochange ) {
    cout << tochange << " open";

    MCFBlock::Vec_Index nms( tochange );
    for( int i = 0 ; i < tochange ; i++ ) {
     MCFBlock::Index pos = drand48() * ( m - opened );
     if( opened + pos >= m )
      pos = m - opened - 1;
     MCFBlock::Index arc = open[ opened + pos ];
     open[ opened + pos ] = open[ opened ];
     open[ opened++ ] = arc;
     nms[ i ] = arc;
     } 

    if( ( mode & 16 ) && ( drand48() < 0.5 ) ) {
     // change via abstract representation
     cout << "(a)";
     auto x = boost::any_cast<std::vector<ColVariable> *>(
				      (oMCFB->get_static_variables())[ 0 ] );
     assert( x );
     for( MCFBlock::Index i = 0 ; i < tochange ; ++i )
      (*x)[ nms[ i ] ].is_fixed( false );
     }
    else  // change via call to chg_* method
     oMCFB->open_arcs( std::move( nms ) );

    cout << " - ";
    }
   }

  // finally, re-solve the problems- - - - - - - - - - - - - - - - - - - - -
  // yet, if the problem is either unfeasible or unbounded, re-load it in
  // both MCFClass and MCFBlock

  if( ! SolveMCF() )
   load( argv[ 1 ] );

  }  // end( main loop )- - - - - - - - - - - - - - - - - - - - - - - - - - -
     // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 // destroy objects and vectors - - - - - - - - - - - - - - - - - - - - - - - 
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 delete dMCFB;
 delete oMCFB;

 // terminate - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 return( 0 );

 }  // end( main )

/*--------------------------------------------------------------------------*/
/*------------------------ End File test.cpp -------------------------------*/
/*--------------------------------------------------------------------------*/
