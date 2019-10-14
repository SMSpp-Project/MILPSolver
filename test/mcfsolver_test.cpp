/*--------------------------------------------------------------------------*/
/*-------------------------- File test.cpp ---------------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Main for testing MILPSolver and MCFSolver on a MCFBlock.
 *
 * \version 3.00
 *
 * \date 20 - 08 - 2019
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
 * - HAVE_CPLEX      for the MCFCplex class
 *
 * - HAVE_SPTRE      for the SPTree class; note that SPTree cannot solve
 *                   most MCF instances, except those with SPT structure
 *
 * Thus, the choice of the specific :MCFClass solver can be done in the
 * makefile with a simple -DHAVE_* argument to the compiler.
 */

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include <fstream>
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
#include "CPXMILPSolver.h"

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

// FIXME: Avoid these declarations
template<> const std::vector<int> MCFSolver<MCFC>::Solver_2_MCFClass_int;
template<> const std::vector<int> MCFSolver<MCFC>::Solver_2_MCFClass_dbl;

/*--------------------------------------------------------------------------*/
/*------------------------------- GLOBALS ----------------------------------*/
/*--------------------------------------------------------------------------*/
unsigned int mode = 0;
MCFBlock * mcfb = nullptr;
bool isnc4 = false;
int file_counter = 0;

/*--------------------------------------------------------------------------*/
/*------------------------------ FUNCTIONS ---------------------------------*/
/*--------------------------------------------------------------------------*/

template< class T >
static inline void Str2Sthg( const char * const str, T & sthg ) {
 istringstream( str ) >> sthg;
}

/*--------------------------------------------------------------------------*/

static inline double rndfctr() {
 // return a random number between 0.5 and 2, with 50% probability of being
 // < 1
 double fctr = drand48() - 0.5;
 return ( fctr < 0 ? -fctr : fctr * 4 );
}

/*--------------------------------------------------------------------------*/

static inline void load( char * fn ) {
 try {

  if( isnc4 ) {
   netCDF::NcFile f( fn, netCDF::NcFile::read );
   if( f.isNull() ) {
    std::cerr << "cannot open nc4 file " << fn << std::endl;
    exit( 1 );
   }

   netCDF::NcGroupAtt gtype = f.getAtt( "SMS++_file_type" );
   if( gtype.isNull() ) {
    std::cerr << fn << " is not an SMS++ nc4 file" << std::endl;
    exit( 1 );
   }

   int type;
   gtype.getValues( &type );

   if( type != eBlockFile ) {
    std::cerr << fn << " is not an SMS++ nc4 Block file" << std::endl;
    exit( 1 );
   }

   netCDF::NcGroup bg = f.getGroup( "Block_0" );
   if( bg.isNull() ) {
    std::cerr << "Block_0 empty or undefined in " << fn << std::endl;
    exit( 1 );
   }

   mcfb->deserialize( bg );

  } else {
   ifstream iFile( fn );
   if( !iFile ) {
    cerr << "Can't open dmx file " << fn << endl;
    exit( 1 );
   }

   iFile.clear();
   iFile.seekg( 0 );
   iFile >> *mcfb;
  }

  mcfb->generate_abstract_constraints();
  mcfb->generate_objective();
 }

 catch( exception & e ) {
  cerr << "MCFClass: " << e.what() << endl;
  exit( 1 );
 }
 catch( ... ) {
  cerr << "Error: unknown exception thrown" << endl;
  exit( 1 );
 }
}

/*--------------------------------------------------------------------------*/

static inline bool SolveMCF() {
 try {
  Solver * milpsolver = mcfb->get_registered_solvers().front();
  Solver * mcfsolver = mcfb->get_registered_solvers().back();
  int milp_status = milpsolver->compute( false );
  int mcf_status = mcfsolver->compute( false );

  // Print problem on file
  // ofstream output_stream;
  // std::ostringstream ss;
  // ss << std::setw( 2 ) << std::setfill( '0' ) << file_counter;
  // output_stream.open( ss.str() + ".dmx" );
  // dynamic_cast<MCFSolver< MCFC > *>(mcfsolver)->WriteMCF( output_stream );
  // output_stream.close();
  // dynamic_cast<CPXMILPSolver *>(milpsolver)->write_lp( ss.str() + ".lp" );
  // cout << ss.str() << " - ";
  // file_counter++;

  // EQUAL STATUS
  if( milp_status == mcf_status ) {
   if( milp_status >= Solver::kOK && milp_status < Solver::kError ) {
    auto fo1 = milpsolver->get_ub();
    auto fo2 = mcfsolver->get_ub();
    if( abs( fo1 - fo2 ) <= 1e-9 * max( double( 1 ), abs( max( fo1, fo2 ) ) ) ) {
     cout << "\033[1;32mOK(f)\033[0m: MILPSolver = " << fo1;
     cout << ", MCFSolver = " << fo2 << endl;
     return true;
    }
   }

   if( milp_status == Solver::kInfeasible ) {
    cout << "\033[1;33mOK(e)\033[0m" << endl;
    return false;
   }

   if( milp_status == Solver::kUnbounded ) {
    cout << "\033[1;33mOK(u)\033[0m" << endl;
    return false;
   }
  }

  // STATUS DIFFERS
  cout << "\033[1;31mNot OK\033[0m: MILPSolver = ";
  switch( milp_status ) {
   case Solver::kInfeasible :
    cout << "Infeasible";
    break;
   case Solver::kUnbounded :
    cout << "Unbounded";
    break;
   default:
    cout << milpsolver->get_ub();
  }

  cout << ", MCFSolver = ";
  switch( mcf_status ) {
   case Solver::kInfeasible :
    cout << "Infeasible";
    break;
   case Solver::kUnbounded :
    cout << "Unbounded";
    break;
   default:
    cout << mcfsolver->get_ub();
  }
  cout << endl;
  return false;

 }
 catch( exception & e ) {
  cerr << e.what() << endl;
  exit( 1 );
 }
 catch( ... ) {
  cerr << "Error: unknown exception thrown" << endl;
  exit( 1 );
 }
 return true;
}

/*--------------------------------------------------------------------------*/

int main( int argc, char ** argv ) {

 long int seed = 1;
 unsigned int wchg = 127;
 MCFClass::Index n_change = 10;
 MCFClass::Index n_repeat = 40;
 int options = 1;

 switch( argc ) {
  case 8:
   Str2Sthg( argv[ 7 ], options );
  case 7:
   Str2Sthg( argv[ 6 ], n_change );
  case 6:
   Str2Sthg( argv[ 5 ], n_repeat );
  case 5:
   Str2Sthg( argv[ 4 ], wchg );
  case 4:
   Str2Sthg( argv[ 3 ], mode );
  case 3:
   Str2Sthg( argv[ 2 ], seed );
  case 2:
   break;
  default:
   cerr << "Usage: " << argv[ 0 ] <<
        " <dmx file> [seed how what #rounds #chng options]"
        << endl <<
        "       how: how to change, coded bit-wise "
        << endl <<
        "             0 = abstract (0) or physical (1) "
        << endl <<
        "             +2 = use ranged changes instead of sparse"
        << endl <<
        "       what: what to change, coded bit-wise "
        << endl <<
        "             0 = cost, 1 = cap, 2 = dfct, 3 = o.arc, 4 = c.arc"
        << endl <<
        "             5 = add arc, 6 = delete arc"
        << endl <<
        "       options: bit 0 = re-optimize, other bits MCF-specific"
        << endl <<
        "              Relax   : > 0 uses Auction"
        << endl <<
        "              Cplex   : network pricing parameter"
        << endl <<
        "              ZIB     : 1st bit == 1 ==> primal +"
        << endl <<
        "                        0 = Dantzig, 2 = First Eligible, 4 = MPP"
        << endl <<
        "              Simplex : 1st bit == 1 ==> primal +"
        << endl <<
        "                        0 = Dantzig, 2 = First Eligible, 4 = MPP"
        << endl;
   return 1;
 }

 bool abstract = false;
 bool ranged = false;

 if( mode & 1u ) {
  std::cout << "Changing physical representation ";
  abstract = false;
 } else {
  std::cout << "Changing abstract representation ";
  abstract = true;
 }

 if( mode & 2u ) {
  std::cout << "using ranged modifications" << std::endl;
  ranged = true;
 } else {
  std::cout << "using sparse modifications" << std::endl;
  ranged = false;
 }

 mcfb = dynamic_cast<MCFBlock *>( Block::new_Block( "MCFBlock" ));
 assert( mcfb );

 // Check if the file is a .dmx (default) or a .nc4 one
 std::string filename( argv[ 1 ] );
 if( filename.size() > 4 ) {
  std::string sffx = filename.substr( filename.size() - 4, 4 );
  std::string nc4( ".nc4" );

  isnc4 = std::equal( sffx.begin(), sffx.end(), nc4.begin(),
                      []( auto a, auto b ) {
                       return ( std::tolower( a ) == std::tolower( b ) );
                      } );
 }

 load( argv[ 1 ] );
 // mcfb->register_Solver( Solver::new_Solver( "CPXMILPSolver" ) );
 // mcfb->register_Solver( Solver::new_Solver( solver_name( MCFC ) ) );
 mcfb->register_Solver( new CPXMILPSolver() );
 mcfb->register_Solver( new MCFSolver<MCFC>() );

 // Compute min/max cost & max deficit
 MCFClass::Index m = mcfb->get_NArcs();
 MCFClass::Index n = mcfb->get_NNodes();

 cout << "n = " << n << ", m = " << m << endl;
 if( n_change > m ) n_change = m;
 MCFClass::CNumber c_max = -OPTtypes_di_unipi_it::Inf< MCFClass::CNumber >();
 MCFClass::CNumber c_min = -c_max;
 MCFClass::FNumber u_avg = 0;
 MCFClass::FNumber u_min = OPTtypes_di_unipi_it::Inf< MCFClass::FNumber >();

 for( MCFClass::Index i = 0; i < m; i++ ) {
  MCFClass::cCNumber ci = mcfb->get_C( i );
  if( ci < c_min ) c_min = ci;
  if( ci > c_max ) c_max = ci;
  MCFClass::cFNumber ui = mcfb->get_U( i );
  u_avg += ui;
  if( ui < u_min ) u_min = ui;
 }

 u_avg /= m;
 bool nz_deficits = false;

 for( MCFClass::Index i = 0; i < n; ) {
  if( mcfb->get_B( i++ ) > 0 ) {
   nz_deficits = true;
   break;
  }
 }

 Solver * mcfsolver = mcfb->get_registered_solvers().back();
 dynamic_cast<MCFSolver< MCFC > *>(mcfsolver)->set_par( Solver::dblAbsAcc, u_avg * 1e-8 );
 // mcfsolver->set_par( Solver::dblAbsAcc , u_avg * 1e-10 );
 cout << "First call: " << endl;
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

 srand48( seed );
 while( n_repeat-- ) {
  cout << endl << "Changing: ";

  // change costs - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  if( wchg & 1u ) {
   MCFBlock::Index tochange = max( double( 1 ), drand48() * n_change );
   cout << tochange << " cost";

   if( tochange == 1 ) {
    MCFBlock::CNumber newcst = c_min + MCFBlock::CNumber( drand48() * ( c_max - c_min ) );
    auto arc = MCFBlock::Index( drand48() * ( m - 1 ) );

    if( abstract ) {
     cout << "(a)" << std::endl;
     auto obj = dynamic_cast<FRealObjective *>( mcfb->get_objective() );
     auto lf = dynamic_cast<LinearFunction *>( obj->get_function() );
     assert( lf );
     LinearFunction::v_coeff nc = { newcst };
     lf->modify_coefficient( arc, nc.front() );
    } else {
     mcfb->chg_cost( newcst, arc );
     cout << "(s)" << std::endl;
    }

   } else {
    MCFBlock::Vec_CNumber newcsts( tochange );
    for( MCFBlock::Index i = 0; i < tochange; i++ ) {
     newcsts[ i ] = c_min + MCFBlock::CNumber( drand48() * ( c_max - c_min ) );
    }

    if( ranged ) {
     MCFBlock::Index strt = drand48() * ( m - tochange );
     MCFBlock::Index stp = strt + tochange;

     if( abstract ) {

      cout << "s(r,a)" << std::endl;
      auto obj = dynamic_cast<FRealObjective *>( mcfb->get_objective() );
      auto lf = dynamic_cast<LinearFunction *>( obj->get_function() );
      assert( lf );
      lf->modify_coefficients( std::move( newcsts ) ,
                               Function::Range( strt , stp ) );

     } else {
      mcfb->chg_costs( newcsts.begin(), Block::Range( strt, stp ) );
      cout << "s(r)" << std::endl;
     }
    } else {

     // Sparse
     Block::Subset nms( m );
     std::iota( nms.begin(), nms.end(), 0 );

     for( Block::Index i = 0; i < tochange; i++ )
      swap( nms[ i ], nms[ i + drand48() * ( m - i ) ] );

     auto end = nms.begin() + tochange;
     sort( nms.begin(), end );
     nms.resize( tochange );

     if( abstract ) {

      // Change via abstract representation
      cout << "s(s,a)" << std::endl;
      auto obj = dynamic_cast<FRealObjective *>( mcfb->get_objective() );
      auto lf = dynamic_cast<LinearFunction *>( obj->get_function() );
      assert( lf );
      lf->modify_coefficients( std::move( newcsts ),
                               std::move( nms ),
                               true );
     } else {
      mcfb->chg_costs( newcsts.begin(), std::move( nms ), true );
      cout << "s(s)" << std::endl;
     }
    }
   }
  }

  // change capacities- - - - - - - - - - - - - - - - - - - - - - - - - - - -

  if( wchg & 2u ) {
   MCFBlock::Index tochange = max( double( 1 ), drand48() * n_change );
   cout << tochange << " capacit";

   if( tochange == 1 ) {
    auto arc = MCFBlock::Index( drand48() * ( m - 1 ) );
    MCFBlock::CNumber newcap = mcfb->get_U( arc ) * rndfctr();

    if( abstract ) {

     cout << "y(a)" << std::endl;
     mcfb->i2p_ub( arc )->set_rhs( newcap );
    } else {

     mcfb->chg_ucap( newcap, arc );
     cout << "y" << std::endl;
    }
   } else {
    MCFBlock::Vec_FNumber newcaps( tochange );

    if( ranged ) {

     MCFBlock::Index strt = drand48() * ( m - tochange );
     MCFBlock::Index stp = strt + tochange;
     for( MCFBlock::Index i = 0; i < tochange; ++i )
      newcaps[ i ] = mcfb->get_U( i + strt ) * rndfctr();

     if( abstract ) {

      cout << "ies(a,r)" << std::endl;
      for( MCFBlock::Index i = 0; i < tochange; ++i )
       mcfb->i2p_ub( i + strt )->set_rhs( newcaps[ i ] );
     } else {

      mcfb->chg_ucaps( newcaps.begin(), Block::Range( strt, stp ) );
      cout << "ies(r)" << std::endl;
     }
    } else {

     // Sparse modification
     Block::Subset nms( m );
     std::iota( nms.begin(), nms.end(), 0 );

     for( Block::Index i = 0; i < tochange; i++ ) {
      swap( nms[ i ], nms[ i + drand48() * ( m - i ) ] );
      newcaps[ i ] = mcfb->get_U( nms[ i ] ) * rndfctr();
     }

     auto end = nms.begin() + tochange;
     sort( nms.begin(), end );
     nms.resize( tochange );

     if( abstract ) {

      cout << "ies(a,s)" << std::endl;
      for( MCFBlock::Index i = 0; i < tochange; ++i )
       mcfb->i2p_ub( nms[ i ] )->set_rhs( newcaps[ i ] );
     } else {

      mcfb->chg_ucaps( newcaps.begin(), std::move( nms ), true );
      cout << "ies(s)" << std::endl;
     }
    }
   }
  }

  // change deficits- - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  if( wchg & 4u ) {
   cout << "2 deficits";

   MCFClass::Index posn;
   MCFClass::Index negn;
   MCFClass::FNumber posd;
   MCFClass::FNumber negd;

   if( nz_deficits ) {
    MCFBlock::Vec_FNumber dfcts( n );

    do
     posn = MCFClass::Index( drand48() * n ); // select node with positive
    while( mcfb->get_B( posn ) <= 0 );         // deficit (one must exist)
    posd = mcfb->get_B( posn );

    do
     negn = MCFClass::Index( drand48() * n ); // select node with negative
    while( mcfb->get_B( negn ) >= 0 );         // deficit (one must exist)
    negd = mcfb->get_B( negn );
   } else {
    posn = MCFClass::Index( drand48() * n );  // just select at random
    negn = MCFClass::Index( drand48() * n );
    posd = negd = 0;
   }

   MCFClass::FNumber Dlt = std::ceil( u_avg * 2 * drand48() );
   if( drand48() <= 0.5 ) {  // in 50% of cases up, in 50% of cases down
    posd += Dlt;
    negd -= Dlt;
   } else {
    Dlt = min( Dlt, max( max( posd, -negd ) / 2, double( 1 ) ) );
    posd -= Dlt;
    negd += Dlt;
   }

   if( abstract ) {

    // Change via abstract representation
    cout << "(a)";
    mcfb->i2p_e( posn )->set_both( posd );
    mcfb->i2p_e( negn )->set_both( negd );
   } else {

    // Change via call to chg_* method
    mcfb->chg_dfct( posd, posn );
    mcfb->chg_dfct( negd, negn );
   }
   cout << "" << std::endl;

  }

  // closing arcs- - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  if( wchg & 8u ) {
   MCFBlock::Index changed = 0;

   MCFBlock::Subset nms( n_change );
   for( MCFBlock::Index i = mcfb->get_NStaticArcs();
        i < mcfb->get_NArcs(); ++i ) {
    if( mcfb->is_deleted( i ) )
     continue;
    if( mcfb->is_closed( i ) )
     continue;
    if( drand48() <= 0.5 )
     continue;

    nms[ changed++ ] = i;

    if( changed >= n_change )
     break;
   }

   if( changed ) {
    nms.resize( changed );
    cout << changed << " close";

    if( abstract ) {

     cout << "(a)";
     for( auto i : nms ) {
      auto x = mcfb->i2p_x( i );
      x->set_value( 0 );
      x->is_fixed( true );
     }
    } else {

     mcfb->close_arcs( std::move( nms ) );
    }
    cout << "" << std::endl;
   }
  }

  // re-opening arcs - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  if( wchg & 16u ) {
   MCFBlock::Index changed = 0;

   MCFBlock::Subset nms( n_change );
   for( MCFBlock::Index i = mcfb->get_NStaticArcs();
        i < mcfb->get_NArcs(); ++i ) {
    if( mcfb->is_deleted( i ) )
     continue;
    if( !mcfb->is_closed( i ) )
     continue;
    if( drand48() <= 0.5 )
     continue;

    nms[ changed++ ] = i;

    if( changed >= n_change )
     break;
   }

   if( changed ) {
    nms.resize( changed );
    cout << changed << " open";

    if( abstract ) {

     cout << "(a)";
     for( auto i : nms )
      mcfb->i2p_x( i )->is_fixed( false );
    } else {

     mcfb->open_arcs( std::move( nms ) );

    }

    cout << "" << std::endl;
   }
  }

  // deleting arcs - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  if( wchg & 32u ) {
   MCFBlock::Index changed = 0;

   if( drand48() < 0.5 ) {
    // delete somewhere in the middle

    for( MCFBlock::Index i = mcfb->get_NStaticArcs();
         i < mcfb->get_NArcs(); ++i ) {
     if( mcfb->is_deleted( i ) )
      continue;
     if( drand48() <= 0.75 )
      continue;

     mcfb->remove_arc( i );
     if( ++changed >= n_change )
      break;
    }

    if( changed )
     cout << changed << " delete(m)" << std::endl;
   } else {
    for( MCFBlock::Index i = mcfb->get_NArcs();
         --i >= mcfb->get_NStaticArcs(); ) {
     if( mcfb->is_deleted( i ) )
      continue;
     if( drand48() <= 0.13 )
      break;

     mcfb->remove_arc( i );
     if( ++changed >= n_change )
      break;
    }

    if( changed )
     cout << changed << " delete(e)" << std::endl;
   }
  }

  // creating new arcs - - - - - - - - - - - - - - - - - - - - - - - - - - -

  if( wchg & 64u ) {

   MCFBlock::Index changed = 0;
   MCFBlock::Index afterend = 0;
   while( changed < n_change ) {
    if( drand48() <= 0.13 )
     break;

    ++changed;

    MCFBlock::Index sn, en;
    do {
     sn = drand48() * mcfb->get_NNodes() + 1;
     en = drand48() * mcfb->get_NNodes() + 1;
    } while( sn == en );

    // random cost in [ - c_max , c_max ]
    auto cst = c_max * ( 1 - 2 * drand48() );

    // random capacity <= 0.75 u_avg
    auto cap = 1.5 * ( u_avg - u_min ) * drand48() + u_min;

    auto arc = mcfb->add_arc( sn, en, cst, cap );

    if( arc >= m )
     ++afterend;

    if( mcfb->get_NArcs() >= mcfb->get_MaxNArcs() )
     break;
   }

   if( changed ) {
    cout << "create " << changed << "(" << afterend << ")";
    // if( diffarcs )
    //  cout << "[d]";
    cout << "" << std::endl;
   }
  }

  // check that the status of the arcs is the same- - - - - - - - - - - - - -

  if( wchg & 120u ) {  // ... if it can ever change

   for( MCFBlock::Index i = 0; i < m; ++i ) {
    if( mcfb->is_deleted( i ) ) {
     if( !mcfb->is_deleted( i ) ) {
      std::cerr << "inconsistent del status for arc " << i << std::endl;
      exit( 1 );
     }
     continue;
    }

    if( mcfb->is_closed( i ) )
     if( !mcfb->is_closed( i ) ) {
      std::cerr << "inconsistent cls status for arc " << i << std::endl;
      exit( 1 );
     }
   }
  }

  // finally, re-solve the problems- - - - - - - - - - - - - - - - - - - - -
  // yet, if the problem is either unfeasible or unbounded, or something has
  // gone awry with the arcs names, re-load it in both MCFClass and MCFBlock

  if( !SolveMCF() ) {
   load( argv[ 1 ] );
   n = mcfb->get_NNodes();
   m = mcfb->get_NArcs();
  }
 }

 delete mcfb;
 return ( 0 );
}
