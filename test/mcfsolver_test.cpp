/*--------------------------------------------------------------------------*/
/*-------------------------- File test.cpp ---------------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Main for testing MCFBlock and MCFSolver.
 *
 * An instance of a MCF in DIMACS format is read from file in both an object
 * of class MCFC derived from MCFClass, and a MCFBlock to which a
 * MCFSolver<MCFC> is attached. The MCF problem is then repeatedly solved
 * with several changes in costs / capacities / deficits, arcs openings /
 * closures and arcs additions / deletions. The same operations are performed
 * on the two solvers, and the results are compared.
 *
 * \version 2.00
 *
 * \date 08 - 05 - 2019
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
 * - HAVE_CPLEX      for the MCFCplex class
 *
 * - HAVE_SPTRE      for the SPTree class; note that SPTree cannot solve
 *                   most MCF instances, except those with SPT structure
 *
 * Thus, the choice of the specific :MCFClass solver can be done in the
 * makefile with a simple -DHAVE_* argument to the compiler.
 */

#define NMS_IS_USED 0

// if NMS_IS_USED > 0, then the Chg****() routines are fed with a
// non-consecutive set of names; otherwise, all the involved arcs are
// consecutive

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
#include "CPXMILPSolver.h"
#include "FakeSolver.h"

/*--------------------------------------------------------------------------*/
/*-------------------------------- MACROS ----------------------------------*/
/*--------------------------------------------------------------------------*/

#define arg_2_str(s) #s

#define solver_name(snm) "MCFSolver<" arg_2_str( snm ) ">"

/*--------------------------------------------------------------------------*/
/*-------------------------------- USING -----------------------------------*/
/*--------------------------------------------------------------------------*/

#if(OPT_USE_NAMESPACES)
using namespace MCFClass_di_unipi_it;
#else
using namespace std;
#endif

using namespace SMSpp_di_unipi_it;

/*--------------------------------------------------------------------------*/
/*------------------------------- GLOBALS ----------------------------------*/
/*--------------------------------------------------------------------------*/

int mode = 0;                  // what is modified, what is solved

MCFBlock* oMCFB = nullptr;    // original MCFBlock
// MCFBlock * dMCFB = nullptr;    // "derived" (i.e., R3B) MCFBlock
MCFBlock* sMCFB = nullptr;    // MCFBlock that is solved
MCFBlock* mMCFB = nullptr;    // MCFBlock that is modified

MCFClass* mcf;                // the MCFClass object

bool isnc4 = false;            // true if the file is a ntCDF one

/*--------------------------------------------------------------------------*/
/*------------------------------ FUNCTIONS ---------------------------------*/
/*--------------------------------------------------------------------------*/

template<class T>
static inline void Str2Sthg(const char* const str, T& sthg) {
 istringstream(str) >> sthg;
}

/*--------------------------------------------------------------------------*/

static inline double rndfctr(void) {
 // return a random number between 0.5 and 2, with 50% probability of being
 // < 1
 double fctr = drand48() - 0.5;
 return (fctr < 0 ? -fctr : fctr*4);
}

/*--------------------------------------------------------------------------*/

static inline void CreateProb(int Optns) {
 bool reoptmz = Optns & 1;
 Optns /= 2;

 mcf = nullptr;  // unknown solver, or the required solver is not
 // available due to the macroes settings

#ifdef HAVE_RELAX  // - - - - - - - - - - - - - - - - - - - - - - - - - - -
 RelaxIV *rlx = new RelaxIV();
#if( AUCTION )
  if( Optns )
   rlx->SetPar( RelaxIV::kAuction , MCFClass::kYes );
#endif
 mcf = rlx;
 cout << "RelaxIV";
#endif
#ifdef HAVE_SPTRE  // - - - - - - - - - - - - - - - - - - - - - - - - - - -
 mcf = new SPTree();
 cout << "SPTree";
 assert( false );  // SPTree not fully supported yet
#endif
#ifdef HAVE_CPLEX  // - - - - - - - - - - - - - - - - - - - - - - - - - - -
 MCFCplex *cpx = new MCFCplex();
 if( Optns >= 0 )
  cpx->SetPar( CPX_PARAM_NETPPRIIND , Optns );
 mcf = cpx;
 cout << "MCFCplex";
 assert( false );  // MCFCplex not fully supported yet
#endif
#ifdef HAVE_MFZIB  // - - - - - - - - - - - - - - - - - - - - - - - - - - -
 MCFZIB *zib = new MCFZIB();
 bool PrmlSmplx = Optns & 1;
 char Prcng;
 switch( Optns / 2 ) {
  case( 0 ): Prcng = char( MCFZIB::kDantzig ); break;
  case( 1 ): Prcng = char( MCFZIB::kFrstElA ); break;
  default:   Prcng = char( MCFZIB::kMltPrPr );
  }
 if( ( ! PrmlSmplx ) && ( Prcng == MCFZIB::kDantzig ) )
  Prcng = char( MCFZIB::kMltPrPr );
 zib->SetAlg( PrmlSmplx , Prcng );
 mcf = zib;
 cout << "MCFZIB";
 assert( false );  // MCFZIB not fully supported yet
#endif
#ifdef HAVE_CSCL2  // - - - - - - - - - - - - - - - - - - - - - - - - - - -
 mcf = new CS2();
 cout << "CS2";
 assert( false );  // CS2 not fully supported yet
#endif
#ifdef HAVE_MFSMX  // - - - - - - - - - - - - - - - - - - - - - - - - - - -
 MCFSimplex* mcfs = new MCFSimplex();
 bool PrmlSmplx = Optns & 1;
 char Prcng;
 switch (Optns/2) {
  case (0):
   Prcng = char(MCFSimplex::kDantzig);
   break;
  case (1):
   Prcng = char(MCFSimplex::kFirstEligibleArc);
   break;
  default:
   Prcng = char(MCFSimplex::kCandidateListPivot);
 }
 if ((!PrmlSmplx) && (Prcng == MCFSimplex::kDantzig))
  Prcng = char(MCFSimplex::kCandidateListPivot);
 mcf = mcfs;
 cout << "MCFSimplex";
#endif

 if (!reoptmz)
  mcf->SetPar(MCFClass::kReopt, MCFClass::kNo);

 if (!mcf)
  throw (std::logic_error("No MCFClass defined"));
}

/*--------------------------------------------------------------------------*/

static inline void load(char* fn) {
 try {
  // note that usually the "original" MCFBlock is loaded, unless mode == 2
  // (==> original solved, R3 modified) *and* the R3 as been constructed
  // already, in which case the R3 is loaded; this is perhaps stretching
  // the concept of "R3Block" close to the breaking point, but in this case
  // it works because the R3B is a copy *and* always the same instance is
  // loaded

  if (isnc4) {
   netCDF::NcFile f(fn, netCDF::NcFile::read);
   if (f.isNull()) {
    std::cerr << "cannot open nc4 file " << fn << std::endl;
    exit(1);
   }

   netCDF::NcGroupAtt gtype = f.getAtt("SMS++_file_type");
   if (gtype.isNull()) {
    std::cerr << fn << " is not an SNS++ nc4 file" << std::endl;
    exit(1);
   }

   int type;
   gtype.getValues(&type);

   if (type != eBlockFile) {
    std::cerr << fn << " is not an SNS++ nc4 Block file" << std::endl;
    exit(1);
   }

   netCDF::NcGroup bg = f.getGroup("Block_0");
   if (bg.isNull()) {
    std::cerr << "Block_0 empty or undefined in " << fn << std::endl;
    exit(1);
   }

   // if( ( ( mode & 3 ) == 2 ) && dMCFB ) {
   //  dMCFB->deserialize( bg );  // load the (derived) MCFBlock
   //
   //   // load the MCFClass out of the MCFBlock using the in-memory interface
   //   mcf->LoadNet( dMCFB->get_MaxNNodes() , dMCFB->get_MaxNArcs() ,
   //   dMCFB->get_NNodes() , dMCFB->get_NArcs() ,
   //   dMCFB->get_U().empty() ? nullptr : dMCFB->get_U().data() ,
   //   dMCFB->get_C().empty() ? nullptr : dMCFB->get_C().data() ,
   //   dMCFB->get_B().empty() ? nullptr : dMCFB->get_B().data() ,
   //   dMCFB->get_SN().data() , dMCFB->get_EN().data() );
   //   }
   //  else {
   oMCFB->deserialize(bg);  // load the (original) MCFBlock

   // load the MCFClass out of the MCFBlock using the in-memory interface
   mcf->LoadNet(oMCFB->get_MaxNNodes(), oMCFB->get_MaxNArcs(),
                oMCFB->get_NNodes(), oMCFB->get_NArcs(),
                oMCFB->get_U().empty() ? nullptr : oMCFB->get_U().data(),
                oMCFB->get_C().empty() ? nullptr : oMCFB->get_C().data(),
                oMCFB->get_B().empty() ? nullptr : oMCFB->get_B().data(),
                oMCFB->get_SN().data(), oMCFB->get_EN().data());
   // }
  } else {
   ifstream iFile(fn);
   if (!iFile) {
    cerr << "Can't open dmx file " << fn << endl;
    exit(1);
   }

   mcf->LoadDMX(iFile);  // load the MCFClass

   iFile.clear();
   iFile.seekg(0);       // rewind the file

   // if( ( ( mode & 3 ) == 2 ) && dMCFB )
   //  iFile >> *dMCFB;       // load the (derived) MCFBlock
   // else
   iFile >> *oMCFB;       // load the (original) MCFBlock
  }

  // immediately run the Modification mapping: this should only map the
  // single NBModification generated by the loading, and it ensures that
  // the MCFBlock that is *not* directly loaded still receives the state
  // of a "freshly minted" MCFBlock before any other Modification, which
  // is necessary for the MCFSolver (see the comments to
  // MCFSolver::add_Modification()

  // if( ( mode & 3 ) && dMCFB ) {
  //  FakeSolver * fs = dynamic_cast<FakeSolver *>(
  // 		  (mMCFB->get_registered_solvers()).front() );
  //  assert( fs );
  //  Lst_sp_Mod & modlist = fs->get_Modification_list();
  //
  //  if( mMCFB == oMCFB ) {  // modify the original, solve the R3
  //   for( auto mod : modlist )
  //    oMCFB->map_forward_Modification( dMCFB , mod );
  //   }
  //  else {                  // modify the R3, solve the original
  //   for( auto mod : modlist )
  //    oMCFB->map_back_Modification( dMCFB , mod );
  //   }
  //
  //  modlist.clear();  // clear the processed Modification
  //  }

  if (mode & 4) {
   // if so instructed, generate abstract representation for oMCFB
   oMCFB->generate_abstract_constraints();
   oMCFB->generate_objective();
  }

  // if( ( mode & 8 ) && dMCFB ) {
  //  // if so instructed, generate abstract representation for dMCFB (if any)
  //  dMCFB->generate_abstract_constraints();
  //  dMCFB->generate_objective();
  //  }
 }
 catch (exception& e) {
  cerr << "MCFClass: " << e.what() << endl;
  exit(1);
 }
 catch (...) {
  cerr << "Error: unknown exception thrown" << endl;
  exit(1);
 }
}  // end( load )

/*--------------------------------------------------------------------------*/

static inline bool SolveMCF(void) {
 try {
  // solve the MCFClass- - - - - - - - - - - - - - - - - - - - - - - - - - - -
  mcf->SolveMCF();
  auto stat = mcf->MCFGetStatus();

  // when the modified MCFBlock is not the solved one, pass the - - - - - - -
  // Modification from one to the other
  // note that also "abstract" Modification (if any) are issued with the
  // "eNoBlck" setting to avoid that they generate a "physical" Modification;
  // this is clearly useless, as the "physical" Modification corresponding to
  // the "abstract" one must already be in the Modification queue of the
  // FakeSolver

  // if( sMCFB != mMCFB ) {
  //  FakeSolver * fs = dynamic_cast<FakeSolver *>(
  // 		  (mMCFB->get_registered_solvers()).front() );
  //  assert( fs );
  //  Lst_sp_Mod & modlist = fs->get_Modification_list();
  //
  //  if( mMCFB == oMCFB ) {  // modify the original, solve the R3
  //   for( auto mod : modlist ) {
  //    //!! std::cout << *mod << std::endl;
  //    oMCFB->map_forward_Modification( dMCFB , mod );
  //    }
  //   }
  //  else {                  // modify the R3, solve the original
  //   for( auto mod : modlist ) {
  //    //!! std::cout << *mod << std::endl;
  //    oMCFB->map_back_Modification( dMCFB , mod );
  //    }
  //   }
  //
  //  modlist.clear();  // clear the processed Modification
  //  }

  // solve the MCFBlock- - - - - - - - - - - - - - - - - - - - - - - - - - - -
  Solver* milpsolver = sMCFB->get_registered_solvers().front();
  Solver* mcfsolver = sMCFB->get_registered_solvers().back();
  int milp_status = milpsolver->compute(false);
  int mcf_status = mcfsolver->compute(false);

  cout << endl;

  // EQUAL STATUS
  if (milp_status == mcf_status) {
   if (milp_status >= Solver::kOK && milp_status < Solver::kError) {
    auto fo1 = milpsolver->get_ub();
    auto fo2 = mcfsolver->get_ub();
    if (abs(fo1 - fo2) <= 1e-9*max(double(1), abs(max(fo1, fo2)))) {
     cout << "OK(f): MILPSolver = " << fo1;
     cout << ", MCFSolver = " << fo2 << endl;
     return (true);
    }
   }

   if (milp_status == Solver::kInfeasible) {
    cout << "OK(e)" << endl;
    return false;
   }

   if (milp_status == Solver::kUnbounded) {
    cout << "OK(u)" << endl;
    return false;
   }
  }

  // STATUS DIFFERS




  // if ((stat == MCFClass::kOK) &&
  //     (rtrn1 >= Solver::kOK) && (rtrn1 < Solver::kError)) {
  //  auto fo1 = mcf->MCFGetFO();
  //  auto fo2 = slvr1->get_ub();
  //  if (abs(fo1 - fo2)
  //      <= 1e-9*max(double(1), abs(max(fo1, fo2)))) {
  //   cout << "OK(f)" << endl;
  //   return (true);
  //  }
  // }

  // if ((stat == MCFClass::kUnfeasible) &&
  //     (rtrn1 == Solver::kInfeasible)) {
  //  cout << "OK(e)" << endl;
  //  return (false);
  // }

  // if ((stat == MCFClass::kUnbounded) &&
  //     (rtrn1 == Solver::kUnbounded)) {
  //  cout << "OK(u)" << endl;
  //  return (false);
  // }

  // cout << "MCFClass = ";
  // switch (mcf->MCFGetStatus()) {
  //  case (MCFClass::kOK):
  //   cout << mcf->MCFGetFO();
  //   break;
  //  case (MCFClass::kUnfeasible):
  //   cout << "        +INF";
  //   break;
  //  case (MCFClass::kUnbounded):
  //   cout << "        -INF";
  //   break;
  //  default:
  //   cout << "      Error!";
  // }

  // if (milp_status >= Solver::kOK && milp_status < Solver::kError) {
   cout << "Not OK: MILPSolver = " << milpsolver->get_ub();
  // } else {
  //  if (milp_status == Solver::kInfeasible) {
  //   cout << "MILPSolver =         +INF" << endl;
  //  } else if (milp_status == Solver::kUnbounded) {
  //   cout << "MILPSolver =         -INF" << endl;
  //  } else {
  //   cout << "MILPSolver =       Error!" << endl;
  //  }
  // }

  // if (mcf_status >= Solver::kOK && mcf_status < Solver::kError) {
   cout << ", MCFSolver = " << mcfsolver->get_ub() << endl;
  // } else {
  //  if (mcf_status == Solver::kInfeasible) {
  //   cout << "MCFSolver =         +INF" << endl;
  //  } else if (mcf_status == Solver::kUnbounded) {
  //   cout << "MCFSolver =         -INF" << endl;
  //  } else {
  //   cout << "MCFSolver =       Error!" << endl;
  //  }
  // }
 }
 catch (exception& e) {
  cerr << e.what() << endl;
  exit(1);
 }
 catch (...) {
  cerr << "Error: unknown exception thrown" << endl;
  exit(1);
 }
}

/*--------------------------------------------------------------------------*/

int main(int argc, char** argv) {
 // reading command line parameters - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 long int seed = 1;
 int wchg = 127;
 double p_change = 0.5;
 MCFClass::Index n_change = 10;
 MCFClass::Index n_repeat = 40;
 int optns = 1;

 switch (argc) {
  case (9):
   Str2Sthg(argv[8], optns);
  case (8):
   Str2Sthg(argv[7], p_change);
  case (7):
   Str2Sthg(argv[6], n_change);
  case (6):
   Str2Sthg(argv[5], n_repeat);
  case (5):
   Str2Sthg(argv[4], wchg);
  case (4):
   Str2Sthg(argv[3], mode);
  case (3):
   Str2Sthg(argv[2], seed);
  case (2):
   break;
  default:
   cerr << "Usage: " << argv[0] <<
        " <dmx file> [seed mode wchg #rounds #chng %chng optns]"
        << endl <<
        "       mode: 0 = only one, 1 = orig -> solve, 2 = solve -> orig"
        << endl <<
        "             +4 = abstract orig, +8 = abstract solve,"
        << endl <<
        "             +16 = also change abstract representation"
        << endl <<
        "       wchg: what to change, coded bit-wise "
        << endl <<
        "             0 = cost, 1 = cap, 2 = dfct, 3 = o.arc, 4 = c.arc"
        << endl <<
        "             5 = add arc, 6 = delete arc"
        << endl <<
        "       optns: bit 0 = re-optimize, other bits MCF-specific"
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
   return (1);
 }

 mode = 16; // SORRYNOTSORRY

 if ((mode & 3) == 3) {
  std::cerr << "Wrong mode (must be 0, 1 or 2)" << std::endl;
  exit(1);
 }

 if (mode & 16)  // if the abstract representations are changed
  mode |= 12;     // ensure they exist in the first place

 // construction and loading of the objects - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 // construct MCFClass- - - - - - - - - - - - - - - - - - - - - - - - - - - -

 CreateProb(optns);

 // construct the "original" MCFBlock - - - - - - - - - - - - - - - - - - - -

 oMCFB = dynamic_cast<MCFBlock*>( Block::new_Block("MCFBlock"));
 assert(oMCFB);

 // load the instance - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // first of check if the file is a .dmx (default) or a .nc4 one

 std::string name(argv[1]);
 if (name.size() > 4) {
  std::string sffx = name.substr(name.size() - 4, 4);
  std::string nc4(".nc4");

  isnc4 = std::equal(sffx.begin(), sffx.end(), nc4.begin(),
                     [](auto a, auto b) {
                      return (std::tolower(a) == std::tolower(b));
                     });
 }

 load(argv[1]);

 // if so instructed, construct the R3 MCFBlock = copy- - - - - - - - - - - -

 // if( mode & 3 ) {
 //  dMCFB = dynamic_cast<MCFBlock *>( oMCFB->get_R3_Block() );
 //  assert( dMCFB );           // excess of caution (we know it is)
 //
 //  if( ( mode & 8 ) && dMCFB ) {
 //   // if so instructed, also generate abstract representation for dMCFB
 //   dMCFB->generate_abstract_constraints();
 //   dMCFB->generate_objective();
 //   }
 //
 //  if( ( mode & 3 ) == 1 ) {  // modify the original, solve the R3
 //   mMCFB = oMCFB;
 //   sMCFB = dMCFB;
 //   }
 //  else {                     // modify the R3, solve the original
 //   mMCFB = dMCFB;
 //   sMCFB = oMCFB;
 //   }
 //
 //  // attach a FakeSolver to the modified one to syphoon off Modification
 //  mMCFB->register_Solver( Solver::new_Solver( "FakeSolver" ) );
 //  }
 // else                        // just use one MCFBlock
 sMCFB = mMCFB = oMCFB;

 //  attach a "true" MCFSolver to the one that is actually solved
 sMCFB->register_Solver(Solver::new_Solver("CPXMILPSolver"));
 sMCFB->register_Solver(Solver::new_Solver(solver_name(MCFC)));

 // compute min/max cost & max deficit- - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 MCFClass::Index n = mcf->MCFn();
 MCFClass::Index m = mcf->MCFm();

 cout << ", n = " << n << ", m = " << m << endl;
 if (n_change > m)
  n_change = m;

 MCFClass::CNumber c_max = -OPTtypes_di_unipi_it::Inf<MCFClass::CNumber>();
 // max cost
 MCFClass::CNumber c_min = -c_max;                     // min cost
 MCFClass::FNumber u_avg = 0;                           // average capacity
 MCFClass::FNumber u_min = OPTtypes_di_unipi_it::Inf<MCFClass::FNumber>();

 for (MCFClass::Index i = 0; i < m; i++) {
  MCFClass::cCNumber ci = mcf->MCFCost(i);
  if (ci < c_min)
   c_min = ci;

  if (ci > c_max)
   c_max = ci;

  MCFClass::cFNumber ui = mcf->MCFUCap(i);
  u_avg += ui;
  if (ui < u_min)
   u_min = ui;
 }

 u_avg /= m;
 bool nzdfct = false;

 for (MCFClass::Index i = 0; i < n;)
  if (mcf->MCFDfct(i++) > 0) {
   nzdfct = true;
   break;
  }

 // first solver call - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 cout << "First call: ";
 cout.setf(ios::scientific, ios::floatfield);
 cout << setprecision(6);

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

 srand48(seed);  // seed the pseudo-random number generator

 bool diffarcs = false;  // whether added arcs ended up with different names

 while (n_repeat--) {

  cout << "Changing: ";

  // change costs - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  if ((wchg & 1) && (drand48() <= p_change)) {
   MCFBlock::Index tochange = max(double(1), drand48()*n_change);
   cout << tochange << " cost";

   if (tochange == 1) {
    MCFBlock::CNumber newcst = c_min +
                               MCFBlock::CNumber(drand48()*(c_max - c_min));

    MCFBlock::Index arc = MCFBlock::Index(drand48()*(m - 1));

    mcf->ChgCost(arc, newcst);

    if ((mode & 16) && (drand48() < 0.5)) {
     // change via abstract representation
     cout << "(a)";
     auto obj = boost::any_cast<FRealObjective*>(mMCFB->get_objective());
     assert(obj);
     auto lf = dynamic_cast<LinearFunction*>( obj->get_function());
     assert(lf);
     LinearFunction::v_coeff nc = {newcst};
     lf->modify_coefficient(mMCFB->i2p_x(arc), nc.front());
    } else  // change via call to chg_* method
     mMCFB->chg_cost(newcst, arc);

    cout << " - ";
   } else {
    MCFBlock::Vec_CNumber newcsts(tochange);
    for (MCFBlock::Index i = 0; i < tochange; i++)
     newcsts[i] = c_min +
                  MCFBlock::CNumber(drand48()*(c_max - c_min));

    // in 50% of the cases do a ranged change, in the others a sparse change
    if (drand48() <= 0.5) {
     MCFBlock::Index strt = drand48()*(m - tochange);
     MCFBlock::Index stp = strt + tochange;
     mcf->ChgCosts(newcsts.data(), nullptr, strt, stp);

     if ((mode & 16) && (drand48() < 0.5)) {
      // change via abstract representation
      cout << "s(r,a) - ";
      auto obj = boost::any_cast<FRealObjective*>(mMCFB->get_objective());
      assert(obj);
      auto lf = dynamic_cast<LinearFunction*>( obj->get_function());
      assert(lf);
      if (mMCFB->HasDynamicX()) {
       LinearFunction::v_coeff_pair chg(tochange);
       for (MCFBlock::Index i = 0; i < tochange; ++i) {
        chg[i].first = mMCFB->i2p_x(i + strt);
        chg[i].second = newcsts[i];
       }
       // chg is ordered only if all arcs are static
       lf->modify_coefficients(chg, stp <= mMCFB->get_NStaticArcs());
      } else {
       // note that all static Variable are consecutive, but not
       // necessarily the first ones as dynamic Variable may come first
       auto dlt = lf->is_active(mMCFB->i2p_x(0));
       strt += dlt;
       stp += dlt;
       lf->modify_coefficients(newcsts.begin(), strt, stp);
      }
     } else {  // change via call to chg_* method
      mMCFB->chg_costs(newcsts.begin(), Block::Range(strt, stp));
      cout << "s(r) - ";
     }
    } else {
     MCFBlock::Subset nms(m + 1);
     for (MCFBlock::Index i = 0; i < m; i++)
      nms[i] = i;

     for (MCFBlock::Index i = 0; i < tochange; i++)
      swap(nms[i], nms[i + drand48()*(m - i)]);

     auto end = nms.begin() + tochange;
     sort(nms.begin(), end);
     *end = OPTtypes_di_unipi_it::Inf<MCFClass::Index>();
     mcf->ChgCosts(newcsts.data(), nms.data());
     nms.resize(tochange);

     if ((mode & 16) && (drand48() < 0.5)) {
      // change via abstract representation
      cout << "s(s,a) - ";
      auto obj = boost::any_cast<FRealObjective*>(mMCFB->get_objective());
      assert(obj);
      auto lf = dynamic_cast<LinearFunction*>( obj->get_function());
      assert(lf);
      if (mMCFB->HasDynamicX()) {
       LinearFunction::v_coeff_pair chg(tochange);
       for (MCFBlock::Index i = 0; i < tochange; ++i) {
        chg[i].first = mMCFB->i2p_x(nms[i]);
        chg[i].second = newcsts[i];
       }
       // chg is ordered only if all arcs are static
       lf->modify_coefficients(chg, nms.back() < mMCFB->get_NStaticArcs());
      } else
       lf->modify_coefficients(newcsts.begin(), nms);
     } else {  // change via call to chg_* method
      mMCFB->chg_costs(newcsts.begin(), std::move(nms), true);
      cout << "s(s) - ";
     }
    }
   }
  }  // end( if( change costs ) )

  // change capacities- - - - - - - - - - - - - - - - - - - - - - - - - - - -

  if ((wchg & 2) && (drand48() <= p_change)) {
   MCFBlock::Index tochange = max(double(1), drand48()*n_change);
   cout << tochange << " capacit";

   if (tochange == 1) {
    MCFBlock::Index arc = MCFBlock::Index(drand48()*(m - 1));
    MCFBlock::CNumber newcap = mcf->MCFUCap(arc)*rndfctr();
    mcf->ChgUCap(arc, newcap);

    if ((mode & 16) && (drand48() < 0.5)) {
     // change via abstract representation
     cout << "y(a) - ";
     mMCFB->i2p_ub(arc)->set_rhs(newcap);
    } else {  // change via call to chg_* method
     mMCFB->chg_ucap(newcap, arc);
     cout << "y - ";
    }
   } else {
    MCFBlock::Vec_FNumber newcaps(tochange);

    // in 50% of the cases do a ranged change, in the others a sparse change
    if (drand48() <= 0.5) {
     MCFBlock::Index strt = drand48()*(m - tochange);
     MCFBlock::Index stp = strt + tochange;
     for (MCFBlock::Index i = 0; i < tochange; ++i)
      newcaps[i] = mcf->MCFUCap(i + strt)*rndfctr();
     mcf->ChgUCaps(newcaps.data(), nullptr, strt, stp);

     if ((mode & 16) && (drand48() < 0.5)) {
      // change via abstract representation
      cout << "ies(a,r) - ";
      for (MCFBlock::Index i = 0; i < tochange; ++i)
       mMCFB->i2p_ub(i + strt)->set_rhs(newcaps[i]);
     } else {  // change via call to chg_* method
      mMCFB->chg_ucaps(newcaps.begin(), Block::Range(strt, stp));
      cout << "ies(r) - ";
     }
    } else {
     MCFBlock::Subset nms(m + 1);
     for (MCFBlock::Index i = 0; i < m; i++)
      nms[i] = i;

     for (MCFBlock::Index i = 0; i < tochange; i++) {
      swap(nms[i], nms[i + drand48()*(m - i)]);
      newcaps[i] = mcf->MCFUCap(nms[i])*rndfctr();
     }

     auto end = nms.begin() + tochange;
     sort(nms.begin(), end);
     *end = OPTtypes_di_unipi_it::Inf<MCFClass::Index>();
     mcf->ChgUCaps(newcaps.data(), nms.data());
     nms.resize(tochange);

     if ((mode & 16) && (drand48() < 0.5)) {
      // change via abstract representation
      cout << "ies(a,s) - ";
      for (MCFBlock::Index i = 0; i < tochange; ++i)
       mMCFB->i2p_ub(nms[i])->set_rhs(newcaps[i]);
     } else {  // change via call to chg_* method
      mMCFB->chg_ucaps(newcaps.begin(), std::move(nms), true);
      cout << "ies(s) - ";
     }
    }
   }
  }  // end( if( change capacities ) )

  // change deficits- - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  if ((wchg & 4) && (drand48() <= p_change)) {
   cout << "2 deficits";

   MCFClass::Index posn;
   MCFClass::Index negn;
   MCFClass::FNumber posd;
   MCFClass::FNumber negd;

   if (nzdfct) {  // if there are nonzero deficits
    MCFBlock::Vec_FNumber dfcts(n);
    mcf->MCFDfcts(dfcts.data());

    do
     posn = MCFClass::Index(drand48()*n);  // select node with positive
    while (dfcts[posn] <= 0);               // deficit (one must exist)
    posd = dfcts[posn];

    do
     negn = MCFClass::Index(drand48()*n);  // select node with negative
    while (dfcts[negn] >= 0);               // deficit (one must exist)
    negd = dfcts[negn];
   } else {
    posn = MCFClass::Index(drand48()*n);   // just select at random
    negn = MCFClass::Index(drand48()*n);
    posd = negd = 0;
   }

   MCFClass::FNumber Dlt = u_avg*2*drand48();
   if (drand48() <= 0.5) {  // in 50% of cases up, in 50% of cases down
    posd += Dlt;
    negd -= Dlt;
   } else {
    Dlt = min(Dlt, max(max(posd, -negd)/2, double(1)));
    posd -= Dlt;
    negd += Dlt;
   }

   mcf->ChgDfct(posn, posd);
   mcf->ChgDfct(negn, negd);

   if ((mode & 16) && (drand48() < 0.5)) {
    // change via abstract representation
    cout << "(a)";
    mMCFB->i2p_e(posn)->set_both(posd);
    mMCFB->i2p_e(negn)->set_both(negd);
   } else {  // change via call to chg_* method
    mMCFB->chg_dfct(posd, posn);
    mMCFB->chg_dfct(negd, negn);
   }

   cout << " - ";

  }  // end( change deficits )

  // closing arcs- - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  if ((wchg & 8) && (drand48() <= p_change)) {
   MCFBlock::Index changed = 0;

   MCFBlock::Subset nms(n_change);
   for (MCFBlock::Index i = mMCFB->get_NStaticArcs();
        i < mMCFB->get_NArcs(); ++i) {
    if (mcf->IsDeletedArc(i))
     continue;
    if (mcf->IsClosedArc(i))
     continue;
    if (drand48() <= 0.5)
     continue;

    nms[changed++] = i;
    mcf->CloseArc(i);

    if (changed >= n_change)
     break;
   }

   if (changed) {
    nms.resize(changed);
    cout << changed << " close";

    if ((mode & 16) && (drand48() < 0.5)) {
     // change via abstract representation
     cout << "(a)";
     for (auto i : nms) {
      auto x = mMCFB->i2p_x(i);
      x->set_value(0);
      x->is_fixed(true);
     }
    } else  // change via call to chg_* method
     mMCFB->close_arcs(std::move(nms));

    cout << " - ";
   }
  }

  // re-opening arcs - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  if ((wchg & 16) && (drand48() <= p_change)) {
   MCFBlock::Index changed = 0;

   MCFBlock::Subset nms(n_change);
   for (MCFBlock::Index i = mMCFB->get_NStaticArcs();
        i < mMCFB->get_NArcs(); ++i) {
    if (mcf->IsDeletedArc(i))
     continue;
    if (!mcf->IsClosedArc(i))
     continue;
    if (drand48() <= 0.5)
     continue;

    nms[changed++] = i;
    mcf->OpenArc(i);

    if (changed >= n_change)
     break;
   }

   if (changed) {
    nms.resize(changed);
    cout << changed << " open";

    if ((mode & 16) && (drand48() < 0.5)) {
     // change via abstract representation
     cout << "(a)";
     for (auto i : nms)
      mMCFB->i2p_x(i)->is_fixed(false);
    } else  // change via call to chg_* method
     mMCFB->open_arcs(std::move(nms));

    cout << " - ";
   }
  }

  // deleting arcs - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  if ((wchg & 32) && (drand48() <= p_change)) {
   MCFBlock::Index changed = 0;

   if (drand48() < 0.5) {
    // delete somewhere in the middle

    for (MCFBlock::Index i = mMCFB->get_NStaticArcs();
         i < mMCFB->get_NArcs(); ++i) {
     if (mcf->IsDeletedArc(i))
      continue;
     if (drand48() <= 0.75)
      continue;

     mcf->DelArc(i);
     mMCFB->remove_arc(i);
     if (++changed >= n_change)
      break;
    }

    if (changed)
     cout << changed << " delete(m) - ";
   } else {
    for (MCFBlock::Index i = mMCFB->get_NArcs();
         --i >= mMCFB->get_NStaticArcs();) {
     if (mcf->IsDeletedArc(i))
      continue;
     if (drand48() <= 0.13)
      break;

     mcf->DelArc(i);
     mMCFB->remove_arc(i);
     if (++changed >= n_change)
      break;
    }

    if (changed)
     cout << changed << " delete(e) - ";
   }
  }

  // creating new arcs - - - - - - - - - - - - - - - - - - - - - - - - - - -

  if ((wchg & 64) && (drand48() <= p_change)) {

   MCFBlock::Index changed = 0;
   MCFBlock::Index afterend = 0;
   while (changed < n_change) {
    if (drand48() <= 0.13)
     break;

    ++changed;

    // random sn != en
    MCFBlock::Index sn, en;
    do {
     sn = drand48()*mMCFB->get_NNodes() + 1;
     en = drand48()*mMCFB->get_NNodes() + 1;
    } while (sn == en);

    // random cost in [ - c_max , c_max ]
    auto cst = c_max*(1 - 2*drand48());

    // random capacity <= 0.75 u_avg
    auto cap = 1.5*(u_avg - u_min)*drand48() + u_min;

    auto arc = mMCFB->add_arc(sn, en, cst, cap);
    if (arc != mcf->AddArc(sn, en, cap, cst))
     diffarcs = false;

    if (arc >= m)
     ++afterend;

    if (mMCFB->get_NArcs() >= mMCFB->get_MaxNArcs())
     break;
   }

   if (changed) {
    cout << "create " << changed << "(" << afterend << ")";
    if (diffarcs)
     cout << "[d]";
    cout << " - ";
   }
  }

  // check that the status of the arcs is the same- - - - - - - - - - - - - -

  if (wchg & 120) {  // ... if it can ever change

   n = mcf->MCFn();
   if (n != mMCFB->get_NNodes()) {
    cerr << endl << "error: different number of nodes" << endl;
    exit(1);
   }

   m = mcf->MCFm();
   if (m != mMCFB->get_NArcs()) {
    cerr << endl << "error: different number of arcs" << endl;
    exit(1);
   }

   if (!diffarcs) {
    for (MCFBlock::Index i = 0; i < m; ++i) {
     if (mcf->IsDeletedArc(i)) {
      if (!mMCFB->is_deleted(i)) {
       std::cerr << "inconsistent del status for arc " << i << std::endl;
       exit(1);
      }
      continue;
     }

     if (mcf->IsClosedArc(i))
      if (!mMCFB->is_closed(i)) {
       std::cerr << "inconsistent cls status for arc " << i << std::endl;
       exit(1);
      }
    }
   }
  }

  // finally, re-solve the problems- - - - - - - - - - - - - - - - - - - - -
  // yet, if the problem is either unfeasible or unbounded, or something has
  // gone awry with the arcs names, re-load it in both MCFClass and MCFBlock

  if ((!SolveMCF()) || diffarcs) {
   load(argv[1]);
   n = mcf->MCFn();
   m = mcf->MCFm();
   diffarcs = false;
  }
 }  // end( main loop )- - - - - - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 // destroy objects and vectors - - - - - - - - - - - - - - - - - - - - - - - 
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 // delete dMCFB;
 delete oMCFB;
 delete mcf;

 // terminate - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 return (0);

}  // end( main )

/*--------------------------------------------------------------------------*/
/*------------------------ End File test.cpp -------------------------------*/
/*--------------------------------------------------------------------------*/
