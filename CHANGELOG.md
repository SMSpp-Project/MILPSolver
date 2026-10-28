# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

- `MILPSolver::previous_column_index()`, `previous_row_index()` and
  `previous_sizes()`: each `load_problem()` but the first finds, by the
  address of its `ColVariable` or `FRowConstraint`, where each column and
  row of the problem it loads was in the one loaded before, so that a
  back-end that loads the problem again whole at a Modification
  (`reload_on_modification()`) can start from the previous solve even when
  columns or rows have been added or removed; `test_qp` checks it after a
  row has been removed

- `dblUpCutOff` and `dblLwCutOff` act on a continuous model too in
  `CPXMILPSolver`, `GRBMILPSolver` and `HiGHSMILPSolver`: the cutoff on the
  bound that the dual simplex certifies (the lower one of a minimization,
  the upper one of a maximization) becomes the objective limit of the dual
  simplex of the back-end, which stops as soon as its dual objective
  passes it and returns `kCutOff`, with a certified bound at least as good
  as the cutoff and a feasible dual solution, i.e., one that gives a valid
  cut, and no bound on the other side, since the solution it stops at is
  not primal feasible (the presolve is turned off for this solve in `GRBMILPSolver` and
  `HiGHSMILPSolver`, whose presolved solve gives no such dual solution or
  ignores the limit, and HiGHS acts on it only when minimizing). The
  parameters changed for it are given back after the solve.
  `test_cutoff_status` checks it on an LP with 60 columns and 40 rows, in
  either sense, with a cutoff that is reached and with one that is not

- `MILPSolver::get_elapsed_iterations()` gives the iterations of the
  back-end in the last `compute()`, all its solves summed: those of DAQP in
  `DAQPMILPSolver` and the working set recalculations of qpOASES in
  `QPOASESMILPSolver`, those of a hot start repeated from scratch included;
  0 for the back-ends that do not count them

- `get_var_solution()`, `get_dual_solution()`, `get_var_direction()` and
  `get_dual_direction()` of every :MILPSolver write only a part of the
  solution when they are given a `SimpleConfiguration< std::vector< int > >`:
  the sub-Blocks of the `Block` it lists, each with all its sub-Blocks, and
  the `Block` itself without them for -1; the columns and rows not asked for
  are skipped a group at a time, and any other `Configuration` gives the
  whole solution as before. `test_partial_solution` checks it on every
  back-end, dynamic rows and columns included

- the makefiles build the C++20 module `smspp.milpsolver` when `SMSppMOD`
  is not empty (e.g., `make SMSppMOD=1`), with clang++ 16 or later [see
  `lib/makefile-mod` of SMS++]

- `GRBMILPSolver`, `CPXMILPSolver` and `SCIPMILPSolver` free after each
  solve of an LP ending optimal the memory that the back-end keeps and the
  next solve does not need: the solution and the duals are kept by the
  `Solver`, and the problem is replaced by a copy of itself (GUROBI,
  CPLEX), whose next solve starts from the kept basis, or the transformed
  problem of SCIP is freed (SCIP starts from scratch after any change in
  any case). On 16 copies of a subproblem of 129676 rows and 124492
  columns of a Benders decomposition of a PyPSA-Eur instance the memory
  after the solves goes from 2.0 to 0.4 GB with GUROBI, from 1.7 to 0.4 GB
  with CPLEX and from 11.1 to 5.9 GB with SCIP, with the same simplex
  iterations, and with the 200 subproblems of the instance with 200
  scenarios the peak of `GRBMILPSolver` goes from 50 to 30 GB, with the
  same rounds. `MOSEKMILPSolver` keeps its task, since a copy of it does
  not lower the memory (what the deleted task had stays with the
  environment of MOSEK), and `HiGHSMILPSolver` its instance: a new one with
  the same options, model and basis frees 2.0 of the 2.3 GB, but its next
  solve does not start from the state the dual simplex had, and on the
  badly scaled master of a Benders decomposition HiGHS then ends in
  `Unknown` and, from scratch, in a wrong `Unbounded`, where the instance
  that is kept solves it.
  `test_free_workspace` checks every `:MILPSolver` of `LPPar-all.txt`
  against the others across warm re-solves after changes of sides, bounds
  and costs, an infeasible solve and a solve with no change

- `HPIPMMILPSolver`, the interface with the dense interior-point solver
  of HPIPM (0.1.2 or later, on BLASFEO) for convex QPs and LPs with linear
  constraints: the rows with equal sides and the fixed variables are its
  equalities, whose dependent ones it removes, the other bounds and rows
  its box and general constraints, an infinite side masked out; the duals
  follow the convention of the other backends, `intHPIPMMode` picks the
  mode of HPIPM (ROBUST by default) and `intHPIPMWarmStart` starts a solve
  from the previous solution; it is built when `MILPSolver_USE_HPIPM` is
  ON and `FindHPIPM` finds it (CMake) or HPIPM is in `HPIPM_ROOT`
  (makefile), and the tests of the QP backends run it

- `dblKKTTol`, off by default: when positive, `compute()` writes the
  solution in the `Block` after a solve that returns `kOK` and measures
  its primal residual, its dual residual and its gap, with the duals in the
  convention of `write_dual_solution()` and the multipliers of the bounds
  taken as the reduced costs that the duals of the rows imply, and returns
  `kError` if one of them exceeds the value; `get_kkt_residuals()` gives
  them, and `test_qp` checks every back-end and a sloppy `OSQPMILPSolver`

- `ClarabelMILPSolver` and `OSQPMILPSolver` give the unbounded ray of a
  dual infeasible problem through `get_var_direction()` and the Farkas
  certificate of a primal infeasible one through `get_dual_direction()`,
  with the sign of those of the other :MILPSolver and with
  `intHomogeneousDirection`, and `ProxQPMILPSolver` gives the ray;
  `test_farkas` checks them, and checks on an unbounded convex QP with a
  singular Hessian that the ray of every :MILPSolver that gives one is a
  ray

- `MOSEKMILPSolver`, the interface with the commercial solver MOSEK (10 or
  later) for LPs, MILPs, convex QPs and MIQPs with linear constraints:
  bounds, sides, costs, coefficients, integrality, rows, columns and terms
  of the quadratic objective are changed in place in the task, which keeps
  its basis for the simplex; the duals and reduced costs follow the
  convention of the other backends, and the certificates of infeasibility
  and unboundedness come with the solve that finds them; `dblUpCutOff` and
  `dblLwCutOff` stop a MIP from the callback of MOSEK, quadratic
  constraints and the separation of user cuts or lazy constraints in a MIP
  are refused, and every parameter of MOSEK can be given with its symbolic
  name, listed for the installed version by `tools/mosek_pars`; it is
  built when `MILPSolver_USE_MOSEK` is ON (CMake) or MOSEK is in
  `MOSEK_ROOT` (makefile)

- `QPOASESMILPSolver`, the interface with the open-source online active-set
  solver qpOASES for convex QPs (and LPs) with linear constraints, dense and
  with no integer variables, whose Hessian may be singular: a direction of
  zero curvature is followed until a constraint blocks it, and an unbounded
  problem is declared so. The problem is loaded again whole at every
  Modification, and a solve is a hot start from the solution of the
  previous one if the size has not changed (`intQPOASESHotStart`), repeated
  from scratch if it does not end optimal; it is built when
  `MILPSolver_USE_qpOASES` is ON (CMake) or qpOASES is in `qpOASES_ROOT`
  (makefile)

- `PIPSMILPSolver` solves convex QPs: the Hessian of the objective is given
  to PIPS-IPM++ node by node through the callbacks that had it empty, a
  term between the Variable of two nodes making `load_problem()` throw, and
  a solve with a Hessian runs with the presolve of PIPS-IPM++ off and its
  outer BiCGStab given no iteration, since the presolve treats the
  objective as linear and the BiCGStab refines on Q + diag( Q ), which
  halves the residual at each step and ends in a crash of the
  factorization; a Modification of the quadratic part loads the model
  again, as any other one does

- `test_qp` checks every `:MILPSolver` of the build on small quadratic
  models with a known optimum: the value, the solution and the duals (with
  the stationarity of the KKT conditions) of a convex QP with a separable
  and a non-separable objective, minimized and maximized, the value of an
  MIQP, and the value after each of a sequence of Modifications against a
  Solver loaded anew, a Solver that is not meant to handle a model having
  to refuse it rather than answer

- the Objective of the Block given to `Solver::set_excluded_objectives()` is
  not loaded, nor is its sense, and its Modification are dropped, a
  GroupModification having one of them inside being processed piece by piece
  rather than as one operation of the back-end

- a dynamic `FRowConstraint` whose Function is not a `LinearFunction` (e.g.,
  a `DQuadFunction`) added to a loaded problem has the problem loaded again,
  which reads it as the rows of the problem as first loaded are, rather than
  have the back-end throw, which adds only linear rows; so do the dynamic
  rows and columns that come or go in a problem with quadratic rows, whose
  duals `GRBMILPSolver` read at the wrong index once a linear row had been
  added in place

- `DAQPMILPSolver`, the interface with the open-source dual active-set
  solver DAQP for convex QPs (and LPs) with linear constraints, dense and
  with no integer variables: the problem is loaded again whole at every
  Modification, and a solve starts from the working set of the previous one
  if the size has not changed (`intDAQPWarmStart`), but for the constraints
  that either problem fixes, which are active only if the new one fixes
  them; it is built when
  `MILPSolver_USE_DAQP` is ON (CMake) or DAQP is in `DAQP_ROOT` (makefile)

- `ClarabelMILPSolver`, the interface with the open-source interior-point
  solver Clarabel for convex QPs (and LPs) with linear constraints, sparse
  and with no integer variables: an unbounded problem is detected as such
  (`kUnbounded`) and a singular Hessian needs no regularization, as in the
  master problems of bundle methods with vertical cuts; Clarabel has no
  warm start, hence the problem is loaded again whole at every Modification
  and each solve starts from scratch; it is built when
  `MILPSolver_USE_Clarabel` is ON (CMake) or Clarabel is in `Clarabel_ROOT`
  (makefile)

- `ProxQPMILPSolver`, the interface with the open-source proximal augmented
  Lagrangian solver ProxQP of ProxSuite for convex QPs (and LPs) with linear
  constraints and no integer variables, with its dense backend (the bounds
  as box constraints) or its sparse one (`intProxQPBackend`): the problem is
  loaded again whole at every Modification, and a solve starts from the
  primal and dual solution of the previous one if the size has not changed
  (`intProxQPWarmStart`); it is built when `MILPSolver_USE_ProxQP` is ON and
  ProxSuite is found (CMake, `proxsuite_ROOT` as a hint) or ProxSuite is in
  `proxsuite_ROOT` (makefile)

- `PIQPMILPSolver`, the interface with the open-source proximal interior
  point solver PIQP for convex QPs (and LPs) with linear constraints and no
  integer variables, with its dense backend or its sparse one
  (`intPIQPBackend`) and the bounds as box constraints of PIQP: the problem
  is loaded again whole at every Modification, and each solve starts from
  scratch, PIQP having no warm start; it is built when
  `MILPSolver_USE_PIQP` is ON and PIQP 0.6 or later is found (CMake,
  `piqp_ROOT` as a hint) or PIQP is in `piqp_ROOT` (makefile)

- `OSQPMILPSolver`, the interface with the open-source operator splitting
  solver OSQP for convex QPs (and LPs) with linear constraints and no
  integer variables, with the bounds as unit rows of the constraints of OSQP
  (the solution is projected on them), the polishing of the solution on and
  the tolerances of the certificates of infeasibility at 1e-8 by default:
  the problem is loaded again whole at every Modification, but the data are
  given to the solver already set up if their sparsity pattern has not
  changed, and each solve starts from the primal and dual solution of the
  previous one (`intOSQPWarmStart`); it is built when `MILPSolver_USE_OSQP`
  is ON and OSQP 1.0 or later is found (CMake, `osqp_ROOT` as a hint) or
  OSQP is in `osqp_ROOT` (makefile)

- `QPMADMILPSolver`, the interface with the open-source dual active-set
  solver qpmad, of Goldfarb and Idnani, for strictly convex QPs with linear
  constraints, dense and with no integer variables: H has to be positive
  definite, and a solve whose Cholesky factorization of H fails ends with
  `kError` unless H is regularized (`dblQPMADReg`); the problem is loaded
  again whole at every Modification and each solve starts from the
  unconstrained minimum, qpmad having no warm start of the active set, but
  the Cholesky factor of H is kept from a solve to the next as long as H
  has not changed (`intQPMADReuseFactor`); the multipliers, of which qpmad
  only gives those of the active inequalities, are computed from its active
  set; it is built when `MILPSolver_USE_qpmad` is ON and qpmad is found
  (CMake, `qpmad_ROOT` as a hint) or qpmad is in `qpmad_ROOT` (makefile)

- `QuadProgMILPSolver`, the interface with the open-source dual active-set
  solver quadprog, the original implementation of the method of Goldfarb
  and Idnani, for the same QPs as `QPMADMILPSolver`: each finite bound and
  side is a constraint of quadprog and an equality is one, the multipliers
  are those of quadprog, the inverse of the Cholesky factor of H, regularized
  if so asked (`dblQuadProgReg`), is computed once and kept as long as H has
  not changed (`intQuadProgReuseFactor`), and quadprog has no limit on the
  iterations or the time; quadprog is under the GPL, and it is built when
  `MILPSolver_USE_quadprog` is ON and `libquadprog` is found (CMake,
  `quadprog_ROOT` as a hint) or is in `quadprog_ROOT/lib` (makefile)

- the C++20 module `smspp.milpsolver` (`src/MILPSolver.cppm`), which exports
  every name that the headers declare at namespace scope and re-exports
  `smspp`, so that a source file can write `import smspp.milpsolver;` in
  place of the `#include` of the headers; it is built with the CMake option
  `SMSpp_BUILD_MODULE`, as the target `SMS++::MILPSolver_module`. The names
  of a backend are there only if the backend is

### Changed

- `DAQPMILPSolver`, `QPOASESMILPSolver`, `ProxQPMILPSolver`,
  `OSQPMILPSolver` and `HPIPMMILPSolver` start a solve from the previous one
  also after columns or rows have been added or removed, on those that are
  still there [see `MILPSolver::previous_column_index()`], instead of only
  when the size of the problem has not changed: DAQP from the working set,
  ProxQP from the primal solution and the multipliers of the equalities,
  those of the inequalities at 0, OSQP from the primal and dual solution,
  HPIPM from the primal solution, repeating from its initial point a solve
  that does not end, and qpOASES with a hot start if no column or row has
  moved, otherwise with a new problem initialized from the working set,
  solution and multipliers, given up for a solve from scratch after half
  as many working set recalculations as variables and rows (at least 10)

- `MILPSolver::get_elapsed_iterations()` also gives the iterations of
  ProxQP, OSQP and HPIPM

- `GRBMILPSolver` starts the environment of GUROBI with `OutputFlag` off,
  which it turns on again right after: the license parameters and banner
  are no longer printed, once per process, and the log of the tests on the
  CI, which start thousands of processes, stays within the limit of
  GitLab

- the pipeline builds the module also with its makefiles, after CMake, and
  runs on a change of the makefiles too [see `ci/makefile` of the umbrella]

- `SMSpp_BUILD_MODULE` no longer stops the configuration when the C++20
  module cannot be built, i.e., with a CMake older than 3.28, a generator
  other than Ninja or Visual Studio, a compiler other than Clang 16 (with
  `clang-scan-deps`) or MSVC 19.34, or a module it imports that is not
  there: the option is turned off with a warning that says why, and the
  build goes on with the headers [see `cmake/CheckModule.cmake` of SMS++]

- `write_dual_solution()` finds the bounds of each column in constant time,
  from a bound per column kept up to date when Variable and bounds are
  added or removed, with the few columns having more than one in a map,
  instead of looking for them with a `dynamic_cast` among the active stuff
  of every Variable at each call, which on an `InvestmentBlock` over a
  `TwoStageStochasticBlock` solved as a MILP took more than half of the
  time; the duals written are the same

- `MILPSolver::clear_problem()` gives back the memory of the arrays the
  model is loaded from, which `clear()` kept, so that a `Solver` no longer
  holds a second copy of its model next to the one of the back-end

- the minimum versions of the QP back-ends are those their code builds and
  passes the tests with, checked by CMake: MOSEK 11.0.1 (was 10.0, which
  lacks some error codes), ProxSuite 0.4 (was 0.6), qpmad 1.5, DAQP 0.8.2
  and Clarabel.cpp 0.11.1, the last two from their headers, which carry no
  version number, as does qpOASES when it has no CMake configuration;
  the configuration of the installed MILPSolver asks for the same
  versions

- the CI runs the batteries of the QPs of qpbenchmark (`batchMPC_C`,
  `batchMM_C` and their dense versions) in a stage `test-qpb` of their
  own, after `test-miqp`, and the job `test` leaves them out

- `load_problem()` scans the columns with buffers kept for the whole load,
  rather than with a std::set, a std::map and two vectors of its own for
  each column, the entries of a column being the same in the same order,
  and `GRBMILPSolver` gathers the nonzero diagonal of a quadratic objective
  in one pass rather than erasing the zeros one at a time

- `HiGHSMILPSolver` leaves `qp_regularization_value` at the 1e-7 of HiGHS,
  the 1e-5 of 0.10.0 making many QPs less precise and slower (the Lagrangian
  dual of a UC master and of an SDDP stage); a QP that the active set QP
  solver of HiGHS does not solve is solved again from scratch and then with
  the Hessian regularized by 1e-5 and by 1e-4, and if none of these solves
  it the Solver throws, saying which Solver or stabilization to use instead,
  rather than returning a solution that is not there

### Fixed

- `HiGHSMILPSolver` returns `kError` on a continuous QP that the active set
  of HiGHS declares optimal with a primal-dual objective error above
  `QPObjectiveError` (1e-3): on the unbounded master problems of a bundle
  with vertical linearizations, whose Hessian is singular, HiGHS 1.15.1
  ends at a finite point with the status optimal and an error close to 1,
  and so it did on 16 of 190 QPs (the master problems of `PFB`, MMCF and UC
  and the QPs of `qpbenchmark`), 11 of them unbounded, all with an error
  of at least 3.6e-3, while 8 of the 174 it solves right have an error
  above 1e-3, which are now `kError` as well

- `ProxQPMILPSolver` and `HPIPMMILPSolver` start from the previous solve
  also after a Modification: the `clear_problem()` that
  `MILPSolver::load_problem()` calls threw away the solution of ProxQP and
  the memory of HPIPM, hence every solve after a Modification started from
  the initial guess

- `SCIPMILPSolver` with `intComputeDuals` gave no duals (all of them
  `SCIP_INVALID`, i.e., 1e+99) after a change of the costs that left the
  previous solution optimal: SCIP keeps the solutions of a solve when the
  problem changes and gives them to the next one as candidates, one of
  them became the best solution in place of that of the LP, and the duals
  belong to the latter. `intComputeDuals` now sets `limits/maxorigsol` to
  0, so that no solution is carried over

- `MOSEKMILPSolver` gives the diagonal of the Hessian the value of the
  quadratic coefficient of the objective, not the sum of the deltas of
  its changes: a coefficient that went back to its value after a few
  changes kept a rounding error, which on a singular Hessian (e.g., the
  Wolfe dual of a SVM with a polynomial kernel) made MOSEK refuse the
  objective as not convex

- `FindHiGHS.cmake` links the extras library of HiGHS (`libhighs_extras.a`)
  after a static `libhighs.a`, which calls its functions from HiGHS 1.15
  on: a static build of HiGHS, as the one of vcpkg, left the tools and the
  consumers of MILPSolver with undefined references to
  `HighsExtras_getApi`

- the duals of a model with quadratic constraints in `GRBMILPSolver` and
  `CPXMILPSolver`, and the index of a linear row among them, no longer
  read `q_part`, which `load_problem()` empties: they tell the quadratic
  rows by the indices the back-end keeps [see `is_quad_row()`]. Since
  `clear_problem()` gives the memory back, the reading crashed, e.g., in
  `linear_constraints_test` of BundleSolver

- `write_dual_solution()` writes in the quadratic rows the dual that the
  back-end gives them, which it skipped, leaving there whatever was there
  before: `CPXMILPSolver` and `GRBMILPSolver` computed it, and it was lost;
  `GRBMILPSolver` gives a quadratic row with a linear part the dual of its
  linear row, as `CPXMILPSolver` does, and no longer that of the quadratic
  constraint Gurobi gets it with, which has the opposite sign

- the dual that `CPXMILPSolver` writes in a quadratic row was always -1, and
  the duals of the linear rows that follow a quadratic row with no linear
  part were shifted by one: a quadratic row with a linear part takes the
  dual of its linear row, and one with no linear part the multiplier that
  its dual slack gives; `test_qcp_dual` checks them on every back-end that
  takes quadratic rows

- the makefiles build `CPXMILPSolver`, `GRBMILPSolver`, `SCIPMILPSolver`
  and `HiGHSMILPSolver` only if CPLEX, Gurobi, SCIP and HiGHS are there,
  as `CPLEX_PRESENT`, `GUROBI_PRESENT`, `SCIP_PRESENT` and `HiGHS_PRESENT`
  of `extlib` say, as for the optional back-ends, and add their -I and
  libraries only then, and `tools/makefile` builds the tool of a solver
  only if the solver is there: a build without CPLEX (e.g., the image of
  the pipelines) failed on `ilcplex/cplex.h`

- the makefiles of the tests and of the tools build also with g++: they
  take the option that stops at the first error and the one of the
  debugging information, which clang++ and g++ spell differently
  (`-ferror-limit=1` and `-glldb` against `-fmax-errors=1` and `-ggdb`),
  from `SMSppSWERR` and `SMSppSWDBG` of `lib/makefile-OS` of SMS++

- the makefile finds `tools/` also when `MILPSSDR` is an absolute path:
  `./` was put in front of it, which turned `/home/...` into `./home/...`

- `test_farkas` takes the values of `MSK_IPAR_OPTIMIZER` from the `mosek.h`
  of the build and `test/LPPar-ipm.txt` gives the interior point of MOSEK 12:
  MOSEK 12 numbers the optimizers differently, and with it the test asked
  for the deprecated interior point and got no ray on the unbounded QP

- an Objective whose Function is a `PolyhedralFunction` makes the
  exception say that the `PolyhedralFunctionBlock` has to be in its
  linearized representation, rather than that the type of the Function is
  unknown
- `variable_with_index()` and `constraint_with_index()` of a dynamic column
  or row check the index against the dynamic ones, rather than reading past
  them
- `SCIPMILPSolver` gives no unbounded ray when the problem has a quadratic
  objective or quadratic rows: the ray of SCIP is then that of its LP
  relaxation, where they are outer approximated, and on an unbounded QP it
  was not a ray of the problem

- `QuadProgMILPSolver` projects the solution on the bounds of the
  variables, which quadprog has as constraints and satisfies only up to the
  accuracy of the dual method on the regularized Hessian: a fixed variable
  came back up to 4e-3 away from its value, and `ColVariable::set_value()`
  refused it (10 of the 60 cases of the PolyhedralFunctionBlock battery
  with QuadProgMILPSolver as the solver of the master of the bundle)

- the HiGHS that `HiGHSMILPSolver` asks for is 1.8.0 or later: it calls the
  callbacks of the C API of HiGHS, which came with 1.6.0, and
  `Highs_setSparseSolution`, which came with 1.8.0, so that with a HiGHS
  between 1.5.3 and 1.7.2 the configure accepted it and the build failed

- `SCIPMILPSolver` tells whether a quadratic objective is convex from the
  pivots of a sparse LDL^T factorization of Q + tol I, rather than from the
  eigenvalues of a dense Q, which on the QPLIB instances of the MIQP
  battery took longer than the solve itself and outside its time limit

- `SCIPMILPSolver` solved a QP whose objective is not separable by spatial
  branch-and-bound when SCIP is built without LAPACK (through Ipopt), as
  it then cannot prove the auxiliary constraint z >= x^T Q x convex; with
  free variables the search need not end, and it ran forever on the
  re-solve after a Modification and on any solve with `intComputeDuals`.
  SCIP is now told that the problem is convex
  (`constraints/nonlinear/assumeconvex`) when the only nonlinear
  constraint is that of the objective and its matrix is positive
  semidefinite (negative in a maximization), as checked before each solve
  that follows a change of the objective; the first quadratic term given
  to an objective that had none built that constraint with null variables

- `PIPSMILPSolver` threw on a non-separable quadratic objective but solved
  a separable one (a `DQuadFunction`) as if its quadratic terms were not
  there, returning the optimum of another problem

- a MIP solve that `CPXMILPSolver`, `GRBMILPSolver`, `HiGHSMILPSolver` or
  `SCIPMILPSolver` interrupted because a cutoff (`dblUpCutOff`,
  `dblLwCutOff`) was reached returned `kError`, with infinite bounds: it
  now returns `kCutOff`, which `Solver.h` defines for this case, with the
  value of the incumbent and the bound of the search from `get_lb()` and
  `get_ub()` as after `kStopTime`, and the incumbent as the solution
  (`SCIPMILPSolver` gives it with `kCutOff` as well), and `GRBMILPSolver`
  no longer warns of an interrupted run for it. A `LagBFunction` or
  `BendersBFunction` passes the status of its inner Solver on, and a
  `BundleSolver` that gives them targets (`intTrgtMng`) stopped with
  `kError` at the first component that reached one, while it reads
  `kCutOff` as an evaluation stopped early; `test_cutoff_status` expects
  `kCutOff`, a solution beyond the cutoff and a valid bound

- the cutoffs `dblUpCutOff` and `dblLwCutOff` of a `ComputeConfig` never
  reached a solve: `CPXMILPSolver`, `GRBMILPSolver`, `HiGHSMILPSolver` and
  `SCIPMILPSolver` set them back to infinite whenever they loaded the
  problem, in `set_Block()` (which a `BlockSolverConfig` calls after giving
  the `ComputeConfig`) and at any reload, and now keep them as any other
  parameter. A lower cutoff alone did not install the callback that checks
  them, which compared it with + infinity, and the callbacks compared the
  cutoffs with the values of the solver, which lack the constant of the
  objective. `GRBMILPSolver` also checks the cutoffs at each new
  incumbent, unless it separates lazy constraints, since a short run may
  find its solutions with no call of the callback in between

- the callback of `HiGHSMILPSolver` did not interrupt a minimization that
  reached a cutoff, and took the incumbent for the bound, so that the
  lower cutoff of a maximization stopped it at an incumbent below the
  cutoff with nothing proved; it now reads the bound of HiGHS and
  interrupts a minimization as a maximization

- `SCIPMILPSolver` did not stop at a cutoff: the constraint handler it
  installed for one did nothing with it, and added to a solved problem, as
  a cutoff given after a solve made it, returned `SCIP_INVALIDCALL` (-8) as
  the status of `compute()`. An event handler now interrupts SCIP at a new
  incumbent or at a solved node that reaches a cutoff, which ends the solve
  with `kCutOff` as for the other `:MILPSolver`; the constraint handler of
  the separation is added once, to the original problem with the
  transformed one gone, and again after the problem is loaded again.
  `test_cutoff_status` attaches each `:MILPSolver`, `SCIPMILPSolver`
  included, with an upper cutoff, a lower one and both in its
  `ComputeConfig`, to a maximization and a minimization MIP with a constant
  in the objective, and checks that the first solve stops at the first
  solution if that reaches a cutoff and gives the optimal value otherwise

- `HiGHSMILPSolver` wrote (`write_lp()`, `strOutputFile`) the model with
  the slack columns of the free rows that an interior point solve had added
  to it: it now writes a copy without them, the rows free again, which
  leaves the model and the solution of HiGHS as they are;
  `test_ipm_free_rows` checks that each `:MILPSolver` writes the same model
  before and after a solve

- `HiGHSMILPSolver` threw "unknown status" on a run that HiGHS ended with
  the model status Interrupt, which its callback asks for when a cutoff is
  reached on a maximization MIP, or with Postsolve error: both are now
  `kError`, as an aborted or interrupted run of `CPXMILPSolver` and
  `GRBMILPSolver` and a presolve or solve error of HiGHS are;
  `test_cutoff_status` solves a maximization MIP with an upper cutoff that
  its first integer solution reaches with each `:MILPSolver`

- the interior point of HiGHS up to 1.15.1, asked for on an LP without
  presolve, solved another LP whenever a row with no finite side came
  before a row with two finite sides: the LP that HiGHS builds for IPX
  drops the free rows and renumbers the others, but puts the slack of a
  row with two finite sides in the row with its original index
  (ERGO-Code/HiGHS#3319), so that a feasible LP could be declared
  infeasible with no error to tell. `HiGHSMILPSolver` now gives each free
  row, before the interior point solves the LP (also as the second attempt
  after a failed simplex), a free slack column s, the row becoming
  a'x + s = 0: the rows keep their index and the free ones their dual value
  of 0, the slacks are out of every solution read and go at the next
  Modification, which loads the problem again. The remedy is compiled only
  with a HiGHS up to 1.15.1; `test_ipm_free_rows` checks the interior point
  of each `:MILPSolver` on such an LP and on its changes

- an LP on which HiGHS ended with model status Solve error got no second
  attempt, and the second attempt after an Unknown was always the interior
  point, even when the interior point was the algorithm that had failed:
  `HiGHSMILPSolver` now solves again with the other algorithm after either
  status, the simplex if the `solver` option asks for an interior point and
  the interior point otherwise. The interior point of HiGHS (up to 1.15.1
  at least) fails at random with a Solve error on an LP with a free row
  that is not the last one, reading past the end of its list of free rows
  (ERGO-Code/HiGHS#3319), which stopped the Solver of the inner Block of
  some `InvestmentBlock` instances with `kError`

- `GRBMILPSolver` gave a row with no finite side, and a relaxed one, to
  GUROBI as a'x <= GRB_INFINITY, a row with right-hand side 1e100 on which
  a warm re-solve after a change of the coefficients of other rows ends in
  numerical trouble (status 12, from a right-hand side of 1e33 on): such a
  row is now the equality a'x + s = 0 with a free slack s, which keeps the
  index of the row, as the ranged rows keep theirs, and the slack gets its
  bounds back when a side becomes finite or the row is enforced; a row that
  is still relaxed stays free when a side of it changes

- the bound of `CPXMILPSolver` on the side of the optimum it has not proved,
  on a continuous problem stopped by a limit, was the value of whatever
  point CPLEX had (0 with no point at all), which is no bound: it is now
  that value only if the point is a dual feasible basis, and -INF / +INF
  otherwise; on a problem with integer variables a best bound that CPLEX
  does not have is now infinite as well. The same for `SCIPMILPSolver`,
  whose bounds were SCIP's infinity (1e20) where SCIP has none

- a double parameter set out of the range that the back-end admits (an
  infinite or a negative `dblMaxTime`, say) was refused by CPLEX, GUROBI
  and HiGHS, which kept the previous value, so that a Solver solved once
  with `dblMaxTime` 0 did not reach the optimum after the time limit was
  lifted with `dblMaxTime` = INF; SCIP refused it with an error. It is now
  taken as the nearest end of the range

- the bound of `GRBMILPSolver` on the side of the optimum it has not proved
  (`get_lb()` when minimizing, `get_ub()` when maximizing) on a continuous
  problem whose solve does not end optimal, e.g., stopped by the time limit
  or with a sub-optimal solution, unless `NonConvex` is 0: it was GUROBI's
  infinity (1e100) or its "undefined" (1e101), finite numbers that the
  callers took for a bound. A `LagBFunction` whose inner Block is such a
  problem built from it a linearization with a constant of the order of
  1e100, which made the master problem of the `BundleSolver` above it fail
  at every call thereafter. The bound is now infinite

- the rows that a Block separates in the callback of a solve became rows of
  the model at the following `compute()`, all of them and for good, whether
  the back-end had used them or not: on a unit commitment re-solved after
  changes of the costs the perspective cuts of the first solve took the
  model from some 1500 rows to more than 70000, half of them copies of one
  another and most nearly parallel, and with them GUROBI declared optimal
  a solution with a value worse than the optimum by 18 out of 34800. Those
  rows are now given to the back-end as user cuts or lazy constraints of
  the solve they are separated in, and are kept in a pool, out of the model,
  that the callback gives back whenever they are violated (the most
  violated one on each set of variables), so that the model keeps the rows
  of the formulation alone; the pool becomes rows of the model, as before,
  only before a solve that does not separate lazy constraints (with HiGHS,
  whose callback cannot give them, always), and a Block that removes one of
  its rows from the pool is followed [see `separate_in_callback()`]

- the coefficients of the quadratic part of a constraint were stored as
  `float` between the scan of the Block and the loading into the backend,
  so that each of them reached the solver rounded to some 7 significant
  digits (those of the envelope of the cosine of `ACNetworkBlock`, for
  instance); they are now kept as `double`, as all the other coefficients
  of the model

## [0.10.0] - 2026-10-09

### Added

- the tests of this directory carry the label of the module, so that the
  pipeline, which selects with `ctest -L <module>`, runs them: they were built
  and never run. The configurations they read name the four `:MILPSolver`
  backends, each with the `ComputeConfig` that goes with it, and
  `keep_available_Solvers()` of the new `test_common.h` keeps the first ones
  this build has before the configuration is applied, so that the same file
  works wherever it is run; a test needing more backends than there are gives
  up with the 77 that `ctest` reads as "skipped" rather than failing, which is
  what `test_dynamic`, which needs two of them, does on a machine with one

- `test/test_farkas.cpp` and `test/test_groups.cpp`, i.e., the two testers
  that used to live in the suite of this module in `tests/`: the first asks
  for the certificate of an infeasible and of an unbounded model and checks
  its sign and its scale, the second builds a model out of the groups of a
  Block and changes it by groups, and neither of the two needs a Block of
  another module to exist, so they belong here

- `process_group_modification()`, which is asked whether a whole
  GroupModification is one single operation of the back-end before the group
  is taken apart: the shape it recognises is the column, declared by a
  `VariableGroupMod`, and it is handed to `add_columns()` / `remove_columns()`,
  virtual and answering false unless a back-end implements them, so that a
  back-end that does not keeps exactly the behaviour it had. GRBMILPSolver
  implements both: a batch of columns is one `GRBaddvar` each plus one
  `GRBchgcoeffs` and one `GRBsetdblattrlist` for all of them, and removing a
  column is one `GRBdelvars` with none of the coefficient changes, which the
  deletion does by itself

- `for_each_row_addition()`, which walks a Modification and the groups inside
  it looking for additions of FRowConstraint

- the shapes that a group need not declare, recognised from what it holds:
  changes of the linear coefficients, of the bounds of the columns and of the
  sides of the rows, each batched on its own and in the order they were
  issued, with `change_coefficients()`, `change_bounds()` and
  `change_sides()` to execute them (virtual, answering false unless a
  back-end implements them, GRBMILPSolver implementing all three). Groups
  nest, so the leaves of a group inside a group belong to the same batch,
  while a structural change, or a group declaring a shape of its own, closes
  the batch, is executed on its own and opens the next one

- the same three shapes executed by CPXMILPSolver, with one `CPXchgcoeflist`
  and one `CPXchgobj` for the coefficients, one `CPXchgbds` for the bounds
  and one `CPXchgrhs` plus one `CPXchgsense` for the sides, and by
  HiGHSMILPSolver, with the `*BySet` calls, whose sets have to be ordered and
  without repetitions; `Highs_changeCoeff` has no batched form, so on HiGHS a
  group of coefficients is executed one change at a time. SCIPMILPSolver
  implements none of the three and keeps the behaviour it had

- the column, i.e., `add_columns()` and `remove_columns()`, also on
  CPXMILPSolver and HiGHSMILPSolver: what a group of the column shape holds
  is read once and for all by `read_column_group()` in the base class, each
  back-end writing it with the calls it has

- the shape of a batch of fixings, i.e., of VariableMod, executed by
  `change_variables()`: fixing a column is writing its two bounds, hence a
  whole set of them is two calls on GUROBI, one on CPLEX and one on HiGHS.
  A batch of fixings and one of changes of the bounds are executed in the
  order they were issued, both being about the same attribute of the model,
  and a batch carrying a change of integrality is refused

### Changed

- a batch of changes of the sides that holds a row GUROBI keeps as a ranged
  one is executed whole rather than refused: such a row is an equality whose
  range lives in the upper bound of the auxiliary column standing for its
  slack, and those are three lists like any other. Only the row that becomes
  ranged for the first time goes on its own, a column having to be created
  for it

- a Modification that the Solver does not execute no longer closes the batch
  of a group: a PolyhedralFunctionMod is one, the PolyhedralFunctionBlock
  answering it with the equivalent changes of the abstract representation,
  and closing the batch on it had the cascade of a whole bundle of cuts come
  out one cut at a time

- the coefficients of the Objective a change has to be read back from are
  read in one call rather than one column at a time, which also brings the
  model up to date once instead of once per column

- GRBMILPSolver maps GRB_SUBOPTIMAL to kLowPrecision rather than to kOK,
  since a solution that does not satisfy the optimality tolerances carries
  no accuracy promise; the callers that read the status of a component,
  BundleSolver among them, take it as inexact information

- the Gurobi environment is one for the whole process, created by the first
  GRBMILPSolver and released by the last one, since every environment is a
  session of the license and one per Solver does not scale: a decomposition
  with one component per Solver has the license service refuse the sessions
  it asks for. What a Solver keeps to itself, i.e., the parameters, lives in
  the environment of its own model, which is where `set_par()` writes them
  and where they are read back from

- the errors that have nothing to do with the model being solved, i.e., the
  license service unreachable or refusing the request, are returned as
  kError by GRBMILPSolver instead of being thrown, so that whoever asked
  decides what to do with a computation that may have been running for hours

- `grb_pars` probes by name the parameters that the enumeration of Gurobi
  does not return, so that the table it writes carries them as well

- the makefile asks for `-O3 -DNDEBUG` and nothing else, the macro of the
  patch for `boost::any` on macOS having no reason to be there since there is
  no `boost::any` left in the core

- whoever links the module keeps it: the classes of a module register
  themselves in the factory from a static initialiser, and a linker that
  drops what looks unused takes the registration away with it, so the target
  now tells whoever links it to keep the symbol that forces the module in,
  and on ELF, where naming the symbol is not enough, the library as a whole

- GRBMILPSolver says out loud which status stopped it when that status is a
  failure of the run and not an answer about the model: an interrupt,
  numerical difficulties it cannot recover from, an asynchronous run that is
  not over and the work and memory limits all reached the caller as `kError`
  and nothing else, so that a batch showed a failure indistinguishable from a
  wrong answer, and from the refusal of the license below

- GRBMILPSolver says out loud when the license service refuses the request:
  the errors that have nothing to do with the model that was being solved
  (`GRB_ERROR_NO_LICENSE`, `GRB_ERROR_NETWORK`, `GRB_ERROR_JOB_REJECTED`,
  `GRB_ERROR_CSWORKER`, `GRB_ERROR_CLOUD` and `GRB_ERROR_SECURITY`) were
  mapped to `kError` in silence, so that a batch running into one of them,
  which is what happens when too many sessions of the same license are open
  at the same moment, showed a failure indistinguishable from a wrong
  answer; the status is unchanged, a line on `cerr` now names the reason

- GRBMILPSolver asks again when the license service refuses an optimization,
  as it already did for the start of the environment: the refusal comes in
  the middle of a run, when the token of the session is renewed while the
  license has as many sessions open on other machines as it allows, and a
  single one made the Solver fail, and with it, e.g., a nested Lagrangian
  dual 30 minutes into its run. `GRBoptimize()` is now attempted up to 16
  times, the wait doubling up to 30 seconds (about 6 minutes in all), with a
  line on `cerr` at the first refusal; after the last one the status is the
  `kError` of before

### Fixed

- a batch of fixings that names the same column more than once, e.g., a
  ColVariable fixed and unfixed on the same channel, gave GRBMILPSolver and
  CPXMILPSolver one entry per change in the same `GRBsetdblattrlist` /
  `CPXchgbds`, neither of which says which entry wins, and on Gurobi the
  column stayed fixed; each column is now sent once, with the last state of
  the batch, as HiGHSMILPSolver already did and as the changes do one at a
  time. The new `test/test_batch_fix.cpp` checks it on every registered
  `:MILPSolver`

- SCIPMILPSolver ignores a change of the bounds of a fixed ColVariable, as
  the other back-ends do: it wrote the new bound over the fixing, so that a
  Variable fixed and then given a new bound came back from SCIP with a value
  other than the one it was fixed at

- the documentation of `intThrowReducedCostException` in `MILPSolver.h` and
  `GRBMILPSolver.h` gives its default as 1, which is what
  `get_dflt_int_par()` returns, and says that a Solver to which no
  `ComputeConfig` has been applied starts with 0

- the reduced costs go to the bounds of the column they belong to also when
  the model has dynamic columns: they were handed out walking the static and
  dynamic Variable of each Block in turn, while the dynamic columns all come
  after the static ones, so that each dynamic column of a Block shifted the
  reduced costs of all the columns of the Block that follow, e.g., the duals
  of the bounds of an easy component in the master problem of a bundle

- the active set QP solver of HiGHS, which declares some convex QPs
  non-convex and stalls on others (ERGO-Code/HiGHS#3322, fixed by #3325,
  which no release of HiGHS has yet), is stopped by HiGHSMILPSolver after
  an iteration limit (the largest of 10000 and 10 times the rows and columns
  of the QP), unless the configuration sets one, and a QP that it has not
  solved is solved again from scratch, which gets it out of a stall, and if
  still not solved once more with the Hessian regularized by 1e-4; the QPs
  solved at the first attempt are solved as before

- `SCIPMILPSolver` keeps finite solutions only in the store of SCIP
  (`misc/finitesolutionstore`, now 1 by default): the presolve could fix at
  an infinite bound a free Variable that the deleted rows had left with no
  constraint, and once rows on it were added back that stored solution,
  whose rows have infinite activity, still passed the check and was taken
  for the optimal one; `batch-dynamic` and `batch-dynamic-L` compared it
  with Gurobi and failed in every pipeline

- `test_cuts` compares the optimal value with a relative tolerance of 1e-6,
  that of the feasibility of the lazy constraints in the Solver, instead of
  1e-7: Gurobi closed `batch-cuts` 1.2e-7 below the optimum on the 2 vCPU
  runners of the pipeline, a value within its tolerance that the test took
  for a wrong one

- the default and the current value of a string parameter of
  `CPXMILPSolver` and `HiGHSMILPSolver`, and the name of a parameter of
  `CPXMILPSolver`, were empty: the C library wrote them into a `std::string`
  only reserved, whose length stayed 0, so that a non-differential
  `ComputeConfig` set every string parameter it did not give to "", which
  CPLEX refuses for `CPXPARAM_CPUmask` (error 3700, nothing solved) and HiGHS
  for its enumerated options; `HiGHSMILPSolver` also wrote the current value
  of the option past the end of an empty string

- `HiGHSMILPSolver::get_dual_direction()` wrote the reduced costs, `numcols`
  of them, where `Highs_getSolution()` puts the `numrows` row values, which
  corrupted the heap when the rows are more than the columns, e.g., in the
  subproblems of a `BendersDecompositionSolver`; the multipliers of the
  columns are now those of the ray, `c - A' y`, or `- A' y` with
  `intHomogeneousDirection`, as in `CPXMILPSolver`, rather than the reduced
  costs of the current solution

- a dynamic `ColVariable` added already fixed, such as the one of a hard
  clause of `SATBlock::add_clauses()`, entered the problem with its bounds
  rather than with its value, in all the four backends, so that the solution
  could move it; `fix_bounds()` restricts the bounds of the new column to
  the value, as the loading of the problem already did

- `HiGHSMILPSolver::get_lb()` of a minimization, and `get_ub()` of a
  maximization, stopped by the time or by the iterations: they gave the value
  of the solution found so far, which is not a bound on the optimum, and a
  cross-check took it for the optimum. A MIP gives now the dual bound of the
  search (the information `mip_dual_bound` of HiGHS, which its C API exposes),
  and an LP, which has no such bound when stopped, gives -INF (+INF)

- `:HiGHSMILPSolver` loads the problem again whole, as for an
  NBModification, whenever a Modification is pending and the algorithm is
  the interior point (`solver ipm`): on a problem changed in place, HiGHS
  1.15.1 can end in a segmentation fault when the crossover is imprecise and
  the simplex cleans the solution up, which it does not on the same problem
  loaded whole. `MILPSolver::reload_on_modification()` is the hook, false
  for every other back-end

- `GRBMILPSolver::has_dual_solution()` probes the reduced costs of a model
  with no rows, which has no multiplier to probe: it used to ask Gurobi for
  the first one, which fails, and answer that there is no dual solution

- `:HiGHSMILPSolver` sizes the scheduler of threads of HiGHS per calling
  thread, since that is how HiGHS keeps it: with a count shared by the whole
  process, a thread whose scheduler had started with the default number of
  threads could be asked to run with another one (e.g., the one of a master
  with `intMaxThread` 1) while other threads were running, and HiGHS refused
  the run with an error, which is what made the subproblems of a
  `ParallelSDDPSolver` fail when all of them used HiGHS

- `CPXMILPSolver::get_lb()` and `get_ub()` give a bound when CPLEX ends in
  `kLowPrecision`, instead of an infinite one on both sides: on the side of
  the value of a solution the value is given if CPLEX says the solution is
  primal feasible, and on the side of the bound only for an LP whose basis
  CPLEX says is both primal and dual feasible, hence optimal. A
  `LagBFunction` whose inner Block is such an LP (e.g., in
  `CPX_STAT_OPTIMAL_INFEAS`) gave the bundle two infinite bounds, and the
  bundle stopped without cuts

- `GRBMILPSolver::get_lb()` for a minimization problem, and `get_ub()` for a
  maximization one, give no longer the value of a `GRB_SUBOPTIMAL` solution
  of a continuous problem as the bound, which it is not: the bound is
  infinite there, while the bound of the branch-and-bound of an integer
  problem stays

- on macOS a program linking the module lost the classes the module
  registers in the factories when the linker dropped the library, as it
  does under `-dead_strip_dylibs`, which conda sets: the target now asks the
  linker for the symbol that forces the module in (`-u`), which ld64,
  unlike the ELF linker, counts as a use of the library

- the active set QP solver of HiGHS declared some convex QPs non-convex,
  e.g., the master problem of the proximal BundleSolver, which then failed
  at most iterations: whether it did depended only on the order of the
  columns, the same model written in the order of a `.lp` file being solved.
  `qp_regularization_value` now defaults to 1e-5 instead of the 1e-7 of
  HiGHS, which solves them and moves the optimum by less than 1e-10; it only
  acts on QPs, and a configuration can still set it

- the two configurations of the testers of the duals asked for
  `CPXPARAM_MIP_Tolerances_Integrality`, a parameter of CPLEX alone and an
  integrality tolerance in a configuration for an LP: a build without CPLEX,
  where the :MILPSolver that reads them is another one, died on the name.
  It is commented out as the other parameters of CPLEX in the same files
  already are

- the testers that try each :MILPSolver in turn skipped the ones the build
  does not have by constructing them, while `Solver::new_Solver()` throws
  rather than returning `nullptr` on a name the factory does not hold, so
  `test_dual`, `test_farkas` and `test_groups` died with `CPXMILPSolver not
  present in Solver factory` wherever CPLEX is not installed; they ask
  `Solver::has_Solver()` first, and the BlockSolverConfig read by `test_dual`,
  which names one :MILPSolver, has a name the factory does not hold replaced
  by the first one it does, saying which took its place

- PIPSMILPSolver takes dynamic Constraint and Variable: the index of a static
  one, an `int` whose absence is `Inf< int >()`, was compared with
  `Inf< Index >()`, which is larger, so that every dynamic element was given
  the row or column `Inf< int >()` and the callbacks threw; the Benders form
  of `tssb_solver -k`, whose coupling is dynamic, is now solved, with the
  here-and-now Variable in the first stage of PIPS-IPM++

- the dual of a dynamic Constraint was the dual of some other row. The rows of
  the model are numbered with the static ones of every Block first and the
  dynamic ones after them, while the duals were written back by walking the
  groups of each Block in turn with a running counter, so that a dynamic row
  was handed the dual of a static one. The static rows are now written by that
  counter alone and the dynamic ones through `idx_to_dcon`, the map kept in
  the actual solver-row order, as `write_var_solution()` already did for the
  dynamic columns with `idx_to_dvar`. A Benders decomposition, whose coupling
  rows are dynamic Constraint of the Block that wraps a subproblem, read from
  them the prices of unrelated rows and built cuts that were not valid, which
  is how this came out; `test/test_dual.cpp` has the case

- GRBMILPSolver set `InfUnbdInfo` on every model, so that a certificate of
  infeasibility or of unboundedness was always there as it is in CPLEX; on
  GUROBI 13.0 a model that the dual simplex solves with that parameter on is
  answered "optimal" although it is unbounded, and the Solver returned a
  finite value for a problem that has none (12.0 answers it correctly, with
  every method). The status is now decided by a solve that leaves the
  parameter alone, and a model that turns out to be infeasible or unbounded
  is solved again with the parameter on, and with the primal simplex when it
  is unbounded, which is the method the ray comes out of; a model that is
  neither no longer gives up the reductions that parameter disables

- a group of static Variable or Constraint made of several arrays, i.e., a
  `std::vector` of `std::vector` or a `boost::multi_array` of `std::vector`,
  was mapped to the rows and columns of the matrix as if it were one array,
  from the address of its first element, and with more than one array the
  back-end got a broken problem (CPLEX crashing, HiGHS refusing it): each run
  of contiguous elements of a group is now mapped on its own

- CPXMILPSolver passed an infinite side to CPLEX when a row changed so as to
  have no bound at all (e.g., the lower bound of a `>=` row removed), and
  CPLEX answered NaN; a side changed by a Modification, or of a dynamic row,
  is now CPX_INFBOUND when infinite, as it already was when loading

- the scan of `perform_separation()`, in all four back-ends, which looked for
  the added rows only at the first level of the Modification list: a Block
  generating its dynamic Constraint inside an open channel had them arrive
  inside a GroupModification, where the scan did not see them, and the cuts
  were silently lost

- native PolyhedralFunction path in all four backends, with the helper
  `scatter_lf_to_csr` factoring CSR scatter logic across them

- batched `add_dynamic_constraints` taking a CSR matrix, with the
  base loop falling back to the per-row single-constraint API

- ignore-sub-Blocks filter promoted from this module to the
  `Solver` base class API (`Solver::set_excluded_blocks` /
  `is_excluded`), so the legacy `vstrMILPIgnSBlks` parameter has
  been retired and any Solver can now be told to skip a subset of
  the Block tree

## [0.9.1] - 2026-09-13

### Changed

- the start of the Gurobi environment is retried, with a growing wait,
  when the license service refuses it for a transient reason, so that a
  long computation is not lost to a momentary refusal

## [0.9.0] - 2026-09-12

### Added

- PIPSMILPSolver, providing the interface with the parallel
  interior-point solver PIPS-IPM++ (LP problems, Linux only),
  plus the pips\_pars header-generator tool

### Changed

- relaxed-integer LP cut-separation loop (`intRelaxIntVars == 2`) is
  now driven by the base `compute()`

- the version of the module is the git tag of its repository, or the
  VERSION.txt of a release tarball, and the shared library carries it: its
  SONAME is major.minor while the major is 0, and it is installed with an
  RPATH relative to itself, so that an installed tree keeps working wherever
  it is moved

### Fixed

- CPXMILPSolver read `cpx_idx_aux_qvar` past its end in
  `get_var_solution()` and in `cpx_index_of_dynamic_variable()`, which
  crashes on a QCP model whose quadratic constraints have no linear part,
  since those are handed to CPLEX as they are and no auxiliary variable is
  built for them, leaving that vector empty

- HiGHSMILPSolver: `lhs = rhs;` → `lhs = con_lhs;` in batched row
  addition (`add_dynamic_constraints`), which previously set the LHS
  to zero on equality rows

- GRBMILPSolver::get_lb() returned minus infinity for a continuous problem
  that GUROBI solved without a branch-and-bound, OBJBOUND being undefined
  there: where the solve ended optimal, the optimal value is the bound

- the callback of GRBMILPSolver read the solution into a buffer as long as
  the columns of the Block, while the model also has the auxiliary variables
  of the ranged constraints and of the quadratic terms and GUROBI writes all
  of them: the heap was corrupted past the end of the buffer, which then
  crashed inside GUROBI itself. It only showed with the integer variables
  on, the callback being called only in a MIP context

- FindCPLEX, FindGUROBI and FindSCIP compute ARCH when it is not given, as for
  a module built on its own or a project using the installed one, where the
  umbrella does not set it

## [0.8.0] - 2025-12-12

### Added

- possibility of redirecting log stream into specific file in
  all *MILPSolvers

- support for static Variable/Constraint to be contained in
  std::vector< std::vector > and multi\_array< std::vector >

- [big] add\_mip\_starts() in base class to provide multiple
  starting solution to the solver, then implemented in all
  *MILPSolver that support it

- callback function in HiGHSMILPSolver

- [big] has\_var\_direction() and get\_var\_direction() in all
  *MILPSolver

- option in SCIPMILPSilver to retrieve dual values by disabling
  presolve operations

- options to retrieve information during computation (callbacks)

- function to evaluate dual values for quadratic constraint
  starting from slack values

- tester for changes involving quadratic [FRow]Constraints and
  [FReal]Objective

- [huge] *MILPSolver class handle [MI-]QP and [MI-]QCQP models

- support new SMS++ Modifications (LinearFunctionModVarsAddd
  and DQuadFunctionModVarsAddd)

- added write_[var/dual]_solution()

- option IntSingleBound in the BlockSolverConfig of test_dual

### Changed

- include files *\_defs.h and *_maps.h are now automatically
  created when needed by both cmake and makefile

- using new un\_any\_thing\_count\_*()

- adapted to new CMake / makefile organisation

- promoted integer parameter intThrowReducedCostException to
  base class MILPSolver

- MILPSolver::set\_par( double ) does nothing and it was not
  defined, but this meant that set\_par( int ) was accidentally
  called in its stead

- improved CPXMILPSolver::objective\_function\_modification: the
  new coefficients are obtained by adding modl->delta() to the
  old ones

- improved query of variables bound in MILPSolver

### Fixed

- several pesky memory leaks

- set\_par( intMaxTime ) in GRBMILPSolver

- writing of dual_values in the Block

- error in SCIPMILPSolver::get\_dual\_solution()

- issues in get\_lb() and get\_ub()

- catastrophic blunder whereby the semantic of "the Objective of
  a Block is its own Objective plus the Objective of all the
  sub-Block, recursively" was *not at all* correctly implemented
  in *MILPSolver

- several minor ones

## [0.7.1] - 2024-02-01

### Added

- Separation of user cuts and lazy constraints is now possible
  with SCIP

- Support for modifying a quadratic objective function in SCIP.

## [0.7.0] - 2023-08-01

### Added

- HiGHS interface

## [0.6.0] - 2023-07-03

### Added

- Gurobi interface

## [0.5.2] - 2023-05-17

### Added

- Support for constant term in the objective function.

- Support for all-important SCIP 8.0.3.

### Fixed

- Destroy the problem before constructing a new one in
  load\_problem() (in SCIPMILPSolver and CPXMILPSolver).

## [0.5.1] - 2022-07-01

### Added

- Support to SCIP 8.0.0 and Cplex 22.1.

- Separation of user cuts and lazy constraints.

### Fixed

- Invert the sign of the dual solution in CPXMILPSolver to follow the
  RowConstraint conventions.

- Scan of ColVariable in CPXMILPSolver.

- Blunder in bound changes in CPXMILPSolver.

- Callback for Cplex versions prior to 12.10.

## [0.5.0] - 2021-12-08

### Added

- MPS support in tool and unit test.

- set\_par( idx_type , std::string && ).

### Fixed

- Bug in SCIPMILPSolver with fixed binary variables.

- Block ownership check while scanning active Constraints.

- SMS++ to CPLEX conversion of Inf values in RHS vector.

- Variable type change handling.

## [0.4.0] - 2021-05-02

### Added

- Full support for CPLEX and SCIP parameters.

- Support for dual solution and dual direction.

- Tools.

### Fixed

- Too many individual fixes to list.

## [0.3.0] - 2020-09-16

### Added

- Support for concurrency.

## [0.2.0] - 2020-03-06

### Added

- SCIP interface.

### Fixed

- Minor bugs.

## [0.1.1] - 2020-02-10

### Fixed

- Minor fix in makefile support.

## [0.1.0] - 2020-01-30

### Added

- First test release.

[Unreleased]: https://gitlab.com/smspp/milpsolver/-/compare/0.10.0...develop
[0.10.0]: https://gitlab.com/smspp/milpsolver/-/compare/0.9.1...0.10.0
[0.9.1]: https://gitlab.com/smspp/milpsolver/-/compare/0.9.0...0.9.1
[0.9.0]: https://gitlab.com/smspp/milpsolver/-/compare/0.8.0...0.9.0
[0.8.0]: https://gitlab.com/smspp/milpsolver/-/compare/0.6.0...0.8.0
[0.6.0]: https://gitlab.com/smspp/milpsolver/-/compare/0.5.2...0.6.0
[0.5.2]: https://gitlab.com/smspp/milpsolver/-/compare/0.5.1...0.5.2
[0.5.1]: https://gitlab.com/smspp/milpsolver/-/compare/0.5.0...0.5.1
[0.5.0]: https://gitlab.com/smspp/milpsolver/-/compare/0.4.0...0.5.0
[0.4.0]: https://gitlab.com/smspp/milpsolver/-/compare/0.3.0...0.4.0
[0.3.0]: https://gitlab.com/smspp/milpsolver/-/compare/0.2.0...0.3.0
[0.2.0]: https://gitlab.com/smspp/milpsolver/-/compare/0.1.1...0.2.0
[0.1.1]: https://gitlab.com/smspp/milpsolver/-/compare/0.1.0...0.1.1
[0.1.0]: https://gitlab.com/smspp/milpsolver/-/tags/0.1.0
