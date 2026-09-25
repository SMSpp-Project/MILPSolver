# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

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

[Unreleased]: https://gitlab.com/smspp/milpsolver/-/compare/0.9.1...develop
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
