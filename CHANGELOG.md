# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

- native PolyhedralFunction path in all four backends, with the helper
  `scatter_lf_to_csr` factoring CSR scatter logic across them

- batched `add_dynamic_constraints` taking a CSR matrix, with the
  base loop falling back to the per-row single-constraint API

- ignore-sub-Blocks filter promoted from this module to the
  `Solver` base class API (`Solver::set_excluded_blocks` /
  `is_excluded`), so the legacy `vstrMILPIgnSBlks` parameter has
  been retired and any Solver can now be told to skip a subset of
  the Block tree

### Changed

- relaxed-integer LP cut-separation loop (`intRelaxIntVars == 2`) is
  now driven by the base `compute()`

### Fixed

- HiGHSMILPSolver: `lhs = rhs;` → `lhs = con_lhs;` in batched row
  addition (`add_dynamic_constraints`), which previously set the LHS
  to zero on equality rows


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

[Unreleased]: https://gitlab.com/smspp/milpsolver/-/compare/0.8.0...develop
[0.8.0]: https://gitlab.com/smspp/milpsolver/-/compare/0.7.1...0.8.0
[0.7.1]: https://gitlab.com/smspp/milpsolver/-/compare/0.7.0...0.7.1
[0.7.0]: https://gitlab.com/smspp/milpsolver/-/compare/0.6.0...0.7.0
[0.6.0]: https://gitlab.com/smspp/milpsolver/-/compare/0.5.2...0.6.0
[0.5.2]: https://gitlab.com/smspp/milpsolver/-/compare/0.5.1...0.5.2
[0.5.1]: https://gitlab.com/smspp/milpsolver/-/compare/0.5.0...0.5.1
[0.5.0]: https://gitlab.com/smspp/milpsolver/-/compare/0.4.0...0.5.0
[0.4.0]: https://gitlab.com/smspp/milpsolver/-/compare/0.3.0...0.4.0
[0.3.0]: https://gitlab.com/smspp/milpsolver/-/compare/0.2.0...0.3.0
[0.2.0]: https://gitlab.com/smspp/milpsolver/-/compare/0.1.0...0.2.0
[0.1.1]: https://gitlab.com/smspp/milpsolver/-/compare/0.1.0...0.1.1
[0.1.0]: https://gitlab.com/smspp/milpsolver/-/tags/0.1.0
