# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [0.7.1] - 2024-02-01

### Added

- Separation of user cuts and lazy constraints is now possible with SCIP
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

- Destroy the problem before constructing a new one in load_problem() (in
  SCIPMILPSolver and CPXMILPSolver).

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
- set_par( idx_type , std::string && ).

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

[Unreleased]: https://gitlab.com/smspp/milpsolver/-/compare/0.5.2...develop
[0.5.2]: https://gitlab.com/smspp/milpsolver/-/compare/0.5.1...0.5.2
[0.5.1]: https://gitlab.com/smspp/milpsolver/-/compare/0.5.0...0.5.1
[0.5.0]: https://gitlab.com/smspp/milpsolver/-/compare/0.4.0...0.5.0
[0.4.0]: https://gitlab.com/smspp/milpsolver/-/compare/0.3.0...0.4.0
[0.3.0]: https://gitlab.com/smspp/milpsolver/-/compare/0.2.0...0.3.0
[0.2.0]: https://gitlab.com/smspp/milpsolver/-/compare/0.1.0...0.2.0
[0.1.1]: https://gitlab.com/smspp/milpsolver/-/compare/0.1.0...0.1.1
[0.1.0]: https://gitlab.com/smspp/milpsolver/-/tags/0.1.0
