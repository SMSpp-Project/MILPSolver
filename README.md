# MILPSolver

[![Pipeline](https://gitlab.com/smspp/milpsolver/badges/develop/pipeline.svg)](https://gitlab.com/smspp/milpsolver/-/pipelines?ref=develop)
[![Release](https://gitlab.com/smspp/milpsolver/-/badges/release.svg)](https://gitlab.com/smspp/milpsolver/-/releases)
[![Documentation](https://img.shields.io/badge/docs-Doxygen-blue)](https://smspp.gitlab.io/smspp-project/)
[![C++ Standard](https://img.shields.io/badge/C%2B%2B-20-blue.svg)](https://en.cppreference.com/w/cpp/20)
[![License: LGPL v3](https://img.shields.io/badge/License-LGPL_v3-blue.svg)](https://www.gnu.org/licenses/lgpl-3.0)

A generic MILP Solver meta-interface for SMS++, with modules for interfacing
with some actual solvers. Supports `LinearFunction`, `DQuadFunction` (diagonal
quadratic) and `QuadFunction` (general quadratic, the last two cases not
necessarily convex) both in the `FRealObjective` and in the `FRowConstraint`,
so it is actually an interface for MI-QCQP in its full generality.

The `MILPSolver` base class (deriving from `CDASolver` for the case when the
MILP is actually an LP/QP/QCQP and therefore dual solutions are available)
provides a "meta" interface between any SMS++ `Block` whose abstract
representation is a Mixed-Integer Quadratically-Constrained Quadratic Program,
i.e.:

- all `Variable` need be `ColVariable`

- all `Objective` need be a `FRealObjective`

- the `Constraint` need all be either `FRowConstraint` or `OneVarConstraint`

- all inner `Function` in the `FRealObjective` and `FRowConstraint` are
  `LinearFunction` or `DQuadFunction` (diagonal quadratic) or `QuadFunction`
  (general quadratic, the last two cases not necessarily convex)

However, `MILPSolver` only reads the abstract representation and prepares
data structures representing the classic (sparse) coefficient matrix of
a MILP and the accompanying vectors (objective, LHS, RHS, LB, UB), plus
the other data structures for representing quadratic `Objective` /
`Constraint`, with the idea that derived classes will then use it to
interface with actual solvers. Indeed, `MILPSolver` also provides a handy
system for dealing with the `Modification` coming from the `Block`, where it
handles the changes in the internal data structures with a call to a number
of protected virtual functions that derived classes can redefine in order to
"communicate" the changes to the underlying actual MILP solvers. Yet,
`MILPSolver` can also be used directly (it is not pure virtual) in case one
just wants to read the coefficient matrix representation of a `Block`.

Currently available derived classes are:

- `CPXMILPSolver`, providing the interface with the commercial
  [IBM ILOG CPLEX](https://www.ibm.com/products/ilog-cplex-optimization-studio)

- `SCIPMILPSolver`, providing the interface with the open-source
  [SCIP](https://www.scipopt.org) (note that since version 8.0.3 SCIP is
  "truly" FOSS by dint of being distributed under the Apache 2.0 License as
  opposed to the previous academic license preventing royalty-free commercial
  use)

- `GRBMILPSolver`, providing the interface with the commercial
  [GUROBI Optimizer](https://www.gurobi.com/solutions/gurobi-optimizer)

- `HiGHSMILPSolver`, providing the interface with the open-source
  [HiGHS](https://highs.dev)

- `MOSEKMILPSolver`, providing the interface with the commercial
  [MOSEK](https://www.mosek.com) (LP, MILP, convex QP and MIQP problems
  with linear constraints)

- `PIPSMILPSolver`, providing the interface with the open-source parallel
  interior-point solver
  [PIPS-IPM++](https://gitlab.com/pips-ipmpp/pips-ipmpp) (LP problems
  only, Linux only; see also its
  [website](https://pips-ipmpp.gitlab.io/index.html))

- `DAQPMILPSolver`, providing the interface with the open-source dual
  active-set solver [DAQP](https://github.com/darnstrom/daqp) (convex QP
  and LP problems with linear constraints, dense, no integer variables)

- `QPOASESMILPSolver`, providing the interface with the open-source online
  active-set solver [qpOASES](https://github.com/coin-or/qpOASES) (convex QP
  and LP problems with linear constraints, dense, no integer variables, the
  Hessian possibly singular, hot starts between similar problems)

- `ClarabelMILPSolver`, providing the interface with the open-source
  interior-point solver [Clarabel](https://clarabel.org) (convex QP and LP
  problems with linear constraints, sparse, no integer variables, with the
  unbounded rays and the Farkas certificates of infeasible problems)

- `ProxQPMILPSolver`, providing the interface with the open-source proximal
  augmented Lagrangian solver ProxQP of
  [ProxSuite](https://github.com/Simple-Robotics/proxsuite) (convex QP and
  LP problems with linear constraints, dense or sparse, no integer
  variables, with the unbounded rays)

- `PIQPMILPSolver`, providing the interface with the open-source proximal
  interior point solver [PIQP](https://github.com/PREDICT-EPFL/piqp)
  (convex QP and LP problems with linear constraints, dense or sparse, no
  integer variables)

- `OSQPMILPSolver`, providing the interface with the open-source operator
  splitting solver [OSQP](https://osqp.org) (convex QP and LP problems with
  linear constraints, sparse, no integer variables, with the unbounded rays
  and the Farkas certificates of infeasible problems)

- `QPMADMILPSolver`, providing the interface with the open-source dual
  active-set solver [qpmad](https://github.com/asherikov/qpmad), of
  Goldfarb and Idnani (strictly convex QP problems with linear constraints,
  dense, the Hessian positive definite, no integer variables)

- `QuadProgMILPSolver`, providing the interface with the open-source dual
  active-set solver [quadprog](https://github.com/quadprog/quadprog), the
  original implementation of the method of Goldfarb and Idnani (same class
  of problems as `QPMADMILPSolver`); quadprog is under the GPL, hence a
  library built with it is too

Whatever the back-end, the solution of a continuous problem can be checked
by `MILPSolver` itself against the optimality conditions: with the parameter
`dblKKTTol` positive, the primal residual, the dual residual and the gap of
the solution are measured on the `Block` after each solve, and a solution
that exceeds the tolerance is reported as `kError` (only the primal residual
is measured on a MIP, or when the back-end gives no duals).

Recent versions of every underlying solver are supported via a mechanism
that auto-generates the `*_defs.h` and `*_maps.h` files for the version
found on the system, both with CMake and with the makefiles (see below
for details). Older versions may fail if the upstream API changed: should
that happen, please upgrade.


## Getting started

These instructions will let you build MILPSolver on your system.

The module also comes ready-made, in any of

```sh
sudo add-apt-repository ppa:smspp-project/smspp   # Ubuntu
sudo apt install libsmspp-milp-dev

conda install -c conda-forge smspp-project        # Linux, macOS, Windows

brew tap SMSpp-Project/smspp                      # macOS, Linux
brew install smspp

vcpkg install "smspp[core,milp]"                  # from the sources
```

where apt and the port give the module alone, while conda and the tap carry
the whole framework. What follows is about building it yourself.

### Requirements

- [SMS++ core library](https://gitlab.com/smspp/smspp)

- for `CPXMILPSolver` you will need
  [IBM ILOG CPLEX](https://www.ibm.com/products/ilog-cplex-optimization-studio)

- for `SCIPMILPSolver` you will need
  [SCIP](https://www.scipopt.org)

- for `GRBMILPSolver` you will need
  [GUROBI Optimizer](https://www.gurobi.com/solutions/gurobi-optimizer)

- for `HiGHSMILPSolver` you will need [HiGHS](https://highs.dev) 1.8.0 or
  later

- for `MOSEKMILPSolver` you will need [MOSEK](https://www.mosek.com) 11.0.1
  or later, in `MOSEK_ROOT` (the directory of a version, e.g.
  `/opt/mosek/11.0`, or the one above it), and a license of it in
  `~/mosek/mosek.lic` or wherever `MOSEKLM_LICENSE_FILE` says (MOSEK gives
  [free academic licenses](https://www.mosek.com/products/academic-licenses/))

- for `PIPSMILPSolver` you will need
  [PIPS-IPM++](https://gitlab.com/pips-ipmpp/pips-ipmpp), which is only
  available under Linux and requires an MPI implementation (note that
  PIPS-IPM++ also needs a sparse linear solver: the fully open-source
  route is MUMPS, while other backends, e.g. PARDISO or MA57, are
  subject to licensing)

- for `DAQPMILPSolver` you will need [DAQP](https://github.com/darnstrom/daqp)
  0.8.2 or later, installed with its CMake (`cmake --install`), in
  `DAQP_ROOT`

- for `QPOASESMILPSolver` you will need
  [qpOASES](https://github.com/coin-or/qpOASES) 3.2 or later, installed
  with its CMake (`cmake --install`), in `qpOASES_ROOT`

- for `ClarabelMILPSolver` you will need the C interface of Clarabel,
  [Clarabel.cpp](https://github.com/oxfordcontrol/Clarabel.cpp) 0.11.1 or
  later, whose build needs a Rust toolchain (`cargo`): since its CMake
  installs nothing, copy its `include` directory and the `libclarabel_c`
  built by `cargo build --release` in `rust_wrapper` to the `include` and
  `lib` directories of `Clarabel_ROOT`

- for `ProxQPMILPSolver` you will need
  [ProxSuite](https://github.com/Simple-Robotics/proxsuite) 0.4 or later,
  which is header-only and needs Eigen (and Simde, unless it is configured
  with `BUILD_WITH_VECTORIZATION_SUPPORT=OFF`), installed with its CMake
  (`cmake --install`) in `proxsuite_ROOT`

- for `PIQPMILPSolver` you will need
  [PIQP](https://github.com/PREDICT-EPFL/piqp) 0.6 or later, which needs
  Eigen, installed with its CMake (`cmake --install`) in `piqp_ROOT`; its
  C++ interface is enough (`BUILD_C_INTERFACE=OFF`), and its compiled
  library (the default `BUILD_WITH_TEMPLATE_INSTANTIATION=ON`) is used if
  it has been built, its headers alone otherwise

- for `OSQPMILPSolver` you will need [OSQP](https://github.com/osqp/osqp)
  1.0 or later, installed with its CMake (`cmake --install`) in
  `osqp_ROOT`, with its builtin linear algebra and in double precision (the
  defaults); its static library is used if it has been built (the
  default), its shared one otherwise. OSQP exports its internal functions,
  and one of them (`csc_copy`) is also in HiGHS: on Linux, its static
  library is linked with its symbols kept local, while with its shared one
  alone `OSQPMILPSolver` crashes in a program that also loads HiGHS

- for `QPMADMILPSolver` you will need
  [qpmad](https://github.com/asherikov/qpmad) 1.5 or later, which is
  header-only and needs Eigen, installed with its CMake (`cmake --install`)
  in `qpmad_ROOT`

- for `QuadProgMILPSolver` you will need
  [quadprog](https://github.com/quadprog/quadprog), which has no library of
  its own: its three C sources (`solve.QP.c`, `linear-algebra.c`,
  `qr-update.c`) are compiled into `libquadprog` in `quadprog_ROOT/lib`,
  exporting only `qpgen2_`, as `INSTALL.sh --with-quadprog` of the umbrella
  does

All actual `:MILPSolver` are optional but you will need at least one of them to
actually solve MILP/LP problems. Without any of them, you can still build a
`MILPSolver` that loads the problem from the SMS++ `Block` and makes it
available as a (sparse) coefficient matrix + accompanying vectors (objective,
LHS, RHS, LB, UB) + other data structures for representing quadratic
`Objective` / `Constraint`.

### Build and install with CMake

Configure and build the library with:

```sh
mkdir build
cd build
cmake ..
cmake --build .
```

The library has the same configuration options of
[SMS++](https://gitlab.com/smspp/smspp-project/-/wikis/Customize-the-configuration).
Moreover, you can use the following configuration options:

| Variable                  | Description    | Default value |
| ------------------------- | -------------- | ------------- |
| `MILPSolver_USE_CPLEX`    | Use CPLEX      | ON            |
| `MILPSolver_USE_SCIP`     | Use SCIP       | ON            |
| `MILPSolver_USE_GUROBI`   | Use GUROBI     | ON            |
| `MILPSolver_USE_HiGHS`    | Use HiGHS      | ON            |
| `MILPSolver_USE_PIPS`     | Use PIPS-IPM++ | ON            |
| `MILPSolver_USE_DAQP`     | Use DAQP       | ON            |
| `MILPSolver_USE_qpOASES`  | Use qpOASES    | ON            |
| `MILPSolver_USE_Clarabel` | Use Clarabel   | ON            |
| `MILPSolver_USE_MOSEK`    | Use MOSEK      | ON            |
| `MILPSolver_USE_ProxQP`   | Use ProxQP     | ON            |
| `MILPSolver_USE_PIQP`     | Use PIQP       | ON            |
| `MILPSolver_USE_OSQP`     | Use OSQP       | ON            |
| `MILPSolver_USE_qpmad`    | Use qpmad      | ON            |
| `MILPSolver_USE_quadprog` | Use quadprog   | ON            |

At least one of the backends must be enabled at build time for the
library to be able to actually solve problems.

Optionally, install the library in the system with:

```sh
cmake --install .
```

### Usage with CMake

After the module is built, you can use it in your CMake project with:
```cmake
find_package(MILPSolver)
target_link_libraries(<my_target> SMS++::MILPSolver)
```

### Running the tests with CMake

Some unit tests will be built with the library.
Launch `ctest` from the build directory to run them.
To disable them, set the option `BUILD_TESTING` to `OFF`.


### Build and install with makefiles

Carefully hand-crafted makefiles have also been developed for those unwilling
to use CMake. Makefiles build the executable in-source (in the same directory
tree where the code is) as opposed to out-of-source (in the copy of the
directory tree constructed in the build/ folder) and therefore it is more
convenient when having to recompile often, such as when developing/debugging
a new module, as opposed to the compile-and-forget usage envisioned by CMake.

Each executable using `MILPSolver` and/or one of the derived `:MILPSolver` has
to include a "main makefile" of the module, which typically is either
[makefile-c](makefile-c) including all necessary libraries comprised the
"core SMS++" one, or [makefile-s](makefile-s) including all necessary
libraries but not the "core SMS++" one (for the common case in which this is
used together with other modules that already include them). If you want to
exclude some specific `:MILPSolver` from being compiled you have to go in
[makefile](makefile), [makefile-c](makefile-c) and [makefile-s](makefile-s)
(depending on which one of the latter two is used) and comment out all the
lines mentioning it. Don't bother about the `$(*H)`, `$(*INC)` etc. variables
(assuming you would) since if they are not defined they are empty and
therefore do no harm. Relevant cases are the testers described below.
The makefiles in turn recursively include all the required other makefiles,
hence one should only need to edit the "main makefile" for compilation type
(C++ compiler and its options) and it all should be good to go. In case some
of the external libraries are not at their default location, it should only be
necessary to create the `../extlib/makefile-paths` out of the
`extlib/makefile-default-paths-*` for your OS `*` and edit the relevant bits
(commenting out all the rest).

Check the [SMS++ installation wiki](https://gitlab.com/smspp/smspp-project/-/wikis/Customize-the-configuration#location-of-required-libraries)
for further details.


## Tools

All the available `:MILPSolver`s support the parameter names of the
corresponding solver in the `Configuration` files. To do so, they need
header files *\_defs.h and *\_maps.h that depend on the versions of the
solver currently installed on the system. Such headers can be generated
with the `cpx_pars` , `scip_pars` , `grb_pars` , `highs_pars` and
`pips_pars` executables in the [tools](tools) folder. This is done
automatically by cmake / make, so you should not bother about it. However,
if you change the version of the underlying solver you may want to delete
the corresponding *\_defs.h and *\_maps.h header files so that they are
rebuilt for the new one.


## Testers

The [test](test/README.md) folder includes some testers that may be useful
for someone willing to write other `:MILPSolver`:

- `test_cuts` tests dynamic generation of constraints

- `test_dual` tests correct signs of dual variables (never to be given for
  granted, every LP solver seems to have a different idea about it)

- `test_dynamic` compares two `:MILPSolver` for the repeated solution of LPs
  changing everything that can be changed, useful to test a new `:MILPSolver`
  against an old and hopefully reliable one


## Known issues

The active set QP solver of HiGHS (up to 1.15.1 at least) may declare
non-convex a convex QP, depending on the order of its columns, or stall on
it ([ERGO-Code/HiGHS#3322](https://github.com/ERGO-Code/HiGHS/issues/3322),
fixed by [#3325](https://github.com/ERGO-Code/HiGHS/pull/3325), which no
release of HiGHS contains yet). `HiGHSMILPSolver` works around it for
continuous QPs: the Hessian is regularized by 1e-5 instead of the 1e-7 of
HiGHS, the active set is stopped after an iteration limit, and a QP that it
has not solved is solved again from scratch and, if needed, once more with
the Hessian regularized by 1e-4. The QPs that HiGHS solves at the first
attempt are not affected; all the details are in the documentation of
`HiGHSMILPSolver`, and the workaround will be removed once a release of HiGHS
with the fix has been checked.


## Getting help

If you need support, you want to submit bugs or propose a new feature, you can
[open a new issue](https://gitlab.com/smspp/milpsolver/-/issues/new).


## Contributing

Please read [CONTRIBUTING.md](CONTRIBUTING.md) for details on our code of
conduct, and the process for submitting merge requests to us.


## Authors

### Current Lead Authors

- **Enrico Calandrini**  
  Dipartimento di Informatica  
  Universita' di Pisa

- **Antonio Frangioni**  
  Dipartimento di Informatica  
  Università di Pisa

- **Donato Meoli**  
  Dipartimento di Informatica  
  Università di Pisa

### Contributors

- **Rafael Durbano Lobato**  
  Dipartimento di Informatica  
  Università di Pisa

- **Niccolò Iardella**  
  Dipartimento di Informatica  
  Università di Pisa


## License

This code is provided free of charge under the [GNU Lesser General Public
License version 3.0](https://opensource.org/licenses/lgpl-3.0.html) -
see the [LICENSE](LICENSE) file for details.


## Disclaimer

The code is currently provided free of charge under an open-source license.
As such, it is provided "*as is*", without any explicit or implicit warranty
that it will properly behave or it will suit your needs. The Authors of
the code cannot be considered liable, either directly or indirectly, for
any damage or loss that anybody could suffer for having used it. More
details about the non-warranty attached to this code are available in the
license description file.

The authors are not affiliated, associated, authorized, endorsed by, or in
any way officially connected with IBM, Gurobi or MOSEK, or any of their
subsidiaries or affiliates. The names IBM, ILOG, CPLEX, Gurobi and MOSEK as
well as related names, marks, emblems and images are registered trademarks
of their respective owners. The authors also do not claim any affiliation with the
open-source projects developing SCIP and HiGHS: any fault in the code
interfacing SMS++ with those is all and entirely ours.

