# MILPSolver

A generic MILP Solver meta-interface for SMS++, with modules for interfacing
with some actual solvers.

The `MILPSolver` base class (deriving from `CDASolver` for the case when the
MILP is actually an LP and therefore dual solutions are available) provides
a "meta" interface between any SMS++ `Block` whose abstract representation is
a Mixed-Integer Linear Program (all `Variable` need be `ColVariable`, all
`Objective` need be a `FRealObjective` whose inner `Function` is a
`LinearFunction`, the `Constraint` need all be either `FRowConstraint` whose
inner `Function` is a `LinearFunction` or `OneVarConstraint`). In fact, the
class is slightly extended to Mixed-Integer Quadratic Program with separable
objective  (the `Function` inside the `FRealObjective` can also be a
`DQuadFunction`). However, `MILPSolver` only reads the abstract representation
and prepares data structures representing the classic (sparse) coefficient
matrix + accompanying vectors (objective, LHS, RHS, LB, UB) with the idea that
derived classes will then use it to interface with actual solvers. Indeed,
`MILPSolver` also provides a handy system for dealing with the `Modification`
coming from the `Block`, where it handles the changes in the internal data
structures with a call to a number of protected virtual functions that
derived classes can redefine in order to "communicate" the changes to the
underlying actual MILP solvers. Yet, `MILPSolver` can also be used directly
(it is not pure virtual) in case one just wants to read the coefficient matrix
representation of a `Block`.

Currently available derived classes are:

- `CPXMILPSolver`, providing the interface with the commercial
  [IBM ILOG CPLEX](https://www.ibm.com/products/ilog-cplex-optimization-studio)

- `SCIPMILPSolver`, providing the interface with the open-source
  [SCIP](https://www.scipopt.org) (note that since version 8.0.3 SCIP is
  "truly" FOSS by dint of being distributed under the Apache 2.0 License as
  opposed to the previous academic license preventing roialty-free commercial
  use)

- `GRBMILPSolver`, providing the interface with the commercial
  [GUROBI Optimizer](https://www.gurobi.com/solutions/gurobi-optimizer)

- `HiGHSMILPSolver`, providing the interface with the open-source
  [HiGHS](https://highs.dev)

## Getting started

These instructions will let you build MILPSolver on your system.

### Requirements

- [SMS++ core library](https://gitlab.com/smspp/smspp)

- for `CPXMILPSolver` you will need
  [IBM ILOG CPLEX](https://www.ibm.com/products/ilog-cplex-optimization-studio)
  (currently supported versions: 12.8, 12.10, 20.10, 22.01)

- for `SCIPMILPSolver` you will need
  [SCIP](https://www.scipopt.org) (currently supported versions: 7.0.0, 7.0.1,
  7.0.2, 7.0.3, 8.0.0, 8.0.3, , 8.1.0)

- for `GRBMILPSolver` you will need
  [GUROBI Optimizer](https://www.gurobi.com/solutions/gurobi-optimizer)
  (currently supported versions: 10.0.1, 10.0.2)

- for `HiGHSMILPSolver` you will need
  [HiGHS](https://highs.dev) (currently supported versions: 1.5.3)

All actual `:MILPSolver` are optional but you will need at least one of them to
actually solve MILP/LP problems. Without any of them, you can still build a
`MILPSolver` that loads the problem from the SMS++ `Block` and makes it
available as a (sparse) coefficient matrix + accompanying vectors (objective,
LHS, RHS, LB, UB).

### Build and install with CMake

Configure and build the library with:

```sh
mkdir build
cd build
cmake ..
make
```

The library has the same configuration options of
[SMS++](https://gitlab.com/smspp/smspp-project/-/wikis/Customize-the-configuration).
Moreover, you can use the following configuration options:

| Variable                 | Description | Default value |
| ------------------------ | ----------- | ------------- |
| `MILPSolver_USE_CPLEX`   | Use CPLEX   | ON            |
| `MILPSolver_USE_SCIP`    | Use SCIP    | ON            |
| `MILPSolver_USE_GUROBI`  | Use GUROBI  | ON            |
| `MILPSolver_USE_HiGHS`   | Use HiGHS   | ON            |

Optionally, install the library in the system with:

```sh
sudo make install
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
[makefile](makefile), [makefile-c](makefile-c) and [makefile-](makefile-s)
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

`CPXMILPSolver` , `SCIPMILPSolver` , `GRBMILPSolver` and `HiGHSMILPSolver` 
support, respectively, CPLEX , SCIP , GUROBI and HiGHS parameter names in 
the `Configuration` files. To do so, they need header files that depend on 
the versions of CPLEX , SCIP , GUROBI and HiGHS currently installed on the
system; such headers can be generated with the `cpx_pars` , `scip_pars` ,
`grb_pars` and `high_pars` executables in the [tools](tools) folder.

> **Note:**
> We provide header files for the versions we already support, so you will
> need these tools only if you have an unsupported version of either CPLEX ,
> SCIP , GUROBI or HiGHS.


## Testers

The repo includes some testers that may be useful for someone willing to
write other `:MILPSolver`:

- [test_cuts](test_cuts/README.md) tests dynamic generation of constraints

- [test_dual](test_dual/README.md) tests correct signs of dual variables
  (never to be given for granted, every LP solver seems to have a different
  idea about it)

- [test_dynamic](test_dynamic/README.md) compares two `:MILPSolver` for
  the repeted solution of LPs changing everything that can be changed,
  useful to test a new `:MILPSolver` against an old an hopefully reliable
  one


## Getting help

If you need support, you want to submit bugs or propose a new feature, you can
[open a new issue](https://gitlab.com/smspp/milpsolver/-/issues/new).


## Contributing

Please read [CONTRIBUTING.md](CONTRIBUTING.md) for details on our code of
conduct, and the process for submitting merge requests to us.


## Authors

### Current Lead Authors

- **Antonio Frangioni**  
  Dipartimento di Informatica  
  Università di Pisa

- **Enrico Calandrini**  
  Dipartimento di Informatica  
  Universita' di Pisa

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
any way officially connected with IBM or Gurobi, or any of its subsidiaries
or its affiliates. The names IBM, ILOG, CPLEX and Gurobi as well as related
names, marks, emblems and images are registered trademarks of their
respective owners.

