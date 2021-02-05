# MILPSolver

A MILP Solver for SMS++ that uses [IBM® ILOG® CPLEX® Optimization Studio] and
[SCIP].

## Getting started

These instructions will let you build MILPSolver on your system.

### Requirements

- [SMS++ core library](https://gitlab.com/smspp/smspp)
- [IBM® ILOG® CPLEX® Optimization Studio] (supported versions: 12.8 - 12.10)
- [SCIP] (supported versions: 7.0.0 - 7.0.1)

Both CPLEX and SCIP are optional but you will need at least one of them to
build a solver that actually solves problems.
Without either of them, you can still build a solver that loads the problem
from the SMS++ blocks and makes it available as a set of vectors.

### Build and install with CMake

Configure and build the library with:

```sh
mkdir build
cd build
cmake ..
make
```

The library has the same configuration options of
[SMS++](https://gitlab.com/smspp/smspp/wikis/custom).
Moreover, you can use the following configuration options:

| Variable               | Description | Default value |
| ---------------------- | ----------- | ------------- |
| `MILPSolver_USE_CPLEX` | Use CPLEX   | ON            |
| `MILPSolver_USE_SCIP`  | Use SCIP    | ON            |

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

> **Note:**
> The tests use [Google Test](https://github.com/google/googletest).
> CMake will fetch and build it automatically.

## Tools

The repository contains some tools that are built with the library.

### Solver

The `milp_solver` tool reads MILP problems from text files and solves them.
Optionally, it writes back the problem in netCDF format.

You can find some example input files in the [`test`](test) directory.
For the block and solver configuration file format,
see the [SMS++ core library](https://gitlab.com/smspp/smspp).

```sh
Usage:
milp_solver [options] <file>
milp_solver -h | --help

Options:
-B, --blockcfg <file>    Block configuration.
-S, --solvercfg <file>   Solver configuration.
-n, --nc4problem <file>  Write nc4 problem on file.
-v, --verbose            Make the solver verbose.
-h, --help               Print this help.
```

### Parameter generators

`CPXMILPSolver` and `SCIPMILPSolver` support,
respectively, CPLEX and SCIP parameter names in the configuration files.
To do so, they need header files that depend on the versions of CPLEX
and SCIP currently installed on the system;
such headers can be generated with the `cpx_pars` and `scip_pars` tools.

> **Note:**
> We provide header files for the versions we already support, so you will
> need these tools only if you have an unsupported version of either CPLEX or SCIP.

## Getting help

If you need support, you want to submit bugs or propose a new feature, you can
[open a new issue](https://gitlab.com/smspp/milpsolver/-/issues/new).

## Contributing

Please read [CONTRIBUTING.md](CONTRIBUTING.md) for details on our code of
conduct, and the process for submitting merge requests to us.

## Authors

### Current Lead Authors

- **Antonio Frangioni**  
  *Operations Research Group*  
  Dipartimento di Informatica  
  Università di Pisa

- **Niccolò Iardella**  
  *Operations Research Group*  
  Dipartimento di Informatica  
  Università di Pisa

### Contributors

- **Rafael Durbano Lobato**  
  Department of Applied Mathematics  
  State University of Campinas, Brazil

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

The authors are not affiliated, associated, authorized, endorsed by,
or in any way officially connected with IBM, or any of its subsidiaries or its affiliates.
The names IBM, ILOG and CPLEX as well as related names, marks, emblems and
images are registered trademarks of their respective owners.

[IBM® ILOG® CPLEX® Optimization Studio]: https://www.ibm.com/products/ilog-cplex-optimization-studio
[SCIP]: https://scipopt.org/index.php
