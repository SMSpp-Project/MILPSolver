# /MILPSolver/test_cuts

A tester which provides tests for any `:MILPSolver` dynamically separating
cuts via `generate_dynamic_constraint()`.

Given one input parameter n, a `NCoCubeBlock` with that size is defined.
The `NCoCubeBlock` is a `Block` that encodes for the N-Co-Cube, i.e., the
unitary ball in the L~1~ norm, with a simple linear `Objective`. Optimizing
over this is trivial, since the N-Co-Cube only has 2n vertices of the form
[ 0 , 0 , ... +/- 1 , ... , 0 ], and therefore also has the integrality
property. However, NCoCubeBlock implements the "crazy" formulation in the
original variable space which has 2^n^ constraints of the form

   ( +/- 1 ) x~1~ + ( +/- 1 ) x~2~ + ... + ( +/- 1 ) x~n~ <= 1

Of course these are implemented as dynamic constraint, whose separation
(trivial and linear) is performed in `generate_dynamic_constraint()`.

The problem is repeatedly solved a number of times randomly changing the
`Objective` (which is the only thing that can change), and the optimal
value is compared with that obtained by the trivial optimization.

The specific `:MILPSolver` used for the tests and is algorithmic parameters
are specified in the [MILPPar.txt](MILPPar.txt) `BlockSolverConfig` file.

The usage of the executable is the following:

       ./test_cuts seed [nvar n_repeat]
       nvar: number of variables [10]
       n_repeat: number of repetitions [10]

A batch file is provided that runs a not-so-large set of tests with
different sizes and seeds of the random generator; all these passing is a
good sign that the given `:MILPSolver` correctly handles dynamic generation
of constraints.

A makefile is also provided that builds the executable including the
`MILPSolver` module and all its dependencies (and, obviously, the 
core SMS++ library).


## Authors

- **Antonio Frangioni**  
  Dipartimento di Informatica  
  Università di Pisa


## License

This code is provided free of charge under the [GNU Lesser General Public
License version 3.0](https://opensource.org/licenses/lgpl-3.0.html) -
see the [LICENSE](../LICENSE) file for details.
