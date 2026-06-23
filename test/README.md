# test

Three testers for the `MILPSolver` module, each registered as a separate
`ctest` and built by the shared `makefile`.

- `test_cuts` checks that a `:MILPSolver` correctly separates dynamic
  constraints via `generate_dynamic_constraint()`. It builds an
  `NCoCubeBlock` (the unit ball in the L1 norm, encoded with its
  exponentially many "crazy" constraints, separated trivially) and compares
  the repeatedly re-solved optimal value against the trivial optimization.
  The solver and its parameters are read from `MILPPar.txt`; see `batch-cuts`.

- `test_dual` verifies that the dual solution returned by a `:MILPSolver`
  respects the convention established by `RowConstraint`. The solver is
  configured by `LPPar-dual.txt`.

- `test_dynamic` provides comprehensive tests for any `CDASolver` able to
  handle Linear Programs: it builds a random LP in an `AbstractBlock`, solves
  it with two attached solvers, then repeatedly modifies the LP (rows,
  variables, bounds, constants) and re-solves, comparing the results each
  time. The solvers are configured by `LPPar-dynamic.txt`; see `batch-dynamic`
  and `batch-dynamic-L`.

Run a tester directly (for instance `./test_dynamic seed [...]`, pass no
argument to see its usage), or run a batch over a set of sizes and seeds.

The `makefile` builds all three executables including the `MILPSolver`
module and all its dependencies, and the core SMS++ library.


## Authors

- **Antonio Frangioni**  
  Dipartimento di Informatica  
  Università di Pisa

- **Rafael Durbano Lobato**  
  Dipartimento di Informatica  
  Università di Pisa

- **Enrico Calandrini**  
  Dipartimento di Informatica  
  Università di Pisa


## License

This code is provided free of charge under the [GNU Lesser General Public
License version 3.0](https://opensource.org/licenses/lgpl-3.0.html),
see the [LICENSE](../LICENSE) file for details.
