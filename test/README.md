# test

Three testers for the `MILPSolver` module, each registered as a separate
`ctest` and built by the shared `makefile`.

- `test_cuts` checks that a `:MILPSolver` correctly separates dynamic
  constraints via `generate_dynamic_constraint()`. It builds an
  `NCoCubeBlock` (the unit ball in the L1 norm, encoded with its
  exponentially many "crazy" constraints, separated trivially) and compares
  the repeatedly re-solved optimal value against the trivial optimization.
  The solver and its parameters are read from `MILPPar.txt`; see `batch-cuts`.

- `test_batch_fix` issues on one channel changes that name the same column
  more than once (a `ColVariable` fixed and unfixed, the bound of a
  `BoxConstraint` changed twice, the two mixed) and checks that every
  registered `:MILPSolver` ends up with the last state, as it does when the
  same changes arrive one at a time.

- `test_dual` verifies that the dual solution returned by a `:MILPSolver`
  respects the convention established by `RowConstraint`. The solver is
  configured by `LPPar-dual.txt`.

- `test_dynamic` provides comprehensive tests for any `CDASolver` able to
  handle Linear Programs: it builds a random LP in an `AbstractBlock`, solves
  it with two attached solvers, then repeatedly modifies the LP (rows,
  variables, bounds, constants) and re-solves, comparing the results each
  time. The solvers are configured by `LPPar-dynamic.txt`; see `batch-dynamic`
  and `batch-dynamic-L`.

- `test_reduced_costs` verifies that the reduced costs become the duals of
  the bounds of the right Variable when the model has dynamic Variable in a
  Block that comes before the static ones, every registered `:MILPSolver`
  solving a small LP whose duals are known.

- `test_stopped_bounds` verifies that each `:MILPSolver` of the build,
  stopped by its time limit on a continuous problem, gives on the side of
  the optimum it has not proved an infinite or a valid bound, and not zero
  or a large finite number, and that once the time limit is lifted on the
  same Solver it solves the problem, with the optimal value as both
  bounds.

- `test_free_rows` verifies that each `:MILPSolver` of the build solves
  again, warm and right, a small LP with rows free on both sides after a
  change of the coefficients of other rows, after rows are relaxed and
  enforced again, and after the free rows get a finite side and lose it,
  each value being also that of a Solver loaded from scratch.

- `test_ipm_free_rows` verifies that the interior point of each
  `:MILPSolver` of the build that has one, without presolve, solves right
  a small LP with a free row followed by a row with two finite sides, and
  solves it again after the free row gets a finite side and loses it and
  after a row is relaxed and enforced again; the Solver and their options
  are those of `LPPar-ipm.txt`.

- `test_cutoff_status` verifies that each `:MILPSolver` of the build,
  given an upper cutoff that the first integer solution of a maximization
  MIP reaches, ends with `kOK` and the optimal value or with `kError`, and
  does not throw; the Solver are those of `MILPPar-cutoff.txt`, all
  but `SCIPMILPSolver`, which fails on the second solve with a cutoff.

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
