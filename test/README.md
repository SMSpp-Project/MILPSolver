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
  after a row is relaxed and enforced again, and that the model it writes
  after the first solve is the one it writes before; the Solver and their
  options are those of `LPPar-ipm.txt`.

- `test_cutoff_status` verifies the cutoffs of each `:MILPSolver` of the
  build on a maximization and a minimization MIP with a constant in the
  objective: the cutoffs of the `ComputeConfig` of `MILPPar-cutoff.txt`
  (an upper one, a lower one, both) have to survive loading the problem,
  the first solve has to stop (`kCutOff`, with the value of its solution
  beyond the cutoff and a valid bound) at a cutoff that the first solution
  reaches and to give the optimal value otherwise, and a second solve with
  a cutoff after a first one has to give either, and not throw.

- `test_free_workspace` verifies that every `:MILPSolver` of
  `LPPar-all.txt` that the build has, each of which may free after a solve
  the memory that its back-end keeps and the next solve does not need,
  gives the same status, value, primal values and duals as the others,
  across warm re-solves after changes of sides, bounds and costs, an
  infeasible solve and a solve with no change.

- `test_partial_solution` verifies that each `:MILPSolver` of the build
  writes, with a `SimpleConfiguration< std::vector< int > >` given to
  `get_var_solution()` and `get_dual_solution()`, only the part of the
  solution of the sub-Blocks it lists (and of the Block itself for -1),
  dynamic rows and columns included, that any other `Configuration` gives
  the whole solution, and that an index of no sub-Block throws.

- `test_qcp_dual` verifies the duals that each `:MILPSolver` of the build
  solving a convex problem with quadratic rows writes in them, one with and
  one without a linear part, and in the linear rows that follow them, on a
  small problem whose multipliers are known; the `:MILPSolver` are those
  that take quadratic rows (`LPPar-qcp.txt`), CPLEX and Gurobi have to solve
  it, SCIP is skipped if it gives no duals of them.

- `test_qp` verifies what each `:MILPSolver` of the build does with small
  quadratic models whose optimum is known: a convex QP with a separable and
  a non-separable objective, minimized and (its opposite) maximized, with
  rows of every kind and finite and infinite bounds, whose value, solution
  and duals are checked (the duals through the stationarity of the KKT
  conditions too), an MIQP, and a sequence of Modifications of the QP
  (coefficients of the objective, a side, a row added and removed, a
  variable fixed), each solved again and compared with a fresh Solver; the
  QP is also solved with `dblKKTTol` on, every Solver having to pass the
  check, and a sloppy `OSQPMILPSolver` having to be caught by it. A
  Solver that is not meant to handle a model has to refuse it with an
  exception or an error status rather than answer.

Run a tester directly (for instance `./test_dynamic seed [...]`, pass no
argument to see its usage), or run a batch over a set of sizes and seeds.

The `makefile` builds `test_cuts`, `test_dual`, `test_dynamic` and
`test_qp`, including the `MILPSolver` module and all its dependencies, and
the core SMS++ library.


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
