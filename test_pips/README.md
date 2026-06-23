# test_pips

Small standalone driver for checking `PIPSMILPSolver` on one `.nc4` instance.

The test:

- deserializes an SMS++ `Block` from an `.nc4` file;
- applies a `BlockSolverConfig`;
- runs the first solver registered on the root block;
- prints the solver log and final status information to standard output;
- asks the solver to write the primal solution back into the SMS++ variables;
- dumps all `ColVariable` values in the block tree to a text file;
- asks the solver to write the dual solution back into the SMS++ constraints;
- dumps all `FRowConstraint` dual values in the block tree to a text file.

## Build

From this directory:

```bash
make release -j2
```

The makefile follows the structure used by `tests/TwoStageStochasticBlock`: it
includes SMS++, MILPSolver and UCBlock as normal SMS++ modules.

Plain `make` uses the default target from `../../tests/makefile_common`, which is
currently `debug`. Use `make release` when measuring timings.

## Run

General syntax:

```bash
./test_pips [instance.nc4] [BlockSolverConfig.txt] [primal.txt] [dual.txt]
```

Defaults:

```bash
./test_pips EC_CO_Test.nc4 BSCfg1-PIPS.txt primal_solution.txt dual_solution.txt
```

Example using the MUMPS-based PIPS configuration:

```bash
./test_pips EC_CO_Test.nc4 BSCfg1-PIPS-MUMPS.txt \
  primal_solution_mumps.txt dual_solution_mumps.txt \
> test_pips_mumps_run.log 2>&1
```

The `LD_LIBRARY_PATH` is needed because the executable is linked against the PIPS
PARDISO wrapper library, even when the selected PIPS options use another linear
solver.

## Solver Configurations

- `BSCfg1-PIPS.txt` uses the default PIPS linear solver selection. On the current
  local setup this selects Panua PARDISO, which requires a valid Panua license.
- `BSCfg1-PIPS-MA57.txt` forces MA57, but it requires `libhsl.so` to be available.
- `BSCfg1-PIPS-MUMPS.txt` forces MUMPS and disables `PARDISO_FOR_GLOBAL_SC`; this
  avoids the Panua license path and is useful for checking that the driver works.

## Output

The console output contains the solver log plus a final summary:

```text
Status
Time_s
Iterations
Lower_bound
Upper_bound
Objective
```

Both solution files are tab-separated. The primal solution file has the following columns:

```text
block_path  kind  group  index  address  value
```

The dual solution file has the following columns:

```text
block_path  kind  group  index  address  dual
```

`block_path` identifies the block in the nested SMS++ block tree. `kind` is
`static` or `dynamic`, and `group`/`index` identify the variable or constraint
location within that block.
