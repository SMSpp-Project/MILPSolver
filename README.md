# MILPSolver
> A MILP Solver for SMS++ using IBM® ILOG® CPLEX® Optimization Studio

## Requirements

The project requires the [core SMS++](https://gitlab.com/frangio68/sms_plus_plus).

It has the same requirements of the core SMS++ project (Boost, Eigen3 and netCDF)
and also IBM® ILOG® CPLEX® Optimization Studio.

## Download and build

- Clone the project from the repository:

      $ git clone https://gitlab.com/niccolo/milpsolver

- Navigate inside the project directory:

      $ cd milpsolver

- Then configure and build the library:

      $ mkdir build
      $ cd build
      $ cmake ..
      $ make

- Optionally, you can install the library with:

      $ make install

- After the library is configured, you can use it in your CMake project with:

      find_package(MILPSolver)
      target_link_libraries(<my_target> SMS++::MILPSolver)

## Test

The test application takes a nc4 or dmx file containing a Min-Cost Flow problem
and solves it using both MILPSolver and MCFSolver, optionally changing the
problem multiple times.
By default, the test executable is:

    build/MILPSolver/test/MCFSolver_test

The usage is:

    MCFSolver_test <nc4 file> [seed how what #rounds #chng options]
           how: how to change, coded bit-wise 
                 0 = abstract (0) or physical (1) 
                 +2 = use ranged changes instead of sparse
           what: what to change, coded bit-wise 
                 0 = cost, 1 = cap, 2 = dfct, 3 = o.arc, 4 = c.arc
                 5 = add arc, 6 = delete arc
           options: bit 0 = re-optimize, other bits MCF-specific
                  Relax   : > 0 uses Auction
                  Cplex   : network pricing parameter
                  ZIB     : 1st bit == 1 ==> primal +
                            0 = Dantzig, 2 = First Eligible, 4 = MPP
                  Simplex : 1st bit == 1 ==> primal +
                            0 = Dantzig, 2 = First Eligible, 4 = MPP
