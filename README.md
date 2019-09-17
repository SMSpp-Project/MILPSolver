# MILPSolver
> A MILP Solver for SMS++ using IBM® ILOG® CPLEX® Optimization Studio

## Requirements

The project requires the core SMS++ and it is assumed to be in MILPSolver
directory besides SMS++.
It has the same requirements of the core SMS++ project (Boost, Eigen3 and netCDF)
and also IBM® ILOG® CPLEX® Optimization Studio.

## Build

You can use CMake to build the project, go
in the main SMS++ directory and type:

     mkdir build
     cd build
     cmake
     make

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

