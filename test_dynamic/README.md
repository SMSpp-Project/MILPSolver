# /MILPSolver/test_dynamic

A tester which provides very comprehensive tests for any `CDASolver` 
able to handle Linear Programs (such as `MILPSolver` and its derived 
classes `CPXMILPSolver` , `SCIPMILPSolver` and `GRBMILPSolver`), as 
well as for some of the mechanics of the "core" SMS++ library.

This executable, given the input parameter n, constructs a "random"
Linear Program with n `ColVariable`, a "linear objective" 
(`FRealObjective` with a `LinearFunction` inside) and "linear constraints"
(`FRowConstraint` with a `LinearFunction` inside) and represent it in an 
`AbstractBlock` (LPBlock). Moreover, the built `ColVariable` can have simple
bound constraints imposed on them, which are implemented as ranged 
constraints if both rhs and lhs are finite.

Two appropriate `CDASolver` are attached to the LPBlock, which can be any
`Solver` capable of handling Linear Programs (say, some derived class 
of `MILPSolver` such as `CPXMILPSolver` , `SCIPMILPSolver` or 
`GRBMILPSolver`).

After all this is done, the LPBlock is solved with registered `Solvers` 
and the results (termination status and objective value, if applicable) 
are compared.

The LP is then repeatedly randomly modified, and re-solved several times; 
each time the results of the two `Solver` are compared.

The usage of the executable is the following:

       ./DynamicLP_test seed [wchg nvar dens #rounds #chng %chng]
       wchg: what to change, coded bit-wise [255]
             0 = add rows, 1 = delete rows 
             2 = modify rows, 3 = modify constants
             4 = change global lower/upper bound
             5 = add variables, 6 = delete variables
             7 = change variables bounds
       nvar: number of variables [10]
       dens: rows / variables [4]
       #rounds: how many iterations [40]
       #chng: number changes [10]
       %chng: probability of changing [0.5]

A batch file is provided that runs a not-so-large set of tests with
different sizes and seeds of the random generator; all these passing is a
good sign that no regressions have been done for the tested modules.

A makefile is also provided that builds the executable including the
`MILPSolver` module and all its dependencies (and, obviously, the 
core SMS++ library).


## Authors

- **Antonio Frangioni**  
  Dipartimento di Informatica  
  Università di Pisa

- **Enrico Calandrini**  
  Dipartimento di Informatica  
  Università di Pisa


## License

This code is provided free of charge under the [GNU Lesser General Public
License version 3.0](https://opensource.org/licenses/lgpl-3.0.html) -
see the [LICENSE](../LICENSE) file for details.
