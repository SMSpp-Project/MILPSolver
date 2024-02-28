# /MILPSolver/test_dual

A tester which verifies the dual solution returned by `*MILPSolver`.
These tests are intended to ensure that the convention established by
`RowConstraint` regarding dual solutions is respected by the
`*MILPSolver`.

A makefile is also provided that builds the executable including the
`MILPSolver` module and all its dependencies (and, obviously, the 
core SMS++ library).


## Authors

- **Rafael Durbano Lobato**  
  Dipartimento di Informatica  
  Università di Pisa


## License

This code is provided free of charge under the [GNU Lesser General Public
License version 3.0](https://opensource.org/licenses/lgpl-3.0.html) -
see the [LICENSE](../LICENSE) file for details.
