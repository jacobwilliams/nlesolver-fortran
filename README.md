# nlesolver-fortran

Nonlinear Equation Solver with Modern Fortran.

A basic Newton-Raphson type nonlinear equation solver for dense systems. A work in progress.

### Features

  * Works with square, under-determined, or over-determined systems.
  * Uses LAPACK routines (`dgesv` or `dgels`) to solve the linear system.
  * Has a Broyden update option.
  * Has various line search options.
  * Is object-oriented.

### License

  * BSD-3

### See also

  * [MINPACK](https://github.com/jacobwilliams/minpack)
