![nlesolver-fortran](media/logo.png)
============

Nonlinear Equation Solver with Modern Fortran.

A basic Newton-Raphson type nonlinear equation solver for dense systems with `m` functions of `n` input variables.

A work in progress.

### Features

  * Is object-oriented.
  * Works with square, under-determined, or over-determined systems.
  * Uses LAPACK routines (`dgesv` or `dgels`) to solve the linear system.
     * if `n=m`,  uses `dgesv` (LU decomposition)
     * if `n/=m`, uses `dgels` (if `m>n` uses QR factorization,
       if `m<n` uses LQ factorization)
  * Has a Broyden update option.
  * Has various line search options.
     * 1 -- use a specified constant step size (0,1]
     * 2 -- backtracking linesearch method
     * 3 -- exact linesearch method using `fmin` minimizer
     * 4 -- evaluate function at specified fixed points

### License

  * BSD-3

### See also

  * [MINPACK](https://github.com/jacobwilliams/minpack)
