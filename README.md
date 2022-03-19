![nlesolver-fortran](media/logo.png)
============

Nonlinear Equation Solver with Modern Fortran.

![Build Status](https://github.com/jacobwilliams/nlesolver-fortran/actions/workflows/CI.yml/badge.svg)

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
     * use a specified constant step size (0,1]
     * backtracking linesearch method
     * exact linesearch method using `fmin` minimizer
     * evaluate function at specified fixed points

### Compiling

* A [Fortran Package Manager](https://github.com/fortran-lang/fpm) file is also included, so that the library and tests cases can be compiled with FPM. For example:

```
fpm build --profile release
fpm test --profile release
```

To use `nlesolver` within your fpm project, add the following to your `fpm.toml` file:
```toml
[dependencies]
fmin = { git="https://github.com/jacobwilliams/nlesolver-fortran.git" }

Note that LAPACK is required to build. The [fmin](https://github.com/jacobwilliams/fmin) library is also a dependency (which will be automatically downloaded by fpm).

### Documentation

 * The API documentation for the current ```master``` branch can be found [here](https://jacobwilliams.github.io/nlesolver-fortran/).  This is generated by processing the source files with [FORD](https://github.com/Fortran-FOSS-Programmers/ford).

### License

 * The NLESolver-Fortran source code and related files and documentation are distributed under a permissive free software [license](https://github.com/jacobwilliams/nlesolver-fortran/blob/master/LICENSE) (BSD-3).

### See also

  * [MINPACK](https://github.com/jacobwilliams/minpack)
