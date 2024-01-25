!******************************************************************************************************
!>
!  A basic multidimensional nonlinear equation solver, using a Newton-Raphson type direct method.
!
!### Features
!  * Works with square, under-determined, or over-determined systems.
!  * Uses LAPACK routines (`dgesv` or `dgels`) to solve the linear system.
!  * Has a Broyden update option.
!  * Has various line search options.
!
!### References
!  * https://en.wikipedia.org/wiki/Backtracking_line_search
!  * http://projecteuclid.org/download/pdf_1/euclid.pjm/1102995080
!  * http://help.agi.com/stk/index.htm#gator/eq-diffcorr.htm
!  * http://gmat.sourceforge.net/doc/R2015a/html/DifferentialCorrector.html
!
!### Author
!  * Jacob Williams
!
!### License
!  * BSD-3
!
!@todo add an `istat` output to func and grad, for user stopping
!      or to take a smaller stop (if istat>0 take a smaller step, if istat<0 abort)
!
!@note The default real kind (`wp`) can be
!      changed using optional preprocessor flags.
!      This library was built with real kind:
#ifdef REAL32
!      `real(kind=real32)` [4 bytes]
#elif REAL64
!      `real(kind=real64)` [8 bytes]
#elif REAL128
!      `real(kind=real128)` [16 bytes]
#else
!      `real(kind=real64)` [8 bytes]
#endif

    module nlesolver_module

    use iso_fortran_env
    use fmin_module,     only: fmin
    use lsqr_module,     only: lsqr_solver_ez
    use lusol_ez_module, only: solve, lusol_settings
    use lsmrModule,      only: lsmr_ez

    implicit none

    private

#ifdef REAL32
    integer,parameter,public :: nlesolver_rk = real32   !! real kind used by this module [4 bytes]
#elif REAL64
    integer,parameter,public :: nlesolver_rk = real64   !! real kind used by this module [8 bytes]
#elif REAL128
    integer,parameter,public :: nlesolver_rk = real128  !! real kind used by this module [16 bytes]
#else
    integer,parameter,public :: nlesolver_rk = real64   !! real kind used by this module [8 bytes]
#endif

    integer,parameter :: wp = nlesolver_rk  !! local copy of `nlesolver_rk` with a shorter name

    real(wp),parameter :: zero = 0.0_wp
    real(wp),parameter :: one  = 1.0_wp
    real(wp),parameter :: two  = 2.0_wp
    real(wp),parameter :: eps  = epsilon(one)  !! machine \( \epsilon \)
    real(wp),parameter :: big  = huge(one)

    !*********************************************************
        type,public :: nlesolver_type

        !! Nonlinear equations solver class.

        private

        integer     :: n                = 0             !! number of opt vars
        integer     :: m                = 0             !! number of constraints
        integer     :: max_iter         = 100           !! maximum number of iterations
        real(wp)    :: tol              = 1.0e-6_wp     !! convergence tolerance for function values
        real(wp)    :: alpha            = 1.0_wp        !! step length (when specified constant)
        real(wp)    :: alpha_min        = 0.1_wp        !! minimum step length (when allowed to vary)
        real(wp)    :: alpha_max        = 1.0_wp        !! maximum step length (when allowed to vary)
        real(wp)    :: tolx             = 1.0e-8_wp     !! convergence tolerance for `x`
        real(wp)    :: c                = 0.5_wp        !! backtracking linesearch parameter (0,1)
        real(wp)    :: tau              = 0.5_wp        !! backtracking linesearch parameter (0,1)
        real(wp)    :: fmin_tol         = 1.0e-5_wp     !! tolerance for "exact" linesearch
        integer     :: n_intervals      = 2             !! number of intervals for fixed point linesearch
        logical     :: use_broyden      = .false.       !! if true, a Broyden update is used
                                                        !! rather than computing the Jacobian
                                                        !! at every step. The `grad` function is
                                                        !! only called for the initial evaluation.
        integer     :: broyden_update_n = 4             !! if this value is `>0`, the Broyden update
                                                        !! is computed at most this many times before
                                                        !! the full Jacobian is recomputed.
        integer     :: n_uphill_max     = 5             !! maximum number of consecutive steps
                                                        !! to allow where the value of `f` increases
        logical     :: verbose          = .false.       !! verbose output printing
        integer     :: iunit            = output_unit   !! output unit for printing (assumed to be open).
        character(len=:),allocatable :: message         !! latest status message
        integer     :: istat            = -999          !! latest status message

        procedure(func_func),pointer       :: func             => null() !! user-supplied routine to compute the function
        procedure(export_func),pointer     :: export_iteration => null() !! user-supplied routine to export iterations
        procedure(wait_func),pointer       :: user_input_check => null() !! user-supplied routine to enable user to stop iterations
        procedure(linesearch_func),pointer :: linesearch       => null() !! line search method (determined by `step_mode` user input in [[nlesolver_type:initialize]])

        ! sparsity options:
        integer :: sparsity_mode = 1  !! sparsity mode:
                                      !!
                                      !! * 1 - assume dense (use dense solver)
                                      !! * 2 - assume sparse (use LSQR sparse solver).
                                      !! * 3 - assume sparse (use LUSOL sparse solver).
                                      !! * 4 - assume sparse (use LSMR sparse solver).
        integer :: n_nonzeros = -1 !! number of nonzero Jacobian elements (used for `sparsity_mode > 1`)
        integer,dimension(:),allocatable :: irow !! sparsity pattern nonzero elements row indices.
        integer,dimension(:),allocatable :: icol !! sparsity pattern nonzero elements column indices

        ! LSQR parameters:
        real(wp) :: atol    = zero  !! relative error in definition of `A`
        real(wp) :: btol    = zero  !! relative error in definition of `b`
        real(wp) :: conlim  = zero  !! An upper limit on `cond(Abar)`, the apparent
                                    !! condition number of the matrix `Abar`.
        integer  :: itnlim  = 100   !! max iterations
        integer  :: nout    = 0     !! output unit for printing
        real(wp) :: damp    = zero  !! damp parameter for LSQR

        ! LUSOL parameters:
        integer :: lusol_method = 0 !! * 0 => TPP: Threshold Partial   Pivoting.
                                    !! * 1 => TRP: Threshold Rook      Pivoting.
                                    !! * 2 => TCP: Threshold Complete  Pivoting.

        ! dense version:
        procedure(grad_func),pointer :: grad => null() !! user-supplied routine to compute the gradient of the function (dense version)

        ! sparse version:
        procedure(grad_func_sparse),pointer :: grad_sparse => null() !! user-supplied routine to compute the gradient of the function (sparse version)

        contains

        private

        procedure,public :: initialize => initialize_nlesolver_variables
        procedure,public :: solve      => nlesolver_solver
        procedure,public :: destroy    => destroy_nlesolver_variables
        procedure,public :: status     => get_status

        procedure :: set_status

        end type nlesolver_type
    !*********************************************************

    abstract interface

        subroutine func_func(me,x,f)
            !! compute the function
            import :: wp,nlesolver_type
            implicit none
            class(nlesolver_type),intent(inout) :: me
            real(wp),dimension(:),intent(in)    :: x
            real(wp),dimension(:),intent(out)   :: f
        end subroutine func_func

        subroutine grad_func(me,x,g)
            !! compute the gradient of the function (Jacobian). Dense version.
            import :: wp,nlesolver_type
            implicit none
            class(nlesolver_type),intent(inout) :: me
            real(wp),dimension(:),intent(in)    :: x
            real(wp),dimension(:,:),intent(out) :: g
        end subroutine grad_func

        subroutine grad_func_sparse(me,x,g)
            !! compute the gradient of the function (Jacobian). Sparse version.
            import :: wp,nlesolver_type
            implicit none
            class(nlesolver_type),intent(inout) :: me
            real(wp),dimension(:),intent(in)  :: x
            real(wp),dimension(:),intent(out) :: g !! sparse jacobian. length is `n_nonzeros`
        end subroutine grad_func_sparse

        subroutine export_func(me,x,f,iter)
            !! export an iteration:
            import :: wp,nlesolver_type
            implicit none
            class(nlesolver_type),intent(inout) :: me
            real(wp),dimension(:),intent(in)    :: x
            real(wp),dimension(:),intent(in)    :: f
            integer,intent(in)                  :: iter !! iteration number
        end subroutine export_func

        subroutine wait_func(me,user_stop)
            !! enable a user-triggered stop of the iterations:
            import :: wp,nlesolver_type
            implicit none
            class(nlesolver_type),intent(inout) :: me
            logical,intent(out)                 :: user_stop
        end subroutine wait_func

        subroutine linesearch_func(me,xold,p,x,f,fvec,fjac,fjac_sparse)
            !! line search method. Note that not all inputs/outputs are
            !! used by all methods.
            import :: wp,nlesolver_type
            implicit none
            class(nlesolver_type),intent(inout) :: me
            real(wp),dimension(me%n),intent(in) :: xold      !! previous value of `x`
            real(wp),dimension(me%n),intent(in) :: p         !! search direction
            real(wp),dimension(me%n),intent(out) :: x        !! new `x`
            real(wp),intent(inout) :: f                      !! * input: current magnitude of `fvec`,
                                                             !! * output: new value of `f`
            real(wp),dimension(me%m),intent(inout) :: fvec   !! * input: current function vector,
                                                             !! * output: new function vector
            real(wp),dimension(:,:),intent(in),optional :: fjac !! jacobian matrix [dense]
            real(wp),dimension(:),intent(in),optional :: fjac_sparse !! jacobian matrix [sparse]
        end subroutine linesearch_func

    end interface

    contains
!*******************************************************************************************************

!*****************************************************************************************
!>
!  Set status flag and message.

    subroutine set_status(me,istat,string,i,r)

    implicit none

    class(nlesolver_type),intent(inout) :: me
    integer,intent(in)                  :: istat   !! status code
    character(len=*),intent(in)         :: string  !! status message
    integer,intent(in),optional         :: i       !! an integer value to append
    real(wp),intent(in),optional        :: r       !! a real value to append

    character(len=256) :: numstr !! for number fo string conversion
    character(len=:),allocatable :: message !! the full message to log
    integer :: iostat !! write `iostat` code

    message = trim(string)
    if (present(i)) then
        numstr = ''
        write(numstr,fmt=*,iostat=iostat) i
        if (iostat/=0) numstr = '****'
        message = message//' '//trim(adjustl(numstr))
    end if
    if (present(r)) then
        numstr = ''
        write(numstr,fmt=*,iostat=iostat) r
        if (iostat/=0) numstr = '****'
        message = message//' '//trim(adjustl(numstr))
    end if

    if (me%verbose) then
        write(me%iunit,'(A)',iostat=iostat) message
    end if

    ! store in the class:
    me%istat = istat
    me%message = message

    end subroutine set_status
!*****************************************************************************************

!*****************************************************************************************
!>
!  Return the status code and message from the [[nlesolver_type]] class.
!
!### Status codes
!
!  * -1   -- Error: Invalid alpha
!  * -2   -- Error: Invalid alpha_min
!  * -3   -- Error: Invalid alpha_max
!  * -4   -- Error: Alpha_min must be < alpha_max
!  * -5   -- Error: Invalid step_mode
!  * -6   -- Error solving linear system
!  * -7   -- Error: More than 5 steps in the uphill direction
!  * -8   -- Error: Divide by zero when computing Broyden update
!  * -9   -- Error: Out of memory
!  * -10  -- Error: function routine is not associated
!  * -11  -- Error: gradient routine is not associated
!  * -12  -- Error: backtracking linesearch c must be in range (0, 1)
!  * -13  -- Error: backtracking linesearch tau must be in range (0, 1)
!  * -14  -- Error: must specify grad_sparse, irow, and icol for sparsity_mode > 1
!  * -15  -- Error: irow and icol must be the same length
!  * -999 -- Error: class has not been initialized
!  * 0    -- Class successfully initialized in [[nlesolver_type:initialize]]
!  * 1    -- Required accuracy achieved
!  * 2    -- Solution cannot be improved
!  * 3    -- Maximum number of iterations reached
!  * 4    -- Stopped by the user

    subroutine get_status(me, istat, message)

    implicit none

    class(nlesolver_type),intent(inout) :: me
    integer,intent(out),optional :: istat  !! Integer status code.
    character(len=:),allocatable,intent(out),optional :: message  !! Text status message

    if (present(istat)) istat = me%istat

    if (present(message)) then
        if (allocated(me%message)) then
            message = trim(me%message)
        else
            message = 'Error: class has not been initialized'
        end if
    end if

    end subroutine get_status
!*****************************************************************************************

!*****************************************************************************************
!>
!  Constructor for the class.

    subroutine initialize_nlesolver_variables(me,&
                    n,m,max_iter,tol,alpha,alpha_min,alpha_max,tolx,fmin_tol,&
                    backtrack_c,backtrack_tau,&
                    use_broyden,broyden_update_n,step_mode,&
                    func,grad,grad_sparse,&
                    export_iteration,user_input_check,&
                    verbose,iunit,n_uphill_max,n_intervals,&
                    sparsity_mode,irow,icol,&
                    atol,btol,conlim,damp,itnlim,nout,&
                    lusol_method )

    implicit none

    class(nlesolver_type),intent(inout) :: me
    integer,intent(in)                :: n                  !! number of optimization variables
    integer,intent(in)                :: m                  !! number of constraints
    integer,intent(in)                :: max_iter           !! maximum number of iterations
    real(wp),intent(in)               :: tol                !! function convergence tolerance
    procedure(func_func)              :: func               !! computes the function vector
    procedure(grad_func),optional        :: grad            !! computes the jacobian [required for dense mode: `sparsity_mode=1`]
    procedure(grad_func_sparse),optional :: grad_sparse     !! computes the jacobian [required for sparse mode: `sparsity_mode>1`]
    integer,intent(in),optional       :: step_mode          !! step mode:
                                                            !!
                                                            !!  * 1 = use the specified `alpha` (0,1]
                                                            !!  * 2 = backtracking linesearch (between `alpha_min` and `alpha_max`)
                                                            !!  * 3 = exact linesearch (between `alpha_min` and `alpha_max`)
                                                            !!  * 4 = evaluate function at specified fixed points  (between `alpha_min` and `alpha_max`)
    real(wp),intent(in),optional      :: alpha              !! constant step length for `step_mode=1` (0,1]
    real(wp),intent(in),optional      :: alpha_min          !! minimum step length (0,1]
    real(wp),intent(in),optional      :: alpha_max          !! maximum step length (0,1]
    real(wp),intent(in),optional      :: tolx               !! convergence tolerance for changes in `x`
    real(wp),intent(in),optional      :: fmin_tol           !! convergence tolerance for [[fmin]] (used when `step_mode=3`)
    real(wp),intent(in),optional      :: backtrack_c        !! backtracking linesearch parameter (0,1)
    real(wp),intent(in),optional      :: backtrack_tau      !! backtracking linesearch parameter (0,1)
    logical,intent(in),optional       :: use_broyden        !! use a Broyden update (default is False)
    integer,intent(in),optional       :: broyden_update_n   !! For Broyden mode, update the full Jacobian
                                                            !! at most every this many iterations (must be >1)
                                                            !! If <=1 then Jacobian is only computed on the
                                                            !! first iteration.
    procedure(export_func),optional   :: export_iteration   !! function to export each iteration
    procedure(wait_func),optional     :: user_input_check   !! check for user input (to stop solver if necessary)
    logical,intent(in),optional       :: verbose            !! for verbose status printing
    integer,intent(in),optional       :: iunit              !! unit for verbose printing (assumed to be open).
                                                            !! by default this is `output_unit`.
    integer,intent(in),optional       :: n_uphill_max       !! maximum number of consecutive steps
                                                            !! to allow where the value of `f` increases
    integer,intent(in),optional       :: n_intervals        !! number of intervals for fixed point linesearch

    integer,intent(in),optional :: sparsity_mode !! sparsity mode:
                                                 !!
                                                 !! * 1 - assume dense (use dense solver)
                                                 !!   Must specify `grad` for this mode.
                                                 !! * 2 - assume sparse (use LSQR sparse solver).
                                                 !! * 3 - assume sparse (use LUSOL sparse solver).
                                                 !! * 4 - assume sparse (use LSMR sparse solver).
                                                 !! Must also specify `grad_sparse` and the `irow` and `icol`
                                                 !! sparsity pattern for `mode>1`.
    integer,dimension(:),intent(in),optional :: irow !! sparsity pattern nonzero elements row indices.
                                                     !! must be specified with `icol` and be the same length (`n_nonzeros`).
    integer,dimension(:),intent(in),optional :: icol !! sparsity pattern nonzero elements column indices
                                                     !! must be specified with `icol` and be the same length (`n_nonzeros`).
    real(wp),intent(in),optional      :: atol    !! `LSQR`: relative error in definition of `A`
    real(wp),intent(in),optional      :: btol    !! `LSQR`: relative error in definition of `b`
    real(wp),intent(in),optional      :: conlim  !! `LSQR`: An upper limit on `cond(Abar)`, the apparent
                                                 !! condition number of the matrix `Abar`.
    real(wp),intent(in),optional      :: damp    !! `LSQR`: damp factor
    integer,intent(in),optional       :: itnlim  !! `LSQR`: max iterations
    integer,intent(in),optional       :: nout    !! `LSQR`: output unit for printing
    integer,intent(in),optional       :: lusol_method !! `LUSOL` method:
                                                      !!
                                                      !! * 0 => TPP: Threshold Partial   Pivoting.
                                                      !! * 1 => TRP: Threshold Rook      Pivoting.
                                                      !! * 2 => TCP: Threshold Complete  Pivoting.

    logical :: status_ok !! true if there were no errors

    status_ok = .true.

    !required:
    me%n        = abs(n)
    me%m        = abs(m)
    me%max_iter = abs(max_iter)
    me%tol      = abs(tol)
    me%func     => func
    if (present(grad)) me%grad               => grad
    if (present(grad_sparse)) me%grad_sparse => grad_sparse

    !optional:

    if (present(step_mode)) then
        select case (step_mode)
        case(1) ! = use the specified `alpha` (0,1]
            me%linesearch => simple_step
        case(2) ! = backtracking linesearch (between `alpha_min` and `alpha_max`)
            me%linesearch => backtracking_linesearch
        case(3) ! = exact linesearch (between `alpha_min` and `alpha_max`)
            me%linesearch => exact_linesearch
        case(4) ! = evaluate function at specified fixed points (between `alpha_min` and `alpha_max`)
            me%linesearch => fixed_point_linesearch
        case default
            status_ok = .false.
            call me%set_status(istat = -5, string = 'Error: invalid step_mode:',i=step_mode)
            return
        end select
    else
        me%linesearch => simple_step
    end if

    if (present(alpha))             me%alpha             = abs(alpha)
    if (present(alpha_min))         me%alpha_min         = abs(alpha_min)
    if (present(alpha_max))         me%alpha_max         = abs(alpha_max)
    if (present(tolx))              me%tolx              = abs(tolx)
    if (present(backtrack_c))       me%c                 = abs(backtrack_c)
    if (present(backtrack_tau))     me%tau               = abs(backtrack_tau)
    if (present(use_broyden))       me%use_broyden       = use_broyden
    if (present(broyden_update_n))  me%broyden_update_n  = abs(broyden_update_n)
    if (present(verbose))           me%verbose           = verbose
    if (present(iunit))             me%iunit             = iunit
    if (present(n_uphill_max))      me%n_uphill_max      = abs(n_uphill_max)
    if (present(n_intervals))       me%n_intervals       = max(abs(n_intervals),1)
    if (present(fmin_tol))          me%fmin_tol          = abs(fmin_tol)

    if (present(export_iteration)) me%export_iteration  => export_iteration
    if (present(user_input_check)) me%user_input_check  => user_input_check

    ! error checks:
    if (me%alpha<zero .or. me%alpha>one) then
        status_ok = .false.
        call me%set_status(istat = -1, string = 'Error: invalid alpha:',r=me%alpha)
        return
    end if
    if (me%alpha_min<zero .or. me%alpha_min>one) then
        status_ok = .false.
        call me%set_status(istat = -2, string = 'Error: invalid alpha_min:',r=me%alpha_min)
        return
    end if
    if (me%alpha_max<zero .or. me%alpha_max>one) then
        status_ok = .false.
        call me%set_status(istat = -3, string = 'Error: invalid alpha_max:',r=me%alpha_max)
        return
    end if
    if (me%alpha_max<=me%alpha_min) then
        status_ok = .false.
        call me%set_status(istat = -4, string = 'Error: alpha_min must be < alpha_max')
        return
    end if
    if (me%c<zero .or. me%c>one) then
        status_ok = .false.
        call me%set_status(istat = -12, string = 'Error: backtracking linesearch c must be in range (0, 1):',r=me%c)
        return
    end if
    if (me%tau<zero .or. me%tau>one) then
        status_ok = .false.
        call me%set_status(istat = -13, string = 'Error: backtracking linesearch tau must be in range (0, 1):',r=me%tau)
        return
    end if

    ! initialize:
    me%n_nonzeros = -1 ! not used
    if (allocated(me%irow)) deallocate(me%irow)
    if (allocated(me%icol)) deallocate(me%icol)

    if (present(sparsity_mode)) then
        me%sparsity_mode = sparsity_mode
        if (sparsity_mode>1) then ! sparse solver method
            if (present(irow) .and. present(icol) .and. present(grad_sparse)) then
                if (size(irow) == size(icol)) then
                    me%n_nonzeros = size(irow)
                    me%irow = irow
                    me%icol = icol
                else
                    call me%set_status(istat = -15, string = 'Error: irow and icol must be the same length')
                    return
                end if
            else
                call me%set_status(istat = -14, string = 'Error: must specify grad_sparse, irow, and icol for sparsity_mode > 1')
                return
            end if
            ! LSQR optional inputs:
            if (present(atol))   me%atol   = atol
            if (present(btol))   me%btol   = btol
            if (present(conlim)) me%conlim = conlim
            if (present(damp))   me%damp   = damp
            if (present(itnlim)) me%itnlim = itnlim
            if (present(nout))   me%nout   = nout

            ! LUSOL method
            if (present(nout)) me%lusol_method = lusol_method

        end if
    end if

    if (status_ok) then
        call me%set_status(istat = 0, string = 'Class successfully initialized')
    end if

    end subroutine initialize_nlesolver_variables
!*****************************************************************************************

!*****************************************************************************************
!>
!  Main solver.

    subroutine nlesolver_solver(me,x)

    implicit none

    class(nlesolver_type),intent(inout) :: me
    real(wp),dimension(:),intent(inout) :: x

    real(wp),dimension(:)  ,allocatable :: fvec             !! function vector
    real(wp),dimension(:,:),allocatable :: fjac             !! jacobian matrix [dense]
    real(wp),dimension(:),  allocatable :: fjac_sparse      !! jacobian matrix [sparse]
    real(wp),dimension(:),  allocatable :: prev_fjac_sparse !! previous jacobian matrix [sparse]
    real(wp),dimension(:)  ,allocatable :: rhs              !! linear system right-hand side
    real(wp),dimension(:)  ,allocatable :: p                !! search direction
    real(wp),dimension(:)  ,allocatable :: xold             !! previous value of `x`
    real(wp),dimension(:)  ,allocatable :: prev_fvec        !! previous function vector
    real(wp),dimension(:,:),allocatable :: prev_fjac        !! previous jacobian matrix
    real(wp),dimension(:,:),allocatable :: delf             !! used for Broyden (rank 2 for `matmul`)
    real(wp),dimension(:,:),allocatable :: delx             !! used for Broyden (rank 2 for `matmul`)
    logical                             :: user_stop        !! user stop button flag
    integer                             :: info             !! status flag from the [[linear_solver]]
    integer                             :: iter             !! iteration counter
    real(wp)                            :: f                !! magnitude of `fvec`
    real(wp)                            :: fold             !! previous value of `f`
    integer                             :: n_uphill         !! number of steps taken in the "uphill" direction
                                                            !! (where `f` is increasing)
    real(wp)                            :: delxmag2         !! used for Broyden
    logical                             :: recompute_jac    !! if using Broyden, and we want to call the user
                                                            !! jacobian routine instead
    integer                             :: broyden_counter  !! number of times the broyden update has been used
    integer                             :: alloc_stat       !! allocation status flag
    type(lsqr_solver_ez) :: sparse_solver  !! sparse LSQR solver class
    type(lusol_settings) :: lusol_options
    integer :: i !! counter
    integer,dimension(:),allocatable :: idx, index_array !! for sparse indexing
    character(len=10) :: i_str !! string version of `i` for row string
    real(wp) :: normA, condA, normr, normAr, normx !! for LSMR
    integer :: localSize, itn !! for LSMR

    if (me%istat<0) return ! class was not initialized properly

    if (.not. associated(me%func)) then
        call me%set_status(istat = -10, string = 'Error: function routine is not associated')
        return
    end if
    if (me%sparsity_mode==1 .and. .not. associated(me%grad)) then
        call me%set_status(istat = -11, string = 'Error: gradient routine is not associated')
        return
    end if
    if (me%sparsity_mode>1 .and. .not. associated(me%grad_sparse)) then
        call me%set_status(istat = -11, string = 'Error: gradient routine is not associated')
        return
    end if

    ! initialize:
    iter            = 0
    n_uphill        = 0
    recompute_jac   = .false.
    broyden_counter = 0

    ! allocate the arrays:
    alloc_stat = 0
    if (alloc_stat==0) allocate(fvec     (me%m)      , stat=alloc_stat)

    if (me%sparsity_mode==1) then
        ! dense
        if (alloc_stat==0) allocate(fjac (me%m,me%n) , stat=alloc_stat)
    else
        ! sparse
        if (alloc_stat==0) allocate(fjac_sparse (me%n_nonzeros) , stat=alloc_stat)
        if (me%use_broyden .and. alloc_stat==0) allocate(prev_fjac_sparse(me%n_nonzeros) , stat=alloc_stat)
    end if

    if (alloc_stat==0) allocate(rhs      (me%m)      , stat=alloc_stat)
    if (alloc_stat==0) allocate(p        (me%n)      , stat=alloc_stat)
    if (alloc_stat==0) allocate(xold     (me%n)      , stat=alloc_stat)
    if (me%use_broyden) then
        ! only need these for broyden:
        if (alloc_stat==0) allocate(prev_fvec(me%m)      , stat=alloc_stat)
        if (me%sparsity_mode/=1 .and. alloc_stat==0) then
            allocate(prev_fjac(me%m,me%n) , stat=alloc_stat)
            index_array = [(i, i=1,me%n_nonzeros)] ! an array to index the sparse jacobian elements
        end if
        if (alloc_stat==0) allocate(delf     (me%m,1)    , stat=alloc_stat)
        if (alloc_stat==0) allocate(delx     (me%n,1)    , stat=alloc_stat)
    end if
    if (alloc_stat/=0) then
        call me%set_status(istat = -9, string = 'Error: Out of memory')
        return
    else
        me%istat = -998
        me%message = 'Unknown error'
    end if

    ! evaluate the function:
    call me%func(x,fvec)
    f = norm2(fvec)

    ! check to see if initial guess is a root:
    if (f <= me%tol) then

        call me%set_status(istat = 1, string = 'Required accuracy achieved')

    else

        ! main iteration loop:
        do iter = 1,me%max_iter

            ! Export the current iteration:
            if (associated(me%export_iteration)) call me%export_iteration(x,fvec,iter)

            ! Check for user stop:
            if (associated(me%user_input_check)) then
                call me%user_input_check(user_stop)
                if (user_stop) then
                    call me%set_status(istat = 4, string = 'Stopped by the user')
                    exit
                end if
            end if

            if (me%use_broyden .and. .not. recompute_jac) then
                if (iter==1) then
                    ! always compute Jacobian on the first iteration
                    !call me%grad(x,fjac)
                    select case (me%sparsity_mode)
                    case (1)
                        call me%grad(x,fjac)
                    case (2:)
                        call me%grad_sparse(x,fjac_sparse)
                    end select
                    broyden_counter = 0
                else
                    ! and use Broyden update to estimate Jacobian
                    ! for subsequent iterations.

                    ! note: fvec was computed the last iteration
                    delx(:,1) = x - xold
                    delf(:,1) = fvec - prev_fvec

                    if (me%sparsity_mode==1) then ! dense

                        delxmag2 = dot_product(delx(:,1),delx(:,1))
                        if (delxmag2 < eps) then
                            call me%set_status(istat = -8, &
                                    string = 'Error: Divide by zero when computing Broyden update')
                            exit
                        end if

                        ! Jacobian estimate:
                        ! This is the equation from wikipedia : https://en.wikipedia.org/wiki/Broyden%27s_method
                        fjac = prev_fjac + matmul((delf-matmul(prev_fjac,delx)), transpose(delx)) / delxmag2

                        ! see also: eqn 4 in Schubert, which is different ...

                    else ! using a sparse option

                        ! Just a row-by-row version of the fjac equation above.
                        ! see L.K. Schubert, 1970
                        do i = 1, me%m

                            idx = pack(index_array, mask = me%irow==i) ! the nonzero elements in this row
                            if (size(idx)==0) cycle ! no non-zeros in this row

                            associate (dx => delx(me%icol(idx),:)) ! nonzero x vec for this row
                                delxmag2 = dot_product(dx(:,1),dx(:,1)) ! only those x's for this row
                                if (delxmag2 < eps) then
                                    write(i_str,'(I10)') i
                                    call me%set_status(istat = -8, &
                                            string = 'Error: Divide by zero when computing sparse Broyden update for row '//&
                                            trim(adjustl(i_str)))
                                    exit
                                end if
                                fjac_sparse(idx) = prev_fjac_sparse(idx) + &
                                        matmul((delf(i,1)-matmul(prev_fjac_sparse(idx),dx)), transpose(dx)) / delxmag2
                            end associate

                        end do

                    end if
                    broyden_counter = broyden_counter + 1
                end if
                prev_fvec = fvec
                if (me%sparsity_mode==1) then
                    prev_fjac = fjac
                else
                    prev_fjac_sparse = fjac_sparse
                end if
            else
                ! compute the jacobian:
                select case (me%sparsity_mode)
                case (1)
                    call me%grad(x,fjac)
                case (2:)
                    call me%grad_sparse(x,fjac_sparse)
                end select
                recompute_jac = .false.  ! for broyden
                broyden_counter = 0
            end if

            xold = x
            fold = f

            ! compute the search direction p by solving linear system:
            rhs = -fvec ! RHS of the linear system
            select case (me%sparsity_mode)
            case (1)
                ! use dense solver
                call linear_solver(me%m,me%n,fjac,rhs,p,info)
            case (2)
                ! initialize the LSQR sparse solver
                call sparse_solver%initialize(me%m,me%n,fjac_sparse,me%irow,me%icol,&
                                              atol   = me%atol, &
                                              btol   = me%btol, &
                                              conlim = me%conlim, &
                                              itnlim = me%itnlim, &
                                              nout   = me%nout)
                call sparse_solver%solve(rhs,me%damp,p,info) ! solve the linear system
                ! check convergence:
                select case (info)
                case(4)
                    call me%set_status(istat = -1004, &
                                string = 'LSQR Error: The system appears to be ill-conditioned. istop =', i=info)
                    exit
                case(5)
                    call me%set_status(istat = -1005, &
                                string = 'LSQR Error: The iteration limit was reached. istop =', i=info)
                    exit
                case default
                    info = 0
                end select
            case (3)
                ! use lusol solver
                lusol_options%method = me%lusol_method
                call solve(me%n,me%m,me%n_nonzeros,me%irow,me%icol,fjac_sparse,rhs,p,info,&
                           settings=lusol_options)
            case (4)

                ! use LSMR solver

                ! TODO this should be an input
                localSize = 0
                !localSize = min(me%m, me%n)

                me%conlim = 1.0_wp/(100*epsilon(1.0_wp))
                call lsmr_ez  ( me%m, me%n, me%irow, me%icol, fjac_sparse, rhs, me%damp, &
                                me%atol, me%btol, me%conlim, me%itnlim, localSize, me%nout, &
                                p, info, itn, normA, condA, normr, normAr, normx )

                ! check convergence:
                select case (info)
                case(4)
                    call me%set_status(istat = -1004, &
                                string = 'LSMR Error: The system appears to be ill-conditioned. istop =', i=info)
                    exit
                case(5)
                    call me%set_status(istat = -1005, &
                                string = 'LSMR Error: The iteration limit was reached. istop =', i=info)
                    exit
                case default
                    info = 0
                end select

            end select

            ! check for errors:
            if (info /= 0) then

                call me%set_status(istat = -6, string = 'Error solving linear system. info =', i=info)
                exit

            else

                ! next step, using the specified method:
                call me%linesearch(xold,p,x,f,fvec,fjac,fjac_sparse)

                ! keep track of the number of steps in the "uphill" direction:
                if (f>fold) then
                    n_uphill = n_uphill + 1
                else
                    n_uphill = 0
                end if

                ! check for stopping conditions
                if (f <= me%tol) then

                    call me%set_status(istat = 1, string = 'Required accuracy achieved')
                    exit

                elseif ( maxval(abs(x-xold)) <= me%tolx ) then

                    call me%set_status(istat = 2, string = 'Solution cannot be improved')
                    exit

                elseif (iter == me%max_iter) then

                    call me%set_status(istat = 3, string = 'Maximum number of iterations reached')
                    exit

                elseif (n_uphill > me%n_uphill_max) then

                    call me%set_status(istat = 5, string = 'Too many steps in the uphill direction')
                    exit

                elseif (me%use_broyden) then

                    ! If delxmag2 is too small when using broyden, just
                    ! call the user-supplied jacobian function to avoid
                    ! a divide by zero on the next step. This should
                    ! normally only happen when the solution is almost converged.
                    if (norm2(x-xold)**2 <= eps) then
                        recompute_jac = .true.
                    else
                        if (me%broyden_update_n>0) then
                            ! Note that we also recompute if we have taken an uphill step
                            if (broyden_counter==me%broyden_update_n .or. n_uphill>0) then
                                ! time to recompute the full jacobian
                                recompute_jac = .true.
                            end if
                        end if
                    end if

                endif

            end if

        end do   !end of iterations loop

    end if

    !Export the last iteration:
    iter = iter + 1
    if (associated(me%export_iteration)) call me%export_iteration(x,fvec,iter)

    end subroutine nlesolver_solver
!*****************************************************************************************

!*****************************************************************************************
!>
!   Destructor

    subroutine destroy_nlesolver_variables(me)

    implicit none

    class(nlesolver_type),intent(out) :: me

    me%message = 'Error: class has not been initialized'
    me%istat   = -999

    end subroutine destroy_nlesolver_variables
!*****************************************************************************************

!*****************************************************************************************
!>
!  Solve the linear system: \( Ax = b \), using a dense, direct method.
!
!  * if `n=m` : use LAPACK `dgesv` (LU decomposition)
!  * if `n/=m` : use LAPACK `dgels` (if m>n uses QR factorization,
!    if m<n uses LQ factorization)

    subroutine linear_solver(m,n,a,b,x,info)

    implicit none

    integer,intent(in)                 :: n          !! number of columns in `a`
    integer,intent(in)                 :: m          !! number of rows in `a`
    real(wp),dimension(m,n),intent(in) :: a          !! `A` matrix of the linear system
    real(wp),dimension(m),intent(in)   :: b          !! RHS of the linear system
    real(wp),dimension(n),intent(out)  :: x          !! the solution of the linear system.
    integer,intent(out)                :: info       !! output status flag (`=0` if success)

    ! LAPACK routine interfaces:
    interface
        subroutine dgesv(n,nrhs,a,lda,ipiv,b,ldb,info)
            !! See: [?gesv](https://www.intel.com/content/www/us/en/develop/documentation/onemkl-developer-reference-fortran/top/lapack-routines/lapack-linear-equation-routines/lapack-linear-equation-driver-routines/gesv.html)
            import :: wp
            implicit none
            integer :: info
            integer :: lda
            integer :: ldb
            integer :: n
            integer :: nrhs
            integer :: ipiv(*)
            real(wp) :: a(lda,*)
            real(wp) :: b(ldb,*)
        end subroutine dgesv
        subroutine dgels(trans,m,n,nrhs,a,lda,b,ldb,work,lwork,info)
            !! See: [?gels](https://www.intel.com/content/www/us/en/develop/documentation/onemkl-developer-reference-fortran/top/lapack-routines/lapack-least-squares-and-eigenvalue-problem/lapack-least-squares-eigenvalue-problem-driver/linear-least-squares-lls-problems-lapack-driver/gels.html)
            import :: wp
            implicit none
            character :: trans
            integer :: info
            integer :: lda
            integer :: ldb
            integer :: lwork
            integer :: m
            integer :: n
            integer :: nrhs
            real(wp) :: a(lda,*)
            real(wp) :: b(ldb,*)
            real(wp) :: work(*)
        end subroutine dgels
    end interface

    integer,dimension(:),allocatable    :: ipiv   !! pivot indices array
    real(wp),dimension(:,:),allocatable :: bmat   !! copy of `b` so it won't be overwritten
    real(wp),dimension(:),allocatable   :: work   !! work array for `dgels`
    real(wp),dimension(:,:),allocatable :: amat   !! copy of `a` so it won't be overwritten
    integer                             :: lwork  !! size of `work`

    allocate(amat(m,n))
    allocate(bmat(max(1,m,n),1))

    if (n==m) then  !normal inverse

        allocate(ipiv(n))

        amat = a
        bmat(1:n,1) = b
        call dgesv(n, 1, amat, n, ipiv, bmat, n, info)
        x = bmat(1:n,1)

    else

        amat = a
        bmat = zero
        bmat(1:m,1) = b
        lwork = min(m,n)+max(1,m,n)
        allocate(work(lwork))
        call dgels('N', m, n, 1, amat, m, bmat, max(1,m,n), work, lwork, info)
        x = bmat(1:n,1)

    end if

    end subroutine linear_solver
!*****************************************************************************************

!*****************************************************************************************
!>
!  Take a simple step in the search direction of `p * alpha`.

    subroutine simple_step(me,xold,p,x,f,fvec,fjac,fjac_sparse)

    implicit none

    class(nlesolver_type),intent(inout) :: me
    real(wp),dimension(me%n),intent(in) :: xold      !! previous value of `x`
    real(wp),dimension(me%n),intent(in) :: p         !! search direction
    real(wp),dimension(me%n),intent(out) :: x        !! new `x`
    real(wp),intent(inout) :: f                      !! magnitude of `fvec`
    real(wp),dimension(me%m),intent(inout) :: fvec   !! function vector
    real(wp),dimension(:,:),intent(in),optional :: fjac !! jacobian matrix [dense]
    real(wp),dimension(:),intent(in),optional :: fjac_sparse !! jacobian matrix [sparse]

    x = xold + p * me%alpha

    !evaluate the function at the new point:
    call me%func(x,fvec)
    f = norm2(fvec)

    end subroutine simple_step
!*****************************************************************************************

!*****************************************************************************************
!>
!  Backtracking line search.
!
!### See also
!  * [Backtracking line search](https://en.wikipedia.org/wiki/Backtracking_line_search)
!
!@note Either `fjac` or `fjac_sparse` should be present.

    subroutine backtracking_linesearch(me,xold,p,x,f,fvec,fjac,fjac_sparse)

    implicit none

    class(nlesolver_type),intent(inout) :: me
    real(wp),dimension(me%n),intent(in) :: xold      !! previous value of `x`
    real(wp),dimension(me%n),intent(in) :: p         !! search direction
    real(wp),dimension(me%n),intent(out) :: x        !! new `x`
    real(wp),intent(inout) :: f                      !! magnitude of `fvec`
    real(wp),dimension(me%m),intent(inout) :: fvec   !! function vector
    real(wp),dimension(:,:),intent(in),optional :: fjac !! jacobian matrix [dense]
    real(wp),dimension(:),intent(in),optional :: fjac_sparse !! jacobian matrix [sparse]

    integer    :: i                  !! counter
    real(wp)   :: slope              !! local slope of the function of `alpha` along the search direction used for line search
    logical    :: min_alpha_reached  !! if the minimum step size is reached during the line search
    real(wp)   :: alpha              !! `alpha` for the line search
    real(wp)   :: ftmp               !! `f` value for linesearch
    real(wp)   :: t                  !! used for line search
    real(wp),dimension(:),allocatable :: gradf      !! line search objective function gradient vector
    real(wp),dimension(:),allocatable :: xtmp       !! `x` value for linesearch
    real(wp),dimension(:),allocatable :: fvectmp    !! `fvec` value for linesearch

    ! allocate arrays:
    allocate(gradf(me%n))
    allocate(xtmp(me%n))
    allocate(fvectmp(me%m))

    ! compute the gradient of the function to be minimized
    ! (which in this case is 1/2 the norm of fvec). Use the chain
    ! rule and the Jacobian matrix already computed.
    if (present(fjac)) then
        ! dense
        do i=1,me%n
            gradf(i) = dot_product(fvec,fjac(:,i))
        end do
    else
        ! sparse
        do i=1,me%n
            gradf(i) = dot_product(fvec,pack(fjac_sparse,mask=me%icol==i))
        end do
    end if
    slope = dot_product(p, gradf)
    t = -me%c * slope

    if (me%verbose) then
        write(me%iunit,'(1P,*(A,1X,E16.6))') '        slope    = ', slope
        write(me%iunit,'(1P,*(A,1X,E16.6))') '        t        = ', t
    end if

    ! perform the line search:
    min_alpha_reached = .false.
    alpha = me%alpha_max  ! start with the largest step
    do

        xtmp = xold + p * alpha
        call me%func(xtmp,fvectmp)
        ftmp = norm2(fvectmp)

        if (me%verbose) then
            write(me%iunit,'(1P,*(A,1X,E16.6))')          '        alpha    = ', alpha,    ' f       = ', ftmp
            if (f - ftmp >= alpha*t) then
                write(me%iunit,'(1P,2(A,1X,E16.6),1X,A)') '        f - ftmp = ', f - ftmp, ' alpha*t = ', alpha*t, ' [ACCEPTED]'
            else
                write(me%iunit,'(1P,*(A,1X,E16.6))')      '        f - ftmp = ', f - ftmp, ' alpha*t = ', alpha*t
            end if
        end if

        if (((f - ftmp) / 2.0_wp >= alpha*t) .or. min_alpha_reached) then
            if (me%verbose .and. min_alpha_reached) then
                write(me%iunit,'(A)') '        Minimum alpha reached'
            end if
            ! Armijo-Goldstein condition is satisfied
            ! (or the min step has been reached)
            x    = xtmp
            fvec = fvectmp
            f    = ftmp
            exit
        end if
        alpha = alpha * me%tau ! reduce step size

        if (alpha<=me%alpha_min) then
            alpha = me%alpha_min
            min_alpha_reached = .true. ! will stop on the next step
        end if

    end do

    end subroutine backtracking_linesearch
!*****************************************************************************************

!*****************************************************************************************
!>
!  An exact linesearch that uses a derivative-free minimizer to
!  find the minimum value of `f(x)` between
!  `x = xold + p * alpha_min` and
!  `x = xold + p * alpha_max`.
!
!  Usually this is overkill and not necessary, but is here as an option for testing.

    subroutine exact_linesearch(me,xold,p,x,f,fvec,fjac,fjac_sparse)

    implicit none

    class(nlesolver_type),intent(inout) :: me
    real(wp),dimension(me%n),intent(in) :: xold      !! previous value of `x`
    real(wp),dimension(me%n),intent(in) :: p         !! search direction
    real(wp),dimension(me%n),intent(out) :: x        !! new `x`
    real(wp),intent(inout) :: f                      !! magnitude of `fvec`
    real(wp),dimension(me%m),intent(inout) :: fvec   !! function vector
    real(wp),dimension(:,:),intent(in),optional :: fjac !! jacobian matrix [dense]
    real(wp),dimension(:),intent(in),optional :: fjac_sparse !! jacobian matrix [sparse]

    real(wp),dimension(:),allocatable :: xnew !! used in [[func_for_fmin]]
    real(wp) :: alpha_min

    allocate(xnew(me%n))

    ! find the minimum value of f in the range of alphas:
    alpha_min = fmin(func_for_fmin,me%alpha_min,me%alpha_max,me%fmin_tol)

    if (me%verbose) write(me%iunit,'(1P,*(A,1X,E16.6))') '        alpha_min = ', alpha_min

    x = xold + p * alpha_min
    if (all(x==xnew)) then
        ! already computed in the func
    else
        call me%func(x,fvec)
        f = norm2(fvec)
    end if

contains

    real(wp) function func_for_fmin(alpha)
    !! function for [[fmin]]
    implicit none
    real(wp),intent(in) :: alpha !! indep variable

    xnew = xold + p * alpha
    call me%func(xnew,fvec)
    func_for_fmin = norm2(fvec) ! return result

    f = func_for_fmin ! just in case this is the solution

    end function func_for_fmin

    end subroutine exact_linesearch
!*****************************************************************************************

!*****************************************************************************************
!>
!  A simple search that just evaluates the function at a specified
!  number of points and picks the one with the minimum function value.

    subroutine fixed_point_linesearch(me,xold,p,x,f,fvec,fjac,fjac_sparse)

    implicit none

    class(nlesolver_type),intent(inout) :: me
    real(wp),dimension(me%n),intent(in) :: xold      !! previous value of `x`
    real(wp),dimension(me%n),intent(in) :: p         !! search direction
    real(wp),dimension(me%n),intent(out) :: x        !! new `x`
    real(wp),intent(inout) :: f                      !! magnitude of `fvec`
    real(wp),dimension(me%m),intent(inout) :: fvec   !! function vector
    real(wp),dimension(:,:),intent(in),optional :: fjac !! jacobian matrix [dense]
    real(wp),dimension(:),intent(in),optional :: fjac_sparse !! jacobian matrix [sparse]

    integer :: i !! counter
    integer :: n_points !! number of points to compute
    real(wp),dimension(:),allocatable :: alphas_to_try !! set of `alpha` values to try
    real(wp),dimension(:),allocatable :: x_tmp !! temp `x`
    real(wp),dimension(:),allocatable :: fvec_tmp !! temp `fvec`
    real(wp) :: f_tmp !! temp `f`
    real(wp) :: step_size !! step size for `alpha`
    integer :: n !! number of steps to divide the interval

    ! 1 o-----------o
    ! 2 o-----o-----o
    ! 3 o---o---o---o

    n = me%n_intervals
    n_points = n + 1

    allocate(alphas_to_try(n_points))
    allocate(x_tmp(me%n))
    allocate(fvec_tmp(me%m))

    step_size = (me%alpha_max - me%alpha_min) / real(n,wp)

    ! compute the alphas:
    alphas_to_try(1) = me%alpha_min
    do i = 2, n
        alphas_to_try(i) = alphas_to_try(i-1) + step_size
    end do
    alphas_to_try(n_points) = me%alpha_max

    ! now compute the functions at these alphas:
    f = big
    do i = 1, n_points

        x_tmp = xold + p * alphas_to_try(i)

        ! evaluate the function at tthis point:
        call me%func(x_tmp,fvec_tmp)
        f_tmp = norm2(fvec_tmp)

        if (f_tmp<=f) then ! new best point
            x = x_tmp
            f = f_tmp
            fvec = fvec_tmp
        end if

    end do

    end subroutine fixed_point_linesearch
!*****************************************************************************************

!******************************************************************************************************
    end module nlesolver_module
!******************************************************************************************************