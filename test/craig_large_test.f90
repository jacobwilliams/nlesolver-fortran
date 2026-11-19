!******************************************************************************************************
!>
!  Large, sparse, underdetermined example for the minimum-norm solvers
!  (CRAIG: `sparsity_mode=6`, and the sparse direct Cholesky solver: `sparsity_mode=7`).
!
!  The system has `m` equations and `n > m` variables. Each equation `i` couples a
!  band of `2*hb+1` consecutive variables centered at column `c_i ~ i*n/m`, so the
!  Jacobian is a "stretched" banded matrix:
!
!  \( f_i(x) = \sum_{k=-hb}^{hb} a_{ik} x_{c_i+k} + 0.1 x_{c_i}^3 + 0.2 \sin(x_{c_i-hb}) - b_i \)
!
!  The right-hand side `b` is computed from a known point `x*`, so the system has
!  a solution (actually an (n-m)-dimensional manifold of solutions). Each Newton
!  step is the minimum-norm solution of the linearized system, so starting from
!  `x0` the solver moves to a nearby root.
!
!  For the default sizes (n=50,000, m=30,000), the dense Jacobian would need
!  ~12 GB of memory. The sparse Jacobian has 270,000 nonzeros (~4 MB including indices).
!
!  The same problem is also solved with the sparse direct Cholesky solver
!  (`sparsity_mode=7`), and with LSQR and LSMR for comparison
!  (both also give minimum-norm steps for compatible systems when started from zero).

    program craig_large_test

    use nlesolver_module, wp => nlesolver_rk
    use iso_fortran_env, only: int64

    implicit none

    integer,parameter :: n  = 50000     !! number of variables
    integer,parameter :: m  = 30000     !! number of equations
    integer,parameter :: hb = 4         !! half bandwidth (each row has `2*hb+1` nonzeros)
    integer,parameter :: nb = 2*hb+1    !! nonzeros per row
    integer,parameter :: nnz = m*nb     !! total number of nonzeros
    real(wp),parameter :: tol = 1.0e-10_wp  !! convergence tolerance for `||f||`

    integer,dimension(nnz) :: irow, icol  !! sparsity pattern
    integer,dimension(m)   :: jstart      !! first column in the band for each row
    real(wp),dimension(nnz) :: a          !! constant (linear) part of the Jacobian
    real(wp),dimension(m)   :: b          !! right-hand side
    real(wp),dimension(n)   :: xstar      !! point used to construct `b`
    integer :: n_func, n_grad             !! function and Jacobian evaluation counters
    integer :: iter                       !! last exported iteration number

    call setup()

    write(*,*) ''
    write(*,*) '*******************************************'
    write(*,*) '* craig_large_test'
    write(*,*) '*******************************************'
    write(*,'(A,I0,A,I0,A,I0,A,I0)') ' n = ', n, ', m = ', m, ', nnz = ', nnz, ', bandwidth = ', nb
    write(*,'(A,F8.3,A)') ' dense jacobian would be : ', real(m,wp)*real(n,wp)*storage_size(1.0_wp)/8/1.0e9_wp, ' GB'
    write(*,'(A,F8.3,A)') ' sparse jacobian is      : ', &
                           real(nnz,wp)*(storage_size(1.0_wp)+2*storage_size(1))/8/1.0e6_wp, ' MB (values + indices)'
    write(*,'(A,1P,E12.4,A)') ' ||x* - x0||             : ', norm2(xstar), &
                              '  (x* = point used to construct b; x0 = 0)'
    write(*,*) ''
    write(*,'(A)') ' solver  istat   newton  f-evals  j-evals      ||f||        ||x-x0||     time (s)'
    write(*,'(A)') ' ------  -----   ------  -------  -------   -----------   -----------   --------'

    call go(NLESOLVER_SPARSITY_CRAIG,    'CRAIG', check=.true.)
    call go(NLESOLVER_SPARSITY_CHOLESKY, 'CHOL',  check=.true.)
    call go(NLESOLVER_SPARSITY_LSQR,  'LSQR',  check=.false.)
    call go(NLESOLVER_SPARSITY_LSMR,  'LSMR',  check=.false.)
    write(*,*) ''

    contains

!******************************************************************************************************
    subroutine setup()
        !! construct the sparsity pattern, the linear coefficients, and the right-hand side.

        integer :: i, k, idx, c
        real(wp),dimension(m) :: f
        real(wp),parameter :: pi = acos(-1.0_wp)

        do i = 1, m
            ! center column of the band for this row, kept away from the edges:
            c = 1 + nint(real(i-1,wp) * real(n-1,wp) / real(m-1,wp))
            jstart(i) = max(1, min(n-nb+1, c-hb))
            do k = 1, nb
                idx = (i-1)*nb + k
                irow(idx) = i
                icol(idx) = jstart(i) + k - 1
                ! diagonally weighted coefficients, varying with the row:
                if (k == hb+1) then
                    a(idx) = 4.0_wp + sin(real(i,wp))
                else
                    a(idx) = -1.0_wp / real(abs(k-hb-1),wp) * (1.0_wp + 0.25_wp*cos(real(i+k,wp)))
                end if
            end do
        end do

        ! a known point on the solution manifold:
        do i = 1, n
            xstar(i) = sin(2.0_wp*pi*real(i,wp)/real(n,wp)) + 0.5_wp*cos(6.0_wp*pi*real(i,wp)/real(n,wp))
        end do
        b = 0.0_wp
        call residual(xstar, f)
        b = f

    end subroutine setup
!******************************************************************************************************

!******************************************************************************************************
    subroutine residual(x, f)
        !! f(x) without counting function evaluations
        real(wp),dimension(:),intent(in)  :: x
        real(wp),dimension(:),intent(out) :: f
        integer :: i, j1, jc, idx
        do i = 1, m
            j1  = jstart(i)
            jc  = j1 + hb
            idx = (i-1)*nb
            f(i) = dot_product(a(idx+1:idx+nb), x(j1:j1+nb-1)) &
                   + 0.1_wp*x(jc)**3 + 0.2_wp*sin(x(j1)) - b(i)
        end do
    end subroutine residual
!******************************************************************************************************

!******************************************************************************************************
    subroutine go(sparsity_mode, name, check)

        integer,intent(in)          :: sparsity_mode
        character(len=*),intent(in) :: name
        logical,intent(in)          :: check  !! if true, stop if not converged

        type(nlesolver_type) :: solver
        real(wp),dimension(n) :: x, x0
        real(wp),dimension(m) :: f
        integer :: istat
        character(len=:),allocatable :: message
        integer(int64) :: count_start, count_end, count_rate

        n_func = 0
        n_grad = 0
        iter   = 0

        call solver%initialize(n = n, m = m, max_iter = 50, tol = tol, &
                               func = func, grad_sparse = grad_sparse, &
                               step_mode = NLESOLVER_LINESEARCH_BACKTRACKING, &
                               export_iteration = export, &
                               sparsity_mode = sparsity_mode, irow = irow, icol = icol, &
                               atol = 1.0e-12_wp, btol = 1.0e-12_wp, itnlim = 20000)
        call solver%status(istat, message)
        if (istat /= 0) then
            write(*,*) message
            error stop 'initialize failed'
        end if

        x0 = 0.0_wp
        x  = x0
        call system_clock(count_start, count_rate)
        call solver%solve(x)
        call system_clock(count_end)
        call solver%status(istat, message)

        call residual(x, f)
        write(*,'(1X,A6,I6,3I9,1P,2E14.4,0P,F11.3)') name, istat, iter-1, n_func, n_grad, &
                                                    norm2(f), norm2(x-x0), &
                                                    real(count_end-count_start,wp)/real(count_rate,wp)
        if (istat /= 1) write(*,'(10X,A)') message
        if (check .and. istat /= 1) error stop 'large test failed to converge'

    end subroutine go
!******************************************************************************************************

!******************************************************************************************************
    subroutine export(me,x,f,it)
        !! just keeps track of the number of iterations
        class(nlesolver_type),intent(inout) :: me
        real(wp),dimension(:),intent(in)    :: x
        real(wp),dimension(:),intent(in)    :: f
        integer,intent(in)                  :: it
        iter = it
    end subroutine export
!******************************************************************************************************

!******************************************************************************************************
    subroutine func(me,x,f)
        class(nlesolver_type),intent(inout) :: me
        real(wp),dimension(:),intent(in)    :: x
        real(wp),dimension(:),intent(out)   :: f
        n_func = n_func + 1
        call residual(x, f)
    end subroutine func
!******************************************************************************************************

!******************************************************************************************************
    subroutine grad_sparse(me,x,g)
        !! sparse Jacobian, in the same order as `irow`, `icol`
        class(nlesolver_type),intent(inout) :: me
        real(wp),dimension(:),intent(in)    :: x
        real(wp),dimension(:),intent(out)   :: g
        integer :: i, j1, idx
        n_grad = n_grad + 1
        g = a
        do i = 1, m
            j1  = jstart(i)
            idx = (i-1)*nb
            g(idx+1)    = g(idx+1)    + 0.2_wp*cos(x(j1))          ! d/dx of sin(x_{j1})
            g(idx+hb+1) = g(idx+hb+1) + 0.3_wp*x(j1+hb)**2          ! d/dx of 0.1*x_c^3
        end do
    end subroutine grad_sparse
!******************************************************************************************************

    end program craig_large_test
!******************************************************************************************************
