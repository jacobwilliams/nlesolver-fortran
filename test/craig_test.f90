!******************************************************************************************************
!>
!  Test of the CRAIG minimum-norm sparse solver (`sparsity_mode=6`)
!  for underdetermined systems (`n > m`).
!
!  1. A linear problem `f(x) = A x - b`. Starting from `x=0` with a full step,
!     one Newton iteration gives the minimum-norm solution, which is compared to
!     the dense LAPACK (`dgels`, LQ factorization) minimum-norm solution.
!  2. A larger nonlinear sparse problem solved with each line search option.

    program craig_test

    use nlesolver_module, wp => nlesolver_rk

    implicit none

    integer,parameter :: m_lin = 30   !! linear problem: number of equations
    integer,parameter :: n_lin = 70   !! linear problem: number of variables
    integer,parameter :: m_nl = 2000  !! nonlinear problem: number of equations
    integer,parameter :: n_nl = 5000  !! nonlinear problem: number of variables

    real(wp),dimension(m_lin,n_lin) :: a_lin  !! dense version of the linear problem matrix
    real(wp),dimension(m_lin) :: b_lin
    integer,dimension(:),allocatable :: irow, icol

    call linear_test()
    call nonlinear_test()

    contains

!******************************************************************************************************
    subroutine linear_test()

    type(nlesolver_type) :: solver
    real(wp),dimension(n_lin) :: x_dense, x_craig
    integer :: i, j, k, istat
    character(len=:),allocatable :: message
    real(wp) :: err

    write(*,*) ''
    write(*,*) '***********************'
    write(*,*) '* craig_test : linear'
    write(*,*) '***********************'

    ! a sparse banded matrix with a few extra "long range" entries:
    a_lin = 0.0_wp
    do i = 1, m_lin
        do j = 1, n_lin
            if (abs(2*i-j) <= 2 .or. mod(i*7+j*3, 23) == 0) &
                a_lin(i,j) = 1.0_wp + 0.1_wp*i - 0.05_wp*j + sin(real(i*j,wp))
        end do
        b_lin(i) = cos(real(i,wp))
    end do
    k = count(a_lin /= 0.0_wp)
    allocate(irow(k), icol(k))
    k = 0
    do j = 1, n_lin
        do i = 1, m_lin
            if (a_lin(i,j) /= 0.0_wp) then
                k = k + 1
                irow(k) = i
                icol(k) = j
            end if
        end do
    end do
    write(*,'(A,I0,A,I0,A,I0)') ' m = ', m_lin, ', n = ', n_lin, ', nnz = ', k

    ! dense min-norm solution:
    call solver%initialize(n = n_lin, m = m_lin, max_iter = 1, tol = 1.0e-10_wp, &
                           func = func_lin, grad = grad_lin, sparsity_mode = NLESOLVER_SPARSITY_DENSE)
    x_dense = 0.0_wp
    call solver%solve(x_dense)
    call solver%status(istat, message)
    write(*,'(A,I3,1X,A)') ' dense : ', istat, message

    ! CRAIG min-norm solution:
    call solver%initialize(n = n_lin, m = m_lin, max_iter = 1, tol = 1.0e-10_wp, &
                           func = func_lin, grad_sparse = grad_lin_sparse, &
                           sparsity_mode = NLESOLVER_SPARSITY_CRAIG, irow = irow, icol = icol, &
                           atol = 1.0e-14_wp, btol = 1.0e-14_wp, itnlim = 1000)
    x_craig = 0.0_wp
    call solver%solve(x_craig)
    call solver%status(istat, message)
    write(*,'(A,I3,1X,A)') ' CRAIG : ', istat, message
    if (istat /= 1) error stop 'CRAIG linear test failed to converge'

    err = norm2(x_craig - x_dense) / norm2(x_dense)
    write(*,'(A,1P,E12.4)') ' ||x||   (dense) = ', norm2(x_dense)
    write(*,'(A,1P,E12.4)') ' ||x||   (CRAIG) = ', norm2(x_craig)
    write(*,'(A,1P,E12.4)') ' rel diff        = ', err
    if (err > 1.0e-8_wp) error stop 'CRAIG solution is not the minimum-norm solution'

    deallocate(irow, icol)

    end subroutine linear_test
!******************************************************************************************************

!******************************************************************************************************
    subroutine nonlinear_test()

    !! f_i(x) = x_i^2 + x_{i+m} + 0.5*x_{i+1} + 0.1*x_{n-i+1} - 1,  i = 1..m
    !! (4 nonzeros per row)

    type(nlesolver_type) :: solver
    real(wp),dimension(n_nl) :: x
    real(wp),dimension(m_nl) :: f
    integer :: i, step_mode, istat, ib
    logical :: use_broyden
    character(len=:),allocatable :: message

    write(*,*) ''
    write(*,*) '***********************'
    write(*,*) '* craig_test : nonlinear'
    write(*,*) '***********************'

    allocate(irow(4*m_nl), icol(4*m_nl))
    do i = 1, m_nl
        irow(4*i-3:4*i) = i
        icol(4*i-3:4*i) = [i, i+m_nl, i+1, n_nl-i+1]
    end do
    write(*,'(A,I0,A,I0,A,I0)') ' m = ', m_nl, ', n = ', n_nl, ', nnz = ', size(irow)

    do ib = 1, 2
    use_broyden = ib == 2
    do step_mode = 1, 4
        call solver%initialize(n = n_nl, m = m_nl, max_iter = 50, tol = 1.0e-10_wp, &
                               func = func_nl, grad_sparse = grad_nl_sparse, &
                               step_mode = step_mode, use_broyden = use_broyden, n_intervals = 4, fmin_tol = 1.0e-4_wp, &
                               sparsity_mode = NLESOLVER_SPARSITY_CRAIG, irow = irow, icol = icol, &
                               atol = 1.0e-12_wp, btol = 1.0e-12_wp, itnlim = 1000)
        call solver%status(istat, message)
        if (istat /= 0) error stop 'initialize failed'
        x = 2.0_wp
        call solver%solve(x)
        call solver%status(istat, message)
        call func_nl(solver, x, f)
        write(*,'(A,I0,A,L1,A,I3,1X,A,A,1P,E12.4)') ' step_mode ', step_mode, ' broyden=', use_broyden, &
                                                ' : ', istat, message, '  ||f|| = ', norm2(f)
        ! note: Broyden with constant full steps (step_mode=1) is not globally convergent,
        ! and diverges from this starting point (the same is true with LSQR/LSMR), so it is not checked.
        if (step_mode == 1 .and. use_broyden) cycle
        ! with Broyden, the approximate Jacobian can make the steps smaller than `tolx`
        ! just before reaching `tol` (istat=2), so in that case just require a small residual.
        if (istat == 1 .or. (use_broyden .and. istat == 2 .and. norm2(f) < 1.0e-8_wp)) cycle
        error stop 'CRAIG nonlinear test failed to converge'
    end do
    end do

    end subroutine nonlinear_test
!******************************************************************************************************

    subroutine func_lin(me,x,f)
        class(nlesolver_type),intent(inout) :: me
        real(wp),dimension(:),intent(in)    :: x
        real(wp),dimension(:),intent(out)   :: f
        f = matmul(a_lin, x) - b_lin
    end subroutine func_lin

    subroutine grad_lin(me,x,g)
        class(nlesolver_type),intent(inout) :: me
        real(wp),dimension(:),intent(in)    :: x
        real(wp),dimension(:,:),intent(out) :: g
        g = a_lin
    end subroutine grad_lin

    subroutine grad_lin_sparse(me,x,g)
        class(nlesolver_type),intent(inout) :: me
        real(wp),dimension(:),intent(in)    :: x
        real(wp),dimension(:),intent(out)   :: g
        integer :: k
        do k = 1, size(irow)
            g(k) = a_lin(irow(k), icol(k))
        end do
    end subroutine grad_lin_sparse

    subroutine func_nl(me,x,f)
        class(nlesolver_type),intent(inout) :: me
        real(wp),dimension(:),intent(in)    :: x
        real(wp),dimension(:),intent(out)   :: f
        integer :: i
        do i = 1, m_nl
            f(i) = x(i)**2 + x(i+m_nl) + 0.5_wp*x(i+1) + 0.1_wp*x(n_nl-i+1) - 1.0_wp
        end do
    end subroutine func_nl

    subroutine grad_nl_sparse(me,x,g)
        class(nlesolver_type),intent(inout) :: me
        real(wp),dimension(:),intent(in)    :: x
        real(wp),dimension(:),intent(out)   :: g
        integer :: i
        do i = 1, m_nl
            g(4*i-3:4*i) = [2.0_wp*x(i), 1.0_wp, 0.5_wp, 0.1_wp]
        end do
    end subroutine grad_nl_sparse

!******************************************************************************************************
    end program craig_test
!******************************************************************************************************
