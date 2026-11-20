!******************************************************************************************************
!>
!  Tests of the bounds modes. When a bounds mode is enabled, the bounds must never be
!  ignored: the function must never be evaluated outside of them.
!
!  1. `f(x) = x - 5` with `0 <= x <= 2`, starting at `x = 0`: one Newton step (`p = 5`)
!     must stop at the upper bound `x = 2` in every bounds mode:
!     * scalar: the component is clipped to the bound.
!     * vector: the step is shortened along the Newton direction (`t = 2/5`).
!     * wall:   the component is put on the bound (and held there during the line search).
!  2. Same problem, starting on the upper bound (`x = 2`), so the Newton step points outward:
!     no step is possible, and `x` must stay at 2 in every mode.
!  3. Two variables. Wall mode with a half step (`alpha = 0.5`): the bounded component
!     is held on the bound, while the other one takes the half step. Vector mode, starting
!     with one variable on its bound and the step pointing outward: the step direction is
!     kept, so no step is possible (the other variable must not move either).
!  4. Randomized: many random problems, bounds, and starting points, for every bounds mode and
!     step mode. The function routine itself checks that it is never called outside the bounds
!     (including by rounding, e.g. `x + (xupp - x) > xupp`).
!  5. Input errors: `bounds_mode > 0` without `xlow` and `xupp`, and an invalid `bounds_mode`.
!  6. The bounds are checked on the largest step actually attempted (`alpha`, or `alpha_max`),
!     not on the full Newton step: with `0 <= x <= 2`, `x0 = 0`, and a half step:
!     * `f = x - 3`: the full step (to 3) would violate the bound, but the half step (to 1.5)
!       does not, so no variable may be put on the bound: `x = 1.5` in every mode.
!     * `f = x - 5`: the half step (to 2.5) violates the bound: `x = 2` in every mode.

    program bounds_test

    use nlesolver_module, wp => nlesolver_rk

    implicit none

    character(len=6),dimension(3),parameter :: names = ['scalar','vector','wall  ']

    ! for the randomized test:
    integer,parameter :: nr = 20                 !! number of variables (and equations)
    real(wp),dimension(nr) :: xlo, xhi, target   !! bounds and the (unconstrained) root
    integer :: n_evals = 0, n_violations = 0
    integer :: seed = 12345

    write(*,*) ''
    write(*,*) '***********************'
    write(*,*) '* bounds_test'
    write(*,*) '***********************'

    call test_one_step()
    call test_start_on_bound()
    call test_wall_two_variables()
    call test_random()
    call test_errors()
    call test_attempted_step()

    write(*,'(/A)') ' all bounds tests passed'

    contains

!******************************************************************************************************
    subroutine test_one_step()
        type(nlesolver_type) :: solver
        real(wp),dimension(1) :: x
        integer :: mode, step_mode
        write(*,'(/A)') ' 1. one step from x=0 to the bound x=2:'
        do step_mode = 1, 4
            do mode = NLESOLVER_SCALAR_BOUNDS, NLESOLVER_WALL_BOUNDS
                call solver%initialize(n=1, m=1, max_iter=1, tol=1.0e-10_wp, func=f1, grad=g1, &
                                       step_mode=step_mode, alpha_min=0.5_wp, &
                                       bounds_mode=mode, xlow=[0.0_wp], xupp=[2.0_wp])
                x = 0.0_wp
                call solver%solve(x)
                write(*,'(A,I0,A,A6,A,F10.6)') '   step_mode=', step_mode, '  bounds_mode=', names(mode), &
                                               '  x = ', x(1)
                if (x(1) < 0.0_wp .or. x(1) > 2.0_wp) error stop 'bounds violated'
                ! the full step to the bound is the best step for all line searches here:
                if (abs(x(1) - 2.0_wp) > 1.0e-6_wp) error stop 'step did not reach the bound'
            end do
        end do
    end subroutine test_one_step

!******************************************************************************************************
    subroutine test_start_on_bound()
        type(nlesolver_type) :: solver
        real(wp),dimension(1) :: x
        integer :: mode, step_mode, istat
        character(len=:),allocatable :: message
        write(*,'(/A)') ' 2. starting on the bound x=2, step pointing outward:'
        do step_mode = 1, 4
            do mode = NLESOLVER_SCALAR_BOUNDS, NLESOLVER_WALL_BOUNDS
                call solver%initialize(n=1, m=1, max_iter=5, tol=1.0e-10_wp, func=f1, grad=g1, &
                                       step_mode=step_mode, &
                                       bounds_mode=mode, xlow=[0.0_wp], xupp=[2.0_wp])
                x = 2.0_wp
                call solver%solve(x)
                call solver%status(istat, message)
                write(*,'(A,I0,A,A6,A,F10.6,A,I0,1X,A)') '   step_mode=', step_mode, '  bounds_mode=', names(mode), &
                                                         '  x = ', x(1), '  status ', istat, message
                if (x(1) /= 2.0_wp) error stop 'x left the bound (bounds ignored)'
                if (istat == 1) error stop 'reported convergence outside the bounds'
            end do
        end do
    end subroutine test_start_on_bound

!******************************************************************************************************
    subroutine test_wall_two_variables()
        !! f = [x1 - 5, x2 - 1], 0 <= x1 <= 2, 0 <= x2 <= 10, x0 = 0, alpha = 0.5
        type(nlesolver_type) :: solver
        real(wp),dimension(2) :: x
        write(*,'(/A)') ' 3. two variables:'
        write(*,'(A)') '   wall mode, alpha = 0.5:'
        call solver%initialize(n=2, m=2, max_iter=1, tol=1.0e-10_wp, func=f2, grad=g2, alpha=0.5_wp, &
                               bounds_mode=NLESOLVER_WALL_BOUNDS, xlow=[0.0_wp,0.0_wp], xupp=[2.0_wp,10.0_wp])
        x = 0.0_wp
        call solver%solve(x)
        write(*,'(A,2F10.6)') '   x = ', x
        if (x(1) /= 2.0_wp) error stop 'wall: the bounded component was not held on the bound'
        if (abs(x(2) - 0.5_wp) > 1.0e-12_wp) error stop 'wall: the free component did not take the half step'

        write(*,'(A)') '   vector mode, starting with x1 on its bound (step pointing outward):'
        call solver%initialize(n=2, m=2, max_iter=1, tol=1.0e-10_wp, func=f2, grad=g2, &
                               bounds_mode=NLESOLVER_VECTOR_BOUNDS, xlow=[0.0_wp,0.0_wp], xupp=[2.0_wp,10.0_wp])
        x = [2.0_wp, 0.0_wp]
        call solver%solve(x)
        write(*,'(A,2F10.6)') '   x = ', x
        if (any(x /= [2.0_wp, 0.0_wp])) error stop 'vector: a step was taken that changes the Newton direction'
    end subroutine test_wall_two_variables

!******************************************************************************************************
    subroutine test_random()
        type(nlesolver_type) :: solver
        real(wp),dimension(nr) :: x
        integer :: itest, mode, step_mode, i, istat
        integer,parameter :: n_tests = 200
        write(*,'(/A)') ' 4. randomized test (the function is never evaluated outside the bounds):'
        n_evals = 0
        n_violations = 0
        do itest = 1, n_tests
            do i = 1, nr
                xlo(i)    = -1.0_wp - 10.0_wp*rand()
                xhi(i)    =  1.0_wp + 10.0_wp*rand()
                target(i) = 30.0_wp*(rand() - 0.5_wp)   ! often outside the bounds
            end do
            do step_mode = 1, 4
                do mode = NLESOLVER_SCALAR_BOUNDS, NLESOLVER_WALL_BOUNDS
                    call solver%initialize(n=nr, m=nr, max_iter=10, tol=1.0e-10_wp, func=fr, grad=gr, &
                                           step_mode=step_mode, bounds_mode=mode, xlow=xlo, xupp=xhi)
                    call solver%status(istat)
                    if (istat /= 0) error stop 'initialize failed'
                    ! random starting point, sometimes outside of the bounds
                    ! (it is projected onto the bounds before the first evaluation):
                    do i = 1, nr
                        x(i) = xlo(i) + (xhi(i)-xlo(i))*(1.4_wp*rand() - 0.2_wp)
                    end do
                    call solver%solve(x)
                    if (any(x < xlo) .or. any(x > xhi)) error stop 'final x is outside the bounds'
                end do
            end do
        end do
        write(*,'(A,I0,A,I0,A)') '   ', n_evals, ' function evaluations, ', n_violations, ' outside the bounds'
        if (n_violations > 0) error stop 'the function was evaluated outside the bounds'
    end subroutine test_random

!******************************************************************************************************
    subroutine test_errors()
        type(nlesolver_type) :: solver
        integer :: istat
        character(len=:),allocatable :: message
        write(*,'(/A)') ' 5. input errors:'
        call solver%initialize(n=1, m=1, max_iter=1, tol=1.0e-10_wp, func=f1, grad=g1, &
                               bounds_mode=NLESOLVER_VECTOR_BOUNDS)
        call solver%status(istat, message)
        write(*,'(A,I0,1X,A)') '   bounds_mode without xlow/xupp: ', istat, message
        if (istat /= -19) error stop 'bounds_mode > 0 without xlow and xupp was not rejected'
        call solver%initialize(n=1, m=1, max_iter=1, tol=1.0e-10_wp, func=f1, grad=g1, &
                               bounds_mode=7, xlow=[0.0_wp], xupp=[2.0_wp])
        call solver%status(istat, message)
        write(*,'(A,I0,1X,A)') '   invalid bounds_mode:            ', istat, message
        if (istat /= -20) error stop 'invalid bounds_mode was not rejected'
    end subroutine test_errors

!******************************************************************************************************
    subroutine test_attempted_step()
        type(nlesolver_type) :: solver
        real(wp),dimension(1) :: x
        integer :: mode, icase, istep
        real(wp) :: expected
        integer,dimension(2),parameter :: step_modes = [NLESOLVER_LINESEARCH_SIMPLE, NLESOLVER_LINESEARCH_BACKTRACKING]
        write(*,'(/A)') ' 6. bounds checked on the attempted step (alpha = alpha_max = 0.5):'
        do icase = 1, 2
            target(1) = merge(3.0_wp, 5.0_wp, icase==1)
            expected  = merge(1.5_wp, 2.0_wp, icase==1)
            do istep = 1, 2
                do mode = NLESOLVER_SCALAR_BOUNDS, NLESOLVER_WALL_BOUNDS
                    call solver%initialize(n=1, m=1, max_iter=1, tol=1.0e-10_wp, func=ft, grad=g1, &
                                           step_mode=step_modes(istep), alpha=0.5_wp, &
                                           alpha_min=0.1_wp, alpha_max=0.5_wp, &
                                           bounds_mode=mode, xlow=[0.0_wp], xupp=[2.0_wp])
                    x = 0.0_wp
                    call solver%solve(x)
                    write(*,'(A,F4.1,A,I0,A,A6,A,F10.6)') '   f = x - ', target(1), '  step_mode=', step_modes(istep), &
                                                          '  bounds_mode=', names(mode), '  x = ', x(1)
                    if (abs(x(1) - expected) > 1.0e-12_wp) error stop 'bounds not checked on the attempted step'
                end do
            end do
        end do
    end subroutine test_attempted_step

    subroutine ft(me,x,f)
        !! f = x - target(1)
        class(nlesolver_type),intent(inout) :: me
        real(wp),dimension(:),intent(in)    :: x
        real(wp),dimension(:),intent(out)   :: f
        f(1) = x(1) - target(1)
    end subroutine ft

!******************************************************************************************************
    real(wp) function rand()
        !! simple deterministic random number generator in [0,1)
        seed = int(mod(int(seed,8) * 16807_8, 2147483647_8))
        rand = real(seed,wp) / 2147483647.0_wp
    end function rand

    subroutine f1(me,x,f)
        class(nlesolver_type),intent(inout) :: me
        real(wp),dimension(:),intent(in)    :: x
        real(wp),dimension(:),intent(out)   :: f
        f(1) = x(1) - 5.0_wp
    end subroutine f1

    subroutine g1(me,x,g)
        class(nlesolver_type),intent(inout) :: me
        real(wp),dimension(:),intent(in)    :: x
        real(wp),dimension(:,:),intent(out) :: g
        g = 1.0_wp
    end subroutine g1

    subroutine f2(me,x,f)
        class(nlesolver_type),intent(inout) :: me
        real(wp),dimension(:),intent(in)    :: x
        real(wp),dimension(:),intent(out)   :: f
        f = [x(1) - 5.0_wp, x(2) - 1.0_wp]
    end subroutine f2

    subroutine g2(me,x,g)
        class(nlesolver_type),intent(inout) :: me
        real(wp),dimension(:),intent(in)    :: x
        real(wp),dimension(:,:),intent(out) :: g
        g = 0.0_wp
        g(1,1) = 1.0_wp
        g(2,2) = 1.0_wp
    end subroutine g2

    subroutine fr(me,x,f)
        !! f_i = x_i - target_i + 0.1*x_{i+1}**2  (coupled, root often outside the bounds).
        !! counts the evaluations outside the bounds.
        class(nlesolver_type),intent(inout) :: me
        real(wp),dimension(:),intent(in)    :: x
        real(wp),dimension(:),intent(out)   :: f
        integer :: i
        n_evals = n_evals + 1
        if (any(x < xlo) .or. any(x > xhi)) n_violations = n_violations + 1
        do i = 1, nr
            f(i) = x(i) - target(i)
            if (i < nr) f(i) = f(i) + 0.1_wp*x(i+1)**2
        end do
    end subroutine fr

    subroutine gr(me,x,g)
        class(nlesolver_type),intent(inout) :: me
        real(wp),dimension(:),intent(in)    :: x
        real(wp),dimension(:,:),intent(out) :: g
        integer :: i
        g = 0.0_wp
        do i = 1, nr
            g(i,i) = 1.0_wp
            if (i < nr) g(i,i+1) = 0.2_wp*x(i+1)
        end do
    end subroutine gr

    end program bounds_test
!******************************************************************************************************
