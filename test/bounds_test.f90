!******************************************************************************************************
!>
!  Test of the bounds modes: for f(x) = x - 5 with 0 <= x <= 2, starting at x = 0,
!  one Newton step (p = 5) must stop at the upper bound x = 2 in every bounds mode:
!
!  * scalar: the component is clipped to the bound.
!  * vector: the step is shortened along the Newton direction (t = 2/5).
!  * wall:   the component is put on the bound (and held there during the line search).

    program bounds_test

    use nlesolver_module, wp => nlesolver_rk

    implicit none

    type(nlesolver_type) :: solver
    real(wp),dimension(1) :: x
    integer :: mode, step_mode
    character(len=6),dimension(3),parameter :: names = ['scalar','vector','wall  ']

    write(*,*) ''
    write(*,*) '***********************'
    write(*,*) '* bounds_test'
    write(*,*) '***********************'

    do step_mode = 1, 4
        do mode = NLESOLVER_SCALAR_BOUNDS, NLESOLVER_WALL_BOUNDS
            call solver%initialize(n=1, m=1, max_iter=1, tol=1.0e-10_wp, func=func, grad=grad, &
                                   step_mode=step_mode, alpha_min=0.5_wp, &
                                   bounds_mode=mode, xlow=[0.0_wp], xupp=[2.0_wp])
            x = 0.0_wp
            call solver%solve(x)
            write(*,'(A,I0,A,A6,A,F10.6)') ' step_mode=', step_mode, '  bounds_mode=', names(mode), &
                                           '  x = ', x(1)
            if (x(1) < 0.0_wp .or. x(1) > 2.0_wp) error stop 'bounds violated'
            ! the full step to the bound is the best step for all line searches here:
            if (abs(x(1) - 2.0_wp) > 1.0e-6_wp) error stop 'step did not reach the bound'
        end do
    end do

    contains

    subroutine func(me,x,f)
        class(nlesolver_type),intent(inout) :: me
        real(wp),dimension(:),intent(in)    :: x
        real(wp),dimension(:),intent(out)   :: f
        f(1) = x(1) - 5.0_wp
    end subroutine func

    subroutine grad(me,x,g)
        class(nlesolver_type),intent(inout) :: me
        real(wp),dimension(:),intent(in)    :: x
        real(wp),dimension(:,:),intent(out) :: g
        g = 1.0_wp
    end subroutine grad

    end program bounds_test
!******************************************************************************************************
