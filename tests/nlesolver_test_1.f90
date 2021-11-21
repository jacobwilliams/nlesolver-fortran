!******************************************************************************************************
!>
!  Tests

    program nlesolver_test_1

    use nlesolver_module
    use iso_fortran_env, only: wp => real64, output_unit

    implicit none

    integer,parameter :: n = 2
    integer,parameter :: m = 2
    integer,parameter :: max_iter = 100
    real(wp),parameter :: tol = 1.0e-8_wp
    logical,parameter :: verbose = .false.

    type(nlesolver_type) :: solver
    real(wp) :: alpha
    logical :: use_broyden
    integer :: step_mode
    integer :: n_intervals
    integer :: istat !! Integer status code.
    character(len=:),allocatable :: message  !! Text status message
    real(wp),dimension(n) :: x
    integer :: f_evals
    integer :: i
    character(len=:),allocatable :: description

    write(*,*) ''
    do i = 1, 8

        select case (i)
        case(1)
            step_mode = 1
            alpha = 1.0_wp
            use_broyden = .false.
            f_evals = 0
            n_intervals = 2
            description = 'Constant alpha'
        case(2)
            step_mode = 1
            alpha = 1.0_wp
            use_broyden = .true.
            f_evals = 0
            n_intervals = 2
            description = 'Constant alpha + broyden'
        case(3)
            step_mode = 2
            alpha = 1.0_wp
            use_broyden = .false.
            f_evals = 0
            n_intervals = 2
            description = 'Backtracking line search'
        case(4)
            step_mode = 2
            alpha = 1.0_wp
            use_broyden = .true.
            f_evals = 0
            n_intervals = 2
            description = 'Backtracking line search + broyden'
        case(5)
            step_mode = 3
            alpha = 1.0_wp
            use_broyden = .false.
            f_evals = 0
            n_intervals = 2
            description = 'Exact line search'
        case(6)
            step_mode = 3
            alpha = 1.0_wp
            use_broyden = .true.
            f_evals = 0
            description = 'Exact line search + broyden'
        case(7)
            step_mode = 4
            alpha = 1.0_wp
            use_broyden = .false.
            f_evals = 0
            n_intervals = 2
            description = 'Fixed point search'
        case(8)
            step_mode = 4
            alpha = 1.0_wp
            use_broyden = .true.
            f_evals = 0
            n_intervals = 2
            description = 'Fixed point search + broyden'
        case default
            error stop 'invalid case'
        end select

        write(*,*) '-------------------------------------------------------'
        write(*,'(A,I3,A,A)') 'Case ', i, ' : ', description
        write(*,*) ''

        call solver%initialize( n = n, &
                                m = m, &
                                max_iter = max_iter, &
                                tol = tol, &
                                func = func, &
                                grad = grad, &
                                step_mode = step_mode,&
                                use_broyden = use_broyden,&
                                export_iteration = export,&
                                n_intervals = n_intervals, &
                                verbose = verbose)
        call solver%status(istat, message)
        write(*,'(I3,1X,A)') istat, message
        if (istat /= 0) error stop

        x = [1.0_wp, 2.0_wp]
        call solver%solve(x)

        call solver%status(istat, message)
        write(*,'(I3,1X,A)') istat, message
        write(*,*) ''

    end do

    contains

    subroutine func(me,x,f)
        !! compute the function
        implicit none
        class(nlesolver_type),intent(inout) :: me
        real(wp),dimension(:),intent(in)    :: x
        real(wp),dimension(:),intent(out)   :: f

        f_evals = f_evals + 1

        f(1) = x(1)**2 + x(2) - 0.1_wp
        f(2) = x(2) + 0.2_wp

        ! root is [0.3, -0.2]

    end subroutine func

    subroutine grad(me,x,g)
        !! compute the gradient of the function (Jacobian):
        implicit none
        class(nlesolver_type),intent(inout) :: me
        real(wp),dimension(:),intent(in)    :: x
        real(wp),dimension(:,:),intent(out) :: g

        f_evals = f_evals + 2   ! to approximate forward diff derivatives

        g(1,1) = 2.0_wp * x(1)  !df(1)/dx
        g(2,1) = 0.0_wp  !df(2)/dx

        g(1,2) = 1.0_wp  !df(1)/dy
        g(2,2) = 1.0_wp  !df(2)/dy

    end subroutine grad

    subroutine export(me,x,f,iter)
        !! export an iteration:
        implicit none
        class(nlesolver_type),intent(inout) :: me
        real(wp),dimension(:),intent(in)    :: x
        real(wp),dimension(:),intent(in)    :: f
        integer,intent(in)                  :: iter !! iteration number

        write(*,'(1P,I3,1X,A,I3,A,*(E15.6))') iter, '(',f_evals,')', x, norm2(f)

    end subroutine export

!******************************************************************************************************
    end program nlesolver_test_1
!******************************************************************************************************