!******************************************************************************************************
!>
!  Unit tests of the sparse direct minimum-norm solver class ([[cholesky_minnorm_solver]]).
!
!  1. Random sparse matrix (rows in random order): compare with the dense LAPACK
!     minimum-norm solution (`dgels`).
!  2. `damp > 0`: compare with the dense solution of \( (A A^T + \delta^2 I) y = b \), \( x = A^T y \).
!  3. Rank deficient matrix: must report `CHOLESKY_ERR_NOT_PD`, and succeed with `damp > 0`.
!  4. Banded matrix with the rows shuffled: the RCM ordering must recover a small bandwidth.
!  5. Ill-conditioned matrix: iterative refinement must reduce the residual.
!  6. Reuse of the analysis with new values (same sparsity pattern).
!  7. Error codes.

    program cholesky_test

    use nlesolver_module, only: wp => nlesolver_rk
    use nlesolver_cholesky_module
    use iso_fortran_env, only: int64

    implicit none

    integer :: seed = 12345 !! for the random number generator
    integer :: n_failed = 0

    write(*,*) ''
    write(*,*) '***********************'
    write(*,*) '* cholesky_test'
    write(*,*) '***********************'

    call test_random()
    call test_damp()
    call test_rank_deficient()
    call test_rcm()
    call test_refinement()
    call test_reuse()
    call test_errors()

    write(*,*) ''
    if (n_failed > 0) then
        write(*,'(A,I0,A)') ' ', n_failed, ' test(s) FAILED'
        error stop 1
    else
        write(*,'(A)') ' all tests passed'
    end if

    contains

!******************************************************************************************************
    real(wp) function rand()
        !! simple deterministic random number generator in [0,1)
        seed = mod(seed * 16807, 2147483647)
        if (seed < 0) seed = seed + 2147483647
        rand = real(seed,wp) / 2147483647.0_wp
    end function rand

    subroutine check(ok, name, val)
        logical,intent(in) :: ok
        character(len=*),intent(in) :: name
        real(wp),intent(in),optional :: val
        if (present(val)) then
            write(*,'(A,A,A,1P,E12.4)') merge(' [pass] ', ' [FAIL] ', ok), name, ' : ', val
        else
            write(*,'(A,A)') merge(' [pass] ', ' [FAIL] ', ok), name
        end if
        if (.not. ok) n_failed = n_failed + 1
    end subroutine check

    subroutine random_sparse(m, n, nper, irow, icol, a, dense)
        !! random sparse m x n matrix with `nper` nonzeros per row (and a "diagonal"
        !! element to make it full rank), rows generated in random order.
        integer,intent(in) :: m, n, nper
        integer,dimension(:),allocatable,intent(out) :: irow, icol
        real(wp),dimension(:),allocatable,intent(out) :: a
        real(wp),dimension(m,n),intent(out) :: dense
        integer :: i, k, j, idx
        allocate(irow(m*(nper+1)), icol(m*(nper+1)), a(m*(nper+1)))
        dense = 0.0_wp
        idx = 0
        do i = 1, m
            do k = 0, nper
                idx = idx + 1
                if (k == 0) then
                    j = 1 + ((i-1)*n)/m  ! a "diagonal" element
                    a(idx) = 3.0_wp + rand()
                else
                    j = 1 + int(rand()*n)
                    a(idx) = rand() - 0.5_wp
                end if
                irow(idx) = i
                icol(idx) = j
                dense(i,j) = dense(i,j) + a(idx)  ! duplicates are summed
            end do
        end do
    end subroutine random_sparse

    function dense_minnorm(m, n, a, b) result(x)
        !! dense minimum-norm solution with LAPACK dgels
        integer,intent(in) :: m, n
        real(wp),dimension(m,n),intent(in) :: a
        real(wp),dimension(m),intent(in) :: b
        real(wp),dimension(n) :: x
        interface
            subroutine dgels(trans,m,n,nrhs,a,lda,b,ldb,work,lwork,info)
                import :: wp
                character :: trans
                integer :: info, lda, ldb, lwork, m, n, nrhs
                real(wp) :: a(lda,*), b(ldb,*), work(*)
            end subroutine dgels
        end interface
        real(wp),dimension(m,n) :: acopy
        real(wp),dimension(max(m,n),1) :: bcopy
        real(wp),dimension(:),allocatable :: work
        integer :: info
        acopy = a
        bcopy = 0.0_wp
        bcopy(1:m,1) = b
        allocate(work(64*(m+n)))
        call dgels('N', m, n, 1, acopy, m, bcopy, max(m,n), work, size(work), info)
        if (info /= 0) error stop 'dgels failed'
        x = bcopy(1:n,1)
    end function dense_minnorm

    function dense_damped(m, n, a, b, damp) result(x)
        !! x = A^T y, (A A^T + damp^2 I) y = b, with LAPACK dgesv
        integer,intent(in) :: m, n
        real(wp),dimension(m,n),intent(in) :: a
        real(wp),dimension(m),intent(in) :: b
        real(wp),intent(in) :: damp
        real(wp),dimension(n) :: x
        interface
            subroutine dgesv(n,nrhs,a,lda,ipiv,b,ldb,info)
                import :: wp
                integer :: info, lda, ldb, n, nrhs, ipiv(*)
                real(wp) :: a(lda,*), b(ldb,*)
            end subroutine dgesv
        end interface
        real(wp),dimension(m,m) :: aat
        real(wp),dimension(m,1) :: y
        integer,dimension(m) :: ipiv
        integer :: i, info
        aat = matmul(a, transpose(a))
        do i = 1, m
            aat(i,i) = aat(i,i) + damp**2
        end do
        y(:,1) = b
        call dgesv(m, 1, aat, m, ipiv, y, m, info)
        if (info /= 0) error stop 'dgesv failed'
        x = matmul(transpose(a), y(:,1))
    end function dense_damped

!******************************************************************************************************
    subroutine test_random()
        integer,parameter :: m = 200, n = 500
        integer,dimension(:),allocatable :: irow, icol
        real(wp),dimension(:),allocatable :: a
        real(wp),dimension(m,n) :: dense
        real(wp),dimension(m) :: b
        real(wp),dimension(n) :: x, xd
        type(cholesky_minnorm_solver) :: solver
        integer :: istat, i
        real(wp) :: rnorm

        call random_sparse(m, n, 6, irow, icol, a, dense)
        do i = 1, m
            b(i) = rand() - 0.5_wp
        end do
        call solver%analyze(m, n, irow, icol, istat)
        if (istat == 0) call solver%factor(a, 0.0_wp, istat)
        if (istat == 0) call solver%solve(b, x, istat, rnorm=rnorm)
        call check(istat == 0, 'random 200x500: status')
        xd = dense_minnorm(m, n, dense, b)
        call check(norm2(x-xd)/norm2(xd) < 1.0e-12_wp, 'random 200x500: rel diff vs dgels', norm2(x-xd)/norm2(xd))
        call check(norm2(matmul(dense,x)-b) < 1.0e-12_wp, 'random 200x500: ||Ax-b||', norm2(matmul(dense,x)-b))
    end subroutine test_random

!******************************************************************************************************
    subroutine test_damp()
        integer,parameter :: m = 100, n = 250
        real(wp),parameter :: damp = 0.7_wp
        integer,dimension(:),allocatable :: irow, icol
        real(wp),dimension(:),allocatable :: a
        real(wp),dimension(m,n) :: dense
        real(wp),dimension(m) :: b
        real(wp),dimension(n) :: x, xd
        type(cholesky_minnorm_solver) :: solver
        integer :: istat, i

        call random_sparse(m, n, 4, irow, icol, a, dense)
        do i = 1, m
            b(i) = rand() - 0.5_wp
        end do
        call solver%analyze(m, n, irow, icol, istat)
        if (istat == 0) call solver%factor(a, damp, istat)
        if (istat == 0) call solver%solve(b, x, istat)
        call check(istat == 0, 'damp: status')
        xd = dense_damped(m, n, dense, b, damp)
        call check(norm2(x-xd)/norm2(xd) < 1.0e-12_wp, 'damp: rel diff vs dense', norm2(x-xd)/norm2(xd))
    end subroutine test_damp

!******************************************************************************************************
    subroutine test_rank_deficient()
        !! row 2 is a multiple of row 1
        integer,dimension(6),parameter :: irow = [1,1,2,2,3,3]
        integer,dimension(6),parameter :: icol = [1,2,1,2,3,4]
        real(wp),dimension(6),parameter :: a = [1.0_wp,2.0_wp,2.0_wp,4.0_wp,1.0_wp,1.0_wp]
        real(wp),dimension(3) :: b = [1.0_wp, 2.0_wp, 1.0_wp]
        real(wp),dimension(4) :: x
        type(cholesky_minnorm_solver) :: solver
        integer :: istat
        call solver%analyze(3, 4, irow, icol, istat)
        call solver%factor(a, 0.0_wp, istat)
        call check(istat == CHOLESKY_ERR_NOT_PD, 'rank deficient: detected (damp=0)')
        call solver%factor(a, 1.0e-3_wp, istat)
        if (istat == 0) call solver%solve(b, x, istat)
        call check(istat == 0, 'rank deficient: solved with damp>0')
    end subroutine test_rank_deficient

!******************************************************************************************************
    subroutine test_rcm()
        !! banded matrix (bandwidth of A*A^T = 2) with the rows randomly shuffled.
        !! the natural order of the shuffled matrix has a large bandwidth.
        integer,parameter :: m = 5000, n = 2*m+1
        integer,dimension(:),allocatable :: irow, icol, p
        real(wp),dimension(:),allocatable :: a
        real(wp),dimension(m) :: b, r
        real(wp),dimension(n) :: x
        type(cholesky_minnorm_solver) :: solver
        integer :: i, k, j, tmp, istat, bw_natural
        allocate(irow(3*m), icol(3*m), a(3*m), p(m))
        p = [(i, i=1,m)]
        do i = m, 2, -1  ! random permutation
            j = 1 + int(rand()*i)
            tmp = p(i); p(i) = p(j); p(j) = tmp
        end do
        k = 0
        do i = 1, m
            do j = 0, 2
                k = k + 1
                irow(k) = p(i)        ! shuffled row index
                icol(k) = 2*i - 1 + j
                a(k) = merge(4.0_wp, -1.0_wp, j==1)
            end do
            b(p(i)) = 1.0_wp
        end do
        ! bandwidth of A*A^T in the shuffled (natural) order:
        bw_natural = 0
        do i = 2, m
            bw_natural = max(bw_natural, abs(p(i)-p(i-1)))
        end do
        call solver%analyze(m, n, irow, icol, istat)
        if (istat == 0) call solver%factor(a, 0.0_wp, istat)
        if (istat == 0) call solver%solve(b, x, istat)
        call check(istat == 0, 'rcm: status')
        write(*,'(A,I0,A,I0,A,I0)') '        bandwidth: shuffled = ', bw_natural, ', after RCM = ', &
                                    solver%bandwidth(), ', envelope size = ', solver%envelope_size()
        call check(solver%bandwidth() >= 0 .and. solver%bandwidth() <= 2, 'rcm: recovered bandwidth')
        r = b
        do k = 1, size(a)
            r(irow(k)) = r(irow(k)) - a(k)*x(icol(k))
        end do
        call check(norm2(r) < 1.0e-12_wp, 'rcm: ||Ax-b||', norm2(r))
    end subroutine test_rcm

!******************************************************************************************************
    subroutine test_refinement()
        !! ill-conditioned matrix: rows 2k-1 and 2k are identical except for a small perturbation.
        integer,parameter :: m = 50, n = 120
        real(wp),parameter :: delta = 1.0e-6_wp  !! cond(A) ~ 1/delta, cond(A*A^T) ~ 1/delta^2
        integer,dimension(:),allocatable :: irow, icol
        real(wp),dimension(:),allocatable :: a
        real(wp),dimension(m,n) :: dense
        real(wp),dimension(m) :: b
        real(wp),dimension(n) :: x0, x2, xd
        type(cholesky_minnorm_solver) :: solver
        integer :: i, k, ip, istat
        real(wp) :: r0, r2
        allocate(irow(4*m), icol(4*m), a(4*m))
        dense = 0.0_wp
        k = 0
        do i = 1, m
            ip = (i+1)/2  ! pair number
            irow(k+1:k+4) = i
            icol(k+1:k+4) = [2*ip-1, 2*ip, 2*ip+1, 60+i]
            a(k+1:k+4)    = [1.0_wp, 2.0_wp, 1.0_wp, delta*(1.0_wp+rand())]
            dense(i, icol(k+1:k+4)) = a(k+1:k+4)
            k = k + 4
            b(i) = rand()
        end do
        call solver%analyze(m, n, irow, icol, istat)
        if (istat == 0) call solver%factor(a, 0.0_wp, istat)
        if (istat == 0) call solver%solve(b, x0, istat, n_refine=0, rnorm=r0)
        if (istat == 0) call solver%solve(b, x2, istat, rnorm=r2)
        call check(istat == 0, 'refinement: status')
        write(*,'(A,1P,E12.4,A,E12.4)') '        ||r|| no refinement = ', r0, ',  with refinement = ', r2
        call check(r2 < r0, 'refinement: residual reduced')
        xd = dense_minnorm(m, n, dense, b)
        call check(norm2(x2-xd)/norm2(xd) < 1.0e-8_wp, 'refinement: rel diff vs dgels', norm2(x2-xd)/norm2(xd))
    end subroutine test_refinement

!******************************************************************************************************
    subroutine test_reuse()
        !! analyze once, factor with two different sets of values
        integer,parameter :: m = 80, n = 150
        integer,dimension(:),allocatable :: irow, icol
        real(wp),dimension(:),allocatable :: a
        real(wp),dimension(m,n) :: dense
        real(wp),dimension(m) :: b
        real(wp),dimension(n) :: x, xd
        type(cholesky_minnorm_solver) :: solver
        integer :: istat, i, k
        call random_sparse(m, n, 5, irow, icol, a, dense)
        b = 1.0_wp
        call solver%analyze(m, n, irow, icol, istat)
        call solver%factor(a, 0.0_wp, istat)
        call solver%solve(b, x, istat)
        ! new values, same pattern:
        a = a * 2.0_wp + 0.1_wp
        dense = 0.0_wp
        do k = 1, size(a)
            dense(irow(k),icol(k)) = dense(irow(k),icol(k)) + a(k)
        end do
        call solver%factor(a, 0.0_wp, istat)
        if (istat == 0) call solver%solve(b, x, istat)
        xd = dense_minnorm(m, n, dense, b)
        call check(istat == 0 .and. norm2(x-xd)/norm2(xd) < 1.0e-12_wp, 'reuse analysis: rel diff vs dgels', &
                   norm2(x-xd)/norm2(xd))
    end subroutine test_reuse

!******************************************************************************************************
    subroutine test_errors()
        type(cholesky_minnorm_solver) :: solver
        real(wp),dimension(2) :: b = 1.0_wp
        real(wp),dimension(3) :: x
        integer :: istat
        call solver%solve(b, x, istat)
        call check(istat == CHOLESKY_ERR_NOT_READY, 'errors: solve before factor')
        call solver%factor([1.0_wp], 0.0_wp, istat)
        call check(istat == CHOLESKY_ERR_NOT_READY, 'errors: factor before analyze')
        call solver%analyze(2, 3, [1,3], [1,2], istat)
        call check(istat == CHOLESKY_ERR_INPUT, 'errors: row index out of range')
        call solver%analyze(2, 3, [1,2], [1,2,3], istat)
        call check(istat == CHOLESKY_ERR_INPUT, 'errors: irow/icol size mismatch')
    end subroutine test_errors

!******************************************************************************************************
    end program cholesky_test
!******************************************************************************************************
