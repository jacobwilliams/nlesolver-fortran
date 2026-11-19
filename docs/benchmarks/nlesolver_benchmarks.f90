!******************************************************************************************************
!>
!  Benchmarks for the nlesolver-fortran technical note (`docs/nlesolver_technote.tex`).
!
!  Writes CSV files to the directory given as the first command line argument
!  (default: `docs/data`). An optional second argument is a comma separated list of the
!  studies to run (`spy,scaling,square,steps,cond,sparsity`; default: all).
!  Run from the repository root with:
!
!```
!  fpm run --example nlesolver_benchmarks --profile release -- docs/data
!```
!
!  Studies:
!
!  1. `scaling_*.csv`    : banded underdetermined problem (m = 0.6 n), all linear solvers, increasing n.
!  2. `square_*.csv`     : banded square problem (m = n), all linear solvers.
!  3. `steps_*.csv`, `hist_*.csv` : step modes and Broyden updates (n = 50,000).
!  4. `cond_*.csv`       : effect of the condition number (CRAIG iterations; Cholesky accuracy).
!  5. `sparsity_*.csv`   : banded, shuffled-banded, and 2-D grid-like sparsity (Cholesky envelope growth).
!  6. `spy_*.csv`        : small sparsity patterns for the figures.

    module bench_problem_module

    use nlesolver_module, wp => nlesolver_rk

    implicit none

    public

    ! generic problem definition:
    !   equation i (stored at f(rl(i))):
    !   f = sum_k coef(k)*x(icol(k)) + ccoef*x(icol(icub(i)))**3 + scoef*sin(x(icol(isin(i)))) - b
    integer :: n = 0, m = 0, nnz = 0
    integer,dimension(:),allocatable :: irow, icol    !! sparsity pattern (irow = rl(row))
    integer,dimension(:),allocatable :: k1, k2        !! COO range of each equation
    integer,dimension(:),allocatable :: rl            !! row label of each equation
    integer,dimension(:),allocatable :: icub, isin    !! COO index of the nonlinear terms (0 = none)
    real(wp),dimension(:),allocatable :: coef         !! linear coefficients
    real(wp),dimension(:),allocatable :: b            !! right-hand side (by row label)
    real(wp) :: ccoef = 0.0_wp, scoef = 0.0_wp

    integer :: nf = 0, ng = 0, last_iter = 0
    logical :: record_history = .false.
    integer :: nhist = 0
    real(wp),dimension(1000) :: hist

    contains

!******************************************************************************************************
    subroutine eval_f(x, f)
        real(wp),dimension(:),intent(in)  :: x
        real(wp),dimension(:),intent(out) :: f
        integer :: i, k
        real(wp) :: s
        do i = 1, m
            s = -b(rl(i))
            do k = k1(i), k2(i)
                s = s + coef(k)*x(icol(k))
            end do
            if (icub(i) > 0) s = s + ccoef*x(icol(icub(i)))**3
            if (isin(i) > 0) s = s + scoef*sin(x(icol(isin(i))))
            f(rl(i)) = s
        end do
    end subroutine eval_f

    subroutine eval_g(x, g)
        real(wp),dimension(:),intent(in)  :: x
        real(wp),dimension(:),intent(out) :: g
        integer :: i
        g = coef
        do i = 1, m
            if (icub(i) > 0) g(icub(i)) = g(icub(i)) + 3.0_wp*ccoef*x(icol(icub(i)))**2
            if (isin(i) > 0) g(isin(i)) = g(isin(i)) + scoef*cos(x(icol(isin(i))))
        end do
    end subroutine eval_g

    subroutine func(me, x, f)
        class(nlesolver_type),intent(inout) :: me
        real(wp),dimension(:),intent(in)    :: x
        real(wp),dimension(:),intent(out)   :: f
        nf = nf + 1
        call eval_f(x, f)
    end subroutine func

    subroutine grad_sparse(me, x, g)
        class(nlesolver_type),intent(inout) :: me
        real(wp),dimension(:),intent(in)    :: x
        real(wp),dimension(:),intent(out)   :: g
        ng = ng + 1
        call eval_g(x, g)
    end subroutine grad_sparse

    subroutine grad_dense(me, x, g)
        class(nlesolver_type),intent(inout) :: me
        real(wp),dimension(:),intent(in)    :: x
        real(wp),dimension(:,:),intent(out) :: g
        real(wp),dimension(nnz) :: gs
        integer :: k
        ng = ng + 1
        call eval_g(x, gs)
        g = 0.0_wp
        do k = 1, nnz
            g(irow(k),icol(k)) = g(irow(k),icol(k)) + gs(k)
        end do
    end subroutine grad_dense

    subroutine export(me, x, f, iter)
        class(nlesolver_type),intent(inout) :: me
        real(wp),dimension(:),intent(in)    :: x
        real(wp),dimension(:),intent(in)    :: f
        integer,intent(in)                  :: iter
        last_iter = iter
        if (record_history .and. iter <= size(hist)) then
            nhist = iter
            hist(iter) = norm2(f)
        end if
    end subroutine export

!******************************************************************************************************
    subroutine set_rhs_from(xstar)
        !! b = (f without b) at xstar, so that xstar is a solution
        real(wp),dimension(n),intent(in) :: xstar
        real(wp),dimension(m) :: f
        b = 0.0_wp
        call eval_f(xstar, f)
        b = f
    end subroutine set_rhs_from

!******************************************************************************************************
    subroutine make_banded(nn, mm, hb, d0, shuffle, nonlinear)
        !! the banded problem: equation i couples 2*hb+1 consecutive variables
        !! centered at column ~ i*n/m (see test/craig_large_test.f90).
        integer,intent(in)  :: nn, mm, hb
        real(wp),intent(in) :: d0       !! diagonal weight (smaller = worse conditioned)
        logical,intent(in)  :: shuffle  !! randomly permute the row labels
        logical,intent(in)  :: nonlinear
        integer :: i, k, idx, c, j1, nb, j, tmp, sd, d
        real(wp),dimension(:),allocatable :: xstar
        real(wp),parameter :: pi = acos(-1.0_wp)
        n = nn; m = mm; nb = 2*hb+1; nnz = m*nb
        call alloc_problem()
        rl = [(i, i=1,m)]
        if (shuffle) then
            sd = 12345
            do i = m, 2, -1
                sd = lcg(sd); j = 1 + mod(sd, i)
                tmp = rl(i); rl(i) = rl(j); rl(j) = tmp
            end do
        end if
        idx = 0
        do i = 1, m
            if (m == n) then
                c = i
            else
                c = 1 + nint(real(i-1,wp) * real(n-1,wp) / real(m-1,wp))
            end if
            j1 = max(1, min(n-nb+1, c-hb))
            k1(i) = idx + 1
            do k = 1, nb
                idx = idx + 1
                irow(idx) = rl(i)
                icol(idx) = j1 + k - 1
                ! the large coefficient is at the center column c (near the ends of the
                ! matrix the band is shifted to stay inside it, so c is not always the middle):
                d = abs(icol(idx) - c)
                if (d == 0) then
                    coef(idx) = d0 + sin(real(i,wp))
                    icub(i) = idx
                else
                    coef(idx) = -1.0_wp / real(d,wp) * (1.0_wp + 0.25_wp*cos(real(i+k,wp)))
                end if
            end do
            k2(i) = idx
            if (nonlinear) then
                isin(i) = k1(i)        ! first in band
            else
                icub(i) = 0; isin(i) = 0
            end if
        end do
        ccoef = merge(0.1_wp, 0.0_wp, nonlinear)
        scoef = merge(0.2_wp, 0.0_wp, nonlinear)
        allocate(xstar(n))
        do i = 1, n
            xstar(i) = sin(2.0_wp*pi*real(i,wp)/real(n,wp)) + 0.5_wp*cos(6.0_wp*pi*real(i,wp)/real(n,wp))
        end do
        call set_rhs_from(xstar)
    end subroutine make_banded

!******************************************************************************************************
    subroutine make_grid(k)
        !! 2-D grid-like problem: m = k*k equations (one per grid node), n = 2m variables
        !! (a state u and a control v at each node):
        !! f_e = 4 u_e - sum(neighbors u) + 0.1 u_e**3 + v_e - b_e   (5-point stencil)
        integer,intent(in) :: k
        integer :: i, j, e, idx
        real(wp),dimension(:),allocatable :: xstar
        m = k*k; n = 2*m; nnz = 6*m   ! upper bound
        call alloc_problem()
        rl = [(i, i=1,m)]
        idx = 0
        do j = 1, k
            do i = 1, k
                e = (j-1)*k + i
                k1(e) = idx + 1
                call add(e, e, 4.0_wp); icub(e) = idx
                if (i > 1) call add(e, e-1, -1.0_wp)
                if (i < k) call add(e, e+1, -1.0_wp)
                if (j > 1) call add(e, e-k, -1.0_wp)
                if (j < k) call add(e, e+k, -1.0_wp)
                call add(e, m+e, 1.0_wp)
                k2(e) = idx
                isin(e) = 0
            end do
        end do
        nnz = idx
        irow = irow(1:nnz); icol = icol(1:nnz); coef = coef(1:nnz)
        ccoef = 0.1_wp; scoef = 0.0_wp
        allocate(xstar(n))
        do i = 1, n
            xstar(i) = sin(real(i,wp)*0.01_wp)
        end do
        call set_rhs_from(xstar)
    contains
        subroutine add(row, col, v)
            integer,intent(in) :: row, col
            real(wp),intent(in) :: v
            idx = idx + 1
            irow(idx) = row; icol(idx) = col; coef(idx) = v
        end subroutine add
    end subroutine make_grid

!******************************************************************************************************
    subroutine make_pairs(mm, delta)
        !! ill-conditioned linear problem: rows 2p-1 and 2p are identical except for a
        !! perturbation delta in a separate column, so cond(A) ~ 1/delta.
        !! m = mm (even), n = 3*mm/2 + mm.
        integer,intent(in) :: mm
        real(wp),intent(in) :: delta
        integer :: i, p, idx
        m = mm; n = 3*(mm/2) + mm; nnz = 4*m
        call alloc_problem()
        rl = [(i, i=1,m)]
        idx = 0
        do i = 1, m
            p = (i+1)/2
            k1(i) = idx + 1
            irow(idx+1:idx+4) = i
            icol(idx+1:idx+4) = [3*p-2, 3*p-1, 3*p, 3*(m/2)+i]
            coef(idx+1:idx+4) = [1.0_wp+0.1_wp*sin(real(p,wp)), 2.0_wp, 1.0_wp, &
                                 merge(delta, 0.5_wp*delta, mod(i,2)==0)]
            idx = idx + 4
            k2(i) = idx
            icub(i) = 0; isin(i) = 0
            b(i) = cos(real(i,wp))
        end do
        ccoef = 0.0_wp; scoef = 0.0_wp
    end subroutine make_pairs

!******************************************************************************************************
    subroutine alloc_problem()
        if (allocated(irow)) deallocate(irow, icol, coef, k1, k2, rl, icub, isin, b)
        allocate(irow(nnz), icol(nnz), coef(nnz), k1(m), k2(m), rl(m), icub(m), isin(m), b(m))
    end subroutine alloc_problem

    integer function lcg(s)
        integer,intent(in) :: s
        integer(8) :: t
        t = mod(int(s,8)*16807_8, 2147483647_8)
        lcg = int(t)
    end function lcg

    end module bench_problem_module
!******************************************************************************************************

!******************************************************************************************************
    program nlesolver_benchmarks

    use nlesolver_module, wp => nlesolver_rk
    use nlesolver_cholesky_module
    use bench_problem_module
    use iso_fortran_env, only: int64, output_unit

    implicit none

    type :: result
        integer :: istat = 0, iters = 0, nfe = 0, nje = 0, inner = -1
        real(wp) :: fnorm = 0.0_wp, dist = 0.0_wp, time = 0.0_wp
    end type result

    integer,parameter :: n_solvers = 6
    integer,dimension(n_solvers),parameter :: solver_modes = [NLESOLVER_SPARSITY_DENSE, NLESOLVER_SPARSITY_LUSOL, &
                                                              NLESOLVER_SPARSITY_LSQR,  NLESOLVER_SPARSITY_LSMR,  &
                                                              NLESOLVER_SPARSITY_CRAIG, NLESOLVER_SPARSITY_CHOLESKY]
    character(len=8),dimension(n_solvers),parameter :: solver_names = ['dense   ','lusol   ','lsqr    ', &
                                                                        'lsmr    ','craig   ','cholesky']
    character(len=:),allocatable :: outdir
    character(len=:),allocatable :: studies  !! comma separated list of studies to run (default: all)
    real(wp),dimension(:),allocatable :: last_x  !! solution of the last single-iteration (linear) run
    integer :: nargs
    character(len=256) :: arg

    nargs = command_argument_count()
    if (nargs >= 1) then
        call get_command_argument(1, arg)
        outdir = trim(arg)
    else
        outdir = 'docs/data'
    end if
    if (nargs >= 2) then
        call get_command_argument(2, arg)
        studies = ','//trim(arg)//','
    else
        studies = ',spy,scaling,square,steps,cond,sparsity,'
    end if
    call execute_command_line('mkdir -p '//outdir)
    write(*,'(A)') 'writing results to: '//outdir

    if (index(studies,',spy,')      > 0) call study_spy()
    if (index(studies,',scaling,')  > 0) call study_scaling()
    if (index(studies,',square,')   > 0) call study_square()
    if (index(studies,',steps,')    > 0) call study_steps()
    if (index(studies,',cond,')     > 0) call study_cond()
    if (index(studies,',sparsity,') > 0) call study_sparsity()

    write(*,'(A)') 'done.'

    contains

!******************************************************************************************************
    function run(mode, step_mode, use_broyden, max_iter, tol, nrep, refine, x0) result(r)
        !! solve the current problem (best time of `nrep` runs)
        integer,intent(in) :: mode, step_mode, max_iter, nrep
        logical,intent(in) :: use_broyden
        real(wp),intent(in) :: tol
        integer,intent(in),optional :: refine
        real(wp),dimension(:),intent(in),optional :: x0
        type(result) :: r
        type(nlesolver_type) :: solver
        real(wp),dimension(:),allocatable :: x, xs, f
        integer(int64) :: c0, c1, cr
        integer :: irep, u, ios, nlines
        character(len=:),allocatable :: msg
        character(len=200) :: line
        logical :: count_inner

        allocate(x(n), xs(n), f(m))
        xs = 0.0_wp
        if (present(x0)) xs = x0
        r%time = huge(1.0_wp)
        count_inner = (mode == NLESOLVER_SPARSITY_CRAIG) .and. max_iter == 1
        do irep = 1, nrep
            nf = 0; ng = 0; last_iter = 0; nhist = 0
            u = 0
            if (count_inner) then
                u = 77
                open(unit=u, status='scratch')
            end if
            if (mode == NLESOLVER_SPARSITY_DENSE) then
                call solver%initialize(n=n, m=m, max_iter=max_iter, tol=tol, func=func, grad=grad_dense, &
                                       step_mode=step_mode, use_broyden=use_broyden, export_iteration=export, &
                                       n_intervals=4, fmin_tol=1.0e-4_wp)
            else
                call solver%initialize(n=n, m=m, max_iter=max_iter, tol=tol, func=func, grad_sparse=grad_sparse, &
                                       step_mode=step_mode, use_broyden=use_broyden, export_iteration=export, &
                                       n_intervals=4, fmin_tol=1.0e-4_wp, &
                                       sparsity_mode=mode, irow=irow, icol=icol, &
                                       atol=1.0e-12_wp, btol=1.0e-12_wp, itnlim=100000, nout=u, &
                                       cholesky_refine=refine)
            end if
            x = xs
            call system_clock(c0, cr)
            call solver%solve(x)
            call system_clock(c1)
            r%time = min(r%time, real(c1-c0,wp)/real(cr,wp))
            if (count_inner) then
                rewind(u); nlines = -1
                do
                    read(u,'(A)',iostat=ios) line
                    if (ios /= 0) exit
                    nlines = nlines + 1
                end do
                close(u)
                r%inner = nlines
            end if
        end do
        call solver%status(r%istat, msg)
        call eval_f(x, f)
        r%iters = last_iter - 1
        r%nfe = nf
        r%nje = ng
        r%fnorm = norm2(f)
        r%dist = norm2(x - xs)
        if (max_iter == 1) last_x = x   ! for the linear studies
    end function run

!******************************************************************************************************
    subroutine write_header(u, cols)
        integer,intent(in) :: u
        character(len=*),intent(in) :: cols
        write(u,'(A)') cols
    end subroutine write_header

    subroutine write_result(u, prefix, r)
        integer,intent(in) :: u
        character(len=*),intent(in) :: prefix
        type(result),intent(in) :: r
        write(u,'(A,I0,",",I0,",",I0,",",I0,",",ES12.4,",",ES14.6,",",ES12.4)') &
            prefix, r%istat, r%iters, r%nfe, r%nje, r%fnorm, r%dist, r%time
        flush(u)
        flush(output_unit)
    end subroutine write_result

!******************************************************************************************************
    subroutine study_scaling()
        !! banded underdetermined (m = 0.6 n): time vs n for each solver
        integer,dimension(9),parameter :: ns = [1000, 2000, 5000, 10000, 20000, 50000, 100000, 200000, 500000]
        integer :: in, is, ud, uw
        type(result) :: r
        real(wp),dimension(n_solvers) :: t
        character(len=32) :: pre
        write(*,'(/A)') '--- scaling (underdetermined banded, m = 0.6 n) ---'
        open(newunit=ud, file=outdir//'/scaling_detail.csv', status='replace')
        open(newunit=uw, file=outdir//'/scaling_time.csv', status='replace')
        call write_header(ud, 'n,m,solver,istat,iters,fevals,jevals,fnorm,dist,time')
        call write_header(uw, 'n,dense,lusol,lsqr,lsmr,craig,cholesky')
        do in = 1, size(ns)
            call make_banded(ns(in), 3*ns(in)/5, 4, 4.0_wp, .false., .true.)
            t = -1.0_wp
            do is = 1, n_solvers
                if (solver_modes(is) == NLESOLVER_SPARSITY_DENSE .and. n > 5000) cycle
                ! LUSOL does not handle m < n (it reports the system as singular),
                ! so it is only run at the smallest size, to document that:
                if (solver_modes(is) == NLESOLVER_SPARSITY_LUSOL .and. n > 1000) cycle
                r = run(solver_modes(is), NLESOLVER_LINESEARCH_BACKTRACKING, .false., 50, 1.0e-10_wp, &
                        merge(3, 1, n <= 50000))
                if (r%istat == 1) t(is) = r%time  ! (failed runs are not plotted)
                write(pre,'(I0,",",I0,",",A,",")') n, m, trim(solver_names(is))
                call write_result(ud, trim(pre), r)
                write(*,'(I8,1X,A8,I4,I4,ES12.4,ES12.4,F10.4)') n, solver_names(is), r%istat, r%iters, &
                                                                 r%fnorm, r%dist, r%time
            end do
            write(uw,'(I0,*(:,",",A))') n, (trim(tstr(t(is))), is=1,n_solvers)
            flush(uw)
        end do
        close(ud); close(uw)
    end subroutine study_scaling

    function tstr(t) result(s)
        !! time as a string (nan if not run)
        real(wp),intent(in) :: t
        character(len=:),allocatable :: s
        character(len=20) :: buf
        if (t < 0.0_wp) then
            s = 'nan'
        else
            write(buf,'(ES12.4)') t
            s = trim(adjustl(buf))
        end if
    end function tstr

!******************************************************************************************************
    subroutine study_square()
        !! banded square (m = n)
        integer,dimension(5),parameter :: ns = [1000, 5000, 20000, 100000, 500000]
        integer :: in, is, ud, uw
        type(result) :: r
        real(wp),dimension(n_solvers) :: t
        character(len=32) :: pre
        write(*,'(/A)') '--- square banded (m = n) ---'
        open(newunit=ud, file=outdir//'/square_detail.csv', status='replace')
        open(newunit=uw, file=outdir//'/square_time.csv', status='replace')
        call write_header(ud, 'n,m,solver,istat,iters,fevals,jevals,fnorm,dist,time')
        call write_header(uw, 'n,dense,lusol,lsqr,lsmr,craig,cholesky')
        do in = 1, size(ns)
            ! (diagonally dominant: d0 = 8 > sum of the off-diagonals)
            call make_banded(ns(in), ns(in), 4, 8.0_wp, .false., .true.)
            t = -1.0_wp
            do is = 1, n_solvers
                if (solver_modes(is) == NLESOLVER_SPARSITY_DENSE .and. n > 5000) cycle
                r = run(solver_modes(is), NLESOLVER_LINESEARCH_BACKTRACKING, .false., 50, 1.0e-10_wp, &
                        merge(3, 1, n <= 20000))
                if (r%istat == 1) t(is) = r%time  ! (failed runs are not plotted)
                write(pre,'(I0,",",I0,",",A,",")') n, m, trim(solver_names(is))
                call write_result(ud, trim(pre), r)
                write(*,'(I8,1X,A8,I4,I4,ES12.4,F10.4)') n, solver_names(is), r%istat, r%iters, r%fnorm, r%time
            end do
            write(uw,'(I0,*(:,",",A))') n, (trim(tstr(t(is))), is=1,n_solvers)
            flush(uw)
        end do
        close(ud); close(uw)
    end subroutine study_square

!******************************************************************************************************
    subroutine study_steps()
        !! step modes x Broyden on the n = 50,000 banded problem (CHOLESKY and CRAIG)
        integer :: ud, uh, sm, ib, is, i
        type(result) :: r
        logical :: br
        integer,dimension(2),parameter :: modes = [NLESOLVER_SPARSITY_CHOLESKY, NLESOLVER_SPARSITY_CRAIG]
        character(len=8),dimension(2),parameter :: names = ['cholesky','craig   ']
        character(len=64) :: pre, fname
        write(*,'(/A)') '--- step modes and Broyden (n = 50,000) ---'
        call make_banded(50000, 30000, 4, 4.0_wp, .false., .true.)
        open(newunit=ud, file=outdir//'/steps_detail.csv', status='replace')
        call write_header(ud, 'solver,step_mode,broyden,istat,iters,fevals,jevals,fnorm,dist,time')
        do is = 1, 2
            do ib = 0, 1
                br = ib == 1
                do sm = 1, 4
                    record_history = .true.
                    r = run(modes(is), sm, br, 50, 1.0e-10_wp, 3)
                    record_history = .false.
                    write(pre,'(A,",",I0,",",I0,",")') trim(names(is)), sm, ib
                    call write_result(ud, trim(pre), r)
                    write(*,'(A8,I3,L3,I4,I4,I6,I4,ES12.4,F10.4)') names(is), sm, br, r%istat, r%iters, &
                                                                    r%nfe, r%nje, r%fnorm, r%time
                    if (is == 1) then
                        write(fname,'(A,I0,A,I0,A)') '/hist_sm', sm, '_b', ib, '.csv'
                        open(newunit=uh, file=outdir//trim(fname), status='replace')
                        write(uh,'(A)') 'iter,fnorm'
                        do i = 1, nhist
                            write(uh,'(I0,",",ES12.4)') i-1, max(hist(i), 1.0e-16_wp)
                        end do
                        close(uh)
                    end if
                end do
            end do
        end do
        close(ud)
    end subroutine study_steps

!******************************************************************************************************
    subroutine study_cond()
        !! (a) CRAIG iterations vs. conditioning (banded, d0 knob), n = 50,000, one linear solve.
        !!     cond(A) is computed with the SVD of a small instance with the same structure.
        !! (b) Cholesky accuracy vs. cond(A) (row pairs, delta knob), with and without refinement.
        real(wp),dimension(7),parameter :: d0s = [8.0_wp, 4.0_wp, 2.0_wp, 1.0_wp, 0.5_wp, 0.3_wp, 0.0_wp]
        real(wp),dimension(9),parameter :: deltas = [1.0e-1_wp,1.0e-2_wp,1.0e-3_wp,1.0e-4_wp,1.0e-5_wp, &
                                                     1.0e-6_wp,1.0e-7_wp,1.0e-8_wp,1.0e-9_wp]
        integer :: i, ua, ub
        type(result) :: rc, rh
        real(wp) :: cnd
        real(wp),dimension(:),allocatable :: xd, xc, x0, x2
        write(*,'(/A)') '--- conditioning ---'
        open(newunit=ua, file=outdir//'/cond_craig.csv', status='replace')
        write(ua,'(A)') 'd0,cond,craig_iters,craig_time,chol_time,craig_res,chol_res'
        do i = 1, size(d0s)
            call make_banded(1000, 600, 4, d0s(i), .false., .false.)
            cnd = cond_dense()
            call make_banded(50000, 30000, 4, d0s(i), .false., .false.)
            rc = run(NLESOLVER_SPARSITY_CRAIG,    1, .false., 1, 1.0e-30_wp, 3)
            rh = run(NLESOLVER_SPARSITY_CHOLESKY, 1, .false., 1, 1.0e-30_wp, 3)
            write(ua,'(F6.2,",",ES12.4,",",I0,",",ES12.4,",",ES12.4,",",ES12.4,",",ES12.4)') &
                d0s(i), cnd, rc%inner, rc%time, rh%time, rc%fnorm, rh%fnorm
            write(*,'(A,F6.2,A,ES10.2,A,I6,2F9.4,2ES11.3)') ' d0=', d0s(i), ' cond=', cnd, ' craig its=', &
                rc%inner, rc%time, rh%time, rc%fnorm, rh%fnorm
        end do
        close(ua)

        open(newunit=ub, file=outdir//'/cond_accuracy.csv', status='replace')
        write(ub,'(A)') 'delta,cond,err_craig,err_chol0,err_chol2,res_craig,res_chol0,res_chol2'
        do i = 1, size(deltas)
            call make_pairs(400, deltas(i))
            cnd = cond_dense()
            xd = minnorm_dense()
            rc = run(NLESOLVER_SPARSITY_CRAIG,    1, .false., 1, 1.0e-30_wp, 1); xc = last_x
            rh = run(NLESOLVER_SPARSITY_CHOLESKY, 1, .false., 1, 1.0e-30_wp, 1, refine=0); x0 = last_x
            rh = run(NLESOLVER_SPARSITY_CHOLESKY, 1, .false., 1, 1.0e-30_wp, 1, refine=5); x2 = last_x
            write(ub,'(ES10.2,",",ES12.4,6(",",ES12.4))') deltas(i), cnd, &
                rel(xc,xd), rel(x0,xd), rel(x2,xd), resid(xc), resid(x0), resid(x2)
            write(*,'(A,ES9.1,A,ES10.2,A,3ES11.3)') ' delta=', deltas(i), ' cond=', cnd, &
                '  err craig/chol0/chol5 =', rel(xc,xd), rel(x0,xd), rel(x2,xd)
        end do
        close(ub)
    end subroutine study_cond

    real(wp) function rel(x, xr)
        real(wp),dimension(:),intent(in) :: x, xr
        rel = norm2(x-xr)/norm2(xr)
    end function rel

    real(wp) function resid(x)
        real(wp),dimension(:),intent(in) :: x
        real(wp),dimension(m) :: f
        call eval_f(x, f)
        resid = norm2(f)/norm2(b)
    end function resid

    function dense_a() result(a)
        real(wp),dimension(m,n) :: a
        integer :: k
        a = 0.0_wp
        do k = 1, nnz
            a(irow(k),icol(k)) = a(irow(k),icol(k)) + coef(k)
        end do
    end function dense_a

    real(wp) function cond_dense()
        !! 2-norm condition number of the (linear part of the) current A, via the SVD
        interface
            subroutine dgesvd(jobu,jobvt,m,n,a,lda,s,u,ldu,vt,ldvt,work,lwork,info)
                import :: wp
                character :: jobu, jobvt
                integer :: m, n, lda, ldu, ldvt, lwork, info
                real(wp) :: a(lda,*), s(*), u(ldu,*), vt(ldvt,*), work(*)
            end subroutine dgesvd
        end interface
        real(wp),dimension(:,:),allocatable :: a
        real(wp),dimension(:),allocatable :: s, work
        real(wp) :: du(1,1), dvt(1,1)
        integer :: info
        a = dense_a()
        allocate(s(min(m,n)), work(10*(m+n)))
        call dgesvd('N','N', m, n, a, m, s, du, 1, dvt, 1, work, size(work), info)
        cond_dense = s(1)/s(min(m,n))
    end function cond_dense

    function minnorm_dense() result(x)
        !! dense minimum-norm solution of A x = b (LAPACK dgels)
        interface
            subroutine dgels(trans,m,n,nrhs,a,lda,b,ldb,work,lwork,info)
                import :: wp
                character :: trans
                integer :: info, lda, ldb, lwork, m, n, nrhs
                real(wp) :: a(lda,*), b(ldb,*), work(*)
            end subroutine dgels
        end interface
        real(wp),dimension(n) :: x
        real(wp),dimension(:,:),allocatable :: a, bb
        real(wp),dimension(:),allocatable :: work
        integer :: info
        a = dense_a()
        allocate(bb(max(m,n),1), work(64*(m+n)))
        bb = 0.0_wp
        bb(1:m,1) = b
        call dgels('N', m, n, 1, a, m, bb, max(m,n), work, size(work), info)
        x = bb(1:n,1)
    end function minnorm_dense

!******************************************************************************************************
    subroutine study_sparsity()
        !! banded (natural and shuffled) and 2-D grid-like problems:
        !! Cholesky envelope size and time, compared with CRAIG.
        integer,dimension(5),parameter :: ks = [25, 50, 100, 150, 200]
        integer :: u, i
        write(*,'(/A)') '--- sparsity ---'
        open(newunit=u, file=outdir//'/sparsity.csv', status='replace')
        write(u,'(A)') 'case,param,n,m,nnz,bw_natural,bw_rcm,envelope,chol_time,craig_time,chol_iters,craig_iters'
        call make_banded(50000, 30000, 4, 4.0_wp, .false., .true.); call sparsity_case(u, 'banded', 0)
        call make_banded(50000, 30000, 4, 4.0_wp, .true., .true.);  call sparsity_case(u, 'shuffled', 0)
        do i = 1, size(ks)
            call make_grid(ks(i)); call sparsity_case(u, 'grid', ks(i))
        end do
        close(u)
    end subroutine study_sparsity

        subroutine sparsity_case(u, name, param)
            integer,intent(in) :: u
            character(len=*),intent(in) :: name
            integer,intent(in) :: param
            type(cholesky_minnorm_solver) :: ch
            type(result) :: rh, rc
            integer :: istat, bwn
            call ch%analyze(m, n, irow, icol, istat)
            bwn = natural_bandwidth()
            rh = run(NLESOLVER_SPARSITY_CHOLESKY, NLESOLVER_LINESEARCH_BACKTRACKING, .false., 50, 1.0e-10_wp, 1)
            rc = run(NLESOLVER_SPARSITY_CRAIG,    NLESOLVER_LINESEARCH_BACKTRACKING, .false., 50, 1.0e-10_wp, 1)
            write(u,'(A,",",I0,",",I0,",",I0,",",I0,",",I0,",",I0,",",I0,",",ES12.4,",",ES12.4,",",I0,",",I0)') &
                name, param, n, m, nnz, bwn, ch%bandwidth(), ch%envelope_size(), rh%time, rc%time, rh%iters, rc%iters
            write(*,'(A10,I5,I8,I8,I8,I8,I12,2F10.4,2I4)') name, param, n, m, bwn, ch%bandwidth(), &
                ch%envelope_size(), rh%time, rc%time, rh%iters, rc%iters
        end subroutine sparsity_case

    integer function natural_bandwidth()
        !! bandwidth of A*A^T in the natural (given) row order
        integer,dimension(:),allocatable :: cmin, cmax
        integer :: k, i
        ! rows sharing a column: bandwidth = max over columns of (max row - min row)
        allocate(cmin(n), cmax(n))
        cmin = huge(1); cmax = 0
        do k = 1, nnz
            cmin(icol(k)) = min(cmin(icol(k)), irow(k))
            cmax(icol(k)) = max(cmax(icol(k)), irow(k))
        end do
        natural_bandwidth = 0
        do i = 1, n
            if (cmax(i) > 0) natural_bandwidth = max(natural_bandwidth, cmax(i)-cmin(i))
        end do
    end function natural_bandwidth

!******************************************************************************************************
    subroutine study_spy()
        !! small sparsity patterns for the figures
        write(*,'(/A)') '--- sparsity patterns ---'
        call make_banded(60, 36, 4, 4.0_wp, .false., .true.)
        call spy('band', .false.)
        call make_banded(60, 36, 4, 4.0_wp, .true., .true.)
        call spy('shuffled', .true.)
        call make_grid(6)
        call spy('grid', .true.)
    end subroutine study_spy

    subroutine spy(name, with_rcm)
        !! write the patterns of A, A*A^T, and (with_rcm) the RCM-reordered A*A^T and its envelope
        character(len=*),intent(in) :: name
        logical,intent(in) :: with_rcm
        logical,dimension(:,:),allocatable :: pa, paat
        integer,dimension(:),allocatable :: perm, iperm
        type(cholesky_minnorm_solver) :: ch
        integer :: k, i, j, u, istat, r
        allocate(pa(m,n), paat(m,m))
        pa = .false.
        do k = 1, nnz
            pa(irow(k),icol(k)) = .true.
        end do
        do i = 1, m
            do j = 1, m
                paat(i,j) = any(pa(i,:) .and. pa(j,:))
            end do
        end do
        open(newunit=u, file=outdir//'/spy_'//name//'_A.csv', status='replace')
        write(u,'(A)') 'row,col'
        do i = 1, m
            do j = 1, n
                if (pa(i,j)) write(u,'(I0,",",I0)') i, j
            end do
        end do
        close(u)
        open(newunit=u, file=outdir//'/spy_'//name//'_AAt.csv', status='replace')
        write(u,'(A)') 'row,col'
        do i = 1, m
            do j = 1, m
                if (paat(i,j)) write(u,'(I0,",",I0)') i, j
            end do
        end do
        close(u)
        if (with_rcm) then
            call ch%analyze(m, n, irow, icol, istat)
            call ch%get_permutation(perm)
            allocate(iperm(m))
            do r = 1, m
                iperm(perm(r)) = r
            end do
            open(newunit=u, file=outdir//'/spy_'//name//'_AAt_rcm.csv', status='replace')
            write(u,'(A)') 'row,col'
            do i = 1, m
                do j = 1, m
                    if (paat(perm(i),perm(j))) write(u,'(I0,",",I0)') i, j
                end do
            end do
            close(u)
            ! envelope boundary (lower triangle): first column in each row
            open(newunit=u, file=outdir//'/spy_'//name//'_env.csv', status='replace')
            write(u,'(A)') 'row,first'
            do i = 1, m
                do j = 1, i
                    if (paat(perm(i),perm(j))) exit
                end do
                write(u,'(I0,",",I0)') i, j
            end do
            close(u)
        end if
        write(*,'(A,A,I4,A,I4,A,I5)') name, ': m=', m, ' n=', n, ' nnz=', nnz
    end subroutine spy

    end program nlesolver_benchmarks
!******************************************************************************************************
