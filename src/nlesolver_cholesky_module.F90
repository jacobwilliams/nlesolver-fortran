!*****************************************************************************************
!>
!  Sparse direct minimum-norm linear solver, using a Cholesky factorization
!  of \( A A^T \) (no external dependencies).
!
!  Computes the minimum-norm solution of the compatible system \( A x = b \),
!  where `A` is an `m x n` sparse matrix (normally `n >= m`) in coordinate (triplet) format:
!
!  \( x = A^T y, \quad (A A^T + \delta^2 I) y = b \)
!
!  With `damp` \( \delta > 0 \), this is the solution of the regularized problem
!  \( \min \| x \|^2 + \| s \|^2 \) subject to \( A x + \delta s = b \)
!  (useful if `A` is rank deficient).
!
!### Method
!
!  1. **Analysis** (done once for a sparsity pattern, see [[cholesky_minnorm_solver:analyze]]):
!     the graph of \( A A^T \) is formed, and its rows are reordered with the reverse
!     Cuthill-McKee (RCM) algorithm to reduce the envelope (profile) of the matrix.
!     The envelope storage of the Cholesky factor is then allocated. Since all fill-in
!     of a Cholesky factorization is confined to the envelope, no further symbolic work
!     is needed.
!  2. **Factorization** (each time the values of `A` change, see [[cholesky_minnorm_solver:factor]]):
!     \( A A^T + \delta^2 I \) is formed directly in the envelope storage from the sparse
!     rows of `A` (the dense matrix is never formed), and factored as \( L L^T \)
!     using a row-oriented envelope Cholesky factorization.
!  3. **Solve** (see [[cholesky_minnorm_solver:solve]]): triangular solves with `L`, followed
!     by `x = A^T y`. Since forming \( A A^T \) squares the condition number, a few steps
!     of iterative refinement are performed using the residual of the original system.
!     Each correction stays in the range of \( A^T \), so the minimum-norm property is kept.
!
!  The cost of the factorization is \( O(m \, b^2) \), where `b` is the (average) bandwidth
!  of the reordered \( A A^T \). This is very efficient for banded or "mostly banded"
!  Jacobians (even if the rows are not numbered in banded order). For general sparse
!  matrices with large fill-in (e.g., 2D/3D meshes), a nested-dissection based solver
!  (such as qr_mumps) will be faster.
!
!### References
!  * A. George and J. W. H. Liu, "Computer Solution of Large Sparse Positive
!    Definite Systems", Prentice-Hall, 1981.
!  * E. Cuthill and J. McKee, "Reducing the bandwidth of sparse symmetric matrices",
!    Proc. 24th ACM National Conference, 1969.
!  * A. Bjorck, "Numerical Methods for Least Squares Problems", SIAM, 1996.

    module nlesolver_cholesky_module

    use iso_fortran_env

    implicit none

    private

#ifdef REAL32
    integer,parameter :: wp = real32   !! real kind used by this module [4 bytes]
#elif REAL64
    integer,parameter :: wp = real64   !! real kind used by this module [8 bytes]
#elif REAL128
    integer,parameter :: wp = real128  !! real kind used by this module [16 bytes]
#else
    integer,parameter :: wp = real64   !! real kind used by this module [8 bytes]
#endif

    integer,parameter :: ip = int64  !! integer kind for indexing the envelope storage

    real(wp),parameter :: zero = 0.0_wp
    real(wp),parameter :: eps  = epsilon(1.0_wp)

    ! status codes:
    integer,parameter,public :: CHOLESKY_SUCCESS        =  0 !! success
    integer,parameter,public :: CHOLESKY_ERR_INPUT      = -1 !! invalid inputs (e.g., `irow`, `icol` out of range)
    integer,parameter,public :: CHOLESKY_ERR_MEMORY     = -2 !! out of memory (envelope too large)
    integer,parameter,public :: CHOLESKY_ERR_NOT_PD     = -3 !! \( A A^T \) is not positive definite (`A` is rank deficient). Try `damp>0`.
    integer,parameter,public :: CHOLESKY_ERR_NOT_READY  = -4 !! `analyze` or `factor` has not been called

    type,public :: cholesky_minnorm_solver
        !! Sparse direct minimum-norm solver class.
        !!
        !!### Example
        !!```fortran
        !!  type(cholesky_minnorm_solver) :: solver
        !!  call solver%analyze(m, n, irow, icol, istat) ! once for the sparsity pattern
        !!  call solver%factor(a, damp, istat)           ! each time the values of A change
        !!  call solver%solve(b, x, istat)               ! for each right-hand side
        !!```
        private
        logical :: analyzed = .false. !! if `analyze` has been called
        logical :: factored = .false. !! if `factor` has been called
        integer :: m   = 0 !! number of rows in `A`
        integer :: n   = 0 !! number of columns in `A`
        integer :: nnz = 0 !! number of nonzeros in `A`
        real(wp) :: damp = zero !! damping parameter used in the factorization

        ! `A` stored by rows (CSR):
        integer,dimension(:),allocatable  :: rowptr !! row pointers (size `m+1`)
        integer,dimension(:),allocatable  :: colidx !! column indices (size `nnz`)
        integer,dimension(:),allocatable  :: map    !! index into the user's `a` array for each CSR element (size `nnz`)
        real(wp),dimension(:),allocatable :: val    !! values of `A` in CSR order (size `nnz`)

        ! graph of `A*A^T` (off-diagonal nonzeros), CSR:
        integer,dimension(:),allocatable :: adjptr !! (size `m+1`)
        integer,dimension(:),allocatable :: adj    !! (size `adjptr(m+1)-1`)

        ! row ordering:
        integer,dimension(:),allocatable :: perm  !! `perm(new) = old`
        integer,dimension(:),allocatable :: iperm !! `iperm(old) = new`

        ! envelope storage of the lower triangle of L (in the new ordering):
        ! row `r` holds columns `first(r):r`, stored contiguously starting at `ptr(r)`.
        integer,dimension(:),allocatable     :: first !! first column in each row (size `m`)
        integer(ip),dimension(:),allocatable :: ptr   !! start of each row in `env` (size `m+1`)
        real(wp),dimension(:),allocatable    :: env   !! envelope values

    contains

        private

        procedure,public :: analyze
        procedure,public :: factor
        procedure,public :: solve
        procedure,public :: destroy
        procedure,public :: is_analyzed
        procedure,public :: envelope_size
        procedure,public :: bandwidth
        procedure,public :: get_permutation

        procedure :: lower_solve

    end type cholesky_minnorm_solver

    contains
!*****************************************************************************************

!*****************************************************************************************
!>
!  Destroy the solver (deallocate all the arrays).

    subroutine destroy(me)

    class(cholesky_minnorm_solver),intent(out) :: me

    end subroutine destroy
!*****************************************************************************************

!*****************************************************************************************
!>
!  Returns true if [[cholesky_minnorm_solver:analyze]] has been called successfully.

    pure logical function is_analyzed(me)

    class(cholesky_minnorm_solver),intent(in) :: me

    is_analyzed = me%analyzed

    end function is_analyzed
!*****************************************************************************************

!*****************************************************************************************
!>
!  Number of elements in the envelope of the Cholesky factor
!  (i.e., the storage required for `L`, in reals). Returns `0` if not analyzed.

    pure integer(int64) function envelope_size(me)

    class(cholesky_minnorm_solver),intent(in) :: me

    if (me%analyzed) then
        envelope_size = me%ptr(me%m+1) - 1
    else
        envelope_size = 0
    end if

    end function envelope_size
!*****************************************************************************************

!*****************************************************************************************
!>
!  Maximum bandwidth of the reordered \( A A^T \). Returns `-1` if not analyzed.

    pure integer function bandwidth(me)

    class(cholesky_minnorm_solver),intent(in) :: me

    integer :: r

    if (me%analyzed) then
        bandwidth = 0
        do r = 1, me%m
            bandwidth = max(bandwidth, r - me%first(r))
        end do
    else
        bandwidth = -1
    end if

    end function bandwidth
!*****************************************************************************************

!*****************************************************************************************
!>
!  Returns the row ordering computed by [[cholesky_minnorm_solver:analyze]]:
!  `perm(new) = old`. Returns an unallocated array if not analyzed.

    subroutine get_permutation(me, perm)

    class(cholesky_minnorm_solver),intent(in) :: me
    integer,dimension(:),allocatable,intent(out) :: perm

    if (me%analyzed) perm = me%perm

    end subroutine get_permutation
!*****************************************************************************************

!*****************************************************************************************
!>
!  Symbolic analysis of the sparsity pattern of `A`.
!  This only has to be called once for a given sparsity pattern.
!
!  Forms the graph of \( A A^T \), computes a reverse Cuthill-McKee ordering
!  of the rows, and allocates the envelope storage for the Cholesky factor.
!
!@note Duplicate `(irow,icol)` entries are allowed (they are summed).

    subroutine analyze(me, m, n, irow, icol, istat)

    class(cholesky_minnorm_solver),intent(inout) :: me
    integer,intent(in)               :: m      !! number of rows in `A`
    integer,intent(in)               :: n      !! number of columns in `A`
    integer,dimension(:),intent(in)  :: irow   !! row indices of the nonzero elements of `A`
    integer,dimension(:),intent(in)  :: icol   !! column indices of the nonzero elements of `A`
                                               !! (same size as `irow`)
    integer,intent(out)              :: istat  !! status code (`CHOLESKY_SUCCESS` if no errors)

    integer,dimension(:),allocatable :: colptr  !! `A` stored by columns: column pointers
    integer,dimension(:),allocatable :: rowidx  !! `A` stored by columns: row indices
    integer,dimension(:),allocatable :: mark    !! marker array
    integer,dimension(:),allocatable :: cnt     !! counter array
    integer :: i, j, k, kk, t, c, r, s, nadj, alloc_stat

    call me%destroy()
    istat = CHOLESKY_SUCCESS

    ! check inputs:
    if (m < 1 .or. n < 1 .or. size(irow) /= size(icol) .or. size(irow) < 1) then
        istat = CHOLESKY_ERR_INPUT; return
    end if
    if (any(irow < 1) .or. any(irow > m) .or. any(icol < 1) .or. any(icol > n)) then
        istat = CHOLESKY_ERR_INPUT; return
    end if

    me%m   = m
    me%n   = n
    me%nnz = size(irow)

    ! A by rows (CSR). the map array keeps track of the original element index:
    allocate(me%rowptr(m+1), me%colidx(me%nnz), me%map(me%nnz), me%val(me%nnz), cnt(max(m,n)))
    cnt(1:m) = 0
    do k = 1, me%nnz
        cnt(irow(k)) = cnt(irow(k)) + 1
    end do
    me%rowptr(1) = 1
    do i = 1, m
        me%rowptr(i+1) = me%rowptr(i) + cnt(i)
    end do
    cnt(1:m) = me%rowptr(1:m)  ! insertion points
    do k = 1, me%nnz
        i = irow(k)
        me%colidx(cnt(i)) = icol(k)
        me%map(cnt(i))    = k
        cnt(i) = cnt(i) + 1
    end do

    ! A by columns (only the row indices are needed):
    allocate(colptr(n+1), rowidx(me%nnz))
    cnt(1:n) = 0
    do k = 1, me%nnz
        cnt(icol(k)) = cnt(icol(k)) + 1
    end do
    colptr(1) = 1
    do j = 1, n
        colptr(j+1) = colptr(j) + cnt(j)
    end do
    cnt(1:n) = colptr(1:n)
    do k = 1, me%nnz
        j = icol(k)
        rowidx(cnt(j)) = irow(k)
        cnt(j) = cnt(j) + 1
    end do

    ! graph of A*A^T: rows i and j are adjacent if they share a column.
    ! two passes: count, then fill.
    allocate(mark(m), me%adjptr(m+1))
    mark = 0
    me%adjptr(1) = 1
    do i = 1, m
        mark(i) = i  ! exclude the diagonal
        nadj = 0
        do kk = me%rowptr(i), me%rowptr(i+1)-1
            c = me%colidx(kk)
            do t = colptr(c), colptr(c+1)-1
                j = rowidx(t)
                if (mark(j) /= i) then
                    mark(j) = i
                    nadj = nadj + 1
                end if
            end do
        end do
        me%adjptr(i+1) = me%adjptr(i) + nadj
    end do
    allocate(me%adj(me%adjptr(m+1)-1), stat=alloc_stat)
    if (alloc_stat /= 0) then
        call me%destroy(); istat = CHOLESKY_ERR_MEMORY; return
    end if
    mark = 0
    do i = 1, m
        mark(i) = i
        s = me%adjptr(i)
        do kk = me%rowptr(i), me%rowptr(i+1)-1
            c = me%colidx(kk)
            do t = colptr(c), colptr(c+1)-1
                j = rowidx(t)
                if (mark(j) /= i) then
                    mark(j) = i
                    me%adj(s) = j
                    s = s + 1
                end if
            end do
        end do
    end do
    deallocate(colptr, rowidx, mark, cnt)

    ! reverse Cuthill-McKee ordering:
    allocate(me%perm(m), me%iperm(m))
    call rcm(m, me%adjptr, me%adj, me%perm)
    do r = 1, m
        me%iperm(me%perm(r)) = r
    end do

    ! envelope of the reordered matrix:
    allocate(me%first(m), me%ptr(m+1))
    do r = 1, m
        me%first(r) = r
    end do
    do i = 1, m
        r = me%iperm(i)
        do t = me%adjptr(i), me%adjptr(i+1)-1
            s = me%iperm(me%adj(t))
            if (s < me%first(r)) me%first(r) = s
        end do
    end do
    me%ptr(1) = 1
    do r = 1, m
        me%ptr(r+1) = me%ptr(r) + int(r - me%first(r) + 1, ip)
    end do
    allocate(me%env(me%ptr(m+1)-1), stat=alloc_stat)
    if (alloc_stat /= 0) then
        call me%destroy(); istat = CHOLESKY_ERR_MEMORY; return
    end if

    me%analyzed = .true.

    end subroutine analyze
!*****************************************************************************************

!*****************************************************************************************
!>
!  Reverse Cuthill-McKee ordering of a symmetric graph.
!
!  Each connected component is ordered by a breadth-first search starting from a
!  pseudo-peripheral node (found with the George-Liu algorithm), visiting the neighbors
!  of each node in order of increasing degree. The final ordering is reversed.

    subroutine rcm(m, adjptr, adj, perm)

    integer,intent(in)                :: m       !! number of nodes
    integer,dimension(:),intent(in)   :: adjptr  !! adjacency pointers (size `m+1`)
    integer,dimension(:),intent(in)   :: adj     !! adjacency lists
    integer,dimension(m),intent(out)  :: perm    !! the ordering: `perm(new) = old`

    integer,dimension(:),allocatable :: deg      !! degree of each node
    logical,dimension(:),allocatable :: visited  !! if the node has been ordered
    integer,dimension(:),allocatable :: order    !! Cuthill-McKee ordering
    integer,dimension(:),allocatable :: q        !! queue for the level structure BFS
    integer,dimension(:),allocatable :: bfsmark  !! marker for the level structure BFS
    integer :: stamp      !! current marker value for `bfsmark`
    integer :: nord       !! number of nodes ordered so far
    integer :: next       !! search pointer for the next unvisited node
    integer :: head, start, v, w, t, i, s, y, h, hy, lstart, lend, iter

    allocate(deg(m), visited(m), order(m), q(m), bfsmark(m))
    do i = 1, m
        deg(i) = adjptr(i+1) - adjptr(i)
    end do
    visited = .false.
    bfsmark = 0
    stamp   = 0
    nord    = 0
    next    = 1

    do  ! loop over the connected components

        do while (next <= m)
            if (.not. visited(next)) exit
            next = next + 1
        end do
        if (next > m) exit
        s = next

        ! find a pseudo-peripheral node in this component:
        call level_structure(s, h, lstart, lend)
        do iter = 1, 50
            ! min degree node in the last level:
            y = q(lstart)
            do i = lstart+1, lend
                if (deg(q(i)) < deg(y)) y = q(i)
            end do
            call level_structure(y, hy, lstart, lend)
            if (hy > h) then
                s = y
                h = hy
            else
                exit
            end if
        end do

        ! Cuthill-McKee ordering of this component:
        nord = nord + 1
        order(nord) = s
        visited(s) = .true.
        head = nord
        do while (head <= nord)
            v = order(head)
            head = head + 1
            start = nord + 1
            do t = adjptr(v), adjptr(v+1)-1
                w = adj(t)
                if (.not. visited(w)) then
                    visited(w) = .true.
                    nord = nord + 1
                    order(nord) = w
                end if
            end do
            call sort_by_degree(order(start:nord))
        end do

    end do

    ! reverse:
    do i = 1, m
        perm(i) = order(m-i+1)
    end do

    contains

        subroutine level_structure(root, nlevels, last_start, last_end)
            !! breadth-first level structure rooted at `root`.
            !! returns the number of levels and the range of the last level in `q`.
            integer,intent(in)  :: root
            integer,intent(out) :: nlevels, last_start, last_end
            integer :: qhead, qtail, lev_end, u, tt, z
            stamp = stamp + 1
            q(1) = root
            bfsmark(root) = stamp
            qhead = 1
            qtail = 1
            nlevels = 0
            do
                ! the current level is q(qhead:qtail)
                nlevels = nlevels + 1
                last_start = qhead
                last_end = qtail
                lev_end = qtail
                do while (qhead <= lev_end)
                    u = q(qhead)
                    qhead = qhead + 1
                    do tt = adjptr(u), adjptr(u+1)-1
                        z = adj(tt)
                        if (bfsmark(z) /= stamp) then
                            bfsmark(z) = stamp
                            qtail = qtail + 1
                            q(qtail) = z
                        end if
                    end do
                end do
                if (qtail == lev_end) exit ! no new level
            end do
        end subroutine level_structure

        subroutine sort_by_degree(list)
            !! insertion sort of a (short) list of nodes by increasing degree
            integer,dimension(:),intent(inout) :: list
            integer :: ii, jj, tmp
            do ii = 2, size(list)
                tmp = list(ii)
                jj = ii - 1
                do while (jj >= 1)
                    if (deg(list(jj)) <= deg(tmp)) exit
                    list(jj+1) = list(jj)
                    jj = jj - 1
                end do
                list(jj+1) = tmp
            end do
        end subroutine sort_by_degree

    end subroutine rcm
!*****************************************************************************************

!*****************************************************************************************
!>
!  Numeric factorization: form \( A A^T + \delta^2 I \) in the envelope storage
!  and compute its Cholesky factorization \( L L^T \).
!
!  [[cholesky_minnorm_solver:analyze]] must be called first. This must be called
!  each time the values of `A` change.

    subroutine factor(me, a, damp, istat)

    class(cholesky_minnorm_solver),intent(inout) :: me
    real(wp),dimension(:),intent(in) :: a      !! nonzero elements of `A` (in the same
                                               !! order as `irow` and `icol` in `analyze`)
    real(wp),intent(in)              :: damp   !! damping parameter \( \delta \) (`0` for none)
    integer,intent(out)              :: istat  !! status code (`CHOLESKY_SUCCESS` if no errors)

    real(wp),dimension(:),allocatable :: w     !! dense work vector (size `n`)
    real(wp),dimension(:),allocatable :: dB    !! diagonal of `A*A^T + damp^2*I` (new ordering)
    integer :: i, j, kk, t, r, q, fr, fq, k0
    integer(ip) :: pr, pq
    real(wp) :: s, d

    istat = CHOLESKY_SUCCESS
    me%factored = .false.
    if (.not. me%analyzed) then
        istat = CHOLESKY_ERR_NOT_READY; return
    end if
    if (size(a) /= me%nnz) then
        istat = CHOLESKY_ERR_INPUT; return
    end if

    me%damp = abs(damp)
    me%val = a(me%map)

    ! form A*A^T + damp^2*I in the envelope (new ordering).
    ! each element (r,q) with q<=r is computed once, from row r.
    allocate(w(me%n), dB(me%m))
    w = zero
    me%env = zero
    do i = 1, me%m
        r  = me%iperm(i)
        fr = me%first(r)
        pr = me%ptr(r)
        ! scatter row i of A into w:
        do kk = me%rowptr(i), me%rowptr(i+1)-1
            w(me%colidx(kk)) = w(me%colidx(kk)) + me%val(kk)
        end do
        ! diagonal:
        s = zero
        do kk = me%rowptr(i), me%rowptr(i+1)-1
            s = s + me%val(kk) * w(me%colidx(kk))
        end do
        s = s + me%damp**2
        me%env(pr + (r-fr)) = s
        dB(r) = s
        ! off-diagonals in the lower triangle:
        do t = me%adjptr(i), me%adjptr(i+1)-1
            j = me%adj(t)
            q = me%iperm(j)
            if (q < r) then
                s = zero
                do kk = me%rowptr(j), me%rowptr(j+1)-1
                    s = s + me%val(kk) * w(me%colidx(kk))
                end do
                me%env(pr + (q-fr)) = s
            end if
        end do
        ! clear w:
        do kk = me%rowptr(i), me%rowptr(i+1)-1
            w(me%colidx(kk)) = zero
        end do
    end do

    ! row-oriented envelope Cholesky factorization (in place):
    do r = 1, me%m
        fr = me%first(r)
        pr = me%ptr(r)
        do q = fr, r-1
            fq = me%first(q)
            pq = me%ptr(q)
            k0 = max(fr, fq)
            ! L(r,q) = ( B(r,q) - sum_k L(r,k)*L(q,k) ) / L(q,q)
            s = me%env(pr + (q-fr))
            if (k0 < q) s = s - dot_product(me%env(pr + (k0-fr) : pr + (q-1-fr)), &
                                            me%env(pq + (k0-fq) : pq + (q-1-fq)))
            me%env(pr + (q-fr)) = s / me%env(pq + (q-fq))
        end do
        ! L(r,r) = sqrt( B(r,r) - sum_k L(r,k)**2 )
        d = me%env(pr + (r-fr))
        if (fr < r) d = d - dot_product(me%env(pr : pr + (r-1-fr)), me%env(pr : pr + (r-1-fr)))
        if (.not. (d > 4.0_wp * eps * dB(r))) then  ! (also catches NaN)
            istat = CHOLESKY_ERR_NOT_PD
            return
        end if
        me%env(pr + (r-fr)) = sqrt(d)
    end do

    me%factored = .true.

    end subroutine factor
!*****************************************************************************************

!*****************************************************************************************
!>
!  Solve \( (A A^T + \delta^2 I) y = c \) using the Cholesky factor.
!  `c` and `y` are in the original row ordering.

    subroutine lower_solve(me, c, y)

    class(cholesky_minnorm_solver),intent(in) :: me
    real(wp),dimension(me%m),intent(in)  :: c
    real(wp),dimension(me%m),intent(out) :: y

    real(wp),dimension(:),allocatable :: z
    integer :: r, fr
    integer(ip) :: pr

    ! permute:
    allocate(z(me%m))
    z = c(me%perm)

    ! forward: L z = c
    do r = 1, me%m
        fr = me%first(r)
        pr = me%ptr(r)
        if (fr < r) z(r) = z(r) - dot_product(me%env(pr : pr + (r-1-fr)), z(fr:r-1))
        z(r) = z(r) / me%env(pr + (r-fr))
    end do

    ! backward: L^T y = z (column oriented, in place)
    do r = me%m, 1, -1
        fr = me%first(r)
        pr = me%ptr(r)
        z(r) = z(r) / me%env(pr + (r-fr))
        if (fr < r) z(fr:r-1) = z(fr:r-1) - z(r) * me%env(pr : pr + (r-1-fr))
    end do

    ! unpermute:
    y(me%perm) = z

    end subroutine lower_solve
!*****************************************************************************************

!*****************************************************************************************
!>
!  Solve for the minimum-norm solution: \( x = A^T y \), \( (A A^T + \delta^2 I) y = b \).
!
!  [[cholesky_minnorm_solver:factor]] must be called first.
!
!  Iterative refinement is used to recover the accuracy lost by forming \( A A^T \):
!  the residual \( r = b - A x - \delta^2 y \) is computed with the original matrix,
!  a correction \( \Delta y \) is computed with the factorization, and
!  \( x \leftarrow x + A^T \Delta y \). Note that `x` is accumulated directly (rather than
!  computed as \( A^T y \)): for an ill-conditioned `A`, `y` can be very large with
!  cancellation in \( A^T y \), while the corrections become small. The refinement stops
!  when the residual no longer decreases. Each correction is in the range of
!  \( A^T \), so the minimum-norm property is kept.
!
!  The refinement converges if \( \epsilon \, \mathrm{cond}(A)^2 \lesssim 1 \)
!  (i.e., \( \mathrm{cond}(A) \lesssim 10^8 \) in double precision).

    subroutine solve(me, b, x, istat, n_refine, rnorm)

    class(cholesky_minnorm_solver),intent(inout) :: me
    real(wp),dimension(:),intent(in)  :: b         !! right-hand side (size `m`)
    real(wp),dimension(:),intent(out) :: x         !! minimum-norm solution (size `n`)
    integer,intent(out)               :: istat     !! status code (`CHOLESKY_SUCCESS` if no errors)
    integer,intent(in),optional       :: n_refine  !! maximum number of iterative refinement steps
                                                   !! (default is 2)
    real(wp),intent(out),optional     :: rnorm     !! final residual norm \( \| b - A x - \delta^2 y \| \)

    real(wp),dimension(:),allocatable :: y, dy, r, xbest
    real(wp) :: rn, rn_best
    integer :: it, nref, i, kk

    istat = CHOLESKY_SUCCESS
    if (.not. me%factored) then
        istat = CHOLESKY_ERR_NOT_READY; return
    end if
    if (size(b) /= me%m .or. size(x) /= me%n) then
        istat = CHOLESKY_ERR_INPUT; return
    end if

    nref = 2
    if (present(n_refine)) nref = max(0, n_refine)

    allocate(y(me%m), dy(me%m), r(me%m))
    y = zero
    x = zero
    r = b
    rn_best = huge(1.0_wp)

    do it = 0, nref

        call me%lower_solve(r, dy)
        y = y + dy

        ! x = x + A^T dy
        do i = 1, me%m
            do kk = me%rowptr(i), me%rowptr(i+1)-1
                x(me%colidx(kk)) = x(me%colidx(kk)) + me%val(kk) * dy(i)
            end do
        end do

        ! r = b - A x - damp^2 y
        do i = 1, me%m
            r(i) = b(i) - me%damp**2 * y(i)
            do kk = me%rowptr(i), me%rowptr(i+1)-1
                r(i) = r(i) - me%val(kk) * x(me%colidx(kk))
            end do
        end do
        rn = norm2(r)

        if (it == 0 .or. rn < rn_best) then
            rn_best = rn
            if (it < nref) xbest = x  ! save in case the next step is worse
        else
            ! no improvement: revert to the previous solution
            x = xbest
            exit
        end if
        if (rn == zero) exit

    end do

    if (present(rnorm)) rnorm = rn_best

    end subroutine solve
!*****************************************************************************************

!*****************************************************************************************
    end module nlesolver_cholesky_module
!*****************************************************************************************
