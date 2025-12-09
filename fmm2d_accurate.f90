module fmm2d_accurate
  !============================================================================
  ! Accurate 2D FMM for 1/r potential
  ! Uses complex Laurent series expansion
  ! Target accuracy: 10^-12 with p_max=15-20
  !
  ! For 1/r in z=0 plane, we use the fact that in 2D:
  ! 1/|z-z'| can be expanded using complex Laurent series
  !
  ! Multipole: M_k = sum_i q_i * (z_i - z_c)^k / k
  ! Local:     L_k = coefficients for (z - z_c)^k expansion
  !
  ! This gives machine precision accuracy for sufficient p
  !============================================================================
  implicit none
  private

  public :: fmm_system, fmm_init, fmm_destroy
  public :: fmm_add_particle, fmm_build_tree, fmm_compute
  public :: compute_direct

  integer, parameter :: dp = selected_real_kind(15, 307)
  real(dp), parameter :: PI = 3.141592653589793238462643383279502884197_dp
  real(dp), parameter :: EPS = 1.0e-14_dp

  type :: particle
    real(dp) :: x, y, q
    real(dp) :: fx, fy, phi
    integer :: box_id
  end type

  type :: box
    real(dp) :: xc, yc, size
    integer :: level, parent, children(4)
    logical :: is_leaf
    integer :: n_part
    integer, allocatable :: parts(:)

    ! Complex moments: M_k for k=0..p
    complex(dp), allocatable :: M(:)  ! Multipole
    complex(dp), allocatable :: L(:)  ! Local

    ! Interaction lists
    integer, allocatable :: near_list(:), far_list(:)
  end type

  type :: fmm_system
    type(particle), allocatable :: particles(:)
    integer :: n_part, max_part

    type(box), allocatable :: boxes(:)
    integer :: n_box, max_box, root, max_lev

    real(dp) :: xmin, xmax, ymin, ymax
    integer :: p_max
    integer :: leaf_size
  end type

contains

  subroutine fmm_init(sys, xmin, xmax, ymin, ymax, p_max, max_part)
    type(fmm_system), intent(out) :: sys
    real(dp), intent(in) :: xmin, xmax, ymin, ymax
    integer, intent(in) :: p_max, max_part

    sys%xmin = xmin
    sys%xmax = xmax
    sys%ymin = ymin
    sys%ymax = ymax
    sys%p_max = p_max
    sys%max_part = max_part
    sys%leaf_size = 20
    sys%n_part = 0
    sys%n_box = 0

    allocate(sys%particles(max_part))
    sys%max_box = max(1000, max_part/10)
    allocate(sys%boxes(sys%max_box))
  end subroutine

  subroutine fmm_destroy(sys)
    type(fmm_system), intent(inout) :: sys
    integer :: i
    if (allocated(sys%particles)) deallocate(sys%particles)
    if (allocated(sys%boxes)) then
      do i = 1, sys%n_box
        if (allocated(sys%boxes(i)%parts)) deallocate(sys%boxes(i)%parts)
        if (allocated(sys%boxes(i)%M)) deallocate(sys%boxes(i)%M)
        if (allocated(sys%boxes(i)%L)) deallocate(sys%boxes(i)%L)
        if (allocated(sys%boxes(i)%near_list)) deallocate(sys%boxes(i)%near_list)
        if (allocated(sys%boxes(i)%far_list)) deallocate(sys%boxes(i)%far_list)
      end do
      deallocate(sys%boxes)
    end if
  end subroutine

  subroutine fmm_add_particle(sys, x, y, q)
    type(fmm_system), intent(inout) :: sys
    real(dp), intent(in) :: x, y, q
    sys%n_part = sys%n_part + 1
    sys%particles(sys%n_part)%x = x
    sys%particles(sys%n_part)%y = y
    sys%particles(sys%n_part)%q = q
    sys%particles(sys%n_part)%fx = 0.0_dp
    sys%particles(sys%n_part)%fy = 0.0_dp
    sys%particles(sys%n_part)%phi = 0.0_dp
  end subroutine

  subroutine fmm_build_tree(sys)
    type(fmm_system), intent(inout) :: sys
    integer, allocatable :: plist(:)
    integer :: i

    sys%n_box = 1
    sys%root = 1
    sys%max_lev = 0

    sys%boxes(1)%xc = 0.5_dp * (sys%xmin + sys%xmax)
    sys%boxes(1)%yc = 0.5_dp * (sys%ymin + sys%ymax)
    sys%boxes(1)%size = max(sys%xmax - sys%xmin, sys%ymax - sys%ymin)
    sys%boxes(1)%level = 0
    sys%boxes(1)%parent = 0
    sys%boxes(1)%children = 0

    allocate(plist(sys%n_part))
    do i = 1, sys%n_part
      plist(i) = i
    end do

    call subdivide(sys, 1, plist, sys%n_part, 0)
    deallocate(plist)
  end subroutine

  recursive subroutine subdivide(sys, bid, plist, np, lev)
    type(fmm_system), intent(inout) :: sys
    integer, intent(in) :: bid, np, lev
    integer, intent(in) :: plist(np)

    integer :: i, ip, cid, quad
    integer :: clists(4, sys%leaf_size*2), ccnts(4)
    real(dp) :: xc, yc, dx, dy

    if (np <= sys%leaf_size .or. lev >= 15) then
      sys%boxes(bid)%is_leaf = .true.
      sys%boxes(bid)%n_part = np
      allocate(sys%boxes(bid)%parts(np))
      sys%boxes(bid)%parts = plist(1:np)
      allocate(sys%boxes(bid)%M(0:sys%p_max))
      allocate(sys%boxes(bid)%L(0:sys%p_max))
      sys%boxes(bid)%M = cmplx(0.0_dp, 0.0_dp, dp)
      sys%boxes(bid)%L = cmplx(0.0_dp, 0.0_dp, dp)
      do i = 1, np
        sys%particles(plist(i))%box_id = bid
      end do
      if (lev > sys%max_lev) sys%max_lev = lev
      return
    end if

    sys%boxes(bid)%is_leaf = .false.
    xc = sys%boxes(bid)%xc
    yc = sys%boxes(bid)%yc

    ccnts = 0
    do i = 1, np
      ip = plist(i)
      dx = sys%particles(ip)%x - xc
      dy = sys%particles(ip)%y - yc
      quad = 1 + merge(0,1,dx<0) + merge(0,2,dy<0)
      ccnts(quad) = ccnts(quad) + 1
      clists(quad, ccnts(quad)) = ip
    end do

    do quad = 1, 4
      if (ccnts(quad) == 0) cycle
      sys%n_box = sys%n_box + 1
      cid = sys%n_box
      sys%boxes(bid)%children(quad) = cid
      sys%boxes(cid)%level = lev + 1
      sys%boxes(cid)%parent = bid
      sys%boxes(cid)%size = 0.5_dp * sys%boxes(bid)%size
      sys%boxes(cid)%children = 0
      dx = 0.25_dp * sys%boxes(bid)%size
      dy = 0.25_dp * sys%boxes(bid)%size
      sys%boxes(cid)%xc = xc + merge(-dx,dx,mod(quad-1,2)==0)
      sys%boxes(cid)%yc = yc + merge(-dy,dy,quad<=2)
      call subdivide(sys, cid, clists(quad,1:ccnts(quad)), ccnts(quad), lev+1)
    end do

    allocate(sys%boxes(bid)%M(0:sys%p_max))
    allocate(sys%boxes(bid)%L(0:sys%p_max))
    sys%boxes(bid)%M = cmplx(0.0_dp, 0.0_dp, dp)
    sys%boxes(bid)%L = cmplx(0.0_dp, 0.0_dp, dp)
  end subroutine

  !============================================================================
  ! Main FMM computation
  !============================================================================
  subroutine fmm_compute(sys)
    type(fmm_system), intent(inout) :: sys
    integer :: i

    ! Reset
    do i = 1, sys%n_part
      sys%particles(i)%fx = 0.0_dp
      sys%particles(i)%fy = 0.0_dp
      sys%particles(i)%phi = 0.0_dp
    end do

    ! Build interaction lists
    call build_lists(sys)

    ! P2M
    do i = 1, sys%n_box
      if (sys%boxes(i)%is_leaf) call p2m(sys, i)
    end do

    ! M2M upward
    call m2m_pass(sys)

    ! M2L interaction
    call m2l_pass(sys)

    ! L2L downward
    call l2l_pass(sys)

    ! L2P + P2P
    do i = 1, sys%n_box
      if (sys%boxes(i)%is_leaf) then
        call l2p(sys, i)
        call p2p_self(sys, i)
        call p2p_near(sys, i)
      end if
    end do
  end subroutine

  !============================================================================
  ! Build interaction lists
  !============================================================================
  subroutine build_lists(sys)
    type(fmm_system), intent(inout) :: sys
    integer :: i, j, nn, nf
    integer :: near_tmp(sys%n_box), far_tmp(sys%n_box)
    real(dp) :: dx, dy, dist, size_i, size_j

    ! Build lists for ALL boxes (not just leaves)
    do i = 1, sys%n_box
      nn = 0
      nf = 0

      do j = 1, sys%n_box
        if (i == j) cycle
        ! Only interact with boxes at same level
        if (sys%boxes(i)%level /= sys%boxes(j)%level) cycle

        dx = sys%boxes(j)%xc - sys%boxes(i)%xc
        dy = sys%boxes(j)%yc - sys%boxes(i)%yc
        dist = sqrt(dx**2 + dy**2)
        size_i = sys%boxes(i)%size
        size_j = sys%boxes(j)%size

        ! Well-separated: dist > size_i + size_j
        if (dist > 1.5_dp * (size_i + size_j)) then
          nf = nf + 1
          far_tmp(nf) = j
        else if (sys%boxes(i)%is_leaf .and. sys%boxes(j)%is_leaf) then
          ! Near list only for leaves
          nn = nn + 1
          near_tmp(nn) = j
        end if
      end do

      if (nn > 0) then
        allocate(sys%boxes(i)%near_list(nn))
        sys%boxes(i)%near_list = near_tmp(1:nn)
      end if

      if (nf > 0) then
        allocate(sys%boxes(i)%far_list(nf))
        sys%boxes(i)%far_list = far_tmp(1:nf)
      end if
    end do
  end subroutine

  !============================================================================
  ! P2M: Particle to Multipole
  ! M_k = sum q_i * z_i^k
  !============================================================================
  subroutine p2m(sys, bid)
    type(fmm_system), intent(inout) :: sys
    integer, intent(in) :: bid
    integer :: i, ip, k
    real(dp) :: dx, dy, q
    complex(dp) :: z, zk

    sys%boxes(bid)%M = cmplx(0.0_dp, 0.0_dp, dp)

    do i = 1, sys%boxes(bid)%n_part
      ip = sys%boxes(bid)%parts(i)
      q = sys%particles(ip)%q
      dx = sys%particles(ip)%x - sys%boxes(bid)%xc
      dy = sys%particles(ip)%y - sys%boxes(bid)%yc
      z = cmplx(dx, dy, dp)

      zk = cmplx(1.0_dp, 0.0_dp, dp)
      do k = 0, sys%p_max
        sys%boxes(bid)%M(k) = sys%boxes(bid)%M(k) + q * zk
        zk = zk * z
      end do
    end do
  end subroutine

  !============================================================================
  ! M2M: Multipole to Multipole (upward pass)
  ! M_parent(k) = sum M_child(k) + binomial translation
  !============================================================================
  subroutine m2m_pass(sys)
    type(fmm_system), intent(inout) :: sys
    integer :: lev, bid, cid, i

    do lev = sys%max_lev, 1, -1
      do bid = 1, sys%n_box
        if (sys%boxes(bid)%level /= lev-1 .or. sys%boxes(bid)%is_leaf) cycle

        do i = 1, 4
          cid = sys%boxes(bid)%children(i)
          if (cid > 0) call m2m(sys, cid, bid)
        end do
      end do
    end do
  end subroutine

  subroutine m2m(sys, cid, pid)
    type(fmm_system), intent(inout) :: sys
    integer, intent(in) :: cid, pid
    integer :: k, j
    real(dp) :: dx, dy
    complex(dp) :: z0, z0k
    complex(dp), allocatable :: M_trans(:)

    dx = sys%boxes(cid)%xc - sys%boxes(pid)%xc
    dy = sys%boxes(cid)%yc - sys%boxes(pid)%yc
    z0 = cmplx(dx, dy, dp)

    allocate(M_trans(0:sys%p_max))
    M_trans = cmplx(0.0_dp, 0.0_dp, dp)

    ! M_parent(k) += sum_{j=0}^k C(k,j) * z0^(k-j) * M_child(j)
    do k = 0, sys%p_max
      z0k = cmplx(1.0_dp, 0.0_dp, dp)
      do j = 0, k
        M_trans(k) = M_trans(k) + binomial(k,j) * z0k * sys%boxes(cid)%M(j)
        z0k = z0k * z0
      end do
    end do

    sys%boxes(pid)%M = sys%boxes(pid)%M + M_trans
    deallocate(M_trans)
  end subroutine

  !============================================================================
  ! M2L: Multipole to Local (interaction pass)
  !============================================================================
  subroutine m2l_pass(sys)
    type(fmm_system), intent(inout) :: sys
    integer :: bid, j, jid
    integer :: lev

    ! Do M2L for all boxes at each level
    do lev = 2, sys%max_lev
      do bid = 1, sys%n_box
        if (sys%boxes(bid)%level /= lev) cycle
        if (.not. allocated(sys%boxes(bid)%far_list)) cycle

        do j = 1, size(sys%boxes(bid)%far_list)
          jid = sys%boxes(bid)%far_list(j)
          call m2l(sys, jid, bid)
        end do
      end do
    end do
  end subroutine

  subroutine m2l(sys, sid, tid)
    type(fmm_system), intent(inout) :: sys
    integer, intent(in) :: sid, tid
    integer :: k, j
    real(dp) :: dx, dy
    complex(dp) :: z0, z0inv, z0k, term
    complex(dp), allocatable :: L_add(:)

    dx = sys%boxes(tid)%xc - sys%boxes(sid)%xc
    dy = sys%boxes(tid)%yc - sys%boxes(sid)%yc
    z0 = cmplx(dx, dy, dp)
    z0inv = 1.0_dp / z0

    allocate(L_add(0:sys%p_max))
    L_add = cmplx(0.0_dp, 0.0_dp, dp)

    ! Correct M2L formula for 2D FMM
    ! L_0 = M_0 * log(-z0) - sum_{j=1}^p M_j / j / z0^j
    ! L_k = -M_0 / k / z0^k - sum_{j=1}^p M_j * C(j+k-1,k) / z0^(j+k)  for k > 0

    ! k = 0 term
    L_add(0) = sys%boxes(sid)%M(0) * log(-z0)
    z0k = z0inv
    do j = 1, sys%p_max
      L_add(0) = L_add(0) - sys%boxes(sid)%M(j) * z0k / real(j, dp)
      z0k = z0k * z0inv
    end do

    ! k > 0 terms
    do k = 1, sys%p_max
      z0k = z0inv**k
      L_add(k) = -sys%boxes(sid)%M(0) * z0k / real(k, dp)

      do j = 1, sys%p_max
        term = -sys%boxes(sid)%M(j) * binomial(j+k-1, k) * z0k * z0inv**j
        L_add(k) = L_add(k) + term
      end do
    end do

    sys%boxes(tid)%L = sys%boxes(tid)%L + L_add
    deallocate(L_add)
  end subroutine

  !============================================================================
  ! L2L: Local to Local (downward pass)
  !============================================================================
  subroutine l2l_pass(sys)
    type(fmm_system), intent(inout) :: sys
    integer :: lev, bid, cid, i

    do lev = 1, sys%max_lev
      do bid = 1, sys%n_box
        if (sys%boxes(bid)%level /= lev-1) cycle

        do i = 1, 4
          cid = sys%boxes(bid)%children(i)
          if (cid > 0) call l2l(sys, bid, cid)
        end do
      end do
    end do
  end subroutine

  subroutine l2l(sys, pid, cid)
    type(fmm_system), intent(inout) :: sys
    integer, intent(in) :: pid, cid
    integer :: k, j
    real(dp) :: dx, dy
    complex(dp) :: z0, z0k
    complex(dp), allocatable :: L_trans(:)

    dx = sys%boxes(cid)%xc - sys%boxes(pid)%xc
    dy = sys%boxes(cid)%yc - sys%boxes(pid)%yc
    z0 = cmplx(dx, dy, dp)

    allocate(L_trans(0:sys%p_max))
    L_trans = sys%boxes(cid)%L

    ! L_child(k) += sum L_parent(j) * z0^(j-k) * C(j,k)
    do k = 0, sys%p_max
      z0k = cmplx(1.0_dp, 0.0_dp, dp)
      do j = k, sys%p_max
        L_trans(k) = L_trans(k) + sys%boxes(pid)%L(j) * z0k * binomial(j,k)
        z0k = z0k * z0
      end do
    end do

    sys%boxes(cid)%L = L_trans
    deallocate(L_trans)
  end subroutine

  !============================================================================
  ! L2P: Local to Particle
  !============================================================================
  subroutine l2p(sys, bid)
    type(fmm_system), intent(inout) :: sys
    integer, intent(in) :: bid
    integer :: i, ip, k
    real(dp) :: dx, dy
    complex(dp) :: z, zk, phi, dphidz

    do i = 1, sys%boxes(bid)%n_part
      ip = sys%boxes(bid)%parts(i)
      dx = sys%particles(ip)%x - sys%boxes(bid)%xc
      dy = sys%particles(ip)%y - sys%boxes(bid)%yc
      z = cmplx(dx, dy, dp)

      phi = cmplx(0.0_dp, 0.0_dp, dp)
      dphidz = cmplx(0.0_dp, 0.0_dp, dp)

      zk = cmplx(1.0_dp, 0.0_dp, dp)
      do k = 0, sys%p_max
        phi = phi + sys%boxes(bid)%L(k) * zk
        if (k > 0) dphidz = dphidz + real(k,dp) * sys%boxes(bid)%L(k) * zk / z
        zk = zk * z
      end do

      sys%particles(ip)%phi = sys%particles(ip)%phi + real(phi, dp)
      sys%particles(ip)%fx = sys%particles(ip)%fx - sys%particles(ip)%q * real(dphidz, dp)
      sys%particles(ip)%fy = sys%particles(ip)%fy - sys%particles(ip)%q * aimag(dphidz)
    end do
  end subroutine

  !============================================================================
  ! P2P: Particle to Particle (direct)
  !============================================================================
  subroutine p2p_self(sys, bid)
    type(fmm_system), intent(inout) :: sys
    integer, intent(in) :: bid
    call p2p_box(sys, bid, bid)
  end subroutine

  subroutine p2p_near(sys, bid)
    type(fmm_system), intent(inout) :: sys
    integer, intent(in) :: bid
    integer :: j, jid

    if (.not. allocated(sys%boxes(bid)%near_list)) return

    do j = 1, size(sys%boxes(bid)%near_list)
      jid = sys%boxes(bid)%near_list(j)
      call p2p_box(sys, bid, jid)
    end do
  end subroutine

  subroutine p2p_box(sys, bid1, bid2)
    type(fmm_system), intent(inout) :: sys
    integer, intent(in) :: bid1, bid2
    integer :: i, j, ip, jp
    real(dp) :: dx, dy, r, r3, qi, qj, fx, fy

    do i = 1, sys%boxes(bid1)%n_part
      ip = sys%boxes(bid1)%parts(i)
      qi = sys%particles(ip)%q

      do j = 1, sys%boxes(bid2)%n_part
        jp = sys%boxes(bid2)%parts(j)
        if (bid1 == bid2 .and. ip >= jp) cycle

        qj = sys%particles(jp)%q
        dx = sys%particles(jp)%x - sys%particles(ip)%x
        dy = sys%particles(jp)%y - sys%particles(ip)%y
        r = sqrt(dx**2 + dy**2)

        if (r < EPS) cycle

        r3 = r**3

        fx = qi * qj * dx / r3
        fy = qi * qj * dy / r3

        sys%particles(ip)%fx = sys%particles(ip)%fx + fx
        sys%particles(ip)%fy = sys%particles(ip)%fy + fy
        sys%particles(ip)%phi = sys%particles(ip)%phi + qj / r

        if (bid1 == bid2) then
          sys%particles(jp)%fx = sys%particles(jp)%fx - fx
          sys%particles(jp)%fy = sys%particles(jp)%fy - fy
          sys%particles(jp)%phi = sys%particles(jp)%phi + qi / r
        end if
      end do
    end do
  end subroutine

  !============================================================================
  ! Direct computation (reference)
  !============================================================================
  subroutine compute_direct(sys)
    type(fmm_system), intent(inout) :: sys
    integer :: i, j
    real(dp) :: dx, dy, r, r3, qi, qj, fx, fy

    do i = 1, sys%n_part
      sys%particles(i)%fx = 0.0_dp
      sys%particles(i)%fy = 0.0_dp
      sys%particles(i)%phi = 0.0_dp
    end do

    do i = 1, sys%n_part
      qi = sys%particles(i)%q
      do j = i+1, sys%n_part
        qj = sys%particles(j)%q
        dx = sys%particles(j)%x - sys%particles(i)%x
        dy = sys%particles(j)%y - sys%particles(i)%y
        r = sqrt(dx**2 + dy**2)

        if (r < EPS) cycle

        r3 = r**3

        fx = qi * qj * dx / r3
        fy = qi * qj * dy / r3

        sys%particles(i)%fx = sys%particles(i)%fx + fx
        sys%particles(i)%fy = sys%particles(i)%fy + fy
        sys%particles(i)%phi = sys%particles(i)%phi + qj / r

        sys%particles(j)%fx = sys%particles(j)%fx - fx
        sys%particles(j)%fy = sys%particles(j)%fy - fy
        sys%particles(j)%phi = sys%particles(j)%phi + qi / r
      end do
    end do
  end subroutine

  !============================================================================
  ! Helper: binomial coefficient
  !============================================================================
  pure function binomial(n, k) result(c)
    integer, intent(in) :: n, k
    real(dp) :: c
    integer :: i
    if (k > n .or. k < 0) then
      c = 0.0_dp
      return
    end if
    c = 1.0_dp
    do i = 0, min(k, n-k) - 1
      c = c * real(n-i, dp) / real(i+1, dp)
    end do
  end function

end module fmm2d_accurate
