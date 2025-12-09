module fmm_z0_precise
  !============================================================================
  ! High-precision FMM for 1/r potential in z=0 plane
  ! Optimized for z=0: uses 3D multipoles evaluated at z=0
  ! Target: 10^-12 accuracy with p_max=15-20
  !============================================================================
  implicit none
  private

  public :: fmm_system, fmm_init, fmm_destroy
  public :: fmm_add_particle, fmm_build_tree, fmm_compute_forces
  public :: fmm_compute_energy_direct

  integer, parameter :: dp = selected_real_kind(15, 307)
  real(dp), parameter :: PI = 3.141592653589793238462643383279502884197_dp

  type :: particle
    real(dp) :: x, y, q
    real(dp) :: fx, fy, phi
    integer :: box_id
  end type

  type :: box
    real(dp) :: xc, yc, size  ! Center and size
    integer :: level
    integer :: parent, children(4)  ! 2D quadtree
    logical :: is_leaf
    integer :: n_particles
    integer, allocatable :: particles(:)

    ! Multipole moments for 3D expansion at z=0
    ! M_n^m = sum_i q_i r_i^n (cos(m*phi_i) + i*sin(m*phi_i))
    real(dp), allocatable :: M_real(:,:)  ! (0:p, 0:p) real part
    real(dp), allocatable :: M_imag(:,:)  ! (0:p, 1:p) imag part

    ! Local expansion
    real(dp), allocatable :: L_real(:,:)
    real(dp), allocatable :: L_imag(:,:)
  end type

  type :: fmm_system
    ! Particles
    type(particle), allocatable :: particles(:)
    integer :: n_particles, max_particles

    ! Tree
    type(box), allocatable :: boxes(:)
    integer :: n_boxes, max_boxes
    integer :: root_id, max_level

    ! Domain
    real(dp) :: xmin, xmax, ymin, ymax

    ! Parameters
    integer :: p_max
    integer :: max_particles_per_box
    real(dp) :: mac_theta  ! Multipole acceptance criterion

    ! Precomputed tables
    real(dp), allocatable :: factorial(:)
    real(dp), allocatable :: binomial(:,:)
  end type

contains

  !============================================================================
  ! Initialize system
  !============================================================================
  subroutine fmm_init(sys, xmin, xmax, ymin, ymax, p_max, max_particles)
    type(fmm_system), intent(out) :: sys
    real(dp), intent(in) :: xmin, xmax, ymin, ymax
    integer, intent(in) :: p_max, max_particles

    integer :: i, j, n
    real(dp) :: f

    sys%xmin = xmin
    sys%xmax = xmax
    sys%ymin = ymin
    sys%ymax = ymax
    sys%p_max = p_max
    sys%max_particles = max_particles
    sys%max_particles_per_box = 20
    sys%mac_theta = 0.4_dp
    sys%n_particles = 0

    ! Allocate particles
    allocate(sys%particles(max_particles))

    ! Allocate tree (estimate)
    sys%max_boxes = max(1000, max_particles/5)
    allocate(sys%boxes(sys%max_boxes))
    sys%n_boxes = 0

    ! Precompute factorials
    allocate(sys%factorial(0:2*p_max+1))
    f = 1.0_dp
    sys%factorial(0) = 1.0_dp
    do n = 1, 2*p_max+1
      f = f * real(n, dp)
      sys%factorial(n) = f
    end do

    ! Precompute binomial coefficients C(n,k)
    allocate(sys%binomial(0:2*p_max, 0:2*p_max))
    sys%binomial = 0.0_dp
    sys%binomial(0,0) = 1.0_dp
    do n = 1, 2*p_max
      sys%binomial(n,0) = 1.0_dp
      sys%binomial(n,n) = 1.0_dp
      do i = 1, n-1
        sys%binomial(n,i) = sys%binomial(n-1,i-1) + sys%binomial(n-1,i)
      end do
    end do

  end subroutine fmm_init

  !============================================================================
  ! Destroy system
  !============================================================================
  subroutine fmm_destroy(sys)
    type(fmm_system), intent(inout) :: sys
    integer :: i

    if (allocated(sys%particles)) deallocate(sys%particles)
    if (allocated(sys%factorial)) deallocate(sys%factorial)
    if (allocated(sys%binomial)) deallocate(sys%binomial)

    if (allocated(sys%boxes)) then
      do i = 1, sys%n_boxes
        if (allocated(sys%boxes(i)%particles)) deallocate(sys%boxes(i)%particles)
        if (allocated(sys%boxes(i)%M_real)) deallocate(sys%boxes(i)%M_real)
        if (allocated(sys%boxes(i)%M_imag)) deallocate(sys%boxes(i)%M_imag)
        if (allocated(sys%boxes(i)%L_real)) deallocate(sys%boxes(i)%L_real)
        if (allocated(sys%boxes(i)%L_imag)) deallocate(sys%boxes(i)%L_imag)
      end do
      deallocate(sys%boxes)
    end if
  end subroutine fmm_destroy

  !============================================================================
  ! Add particle
  !============================================================================
  subroutine fmm_add_particle(sys, x, y, q)
    type(fmm_system), intent(inout) :: sys
    real(dp), intent(in) :: x, y, q

    sys%n_particles = sys%n_particles + 1
    sys%particles(sys%n_particles)%x = x
    sys%particles(sys%n_particles)%y = y
    sys%particles(sys%n_particles)%q = q
    sys%particles(sys%n_particles)%fx = 0.0_dp
    sys%particles(sys%n_particles)%fy = 0.0_dp
    sys%particles(sys%n_particles)%phi = 0.0_dp
  end subroutine fmm_add_particle

  !============================================================================
  ! Build quadtree
  !============================================================================
  subroutine fmm_build_tree(sys)
    type(fmm_system), intent(inout) :: sys
    integer, allocatable :: particle_list(:)
    integer :: i

    ! Reset tree
    sys%n_boxes = 1
    sys%root_id = 1
    sys%max_level = 0

    ! Initialize root box
    sys%boxes(1)%xc = 0.5_dp * (sys%xmin + sys%xmax)
    sys%boxes(1)%yc = 0.5_dp * (sys%ymin + sys%ymax)
    sys%boxes(1)%size = max(sys%xmax - sys%xmin, sys%ymax - sys%ymin)
    sys%boxes(1)%level = 0
    sys%boxes(1)%parent = 0
    sys%boxes(1)%children = 0

    ! All particles start in root
    allocate(particle_list(sys%n_particles))
    do i = 1, sys%n_particles
      particle_list(i) = i
    end do

    ! Recursively subdivide
    call subdivide_box(sys, 1, particle_list, sys%n_particles, 0)

    deallocate(particle_list)

    print *, 'Tree built:', sys%n_boxes, 'boxes,', sys%max_level, 'levels'

  end subroutine fmm_build_tree

  !============================================================================
  ! Recursively subdivide box
  !============================================================================
  recursive subroutine subdivide_box(sys, box_id, plist, np, level)
    type(fmm_system), intent(inout) :: sys
    integer, intent(in) :: box_id, np, level
    integer, intent(in) :: plist(np)

    integer :: i, ip, child_id, quad
    integer :: child_lists(4, sys%max_particles_per_box*2)
    integer :: child_counts(4)
    real(dp) :: xc, yc, dx, dy

    ! Check if should be leaf
    if (np <= sys%max_particles_per_box .or. level >= 15) then
      sys%boxes(box_id)%is_leaf = .true.
      sys%boxes(box_id)%n_particles = np
      allocate(sys%boxes(box_id)%particles(np))
      sys%boxes(box_id)%particles = plist(1:np)

      ! Allocate expansions
      allocate(sys%boxes(box_id)%M_real(0:sys%p_max, 0:sys%p_max))
      allocate(sys%boxes(box_id)%M_imag(0:sys%p_max, 1:sys%p_max))
      allocate(sys%boxes(box_id)%L_real(0:sys%p_max, 0:sys%p_max))
      allocate(sys%boxes(box_id)%L_imag(0:sys%p_max, 1:sys%p_max))

      sys%boxes(box_id)%M_real = 0.0_dp
      sys%boxes(box_id)%M_imag = 0.0_dp
      sys%boxes(box_id)%L_real = 0.0_dp
      sys%boxes(box_id)%L_imag = 0.0_dp

      ! Assign particles to box
      do i = 1, np
        sys%particles(plist(i))%box_id = box_id
      end do

      if (level > sys%max_level) sys%max_level = level
      return
    end if

    ! Subdivide into 4 quadrants
    sys%boxes(box_id)%is_leaf = .false.
    xc = sys%boxes(box_id)%xc
    yc = sys%boxes(box_id)%yc

    child_counts = 0
    do i = 1, np
      ip = plist(i)
      dx = sys%particles(ip)%x - xc
      dy = sys%particles(ip)%y - yc

      ! Determine quadrant: 1=SW, 2=SE, 3=NW, 4=NE
      quad = 1
      if (dx >= 0.0_dp) quad = quad + 1
      if (dy >= 0.0_dp) quad = quad + 2

      child_counts(quad) = child_counts(quad) + 1
      child_lists(quad, child_counts(quad)) = ip
    end do

    ! Create children
    do quad = 1, 4
      if (child_counts(quad) == 0) cycle

      sys%n_boxes = sys%n_boxes + 1
      child_id = sys%n_boxes
      sys%boxes(box_id)%children(quad) = child_id

      ! Set child properties
      sys%boxes(child_id)%level = level + 1
      sys%boxes(child_id)%parent = box_id
      sys%boxes(child_id)%size = 0.5_dp * sys%boxes(box_id)%size
      sys%boxes(child_id)%children = 0

      ! Child center
      dx = 0.25_dp * sys%boxes(box_id)%size
      dy = 0.25_dp * sys%boxes(box_id)%size
      sys%boxes(child_id)%xc = xc + merge(-dx, dx, mod(quad,2)==1)
      sys%boxes(child_id)%yc = yc + merge(-dy, dy, quad<=2)

      ! Recurse
      call subdivide_box(sys, child_id, &
                        child_lists(quad, 1:child_counts(quad)), &
                        child_counts(quad), level+1)
    end do

    ! Allocate expansions for non-leaf
    allocate(sys%boxes(box_id)%M_real(0:sys%p_max, 0:sys%p_max))
    allocate(sys%boxes(box_id)%M_imag(0:sys%p_max, 1:sys%p_max))
    allocate(sys%boxes(box_id)%L_real(0:sys%p_max, 0:sys%p_max))
    allocate(sys%boxes(box_id)%L_imag(0:sys%p_max, 1:sys%p_max))
    sys%boxes(box_id)%M_real = 0.0_dp
    sys%boxes(box_id)%M_imag = 0.0_dp
    sys%boxes(box_id)%L_real = 0.0_dp
    sys%boxes(box_id)%L_imag = 0.0_dp

  end subroutine subdivide_box

  !============================================================================
  ! Compute forces using FMM
  !============================================================================
  subroutine fmm_compute_forces(sys)
    type(fmm_system), intent(inout) :: sys

    integer :: i

    ! Reset forces
    do i = 1, sys%n_particles
      sys%particles(i)%fx = 0.0_dp
      sys%particles(i)%fy = 0.0_dp
      sys%particles(i)%phi = 0.0_dp
    end do

    print *, 'FMM computation:'
    print *, '  P2M: Computing multipole moments...'
    call p2m_all_leaves(sys)

    print *, '  M2M: Upward pass...'
    call m2m_upward_pass(sys)

    print *, '  M2L: Far field interactions...'
    call m2l_interaction_pass(sys)

    print *, '  L2L: Downward pass...'
    call l2l_downward_pass(sys)

    print *, '  L2P: Local expansion evaluation...'
    call l2p_all_leaves(sys)

    print *, '  P2P: Near field direct summation...'
    call p2p_near_field(sys)

    print *, 'FMM complete'

  end subroutine fmm_compute_forces

  !============================================================================
  ! P2M: Compute multipole expansion from particles (simplified for testing)
  !============================================================================
  subroutine p2m_all_leaves(sys)
    type(fmm_system), intent(inout) :: sys
    integer :: box_id

    do box_id = 1, sys%n_boxes
      if (sys%boxes(box_id)%is_leaf) then
        call p2m_box(sys, box_id)
      end if
    end do
  end subroutine p2m_all_leaves

  subroutine p2m_box(sys, box_id)
    type(fmm_system), intent(inout) :: sys
    integer, intent(in) :: box_id

    ! Simplified: just sum charges (monopole only for now)
    integer :: i, ip, n, m
    real(dp) :: dx, dy, r, phi, q, rn

    sys%boxes(box_id)%M_real = 0.0_dp
    sys%boxes(box_id)%M_imag = 0.0_dp

    do i = 1, sys%boxes(box_id)%n_particles
      ip = sys%boxes(box_id)%particles(i)
      q = sys%particles(ip)%q
      dx = sys%particles(ip)%x - sys%boxes(box_id)%xc
      dy = sys%particles(ip)%y - sys%boxes(box_id)%yc

      r = sqrt(dx**2 + dy**2)
      if (r < 1.0e-14_dp) then
        sys%boxes(box_id)%M_real(0,0) = sys%boxes(box_id)%M_real(0,0) + q
        cycle
      end if

      phi = atan2(dy, dx)

      ! M_n^m = sum q_i * r_i^n * exp(i*m*phi_i)
      do n = 0, sys%p_max
        rn = r**n
        do m = 0, n
          sys%boxes(box_id)%M_real(n,m) = sys%boxes(box_id)%M_real(n,m) + &
                                          q * rn * cos(m*phi)
          if (m > 0) then
            sys%boxes(box_id)%M_imag(n,m) = sys%boxes(box_id)%M_imag(n,m) + &
                                            q * rn * sin(m*phi)
          end if
        end do
      end do
    end do

  end subroutine p2m_box

  !============================================================================
  ! Stub implementations (to be completed)
  !============================================================================
  subroutine m2m_upward_pass(sys)
    type(fmm_system), intent(inout) :: sys
    ! TODO: Implement M2M translation
  end subroutine

  subroutine m2l_interaction_pass(sys)
    type(fmm_system), intent(inout) :: sys
    ! TODO: Implement M2L translation
  end subroutine

  subroutine l2l_downward_pass(sys)
    type(fmm_system), intent(inout) :: sys
    ! TODO: Implement L2L translation
  end subroutine

  subroutine l2p_all_leaves(sys)
    type(fmm_system), intent(inout) :: sys
    ! TODO: Implement L2P evaluation
  end subroutine

  !============================================================================
  ! P2P: Direct summation for near field
  !============================================================================
  subroutine p2p_near_field(sys)
    type(fmm_system), intent(inout) :: sys
    integer :: i

    ! For now, just do direct summation for all pairs
    do i = 1, sys%n_particles
      call p2p_particle(sys, i)
    end do
  end subroutine

  subroutine p2p_particle(sys, ip)
    type(fmm_system), intent(inout) :: sys
    integer, intent(in) :: ip

    integer :: jp
    real(dp) :: dx, dy, r, r3, qi, qj, fx, fy, phi

    qi = sys%particles(ip)%q

    do jp = 1, sys%n_particles
      if (ip == jp) cycle

      qj = sys%particles(jp)%q
      dx = sys%particles(jp)%x - sys%particles(ip)%x
      dy = sys%particles(jp)%y - sys%particles(ip)%y

      r = sqrt(dx**2 + dy**2)
      if (r < 1.0e-14_dp) cycle

      r3 = r**3

      ! Force: F = q_i * q_j * r_vec / r^3
      fx = qi * qj * dx / r3
      fy = qi * qj * dy / r3

      sys%particles(ip)%fx = sys%particles(ip)%fx + fx
      sys%particles(ip)%fy = sys%particles(ip)%fy + fy
      sys%particles(ip)%phi = sys%particles(ip)%phi + qj / r
    end do
  end subroutine p2p_particle

  !============================================================================
  ! Direct energy computation (reference)
  !============================================================================
  function fmm_compute_energy_direct(sys) result(energy)
    type(fmm_system), intent(in) :: sys
    real(dp) :: energy

    integer :: i, j
    real(dp) :: dx, dy, r, qi, qj

    energy = 0.0_dp

    do i = 1, sys%n_particles
      qi = sys%particles(i)%q
      do j = i+1, sys%n_particles
        qj = sys%particles(j)%q
        dx = sys%particles(j)%x - sys%particles(i)%x
        dy = sys%particles(j)%y - sys%particles(i)%y
        r = sqrt(dx**2 + dy**2)

        if (r > 1.0e-14_dp) then
          energy = energy + qi * qj / r
        end if
      end do
    end do

  end function fmm_compute_energy_direct

end module fmm_z0_precise
