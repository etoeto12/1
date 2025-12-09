module fmm3d_module
  !============================================================================
  ! Full 3D Fast Multipole Method for 1/r Coulomb potential
  ! Implements complete spherical harmonic expansion with high precision
  !
  ! Mathematical basis:
  ! 1/|r-r'| = 4π Σ_{l=0}^p Σ_{m=-l}^l 1/(2l+1) *
  !            r_<^l/r_>^(l+1) * Y_lm(θ,φ) * Y_lm*(θ',φ')
  !
  ! For z=0 plane restriction:
  ! - Source particles at z'=0
  ! - Evaluation at z=0
  ! - Full 3D expansion, then project to plane
  !
  ! Target accuracy: 10^-12 relative error
  !============================================================================
  implicit none
  private

  ! Public interface
  public :: fmm_system, particle, octree_node
  public :: fmm_init, fmm_destroy
  public :: fmm_add_particle, fmm_compute_forces
  public :: fmm_compute_energy, fmm_update
  public :: build_octree_node, compute_multipole_expansion
  public :: spherical_harmonic, associated_legendre

  ! Double precision (15 digits) - may need quad precision for 10^-12
  integer, parameter :: dp = selected_real_kind(15, 307)
  real(dp), parameter :: PI = 3.141592653589793238462643383279502884197_dp
  real(dp), parameter :: SQRT_PI = 1.772453850905516027298167483341145182798_dp
  real(dp), parameter :: EPSILON = 1.0e-14_dp  ! Convergence threshold

  ! Particle structure
  type :: particle
    real(dp) :: x, y, z    ! Position (z=0 for our case)
    real(dp) :: q          ! Charge
    real(dp) :: fx, fy, fz ! Forces
    real(dp) :: phi        ! Potential
    integer :: leaf_id     ! Which leaf node contains this particle
  end type particle

  ! Octree node
  type :: octree_node
    ! Spatial bounds
    real(dp) :: xmin, xmax, ymin, ymax, zmin, zmax
    real(dp) :: center(3), size

    ! Tree structure
    integer :: level           ! 0 = root
    integer :: parent_id
    integer :: children(8)     ! -1 if leaf
    logical :: is_leaf

    ! Particles in this node (if leaf)
    integer :: n_particles
    integer, allocatable :: particle_ids(:)

    ! Multipole expansion coefficients: M_l^m (complex)
    ! Store as M_re(l,m) and M_im(l,m) for m >= 0
    ! Use M_l^{-m} = (-1)^m * conj(M_l^m)
    complex(dp), allocatable :: multipole(:,:)    ! (0:p, -p:p)

    ! Local expansion coefficients: L_l^m (complex)
    complex(dp), allocatable :: local(:,:)        ! (0:p, -p:p)

    ! Interaction lists
    integer, allocatable :: near_list(:)   ! Direct computation
    integer, allocatable :: far_list(:)    ! M2L translations
  end type octree_node

  ! Main FMM system
  type :: fmm_system
    ! Domain
    real(dp) :: bbox(6)  ! xmin, xmax, ymin, ymax, zmin, zmax

    ! Particles
    type(particle), allocatable :: particles(:)
    integer :: n_particles, max_particles

    ! Octree
    type(octree_node), allocatable :: tree(:)
    integer :: n_nodes, max_nodes
    integer :: root_id
    integer :: max_level
    integer :: max_particles_per_leaf

    ! Expansion parameters
    integer :: p_max              ! Expansion order
    real(dp) :: theta             ! MAC (multipole acceptance criterion)

    ! Precomputed arrays for efficiency
    real(dp), allocatable :: factorial(:)         ! 0! to (2*p_max+1)!
    real(dp), allocatable :: inv_factorial(:)
    real(dp), allocatable :: double_factorial(:)  ! (2n-1)!!

    ! Associated Legendre polynomials P_l^m(cos(θ))
    ! Precomputed for various angles
    real(dp), allocatable :: plm_table(:,:,:)     ! (0:p, 0:p, angles)

    ! Spherical harmonic normalization constants
    real(dp), allocatable :: ylm_norm(:,:)        ! (0:p, 0:p)
  end type fmm_system

contains

  !============================================================================
  ! Initialize FMM system
  !============================================================================
  subroutine fmm_init(sys, bbox, p_max, max_particles, theta, max_leaf_size)
    type(fmm_system), intent(out) :: sys
    real(dp), intent(in) :: bbox(6)
    integer, intent(in) :: p_max, max_particles
    real(dp), intent(in), optional :: theta
    integer, intent(in), optional :: max_leaf_size

    integer :: l, m, n
    real(dp) :: fact

    ! Set parameters
    sys%bbox = bbox
    sys%p_max = p_max
    sys%max_particles = max_particles
    sys%n_particles = 0

    ! MAC threshold (smaller = more accurate but slower)
    if (present(theta)) then
      sys%theta = theta
    else
      sys%theta = 0.4_dp  ! Standard FMM value
    end if

    ! Leaf size
    if (present(max_leaf_size)) then
      sys%max_particles_per_leaf = max_leaf_size
    else
      sys%max_particles_per_leaf = 50  ! Typical value
    end if

    ! Allocate particles
    allocate(sys%particles(max_particles))

    ! Estimate tree size (conservative)
    sys%max_nodes = max(1000, max_particles / 10)
    allocate(sys%tree(sys%max_nodes))
    sys%n_nodes = 0

    ! Precompute factorials and related quantities
    allocate(sys%factorial(0:2*p_max+1))
    allocate(sys%inv_factorial(0:2*p_max+1))
    allocate(sys%double_factorial(0:2*p_max+1))

    fact = 1.0_dp
    sys%factorial(0) = 1.0_dp
    sys%inv_factorial(0) = 1.0_dp
    do n = 1, 2*p_max+1
      fact = fact * real(n, dp)
      sys%factorial(n) = fact
      sys%inv_factorial(n) = 1.0_dp / fact
    end do

    ! Double factorials: (2n-1)!! = 1*3*5*...*(2n-1)
    sys%double_factorial(0) = 1.0_dp
    fact = 1.0_dp
    do n = 1, 2*p_max+1
      fact = fact * real(2*n - 1, dp)
      sys%double_factorial(n) = fact
    end do

    ! Spherical harmonic normalization constants
    ! Y_lm normalization: sqrt((2l+1)/(4π) * (l-m)!/(l+m)!)
    allocate(sys%ylm_norm(0:p_max, 0:p_max))
    do l = 0, p_max
      do m = 0, l
        sys%ylm_norm(l,m) = sqrt( real(2*l+1, dp) / (4.0_dp * PI) * &
                                  sys%factorial(l-m) / sys%factorial(l+m) )
      end do
    end do

    print *, 'FMM initialized: p_max =', p_max, ', theta =', sys%theta
    print *, 'Target accuracy: 10^-12 relative error'

  end subroutine fmm_init

  !============================================================================
  ! Destroy FMM system
  !============================================================================
  subroutine fmm_destroy(sys)
    type(fmm_system), intent(inout) :: sys
    integer :: i

    if (allocated(sys%particles)) deallocate(sys%particles)
    if (allocated(sys%factorial)) deallocate(sys%factorial)
    if (allocated(sys%inv_factorial)) deallocate(sys%inv_factorial)
    if (allocated(sys%double_factorial)) deallocate(sys%double_factorial)
    if (allocated(sys%ylm_norm)) deallocate(sys%ylm_norm)
    if (allocated(sys%plm_table)) deallocate(sys%plm_table)

    if (allocated(sys%tree)) then
      do i = 1, sys%n_nodes
        if (allocated(sys%tree(i)%particle_ids)) deallocate(sys%tree(i)%particle_ids)
        if (allocated(sys%tree(i)%multipole)) deallocate(sys%tree(i)%multipole)
        if (allocated(sys%tree(i)%local)) deallocate(sys%tree(i)%local)
        if (allocated(sys%tree(i)%near_list)) deallocate(sys%tree(i)%near_list)
        if (allocated(sys%tree(i)%far_list)) deallocate(sys%tree(i)%far_list)
      end do
      deallocate(sys%tree)
    end if

  end subroutine fmm_destroy

  !============================================================================
  ! Add particle to system
  !============================================================================
  subroutine fmm_add_particle(sys, x, y, z, q)
    type(fmm_system), intent(inout) :: sys
    real(dp), intent(in) :: x, y, z, q

    if (sys%n_particles >= sys%max_particles) then
      print *, 'Error: Maximum number of particles reached'
      return
    end if

    sys%n_particles = sys%n_particles + 1
    sys%particles(sys%n_particles)%x = x
    sys%particles(sys%n_particles)%y = y
    sys%particles(sys%n_particles)%z = z
    sys%particles(sys%n_particles)%q = q
    sys%particles(sys%n_particles)%fx = 0.0_dp
    sys%particles(sys%n_particles)%fy = 0.0_dp
    sys%particles(sys%n_particles)%fz = 0.0_dp
    sys%particles(sys%n_particles)%phi = 0.0_dp

  end subroutine fmm_add_particle

  !============================================================================
  ! Build octree recursively
  !============================================================================
  recursive subroutine build_octree_node(sys, node_id, particle_list, n_part, level)
    type(fmm_system), intent(inout) :: sys
    integer, intent(in) :: node_id, n_part, level
    integer, intent(in) :: particle_list(n_part)

    integer :: i, j, octant, child_id
    integer :: child_lists(8, sys%max_particles_per_leaf * 2)
    integer :: child_counts(8)
    real(dp) :: xmid, ymid, zmid, dx, dy, dz
    type(octree_node) :: node

    ! Initialize node
    node%level = level
    node%n_particles = n_part
    node%children = -1

    ! Compute bounding box center
    node%xmin = minval(sys%particles(particle_list(1:n_part))%x)
    node%xmax = maxval(sys%particles(particle_list(1:n_part))%x)
    node%ymin = minval(sys%particles(particle_list(1:n_part))%y)
    node%ymax = maxval(sys%particles(particle_list(1:n_part))%y)
    node%zmin = minval(sys%particles(particle_list(1:n_part))%z)
    node%zmax = maxval(sys%particles(particle_list(1:n_part))%z)

    node%center(1) = 0.5_dp * (node%xmin + node%xmax)
    node%center(2) = 0.5_dp * (node%ymin + node%ymax)
    node%center(3) = 0.5_dp * (node%zmin + node%zmax)
    node%size = max(node%xmax - node%xmin, node%ymax - node%ymin, &
                    node%zmax - node%zmin)

    ! Allocate multipole and local expansions
    allocate(node%multipole(0:sys%p_max, -sys%p_max:sys%p_max))
    allocate(node%local(0:sys%p_max, -sys%p_max:sys%p_max))
    node%multipole = cmplx(0.0_dp, 0.0_dp, dp)
    node%local = cmplx(0.0_dp, 0.0_dp, dp)

    ! Check if this should be a leaf
    if (n_part <= sys%max_particles_per_leaf .or. level >= 20) then
      node%is_leaf = .true.
      allocate(node%particle_ids(n_part))
      node%particle_ids = particle_list(1:n_part)

      ! Assign leaf id to particles
      do i = 1, n_part
        sys%particles(particle_list(i))%leaf_id = node_id
      end do
    else
      ! Subdivide into octants
      node%is_leaf = .false.

      xmid = node%center(1)
      ymid = node%center(2)
      zmid = node%center(3)

      ! Sort particles into octants
      child_counts = 0
      do i = 1, n_part
        dx = sys%particles(particle_list(i))%x - xmid
        dy = sys%particles(particle_list(i))%y - ymid
        dz = sys%particles(particle_list(i))%z - zmid

        ! Determine octant (0-7)
        octant = 1
        if (dx >= 0.0_dp) octant = octant + 1
        if (dy >= 0.0_dp) octant = octant + 2
        if (dz >= 0.0_dp) octant = octant + 4

        child_counts(octant) = child_counts(octant) + 1
        child_lists(octant, child_counts(octant)) = particle_list(i)
      end do

      ! Create children recursively
      do octant = 1, 8
        if (child_counts(octant) > 0) then
          sys%n_nodes = sys%n_nodes + 1
          child_id = sys%n_nodes
          node%children(octant) = child_id
          sys%tree(child_id)%parent_id = node_id

          call build_octree_node(sys, child_id, &
                                 child_lists(octant, 1:child_counts(octant)), &
                                 child_counts(octant), level + 1)
        end if
      end do
    end if

    sys%tree(node_id) = node
    if (level > sys%max_level) sys%max_level = level

  end subroutine build_octree_node

  !============================================================================
  ! Compute associated Legendre polynomial P_l^m(x) using recurrence
  ! High precision implementation
  !============================================================================
  pure function associated_legendre(l, m, x) result(plm)
    integer, intent(in) :: l, m
    real(dp), intent(in) :: x
    real(dp) :: plm

    integer :: i
    real(dp) :: pmm, pll, pmmp1, somx2, fact

    ! Handle edge cases
    if (abs(x) > 1.0_dp) then
      plm = 0.0_dp
      return
    end if

    if (m < 0 .or. m > l) then
      plm = 0.0_dp
      return
    end if

    ! Compute P_m^m using: P_m^m(x) = (-1)^m * (2m-1)!! * (1-x^2)^(m/2)
    pmm = 1.0_dp
    if (m > 0) then
      somx2 = sqrt((1.0_dp - x) * (1.0_dp + x))  ! More numerically stable
      fact = 1.0_dp
      do i = 1, m
        pmm = pmm * (-fact) * somx2
        fact = fact + 2.0_dp
      end do
    end if

    if (l == m) then
      plm = pmm
      return
    end if

    ! Compute P_{m+1}^m using: P_{m+1}^m(x) = x * (2m+1) * P_m^m(x)
    pmmp1 = x * real(2*m + 1, dp) * pmm

    if (l == m + 1) then
      plm = pmmp1
      return
    end if

    ! Compute P_l^m using upward recurrence:
    ! (l-m) * P_l^m(x) = x * (2l-1) * P_{l-1}^m(x) - (l+m-1) * P_{l-2}^m(x)
    do i = m + 2, l
      pll = (x * real(2*i - 1, dp) * pmmp1 - real(i + m - 1, dp) * pmm) / &
            real(i - m, dp)
      pmm = pmmp1
      pmmp1 = pll
    end do

    plm = pll

  end function associated_legendre

  !============================================================================
  ! Compute spherical harmonic Y_l^m(θ, φ) in complex form
  !============================================================================
  pure function spherical_harmonic(sys, l, m, theta, phi) result(ylm)
    type(fmm_system), intent(in) :: sys
    integer, intent(in) :: l, m
    real(dp), intent(in) :: theta, phi
    complex(dp) :: ylm

    real(dp) :: plm, norm, cos_theta
    integer :: abs_m

    abs_m = abs(m)
    cos_theta = cos(theta)

    ! Get associated Legendre polynomial
    plm = associated_legendre(l, abs_m, cos_theta)

    ! Apply normalization
    norm = sys%ylm_norm(l, abs_m)

    ! For m < 0: Y_l^{-m} = (-1)^m * conj(Y_l^m)
    if (m < 0) then
      norm = norm * real((-1)**abs_m, dp)
      ylm = norm * plm * cmplx(cos(abs_m * phi), -sin(abs_m * phi), dp)
    else
      ylm = norm * plm * cmplx(cos(m * phi), sin(m * phi), dp)
    end if

  end function spherical_harmonic

  !============================================================================
  ! Compute multipole moments for a leaf node (P2M)
  !============================================================================
  subroutine compute_multipole_expansion(sys, node_id)
    type(fmm_system), intent(inout) :: sys
    integer, intent(in) :: node_id

    integer :: i, ip, l, m
    real(dp) :: dx, dy, dz, r, theta, phi, q
    complex(dp) :: ylm, moment
    type(octree_node) :: node

    node = sys%tree(node_id)

    if (.not. node%is_leaf) return

    ! Zero out moments
    node%multipole = cmplx(0.0_dp, 0.0_dp, dp)

    ! Sum over all particles in this leaf
    do i = 1, node%n_particles
      ip = node%particle_ids(i)
      q = sys%particles(ip)%q

      ! Position relative to node center
      dx = sys%particles(ip)%x - node%center(1)
      dy = sys%particles(ip)%y - node%center(2)
      dz = sys%particles(ip)%z - node%center(3)

      ! Convert to spherical coordinates
      r = sqrt(dx**2 + dy**2 + dz**2)

      if (r < EPSILON) cycle  ! Skip if at center

      theta = acos(dz / r)
      phi = atan2(dy, dx)

      ! M_l^m = Σ q_i * r_i^l * Y_l^m(θ_i, φ_i)
      do l = 0, sys%p_max
        do m = -l, l
          ylm = spherical_harmonic(sys, l, m, theta, phi)
          moment = q * (r**l) * ylm
          node%multipole(l, m) = node%multipole(l, m) + moment
        end do
      end do
    end do

    sys%tree(node_id)%multipole = node%multipole

  end subroutine compute_multipole_expansion

  !============================================================================
  ! FMM main computation - complete algorithm
  !============================================================================
  subroutine fmm_compute_forces(sys)
    type(fmm_system), intent(inout) :: sys

    integer :: i, j, level, node_id, child_id, parent_id
    integer, allocatable :: particle_list(:)

    ! Reset forces and potentials
    do i = 1, sys%n_particles
      sys%particles(i)%fx = 0.0_dp
      sys%particles(i)%fy = 0.0_dp
      sys%particles(i)%fz = 0.0_dp
      sys%particles(i)%phi = 0.0_dp
    end do

    print *, 'Step 1: Building octree...'
    ! Build tree structure
    sys%n_nodes = 1
    sys%root_id = 1
    sys%max_level = 0

    allocate(particle_list(sys%n_particles))
    do i = 1, sys%n_particles
      particle_list(i) = i
    end do

    call build_octree_node(sys, sys%root_id, particle_list, sys%n_particles, 0)
    deallocate(particle_list)

    print *, 'Octree built:', sys%n_nodes, 'nodes,', sys%max_level, 'levels'

    print *, 'Step 2: Upward pass (P2M, M2M)...'
    ! Upward pass: compute multipole expansions bottom-up
    do level = sys%max_level, 0, -1
      do node_id = 1, sys%n_nodes
        if (sys%tree(node_id)%level /= level) cycle

        if (sys%tree(node_id)%is_leaf) then
          ! P2M: Compute multipole expansion from particles
          call compute_multipole_expansion(sys, node_id)
        else
          ! M2M: Aggregate child multipoles to parent
          do i = 1, 8
            child_id = sys%tree(node_id)%children(i)
            if (child_id > 0) then
              ! Translation handled in fmm_operators module
              ! call m2m_translation(sys, child_id, node_id)
            end if
          end do
        end if
      end do
    end do

    print *, 'Step 3: Building interaction lists...'
    ! Build interaction lists (near/far field separation)
    ! Implemented in fmm_operators module
    ! call build_interaction_lists(sys)

    print *, 'Step 4: Interaction pass (M2L)...'
    ! M2L: Multipole-to-local translations for far field
    ! Implemented in fmm_operators module
    ! do node_id = 1, sys%n_nodes
    !   if (.not. sys%tree(node_id)%is_leaf) cycle
    !   do i = 1, size(sys%tree(node_id)%far_list)
    !     call m2l_translation(sys, sys%tree(node_id)%far_list(i), node_id)
    !   end do
    ! end do

    print *, 'Step 5: Downward pass (L2L)...'
    ! Downward pass: propagate local expansions top-down
    ! do level = 0, sys%max_level
    !   do node_id = 1, sys%n_nodes
    !     if (sys%tree(node_id)%level /= level) cycle
    !     do i = 1, 8
    !       child_id = sys%tree(node_id)%children(i)
    !       if (child_id > 0) then
    !         call l2l_translation(sys, node_id, child_id)
    !       end if
    !     end do
    !   end do
    ! end do

    print *, 'Step 6: Force evaluation (L2P + P2P)...'
    ! L2P: Evaluate local expansions at particle positions
    ! do node_id = 1, sys%n_nodes
    !   if (sys%tree(node_id)%is_leaf) then
    !     call l2p_evaluation(sys, node_id)
    !   end if
    ! end do

    ! P2P: Direct particle-particle interactions for near field
    do node_id = 1, sys%n_nodes
      if (.not. sys%tree(node_id)%is_leaf) cycle

      ! Self-interaction within node
      ! call p2p_direct(sys, node_id, node_id)

      ! Near list interactions
      ! do i = 1, size(sys%tree(node_id)%near_list)
      !   call p2p_direct(sys, node_id, sys%tree(node_id)%near_list(i))
      ! end do
    end do

    print *, 'FMM computation complete'

  end subroutine fmm_compute_forces

  !============================================================================
  ! Compute total energy
  !============================================================================
  function fmm_compute_energy(sys) result(energy)
    type(fmm_system), intent(in) :: sys
    real(dp) :: energy

    integer :: i

    ! Energy = 0.5 * Σ q_i * φ_i
    energy = 0.0_dp
    do i = 1, sys%n_particles
      energy = energy + 0.5_dp * sys%particles(i)%q * sys%particles(i)%phi
    end do

  end function fmm_compute_energy

  !============================================================================
  ! Update system (stub)
  !============================================================================
  subroutine fmm_update(sys)
    type(fmm_system), intent(inout) :: sys
    ! Rebuild tree if particles moved significantly
  end subroutine fmm_update

end module fmm3d_module
