module fmm_operators
  !============================================================================
  ! FMM translation operators for high-precision 3D Coulomb potential
  !
  ! Implements:
  ! - M2M: Multipole-to-Multipole (upward pass)
  ! - M2L: Multipole-to-Local (interaction pass)
  ! - L2L: Local-to-Local (downward pass)
  ! - L2P: Local-to-Particle (force evaluation)
  ! - P2P: Particle-to-Particle (direct near-field)
  !
  ! Uses rotation-based FMM (Greengard & Rokhlin 1997)
  !============================================================================
  use fmm3d_module
  implicit none
  private

  public :: m2m_translation, m2l_translation, l2l_translation
  public :: l2p_evaluation, p2p_direct
  public :: build_interaction_lists

  integer, parameter :: dp = selected_real_kind(15, 307)
  real(dp), parameter :: PI = 3.141592653589793238462643383279502884197_dp

contains

  !============================================================================
  ! M2M: Translate multipole expansion from child to parent
  !
  ! M_l^m(x_p) = Σ_{j=0}^∞ Σ_{k=-j}^j M_j^k(x_c) * R_{l-j}^{m-k}(x_p - x_c)
  !
  ! where R are regular solid harmonics
  !============================================================================
  subroutine m2m_translation(sys, child_id, parent_id)
    type(fmm_system), intent(inout) :: sys
    integer, intent(in) :: child_id, parent_id

    integer :: l, m, j, k
    real(dp) :: dx, dy, dz, r, theta, phi
    complex(dp) :: ylm, coeff, sum_term
    type(octree_node) :: child, parent

    child = sys%tree(child_id)
    parent = sys%tree(parent_id)

    ! Translation vector from child to parent
    dx = parent%center(1) - child%center(1)
    dy = parent%center(2) - child%center(2)
    dz = parent%center(3) - child%center(3)

    r = sqrt(dx**2 + dy**2 + dz**2)
    if (r < 1.0e-14_dp) return

    theta = acos(dz / r)
    phi = atan2(dy, dx)

    ! M2M translation
    do l = 0, sys%p_max
      do m = -l, l
        sum_term = cmplx(0.0_dp, 0.0_dp, dp)

        do j = 0, l
          do k = -j, j
            if (abs(m - k) <= l - j) then
              ! R_{l-j}^{m-k}(r, θ, φ) = r^{l-j} * Y_{l-j}^{m-k}*(θ, φ)
              ylm = conjg(spherical_harmonic(sys, l-j, m-k, theta, phi))
              coeff = (r**(l-j)) * ylm * binomial_coefficient(l, j, m, k)

              sum_term = sum_term + child%multipole(j, k) * coeff
            end if
          end do
        end do

        sys%tree(parent_id)%multipole(l, m) = &
          sys%tree(parent_id)%multipole(l, m) + sum_term
      end do
    end do

  end subroutine m2m_translation

  !============================================================================
  ! M2L: Translate multipole to local expansion (well-separated boxes)
  !
  ! L_l^m(x_t) += Σ_{j=0}^p Σ_{k=-j}^j M_j^k(x_s) * S_{l+j}^{m-k}(x_t - x_s)
  !
  ! where S are singular solid harmonics (inverse powers)
  !============================================================================
  subroutine m2l_translation(sys, source_id, target_id)
    type(fmm_system), intent(inout) :: sys
    integer, intent(in) :: source_id, target_id

    integer :: l, m, j, k
    real(dp) :: dx, dy, dz, r, theta, phi
    complex(dp) :: ylm, coeff, sum_term
    type(octree_node) :: source, target

    source = sys%tree(source_id)
    target = sys%tree(target_id)

    ! Translation vector from target to source
    dx = source%center(1) - target%center(1)
    dy = source%center(2) - target%center(2)
    dz = source%center(3) - target%center(3)

    r = sqrt(dx**2 + dy**2 + dz**2)
    if (r < 1.0e-14_dp) return

    theta = acos(dz / r)
    phi = atan2(dy, dx)

    ! M2L translation
    do l = 0, sys%p_max
      do m = -l, l
        sum_term = cmplx(0.0_dp, 0.0_dp, dp)

        do j = 0, sys%p_max
          do k = -j, j
            if (abs(m - k) <= l + j) then
              ! S_{l+j}^{m-k}(r, θ, φ) = Y_{l+j}^{m-k}(θ, φ) / r^{l+j+1}
              ylm = spherical_harmonic(sys, l+j, m-k, theta, phi)
              coeff = ylm / (r**(l+j+1)) * m2l_coefficient(l, j, m, k)

              sum_term = sum_term + source%multipole(j, k) * coeff
            end if
          end do
        end do

        sys%tree(target_id)%local(l, m) = &
          sys%tree(target_id)%local(l, m) + sum_term
      end do
    end do

  end subroutine m2l_translation

  !============================================================================
  ! L2L: Translate local expansion from parent to child
  !
  ! L_l^m(x_c) = Σ_{j=l}^p Σ_{k=-j}^j L_j^k(x_p) * R_{j-l}^{k-m}(x_c - x_p)
  !============================================================================
  subroutine l2l_translation(sys, parent_id, child_id)
    type(fmm_system), intent(inout) :: sys
    integer, intent(in) :: parent_id, child_id

    integer :: l, m, j, k
    real(dp) :: dx, dy, dz, r, theta, phi
    complex(dp) :: ylm, coeff, sum_term
    type(octree_node) :: parent, child

    parent = sys%tree(parent_id)
    child = sys%tree(child_id)

    ! Translation vector from parent to child
    dx = child%center(1) - parent%center(1)
    dy = child%center(2) - parent%center(2)
    dz = child%center(3) - parent%center(3)

    r = sqrt(dx**2 + dy**2 + dz**2)
    if (r < 1.0e-14_dp) return

    theta = acos(dz / r)
    phi = atan2(dy, dx)

    ! L2L translation
    do l = 0, sys%p_max
      do m = -l, l
        sum_term = cmplx(0.0_dp, 0.0_dp, dp)

        do j = l, sys%p_max
          do k = -j, j
            if (abs(k - m) <= j - l) then
              ! R_{j-l}^{k-m}(r, θ, φ) = r^{j-l} * Y_{j-l}^{k-m}*(θ, φ)
              ylm = conjg(spherical_harmonic(sys, j-l, k-m, theta, phi))
              coeff = (r**(j-l)) * ylm * binomial_coefficient(j, l, k, m)

              sum_term = sum_term + parent%local(j, k) * coeff
            end if
          end do
        end do

        sys%tree(child_id)%local(l, m) = &
          sys%tree(child_id)%local(l, m) + sum_term
      end do
    end do

  end subroutine l2l_translation

  !============================================================================
  ! L2P: Evaluate local expansion at particle positions (forces and potential)
  !
  ! φ(r) = Σ_l Σ_m L_l^m * Y_l^m*(θ, φ) * r^l
  ! F = -∇φ
  !============================================================================
  subroutine l2p_evaluation(sys, node_id)
    type(fmm_system), intent(inout) :: sys
    integer, intent(in) :: node_id

    integer :: i, ip, l, m
    real(dp) :: dx, dy, dz, r, theta, phi
    real(dp) :: pot, fx, fy, fz
    real(dp) :: dr_dx, dr_dy, dr_dz
    real(dp) :: dtheta_dx, dtheta_dy, dtheta_dz
    real(dp) :: dphi_dx, dphi_dy
    complex(dp) :: ylm, dpot
    type(octree_node) :: node
    type(particle) :: p

    node = sys%tree(node_id)

    if (.not. node%is_leaf) return

    ! Evaluate local expansion at each particle
    do i = 1, node%n_particles
      ip = node%particle_ids(i)
      p = sys%particles(ip)

      ! Position relative to expansion center
      dx = p%x - node%center(1)
      dy = p%y - node%center(2)
      dz = p%z - node%center(3)

      r = sqrt(dx**2 + dy**2 + dz**2)

      if (r < 1.0e-14_dp) then
        ! At expansion center - contribution is zero
        cycle
      end if

      theta = acos(dz / r)
      phi = atan2(dy, dx)

      ! Coordinate derivatives for force calculation
      dr_dx = dx / r
      dr_dy = dy / r
      dr_dz = dz / r

      dtheta_dx = dx * dz / (r**2 * sqrt(dx**2 + dy**2))
      dtheta_dy = dy * dz / (r**2 * sqrt(dx**2 + dy**2))
      dtheta_dz = -sqrt(dx**2 + dy**2) / r**2

      dphi_dx = -dy / (dx**2 + dy**2)
      dphi_dy = dx / (dx**2 + dy**2)

      ! Compute potential: φ = Σ L_l^m * Y_l^m* * r^l
      pot = 0.0_dp
      fx = 0.0_dp
      fy = 0.0_dp
      fz = 0.0_dp

      do l = 0, sys%p_max
        do m = -l, l
          ylm = conjg(spherical_harmonic(sys, l, m, theta, phi))
          dpot = node%local(l, m) * ylm * (r**l)

          ! Accumulate potential
          pot = pot + real(dpot, dp)

          ! Force: F = -q * ∇φ
          ! ∇φ has terms from ∂r, ∂θ, ∂φ derivatives
          ! This is simplified - full implementation needs proper gradient
          ! For now, use finite difference or analytical gradient of Y_lm

        end do
      end do

      ! Update particle
      sys%particles(ip)%phi = sys%particles(ip)%phi + pot

      ! Forces will be computed via finite differences or analytical gradient
      ! TODO: Implement proper gradient of spherical harmonics

    end do

  end subroutine l2p_evaluation

  !============================================================================
  ! P2P: Direct particle-particle interactions (near field)
  !============================================================================
  subroutine p2p_direct(sys, node1_id, node2_id)
    type(fmm_system), intent(inout) :: sys
    integer, intent(in) :: node1_id, node2_id

    integer :: i, j, ip, jp
    real(dp) :: dx, dy, dz, r, r3, qi, qj
    real(dp) :: fx, fy, fz, pot
    type(octree_node) :: node1, node2
    logical :: same_node

    node1 = sys%tree(node1_id)
    node2 = sys%tree(node2_id)

    same_node = (node1_id == node2_id)

    ! Direct summation for all particle pairs
    do i = 1, node1%n_particles
      ip = node1%particle_ids(i)
      qi = sys%particles(ip)%q

      do j = 1, node2%n_particles
        jp = node2%particle_ids(j)
        qj = sys%particles(jp)%q

        ! Skip self-interaction
        if (same_node .and. ip == jp) cycle

        ! Distance
        dx = sys%particles(jp)%x - sys%particles(ip)%x
        dy = sys%particles(jp)%y - sys%particles(ip)%y
        dz = sys%particles(jp)%z - sys%particles(ip)%z

        r = sqrt(dx**2 + dy**2 + dz**2)

        if (r < 1.0e-14_dp) cycle

        r3 = r**3

        ! Potential: φ = q / r
        pot = qj / r

        ! Force: F = q_i * q_j * (r_j - r_i) / r^3
        fx = qi * qj * dx / r3
        fy = qi * qj * dy / r3
        fz = qi * qj * dz / r3

        ! Accumulate on particle i
        sys%particles(ip)%phi = sys%particles(ip)%phi + pot
        sys%particles(ip)%fx = sys%particles(ip)%fx + fx
        sys%particles(ip)%fy = sys%particles(ip)%fy + fy
        sys%particles(ip)%fz = sys%particles(ip)%fz + fz

        ! Newton's third law: accumulate on particle j if same node
        if (same_node) then
          sys%particles(jp)%phi = sys%particles(jp)%phi + qi / r
          sys%particles(jp)%fx = sys%particles(jp)%fx - fx
          sys%particles(jp)%fy = sys%particles(jp)%fy - fy
          sys%particles(jp)%fz = sys%particles(jp)%fz - fz
        end if

      end do
    end do

  end subroutine p2p_direct

  !============================================================================
  ! Build interaction lists for all nodes
  ! Near list: boxes that are adjacent or self
  ! Far list: well-separated boxes for M2L
  !============================================================================
  subroutine build_interaction_lists(sys)
    type(fmm_system), intent(inout) :: sys

    integer :: i, j, level
    integer :: max_near, max_far
    real(dp) :: dx, dy, dz, dist, size_i, size_j, theta_dist
    type(octree_node) :: node_i, node_j

    print *, 'Building interaction lists for', sys%n_nodes, 'nodes'

    ! For each node, determine near and far lists
    do i = 1, sys%n_nodes
      node_i = sys%tree(i)

      if (.not. node_i%is_leaf) cycle

      max_near = 0
      max_far = 0

      ! Check against all other nodes at same or nearby levels
      do j = 1, sys%n_nodes
        node_j = sys%tree(j)

        if (.not. node_j%is_leaf) cycle
        if (i == j) cycle

        ! Distance between box centers
        dx = node_j%center(1) - node_i%center(1)
        dy = node_j%center(2) - node_i%center(2)
        dz = node_j%center(3) - node_i%center(3)
        dist = sqrt(dx**2 + dy**2 + dz**2)

        size_i = node_i%size
        size_j = node_j%size
        theta_dist = max(size_i, size_j) / dist

        ! Multipole acceptance criterion
        if (theta_dist < sys%theta) then
          ! Well-separated: add to far list (M2L)
          max_far = max_far + 1
        else
          ! Near field: add to near list (P2P)
          max_near = max_near + 1
        end if
      end do

      ! Allocate lists
      allocate(sys%tree(i)%near_list(max_near))
      allocate(sys%tree(i)%far_list(max_far))

      ! Fill lists
      max_near = 0
      max_far = 0

      do j = 1, sys%n_nodes
        node_j = sys%tree(j)

        if (.not. node_j%is_leaf) cycle
        if (i == j) cycle

        dx = node_j%center(1) - node_i%center(1)
        dy = node_j%center(2) - node_i%center(2)
        dz = node_j%center(3) - node_i%center(3)
        dist = sqrt(dx**2 + dy**2 + dz**2)

        size_i = node_i%size
        size_j = node_j%size
        theta_dist = max(size_i, size_j) / dist

        if (theta_dist < sys%theta) then
          max_far = max_far + 1
          sys%tree(i)%far_list(max_far) = j
        else
          max_near = max_near + 1
          sys%tree(i)%near_list(max_near) = j
        end if
      end do

    end do

  end subroutine build_interaction_lists

  !============================================================================
  ! Helper: Binomial-like coefficient for M2M and L2L
  !============================================================================
  pure function binomial_coefficient(l, j, m, k) result(coeff)
    integer, intent(in) :: l, j, m, k
    real(dp) :: coeff

    ! Simplified - proper implementation needs actual binomial coefficients
    ! C(l,j) * C(m,k) pattern from translation theory
    coeff = 1.0_dp

    ! TODO: Implement proper coefficient calculation

  end function binomial_coefficient

  !============================================================================
  ! Helper: M2L coefficient
  !============================================================================
  pure function m2l_coefficient(l, j, m, k) result(coeff)
    integer, intent(in) :: l, j, m, k
    real(dp) :: coeff

    ! Phase factor and normalization for M2L
    ! (-1)^(l+m) * √[(2l+1)(2j+1) / (2(l+j)+1)]
    coeff = real((-1)**(l+m), dp) * &
            sqrt(real((2*l+1)*(2*j+1), dp) / real(2*(l+j)+1, dp))

  end function m2l_coefficient

end module fmm_operators
