module multipole_module
  !============================================================================
  ! Multipole expansion for 1/r potential in z=0 plane
  ! Using real spherical harmonics restricted to equatorial plane
  !
  ! For particles in z=0 plane:
  ! Phi(r) = sum_i q_i / |r - r_i|
  !
  ! Multipole expansion:
  ! 1/|r-r'| = sum_{n=0}^{p} sum_{m=-n}^{n} M_n^m(r') * R_n^m(r)
  !
  ! where M_n^m are multipole moments and R_n^m are regular solid harmonics
  !============================================================================
  implicit none
  private

  ! Public interface
  public :: multipole_system
  public :: init_system, destroy_system
  public :: add_particle, compute_energy, compute_force
  public :: compute_all_forces, move_particles, update_cells
  public :: write_frame

  ! Precision
  integer, parameter :: dp = selected_real_kind(15, 307)
  real(dp), parameter :: PI = 3.141592653589793238462643383279502884197_dp

  ! Particle type
  type :: particle
    real(dp) :: q              ! Charge
    real(dp) :: x, y           ! Position (z=0 implied)
    real(dp) :: fx, fy         ! Force components
    integer :: cell_i, cell_j  ! Cell indices
  end type particle

  ! Cell type
  type :: cell
    integer :: n_particles
    integer, allocatable :: particle_ids(:)
    integer :: capacity

    ! Multipole moments: M_n^m = sum_i q_i * r_i^n * Y_n^m(theta_i, phi_i)
    ! For z=0: only m=0 terms survive (axisymmetric in simplified case)
    ! Full 2D: we need cos(m*phi) and sin(m*phi) terms
    real(dp), allocatable :: M_cos(:,:)  ! (0:p_max, 0:p_max) - cosine moments
    real(dp), allocatable :: M_sin(:,:)  ! (0:p_max, 1:p_max) - sine moments

    real(dp) :: center_x, center_y
  end type cell

  ! Main system type
  type :: multipole_system
    ! Domain
    real(dp) :: xmax, ymax
    integer :: nx_cells, ny_cells
    real(dp) :: dx_cell, dy_cell

    ! Grid
    type(cell), allocatable :: cells(:,:)

    ! Particles
    type(particle), allocatable :: particles(:)
    integer :: n_particles
    integer :: max_particles

    ! Expansion order
    integer :: p_max

    ! Precomputed coefficients
    real(dp), allocatable :: factorial(:)          ! 0! to (2*p_max)!
    real(dp), allocatable :: inv_factorial(:)      ! 1/n!
    real(dp), allocatable :: binomial(:,:)         ! C(n,k)
    real(dp), allocatable :: legendre_00(:)        ! P_n(0) for n=0..p_max
  end type multipole_system

contains

  !============================================================================
  ! Initialize system
  !============================================================================
  subroutine init_system(sys, xmax, ymax, nx, ny, p_max, max_part)
    type(multipole_system), intent(out) :: sys
    real(dp), intent(in) :: xmax, ymax
    integer, intent(in) :: nx, ny, p_max, max_part
    integer :: i, j, n, k
    real(dp) :: fact_n

    ! Domain parameters
    sys%xmax = xmax
    sys%ymax = ymax
    sys%nx_cells = nx
    sys%ny_cells = ny
    sys%dx_cell = 2.0_dp * xmax / real(nx, dp)
    sys%dy_cell = 2.0_dp * ymax / real(ny, dp)
    sys%p_max = p_max

    ! Allocate particles
    sys%max_particles = max_part
    sys%n_particles = 0
    allocate(sys%particles(max_part))

    ! Allocate cells
    allocate(sys%cells(nx, ny))
    do j = 1, ny
      do i = 1, nx
        sys%cells(i,j)%n_particles = 0
        sys%cells(i,j)%capacity = 10
        allocate(sys%cells(i,j)%particle_ids(10))
        allocate(sys%cells(i,j)%M_cos(0:p_max, 0:p_max))
        allocate(sys%cells(i,j)%M_sin(0:p_max, 1:p_max))
        sys%cells(i,j)%M_cos = 0.0_dp
        sys%cells(i,j)%M_sin = 0.0_dp

        ! Cell center
        sys%cells(i,j)%center_x = -xmax + (real(i, dp) - 0.5_dp) * sys%dx_cell
        sys%cells(i,j)%center_y = -ymax + (real(j, dp) - 0.5_dp) * sys%dy_cell
      end do
    end do

    ! Precompute mathematical constants
    allocate(sys%factorial(0:2*p_max+5))
    allocate(sys%inv_factorial(0:2*p_max+5))
    allocate(sys%binomial(0:2*p_max, 0:2*p_max))
    allocate(sys%legendre_00(0:p_max))

    ! Factorials
    fact_n = 1.0_dp
    sys%factorial(0) = 1.0_dp
    sys%inv_factorial(0) = 1.0_dp
    do n = 1, 2*p_max+5
      fact_n = fact_n * real(n, dp)
      sys%factorial(n) = fact_n
      sys%inv_factorial(n) = 1.0_dp / fact_n
    end do

    ! Binomial coefficients
    do n = 0, 2*p_max
      do k = 0, n
        sys%binomial(n,k) = sys%factorial(n) / (sys%factorial(k) * sys%factorial(n-k))
      end do
    end do

    ! Legendre P_n(0) values (many are zero!)
    do n = 0, p_max
      if (mod(n, 2) == 0) then
        ! P_n(0) = (-1)^(n/2) * (n)! / (2^n * ((n/2)!)^2)
        sys%legendre_00(n) = (-1.0_dp)**(n/2) * sys%factorial(n) / &
          (2.0_dp**n * sys%factorial(n/2)**2)
      else
        sys%legendre_00(n) = 0.0_dp
      end if
    end do

    print *, "System initialized:"
    print *, "  Domain: [", -xmax, ",", xmax, "] x [", -ymax, ",", ymax, "]"
    print *, "  Grid: ", nx, "x", ny
    print *, "  p_max: ", p_max
    print *, "  max_particles: ", max_part
  end subroutine init_system

  !============================================================================
  ! Destroy system
  !============================================================================
  subroutine destroy_system(sys)
    type(multipole_system), intent(inout) :: sys
    integer :: i, j

    if (allocated(sys%particles)) deallocate(sys%particles)
    if (allocated(sys%factorial)) deallocate(sys%factorial)
    if (allocated(sys%inv_factorial)) deallocate(sys%inv_factorial)
    if (allocated(sys%binomial)) deallocate(sys%binomial)
    if (allocated(sys%legendre_00)) deallocate(sys%legendre_00)

    if (allocated(sys%cells)) then
      do j = 1, sys%ny_cells
        do i = 1, sys%nx_cells
          if (allocated(sys%cells(i,j)%particle_ids)) deallocate(sys%cells(i,j)%particle_ids)
          if (allocated(sys%cells(i,j)%M_cos)) deallocate(sys%cells(i,j)%M_cos)
          if (allocated(sys%cells(i,j)%M_sin)) deallocate(sys%cells(i,j)%M_sin)
        end do
      end do
      deallocate(sys%cells)
    end if
  end subroutine destroy_system

  !============================================================================
  ! Add particle to system
  !============================================================================
  subroutine add_particle(sys, q, x, y)
    type(multipole_system), intent(inout) :: sys
    real(dp), intent(in) :: q, x, y
    integer :: ci, cj
    type(particle), allocatable :: temp(:)

    ! Resize if needed
    if (sys%n_particles >= sys%max_particles) then
      sys%max_particles = sys%max_particles * 2
      allocate(temp(sys%max_particles))
      temp(1:sys%n_particles) = sys%particles(1:sys%n_particles)
      call move_alloc(temp, sys%particles)
    end if

    sys%n_particles = sys%n_particles + 1
    sys%particles(sys%n_particles)%q = q
    sys%particles(sys%n_particles)%x = x
    sys%particles(sys%n_particles)%y = y
    sys%particles(sys%n_particles)%fx = 0.0_dp
    sys%particles(sys%n_particles)%fy = 0.0_dp

    ! Assign to cell
    call get_cell_indices(sys, x, y, ci, cj)
    sys%particles(sys%n_particles)%cell_i = ci
    sys%particles(sys%n_particles)%cell_j = cj
    call add_to_cell(sys%cells(ci, cj), sys%n_particles)
  end subroutine add_particle

  !============================================================================
  ! Get cell indices for position (x,y)
  !============================================================================
  subroutine get_cell_indices(sys, x, y, ci, cj)
    type(multipole_system), intent(in) :: sys
    real(dp), intent(in) :: x, y
    integer, intent(out) :: ci, cj

    ci = int((x + sys%xmax) / sys%dx_cell) + 1
    cj = int((y + sys%ymax) / sys%dy_cell) + 1
    ci = max(1, min(sys%nx_cells, ci))
    cj = max(1, min(sys%ny_cells, cj))
  end subroutine get_cell_indices

  !============================================================================
  ! Add particle to cell
  !============================================================================
  subroutine add_to_cell(c, pid)
    type(cell), intent(inout) :: c
    integer, intent(in) :: pid
    integer, allocatable :: temp(:)

    if (c%n_particles >= c%capacity) then
      c%capacity = c%capacity * 2
      allocate(temp(c%capacity))
      temp(1:c%n_particles) = c%particle_ids(1:c%n_particles)
      call move_alloc(temp, c%particle_ids)
    end if

    c%n_particles = c%n_particles + 1
    c%particle_ids(c%n_particles) = pid
  end subroutine add_to_cell

  !============================================================================
  ! Compute multipole moments for all cells
  ! M_n^m = sum_i q_i * r_i^n * [cos(m*phi_i), sin(m*phi_i)]
  !============================================================================
  subroutine compute_multipole_moments(sys)
    type(multipole_system), intent(inout) :: sys
    integer :: i, j, k, pid, n, m
    real(dp) :: dx, dy, r, r_n, phi, q
    real(dp) :: cos_m_phi, sin_m_phi

    ! Reset moments
    do j = 1, sys%ny_cells
      do i = 1, sys%nx_cells
        sys%cells(i,j)%M_cos = 0.0_dp
        sys%cells(i,j)%M_sin = 0.0_dp
      end do
    end do

    ! Compute moments
    do j = 1, sys%ny_cells
      do i = 1, sys%nx_cells
        do k = 1, sys%cells(i,j)%n_particles
          pid = sys%cells(i,j)%particle_ids(k)
          q = sys%particles(pid)%q

          ! Position relative to cell center
          dx = sys%particles(pid)%x - sys%cells(i,j)%center_x
          dy = sys%particles(pid)%y - sys%cells(i,j)%center_y
          r = sqrt(dx*dx + dy*dy)

          if (r < 1.0e-12_dp) then
            ! Monopole only
            sys%cells(i,j)%M_cos(0,0) = sys%cells(i,j)%M_cos(0,0) + q
            cycle
          end if

          phi = atan2(dy, dx)

          ! Accumulate moments
          r_n = 1.0_dp
          do n = 0, sys%p_max
            do m = 0, n
              if (m == 0) then
                ! m=0: only cosine term (azimuthally symmetric part)
                sys%cells(i,j)%M_cos(n,0) = sys%cells(i,j)%M_cos(n,0) + &
                  q * r_n * sys%legendre_00(n)
              else
                ! m>0: both sine and cosine
                cos_m_phi = cos(real(m, dp) * phi)
                sin_m_phi = sin(real(m, dp) * phi)
                sys%cells(i,j)%M_cos(n,m) = sys%cells(i,j)%M_cos(n,m) + q * r_n * cos_m_phi
                sys%cells(i,j)%M_sin(n,m) = sys%cells(i,j)%M_sin(n,m) + q * r_n * sin_m_phi
              end if
            end do
            r_n = r_n * r
          end do
        end do
      end do
    end do
  end subroutine compute_multipole_moments

  !============================================================================
  ! Compute total energy
  !============================================================================
  function compute_energy(sys, use_multipole) result(energy)
    type(multipole_system), intent(inout) :: sys
    logical, intent(in) :: use_multipole
    real(dp) :: energy

    if (use_multipole) then
      call compute_multipole_moments(sys)
      energy = compute_energy_multipole(sys)
    else
      energy = compute_energy_direct(sys)
    end if
  end function compute_energy

  !============================================================================
  ! Direct summation energy
  !============================================================================
  function compute_energy_direct(sys) result(energy)
    type(multipole_system), intent(in) :: sys
    real(dp) :: energy
    integer :: i, j
    real(dp) :: dx, dy, r

    energy = 0.0_dp
    do i = 1, sys%n_particles - 1
      do j = i + 1, sys%n_particles
        dx = sys%particles(i)%x - sys%particles(j)%x
        dy = sys%particles(i)%y - sys%particles(j)%y
        r = sqrt(dx*dx + dy*dy)
        if (r > 1.0e-12_dp) then
          energy = energy + sys%particles(i)%q * sys%particles(j)%q / r
        end if
      end do
    end do
  end function compute_energy_direct

  !============================================================================
  ! Multipole energy (simplified: cell-cell interactions)
  !============================================================================
  function compute_energy_multipole(sys) result(energy)
    type(multipole_system), intent(in) :: sys
    real(dp) :: energy
    integer :: i1, j1, i2, j2

    energy = 0.0_dp

    ! Self energy (within cells) - use direct
    do j1 = 1, sys%ny_cells
      do i1 = 1, sys%nx_cells
        energy = energy + cell_self_energy(sys, i1, j1)
      end do
    end do

    ! Cell-cell interactions
    do j1 = 1, sys%ny_cells
      do i1 = 1, sys%nx_cells
        do j2 = j1, sys%ny_cells
          do i2 = 1, sys%nx_cells
            if (i2 == i1 .and. j2 == j1) cycle
            if (j2 == j1 .and. i2 < i1) cycle

            ! Use multipole if well-separated
            if (abs(i2-i1) > 1 .or. abs(j2-j1) > 1) then
              energy = energy + cell_cell_energy_multipole(sys, i1, j1, i2, j2)
            else
              energy = energy + 0.5_dp * cell_cell_energy_direct(sys, i1, j1, i2, j2)
            end if
          end do
        end do
      end do
    end do
  end function compute_energy_multipole

  !============================================================================
  ! Cell self energy (direct)
  !============================================================================
  function cell_self_energy(sys, ci, cj) result(energy)
    type(multipole_system), intent(in) :: sys
    integer, intent(in) :: ci, cj
    real(dp) :: energy
    integer :: i, j, pi, pj
    real(dp) :: dx, dy, r

    energy = 0.0_dp
    do i = 1, sys%cells(ci,cj)%n_particles - 1
      pi = sys%cells(ci,cj)%particle_ids(i)
      do j = i+1, sys%cells(ci,cj)%n_particles
        pj = sys%cells(ci,cj)%particle_ids(j)
        dx = sys%particles(pi)%x - sys%particles(pj)%x
        dy = sys%particles(pi)%y - sys%particles(pj)%y
        r = sqrt(dx*dx + dy*dy)
        if (r > 1.0e-12_dp) then
          energy = energy + sys%particles(pi)%q * sys%particles(pj)%q / r
        end if
      end do
    end do
  end function cell_self_energy

  !============================================================================
  ! Cell-cell energy (direct)
  !============================================================================
  function cell_cell_energy_direct(sys, ci1, cj1, ci2, cj2) result(energy)
    type(multipole_system), intent(in) :: sys
    integer, intent(in) :: ci1, cj1, ci2, cj2
    real(dp) :: energy
    integer :: i, j, pi, pj
    real(dp) :: dx, dy, r

    energy = 0.0_dp
    do i = 1, sys%cells(ci1,cj1)%n_particles
      pi = sys%cells(ci1,cj1)%particle_ids(i)
      do j = 1, sys%cells(ci2,cj2)%n_particles
        pj = sys%cells(ci2,cj2)%particle_ids(j)
        dx = sys%particles(pi)%x - sys%particles(pj)%x
        dy = sys%particles(pi)%y - sys%particles(pj)%y
        r = sqrt(dx*dx + dy*dy)
        if (r > 1.0e-12_dp) then
          energy = energy + sys%particles(pi)%q * sys%particles(pj)%q / r
        end if
      end do
    end do
  end function cell_cell_energy_direct

  !============================================================================
  ! Cell-cell energy using multipole expansion
  ! Phi(R) ~= sum_n sum_m M_n^m / R^(n+1) * Y_n^m(Theta, Phi)
  !============================================================================
  function cell_cell_energy_multipole(sys, ci1, cj1, ci2, cj2) result(energy)
    type(multipole_system), intent(in) :: sys
    integer, intent(in) :: ci1, cj1, ci2, cj2
    real(dp) :: energy
    real(dp) :: dx, dy, R, phi_R, inv_R, inv_R_n
    integer :: n, m
    real(dp) :: cos_m_phi, sin_m_phi, phi_contrib

    energy = 0.0_dp

    ! Distance between cell centers
    dx = sys%cells(ci1,cj1)%center_x - sys%cells(ci2,cj2)%center_x
    dy = sys%cells(ci1,cj1)%center_y - sys%cells(ci2,cj2)%center_y
    R = sqrt(dx*dx + dy*dy)

    if (R < 1.0e-12_dp) return

    phi_R = atan2(dy, dx)
    inv_R = 1.0_dp / R

    ! Potential from cell2 at cell1's center
    phi_contrib = 0.0_dp
    inv_R_n = inv_R

    do n = 0, sys%p_max
      do m = 0, n
        if (m == 0) then
          phi_contrib = phi_contrib + sys%cells(ci2,cj2)%M_cos(n,0) * inv_R_n
        else
          cos_m_phi = cos(real(m,dp) * phi_R)
          sin_m_phi = sin(real(m,dp) * phi_R)
          phi_contrib = phi_contrib + inv_R_n * &
            (sys%cells(ci2,cj2)%M_cos(n,m) * cos_m_phi + &
             sys%cells(ci2,cj2)%M_sin(n,m) * sin_m_phi)
        end if
      end do
      inv_R_n = inv_R_n * inv_R
    end do

    ! Energy = Q1 * Phi(R1)
    energy = sys%cells(ci1,cj1)%M_cos(0,0) * phi_contrib
  end function cell_cell_energy_multipole

  !============================================================================
  ! Compute force on particle
  !============================================================================
  subroutine compute_force(sys, pid, fx, fy, use_multipole)
    type(multipole_system), intent(inout) :: sys
    integer, intent(in) :: pid
    real(dp), intent(out) :: fx, fy
    logical, intent(in) :: use_multipole

    if (use_multipole) then
      call compute_multipole_moments(sys)
      call compute_force_multipole(sys, pid, fx, fy)
    else
      call compute_force_direct(sys, pid, fx, fy)
    end if
  end subroutine compute_force

  !============================================================================
  ! Direct force calculation
  !============================================================================
  subroutine compute_force_direct(sys, pid, fx, fy)
    type(multipole_system), intent(in) :: sys
    integer, intent(in) :: pid
    real(dp), intent(out) :: fx, fy
    integer :: j
    real(dp) :: dx, dy, r, r3, f_mag

    fx = 0.0_dp
    fy = 0.0_dp

    do j = 1, sys%n_particles
      if (j == pid) cycle

      dx = sys%particles(pid)%x - sys%particles(j)%x
      dy = sys%particles(pid)%y - sys%particles(j)%y
      r = sqrt(dx*dx + dy*dy)

      if (r > 1.0e-12_dp) then
        r3 = r * r * r
        f_mag = sys%particles(pid)%q * sys%particles(j)%q / r3
        fx = fx + f_mag * dx
        fy = fy + f_mag * dy
      end if
    end do
  end subroutine compute_force_direct

  !============================================================================
  ! Multipole force calculation
  ! F = -q * grad(Phi)
  ! grad(Phi) in 2D using multipoles
  !============================================================================
  subroutine compute_force_multipole(sys, pid, fx, fy)
    type(multipole_system), intent(in) :: sys
    integer, intent(in) :: pid
    real(dp), intent(out) :: fx, fy
    integer :: ci, cj, i, j, k, p_other
    real(dp) :: dx, dy, r, r3, f_mag

    fx = 0.0_dp
    fy = 0.0_dp

    ci = sys%particles(pid)%cell_i
    cj = sys%particles(pid)%cell_j

    ! Near-field: direct calculation for same and neighbor cells
    do j = max(1, cj-1), min(sys%ny_cells, cj+1)
      do i = max(1, ci-1), min(sys%nx_cells, ci+1)
        do k = 1, sys%cells(i,j)%n_particles
          p_other = sys%cells(i,j)%particle_ids(k)
          if (p_other == pid) cycle

          dx = sys%particles(pid)%x - sys%particles(p_other)%x
          dy = sys%particles(pid)%y - sys%particles(p_other)%y
          r = sqrt(dx*dx + dy*dy)

          if (r > 1.0e-12_dp) then
            r3 = r * r * r
            f_mag = sys%particles(pid)%q * sys%particles(p_other)%q / r3
            fx = fx + f_mag * dx
            fy = fy + f_mag * dy
          end if
        end do
      end do
    end do

    ! Far-field: multipole expansion
    do j = 1, sys%ny_cells
      do i = 1, sys%nx_cells
        if (abs(i-ci) > 1 .or. abs(j-cj) > 1) then
          call add_far_field_force(sys, pid, i, j, fx, fy)
        end if
      end do
    end do
  end subroutine compute_force_multipole

  !============================================================================
  ! Add far-field force contribution from cell (ci,cj)
  !============================================================================
  subroutine add_far_field_force(sys, pid, ci, cj, fx, fy)
    type(multipole_system), intent(in) :: sys
    integer, intent(in) :: pid, ci, cj
    real(dp), intent(inout) :: fx, fy
    real(dp) :: dx, dy, r, phi, inv_r_n, inv_r_np1
    integer :: n, m
    real(dp) :: cos_m_phi, sin_m_phi
    real(dp) :: dphi_dx, dphi_dy, q

    dx = sys%particles(pid)%x - sys%cells(ci,cj)%center_x
    dy = sys%particles(pid)%y - sys%cells(ci,cj)%center_y
    r = sqrt(dx*dx + dy*dy)

    if (r < 1.0e-12_dp) return

    phi = atan2(dy, dx)
    q = sys%particles(pid)%q

    ! grad(Phi) = sum_n sum_m -(n+1) * M_n^m / r^(n+2) * grad(r^(-1) * Y_n^m)
    dphi_dx = 0.0_dp
    dphi_dy = 0.0_dp

    inv_r_n = 1.0_dp / r
    do n = 0, sys%p_max
      inv_r_np1 = inv_r_n / r

      do m = 0, n
        cos_m_phi = cos(real(m,dp) * phi)
        sin_m_phi = sin(real(m,dp) * phi)

        ! Simplified gradient (monopole approximation for far field)
        if (m == 0) then
          dphi_dx = dphi_dx - sys%cells(ci,cj)%M_cos(n,0) * real(n+1,dp) * inv_r_np1 * (dx/r)
          dphi_dy = dphi_dy - sys%cells(ci,cj)%M_cos(n,0) * real(n+1,dp) * inv_r_np1 * (dy/r)
        end if
      end do

      inv_r_n = inv_r_np1
    end do

    fx = fx - q * dphi_dx
    fy = fy - q * dphi_dy
  end subroutine add_far_field_force

  !============================================================================
  ! Compute forces for all particles
  !============================================================================
  subroutine compute_all_forces(sys, use_multipole)
    type(multipole_system), intent(inout) :: sys
    logical, intent(in) :: use_multipole
    integer :: i

    if (use_multipole) call compute_multipole_moments(sys)

    do i = 1, sys%n_particles
      if (use_multipole) then
        call compute_force_multipole(sys, i, sys%particles(i)%fx, sys%particles(i)%fy)
      else
        call compute_force_direct(sys, i, sys%particles(i)%fx, sys%particles(i)%fy)
      end if
    end do
  end subroutine compute_all_forces

  !============================================================================
  ! Move particles: dr/dt = eta * F
  !============================================================================
  subroutine move_particles(sys, eta, dt)
    type(multipole_system), intent(inout) :: sys
    real(dp), intent(in) :: eta, dt
    integer :: i
    real(dp) :: dx, dy, displacement
    real(dp), parameter :: MAX_DISP = 0.01_dp

    do i = 1, sys%n_particles
      dx = eta * sys%particles(i)%fx * dt
      dy = eta * sys%particles(i)%fy * dt

      ! Limit displacement
      displacement = sqrt(dx*dx + dy*dy)
      if (displacement > MAX_DISP) then
        dx = dx * MAX_DISP / displacement
        dy = dy * MAX_DISP / displacement
      end if

      sys%particles(i)%x = sys%particles(i)%x + dx
      sys%particles(i)%y = sys%particles(i)%y + dy

      ! Reflective boundaries
      if (sys%particles(i)%x < -sys%xmax) then
        sys%particles(i)%x = -2.0_dp*sys%xmax - sys%particles(i)%x
      else if (sys%particles(i)%x > sys%xmax) then
        sys%particles(i)%x = 2.0_dp*sys%xmax - sys%particles(i)%x
      end if

      if (sys%particles(i)%y < -sys%ymax) then
        sys%particles(i)%y = -2.0_dp*sys%ymax - sys%particles(i)%y
      else if (sys%particles(i)%y > sys%ymax) then
        sys%particles(i)%y = 2.0_dp*sys%ymax - sys%particles(i)%y
      end if
    end do
  end subroutine move_particles

  !============================================================================
  ! Update cell assignments
  !============================================================================
  subroutine update_cells(sys)
    type(multipole_system), intent(inout) :: sys
    integer :: i, ci, cj

    ! Clear cells
    do cj = 1, sys%ny_cells
      do ci = 1, sys%nx_cells
        sys%cells(ci,cj)%n_particles = 0
      end do
    end do

    ! Reassign particles
    do i = 1, sys%n_particles
      call get_cell_indices(sys, sys%particles(i)%x, sys%particles(i)%y, ci, cj)
      sys%particles(i)%cell_i = ci
      sys%particles(i)%cell_j = cj
      call add_to_cell(sys%cells(ci,cj), i)
    end do
  end subroutine update_cells

  !============================================================================
  ! Write frame to CSV
  !============================================================================
  subroutine write_frame(sys, frame_num, output_dir)
    type(multipole_system), intent(in) :: sys
    integer, intent(in) :: frame_num
    character(len=*), intent(in) :: output_dir
    character(len=512) :: filename
    integer :: unit, i

    write(filename, '(A,A,I0.6,A)') trim(output_dir), '/frame_', frame_num, '.csv'
    open(newunit=unit, file=trim(filename), status='replace')
    write(unit, '(A)') 'x,y,charge'
    do i = 1, sys%n_particles
      write(unit, '(F12.6,A,F12.6,A,F12.6)') &
        sys%particles(i)%x, ',', sys%particles(i)%y, ',', sys%particles(i)%q
    end do
    close(unit)
  end subroutine write_frame

end module multipole_module
