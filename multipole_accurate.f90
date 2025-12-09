module multipole_accurate
  !============================================================================
  ! Accurate multipole expansion for 1/r potential in z=0 plane
  ! Using 2D complex moments (Laurent series expansion)
  ! Uniform grid implementation
  ! Target: 10^-12 accuracy
  !============================================================================
  implicit none
  private

  public :: multipole_system
  public :: init_system, destroy_system
  public :: add_particle, compute_all_forces
  public :: compute_energy_direct, compute_energy_multipole
  public :: move_particles, write_frame

  integer, parameter :: dp = selected_real_kind(15, 307)
  real(dp), parameter :: PI = 3.141592653589793238462643383279502884197_dp
  real(dp), parameter :: EPS = 1.0e-14_dp

  type :: particle
    real(dp) :: q, x, y
    real(dp) :: fx, fy, phi
    integer :: cell_i, cell_j
  end type

  type :: cell
    integer :: n_particles, capacity
    integer, allocatable :: particle_ids(:)
    real(dp) :: center_x, center_y

    ! Complex multipole moments M_k for k=0..p
    complex(dp), allocatable :: M(:)
    ! Complex local expansion L_k for k=0..p
    complex(dp), allocatable :: L(:)
  end type

  type :: multipole_system
    real(dp) :: xmax, ymax
    integer :: nx_cells, ny_cells
    real(dp) :: dx_cell, dy_cell
    integer :: p_max

    type(cell), allocatable :: cells(:,:)
    type(particle), allocatable :: particles(:)
    integer :: n_particles, max_particles

    ! Precomputed binomial coefficients
    real(dp), allocatable :: binomial(:,:)
  end type

contains

  subroutine init_system(sys, xmax, ymax, nx, ny, p_max, max_part)
    type(multipole_system), intent(out) :: sys
    real(dp), intent(in) :: xmax, ymax
    integer, intent(in) :: nx, ny, p_max, max_part
    integer :: i, j, n, k

    sys%xmax = xmax
    sys%ymax = ymax
    sys%nx_cells = nx
    sys%ny_cells = ny
    sys%dx_cell = 2.0_dp * xmax / real(nx, dp)
    sys%dy_cell = 2.0_dp * ymax / real(ny, dp)
    sys%p_max = p_max
    sys%max_particles = max_part
    sys%n_particles = 0

    allocate(sys%particles(max_part))
    allocate(sys%cells(nx, ny))

    ! Initialize cells
    do j = 1, ny
      do i = 1, nx
        sys%cells(i,j)%center_x = -xmax + (i - 0.5_dp) * sys%dx_cell
        sys%cells(i,j)%center_y = -ymax + (j - 0.5_dp) * sys%dy_cell
        sys%cells(i,j)%n_particles = 0
        sys%cells(i,j)%capacity = 100
        allocate(sys%cells(i,j)%particle_ids(100))
        allocate(sys%cells(i,j)%M(0:p_max))
        allocate(sys%cells(i,j)%L(0:p_max))
        sys%cells(i,j)%M = cmplx(0.0_dp, 0.0_dp, dp)
        sys%cells(i,j)%L = cmplx(0.0_dp, 0.0_dp, dp)
      end do
    end do

    ! Precompute binomial coefficients
    allocate(sys%binomial(0:2*p_max, 0:2*p_max))
    sys%binomial = 0.0_dp
    sys%binomial(0,0) = 1.0_dp
    do n = 1, 2*p_max
      sys%binomial(n,0) = 1.0_dp
      sys%binomial(n,n) = 1.0_dp
      do k = 1, n-1
        sys%binomial(n,k) = sys%binomial(n-1,k-1) + sys%binomial(n-1,k)
      end do
    end do

  end subroutine

  subroutine destroy_system(sys)
    type(multipole_system), intent(inout) :: sys
    integer :: i, j

    if (allocated(sys%particles)) deallocate(sys%particles)
    if (allocated(sys%binomial)) deallocate(sys%binomial)

    if (allocated(sys%cells)) then
      do j = 1, sys%ny_cells
        do i = 1, sys%nx_cells
          if (allocated(sys%cells(i,j)%particle_ids)) &
            deallocate(sys%cells(i,j)%particle_ids)
          if (allocated(sys%cells(i,j)%M)) deallocate(sys%cells(i,j)%M)
          if (allocated(sys%cells(i,j)%L)) deallocate(sys%cells(i,j)%L)
        end do
      end do
      deallocate(sys%cells)
    end if
  end subroutine

  subroutine add_particle(sys, x, y, q)
    type(multipole_system), intent(inout) :: sys
    real(dp), intent(in) :: x, y, q
    integer :: ci, cj
    integer, allocatable :: tmp(:)
    integer :: new_cap

    sys%n_particles = sys%n_particles + 1
    sys%particles(sys%n_particles)%x = x
    sys%particles(sys%n_particles)%y = y
    sys%particles(sys%n_particles)%q = q
    sys%particles(sys%n_particles)%fx = 0.0_dp
    sys%particles(sys%n_particles)%fy = 0.0_dp
    sys%particles(sys%n_particles)%phi = 0.0_dp

    ! Assign to cell
    ci = int((x + sys%xmax) / sys%dx_cell) + 1
    cj = int((y + sys%ymax) / sys%dy_cell) + 1
    ci = max(1, min(sys%nx_cells, ci))
    cj = max(1, min(sys%ny_cells, cj))

    sys%particles(sys%n_particles)%cell_i = ci
    sys%particles(sys%n_particles)%cell_j = cj

    ! Add to cell (resize if needed)
    if (sys%cells(ci,cj)%n_particles >= sys%cells(ci,cj)%capacity) then
      new_cap = sys%cells(ci,cj)%capacity * 2
      allocate(tmp(new_cap))
      tmp(1:sys%cells(ci,cj)%n_particles) = sys%cells(ci,cj)%particle_ids(1:sys%cells(ci,cj)%n_particles)
      deallocate(sys%cells(ci,cj)%particle_ids)
      allocate(sys%cells(ci,cj)%particle_ids(new_cap))
      sys%cells(ci,cj)%particle_ids(1:sys%cells(ci,cj)%n_particles) = tmp(1:sys%cells(ci,cj)%n_particles)
      sys%cells(ci,cj)%capacity = new_cap
      deallocate(tmp)
    end if

    sys%cells(ci,cj)%n_particles = sys%cells(ci,cj)%n_particles + 1
    sys%cells(ci,cj)%particle_ids(sys%cells(ci,cj)%n_particles) = sys%n_particles
  end subroutine

  !============================================================================
  ! Compute all forces using FMM
  !============================================================================
  subroutine compute_all_forces(sys)
    type(multipole_system), intent(inout) :: sys
    integer :: i

    ! Reset forces
    do i = 1, sys%n_particles
      sys%particles(i)%fx = 0.0_dp
      sys%particles(i)%fy = 0.0_dp
      sys%particles(i)%phi = 0.0_dp
    end do

    ! P2M: Compute multipole moments
    call p2m_all(sys)

    ! M2L: Far field interactions
    call m2l_all(sys)

    ! L2P: Evaluate local expansions
    call l2p_all(sys)

    ! P2P: Near field direct
    call p2p_all(sys)

  end subroutine

  !============================================================================
  ! P2M: Particle to Multipole
  !============================================================================
  subroutine p2m_all(sys)
    type(multipole_system), intent(inout) :: sys
    integer :: ci, cj, i, ip, k
    real(dp) :: dx, dy, q
    complex(dp) :: z, zk

    ! Reset moments
    do cj = 1, sys%ny_cells
      do ci = 1, sys%nx_cells
        sys%cells(ci,cj)%M = cmplx(0.0_dp, 0.0_dp, dp)
      end do
    end do

    ! Compute moments
    do cj = 1, sys%ny_cells
      do ci = 1, sys%nx_cells
        do i = 1, sys%cells(ci,cj)%n_particles
          ip = sys%cells(ci,cj)%particle_ids(i)
          q = sys%particles(ip)%q
          dx = sys%particles(ip)%x - sys%cells(ci,cj)%center_x
          dy = sys%particles(ip)%y - sys%cells(ci,cj)%center_y
          z = cmplx(dx, dy, dp)

          zk = cmplx(1.0_dp, 0.0_dp, dp)
          do k = 0, sys%p_max
            sys%cells(ci,cj)%M(k) = sys%cells(ci,cj)%M(k) + q * zk
            zk = zk * z
          end do
        end do
      end do
    end do
  end subroutine

  !============================================================================
  ! M2L: Multipole to Local (far field)
  !============================================================================
  subroutine m2l_all(sys)
    type(multipole_system), intent(inout) :: sys
    integer :: ci, cj, si, sj

    ! Reset local expansions
    do cj = 1, sys%ny_cells
      do ci = 1, sys%nx_cells
        sys%cells(ci,cj)%L = cmplx(0.0_dp, 0.0_dp, dp)
      end do
    end do

    ! Compute M2L for far cells
    do cj = 1, sys%ny_cells
      do ci = 1, sys%nx_cells
        do sj = 1, sys%ny_cells
          do si = 1, sys%nx_cells
            ! Skip near neighbors (will be handled by P2P)
            if (abs(si-ci) <= 1 .and. abs(sj-cj) <= 1) cycle

            call m2l_single(sys, si, sj, ci, cj)
          end do
        end do
      end do
    end do
  end subroutine

  subroutine m2l_single(sys, si, sj, ti, tj)
    type(multipole_system), intent(inout) :: sys
    integer, intent(in) :: si, sj, ti, tj
    integer :: k, j
    real(dp) :: dx, dy
    complex(dp) :: z0, z0inv, z0k, term

    dx = sys%cells(ti,tj)%center_x - sys%cells(si,sj)%center_x
    dy = sys%cells(ti,tj)%center_y - sys%cells(si,sj)%center_y
    z0 = cmplx(dx, dy, dp)

    if (abs(z0) < EPS) return

    z0inv = 1.0_dp / z0

    ! k=0 term
    sys%cells(ti,tj)%L(0) = sys%cells(ti,tj)%L(0) + &
                            sys%cells(si,sj)%M(0) * log(-z0)

    z0k = z0inv
    do j = 1, sys%p_max
      sys%cells(ti,tj)%L(0) = sys%cells(ti,tj)%L(0) - &
                              sys%cells(si,sj)%M(j) * z0k / real(j, dp)
      z0k = z0k * z0inv
    end do

    ! k>0 terms
    do k = 1, sys%p_max
      z0k = z0inv**k
      term = -sys%cells(si,sj)%M(0) * z0k / real(k, dp)

      do j = 1, sys%p_max
        term = term - sys%cells(si,sj)%M(j) * sys%binomial(j+k-1,k) * &
               z0k * z0inv**j
      end do

      sys%cells(ti,tj)%L(k) = sys%cells(ti,tj)%L(k) + term
    end do
  end subroutine

  !============================================================================
  ! L2P: Local to Particle
  !============================================================================
  subroutine l2p_all(sys)
    type(multipole_system), intent(inout) :: sys
    integer :: ci, cj, i, ip, k
    real(dp) :: dx, dy
    complex(dp) :: z, zk, phi, dphidz

    do cj = 1, sys%ny_cells
      do ci = 1, sys%nx_cells
        do i = 1, sys%cells(ci,cj)%n_particles
          ip = sys%cells(ci,cj)%particle_ids(i)
          dx = sys%particles(ip)%x - sys%cells(ci,cj)%center_x
          dy = sys%particles(ip)%y - sys%cells(ci,cj)%center_y
          z = cmplx(dx, dy, dp)

          phi = cmplx(0.0_dp, 0.0_dp, dp)
          dphidz = cmplx(0.0_dp, 0.0_dp, dp)

          zk = cmplx(1.0_dp, 0.0_dp, dp)
          do k = 0, sys%p_max
            phi = phi + sys%cells(ci,cj)%L(k) * zk
            if (k > 0) then
              dphidz = dphidz + real(k,dp) * sys%cells(ci,cj)%L(k) * zk / z
            end if
            zk = zk * z
          end do

          sys%particles(ip)%phi = sys%particles(ip)%phi + real(phi, dp)
          sys%particles(ip)%fx = sys%particles(ip)%fx - &
                                 sys%particles(ip)%q * real(dphidz, dp)
          sys%particles(ip)%fy = sys%particles(ip)%fy - &
                                 sys%particles(ip)%q * aimag(dphidz)
        end do
      end do
    end do
  end subroutine

  !============================================================================
  ! P2P: Particle to Particle (near field)
  !============================================================================
  subroutine p2p_all(sys)
    type(multipole_system), intent(inout) :: sys
    integer :: ci, cj, si, sj

    do cj = 1, sys%ny_cells
      do ci = 1, sys%nx_cells
        ! Self and immediate neighbors
        do sj = max(1,cj-1), min(sys%ny_cells,cj+1)
          do si = max(1,ci-1), min(sys%nx_cells,ci+1)
            call p2p_cells(sys, ci, cj, si, sj)
          end do
        end do
      end do
    end do
  end subroutine

  subroutine p2p_cells(sys, ci, cj, si, sj)
    type(multipole_system), intent(inout) :: sys
    integer, intent(in) :: ci, cj, si, sj
    integer :: i, j, ip, jp
    real(dp) :: dx, dy, r, r3, qi, qj, fx, fy
    logical :: same_cell

    same_cell = (ci == si .and. cj == sj)

    do i = 1, sys%cells(ci,cj)%n_particles
      ip = sys%cells(ci,cj)%particle_ids(i)
      qi = sys%particles(ip)%q

      do j = 1, sys%cells(si,sj)%n_particles
        jp = sys%cells(si,sj)%particle_ids(j)
        if (same_cell .and. ip >= jp) cycle

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

        if (same_cell) then
          sys%particles(jp)%fx = sys%particles(jp)%fx - fx
          sys%particles(jp)%fy = sys%particles(jp)%fy - fy
          sys%particles(jp)%phi = sys%particles(jp)%phi + qi / r
        end if
      end do
    end do
  end subroutine

  !============================================================================
  ! Direct energy computation (reference)
  !============================================================================
  function compute_energy_direct(sys) result(energy)
    type(multipole_system), intent(in) :: sys
    real(dp) :: energy
    integer :: i, j
    real(dp) :: dx, dy, r

    energy = 0.0_dp
    do i = 1, sys%n_particles
      do j = i+1, sys%n_particles
        dx = sys%particles(j)%x - sys%particles(i)%x
        dy = sys%particles(j)%y - sys%particles(i)%y
        r = sqrt(dx**2 + dy**2)
        if (r > EPS) then
          energy = energy + sys%particles(i)%q * sys%particles(j)%q / r
        end if
      end do
    end do
  end function

  function compute_energy_multipole(sys) result(energy)
    type(multipole_system), intent(in) :: sys
    real(dp) :: energy
    integer :: i

    energy = 0.0_dp
    do i = 1, sys%n_particles
      energy = energy + 0.5_dp * sys%particles(i)%q * sys%particles(i)%phi
    end do
  end function

  !============================================================================
  ! Move particles (simple Euler)
  !============================================================================
  subroutine move_particles(sys, dt, eta)
    type(multipole_system), intent(inout) :: sys
    real(dp), intent(in) :: dt, eta
    integer :: i, old_ci, old_cj, new_ci, new_cj, k
    real(dp) :: vx, vy, new_x, new_y, max_disp
    integer, allocatable :: tmp_ids(:)
    integer :: tmp_cap

    max_disp = 0.1_dp * min(sys%dx_cell, sys%dy_cell)

    do i = 1, sys%n_particles
      vx = eta * sys%particles(i)%fx
      vy = eta * sys%particles(i)%fy

      ! Limit displacement
      if (sqrt(vx**2 + vy**2) * dt > max_disp) then
        vx = vx * max_disp / (sqrt(vx**2 + vy**2) * dt)
        vy = vy * max_disp / (sqrt(vx**2 + vy**2) * dt)
      end if

      new_x = sys%particles(i)%x + vx * dt
      new_y = sys%particles(i)%y + vy * dt

      ! Reflecting boundaries
      if (new_x < -sys%xmax) new_x = -2.0_dp * sys%xmax - new_x
      if (new_x > sys%xmax) new_x = 2.0_dp * sys%xmax - new_x
      if (new_y < -sys%ymax) new_y = -2.0_dp * sys%ymax - new_y
      if (new_y > sys%ymax) new_y = 2.0_dp * sys%ymax - new_y

      old_ci = sys%particles(i)%cell_i
      old_cj = sys%particles(i)%cell_j

      sys%particles(i)%x = new_x
      sys%particles(i)%y = new_y

      ! Update cell assignment
      new_ci = int((new_x + sys%xmax) / sys%dx_cell) + 1
      new_cj = int((new_y + sys%ymax) / sys%dy_cell) + 1
      new_ci = max(1, min(sys%nx_cells, new_ci))
      new_cj = max(1, min(sys%ny_cells, new_cj))

      if (new_ci /= old_ci .or. new_cj /= old_cj) then
        ! Remove from old cell
        do k = 1, sys%cells(old_ci,old_cj)%n_particles
          if (sys%cells(old_ci,old_cj)%particle_ids(k) == i) then
            sys%cells(old_ci,old_cj)%particle_ids(k) = &
              sys%cells(old_ci,old_cj)%particle_ids(sys%cells(old_ci,old_cj)%n_particles)
            sys%cells(old_ci,old_cj)%n_particles = sys%cells(old_ci,old_cj)%n_particles - 1
            exit
          end if
        end do

        ! Add to new cell (resize if needed)
        if (sys%cells(new_ci,new_cj)%n_particles >= sys%cells(new_ci,new_cj)%capacity) then
          tmp_cap = sys%cells(new_ci,new_cj)%capacity * 2
          allocate(tmp_ids(tmp_cap))
          tmp_ids(1:sys%cells(new_ci,new_cj)%n_particles) = &
            sys%cells(new_ci,new_cj)%particle_ids(1:sys%cells(new_ci,new_cj)%n_particles)
          deallocate(sys%cells(new_ci,new_cj)%particle_ids)
          allocate(sys%cells(new_ci,new_cj)%particle_ids(tmp_cap))
          sys%cells(new_ci,new_cj)%particle_ids(1:sys%cells(new_ci,new_cj)%n_particles) = &
            tmp_ids(1:sys%cells(new_ci,new_cj)%n_particles)
          sys%cells(new_ci,new_cj)%capacity = tmp_cap
          deallocate(tmp_ids)
        end if
        sys%cells(new_ci,new_cj)%n_particles = sys%cells(new_ci,new_cj)%n_particles + 1
        sys%cells(new_ci,new_cj)%particle_ids(sys%cells(new_ci,new_cj)%n_particles) = i

        sys%particles(i)%cell_i = new_ci
        sys%particles(i)%cell_j = new_cj
      end if
    end do
  end subroutine

  !============================================================================
  ! Write frame
  !============================================================================
  subroutine write_frame(sys, filename)
    type(multipole_system), intent(in) :: sys
    character(len=*), intent(in) :: filename
    integer :: unit_num, i

    open(newunit=unit_num, file=trim(filename), status='replace')
    write(unit_num, '(A)') 'x,y,charge'
    do i = 1, sys%n_particles
      write(unit_num, '(F12.6,A,F12.6,A,F12.6)') &
        sys%particles(i)%x, ',', sys%particles(i)%y, ',', sys%particles(i)%q
    end do
    close(unit_num)
  end subroutine

end module multipole_accurate
