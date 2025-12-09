program multipole_simulation
  !============================================================================
  ! Multipole simulation for 1/r potential in z=0 plane
  ! Using uniform grid FMM with complex moments
  !============================================================================
  use multipole_accurate
  implicit none

  integer, parameter :: dp = selected_real_kind(15, 307)

  ! System parameters
  type(multipole_system) :: sys
  real(dp) :: box_size
  integer :: n_cells, p_order
  integer :: n_particles, max_particles

  ! Simulation parameters
  integer :: n_steps, write_interval
  real(dp) :: dt, eta

  ! Energy and diagnostics
  real(dp) :: energy_direct, energy_fmm, rel_error

  ! Loop variables
  integer :: i, step
  character(len=100) :: filename
  real(dp) :: x, y, charge
  real(dp) :: start_time, end_time

  !----------------------------------------------------------------------------
  ! Configuration
  !----------------------------------------------------------------------------

  ! Domain: [-box_size, box_size] x [-box_size, box_size]
  box_size = 10.0_dp

  ! Grid parameters
  n_cells = 5               ! 5x5 uniform grid
  p_order = 12              ! Multipole expansion order (higher = more accurate)

  ! Particles
  n_particles = 100
  max_particles = 1000

  ! Time integration
  n_steps = 100
  dt = 0.01_dp
  eta = 0.1_dp              ! Mobility coefficient
  write_interval = 10

  !----------------------------------------------------------------------------
  ! Initialize system
  !----------------------------------------------------------------------------

  print *, '============================================'
  print *, 'Multipole Simulation (Uniform Grid FMM)'
  print *, '============================================'
  print *, 'Box size:        ', box_size
  print *, 'Grid:            ', n_cells, 'x', n_cells
  print *, 'Multipole order: ', p_order
  print *, 'Particles:       ', n_particles
  print *, 'Time steps:      ', n_steps
  print *, '============================================'
  print *

  call init_system(sys, box_size, box_size, n_cells, n_cells, p_order, max_particles)

  ! Add particles with random positions and charges
  call random_seed()
  do i = 1, n_particles
    call random_number(x)
    call random_number(y)
    call random_number(charge)

    x = (x - 0.5_dp) * 2.0_dp * box_size * 0.9_dp  ! [-0.9*box_size, 0.9*box_size]
    y = (y - 0.5_dp) * 2.0_dp * box_size * 0.9_dp
    charge = (charge - 0.5_dp) * 2.0_dp            ! [-1, 1]

    call add_particle(sys, x, y, charge)
  end do

  print *, 'Added', sys%n_particles, 'particles'
  print *

  !----------------------------------------------------------------------------
  ! Initial energy check
  !----------------------------------------------------------------------------

  print *, 'Computing initial forces and energy...'
  call cpu_time(start_time)
  call compute_all_forces(sys)
  call cpu_time(end_time)

  energy_fmm = compute_energy_multipole(sys)
  energy_direct = compute_energy_direct(sys)
  rel_error = abs(energy_fmm - energy_direct) / abs(energy_direct)

  print *, 'FMM computation time:    ', end_time - start_time, 's'
  print *, 'Energy (direct):         ', energy_direct
  print *, 'Energy (FMM):            ', energy_fmm
  print *, 'Relative error:          ', rel_error
  print *

  !----------------------------------------------------------------------------
  ! Time integration loop
  !----------------------------------------------------------------------------

  print *, 'Starting time integration...'
  print *

  do step = 1, n_steps

    ! Compute forces
    call compute_all_forces(sys)

    ! Move particles
    call move_particles(sys, dt, eta)

    ! Write output
    if (mod(step, write_interval) == 0) then
      write(filename, '(A,I6.6,A)') 'particles_', step, '.csv'
      call write_frame(sys, filename)

      energy_fmm = compute_energy_multipole(sys)

      print '(A,I6,A,E14.6)', 'Step ', step, '  Energy: ', energy_fmm
    end if

  end do

  print *
  print *, '============================================'
  print *, 'Simulation complete!'
  print *, '============================================'
  print *

  ! Final diagnostics
  print *, 'Computing final accuracy...'
  call compute_all_forces(sys)

  energy_fmm = compute_energy_multipole(sys)
  energy_direct = compute_energy_direct(sys)
  rel_error = abs(energy_fmm - energy_direct) / abs(energy_direct)

  print *, 'Final energy (direct):   ', energy_direct
  print *, 'Final energy (FMM):      ', energy_fmm
  print *, 'Relative error:          ', rel_error
  print *

  ! Cleanup
  call destroy_system(sys)

  print *, 'Done!'

end program multipole_simulation
