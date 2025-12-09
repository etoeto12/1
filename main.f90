program main
  !============================================================================
  ! Multipole simulation for 1/r potential in z=0 plane
  ! Simplified and optimized version
  !============================================================================
  use multipole_module
  implicit none

  type(multipole_system) :: sys
  integer :: i, frame
  real(8) :: t_start, t_end, time_direct, time_multipole
  real(8) :: energy_direct, energy_multipole, rel_error
  real(8) :: force_direct(2), force_multipole(2), force_error

  !============================================================================
  ! SIMULATION PARAMETERS - All controlled from here
  !============================================================================

  ! Domain size
  real(8), parameter :: XMAX = 1.0d0
  real(8), parameter :: YMAX = 1.0d0

  ! Grid parameters
  integer, parameter :: NX_CELLS = 10
  integer, parameter :: NY_CELLS = 10

  ! Multipole expansion order (adjust for accuracy/speed tradeoff)
  ! p=5:  fast, ~10^-3 accuracy
  ! p=10: balanced, ~10^-5 accuracy
  ! p=20: accurate, ~10^-8 accuracy
  integer, parameter :: P_MAX = 10

  ! Number of particles
  integer, parameter :: N_PARTICLES = 10000

  ! Dynamics parameters
  real(8), parameter :: ETA = 1.0d0        ! Mobility
  real(8), parameter :: DT = 0.001d0       ! Time step
  integer, parameter :: NSTEPS = 100       ! Number of steps

  ! Simulation control
  logical, parameter :: TEST_ACCURACY = .true.   ! Run accuracy test
  logical, parameter :: RUN_DYNAMICS = .true.    ! Run dynamics
  integer, parameter :: FRAME_SKIP = 5           ! Write every N frames

  ! Output directory - will be created based on parameters
  character(len=256) :: output_dir

  !============================================================================
  ! INITIALIZATION
  !============================================================================

  print *, "======================================================================"
  print *, "Multipole Simulation: 1/r Potential (z=0 plane)"
  print *, "======================================================================"
  print *, ""
  print *, "Parameters:"
  print *, "  Domain:     [", -XMAX, ",", XMAX, "] x [", -YMAX, ",", YMAX, "]"
  print *, "  Grid:       ", NX_CELLS, "x", NY_CELLS
  print *, "  p_max:      ", P_MAX
  print *, "  Particles:  ", N_PARTICLES
  print *, "  Timestep:   ", DT
  print *, "  Steps:      ", NSTEPS
  print *, ""

  ! Initialize system
  call init_system(sys, XMAX, YMAX, NX_CELLS, NY_CELLS, P_MAX, N_PARTICLES*2)

  !============================================================================
  ! GENERATE PARTICLES
  !============================================================================

  print *, "Generating particles..."
  call generate_random_particles(sys, N_PARTICLES)
  print *, "  Total particles: ", sys%n_particles
  print *, ""

  !============================================================================
  ! ACCURACY TEST
  !============================================================================

  if (TEST_ACCURACY) then
    print *, "======================================================================"
    print *, "ACCURACY TEST"
    print *, "======================================================================"

    ! Energy test
    if (sys%n_particles <= 5000) then
      print *, ""
      print *, "Energy comparison:"
      call cpu_time(t_start)
      energy_direct = compute_energy(sys, .false.)
      call cpu_time(t_end)
      time_direct = t_end - t_start
      print *, "  Direct:    ", energy_direct, " (", time_direct*1000, " ms)"

      call cpu_time(t_start)
      energy_multipole = compute_energy(sys, .true.)
      call cpu_time(t_end)
      time_multipole = t_end - t_start
      print *, "  Multipole: ", energy_multipole, " (", time_multipole*1000, " ms)"

      if (abs(energy_direct) > 1.0d-10) then
        rel_error = abs(energy_direct - energy_multipole) / abs(energy_direct)
        print *, "  Relative error: ", rel_error
        print *, "  Speedup:        ", time_direct / time_multipole, "x"
      end if
    else
      print *, "Skipping energy test (N too large for direct O(N^2) calculation)"
    end if

    ! Force test (sample 5 particles)
    print *, ""
    print *, "Force comparison (first 5 particles):"
    print *, "  ID  |    Direct Force    |  Multipole Force   | Rel. Error"
    print *, "  ---------------------------------------------------------------"
    do i = 1, min(5, sys%n_particles)
      call compute_force(sys, i, force_direct(1), force_direct(2), .false.)
      call compute_force(sys, i, force_multipole(1), force_multipole(2), .true.)

      force_error = sqrt((force_direct(1) - force_multipole(1))**2 + &
                        (force_direct(2) - force_multipole(2))**2) / &
                   max(sqrt(force_direct(1)**2 + force_direct(2)**2), 1.0d-10)

      write(*, '(I5,A,2E11.3,A,2E11.3,A,E10.2)') &
        i, '  | ', force_direct, '  | ', force_multipole, '  | ', force_error
    end do
    print *, ""
  end if

  !============================================================================
  ! DYNAMICS SIMULATION
  !============================================================================

  if (RUN_DYNAMICS) then
    print *, "======================================================================"
    print *, "DYNAMICS SIMULATION"
    print *, "======================================================================"
    print *, ""

    ! Create output directory with informative name
    write(output_dir, '(A,I0,A,I0,A,I0,A,I0)') &
      'output_N', N_PARTICLES, '_grid', NX_CELLS, '_p', P_MAX, '_steps', NSTEPS

    ! Create directory (Fortran-compatible way)
    call execute_command_line('mkdir ' // trim(output_dir), wait=.true.)
    print *, "Output directory: ", trim(output_dir)
    print *, ""

    ! Initial frame
    call write_frame(sys, 0, output_dir)

    print *, "Running simulation..."
    do frame = 1, NSTEPS
      ! Compute forces
      call compute_all_forces(sys, .true.)

      ! Update positions
      call move_particles(sys, ETA, DT)

      ! Update cell assignments
      call update_cells(sys)

      ! Write output
      if (mod(frame, FRAME_SKIP) == 0) then
        call write_frame(sys, frame, OUTPUT_DIR)
        write(*, '(A,I6,A,I6,A,F8.3,A)') &
          '  Frame ', frame, ' / ', NSTEPS, ' (', &
          100.0d0*real(frame)/real(NSTEPS), '%)'
      end if
    end do

    print *, ""
    print *, "Simulation complete!"
    print *, "Output written to: ", trim(output_dir)
    print *, ""
    print *, "Visualize with: python3 visualize.py --dir ", trim(output_dir)
  end if

  !============================================================================
  ! CLEANUP
  !============================================================================

  call destroy_system(sys)

  print *, ""
  print *, "======================================================================"
  print *, "Done!"
  print *, "======================================================================"

contains

  !============================================================================
  ! Generate random particles
  !============================================================================
  subroutine generate_random_particles(sys, n)
    type(multipole_system), intent(inout) :: sys
    integer, intent(in) :: n
    integer :: i
    real(8) :: x, y, q

    call random_seed()

    do i = 1, n
      call random_number(x)
      call random_number(y)
      call random_number(q)

      ! Position: uniform in domain
      x = (x - 0.5d0) * 1.8d0 * XMAX
      y = (y - 0.5d0) * 1.8d0 * YMAX

      ! Charge: all positive (superconductor vortices)
      q = 1.0d0

      call add_particle(sys, q, x, y)
    end do
  end subroutine generate_random_particles

end program main
