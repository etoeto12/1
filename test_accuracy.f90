program test_accuracy
  !============================================================================
  ! Test program for multipole expansion accuracy
  ! Tests convergence with increasing p_max
  !============================================================================
  use multipole_module
  implicit none

  type(multipole_system) :: sys
  integer :: p, i, n_test
  integer, parameter :: N_PARTICLES = 50
  integer, parameter :: MAX_P = 30
  real(8) :: energy_direct, energy_multipole, rel_error
  real(8) :: force_direct(2), force_multipole(2), force_error
  real(8) :: max_force_error, avg_force_error
  real(8) :: t_start, t_end, time_direct, time_multipole
  real(8) :: x, y

  print *, "========================================================================"
  print *, "MULTIPOLE ACCURACY TEST"
  print *, "========================================================================"
  print *, ""
  print *, "Testing with N =", N_PARTICLES, "particles"
  print *, "Grid: 5x5 cells"
  print *, ""

  ! Test convergence with increasing p_max
  print *, "P_MAX | Energy Rel.Error | Max Force Error | Avg Force Error | Time Ratio"
  print *, "----------------------------------------------------------------------"

  do p = 2, MAX_P, 2
    ! Initialize system with current p_max
    call init_system(sys, 1.0d0, 1.0d0, 5, 5, p, N_PARTICLES*2)

    ! Generate fixed random configuration
    call random_seed(put=[(12345, i=1,8)])
    do i = 1, N_PARTICLES
      call random_number(x)
      call random_number(y)
      x = (x - 0.5d0) * 1.6d0
      y = (y - 0.5d0) * 1.6d0
      call add_particle(sys, 1.0d0, x, y)
    end do

    ! Energy test
    call cpu_time(t_start)
    energy_direct = compute_energy(sys, .false.)
    call cpu_time(t_end)
    time_direct = t_end - t_start

    call cpu_time(t_start)
    energy_multipole = compute_energy(sys, .true.)
    call cpu_time(t_end)
    time_multipole = t_end - t_start

    rel_error = abs(energy_direct - energy_multipole) / abs(energy_direct)

    ! Force test on sample particles
    max_force_error = 0.0d0
    avg_force_error = 0.0d0
    n_test = min(10, sys%n_particles)

    do i = 1, n_test
      call compute_force(sys, i, force_direct(1), force_direct(2), .false.)
      call compute_force(sys, i, force_multipole(1), force_multipole(2), .true.)

      force_error = sqrt((force_direct(1) - force_multipole(1))**2 + &
                        (force_direct(2) - force_multipole(2))**2) / &
                   max(sqrt(force_direct(1)**2 + force_direct(2)**2), 1.0d-10)

      max_force_error = max(max_force_error, force_error)
      avg_force_error = avg_force_error + force_error
    end do
    avg_force_error = avg_force_error / real(n_test, 8)

    write(*, '(I5,A,E13.5,A,E13.5,A,E13.5,A,F8.3)') &
      p, '  | ', rel_error, '  | ', max_force_error, '  | ', avg_force_error, &
      '  | ', time_multipole / time_direct

    call destroy_system(sys)

    ! Stop if we reached target accuracy
    if (rel_error < 1.0d-12 .and. max_force_error < 1.0d-12) then
      print *, ""
      print *, "*** TARGET ACCURACY 10^-12 ACHIEVED AT P_MAX =", p, " ***"
      exit
    end if
  end do

  print *, ""
  print *, "========================================================================"

end program test_accuracy
