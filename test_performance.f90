program test_performance
  !============================================================================
  ! Performance benchmark for multipole method
  !============================================================================
  use multipole_module
  implicit none

  type(multipole_system) :: sys
  integer :: i, n, p
  integer, dimension(6) :: n_values = [100, 500, 1000, 5000, 10000, 50000]
  integer, parameter :: P_MAX = 10
  real(8) :: t_start, t_end, time_direct, time_multipole
  real(8) :: energy_direct, energy_multipole
  real(8) :: x, y
  integer :: ncells

  print *, "========================================================================"
  print *, "PERFORMANCE BENCHMARK"
  print *, "========================================================================"
  print *, ""
  print *, "P_MAX =", P_MAX
  print *, ""
  print *, "    N  | Grid |  Direct (ms) | Multipole (ms) |  Speedup  | Particles/sec"
  print *, "------------------------------------------------------------------------"

  do i = 1, size(n_values)
    n = n_values(i)

    ! Choose grid size based on N
    if (n <= 500) then
      ncells = 5
    else if (n <= 5000) then
      ncells = 10
    else if (n <= 20000) then
      ncells = 20
    else
      ncells = 30
    end if

    ! Initialize
    call init_system(sys, 1.0d0, 1.0d0, ncells, ncells, P_MAX, n*2)

    ! Generate particles
    call random_seed()
    do while (sys%n_particles < n)
      call random_number(x)
      call random_number(y)
      x = (x - 0.5d0) * 1.8d0
      y = (y - 0.5d0) * 1.8d0
      call add_particle(sys, 1.0d0, x, y)
    end do

    ! Benchmark multipole (always)
    call cpu_time(t_start)
    energy_multipole = compute_energy(sys, .true.)
    call cpu_time(t_end)
    time_multipole = (t_end - t_start) * 1000.0d0  ! Convert to ms

    ! Benchmark direct (only for N <= 5000)
    if (n <= 5000) then
      call cpu_time(t_start)
      energy_direct = compute_energy(sys, .false.)
      call cpu_time(t_end)
      time_direct = (t_end - t_start) * 1000.0d0

      write(*, '(I6,A,I4,A,F12.2,A,F14.2,A,F10.2,A,F12.0)') &
        n, '  | ', ncells, ' | ', time_direct, '  | ', time_multipole, &
        '  | ', time_direct/time_multipole, 'x | ', &
        real(n, 8) / (time_multipole / 1000.0d0)
    else
      write(*, '(I6,A,I4,A,A14,A,F14.2,A,A10,A,F12.0)') &
        n, '  | ', ncells, ' | ', '      -      ', '  | ', time_multipole, &
        '  | ', '    -     ', ' | ', real(n, 8) / (time_multipole / 1000.0d0)
    end if

    call destroy_system(sys)
  end do

  print *, ""
  print *, "========================================================================"

end program test_performance
