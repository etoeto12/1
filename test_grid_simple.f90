program test_grid_simple
  use multipole_accurate
  implicit none

  integer, parameter :: dp = selected_real_kind(15, 307)
  type(multipole_system) :: sys
  integer :: i, n
  real(dp) :: x, y
  real(dp) :: energy_fmm
  integer :: seed_size
  integer, allocatable :: seed(:)

  print *, '============================================================'
  print *, 'Simple Uniform Grid FMM Test'
  print *, '============================================================'

  call random_seed(size=seed_size)
  allocate(seed(seed_size))
  seed = 12345
  call random_seed(put=seed)

  n = 100
  print *, 'N =', n
  print *, 'Grid = 10x10'
  print *, 'p_max = 15'
  print *, ''

  call init_system(sys, 1.0_dp, 1.0_dp, 10, 10, p_max=15, max_part=n)

  ! Generate particles
  do i = 1, n
    call random_number(x)
    call random_number(y)
    x = 1.8_dp * x - 0.9_dp
    y = 1.8_dp * y - 0.9_dp
    call add_particle(sys, x, y, 1.0_dp)
  end do

  print *, 'Particles generated'
  print *, ''

  ! Compute forces
  call compute_all_forces(sys)

  ! Show first 10 particles
  print *, 'Forces on first 10 particles:'
  do i = 1, min(10, n)
    write(*, '(A,I3,A,2F12.4)') '  Particle ', i, ': F = (', &
          sys%particles(i)%fx, sys%particles(i)%fy
  end do

  energy_fmm = compute_energy_multipole(sys)
  write(*, '(A,ES15.6)') 'Energy (FMM): ', energy_fmm

  ! Compare with direct
  energy_fmm = compute_energy_direct(sys)
  write(*, '(A,ES15.6)') 'Energy (direct): ', energy_fmm

  call destroy_system(sys)
  deallocate(seed)

  print *, ''
  print *, 'Test complete'

end program test_grid_simple
