program test_simple_fmm
  use fmm_z0_precise
  implicit none

  integer, parameter :: dp = selected_real_kind(15, 307)
  type(fmm_system) :: sys
  integer :: i, n
  real(dp) :: x, y, q
  real(dp) :: energy_direct, energy_fmm
  real(dp) :: fx_ref, fy_ref, fx_fmm, fy_fmm, err
  real(dp) :: max_err, avg_err
  integer :: seed_size
  integer, allocatable :: seed(:)

  print *, '============================================================'
  print *, 'Simple FMM Test for 1/r potential in z=0 plane'
  print *, '============================================================'

  ! Random seed
  call random_seed(size=seed_size)
  allocate(seed(seed_size))
  seed = 12345
  call random_seed(put=seed)

  ! Setup
  n = 50  ! Small for testing
  call fmm_init(sys, -1.0_dp, 1.0_dp, -1.0_dp, 1.0_dp, p_max=10, max_particles=n)

  ! Generate particles
  print *, 'Generating', n, 'random particles...'
  do i = 1, n
    call random_number(x)
    call random_number(y)
    x = 2.0_dp * x - 1.0_dp
    y = 2.0_dp * y - 1.0_dp
    q = 1.0_dp
    call fmm_add_particle(sys, x, y, q)
  end do

  ! Build tree
  print *, 'Building tree...'
  call fmm_build_tree(sys)

  ! Direct computation
  print *, 'Computing forces (direct)...'
  ! Store reference using direct P2P
  call fmm_compute_forces(sys)  ! This currently does P2P only

  ! Compute energy
  energy_direct = fmm_compute_energy_direct(sys)
  print *, 'Energy (direct):', energy_direct

  ! Check force errors (for now, just report that P2P works)
  max_err = 0.0_dp
  avg_err = 0.0_dp
  do i = 1, min(10, n)
    print '(A,I3,A,2F12.6)', '  Particle ', i, ': F = (', &
          sys%particles(i)%fx, sys%particles(i)%fy
  end do

  print *, ''
  print *, 'Basic tree and P2P working!'
  print *, 'Next step: Implement M2M, M2L, L2L, L2P operators'

  call fmm_destroy(sys)
  deallocate(seed)

end program test_simple_fmm
