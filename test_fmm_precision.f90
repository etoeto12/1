program test_fmm_precision
  !============================================================================
  ! Test FMM precision for 1/r potential in z=0 plane
  ! Target: achieve 10^-12 relative error in forces and energy
  !============================================================================
  use fmm3d_module
  implicit none

  integer, parameter :: dp = selected_real_kind(15, 307)
  type(fmm_system) :: sys
  real(dp) :: bbox(6)
  integer :: n_particles, p_max, i, j
  real(dp) :: x, y, q
  real(dp) :: energy_direct, energy_fmm
  real(dp) :: fx_direct, fy_direct, fx_fmm, fy_fmm
  real(dp) :: max_force_error, avg_force_error, energy_error
  real(dp) :: dx, dy, dz, r, r3, qi, qj
  real(dp), allocatable :: fx_ref(:), fy_ref(:), phi_ref(:)
  integer :: seed_size
  integer, allocatable :: seed(:)

  print *, '============================================================'
  print *, 'FMM High-Precision Test for 1/r Potential (z=0 plane)'
  print *, 'Target: 10^-12 relative error'
  print *, '============================================================'
  print *, ''

  ! Random seed for reproducibility
  call random_seed(size=seed_size)
  allocate(seed(seed_size))
  seed = 12345
  call random_seed(put=seed)

  ! Test configuration
  n_particles = 100  ! Small number for direct comparison
  p_max = 20         ! High order for precision

  print *, 'Configuration:'
  print *, '  N particles:', n_particles
  print *, '  p_max:', p_max
  print *, ''

  ! Domain: [-1, 1]^3 but particles only in z=0 plane
  bbox = [-1.0_dp, 1.0_dp, -1.0_dp, 1.0_dp, -0.1_dp, 0.1_dp]

  ! Initialize FMM with strict accuracy
  call fmm_init(sys, bbox, p_max, n_particles, theta=0.3_dp, max_leaf_size=20)

  ! Generate random particles in z=0 plane
  print *, 'Generating', n_particles, 'random particles in z=0 plane...'
  do i = 1, n_particles
    call random_number(x)
    call random_number(y)
    x = 2.0_dp * x - 1.0_dp  ! [-1, 1]
    y = 2.0_dp * y - 1.0_dp
    q = 1.0_dp  ! All positive charges

    call fmm_add_particle(sys, x, y, 0.0_dp, q)
  end do

  print *, 'Particles generated'
  print *, ''

  ! Allocate reference arrays
  allocate(fx_ref(n_particles), fy_ref(n_particles), phi_ref(n_particles))

  !--------------------------------------------------------------------------
  ! DIRECT SUMMATION (reference)
  !--------------------------------------------------------------------------
  print *, 'Computing reference solution (direct O(N^2) summation)...'

  fx_ref = 0.0_dp
  fy_ref = 0.0_dp
  phi_ref = 0.0_dp
  energy_direct = 0.0_dp

  do i = 1, n_particles
    qi = sys%particles(i)%q

    do j = 1, n_particles
      if (i == j) cycle

      qj = sys%particles(j)%q

      dx = sys%particles(j)%x - sys%particles(i)%x
      dy = sys%particles(j)%y - sys%particles(i)%y
      dz = sys%particles(j)%z - sys%particles(i)%z

      r = sqrt(dx**2 + dy**2 + dz**2)

      if (r < 1.0e-14_dp) cycle

      r3 = r**3

      ! Potential
      phi_ref(i) = phi_ref(i) + qj / r

      ! Force: F_i = q_i * Σ_j q_j * (r_j - r_i) / |r_j - r_i|^3
      fx_ref(i) = fx_ref(i) + qi * qj * dx / r3
      fy_ref(i) = fy_ref(i) + qi * qj * dy / r3
    end do

    ! Energy: E = 0.5 * Σ_i q_i * φ_i
    energy_direct = energy_direct + 0.5_dp * qi * phi_ref(i)
  end do

  print *, 'Direct summation complete'
  print *, '  Energy (direct):', energy_direct
  print *, ''

  !--------------------------------------------------------------------------
  ! FMM COMPUTATION
  !--------------------------------------------------------------------------
  print *, 'Computing FMM solution...'
  print *, ''

  call fmm_compute_forces(sys)

  print *, ''
  print *, 'FMM computation complete'

  ! Compute FMM energy
  energy_fmm = fmm_compute_energy(sys)
  print *, '  Energy (FMM):', energy_fmm
  print *, ''

  !--------------------------------------------------------------------------
  ! ERROR ANALYSIS
  !--------------------------------------------------------------------------
  print *, '============================================================'
  print *, 'ERROR ANALYSIS'
  print *, '============================================================'
  print *, ''

  ! Energy error
  energy_error = abs(energy_fmm - energy_direct) / abs(energy_direct)
  print *, 'Energy:'
  print *, '  Direct:', energy_direct
  print *, '  FMM:   ', energy_fmm
  print *, '  Relative error:', energy_error
  print *, ''

  ! Force errors
  max_force_error = 0.0_dp
  avg_force_error = 0.0_dp

  print *, 'Force errors (first 10 particles):'
  print *, '   i       Fx_direct       Fx_FMM        Fy_direct       Fy_FMM      Rel.Error'

  do i = 1, min(10, n_particles)
    fx_fmm = sys%particles(i)%fx
    fy_fmm = sys%particles(i)%fy
    fx_direct = fx_ref(i)
    fy_direct = fy_ref(i)

    ! Relative error in force magnitude
    r = sqrt((fx_fmm - fx_direct)**2 + (fy_fmm - fy_direct)**2) / &
        (sqrt(fx_direct**2 + fy_direct**2) + 1.0e-30_dp)

    write(*, '(I4, 4F16.8, ES14.4)') i, fx_direct, fx_fmm, fy_direct, fy_fmm, r

    max_force_error = max(max_force_error, r)
    avg_force_error = avg_force_error + r
  end do

  ! Compute full statistics
  do i = 1, n_particles
    fx_fmm = sys%particles(i)%fx
    fy_fmm = sys%particles(i)%fy
    fx_direct = fx_ref(i)
    fy_direct = fy_ref(i)

    r = sqrt((fx_fmm - fx_direct)**2 + (fy_fmm - fy_direct)**2) / &
        (sqrt(fx_direct**2 + fy_direct**2) + 1.0e-30_dp)

    if (i > 10) then  ! Continue accumulating for particles > 10
      max_force_error = max(max_force_error, r)
      avg_force_error = avg_force_error + r
    end if
  end do

  avg_force_error = avg_force_error / real(n_particles, dp)

  print *, ''
  print *, 'Force statistics:'
  print *, '  Max relative error:', max_force_error
  print *, '  Avg relative error:', avg_force_error
  print *, ''

  !--------------------------------------------------------------------------
  ! VERDICT
  !--------------------------------------------------------------------------
  print *, '============================================================'
  print *, 'VERDICT'
  print *, '============================================================'
  print *, ''

  if (energy_error < 1.0e-12_dp .and. max_force_error < 1.0e-12_dp) then
    print *, '✓ SUCCESS: Achieved 10^-12 precision!'
    print *, '  Energy error: ', energy_error, ' < 10^-12'
    print *, '  Force error:  ', max_force_error, ' < 10^-12'
  else if (energy_error < 1.0e-10_dp .and. max_force_error < 1.0e-10_dp) then
    print *, '⚠ PARTIAL: Achieved 10^-10 precision (close to target)'
    print *, '  Energy error: ', energy_error
    print *, '  Force error:  ', max_force_error
    print *, ''
    print *, 'Recommendations:'
    print *, '  - Increase p_max beyond', p_max
    print *, '  - Reduce theta (currently', sys%theta, ')'
    print *, '  - Use quad precision (real(16)) instead of double'
  else if (energy_error < 1.0e-6_dp .and. max_force_error < 1.0e-6_dp) then
    print *, '⚠ MODERATE: Achieved ~10^-6 precision (needs improvement)'
    print *, '  Energy error: ', energy_error
    print *, '  Force error:  ', max_force_error
    print *, ''
    print *, 'Issues detected:'
    print *, '  - Current implementation may have bugs in M2L/L2L operators'
    print *, '  - Spherical harmonic normalization may be incorrect'
    print *, '  - Translation coefficients need verification'
  else
    print *, '✗ FAILED: Accuracy is poor (>', max_force_error, ')'
    print *, '  Energy error: ', energy_error
    print *, '  Force error:  ', max_force_error
    print *, ''
    print *, 'Critical issues:'
    print *, '  - Fundamental problem with FMM implementation'
    print *, '  - Need to verify all mathematical formulas'
    print *, '  - Check octree construction and list building'
  end if

  print *, ''
  print *, 'NOTE: Current implementation is incomplete (operators commented out)'
  print *, 'This is a framework demonstration. Full implementation requires:'
  print *, '  1. Complete M2M, M2L, L2L, L2P operators'
  print *, '  2. Proper rotation-based FMM or translation operators'
  print *, '  3. Careful numerical stability analysis'
  print *, '  4. Extensive testing and validation'
  print *, ''

  ! Cleanup
  deallocate(fx_ref, fy_ref, phi_ref, seed)
  call fmm_destroy(sys)

  print *, 'Test complete'

end program test_fmm_precision
