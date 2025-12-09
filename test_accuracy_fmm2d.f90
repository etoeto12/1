program test_accuracy_fmm2d
  use fmm2d_accurate
  implicit none

  integer, parameter :: dp = selected_real_kind(15, 307)
  type(fmm_system) :: sys
  integer :: i, n, p
  real(dp) :: x, y, q
  real(dp) :: fx_ref(1000), fy_ref(1000), phi_ref(1000)
  real(dp) :: err_fx, err_fy, err_phi, max_err_f, avg_err_f
  real(dp) :: energy_ref, energy_fmm
  integer :: seed_size
  integer, allocatable :: seed(:)

  print *, '============================================================'
  print *, 'Accurate 2D FMM Test for 1/r potential (target: 10^-12)'
  print *, '============================================================'
  print *, ''

  ! Random seed
  call random_seed(size=seed_size)
  allocate(seed(seed_size))
  seed = 12345
  call random_seed(put=seed)

  ! Test different parameters
  do n = 50, 100, 50
    do p = 10, 20, 5
      print *, '------------------------------------------------------------'
      write(*, '(A,I4,A,I3)') ' Test: N =', n, ', p_max =', p
      print *, '------------------------------------------------------------'

      ! Initialize
      call fmm_init(sys, -1.0_dp, 1.0_dp, -1.0_dp, 1.0_dp, p_max=p, max_part=n)

      ! Generate particles
      do i = 1, n
        call random_number(x)
        call random_number(y)
        x = 1.8_dp * x - 0.9_dp
        y = 1.8_dp * y - 0.9_dp
        q = 1.0_dp
        call fmm_add_particle(sys, x, y, q)
      end do

      ! Build tree
      call fmm_build_tree(sys)

      ! Reference (direct)
      call compute_direct(sys)
      do i = 1, n
        fx_ref(i) = sys%particles(i)%fx
        fy_ref(i) = sys%particles(i)%fy
        phi_ref(i) = sys%particles(i)%phi
      end do

      energy_ref = 0.0_dp
      do i = 1, n
        energy_ref = energy_ref + 0.5_dp * sys%particles(i)%q * phi_ref(i)
      end do

      ! FMM
      call fmm_compute(sys)

      energy_fmm = 0.0_dp
      do i = 1, n
        energy_fmm = energy_fmm + 0.5_dp * sys%particles(i)%q * sys%particles(i)%phi
      end do

      ! Compute errors
      max_err_f = 0.0_dp
      avg_err_f = 0.0_dp

      do i = 1, n
        err_fx = abs(sys%particles(i)%fx - fx_ref(i))
        err_fy = abs(sys%particles(i)%fy - fy_ref(i))
        err_phi = abs(sys%particles(i)%phi - phi_ref(i))

        ! Relative force error
        err_fx = err_fx / (abs(fx_ref(i)) + 1.0e-10_dp)
        err_fy = err_fy / (abs(fy_ref(i)) + 1.0e-10_dp)

        max_err_f = max(max_err_f, max(err_fx, err_fy))
        avg_err_f = avg_err_f + 0.5_dp * (err_fx + err_fy)
      end do

      avg_err_f = avg_err_f / real(n, dp)

      ! Results
      write(*, '(A,ES12.4)') '  Energy (direct): ', energy_ref
      write(*, '(A,ES12.4)') '  Energy (FMM):    ', energy_fmm
      write(*, '(A,ES12.4)') '  Energy error:    ', abs(energy_fmm - energy_ref) / abs(energy_ref)
      write(*, '(A,ES12.4)') '  Max force error: ', max_err_f
      write(*, '(A,ES12.4)') '  Avg force error: ', avg_err_f
      print *, ''

      ! Verdict
      if (max_err_f < 1.0e-12_dp) then
        print *, '  ✓ SUCCESS: Achieved 10^-12 precision!'
      else if (max_err_f < 1.0e-10_dp) then
        print *, '  ⚠ GOOD: Achieved 10^-10 precision (close to target)'
      else if (max_err_f < 1.0e-6_dp) then
        print *, '  ⚠ MODERATE: ~10^-6 precision (needs improvement)'
      else
        print *, '  ✗ FAILED: Poor accuracy'
      end if

      print *, ''

      ! Cleanup
      call fmm_destroy(sys)
    end do
  end do

  deallocate(seed)

  print *, '============================================================'
  print *, 'Testing complete'
  print *, '============================================================'

end program test_accuracy_fmm2d
