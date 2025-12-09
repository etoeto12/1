program test_grid_accuracy
  use multipole_accurate
  implicit none

  integer, parameter :: dp = selected_real_kind(15, 307)
  type(multipole_system) :: sys
  integer :: i, n, grid_size, p
  real(dp) :: x, y, q
  real(dp), allocatable :: fx_ref(:), fy_ref(:)
  real(dp) :: err_fx, err_fy, max_err, avg_err
  real(dp) :: energy_ref, energy_fmm, energy_err
  integer :: seed_size
  integer, allocatable :: seed(:)

  print *, '============================================================'
  print *, 'Uniform Grid FMM Test for 1/r potential'
  print *, 'Target: 10^-12 accuracy'
  print *, '============================================================'
  print *, ''

  call random_seed(size=seed_size)
  allocate(seed(seed_size))
  seed = 12345
  call random_seed(put=seed)

  ! Test different configurations
  do n = 50, 200, 50
    do grid_size = 5, 10, 5
      do p = 10, 15, 5
        print *, '------------------------------------------------------------'
        write(*, '(A,I4,A,I3,A,I3)') ' N =', n, ', grid =', grid_size, 'x', grid_size, ', p_max =', p
        print *, '------------------------------------------------------------'

        ! Initialize
        call init_system(sys, 1.0_dp, 1.0_dp, grid_size, grid_size, p_max=p, max_part=n)

        ! Generate particles
        do i = 1, n
          call random_number(x)
          call random_number(y)
          x = 1.8_dp * x - 0.9_dp
          y = 1.8_dp * y - 0.9_dp
          q = 1.0_dp
          call add_particle(sys, x, y, q)
        end do

        ! Allocate reference arrays
        allocate(fx_ref(n), fy_ref(n))

        ! Direct computation
        call compute_all_forces(sys)
        do i = 1, n
          fx_ref(i) = sys%particles(i)%fx
          fy_ref(i) = sys%particles(i)%fy
        end do
        energy_ref = compute_energy_direct(sys)

        ! FMM computation
        call compute_all_forces(sys)
        energy_fmm = compute_energy_multipole(sys)

        ! Compute errors
        max_err = 0.0_dp
        avg_err = 0.0_dp

        do i = 1, n
          err_fx = abs(sys%particles(i)%fx - fx_ref(i))
          err_fy = abs(sys%particles(i)%fy - fy_ref(i))

          err_fx = err_fx / (abs(fx_ref(i)) + 1.0e-10_dp)
          err_fy = err_fy / (abs(fy_ref(i)) + 1.0e-10_dp)

          max_err = max(max_err, max(err_fx, err_fy))
          avg_err = avg_err + 0.5_dp * (err_fx + err_fy)
        end do

        avg_err = avg_err / real(n, dp)
        energy_err = abs(energy_fmm - energy_ref) / abs(energy_ref)

        write(*, '(A,ES12.4)') '  Energy (direct): ', energy_ref
        write(*, '(A,ES12.4)') '  Energy (FMM):    ', energy_fmm
        write(*, '(A,ES12.4)') '  Energy error:    ', energy_err
        write(*, '(A,ES12.4)') '  Max force error: ', max_err
        write(*, '(A,ES12.4)') '  Avg force error: ', avg_err

        if (max_err < 1.0e-12_dp) then
          print *, '  ✓ SUCCESS: 10^-12 achieved!'
        else if (max_err < 1.0e-10_dp) then
          print *, '  ⚠ GOOD: 10^-10 (close)'
        else if (max_err < 1.0e-6_dp) then
          print *, '  ⚠ MODERATE: ~10^-6'
        else
          print *, '  ✗ FAILED: Poor accuracy'
        end if

        print *, ''

        deallocate(fx_ref, fy_ref)
        call destroy_system(sys)
      end do
    end do
  end do

  deallocate(seed)

  print *, '============================================================'
  print *, 'Testing complete'
  print *, '============================================================'

end program test_grid_accuracy
