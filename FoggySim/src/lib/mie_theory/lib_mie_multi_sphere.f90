module lib_mie_multi_sphere
    use libmath
    use lib_field
    use lib_mie_type
    use lib_mie_type_functions
    use lib_mie_ms_helper_functions
    use lib_mie_ms_solver_gmres
    use lib_mie_single_sphere
    implicit none

    private

    public :: lib_mie_multi_sphere_test_functions

    contains

!        ! Arguments
!        ! ----
!        !   k: type(spherical_coordinate_real_type)
!        !       ?
!        !   lambda: double precision
!        !       vacuum wave length
!        !   n_medium: double precision
!        !       refractive index of the medium
!        !   sphere: type(sphere_type), dimension(:)
!        !       list of spheres
!        !   sphere_parameter: type(sphere_parameter_type), dimension(:)
!        !       list of shared sphere parameters
!        !   z_selector: integer
!        !       parameter of the spherical harmonics
!        !       values:
!        !           1: spherical Bessel function first kind   j_n
!        !           2: spherical Bessel function second kind  y_n
!        !           3: spherical Hankel function first kind   h^(1)_n
!        !           4: spherical Hankel function second kind  h^(2)_n
!        !   a_old, b_old: type(list_list_cmplx), optional
!        !       coefficient of the previous calculation, if f < 1
!        !   f: double precision, optional (std: 1)
!        !       numerical factor (0, 1]
!        !       "In our actual calculations, some multi-
!        !        sphere systems do not converge if f 5 1, but they do
!        !        converge when the value of f is reduced to, say, 0.7."[1]
!        ! Returns
!        ! ----
!        !
!        !
!        !
!        !
!        ! Reference: [1] Electromagnetic scattering by an aggregate of spheres, Yu-lin Xu, eq. 30
!        subroutine lib_mie_ms_get_interactive_scattering_coefficients(k, lambda, n_medium, &
!                                                                   sphere, sphere_parameter, sphere_j, &
!                                                                   z_selector, &
!                                                                   a, b)
!            implicit none
!            ! dummy
!            type(spherical_coordinate_real_type), intent(in) :: k
!            double precision, intent(in) :: lambda
!            double precision, intent(in) :: n_medium
!            type(lib_mie_sphere_type), dimension(:), intent(in) :: sphere
!            type(lib_mie_sphere_parameter_type), dimension(:), intent(in) :: sphere_parameter
!            integer :: sphere_j
!            integer(kind=1) :: z_selector
!
!            type(list_list_cmplx), intent(inout) :: a
!            type(list_list_cmplx), intent(inout) :: b
!
!            ! auxiliary
!            integer :: mu
!            integer :: nu
!
!
!        end subroutine lib_mie_ms_get_interactive_scattering_coefficients

        subroutine lib_mie_ms_calculate_scattering_coefficients_ab_nm(simulation)
            use file_io
            implicit none
            ! dummy
            type(lib_mie_simulation_parameter_type), intent(inout) :: simulation

            ! auxiliary
            double complex, dimension(:,:), allocatable :: matrix_a
            double complex, dimension(:), allocatable :: vector_x
            double complex, dimension(:), allocatable :: vector_b

            ! test
            logical :: test

            ! initial values
            call lib_mie_ss_calculate_scattering_coefficients_ab_nm(simulation)

!            ! test
!            simulation%sphere_parameter_list(1)%n_range = (/ 1, 3 /)
!            simulation%sphere_parameter_list(2)%n_range = (/ 1, 3 /)
!            ! ~~~ test ~~~

            call lib_mie_ms_solver_gmres_run_without_ml_fmm(simulation)

        end subroutine lib_mie_ms_calculate_scattering_coefficients_ab_nm

        ! Argument
        ! ----
        !   simulation: type(lib_mie_simulation_parameter_type)
        !       simulation data set
        !   x_0: cartesian_coordinate_real_type
        !
        function lib_mie_ms_get_field(simulation, x_0) result(field)
            implicit none
            ! dummy
            type(lib_mie_simulation_parameter_type), intent(in) :: simulation
            type(cartesian_coordinate_real_type), intent(in) :: x_0

            type(cartesian_coordinate_cmplx_type), dimension(2) :: field

            ! auxiliaray
            integer :: sphere_no
            integer :: parameter_no

            double precision :: e_field_0
            double precision :: k_0
            double precision :: n_medium

            type(cartesian_coordinate_real_type) :: x_j

            type(lib_mie_sphere_type) :: sphere
            type(lib_mie_sphere_parameter_type) :: sphere_parameter

            type(spherical_coordinate_cmplx_type), dimension(:, :), allocatable :: buffer_field

            logical :: inside_sphere

            inside_sphere = .true.

            e_field_0 = simulation%illumination%e_field_0
            k_0 = 2 * PI / simulation%illumination%lambda_0
            n_medium = simulation%refractive_index_medium

            allocate(buffer_field(lbound(simulation%sphere_list, 1):ubound(simulation%sphere_list, 1), 2))

            field(:)%x = dcmplx(0,0)
            field(:)%y = dcmplx(0,0)
            field(:)%z = dcmplx(0,0)

            do sphere_no = lbound(simulation%sphere_list, 1), ubound(simulation%sphere_list, 1)
                if (inside_sphere) then
                    sphere = simulation%sphere_list(sphere_no)
                    parameter_no = sphere%sphere_parameter_index
                    sphere_parameter = simulation%sphere_parameter_list(parameter_no)

                    x_j = sphere%d_0_j - x_0

                    if (abs(x_j) .gt. sphere_parameter%radius) then
                        buffer_field(sphere_no, 1:2) = lib_mie_ss_get_field(x_0, &
                                                                         e_field_0, k_0, n_medium, &
                                                                         sphere, sphere_parameter)

                    else
                        inside_sphere = .false. ! OMP critical ?!
                        exit
                    end if
                end if
            end do

            do sphere_no = lbound(simulation%sphere_list, 1), ubound(simulation%sphere_list, 1)
                field(1) = field(1) + make_cartesian(buffer_field(sphere_no, 1), x_j)
                field(2) = field(2) + make_cartesian(buffer_field(sphere_no, 2), x_j)
            end do
        end function

        function lib_mie_multi_sphere_test_functions() result(rv)
            implicit none
            ! dummy
            integer :: rv

            ! auxiliaray
            double precision, parameter :: ground_truth_e = 10.0_8**(-6.0_8)
            ! CPU-time
            real :: test_start, test_finish
            ! WALL-time
            INTEGER :: test_count_start, test_count_finish, test_count_rate

            rv = 0

            call system_clock(test_count_start, test_count_rate)
            call cpu_time(test_start)

            if (.not. test_lib_mie_ms_get_field()) then
                rv = rv + 1
            end if

            call cpu_time(test_finish)
            call system_clock(test_count_finish, test_count_rate)

            print *, ""
            print *, "------lib_mie_multi_sphere_test_functions------"
            print '("  CPU-Time = ",f10.3," seconds.")',test_finish-test_start
            print '("  WALL-Time = ",f10.3," seconds.")',(test_count_finish-test_count_start) / real(test_count_rate)
            print *, ""
            if (rv == 0) then
                print *, "lib_mie_multi_sphere_test_functions tests: OK"
            else
                print *, rv,"lib_mie_multi_sphere_test_functions test(s) FAILED"
            end if
            print *, "------------------------------------------------------------"
            print *, ""

            contains

            function test_lib_mie_ms_get_field() result(rv)
                implicit none
                ! dummy
                logical :: rv

                ! parameter
                integer, parameter :: number_of_waves = 1
                ! auxiliaray
                type(lib_mie_simulation_parameter_type) :: simulation

                ! illumination parameter
                double precision :: lambda_0
                double precision :: k_0
                double precision :: e_field_0

                double precision, dimension(number_of_waves) :: plane_wave_g
                type(cartesian_coordinate_real_type), dimension(number_of_waves) :: plane_wave_k
                type(cartesian_coordinate_real_type), dimension(number_of_waves) :: plane_wave_d_0_i

                double precision :: n_medium

                type(cartesian_coordinate_real_type) :: sphere_d_0_j
                double precision :: r_particle
                double complex :: n_particle
                integer, dimension(2) :: n_range

                type(cartesian_coordinate_real_type) :: buffer_car

                type(cartesian_coordinate_cmplx_type), dimension(2) :: field

                integer :: i
                integer :: ii
                double precision :: x
                double precision :: y
                double precision :: z
                type(cartesian_coordinate_real_type) :: x_0
                type(cartesian_coordinate_cmplx_type), dimension(:, :), allocatable :: e_field_s
                type(cartesian_coordinate_cmplx_type), dimension(:, :), allocatable :: h_field_s
                double precision, dimension(2) :: x_range
                double precision, dimension(2) :: z_range
                real(kind=8) :: step_size

                integer :: no_x_values
                integer :: no_z_values

                simulation%spherical_harmonics%z_selector_incident_wave = 1
                simulation%spherical_harmonics%z_selector_scatterd_wave = 3
                simulation%spherical_harmonics%z_selector_translation = 1

                ! set illumination parameter
                e_field_0 = 1
                lambda_0 = 1 * unit_mu

                n_medium = 1

                k_0 = 2 * PI / lambda_0

                simulation%illumination%lambda_0 = lambda_0

                plane_wave_g(:) = 1

                buffer_car%x = 0
                buffer_car%y = 0
                buffer_car%z = 0
                plane_wave_d_0_i(:) = buffer_car

                buffer_car%x = 0
                buffer_car%y = 0
                buffer_car%z = 1
                buffer_car = buffer_car / abs(buffer_car) / lambda_0
                plane_wave_k(1) = buffer_car

!                buffer_car%x = -1
!                buffer_car%y = 0
!                buffer_car%z = 1
!                buffer_car = buffer_car / abs(buffer_car) / lambda
!                plane_wave_k(2) = buffer_car

                simulation%illumination = lib_mie_type_func_get_plane_wave_illumination(lambda_0, e_field_0, &
                                                                                        plane_wave_g, &
                                                                                        plane_wave_k, &
                                                                                        plane_wave_d_0_i)

                simulation%refractive_index_medium = n_medium

                ! set spheres
                allocate(simulation%sphere_list(2))
                simulation%sphere_list(1)%sphere_parameter_index = 1
                simulation%sphere_list(2)%sphere_parameter_index = 1

                sphere_d_0_j%x = -2 * unit_mu
                sphere_d_0_j%y = 0
                sphere_d_0_j%z = 0
                simulation%sphere_list(1)%d_0_j = sphere_d_0_j

                sphere_d_0_j%x = 2 * unit_mu
                sphere_d_0_j%y = 0
                sphere_d_0_j%z = 1
                simulation%sphere_list(2)%d_0_j = sphere_d_0_j

                ! set sphere parameter
                allocate(simulation%sphere_parameter_list(2))

                ! set 1
                r_particle = 1 * unit_mu
                n_particle = dcmplx(1.33_8, 0)

                n_range = lib_mie_ss_test_convergence_plane_wave(lambda_0, n_medium, r_particle, n_particle)
                if (n_range(1) .gt. 0) then
                    n_range(2) = n_range(1)
                    n_range(1) = 1
                else
                    print *, "test_lib_mie_ms_get_field: ERROR"
                    print *, "  lib_mie_ss_test_convergence_plane_wave"
                    print *, "  rv(1): ", n_range(1)
                    print *, "  rv(2): ", n_range(2)
                end if

                simulation%spherical_harmonics%n_range = n_range

                simulation%sphere_parameter_list(1) = lib_mie_type_func_get_sphere_parameter(lambda_0, n_medium, &
                                                                                          r_particle, n_particle, &
                                                                                          n_range)

                ! set 2
                r_particle = 0.5 * unit_mu
                n_particle = dcmplx(1.33_8, 0)

                n_range = lib_mie_ss_test_convergence_plane_wave(lambda_0, n_medium, r_particle, n_particle)
                if (n_range(1) .gt. 0) then
                    n_range(2) = n_range(1)
                    n_range(1) = 1
                else
                    print *, "test_lib_mie_ms_get_field: ERROR"
                    print *, "  lib_mie_ss_test_convergence_plane_wave"
                    print *, "  rv(1): ", n_range(1)
                    print *, "  rv(2): ", n_range(2)
                end if

                simulation%sphere_parameter_list(2) = lib_mie_type_func_get_sphere_parameter(lambda_0, n_medium, &
                                                                                          r_particle, n_particle, &
                                                                                          n_range)

                call lib_mie_ss_constructor(simulation)


                call lib_mie_ms_calculate_scattering_coefficients_ab_nm(simulation)

                ! evaluate and export
                x_range = (/ -5.0_8 * unit_mu, 5.0_8 * unit_mu /)
                z_range = (/ -5_8 * unit_mu, 10.0_8 * unit_mu /)
                step_size = 0.05_8 * unit_mu

                no_x_values = abs(int(floor((x_range(2)-x_range(1))/step_size)))
                no_z_values = abs(int(floor((z_range(2)-z_range(1))/step_size)))

                allocate(e_field_s(no_x_values, no_z_values))
                allocate(h_field_s(no_x_values, no_z_values))

                x = 0
                y = 0
                z = 0

                do i=1, no_x_values
                    x = x_range(1) + (i-1) * step_size
                    do ii=1, no_z_values
                        z = z_range(1) + (ii-1) * step_size

                        x_0%x = x
                        x_0%y = y
                        x_0%z = z

                        field = lib_mie_ms_get_field(simulation, x_0)
                        e_field_s(i,ii) = field(1)
                        h_field_s(i,ii) = field(2)
                    end do
                 end do

                rv = lib_field_export(e_field_s, h_field_s, "temp/real/")


            end function
        end function
end module lib_mie_multi_sphere
