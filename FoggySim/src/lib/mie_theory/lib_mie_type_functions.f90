module lib_mie_type_functions
    use libmath
    use lib_mie_type
    use lib_mie_ss_helper_functions
    use lib_mie_vector_spherical_harmonics
    implicit none

    contains

        ! Argument
        ! ----
        !   lambda_0: double precision
        !       vacuum wave length
        !   n_medium: double precision
        !       refractive index of the medium
        !   r_particle: double precision
        !       radius of the particel
        !   n_particle: double complex
        !       refractive index of the particle
        !   n_range: integer, dimension(2)
        !       first and last index (degree) of the sum to calculate the electical field
        !       e.g. n = (/ 1, 45 /)
        !
        ! Returns
        ! ----
        !   sphere_parameter: type(lib_mie_sshere_parameter_type)
        function lib_mie_type_func_get_sphere_parameter(lambda_0, n_medium, &
                                              r_particle, n_particle,&
                                              n_range) &
                                            result (sphere_parameter)
            implicit none
            ! dummy
            double precision, intent(in) :: lambda_0
            double precision, intent(in) :: n_medium
            double precision, intent(in) :: r_particle
            double complex, intent(in) :: n_particle
            integer, dimension(2) :: n_range

            type(lib_mie_sphere_parameter_type) :: sphere_parameter

            ! auxiliary
            double precision :: size_parameter

            size_parameter = 2 * PI * n_medium * r_particle / lambda_0

            sphere_parameter%n_range = n_range
            sphere_parameter%size_parameter = size_parameter
            sphere_parameter%radius = r_particle
            sphere_parameter%refractive_index = n_particle

            allocate(sphere_parameter%a_n%item(n_range(1):n_range(2)))
            allocate(sphere_parameter%b_n%item(n_range(1):n_range(2)))

            if (aimag(n_particle) .eq. 0d0) then
                call lib_mie_ss_hf_get_coefficients_a_n_b_n(size_parameter, real(n_particle)/n_medium, n_range, &
                                                       sphere_parameter%a_n%item, sphere_parameter%b_n%item)
            else
                call lib_mie_ss_hf_get_coefficients_a_n_b_n(size_parameter, n_particle/n_medium, n_range, &
                                                           sphere_parameter%a_n%item, sphere_parameter%b_n%item)
            end if
        end function lib_mie_type_func_get_sphere_parameter

        ! Skatch
        ! ----
        !             _________
        !             ___k^____
        !             ____|____
        !             _________ plane wave
        !                 z
        !                 ^
        !             K_i |
        !                 --> x
        !                ^
        !               /
        !           z  /d_0_i
        !           ^ /
        !       K_0 |/
        !           --> x
        !
        ! K_0: world coordinate system
        ! K_i: illumination coordinate system
        !
        !
        ! Argument
        ! ----
        !   type: integer
        !       illumination type
        !       - 1: plane wave
        !   lambda_0: double precision
        !       vacuum wave length
        !   e_field_0: double precision
        !       magnitude of the oscilating electrical field (peak value)
        !   k: type(cartesian_coordinate_real_type)
        !       sets only the direction of the wave vector, the length of this vector is calculated internally as follows:
        !       |k| = 2 Pi / lambda
        !   d_0_i: type(cartesian_coordinate_real_type)
        !       position of the illumination coordinate system respect to the world coordinate system
        !
        ! Returns
        ! ----
        !   illumination: lib_mie_illumination_parameter
        !
        function lib_mie_type_func_get_plane_wave_illumination(lambda_0, e_field_0, &
                                                               g, k, d_0_i) result(illumination)
            implicit none
            ! dummy
            double precision, intent(in) :: lambda_0
            double precision, intent(in) :: e_field_0
            double precision, dimension(:), intent(in) :: g
            type(cartesian_coordinate_real_type), dimension(lbound(g, 1):ubound(g, 1)), intent(in) :: k
            type(cartesian_coordinate_real_type), dimension(lbound(g, 1):ubound(g, 1)), intent(in) :: d_0_i

            type(lib_mie_illumination_parameter) :: illumination

            ! auxiliary
            integer :: i
            double precision :: abs_k

            illumination%type = 1
            illumination%lambda_0 = lambda_0
            illumination%e_field_0 = e_field_0

            allocate(illumination%plane_wave(lbound(g, 1):ubound(g, 1)))

            do i = lbound(g, 1), ubound(g, 1)
                abs_k = abs(k(i))

                if (abs_k .gt. 0) then
                    illumination%plane_wave(i)%g = g(i)
                    illumination%plane_wave(i)%d_0_i = d_0_i(i)
                    illumination%plane_wave(i)%wave_vector_0 = k(i) * 2D0 * PI / (abs_k * lambda_0)
                else
                    print *, "lib_mie_type_func_get_illumination: ERROR"
                    print *, "  abs(k(i=", i, ")) = 0"
                end if
            end do

        end function lib_mie_type_func_get_plane_wave_illumination



        function lib_mie_type_functions_test_functions() result(rv)
            implicit none
            ! dummy
            integer :: rv

            ! auxiliary

            contains
                function test_() result(rv)
                    implicit none
                    ! dummy
                    logical :: rv

                end function
        end function lib_mie_type_functions_test_functions

end module lib_mie_type_functions
