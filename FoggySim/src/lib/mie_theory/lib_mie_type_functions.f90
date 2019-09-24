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
        !       sets only the direction of the wave vector, the length of this vector is calculated as follows:
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

        ! Formats problem of  multi sphere scattering to be able to use a solver.
        !
        ! Formula: A x = b
        !   x: scattering coefficients
        !   b: illumination coefficients
        !   A: mixture of Mie coefficients and vector spherical translation coefficients
        !
        !   -> apply eq. 23 [1] to eq. 30 [2]
        !
        ! Argument
        ! ----
        !   simulation_parameter: type(lib_mie_simulation_parameter_type)
        !
        ! Returns
        ! ----
        !   vector_b: double complex, dimension(:)
        !
        ! Reference: [1] Computation of scattering from clusters of spheres using the fast multipole method, Nail A. Gumerov, and Ramani Duraiswami
        !            [2] Electromagnetic scatteringby an aggregate of spheres, Yu-lin Xu
        subroutine lib_mie_type_func_solver_get_vector_b(simulation_parameter, vector_b)
            implicit none
            ! dummy
            type(lib_mie_simulation_parameter_type), intent(in) :: simulation_parameter

            double complex, dimension(:), allocatable, intent(inout) :: vector_b

            ! auxiliary
            integer :: i
            integer :: ii
            integer :: n
            integer :: m

            integer :: first_sphere
            integer :: last_sphere

            integer :: sphere_parameter_no

            integer, dimension(:), allocatable :: counter
            integer :: counter_sum

            type(cartesian_coordinate_real_type) :: d_0_j
            integer(kind=4), dimension(2) :: n_range

            type(list_list_cmplx), dimension(:), allocatable :: p_nm
            type(list_list_cmplx), dimension(:), allocatable :: q_nm

            first_sphere = lbound(simulation_parameter%sphere_list, 1)
            last_sphere = ubound(simulation_parameter%sphere_list, 1)

            allocate(p_nm(first_sphere:last_sphere))
            allocate(q_nm(first_sphere:last_sphere))

            allocate(counter(first_sphere:last_sphere))

            ! create vector_b
            counter = 0
            !$OMP PARALLEL DO PRIVATE(i, d_0_j, sphere_parameter_no, n_range)
            do i = first_sphere, last_sphere

                d_0_j = simulation_parameter%sphere_list(i)%d_0_j
                sphere_parameter_no = simulation_parameter%sphere_list(i)%sphere_parameter_index
                n_range = simulation_parameter%sphere_parameter_list(sphere_parameter_no)%n_range

                call lib_mie_ss_hf_get_p_q_j_j(simulation_parameter%illumination, &
                                               simulation_parameter%refractive_index_medium, &
                                               d_0_j, n_range, p_nm(i), q_nm(i))
                counter(i) = (1 + n_range(2))**2 - n_range(1)**2
            end do
            !$OMP END PARALLEL DO

            counter_sum = 2 * sum(counter)

            if (allocated(vector_b)) then
                deallocate(vector_b)
            end if
            allocate(vector_b(counter_sum))

            vector_b = 0

            ii = 1
            do i = first_sphere, last_sphere
                sphere_parameter_no = simulation_parameter%sphere_list(i)%sphere_parameter_index
                n_range = simulation_parameter%sphere_parameter_list(sphere_parameter_no)%n_range

                do n = n_range(1), n_range(2)
                    do m = -n, n
                        vector_b(ii) = p_nm(i)%item(n)%item(m)
                        vector_b(ii+1) = q_nm(i)%item(n)%item(m)
                        ii = ii + 2
                    end do
                end do
            end do

        end subroutine lib_mie_type_func_solver_get_vector_b

        ! Formats problem of  multi sphere scattering to be able to use a solver.
        !
        ! Formula: A x = b
        !   x: scattering coefficients
        !   b: illumination coefficients
        !   A: mixture of Mie coefficients and vector spherical translation coefficients
        !
        !   -> apply eq. 23 [1] to eq. 30 [2]
        !
        ! Argument
        ! ----
        !   simulation_parameter: type(lib_mie_simulation_parameter_type)
        !
        ! Returns
        ! ----
        !   vector_b: double complex, dimension(:)
        !
        ! Reference: [1] Computation of scattering from clusters of spheres using the fast multipole method, Nail A. Gumerov, and Ramani Duraiswami
        !            [2] Electromagnetic scatteringby an aggregate of spheres, Yu-lin Xu
        subroutine lib_mie_type_func_solver_get_vector_x(simulation_parameter, vector_x)
            implicit none
            ! dummy
            type(lib_mie_simulation_parameter_type), intent(in) :: simulation_parameter

            double complex, dimension(:), allocatable, intent(inout) :: vector_x

            ! auxiliary
            integer :: i
            integer :: ii
            integer :: n
            integer :: m

            integer :: first_sphere
            integer :: last_sphere

            integer :: sphere_parameter_no

            integer, dimension(:), allocatable :: counter
            integer :: counter_sum

            integer(kind=4), dimension(2) :: n_range

            first_sphere = lbound(simulation_parameter%sphere_list, 1)
            last_sphere = ubound(simulation_parameter%sphere_list, 1)

            allocate(counter(first_sphere:last_sphere))

            ! create vector_x
            counter = 0
            do i = first_sphere, last_sphere
                sphere_parameter_no = simulation_parameter%sphere_list(i)%sphere_parameter_index
                n_range = simulation_parameter%sphere_parameter_list(sphere_parameter_no)%n_range

                counter(i) = (1 + n_range(2))**2 - n_range(1)**2
            end do

            counter_sum = 2 * sum(counter)

            if (allocated(vector_x)) then
                deallocate(vector_x)
            end if
            allocate(vector_x(counter_sum))

            vector_x = 0

            ii = 1
            do i = first_sphere, last_sphere
                sphere_parameter_no = simulation_parameter%sphere_list(i)%sphere_parameter_index
                n_range = simulation_parameter%sphere_parameter_list(sphere_parameter_no)%n_range

                do n = n_range(1), n_range(2)
                    do m = -n, n
                        vector_x(ii) = simulation_parameter%sphere_list(i)%a_nm%item(n)%item(m)
                        vector_x(ii+1) = simulation_parameter%sphere_list(i)%b_nm%item(n)%item(m)
                        ii = ii + 2
                    end do
                end do
            end do

        end subroutine lib_mie_type_func_solver_get_vector_x

        ! Formats problem of multi sphere scattering to be able to use a solver.
        !
        ! Formula: A x = b
        !   x: scattering coefficients
        !   b: illumination coefficients
        !   A: mixture of Mie coefficients and vector spherical translation coefficients
        !
        !   -> apply eq. 23 [1] to eq. 30 [2]
        !
        ! Argument
        ! ----
        !   vector_x: double complex, dimension(:)
        !
        ! Returns
        ! ----
        !   simulation_parameter: type(lib_mie_simulation_parameter_type)
        !
        ! Reference: [1] Computation of scattering from clusters of spheres using the fast multipole method, Nail A. Gumerov, and Ramani Duraiswami
        !            [2] Electromagnetic scatteringby an aggregate of spheres, Yu-lin Xu
        subroutine lib_mie_type_func_solver_set_sphere_parameter_ab_nm(vector_x, simulation_parameter)
            implicit none
            ! dummy
            double complex, dimension(:), intent(in) :: vector_x

            type(lib_mie_simulation_parameter_type), intent(inout) :: simulation_parameter

            ! auxiliary
            integer :: i
            integer :: ii
            integer :: n
            integer :: m

            integer :: first_sphere
            integer :: last_sphere

            integer :: sphere_parameter_no

            integer(kind=4), dimension(2) :: n_range

            first_sphere = lbound(simulation_parameter%sphere_list, 1)
            last_sphere = ubound(simulation_parameter%sphere_list, 1)

            ii = 1
            do i = first_sphere, last_sphere
                sphere_parameter_no = simulation_parameter%sphere_list(i)%sphere_parameter_index
                n_range = simulation_parameter%sphere_parameter_list(sphere_parameter_no)%n_range

                do n = n_range(1), n_range(2)
                    do m = -n, n
                        simulation_parameter%sphere_list(i)%a_nm%item(n)%item(m) = vector_x(ii)
                        simulation_parameter%sphere_list(i)%b_nm%item(n)%item(m) = vector_x(ii+1)
                        ii = ii + 2
                    end do
                end do
            end do

        end subroutine lib_mie_type_func_solver_set_sphere_parameter_ab_nm

        ! Formats problem of  multi sphere scattering to be able to use a solver.
        !
        ! Formula: A x = b
        !   x: scattering coefficients
        !   b: illumination coefficients
        !   A: mixture of Mie coefficients and vector spherical translation coefficients
        !
        !   -> apply eq. 23 [1] to eq. 30 [2]
        !
        ! Argument
        ! ----
        !   simulation_parameter: type(lib_mie_simulation_parameter_type)
        !
        ! Returns
        ! ----
        !   vector_b: double complex, dimension(:)
        !
        ! Reference: [1] Computation of scattering from clusters of spheres using the fast multipole method, Nail A. Gumerov, and Ramani Duraiswami
        !            [2] Electromagnetic scatteringby an aggregate of spheres, Yu-lin Xu
        subroutine lib_mie_type_func_solver_get_matrix_a(simulation_parameter, matrix_a)
            implicit none
            ! dummy
            type(lib_mie_simulation_parameter_type), intent(in) :: simulation_parameter

            double complex, dimension(:,:), allocatable, intent(inout) :: matrix_a

            ! auxiliary
            integer :: i
            integer :: ii
            integer :: m_i
            integer :: m_ii
            integer :: n
            integer :: m
            integer :: nu
            integer :: mu

            integer :: first_sphere
            integer :: last_sphere

            integer :: sphere_parameter_no

            integer, dimension(:,:), allocatable :: counter
            integer :: counter_sum

            integer(kind=1) :: z_selector

            type(cartesian_coordinate_real_type) :: d_0_j
            type(cartesian_coordinate_real_type) :: d_0_l
            type(cartesian_coordinate_real_type) :: d_j_l
            integer(kind=4), dimension(2) :: n_range_j
            integer(kind=4), dimension(2) :: n_range_l

            type(list_cmplx), dimension(:), allocatable :: a_n
            type(list_cmplx), dimension(:), allocatable :: b_n
            type(list_4_cmplx), dimension(:,:), allocatable :: a_nmnumu
            type(list_4_cmplx), dimension(:,:), allocatable :: b_nmnumu

            first_sphere = lbound(simulation_parameter%sphere_list, 1)
            last_sphere = ubound(simulation_parameter%sphere_list, 1)

            allocate(a_n(first_sphere:last_sphere))
            allocate(b_n(first_sphere:last_sphere))

            allocate(a_nmnumu(first_sphere:last_sphere, first_sphere:last_sphere))
            allocate(b_nmnumu(first_sphere:last_sphere, first_sphere:last_sphere))

            allocate(counter(first_sphere:last_sphere, first_sphere:last_sphere))

            z_selector = simulation_parameter%spherical_harmonics%z_selector_translation

            ! calculate elements of matrix_a
            counter = 0
!            !$OMP PARALLEL DO PRIVATE(i, ii, d_0_j, d_0_l, d_j_l, sphere_parameter_no, n_range_l, n_range_j)
            do i = first_sphere, last_sphere
                do ii = first_sphere, last_sphere

                    if (i .eq. ii) then
                        sphere_parameter_no = simulation_parameter%sphere_list(i)%sphere_parameter_index
                        n_range_j = simulation_parameter%sphere_parameter_list(sphere_parameter_no)%n_range

                        a_n(i) = simulation_parameter%sphere_parameter_list(sphere_parameter_no)%a_n
                        b_n(i) = simulation_parameter%sphere_parameter_list(sphere_parameter_no)%b_n

                        counter(i, ii) = (1 + n_range_j(2))**2 - n_range_j(1)**2
!                        counter(i, ii) = n_range_j(2) - n_range_j(1) + 1
                    else
                        d_0_j = simulation_parameter%sphere_list(i)%d_0_j
                        sphere_parameter_no = simulation_parameter%sphere_list(i)%sphere_parameter_index
                        n_range_j = simulation_parameter%sphere_parameter_list(sphere_parameter_no)%n_range

                        d_0_l = simulation_parameter%sphere_list(ii)%d_0_j
                        sphere_parameter_no = simulation_parameter%sphere_list(ii)%sphere_parameter_index
                        n_range_l = simulation_parameter%sphere_parameter_list(sphere_parameter_no)%n_range

                        d_j_l = d_0_l - d_0_j
                        call lib_mie_vector_spherical_harmonics_translation_coefficient(d_j_l, &
                                n_range_j, n_range_l, z_selector, &
                                a_nmnumu(i, ii), b_nmnumu(i, ii))

                        counter(i, ii) = (1 + n_range_j(2))**2 - n_range_j(1)**2
!                        counter(i, ii) = ( (1 + n_range_j(2))**2 - n_range_j(1)**2 ) &
!                                     * ( (1 + n_range_l(2))**2 - n_range_l(1)**2 )
                    end if
                end do
            end do
!            !$OMP END PARALLEL DO

            counter_sum = 2 * sum(counter)

            if (allocated(matrix_a)) then
                deallocate(matrix_a)
            end if
            allocate(matrix_a(counter_sum, counter_sum))

            matrix_a = 0

            m_i = 1
            m_ii = 1
            do i = first_sphere, last_sphere
                do ii = first_sphere, last_sphere

                    if (i .eq. ii) then
                        sphere_parameter_no = simulation_parameter%sphere_list(i)%sphere_parameter_index
                        n_range_j = simulation_parameter%sphere_parameter_list(sphere_parameter_no)%n_range

                        do n = n_range_j(1), n_range_j(2)
                            matrix_a(m_ii, m_i) = 1D0 / a_n(i)%item(n)
                            matrix_a(m_ii, m_i+1) = 1D0 / b_n(i)%item(n)
                            m_i = m_i + 2
                        end do
                    else
                        d_0_j = simulation_parameter%sphere_list(i)%d_0_j
                        sphere_parameter_no = simulation_parameter%sphere_list(i)%sphere_parameter_index
                        n_range_j = simulation_parameter%sphere_parameter_list(sphere_parameter_no)%n_range

                        d_0_l = simulation_parameter%sphere_list(i)%d_0_j
                        sphere_parameter_no = simulation_parameter%sphere_list(i)%sphere_parameter_index
                        n_range_l = simulation_parameter%sphere_parameter_list(sphere_parameter_no)%n_range

                        ! todo: bug fix: freeze m, n per row of matrix_a
                        do n = n_range_j(1), n_range_j(2)
                            do m = -n, n
                                do nu = n_range_l(1), n_range_l(2)
                                    do mu = -nu, nu
                                        matrix_a(m_ii, m_i) = a_nmnumu(i, ii)%item(n)%item(m)%item(nu)%item(mu)
                                        matrix_a(m_ii, m_i+1) = b_nmnumu(i, ii)%item(n)%item(m)%item(nu)%item(mu)
                                        m_i = m_i + 2
                                    end do
                                end do
                            end do
                        end do
                    end if
                end do
                m_ii = m_ii + 1
            end do

        end subroutine lib_mie_type_func_solver_get_matrix_a

        function lib_mie_type_functions_test_functions() result(rv)
            implicit none
            ! dummy
            integer :: rv

            ! auxiliary

            contains
                function test_lib_mie_type_func_solver_get_matrix_a() result(rv)
                    implicit none
                    ! dummy
                    logical :: rv

                    ! auxiliary

                end function test_lib_mie_type_func_solver_get_matrix_a

        end function lib_mie_type_functions_test_functions

end module lib_mie_type_functions
