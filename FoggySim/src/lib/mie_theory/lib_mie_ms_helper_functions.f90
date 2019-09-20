module lib_mie_ms_helper_functions
    use libmath
    use lib_mie_type
    use lib_mie_ss_helper_functions
    use lib_mie_vector_spherical_harmonics
    implicit none

    contains

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
    function lib_mie_ms_hf_get_vector_b(simulation_parameter) result(vector_b)
        implicit none
        ! dummy
        type(lib_mie_simulation_parameter_type), intent(in) :: simulation_parameter

        double complex, dimension(:), allocatable :: vector_b

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

        counter_sum = sum(counter)

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

    end function lib_mie_ms_hf_get_vector_b

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
    function lib_mie_ms_hf_get_matrix_a(simulation_parameter) result(matrix_a)
        implicit none
        ! dummy
        type(lib_mie_simulation_parameter_type), intent(in) :: simulation_parameter

        double complex, dimension(:,:), allocatable :: matrix_a

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
        !$OMP PARALLEL DO PRIVATE(i, ii, d_0_j, d_0_l, d_j_l, sphere_parameter_no, n_range_l, n_range_j)
        do i = first_sphere, last_sphere
            do ii = first_sphere, last_sphere

                if (i .eq. ii) then
                    sphere_parameter_no = simulation_parameter%sphere_list(i)%sphere_parameter_index
                    n_range_j = simulation_parameter%sphere_parameter_list(sphere_parameter_no)%n_range

                    a_n(i) = simulation_parameter%sphere_parameter_list(sphere_parameter_no)%a_n
                    b_n(i) = simulation_parameter%sphere_parameter_list(sphere_parameter_no)%b_n

                    counter(i, ii) = n_range_j(2) - n_range_j(1) + 1
                else
                    d_0_j = simulation_parameter%sphere_list(i)%d_0_j
                    sphere_parameter_no = simulation_parameter%sphere_list(i)%sphere_parameter_index
                    n_range_j = simulation_parameter%sphere_parameter_list(sphere_parameter_no)%n_range

                    d_0_l = simulation_parameter%sphere_list(i)%d_0_j
                    sphere_parameter_no = simulation_parameter%sphere_list(i)%sphere_parameter_index
                    n_range_l = simulation_parameter%sphere_parameter_list(sphere_parameter_no)%n_range

                    d_j_l = d_0_l - d_0_j
                    call lib_mie_vector_spherical_harmonics_translation_coefficient(d_j_l, &
                            n_range_j, n_range_l, z_selector, &
                            a_nmnumu(i, ii), b_nmnumu(i, ii))

                    counter(i, ii) = ( (1 + n_range_j(2))**2 - n_range_j(1)**2 ) &
                                 * ( (1 + n_range_l(2))**2 - n_range_l(1)**2 )
                end if
            end do
        end do
        !$OMP END PARALLEL DO

        counter_sum = sum(counter(first_sphere, :))

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

    end function lib_mie_ms_hf_get_matrix_a

    function lib_mie_ms_helper_functions_test_functions() result(rv)
        implicit none
        ! dummy
        integer :: rv

        ! auxiliary

        contains
            function test_lib_mie_ms_hf_get_matrix_a() result(rv)
                implicit none
                ! dummy
                logical :: rv

                ! auxiliary

            end function test_lib_mie_ms_hf_get_matrix_a

    end function

end module lib_mie_ms_helper_functions
