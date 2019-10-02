module lib_mie_ms_solver
    use libmath
    use lib_mie_vector_spherical_harmonics
    use lib_mie_type
    use lib_mie_type_functions
    implicit none

    private

    public :: lib_mie_ms_solver_get_vector_b
    public :: lib_mie_ms_solver_get_vector_x
    public :: lib_mie_ms_solver_set_sphere_parameter_ab_nm
    public :: lib_mie_ms_solver_calculate_vector_b

    public :: lib_mie_ms_solver_test_functions

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
        subroutine lib_mie_ms_solver_get_vector_b(simulation_parameter, vector_b)
            implicit none
            ! dummy
            type(lib_mie_simulation_parameter_type), intent(in) :: simulation_parameter

            double complex, dimension(:), allocatable, intent(inout) :: vector_b

            ! auxiliary
            integer :: i

            integer :: first_sphere
            integer :: last_sphere

            integer :: first
            integer :: last

            integer :: sphere_parameter_no

            integer, dimension(:), allocatable :: counter
            integer :: counter_sum

            type(cartesian_coordinate_real_type) :: d_0_j
            integer(kind=4), dimension(2) :: n_range

            type(list_list_cmplx), dimension(:), allocatable :: p_nm
            type(list_list_cmplx), dimension(:), allocatable :: q_nm

            double complex, dimension(:), allocatable :: buffer_array


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
                if (size(vector_b, 1) .ne. counter_sum) then
                    deallocate(vector_b)
                    allocate(vector_b(counter_sum))
                end if
            else
                allocate(vector_b(counter_sum))
            end if

            !$OMP PARALLEL DO PRIVATE(i, first, last, buffer_array)
            do i = first_sphere, last_sphere
                call make_array(p_nm(i), buffer_array)
                first = 2 * (sum(counter(first_sphere:i)) - counter(i)) + 1
                last = first + counter(i) - 1
                vector_b(first:last) = buffer_array

                call make_array(q_nm(i), buffer_array)
                first = last + 1
                last = first + counter(i) - 1
                vector_b(first:last) = buffer_array
            end do
            !$OMP END PARALLEL DO

        end subroutine lib_mie_ms_solver_get_vector_b

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
        !   vector_x: double complex, dimension(:)
        !       a_nm, b_nm coefficients of all spheres
        !
        ! Reference: [1] Computation of scattering from clusters of spheres using the fast multipole method, Nail A. Gumerov, and Ramani Duraiswami
        !            [2] Electromagnetic scatteringby an aggregate of spheres, Yu-lin Xu
        subroutine lib_mie_ms_solver_get_vector_x(simulation_parameter, vector_x)
            implicit none
            ! dummy
            type(lib_mie_simulation_parameter_type), intent(in) :: simulation_parameter

            double complex, dimension(:), allocatable, intent(inout) :: vector_x

            ! auxiliary
            integer :: i

            integer :: first_sphere
            integer :: last_sphere

            integer :: first
            integer :: last

            integer :: sphere_parameter_no

            integer, dimension(:), allocatable :: counter
            integer :: counter_sum

            integer(kind=4), dimension(2) :: n_range

            double complex, dimension(:), allocatable :: buffer_array

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
                if (size(vector_x, 1) .ne. counter_sum) then
                    deallocate(vector_x)
                    allocate(vector_x(counter_sum))
                end if
            else
                allocate(vector_x(counter_sum))
            end if

            !$OMP PARALLEL DO PRIVATE(i, first, last, buffer_array)
            do i = first_sphere, last_sphere
                call make_array(simulation_parameter%sphere_list(i)%a_nm, buffer_array)
                first = 2 * (sum(counter(first_sphere:i)) - counter(i)) + 1
                last = first + counter(i) - 1
                vector_x(first:last) = buffer_array

                call make_array(simulation_parameter%sphere_list(i)%b_nm, buffer_array)
                first = last + 1
                last = first + counter(i) - 1
                vector_x(first:last) = buffer_array
            end do
            !$OMP END PARALLEL DO

        end subroutine lib_mie_ms_solver_get_vector_x

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
        subroutine lib_mie_ms_solver_set_sphere_parameter_ab_nm(vector_x, simulation_parameter)
            implicit none
            ! dummy
            double complex, dimension(:), intent(in) :: vector_x

            type(lib_mie_simulation_parameter_type), intent(inout) :: simulation_parameter

            ! auxiliary
            integer :: i

            integer :: first
            integer :: last

            integer :: first_sphere
            integer :: last_sphere

            integer :: sphere_parameter_no

            integer, dimension(:,:), allocatable :: n_range
            integer, dimension(:), allocatable :: counter

            type(list_list_cmplx) :: buffer_list

            first_sphere = lbound(simulation_parameter%sphere_list, 1)
            last_sphere = ubound(simulation_parameter%sphere_list, 1)

            allocate(n_range(first_sphere:last_sphere, 2))
            allocate(counter(first_sphere:last_sphere))

            counter = 0
            do i = first_sphere, last_sphere
                sphere_parameter_no = simulation_parameter%sphere_list(i)%sphere_parameter_index
                n_range(i, :) = simulation_parameter%sphere_parameter_list(sphere_parameter_no)%n_range

                counter(i) = (1 + n_range(i, 2))**2 - n_range(i, 1)**2
            end do

            !$OMP PARALLEL DO PRIVATE(i, first, last, buffer_list)
            do i = first_sphere, last_sphere
                first = 2 * (sum(counter(first_sphere:i)) - counter(i)) + 1
                last = first + counter(i) - 1
                call make_list(vector_x(first:last), &
                               n_range(i, 1), n_range(i, 2) - n_range(i, 1) + 1, &
                               buffer_list)
                call move_alloc(buffer_list%item, simulation_parameter%sphere_list(i)%a_nm%item)

                first = last + 1
                last = first + counter(i) - 1
                call make_list(vector_x(first:last), &
                               n_range(i, 1), n_range(i, 2) - n_range(i, 1) + 1, &
                               buffer_list)
                call move_alloc(buffer_list%item, simulation_parameter%sphere_list(i)%b_nm%item)
            end do
            !$OMP END PARALLEL DO

        end subroutine lib_mie_ms_solver_set_sphere_parameter_ab_nm

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
        !       result of the matrix vector multiplication
        !
        ! Reference: [1] Computation of scattering from clusters of spheres using the fast multipole method, Nail A. Gumerov, and Ramani Duraiswami
        !            [2] Electromagnetic scatteringby an aggregate of spheres, Yu-lin Xu
        subroutine lib_mie_ms_solver_calculate_vector_b(simulation_parameter, vector_x, vector_b, vector_b_t)
            implicit none
            ! dummy
            type(lib_mie_simulation_parameter_type), intent(in) :: simulation_parameter
            double complex, dimension(:), intent(in) :: vector_x

            double complex, dimension(:), allocatable, intent(inout) :: vector_b
            double complex, dimension(:), allocatable, intent(inout) :: vector_b_t

            ! auxiliary
            integer :: i
            integer :: j
            integer :: l
            integer :: n
            integer :: m

            integer :: first
            integer :: last

            integer :: first_sphere
            integer :: last_sphere

            integer :: sphere_parameter_no

            integer, dimension(:), allocatable :: counter

            integer(kind=1) :: z_selector

            type(cartesian_coordinate_real_type) :: d_0_j
            type(cartesian_coordinate_real_type) :: d_0_l
            type(cartesian_coordinate_real_type) :: d_j_l
            integer(kind=4), dimension(2) :: n_range_j
            integer(kind=4), dimension(2) :: n_range_l

            type(list_cmplx) :: a_n
            type(list_cmplx) :: b_n
            type(list_4_cmplx) :: a_nmnumu
            type(list_4_cmplx) :: b_nmnumu

            type(list_list_cmplx) :: buffer_1_nm
            type(list_list_cmplx) :: buffer_2_nm

            type(list_list_cmplx) :: buffer_x_1_nm
            type(list_list_cmplx) :: buffer_x_2_nm

            type(list_list_cmplx) :: buffer_b_1_nm
            type(list_list_cmplx) :: buffer_b_2_nm

            double complex, dimension(:), allocatable :: buffer_array

            first_sphere = lbound(simulation_parameter%sphere_list, 1)
            last_sphere = ubound(simulation_parameter%sphere_list, 1)

            allocate(counter(first_sphere:last_sphere))

            z_selector = simulation_parameter%spherical_harmonics%z_selector_translation

            counter = 0
            do i = first_sphere, last_sphere
                sphere_parameter_no = simulation_parameter%sphere_list(i)%sphere_parameter_index
                n_range_j = simulation_parameter%sphere_parameter_list(sphere_parameter_no)%n_range

                counter(i) = (1 + n_range_j(2))**2 - n_range_j(1)**2
            end do

!            !$OMP PARALLEL DO PRIVATE(i, ii, d_0_j, d_0_l, d_j_l, sphere_parameter_no, n_range_l, n_range_j)
            do j = first_sphere, last_sphere
                !
                !
                !     Sphere 1     Sphere 2             Sphere 2
                !     v            v                    v
                !  | 1 0  a^1_n*A_nmnumu(2,1)  a^1_n*B_nmnumu(2,1)  ... | | a^1_nm |   | a^1_n p^11_nm | < Sphere 1   <- first
                !  | 0 1  b^1_n*B_nmnumu(2,1)  b^1_n*A_nmnumu(2,1)  ... | | b^1_nm |   | b^1_n p^11_nm | < Sphere 1   <- last
                !  | ...           1                    0           ... |*|  ...   | = | a^2_n p^22_nm |
                !  | ...           0                    1           ... | |  ...   |   | b^2_n p^22_nm |
                !                         ^                                   ^            ^
                !                         Matrix A                            vector_x     vector_b
                !
                !
                !  first: first element of sphere j (of vector x and b: A is symmetrical)
                !  last : last element of sphere j (of vector x and b: A is symmetrical)
                !
                !  buffer_x_1 = a^l_nm
                !  buffer_x_2 = b^l_nm
                !
                !  buffer_b_1 = a^l_n p^ll_nm
                !  buffer_b_2 = b^l_n p^ll_nm
                !
                first = 2 * (sum(counter(first_sphere:j)) - counter(j)) + 1

                sphere_parameter_no = simulation_parameter%sphere_list(j)%sphere_parameter_index
                n_range_j = simulation_parameter%sphere_parameter_list(sphere_parameter_no)%n_range

                a_n = simulation_parameter%sphere_parameter_list(sphere_parameter_no)%a_n
                b_n = simulation_parameter%sphere_parameter_list(sphere_parameter_no)%b_n

                do l = first_sphere, last_sphere
                    last = first + 2*counter(l) - 1

                    if (j .eq. l) then
                        vector_b(first:last) = vector_x(first:last)
                    else
                        d_0_j = simulation_parameter%sphere_list(j)%d_0_j
                        sphere_parameter_no = simulation_parameter%sphere_list(j)%sphere_parameter_index
                        n_range_j = simulation_parameter%sphere_parameter_list(sphere_parameter_no)%n_range

                        d_0_l = simulation_parameter%sphere_list(l)%d_0_j
                        sphere_parameter_no = simulation_parameter%sphere_list(l)%sphere_parameter_index
                        n_range_l = simulation_parameter%sphere_parameter_list(sphere_parameter_no)%n_range

                        d_j_l = d_0_l - d_0_j
                        call lib_math_vector_spherical_harmonics_translation_coefficient(d_j_l, &
                                n_range_j, n_range_l, z_selector, &
                                a_nmnumu, b_nmnumu)

                        call init_list(buffer_1_nm, n_range_l(1), n_range_l(2)-n_range_l(1)+1, dcmplx(0,0))
                        call init_list(buffer_2_nm, n_range_l(1), n_range_l(2)-n_range_l(1)+1, dcmplx(0,0))

                        call init_list(buffer_b_1_nm, n_range_l(1), n_range_l(2)-n_range_l(1)+1, dcmplx(0,0))
                        call init_list(buffer_b_2_nm, n_range_l(1), n_range_l(2)-n_range_l(1)+1, dcmplx(0,0))

                        do n = n_range_j(1), n_range_j(2)
                            do m = -n, n
                                buffer_1_nm = a_n%item(n) * a_nmnumu%item(n)%item(m)
                                buffer_2_nm = a_n%item(n) * b_nmnumu%item(n)%item(m)

                                call make_list(vector_x(first:first+counter(l)-1), &
                                               n_range_l(1), n_range_l(2)-n_range_l(1)+1, &
                                               buffer_x_1_nm)

                                call make_list(vector_x(first+counter(l):first+2*counter(l)-1), &
                                               n_range_l(1), n_range_l(2)-n_range_l(1)+1, &
                                               buffer_x_2_nm)

                                buffer_b_1_nm%item(n)%item(m) = sum(buffer_1_nm * buffer_x_1_nm)
                                buffer_b_2_nm%item(n)%item(m) = sum(buffer_2_nm * buffer_x_2_nm)
                            end do
                        end do

                        call make_array(buffer_b_1_nm, buffer_array)
                        vector_b(first:first+counter(l)-1) = buffer_array

                        call make_array(buffer_b_2_nm, buffer_array)
                        vector_b(first+counter(l):first+2*counter(l)-1) = buffer_array
                    end if
                end do
            end do
!            !$OMP END PARALLEL DO

        end subroutine lib_mie_ms_solver_calculate_vector_b

        function lib_mie_ms_solver_test_functions() result(rv)
            implicit none
            ! dummy
            integer :: rv

            rv = 0

            if (.not. test_lib_mie_ms_solver_calculate_vector_b()) rv = rv + 1

            contains

                function test_lib_mie_ms_solver_calculate_vector_b() result(rv)
                    implicit none
                    ! dummy
                    logical :: rv

                    ! auxiliary
                    integer :: i
                    type(lib_mie_simulation_parameter_type) :: simulation
                    double complex, dimension(:), allocatable :: vector_x

                    double complex, dimension(:), allocatable :: vector_b
                    double complex, dimension(:), allocatable :: vector_b_t

                    integer :: no_of_elements

                    double precision :: n_medium

                    ! illumination
                    double precision :: lambda_0
                    double precision :: e_field_0

                    double precision, dimension(1) :: g
                    type(cartesian_coordinate_real_type), dimension(1) :: k
                    type(cartesian_coordinate_real_type), dimension(1) :: d_0_i

                    ! sphere parameter
                    type(lib_mie_sphere_parameter_type) :: sphere_para
                    double precision :: r_particle
                    double complex :: n_particle
                    integer, dimension(2) :: n_range

                    type(cartesian_coordinate_real_type) :: d_0_j

                    lambda_0 = 1 * unit_mu
                    e_field_0 = 1
                    n_medium = 1

                    g(:) = 1
                    k(1) = make_cartesian(0d0,0d0,1d0)
                    d_0_i(1) = make_cartesian(0d0,0d0,0d0) * unit_mu

                    simulation%illumination = lib_mie_type_func_get_plane_wave_illumination(lambda_0, e_field_0, &
                                                                                            g, k, &
                                                                                            d_0_i)

                    simulation%refractive_index_medium = n_medium

                    r_particle = 2 * unit_mu
                    n_particle = 1.33
                    n_range = (/ 1, 3 /)

                    sphere_para = lib_mie_type_func_get_sphere_parameter(lambda_0, n_medium, &
                                                                         r_particle, n_particle, n_range)

                    allocate(simulation%sphere_parameter_list(1))
                    simulation%sphere_parameter_list(1) = sphere_para


                    allocate(simulation%sphere_list(2))
                    simulation%sphere_list(:)%sphere_parameter_index = 1

                    d_0_j = make_cartesian(-2d0,0d0,1d0) * unit_mu
                    simulation%sphere_list(1)%d_0_j = d_0_j

                    d_0_j = make_cartesian(2d0,0d0,1d0) * unit_mu
                    simulation%sphere_list(2)%d_0_j = d_0_j

                    simulation%spherical_harmonics%z_selector_incident_wave = 1
                    simulation%spherical_harmonics%z_selector_scatterd_wave = 3
                    simulation%spherical_harmonics%z_selector_translation = 1

                    no_of_elements = size(simulation%sphere_list) * 2 * ( (1+n_range(2))**2 - n_range(1)**2 )
                    allocate(vector_x(no_of_elements))
                    allocate(vector_b(no_of_elements))

                    do i = 1, no_of_elements
                        vector_x(i) = i
                    end do

                    call lib_mie_ms_solver_calculate_vector_b(simulation, vector_x, vector_b, vector_b_t)

                end function


        end function
end module lib_mie_ms_solver
