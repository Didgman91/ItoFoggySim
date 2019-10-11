module lib_mie_ms_solver_interface
    use libmath
    use lib_mie_type
    use lib_mie_type_functions
    implicit none

    private

    public :: lib_mie_ms_solver_get_vector_b
    public :: lib_mie_ms_solver_get_vector_x
    public :: lib_mie_ms_solver_set_sphere_parameter_ab_nm
    public :: lib_mie_ms_solver_calculate_vector_b
    public :: lib_mie_ms_solver_allocate_vector_x_b

    public :: lib_mie_ms_solver_interface_test_functions

    ! --- interface ---
    interface lib_mie_ms_solver_calculate_vector_b
        module procedure lib_mie_ms_solver_calculate_vector_b_with_transposed
        module procedure lib_mie_ms_solver_calculate_vector_b_without_transposed_all
        module procedure lib_mie_ms_solver_calculate_vector_b_without_transposed_select
    end interface

    interface lib_mie_ms_solver_get_vector_b
        module procedure lib_mie_ms_solver_get_vector_b_all
        module procedure lib_mie_ms_solver_get_vector_b_selection
    end interface

    interface lib_mie_ms_solver_get_vector_x
        module procedure lib_mie_ms_solver_get_vector_x_all
        module procedure lib_mie_ms_solver_get_vector_x_selection
    end interface

    interface lib_mie_ms_solver_set_sphere_parameter_ab_nm
        module procedure lib_mie_ms_solver_set_sphere_parameter_ab_nm_all
        module procedure lib_mie_ms_solver_set_sphere_parameter_ab_nm_selection
    end interface

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
        subroutine lib_mie_ms_solver_get_vector_b_all(simulation_parameter, vector_b)
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

            integer :: counter
            integer :: counter_sum

            type(cartesian_coordinate_real_type) :: d_0_j
            integer(kind=4), dimension(2) :: n_range

            type(list_cmplx) :: a_n
            type(list_cmplx) :: b_n

            type(list_list_cmplx), dimension(:), allocatable :: p_nm
            type(list_list_cmplx), dimension(:), allocatable :: q_nm

            double complex, dimension(:), allocatable :: buffer_array


            first_sphere = lbound(simulation_parameter%sphere_list, 1)
            last_sphere = ubound(simulation_parameter%sphere_list, 1)

            allocate(p_nm(first_sphere:last_sphere))
            allocate(q_nm(first_sphere:last_sphere))

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
            end do
            !$OMP END PARALLEL DO

            n_range = simulation_parameter%spherical_harmonics%n_range
            counter = (1 + n_range(2))**2 - n_range(1)**2
            counter_sum = (last_sphere - first_sphere + 1) * (2 * counter)

            if (allocated(vector_b)) then
                if (size(vector_b, 1) .ne. counter_sum) then
                    deallocate(vector_b)
                    allocate(vector_b(counter_sum))
                end if
            else
                allocate(vector_b(counter_sum))
            end if

            !$OMP PARALLEL DO PRIVATE(i, first, last, buffer_array, a_n, b_n)
            do i = first_sphere, last_sphere
                sphere_parameter_no = simulation_parameter%sphere_list(i)%sphere_parameter_index
                a_n = simulation_parameter%sphere_parameter_list(sphere_parameter_no)%a_n
                b_n = simulation_parameter%sphere_parameter_list(sphere_parameter_no)%b_n

                p_nm(i) = a_n * p_nm(i)
                q_nm(i) = b_n * p_nm(i)

                call make_array(p_nm(i), buffer_array, n_range(1), n_range(2) - n_range(1) + 1)
                first = 2 * (i - 1) * counter + 1
                last = first + counter - 1
                vector_b(first:last) = buffer_array

                call make_array(q_nm(i), buffer_array, n_range(1), n_range(2) - n_range(1) + 1)
                first = last + 1
                last = first + counter - 1
                vector_b(first:last) = buffer_array
            end do
            !$OMP END PARALLEL DO

        end subroutine lib_mie_ms_solver_get_vector_b_all

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
        !   sphere_no_list: integer, dimension(:)
        !       list of sphere numbers (index number of simulation_parameter%sphere_list)
        !
        ! Returns
        ! ----
        !   vector_b: double complex, dimension(:)
        !
        ! Reference: [1] Computation of scattering from clusters of spheres using the fast multipole method, Nail A. Gumerov, and Ramani Duraiswami
        !            [2] Electromagnetic scatteringby an aggregate of spheres, Yu-lin Xu
        subroutine lib_mie_ms_solver_get_vector_b_selection(simulation_parameter, sphere_no_list, vector_b)
            implicit none
            ! dummy
            type(lib_mie_simulation_parameter_type), intent(in) :: simulation_parameter
            integer, dimension(:), intent(in) :: sphere_no_list

            double complex, dimension(:), allocatable, intent(inout) :: vector_b

            ! auxiliary
            integer :: i
            integer :: ii

            integer :: first_sphere
            integer :: last_sphere

            integer :: first
            integer :: last

            integer :: sphere_parameter_no

            integer :: counter

            type(cartesian_coordinate_real_type) :: d_0_j
            integer(kind=4), dimension(2) :: n_range

            type(list_cmplx) :: a_n
            type(list_cmplx) :: b_n

            type(list_list_cmplx), dimension(:), allocatable :: p_nm
            type(list_list_cmplx), dimension(:), allocatable :: q_nm

            double complex, dimension(:), allocatable :: buffer_array


            first_sphere = lbound(simulation_parameter%sphere_list, 1)
            last_sphere = ubound(simulation_parameter%sphere_list, 1)

            allocate(p_nm(lbound(sphere_no_list, 1):ubound(sphere_no_list, 1)))
            allocate(q_nm(lbound(sphere_no_list, 1):ubound(sphere_no_list, 1)))

            ! create vector_b
            counter = 0
            !$OMP PARALLEL DO PRIVATE(i, ii, d_0_j, sphere_parameter_no, n_range)
            do i = lbound(sphere_no_list, 1), ubound(sphere_no_list, 1)
                ii = sphere_no_list(i)
                d_0_j = simulation_parameter%sphere_list(ii)%d_0_j
                sphere_parameter_no = simulation_parameter%sphere_list(ii)%sphere_parameter_index
                n_range = simulation_parameter%sphere_parameter_list(sphere_parameter_no)%n_range

                call lib_mie_ss_hf_get_p_q_j_j(simulation_parameter%illumination, &
                                               simulation_parameter%refractive_index_medium, &
                                               d_0_j, n_range, p_nm(i), q_nm(i))
            end do
            !$OMP END PARALLEL DO

            n_range = simulation_parameter%spherical_harmonics%n_range
            counter = (1 + n_range(2))**2 - n_range(1)**2

            !$OMP PARALLEL DO PRIVATE(i, ii, first, last, buffer_array, a_n, b_n)
            do i = lbound(sphere_no_list, 1), ubound(sphere_no_list, 1)
                ii = sphere_no_list(i)
                sphere_parameter_no = simulation_parameter%sphere_list(ii)%sphere_parameter_index
                a_n = simulation_parameter%sphere_parameter_list(sphere_parameter_no)%a_n
                b_n = simulation_parameter%sphere_parameter_list(sphere_parameter_no)%b_n

                p_nm(i) = a_n * p_nm(i)
                q_nm(i) = b_n * p_nm(i)

                call make_array(p_nm(i), buffer_array, n_range(1), n_range(2) - n_range(1) + 1)
                first = 2 * (ii - 1) * counter + 1
                last = first + counter - 1
                vector_b(first:last) = buffer_array

                call make_array(q_nm(i), buffer_array, n_range(1), n_range(2) - n_range(1) + 1)
                first = last + 1
                last = first + counter - 1
                vector_b(first:last) = buffer_array
            end do
            !$OMP END PARALLEL DO

        end subroutine lib_mie_ms_solver_get_vector_b_selection

        subroutine lib_mie_ms_solver_allocate_vector_x_b(simulation_parameter, vector)
            implicit none
            ! dummy
            type(lib_mie_simulation_parameter_type), intent(in) :: simulation_parameter

            double complex, dimension(:), allocatable, intent(inout) :: vector

            ! auxiliary
            integer, dimension(2) :: n_range
            integer(kind=8) :: counter
            integer(kind=8) :: counter_sum

            n_range = simulation_parameter%spherical_harmonics%n_range
            counter = (1 + n_range(2))**2 - n_range(1)**2
            counter_sum = size(simulation_parameter%sphere_list) * (2 * counter)

            if (allocated(vector)) then
                if (size(vector, 1) .ne. counter_sum) then
                    deallocate(vector)
                    allocate(vector(counter_sum))
                end if
            else
                allocate(vector(counter_sum))
            end if
        end subroutine lib_mie_ms_solver_allocate_vector_x_b

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
        subroutine lib_mie_ms_solver_get_vector_x_all(simulation_parameter, vector_x)
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

            integer :: counter
            integer :: counter_sum

            integer(kind=4), dimension(2) :: n_range

            double complex, dimension(:), allocatable :: buffer_array

            first_sphere = lbound(simulation_parameter%sphere_list, 1)
            last_sphere = ubound(simulation_parameter%sphere_list, 1)

            ! create vector_x
            n_range = simulation_parameter%spherical_harmonics%n_range

            counter = (1 + n_range(2))**2 - n_range(1)**2
            counter_sum = (last_sphere - first_sphere + 1) * 2 * counter

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
                call make_array(simulation_parameter%sphere_list(i)%a_nm, buffer_array, &
                                n_range(1), n_range(2) - n_range(1) + 1)
                first = 2 * (i - first_sphere) * counter + 1
                last = first + counter - 1
                vector_x(first:last) = buffer_array

                call make_array(simulation_parameter%sphere_list(i)%b_nm, buffer_array, &
                                n_range(1), n_range(2) - n_range(1) + 1)
                first = last + 1
                last = first + counter - 1
                vector_x(first:last) = buffer_array
            end do
            !$OMP END PARALLEL DO

        end subroutine lib_mie_ms_solver_get_vector_x_all

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
        !   sphere_no_list: integer, dimension(:)
        !       list of sphere numbers (index number of simulation_parameter%sphere_list)
        !
        ! Returns
        ! ----
        !   vector_x: double complex, dimension(:)
        !       a_nm, b_nm coefficients of all spheres
        !
        ! Reference: [1] Computation of scattering from clusters of spheres using the fast multipole method, Nail A. Gumerov, and Ramani Duraiswami
        !            [2] Electromagnetic scatteringby an aggregate of spheres, Yu-lin Xu
        subroutine lib_mie_ms_solver_get_vector_x_selection(simulation_parameter, sphere_no_list, vector_x)
            implicit none
            ! dummy
            type(lib_mie_simulation_parameter_type), intent(in) :: simulation_parameter
            integer, dimension(:), intent(in) :: sphere_no_list

            double complex, dimension(:), allocatable, intent(inout) :: vector_x

            ! auxiliary
            integer :: i
            integer :: ii

            integer :: first_sphere
            integer :: last_sphere

            integer :: first
            integer :: last

            integer :: counter

            integer(kind=4), dimension(2) :: n_range

            double complex, dimension(:), allocatable :: buffer_array

            first_sphere = lbound(simulation_parameter%sphere_list, 1)
            last_sphere = ubound(simulation_parameter%sphere_list, 1)

            ! create vector_x
            n_range = simulation_parameter%spherical_harmonics%n_range

            counter = (1 + n_range(2))**2 - n_range(1)**2

            !$OMP PARALLEL DO PRIVATE(i, ii, first, last, buffer_array)
            do i = lbound(sphere_no_list, 1), ubound(sphere_no_list, 1)
                ii = sphere_no_list(i)

                call make_array(simulation_parameter%sphere_list(ii)%a_nm, buffer_array, &
                                n_range(1), n_range(2) - n_range(1) + 1)
                first = 2 * (ii - first_sphere) * counter + 1
                last = first + counter - 1
                vector_x(first:last) = buffer_array

                call make_array(simulation_parameter%sphere_list(ii)%b_nm, buffer_array, &
                                n_range(1), n_range(2) - n_range(1) + 1)
                first = last + 1
                last = first + counter - 1
                vector_x(first:last) = buffer_array
            end do
            !$OMP END PARALLEL DO

        end subroutine lib_mie_ms_solver_get_vector_x_selection

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
        subroutine lib_mie_ms_solver_set_sphere_parameter_ab_nm_all(vector_x, simulation_parameter)
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

            integer, dimension(2) :: n_range
            integer :: counter

            type(list_list_cmplx) :: buffer_list

            first_sphere = lbound(simulation_parameter%sphere_list, 1)
            last_sphere = ubound(simulation_parameter%sphere_list, 1)

            n_range = simulation_parameter%spherical_harmonics%n_range
            counter = (1 + n_range(2))**2 - n_range(1)**2

            !$OMP PARALLEL DO PRIVATE(i, first, last, buffer_list)
            do i = first_sphere, last_sphere
                first = 2 * (i - first_sphere) * counter + 1
                last = first + counter - 1
                call make_list(vector_x(first:last), &
                               n_range(1), n_range(2) - n_range(1) + 1, &
                               buffer_list)
                call remove_zeros(buffer_list)
                call move_alloc(buffer_list%item, simulation_parameter%sphere_list(i)%a_nm%item)

                first = last + 1
                last = first + counter - 1
                call make_list(vector_x(first:last), &
                               n_range(1), n_range(2) - n_range(1) + 1, &
                               buffer_list)
                call remove_zeros(buffer_list)
                call move_alloc(buffer_list%item, simulation_parameter%sphere_list(i)%b_nm%item)
            end do
            !$OMP END PARALLEL DO

        end subroutine lib_mie_ms_solver_set_sphere_parameter_ab_nm_all

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
        !   sphere_no_list: integer, dimension(:)
        !       list of sphere numbers (index number of simulation_parameter%sphere_list)

        !
        ! Returns
        ! ----
        !   simulation_parameter: type(lib_mie_simulation_parameter_type)
        !
        ! Reference: [1] Computation of scattering from clusters of spheres using the fast multipole method, Nail A. Gumerov, and Ramani Duraiswami
        !            [2] Electromagnetic scatteringby an aggregate of spheres, Yu-lin Xu
        subroutine lib_mie_ms_solver_set_sphere_parameter_ab_nm_selection(vector_x, sphere_no_list, simulation_parameter)
            implicit none
            ! dummy
            double complex, dimension(:), intent(in) :: vector_x
            integer, dimension(:), intent(in) :: sphere_no_list

            type(lib_mie_simulation_parameter_type), intent(inout) :: simulation_parameter


            ! auxiliary
            integer :: i
            integer :: ii

            integer :: first
            integer :: last

            integer :: first_sphere
            integer :: last_sphere

            integer, dimension(2) :: n_range
            integer :: counter

            type(list_list_cmplx) :: buffer_list

            first_sphere = lbound(simulation_parameter%sphere_list, 1)
            last_sphere = ubound(simulation_parameter%sphere_list, 1)

            n_range = simulation_parameter%spherical_harmonics%n_range
            counter = (1 + n_range(2))**2 - n_range(1)**2

            !$OMP PARALLEL DO PRIVATE(i, ii, first, last, buffer_list)
            do i = lbound(sphere_no_list, 1), ubound(sphere_no_list, 1)
                ii = sphere_no_list(i)

                first = 2 * (ii - first_sphere) * counter + 1
                last = first + counter - 1
                call make_list(vector_x(first:last), &
                               n_range(1), n_range(2) - n_range(1) + 1, &
                               buffer_list)
                call remove_zeros(buffer_list)
                call move_alloc(buffer_list%item, simulation_parameter%sphere_list(ii)%a_nm%item)

                first = last + 1
                last = first + counter - 1
                call make_list(vector_x(first:last), &
                               n_range(1), n_range(2) - n_range(1) + 1, &
                               buffer_list)
                call remove_zeros(buffer_list)
                call move_alloc(buffer_list%item, simulation_parameter%sphere_list(ii)%b_nm%item)
            end do
            !$OMP END PARALLEL DO

        end subroutine lib_mie_ms_solver_set_sphere_parameter_ab_nm_selection

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
        subroutine lib_mie_ms_solver_calculate_vector_b_without_transposed_all(simulation_parameter, vector_x, &
                                                        vector_b)
            implicit none
            ! dummy
            type(lib_mie_simulation_parameter_type), intent(in) :: simulation_parameter
            double complex, dimension(:), intent(in) :: vector_x

            double complex, dimension(:), allocatable, intent(inout) :: vector_b

            ! auxiliary
            integer :: j
            integer :: l
            integer :: n
            integer :: m

            integer :: first
            integer :: last

            integer :: first_sphere
            integer :: last_sphere

            integer :: sphere_parameter_no

            integer :: counter
            integer :: counter_sum

            integer(kind=1) :: z_selector

            type(cartesian_coordinate_real_type) :: d_0_j
            type(cartesian_coordinate_real_type) :: d_0_l
            type(cartesian_coordinate_real_type) :: d_j_l
            integer(kind=4), dimension(2) :: n_range
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

            !$  logical :: thread_first_run


            first_sphere = lbound(simulation_parameter%sphere_list, 1)
            last_sphere = ubound(simulation_parameter%sphere_list, 1)

            z_selector = simulation_parameter%spherical_harmonics%z_selector_translation

            n_range = simulation_parameter%spherical_harmonics%n_range

            counter = (1 + n_range(2))**2 - n_range(1)**2
            counter_sum = 2 * (last_sphere - first_sphere + 1) * counter

            if (allocated(vector_b)) then
                if (size(vector_b, 1) .ne. counter_sum) then
                    deallocate(vector_b)
                    allocate(vector_b(counter_sum))
                end if
            else
                allocate(vector_b(counter_sum))
            end if
            vector_b = 0


            !$  thread_first_run = .true.
            !$OMP PARALLEL DO PRIVATE(j, l, first, last, d_0_j, d_0_l, d_j_l, sphere_parameter_no, n_range_l, n_range_j) &
            !$OMP  PRIVATE(a_n, b_n, a_nmnumu, b_nmnumu, buffer_1_nm, buffer_2_nm, buffer_x_1_nm, buffer_x_2_nm) &
            !$OMP  PRIVATE(buffer_b_1_nm, buffer_b_2_nm, n, m, buffer_array) &
            !$OMP  FIRSTPRIVATE(thread_first_run)
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
                sphere_parameter_no = simulation_parameter%sphere_list(j)%sphere_parameter_index

                a_n = simulation_parameter%sphere_parameter_list(sphere_parameter_no)%a_n
                b_n = simulation_parameter%sphere_parameter_list(sphere_parameter_no)%b_n


                do l = first_sphere, last_sphere

                    first = 2 * (j - first_sphere) * counter + 1
                    if (j .eq. l) then
                        last = first + 2*counter - 1
                        vector_b(first:last) = vector_b(first:last) + vector_x(first:last)
                    else
                        d_0_j = simulation_parameter%sphere_list(j)%d_0_j
                        sphere_parameter_no = simulation_parameter%sphere_list(j)%sphere_parameter_index
                        n_range_j = simulation_parameter%sphere_parameter_list(sphere_parameter_no)%n_range

                        d_0_l = simulation_parameter%sphere_list(l)%d_0_j
                        sphere_parameter_no = simulation_parameter%sphere_list(l)%sphere_parameter_index
                        n_range_l = simulation_parameter%sphere_parameter_list(sphere_parameter_no)%n_range

                        n_range_j(1) = max(n_range_j(1), n_range(1))
                        n_range_j(2) = min(n_range_j(2), n_range(2))

                        n_range_l(1) = max(n_range_l(1), n_range(1))
                        n_range_l(2) = min(n_range_l(2), n_range(2))

                        !$  if (thread_first_run) then
                        !$    call fwig_thread_temp_init(4 * max(n_range_l(2), n_range_l(2)))     ! multi threaded
                        !$    thread_first_run = .false.
                        !$  endif

                        d_j_l = d_0_l - d_0_j
                        call lib_math_vector_spherical_harmonics_translation_coefficient(d_j_l, &
                                n_range_j, n_range_l, z_selector, &
                                a_nmnumu, b_nmnumu)

                        call init_list(buffer_1_nm, n_range_l(1), n_range_l(2)-n_range_l(1)+1, dcmplx(0,0))
                        call init_list(buffer_2_nm, n_range_l(1), n_range_l(2)-n_range_l(1)+1, dcmplx(0,0))

                        call init_list(buffer_b_1_nm, n_range(1), n_range(2)-n_range(1)+1, dcmplx(0,0))
                        call init_list(buffer_b_2_nm, n_range(1), n_range(2)-n_range(1)+1, dcmplx(0,0))

                        ! calculate vector b
                        call make_list(vector_x(first:first+counter-1), &
                                       n_range(1), n_range(2)-n_range(1)+1, &
                                       buffer_x_1_nm)
                        call remove_zeros(buffer_x_1_nm)

                        call make_list(vector_x(first+counter:first+2*counter-1), &
                                       n_range(1), n_range(2)-n_range(1)+1, &
                                       buffer_x_2_nm)
                        call remove_zeros(buffer_x_2_nm)

                        !$OMP PARALLEL DO PRIVATE(n, m) &
                        !$OMP  PRIVATE(buffer_1_nm, buffer_2_nm)
                        do n = n_range(1), n_range(2)
                            do m = -n, n
                                if (n_range(2) .le. n_range_j(2) &
                                    .and. n_range(1) .ge. n_range_j(1)) then
                                    buffer_1_nm = a_n%item(n) * a_nmnumu%item(n)%item(m)
                                    buffer_2_nm = a_n%item(n) * b_nmnumu%item(n)%item(m)

                                    buffer_b_1_nm%item(n)%item(m) = sum(buffer_1_nm * buffer_x_1_nm &
                                                                        + buffer_2_nm * buffer_x_2_nm)

                                    buffer_1_nm = b_n%item(n) * b_nmnumu%item(n)%item(m)
                                    buffer_2_nm = b_n%item(n) * a_nmnumu%item(n)%item(m)

                                    buffer_b_2_nm%item(n)%item(m) = sum(buffer_1_nm * buffer_x_1_nm &
                                                                        + buffer_2_nm * buffer_x_2_nm)
                                else
                                    buffer_b_1_nm%item(n)%item(m) = dcmplx(0,0)
                                    buffer_b_2_nm%item(n)%item(m) = dcmplx(0,0)
                                end if
                            end do
                        end do
                        !$OMP END PARALLEL DO

                        n_range = simulation_parameter%spherical_harmonics%n_range

                        call make_array(buffer_b_1_nm, buffer_array, n_range(1) , n_range(2) - n_range(1) + 1)
                        last = first + counter - 1
                        vector_b(first:last) = vector_b(first:last) + buffer_array

                        call make_array(buffer_b_2_nm, buffer_array, n_range(1) , n_range(2) - n_range(1) + 1)
                        first = last + 1
                        last = first + counter - 1
                        vector_b(first:last) = vector_b(first:last) + buffer_array
                    end if
                end do
            end do
            !$OMP END PARALLEL DO
        end subroutine lib_mie_ms_solver_calculate_vector_b_without_transposed_all

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
        subroutine lib_mie_ms_solver_calculate_vector_b_without_transposed_select(simulation_parameter, sphere_no_list, &
                                                                                     vector_x, vector_b)
            implicit none
            ! dummy
            type(lib_mie_simulation_parameter_type), intent(in) :: simulation_parameter
            integer, dimension(:), intent(in) :: sphere_no_list
            double complex, dimension(:), intent(in) :: vector_x

            double complex, dimension(:), allocatable, intent(inout) :: vector_b

            ! auxiliary
            integer :: i
            integer :: ii

            integer :: j
            integer :: l
            integer :: n
            integer :: m

            integer :: first
            integer :: last

            integer :: first_sphere
            integer :: last_sphere

            integer :: sphere_parameter_no

            integer :: counter
            integer :: counter_sum

            integer(kind=1) :: z_selector

            type(cartesian_coordinate_real_type) :: d_0_j
            type(cartesian_coordinate_real_type) :: d_0_l
            type(cartesian_coordinate_real_type) :: d_j_l
            integer(kind=4), dimension(2) :: n_range
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

            !$  logical :: thread_first_run


            first_sphere = lbound(simulation_parameter%sphere_list, 1)
            last_sphere = ubound(simulation_parameter%sphere_list, 1)

            z_selector = simulation_parameter%spherical_harmonics%z_selector_translation

            n_range = simulation_parameter%spherical_harmonics%n_range

            counter = (1 + n_range(2))**2 - n_range(1)**2


            !$  thread_first_run = .true.
            !$OMP PARALLEL DO PRIVATE(j, l, first, last, d_0_j, d_0_l, d_j_l, sphere_parameter_no, n_range_l, n_range_j) &
            !$OMP  PRIVATE(a_n, b_n, a_nmnumu, b_nmnumu, buffer_1_nm, buffer_2_nm, buffer_x_1_nm, buffer_x_2_nm) &
            !$OMP  PRIVATE(buffer_b_1_nm, buffer_b_2_nm, n, m, buffer_array) &
            !$OMP  FIRSTPRIVATE(thread_first_run)
            do i = lbound(sphere_no_list, 1), ubound(sphere_no_list, 1)
                j = sphere_no_list(i)
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
                sphere_parameter_no = simulation_parameter%sphere_list(j)%sphere_parameter_index

                a_n = simulation_parameter%sphere_parameter_list(sphere_parameter_no)%a_n
                b_n = simulation_parameter%sphere_parameter_list(sphere_parameter_no)%b_n


                do ii = lbound(sphere_no_list, 1), ubound(sphere_no_list, 1)
                    l = sphere_no_list(ii)

                    first = 2 * (j - first_sphere) * counter + 1
                    if (j .eq. l) then
                        last = first + 2*counter - 1
                        vector_b(first:last) = vector_b(first:last) + vector_x(first:last)
                    else
                        d_0_j = simulation_parameter%sphere_list(j)%d_0_j
                        sphere_parameter_no = simulation_parameter%sphere_list(j)%sphere_parameter_index
                        n_range_j = simulation_parameter%sphere_parameter_list(sphere_parameter_no)%n_range

                        d_0_l = simulation_parameter%sphere_list(l)%d_0_j
                        sphere_parameter_no = simulation_parameter%sphere_list(l)%sphere_parameter_index
                        n_range_l = simulation_parameter%sphere_parameter_list(sphere_parameter_no)%n_range

                        n_range_j(1) = max(n_range_j(1), n_range(1))
                        n_range_j(2) = min(n_range_j(2), n_range(2))

                        n_range_l(1) = max(n_range_l(1), n_range(1))
                        n_range_l(2) = min(n_range_l(2), n_range(2))

                        !$  if (thread_first_run) then
                        !$    call fwig_thread_temp_init(4 * max(n_range_l(2), n_range_l(2)))     ! multi threaded
                        !$    thread_first_run = .false.
                        !$  endif

                        d_j_l = d_0_l - d_0_j
                        call lib_math_vector_spherical_harmonics_translation_coefficient(d_j_l, &
                                n_range_j, n_range_l, z_selector, &
                                a_nmnumu, b_nmnumu)

                        call init_list(buffer_1_nm, n_range_l(1), n_range_l(2)-n_range_l(1)+1, dcmplx(0,0))
                        call init_list(buffer_2_nm, n_range_l(1), n_range_l(2)-n_range_l(1)+1, dcmplx(0,0))

                        call init_list(buffer_b_1_nm, n_range(1), n_range(2)-n_range(1)+1, dcmplx(0,0))
                        call init_list(buffer_b_2_nm, n_range(1), n_range(2)-n_range(1)+1, dcmplx(0,0))

                        ! calculate vector b
                        call make_list(vector_x(first:first+counter-1), &
                                       n_range(1), n_range(2)-n_range(1)+1, &
                                       buffer_x_1_nm)
                        call remove_zeros(buffer_x_1_nm)

                        call make_list(vector_x(first+counter:first+2*counter-1), &
                                       n_range(1), n_range(2)-n_range(1)+1, &
                                       buffer_x_2_nm)
                        call remove_zeros(buffer_x_2_nm)

                        !$OMP PARALLEL DO PRIVATE(n, m) &
                        !$OMP  PRIVATE(buffer_1_nm, buffer_2_nm)
                        do n = n_range(1), n_range(2)
                            do m = -n, n
                                if (n_range(2) .le. n_range_j(2) &
                                    .and. n_range(1) .ge. n_range_j(1)) then
                                    buffer_1_nm = a_n%item(n) * a_nmnumu%item(n)%item(m)
                                    buffer_2_nm = a_n%item(n) * b_nmnumu%item(n)%item(m)

                                    buffer_b_1_nm%item(n)%item(m) = sum(buffer_1_nm * buffer_x_1_nm &
                                                                        + buffer_2_nm * buffer_x_2_nm)

                                    buffer_1_nm = b_n%item(n) * b_nmnumu%item(n)%item(m)
                                    buffer_2_nm = b_n%item(n) * a_nmnumu%item(n)%item(m)

                                    buffer_b_2_nm%item(n)%item(m) = sum(buffer_1_nm * buffer_x_1_nm &
                                                                        + buffer_2_nm * buffer_x_2_nm)
                                else
                                    buffer_b_1_nm%item(n)%item(m) = dcmplx(0,0)
                                    buffer_b_2_nm%item(n)%item(m) = dcmplx(0,0)
                                end if
                            end do
                        end do
                        !$OMP END PARALLEL DO

                        n_range = simulation_parameter%spherical_harmonics%n_range

                        call make_array(buffer_b_1_nm, buffer_array, n_range(1) , n_range(2) - n_range(1) + 1)
                        last = first + counter - 1
                        vector_b(first:last) = vector_b(first:last) + buffer_array

                        call make_array(buffer_b_2_nm, buffer_array, n_range(1) , n_range(2) - n_range(1) + 1)
                        first = last + 1
                        last = first + counter - 1
                        vector_b(first:last) = vector_b(first:last) + buffer_array
                    end if
                end do
            end do
            !$OMP END PARALLEL DO
        end subroutine lib_mie_ms_solver_calculate_vector_b_without_transposed_select

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
        subroutine lib_mie_ms_solver_calculate_vector_b_with_transposed(simulation_parameter, vector_x, &
                                                        vector_b, vector_b_t, &
                                                        calc_vector_b, calc_vector_b_t)
            implicit none
            ! dummy
            type(lib_mie_simulation_parameter_type), intent(in) :: simulation_parameter
            double complex, dimension(:), intent(in) :: vector_x

            double complex, dimension(:), allocatable, intent(inout) :: vector_b
            double complex, dimension(:), allocatable, intent(inout) :: vector_b_t

            logical, intent(in), optional :: calc_vector_b
            logical, intent(in), optional :: calc_vector_b_t

            ! auxiliary
            integer :: j
            integer :: l
            integer :: n
            integer :: m
            integer :: nu
            integer :: mu

            integer :: first
            integer :: last

            integer :: first_sphere
            integer :: last_sphere

            integer :: sphere_parameter_no

            integer :: counter
            integer :: counter_sum

            integer(kind=1) :: z_selector

            type(cartesian_coordinate_real_type) :: d_0_j
            type(cartesian_coordinate_real_type) :: d_0_l
            type(cartesian_coordinate_real_type) :: d_j_l
            integer(kind=4), dimension(2) :: n_range
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

            logical :: m_calc_vector_b
            logical :: m_calc_vector_b_t

            ! set std values
            if (present(calc_vector_b)) then
                m_calc_vector_b = calc_vector_b
            else
                m_calc_vector_b = .true.
            end if

            if (present(calc_vector_b_t)) then
                m_calc_vector_b_t = calc_vector_b_t
            else
                m_calc_vector_b_t = .false.
            end if

            first_sphere = lbound(simulation_parameter%sphere_list, 1)
            last_sphere = ubound(simulation_parameter%sphere_list, 1)

            z_selector = simulation_parameter%spherical_harmonics%z_selector_translation

            n_range = simulation_parameter%spherical_harmonics%n_range

            counter = (1 + n_range(2))**2 - n_range(1)**2
            counter_sum = 2 * (last_sphere - first_sphere + 1) * counter

            if (m_calc_vector_b) then
                if (allocated(vector_b)) then
                    if (size(vector_b, 1) .ne. counter_sum) then
                        deallocate(vector_b)
                        allocate(vector_b(counter_sum))
                    end if
                else
                    allocate(vector_b(counter_sum))
                end if
                vector_b = 0
            end if

            if (m_calc_vector_b_t) then
                if (allocated(vector_b_t)) then
                if (size(vector_b_t, 1) .ne. counter_sum) then
                    deallocate(vector_b_t)
                    allocate(vector_b_t(counter_sum))
                end if
            else
                allocate(vector_b_t(counter_sum))
            end if
                vector_b_t = 0
            end if

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
                sphere_parameter_no = simulation_parameter%sphere_list(j)%sphere_parameter_index

                a_n = simulation_parameter%sphere_parameter_list(sphere_parameter_no)%a_n
                b_n = simulation_parameter%sphere_parameter_list(sphere_parameter_no)%b_n


                do l = first_sphere, last_sphere


                    first = 2 * (j - first_sphere) * counter + 1
                    if (j .eq. l) then
                        last = first + 2*counter - 1
                        if (calc_vector_b_t) then
                            vector_b(first:last) = vector_b(first:last) + vector_x(first:last)
                        end if

                        if (calc_vector_b_t) then
                            vector_b_t(first:last) = vector_b_t(first:last) + vector_x(first:last)
                        end if
                    else
                        d_0_j = simulation_parameter%sphere_list(j)%d_0_j
                        sphere_parameter_no = simulation_parameter%sphere_list(j)%sphere_parameter_index
                        n_range_j = simulation_parameter%sphere_parameter_list(sphere_parameter_no)%n_range

                        d_0_l = simulation_parameter%sphere_list(l)%d_0_j
                        sphere_parameter_no = simulation_parameter%sphere_list(l)%sphere_parameter_index
                        n_range_l = simulation_parameter%sphere_parameter_list(sphere_parameter_no)%n_range

                        n_range_j(1) = max(n_range_j(1), n_range(1))
                        n_range_j(2) = min(n_range_j(2), n_range(2))

                        n_range_l(1) = max(n_range_l(1), n_range(1))
                        n_range_l(2) = min(n_range_l(2), n_range(2))

                        d_j_l = d_0_l - d_0_j
                        call lib_math_vector_spherical_harmonics_translation_coefficient(d_j_l, &
                                n_range_j, n_range_l, z_selector, &
                                a_nmnumu, b_nmnumu)

                        call init_list(buffer_1_nm, n_range_l(1), n_range_l(2)-n_range_l(1)+1, dcmplx(0,0))
                        call init_list(buffer_2_nm, n_range_l(1), n_range_l(2)-n_range_l(1)+1, dcmplx(0,0))

                        call init_list(buffer_b_1_nm, n_range_l(1), n_range_l(2)-n_range_l(1)+1, dcmplx(0,0))
                        call init_list(buffer_b_2_nm, n_range_l(1), n_range_l(2)-n_range_l(1)+1, dcmplx(0,0))


                        if (m_calc_vector_b) then
                            ! calculate vector b
                            call make_list(vector_x(first:first+counter-1), &
                                           n_range(1), n_range(2)-n_range(1)+1, &
                                           buffer_x_1_nm)
                            call remove_zeros(buffer_x_1_nm)

                            call make_list(vector_x(first+counter:first+2*counter-1), &
                                           n_range(1), n_range(2)-n_range(1)+1, &
                                           buffer_x_2_nm)
                            call remove_zeros(buffer_x_2_nm)

                            do n = n_range_j(1), n_range_j(2)
                                do m = -n, n
                                    buffer_1_nm = a_n%item(n) * a_nmnumu%item(n)%item(m)
                                    buffer_2_nm = a_n%item(n) * b_nmnumu%item(n)%item(m)

                                    buffer_b_1_nm%item(n)%item(m) = sum(buffer_1_nm * buffer_x_1_nm &
                                                                        + buffer_2_nm * buffer_x_2_nm)

                                    buffer_1_nm = b_n%item(n) * b_nmnumu%item(n)%item(m)
                                    buffer_2_nm = b_n%item(n) * a_nmnumu%item(n)%item(m)

                                    buffer_b_2_nm%item(n)%item(m) = sum(buffer_1_nm * buffer_x_1_nm &
                                                                        + buffer_2_nm * buffer_x_2_nm)
                                end do
                            end do

                            call make_array(buffer_b_1_nm, buffer_array, n_range(1) , n_range(2) - n_range(1) + 1)
                            last = first + counter - 1
                            vector_b(first:last) = vector_b(first:last) + buffer_array

                            call make_array(buffer_b_2_nm, buffer_array, n_range(1) , n_range(2) - n_range(1) + 1)
                            first = last + 1
                            last = first + counter - 1
                            vector_b(first:last) = vector_b(first:last) + buffer_array
                        end if

                        if (m_calc_vector_b_t) then
                            ! calculate vector_b_t
                            ! -> Matrix is transposed
                            first = 2 * (j - first_sphere) * counter + 1
                            last = first

                            call make_list(vector_x(first:first+counter-1), &
                                           n_range(1), n_range(2)-n_range(1)+1, &
                                           buffer_x_1_nm)

                            call make_list(vector_x(first+counter:first+2*counter-1), &
                                           n_range(1), n_range(2)-n_range(1)+1, &
                                           buffer_x_2_nm)
                            call init_list(buffer_b_1_nm, n_range(1), n_range(2)-n_range(1)+1, dcmplx(0,0))
                            call init_list(buffer_b_2_nm, n_range(1), n_range(2)-n_range(1)+1, dcmplx(0,0))

                            do n = n_range_j(1), n_range_j(2)
                                do m = -n, n
                                    do nu = n_range_l(1), n_range_l(2)
                                        do mu = -nu, nu
                                            buffer_b_1_nm%item(n)%item(m) = buffer_b_1_nm%item(n)%item(m) &
                                                + a_n%item(nu) * a_nmnumu%item(nu)%item(mu)%item(n)%item(m) &
                                                  * buffer_x_1_nm%item(nu)%item(mu)
                                            buffer_b_1_nm%item(n)%item(m) = buffer_b_1_nm%item(n)%item(m) &
                                                + a_n%item(nu) * b_nmnumu%item(nu)%item(mu)%item(n)%item(m) &
                                                  * buffer_x_2_nm%item(nu)%item(mu)

                                            buffer_b_2_nm%item(n)%item(m) = buffer_b_2_nm%item(n)%item(m) &
                                                + a_n%item(nu) * b_nmnumu%item(nu)%item(mu)%item(n)%item(m) &
                                                  * buffer_x_1_nm%item(nu)%item(mu)
                                            buffer_b_2_nm%item(n)%item(m) = buffer_b_2_nm%item(n)%item(m) &
                                                + b_n%item(nu) * a_nmnumu%item(nu)%item(mu)%item(n)%item(m) &
                                                  * buffer_x_2_nm%item(nu)%item(mu)
                                        end do
                                    end do
                                end do
                            end do

                            ! insert buffer_b into vector_b_t
                            first = 2 * (l - first_sphere) * counter + 1
                            call make_array(buffer_b_1_nm, buffer_array, n_range(1) , n_range(2) - n_range(1) + 1)
                            last = first + counter - 1
                            vector_b_t(first:last) =vector_b_t(first:last) + buffer_array

                            call make_array(buffer_b_2_nm, buffer_array, n_range(1) , n_range(2) - n_range(1) + 1)
                            first = last + 1
                            last = first + counter - 1
                            vector_b_t(first:last) = vector_b_t(first:last) + buffer_array
                        end if

                    end if
                end do
            end do
!            !$OMP END PARALLEL DO
        end subroutine lib_mie_ms_solver_calculate_vector_b_with_transposed

         ! todo: optimise for ML FMM
         !
!        ! Formats problem of  multi sphere scattering to be able to use a solver.
!        !
!        ! Formula: A x = b
!        !   x: scattering coefficients
!        !   b: illumination coefficients
!        !   A: mixture of Mie coefficients and vector spherical translation coefficients
!        !
!        !   -> apply eq. 23 [1] to eq. 30 [2]
!        !
!        ! Argument
!        ! ----
!        !   simulation_parameter: type(lib_mie_simulation_parameter_type)
!        !   vector_x:
!        !
!        ! Returns
!        ! ----
!        !   vector_b: double complex, dimension(:)
!        !       result of the matrix vector multiplication
!        !
!        ! Reference: [1] Computation of scattering from clusters of spheres using the fast multipole method, Nail A. Gumerov, and Ramani Duraiswami
!        !            [2] Electromagnetic scatteringby an aggregate of spheres, Yu-lin Xu
!        subroutine lib_mie_ms_solver_calculate_vector_b(simulation_parameter, vector_x, &
!                                                        vector_b, vector_b_t, &
!                                                        calc_vector_b, calc_vector_b_t)
!            implicit none
!            ! dummy
!            type(lib_mie_simulation_parameter_type), intent(inout) :: simulation_parameter
!

!        subroutine lib_mie_ms_solver_calculate_vector_b_elements()
!            implicit none
!            ! dummy
!
!        end subroutine

        function lib_mie_ms_solver_interface_test_functions() result(rv)
            implicit none
            ! dummy
            integer :: rv

            ! auxiliaray
            double precision, parameter :: ground_truth_e = 10.0_8**(-12.0_8)

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

                    double complex, dimension(:), allocatable :: ground_truth_vector_b
                    double complex, dimension(:), allocatable :: ground_truth_vector_b_t

                    integer :: n
                    integer :: m
                    integer :: nu
                    integer :: mu

                    integer :: first
                    integer :: last
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
                    type(cartesian_coordinate_real_type) :: d_0_l
                    type(cartesian_coordinate_real_type) :: d_j_l

                    type(list_cmplx) :: a_n
                    type(list_cmplx) :: b_n
                    type(list_list_cmplx) :: a_numu
                    type(list_list_cmplx) :: b_numu
                    type(list_4_cmplx) :: a_nmnumu
                    type(list_4_cmplx) :: b_nmnumu

                    type(list_list_cmplx) :: buffer_vector_b_1
                    type(list_list_cmplx) :: buffer_vector_b_2

                    double complex, dimension(:), allocatable :: buffer_array

                    double complex :: buffer_cmplx
                    double precision :: buffer

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
                    n_range = (/ 1, 5 /)

                    sphere_para = lib_mie_type_func_get_sphere_parameter(lambda_0, n_medium, &
                                                                         r_particle, n_particle, n_range)

                    allocate(simulation%sphere_parameter_list(2))
                    sphere_para%a_n%item = dcmplx(0,0)
                    sphere_para%b_n%item = dcmplx(0,0)
                    simulation%sphere_parameter_list(1) = sphere_para
                    sphere_para%a_n%item = dcmplx(1,0)
                    sphere_para%b_n%item = dcmplx(1,0)
                    simulation%sphere_parameter_list(2) = sphere_para


                    allocate(simulation%sphere_list(2))
                    simulation%sphere_list(:)%sphere_parameter_index = 1
                    simulation%sphere_list(1)%sphere_parameter_index = 2

                    d_0_j = make_cartesian(-2d0,0d0,1d0) * unit_mu
                    simulation%sphere_list(1)%d_0_j = d_0_j

                    d_0_j = make_cartesian(2d0,0d0,1d0) * unit_mu
                    simulation%sphere_list(2)%d_0_j = d_0_j

!                    d_0_j = make_cartesian(0d0,0d0,3d0) * unit_mu
!                    simulation%sphere_list(3)%d_0_j = d_0_j

                    simulation%spherical_harmonics%z_selector_incident_wave = 1
                    simulation%spherical_harmonics%z_selector_scatterd_wave = 3
                    simulation%spherical_harmonics%z_selector_translation = 1
                    n_range = (/ 1, 5 /)
                    simulation%spherical_harmonics%n_range = n_range

                    no_of_elements = 2 * ( (1+n_range(2))**2 - n_range(1)**2 )
                    allocate(vector_x(size(simulation%sphere_list) * no_of_elements))

                    do i = 1, size(simulation%sphere_list)
                        first = (i - 1) * no_of_elements + 1
                        last = i * no_of_elements
                        vector_x(first:last) = i
                    end do

                    call lib_mie_ms_solver_calculate_vector_b(simulation, vector_x, &
                                                              vector_b, vector_b_t, &
                                                              calc_vector_b = .true., &
                                                              calc_vector_b_t = .true.)

                    ! generate ground truth
                    ! > vector_b
                    allocate(ground_truth_vector_b(size(simulation%sphere_list) * no_of_elements))
                    allocate(ground_truth_vector_b_t(size(simulation%sphere_list) * no_of_elements))

                    do i = 1, size(simulation%sphere_list)
                        first = (i - 1) * no_of_elements + 1
                        last = i * no_of_elements
                        ground_truth_vector_b(first:last) = i
                        ground_truth_vector_b_t(first:last) = i
                    end do

                    d_0_j = simulation%sphere_list(1)%d_0_j
                    d_0_l = simulation%sphere_list(2)%d_0_j
                    d_j_l = d_0_l - d_0_j

                    call lib_math_vector_spherical_harmonics_translation_coefficient(d_j_l, &
                            n_range, n_range, simulation%spherical_harmonics%z_selector_translation, &
                            a_nmnumu, b_nmnumu)
                    i = simulation%sphere_list(1)%sphere_parameter_index
                    a_n = simulation%sphere_parameter_list(i)%a_n
                    b_n = simulation%sphere_parameter_list(i)%b_n

                    call init_list(buffer_vector_b_1, n_range(1), n_range(2) - n_range(1) + 1)
                    call init_list(buffer_vector_b_2, n_range(1), n_range(2) - n_range(1) + 1)

                    do n = n_range(1), n_range(2)
                        do m = -n, n
                            a_numu = a_nmnumu%item(n)%item(m)
                            b_numu = b_nmnumu%item(n)%item(m)
                            buffer_vector_b_1%item(n)%item(m) = sum(a_n * a_numu + a_n * b_numu)
                            buffer_vector_b_2%item(n)%item(m) = sum(b_n * b_numu + b_n * a_numu)
                        end do
                    end do

                    call make_array(buffer_vector_b_1, buffer_array)
                    ground_truth_vector_b(1:no_of_elements/2) = ground_truth_vector_b(1:no_of_elements/2) &
                                                                + buffer_array

                    call make_array(buffer_vector_b_2, buffer_array)
                    ground_truth_vector_b(no_of_elements/2+1:no_of_elements) = &
                            ground_truth_vector_b(no_of_elements/2+1:no_of_elements) + buffer_array

                    ! > vector_b_t
                    call init_list(buffer_vector_b_1, n_range(1), n_range(2) - n_range(1) + 1, dcmplx(0,0))
                    call init_list(buffer_vector_b_2, n_range(1), n_range(2) - n_range(1) + 1, dcmplx(0,0))

                    do n = n_range(1), n_range(2)
                        do m = -n, n
                            do nu = n_range(1), n_range(2)
                                do mu = -nu, nu
                                    buffer_vector_b_1%item(n)%item(m) = buffer_vector_b_1%item(n)%item(m) &
                                        + a_n%item(nu) * a_nmnumu%item(nu)%item(mu)%item(n)%item(m) &
                                        + b_n%item(nu) * b_nmnumu%item(nu)%item(mu)%item(n)%item(m)

                                    buffer_vector_b_2%item(n)%item(m) = buffer_vector_b_2%item(n)%item(m) &
                                        + a_n%item(nu) * b_nmnumu%item(nu)%item(mu)%item(n)%item(m) &
                                        + b_n%item(nu) * a_nmnumu%item(nu)%item(mu)%item(n)%item(m)

                                end do
                            end do
                        end do
                    end do

                    call make_array(buffer_vector_b_1, buffer_array)
                    first = no_of_elements + 1
                    last = first + no_of_elements / 2 - 1
                    ground_truth_vector_b_t(first:last) = ground_truth_vector_b_t(first:last) &
                                                            + buffer_array

                    call make_array(buffer_vector_b_2, buffer_array)
                    first = last + 1
                    last = first + no_of_elements / 2 - 1
                    ground_truth_vector_b_t(first:last) = ground_truth_vector_b_t(first:last) &
                                                            +  buffer_array



                    rv = .true.
                    print *, "test_lib_mie_ms_solver_calculate_vector_b:"
                    do i=lbound(vector_b, 1), ubound(vector_b, 1)
                        buffer = abs(vector_b(i) - ground_truth_vector_b(i))
                        if (buffer .gt. ground_truth_e) then
                            print *, "  vector_b: ", i , "difference: ", buffer, " : FAILED"
                            rv = .false.
                        else
                            print *, "  vector_b: ", i, ": OK"
                        end if
                    end do

                    do i=lbound(vector_b_t, 1), ubound(vector_b_t, 1)
                        buffer = abs(vector_b_t(i) - ground_truth_vector_b_t(i))
                        if (buffer .gt. ground_truth_e) then
                            print *, "  vector_b_t: ", i , "difference: ", buffer, " : FAILED"
                            rv = .false.
                        else
                            print *, "  vector_b_t: ", i, ": OK"
                        end if
                    end do

                end function test_lib_mie_ms_solver_calculate_vector_b


        end function lib_mie_ms_solver_interface_test_functions
end module lib_mie_ms_solver_interface
