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
            integer :: ii
            integer :: m_i
            integer :: m_ii
            integer :: m_i_start
            integer :: m_ii_start
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

            type(list_cmplx) :: a_n
            type(list_cmplx) :: b_n
            type(list_4_cmplx) :: a_nmnumu
            type(list_4_cmplx) :: b_nmnumu

            type(list_list_cmplx) :: b1_nm
            type(list_list_cmplx) :: b2_nm

            first_sphere = lbound(simulation_parameter%sphere_list, 1)
            last_sphere = ubound(simulation_parameter%sphere_list, 1)

            allocate(counter(first_sphere:last_sphere, first_sphere:last_sphere))

            z_selector = simulation_parameter%spherical_harmonics%z_selector_translation

            ! calculate elements of matrix_a
            counter = 0
!            !$OMP PARALLEL DO PRIVATE(i, ii, d_0_j, d_0_l, d_j_l, sphere_parameter_no, n_range_l, n_range_j)
            do i = first_sphere, last_sphere
                sphere_parameter_no = simulation_parameter%sphere_list(i)%sphere_parameter_index
                n_range_j = simulation_parameter%sphere_parameter_list(sphere_parameter_no)%n_range

                a_n = simulation_parameter%sphere_parameter_list(sphere_parameter_no)%a_n
                b_n = simulation_parameter%sphere_parameter_list(sphere_parameter_no)%b_n

                do ii = first_sphere, last_sphere

                    if (i .eq. ii) then

                    else
                        d_0_j = simulation_parameter%sphere_list(i)%d_0_j
                        sphere_parameter_no = simulation_parameter%sphere_list(i)%sphere_parameter_index
                        n_range_j = simulation_parameter%sphere_parameter_list(sphere_parameter_no)%n_range

                        d_0_l = simulation_parameter%sphere_list(ii)%d_0_j
                        sphere_parameter_no = simulation_parameter%sphere_list(ii)%sphere_parameter_index
                        n_range_l = simulation_parameter%sphere_parameter_list(sphere_parameter_no)%n_range

                        d_j_l = d_0_l - d_0_j
                        call lib_math_vector_spherical_harmonics_translation_coefficient(d_j_l, &
                                n_range_j, n_range_l, z_selector, &
                                a_nmnumu, b_nmnumu)

                        counter(i, ii) = (1 + n_range_j(2))**2 - n_range_j(1)**2
!                        counter(i, ii) = ( (1 + n_range_j(2))**2 - n_range_j(1)**2 ) &
!                                     * ( (1 + n_range_l(2))**2 - n_range_l(1)**2 )
                    end if
                end do
            end do
!            !$OMP END PARALLEL DO

            counter_sum = 2 * sum(counter)

!            if (allocated(matrix_a)) then
!                deallocate(matrix_a)
!            end if
!            allocate(matrix_a(counter_sum, counter_sum))
!
!            matrix_a = 0
!
!!            m_i = 1
!!            m_ii = 1
!            do i = first_sphere, last_sphere
!                do ii = first_sphere, last_sphere
!
!                    if (i .eq. ii) then
!                        sphere_parameter_no = simulation_parameter%sphere_list(i)%sphere_parameter_index
!                        n_range_j = simulation_parameter%sphere_parameter_list(sphere_parameter_no)%n_range
!
!                        m_i_start = sum(counter(1:i,ii)) - counter(i,ii) + 1
!                        m_ii_start = sum(counter(i,1:ii)) - counter(i,ii) + 1
!
!                        m_i = m_i_start
!                        m_ii = m_ii_start
!                        do n = n_range_j(1), n_range_j(2)
!                            do m = -n, n
!    !                            matrix_a(m_i, m_ii) = 1D0
!    !                            matrix_a(m_i+1, m_ii+1) = 1D0
!
!                                ! test
!                                matrix_a(m_i, m_ii) = a_n(i)%item(n)
!                                matrix_a(m_i+1, m_ii+1) = b_n(i)%item(n)
!                                ! ~~~ test ~~~
!
!                                m_i = m_i + 2
!                                m_ii = m_ii + 2
!                            end do
!                        end do
!                    else
!!                        d_0_j = simulation_parameter%sphere_list(i)%d_0_j
!!                        sphere_parameter_no = simulation_parameter%sphere_list(i)%sphere_parameter_index
!!                        n_range_j = simulation_parameter%sphere_parameter_list(sphere_parameter_no)%n_range
!!
!!                        d_0_l = simulation_parameter%sphere_list(i)%d_0_j
!!                        sphere_parameter_no = simulation_parameter%sphere_list(i)%sphere_parameter_index
!!                        n_range_l = simulation_parameter%sphere_parameter_list(sphere_parameter_no)%n_range
!!
!!                        m_i_start = sum(counter(1:i,ii)) - counter(i,ii) + 1
!!                        m_ii_start = sum(counter(i,1:ii)) - counter(i,ii) + 1
!!
!!                        m_i = m_i_start
!!                        do n = n_range_l(1), n_range_l(2)
!!                            do m = -n, n
!!                                m_ii = m_ii_start
!!                                do nu = n_range_l(1), n_range_l(2)
!!                                    do mu = -nu, nu
!!                                        matrix_a(m_i, m_ii) = a_nmnumu(i, ii)%item(n)%item(m)%item(nu)%item(mu)
!!                                        matrix_a(m_i, m_ii+1) = b_nmnumu(i, ii)%item(n)%item(m)%item(nu)%item(mu)
!!
!!                                        matrix_a(m_i+1, m_ii+1) = a_nmnumu(i, ii)%item(n)%item(m)%item(nu)%item(mu)
!!                                        matrix_a(m_i+1, m_ii) = b_nmnumu(i, ii)%item(n)%item(m)%item(nu)%item(mu)
!!
!!                                        m_ii = m_ii + 2
!!                                    end do
!!                                end do
!!                            end do
!!                            m_i = m_i + 2
!!                        end do
!                    end if
!                end do
!            end do

        end subroutine lib_mie_ms_solver_calculate_vector_b
end module lib_mie_ms_solver
