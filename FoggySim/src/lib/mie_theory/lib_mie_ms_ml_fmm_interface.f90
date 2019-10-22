module lib_mie_ms_ml_fmm_interface
    use libmath

    use lib_tree_public
    use ml_fmm_type
    use lib_ml_fmm

    use lib_mie_type
    use lib_mie_type_functions
    use lib_mie_ms_solver_interface_helper_functions

    implicit none

    private

    public :: lib_mie_ms_ml_fmm_constructor
    public :: lib_mie_ms_ml_fmm_destructor
    public :: lib_mie_ms_ml_fmm_calculate_vector_b

    contains

        ! HINT: simulation_data @ lib_mie_ms_data_container has to be initialised
        !       - sphere_list
        subroutine lib_mie_ms_ml_fmm_constructor()
            use lib_ml_fmm_type_operator
            use ml_fmm_math
            use lib_mie_ms_data_container
            implicit none
            ! dummy

            ! auxiliary
            integer :: i
            integer :: no

            type(ml_fmm_type_operator_procedures) :: operator_procedures
            type(lib_ml_fmm_procedure_handles) :: ml_fmm_procedures
            type(lib_ml_fmm_data) :: data_elements

            allocate(data_elements%XY(size(simulation_data%sphere_list)))

            do i = lbound(simulation_data%sphere_list, 1), ubound(simulation_data%sphere_list, 1)
                no = i - lbound(simulation_data%sphere_list, 1) + 1

                data_elements%XY(no)%element_type = 1 ! todo: is this important for the ML FMM?
                data_elements%XY(no)%hierarchy = HIERARCHY_XY

                data_elements%XY(no)%point_x%x(1) = simulation_data%sphere_list(i)%d_0_j%x
                data_elements%XY(no)%point_x%x(2) = simulation_data%sphere_list(i)%d_0_j%y
                data_elements%XY(no)%point_x%x(3) = simulation_data%sphere_list(i)%d_0_j%z
            end do

            operator_procedures = ml_fmm_type_operator_get_procedures(1)
            ml_fmm_procedures = lib_mie_ms_ml_fmm_get_procedures()

            call lib_ml_fmm_constructor(data_elements, operator_procedures, ml_fmm_procedures)

        end subroutine lib_mie_ms_ml_fmm_constructor

        subroutine lib_mie_ms_ml_fmm_destructor
            implicit none

            call lib_ml_fmm_destructor()

        end subroutine lib_mie_ms_ml_fmm_destructor

        subroutine lib_mie_ms_ml_fmm_calculate_vector_b(vector_x, vector_b)
            use lib_mie_ms_data_container
            implicit none
            ! dummy
            double complex, dimension(:), intent(in) :: vector_x

            double complex, dimension(:), allocatable, intent(inout) :: vector_b

            ! auxiliary
            integer :: i
            integer :: j

            integer :: counter

            integer :: first_sphere
            integer :: last_sphere

            integer(kind=4), dimension(2) :: n_range

            type(list_list_cmplx) :: a_nm
            type(list_list_cmplx) :: b_nm

            type(lib_ml_fmm_v), dimension(:), allocatable :: b
            type(lib_ml_fmm_v), dimension(:), allocatable :: x

            first_sphere = lbound(simulation_data%sphere_list, 1)
            last_sphere = ubound(simulation_data%sphere_list, 1)

            n_range = simulation_data%spherical_harmonics%n_range

            allocate ( x(size(simulation_data%sphere_list)) )

            !$OMP PARALLEL DO PRIVATE(j, i, a_nm, b_nm)
            do j = first_sphere, last_sphere
                i = j - first_sphere + 1

                call lib_mie_ms_solver_hf_get_list_list_cmplx_from_array(vector_x, i, n_range, &
                                                                         a_nm, b_nm)

                call move_alloc(a_nm%item, x(i)%a_nm%item)
                call move_alloc(b_nm%item, x(i)%b_nm%item)
            end do
            !$OMP END PARALLEL DO

            call lib_ml_fmm_run(x, b)

            counter = (1 + n_range(2))**2 - n_range(1)**2
            if (allocated(vector_b)) deallocate(vector_b)
            allocate(vector_b(2 * counter * size(simulation_data%sphere_list)))

            !$OMP PARALLEL DO PRIVATE(i, a_nm, b_nm)
            do i = 1, size(simulation_data%sphere_list)
                a_nm = b(i)%a_nm
                b_nm = b(i)%b_nm
                call lib_mie_ms_solver_hf_insert_list_list_cmplx_into_array(a_nm, b_nm, i, n_range, vector_b)
            end do
            !$OMP END PARALLEL DO

        end subroutine lib_mie_ms_ml_fmm_calculate_vector_b

        function lib_mie_ms_ml_fmm_get_procedures() result (handle)
!            use ml_fmm_type
            use lib_ml_fmm_type_operator
            implicit none
            ! dummy
            type(lib_ml_fmm_procedure_handles) :: handle

            handle%get_u_phi_i_j => lib_mie_ms_ml_fmm_get_u_phi_i_j
            handle%get_u_B_i => lib_mie_ms_ml_fmm_get_u_B_i
            handle%get_translation_SS => lib_mie_ms_ml_fmm_translation_SS
            handle%get_translation_SR => lib_mie_ms_ml_fmm_translation_SR
            handle%get_translation_RR => lib_mie_ms_ml_fmm_translation_RR
            handle%dor => lib_mie_ms_ml_fmm_dor

        end function lib_mie_ms_ml_fmm_get_procedures

        function lib_mie_ms_ml_fmm_get_u_B_i(x, data_element, element_number) result(u_B_i)
            use libmath
            use lib_tree_public
            use lib_ml_fmm_data_container
            use lib_mie_ms_data_container
            implicit none
            ! dummy
            type(lib_tree_spatial_point), intent(in) :: x
            type(lib_tree_data_element), intent(in) :: data_element
            integer(kind=4), intent(in) :: element_number

            type(lib_ml_fmm_coefficient) :: u_B_i

            ! auxiliary
            integer :: i
            integer :: no

            integer :: n
            integer :: m

            integer(kind=1) :: z_selector

            TYPE(lib_tree_spatial_point) :: x_unscaled
            type(cartesian_coordinate_real_type) :: d_0_j
            type(cartesian_coordinate_real_type) :: d_0_l
            type(cartesian_coordinate_real_type) :: d_j_l
            integer(kind=4), dimension(2) :: n_range
            integer(kind=4), dimension(2) :: n_range_j

            type(list_list_cmplx) :: buffer_1_nm
            type(list_list_cmplx) :: buffer_2_nm

            type(list_list_cmplx) :: buffer_x_1_nm
            type(list_list_cmplx) :: buffer_x_2_nm

            type(list_4_cmplx) :: a_nmnumu
            type(list_4_cmplx) :: b_nmnumu

            i = element_number - lbound(simulation_data%sphere_list, 1) + 1
            d_0_j = simulation_data%sphere_list(i)%d_0_j
            no = simulation_data%sphere_list(i)%sphere_parameter_index
            n_range_j = simulation_data%sphere_parameter_list(no)%n_range

            n_range = simulation_data%spherical_harmonics%n_range

            n_range_j(1) = max(n_range_j(1), n_range(1))
            n_range_j(2) = min(n_range_j(2), n_range(2))

            x_unscaled = lib_tree_get_unscaled_point(x)
            d_0_l = make_cartesian(x_unscaled%x(1), &
                                   x_unscaled%x(2), &
                                   x_unscaled%x(3))

            d_j_l = d_0_l - d_0_j

            if (abs(d_j_l) .eq. 0d0) then
!                u_B_i%a_nm = simulation_data%sphere_list(i)%a_nm
!                u_B_i%b_nm = simulation_data%sphere_list(i)%b_nm
                u_B_i%a_nm = m_ml_fmm_u(i)%a_nm
                u_B_i%b_nm = m_ml_fmm_u(i)%b_nm

                return

            else if (abs(d_j_l) .gt. simulation_data%sphere_parameter_list(no)%radius) then
                z_selector = simulation_data%spherical_harmonics%z_selector_translation_gt_r
            else
                z_selector = simulation_data%spherical_harmonics%z_selector_translation_le_r
            end if

            call lib_math_vector_spherical_harmonics_translation_coefficient(d_j_l, &
                    n_range, n_range_j, z_selector, &
                    a_nmnumu, b_nmnumu)

!            buffer_x_1_nm = simulation_data%sphere_list(i)%a_nm
!            buffer_x_2_nm = simulation_data%sphere_list(i)%b_nm

            ! vector u aka vector x (solver_ Mx = b)
            buffer_x_1_nm = m_ml_fmm_u(i)%a_nm
            buffer_x_2_nm = m_ml_fmm_u(i)%b_nm

            call init_list(buffer_1_nm, n_range(1), n_range(2) - n_range(1) + 1)
            call init_list(buffer_2_nm, n_range(1), n_range(2) - n_range(1) + 1)

            !$OMP PARALLEL DO PRIVATE(n, m)
            do n = n_range(1), n_range(2)
                do m = -n, n
                    if (n_range(2) .le. n_range_j(2) &
                        .and. n_range(1) .ge. n_range_j(1)) then
                        buffer_1_nm%item(n)%item(m) = sum(a_nmnumu%item(n)%item(m) * buffer_x_1_nm &
                                                         + b_nmnumu%item(n)%item(m) * buffer_x_2_nm)

                        buffer_2_nm%item(n)%item(m) = sum(b_nmnumu%item(n)%item(m) * buffer_x_1_nm &
                                                         + a_nmnumu%item(n)%item(m) * buffer_x_2_nm)
                    else
                        buffer_1_nm%item(n)%item(m) = dcmplx(0,0)
                        buffer_2_nm%item(n)%item(m) = dcmplx(0,0)
                    end if
                end do
            end do
            !$OMP END PARALLEL DO

            call move_alloc(buffer_1_nm%item, u_B_i%a_nm%item)
            call move_alloc(buffer_2_nm%item, u_B_i%b_nm%item)


        end function lib_mie_ms_ml_fmm_get_u_B_i

        function lib_mie_ms_ml_fmm_translation_SS(B_i_1, x_1, x_2) result(B_i_2)
            use lib_tree_public
            use ml_fmm_type
            use lib_mie_ms_data_container
            implicit none
            ! dummy
            type(lib_ml_fmm_coefficient), intent(in) :: B_i_1
            type(lib_tree_spatial_point), intent(in) :: x_1
            type(lib_tree_spatial_point), intent(in) :: x_2
            type(lib_ml_fmm_coefficient) :: B_i_2

            ! auxiliaray
            integer :: n
            integer :: m

            TYPE(lib_tree_spatial_point) :: x_unscaled
            type(cartesian_coordinate_real_type) :: d_0_j
            type(cartesian_coordinate_real_type) :: d_0_l
            type(cartesian_coordinate_real_type) :: d_j_l

            integer, dimension(2) :: n_range
            integer(kind=1) :: z_selector

            type(list_4_cmplx) :: a_nmnumu
            type(list_4_cmplx) :: b_nmnumu

            z_selector = simulation_data%spherical_harmonics%z_selector_translation_le_r
            n_range = simulation_data%spherical_harmonics%n_range


            x_unscaled = lib_tree_get_unscaled_point(x_1)
            d_0_l = make_cartesian(x_unscaled%x(1), &
                                   x_unscaled%x(2), &
                                   x_unscaled%x(3))

            x_unscaled = lib_tree_get_unscaled_point(x_2)
            d_0_j = make_cartesian(x_unscaled%x(1), &
                                   x_unscaled%x(2), &
                                   x_unscaled%x(3))

            d_j_l = d_0_l - d_0_j
            call lib_math_vector_spherical_harmonics_translation_coefficient(d_j_l, &
                    n_range, n_range, z_selector, &
                    a_nmnumu, b_nmnumu)

            !$OMP PARALLEL DO PRIVATE(n, m)
            do n = n_range(1), n_range(2)
                do m = -n, n
                    B_i_2%a_nm%item(n)%item(m) = sum(a_nmnumu%item(n)%item(m) * B_i_1%a_nm &
                                                     + b_nmnumu%item(n)%item(m) * B_i_1%b_nm)

                    B_i_2%b_nm%item(n)%item(m) = sum(b_nmnumu%item(n)%item(m) * B_i_1%a_nm &
                                                     + a_nmnumu%item(n)%item(m) * B_i_1%b_nm)
                end do
            end do
            !$OMP END PARALLEL DO

        end function lib_mie_ms_ml_fmm_translation_SS

        function lib_mie_ms_ml_fmm_translation_SR(B_i_1, x_1, x_2) result(A_i_2)
            use lib_tree_public
            use ml_fmm_type
            use lib_mie_ms_data_container
            implicit none
            ! dummy
            type(lib_ml_fmm_coefficient), intent(in) :: B_i_1
            type(lib_tree_spatial_point), intent(in) :: x_1
            type(lib_tree_spatial_point), intent(in) :: x_2
            type(lib_ml_fmm_coefficient) :: A_i_2

            ! auxiliaray
            integer :: n
            integer :: m

            TYPE(lib_tree_spatial_point) :: x_unscaled
            type(cartesian_coordinate_real_type) :: d_0_j
            type(cartesian_coordinate_real_type) :: d_0_l
            type(cartesian_coordinate_real_type) :: d_j_l

            integer, dimension(2) :: n_range
            integer(kind=1) :: z_selector

            type(list_4_cmplx) :: a_nmnumu
            type(list_4_cmplx) :: b_nmnumu

            z_selector = simulation_data%spherical_harmonics%z_selector_translation_gt_r
            n_range = simulation_data%spherical_harmonics%n_range


            x_unscaled = lib_tree_get_unscaled_point(x_1)
            d_0_l = make_cartesian(x_unscaled%x(1), &
                                   x_unscaled%x(2), &
                                   x_unscaled%x(3))

            x_unscaled = lib_tree_get_unscaled_point(x_2)
            d_0_j = make_cartesian(x_unscaled%x(1), &
                                   x_unscaled%x(2), &
                                   x_unscaled%x(3))

            d_j_l = d_0_l - d_0_j
            call lib_math_vector_spherical_harmonics_translation_coefficient(d_j_l, &
                    n_range, n_range, z_selector, &
                    a_nmnumu, b_nmnumu)

            call init_list(A_i_2%a_nm, n_range(1), n_range(2) - n_range(1) +1)
            call init_list(A_i_2%b_nm, n_range(1), n_range(2) - n_range(1) +1)

            !$OMP PARALLEL DO PRIVATE(n, m)
            do n = n_range(1), n_range(2)
                do m = -n, n
                    A_i_2%a_nm%item(n)%item(m) = sum(a_nmnumu%item(n)%item(m) * B_i_1%a_nm &
                                                     + b_nmnumu%item(n)%item(m) * B_i_1%b_nm)

                    A_i_2%b_nm%item(n)%item(m) = sum(b_nmnumu%item(n)%item(m) * B_i_1%a_nm &
                                                     + a_nmnumu%item(n)%item(m) * B_i_1%b_nm)
                end do
            end do
            !$OMP END PARALLEL DO

        end function lib_mie_ms_ml_fmm_translation_SR

        function lib_mie_ms_ml_fmm_translation_RR(A_i_1, x_1, x_2) result(A_i_2)
            use lib_tree_public
            use ml_fmm_type
            use lib_mie_ms_data_container
            implicit none
            ! dummy
            type(lib_ml_fmm_coefficient), intent(in) :: A_i_1
            type(lib_tree_spatial_point), intent(in) :: x_1
            type(lib_tree_spatial_point), intent(in) :: x_2
            type(lib_ml_fmm_coefficient) :: A_i_2

            ! auxiliaray
            integer :: n
            integer :: m

            TYPE(lib_tree_spatial_point) :: x_unscaled
            type(cartesian_coordinate_real_type) :: d_0_j
            type(cartesian_coordinate_real_type) :: d_0_l
            type(cartesian_coordinate_real_type) :: d_j_l

            integer, dimension(2) :: n_range
            integer(kind=1) :: z_selector

            type(list_4_cmplx) :: a_nmnumu
            type(list_4_cmplx) :: b_nmnumu

            z_selector = simulation_data%spherical_harmonics%z_selector_translation_le_r
            n_range = simulation_data%spherical_harmonics%n_range


            x_unscaled = lib_tree_get_unscaled_point(x_1)
            d_0_l = make_cartesian(x_unscaled%x(1), &
                                   x_unscaled%x(2), &
                                   x_unscaled%x(3))

            x_unscaled = lib_tree_get_unscaled_point(x_2)
            d_0_j = make_cartesian(x_unscaled%x(1), &
                                   x_unscaled%x(2), &
                                   x_unscaled%x(3))

            d_j_l = d_0_l - d_0_j
            call lib_math_vector_spherical_harmonics_translation_coefficient(d_j_l, &
                    n_range, n_range, z_selector, &
                    a_nmnumu, b_nmnumu)

            !$OMP PARALLEL DO PRIVATE(n, m)
            do n = n_range(1), n_range(2)
                do m = -n, n
                    A_i_2%a_nm%item(n)%item(m) = sum(a_nmnumu%item(n)%item(m) * A_i_1%a_nm &
                                                     + b_nmnumu%item(n)%item(m) * A_i_1%b_nm)

                    A_i_2%b_nm%item(n)%item(m) = sum(b_nmnumu%item(n)%item(m) * A_i_1%a_nm &
                                                     + a_nmnumu%item(n)%item(m) * A_i_1%b_nm)
                end do
            end do
            !$OMP END PARALLEL DO

        end function lib_mie_ms_ml_fmm_translation_RR

        ! Argument
        ! ----
        !   data_element_i: type(lib_tree_data_element)
        !       tree specific data of the i-th element
        !   element_number_i: integer
        !       element number according to the tree list (1:N)
        !   y_j: type(lib_tree_spatial_point)
        !       scaled coordinated of the j-th element
        !   element_number_j: integer
        !       element number according to the tree list (1:N)
        !
        ! Returns
        ! ----
        !   rv: type(lib_ml_fmm_v)
        !       one summand of the sum, eq. 38 [1]
        !
        ! Reference: [1] Data Structures, Optimal Choice of Parameters, and Complexity Results for
        !                Generalized Multilevel Fast Multipole Methods in d Dimensions,
        !                Nail Gumerov, Ramani Duraiswami, Eugene Borovikov
        function lib_mie_ms_ml_fmm_get_u_phi_i_j(data_element_i, element_number_i, y_j, element_number_j) result(rv)
            use lib_tree_public
            use ml_fmm_type
            use lib_ml_fmm_data_container
            use lib_mie_ms_data_container
            implicit none
            ! dummy
            type(lib_tree_data_element), intent(in) :: data_element_i
            integer(kind=4), intent(in) :: element_number_i
            type(lib_tree_spatial_point), intent(in) :: y_j
            integer(kind=4), intent(in) :: element_number_j
            type(lib_ml_fmm_v) :: rv

            ! auxiliary
            integer :: j
            integer :: l

            integer :: n
            integer :: m

            integer :: no

!            TYPE(lib_tree_spatial_point) :: x_unscaled
            type(cartesian_coordinate_real_type) :: d_0_j
            type(cartesian_coordinate_real_type) :: d_0_l
            type(cartesian_coordinate_real_type) :: d_j_l

            integer, dimension(2) :: n_range
            integer(kind=4), dimension(2) :: n_range_j
            integer(kind=4), dimension(2) :: n_range_l
            integer(kind=1) :: z_selector

            type(list_list_cmplx) :: buffer_1_nm
            type(list_list_cmplx) :: buffer_2_nm

            type(list_list_cmplx) :: buffer_b_1_nm
            type(list_list_cmplx) :: buffer_b_2_nm

            type(list_list_cmplx) :: buffer_x_1_nm
            type(list_list_cmplx) :: buffer_x_2_nm

            type(list_cmplx) :: a_n
            type(list_cmplx) :: b_n

            type(list_4_cmplx) :: a_nmnumu
            type(list_4_cmplx) :: b_nmnumu

            j = element_number_j - lbound(simulation_data%sphere_list, 1) + 1
            d_0_j = simulation_data%sphere_list(j)%d_0_j
            no = simulation_data%sphere_list(j)%sphere_parameter_index
            n_range_j = simulation_data%sphere_parameter_list(no)%n_range

            a_n = simulation_data%sphere_parameter_list(no)%a_n
            b_n = simulation_data%sphere_parameter_list(no)%b_n

            n_range = simulation_data%spherical_harmonics%n_range
            z_selector = simulation_data%spherical_harmonics%z_selector_translation_gt_r

            buffer_x_1_nm = m_ml_fmm_u(j)%a_nm
            buffer_x_2_nm = m_ml_fmm_u(j)%b_nm

            l = element_number_i - lbound(simulation_data%sphere_list, 1) + 1

            if (j .eq. l) then
                buffer_b_1_nm = buffer_x_1_nm
                buffer_b_2_nm = buffer_x_2_nm
            else
                d_0_l = simulation_data%sphere_list(l)%d_0_j
                no = simulation_data%sphere_list(l)%sphere_parameter_index
                n_range_l = simulation_data%sphere_parameter_list(no)%n_range

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

                call init_list(buffer_b_1_nm, n_range(1), n_range(2)-n_range(1)+1, dcmplx(0,0))
                call init_list(buffer_b_2_nm, n_range(1), n_range(2)-n_range(1)+1, dcmplx(0,0))

                ! calculate element j of vector v aka vector b (solver: Mx=b)

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
            end if

            rv%a_nm = buffer_b_1_nm
            rv%b_nm = buffer_b_2_nm

        end function lib_mie_ms_ml_fmm_get_u_phi_i_j

        function lib_mie_ms_ml_fmm_dor(D, x_c, y_j, element_number_j)  result (rv)
            use lib_tree_public
            use ml_fmm_type
            use lib_ml_fmm_data_container
            use lib_mie_ms_data_container
            implicit none
            ! dummy
            type(lib_ml_fmm_coefficient), intent(in) :: D
            type(lib_tree_spatial_point), intent(in) :: x_c
            type(lib_tree_spatial_point), intent(in) :: y_j
            integer(kind=4), intent(in) :: element_number_j
            type(lib_ml_fmm_v) :: rv

            ! auxiliary
            integer :: j

            integer :: n
            integer :: m

            integer :: no

            type(lib_tree_spatial_point) :: x_unscaled
            type(cartesian_coordinate_real_type) :: d_j_l

            integer, dimension(2) :: n_range
            integer(kind=4), dimension(2) :: n_range_j
            integer(kind=4), dimension(2) :: n_range_l
            integer(kind=1) :: z_selector

            type(list_list_cmplx) :: buffer_1_nm
            type(list_list_cmplx) :: buffer_2_nm

            type(list_list_cmplx) :: buffer_b_1_nm
            type(list_list_cmplx) :: buffer_b_2_nm

            type(list_cmplx) :: a_n
            type(list_cmplx) :: b_n

            type(list_4_cmplx) :: a_nmnumu
            type(list_4_cmplx) :: b_nmnumu

            x_unscaled = lib_tree_get_unscaled_point(y_j - x_c)

            d_j_l = make_cartesian(x_unscaled%x(1), &
                                   x_unscaled%x(2), &
                                   x_unscaled%x(3))

            j = element_number_j - lbound(simulation_data%sphere_list, 1) + 1
            no = simulation_data%sphere_list(j)%sphere_parameter_index
            n_range_j = simulation_data%sphere_parameter_list(no)%n_range

            a_n = simulation_data%sphere_parameter_list(no)%a_n
            b_n = simulation_data%sphere_parameter_list(no)%b_n

            n_range_l(1) = lbound(d%a_nm%item, 1)
            n_range_l(2) = ubound(d%a_nm%item, 1)

            n_range = simulation_data%spherical_harmonics%n_range

            z_selector = simulation_data%spherical_harmonics%z_selector_translation_le_r

            n_range_j(1) = max(n_range_j(1), n_range(1))
            n_range_j(2) = min(n_range_j(2), n_range(2))

            n_range_l(1) = max(n_range_l(1), n_range(1))
            n_range_l(2) = min(n_range_l(2), n_range(2))

            if (abs(d_j_l) .eq. 0d0) then
                buffer_b_1_nm = D%a_nm
                buffer_b_2_nm = D%b_nm
            else
                call lib_math_vector_spherical_harmonics_translation_coefficient(d_j_l, &
                        n_range_j, n_range_l, z_selector, &
                        a_nmnumu, b_nmnumu)

                call init_list(buffer_1_nm, n_range_l(1), n_range_l(2)-n_range_l(1)+1, dcmplx(0,0))
                call init_list(buffer_2_nm, n_range_l(1), n_range_l(2)-n_range_l(1)+1, dcmplx(0,0))

                call init_list(buffer_b_1_nm, n_range(1), n_range(2)-n_range(1)+1, dcmplx(0,0))
                call init_list(buffer_b_2_nm, n_range(1), n_range(2)-n_range(1)+1, dcmplx(0,0))

                ! calculate element j of vector v aka vector b (solver: Mx=b)

                !$OMP PARALLEL DO PRIVATE(n, m) &
                !$OMP  PRIVATE(buffer_1_nm, buffer_2_nm)
                do n = n_range(1), n_range(2)
                    do m = -n, n
                        if (n_range(2) .le. n_range_j(2) &
                            .and. n_range(1) .ge. n_range_j(1)) then
                            buffer_1_nm = a_n%item(n) * a_nmnumu%item(n)%item(m)
                            buffer_2_nm = a_n%item(n) * b_nmnumu%item(n)%item(m)

                            buffer_b_1_nm%item(n)%item(m) = sum(buffer_1_nm + buffer_2_nm)

                            buffer_1_nm = b_n%item(n) * b_nmnumu%item(n)%item(m)
                            buffer_2_nm = b_n%item(n) * a_nmnumu%item(n)%item(m)

                            buffer_b_2_nm%item(n)%item(m) = sum(buffer_1_nm + buffer_2_nm)
                        else
                            buffer_b_1_nm%item(n)%item(m) = dcmplx(0,0)
                            buffer_b_2_nm%item(n)%item(m) = dcmplx(0,0)
                        end if
                    end do
                end do
                !$OMP END PARALLEL DO
            end if

            rv%a_nm = buffer_b_1_nm
            rv%b_nm = buffer_b_2_nm

        end function lib_mie_ms_ml_fmm_dor

end module lib_mie_ms_ml_fmm_interface
