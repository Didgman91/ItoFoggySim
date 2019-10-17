module lib_mie_ms_ml_fmm_interface
    use libmath

    use lib_tree_public
    use lib_ml_fmm
    use ml_fmm_type

    use lib_mie_type
    use lib_mie_type_functions

    implicit none

    private

    integer :: m_first_sphere
    integer :: m_last_sphere

    contains

        subroutine lib_mie_ms_ml_fmm_constructor()
            use lib_mie_ms_data_container
            implicit none
            ! dummy

            ! auxiliary
            integer :: i
            integer :: no
            type(lib_ml_fmm_data) :: data_elements

            allocate(data_elements%XY(size(simulation_data%sphere_list)))

            m_first_sphere = lbound(simulation_data%sphere_list, 1)
            m_last_sphere = ubound(simulation_data%sphere_list, 1)
            do i = m_first_sphere, m_last_sphere
                no = i - m_first_sphere + 1

                data_elements%XY(no)%element_type = 1 ! todo: is this important for the ML FMM?
                data_elements%XY(no)%hierarchy = HIERARCHY_XY

                data_elements%XY(no)%point_x%x(1) = simulation_data%sphere_list(i)%d_0_j%x
                data_elements%XY(no)%point_x%x(2) = simulation_data%sphere_list(i)%d_0_j%y
                data_elements%XY(no)%point_x%x(3) = simulation_data%sphere_list(i)%d_0_j%z
            end do

            call lib_ml_fmm_constructor(data_elements)

        end subroutine lib_mie_ms_ml_fmm_constructor

        subroutine lib_mie_ms_ml_fmm_destructor
            implicit none

            call lib_ml_fmm_destructor()

        end subroutine lib_mie_ms_ml_fmm_destructor


        function lib_mie_ms_ml_fmm_get_u_B_i(x, data_element, element_number) result(u_B_i)
            use libmath
            use lib_tree_public
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

            i = element_number - m_first_sphere + 1
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
                u_B_i%a_nm = simulation_data%sphere_list(i)%a_nm
                u_B_i%b_nm = simulation_data%sphere_list(i)%b_nm

                return

            else if (abs(d_j_l) .gt. simulation_data%sphere_parameter_list(no)%radius) then
                z_selector = simulation_data%spherical_harmonics%z_selector_translation_gt_r
            else
                z_selector = simulation_data%spherical_harmonics%z_selector_translation_le_r
            end if

            call lib_math_vector_spherical_harmonics_translation_coefficient(d_j_l, &
                    n_range, n_range_j, z_selector, &
                    a_nmnumu, b_nmnumu)

            buffer_x_1_nm = simulation_data%sphere_list(i)%a_nm
            buffer_x_2_nm = simulation_data%sphere_list(i)%b_nm

            !$OMP PARALLEL DO PRIVATE(n, m) &
            !$OMP  PRIVATE(buffer_1_nm, buffer_2_nm)
            do n = n_range(1), n_range(2)
                do m = -n, n
                    if (n_range(2) .le. n_range_j(2) &
                        .and. n_range(1) .ge. n_range_j(1)) then
                        u_B_i%a_nm%item(n)%item(m) = sum(a_nmnumu%item(n)%item(m) * buffer_x_1_nm &
                                                         + b_nmnumu%item(n)%item(m) * buffer_x_2_nm)

                        u_B_i%b_nm%item(n)%item(m) = sum(b_nmnumu%item(n)%item(m) * buffer_x_1_nm &
                                                         + a_nmnumu%item(n)%item(m) * buffer_x_2_nm)
                    else
                        u_B_i%a_nm%item(n)%item(m) = dcmplx(0,0)
                        u_B_i%b_nm%item(n)%item(m) = dcmplx(0,0)
                    end if
                end do
            end do
            !$OMP END PARALLEL DO

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

end module lib_mie_ms_ml_fmm_interface
