module lib_mie_ms_ml_fmm_interface
    use libmath

    use lib_tree_public
    use lib_ml_fmm
    use ml_fmm_type

    use lib_mie_type
    use lib_mie_type_functions
    implicit none

    private


    type(lib_mie_simulation_parameter_type) :: m_simulation_parameter

    integer :: m_first_sphere
    integer :: m_last_sphere

    contains

        subroutine lib_mie_ms_ml_fmm_constructor(simulation)
            implicit none
            ! dummy
            type(lib_mie_simulation_parameter_type), intent(inout) :: simulation

            ! auxiliary
            integer :: i
            integer :: no
            type(lib_ml_fmm_data) :: data_elements

            m_simulation_parameter = simulation

            allocate(data_elements%XY(size(m_simulation_parameter%sphere_list)))

            m_first_sphere = lbound(m_simulation_parameter%sphere_list, 1)
            m_last_sphere = ubound(m_simulation_parameter%sphere_list, 1)
            do i = m_first_sphere, m_last_sphere
                no = i - m_first_sphere + 1

                data_elements%XY(no)%element_type = 1 ! todo: is this important for the ML FMM?
                data_elements%XY(no)%hierarchy = HIERARCHY_XY

                data_elements%XY(no)%point_x%x(1) = m_simulation_parameter%sphere_list(i)%d_0_j%x
                data_elements%XY(no)%point_x%x(2) = m_simulation_parameter%sphere_list(i)%d_0_j%y
                data_elements%XY(no)%point_x%x(3) = m_simulation_parameter%sphere_list(i)%d_0_j%z
            end do

            call lib_ml_fmm_constructor(data_elements)

        end subroutine lib_mie_ms_ml_fmm_constructor

        subroutine lib_mie_ms_ml_fmm_destructor
            implicit none

            call lib_ml_fmm_destructor()

        end subroutine lib_mie_ms_ml_fmm_destructor


        function lib_mie_ms_ml_fmm_get_u_B_i(x, data_element, element_number) result(u_B_i)
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

            type(cartesian_coordinate_real_type) :: d_0_j
            type(cartesian_coordinate_real_type) :: d_0_l
            type(cartesian_coordinate_real_type) :: d_j_l
            integer(kind=4), dimension(2) :: n_range
            integer(kind=4), dimension(2) :: n_range_j

            type(list_cmplx) :: a_n
            type(list_cmplx) :: b_n

            type(list_4_cmplx) :: a_nmnumu
            type(list_4_cmplx) :: b_nmnumu

            i = element_number - m_first_sphere + 1
            no = m_simulation_parameter%sphere_list(i)%sphere_parameter_index
            n_range_j = m_simulation_parameter%sphere_parameter_list(no)%n_range

            n_range = m_simulation_parameter%spherical_harmonics%n_range

            n_range_j(1) = max(n_range_j(1), n_range(1))
            n_range_j(2) = min(n_range_j(2), n_range(2))

            d_j_l = d_0_l - d_0_j
            call lib_math_vector_spherical_harmonics_translation_coefficient(d_j_l, &
                    n_range, n_range_j, z_selector, &
                    a_nmnumu, b_nmnumu)

!            m_simulation_parameter%
    !        allocate(B_i%dummy, source = (/data_element%uindex%n + abs(x)/))
            allocate(u_B_i%dummy, source = (/real(data_element%uindex%n, kind=LIB_ML_FMM_COEFFICIENT_KIND)/))
        end function

end module lib_mie_ms_ml_fmm_interface
