module lib_mie_ms_solver_interface_helper_functions
    use libmath
    implicit none

    private

    public :: lib_mie_ms_solver_hf_insert_list_list_cmplx_into_array
    public :: lib_mie_ms_solver_hf_add_list_list_cmplx_at_array
    public :: lib_mie_ms_solver_hf_get_list_list_cmplx_from_array

    contains

        ! Argument
        ! ----
        !   a_nm: type(list_list_cmplx)
        !       coefficient of the element
        !   b_nm: type(list_list_cmplx)
        !       coefficient of the element
        !   element_no: integer
        !       number of the element
        !   n_range: integer, dimension(2)
        !
        ! Returns
        ! ----
        !   array: double complex, dimension(:)
        !
        subroutine lib_mie_ms_solver_hf_insert_list_list_cmplx_into_array(a_nm, b_nm, element_no, n_range, array)
            implicit none
            ! dummy
            type(list_list_cmplx), intent(in) :: a_nm
            type(list_list_cmplx), intent(in) :: b_nm
            integer, intent(in) :: element_no
            integer, dimension(2), intent(in) :: n_range

            double complex, dimension(:), intent(inout) :: array

            ! auxiliary
            integer :: counter
            integer :: first
            integer :: last

            double complex, dimension(:), allocatable :: buffer_array

            counter = (1 + n_range(2))**2 - n_range(1)**2

            first = 2 * (element_no - 1) * counter + 1
            call make_array(a_nm, buffer_array, n_range(1) , n_range(2) - n_range(1) + 1)
            last = first + counter - 1
            array(first:last) = buffer_array

            call make_array(b_nm, buffer_array, n_range(1) , n_range(2) - n_range(1) + 1)
            first = last + 1
            last = first + counter - 1
            array(first:last) = buffer_array

        end subroutine lib_mie_ms_solver_hf_insert_list_list_cmplx_into_array

        ! Argument
        ! ----
        !   a_nm: type(list_list_cmplx)
        !       coefficient of the element
        !   b_nm: type(list_list_cmplx)
        !       coefficient of the element
        !   element_no: integer
        !       number of the element
        !   n_range: integer, dimension(2)
        !
        ! Returns
        ! ----
        !   array: double complex, dimension(:)
        !
        subroutine lib_mie_ms_solver_hf_add_list_list_cmplx_at_array(a_nm, b_nm, element_no, n_range, array)
            implicit none
            ! dummy
            type(list_list_cmplx), intent(in) :: a_nm
            type(list_list_cmplx), intent(in) :: b_nm
            integer, intent(in) :: element_no
            integer, dimension(2), intent(in) :: n_range

            double complex, dimension(:), intent(inout) :: array

            ! auxiliary
            integer :: counter
            integer :: first
            integer :: last

            double complex, dimension(:), allocatable :: buffer_array

            counter = (1 + n_range(2))**2 - n_range(1)**2

            first = 2 * (element_no - 1) * counter + 1
            call make_array(a_nm, buffer_array, n_range(1) , n_range(2) - n_range(1) + 1)
            last = first + counter - 1
            array(first:last) =array(first:last) + buffer_array

            call make_array(b_nm, buffer_array, n_range(1) , n_range(2) - n_range(1) + 1)
            first = last + 1
            last = first + counter - 1
            array(first:last) = array(first:last) + buffer_array

        end subroutine lib_mie_ms_solver_hf_add_list_list_cmplx_at_array

        ! Argument
        ! ----
        !   array: double complex, dimension(:)
        !       array of all a_nm and b_nm coefficients
        !   element_no: integer
        !       number of the element
        !   n_range: integer, dimension(2)
        !
        subroutine lib_mie_ms_solver_hf_get_list_list_cmplx_from_array(array, element_no, n_range, a_nm, b_nm)
            implicit none
            ! dummy
            double complex, dimension(:), intent(in) :: array
            integer, intent(in) :: element_no
            integer, dimension(2), intent(in) :: n_range

            type(list_list_cmplx), intent(inout) :: a_nm
            type(list_list_cmplx), intent(inout) :: b_nm

            ! auxiliary
            integer :: counter
            integer :: first
            integer :: last

            double complex, dimension(:), allocatable :: buffer_array

            counter = (1 + n_range(2))**2 - n_range(1)**2

            first = 2 * (element_no - 1) * counter + 1
            last = first + counter - 1
            call make_list(array(first:last), &
                           n_range(1), n_range(2)-n_range(1)+1, &
                           a_nm)
            call remove_zeros(a_nm)

            first = last + 1
            last = first + counter - 1
            call make_list(array(first:last), &
                           n_range(1), n_range(2)-n_range(1)+1, &
                           b_nm)
            call remove_zeros(b_nm)

        end subroutine lib_mie_ms_solver_hf_get_list_list_cmplx_from_array
end module lib_mie_ms_solver_interface_helper_functions
