module ml_fmm_math
    use ml_fmm_type
    use lib_ml_fmm_type_operator
    implicit none

    private

    ! --- public functions ---
    public :: ml_fmm_type_operator_get_procedures

    public :: ml_fmm_get_test_procedures
    public :: lib_ml_fmm_type_operator_test_functions

    contains

    ! Argument
    ! ----
    !   m_procedure_type: integer, optional (std=0)
    !       operators manipulating the vector v, u and the coefficients C, D
    !       0: dummy
    !       1: a_nm and b_nm: type(list_list_cmplx)
    function ml_fmm_type_operator_get_procedures(procedure_type) result (operator_procedures)
        implicit none
        ! dummy
        integer, intent(in), optional :: procedure_type
        type(ml_fmm_type_operator_procedures) :: operator_procedures

        ! auxiliary
        integer :: m_procedure_type

        m_procedure_type = 0
        if (present(procedure_type)) m_procedure_type = procedure_type

        if (m_procedure_type .eq. 0) then
            ! load test procedure functions
            operator_procedures%coefficient_add => test_c_add
            operator_procedures%coefficient_set_zero => test_set_coefficient_zero
!            operator_procedures%u_dot_coefficient => test_u_dot_coefficient
!            operator_procedures%coefficient_eq => test_coefficient_eq
!            operator_procedures%coefficient_ne => test_coefficient_ne
!            operator_procedures%v_add => test_v_operator_add
            operator_procedures%v_add_0D => test_v_operator_add_0D
!            operator_procedures%v_sub => test_v_operator_sub
!            operator_procedures%v_sub_0D => test_v_operator_sub_0D
        else if (m_procedure_type .eq. 1) then
            operator_procedures%coefficient_add => ml_fmm_coefficient_add_operator_list_2_cmplx
            operator_procedures%v_add_0D => ml_fmm_type_operator_v_add_0d_list_2_cmplx
            operator_procedures%coefficient_set_zero => lib_ml_fmm_type_operator_set_coefficient_zero_list_2_cmplx
        end if

    end function ml_fmm_type_operator_get_procedures

    function ml_fmm_get_test_procedures() result (handle)
        implicit none
        ! dummy
        type(lib_ml_fmm_procedure_handles) :: handle

        handle%get_u_B_i => test_get_u_B_i
        handle%get_u_phi_i_j => test_get_u_phi_i_j
        handle%get_translation_RR  => test_translation_RR
        handle%get_translation_SR  => test_translation_SR
        handle%get_translation_SS  => test_translation_SS
        handle%dor => test_dor
    end function

    function ml_fmm_coefficient_add_operator_list_2_cmplx(lhs,rhs) result (rv)
        use libmath
        use ml_fmm_type
        implicit none
        ! dummy
        type (lib_ml_fmm_coefficient), intent(in) :: lhs, rhs
        type (lib_ml_fmm_coefficient) :: rv

        rv%a_nm = lhs%a_nm + rhs%a_nm
        rv%b_nm = lhs%b_nm + rhs%b_nm

    end function ml_fmm_coefficient_add_operator_list_2_cmplx

    function ml_fmm_type_operator_v_add_0d_list_2_cmplx(lhs, rhs) result(rv)
        use libmath
        use ml_fmm_type
        implicit none
        ! dummy
        type (lib_ml_fmm_v), intent(in) :: lhs
        type (lib_ml_fmm_v), intent(in) :: rhs
        type (lib_ml_fmm_v) :: rv

        rv%a_nm = lhs%a_nm + rhs%a_nm
        rv%b_nm = lhs%b_nm + rhs%b_nm

    end function ml_fmm_type_operator_v_add_0d_list_2_cmplx

    subroutine lib_ml_fmm_type_operator_set_coefficient_zero_list_2_cmplx(coefficient)
        use ml_fmm_type
        implicit none
        ! dummy
        type(lib_ml_fmm_coefficient), intent(inout) :: coefficient

        call init_list(coefficient%a_nm, 1 ,1, dcmplx(0,0))
        call init_list(coefficient%b_nm, 1 ,1, dcmplx(0,0))

    end subroutine lib_ml_fmm_type_operator_set_coefficient_zero_list_2_cmplx

    ! ----- test functions ----
    function lib_ml_fmm_type_operator_test_functions() result (error_counter)
        implicit none
        integer :: error_counter

        error_counter = 0

        if (.not. test_lib_ml_fmm_type_operator_constructor()) then
            error_counter = error_counter + 1
        end if

        print *, "-------------lib_ml_fmm_type_operator_test_functions----------------"
        if (error_counter == 0) then
            print *, "lib_ml_fmm_type_operator_test_functions tests: OK"
        else
            print *, error_counter,"lib_ml_fmm_type_operator_test_functions test(s) FAILED"
        end if
        print *, "--------------------------------------------------------------------"

        contains

        function test_lib_ml_fmm_type_operator_constructor() result (rv)
            implicit none
            ! dummy
            logical :: rv

            ! auxiliary
            integer(kind=1), parameter :: length = 3
            type(lib_ml_fmm_v) :: lhs_u
            type(lib_ml_fmm_coefficient) :: lhs_coeff
            type(lib_ml_fmm_coefficient) :: rhs_coeff
            type(lib_ml_fmm_coefficient) :: res_coeff
            type(lib_tree_spatial_point) :: rhs_R
!            type(lib_ml_fmm_v) :: res_v

            type(ml_fmm_type_operator_procedures) :: ml_fmm_operator_procedures

            ml_fmm_operator_procedures%coefficient_add => test_c_add
            ml_fmm_operator_procedures%coefficient_set_zero => test_set_coefficient_zero
!            ml_fmm_operator_procedures%dor => test_dor
!            ml_fmm_operator_procedures%u_dot_coefficient => test_u_dot_coefficient
!            ml_fmm_operator_procedures%v_add => test_v_operator_add
            ml_fmm_operator_procedures%v_add_0D => test_v_operator_add_0D
!            ml_fmm_operator_procedures%v_sub => test_v_operator_sub
!            ml_fmm_operator_procedures%v_sub_0D => test_v_operator_sub_0D

            call lib_ml_fmm_type_operator_constructor(ml_fmm_operator_procedures)
!            call lib_ml_fmm_type_operator_constructor(test_c_add, test_u_dot_coefficient, test_cor, &
!                                                      test_set_coefficient_zero)!, &
!!                                                      test_allocate_coefficient_list, &
!!                                                      test_set_coefficient, test_get_coefficient, &
!!                                                      test_deallocate_coefficient_list)

            allocate (lhs_u%dummy(length))
            allocate (lhs_coeff%dummy(length))
            allocate (rhs_coeff%dummy(length))

            lhs_coeff%dummy = (/5, 4, 3/)
            call lib_ml_fmm_type_operator_set_coefficient_zero(lhs_coeff)
            lhs_coeff%dummy = (/5, 4, 3/)
            rhs_coeff%dummy = (/2, 1, 3/)

            lhs_u%dummy = (/1, 2, 3/)

            res_coeff = lhs_coeff + rhs_coeff

!            res_coeff = lhs_u * res_coeff

            rhs_R%x(:2) = (/1.5, 7.0/)
!            res_v = lhs_coeff .dor. rhs_R

            ! todo: add condition
            rv = .true.

        end function test_lib_ml_fmm_type_operator_constructor

    end function lib_ml_fmm_type_operator_test_functions

    function test_c_add(lhs, rhs) result (rv)
        implicit none
        ! dummy
        type (lib_ml_fmm_coefficient), intent(in) :: lhs, rhs
        type (lib_ml_fmm_coefficient) :: rv

        ! auxiliary
        integer :: i
        integer :: length

        if (size(lhs%dummy) .eq. size(rhs%dummy)) then
            length = size(lhs%dummy)
            allocate (rv%dummy(length))
            do i=1, length
                rv%dummy(i) = lhs%dummy(i) + rhs%dummy(i)
            end do
        end if

    end function

    function test_u_dot_coefficient(lhs, rhs) result (rv)
        implicit none
        ! dummy
        type (lib_ml_fmm_v), intent(in) :: lhs
        type (lib_ml_fmm_coefficient), intent(in) :: rhs
        type (lib_ml_fmm_coefficient) :: rv

        ! auxiliary
        integer :: i
        integer :: length

        if (size(lhs%dummy) .eq. size(rhs%dummy)) then
            length = size(lhs%dummy)
            allocate (rv%dummy(length))
            do i=1, length
                rv%dummy(i) = lhs%dummy(i) * rhs%dummy(i)
            end do
        end if
    end function test_u_dot_coefficient

    function test_dor(D, x_c, y_j, element_number_j)  result (rv)
        implicit none
        ! dummy
        type(lib_ml_fmm_coefficient), intent(in) :: D
        type(lib_tree_spatial_point), intent(in) :: x_c
        type(lib_tree_spatial_point), intent(in) :: y_j
        integer(kind=4), intent(in) :: element_number_j
        type(lib_ml_fmm_v) :: rv

        ! auxiliary
        integer :: i
        integer :: length

        allocate (rv%dummy(1))
        rv%dummy(1) = 0

        length = size(D%dummy)
        do i=1, length
            rv%dummy(1) = rv%dummy(1) + D%dummy(i) * (y_j%x(1) - x_c%x(1))
        end do

    end function test_dor

    function test_coefficient_eq(lhs, rhs) result (rv)
        implicit none
        ! dummy
        type(lib_ml_fmm_coefficient), intent(in) :: lhs
        type(lib_ml_fmm_coefficient), intent(in) :: rhs
        logical :: rv

        if (lhs%dummy(1) .eq. rhs%dummy(1)) then
            rv = .true.
        else
            rv = .false.
        end if

    end function test_coefficient_eq

    function test_coefficient_ne(lhs, rhs) result (rv)
        implicit none
        ! dummy
        type(lib_ml_fmm_coefficient), intent(in) :: lhs
        type(lib_ml_fmm_coefficient), intent(in) :: rhs
        logical :: rv

        if (lhs%dummy(1) .ne. rhs%dummy(1)) then
            rv = .true.
        else
            rv = .false.
        end if

    end function test_coefficient_ne

    subroutine test_set_coefficient_zero(coefficient)
        implicit none
        ! dummy
        type(lib_ml_fmm_coefficient), intent(inout) :: coefficient

        if (allocated(coefficient%dummy)) then
            coefficient%dummy(:) = 0
        else
            allocate(coefficient%dummy(1))
            coefficient%dummy(1) = 0
        end if
    end subroutine

    function test_v_operator_add(lhs, rhs) result(rv)
        use ml_fmm_type
        implicit none
        ! dummy
        type (lib_ml_fmm_v), dimension(:), intent(in) :: lhs
        type (lib_ml_fmm_v), dimension(size(lhs)), intent(in) :: rhs
        type (lib_ml_fmm_v), dimension(size(lhs)) :: rv

        ! auxilary
        integer :: i

        do i=1, size(lhs)
            rv(i)%dummy = lhs(i)%dummy + rhs(i)%dummy
        end do
    end function

    function test_v_operator_add_0d(lhs, rhs) result(rv)
        use ml_fmm_type
        implicit none
        ! dummy
        type (lib_ml_fmm_v), intent(in) :: lhs
        type (lib_ml_fmm_v), intent(in) :: rhs
        type (lib_ml_fmm_v) :: rv

        allocate(rv%dummy, source = lhs%dummy + rhs%dummy)
!        rv%dummy = lhs%dummy + rhs%dummy
    end function

    function test_v_operator_sub(lhs, rhs) result(rv)
        use ml_fmm_type
        implicit none
        ! dummy
        type (lib_ml_fmm_v), dimension(:), intent(in) :: lhs
        type (lib_ml_fmm_v), dimension(size(lhs)), intent(in) :: rhs
        type (lib_ml_fmm_v), dimension(size(lhs)) :: rv

        ! auxilary
        integer :: i

        do i=1, size(lhs)
            allocate(rv(i)%dummy, source = lhs(i)%dummy - rhs(i)%dummy)
!            rv(i)%dummy = lhs(i)%dummy - rhs(i)%dummy
        end do
    end function

    function test_v_operator_sub_0d(lhs, rhs) result(rv)
        use ml_fmm_type
        implicit none
        ! dummy
        type (lib_ml_fmm_v), intent(in) :: lhs
        type (lib_ml_fmm_v), intent(in) :: rhs
        type (lib_ml_fmm_v) :: rv

        allocate(rv%dummy, source = lhs%dummy - rhs%dummy)
!        rv%dummy = lhs%dummy - rhs%dummy
    end function

    function test_get_u_B_i(x, data_element, element_number) result(u_B_i)
        use lib_tree_public
        use ml_fmm_type
        implicit none
        ! dummy
        type(lib_tree_spatial_point), intent(in) :: x
        type(lib_tree_data_element), intent(in) :: data_element
        integer(kind=4), intent(in) :: element_number

        type(lib_ml_fmm_coefficient) :: u_B_i

!        allocate(B_i%dummy, source = (/data_element%uindex%n + abs(x)/))
        allocate(u_B_i%dummy, source = (/real(data_element%uindex%n, kind=LIB_ML_FMM_COEFFICIENT_KIND)/))
    end function

    function test_get_u_phi_i_j(data_element_i, element_number_i, y_j, element_number_j) result(rv)
        use lib_tree_public
        use ml_fmm_type
        implicit none
        ! dummy
        type(lib_tree_data_element), intent(in) :: data_element_i
        integer(kind=4), intent(in) :: element_number_i
        type(lib_tree_spatial_point), intent(in) :: y_j
        integer(kind=4), intent(in) :: element_number_j
        type(lib_ml_fmm_v) :: rv

!        allocate(rv%dummy, source = (/data_element_i%uindex%n + abs(y_j)/))
        allocate(rv%dummy, source = (/real(data_element_i%uindex%n, kind=LIB_ML_FMM_COEFFICIENT_KIND)/))

    end function

    function test_translation_RR(A_i_1, x_1, x_2) result(A_i_2)
        use lib_tree_public
        use ml_fmm_type
        implicit none
        ! dummy
        type(lib_ml_fmm_coefficient), intent(in) :: A_i_1
        type(lib_tree_spatial_point), intent(in) :: x_1
        type(lib_tree_spatial_point), intent(in) :: x_2
        type(lib_ml_fmm_coefficient) :: A_i_2


!        allocate(A_i_2%dummy, source = (/A_i_1%dummy + abs(x_2 - x_1)/))
        allocate(A_i_2%dummy, source = (/A_i_1%dummy(1)/))

    end function

    function test_translation_SR(B_i_1, x_1, x_2) result(A_i_2)
        use lib_tree_public
        use ml_fmm_type
        implicit none
        ! dummy
        type(lib_ml_fmm_coefficient), intent(in) :: B_i_1
        type(lib_tree_spatial_point), intent(in) :: x_1
        type(lib_tree_spatial_point), intent(in) :: x_2
        type(lib_ml_fmm_coefficient) :: A_i_2

!        allocate(A_i_2%dummy, source = (/B_i_1%dummy + abs(x_2 - x_1)/))
        allocate(A_i_2%dummy, source = (/B_i_1%dummy(1)/))

    end function

    function test_translation_SS(B_i_1, x_1, x_2) result(B_i_2)
        use lib_tree_public
        use ml_fmm_type
        implicit none
        ! dummy
        type(lib_ml_fmm_coefficient), intent(in) :: B_i_1
        type(lib_tree_spatial_point), intent(in) :: x_1
        type(lib_tree_spatial_point), intent(in) :: x_2
        type(lib_ml_fmm_coefficient) :: B_i_2

!        allocate(B_i_2%dummy, source = (/B_i_1%dummy(1) + abs(x_2 - x_1)/))
        allocate(B_i_2%dummy, source = (/B_i_1%dummy(1)/))

    end function
end module ml_fmm_math
