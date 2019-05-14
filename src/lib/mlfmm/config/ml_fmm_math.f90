module ml_fmm_math
    use ml_fmm_type
    use lib_ml_fmm_type_operator
    implicit none

    private

    ! --- public functions ---
    public :: ml_fmm_type_operator_get_procedures

    public :: lib_ml_fmm_type_operator_test_functions

    contains

    function ml_fmm_type_operator_get_procedures() result (operator_procedures)
        implicit none
        ! dummy
        type(ml_fmm_type_operator_procedures) :: operator_procedures

        ! load test procedure functions
        operator_procedures%coefficient_add => test_c_add
        operator_procedures%coefficient_set_zero => test_set_coefficient_zero
        operator_procedures%cor => test_cor
        operator_procedures%u_dot_coefficient => test_u_dot_coefficient

    end function ml_fmm_type_operator_get_procedures

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
            type(lib_ml_fmm_v) :: res_v

            type(ml_fmm_type_operator_procedures) :: ml_fmm_operator_procedures

            ml_fmm_operator_procedures%coefficient_add => test_c_add
            ml_fmm_operator_procedures%coefficient_set_zero => test_set_coefficient_zero
            ml_fmm_operator_procedures%cor => test_cor
            ml_fmm_operator_procedures%u_dot_coefficient => test_u_dot_coefficient

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

            res_coeff = lhs_u * res_coeff

            rhs_R%x(:2) = (/1.5, 7.0/)
            res_v = lhs_coeff .cor. rhs_R

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

    function test_cor(lhs, rhs)  result (rv)
        implicit none
        ! dummy
        type(lib_ml_fmm_coefficient), intent(in) :: lhs
        type(lib_tree_spatial_point), intent(in) :: rhs
        type(lib_ml_fmm_v) :: rv

        ! auxiliary
        integer :: i
        integer :: length

        allocate (rv%dummy(1))
        rv%dummy(1) = 0

        length = size(lhs%dummy)
        do i=1, length
            rv%dummy(1) = rv%dummy(1) + lhs%dummy(i) * rhs%x(1)
        end do

    end function test_cor

    subroutine test_set_coefficient_zero(coefficient)
        implicit none
        ! dummy
        type(lib_ml_fmm_coefficient), intent(inout) :: coefficient

        if (allocated(coefficient%dummy)) then
            coefficient%dummy(:) = 0
        end if
    end subroutine
end module ml_fmm_math
