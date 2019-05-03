module lib_ml_fmm_type_operator
    implicit none

    private

    ! ----- test functions -----
    public :: lib_ml_fmm_type_operator_constructor
    public :: lib_ml_fmm_type_operator_test_functions

    public :: operator (*)
    public :: operator (+)
    public :: operator (.CoR.)
    public :: m_coefficient_set_zero

    ! ----- operator -----
    interface operator (+)
        module procedure m_coefficient_add
    end interface

    interface operator (*)
        module procedure m_u_dot_coefficient
    end interface

    ! Sum_q=0_p-1[C_q R_q(y_j − x_∗)] = C o R(y_j − x_∗)
    !                                     ^ cor-operator
    ! Reference: Data_Structures_Optimal_Choice_of_Parameters_and_C, p. 6
    interface operator (.CoR.)
        module procedure m_cor
    end interface

    ! ----- interfaces -----
    interface
        function ml_fmm_coefficient_add_operator(lhs,rhs) result (rv)
            use lib_ml_fmm_type
            implicit none
            ! dummy
            type (lib_ml_fmm_coefficient), intent(in) :: lhs, rhs
            type (lib_ml_fmm_coefficient) :: rv
        end function ml_fmm_coefficient_add_operator

        function ml_fmm_u_dot_coefficient_operator(lhs,rhs) result (rv)
            use lib_ml_fmm_type
            implicit none
            ! dummy
            type (lib_ml_fmm_v), intent(in) :: lhs
            type (lib_ml_fmm_coefficient), intent(in) :: rhs
            type (lib_ml_fmm_coefficient) :: rv
        end function ml_fmm_u_dot_coefficient_operator

        function ml_fmm_cor_operator(lhs, rhs) result(rv)
            use lib_ml_fmm_type
            implicit none
            ! dummy
            type(lib_ml_fmm_coefficient), intent(in) :: lhs
            type(lib_ml_fmm_R), intent(in) :: rhs
            type(lib_ml_fmm_v) :: rv
        end function ml_fmm_cor_operator

        subroutine ml_fmm_coefficient_set_zero(coefficient)
            use lib_ml_fmm_type
            implicit none
            ! dummy
            type(lib_ml_fmm_coefficient), intent(inout) :: coefficient
        end subroutine
    end interface

    ! ----- member procedures -----
    procedure(ml_fmm_coefficient_add_operator), pointer :: m_coefficient_add => null()
    procedure(ml_fmm_u_dot_coefficient_operator), pointer :: m_u_dot_coefficient => null()
    procedure(ml_fmm_cor_operator), pointer :: m_cor => null()
    procedure(ml_fmm_coefficient_set_zero), pointer :: m_coefficient_set_zero => null()

    contains

    subroutine lib_ml_fmm_type_operator_constructor(c_add, u_dot_coefficient, cor, coefficient_set_zero)
        use lib_ml_fmm_type
        implicit none
        procedure(ml_fmm_coefficient_add_operator) :: c_add
        procedure(ml_fmm_u_dot_coefficient_operator) :: u_dot_coefficient
        procedure(ml_fmm_cor_operator) :: cor
        procedure(ml_fmm_coefficient_set_zero) :: coefficient_set_zero


        m_coefficient_add => c_add
        m_u_dot_coefficient => u_dot_coefficient
        m_cor => cor
        m_coefficient_set_zero => coefficient_set_zero

    end subroutine lib_ml_fmm_type_operator_constructor

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
            use lib_ml_fmm_type
            implicit none
            ! dummy
            logical :: rv

            ! auxiliary
            integer(kind=1), parameter :: length = 3
            type(lib_ml_fmm_v) :: lhs_u
            type(lib_ml_fmm_coefficient) :: lhs_coeff
            type(lib_ml_fmm_coefficient) :: rhs_coeff
            type(lib_ml_fmm_coefficient) :: res_coeff
            type(lib_ml_fmm_R) :: rhs_R
            type(lib_ml_fmm_v) :: res_v

            call lib_ml_fmm_type_operator_constructor(test_c_add, test_u_dot_coefficient, test_cor, test_coefficient_set_zero)

            allocate (lhs_u%dummy(length))
            allocate (lhs_coeff%dummy(length))
            allocate (rhs_coeff%dummy(length))
            allocate (rhs_R%dummy(length))

            lhs_coeff%dummy = (/5, 4, 3/)
            call m_coefficient_set_zero(lhs_coeff)
            lhs_coeff%dummy = (/5, 4, 3/)
            rhs_coeff%dummy = (/2, 1, 3/)

            lhs_u%dummy = (/1, 2, 3/)

            res_coeff = lhs_coeff + rhs_coeff

            res_coeff = lhs_u * res_coeff

            rhs_R%dummy = (/0, 7, 1/)
            res_v = lhs_coeff .cor. rhs_R

            rv = .false.

        end function test_lib_ml_fmm_type_operator_constructor

        function test_c_add(lhs, rhs) result (rv)
            use lib_ml_fmm_type
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
            use lib_ml_fmm_type
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
            use lib_ml_fmm_type
            implicit none
            ! dummy
            type(lib_ml_fmm_coefficient), intent(in) :: lhs
            type(lib_ml_fmm_R), intent(in) :: rhs
            type(lib_ml_fmm_v) :: rv

            ! auxiliary
            integer :: i
            integer :: length

            allocate (rv%dummy(1))
            rv%dummy(1) = 0

            if (size(lhs%dummy) .eq. size(rhs%dummy)) then
                length = size(lhs%dummy)
                do i=1, length
                    rv%dummy(1) = rv%dummy(1) + lhs%dummy(i) * rhs%dummy(i)
                end do
            end if

        end function test_cor

        subroutine test_coefficient_set_zero(coefficient)
            use lib_ml_fmm_type
            implicit none
            ! dummy
            type(lib_ml_fmm_coefficient), intent(inout) :: coefficient

            if (allocated(coefficient%dummy)) then
                coefficient%dummy(:) = 0
            end if
        end subroutine

    end function lib_ml_fmm_type_operator_test_functions

end module lib_ml_fmm_type_operator
