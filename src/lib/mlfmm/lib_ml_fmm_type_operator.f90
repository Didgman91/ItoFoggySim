module lib_ml_fmm_type_operator
    implicit none

    private

    ! ----- public functions -----
    public :: lib_ml_fmm_type_operator_constructor

    public :: operator (*)
    public :: operator (+)
    public :: operator (.CoR.)

    public :: lib_ml_fmm_type_operator_set_coefficient_zero
!    public :: lib_ml_fmm_type_operator_allocate_coefficient_list
!    public :: lib_ml_fmm_type_operator_deallocate_coefficient_list
    public :: lib_ml_fmm_type_operator_set_coefficient
    public :: lib_ml_fmm_type_operator_get_coefficient

    ! ---- public type definitions -----
    public :: ml_fmm_type_operator_procedures
!    public :: ml_fmm_coefficient_add_operator

    ! ----- operator -----
    interface operator (+)
        module procedure m_coefficient_add
!        module procedure m_procedures%m_coefficient_add
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
            use ml_fmm_type
            implicit none
            ! dummy
            type (lib_ml_fmm_coefficient), intent(in) :: lhs, rhs
            type (lib_ml_fmm_coefficient) :: rv
        end function ml_fmm_coefficient_add_operator

        function ml_fmm_u_dot_coefficient_operator(lhs,rhs) result (rv)
            use ml_fmm_type
            implicit none
            ! dummy
            type (lib_ml_fmm_v), intent(in) :: lhs
            type (lib_ml_fmm_coefficient), intent(in) :: rhs
            type (lib_ml_fmm_coefficient) :: rv
        end function ml_fmm_u_dot_coefficient_operator

        function ml_fmm_cor_operator(lhs, rhs) result(rv)
            use ml_fmm_type
            implicit none
            ! dummy
            type(lib_ml_fmm_coefficient), intent(in) :: lhs
            type(lib_ml_fmm_R), intent(in) :: rhs
            type(lib_ml_fmm_v) :: rv
        end function ml_fmm_cor_operator

        subroutine ml_fmm_coefficient_set_zero(coefficient)
            use ml_fmm_type
            implicit none
            ! dummy
            type(lib_ml_fmm_coefficient), intent(inout) :: coefficient
        end subroutine

!        subroutine ml_fmm_allocate_coefficient_list(coefficient_list, l_min, l_max)
!            use ml_fmm_type
!            implicit none
!            ! dummy
!            type(lib_ml_fmm_coefficient_list_list) :: coefficient_list
!            integer(kind=1) :: l_min
!            integer(kind=1) :: l_max
!        end subroutine
!
!        subroutine ml_fmm_deallocate_coefficient_list(coefficient_list)
!            use ml_fmm_type
!            implicit none
!            type(lib_ml_fmm_coefficient_list_list), intent(inout) :: coefficient_list
!        end subroutine

        subroutine ml_fmm_set_coefficient(coefficient, uindex, hierarchy)
            use ml_fmm_type
            use lib_tree_type
            implicit none
            ! dummy
            type(lib_tree_universal_index), intent(in) :: uindex
            type(lib_ml_fmm_coefficient), intent(in) :: coefficient
            type(lib_ml_fmm_hierarchy), dimension(:), allocatable, intent(inout) :: hierarchy
        end subroutine

        function ml_fmm_get_coefficient(uindex, hierarchy) result(coefficient)
            use ml_fmm_type
            use lib_tree_type
            implicit none
            ! dummy
            type(lib_tree_universal_index), intent(in) :: uindex
            type(lib_ml_fmm_hierarchy), dimension(:), allocatable, intent(inout) :: hierarchy
            type(lib_ml_fmm_coefficient) :: coefficient
        end function
    end interface

    type ml_fmm_type_operator_procedures
        procedure(ml_fmm_coefficient_add_operator), pointer, nopass :: coefficient_add => null()
        procedure(ml_fmm_u_dot_coefficient_operator), pointer, nopass :: u_dot_coefficient => null()
        procedure(ml_fmm_cor_operator), pointer, nopass :: cor => null()
        procedure(ml_fmm_coefficient_set_zero), pointer, nopass :: coefficient_set_zero => null()
!        procedure(ml_fmm_allocate_coefficient_list), pointer, nopass :: allocate_coefficient_list => null()
        procedure(ml_fmm_set_coefficient), pointer, nopass :: set_coefficient => null()
        procedure(ml_fmm_get_coefficient), pointer, nopass :: get_coefficient => null()
!        procedure(ml_fmm_deallocate_coefficient_list), pointer, nopass :: deallocate_coefficient_list => null()
    end type

    ! ----- member procedures -----
    procedure(ml_fmm_coefficient_add_operator), pointer :: m_coefficient_add => null()
    procedure(ml_fmm_u_dot_coefficient_operator), pointer :: m_u_dot_coefficient => null()
    procedure(ml_fmm_cor_operator), pointer :: m_cor => null()
    procedure(ml_fmm_coefficient_set_zero), pointer :: m_coefficient_set_zero => null()
!    procedure(ml_fmm_allocate_coefficient_list), pointer :: m_allocate_coefficient_list => null()
    procedure(ml_fmm_set_coefficient), pointer :: m_set_coefficient => null()
    procedure(ml_fmm_get_coefficient), pointer :: m_get_coefficient => null()
!    procedure(ml_fmm_deallocate_coefficient_list), pointer :: m_deallocate_coefficient_list => null()

    contains

    subroutine lib_ml_fmm_type_operator_constructor(operator_procedures)
!    subroutine lib_ml_fmm_type_operator_constructor(coefficient_add, u_dot_coefficient, cor, &
!                                                    set_coefficient_zero)!, &
!                                                    allocate_coefficient_list, &
!                                                    set_coefficient, get_coefficient, &
!                                                    deallocate_coefficient_list)
        use ml_fmm_type
        implicit none
        ! dummy
!        procedure(ml_fmm_coefficient_add_operator) :: coefficient_add
!        procedure(ml_fmm_u_dot_coefficient_operator) :: u_dot_coefficient
!        procedure(ml_fmm_cor_operator) :: cor
!        procedure(ml_fmm_coefficient_set_zero) :: set_coefficient_zero
!        procedure(ml_fmm_allocate_coefficient_list) :: allocate_coefficient_list
!        procedure(ml_fmm_set_coefficient) :: set_coefficient
!        procedure(ml_fmm_get_coefficient) :: get_coefficient
!        procedure(ml_fmm_deallocate_coefficient_list) :: deallocate_coefficient_list
!
!        m_coefficient_add => coefficient_add
!        m_u_dot_coefficient => u_dot_coefficient
!        m_cor => cor
!        m_coefficient_set_zero => set_coefficient_zero
!        m_allocate_coefficient_list => allocate_coefficient_list
!        m_set_coefficient => set_coefficient
!        m_get_coefficient => get_coefficient
!        m_deallocate_coefficient_list => deallocate_coefficient_list


        type(ml_fmm_type_operator_procedures) :: operator_procedures
        m_coefficient_add => operator_procedures%coefficient_add
        m_u_dot_coefficient => operator_procedures%u_dot_coefficient
        m_cor => operator_procedures%cor
        m_coefficient_set_zero => operator_procedures%coefficient_set_zero
!        m_allocate_coefficient_list => operator_procedures%allocate_coefficient_list
        m_set_coefficient => operator_procedures%set_coefficient
        m_get_coefficient => operator_procedures%get_coefficient
!        m_deallocate_coefficient_list => operator_procedures%deallocate_coefficient_list

    end subroutine lib_ml_fmm_type_operator_constructor

    subroutine lib_ml_fmm_type_operator_set_coefficient_zero(coefficient)
        use ml_fmm_type
        implicit none
        ! dummy
        type(lib_ml_fmm_coefficient), intent(inout) :: coefficient

        if (associated (m_coefficient_set_zero)) then
            call m_coefficient_set_zero(coefficient)
        else
            print *, "lib_ml_fmm_type_operator_set_coefficient_zero:  ERROR"
            print *, "  m_coefficient_set_zero is not associated"
        end if
    end subroutine lib_ml_fmm_type_operator_set_coefficient_zero

!    subroutine lib_ml_fmm_type_operator_allocate_coefficient_list(coefficient_list, l_min, l_max)
!            use ml_fmm_type
!            implicit none
!            ! dummy
!            type(lib_ml_fmm_coefficient_list_list) :: coefficient_list
!            integer(kind=1) :: l_min
!            integer(kind=1) :: l_max
!
!            if ( associated(m_allocate_coefficient_list) ) then
!                call m_allocate_coefficient_list(coefficient_list, l_min, l_max)
!            else
!                print *, "lib_ml_fmm_type_operator_allocate_coefficient_list:  ERROR"
!                print *, "  m_allocate_coefficient_list is not associated"
!            end if
!    end subroutine lib_ml_fmm_type_operator_allocate_coefficient_list
!
!    subroutine lib_ml_fmm_type_operator_deallocate_coefficient_list(coefficient_list)
!            use ml_fmm_type
!            implicit none
!            ! dummy
!            type(lib_ml_fmm_coefficient_list_list) :: coefficient_list
!            integer(kind=1) :: l_min
!            integer(kind=1) :: l_max
!
!            if ( associated(m_deallocate_coefficient_list) ) then
!                call m_deallocate_coefficient_list(coefficient_list)
!            else
!                print *, "lib_ml_fmm_type_operator_deallocate_coefficient_list:  ERROR"
!                print *, "  m_deallocate_coefficient_list is not associated"
!            end if
!    end subroutine lib_ml_fmm_type_operator_deallocate_coefficient_list

    subroutine lib_ml_fmm_type_operator_set_coefficient(coefficient, uindex, hierarchy)
            use ml_fmm_type
            use lib_tree_type
            implicit none
            ! dummy
            type(lib_tree_universal_index), intent(inout) :: uindex
            type(lib_ml_fmm_coefficient), intent(inout) :: coefficient
            type(lib_ml_fmm_hierarchy), dimension(:), allocatable, intent(inout) :: hierarchy

            if ( associated(m_set_coefficient)) then
                call m_set_coefficient(coefficient, uindex, hierarchy)
            else
                print *, "lib_ml_fmm_type_operator_set_coefficient:  ERROR"
                print *, "  m_set_coefficient is not associated"
            end if
        end subroutine

        function lib_ml_fmm_type_operator_get_coefficient(uindex, hierarchy) result(coefficient)
            use ml_fmm_type
            use lib_tree_type
            implicit none
            ! dummy
            type(lib_tree_universal_index), intent(inout) :: uindex
            type(lib_ml_fmm_hierarchy), dimension(:), allocatable, intent(inout) :: hierarchy
            type(lib_ml_fmm_coefficient) :: coefficient

            if ( associated(m_set_coefficient)) then
                coefficient = m_get_coefficient(uindex, hierarchy)
            else
                print *, "lib_ml_fmm_type_operator_get_coefficient:  ERROR"
                print *, "  m_get_coefficient is not associated"
            end if
        end function

end module lib_ml_fmm_type_operator
