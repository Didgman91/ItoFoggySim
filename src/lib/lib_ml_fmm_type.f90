! LIB: Mulitlevel Fast Multipole Method - type definitions
!
module lib_ml_fmm_type
    use lib_tree_type
    implicit none

    type lib_ml_fmm_u
        integer, dimension(:), allocatable :: u
    end type

!    type lib_ml_fmm_A
!        integer :: dummy
!    end type

    type lib_ml_fmm_A_i
        integer :: dummy
    end type

!    type lib_ml_fmm_B
!        integer :: dummy
!    end type

    type lib_ml_fmm_B_i
        integer :: dummy
    end type

    type lib_ml_fmm_C
        integer :: dummy
    end type

    type lib_ml_fmm_R
        integer :: dummy
    end type

    type lib_ml_fmm_S
        integer :: dummy
    end type

end module lib_ml_fmm_type
