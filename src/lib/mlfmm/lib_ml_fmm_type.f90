! LIB: Mulitlevel Fast Multipole Method - type definitions
!
module lib_ml_fmm_type
    use lib_tree_type
    implicit none

    integer(kind=1), parameter :: LIB_ML_FMM_COEFFICIENT_KIND = 8

    type lib_ml_fmm_v
        real(kind=LIB_ML_FMM_COEFFICIENT_KIND), dimension(:), allocatable :: dummy
    end type

    ! e.g. A, B, C or D coefficient
    ! Reference: Data_Structures_Optimal_Choice_of_Parameters_and_C, p. 6
    type lib_ml_fmm_coefficient
        real(kind=LIB_ML_FMM_COEFFICIENT_KIND), dimension(:), allocatable :: dummy
    end type

    ! values of the local (regular) basis function
    type lib_ml_fmm_R
        real(kind=LIB_ML_FMM_COEFFICIENT_KIND), dimension(:), allocatable :: dummy
    end type

    ! values of the multipole (singular) basis function
    type lib_ml_fmm_S
        real(kind=LIB_ML_FMM_COEFFICIENT_KIND), dimension(:), allocatable :: dummy
    end type

end module lib_ml_fmm_type
