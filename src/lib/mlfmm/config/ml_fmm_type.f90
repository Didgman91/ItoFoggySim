! LIB: Mulitlevel Fast Multipole Method - type definitions
!
module ml_fmm_type
    use lib_tree_type
    implicit none

    integer(kind=1), parameter, public :: LIB_ML_FMM_COEFFICIENT_KIND = 8

    type lib_ml_fmm_v
        real(kind=LIB_ML_FMM_COEFFICIENT_KIND), dimension(:), allocatable :: dummy
    end type

    ! e.g. A, B, C or D coefficient
    ! Reference: Data_Structures_Optimal_Choice_of_Parameters_and_C, p. 6
    !
    !                   -----
    !                   | 0 |                 l = 0
    !                   -----
    !         -----               -----
    !         | 0 |               | 1 |       l = 1
    !         -----               -----
    !    -----     -----     -----     -----
    !    | 0 |     | 1 |     | 2 |     | 3 |  l = 2
    !    -----     -----     -----     -----
    !      ^
    !
    ! coefficients of the box uindex(n,l)
    !   e.g. (0,2): C_i: type(lib_ml_coefficient)
    !
    !
    type lib_ml_fmm_coefficient
        real(kind=LIB_ML_FMM_COEFFICIENT_KIND), dimension(:), allocatable :: dummy
    end type

    ! e.g. A, B, C or D coefficient
    ! Reference: Data_Structures_Optimal_Choice_of_Parameters_and_C, p. 6
    !
    !                   -----
    !                   | 0 |                 l = 0
    !                   -----
    !         -----               -----
    !         | 0 |               | 1 |       l = 1
    !         -----               -----
    !    -----     -----     -----     -----
    !    | 0 |     | 1 |     | 2 |     | 3 |  l = 2    <--
    !    -----     -----     -----     -----
    !
    ! list of coefficients of the boxes at a specific level
    !   e.g. ([0,3], 2): C_list = ( C(0,2), C(1,2), C(2,2), C(3,2) )
    type lib_ml_fmm_coefficient_list
        type(lib_ml_fmm_coefficient), dimension(:), allocatable :: dummy
    end type

    ! e.g. A, B, C or D coefficient
    ! Reference: Data_Structures_Optimal_Choice_of_Parameters_and_C, p. 6
    !
    !                   -----
    !                   | 0 |                 l = 0    <--
    !                   -----
    !         -----               -----
    !         | 0 |               | 1 |       l = 1    <--
    !         -----               -----
    !    -----     -----     -----     -----
    !    | 0 |     | 1 |     | 2 |     | 3 |  l = 2    <--
    !    -----     -----     -----     -----
    !
    ! list of coefficients of all boxes at a level gathered in a list
    !   e.g. (n, [2,0]): C_list_list = ( C([0,3],2), C([0,1],1), C(0,0) )
    type lib_ml_fmm_coefficient_list_list
        type(lib_ml_fmm_coefficient_list), dimension(:), allocatable :: dummy
    end type

    ! values of the local (regular) basis function
    type lib_ml_fmm_R
        real(kind=LIB_ML_FMM_COEFFICIENT_KIND), dimension(:), allocatable :: dummy
    end type

    ! values of the multipole (singular) basis function
    type lib_ml_fmm_S
        real(kind=LIB_ML_FMM_COEFFICIENT_KIND), dimension(:), allocatable :: dummy
    end type

end module ml_fmm_type
