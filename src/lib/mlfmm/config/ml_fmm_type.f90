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

    type lib_ml_fmm_uindex_list
        type(lib_tree_universal_index), dimension(:), allocatable :: uindex
    end type

    type lib_ml_fmm_hashed_coeffcient_index
        integer(kind=LIB_ML_FMM_COEFFICIENT_KIND) :: array_position
        integer(kind=2) :: number_of_hash_runs
    end type

    ! hierarchy of the data structure for the ML FMM
    ! Reference: Reference: Data_Structures_Optimal_Choice_of_Parameters_and_C, chapter 2. Data hierarchies
    !
    ! This type represents a dataset of a hierarchy level.
    !
    ! Example
    ! ----
    !
    !                   -----
    !                   | 0 |                 l = 0
    !                   -----
    !         -----               -----
    !         | 0 |               | 1 |       l = 1
    !         -----            v  -----
    !    -----     -----     -----     -----
    !    | 0 |     | 1 |     | 2 |     | 3 |  l = 2    <---
    !    -----     -----     -----     -----
    !
    !     Variable                    Value
    !    -----------------------------------
    !     is_hashed                     F
    !     coefficient_type              C
    !     number_of_boxes               1
    !     type                          X
    !     coefficient_list_index(2)     1
    !     size(coefficient_list)        1
    !
    !   Only one box contains coefficients at level 1. These coefficients are C expansion coefficients

    type lib_ml_fmm_hierarchy
        type(lib_ml_fmm_coefficient), dimension(:), allocatable :: coefficient_list
        type(lib_ml_fmm_hashed_coeffcient_index), dimension(:), allocatable :: hashed_coefficient_list_index
        integer(kind=LIB_ML_FMM_COEFFICIENT_KIND), dimension(:), allocatable :: coefficient_list_index
        logical :: is_hashed                                            ! true: access coefficient with hashed uindex%n, false: access coefficient with uindex%n
        integer(kind=1) :: coefficient_type                             ! [C, D^~, D]
        integer(kind=LIB_ML_FMM_COEFFICIENT_KIND) :: number_of_boxes    ! number of not empty boxes
        integer(kind=1) :: type                                         ! [X, Y, XY] hierarchy type
    end type
end module ml_fmm_type
