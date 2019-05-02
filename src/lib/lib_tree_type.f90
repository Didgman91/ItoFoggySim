! spatial dimension, value = [2,3]
#define _FMM_DIMENSION_ 2

! 1: true, 0: false (-> spatial point is real)
#define _SPATIAL_POINT_IS_DOUBLE_ 1

! number of bytes of the universal index, value = [4,8,16]
! standard value: 8
!
! Constraint
! ----
!          |  _FMM_DIMENSION_  |
!   value  |    2     |   3    |
!   -----------------------------
!   single | [4,8]    | [8]    |
!   double | [8,16]   | [8,16] |
!
#define _UINDEX_BYTES_ 8

module lib_tree_type
    implicit none

    ! parameter
    integer(kind=1), public, parameter :: TREE_DIMENSIONS = _FMM_DIMENSION_ ! dimensions

#if(_SPATIAL_POINT_IS_DOUBLE_ == 1)
    integer(kind=1), public, parameter :: COORDINATE_BINARY_BYTES = 8
#elif(_SPATIAL_POINT_IS_DOUBLE_ == 0)
    integer(kind=1), public, parameter :: COORDINATE_BINARY_BYTES = 4
#endif

    integer(kind=1), public, parameter :: UINDEX_BYTES = _UINDEX_BYTES_
    ! ~ parameter ~

    ! type definitions
    type lib_tree_spatial_point
#if (_SPATIAL_POINT_IS_DOUBLE_ == 1)
        double precision, dimension(TREE_DIMENSIONS) :: x
#elif (_SPATIAL_POINT_IS_DOUBLE_ == 0)
        real, dimension(TREE_DIMENSIONS) :: x
#endif
    end type lib_tree_spatial_point

    type lib_tree_universal_index
        integer(kind=UINDEX_BYTES) :: n
!        integer(kind=INTERLEAVE_BITS_INTEGER_KIND), &
!            dimension(TREE_DIMENSIONS * COORDINATE_BINARY_BYTES/INTERLEAVE_BITS_INTEGER_KIND) &
!            :: interleaved_bits
        integer(kind=1) :: l
    end type lib_tree_universal_index
    ! ~ type definitions ~

        ! --- type definition ---
    type lib_tree_data_element
        type(lib_tree_spatial_point) :: point_x
        type(lib_tree_universal_index) :: uindex
        integer(kind=1) :: element_type
    end type lib_tree_data_element

end module lib_tree_type
