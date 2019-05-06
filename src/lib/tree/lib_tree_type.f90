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

!    interface operator (+)
!        module procedure lib_tree_type_spatial_point_operator_add
!    end interface
!
!    interface operator (-)
!        module procedure lib_tree_type_spatial_point_operator_sub
!    end interface
!
!    interface operator (*)
!        module procedure lib_tree_type_spatial_point_operator_mul_scalar
!        module procedure lib_tree_type_spatial_point_operator_scalar_product
!    end interface

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

!    contains
!
!    function lib_tree_type_spatial_point_operator_add(lhs, rhs) result(rv)
!        implicit none
!        ! dummy
!        type (lib_tree_spatial_point), intent(in) :: lhs
!        type (lib_tree_spatial_point), intent(in):: rhs
!        type (lib_tree_spatial_point) :: rv
!
!        ! auxilary
!        integer :: i
!
!        do i=1, TREE_DIMENSIONS
!            rv%x(i) = lhs%x(i) + rhs%x(i)
!        end do
!    end function
!
!    function lib_tree_type_spatial_point_operator_sub(lhs, rhs) result(rv)
!        implicit none
!        ! dummy
!        type (lib_tree_spatial_point), intent(in) :: lhs
!        type (lib_tree_spatial_point), intent(in):: rhs
!        type (lib_tree_spatial_point) :: rv
!
!        ! auxilary
!        integer :: i
!
!        do i=1, TREE_DIMENSIONS
!            rv%x(i) = lhs%x(i) - rhs%x(i)
!        end do
!    end function
!
!    function lib_tree_type_spatial_point_operator_mul_scalar(lhs, rhs) result(rv)
!        use lib_tree_type
!        implicit none
!        ! dummy
!        type (lib_tree_spatial_point), intent(in) :: lhs
!#if (_SPATIAL_POINT_IS_DOUBLE_ == 1)
!        double precision, intent(in):: rhs
!#elif (_SPATIAL_POINT_IS_DOUBLE_ == 0)
!        real, intent(in):: rhs
!#endif
!        type (lib_tree_spatial_point) :: rv
!
!        ! auxilary
!        integer :: i
!
!        do i=1, TREE_DIMENSIONS
!            rv%x(i) = lhs%x(i) * rhs%x(i)
!        end do
!    end function
!
!    function lib_tree_type_spatial_point_operator_scalar_mul(lhs, rhs) result(rv)
!        implicit none
!        ! dummy
!#if (_SPATIAL_POINT_IS_DOUBLE_ == 1)
!        double precision, intent(in):: lhs
!#elif (_SPATIAL_POINT_IS_DOUBLE_ == 0)
!        real, intent(in):: lhs
!#endif
!        type (lib_tree_spatial_point), intent(in) :: rhs
!        type (lib_tree_spatial_point) :: rv
!
!        ! auxilary
!        integer :: i
!
!        do i=1, TREE_DIMENSIONS
!            rv%x(i) = lhs%x(i) * rhs%x(i)
!        end do
!    end function
!
!    function lib_tree_type_spatial_point_operator_scalar_product(lhs, rhs) result(rv)
!        implicit none
!        ! dummy
!        type (lib_tree_spatial_point), intent(in) :: lhs
!        type (lib_tree_spatial_point), intent(in) :: rhs
!#if (_SPATIAL_POINT_IS_DOUBLE_ == 1)
!        double precision :: rv
!#elif (_SPATIAL_POINT_IS_DOUBLE_ == 0)
!        real :: rv
!#endif
!        ! auxilary
!        integer :: i
!
!        rv = 0
!        do i=1, TREE_DIMENSIONS
!            rv = rv + rhs(i) * lhs(ii)
!        end do
!    end function
end module lib_tree_type
