module lib_data_types
    implicit none

    integer, parameter :: element_kind = 4

    ! list element
    type octree_element
        integer(kind=element_kind)   :: start
        integer(kind=element_kind)   :: end
        integer(kind=element_kind)   :: parent_element
    end type octree_element

!    type octree_level
!        type(octree_element), dimension() ::
!    end type octree_level

end module lib_data_types
