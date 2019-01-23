module lib_data_types
    implicit none

    integer, parameter :: element_kind = 4

    ! list element
    type lib_octree_element
        integer(kind=element_kind)   :: start
        integer(kind=element_kind)   :: end
        integer(kind=element_kind)   :: parent_element
    end type lib_octree_element

!    type octree_parent
!        type(octree_element), dimension(8) :: children
!        type(ectree_element), dimension() :: neighbor
!    end type octree_level



end module lib_data_types
