module lib_octree
use lib_octree_helper_functions
    implicit none
    ! Data Structures, Optimal Choice of Parameters, and Complexity Results for Generalized Multilevel Fast Multipole Methods in d Dimensions

!    private
!
!    public :: lib_octree_get_parent

    integer, parameter :: octree_integer_kind = 4
    integer, parameter :: fmm_d = 3 ! dimensions

    type lib_octree_box_index
        integer(kind=octree_integer_kind)   :: n
        integer(kind=octree_integer_kind)   :: l
    end type lib_octree_box_index

!    ! a spatial point corresponds to our element, which consists of nodes
!    type lib_octree_spatial_point
!        integer(kind=octree_integer_kind), dimension(fmm_d)   :: x
!    end type lib_octree_spatial_point

    contains

    function lib_octree_get_parent(n,l) result (rv)
    implicit none
        ! Arguments
        ! ----
        !   n
        !       number of the node, with n ranging form 0 to 2^(3*l) in a three-dimensional space
        !   l
        !       number of the level
        !
        ! Returns
        ! ----
        !   the index of the parent box

        ! dummy arguments
        integer(kind=octree_integer_kind), intent (in) :: n
        integer(kind=octree_integer_kind), intent (in) :: l
        type(lib_octree_box_index) :: rv

!        Todo: calculate the id of the parent box

!        parent=
        rv%n = n
        rv%l = l
    end function lib_octree_get_parent

    function lib_octree_get_children(n,l) result (rv)
    implicit none
        ! Arguments
        ! ----
        !   n
        !       number of the node, with n ranging form 0 to 2^(3*l) in a three-dimensional space
        !   l
        !       number of the level
        !
        ! Returns
        ! ----
        !   the index of the parent box

        ! dummy arguments
        integer(kind=octree_integer_kind), intent (in) :: n
        integer(kind=octree_integer_kind), intent (in) :: l
        type(lib_octree_box_index) :: rv

!        Todo: calculate the id's of the children boxes
        rv%n = n
        rv%l = l
    end function lib_octree_get_children

    function lib_octree_get_neighbours(k,n,l) result (rv)
    implicit none
        ! Arguments
        ! ----
        !   k
        !       k-neighbour of the box (n,l)
        !   n
        !       number of the node, with n ranging form 0 to 2^(3*l) in a three-dimensional space
        !   l
        !       number of the level
        !
        ! Returns
        ! ----
        !   the indeces of the k-neighours

        ! dummy arguments
        integer(kind=octree_integer_kind), intent (in)  :: k
        integer(kind=octree_integer_kind), intent (in)  :: n
        integer(kind=octree_integer_kind), intent (in)  :: l
        type(lib_octree_box_index) :: rv

!        Todo: calculate the id's of the children boxes
        rv%n = n
        rv%l = l
    end function lib_octree_get_neighbours

    function lib_octree_get_domain_e1(n,l) result (rv)
    implicit none
        ! Arguments
        ! ----
        !   n
        !       number of the node, with n ranging form 0 to 2^(3*l) in a three-dimensional space
        !   l
        !       number of the level
        !
        ! Returns
        ! ----
        !   the spatial points inside the box (n,l)

        ! dummy arguments
        integer(kind=octree_integer_kind), intent (in)  :: n
        integer(kind=octree_integer_kind), intent (in)  :: l
        type(lib_octree_spatial_point) :: rv

!         rv%x = (/1,2,3/)
    end function lib_octree_get_domain_e1

    function lib_octree_get_domain_e2(k,n,l) result (rv)
    implicit none
        ! Arguments
        ! ----
        !   k
        !       k-neighbour of the box (n,l)
        !   n
        !       number of the node, with n ranging form 0 to 2^(3*l) in a three-dimensional space
        !   l
        !       number of the level
        !
        ! Returns
        ! ----
        !   the spatial points in the k-neighorhood of box (n,l)

        ! dummy arguments
        integer(kind=octree_integer_kind), intent (in) :: k
        integer(kind=octree_integer_kind), intent (in) :: n
        integer(kind=octree_integer_kind), intent (in) :: l
        type(lib_octree_spatial_point) :: rv

!        rv%x = (/1,2,3/)
    end function lib_octree_get_domain_e2

    function lib_octree_get_domain_e3(k,n,l) result (rv)
    implicit none
        ! Arguments
        ! ----
        !   k
        !       k-neighbour of the box (n,l)
        !   n
        !       number of the node, with n ranging form 0 to 2^(3*l) in a three-dimensional space
        !   l
        !       number of the level
        !
        ! Returns
        ! ----
        !   the spatial points outside the k-neighorhood of box (n,l)

        ! dummy arguments
        integer(kind=octree_integer_kind), intent (in) :: k
        integer(kind=octree_integer_kind), intent (in) :: n
        integer(kind=octree_integer_kind), intent (in) :: l
        type(lib_octree_spatial_point) :: rv

!        rv%x = (/1,2,3/)
    end function lib_octree_get_domain_e3

    function lib_octree_get_domain_e4(k,n,l) result (rv)
    implicit none
        ! Arguments
        ! ----
        !   k
        !       k-neighbour of the box (n,l)
        !   n
        !       number of the node, with n ranging form 0 to 2^(3*l) in a three-dimensional space
        !   l
        !       number of the level
        !
        ! Returns
        ! ----
        !   the spatial points in the k-neighorhood of parent box of box (n,l), which do not belong
        !   the k-neighbourhood of the box itself.

        ! dummy arguments
        integer(kind=octree_integer_kind), intent (in) :: k
        integer(kind=octree_integer_kind), intent (in) :: n
        integer(kind=octree_integer_kind), intent (in) :: l
        type(lib_octree_spatial_point) :: rv

!        rv%x = (/1,2,3/)
    end function lib_octree_get_domain_e4

end module lib_octree
