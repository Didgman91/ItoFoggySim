module lib_tree
use lib_tree_helper_functions
    implicit none
    ! Data Structures, Optimal Choice of Parameters, and Complexity Results for Generalized Multilevel Fast Multipole Methods in d Dimensions

    private

    ! parameter
    public :: TREE_DIMENSIONS

    ! test functions
    public :: lib_tree_hf_test_functions
    public :: lib_tree_hf_benchmark
    public :: lib_tree_hf_destructor

!    ! a spatial point corresponds to our element, which consists of nodes
!    type lib_tree_spatial_point
!        integer(kind=tree_integer_kind), dimension(fmm_d)   :: x
!    end type lib_tree_spatial_point

    contains

    function lib_tree_get_parent(uindex) result (rv)
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
        type(lib_tree_universal_index), intent (in) :: uindex
        type(lib_tree_universal_index) :: rv

        rv%n = lib_tree_hf_get_parent(uindex%n)
        rv%l = uindex%l - int(1,1)
    end function lib_tree_get_parent

    function lib_tree_get_children(uindex) result (rv)
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
        type(lib_tree_universal_index), intent (in) :: uindex
        integer(kind=COORDINATE_BINARY_BYTES), dimension(2**TREE_DIMENSIONS):: rv

        rv = lib_tree_hf_get_children_all(uindex%n)

    end function lib_tree_get_children

    function lib_tree_get_neighbours(k,uindex) result (rv)
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
        type(lib_tree_universal_index), intent (in) :: uindex
        integer(kind=COORDINATE_BINARY_BYTES) :: k
        integer(kind=COORDINATE_BINARY_BYTES), dimension(3**TREE_DIMENSIONS-1) :: rv

        rv = lib_tree_hf_get_neighbour_all_xD(k, uindex%n, uindex%l)

    end function lib_tree_get_neighbours

    function lib_tree_get_domain_e1(uindex) result (rv)
    implicit none
        ! Arguments
        ! ----
        !   uindex: type(lib_tree_universal_index)
        !       uindex%n
        !           number of the node, with n ranging form 0 to 2^(3*l) in a three-dimensional space
        !       uindex%l
        !           number of the level
        !
        ! Returns
        ! ----
        !   the spatial points inside the box (n,l)

        ! dummy arguments
        type(lib_tree_universal_index) :: uindex
        type(lib_tree_spatial_point) :: rv

        ! auxiliar
        type(lib_tree_universal_index), dimension(1) :: domain_box

        domain_box = uindex



    end function lib_tree_get_domain_e1

    function lib_tree_get_domain_e2(k,uindex) result (rv)
    implicit none
        ! Arguments
        ! ----
        !   k
        !       k-neighbour of the box (n,l)
        !   uindex: type(lib_tree_universal_index)
        !       uindex%n
        !           number of the node, with n ranging form 0 to 2^(3*l) in a three-dimensional space
        !       uindex%l
        !           number of the level
        !
        ! Returns
        ! ----
        !   the spatial points in the k-neighorhood of box (n,l)

        ! dummy arguments
        integer(kind=COORDINATE_BINARY_BYTES), intent (in) :: k
        type(lib_tree_universal_index), intent(in) :: uindex
        type(lib_tree_spatial_point) :: rv

        ! auxiliar
!        type(lib_tree_universal_index), dimension(:), allocatable :: domain_box
        integer(kind=COORDINATE_BINARY_BYTES) :: number_of_boxes
        integer(kind=COORDINATE_BINARY_BYTES), dimension(3**TREE_DIMENSIONS-1, k) :: buffer

        integer(kind=COORDINATE_BINARY_BYTES) :: i

        do i=1, k
            buffer(:, i) = lib_tree_hf_get_neighbour_all_xD(k, uindex%n, uindex%l)
        end do

!        rv%x = (/1,2,3/)
    end function lib_tree_get_domain_e2

    function lib_tree_get_domain_e3(k,n,l) result (rv)
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
        integer(kind=COORDINATE_BINARY_BYTES), intent (in) :: k
        integer(kind=COORDINATE_BINARY_BYTES), intent (in) :: n
        integer(kind=tree_integer_kind), intent (in) :: l
        type(lib_tree_spatial_point) :: rv

!        rv%x = (/1,2,3/)
    end function lib_tree_get_domain_e3

    function lib_tree_get_domain_e4(k,n,l) result (rv)
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
        integer(kind=COORDINATE_BINARY_BYTES), intent (in) :: k
        integer(kind=COORDINATE_BINARY_BYTES), intent (in) :: n
        integer(kind=tree_integer_kind), intent (in) :: l
        type(lib_tree_spatial_point) :: rv

!        rv%x = (/1,2,3/)
    end function lib_tree_get_domain_e4

end module lib_tree
