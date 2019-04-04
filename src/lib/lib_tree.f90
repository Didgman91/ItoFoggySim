
module lib_tree
use lib_tree_helper_functions
    implicit none
    ! Data Structures, Optimal Choice of Parameters, and Complexity Results for Generalized Multilevel Fast Multipole Methods in d Dimensions

    private

    ! parameter
    public :: TREE_DIMENSIONS

    public :: lib_tree_destructor

    ! test functions
    public :: lib_tree_hf_test_functions
    public :: lib_tree_hf_benchmark
    public :: lib_tree_hf_destructor

    public :: lib_tree_test_functions

    type lib_tree_data_element
        type(lib_tree_spatial_point) :: point_x
        integer(kind=1) :: element_type
    end type lib_tree_data_element

!    ! a spatial point corresponds to our element, which consists of nodes
    type(lib_tree_data_element), dimension (:), allocatable :: lib_tree_data

    contains

    ! cleans up the memory
    subroutine lib_tree_destructor()
        call lib_tree_hf_destructor()
    end subroutine lib_tree_destructor

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
    function lib_tree_get_parent(uindex) result (rv)
    implicit none


        ! dummy arguments
        type(lib_tree_universal_index), intent (in) :: uindex
        type(lib_tree_universal_index) :: rv

        rv%n = lib_tree_hf_get_parent(uindex%n)
        rv%l = uindex%l - int(1,1)
    end function lib_tree_get_parent

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
    function lib_tree_get_children(uindex) result (rv)
    implicit none
        ! dummy arguments
        type(lib_tree_universal_index), intent (in) :: uindex
        integer(kind=COORDINATE_BINARY_BYTES), dimension(2**TREE_DIMENSIONS):: rv

        rv = lib_tree_hf_get_children_all(uindex%n)

    end function lib_tree_get_children

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
    function lib_tree_get_neighbours(k,uindex) result (rv)
    implicit none
        ! dummy arguments
        type(lib_tree_universal_index), intent (in) :: uindex
        integer(kind=COORDINATE_BINARY_BYTES) :: k
        integer(kind=COORDINATE_BINARY_BYTES), dimension(3**TREE_DIMENSIONS-1) :: rv

        rv = lib_tree_hf_get_neighbour_all_xD(k, uindex%n, uindex%l)

    end function lib_tree_get_neighbours

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
    function lib_tree_get_domain_e1(uindex) result (rv)
    implicit none
        ! dummy arguments
        type(lib_tree_universal_index) :: uindex
        type(lib_tree_spatial_point) :: rv

        ! auxiliar
        type(lib_tree_universal_index), dimension(1) :: domain_box

        domain_box = uindex



    end function lib_tree_get_domain_e1

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
    function lib_tree_get_domain_e2(k,uindex) result (rv)
    implicit none
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
    function lib_tree_get_domain_e3(k,n,l) result (rv)
    implicit none
        ! dummy arguments
        integer(kind=COORDINATE_BINARY_BYTES), intent (in) :: k
        integer(kind=COORDINATE_BINARY_BYTES), intent (in) :: n
        integer(kind=tree_integer_kind), intent (in) :: l
        type(lib_tree_spatial_point) :: rv

!        rv%x = (/1,2,3/)
    end function lib_tree_get_domain_e3

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
    function lib_tree_get_domain_e4(k,n,l) result (rv)
    implicit none
        ! dummy arguments
        integer(kind=COORDINATE_BINARY_BYTES), intent (in) :: k
        integer(kind=COORDINATE_BINARY_BYTES), intent (in) :: n
        integer(kind=tree_integer_kind), intent (in) :: l
        type(lib_tree_spatial_point) :: rv

!        rv%x = (/1,2,3/)
    end function lib_tree_get_domain_e4

!    function lib_tree_get_

    ! Sorting data and saving in the module globally variable *lib_tree_data*
    !
    ! Arguments
    ! ----
    !   element_list
    !       list of data elements
    !
    !   treshold_level
    !       Threshold value (l_max) at which two adjacent elements can still be distinguished.
    !
    ! Returns
    !
    !
    function lib_tree_create_correspondence_vector(element_list, treshold_level) result (rv)
        implicit none
        ! dummy
        type(lib_tree_data_element), dimension(:), intent(in) :: element_list
        integer(kind=1) :: treshold_level
        type(lib_tree_data_element), dimension(2**(treshold_level * TREE_DIMENSIONS)) :: rv

        ! auxiliary
        type(lib_tree_universal_index) :: uindex
        integer(kind=4) :: i

        do i=1, size(element_list)
            uindex = lib_tree_hf_get_universal_index(element_list(i)%point_x, treshold_level)

            rv(uindex%n) = element_list(i)
        end do

    end function lib_tree_create_correspondence_vector


    ! ----------------- test functions -----------------
    subroutine lib_tree_test_functions()
        implicit none

        integer :: error_counter

        error_counter = 0

        if (.not. test_lib_tree_create_correspondence_vector()) then
            error_counter = error_counter + 1
        end if

        if (error_counter == 0) then
            print *, "All tests: OK"
        else
            print *, error_counter,"test(s) FAILED"
        end if

        contains

        function test_lib_tree_create_correspondence_vector() result(rv)
            implicit none
            ! dummy
            logical :: rv


            integer(kind=1), parameter :: l_th = 3 ! threshold level
            type(lib_tree_data_element), dimension(3) :: element_list
            type(lib_tree_data_element), dimension(2**(l_th * TREE_DIMENSIONS)) :: correspondece_vector


#if (FMM_DIMENSION == 2)
            element_list(1)%point_x%x(1) = 0.1
            element_list(1)%point_x%x(2) = 0.256

            element_list(2)%point_x%x(1) = 0.4
            element_list(2)%point_x%x(2) = 0.256

            element_list(3)%point_x%x(1) = 0.7
            element_list(3)%point_x%x(2) = 0.56
#endif
            correspondece_vector = lib_tree_create_correspondence_vector(element_list, l_th)

        end function test_lib_tree_create_correspondence_vector

    end subroutine lib_tree_test_functions
end module lib_tree
