
module lib_tree
use lib_tree_helper_functions
use lib_hash_function
    implicit none
    ! Data Structures, Optimal Choice of Parameters, and Complexity Results for Generalized Multilevel Fast Multipole Methods in d Dimensions

    private

    ! parameter
    public :: TREE_DIMENSIONS

    public :: lib_tree_destructor

    ! test functions
    public :: lib_tree_test_functions
    public :: lib_tree_benchmark

    public :: lib_tree_hf_test_functions
    public :: lib_tree_hf_benchmark
    public :: lib_tree_hf_destructor




    ! member
    integer(kind=1), parameter :: CORRESPONDENCE_VECTOR_KIND = 8


    ! type definition
    type lib_tree_data_element
        type(lib_tree_spatial_point) :: point_x
        integer(kind=1) :: element_type
    end type lib_tree_data_element

    type lib_tree_correspondece_vector_element
        integer(kind=CORRESPONDENCE_VECTOR_KIND) :: data_element_number
        integer(kind=1) :: number_of_hash_runs
    end type lib_tree_correspondece_vector_element


    ! module global
    type(lib_tree_data_element), dimension (:), allocatable :: lib_tree_data_element_list
    type(lib_tree_correspondece_vector_element), dimension(:), allocatable :: lib_tree_correspondence_vector

    contains

    ! cleans up the memory
    subroutine lib_tree_destructor()
        if (allocated (lib_tree_correspondence_vector)) then
            deallocate (lib_tree_correspondence_vector)
        end if

        if (allocated (lib_tree_data_element_list)) then
            deallocate (lib_tree_data_element_list)
        end if

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

    function lib_tree_get_threshold_level() result(l_max)
        implicit none

        integer :: i, j
        integer :: s
        integer :: m
        integer :: Bit_max
        integer :: N
        integer :: l_max
        type(lib_tree_universal_index) :: a
        type(lib_tree_universal_index) :: b

        i = 0
        m = s
        do
            i = i + 1
            m = m + 1
!            a = Interleaved(v(ind(i))
!            b = Interleaved(v(ind(m))
            j = Bit_max + 1

            do
                j = j - 1
                a = lib_tree_get_parent(a)
                b = lib_tree_get_parent(b)
                l_max = max(l_max , j)
                if (a%n .eq. b%n) then
                    exit
                end if
            end do

            if (m > N) then
                exit
            end if
        end do

    end function lib_tree_get_threshold_level

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
    !   margin
    !       length of the correspondece vector compared with the length of the element_list
    !       recommented values: 110-125
    !
    ! Returns
    ! ----
    !       the reduced correspondece vector
    !
    !
    ! Correspondence vector
    ! ----
    !   The i-th elements of this vector represents the data point of the i-th box.
    !
    !   Reference: Data_Structures_Optimal_Choice_of_Parameters_and_C: equation (97) and (98)
    !
    !   element_list
    !       d=2, l_th=2
    !       -----------------
    !       | 3 | 5 | 7 | 1 |   // number = universal index
    !       -----------------
    !
    !   complete correspondence vector (can be sparse)
    !       -> dimension = 2**(d*l_th) = 2**(2*2) = 16
    !       -----------------------------------------------------------------
    !       |   | 3 |   | 0 |   | 1 |   | 2 |   |   |   |   |   |   |   |   |   // number = element_list index
    !       -----------------------------------------------------------------
    !         ^                           ^                               ^
    !  index: 0                           7                              16
    !
    !        |  1. calculate the hash of the universal index
    !        |  2. save elements at the position corresponding to the hash value
    !        V
    !   correspondence vector
    !       -------------------------
    !       |   |   |   |   |   |   |
    !       -------------------------
    !         ^                   ^
    !  index: 0                   5     // hash(universal index)
    !
    subroutine lib_tree_create_correspondence_vector(element_list, threshold_level, margin)
        implicit none
        ! dummy
        type(lib_tree_data_element), dimension(:), intent(in) :: element_list
        integer(kind=1), intent(in) :: threshold_level
        integer(kind=2), intent(in) :: margin

        ! auxiliary
        integer(kind=4) :: correspondence_vector_dimension
        integer(kind=1) :: number_of_bits
        integer(kind=CORRESPONDENCE_VECTOR_KIND) :: hashed_uindex
        type(lib_tree_universal_index) :: uindex
        integer(kind=CORRESPONDENCE_VECTOR_KIND) :: i
        integer(kind=1) :: ii
        integer(kind=CORRESPONDENCE_VECTOR_KIND) :: hash_max
        integer(kind=CORRESPONDENCE_VECTOR_KIND) :: hash_idum
        integer(kind=2) :: hash_i

        ! copy element list to the module global *lib_tree_data_element_list* list.
        if (allocated (lib_tree_data_element_list)) then
            deallocate (lib_tree_data_element_list)
        end if
        lib_tree_data_element_list = element_list

        ! check if the significant bits doesn't exceed the number of bits of the universal index
        number_of_bits = threshold_level * TREE_DIMENSIONS
        if (number_of_bits .le. (COORDINATE_BINARY_BYTES * NUMBER_OF_BITS_PER_BYTE)) then
            correspondence_vector_dimension = ceiling(size(element_list) * margin / 100.0)

            hash_max = correspondence_vector_dimension

            ! initiate lib_tree_correspondence_vector
            if (allocated (lib_tree_correspondence_vector)) then
                deallocate (lib_tree_correspondence_vector)
            end if
            allocate( lib_tree_correspondence_vector(correspondence_vector_dimension) )
            lib_tree_correspondence_vector(:)%number_of_hash_runs = 0

            ! iterate element_list
            do i=1, size(lib_tree_data_element_list)
                uindex = lib_tree_hf_get_universal_index(element_list(i)%point_x, threshold_level)

                hash_idum = uindex%n
                do ii=1,1+2,1  !2=Offset zum einschwingen !
                    hash_idum=int(modulo(real(16807.0D0*hash_idum,kind=16),2147483647.0D0),kind=8)
                end do
                hashed_uindex=hashpp(hash_max, hash_idum)+1  ! no zero
                ! find unique hashed universal index
                do ii=1, huge(lib_tree_correspondence_vector(1)%number_of_hash_runs)-1
                    ! save
                    if (lib_tree_correspondence_vector(hashed_uindex)%number_of_hash_runs .eq. 0) then
                        lib_tree_correspondence_vector(hashed_uindex)%data_element_number = i
                        lib_tree_correspondence_vector(hashed_uindex)%number_of_hash_runs = ii
                        ! unique hash found -> terminate the inner do loop immediatly
                        exit
                    end if
                    do hash_i=int(1,2), ii*int(5, 2)
                        hash_idum=int(modulo(real(16807.0D0*hash_idum,kind=16),2147483647.0D0),kind=8)
                    end do
                    hashed_uindex=hashpp(hash_max, hash_idum)+1  ! no zero
                end do
            end do
        else
            print *, "lib_tree_create_correspondence_vector: tree is too deep"
            print *, "   the biggest universal index would exceed (tree_dimension * threshold level * number_of_bits_per_byte)"
            print *, "   number_of_bits = ", number_of_bits
        end if

    end subroutine lib_tree_create_correspondence_vector


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
            integer(kind=2) :: margin

            margin = 125


#if (FMM_DIMENSION == 2)
            element_list(1)%point_x%x(1) = 0.1
            element_list(1)%point_x%x(2) = 0.256

            element_list(2)%point_x%x(1) = 0.4
            element_list(2)%point_x%x(2) = 0.256

            element_list(3)%point_x%x(1) = 0.7
            element_list(3)%point_x%x(2) = 0.56
#endif
            call lib_tree_create_correspondence_vector(element_list, l_th, margin)

            rv = .true.
        end function test_lib_tree_create_correspondence_vector

    end subroutine lib_tree_test_functions

    subroutine lib_tree_benchmark
        implicit none

        call benchmark_lib_tree_create_correspondece_vector()

        contains

        subroutine benchmark_lib_tree_create_correspondece_vector()
            implicit none

            real :: start, finish
            double precision :: delta

            integer(kind=4), parameter :: list_length = 10**5

            integer(kind=1), parameter :: l_th = 16 ! threshold level
            type(lib_tree_data_element), dimension(list_length) :: element_list
            integer(kind=2) :: margin

            integer(kind=4) :: i
            integer(kind=4) :: number

            margin = 110


#if (FMM_DIMENSION == 2)

            do i=1, list_length
                element_list(i)%point_x%x(1) = (0.9 * i)/(1.0*list_length)
                element_list(i)%point_x%x(2) = 0.9 + (-0.9 * i)/list_length
            end do
#endif

            print *, "benchmark_lib_tree_create_correspondece_vector"
            call cpu_time(start)
            call lib_tree_create_correspondence_vector(element_list, l_th, margin)
            call cpu_time(finish)
            print *, "  Time = ", finish-start, " seconds."

            number = 0
            do i=1, list_length
                if (lib_tree_correspondence_vector(i)%number_of_hash_runs .ne. 0) then
                    number = number + 1
                end if
            end do

            if (number .eq. list_length) then
                print *, "  OK"
            else
                print *, "  number is NOT equal to list_length "
                print *, "  number: ", number
                print *, "  list_length: ", list_length
                print *, "    -> ", 1.0*number / list_length, "%"
            end if

        end subroutine benchmark_lib_tree_create_correspondece_vector
    end subroutine lib_tree_benchmark
end module lib_tree
