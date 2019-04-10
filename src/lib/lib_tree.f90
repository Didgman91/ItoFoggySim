#define _FMM_DIMENSION_ 3

! number of bytes of the universal index, value = [4,8,16]
! standard value: 8
#define _UINDEX_BYTES_ 8

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
    integer(kind=1), parameter :: CORRESPONDENCE_VECTOR_KIND = 4    ! limited by the total number of elements -> 2**(8*4) = 4,294,967,296 elements
    integer(kind=4), parameter :: LIB_TREE_MAX_HASH_RUNS = 400
    integer(kind=8), parameter :: LIB_TREE_HASH_I= 10



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
    integer(kind=1) :: lib_tree_l_th    ! threshold level
    integer(kind=8) :: hash_max
    integer(kind=2) :: lib_tree_max_number_of_hash_runs

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
        integer(kind=UINDEX_BYTES) :: k
        integer(kind=UINDEX_BYTES), dimension(3**TREE_DIMENSIONS-1) :: rv

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
        integer(kind=UINDEX_BYTES), intent (in) :: k
        type(lib_tree_universal_index), intent(in) :: uindex
        type(lib_tree_spatial_point) :: rv

        ! auxiliar
!        type(lib_tree_universal_index), dimension(:), allocatable :: domain_box
        integer(kind=UINDEX_BYTES) :: number_of_boxes
        integer(kind=UINDEX_BYTES), dimension(3**TREE_DIMENSIONS-1, k) :: buffer

        integer(kind=UINDEX_BYTES) :: i

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
        integer(kind=4) :: hashed_uindex
        type(lib_tree_universal_index) :: uindex
        integer(kind=4) :: i
        integer(kind=2) :: ii
        integer(kind=8) :: hash_max
#if (_UINDEX_BYTES_ == 16)
        integer(kind=8), dimension(2) :: hash_idum
        integer(kind=2) :: number_of_bits
#else
        integer(kind=8) :: hash_idum
        integer(kind=1) :: number_of_bits
#endif
        integer(kind=4) :: hash_overflow_counter
        logical :: element_saved

        ! copy element list to the module global *lib_tree_data_element_list* list.

        if ( .not. allocated(lib_tree_data_element_list) ) then
            allocate( lib_tree_data_element_list(size(element_list)) )
        else
            if ( size(lib_tree_data_element_list) .ne. size(element_list)) then
                deallocate (lib_tree_data_element_list)
                allocate( lib_tree_data_element_list(size(element_list)) )
            end if
        end if
        lib_tree_data_element_list = element_list
        lib_tree_l_th = threshold_level

        ! check if the significant bits doesn't exceed the number of bits of the universal index
        number_of_bits = threshold_level * TREE_DIMENSIONS
#if (_UINDEX_BYTES_ == 16)
        if (number_of_bits .le. (int(UINDEX_BYTES,2) * int(NUMBER_OF_BITS_PER_BYTE,2))) then
#else
        if (number_of_bits .le. (UINDEX_BYTES * NUMBER_OF_BITS_PER_BYTE)) then
#endif
            correspondence_vector_dimension = ceiling(size(element_list) * margin / 100.0)

            hash_max = correspondence_vector_dimension

            ! initiate lib_tree_correspondence_vector
            if (.not. allocated (lib_tree_correspondence_vector)) then
                allocate( lib_tree_correspondence_vector(correspondence_vector_dimension) )
            else
                if (size(lib_tree_correspondence_vector) .ne. correspondence_vector_dimension) then
                    deallocate (lib_tree_correspondence_vector)
                    allocate( lib_tree_correspondence_vector(correspondence_vector_dimension) )
                end if
            end if
            lib_tree_correspondence_vector(:)%number_of_hash_runs = 0

            ii = 0
            do i=1, correspondence_vector_dimension
                if (lib_tree_correspondence_vector(i)%number_of_hash_runs .ne. 0) then
                    ii = ii + int(1,2)
                end if
            end do

            ! iterate element_list
            lib_tree_max_number_of_hash_runs = 0
            hash_overflow_counter = 0
            do i=1, size(lib_tree_data_element_list)
                uindex = lib_tree_hf_get_universal_index(element_list(i)%point_x, threshold_level)

                ! find unique hashed universal index
#if (_UINDEX_BYTES_ == 16)
                hashed_uindex = hash_kf_16_byte(int(uindex%n,16),int(4,8),hash_max,LIB_TREE_HASH_I,hash_idum)
#else
                hashed_uindex = hash_kf(int(uindex%n,8),int(4,8),hash_max,LIB_TREE_HASH_I,hash_idum)
#endif

                element_saved = .false.
                do ii=1, LIB_TREE_MAX_HASH_RUNS !huge(lib_tree_correspondence_vector(1)%number_of_hash_runs)-1
                    ! save
                    if (lib_tree_correspondence_vector(hashed_uindex)%number_of_hash_runs .eq. 0) then
                        if ((hashed_uindex .gt. 0) .or. (hashed_uindex .le. hash_max)) then
                            lib_tree_correspondence_vector(hashed_uindex)%data_element_number = i
                            lib_tree_correspondence_vector(hashed_uindex)%number_of_hash_runs = ii

                            if (lib_tree_max_number_of_hash_runs .lt. ii) then
                                lib_tree_max_number_of_hash_runs = ii
                            end if
                            ! unique hash found -> terminate the inner do loop immediatly
                            element_saved = .true.
                            exit
                        end if
                    end if
#if (_UINDEX_BYTES_ == 16)
                    hashed_uindex = hashpp_kf_16_byte(hash_max, hash_idum)
#else
                    hashed_uindex = hashpp_kf(hash_max, hash_idum)
#endif
                end do
                if (.not. element_saved) then
                    hash_overflow_counter = hash_overflow_counter + 1
                end if
            end do
            if (hash_overflow_counter .ne. 0) then
                print *, "lib_tree_create_correspondence_vector  ..ERROR"
                print *, "    number of hash overflows : ", hash_overflow_counter
            end if
        else
            print *, "lib_tree_create_correspondence_vector: tree is too deep"
            print *, "   the biggest universal index would exceed (tree_dimension * threshold level * number_of_bits_per_byte)"
            print *, "   number_of_bits = ", number_of_bits
        end if

    end subroutine lib_tree_create_correspondence_vector

    function lib_tree_get_element_from_correspondence_vector(n) result(rv)
        implicit none
        !dummy
        integer(kind=UINDEX_BYTES), intent(in) :: n
        type(lib_tree_data_element) :: rv

        ! auxiliary
        integer(kind=8) :: hash_idum
        integer(kind=4) :: hashed_uindex
        integer(kind=2) :: i
        integer(kind=CORRESPONDENCE_VECTOR_KIND) :: element_number
        logical :: element_found

        if (allocated(lib_tree_correspondence_vector)) then
            hashed_uindex = hash_kf(int(n,8),int(4,8),int(size(lib_tree_correspondence_vector),8),LIB_TREE_HASH_I,hash_idum)

            element_found = .false.
            do i=1, lib_tree_max_number_of_hash_runs
                if (lib_tree_correspondence_vector(hashed_uindex)%number_of_hash_runs .eq. i) then
                    element_number = lib_tree_correspondence_vector(hashed_uindex)%data_element_number
                    rv = lib_tree_data_element_list(element_number)

                    element_found = .true.
                else
                    hashed_uindex = hashpp_kf(hash_max, hash_idum)
                end if

            end do
            if (.not. element_found) then
                print *, "Element could not be found ..ERROR"
            end if

        end if


    end function lib_tree_get_element_from_correspondence_vector


    ! ----------------- test functions -----------------
    function lib_tree_test_functions() result(error_counter)
        implicit none

        integer :: error_counter

        error_counter = 0

        if (.not. test_lib_tree_create_correspondence_vector()) then
            error_counter = error_counter + 1
        end if
        if (.not. test_lib_tree_get_element_from_correspondence_vector()) then
            error_counter = error_counter + 1
        end if

        print *, "-------------lib_tree_test_functions----------------"
        if (error_counter == 0) then
            print *, "lib_tree_test_functions tests: OK"
        else
            print *, error_counter,"lib_tree_test_functions test(s) FAILED"
        end if
        print *, "----------------------------------------------------"

        contains

        function test_lib_tree_create_correspondence_vector() result(rv)
            implicit none
            ! dummy
            logical :: rv

            integer(kind=4), parameter :: list_length = 10**5

            integer(kind=1), parameter :: l_th = 16 ! threshold level
            type(lib_tree_data_element), dimension(list_length) :: element_list
            integer(kind=2) :: margin

            integer(kind=4) :: i
            integer(kind=4) :: number

            margin = 400


#if (_FMM_DIMENSION_ == 2)
            do i=1, list_length
                element_list(i)%point_x%x(1) = (0.999 * i)/(1.0*list_length)
                element_list(i)%point_x%x(2) = 0.999 + (-0.999 * i)/(1.0*list_length)
            end do
#elif (_FMM_DIMENSION_ == 3)
            do i=1, list_length
                element_list(i)%point_x%x(1) = (0.9 * i)/(1.0*list_length)
                element_list(i)%point_x%x(2) = 0.9 + (-0.9 * i)/list_length
                element_list(i)%point_x%x(3) = 0.9 + (-0.9 * i)/list_length
            end do
#endif
            call lib_tree_create_correspondence_vector(element_list, l_th, margin)

            number = 0
            do i=1, size(lib_tree_correspondence_vector)
                if (lib_tree_correspondence_vector(i)%number_of_hash_runs .ne. 0) then
                    number = number + 1
                end if
            end do


            print *, "test_lib_tree_create_correspondece_vector"
            print *, "  number of data points: ", list_length
            print *, "  max number of hash runs: ", lib_tree_max_number_of_hash_runs

            if (number .eq. list_length) then
                rv = .true.
            else
                rv = .false.
                print *, "  number is NOT equal to list_length "
                print *, "  number: ", number
                print *, "  list_length: ", list_length
                print *, "    -> ", 1.0*number / list_length, "%"
            end if

        end function test_lib_tree_create_correspondence_vector

        function test_lib_tree_get_element_from_correspondence_vector() result (rv)
            implicit none
            ! dummy
            logical :: rv

            ! auxiliary
            type(lib_tree_universal_index) :: uindex
            type(lib_tree_data_element) :: data_element

            ! generate dataset
            integer(kind=4), parameter :: list_length = 10**5

            integer(kind=1), parameter :: l_th = 16 ! threshold level
            type(lib_tree_data_element), dimension(list_length) :: element_list
            integer(kind=2) :: margin

            integer(kind=4) :: i
            integer(kind=4) :: number

            margin = 400


#if (_FMM_DIMENSION_ == 2)
            do i=1, list_length
                element_list(i)%point_x%x(1) = (0.999 * i)/(1.0*list_length)
                element_list(i)%point_x%x(2) = 0.999 + (-0.999 * i)/(1.0*list_length)
            end do
#elif (_FMM_DIMENSION_ == 3)
            do i=1, list_length
                element_list(i)%point_x%x(1) = (0.9 * i)/(1.0*list_length)
                element_list(i)%point_x%x(2) = 0.9 + (-0.9 * i)/list_length
                element_list(i)%point_x%x(3) = 0.9 + (-0.9 * i)/list_length
            end do
#endif
            call lib_tree_create_correspondence_vector(element_list, l_th, margin)


            ! test dataset
            uindex = lib_tree_hf_get_universal_index(element_list(1)%point_x, lib_tree_l_th)
            data_element = lib_tree_get_element_from_correspondence_vector(uindex%n)

            if (element_list(1)%point_x%x(1) .eq. data_element%point_x%x(1)) then
                if (element_list(1)%point_x%x(2) .eq. data_element%point_x%x(2)) then
                    print *, "test_lib_tree_get_element_from_correspondence_vector: ", "OK"
                    rv = .true.
                else
                    print *, "test_lib_tree_get_element_from_correspondence_vector: ", "FAILED"
                    rv = .false.
                end if
            else
                print *, "test_lib_tree_get_element_from_correspondence_vector: ", "FAILED"
                rv = .false.
            end if

        end function test_lib_tree_get_element_from_correspondence_vector

    end function lib_tree_test_functions

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

            margin = 125


#if (_FMM_DIMENSION_ == 2)
            do i=1, list_length
                element_list(i)%point_x%x(1) = (0.999 * i)/(1.0*list_length)
                element_list(i)%point_x%x(2) = 0.999 + (-0.999 * i)/(1.0*list_length)
            end do
#elif (_FMM_DIMENSION_ == 3)
            do i=1, list_length
                element_list(i)%point_x%x(1) = (0.9 * i)/(1.0*list_length)
                element_list(i)%point_x%x(2) = 0.9 + (-0.9 * i)/list_length
                element_list(i)%point_x%x(3) = 0.9 + (-0.9 * i)/list_length
            end do
#endif


            call cpu_time(start)
            call lib_tree_create_correspondence_vector(element_list, l_th, margin)
            call cpu_time(finish)
            print *, "benchmark_lib_tree_create_correspondece_vector:"
            print *, " create correspondence_vector time: ", finish-start, " seconds."

            number = 0
            do i=1, size(lib_tree_correspondence_vector)
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
