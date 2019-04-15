#define _FMM_DIMENSION_ 2

! number of bytes of the universal index, value = [4,8,16]
! standard value: 8
#define _UINDEX_BYTES_ 8

! 1: true, 0: false (-> spatial point is real)
#define _SPATIAL_POINT_IS_DOUBLE_ 1

module lib_tree
    !$  use omp_lib
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

    integer(kind=1), parameter :: LIB_TREE_ELEMENT_TYPE_EMPTY = -1


    ! type definition
    type lib_tree_data_element
        type(lib_tree_spatial_point) :: point_x
        type(lib_tree_universal_index) :: uindex
        integer(kind=1) :: element_type
    end type lib_tree_data_element

    type lib_tree_correspondece_vector_element
        integer(kind=CORRESPONDENCE_VECTOR_KIND) :: data_element_number
        integer(kind=2) :: number_of_hash_runs
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
        !$OMP SINGLE
        if (allocated (lib_tree_correspondence_vector)) then
            deallocate (lib_tree_correspondence_vector)
        end if

        if (allocated (lib_tree_data_element_list)) then
            deallocate (lib_tree_data_element_list)
        end if

        call lib_tree_hf_destructor()
        !$OMP END SINGLE
    end subroutine lib_tree_destructor

    !
    ! D = max_d D_d                        (59)
    !
    ! \overline{x} = (x − x_min) / D       (61)
    !
    ! Modification of equation (59) to avoid ones as a result for the
    ! normalised coordinates \overline{x}:
    !
    ! \overline{x} = ((x − x_min) / D) * (1 + 2**(-BITS_MANTISSA))
    !
    ! BITS_MANTISSA
    !   number of bits of the Mantissa (floating point number)
    !
    function lib_tree_get_scaled_element_list(element_list) result(rv)
        implicit none
        ! dummy
        type(lib_tree_data_element), dimension(:), intent(in) :: element_list
        type(lib_tree_data_element), dimension(size(element_list)) :: rv

        ! auxiliary
        type(lib_tree_spatial_point) :: D
        type(lib_tree_spatial_point) :: D_max_buffer
        type(lib_tree_spatial_point) :: x_min
        type(lib_tree_spatial_point) :: x_max

        integer(kind=COORDINATE_BINARY_BYTES) :: i
        integer(kind=1) :: ii

        x_max%x(:) = 0
        x_min%x(:) = 1

        do ii=1, TREE_DIMENSIONS
            x_max%x(ii) = maxval(element_list(:)%point_x%x(ii))
            x_min%x(ii) = minval(element_list(:)%point_x%x(ii))
#if (_SPATIAL_POINT_IS_DOUBLE_ == 1)
            D%x(ii) = (x_max%x(ii) - x_min%x(ii)) * (1 + 2D0**(-52))    ! double: BITS_MANTISSA = 52
#else
            D%x(ii) = (x_max%x(ii) - x_min%x(ii)) * (1 + 2.0**(-23))    ! single: BITS_MANTISSA = 23
#endif
        end do

        D_max_buffer%x(1) = maxval(D%x)

        !$OMP PARALLEL DO PRIVATE(i, ii)
        do i=1, size(element_list)
            do ii=1, TREE_DIMENSIONS
                rv(i)%point_x%x(ii) = (element_list(i)%point_x%x(ii) - x_min%x(ii)) / D_max_buffer%x(1)
            end do
        end do
        !$OMP END PARALLEL DO

        rv(:)%element_type = element_list(:)%element_type

    end function lib_tree_get_scaled_element_list

    ! Arguments
    ! ----
    !   uindex: type(lib_tree_universal_index)
    !       uindex%n
    !           number of the box, with n ranging form 0 to 2^(3*l) in a three-dimensional space
    !       uindex%l
    !           number of the level
    !
    ! Returns
    ! ----
    !   the universal index of the parent box
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
    !   uindex: type(lib_tree_universal_index)
    !       uindex%n
    !           number of the box, with n ranging form 0 to 2^(3*l) in a three-dimensional space
    !       uindex%l
    !           number of the level
    !
    ! Returns
    ! ----
    !   all the univerval indeces of the children's boxes.
    function lib_tree_get_children(uindex) result (rv)
        implicit none
        ! dummy arguments
        type(lib_tree_universal_index), intent (in) :: uindex
        type(lib_tree_universal_index), dimension(2**TREE_DIMENSIONS) :: rv

        rv(:)%n = lib_tree_hf_get_children_all(uindex%n)
        rv(:)%l = uindex%l + int(1,1)

    end function lib_tree_get_children

    ! Arguments
    ! ----
    !   k
    !       k-neighbour of the box (n,l)
    !   uindex
    !       universal index of the box
    !
    ! Returns
    ! ----
    !   all the universal indeces of the k-neighours.
    function lib_tree_get_neighbours(k,uindex) result (rv)
        implicit none
        ! dummy arguments
        type(lib_tree_universal_index), intent (in) :: uindex
        integer(kind=UINDEX_BYTES) :: k
        type(lib_tree_universal_index), dimension(3**TREE_DIMENSIONS-1) :: rv

        rv(:)%n = lib_tree_hf_get_neighbour_all_xD(k, uindex%n, uindex%l)
        rv(:)%l = uindex%l

    end function lib_tree_get_neighbours

    ! Arguments
    ! ----
    !   uindex: type(lib_tree_universal_index)
    !       uindex%n
    !           number of the box, with n ranging form 0 to 2^(3*l) in a three-dimensional space
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
        type(lib_tree_data_element), dimension(:), allocatable :: rv

        ! auxiliar
        type(lib_tree_universal_index), dimension(:), allocatable :: buffer_1
        type(lib_tree_universal_index), dimension(:), allocatable :: buffer_2
        integer(kind=1) :: l_diff
        integer(kind=UINDEX_BYTES) :: i
        integer(kind=UINDEX_BYTES) :: ii
        integer(kind=UINDEX_BYTES) :: index_start
        integer(kind=UINDEX_BYTES) :: index_end

        l_diff = lib_tree_l_th - uindex%l

        allocate (rv(2**(l_diff*TREE_DIMENSIONS)))
        allocate (buffer_1(2**(l_diff*TREE_DIMENSIONS)))
        allocate (buffer_2(2**(l_diff*TREE_DIMENSIONS)))

        ! get all universal indices at the threshold level
        buffer_1(1) = uindex
        do i=1, l_diff
            do ii=1, 2**((i-1)*TREE_DIMENSIONS)
                index_start = (ii-1)*2**TREE_DIMENSIONS+1
                index_end = index_start + 2**TREE_DIMENSIONS-1
                buffer_2(index_start:index_end) = lib_tree_get_children(buffer_1(ii))
            end do
            buffer_1 = buffer_2
        end do
        deallocate (buffer_2)

        ! get all elements
        do i=1, size(buffer_1)
            rv(i) = lib_tree_get_element_from_correspondence_vector(buffer_1(i))
        end do

        deallocate (buffer_1)


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
        type(lib_tree_data_element) :: rv

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
    !   threshold_level
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
        integer(kind=8) :: hashed_uindex
        type(lib_tree_universal_index) :: uindex
        integer(kind=4) :: i
        integer(kind=2) :: ii
        integer(kind=8) :: hash_max
        integer(kind=2) :: number_of_hash_runs
#if (_UINDEX_BYTES_ == 16)
        integer(kind=8), dimension(2) :: hash_idum
        integer(kind=2) :: number_of_bits
#else
        integer(kind=8) :: hash_idum
        integer(kind=1) :: number_of_bits
#endif
        integer(kind=4) :: hash_overflow_counter
        logical :: element_saved

        ! ---- OMP ----
        ! simple implementation of a semaphore
        ! multiple read, single write
        !$  logical :: semaphore_write_correspondence_vector

        ! read and write is initially possilble
        !$  semaphore_write_correspondence_vector = .true.
        ! ~~~~ OMP ~~~~


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

            ! iterate element_list
            lib_tree_max_number_of_hash_runs = 0
            hash_overflow_counter = 0
            !$OMP PARALLEL DO PRIVATE(uindex, hashed_uindex, element_saved, hash_idum , i, ii)
            do i=1, size(lib_tree_data_element_list)
                uindex = lib_tree_hf_get_universal_index(lib_tree_data_element_list(i)%point_x, threshold_level)
                lib_tree_data_element_list(i)%uindex = uindex
                ! find unique hashed universal index
#if (_UINDEX_BYTES_ == 16)
                hashed_uindex = int(1 + hash_kf_16_byte(int(uindex%n,16),int(4,8),hash_max,LIB_TREE_HASH_I,hash_idum),4)
#else
                hashed_uindex = 1 + hash_kf(int(uindex%n,8),int(4,8),hash_max,LIB_TREE_HASH_I,hash_idum)
!                hashed_uindex = int(1 + hash_kf_8_byte(int(uindex%n,8),int(4,4),hash_max,LIB_TREE_HASH_I,hash_idum),4)
!                hashed_uindex = hash_fnv1a_8_byte(uindex%n)
!                hashed_uindex = 1+mod(hashed_uindex, hash_max/2)
#endif

                element_saved = .false.
                do ii=1, LIB_TREE_MAX_HASH_RUNS !huge(lib_tree_correspondence_vector(1)%number_of_hash_runs)-1
                    ! ---- OMP semaphore ----
                    !$  if (semaphore_write_correspondence_vector .eqv. .false.) then
                    !$      ! wait until write process is finished
                    !$      do
                    !$          if (semaphore_write_correspondence_vector .eqv. .true.) then
                    !$              exit
                    !$          end if
                    !$      end do
                    !$  end if
                    ! ~~~~ OMP semaphore ~~~~
                    ! save
                    if (lib_tree_correspondence_vector(hashed_uindex)%number_of_hash_runs .eq. 0) then
                        !$  semaphore_write_correspondence_vector = .false.
                        if ((hashed_uindex .gt. 0) .or. (hashed_uindex .le. hash_max)) then
                            lib_tree_correspondence_vector(hashed_uindex)%data_element_number = i
                            lib_tree_correspondence_vector(hashed_uindex)%number_of_hash_runs = ii

                            if (lib_tree_max_number_of_hash_runs .lt. ii) then
                                lib_tree_max_number_of_hash_runs = ii
                            end if
                            ! unique hash found -> terminate the inner do loop immediatly
                            !$  semaphore_write_correspondence_vector = .true.
                            element_saved = .true.
                            exit
                        end if
                        !$  semaphore_write_correspondence_vector = .true.
                    end if
#if (_UINDEX_BYTES_ == 16)
                    hashed_uindex = int(1 + hashpp_kf_16_byte(hash_max, hash_idum),4)
#else
                    hashed_uindex = 1 + hashpp_kf(hash_max, hash_idum)
!                    hashed_uindex = int(1 + hashpp_kf_8_byte(hash_max, hash_idum),4)
!                    hashed_uindex = hash_fnv1a(hashed_uindex)
!                    hashed_uindex = 1+mod(hashed_uindex, hash_max/2)
#endif
                end do
                if (.not. element_saved) then
                    hash_overflow_counter = hash_overflow_counter + 1
                end if
            end do
            !$OMP END PARALLEL DO
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

    ! Arguments
    ! ----
    !   uindex: type(lib_tree_universal_index)
    !       uindex%n
    !           number of the box, with n ranging form 0 to 2^(3*l) in a three-dimensional space
    !       uindex%l
    !           number of the level, has to be equal to the threshold level (*lib_tree_l_th*)
    !
    function lib_tree_get_element_from_correspondence_vector(uindex) result(rv)
        implicit none
        !dummy
        type(lib_tree_universal_index), intent(in) :: uindex
        type(lib_tree_data_element) :: rv

        ! auxiliary
#if (_UINDEX_BYTES_ == 16)
        integer(kind=8), dimension(2) :: hash_idum
#else
        integer(kind=8) :: hash_idum
#endif
        integer(kind=4) :: hashed_uindex
        integer(kind=2) :: i
        integer(kind=CORRESPONDENCE_VECTOR_KIND) :: element_number
        logical :: element_found

        !$  logical :: opm_end_do_loop

        if (uindex%l .ne. lib_tree_l_th) then
            print *, "lib_tree_get_element_from_correspondence_vector:"
            print *, "    level is NOT equal to the threshold level"
            print *, "    l: ", uindex%l
            print *, "    l_th: ", lib_tree_l_th

            return
        end if

        if (allocated(lib_tree_correspondence_vector)) then
#if (_UINDEX_BYTES_ == 16)
            hashed_uindex = int(1 + hash_kf_16_byte(int(uindex%n,16),int(4,8), &
                int(size(lib_tree_correspondence_vector),8), &
                LIB_TREE_HASH_I,hash_idum), 4)
#else
             hashed_uindex = 1 + hash_kf(int(uindex%n,8),int(4,8), &
                                         int(size(lib_tree_correspondence_vector),8), &
                                         LIB_TREE_HASH_I,hash_idum)
!            hashed_uindex = int(1 + hash_kf_8_byte(int(uindex%n,8),int(4,4), &
!                                            int(size(lib_tree_correspondence_vector),4), &
!                                            LIB_TREE_HASH_I,hash_idum),4)
!            hashed_uindex = hash_fnv1a_8_byte(uindex%n)
!            hashed_uindex = 1+mod(hashed_uindex, hash_max/2)
#endif


            element_found = .false.
!            !$  opm_end_do_loop = .false.
!            !  $   OMP PARALLEL DO PRIVATE(i, hash_idum)
            do i=1, lib_tree_max_number_of_hash_runs
!                !$  if (.not. opm_end_do_loop) then
                if (lib_tree_correspondence_vector(hashed_uindex)%number_of_hash_runs .eq. i) then
!                    !$  opm_end_do_loop = .true.
!                    print *, "element found"
                    element_found = .true.

                    element_number = lib_tree_correspondence_vector(hashed_uindex)%data_element_number
                    rv = lib_tree_data_element_list(element_number)

                    exit
                else
#if (_UINDEX_BYTES_ == 16)
                    hashed_uindex = int(1 + hashpp_kf_16_byte(hash_max, hash_idum), 4)
#else
                    hashed_uindex = 1 + hashpp_kf(hash_max, hash_idum)
!                    hashed_uindex = int(1 + hashpp_kf_8_byte(hash_max, hash_idum), 4)
!                    hashed_uindex = hash_fnv1a(hashed_uindex)
!                    hashed_uindex = 1+mod(hashed_uindex, hash_max/2)
#endif
                end if
!                !$  end if

            end do
!            !   $   OMP END PARALLEL DO
            if (.not. element_found) then
                rv%element_type = LIB_TREE_ELEMENT_TYPE_EMPTY
!                print *, "Element could not be found: n=", n," ..warning"
            end if

        end if


    end function lib_tree_get_element_from_correspondence_vector


    ! ----------------- test functions -----------------
    function lib_tree_test_functions() result(error_counter)
        implicit none

        integer :: error_counter

        !$  print *, "number of threads: ", omp_get_num_threads()
        !$  print *, "max number of threads: ", omp_get_max_threads()

        error_counter = 0

        if (.not. test_lib_tree_create_correspondence_vector()) then
            error_counter = error_counter + 1
        end if
        if (.not. test_lib_tree_get_element_from_correspondence_vector()) then
            error_counter = error_counter + 1
        end if
        if (.not. test_lib_tree_get_domain_e1()) then
            error_counter = error_counter + 1
        end if
        if (.not. test_lib_tree_get_scaled_element_list()) then
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

        function test_lib_tree_get_scaled_element_list() result(rv)
            implicit none
            ! dummy
            logical :: rv

            ! auxiliary
            integer(kind=4), parameter :: list_length = 10**7

            integer(kind=1), parameter :: l_th = 4 ! threshold level
            type(lib_tree_data_element), dimension(list_length) :: element_list
            type(lib_tree_data_element), dimension(list_length) :: element_list_rv

            integer(kind=4) :: i

            ! create test data
#if (_FMM_DIMENSION_ == 2)
            do i=1, list_length
                element_list(i)%point_x%x(1) = i
                element_list(i)%point_x%x(2) = i
            end do
#elif (_FMM_DIMENSION_ == 3)
            do i=1, list_length
                element_list(i)%point_x%x(1) = i
                element_list(i)%point_x%x(2) = i
                element_list(i)%point_x%x(3) = i
            end do
#endif
            element_list(:)%element_type = 1

            ! normalise element data point

            element_list_rv = lib_tree_get_scaled_element_list(element_list)

            print *, "test_lib_tree_get_scaled_element_list"
            print *, "  last x-value: ", element_list_rv(list_length)%point_x%x(1)

            ! todo: evaluate element_list_rv

            rv = .true.

        end function test_lib_tree_get_scaled_element_list

        function test_lib_tree_get_domain_e1() result(rv)
            implicit none
            ! dummy
            logical :: rv

            integer(kind=4), parameter :: list_length = 10**1

            integer(kind=1), parameter :: l_th = 4 ! threshold level
            type(lib_tree_data_element), dimension(list_length) :: element_list
            integer(kind=2) :: margin

            integer(kind=4) :: i
            integer(kind=4) :: number

            type(lib_tree_universal_index) :: uindex
            type(lib_tree_data_element), dimension(:), allocatable :: domain_element_list

            ! ---- create correspondece vector ----

            margin = 300


#if (_FMM_DIMENSION_ == 2)
            do i=1, list_length
                element_list(i)%point_x%x(1) = (0.999 * i)/(1.0*list_length)
                element_list(i)%point_x%x(2) = (0.999 * i)/(1.0*list_length)
            end do
#elif (_FMM_DIMENSION_ == 3)
            do i=1, list_length
                element_list(i)%point_x%x(1) = (0.9 * i)/(1.0*list_length)
                element_list(i)%point_x%x(2) = (0.9 * i)/(1.0*list_length)
                element_list(i)%point_x%x(3) = (0.9 * i)/(1.0*list_length)
            end do
#endif
            element_list(:)%element_type = 1
            call lib_tree_create_correspondence_vector(element_list, l_th, margin)

            number = 0
            do i=1, size(lib_tree_correspondence_vector)
                if (lib_tree_correspondence_vector(i)%number_of_hash_runs .ne. 0) then
                    number = number + 1
                end if
            end do

            if (number .eq. list_length) then
                rv = .true.
            else
                rv = .false.
                print *, "test_lib_tree_get_domain_e1_create_correspondence_vector"
                print *, "  number of data points: ", list_length
                print *, "  max number of hash runs: ", lib_tree_max_number_of_hash_runs
                print *, "  number is NOT equal to list_length "
                print *, "  number: ", number
                print *, "  list_length: ", list_length
                print *, "    -> ", 1.0*number / list_length, "%"
            end if

            ! ---- test get_domain_e1 ----

            uindex%n = 0
            uindex%l = 2

            domain_element_list = lib_tree_get_domain_e1(uindex)

            deallocate (domain_element_list)

        end function test_lib_tree_get_domain_e1

        function test_lib_tree_create_correspondence_vector() result(rv)
            implicit none
            ! dummy
            logical :: rv

            integer(kind=4), parameter :: list_length = 10**4

            integer(kind=1), parameter :: l_th = 16 ! threshold level
            type(lib_tree_data_element), dimension(list_length) :: element_list
            integer(kind=2) :: margin

            integer(kind=4) :: i
            integer(kind=4) :: number

            margin = 300


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


            if (number .eq. list_length) then
                rv = .true.
                print *, "test_lib_tree_create_correspondece_vector: OK"
                print *, "  number of data points: ", list_length
                print *, "  max number of hash runs: ", lib_tree_max_number_of_hash_runs
            else
                rv = .false.
                print *, "test_lib_tree_create_correspondece_vector: FAILED"
                print *, "  number of data points: ", list_length
                print *, "  max number of hash runs: ", lib_tree_max_number_of_hash_runs
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
            integer(kind=4), parameter :: list_length = 10**1

            integer(kind=1), parameter :: l_th = 16 ! threshold level
            type(lib_tree_data_element), dimension(list_length) :: element_list
            integer(kind=2) :: margin

            integer(kind=4) :: i
            integer(kind=4) :: wrong

            margin = 300


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
            wrong = 0
            do i=1, list_length
                uindex = lib_tree_hf_get_universal_index(element_list(i)%point_x, lib_tree_l_th)
                data_element = lib_tree_get_element_from_correspondence_vector(uindex)
                if ((element_list(i)%point_x%x(1) .ne. data_element%point_x%x(1)) .or. &
                    (element_list(i)%point_x%x(2) .ne. data_element%point_x%x(2))) then
                   wrong = wrong + 1
                end if
            end do

            if (wrong .eq. 0) then
                print *, "test_lib_tree_get_element_from_correspondence_vector: ", "OK"
                rv = .true.
            else
                print *, "test_lib_tree_get_element_from_correspondence_vector: ", "FAILED"
                print *, "  wrong elements: ", wrong
                print *, "  wrong elements [%]: ", 100.0 * wrong / list_length
                rv = .false.
            end if

        end function test_lib_tree_get_element_from_correspondence_vector

    end function lib_tree_test_functions

    subroutine lib_tree_benchmark
        implicit none

        call benchmark_lib_tree_create_correspondece_vector()
        call benchmark_lib_tree_get_element_from_correspondence_vector()

    contains

        subroutine benchmark_lib_tree_create_correspondece_vector()
            implicit none

            real :: start, finish
            double precision :: delta

            integer(kind=4), parameter :: list_length = 10**5
            integer(kind=1), parameter :: number_of_runs = 100

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


            call cpu_time(start)
            do i=1, number_of_runs
                call lib_tree_create_correspondence_vector(element_list, l_th, margin)
            end do
            call cpu_time(finish)
            print *, "benchmark_lib_tree_create_correspondece_vector:"
            delta = (finish-start)/number_of_runs
            print *, " create correspondence_vector time: ", delta, " seconds."
            print *, "   total run time: ", finish-start, " seconds"

        end subroutine benchmark_lib_tree_create_correspondece_vector

        subroutine benchmark_lib_tree_get_element_from_correspondence_vector()
            implicit none
            real :: start, finish
            double precision :: delta

            integer(kind=4), parameter :: list_length = 10**5
            integer(kind=4), parameter :: number_of_runs = 10**7

            integer(kind=1), parameter :: l_th = 16 ! threshold level
            type(lib_tree_data_element), dimension(list_length) :: element_list
            type(lib_tree_data_element) :: data_element
            type(lib_tree_universal_index) :: uindex
            integer(kind=2) :: margin

            integer(kind=4) :: i

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

            uindex = lib_tree_hf_get_universal_index(element_list(1)%point_x, lib_tree_l_th)

            call cpu_time(start)
            do i=1, number_of_runs
                data_element = lib_tree_get_element_from_correspondence_vector(uindex)
            end do
            call cpu_time(finish)
            print *, "benchmark_lib_tree_get_element_from_correspondence_vector:"
            delta = (finish-start)/number_of_runs
            print *, " get element from correspondence vector time: ", delta, " seconds."
            print *, "   total run time: ", finish-start, " seconds"

        end subroutine benchmark_lib_tree_get_element_from_correspondence_vector
    end subroutine lib_tree_benchmark
end module lib_tree
