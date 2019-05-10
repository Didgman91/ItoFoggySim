#define _DEBUG_

! LIB: Mulitlevel Fast Multipole Method - Helper Functions
!
module lib_ml_fmm_helper_functions
    use lib_hash_function
    use lib_tree_public
    use ml_fmm_type
    implicit none

    private

    ! --- public functions ---
    public :: lib_ml_fmm_hf_get_neighbourhood_size
    public :: lib_ml_fmm_hf_create_hierarchy

    public :: lib_ml_fmm_hf_test_functions

    ! --- parameter ---
    integer(kind=2), public, parameter :: LIB_ML_FMM_HF_HIERARCHY_MARGIN = 200

    integer(kind=1), parameter :: IGNORE_ENTRY = -1
    integer(kind=2), parameter :: MAXIMUM_NUMBER_OF_HASH_RUNS = 200

    ! --- type definitions ---
    type hash_list_type
        integer(kind=UINDEX_BYTES) :: value
        integer(kind=2) :: hash_runs
    end type

    contains

        ! "The S-expansion (10) near the center of the nth box at level l for
        !  x_i ∈ E_1 (n,l) is valid for any y in the domain E_3 (n,l)."
        !
        ! Equation
        !   k >= 1/2 (R_c * d^(1/2) - 1)   (25)
        !
        ! Reference: Data_Structures_Optimal_Choice_of_Parameters_and_C
        !
        ! Arguments
        ! ----
        !   R_c: double precision
        !
        !
        ! Returns
        ! ----
        !   k: integer
        !       the minimum number of neighbours
        !
        !
        function lib_ml_fmm_hf_get_neighbourhood_size_S(R_c) result (k)
            implicit none

            double precision, intent (in) :: R_c
            double precision :: buffer
            integer(kind=UINDEX_BYTES) :: k

            buffer = 0.5 * ( R_c * sqrt(real(TREE_DIMENSIONS)) - 1 )

            k = ceiling(buffer)

        end function lib_ml_fmm_hf_get_neighbourhood_size_S

        ! "The R-expansion (8) near the center of the nth box at level l for x i ∈ E_3 (n,l) 
        ! is valid for any y from the domain E_1 (n,l)."
        !
        ! Equation
        !   k >= 1/2 (1/r_c * d^(1/2) - 1)   (26)
        !
        ! Reference: Data_Structures_Optimal_Choice_of_Parameters_and_C
        !
        ! Arguments
        ! ----
        !   R_c: double precision
        !       number of the node, with n ranging form 0 to 2^(3*l) in a three-dimensional space
        !
        ! Returns
        ! ----
        !   k: integer
        !       the minimum number of neighbours
        !
        !
        function lib_ml_fmm_hf_get_neighbourhood_size_R(r_c) result (k)
            implicit none

            double precision, intent (in) :: r_c
            double precision :: buffer
            integer(kind=UINDEX_BYTES) :: k

            buffer = 0.5 * (1/r_c * sqrt(real(TREE_DIMENSIONS)) -1 )

            k = ceiling(buffer)

        end function lib_ml_fmm_hf_get_neighbourhood_size_R

        ! "The R-expansion (8) near the center of the nth box at level l for x i ∈ E_3 (n,l) 
        ! is valid for any y from the domain E_1 (n,l)."
        !
        ! Equation
        !   k >= 1/2 (max(1/r_c, R_c) * d^(1/2) - 1)   (26)
        !
        ! Reference: Data_Structures_Optimal_Choice_of_Parameters_and_C
        !
        ! Arguments
        ! ----
        !   R_c: double precision
        !       number of the node, with n ranging form 0 to 2^(3*l) in a three-dimensional space
        !
        ! Returns
        ! ----
        !   k: integer
        !       the minimum number of neighbours
        !
        !
        function lib_ml_fmm_hf_get_neighbourhood_size(R_c1, r_c2) result (k)
            implicit none

            double precision, intent (in) :: R_c1
            double precision, intent (in) :: r_c2
            double precision :: buffer
            integer(kind=UINDEX_BYTES) :: k

            buffer = 0.5 * (max(1/r_c2, R_c1) * sqrt(real(TREE_DIMENSIONS)) -1 )

            k = ceiling(buffer)

        end function lib_ml_fmm_hf_get_neighbourhood_size

        ! Calculates if the R expansino is valid
        !
        ! Formula
        !   |y − x_∗| .le. r c |x_i − x_∗|       (8)
        !
        !   Reference: Data_Structures_Optimal_Choice_of_Parameters_and_C
        !
        ! Arguments
        ! ----
        !   y
        !   x_i
        !   x
        !   r_c
        !
        ! Returns
        ! ----
        !   rv: logical
        !
        function lib_ml_fmm_hf_check_validity_expansion_R(y, x_i, x, r_c) result(rv)
            implicit none
            ! dummy
            type(lib_tree_spatial_point), intent(in) :: y
            type(lib_tree_spatial_point), intent(in) :: x_i
            type(lib_tree_spatial_point), intent(in) :: x
            real, intent(in) :: r_c
            logical :: rv

            if (abs(y - x) .le. r_c * abs(x_i - x)) then
                rv = .true.
            else
                rv = .false.
            end if

        end function lib_ml_fmm_hf_check_validity_expansion_R

        ! Calculates if the S expansino is valid
        !
        ! Formula
        !   |y − x_∗| .ge. R_c |x_i − x_∗|       (10)
        !
        !   Reference: Data_Structures_Optimal_Choice_of_Parameters_and_C
        !
        ! Arguments
        ! ----
        !   y
        !   x_i
        !   x
        !   R_c
        !
        ! Returns
        ! ----
        !   rv: logical
        !
        function lib_ml_fmm_hf_check_validity_expansion_S(y, x_i, x, R_c) result(rv)
            implicit none
            ! dummy
            type(lib_tree_spatial_point), intent(in) :: y
            type(lib_tree_spatial_point), intent(in) :: x_i
            type(lib_tree_spatial_point), intent(in) :: x
            real, intent(in) :: R_c
            logical :: rv

            if (abs(y - x) .ge. R_c * abs(x_i - x)) then
                rv = .true.
            else
                rv = .false.
            end if

        end function lib_ml_fmm_hf_check_validity_expansion_S

        ! Calculates if the RR translation is valid
        !
        ! Formula
        !   |y − x_∗2 | .le. r_c |x_i − x_∗1 | − |x_∗1 − x_∗2 |       (12)
        !
        !   Reference: Data_Structures_Optimal_Choice_of_Parameters_and_C
        !
        ! Arguments
        ! ----
        !   y
        !   x_i
        !   x_1: spatial point
        !       origin of coordinate system 1
        !   x_2: spatial point
        !       origin of coordinate system 2
        !   r_c
        !
        ! Returns
        ! ----
        !   rv: logical
        !
        function lib_ml_fmm_hf_check_validity_translation_RR(y, x_i, x_1, x_2, r_c) result(rv)
            implicit none
            ! dummy
            type(lib_tree_spatial_point), intent(in) :: y
            type(lib_tree_spatial_point), intent(in) :: x_i
            type(lib_tree_spatial_point), intent(in) :: x_1
            type(lib_tree_spatial_point), intent(in) :: x_2
            real, intent(in) :: r_c
            logical :: rv

            if (abs(y - x_2) .le. (r_c * abs(x_i - x_1) - abs(x_1 - x_2))) then
                rv = .true.
            else
                rv = .false.
            end if

        end function lib_ml_fmm_hf_check_validity_translation_RR

        ! Calculates if the SR translation is valid
        !
        ! Formula
        !   |y − x_∗2 | .le. min(|x_∗2 − x_∗1 | − R_c |x_i − x_∗1 | , r_c |x_i − x_∗2 |)       (14)
        !
        !   Reference: Data_Structures_Optimal_Choice_of_Parameters_and_C
        !
        ! Arguments
        ! ----
        !   y
        !   x_i
        !   x_1: spatial point
        !       origin of coordinate system 1
        !   x_2: spatial point
        !       origin of coordinate system 2
        !   R_c1
        !       (the number doesn't refere to a coordinate system, Fortran is just not case-sensitive)
        !   r_c2
        !       (the number doesn't refere to a coordinate system, Fortran is just not case-sensitive)
        !
        ! Returns
        ! ----
        !   rv: logical
        !
        function lib_ml_fmm_hf_check_validity_translation_SR(y, x_i, x_1, x_2, R_c1, r_c2) result(rv)
            implicit none
            ! dummy
            type(lib_tree_spatial_point), intent(in) :: y
            type(lib_tree_spatial_point), intent(in) :: x_i
            type(lib_tree_spatial_point), intent(in) :: x_1
            type(lib_tree_spatial_point), intent(in) :: x_2
            real, intent(in) :: R_c1
            real, intent(in) :: r_c2
            logical :: rv

            if (abs(y - x_2) .le. min(abs(x_2 - x_1) - R_c1 * abs(x_i - x_1), r_c2 * abs(x_i - x_2))) then
                rv = .true.
            else
                rv = .false.
            end if

        end function lib_ml_fmm_hf_check_validity_translation_SR

        ! Calculates if the SR translation is valid
        !
        ! Formula
        !   |y − x_∗2 | > R_c |x_∗2 − x_∗1 | + R_c |x_i − x_∗1 |       (16)
        !
        !   Reference: Data_Structures_Optimal_Choice_of_Parameters_and_C
        !
        ! Arguments
        ! ----
        !   y
        !   x_i
        !   x_1: spatial point
        !       origin of coordinate system 1
        !   x_2: spatial point
        !       origin of coordinate system 2
        !   R_c
        !
        ! Returns
        ! ----
        !   rv: logical
        !
        function lib_ml_fmm_hf_check_validity_translation_SS(y, x_i, x_1, x_2, R_c) result(rv)
            implicit none
            ! dummy
            type(lib_tree_spatial_point), intent(in) :: y
            type(lib_tree_spatial_point), intent(in) :: x_i
            type(lib_tree_spatial_point), intent(in) :: x_1
            type(lib_tree_spatial_point), intent(in) :: x_2
            real, intent(in) :: R_c
            logical :: rv

            if (abs(y - x_2) .gt. R_c * abs(x_2 - x_1) + R_c * abs(x_i - x_1)) then
                rv = .true.
            else
                rv = .false.
            end if

        end function lib_ml_fmm_hf_check_validity_translation_SS

        ! Argument
        ! ----
        !   data: array<lib_tree_data_element>
        !       list of data points
        !
        ! Result
        ! ----
        !   hierarchy
        !       hierarchy of the X- and Y-hierarchy (sources and targets)
        function lib_ml_fmm_hf_create_hierarchy(data_elements, correspondence_vector, length, l_min, l_max) result(hierarchy)
            implicit none
            ! dummy
            type(lib_tree_data_element), dimension(:), intent(inout) :: data_elements
            type(lib_tree_correspondece_vector_element), dimension(:), allocatable, intent(inout) :: correspondence_vector
            integer(kind=UINDEX_BYTES), dimension(3), intent(in) :: length
            integer(kind=1), intent(in) :: l_max
            integer(kind=1), intent(in) :: l_min
            type(lib_ml_fmm_hierarchy), dimension(:), allocatable :: hierarchy

            ! auxiliary
            integer(kind=1) :: buffer_l
            integer(kind=UINDEX_BYTES) :: i
            type(lib_tree_universal_index) :: buffer_uindex
            type(lib_tree_universal_index), dimension(:), allocatable :: uindex_list_X
            type(lib_tree_universal_index), dimension(:), allocatable :: uindex_list_Y
            type(lib_tree_universal_index), dimension(:), allocatable :: uindex_list_XY
            integer(kind=1) :: l_th
            integer(kind=1) :: hierarchy_type
            integer(kind=UINDEX_BYTES), dimension(3) :: uindex_list_counter
            integer(kind=UINDEX_BYTES) :: number_of_boxes_at_l_max

            integer(kind=CORRESPONDENCE_VECTOR_KIND) :: element_index

            allocate( hierarchy(l_min:l_max))

            allocate( uindex_list_X(length(HIERARCHY_X)) )
            allocate( uindex_list_Y(length(HIERARCHY_Y)) )
            allocate( uindex_list_XY(length(HIERARCHY_XY)) )
            uindex_list_counter(:) = 0

            ! iterate over all correspondence vector elements (l = l_th)
            do i=1, size(correspondence_vector)
                if (correspondence_vector(i)%number_of_hash_runs .gt. 0) then
                    element_index = correspondence_vector(i)%data_element_number(1)
                    hierarchy_type = data_elements(element_index)%hierarchy
                    buffer_uindex = data_elements(element_index)%uindex
                    if (hierarchy_type .eq. HIERARCHY_X) then
                        uindex_list_counter(HIERARCHY_X) = uindex_list_counter(HIERARCHY_X) + 1
                        uindex_list_X(uindex_list_counter(HIERARCHY_X)) = buffer_uindex
                    else if (hierarchy_type .eq. HIERARCHY_Y) then
                        uindex_list_counter(HIERARCHY_Y) = uindex_list_counter(HIERARCHY_Y) + 1
                        uindex_list_Y(uindex_list_counter(HIERARCHY_Y)) = buffer_uindex
                    else if (hierarchy_type .eq. HIERARCHY_XY) then
                        uindex_list_counter(HIERARCHY_XY) = uindex_list_counter(HIERARCHY_XY) + 1
                        uindex_list_XY(uindex_list_counter(HIERARCHY_XY)) = buffer_uindex
                    else
                        print *, "lib_ml_fmm_hf_create_hierarchy: ERROR"
                        print *, "  hierarchy not defined"
                    end if
                end if
            end do

            l_th = buffer_uindex%l

            ! if necessary, get indices at level l_max
            if (l_th .gt. l_max) then
                call lib_ml_fmm_hf_get_parent_uindex_lists(uindex_list_X, l_th-l_max)
                call lib_ml_fmm_hf_get_parent_uindex_lists(uindex_list_Y, l_th-l_max)
                call lib_ml_fmm_hf_get_parent_uindex_lists(uindex_list_XY, l_th-l_max)
                uindex_list_X = lib_ml_fmm_hf_make_uindex_list_unique(uindex_list_X)
                uindex_list_Y = lib_ml_fmm_hf_make_uindex_list_unique(uindex_list_Y)
                uindex_list_XY = lib_ml_fmm_hf_make_uindex_list_unique(uindex_list_XY)
            end if

            ! setup hierarcy at level l_max
            call lib_ml_fmm_hf_add_uindex_to_hierarchy(hierarchy, l_max, uindex_list_X, uindex_list_Y, uindex_list_XY)

            ! setup hierary up to the level l_min
            if (l_max .gt. l_min) then
                do i=l_max-1, l_min, -1
                    call lib_ml_fmm_hf_get_parent_uindex_lists(uindex_list_X)
                    call lib_ml_fmm_hf_get_parent_uindex_lists(uindex_list_Y)
                    call lib_ml_fmm_hf_get_parent_uindex_lists(uindex_list_XY)
                    uindex_list_X = lib_ml_fmm_hf_make_uindex_list_unique(uindex_list_X)
                    uindex_list_Y = lib_ml_fmm_hf_make_uindex_list_unique(uindex_list_Y)
                    uindex_list_XY = lib_ml_fmm_hf_make_uindex_list_unique(uindex_list_XY)
                    call lib_ml_fmm_hf_add_uindex_to_hierarchy(hierarchy, int(i, 1), uindex_list_X, uindex_list_Y, uindex_list_XY)
                end do
            end if
        end function

        ! Adds lists of universal indices of different hierachical type to the hierarchy.
        !
        ! Argument
        ! ----
        !   hierarchy: array<lib_ml_fmm_hierarchy>, inout
        !       reference to the hierarchy object to be processed
        !   level: integer, in
        !       level of the hierarchy to be processed
        !   uindex_list_[X,Y,XY]: array<lib_tree_universal_index>, inout (read-only)
        !       list of uindex' of the type [X,Y,XY] to add to the hierarchy
        !
        subroutine lib_ml_fmm_hf_add_uindex_to_hierarchy(hierarchy, level, uindex_list_X, uindex_list_Y, uindex_list_XY)
            implicit none
            ! dummy
            type(lib_ml_fmm_hierarchy), dimension(:), allocatable, intent(inout) :: hierarchy
            integer(kind=1), intent(in) :: level
            type(lib_tree_universal_index), dimension(:), intent(inout) :: uindex_list_X
            type(lib_tree_universal_index), dimension(:), intent(inout) :: uindex_list_Y
            type(lib_tree_universal_index), dimension(:), intent(inout) :: uindex_list_XY

            !auxiliary
            integer(kind=UINDEX_BYTES) :: number_of_boxes
            integer(kind=UINDEX_BYTES) :: number_of_entries_log_2
            integer(kind=UINDEX_BYTES) :: counter
            integer(kind=UINDEX_BYTES) :: i
            integer(kind=2) :: ii
            integer(kind=UINDEX_BYTES) :: hash
            integer(kind=UINDEX_BYTES) :: max_value
            integer(kind=UINDEX_BYTES) :: n
            logical :: no_hash

            number_of_boxes = size(uindex_list_X) &
                              + size(uindex_list_Y) &
                              + size(uindex_list_XY)
            hierarchy(level)%number_of_boxes = number_of_boxes

            allocate (hierarchy(level)%coefficient_list(number_of_boxes))
            allocate (hierarchy(level)%hierarchy_type(number_of_boxes))
            allocate (hierarchy(level)%coefficient_type(number_of_boxes))

            hierarchy(level)%hierarchy_type(:) = IGNORE_ENTRY

            ! 2**(dl) < 256
            if (log(real(level)) / log(2.0) .lt. 8) then
                no_hash = .true.
            else
                no_hash = .false.
            end if
            no_hash = .false.

            ! calculates whether hash access or direct access via the universal index requires less memory
            !
            ! number_of_boxes * margin/100.0 < 2**(dimension * level)
            number_of_entries_log_2 = ceiling(log(number_of_boxes * LIB_ML_FMM_HF_HIERARCHY_MARGIN/100.0) / log(2.0))
            if ((number_of_entries_log_2 .lt. TREE_DIMENSIONS * level) .and. &
                .not. no_hash) then
                hierarchy(level)%is_hashed = .true.
                max_value = ceiling(number_of_boxes * LIB_ML_FMM_HF_HIERARCHY_MARGIN/100.0)
                allocate (hierarchy(level)%hashed_coefficient_list_index(max_value))
                hierarchy(level)%hashed_coefficient_list_index(:)%number_of_hash_runs = IGNORE_ENTRY
                hierarchy(level)%maximum_number_of_hash_runs = 0
                allocate (hierarchy(level)%coefficient_list_index(number_of_boxes))
                hierarchy(level)%coefficient_list_index(:) = IGNORE_ENTRY

                counter = 0
                call lib_ml_fmm_hf_add_uindex_to_hierarchy_hashed(hierarchy, level, uindex_list_X, HIERARCHY_X, counter)
                call lib_ml_fmm_hf_add_uindex_to_hierarchy_hashed(hierarchy, level, uindex_list_Y, HIERARCHY_Y, counter)
                call lib_ml_fmm_hf_add_uindex_to_hierarchy_hashed(hierarchy, level, uindex_list_XY, HIERARCHY_XY, counter)
             else
                hierarchy(level)%is_hashed = .false.
                i = 2**(TREE_DIMENSIONS * level)
                allocate (hierarchy(level)%coefficient_list_index(i))

                hierarchy(level)%coefficient_list_index(:) = IGNORE_ENTRY

                ! setup the arrays coeeficient_list_index and hierarchy_type
                counter = 0

                call lib_ml_fmm_hf_add_uindex_to_hierarchy_lookup_table(hierarchy, level, uindex_list_X, HIERARCHY_X, counter)
                call lib_ml_fmm_hf_add_uindex_to_hierarchy_lookup_table(hierarchy, level, uindex_list_Y, HIERARCHY_Y, counter)
                call lib_ml_fmm_hf_add_uindex_to_hierarchy_lookup_table(hierarchy, level, uindex_list_XY, HIERARCHY_XY, counter)
            end if

        end subroutine lib_ml_fmm_hf_add_uindex_to_hierarchy

        ! Add the uindex to the hierarchy by using a hash function
        !
        ! Arguments
        ! ----
        !   hierarchy: array<lib_ml_fmm_hierarchy>, inout
        !       reference to the hierarchy object to be processed
        !   level: integer, in
        !       level of the hierarchy to be processed
        !   uindex_list: array<lib_tree_universal_index>, inout (read-only)
        !       list of uindex' to add to the hierarchy
        !   hierarchy_type: integer, in
        !       type of the uindex_list [HIERARCHY_X, HIERARCHY_Y, HIERARCHY_XY]
        !   counter: integer, inout
        !       counter of added uindex' to the hierarchy()
        !
        subroutine lib_ml_fmm_hf_add_uindex_to_hierarchy_hashed(hierarchy, level, uindex_list, hierarchy_type, counter)
            implicit none
            ! dummy
            type(lib_ml_fmm_hierarchy), dimension(:), allocatable, intent(inout) :: hierarchy
            integer(kind=1), intent(in) :: level
            type(lib_tree_universal_index), dimension(:), intent(inout) :: uindex_list
            integer(kind=1), intent(in) :: hierarchy_type
            integer(kind=UINDEX_BYTES) :: counter

            ! auxiliary
            integer(kind=UINDEX_BYTES) :: i
            integer(kind=2) :: ii
            integer(kind=UINDEX_BYTES) :: hash
            integer(kind=UINDEX_BYTES) :: max_value

            max_value = size(hierarchy(level)%hashed_coefficient_list_index)

            ! interate over all indices
            do i=1, size(uindex_list)
                ! determine a valid hash
                hash = hash_fnv1a(uindex_list(i)%n, max_value)
                do ii=1, MAXIMUM_NUMBER_OF_HASH_RUNS
                    if (hierarchy(level)%hashed_coefficient_list_index(hash)%number_of_hash_runs .eq. IGNORE_ENTRY) then
                        ! saving the uindex in the hierarchy
                        counter = counter + 1
                        hierarchy(level)%hashed_coefficient_list_index(hash)%number_of_hash_runs = ii
                        hierarchy(level)%hashed_coefficient_list_index(hash)%array_position = counter
                        hierarchy(level)%coefficient_list_index(counter) = uindex_list(i)%n
                        hierarchy(level)%hierarchy_type(counter) = hierarchy_type

                        if (hierarchy(level)%maximum_number_of_hash_runs .lt. ii) then
                            hierarchy(level)%maximum_number_of_hash_runs = ii
                        end if

                        exit
                    else if (hierarchy(level)%coefficient_list_index(counter+1) .ne. uindex_list(i)%n) then
                        ! hash collision -> re-hash
                        hash = IEOR(hash, int(ii, UINDEX_BYTES))
                        hash = hash_fnv1a(hash, max_value)
                    else
                        ! uindex exists already in the hierarchy
                        if (hierarchy(level)%hierarchy_type(counter+1) .ne. hierarchy_type) then
                            ! hierarchy type differ -> merge hierarchy types
                            hierarchy(level)%hierarchy_type(counter+1) = HIERARCHY_XY
                        end if
#ifdef _DEBUG_
                        print *, "lib_ml_fmm_hf_add_uindex_to_hierarchy_hashed: NOTE"
                        print *, "  uindex bypassed at level: ", level
                        print *, "  hierarchy_type", hierarchy_type
                        print *, "  uindex_list(i)%n: ", uindex_list(i)%n
#endif
                    end if
                end do
                if (ii .gt. MAXIMUM_NUMBER_OF_HASH_RUNS) then
                    print *, "lib_ml_fmm_hf_add_uindex_to_hierarchy_hashed: ERROR"
                    print *, "  uindex bypassed at level: ", level
                    print *, "  hierarchy_type", hierarchy_type
                    print *, "  uindex_list(i)%n: ", uindex_list(i)%n
                end if
            end do

        end subroutine lib_ml_fmm_hf_add_uindex_to_hierarchy_hashed

        ! Add the uindex to the hierarchy by using a lookup table
        !
        ! Arguments
        ! ----
        !   hierarchy: array<lib_ml_fmm_hierarchy>, inout
        !       reference to the hierarchy object to be processed
        !   level: integer, in
        !       level of the hierarchy to be processed
        !   uindex_list: array<lib_tree_universal_index>, inout (read-only)
        !       list of uindex' to add to the hierarchy
        !   hierarchy_type: integer, in
        !       type of the uindex_list [HIERARCHY_X, HIERARCHY_Y, HIERARCHY_XY]
        !   counter: integer, inout
        !       counter of added uindex' to the hierarchy()
        !
        subroutine lib_ml_fmm_hf_add_uindex_to_hierarchy_lookup_table(hierarchy, level, uindex_list, hierarchy_type, counter)
            implicit none
            ! dummy
            type(lib_ml_fmm_hierarchy), dimension(:), allocatable, intent(inout) :: hierarchy
            integer(kind=1), intent(in) :: level
            type(lib_tree_universal_index), dimension(:), intent(inout) :: uindex_list
            integer(kind=1), intent(in) :: hierarchy_type
            integer(kind=UINDEX_BYTES) :: counter

            ! auxiliary
            integer(kind=UINDEX_BYTES) :: i
            integer(kind=UINDEX_BYTES) :: n

            do i=1, size(uindex_list)
                counter = counter + 1
                n = uindex_list(i)%n + 1
                if (hierarchy(level)%coefficient_list_index(n) .eq. IGNORE_ENTRY) then
                    hierarchy(level)%coefficient_list_index(n) = counter
                    hierarchy(level)%hierarchy_type(counter) = hierarchy_type
                else
                    if (hierarchy(level)%hierarchy_type(counter) .ne. hierarchy_type) then
                        hierarchy(level)%hierarchy_type(counter) = HIERARCHY_XY
                    end if
#ifdef _DEBUG_
                    print *, "lib_ml_fmm_hf_add_uindex_to_hierarchy: NOTE"
                    print *, "  uindex bypassed at level: ", level
                    print *, "  hierarchy_type", hierarchy_type
                    print *, "  uindex_list(i)%n: ", n
#endif
                end if
            end do

        end subroutine lib_ml_fmm_hf_add_uindex_to_hierarchy_lookup_table



        function lib_ml_fmm_hf_get_index(hierarchy, uindex) result(rv)
            implicit none
            ! dummy
            type(lib_ml_fmm_hierarchy), dimension(:), allocatable, intent(inout) :: hierarchy
            type(lib_tree_universal_index), intent(inout) :: uindex
            integer(UINDEX_BYTES) :: rv

            if (hierarchy(uindex%l)%is_hashed) then
                ! todo: hash
                continue
            else
                rv = hierarchy(uindex%l)%coefficient_list_index(uindex%n)
            end if

        end function lib_ml_fmm_hf_get_index

        ! Calculates all parent universal indices of a list of universal indices
        !
        ! Argument
        ! ----
        !   uindex_list: array<lib_tree_universal_index>, inout
        !       list of universal indices, inplace replacement with the parent universal indices
        !   step = 1: integer, optional
        !       selects the relative level: [[x-th great ,]grand-]parent
        !
        subroutine lib_ml_fmm_hf_get_parent_uindex_lists(uindex_list, step)
            implicit none
            type(lib_tree_universal_index), dimension(:), allocatable, intent(inout) :: uindex_list
            integer(kind=1), optional, intent(in) :: step

            ! auxiliary
            integer(kind=UINDEX_BYTES) :: i
            integer(kind=1) :: m_step

            if (allocated(uindex_list)) then
                m_step = 1
                if (present(step)) m_step = step

                ! get parent uindex lists
                do i=1, size(uindex_list)
                    uindex_list(i) = lib_tree_get_parent(uindex_list(i), m_step)
                end do
            end if
        end subroutine lib_ml_fmm_hf_get_parent_uindex_lists

        ! Returns a list of unique univeral indices
        !
        ! Argument
        ! ----
        !   uindex_list: array<lib_tree_universal_index>, inout (read-only)
        !       list of universal indices
        !
        ! Returns
        ! ----
        !   rv: array<lib_tree_universal_index>, inout (read-only)
        !       list of unique universal indices
        !
        function lib_ml_fmm_hf_make_uindex_list_unique(uindex_list) result(rv)
            implicit none
            ! dummy
            type(lib_tree_universal_index), dimension(:), allocatable, intent(inout) :: uindex_list
            type(lib_tree_universal_index), dimension(:), allocatable :: rv

            ! auxiliary
            type(hash_list_type), dimension(size(uindex_list) * LIB_ML_FMM_HF_HIERARCHY_MARGIN) :: hash_list
            type(lib_tree_universal_index) :: buffer_uindex
            integer(kind=4) :: max_value
            integer(kind=4) :: hash
            integer(kind=UINDEX_BYTES) :: counter
            integer(kind=UINDEX_BYTES) :: i
            integer(kind=2) :: ii

            if (size(uindex_list) .gt. 0) then

                hash_list(:)%value = -1
                hash_list(:)%hash_runs = 0

                max_value = size(hash_list)

                counter = 0
                do i=1, size(uindex_list)
                    ! hash runs
                    do ii=1, MAXIMUM_NUMBER_OF_HASH_RUNS
                        hash = hash_fnv1a(uindex_list(i)%n, max_value)
                        if (hash_list(hash)%hash_runs .eq. 0) then
                            hash_list(hash)%value = uindex_list(i)%n
                            hash_list(hash)%hash_runs = ii
                            counter = counter + 1
                            exit
                        else if (hash_list(hash)%value .ne. uindex_list(i)%n) then
                            hash = IEOR(hash, int(ii, 4))
                            hash = hash_fnv1a(hash, max_value)
                        end if
                    end do
                end do

                allocate (rv(counter))

                counter = 0
                buffer_uindex%l = uindex_list(1)%l
                do i=1, size(hash_list)
                    if (hash_list(i)%hash_runs .ne. 0) then
                        counter = counter + 1
                        buffer_uindex%n = hash_list(i)%value
                        rv(counter) = buffer_uindex
                    end if
                end do

            else
                rv = uindex_list
            end if
        end function

        ! --- test functions ---
        function lib_ml_fmm_hf_test_functions() result(error_counter)
            implicit none
            ! dummy
            integer :: error_counter

            error_counter = 0

            if (.not. test_lib_ml_fmm_hf_get_parent_uindex_lists()) then
                error_counter = error_counter + 1
            end if
            if (.not. test_lib_ml_fmm_hf_make_uindex_list_unique()) then
                error_counter = error_counter + 1
            end if

            contains

                function test_lib_ml_fmm_hf_get_parent_uindex_lists() result(rv)
                    implicit none
                    ! dummy
                    logical :: rv

                    ! auxiliary
                    type(lib_tree_universal_index), dimension(:), allocatable :: uindex_list
                    type(lib_tree_universal_index), dimension(:), allocatable :: ground_truth_uindex_parent
                    type(lib_tree_universal_index), dimension(:), allocatable :: uindex_parent

                    integer(kind=1) :: level
                    integer(kind=1) :: step
                    integer(kind=1) :: ground_truth_level
                    integer :: length
                    integer :: ground_truth_length
                    integer :: i

                    step = 1
                    level = 2
                    ground_truth_level = level - step

                    length = 2**(TREE_DIMENSIONS*level)
                    ground_truth_length = length

                    allocate (uindex_list(length))
                    allocate (ground_truth_uindex_parent(ground_truth_length))

                    do i=1, length
                        uindex_list(i)%n = i-1
                        uindex_list(i)%l = level
                    end do

                    do i=1, ground_truth_length
                        ground_truth_uindex_parent(i) = lib_tree_get_parent(uindex_list(i), step)
                    end do

                    allocate(uindex_parent, source=uindex_list)
                    call lib_ml_fmm_hf_get_parent_uindex_lists(uindex_parent, step)

                    rv = .true.
                    do i=1, ground_truth_length
                        if (ground_truth_uindex_parent(i) .ne. uindex_parent(i)) then
                            rv = .false.
                        end if
                    end do

                    if (rv) then
                        print *, "test_lib_ml_fmm_hf_get_parent_uindex_lists: OK"
                    else
                        print *, "test_lib_ml_fmm_hf_get_parent_uindex_lists: FAILED"
                    end if

                    rv = .false.
                end function test_lib_ml_fmm_hf_get_parent_uindex_lists

                function test_lib_ml_fmm_hf_make_uindex_list_unique() result(rv)
                    implicit none
                    ! dummy
                    logical :: rv

                    ! auxiliary
                    type(lib_tree_universal_index), dimension(:), allocatable :: uindex_list
                    type(lib_tree_universal_index), dimension(:), allocatable :: ground_truth_uindex_reduced
                    type(lib_tree_universal_index), dimension(:), allocatable :: uindex_reduced

                    integer(kind=UINDEX_BYTES) :: uindex_sum
                    integer(kind=UINDEX_BYTES) :: ground_truth_uindex_sum
                    integer :: length
                    integer :: ground_truth_length
                    integer :: i

                    length = 10
                    ground_truth_length = 4

                    allocate (uindex_list(length))
                    allocate (ground_truth_uindex_reduced(ground_truth_length))

                    do i=1, length
                        uindex_list(i)%n = mod(i, ground_truth_length)
                    end do

                    do i=1, ground_truth_length
                        ground_truth_uindex_reduced(i)%n = i-1
                    end do

                    uindex_reduced = lib_ml_fmm_hf_make_uindex_list_unique(uindex_list)

                    if (size(uindex_reduced) .eq. size(ground_truth_uindex_reduced)) then
                        uindex_sum = 0
                        ground_truth_uindex_sum = 0
                        do i=1, ground_truth_length
                            uindex_sum = uindex_sum + uindex_reduced(i)%n
                            ground_truth_uindex_sum = ground_truth_uindex_sum &
                                                      + ground_truth_uindex_reduced(i)%n
                        end do

                        if (uindex_sum .eq. ground_truth_uindex_sum) then
                            print *, "test_lib_ml_fmm_hf_make_uindex_list_unique: OK"
                        else
                            print *, "test_lib_ml_fmm_hf_make_uindex_list_unique: FAILED"
                        end if

                    else
                        rv = .false.
                        print *, "test_lib_ml_fmm_hf_make_uindex_list_unique: FAILED"
                    end if

                end function test_lib_ml_fmm_hf_make_uindex_list_unique

        end function lib_ml_fmm_hf_test_functions

end module lib_ml_fmm_helper_functions
