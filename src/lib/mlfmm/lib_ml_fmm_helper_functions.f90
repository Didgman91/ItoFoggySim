! LIB: Mulitlevel Fast Multipole Method - Helper Functions
!
module lib_ml_fmm_helper_functions
    use lib_tree_public
    use ml_fmm_type
    implicit none

    private

    ! --- public functions ---
    public :: lib_ml_fmm_hf_get_neighbourhood_size
    public :: lib_ml_fmm_hf_create_hierarchy

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
        function lib_ml_fmm_hf_create_hierarchy(data_elements, length, l_max, l_min) result(hierarchy)
            implicit none
            ! dummy
            type(lib_tree_data_element), dimension(:), intent(in) :: data_elements
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
            integer(kind=1) :: hierarchy_type
            integer(kind=UINDEX_BYTES), dimension(3) :: uindex_list_counter
            integer(kind=UINDEX_BYTES) :: number_of_boxes_at_l_max

            type(lib_tree_correspondece_vector_element), dimension(:), allocatable :: correspondence_vector
            integer(kind=CORRESPONDENCE_VECTOR_KIND) :: element_index

            allocate( hierarchy(l_min:l_max))

            allocate( uindex_list_X(length(HIERARCHY_X)) )
            allocate( uindex_list_Y(length(HIERARCHY_Y)) )
            allocate( uindex_list_XY(length(HIERARCHY_XY)) )
            uindex_list_counter(:) = 0

            ! iterate over all correspondence vector elements (l = l_th)
            correspondence_vector = lib_tree_get_correspondence_vector()
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

!            ! iterate over all data elements of the X-hierarchy
!            if( length(HIERARCHY_X) .gt. 0 ) then
!                allocate( uindex_list(length(HIERARCHY_X)) )
!                ! setup at level l_max
!                do i=1, length(HIERARCHY_X)
!                    ! get box uindex at l=l_max
!                    buffer_uindex = data_elements(i)%uindex
!                    buffer_l = l_max - buffer_uindex%l
!                    if (buffer_l .gt. 0) then
!                        buffer_uindex = lib_tree_get_parent(buffer_uindex, buffer_l)
!                    else if (buffer_l .eq. 0) then
!                        continue
!                    else
!                        print *, "lib_ml_fmm_create_hierarchy: ERROR"
!                        print *, "  l_th < l_max"
!                        print *, "  check for example the Tree constructor (lib_tree_constructor)"
!                    end if
!
!                    uindex_list(i) = buffer_uindex
!                end do
!                hierarchy(l_max)%type = HIERARCHY_X
!
!
!            end if
        end function

        subroutine lib_ml_fmm_hf_add_uindex_to_hierarchy(hierarchy, uindex)
            implicit none
            ! dummy
            type(lib_ml_fmm_hierarchy), dimension(:), intent(inout) :: hierarchy
            type(lib_tree_universal_index), intent(in) :: uindex

!            hierarchy()
        end subroutine

end module lib_ml_fmm_helper_functions
