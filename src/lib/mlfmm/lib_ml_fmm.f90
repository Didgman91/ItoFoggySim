! LIB: Mulitlevel Fast Multipole Method
!
module lib_ml_fmm
    use lib_tree
    use lib_tree_type
    use ml_fmm_type
    use ml_fmm_math
    use lib_ml_fmm_type_operator
    use lib_ml_fmm_helper_functions
    implicit none

    private

    ! --- public functions ---
    public :: lib_ml_fmm_test_functions

    interface
!        ! Basis function: A
!        function lib_ml_fmm_get_A(x) result(A)
!            use lib_tree_type
!            use ml_fmm_type
!            implicit none
!            ! dummy
!            type(lib_tree_spatial_point), intent(in) :: x
!            type(lib_ml_fmm_A) :: A
!        end function

        ! Basis function: A_i
        function lib_ml_fmm_get_A_i(x, data_element) result(A_i)
            use lib_tree_type
            use ml_fmm_type
            implicit none
            ! dummy
            type(lib_tree_spatial_point), intent(in) :: x
            type(lib_tree_data_element) :: data_element
            type(lib_ml_fmm_coefficient) :: A_i
        end function

!        ! Basis function: B
!        function lib_ml_fmm_get_B(x) result(B)
!            use lib_tree_type
!            use ml_fmm_type
!            implicit none
!            ! dummy
!            type(lib_tree_spatial_point), intent(in) :: x
!            type(lib_ml_fmm_B) :: B
!        end function

        ! Basis function: B_i
        function lib_ml_fmm_get_B_i(x, data_element) result(B_i)
            use lib_tree_type
            use ml_fmm_type
            implicit none
            ! dummy
            type(lib_tree_spatial_point), intent(in) :: x
            type(lib_tree_data_element) :: data_element
            type(lib_ml_fmm_coefficient) :: B_i
        end function

        ! Basis function: S
        function lib_ml_fmm_get_S(x) result(S)
            use lib_tree_type
            use ml_fmm_type
            implicit none
            ! dummy
            type(lib_tree_spatial_point), intent(in) :: x
            type(lib_ml_fmm_S) :: S
        end function

        ! Basis function: R
        function lib_ml_fmm_get_R(x) result(R)
            use lib_tree_type
            use ml_fmm_type
            implicit none
            ! dummy
            type(lib_tree_spatial_point), intent(in) :: x
            type(lib_ml_fmm_R) :: R
        end function

        ! Local expansion (inner or Regular expansion)
        !
        ! Restriction
        ! ----
        ! Example calculation:
        !       phi_i(y) = A_i (x_∗) o R(y − x_∗)     (8)
        !
        !   "Here the series is valid in the domain |y − x_∗| .le. r_c |x_i − x_∗| (see Fig. 2),
        !   where 0 < r c < 1 is some real number."
        !
        !   Reference: Data_Structures_Optimal_Choice_of_Parameters_and_C
        !
        ! Arguments
        ! ----
        !   i
        !   x_*
        !   y
        !
        ! Returns
        ! ----
        !   phi_i
        !
        !
        function lib_ml_fmm_expansion_R(i, x, y) result(phi_i)
            use lib_tree_type
            implicit none
            ! dummy
            integer, intent(in) :: i
            type(lib_tree_spatial_point), intent(in) :: x
            type(lib_tree_spatial_point), intent(in) :: y
            integer :: phi_i                   ! todo: define type
        end function lib_ml_fmm_expansion_R

        ! Far field expansion (outer, Singular or multipole expansion)
        !
        ! Restriction
        ! ----
        ! Example calculation:
        !   "Any function phi_i(y) has a complementary expansion valid outside
        !   a d-dimensional sphere centered at y = x_∗ with radius R_c |x_i − x_∗| :
        !       phi_i(y) = B_i (x_∗) o S(y − x_∗), |y - x_*| .ge. R_c |x_i - x_*|,     (10)
        !
        !   where R_c > 1 is a real number similar to r_c".
        !   "Even though for many physical fields, such as the Green’s function
        !   for Laplace’s equation, the function S(y − x_∗) is singular at y = x_∗ ,
        !   this condition is not necessary. In particular we can have S = R."
        !
        !   Reference: Data_Structures_Optimal_Choice_of_Parameters_and_C
        !
        ! Arguments
        ! ----
        !   i
        !   x_*
        !   y
        !
        ! Returns
        ! ----
        !   phi_i
        !
        !
        function lib_ml_fmm_expansion_S(i, x, y) result(phi_i)
            use lib_tree_type
            implicit none
            ! dummy
            integer, intent(in) :: i
            type(lib_tree_spatial_point), intent(in) :: x
            type(lib_tree_spatial_point), intent(in) :: y
            integer :: phi_i                   ! todo: define type
        end function lib_ml_fmm_expansion_S

        ! Translation: local-to-local (Regular-to-Regular)
        !
        ! Arguments
        ! ----
        !   A_i_1
        !   x_1
        !   x_2
        !
        ! Returns
        ! ----
        !   A_i_2
        !
        function lib_ml_fmm_translation_RR(A_i_1, x_1, x_2) result(A_i_2)
            use lib_tree_type
            use ml_fmm_type
            implicit none
            ! dummy
            type(lib_ml_fmm_coefficient), intent(in) :: A_i_1
            type(lib_tree_spatial_point), intent(in) :: x_1
            type(lib_tree_spatial_point), intent(in) :: x_2
            type(lib_ml_fmm_coefficient) :: A_i_2
        end function

        ! Translation: far-to-local (Singular-to-Regular)
        !
        ! Arguments
        ! ----
        !   B_i_1
        !       set of expansion coefficients
        !   x_1: spatial point
        !       origin of coordinate system 1
        !   x_1: spatial point
        !       origin of coordinate system 2
        !
        ! Returns
        ! ----
        !   A_i_2
        !       set of expansion coefficients
        !
        function lib_ml_fmm_translation_SR(B_i_1, x_1, x_2) result(A_i_2)
            use lib_tree_type
            use ml_fmm_type
            implicit none
            ! dummy
            type(lib_ml_fmm_coefficient), intent(in) :: B_i_1
            type(lib_tree_spatial_point), intent(in) :: x_1
            type(lib_tree_spatial_point), intent(in) :: x_2
            type(lib_ml_fmm_coefficient) :: A_i_2
        end function

        ! Translation: far-to-far (Singular-to-Singular)
        !
        ! Arguments
        ! ----
        !   B_i_1
        !       set of expansion coefficients
        !   x_1: spatial point
        !       origin of coordinate system 1
        !   x_2: spatial point
        !       origin of coordinate system 2
        !
        ! Returns
        ! ----
        !   B_i_2
        !       set of expansion coefficients
        !
        function lib_ml_fmm_translation_SS(B_i_1, x_1, x_2) result(B_i_2)
            use lib_tree_type
            use ml_fmm_type
            implicit none
            ! dummy
            type(lib_ml_fmm_coefficient), intent(in) :: B_i_1
            type(lib_tree_spatial_point), intent(in) :: x_1
            type(lib_tree_spatial_point), intent(in) :: x_2
            type(lib_ml_fmm_coefficient) :: B_i_2
        end function
    end interface


    type lib_ml_fmm_procedure_handles
        procedure(lib_ml_fmm_get_A_i), pointer, nopass :: get_A_i => null()
        procedure(lib_ml_fmm_get_B_i), pointer, nopass :: get_B_i => null()
        procedure(lib_ml_fmm_get_S), pointer, nopass :: get_S => null()
        procedure(lib_ml_fmm_get_R), pointer, nopass :: get_R => null()
        procedure(lib_ml_fmm_expansion_R), pointer, nopass :: expansion_R => null()
        procedure(lib_ml_fmm_expansion_S), pointer, nopass :: expansion_S => null()
        procedure(lib_ml_fmm_translation_RR), pointer, nopass :: get_translation_RR => null()
        procedure(lib_ml_fmm_translation_SR), pointer, nopass :: get_translation_SR => null()
        procedure(lib_ml_fmm_translation_SS), pointer, nopass :: get_translation_SS => null()
    end type lib_ml_fmm_procedure_handles

    type lib_ml_fmm_data
        type(lib_tree_data_element), dimension(:), allocatable :: X
        type(lib_tree_data_element), dimension(:), allocatable :: Y
        type(lib_tree_data_element), dimension(:), allocatable :: XY
    end type

    ! --- member ---
    type(lib_ml_fmm_procedure_handles) :: m_ml_fmm_handles
    integer(kind=2) :: m_ml_fmm_p_truncation

    ! e.g. matrix vector product v = u*phi
    type(lib_ml_fmm_v), dimension(:), allocatable :: m_ml_fmm_u
    type(lib_ml_fmm_v), dimension(:), allocatable :: m_ml_fmm_v

    type(lib_ml_fmm_coefficient_list_list) :: m_ml_fmm_expansion_coefficients

    type(lib_ml_fmm_hierarchy) :: m_ml_fmm_hierarchy

    ! Tree parameters
    integer(kind=UINDEX_BYTES) :: m_tree_neighbourhodd_size_k
    integer(kind=4) :: m_tree_s_opt
    integer(kind=1) :: m_tree_l_min
    integer(kind=1) :: m_tree_l_max


    contains

    ! Procedure
    ! 1. call lib_tree_constructor
    !       - determine s_opt
    !       - determine l_min and l_max
    ! 2. create ml fmm data set
    !       - C per box and per level (l_max up to l_min)
    !       - D per box and per level (l_min down to l_max)
    subroutine lib_ml_fmm_constructor(data_elements)
        implicit none
        ! dummy
        type(lib_ml_fmm_data), intent(in) :: data_elements

        ! auxiliaray
        type(lib_tree_data_element), dimension(:), allocatable :: data_concatenated
        integer(kind=UINDEX_BYTES), dimension(3) :: length
        type(ml_fmm_type_operator_procedures) :: operator_procedures

        ! concatenate X-, Y- and XY-hierarchy data
        !   length(1) = size(data_elements%X)
        !   length(2) = size(data_elements%Y)
        !   length(3) = size(data_elements%XY)
        data_concatenated = lib_ml_fmm_concatenate_data_array(data_elements, length)

        ! initiate the Tree
        call lib_tree_constructor(data_concatenated)

        ! initiate the X- and Y-hierarchy
        m_ml_fmm_hierarchy = lib_ml_fmm_hf_create_hierarchy(data_concatenated, length, m_tree_l_max, m_tree_l_min)

        ! initiate the ml fmm type operators
        operator_procedures = ml_fmm_type_operator_get_procedures()
        call lib_ml_fmm_type_operator_constructor(operator_procedures)

        ! adjust Tree parameters
        m_tree_neighbourhodd_size_k = 1!lib_ml_fmm_hf_get_neighbourhood_size(R_c1, r_c2)
        m_tree_l_min = lib_tree_get_level_min(m_tree_neighbourhodd_size_k)
        m_tree_l_max = lib_tree_get_level_max(m_tree_s_opt)

        ! setup the ml fmm data structure
        call lib_ml_fmm_type_operator_allocate_coefficient_list(m_ml_fmm_expansion_coefficients, m_tree_l_min, m_tree_l_max)


    end subroutine lib_ml_fmm_constructor

    ! clean up
    ! - coefficents
    ! - Tree
    subroutine lib_ml_fmm_destructor()
        implicit none

        ! auxilary
        ! example for an indiviual deallocation
!        integer(kind=1) :: i
!        integer(kind=4) :: ii

        call lib_ml_fmm_type_operator_deallocate_coefficient_list(m_ml_fmm_expansion_coefficients)

!        if ( allocated(m_ml_fmm_C) ) then
!            ! HINT: m_ml_fmm_C has a deep structure. If the automatic deallocation fails,
!            !       each subelement should be deallocated individually.
!            deallocate (m_ml_fmm_C)
!            ! example for an indiviual deallocation
!!            do i=1, size(m_ml_fmm_C)
!!                do ii=1, size(m_ml_fmm_C(i)%dummy)
!!                    if (allocated (m_ml_fmm_C(i)%dummy))
!!                        deallocate (m_ml_fmm_C(i)%dummy)
!!                    end if
!!                end do
!!            end do
!        end if
!
!        if ( allocated(m_ml_fmm_D) ) then
!            ! HINT: m_ml_fmm_D has a deep structure. If the automatic deallocation fails,
!            !       each subelement should be deallocated individually.
!            !       (as like m_ml_fmm_C)
!            deallocate(m_ml_fmm_D)
!        end if

        call lib_tree_destructor()

    end subroutine lib_ml_fmm_destructor

    function lib_ml_fmm_concatenate_data_array(data_elements, length) result(data_concatenated)
        implicit none
        type(lib_ml_fmm_data), intent(in) :: data_elements
        integer(kind=UINDEX_BYTES), dimension(3), intent(out) :: length
        type(lib_tree_data_element), dimension(:), allocatable :: data_concatenated

        ! auxiliaray
        integer(kind=1) :: i
        integer(kind=UINDEX_BYTES) :: start, last

        length(:) = 0
        if( allocated(data_elements%X) ) then
            length(HIERARCHY_X) = size(data_elements%X)
        end if

        if( allocated(data_elements%Y) ) then
            length(HIERARCHY_Y) = size(data_elements%Y)
        end if

        if( allocated(data_elements%XY) ) then
            length(HIERARCHY_XY) = size(data_elements%XY)
        end if

        ! concatenate all three datasets
        allocate( data_concatenated(sum(length)))

        do i=1, 3
            if (length(i) .gt. 0) then
                start = sum(length(:i))
                last = start + length(i)
                if (i .eq. HIERARCHY_X) then
                    data_concatenated(start:last) = data_elements%X
                    data_concatenated(start:last)%hierarchy = HIERARCHY_X
                else if (i .eq. HIERARCHY_Y) then
                    data_concatenated(start:last) = data_elements%Y
                    data_concatenated(start:last)%hierarchy = HIERARCHY_Y
                else if (i .eq. HIERARCHY_XY) then
                    data_concatenated(start:last) = data_elements%XY
                    data_concatenated(start:last)%hierarchy = HIERARCHY_XY
                end if
            end if
        end do
    end function



    subroutine lib_ml_fmm_calculate_upward_pass()
        implicit none

        call lib_ml_fmm_calculate_upward_pass_step_1()
        call lib_ml_fmm_calculate_upward_pass_step_2()
    end subroutine

    ! Upward pass - step 1
    !
    ! Procedure
    ! ----
    !   for each box at l_max
    !       1. get all elements
    !       2. calc B coefficients
    !       3. multiply with u_i
    !       4. calculate the sum of all elements of a box
    !
    ! Reference: Data_Structures_Optimal_Choice_of_Parameters_and_C, eq.(32)
    !
    subroutine lib_ml_fmm_calculate_upward_pass_step_1()
        use lib_tree
        implicit none
        ! dummy
        type(lib_ml_fmm_coefficient) :: C

        ! auxiliaray
        integer(kind=UINDEX_BYTES) :: number_of_boxes
        type(lib_tree_universal_index) :: uindex
        logical :: ignore_box

        integer(kind=UINDEX_BYTES) :: i

        number_of_boxes = lib_tree_get_number_of_boxes(m_tree_l_max)

        uindex%l = m_tree_l_max
        do i=1, number_of_boxes
            uindex%n = i

!            C = lib_ml_fmm_get_C_i_from_elements_at_box(uindex, ignore_box)
!            if (.not. ignore_box) then
!                call lib_ml_fmm_type_operator_set_coefficient(C, uindex, m_ml_fmm_expansion_coefficient)
!            end if
        end do

    end subroutine lib_ml_fmm_calculate_upward_pass_step_1

    !
    ! Arguments
    ! ----
    !   uindex: type(lib_tree_universal_index)
    !       universal index of the a Tree-box
    !
    ! Returns
    ! ----
    !   C_i: type(lib_ml_fmm_coefficient)
    !       C coefficient of the Tree-box
    function lib_ml_fmm_get_C_i_from_elements_at_box(uindex, ignore_box) result(C_i)
        implicit none
        ! dummy
        type(lib_tree_universal_index), intent(in) :: uindex
        logical, intent(out) :: ignore_box
        type(lib_ml_fmm_coefficient) :: C_i

        ! auxiliary
        type(lib_ml_fmm_coefficient) :: buffer_C_i

        type(lib_tree_data_element), dimension(:), allocatable :: data_element
        integer(kind=4), dimension(:), allocatable :: element_number

        type(lib_tree_spatial_point) :: x_c

        integer(kind=UINDEX_BYTES) :: i

        call lib_ml_fmm_type_operator_set_coefficient_zero(C_i)

        data_element = lib_tree_get_domain_e1(uindex, element_number)
        if ((allocated (data_element)) &
            .and. (size(data_element) .gt. 0)) then
            ignore_box = .false.
            x_c = lib_tree_get_centre_of_box(uindex)
            do i=1, size(data_element)
                buffer_C_i = m_ml_fmm_u(element_number(i)) * m_ml_fmm_handles%get_B_i(x_c, data_element(i))
                C_i = C_i + buffer_C_i
            end do
        else
            ignore_box = .true.
        end if
    end function lib_ml_fmm_get_C_i_from_elements_at_box

    ! Upward pass - step 2
    !
    ! Procedure
    ! ----
    !   for each box at l_max -1, ..., l_min
    !       1. "reexpanding v^(1)_Children(X;n,l),l+1 (y) near the center of box (n, l) and
    !          summing up the contribution of all the child boxes"
    !
    ! Reference: Data_Structures_Optimal_Choice_of_Parameters_and_C, eq.(33)
    !
    subroutine lib_ml_fmm_calculate_upward_pass_step_2()
        implicit none
        ! auxiliary
        type(lib_ml_fmm_coefficient) :: C
        type(lib_tree_universal_index) :: uindex

        integer(kind=UINDEX_BYTES) :: number_of_boxes
        integer(kind=1) :: l

        integer(kind=UINDEX_BYTES) :: n


        do l=m_tree_l_max-int(1, 1), m_tree_l_min, -int(1, 1)
            uindex%l = l
            number_of_boxes = lib_tree_get_number_of_boxes(l)
            do n=1, number_of_boxes
                uindex%n = n - int(1, UINDEX_BYTES)
                C = lib_ml_fmm_get_C_of_box(uindex)
!                call lib_ml_fmm_type_operator_set_coefficient(C, uindex, m_ml_fmm_C)
            end do
        end do

    end subroutine lib_ml_fmm_calculate_upward_pass_step_2

    ! Argument
    ! ----
    !   uindex: type(lib_tree_universal_index)
    !       universal index of a box
    !
    ! Returns
    ! ----
    !   C: type(lib_ml_fmm_coefficient)
    !       expansion coefficients of the box uindex
    function lib_ml_fmm_get_C_of_box(uindex) result(C)
        implicit none
        ! dummy
        type(lib_tree_universal_index), intent (in) :: uindex
        type(lib_ml_fmm_coefficient) :: C

        ! auxiliary
        type(lib_ml_fmm_coefficient) :: C_child
        type(lib_tree_universal_index), dimension(:), allocatable :: uindex_children
        type(lib_tree_spatial_point) :: x_c
        type(lib_tree_spatial_point) :: x_c_child

        integer(kind=UINDEX_BYTES) :: i

        x_c = lib_tree_get_centre_of_box(uindex)

        ! get child boxes
        uindex_children = lib_tree_get_children(uindex)

        call lib_ml_fmm_type_operator_set_coefficient_zero(C)
        do i=1, size(uindex_children)
            x_c_child = lib_tree_get_centre_of_box(uindex_children(i))
!            C_child = lib_ml_fmm_type_operator_get_coefficient(uindex_children(i), m_ml_fmm_C)
            C = C + m_ml_fmm_handles%get_translation_SS(C_child, x_c_child, x_c)
        end do
    end function lib_ml_fmm_get_C_of_box

    ! Downward pass - step 1
    !
    ! Procedure
    ! ----
    !   for each box at l_min, ..., l_max
    !       - SR-transformation of the C coefficients into the middle of box I_4
    !
    ! Reference: Data_Structures_Optimal_Choice_of_Parameters_and_C, eq.(34)
    !
    function lib_ml_fmm_calculate_downward_pass_step_1() result(dummy)
        implicit none
        integer :: dummy


    end function lib_ml_fmm_calculate_downward_pass_step_1

    ! Downward pass - step 2
    !
    ! Procedure
    ! ----
    !   for boxes at l_min, ..., l_max
    !       1.
    !
    ! Reference: Data_Structures_Optimal_Choice_of_Parameters_and_C, eq.(36)
    !
    function lib_ml_fmm_calculate_downward_pass_step_2() result(dummy)
        implicit none
        integer :: dummy
    end function lib_ml_fmm_calculate_downward_pass_step_2

    ! ----- test functions -----
    function lib_ml_fmm_test_functions() result(error_counter)
        implicit none

        integer :: error_counter

        if (.not. test_calculate_upward_pass_step_1()) then
            error_counter = error_counter + 1
        end if

        print *, "-------------lib_ml_fmm_test_functions----------------"
        if (error_counter == 0) then
            print *, "lib_ml_fmm_test_functions tests: OK"
        else
            print *, error_counter,"lib_ml_fmm_test_functions test(s) FAILED"
        end if
        print *, "------------------------------------------------------"

        contains

        function test_lib_ml_fmm_constructor() result(rv)
            implicit none
            ! dummy
            logical :: rv

            ! auxiliary
            integer(kind=UINDEX_BYTES), parameter :: list_length = 10
            integer(kind=1), parameter :: l_th = 5
            integer(kind=1), parameter :: element_type = 1
            type(lib_tree_data_element), dimension(list_length) :: element_list

            element_list = lib_tree_get_diagonal_test_dataset(list_length, l_th, element_type, HIERARCHY_X)

            rv = .false.

        end function

        function test_calculate_upward_pass_step_1() result(rv)
            implicit none
            ! dummy
            logical :: rv

        end function test_calculate_upward_pass_step_1

    end function lib_ml_fmm_test_functions

end module lib_ml_fmm
