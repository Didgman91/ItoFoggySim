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

    public :: lib_ml_fmm_hf_test_functions

    ! --- member ---
    type(lib_ml_fmm_procedure_handles) :: m_ml_fmm_handles
    integer(kind=2) :: m_ml_fmm_p_truncation

    ! e.g. matrix vector product v = u*phi
    type(lib_ml_fmm_v), dimension(:), allocatable :: m_ml_fmm_u
    type(lib_ml_fmm_v), dimension(:), allocatable :: m_ml_fmm_v

    type(lib_ml_fmm_hierarchy), dimension(:), allocatable :: m_ml_fmm_hierarchy

    ! Tree parameters
    integer(kind=UINDEX_BYTES) :: m_tree_neighbourhood_size_k
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
        type(lib_ml_fmm_data), intent(inout) :: data_elements

        ! auxiliaray
        type(lib_tree_data_element), dimension(:), allocatable :: data_concatenated
        type(lib_tree_correspondece_vector_element), dimension(:), allocatable :: correspondence_vector
        integer(kind=UINDEX_BYTES), dimension(3) :: length
        type(ml_fmm_type_operator_procedures) :: operator_procedures

        ! concatenate X-, Y- and XY-hierarchy data
        !   length(1) = size(data_elements%X)
        !   length(2) = size(data_elements%Y)
        !   length(3) = size(data_elements%XY)
        allocate(data_concatenated, source=lib_ml_fmm_concatenate_data_array(data_elements, length))

        m_tree_s_opt = 1 ! todo: calculate s

        ! initiate the Tree
        call lib_tree_constructor(data_concatenated, m_tree_s_opt)
        data_concatenated = lib_tree_get_element_list()

        ! initiate the ml fmm type operators
        operator_procedures = ml_fmm_type_operator_get_procedures()
        call lib_ml_fmm_type_operator_constructor(operator_procedures)

        ! adjust Tree parameters
        m_tree_neighbourhood_size_k = 1!lib_ml_fmm_hf_get_neighbourhood_size(R_c1, r_c2)
        m_tree_l_min = lib_tree_get_level_min(m_tree_neighbourhood_size_k)
        m_tree_l_max = lib_tree_get_level_max(m_tree_s_opt)

        ! initiate the X- and Y-hierarchy
        correspondence_vector = lib_tree_get_correspondence_vector()
        call lib_ml_fmm_hf_create_hierarchy(data_concatenated, correspondence_vector, &
                                            length, m_tree_l_min, m_tree_l_max,&
                                            m_ml_fmm_hierarchy)

        allocate (m_ml_fmm_u(length(HIERARCHY_X) + length(HIERARCHY_XY)))
        allocate (m_ml_fmm_v(length(HIERARCHY_Y) + length(HIERARCHY_XY)))
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

        if (allocated(m_ml_fmm_hierarchy)) then
            deallocate(m_ml_fmm_hierarchy)
        end if

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

    function lib_ml_fmm_run(vector_u) result(vector_v)
        implicit none
        ! dummy
        type(lib_ml_fmm_v), dimension(:), allocatable, intent(inout) :: vector_u
        type(lib_ml_fmm_v), dimension(:), allocatable :: vector_v

        if (allocated(m_ml_fmm_u)) then
            m_ml_fmm_u = vector_u
        else
            allocate(m_ml_fmm_u, source=vector_u)
        end if

        call lib_ml_fmm_calculate_upward_pass()
        call lib_ml_fmm_calculate_downward_pass()
        call lib_ml_fmm_final_summation()

        allocate(vector_v, source=m_ml_fmm_v)

    end function

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
                if (i .eq. 1) then
                    start = 1
                else
                    start = sum(length(:i-1)) + 1
                end if
                last = start + length(i) - 1
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
    !       1. get all elements of the X-hierarchy
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
        integer(kind=UINDEX_BYTES) :: list_index
        integer(kind=1) :: hierarchy_type

        number_of_boxes = size(m_ml_fmm_hierarchy(m_tree_l_max)%coefficient_list_index)

        uindex%l = m_tree_l_max
        if (m_ml_fmm_hierarchy(m_tree_l_max)%is_hashed) then
            do i=1, number_of_boxes
                uindex%n = m_ml_fmm_hierarchy(uindex%l)%coefficient_list_index(i)
                hierarchy_type = m_ml_fmm_hierarchy(uindex%l)%hierarchy_type(i)
                if ((uindex%n .ge. 0) .and. &
                    ((hierarchy_type .eq. HIERARCHY_X) .or. &
                     (hierarchy_type .eq. HIERARCHY_XY))) then
                    C = lib_ml_fmm_get_C_i_from_elements_at_box(uindex, ignore_box)
                    if (.not. ignore_box) then
                        m_ml_fmm_hierarchy(m_tree_l_max)%coefficient_list(i) = C
                        m_ml_fmm_hierarchy(m_tree_l_max)%coefficient_type(i) = LIB_ML_FMM_COEFFICIENT_TYPE_C
                    end if
                end if
            end do
        else
            do i=1, number_of_boxes
                uindex%n = i - int(1, 1)
                list_index = m_ml_fmm_hierarchy(m_tree_l_max)%coefficient_list_index(i)
                hierarchy_type = m_ml_fmm_hierarchy(uindex%l)%hierarchy_type(i)
                if ((list_index .gt. 0) .and. &
                    ((hierarchy_type .eq. HIERARCHY_X) .or. &
                     (hierarchy_type .eq. HIERARCHY_XY))) then
                    C = lib_ml_fmm_get_C_i_from_elements_at_box(uindex, ignore_box)
                    if (.not. ignore_box) then
                        m_ml_fmm_hierarchy(m_tree_l_max)%coefficient_list(list_index) = C
                        m_ml_fmm_hierarchy(m_tree_l_max)%coefficient_type(list_index) = LIB_ML_FMM_COEFFICIENT_TYPE_C
                    end if
                end if
            end do
        end if

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

        integer(kind=UINDEX_BYTES) :: i
        integer(kind=UINDEX_BYTES) :: list_index
        integer(kind=1) :: hierarchy_type


        do l=m_tree_l_max-int(1, 1), m_tree_l_min, -int(1, 1)
            uindex%l = l
            number_of_boxes = size(m_ml_fmm_hierarchy(l)%coefficient_list_index)

            if (m_ml_fmm_hierarchy(l)%is_hashed) then
                do i=1, number_of_boxes
                    uindex%n = m_ml_fmm_hierarchy(l)%coefficient_list_index(i)
                    hierarchy_type = m_ml_fmm_hierarchy(l)%hierarchy_type(i)
                    if ((uindex%n .ge. 0) .and. &
                        ((hierarchy_type .eq. HIERARCHY_X) .or. &
                         (hierarchy_type .eq. HIERARCHY_XY))) then
                        C = lib_ml_fmm_get_C_of_box(uindex)
                        m_ml_fmm_hierarchy(l)%coefficient_list(i) = C
                        m_ml_fmm_hierarchy(l)%coefficient_type(i) = LIB_ML_FMM_COEFFICIENT_TYPE_C
                    end if
                end do
            else
                do i=1, number_of_boxes
                    uindex%n = i - int(1, 1)
                    list_index = m_ml_fmm_hierarchy(l)%coefficient_list_index(i)
                    hierarchy_type = m_ml_fmm_hierarchy(l)%hierarchy_type(i)
                    if ((list_index .gt. 0) .and. &
                        ((hierarchy_type .eq. HIERARCHY_X) .or. &
                         (hierarchy_type .eq. HIERARCHY_XY))) then
                        C = lib_ml_fmm_get_C_of_box(uindex)
                        m_ml_fmm_hierarchy(l)%coefficient_list(list_index) = C
                        m_ml_fmm_hierarchy(l)%coefficient_type(list_index) = LIB_ML_FMM_COEFFICIENT_TYPE_C
                    end if
                end do
            end if
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
        integer(kind=UINDEX_BYTES) :: list_index

        x_c = lib_tree_get_centre_of_box(uindex)

        ! get child boxes
        uindex_children = lib_tree_get_children(uindex)

        call lib_ml_fmm_type_operator_set_coefficient_zero(C)
        do i=1, size(uindex_children)
            x_c_child = lib_tree_get_centre_of_box(uindex_children(i))
            list_index = lib_ml_fmm_hf_get_hierarchy_index(m_ml_fmm_hierarchy, uindex_children(i))
            if (list_index .gt. 0) then
                C_child = m_ml_fmm_hierarchy(uindex%l)%coefficient_list(list_index)
                C = C + m_ml_fmm_handles%get_translation_SS(C_child, x_c_child, x_c)
            end if
        end do
    end function lib_ml_fmm_get_C_of_box

    subroutine lib_ml_fmm_calculate_downward_pass
        implicit none

        call lib_ml_fmm_calculate_downward_pass_step_1()
        call lib_ml_fmm_calculate_downward_pass_step_2()


    end subroutine

    ! Downward pass - step 1
    !
    ! Procedure
    ! ----
    !   for each box at l_min, ..., l_max
    !       - SR-transformation of the C coefficients into the middle of box I_4
    !
    ! Reference: Data_Structures_Optimal_Choice_of_Parameters_and_C, eq.(34)
    !
    subroutine lib_ml_fmm_calculate_downward_pass_step_1()
        implicit none
        ! auxilary
        type(lib_tree_universal_index) :: uindex
        type(lib_ml_fmm_coefficient) :: D_tilde

        integer(kind=1) :: l
        integer(kind=UINDEX_BYTES) :: i
        integer(kind=UINDEX_BYTES) :: number_of_boxes
        integer(kind=UINDEX_BYTES) :: list_index
        integer(kind=1) :: hierarchy_type

        call lib_ml_fmm_type_operator_set_coefficient_zero(D_tilde)
        do l=m_tree_l_min, m_tree_l_max
            uindex%l = l
            number_of_boxes = size(m_ml_fmm_hierarchy(l)%coefficient_list_index)

            if (m_ml_fmm_hierarchy(l)%is_hashed) then
                do i=1, number_of_boxes
                    uindex%n = m_ml_fmm_hierarchy(l)%coefficient_list_index(i)
                    hierarchy_type = m_ml_fmm_hierarchy(l)%hierarchy_type(i)
                    if ((uindex%n .ge. 0) .and. &
                        ((hierarchy_type .eq. HIERARCHY_Y) .or. &
                         (hierarchy_type .eq. HIERARCHY_XY))) then
                        D_tilde = lib_ml_fmm_get_D_tilde_of_box(uindex)
                        m_ml_fmm_hierarchy(l)%coefficient_list(i) = D_tilde
                        m_ml_fmm_hierarchy(l)%coefficient_type(i) = LIB_ML_FMM_COEFFICIENT_TYPE_D_TILDE
                    end if
                end do
            else
                do i=1, number_of_boxes
                    uindex%n = i - int(1, 1)
                    list_index = m_ml_fmm_hierarchy(l)%coefficient_list_index(i)
                    hierarchy_type = m_ml_fmm_hierarchy(l)%hierarchy_type(i)
                    if ((list_index .gt. 0) .and. &
                        ((hierarchy_type .eq. HIERARCHY_Y) .or. &
                         (hierarchy_type .eq. HIERARCHY_XY))) then
                        D_tilde = lib_ml_fmm_get_D_tilde_of_box(uindex)
                        m_ml_fmm_hierarchy(l)%coefficient_list(i) = D_tilde
                        m_ml_fmm_hierarchy(l)%coefficient_type(i) = LIB_ML_FMM_COEFFICIENT_TYPE_D_TILDE
                    end if
                end do
            end if
        end do


    end subroutine lib_ml_fmm_calculate_downward_pass_step_1

    function lib_ml_fmm_get_D_tilde_of_box(uindex) result(D_tilde)
        implicit none
        ! dummy
        type(lib_tree_universal_index), intent(inout) :: uindex
        type(lib_ml_fmm_coefficient) :: D_tilde

        ! auxiliary
        type(lib_ml_fmm_coefficient) :: C
        type(lib_tree_universal_index), dimension(:), allocatable :: boxes_i4

        type(lib_tree_spatial_point) :: x_c
        type(lib_tree_spatial_point) :: x_c_neighbour

        integer(kind=UINDEX_BYTES) :: i
        integer(kind=UINDEX_BYTES) :: list_index
        integer(kind=1) :: hierarchy_type

        allocate(boxes_i4, source=lib_tree_get_domain_i4(uindex, m_tree_neighbourhood_size_k))

        call lib_ml_fmm_type_operator_set_coefficient_zero(D_tilde)

        x_c = lib_tree_get_centre_of_box(uindex)

        do i=1, size(boxes_i4)
            list_index = lib_ml_fmm_hf_get_hierarchy_index(m_ml_fmm_hierarchy, uindex)
            hierarchy_type = m_ml_fmm_hierarchy(uindex%l)%hierarchy_type(list_index)
            if ((list_index .gt. 0) .and. &
                ((hierarchy_type .eq. HIERARCHY_X) .or. &
                 (hierarchy_type .eq. HIERARCHY_XY))) then
                C = m_ml_fmm_hierarchy(uindex%l)%coefficient_list(list_index)
                x_c_neighbour = lib_tree_get_centre_of_box(boxes_i4(i))

                D_tilde = D_tilde + m_ml_fmm_handles%get_translation_SR(C, x_c, x_c_neighbour)
            end if
        end do

    end function

    ! Downward pass - step 2
    !
    ! Procedure
    ! ----
    !   for boxes at l_min, ..., l_max
    !       1.
    !
    ! Reference: Data_Structures_Optimal_Choice_of_Parameters_and_C, eq.(36)
    !
    subroutine lib_ml_fmm_calculate_downward_pass_step_2()
        implicit none
        ! auxiliaray
        type(lib_ml_fmm_coefficient) :: D
        type(lib_tree_universal_index) :: uindex
        integer(kind=1) :: l
        integer(kind=UINDEX_BYTES) :: i
        integer(kind=UINDEX_BYTES) :: number_of_boxes
        integer(kind=UINDEX_BYTES) :: list_index
        integer(kind=1) :: hierarchy_type

        m_ml_fmm_hierarchy(m_tree_l_min)%coefficient_type(:) = LIB_ML_FMM_COEFFICIENT_TYPE_D

        do l=m_tree_l_min+int(1,1), m_tree_l_max
            uindex%l = l
            number_of_boxes = size(m_ml_fmm_hierarchy(l)%coefficient_list_index)

            if (m_ml_fmm_hierarchy(l)%is_hashed) then
                do i=1, number_of_boxes
                    uindex%n = m_ml_fmm_hierarchy(l)%coefficient_list_index(i)
                    hierarchy_type = m_ml_fmm_hierarchy(l)%hierarchy_type(i)
                    if ((uindex%n .ge. 0) .and. &
                        ((hierarchy_type .eq. HIERARCHY_Y) .or. &
                         (hierarchy_type .eq. HIERARCHY_XY))) then
                        D = lib_ml_fmm_get_D_of_box(uindex)
                        m_ml_fmm_hierarchy(l)%coefficient_list(i) = D
                        m_ml_fmm_hierarchy(l)%coefficient_type(i) = LIB_ML_FMM_COEFFICIENT_TYPE_D
                    end if
                end do
            else
                do i=1, number_of_boxes
                    uindex%n = i - int(1, 1)
                    list_index = m_ml_fmm_hierarchy(l)%coefficient_list_index(i)
                    hierarchy_type = m_ml_fmm_hierarchy(l)%hierarchy_type(i)
                    if ((list_index .gt. 0) .and. &
                        ((hierarchy_type .eq. HIERARCHY_Y) .or. &
                         (hierarchy_type .eq. HIERARCHY_XY))) then
                        D = lib_ml_fmm_get_D_of_box(uindex)
                        m_ml_fmm_hierarchy(l)%coefficient_list(i) = D
                        m_ml_fmm_hierarchy(l)%coefficient_type(i) = LIB_ML_FMM_COEFFICIENT_TYPE_D
                    end if
                end do
            end if
        end do

    end subroutine lib_ml_fmm_calculate_downward_pass_step_2

    function lib_ml_fmm_get_D_of_box(uindex) result(D)
        implicit none
        ! dummy
        type(lib_tree_universal_index), intent(inout) :: uindex
        type(lib_ml_fmm_coefficient) :: D

        ! auxiliary
        type(lib_ml_fmm_coefficient) :: D_parent
        type(lib_ml_fmm_coefficient) :: D_tilde
        type(lib_tree_universal_index) :: uindex_parent
        integer(kind=1) :: coefficient_type

        type(lib_tree_spatial_point) :: x_c
        type(lib_tree_spatial_point) :: x_c_parent

        uindex_parent = lib_tree_get_parent(uindex)

        if (uindex_parent%n .ge. 0) then
            x_c = lib_tree_get_centre_of_box(uindex)
            x_c_parent = lib_tree_get_centre_of_box(uindex_parent)
            D_tilde = lib_ml_fmm_hf_get_hierarchy_coefficient(m_ml_fmm_hierarchy, uindex, coefficient_type)
            D_parent = lib_ml_fmm_hf_get_hierarchy_coefficient(m_ml_fmm_hierarchy, uindex_parent, coefficient_type)
            D = D_tilde + m_ml_fmm_handles%get_translation_RR(D_parent, x_c_parent, x_c)
        else
            print *, "lib_ml_fmm_get_D_of_box: ERROR"
            print *, "  parent box not found"
        end if

    end function lib_ml_fmm_get_D_of_box

    subroutine lib_ml_fmm_final_summation()
        implicit none

        ! auxiliary
        type(lib_tree_universal_index) :: uindex
        integer(kind=UINDEX_BYTES) :: i
        integer(kind=UINDEX_BYTES) :: number_of_boxes
        integer(kind=UINDEX_BYTES) :: list_index
        integer(kind=1) :: hierarchy_type

        uindex%l = m_tree_l_max
        number_of_boxes = size(m_ml_fmm_hierarchy(uindex%l)%coefficient_list_index)

        if (m_ml_fmm_hierarchy(uindex%l)%is_hashed) then
            do i=1, number_of_boxes
                uindex%n = m_ml_fmm_hierarchy(uindex%l)%coefficient_list_index(i)
                hierarchy_type = m_ml_fmm_hierarchy(uindex%l)%hierarchy_type(i)
                if ((uindex%n .ge. 0) .and. &
                    ((hierarchy_type .eq. HIERARCHY_Y) .or. &
                     (hierarchy_type .eq. HIERARCHY_XY))) then

                    call lib_ml_fmm_calculate_all_v_y_j_at_uindex(uindex)
                end if
            end do
        else
            do i=1, number_of_boxes
                uindex%n = i - int(1, 1)
                list_index = m_ml_fmm_hierarchy(uindex%l)%coefficient_list_index(i)
                hierarchy_type = m_ml_fmm_hierarchy(uindex%l)%hierarchy_type(i)
                if ((list_index .gt. 0) .and. &
                    ((hierarchy_type .eq. HIERARCHY_Y) .or. &
                     (hierarchy_type .eq. HIERARCHY_XY))) then

                    call lib_ml_fmm_calculate_all_v_y_j_at_uindex(uindex)
                end if
            end do
        end if

    end subroutine lib_ml_fmm_final_summation

    subroutine lib_ml_fmm_calculate_all_v_y_j_at_uindex(uindex)
        implicit none
        ! dummy
        type(lib_tree_universal_index), intent(inout) :: uindex

        ! auxiliary
        type(lib_ml_fmm_coefficient) :: D
        integer(kind=1) :: coefficient_type
        integer(kind=UINDEX_BYTES) :: i
        type(lib_tree_data_element), dimension(:), allocatable :: data_element_e1
        integer(kind=CORRESPONDENCE_VECTOR_KIND), dimension(:), allocatable :: element_number_e1
        type(lib_tree_data_element), dimension(:), allocatable :: data_element_e2
        integer(kind=CORRESPONDENCE_VECTOR_KIND), dimension(:), allocatable :: element_number_e2


        allocate(data_element_e1, source = lib_tree_get_domain_e1(uindex, &
                                                                  element_number_e1))
        allocate(data_element_e2, source = lib_tree_get_domain_e2(m_tree_neighbourhood_size_k, &
                                                                  uindex, &
                                                                  element_number_e2))

        D = lib_ml_fmm_hf_get_hierarchy_coefficient(m_ml_fmm_hierarchy, uindex, coefficient_type)
        do i=1, size(data_element_e1)
            if ((data_element_e1(i)%hierarchy .eq. HIERARCHY_Y) .or. &
                (data_element_e1(i)%hierarchy .eq. HIERARCHY_XY)) then
                m_ml_fmm_v(element_number_e1(i)) = lib_ml_fmm_calculate_v_y_j(data_element_e1(i), data_element_e2, D)
            end if
        end do
    end subroutine lib_ml_fmm_calculate_all_v_y_j_at_uindex

    function lib_ml_fmm_calculate_v_y_j(data_element_y_j, data_element_e2, D) result(rv)
        implicit none
        ! dummy
        type(lib_tree_data_element), intent(inout) :: data_element_y_j
        type(lib_tree_data_element), dimension(:), allocatable, intent(inout) :: data_element_e2
        type(lib_ml_fmm_coefficient), intent(inout) :: D
        type(lib_ml_fmm_v) :: rv

        ! auxiliary
        integer(kind=UINDEX_BYTES) :: i
        type(lib_tree_spatial_point) :: x_c
        type(lib_tree_spatial_point) :: y_j

        x_c = lib_tree_get_centre_of_box(data_element_y_j%uindex)
        y_j = data_element_y_j%point_x

        rv = D .cor. (y_j - x_c)

        do i=1, size(data_element_e2)
            rv = rv + m_ml_fmm_handles%get_phi_i_j(data_element_e2(i), y_j)
        end do

    end function lib_ml_fmm_calculate_v_y_j

    ! ----- test functions -----
    function lib_ml_fmm_test_functions() result(error_counter)
        implicit none
        ! dummy
        integer :: error_counter

        error_counter = 0

        if (.not. test_lib_ml_fmm_constructor()) then
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
            type(lib_ml_fmm_v), dimension(:), allocatable :: vector_u
            type(lib_ml_fmm_v), dimension(:), allocatable :: vector_v

            integer(kind=UINDEX_BYTES), parameter :: list_length = 10
            integer(kind=1), parameter :: element_type = 1
            type(lib_ml_fmm_data) :: data_elements

            real(kind=LIB_ML_FMM_COEFFICIENT_KIND), dimension(:), allocatable :: dummy
            integer(kind=LIB_ML_FMM_COEFFICIENT_KIND) :: i

            allocate(data_elements%X, source=lib_tree_get_diagonal_test_dataset(list_length, element_type, HIERARCHY_X))
            allocate(data_elements%Y, source=lib_tree_get_diagonal_test_dataset(4_8, element_type, HIERARCHY_Y, .true.))

            call lib_ml_fmm_constructor(data_elements)

            allocate(vector_u(list_length+4_8))
            allocate(dummy(1))
            dummy(1) = 1.0
            do i=1, size(vector_u)
                vector_u(i)%dummy = dummy
            end do

            allocate(vector_v, source = lib_ml_fmm_run(vector_u))

            rv = .false.

        end function

        function test_calculate_upward_pass_step_1() result(rv)
            implicit none
            ! dummy
            logical :: rv

        end function test_calculate_upward_pass_step_1

    end function lib_ml_fmm_test_functions

end module lib_ml_fmm
