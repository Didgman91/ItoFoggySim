! LIB: Mulitlevel Fast Multipole Method
!
module lib_ml_fmm
    use lib_tree
    use lib_ml_fmm_type
    use lib_ml_fmm_helper_functions
    implicit none

    private

    interface
!        ! Basis function: A
!        function lib_ml_fmm_get_A(x) result(A)
!            use lib_tree_type
!            use lib_ml_fmm_type
!            implicit none
!            ! dummy
!            type(lib_tree_spatial_point), intent(in) :: x
!            type(lib_ml_fmm_A) :: A
!        end function

        ! Basis function: A_i
        function lib_ml_fmm_get_A_i(x, data_element) result(A_i)
            use lib_tree_type
            use lib_ml_fmm_type
            implicit none
            ! dummy
            type(lib_tree_spatial_point), intent(in) :: x
            type(lib_tree_data_element) :: data_element
            type(lib_ml_fmm_A_i) :: A_i
        end function

!        ! Basis function: B
!        function lib_ml_fmm_get_B(x) result(B)
!            use lib_tree_type
!            use lib_ml_fmm_type
!            implicit none
!            ! dummy
!            type(lib_tree_spatial_point), intent(in) :: x
!            type(lib_ml_fmm_B) :: B
!        end function

        ! Basis function: B_i
        function lib_ml_fmm_get_B_i(x, data_element) result(B_i)
            use lib_tree_type
            use lib_ml_fmm_type
            implicit none
            ! dummy
            type(lib_tree_spatial_point), intent(in) :: x
            type(lib_tree_data_element) :: data_element
            type(lib_ml_fmm_B_i) :: B_i
        end function

        ! Basis function: S
        function lib_ml_fmm_get_S(x) result(S)
            use lib_tree_type
            use lib_ml_fmm_type
            implicit none
            ! dummy
            type(lib_tree_spatial_point), intent(in) :: x
            type(lib_ml_fmm_S) :: S
        end function

        ! Basis function: R
        function lib_ml_fmm_get_R(x) result(R)
            use lib_tree_type
            use lib_ml_fmm_type
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
            use lib_ml_fmm_type
            implicit none
            ! dummy
            type(lib_ml_fmm_A_i), intent(in) :: A_i_1
            type(lib_tree_spatial_point), intent(in) :: x_1
            type(lib_tree_spatial_point), intent(in) :: x_2
            type(lib_ml_fmm_A_i) :: A_i_2
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
            use lib_ml_fmm_type
            implicit none
            ! dummy
            type(lib_ml_fmm_B_i), intent(in) :: B_i_1
            type(lib_tree_spatial_point), intent(in) :: x_1
            type(lib_tree_spatial_point), intent(in) :: x_2
            type(lib_ml_fmm_A_i) :: A_i_2
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
            use lib_ml_fmm_type
            implicit none
            ! dummy
            type(lib_ml_fmm_B_i), intent(in) :: B_i_1
            type(lib_tree_spatial_point), intent(in) :: x_1
            type(lib_tree_spatial_point), intent(in) :: x_2
            type(lib_ml_fmm_B_i) :: B_i_2
        end function
    end interface



    contains

    ! Procedure
    ! 1. call lib_tree_constructor
    !       - determine s_opt
    !       - determine l_min and l_max
    ! 2. create ml fmm data set
    !       - C per box and per level (l_max up to l_min)
    !       - D per box and per level (l_min down to l_max)
    subroutine lib_ml_fmm_constructor(data)
        implicit none
        type(lib_tree_data_element), dimension(:) :: data



    end subroutine lib_ml_fmm_constructor

    ! clean up
    ! - coefficents
    ! - Tree
    subroutine lib_ml_fmm_destructor()
        implicit none
        call lib_tree_destructor()

    end subroutine lib_ml_fmm_destructor

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
    ! Arguments
    ! ----
    !   h_get_B_i: function handle
    !
    !
    function lib_ml_fmm_calculate_upward_pass_step_1(h_get_B_i) result(C)
        implicit none
        ! dummy
        procedure(lib_ml_fmm_get_B_i) :: h_get_B_i
        type(lib_ml_fmm_C) :: C

        ! auxiliaray
        type(lib_tree_spatial_point) :: x_c
        type(lib_ml_fmm_B_i) :: B_i
        type(lib_tree_data_element) :: data_element

        B_i = h_get_B_i(x_c, data_element)



    end function lib_ml_fmm_calculate_upward_pass_step_1

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
    function lib_ml_fmm_calculate_upward_pass_step_2(C_l_plus_1) result(C_l)
        implicit none
        ! dummy
        type(lib_ml_fmm_C), intent(in) :: C_l_plus_1
        type(lib_ml_fmm_C) :: C_l

    end function lib_ml_fmm_calculate_upward_pass_step_2

    ! Downward pass - step 1
    !
    ! Procedure
    ! ----
    !   for each box at l_min, ..., l_max
    !       1.
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

end module lib_ml_fmm
