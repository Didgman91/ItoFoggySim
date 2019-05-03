! LIB: Mulitlevel Fast Multipole Method - Helper Functions
!
module lib_ml_fmm_helper_functions
    use lib_tree_type
    implicit none

    private

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

        end function lib_ml_fmm_hf_check_validity_translation_SS

end module lib_ml_fmm_helper_functions
