module lib_octree_helper_functions
    implicit none

    private

    integer, parameter :: octree_integer_kind = 4
    integer, parameter :: fmm_d = 3 ! dimensions

    contains

    function lib_octree_hf_get_neighbourhood_size_S(R_c) result (k)
    implicit none
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
        double precision, intent (in) :: R_c
        double precision :: buffer
        integer(kind=octree_integer_kind) :: k

        buffer = 0.5 * (R_c * sqrt(fmm_d) -1 )

        k = ceiling(buffer)

    end function lib_octree_hf_get_neighbourhood_size_S

    function lib_octree_hf_get_neighbourhood_size_R(r_c) result (k)
    implicit none
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
        double precision, intent (in) :: r_c
        double precision :: buffer
        integer(kind=octree_integer_kind) :: k

        buffer = 0.5 * (1/r_c * sqrt(fmm_d) -1 )

        k = ceiling(buffer)

    end function lib_octree_hf_get_neighbourhood_size_S

    function lib_octree_hf_get_neighbourhood_size(R_c, r_c) result (k)
    implicit none
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
        double precision, intent (in) :: r_c
        double precision :: buffer
        integer(kind=octree_integer_kind) :: k

        buffer = 0.5 * (max(1/r_c, R_c) * sqrt(fmm_d) -1 )

        k = ceiling(buffer)

    end function lib_octree_hf_get_neighbourhood_size_S



end module lib_octree_helper_functions
