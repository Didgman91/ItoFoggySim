module lib_math_type_operator
    use lib_math_type
    implicit none

    ! ----- operator -----
    interface operator (+)

    end interface

    contains

        function lib_math_spherical_operator_add(lhs, rhs) result(rv)
            use ml_fmm_type
            implicit none
            ! dummy
            type (spherical_coordinate_cmplx_type), dimension(:), intent(in) :: lhs
            type (spherical_coordinate_cmplx_type), dimension(size(lhs)), intent(in) :: rhs
            type (spherical_coordinate_cmplx_type), dimension(size(lhs)) :: rv



        end function

end module lib_math_type_operator
