module lib_mie_type
    use libmath
    implicit none

    public

    ! sphere
    type sphere_type
        type(cartesian_coordinate_real_type) :: d_0_j
        integer :: sphere_parameter_index
    end type sphere_type

    ! sphere parameter
    type sphere_parameter_type
        ! initial values
        double complex :: refractive_index
        double precision :: radius
        ! calculated values
        double precision :: size_parameter ! = k * x = 2 Pi * n_medium / lambda * x
        integer, dimension(2) :: n_range
        type(list_cmplx) :: a_n
        type(list_cmplx) :: b_n
    end type sphere_parameter_type

end module lib_mie_type
