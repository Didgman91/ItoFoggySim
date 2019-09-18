module lib_mie_type
    use libmath
    implicit none

    public

    ! simulation parameter
    type lib_mie_simulation_parameter_type
        double precision :: refractive_index_medium
        double precision :: lambda ! wave_length_vaccum
        type(cartesian_coordinate_real_type) :: wave_vector ! = 2 Pi * refractive_index_medium / lambda
    end type lib_mie_simulation_parameter_type

    ! sphere
    type lib_mie_sphere_type
        type(cartesian_coordinate_real_type) :: d_0_j
        integer :: sphere_parameter_index
        type(list_list_cmplx) :: a_nm
        type(list_list_cmplx) :: b_nm
    end type lib_mie_sphere_type

    ! sphere parameter
    type lib_mie_sphere_parameter_type
        double complex :: refractive_index
        double precision :: radius
        double precision :: size_parameter ! = k * x = 2 Pi * n_medium / lambda * x
        integer, dimension(2) :: n_range
        type(list_cmplx) :: a_n
        type(list_cmplx) :: b_n
    end type lib_mie_sphere_parameter_type

end module lib_mie_type
