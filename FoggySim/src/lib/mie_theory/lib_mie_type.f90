module lib_mie_type
    use libmath
    implicit none

    public

    ! illumination parameter
    !
    !             _________
    !             ___k^____
    !             ____|____
    !             _________ plane wave
    !                 z
    !                 ^
    !             K_i |
    !                 --> x
    !                ^
    !               /
    !           z  /d_0_i
    !           ^ /
    !       K_0 |/
    !           --> x
    !
    ! K_0: world coordinate system
    ! K_i: illumination coordinate system
    !
    type lib_mie_illumination_parameter
        integer :: illumination_type    ! 1: plane wave
        double precision :: lambda_0 ! wave_length_vaccum
        type(cartesian_coordinate_real_type) :: d_0_i
        type(cartesian_coordinate_real_type) :: wave_vector_0 ! |k| = 2 Pi / lambda
    end type lib_mie_illumination_parameter

    ! sphere
    !                 z
    !                 ^
    !             K_j |
    !                 o--> x
    !                ^
    !               /
    !           z  /d_0_j
    !           ^ /
    !       K_0 |/
    !           --> x
    !
    ! K_0: world coordinate system
    ! K_j: sphere coordinate system
    !   o: j-th sphere
    !
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

    ! simulation parameter
    type lib_mie_simulation_parameter_type
        double precision :: refractive_index_medium
        type(lib_mie_illumination_parameter) :: illumination
        type(lib_mie_sphere_type), dimension(:), allocatable :: sphere_list
        type(lib_mie_sphere_parameter_type), dimension(:), allocatable :: sphere_parameter_list
    end type lib_mie_simulation_parameter_type

end module lib_mie_type
