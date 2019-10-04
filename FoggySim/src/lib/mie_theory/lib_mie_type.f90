module lib_mie_type
    use libmath
    implicit none

    public

    type lib_mie_illumination_plane_wave_parameter
        double precision :: g ! e_x_field_0 / e_field_0
        type(cartesian_coordinate_real_type) :: d_0_i ! [m]
        type(cartesian_coordinate_real_type) :: wave_vector_0 ! |k| = 2 Pi / lambda, wave_vector = [1/m]
    end type lib_mie_illumination_plane_wave_parameter

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
        integer :: type    ! 1: plane wave
        double precision :: e_field_0 ! [V/m]
        double precision :: lambda_0 ! wave_length_vaccum
        type(lib_mie_illumination_plane_wave_parameter), dimension(:), allocatable :: plane_wave
    end type lib_mie_illumination_parameter

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

    ! Depending on the initial definition of the electromagnetic wave,
    ! the spherical harmonics have to be calculated differently.
    !
    ! case: exp(i (k x + omega t) // actual case
    !   z_selector_incident_wave = 1
    !   z_selector_scatterd_wave = 3
    !   z_selector_translation = 3
    !
    ! case: exp(-i (k x + omega t)
    !   z_selector_incident_wave = 2
    !   z_selector_scatterd_wave = 4
    !   z_selector_translation = 4
    !
    ! n_range
    !   Defines the minimum and maximum degree for the entire simulation.
    !
    type lib_mie_vector_spherical_harmonics_type
        integer(kind=1) :: z_selector_incident_wave
        integer(kind=1) :: z_selector_scatterd_wave
        integer(kind=1) :: z_selector_translation
        integer, dimension(2) :: n_range
    end type lib_mie_vector_spherical_harmonics_type

    ! simulation parameter
    type lib_mie_simulation_parameter_type
        double precision :: refractive_index_medium
        type(lib_mie_illumination_parameter) :: illumination
        type(lib_mie_sphere_type), dimension(:), allocatable :: sphere_list
        type(lib_mie_sphere_parameter_type), dimension(:), allocatable :: sphere_parameter_list
        type(lib_mie_vector_spherical_harmonics_type) :: spherical_harmonics
    end type lib_mie_simulation_parameter_type

end module lib_mie_type
