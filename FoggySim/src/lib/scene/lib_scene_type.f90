module lib_scene_type
    use libmath
    implicit none

    public

    !                 z
    !                 ^
    !             K_j |
    !                 o--> x
    !                ^
    !               /
    !           z  /d_0_j
    !           ^ /
    !       K_o |/
    !           --> x
    !
    ! K_o: object coordinate system
    ! K_j: sphere coordinate system
    !   o: j-th sphere
    !
    type lib_scene_object_sphere_type
        type(cartesian_coordinate_real_type) :: d_o_j ! [m]
        integer :: sphere_parameter_index
        double precision :: radius  ! [m]
    end type lib_scene_object_sphere_type


    !
    !           z o o o o
    !           ^  o o o
    !       K_o | o o o o
    !           --> x
    !
    !   K_o: object coordinate system
    !   o: spheres at the hcp_lattice_coordiantes
    !
    type lib_scene_object_hcp_cuboid
        type(cartesian_coordinate_real_type) :: d_o_j
        double precision :: sphere_radius
        type(cartesian_coordinate_real_type), dimension(:,:,:), allocatable :: hcp_lattice_coordiantes
    end type

    type lib_scene_object_hcp_sphere

    end type

    !                 z
    !                 ^
    !             K_o |
    !                 --> x
    !                ^
    !               /
    !           z  /d_0_o
    !           ^ /
    !       K_0 |/
    !           --> x
    !
    ! K_0: world coordinate system
    ! K_o: object coordinate system
    !
    type lib_scene_opject_type
        type(cartesian_coordinate_real_type) :: d_0_o ! [m]
        type(lib_scene_object_sphere_type), dimension(:), allocatable :: sphere
        type(lib_scene_object_hcp_cuboid), dimension(:), allocatable :: hcp_cuboid
    end type

end module lib_scene_type
