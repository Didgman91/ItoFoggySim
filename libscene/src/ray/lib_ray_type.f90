module lib_ray_type
    use libmath
    use lib_scene_object_type
    implicit none

    type ray_type
        type(cartesian_coordinate_real_type) :: p    ! ray origin
        type(cartesian_coordinate_real_type) :: n    ! ray direction
        real length                                  ! ray path length
    end type ray_type

    type ray_object_type
        type(lib_scene_opject_type) :: object
        type(cartesian_coordinate_real_type) :: intersection_point
        logical :: apply_scattering
    end type ray_object_type


end module lib_ray_type
