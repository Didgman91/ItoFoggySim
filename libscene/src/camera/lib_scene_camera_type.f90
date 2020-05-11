module lib_scene_camera_type
    use libmath
    implicit none

    type colour_type
        real :: R
        real :: G
        real :: B
    end type colour_type

    type pixel_type
        type(colour_type) :: filter
        type(colour_type) :: source
        logical :: time_of_flight
    end type pixel_type

end module lib_scene_camera_type
