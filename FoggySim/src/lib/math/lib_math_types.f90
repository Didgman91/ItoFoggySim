module lib_math_types
    implicit none

    type spherical_coordinate_type
        double precision :: rho
        double precision :: theta
        double precision :: phi
    end type spherical_coordinate_type

    type list_spherical_coordinate_type
        type(spherical_coordinate_type), dimension(:), allocatable :: coordinate
    end type list_spherical_coordinate_type
end module lib_math_types
