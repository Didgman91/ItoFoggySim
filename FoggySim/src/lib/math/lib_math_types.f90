module lib_math_types
    implicit none

    type spherical_coordinate_real_type
        double precision :: rho
        double precision :: theta
        double precision :: phi
    end type spherical_coordinate_real_type

    type list_spherical_coordinate_real_type
        type(spherical_coordinate_real_type), dimension(:), allocatable :: coordinate
    end type list_spherical_coordinate_real_type

    type spherical_coordinate_cmplx_type
        complex(kind=8) :: rho
        complex(kind=8) :: theta
        complex(kind=8) :: phi
    end type spherical_coordinate_cmplx_type

    type list_spherical_coordinate_cmplx_type
        type(spherical_coordinate_cmplx_type), dimension(:), allocatable :: coordinate
    end type list_spherical_coordinate_cmplx_type

end module lib_math_types
