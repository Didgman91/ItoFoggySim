module lib_math_type
    implicit none

    integer(kind=1), parameter :: lib_math_type_kind = 8

    type spherical_coordinate_real_type
        real(kind=lib_math_type_kind) :: rho
        real(kind=lib_math_type_kind) :: theta
        real(kind=lib_math_type_kind) :: phi
    end type spherical_coordinate_real_type

    type list_spherical_coordinate_real_type
        type(spherical_coordinate_real_type), dimension(:), allocatable :: coordinate
    end type list_spherical_coordinate_real_type

    type spherical_coordinate_cmplx_type
        complex(kind=lib_math_type_kind) :: rho
        complex(kind=lib_math_type_kind) :: theta
        complex(kind=lib_math_type_kind) :: phi
    end type spherical_coordinate_cmplx_type

    type list_spherical_coordinate_cmplx_type
        type(spherical_coordinate_cmplx_type), dimension(:), allocatable :: coordinate
    end type list_spherical_coordinate_cmplx_type

end module lib_math_type
