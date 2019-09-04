module lib_math_type
    implicit none

    integer(kind=1), parameter :: lib_math_type_kind = 8

    ! cartesian coordinates
    type cartesian_coordinate_real_type
        real(kind=lib_math_type_kind) :: x
        real(kind=lib_math_type_kind) :: y
        real(kind=lib_math_type_kind) :: z
    end type cartesian_coordinate_real_type

    type list_cartesian_coordinate_real_type
        type(cartesian_coordinate_real_type), dimension(:), allocatable :: coordinate
    end type list_cartesian_coordinate_real_type

    type cartesian_coordinate_cmplx_type
        complex(kind=lib_math_type_kind) :: x
        complex(kind=lib_math_type_kind) :: y
        complex(kind=lib_math_type_kind) :: z
    end type cartesian_coordinate_cmplx_type

    type list_cartesian_coordinate_cmplx_type
        type(cartesian_coordinate_cmplx_type), dimension(:), allocatable :: coordinate
    end type list_cartesian_coordinate_cmplx_type

    ! spherical coordinates
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

    ! lists
    type list_logical
        logical, dimension(:), allocatable :: item
    end type

    type list_list_logical
        type(list_logical), dimension(:), allocatable :: item
    end type

    type list_integer_sys
        integer, dimension(:), allocatable :: item
    end type

    type list_list_integer_sys
        type(list_integer_sys), dimension(:), allocatable :: item
    end type

    type list_integer
        integer(kind=lib_math_type_kind), dimension(:), allocatable :: item
    end type

    type list_list_integer
        type(list_integer), dimension(:), allocatable :: item
    end type

    type list_real
        real(kind=lib_math_type_kind), dimension(:), allocatable :: item
    end type

    type list_list_real
        type(list_real), dimension(:), allocatable :: item
    end type

    type list_list_list_real
        type(list_list_real), dimension(:), allocatable :: item
    end type

    type list_4_real
        type(list_list_list_real), dimension(:), allocatable :: item
    end type

    type list_cmplx
        complex(kind=lib_math_type_kind), dimension(:), allocatable :: item
    end type

    type list_list_cmplx
        type(list_cmplx), dimension(:), allocatable :: item
    end type

    type list_list_list_cmplx
        type(list_list_cmplx), dimension(:), allocatable :: item
    end type

    type list_4_cmplx
        type(list_list_list_cmplx), dimension(:), allocatable :: item
    end type

end module lib_math_type
