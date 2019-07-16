module lib_math_type_operator
    use lib_math_type
    implicit none

    private

    ! --- public ---
    public :: operator (+)
    public :: operator (-)
    public :: operator (*)
    public :: operator (/)
    public :: spherical_abs
    public :: init_list

    public :: lib_math_type_operator_test_functions

    ! ---- operator ----
    interface operator (+)
        module procedure lib_math_spherical_operator_add
        module procedure lib_math_spherical_operator_add_array
        module procedure lib_math_spherical_operator_0d_add_array
        module procedure lib_math_spherical_operator_array_add_0d
        module procedure lib_math_list_spherical_operator_add_array
        module procedure lib_math_list_spherical_operator_0d_add_array
        module procedure lib_math_list_spherical_operator_array_add_0d
    end interface

    interface operator (-)
        module procedure lib_math_spherical_operator_sub
        module procedure lib_math_spherical_operator_sub_array
        module procedure lib_math_spherical_operator_array_sub_0d
        module procedure lib_math_list_spherical_operator_sub_array
        module procedure lib_math_list_spherical_operator_array_sub_0d
    end interface

    interface operator (*)
        module procedure lib_math_spherical_operator_scalar_real_mul
        module procedure lib_math_spherical_operator_scalar_cmplx_mul
        module procedure lib_math_spherical_operator_mul_scalar_real
        module procedure lib_math_spherical_operator_mul_scalar_cmplx
        module procedure lib_math_spherical_operator_array_real_mul_array
        module procedure lib_math_spherical_operator_array_cmplx_mul_array
        module procedure lib_math_spherical_operator_0d_real_mul_array
        module procedure lib_math_spherical_operator_0d_cmplx_mul_array
        module procedure lib_math_list_spherical_operator_array_real_mul_array
        module procedure lib_math_list_spherical_operator_array_cmplx_mul_array
        module procedure lib_math_list_spherical_operator_real_mul_array
        module procedure lib_math_list_spherical_operator_cmplx_mul_array
    end interface

    interface operator (/)
        module procedure lib_math_spherical_operator_divide_by_scalar_real
        module procedure lib_math_spherical_operator_divide_by_scalar_cmplx
        module procedure lib_math_spherical_operator_array_divide_by_scalar_array_real
        module procedure lib_math_spherical_operator_array_divide_by_scalar_array_cmplx
        module procedure lib_math_spherical_operator_array_divide_by_0d_real
        module procedure lib_math_spherical_operator_array_divide_by_0d_cmplx
        module procedure lib_math_list_spherical_operator_array_divide_by_real_array
        module procedure lib_math_list_spherical_operator_array_divide_by_cmplx_array
        module procedure lib_math_list_spherical_operator_array_divide_by_real
        module procedure lib_math_list_spherical_operator_array_divide_by_cmplx
    end interface

    interface spherical_abs
        module procedure lib_math_spherical_operator_abs_cmplx
    end interface


    interface init_list
        module procedure lib_math_list_list_real_init
        module procedure lib_math_list_list_cmplx_init
        module procedure lib_math_list_spherical_coordinate_cmplx_type_init
    end interface

    contains

! ---- single ----
        function lib_math_spherical_operator_add(lhs, rhs) result(rv)
            implicit none
            ! dummy
            type (spherical_coordinate_cmplx_type), intent(in) :: lhs
            type (spherical_coordinate_cmplx_type), intent(in) :: rhs
            type (spherical_coordinate_cmplx_type) :: rv

            rv%rho = lhs%rho + rhs%rho
            rv%phi = lhs%phi + rhs%phi
            rv%theta = lhs%theta + rhs%theta

        end function

        function lib_math_spherical_operator_sub(lhs, rhs) result(rv)
            implicit none
            ! dummy
            type (spherical_coordinate_cmplx_type), intent(in) :: lhs
            type (spherical_coordinate_cmplx_type), intent(in) :: rhs
            type (spherical_coordinate_cmplx_type) :: rv

            rv%rho = lhs%rho - rhs%rho
            rv%phi = lhs%phi - rhs%phi
            rv%theta = lhs%theta - rhs%theta

        end function

        function lib_math_spherical_operator_scalar_real_mul(lhs, rhs) result(rv)
            implicit none
            ! dummy
            real(kind=8), intent(in) :: lhs
            type (spherical_coordinate_cmplx_type), intent(in) :: rhs
            type (spherical_coordinate_cmplx_type) :: rv

            rv%rho = lhs * rhs%rho
            rv%phi = lhs * rhs%phi
            rv%theta = lhs * rhs%theta

        end function

        function lib_math_spherical_operator_scalar_cmplx_mul(lhs, rhs) result(rv)
            implicit none
            ! dummy
            complex(kind=8), intent(in) :: lhs
            type (spherical_coordinate_cmplx_type), intent(in) :: rhs
            type (spherical_coordinate_cmplx_type) :: rv

            rv%rho = lhs * rhs%rho
            rv%phi = lhs * rhs%phi
            rv%theta = lhs * rhs%theta

        end function

        function lib_math_spherical_operator_mul_scalar_real(lhs, rhs) result(rv)
            implicit none
            ! dummy
            type (spherical_coordinate_cmplx_type), intent(in) :: lhs
            real(kind=8), intent(in) :: rhs
            type (spherical_coordinate_cmplx_type) :: rv

            rv%rho = lhs%rho * rhs
            rv%phi = lhs%phi * rhs
            rv%theta = lhs%theta * rhs

        end function

        function lib_math_spherical_operator_mul_scalar_cmplx(lhs, rhs) result(rv)
            implicit none
            ! dummy
            type (spherical_coordinate_cmplx_type), intent(in) :: lhs
            complex(kind=8), intent(in) :: rhs
            type (spherical_coordinate_cmplx_type) :: rv

            rv%rho = lhs%rho * rhs
            rv%phi = lhs%phi * rhs
            rv%theta = lhs%theta * rhs

        end function

        function lib_math_spherical_operator_divide_by_scalar_real(lhs, rhs) result(rv)
            implicit none
            ! dummy
            type (spherical_coordinate_cmplx_type), intent(in) :: lhs
            real(kind=8), intent(in) :: rhs
            type (spherical_coordinate_cmplx_type) :: rv

            rv%rho = lhs%rho / rhs
            rv%phi = lhs%phi / rhs
            rv%theta = lhs%theta / rhs

        end function

        function lib_math_spherical_operator_divide_by_scalar_cmplx(lhs, rhs) result(rv)
            implicit none
            ! dummy
            type (spherical_coordinate_cmplx_type), intent(in) :: lhs
            complex(kind=8), intent(in) :: rhs
            type (spherical_coordinate_cmplx_type) :: rv

            rv%rho = lhs%rho / rhs
            rv%phi = lhs%phi / rhs
            rv%theta = lhs%theta / rhs

        end function

        function lib_math_spherical_operator_abs_cmplx(rhs) result(rv)
            implicit none
            ! dummy
            type (spherical_coordinate_cmplx_type), intent(in) :: rhs

            complex(kind=8) :: rv

            rv = sqrt(rhs%rho*rhs%rho + rhs%phi*rhs%phi + rhs%theta*rhs%theta)

        end function

! ---- array ----
        function lib_math_spherical_operator_add_array(lhs, rhs) result(rv)
            implicit none
            ! dummy
            type (spherical_coordinate_cmplx_type), dimension(:), intent(in) :: lhs
            type (spherical_coordinate_cmplx_type), dimension(lbound(lhs,1):ubound(lhs,1)), intent(in) :: rhs
            type (spherical_coordinate_cmplx_type), dimension(lbound(lhs,1):ubound(lhs,1)) :: rv

            ! auxiliary
            integer :: i

            do i=lbound(lhs,1), ubound(lhs,1)
                rv(i) = lhs(i) + rhs(i)
            end do

        end function

        function lib_math_spherical_operator_sub_array(lhs, rhs) result(rv)
            implicit none
            ! dummy
            type (spherical_coordinate_cmplx_type), dimension(:), intent(in) :: lhs
            type (spherical_coordinate_cmplx_type), dimension(lbound(lhs,1):ubound(lhs,1)), intent(in) :: rhs
            type (spherical_coordinate_cmplx_type), dimension(lbound(lhs,1):ubound(lhs,1)) :: rv

            ! auxiliary
            integer :: i

            do i=lbound(lhs,1), ubound(lhs,1)
                rv(i) = lhs(i) - rhs(i)
            end do

        end function

        function lib_math_spherical_operator_0d_add_array(lhs, rhs) result(rv)
            implicit none
            ! dummy
            type (spherical_coordinate_cmplx_type), intent(in) :: lhs
            type (spherical_coordinate_cmplx_type), dimension(:), intent(in) :: rhs
            type (spherical_coordinate_cmplx_type), dimension(lbound(rhs,1):ubound(rhs,1)) :: rv

            ! auxiliary
            integer :: i

            do i=lbound(rhs,1), ubound(rhs,1)
                rv(i) = lhs + rhs(i)
            end do

        end function

        function lib_math_spherical_operator_array_add_0d(lhs, rhs) result(rv)
            implicit none
            ! dummy
            type (spherical_coordinate_cmplx_type), dimension(:), intent(in) :: lhs
            type (spherical_coordinate_cmplx_type), intent(in) :: rhs
            type (spherical_coordinate_cmplx_type), dimension(lbound(lhs,1):ubound(lhs,1)) :: rv

            ! auxiliary
            integer :: i

            do i=lbound(lhs,1), ubound(lhs,1)
                rv(i) = lhs(i) + rhs
            end do

        end function

        function lib_math_spherical_operator_array_sub_0d(lhs, rhs) result(rv)
            implicit none
            ! dummy
            type (spherical_coordinate_cmplx_type), dimension(:), intent(in) :: lhs
            type (spherical_coordinate_cmplx_type), intent(in) :: rhs
            type (spherical_coordinate_cmplx_type), dimension(lbound(lhs,1):ubound(lhs,1)) :: rv

            ! auxiliary
            integer :: i

            do i=lbound(lhs,1), ubound(lhs,1)
                rv(i) = lhs(i) - rhs
            end do

        end function

        function lib_math_spherical_operator_array_real_mul_array(lhs, rhs) result(rv)
            implicit none
            ! dummy
            real(kind=lib_math_type_kind), dimension(:), intent(in) :: lhs
            type (spherical_coordinate_cmplx_type), dimension(lbound(lhs,1):ubound(lhs,1)), intent(in) :: rhs
            type (spherical_coordinate_cmplx_type), dimension(lbound(lhs,1):ubound(lhs,1)) :: rv

            ! auxiliary
            integer :: i

            do i=lbound(lhs,1), ubound(lhs,1)
                rv(i) = lhs(i) * rhs(i)
            end do

        end function

        function lib_math_spherical_operator_array_cmplx_mul_array(lhs, rhs) result(rv)
            implicit none
            ! dummy
            complex(kind=lib_math_type_kind), dimension(:), intent(in) :: lhs
            type (spherical_coordinate_cmplx_type), dimension(lbound(lhs,1):ubound(lhs,1)), intent(in) :: rhs
            type (spherical_coordinate_cmplx_type), dimension(lbound(lhs,1):ubound(lhs,1)) :: rv

            ! auxiliary
            integer :: i

            do i=lbound(lhs,1), ubound(lhs,1)
                rv(i) = lhs(i) * rhs(i)
            end do

        end function

        function lib_math_spherical_operator_0d_real_mul_array(lhs, rhs) result(rv)
            implicit none
            ! dummy
            real(kind=lib_math_type_kind), intent(in) :: lhs
            type (spherical_coordinate_cmplx_type), dimension(:), intent(in) :: rhs
            type (spherical_coordinate_cmplx_type), dimension(lbound(rhs,1):ubound(rhs,1)) :: rv

            ! auxiliary
            integer :: i

            do i=lbound(rhs,1), ubound(rhs,1)
                rv(i) = lhs * rhs(i)
            end do

        end function

        function lib_math_spherical_operator_0d_cmplx_mul_array(lhs, rhs) result(rv)
            implicit none
            ! dummy
            complex(kind=lib_math_type_kind), intent(in) :: lhs
            type (spherical_coordinate_cmplx_type), dimension(:), intent(in) :: rhs
            type (spherical_coordinate_cmplx_type), dimension(lbound(rhs,1):ubound(rhs,1)) :: rv

            ! auxiliary
            integer :: i

            do i=lbound(rhs,1), ubound(rhs,1)
                rv(i) = lhs * rhs(i)
            end do

        end function

        function lib_math_spherical_operator_array_divide_by_scalar_array_real(lhs, rhs) result(rv)
            implicit none
            ! dummy
            type (spherical_coordinate_cmplx_type), dimension(:), intent(in) :: lhs
            real(kind=lib_math_type_kind), dimension(lbound(lhs,1):ubound(lhs,1)), intent(in) :: rhs
            type (spherical_coordinate_cmplx_type), dimension(lbound(lhs,1):ubound(lhs,1)) :: rv

            ! auxiliary
            integer :: i

            do i=lbound(lhs,1), ubound(lhs,1)
                rv(i) = lhs(i) / rhs(i)
            end do

        end function

        function lib_math_spherical_operator_array_divide_by_scalar_array_cmplx(lhs, rhs) result(rv)
            implicit none
            ! dummy
            type (spherical_coordinate_cmplx_type), dimension(:), intent(in) :: lhs
            complex(kind=lib_math_type_kind), dimension(lbound(lhs,1):ubound(lhs,1)), intent(in) :: rhs
            type (spherical_coordinate_cmplx_type), dimension(lbound(lhs,1):ubound(lhs,1)) :: rv

            ! auxiliary
            integer :: i

            do i=lbound(lhs,1), ubound(lhs,1)
                rv(i) = lhs(i) / rhs(i)
            end do

        end function

        function lib_math_spherical_operator_array_divide_by_0d_real(lhs, rhs) result(rv)
            implicit none
            ! dummy
            type (spherical_coordinate_cmplx_type), dimension(:), intent(in) :: lhs
            real(kind=lib_math_type_kind), intent(in) :: rhs
            type (spherical_coordinate_cmplx_type), dimension(lbound(lhs,1):ubound(lhs,1)) :: rv

            ! auxiliary
            integer :: i

            do i=lbound(lhs,1), ubound(lhs,1)
                rv(i) = lhs(i) / rhs
            end do

        end function

        function lib_math_spherical_operator_array_divide_by_0d_cmplx(lhs, rhs) result(rv)
            implicit none
            ! dummy
            type (spherical_coordinate_cmplx_type), dimension(:), intent(in) :: lhs
            complex(kind=lib_math_type_kind), intent(in) :: rhs
            type (spherical_coordinate_cmplx_type), dimension(lbound(lhs,1):ubound(lhs,1)) :: rv

            ! auxiliary
            integer :: i

            do i=lbound(lhs,1), ubound(lhs,1)
                rv(i) = lhs(i) / rhs
            end do

        end function

! ---- list_spherical_coordinate ----
        function lib_math_list_spherical_operator_add_array(lhs, rhs) result(rv)
            implicit none
            ! dummy
            type (list_spherical_coordinate_cmplx_type), intent(in) :: lhs
            type (list_spherical_coordinate_cmplx_type), intent(in) :: rhs
            type (list_spherical_coordinate_cmplx_type) :: rv

            ! auxiliary
            integer :: i

            if (lbound(lhs%coordinate,1) .eq. lbound(rhs%coordinate,1) .and. &
                ubound(lhs%coordinate,1) .eq. ubound(rhs%coordinate,1)) then

                allocate(rv%coordinate(lbound(rhs%coordinate,1):ubound(rhs%coordinate,1)))

                do i=lbound(lhs%coordinate,1), ubound(lhs%coordinate,1)
                    rv%coordinate(i) = lhs%coordinate(i) + rhs%coordinate(i)
                end do
            end if

        end function

        function lib_math_list_spherical_operator_sub_array(lhs, rhs) result(rv)
            implicit none
            ! dummy
            type (list_spherical_coordinate_cmplx_type), intent(in) :: lhs
            type (list_spherical_coordinate_cmplx_type), intent(in) :: rhs
            type (list_spherical_coordinate_cmplx_type) :: rv

            ! auxiliary
            integer :: i

            if (lbound(lhs%coordinate,1) .eq. lbound(rhs%coordinate,1) .and. &
                ubound(lhs%coordinate,1) .eq. ubound(rhs%coordinate,1)) then

                allocate(rv%coordinate(lbound(rhs%coordinate,1):ubound(rhs%coordinate,1)))

                do i=lbound(lhs%coordinate,1), ubound(lhs%coordinate,1)
                    rv%coordinate(i) = lhs%coordinate(i) - rhs%coordinate(i)
                end do
            end if
        end function

        function lib_math_list_spherical_operator_0d_add_array(lhs, rhs) result(rv)
            implicit none
            ! dummy
            type (spherical_coordinate_cmplx_type), intent(in) :: lhs
            type (list_spherical_coordinate_cmplx_type), intent(in) :: rhs
            type (list_spherical_coordinate_cmplx_type) :: rv

            ! auxiliary
            integer :: i

            allocate(rv%coordinate(lbound(rhs%coordinate,1):ubound(rhs%coordinate,1)))

            do i=lbound(rhs%coordinate,1), ubound(rhs%coordinate,1)
                rv%coordinate(i) = lhs + rhs%coordinate(i)
            end do

        end function

        function lib_math_list_spherical_operator_array_add_0d(lhs, rhs) result(rv)
            implicit none
            ! dummy
            type (list_spherical_coordinate_cmplx_type), intent(in) :: lhs
            type (spherical_coordinate_cmplx_type), intent(in) :: rhs
            type (list_spherical_coordinate_cmplx_type) :: rv

            ! auxiliary
            integer :: i

            allocate(rv%coordinate(lbound(lhs%coordinate,1):ubound(lhs%coordinate,1)))

            do i=lbound(lhs%coordinate,1), ubound(lhs%coordinate,1)
                rv%coordinate(i) = lhs%coordinate(i) + rhs
            end do

        end function

        function lib_math_list_spherical_operator_array_sub_0d(lhs, rhs) result(rv)
            implicit none
            ! dummy
            type (list_spherical_coordinate_cmplx_type), intent(in) :: lhs
            type (spherical_coordinate_cmplx_type), intent(in) :: rhs
            type (list_spherical_coordinate_cmplx_type) :: rv

            ! auxiliary
            integer :: i

            allocate(rv%coordinate(lbound(lhs%coordinate,1):ubound(lhs%coordinate,1)))

            do i=lbound(lhs%coordinate,1), ubound(lhs%coordinate,1)
                rv%coordinate(i) = lhs%coordinate(i) - rhs
            end do

        end function

        function lib_math_list_spherical_operator_array_real_mul_array(lhs, rhs) result(rv)
            implicit none
            ! dummy
            real(kind=lib_math_type_kind), dimension(:), intent(in) :: lhs
            type (list_spherical_coordinate_cmplx_type), intent(in) :: rhs
            type (list_spherical_coordinate_cmplx_type) :: rv

            ! auxiliary
            integer :: i

            if (lbound(lhs,1) .eq. lbound(rhs%coordinate,1) .and. &
                ubound(lhs,1) .eq. ubound(rhs%coordinate,1)) then

                allocate(rv%coordinate(lbound(lhs,1):ubound(lhs,1)))

                do i=lbound(lhs,1), ubound(lhs,1)
                    rv%coordinate(i) = lhs(i) * rhs%coordinate(i)
                end do
            else
                print *, "lib_math_list_spherical_operator_array_real_mul_array: ERROR"
            end if

        end function

        function lib_math_list_spherical_operator_array_cmplx_mul_array(lhs, rhs) result(rv)
            implicit none
            ! dummy
            complex(kind=lib_math_type_kind), dimension(:), intent(in) :: lhs
            type (list_spherical_coordinate_cmplx_type), intent(in) :: rhs
            type (list_spherical_coordinate_cmplx_type) :: rv

            ! auxiliary
            integer :: i

            if (lbound(lhs,1) .eq. lbound(rhs%coordinate,1) .and. &
                ubound(lhs,1) .eq. ubound(rhs%coordinate,1)) then

                allocate(rv%coordinate(lbound(rhs%coordinate,1):ubound(rhs%coordinate,1)))

                do i=lbound(lhs,1), ubound(lhs,1)
                    rv%coordinate(i) = lhs(i) * rhs%coordinate(i)
                end do

            end if

        end function

        function lib_math_list_spherical_operator_real_mul_array(lhs, rhs) result(rv)
            implicit none
            ! dummy
            real(kind=lib_math_type_kind), intent(in) :: lhs
            type (list_spherical_coordinate_cmplx_type), intent(in) :: rhs
            type (list_spherical_coordinate_cmplx_type) :: rv

            ! auxiliary
            integer :: i

            allocate(rv%coordinate(lbound(rhs%coordinate,1):ubound(rhs%coordinate,1)))

            do i=lbound(rhs%coordinate,1), ubound(rhs%coordinate,1)
                rv%coordinate(i) = lhs * rhs%coordinate(i)
            end do

        end function

        function lib_math_list_spherical_operator_cmplx_mul_array(lhs, rhs) result(rv)
            implicit none
            ! dummy
            complex(kind=lib_math_type_kind), intent(in) :: lhs
            type (list_spherical_coordinate_cmplx_type), intent(in) :: rhs
            type (list_spherical_coordinate_cmplx_type) :: rv

            ! auxiliary
            integer :: i

            allocate(rv%coordinate(lbound(rhs%coordinate,1):ubound(rhs%coordinate,1)))

            do i=lbound(rhs%coordinate,1), ubound(rhs%coordinate,1)
                rv%coordinate(i) = lhs * rhs%coordinate(i)
            end do

        end function

        function lib_math_list_spherical_operator_array_divide_by_real_array(lhs, rhs) result(rv)
            implicit none
            ! dummy
            type (list_spherical_coordinate_cmplx_type), intent(in) :: lhs
            real(kind=lib_math_type_kind), dimension(:), intent(in) :: rhs
            type (list_spherical_coordinate_cmplx_type) :: rv

            ! auxiliary
            integer :: i

            if (lbound(lhs%coordinate,1) .eq. lbound(rhs,1) .and. &
                ubound(lhs%coordinate,1) .eq. ubound(rhs,1)) then

                allocate(rv%coordinate(lbound(lhs%coordinate,1):ubound(lhs%coordinate,1)))

                do i=lbound(lhs%coordinate,1), ubound(lhs%coordinate,1)
                    rv%coordinate(i) = lhs%coordinate(i) / rhs(i)
                end do

            end if

        end function

        function lib_math_list_spherical_operator_array_divide_by_cmplx_array(lhs, rhs) result(rv)
            implicit none
            ! dummy
            type (list_spherical_coordinate_cmplx_type), intent(in) :: lhs
            complex(kind=lib_math_type_kind), dimension(:), intent(in) :: rhs
            type (list_spherical_coordinate_cmplx_type) :: rv

            ! auxiliary
            integer :: i

            if (lbound(lhs%coordinate,1) .eq. lbound(rhs,1) .and. &
                ubound(lhs%coordinate,1) .eq. ubound(rhs,1)) then

                allocate(rv%coordinate(lbound(lhs%coordinate,1):ubound(lhs%coordinate,1)))

                do i=lbound(lhs%coordinate,1), ubound(lhs%coordinate,1)
                    rv%coordinate(i) = lhs%coordinate(i) / rhs(i)
                end do

            end if

        end function

        function lib_math_list_spherical_operator_array_divide_by_real(lhs, rhs) result(rv)
            implicit none
            ! dummy
            type (list_spherical_coordinate_cmplx_type), intent(in) :: lhs
            real(kind=lib_math_type_kind), intent(in) :: rhs
            type (list_spherical_coordinate_cmplx_type) :: rv

            ! auxiliary
            integer :: i

            allocate(rv%coordinate(lbound(lhs%coordinate,1):ubound(lhs%coordinate,1)))

            do i=lbound(lhs%coordinate,1), ubound(lhs%coordinate,1)
                rv%coordinate(i) = lhs%coordinate(i) / rhs
            end do

        end function

        function lib_math_list_spherical_operator_array_divide_by_cmplx(lhs, rhs) result(rv)
            implicit none
            ! dummy
            type (list_spherical_coordinate_cmplx_type), intent(in) :: lhs
            complex(kind=lib_math_type_kind), intent(in) :: rhs
            type (list_spherical_coordinate_cmplx_type) :: rv

            ! auxiliary
            integer :: i

            allocate(rv%coordinate(lbound(lhs%coordinate,1):ubound(lhs%coordinate,1)))

            do i=lbound(lhs%coordinate,1), ubound(lhs%coordinate,1)
                rv%coordinate(i) = lhs%coordinate(i) / rhs
            end do

        end function

! ---- list_list ----

        ! Arguments
        ! ----
        !   list: type (list_list_real)
        !       derived data type to initialize
        !   fnu: integer
        !       start index
        !   n: integer
        !       number of elements, n .GE. 1
        !
        subroutine lib_math_list_list_real_init(list, fnu, n)
            implicit none
            ! dummy
            type (list_list_real), intent(inout) :: list
            integer, intent(in) :: fnu
            integer, intent(in) :: n

            ! auxiliary
            integer :: i

            if( allocated(list%item) ) then
                deallocate( list%item )
            end if

            allocate( list%item(fnu:fnu+n-1) )
            do i=fnu, fnu+n-1
                allocate (list%item(i)%item(-i:i) )
            end do

        end subroutine lib_math_list_list_real_init

        ! Arguments
        ! ----
        !   list: type (list_list_real)
        !       derived data type to initialize
        !   fnu: integer
        !       start index
        !   n: integer
        !       number of elements, n .GE. 1
        !
        subroutine lib_math_list_list_cmplx_init(list, fnu, n)
            implicit none
            ! dummy
            type (list_list_cmplx), intent(inout) :: list
            integer, intent(in) :: fnu
            integer, intent(in) :: n

            ! auxiliary
            integer :: i

            if( allocated(list%item) ) then
                deallocate( list%item )
            end if

            allocate( list%item(fnu:fnu+n-1) )
            do i=fnu, fnu+n-1
                allocate (list%item(i)%item(-i:i) )
            end do

        end subroutine lib_math_list_list_cmplx_init

        ! Arguments
        ! ----
        !   list: type (list_spherical_coordinate_cmplx_type)
        !       derived data type to initialize
        !   fnu: integer
        !       start index
        !   n: integer
        !       number of elements, n .GE. 1
        !
        subroutine lib_math_list_spherical_coordinate_cmplx_type_init(list, fnu, n)
            implicit none
            ! dummy
            type (list_spherical_coordinate_cmplx_type), dimension(:), allocatable, intent(inout) :: list
            integer, intent(in) :: fnu
            integer, intent(in) :: n

            ! auxiliary
            integer :: i

            if( allocated(list) ) then
                deallocate( list )
            end if

            allocate( list(fnu:fnu+n-1) )
            do i=fnu, fnu+n-1
                allocate (list(i)%coordinate(-i:i))
            end do

        end subroutine lib_math_list_spherical_coordinate_cmplx_type_init

        function lib_math_type_operator_test_functions() result (rv)
            implicit none
            ! dummy
            integer :: rv

            rv = 0

            if (.not. test_lib_math_list_spherical_operator_add_array()) then
                rv = rv + 1
            end if
            if (.not. test_lib_math_list_spherical_operator_sub_array()) then
                rv = rv + 1
            end if
            if (.not. test_lib_math_list_spherical_operator_0d_add_array()) then
                rv = rv + 1
            end if
            if (.not. test_lib_math_list_spherical_operator_array_add_0d()) then
                rv = rv + 1
            end if
            if (.not. test_lib_math_list_spherical_operator_array_sub_0d()) then
                rv = rv + 1
            end if
            if (.not. test_lib_math_list_spherical_operator_array_real_mul_array()) then
                rv = rv + 1
            end if
            if (.not. test_lib_math_list_spherical_operator_array_cmplx_mul_array()) then
                rv = rv + 1
            end if
            if (.not. test_lib_math_list_spherical_operator_real_mul_array()) then
                rv = rv + 1
            end if
            if (.not. test_lib_math_list_spherical_operator_cmplx_mul_array()) then
                rv = rv + 1
            end if
            if (.not. test_lib_math_list_spherical_operator_array_divide_by_r_array()) then
                rv = rv + 1
            end if
            if (.not. test_lib_math_list_spherical_operator_array_divide_by_c_array()) then
                rv = rv + 1
            end if
            if (.not. test_lib_math_list_spherical_operator_array_divide_by_real()) then
                rv = rv + 1
            end if
            if (.not. test_lib_math_list_spherical_operator_array_divide_by_cmplx()) then
                rv = rv + 1
            end if

            print *, "--------lib_math_type_operator_test_functions--------"
            if (rv == 0) then
                print *, "lib_math_type_operator_test_functions tests: OK"
            else
                print *, rv,"lib_math_type_operator_test_functions test(s) FAILED"
            end if
            print *, "-----------------------------------------------------"

            contains

                function get_test_values_1() result(rv)
                    implicit none
                    ! dummy
                    type (spherical_coordinate_cmplx_type) :: rv

                    rv%phi = cmplx(1,2, kind=8)
                    rv%rho = cmplx(3,4, kind=8)
                    rv%theta = cmplx(5,6, kind=8)
                end function

                function get_test_values_2() result(rv)
                    implicit none
                    ! dummy
                    type (spherical_coordinate_cmplx_type) :: rv

                    rv%phi = cmplx(6,5, kind=8)
                    rv%rho = cmplx(4,3, kind=8)
                    rv%theta = cmplx(2,1, kind=8)
                end function

                function evaluate(error, counter, str) result (rv)
                    implicit none
                    ! dummy
                    real(kind=8) :: error
                    integer(kind=4) :: counter
                    character(len=*), optional :: str

                    logical :: rv

                    ! auxiliary
                    double precision, parameter :: ground_truth_e = 10.0_8**(-6.0_8)

                    if (error .gt. ground_truth_e) then
                        if (present(str)) then
                            print *, "  ", str, " ", counter , "error: ", error, " : FAILED"
                        else
                            print *, "  ", counter, "error: ", error, " : FAILED"
                        end if
                        rv = .false.
                    else
                        if (present(str)) then
                            print *, "  ", str, " ", counter , "OK"
                        else
                            print *, "  ", counter, "OK"
                        end if
                        rv = .true.
                    end if
                end function evaluate

                function test_lib_math_list_spherical_operator_add_array() result (rv)
                    implicit none
                    ! dummy
                    logical :: rv

                    ! parameter
                    integer(kind=4), parameter :: d = 2

                    ! auxiliary
                    integer(kind=4) :: i
                    real(kind=8) :: buffer
                    type (list_spherical_coordinate_cmplx_type) :: lhs
                    type (list_spherical_coordinate_cmplx_type) :: rhs
                    type (list_spherical_coordinate_cmplx_type) :: value

                    type (list_spherical_coordinate_cmplx_type) :: ground_truth_value

                    allocate(lhs%coordinate(d))
                    allocate(rhs%coordinate(d))
                    allocate(ground_truth_value%coordinate(d))

                    ! lhs
                    lhs%coordinate(1) = get_test_values_1()
                    lhs%coordinate(2) = get_test_values_1()

                    ! rhs
                    rhs%coordinate(1) = get_test_values_1()
                    rhs%coordinate(2) = get_test_values_2()

                    ! ground truth
                    ground_truth_value%coordinate(1)%phi = cmplx(2,4, kind=8)
                    ground_truth_value%coordinate(1)%rho = cmplx(6,8, kind=8)
                    ground_truth_value%coordinate(1)%theta = cmplx(10,12, kind=8)

                    ground_truth_value%coordinate(2)%phi = cmplx(7,7, kind=8)
                    ground_truth_value%coordinate(2)%rho = cmplx(7,7, kind=8)
                    ground_truth_value%coordinate(2)%theta = cmplx(7,7, kind=8)

                    value = lib_math_list_spherical_operator_add_array(lhs, rhs)

                    rv = .true.
                    print *, "test_lib_math_list_spherical_operator_add_array:"
                    do i=1, d
                        buffer = abs(value%coordinate(i)%phi - ground_truth_value%coordinate(i)%phi)
                        rv = evaluate(buffer, i, 'phi')

                        buffer = abs(value%coordinate(i)%rho - ground_truth_value%coordinate(i)%rho)
                        rv = evaluate(buffer, i, "rho")

                        buffer = abs(value%coordinate(i)%theta - ground_truth_value%coordinate(i)%theta)
                        rv = evaluate(buffer, i, "theta")
                    end do

                end function test_lib_math_list_spherical_operator_add_array

                function test_lib_math_list_spherical_operator_sub_array() result (rv)
                    implicit none
                    ! dummy
                    logical :: rv

                    ! parameter
                    integer(kind=4), parameter :: d = 2

                    ! auxiliary
                    integer(kind=4) :: i
                    real(kind=8) :: buffer
                    type (list_spherical_coordinate_cmplx_type) :: lhs
                    type (list_spherical_coordinate_cmplx_type) :: rhs
                    type (list_spherical_coordinate_cmplx_type) :: value

                    type (list_spherical_coordinate_cmplx_type) :: ground_truth_value

                    allocate(lhs%coordinate(d))
                    allocate(rhs%coordinate(d))
                    allocate(ground_truth_value%coordinate(d))

                    ! lhs
                    lhs%coordinate(1) = get_test_values_1()
                    lhs%coordinate(2) = get_test_values_1()

                    ! rhs
                    rhs%coordinate(1) = get_test_values_1()
                    rhs%coordinate(2) = get_test_values_2()

                    ! ground truth
                    ground_truth_value%coordinate(1)%phi = cmplx(0,0, kind=8)
                    ground_truth_value%coordinate(1)%rho = cmplx(0,0, kind=8)
                    ground_truth_value%coordinate(1)%theta = cmplx(0,0, kind=8)

                    ground_truth_value%coordinate(2)%phi = cmplx(-5,-3, kind=8)
                    ground_truth_value%coordinate(2)%rho = cmplx(-1,1, kind=8)
                    ground_truth_value%coordinate(2)%theta = cmplx(3,5, kind=8)

                    value = lib_math_list_spherical_operator_sub_array(lhs, rhs)

                    rv = .true.
                    print *, "test_lib_math_list_spherical_operator_sub_array:"
                    do i=1, d
                        buffer = abs(value%coordinate(i)%phi - ground_truth_value%coordinate(i)%phi)
                        rv = evaluate(buffer, i, 'phi')

                        buffer = abs(value%coordinate(i)%rho - ground_truth_value%coordinate(i)%rho)
                        rv = evaluate(buffer, i, "rho")

                        buffer = abs(value%coordinate(i)%theta - ground_truth_value%coordinate(i)%theta)
                        rv = evaluate(buffer, i, "theta")
                    end do

                end function test_lib_math_list_spherical_operator_sub_array

                function test_lib_math_list_spherical_operator_0d_add_array() result (rv)
                    implicit none
                    ! dummy
                    logical :: rv

                    ! parameter
                    integer(kind=4), parameter :: d = 2

                    ! auxiliary
                    integer(kind=4) :: i
                    real(kind=8) :: buffer
                    type (spherical_coordinate_cmplx_type) :: lhs
                    type (list_spherical_coordinate_cmplx_type) :: rhs
                    type (list_spherical_coordinate_cmplx_type) :: value

                    type (list_spherical_coordinate_cmplx_type) :: ground_truth_value

                    allocate(rhs%coordinate(d))
                    allocate(ground_truth_value%coordinate(d))

                    ! lhs
                    lhs = get_test_values_1()

                    ! rhs
                    rhs%coordinate(1) = get_test_values_1()
                    rhs%coordinate(2) = get_test_values_2()

                    ! ground truth
                    ground_truth_value%coordinate(1)%phi = cmplx(2,4, kind=8)
                    ground_truth_value%coordinate(1)%rho = cmplx(6,8, kind=8)
                    ground_truth_value%coordinate(1)%theta = cmplx(10,12, kind=8)

                    ground_truth_value%coordinate(2)%phi = cmplx(7,7, kind=8)
                    ground_truth_value%coordinate(2)%rho = cmplx(7,7, kind=8)
                    ground_truth_value%coordinate(2)%theta = cmplx(7,7, kind=8)

                    value = lib_math_list_spherical_operator_0d_add_array(lhs, rhs)

                    rv = .true.
                    print *, "test_lib_math_list_spherical_operator_0d_add_array:"
                    do i=1, d
                        buffer = abs(value%coordinate(i)%phi - ground_truth_value%coordinate(i)%phi)
                        rv = evaluate(buffer, i, 'phi')

                        buffer = abs(value%coordinate(i)%rho - ground_truth_value%coordinate(i)%rho)
                        rv = evaluate(buffer, i, "rho")

                        buffer = abs(value%coordinate(i)%theta - ground_truth_value%coordinate(i)%theta)
                        rv = evaluate(buffer, i, "theta")
                    end do

                end function test_lib_math_list_spherical_operator_0d_add_array

                function test_lib_math_list_spherical_operator_array_add_0d() result (rv)
                    implicit none
                    ! dummy
                    logical :: rv

                    ! parameter
                    integer(kind=4), parameter :: d = 2

                    ! auxiliary
                    integer(kind=4) :: i
                    real(kind=8) :: buffer
                    type (list_spherical_coordinate_cmplx_type) :: lhs
                    type (spherical_coordinate_cmplx_type) :: rhs
                    type (list_spherical_coordinate_cmplx_type) :: value

                    type (list_spherical_coordinate_cmplx_type) :: ground_truth_value

                    allocate(lhs%coordinate(d))
                    allocate(ground_truth_value%coordinate(d))

                    ! lhs
                    lhs%coordinate(1) = get_test_values_1()
                    lhs%coordinate(2) = get_test_values_2()

                    ! rhs
                    rhs = get_test_values_1()

                    ! ground truth
                    ground_truth_value%coordinate(1)%phi = cmplx(2,4, kind=8)
                    ground_truth_value%coordinate(1)%rho = cmplx(6,8, kind=8)
                    ground_truth_value%coordinate(1)%theta = cmplx(10,12, kind=8)

                    ground_truth_value%coordinate(2)%phi = cmplx(7,7, kind=8)
                    ground_truth_value%coordinate(2)%rho = cmplx(7,7, kind=8)
                    ground_truth_value%coordinate(2)%theta = cmplx(7,7, kind=8)

                    value = lib_math_list_spherical_operator_array_add_0d(lhs, rhs)

                    rv = .true.
                    print *, "test_lib_math_list_spherical_operator_array_add_0d:"
                    do i=1, d
                        buffer = abs(value%coordinate(i)%phi - ground_truth_value%coordinate(i)%phi)
                        rv = evaluate(buffer, i, 'phi')

                        buffer = abs(value%coordinate(i)%rho - ground_truth_value%coordinate(i)%rho)
                        rv = evaluate(buffer, i, "rho")

                        buffer = abs(value%coordinate(i)%theta - ground_truth_value%coordinate(i)%theta)
                        rv = evaluate(buffer, i, "theta")
                    end do

                end function test_lib_math_list_spherical_operator_array_add_0d

                function test_lib_math_list_spherical_operator_array_sub_0d() result (rv)
                    implicit none
                    ! dummy
                    logical :: rv

                    ! parameter
                    integer(kind=4), parameter :: d = 2

                    ! auxiliary
                    integer(kind=4) :: i
                    real(kind=8) :: buffer
                    type (list_spherical_coordinate_cmplx_type) :: lhs
                    type (spherical_coordinate_cmplx_type) :: rhs
                    type (list_spherical_coordinate_cmplx_type) :: value

                    type (list_spherical_coordinate_cmplx_type) :: ground_truth_value

                    allocate(lhs%coordinate(d))
                    allocate(ground_truth_value%coordinate(d))

                    ! lhs
                    lhs%coordinate(1) = get_test_values_1()
                    lhs%coordinate(2) = get_test_values_2()

                    ! rhs
                    rhs = get_test_values_1()

                    ! ground truth
                    ground_truth_value%coordinate(1)%phi = cmplx(0,0, kind=8)
                    ground_truth_value%coordinate(1)%rho = cmplx(0,0, kind=8)
                    ground_truth_value%coordinate(1)%theta = cmplx(0,0, kind=8)

                    ground_truth_value%coordinate(2)%phi = cmplx(5,3, kind=8)
                    ground_truth_value%coordinate(2)%rho = cmplx(1,-1, kind=8)
                    ground_truth_value%coordinate(2)%theta = cmplx(-3,-5, kind=8)

                    value = lib_math_list_spherical_operator_array_sub_0d(lhs, rhs)

                    rv = .true.
                    print *, "test_lib_math_list_spherical_operator_array_sub_0d:"
                    do i=1, d
                        buffer = abs(value%coordinate(i)%phi - ground_truth_value%coordinate(i)%phi)
                        rv = evaluate(buffer, i, 'phi')

                        buffer = abs(value%coordinate(i)%rho - ground_truth_value%coordinate(i)%rho)
                        rv = evaluate(buffer, i, "rho")

                        buffer = abs(value%coordinate(i)%theta - ground_truth_value%coordinate(i)%theta)
                        rv = evaluate(buffer, i, "theta")
                    end do

                end function test_lib_math_list_spherical_operator_array_sub_0d

                function test_lib_math_list_spherical_operator_array_real_mul_array() result (rv)
                    implicit none
                    ! dummy
                    logical :: rv

                    ! parameter
                    integer(kind=4), parameter :: d = 2

                    ! auxiliary
                    integer(kind=4) :: i
                    real(kind=8) :: buffer
                    real(kind=lib_math_type_kind), dimension(d) :: lhs
                    type (list_spherical_coordinate_cmplx_type)  :: rhs
                    type (list_spherical_coordinate_cmplx_type) :: value

                    type (list_spherical_coordinate_cmplx_type) :: ground_truth_value

                    allocate(rhs%coordinate(d))
                    allocate(ground_truth_value%coordinate(d))

                    ! lhs
                    lhs(1) = 2.0
                    lhs(2) = 5.0

                    ! rhs
                    rhs%coordinate(1) = get_test_values_1()
                    rhs%coordinate(2) = get_test_values_2()

                    ! ground truth
                    ground_truth_value%coordinate(1)%phi = cmplx(2,4, kind=8)
                    ground_truth_value%coordinate(1)%rho = cmplx(6,8, kind=8)
                    ground_truth_value%coordinate(1)%theta = cmplx(10,12, kind=8)

                    ground_truth_value%coordinate(2)%phi = cmplx(30,25, kind=8)
                    ground_truth_value%coordinate(2)%rho = cmplx(20,15, kind=8)
                    ground_truth_value%coordinate(2)%theta = cmplx(10,5, kind=8)

                    value = lib_math_list_spherical_operator_array_real_mul_array(lhs, rhs)

                    rv = .true.
                    print *, "test_lib_math_list_spherical_operator_array_real_mul_array:"
                    do i=1, d
                        buffer = abs(value%coordinate(i)%phi - ground_truth_value%coordinate(i)%phi)
                        rv = evaluate(buffer, i, 'phi')

                        buffer = abs(value%coordinate(i)%rho - ground_truth_value%coordinate(i)%rho)
                        rv = evaluate(buffer, i, "rho")

                        buffer = abs(value%coordinate(i)%theta - ground_truth_value%coordinate(i)%theta)
                        rv = evaluate(buffer, i, "theta")
                    end do

                end function test_lib_math_list_spherical_operator_array_real_mul_array

                function test_lib_math_list_spherical_operator_array_cmplx_mul_array() result (rv)
                    implicit none
                    ! dummy
                    logical :: rv

                    ! parameter
                    integer(kind=4), parameter :: d = 2

                    ! auxiliary
                    integer(kind=4) :: i
                    real(kind=8) :: buffer
                    complex(kind=lib_math_type_kind), dimension(d) :: lhs
                    type (list_spherical_coordinate_cmplx_type) :: rhs
                    type (list_spherical_coordinate_cmplx_type) :: value

                    type (list_spherical_coordinate_cmplx_type) :: ground_truth_value

                    allocate(rhs%coordinate(d))
                    allocate(ground_truth_value%coordinate(d))

                    ! lhs
                    lhs(1) = cmplx(0, 1, kind=8)
                    lhs(2) = cmplx(2, 0, kind=8)

                    ! rhs
                    rhs%coordinate(1) = get_test_values_1()
                    rhs%coordinate(2) = get_test_values_2()

                    ! ground truth
                    ground_truth_value%coordinate(1)%phi = cmplx(-2,1, kind=8)
                    ground_truth_value%coordinate(1)%rho = cmplx(-4,3, kind=8)
                    ground_truth_value%coordinate(1)%theta = cmplx(-6,5, kind=8)

                    ground_truth_value%coordinate(2)%phi = cmplx(12,10, kind=8)
                    ground_truth_value%coordinate(2)%rho = cmplx(8,6, kind=8)
                    ground_truth_value%coordinate(2)%theta = cmplx(4,2, kind=8)

                    value = lib_math_list_spherical_operator_array_cmplx_mul_array(lhs, rhs)

                    rv = .true.
                    print *, "test_lib_math_list_spherical_operator_array_cmplx_mul_array:"
                    do i=1, d
                        buffer = abs(value%coordinate(i)%phi - ground_truth_value%coordinate(i)%phi)
                        rv = evaluate(buffer, i, 'phi')

                        buffer = abs(value%coordinate(i)%rho - ground_truth_value%coordinate(i)%rho)
                        rv = evaluate(buffer, i, "rho")

                        buffer = abs(value%coordinate(i)%theta - ground_truth_value%coordinate(i)%theta)
                        rv = evaluate(buffer, i, "theta")
                    end do

                end function test_lib_math_list_spherical_operator_array_cmplx_mul_array

                function test_lib_math_list_spherical_operator_real_mul_array() result (rv)
                    implicit none
                    ! dummy
                    logical :: rv

                    ! parameter
                    integer(kind=4), parameter :: d = 2

                    ! auxiliary
                    integer(kind=4) :: i
                    real(kind=8) :: buffer
                    real(kind=lib_math_type_kind) :: lhs
                    type (list_spherical_coordinate_cmplx_type) :: rhs
                    type (list_spherical_coordinate_cmplx_type) :: value

                    type (list_spherical_coordinate_cmplx_type) :: ground_truth_value

                    allocate(rhs%coordinate(d))
                    allocate(ground_truth_value%coordinate(d))

                    ! lhs
                    lhs = 4

                    ! rhs
                    rhs%coordinate(1) = get_test_values_1()
                    rhs%coordinate(2) = get_test_values_2()

                    ! ground truth
                    ground_truth_value%coordinate(1)%phi = cmplx(4,8, kind=8)
                    ground_truth_value%coordinate(1)%rho = cmplx(12,16, kind=8)
                    ground_truth_value%coordinate(1)%theta = cmplx(20,24, kind=8)

                    ground_truth_value%coordinate(2)%phi = cmplx(24,20, kind=8)
                    ground_truth_value%coordinate(2)%rho = cmplx(16,12, kind=8)
                    ground_truth_value%coordinate(2)%theta = cmplx(8,4, kind=8)

                    value = lib_math_list_spherical_operator_real_mul_array(lhs, rhs)

                    rv = .true.
                    print *, "test_lib_math_list_spherical_operator_real_mul_array:"
                    do i=1, d
                        buffer = abs(value%coordinate(i)%phi - ground_truth_value%coordinate(i)%phi)
                        rv = evaluate(buffer, i, 'phi')

                        buffer = abs(value%coordinate(i)%rho - ground_truth_value%coordinate(i)%rho)
                        rv = evaluate(buffer, i, "rho")

                        buffer = abs(value%coordinate(i)%theta - ground_truth_value%coordinate(i)%theta)
                        rv = evaluate(buffer, i, "theta")
                    end do

                end function test_lib_math_list_spherical_operator_real_mul_array

                function test_lib_math_list_spherical_operator_cmplx_mul_array() result (rv)
                    implicit none
                    ! dummy
                    logical :: rv

                    ! parameter
                    integer(kind=4), parameter :: d = 2

                    ! auxiliary
                    integer(kind=4) :: i
                    real(kind=8) :: buffer
                    complex(kind=lib_math_type_kind) :: lhs
                    type (list_spherical_coordinate_cmplx_type) :: rhs
                    type (list_spherical_coordinate_cmplx_type) :: value

                    type (list_spherical_coordinate_cmplx_type) :: ground_truth_value

                    allocate(rhs%coordinate(d))
                    allocate(ground_truth_value%coordinate(d))

                    ! lhs
                    lhs = cmplx(0, 4, kind=8)

                    ! rhs
                    rhs%coordinate(1) = get_test_values_1()
                    rhs%coordinate(2) = get_test_values_2()

                    ! ground truth
                    ground_truth_value%coordinate(1)%phi = cmplx(-8,4, kind=8)
                    ground_truth_value%coordinate(1)%rho = cmplx(-16,12, kind=8)
                    ground_truth_value%coordinate(1)%theta = cmplx(-24,20, kind=8)

                    ground_truth_value%coordinate(2)%phi = cmplx(-20,24, kind=8)
                    ground_truth_value%coordinate(2)%rho = cmplx(-12,16, kind=8)
                    ground_truth_value%coordinate(2)%theta = cmplx(-4,8, kind=8)

                    value = lib_math_list_spherical_operator_cmplx_mul_array(lhs, rhs)

                    rv = .true.
                    print *, "test_lib_math_list_spherical_operator_cmplx_mul_array:"
                    do i=1, d
                        buffer = abs(value%coordinate(i)%phi - ground_truth_value%coordinate(i)%phi)
                        rv = evaluate(buffer, i, 'phi')

                        buffer = abs(value%coordinate(i)%rho - ground_truth_value%coordinate(i)%rho)
                        rv = evaluate(buffer, i, "rho")

                        buffer = abs(value%coordinate(i)%theta - ground_truth_value%coordinate(i)%theta)
                        rv = evaluate(buffer, i, "theta")
                    end do

                end function test_lib_math_list_spherical_operator_cmplx_mul_array

                function test_lib_math_list_spherical_operator_array_divide_by_r_array() result (rv)
                    implicit none
                    ! dummy
                    logical :: rv

                    ! parameter
                    integer(kind=4), parameter :: d = 2

                    ! auxiliary
                    integer(kind=4) :: i
                    real(kind=8) :: buffer
                    type (list_spherical_coordinate_cmplx_type) :: lhs
                    real(kind=lib_math_type_kind), dimension(d) :: rhs
                    type (list_spherical_coordinate_cmplx_type) :: value

                    type (list_spherical_coordinate_cmplx_type) :: ground_truth_value

                    allocate(lhs%coordinate(d))
                    allocate(ground_truth_value%coordinate(d))

                    ! lhs
                    lhs%coordinate(1) = get_test_values_1()
                    lhs%coordinate(2) = get_test_values_2()

                    ! rhs
                    rhs(1) = 2.0
                    rhs(2) = 4.0

                    ! ground truth
                    ground_truth_value%coordinate(1)%phi = cmplx(0.5,1, kind=8)
                    ground_truth_value%coordinate(1)%rho = cmplx(1.5,2, kind=8)
                    ground_truth_value%coordinate(1)%theta = cmplx(2.5,3, kind=8)

                    ground_truth_value%coordinate(2)%phi = cmplx(1.5,1.25, kind=8)
                    ground_truth_value%coordinate(2)%rho = cmplx(1,0.625, kind=8)
                    ground_truth_value%coordinate(2)%theta = cmplx(0.5,0.25, kind=8)

                    value = lib_math_list_spherical_operator_array_divide_by_real_array(lhs, rhs)

                    rv = .true.
                    print *, "test_lib_math_list_spherical_operator_array_divide_by_r_array:"
                    do i=1, d
                        buffer = abs(value%coordinate(i)%phi - ground_truth_value%coordinate(i)%phi)
                        rv = evaluate(buffer, i, 'phi')

                        buffer = abs(value%coordinate(i)%rho - ground_truth_value%coordinate(i)%rho)
                        rv = evaluate(buffer, i, "rho")

                        buffer = abs(value%coordinate(i)%theta - ground_truth_value%coordinate(i)%theta)
                        rv = evaluate(buffer, i, "theta")
                    end do

                end function test_lib_math_list_spherical_operator_array_divide_by_r_array

                function test_lib_math_list_spherical_operator_array_divide_by_c_array() result (rv)
                    implicit none
                    ! dummy
                    logical :: rv

                    ! parameter
                    integer(kind=4), parameter :: d = 2

                    ! auxiliary
                    integer(kind=4) :: i
                    real(kind=8) :: buffer
                    type (list_spherical_coordinate_cmplx_type) :: lhs
                    complex(kind=lib_math_type_kind), dimension(d) :: rhs
                    type (list_spherical_coordinate_cmplx_type) :: value

                    type (list_spherical_coordinate_cmplx_type) :: ground_truth_value

                    allocate(lhs%coordinate(d))
                    allocate(ground_truth_value%coordinate(d))

                    ! lhs
                    lhs%coordinate(1) = get_test_values_1()
                    lhs%coordinate(2) = get_test_values_2()

                    ! rhs
                    rhs(1) = cmplx(0, 2, kind=8)
                    rhs(2) = cmplx(0, 4, kind=8)

                    ! ground truth
                    ground_truth_value%coordinate(1)%phi = cmplx(1, -0.5, kind=8)
                    ground_truth_value%coordinate(1)%rho = cmplx(2, -1.5, kind=8)
                    ground_truth_value%coordinate(1)%theta = cmplx(3, -2.5, kind=8)

                    ground_truth_value%coordinate(2)%phi = cmplx(1.25, -1.5, kind=8)
                    ground_truth_value%coordinate(2)%rho = cmplx(0.625, -1, kind=8)
                    ground_truth_value%coordinate(2)%theta = cmplx(0.25, -0.5, kind=8)

                    value = lib_math_list_spherical_operator_array_divide_by_cmplx_array(lhs, rhs)

                    rv = .true.
                    print *, "test_lib_math_list_spherical_operator_array_divide_by_c_array:"
                    do i=1, d
                        buffer = abs(value%coordinate(i)%phi - ground_truth_value%coordinate(i)%phi)
                        rv = evaluate(buffer, i, 'phi')

                        buffer = abs(value%coordinate(i)%rho - ground_truth_value%coordinate(i)%rho)
                        rv = evaluate(buffer, i, "rho")

                        buffer = abs(value%coordinate(i)%theta - ground_truth_value%coordinate(i)%theta)
                        rv = evaluate(buffer, i, "theta")
                    end do

                end function test_lib_math_list_spherical_operator_array_divide_by_c_array

                function test_lib_math_list_spherical_operator_array_divide_by_real() result (rv)
                    implicit none
                    ! dummy
                    logical :: rv

                    ! parameter
                    integer(kind=4), parameter :: d = 2

                    ! auxiliary
                    integer(kind=4) :: i
                    real(kind=8) :: buffer
                    type (list_spherical_coordinate_cmplx_type) :: lhs
                    real(kind=lib_math_type_kind) :: rhs
                    type (list_spherical_coordinate_cmplx_type) :: value

                    type (list_spherical_coordinate_cmplx_type) :: ground_truth_value

                    allocate(lhs%coordinate(d))
                    allocate(ground_truth_value%coordinate(d))

                    ! lhs
                    lhs%coordinate(1) = get_test_values_1()
                    lhs%coordinate(2) = get_test_values_2()

                    ! rhs
                    rhs = 2

                    ! ground truth
                    ground_truth_value%coordinate(1)%phi = cmplx(0.5, 1, kind=8)
                    ground_truth_value%coordinate(1)%rho = cmplx(1.5, 2, kind=8)
                    ground_truth_value%coordinate(1)%theta = cmplx(2.5, 3, kind=8)

                    ground_truth_value%coordinate(2)%phi = cmplx(3, 2.5, kind=8)
                    ground_truth_value%coordinate(2)%rho = cmplx(2, 1.5, kind=8)
                    ground_truth_value%coordinate(2)%theta = cmplx(1, 0.5, kind=8)

                    value = lib_math_list_spherical_operator_array_divide_by_real(lhs, rhs)

                    rv = .true.
                    print *, "test_lib_math_list_spherical_operator_array_divide_by_real:"
                    do i=1, d
                        buffer = abs(value%coordinate(i)%phi - ground_truth_value%coordinate(i)%phi)
                        rv = evaluate(buffer, i, 'phi')

                        buffer = abs(value%coordinate(i)%rho - ground_truth_value%coordinate(i)%rho)
                        rv = evaluate(buffer, i, "rho")

                        buffer = abs(value%coordinate(i)%theta - ground_truth_value%coordinate(i)%theta)
                        rv = evaluate(buffer, i, "theta")
                    end do

                end function test_lib_math_list_spherical_operator_array_divide_by_real

                function test_lib_math_list_spherical_operator_array_divide_by_cmplx() result (rv)
                    implicit none
                    ! dummy
                    logical :: rv

                    ! parameter
                    integer(kind=4), parameter :: d = 2

                    ! auxiliary
                    integer(kind=4) :: i
                    real(kind=8) :: buffer
                    type (list_spherical_coordinate_cmplx_type) :: lhs
                    complex(kind=lib_math_type_kind) :: rhs
                    type (list_spherical_coordinate_cmplx_type) :: value

                    type (list_spherical_coordinate_cmplx_type) :: ground_truth_value

                    allocate(lhs%coordinate(d))
                    allocate(ground_truth_value%coordinate(d))

                    ! lhs
                    lhs%coordinate(1) = get_test_values_1()
                    lhs%coordinate(2) = get_test_values_2()

                    ! rhs
                    rhs = cmplx(0, 2, kind=8)

                    ! ground truth
                    ground_truth_value%coordinate(1)%phi = cmplx(1, -0.5, kind=8)
                    ground_truth_value%coordinate(1)%rho = cmplx(2, -1.5, kind=8)
                    ground_truth_value%coordinate(1)%theta = cmplx(3, -2.5, kind=8)

                    ground_truth_value%coordinate(2)%phi = cmplx(2.5, -3, kind=8)
                    ground_truth_value%coordinate(2)%rho = cmplx(1.25, -2, kind=8)
                    ground_truth_value%coordinate(2)%theta = cmplx(0.5, -1, kind=8)

                    value = lib_math_list_spherical_operator_array_divide_by_cmplx(lhs, rhs)

                    rv = .true.
                    print *, "test_lib_math_list_spherical_operator_array_divide_by_cmplx:"
                    do i=1, d
                        buffer = abs(value%coordinate(i)%phi - ground_truth_value%coordinate(i)%phi)
                        rv = evaluate(buffer, i, 'phi')

                        buffer = abs(value%coordinate(i)%rho - ground_truth_value%coordinate(i)%rho)
                        rv = evaluate(buffer, i, "rho")

                        buffer = abs(value%coordinate(i)%theta - ground_truth_value%coordinate(i)%theta)
                        rv = evaluate(buffer, i, "theta")
                    end do

                end function test_lib_math_list_spherical_operator_array_divide_by_cmplx

        end function lib_math_type_operator_test_functions

end module lib_math_type_operator
