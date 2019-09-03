module lib_math_type_operator
    use lib_math_type
    implicit none

    private

    ! --- public ---
    public :: operator (+)
    public :: operator (-)
    public :: operator (*)
    public :: operator (/)
    public :: assignment (=)
    public :: abs
    public :: spherical_abs
    public :: cartesian_abs
    public :: init_list

    public :: make_cartesian
    public :: make_spherical

    public :: cross_product
    public :: dot_product

    public :: lib_math_type_operator_test_functions

    ! ---- operator ----
    interface operator (+)
        module procedure lib_math_cartesian_operator_add
        module procedure lib_math_cartesian_operator_add_array
        module procedure lib_math_cartesian_operator_0d_add_array
        module procedure lib_math_cartesian_operator_array_add_0d
        module procedure lib_math_list_cartesian_operator_add_array_cmplx
        module procedure lib_math_list_cartesian_operator_0d_add_array_cmplx
        module procedure lib_math_list_cartesian_operator_array_add_0d_cmplx

        module procedure lib_math_spherical_operator_add
        module procedure lib_math_spherical_operator_add_array
        module procedure lib_math_spherical_operator_0d_add_array
        module procedure lib_math_spherical_operator_array_add_0d
        module procedure lib_math_list_spherical_operator_add_array_cmplx
        module procedure lib_math_list_spherical_operator_0d_add_array_cmplx
        module procedure lib_math_list_spherical_operator_array_add_0d_cmplx
    end interface

    interface operator (-)
        module procedure lib_math_cartesian_operator_sub
        module procedure lib_math_cartesian_operator_sub_array
        module procedure lib_math_cartesian_operator_array_sub_0d
        module procedure lib_math_list_cartesian_operator_sub_array_cmplx
        module procedure lib_math_list_cartesian_operator_array_sub_0d_cmplx

        module procedure lib_math_spherical_operator_sub
        module procedure lib_math_spherical_operator_sub_array
        module procedure lib_math_spherical_operator_array_sub_0d
        module procedure lib_math_list_spherical_operator_sub_array_cmplx
        module procedure lib_math_list_spherical_operator_array_sub_0d_cmplx
    end interface

    interface operator (*)
        module procedure lib_math_cartesian_operator_scalar_real_mul
        module procedure lib_math_cartesian_operator_scalar_cmplx_mul
        module procedure lib_math_cartesian_operator_mul_scalar_real
        module procedure lib_math_cartesian_operator_mul_scalar_cmplx
        module procedure lib_math_cartesian_operator_array_real_mul_array
        module procedure lib_math_cartesian_operator_array_cmplx_mul_array
        module procedure lib_math_cartesian_operator_0d_real_mul_array
        module procedure lib_math_cartesian_operator_0d_cmplx_mul_array
        module procedure lib_math_list_cartesian_operator_array_real_mul_array_c
        module procedure lib_math_list_cartesian_operator_array_c_mul_array_cmplx
        module procedure lib_math_list_cartesian_operator_real_mul_array_cmplx
        module procedure lib_math_list_cartesian_operator_cmplx_mul_array_cmplx

        module procedure lib_math_spherical_operator_scalar_real_mul
        module procedure lib_math_spherical_operator_scalar_cmplx_mul
        module procedure lib_math_spherical_operator_mul_scalar_real
        module procedure lib_math_spherical_operator_mul_scalar_cmplx
        module procedure lib_math_spherical_operator_array_real_mul_array
        module procedure lib_math_spherical_operator_array_cmplx_mul_array
        module procedure lib_math_spherical_operator_0d_real_mul_array
        module procedure lib_math_spherical_operator_0d_cmplx_mul_array
        module procedure lib_math_list_spherical_operator_array_real_mul_array_c
        module procedure lib_math_list_spherical_operator_array_c_mul_array_cmplx
        module procedure lib_math_list_spherical_operator_real_mul_array_cmplx
        module procedure lib_math_list_spherical_operator_cmplx_mul_array_cmplx
    end interface

    interface operator (/)
        module procedure lib_math_cartesian_operator_divide_by_scalar_real
        module procedure lib_math_cartesian_operator_divide_by_scalar_cmplx
        module procedure lib_math_cartesian_operator_array_divide_by_scalar_array_real
        module procedure lib_math_cartesian_operator_array_divide_by_scalar_array_cmplx
        module procedure lib_math_cartesian_operator_array_divide_by_0d_real
        module procedure lib_math_cartesian_operator_array_divide_by_0d_cmplx
        module procedure lib_math_list_cartesian_operator_array_divide_by_real_array
        module procedure lib_math_list_cartesian_operator_c_array_divide_by_c_array
        module procedure lib_math_list_cartesian_operator_array_divide_by_real
        module procedure lib_math_list_cartesian_operator_array_divide_by_cmplx

        module procedure lib_math_spherical_operator_divide_by_scalar_real
        module procedure lib_math_spherical_operator_divide_by_scalar_cmplx
        module procedure lib_math_spherical_operator_array_divide_by_scalar_array_real
        module procedure lib_math_spherical_operator_array_divide_by_scalar_array_cmplx
        module procedure lib_math_spherical_operator_array_divide_by_0d_real
        module procedure lib_math_spherical_operator_array_divide_by_0d_cmplx
        module procedure lib_math_list_spherical_operator_array_divide_by_real_array
        module procedure lib_math_list_spherical_operator_c_array_divide_by_c_array
        module procedure lib_math_list_spherical_operator_array_divide_by_real
        module procedure lib_math_list_spherical_operator_array_divide_by_cmplx
    end interface

    interface assignment (=)
        module procedure lib_math_spherical_point_to_cartesian_point
        module procedure lib_math_cartesian_point_to_spherical_point
    end interface

    interface spherical_abs
        module procedure lib_math_spherical_operator_abs_cmplx
    end interface

    interface cartesian_abs
        module procedure lib_math_cartesian_operator_abs_cmplx
        module procedure lib_math_cartesian_operator_abs_real
    end interface

    interface abs
        module procedure lib_math_cartesian_operator_abs_cmplx
        module procedure lib_math_cartesian_operator_abs_real

        module procedure lib_math_spherical_operator_abs_cmplx
    end interface


    interface init_list
        module procedure lib_math_list_list_logical_init
        module procedure lib_math_list_list_integer_sys_init
        module procedure lib_math_list_list_integer_init
        module procedure lib_math_list_list_real_init
        module procedure lib_math_list_list_cmplx_init
        module procedure lib_math_list_spherical_coordinate_cmplx_type_init
    end interface

    interface make_cartesian
        module procedure lib_math_spherical_components_to_cartesian_components_cmplx_a
        module procedure lib_math_spherical_components_to_cartesian_components_cmplx_c
        module procedure lib_math_spherical_components_to_cartesian_components_cmplx_s
    end interface

    interface make_spherical
        module procedure lib_math_cartesian_components_to_spherical_components_cmplx_a
        module procedure lib_math_cartesian_components_to_spherical_components_cmplx_c
        module procedure lib_math_cartesian_components_to_spherical_components_cmplx_s
    end interface

    interface cross_product
        module procedure lib_math_cartesian_cross_product_real
        module procedure lib_math_cartesian_cross_product_real_list
        module procedure lib_math_cartesian_cross_product_cmplx
        module procedure lib_math_cartesian_cross_product_cmplx_list
    end interface

    interface dot_product
        module procedure lib_math_cartesian_dot_product_real
        module procedure lib_math_cartesian_dot_product_real_list
        module procedure lib_math_cartesian_dot_product_cmplx
        module procedure lib_math_cartesian_dot_product_cmplx_list
    end interface

    contains

! ---- single spherical coordinate ----
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

! ---- array spherical coordinate ----
        function lib_math_spherical_operator_add_array(lhs, rhs) result(rv)
            implicit none
            ! dummy
            type (spherical_coordinate_cmplx_type), dimension(:), intent(in) :: lhs
            type (spherical_coordinate_cmplx_type), dimension(lbound(lhs,1):ubound(lhs,1)), intent(in) :: rhs
            type (spherical_coordinate_cmplx_type), dimension(lbound(lhs,1):ubound(lhs,1)) :: rv

            ! auxiliary
            integer :: i

            !$OMP PARALLEL DO PRIVATE(i)
            do i=lbound(lhs,1), ubound(lhs,1)
                rv(i) = lhs(i) + rhs(i)
            end do
            !$OMP END PARALLEL DO

        end function

        function lib_math_spherical_operator_sub_array(lhs, rhs) result(rv)
            implicit none
            ! dummy
            type (spherical_coordinate_cmplx_type), dimension(:), intent(in) :: lhs
            type (spherical_coordinate_cmplx_type), dimension(lbound(lhs,1):ubound(lhs,1)), intent(in) :: rhs
            type (spherical_coordinate_cmplx_type), dimension(lbound(lhs,1):ubound(lhs,1)) :: rv

            ! auxiliary
            integer :: i

            !$OMP PARALLEL DO PRIVATE(i)
            do i=lbound(lhs,1), ubound(lhs,1)
                rv(i) = lhs(i) - rhs(i)
            end do
            !$OMP END PARALLEL DO

        end function

        function lib_math_spherical_operator_0d_add_array(lhs, rhs) result(rv)
            implicit none
            ! dummy
            type (spherical_coordinate_cmplx_type), intent(in) :: lhs
            type (spherical_coordinate_cmplx_type), dimension(:), intent(in) :: rhs
            type (spherical_coordinate_cmplx_type), dimension(lbound(rhs,1):ubound(rhs,1)) :: rv

            ! auxiliary
            integer :: i

            !$OMP PARALLEL DO PRIVATE(i)
            do i=lbound(rhs,1), ubound(rhs,1)
                rv(i) = lhs + rhs(i)
            end do
            !$OMP END PARALLEL DO

        end function

        function lib_math_spherical_operator_array_add_0d(lhs, rhs) result(rv)
            implicit none
            ! dummy
            type (spherical_coordinate_cmplx_type), dimension(:), intent(in) :: lhs
            type (spherical_coordinate_cmplx_type), intent(in) :: rhs
            type (spherical_coordinate_cmplx_type), dimension(lbound(lhs,1):ubound(lhs,1)) :: rv

            ! auxiliary
            integer :: i

            !$OMP PARALLEL DO PRIVATE(i)
            do i=lbound(lhs,1), ubound(lhs,1)
                rv(i) = lhs(i) + rhs
            end do
            !$OMP END PARALLEL DO

        end function

        function lib_math_spherical_operator_array_sub_0d(lhs, rhs) result(rv)
            implicit none
            ! dummy
            type (spherical_coordinate_cmplx_type), dimension(:), intent(in) :: lhs
            type (spherical_coordinate_cmplx_type), intent(in) :: rhs
            type (spherical_coordinate_cmplx_type), dimension(lbound(lhs,1):ubound(lhs,1)) :: rv

            ! auxiliary
            integer :: i

            !$OMP PARALLEL DO PRIVATE(i)
            do i=lbound(lhs,1), ubound(lhs,1)
                rv(i) = lhs(i) - rhs
            end do
            !$OMP END PARALLEL DO

        end function

        function lib_math_spherical_operator_array_real_mul_array(lhs, rhs) result(rv)
            implicit none
            ! dummy
            real(kind=lib_math_type_kind), dimension(:), intent(in) :: lhs
            type (spherical_coordinate_cmplx_type), dimension(lbound(lhs,1):ubound(lhs,1)), intent(in) :: rhs
            type (spherical_coordinate_cmplx_type), dimension(lbound(lhs,1):ubound(lhs,1)) :: rv

            ! auxiliary
            integer :: i

            !$OMP PARALLEL DO PRIVATE(i)
            do i=lbound(lhs,1), ubound(lhs,1)
                rv(i) = lhs(i) * rhs(i)
            end do
            !$OMP END PARALLEL DO

        end function

        function lib_math_spherical_operator_array_cmplx_mul_array(lhs, rhs) result(rv)
            implicit none
            ! dummy
            complex(kind=lib_math_type_kind), dimension(:), intent(in) :: lhs
            type (spherical_coordinate_cmplx_type), dimension(lbound(lhs,1):ubound(lhs,1)), intent(in) :: rhs
            type (spherical_coordinate_cmplx_type), dimension(lbound(lhs,1):ubound(lhs,1)) :: rv

            ! auxiliary
            integer :: i

            !$OMP PARALLEL DO PRIVATE(i)
            do i=lbound(lhs,1), ubound(lhs,1)
                rv(i) = lhs(i) * rhs(i)
            end do
            !$OMP END PARALLEL DO

        end function

        function lib_math_spherical_operator_0d_real_mul_array(lhs, rhs) result(rv)
            implicit none
            ! dummy
            real(kind=lib_math_type_kind), intent(in) :: lhs
            type (spherical_coordinate_cmplx_type), dimension(:), intent(in) :: rhs
            type (spherical_coordinate_cmplx_type), dimension(lbound(rhs,1):ubound(rhs,1)) :: rv

            ! auxiliary
            integer :: i

            !$OMP PARALLEL DO PRIVATE(i)
            do i=lbound(rhs,1), ubound(rhs,1)
                rv(i) = lhs * rhs(i)
            end do
            !$OMP END PARALLEL DO

        end function

        function lib_math_spherical_operator_0d_cmplx_mul_array(lhs, rhs) result(rv)
            implicit none
            ! dummy
            complex(kind=lib_math_type_kind), intent(in) :: lhs
            type (spherical_coordinate_cmplx_type), dimension(:), intent(in) :: rhs
            type (spherical_coordinate_cmplx_type), dimension(lbound(rhs,1):ubound(rhs,1)) :: rv

            ! auxiliary
            integer :: i

            !$OMP PARALLEL DO PRIVATE(i)
            do i=lbound(rhs,1), ubound(rhs,1)
                rv(i) = lhs * rhs(i)
            end do
            !$OMP END PARALLEL DO

        end function

        function lib_math_spherical_operator_array_divide_by_scalar_array_real(lhs, rhs) result(rv)
            implicit none
            ! dummy
            type (spherical_coordinate_cmplx_type), dimension(:), intent(in) :: lhs
            real(kind=lib_math_type_kind), dimension(lbound(lhs,1):ubound(lhs,1)), intent(in) :: rhs
            type (spherical_coordinate_cmplx_type), dimension(lbound(lhs,1):ubound(lhs,1)) :: rv

            ! auxiliary
            integer :: i

            !$OMP PARALLEL DO PRIVATE(i)
            do i=lbound(lhs,1), ubound(lhs,1)
                rv(i) = lhs(i) / rhs(i)
            end do
            !$OMP END PARALLEL DO

        end function

        function lib_math_spherical_operator_array_divide_by_scalar_array_cmplx(lhs, rhs) result(rv)
            implicit none
            ! dummy
            type (spherical_coordinate_cmplx_type), dimension(:), intent(in) :: lhs
            complex(kind=lib_math_type_kind), dimension(lbound(lhs,1):ubound(lhs,1)), intent(in) :: rhs
            type (spherical_coordinate_cmplx_type), dimension(lbound(lhs,1):ubound(lhs,1)) :: rv

            ! auxiliary
            integer :: i

            !$OMP PARALLEL DO PRIVATE(i)
            do i=lbound(lhs,1), ubound(lhs,1)
                rv(i) = lhs(i) / rhs(i)
            end do
            !$OMP END PARALLEL DO

        end function

        function lib_math_spherical_operator_array_divide_by_0d_real(lhs, rhs) result(rv)
            implicit none
            ! dummy
            type (spherical_coordinate_cmplx_type), dimension(:), intent(in) :: lhs
            real(kind=lib_math_type_kind), intent(in) :: rhs
            type (spherical_coordinate_cmplx_type), dimension(lbound(lhs,1):ubound(lhs,1)) :: rv

            ! auxiliary
            integer :: i

            !$OMP PARALLEL DO PRIVATE(i)
            do i=lbound(lhs,1), ubound(lhs,1)
                rv(i) = lhs(i) / rhs
            end do
            !$OMP END PARALLEL DO

        end function

        function lib_math_spherical_operator_array_divide_by_0d_cmplx(lhs, rhs) result(rv)
            implicit none
            ! dummy
            type (spherical_coordinate_cmplx_type), dimension(:), intent(in) :: lhs
            complex(kind=lib_math_type_kind), intent(in) :: rhs
            type (spherical_coordinate_cmplx_type), dimension(lbound(lhs,1):ubound(lhs,1)) :: rv

            ! auxiliary
            integer :: i

            !$OMP PARALLEL DO PRIVATE(i)
            do i=lbound(lhs,1), ubound(lhs,1)
                rv(i) = lhs(i) / rhs
            end do
            !$OMP END PARALLEL DO

        end function

! ---- list_spherical_coordinate ----
        function lib_math_list_spherical_operator_add_array_cmplx(lhs, rhs) result(rv)
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

                !$OMP PARALLEL DO PRIVATE(i)
                do i=lbound(lhs%coordinate,1), ubound(lhs%coordinate,1)
                    rv%coordinate(i) = lhs%coordinate(i) + rhs%coordinate(i)
                end do
                !$OMP END PARALLEL DO
            end if

        end function

        function lib_math_list_spherical_operator_sub_array_cmplx(lhs, rhs) result(rv)
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

                !$OMP PARALLEL DO PRIVATE(i)
                do i=lbound(lhs%coordinate,1), ubound(lhs%coordinate,1)
                    rv%coordinate(i) = lhs%coordinate(i) - rhs%coordinate(i)
                end do
                !$OMP END PARALLEL DO
            end if
        end function

        function lib_math_list_spherical_operator_0d_add_array_cmplx(lhs, rhs) result(rv)
            implicit none
            ! dummy
            type (spherical_coordinate_cmplx_type), intent(in) :: lhs
            type (list_spherical_coordinate_cmplx_type), intent(in) :: rhs
            type (list_spherical_coordinate_cmplx_type) :: rv

            ! auxiliary
            integer :: i

            allocate(rv%coordinate(lbound(rhs%coordinate,1):ubound(rhs%coordinate,1)))

            !$OMP PARALLEL DO PRIVATE(i)
            do i=lbound(rhs%coordinate,1), ubound(rhs%coordinate,1)
                rv%coordinate(i) = lhs + rhs%coordinate(i)
            end do
            !$OMP END PARALLEL DO

        end function

        function lib_math_list_spherical_operator_array_add_0d_cmplx(lhs, rhs) result(rv)
            implicit none
            ! dummy
            type (list_spherical_coordinate_cmplx_type), intent(in) :: lhs
            type (spherical_coordinate_cmplx_type), intent(in) :: rhs
            type (list_spherical_coordinate_cmplx_type) :: rv

            ! auxiliary
            integer :: i

            allocate(rv%coordinate(lbound(lhs%coordinate,1):ubound(lhs%coordinate,1)))

            !$OMP PARALLEL DO PRIVATE(i)
            do i=lbound(lhs%coordinate,1), ubound(lhs%coordinate,1)
                rv%coordinate(i) = lhs%coordinate(i) + rhs
            end do
            !$OMP END PARALLEL DO

        end function

        function lib_math_list_spherical_operator_array_sub_0d_cmplx(lhs, rhs) result(rv)
            implicit none
            ! dummy
            type (list_spherical_coordinate_cmplx_type), intent(in) :: lhs
            type (spherical_coordinate_cmplx_type), intent(in) :: rhs
            type (list_spherical_coordinate_cmplx_type) :: rv

            ! auxiliary
            integer :: i

            allocate(rv%coordinate(lbound(lhs%coordinate,1):ubound(lhs%coordinate,1)))

            !$OMP PARALLEL DO PRIVATE(i)
            do i=lbound(lhs%coordinate,1), ubound(lhs%coordinate,1)
                rv%coordinate(i) = lhs%coordinate(i) - rhs
            end do
            !$OMP END PARALLEL DO

        end function

        function lib_math_list_spherical_operator_array_real_mul_array_c(lhs, rhs) result(rv)
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

                !$OMP PARALLEL DO PRIVATE(i)
                do i=lbound(lhs,1), ubound(lhs,1)
                    rv%coordinate(i) = lhs(i) * rhs%coordinate(i)
                end do
                !$OMP END PARALLEL DO
            else
                print *, "lib_math_list_spherical_operator_array_real_mul_array_c: ERROR"
            end if

        end function

        function lib_math_list_spherical_operator_array_c_mul_array_cmplx(lhs, rhs) result(rv)
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

                !$OMP PARALLEL DO PRIVATE(i)
                do i=lbound(lhs,1), ubound(lhs,1)
                    rv%coordinate(i) = lhs(i) * rhs%coordinate(i)
                end do
                !$OMP END PARALLEL DO

            end if

        end function

        function lib_math_list_spherical_operator_real_mul_array_cmplx(lhs, rhs) result(rv)
            implicit none
            ! dummy
            real(kind=lib_math_type_kind), intent(in) :: lhs
            type (list_spherical_coordinate_cmplx_type), intent(in) :: rhs
            type (list_spherical_coordinate_cmplx_type) :: rv

            ! auxiliary
            integer :: i

            allocate(rv%coordinate(lbound(rhs%coordinate,1):ubound(rhs%coordinate,1)))

            !$OMP PARALLEL DO PRIVATE(i)
            do i=lbound(rhs%coordinate,1), ubound(rhs%coordinate,1)
                rv%coordinate(i) = lhs * rhs%coordinate(i)
            end do
            !$OMP END PARALLEL DO

        end function

        function lib_math_list_spherical_operator_cmplx_mul_array_cmplx(lhs, rhs) result(rv)
            implicit none
            ! dummy
            complex(kind=lib_math_type_kind), intent(in) :: lhs
            type (list_spherical_coordinate_cmplx_type), intent(in) :: rhs
            type (list_spherical_coordinate_cmplx_type) :: rv

            ! auxiliary
            integer :: i

            allocate(rv%coordinate(lbound(rhs%coordinate,1):ubound(rhs%coordinate,1)))

            !$OMP PARALLEL DO PRIVATE(i)
            do i=lbound(rhs%coordinate,1), ubound(rhs%coordinate,1)
                rv%coordinate(i) = lhs * rhs%coordinate(i)
            end do
            !$OMP END PARALLEL DO

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

                !$OMP PARALLEL DO PRIVATE(i)
                do i=lbound(lhs%coordinate,1), ubound(lhs%coordinate,1)
                    rv%coordinate(i) = lhs%coordinate(i) / rhs(i)
                end do
                !$OMP END PARALLEL DO

            end if

        end function

        function lib_math_list_spherical_operator_c_array_divide_by_c_array(lhs, rhs) result(rv)
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

                !$OMP PARALLEL DO PRIVATE(i)
                do i=lbound(lhs%coordinate,1), ubound(lhs%coordinate,1)
                    rv%coordinate(i) = lhs%coordinate(i) / rhs(i)
                end do
                !$OMP END PARALLEL DO

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

            !$OMP PARALLEL DO PRIVATE(i)
            do i=lbound(lhs%coordinate,1), ubound(lhs%coordinate,1)
                rv%coordinate(i) = lhs%coordinate(i) / rhs
            end do
            !$OMP END PARALLEL DO

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

            !$OMP PARALLEL DO PRIVATE(i)
            do i=lbound(lhs%coordinate,1), ubound(lhs%coordinate,1)
                rv%coordinate(i) = lhs%coordinate(i) / rhs
            end do
            !$OMP END PARALLEL DO

        end function

! ---- single cartesian coordinate ----
        function lib_math_cartesian_operator_add(lhs, rhs) result(rv)
            implicit none
            ! dummy
            type (cartesian_coordinate_cmplx_type), intent(in) :: lhs
            type (cartesian_coordinate_cmplx_type), intent(in) :: rhs
            type (cartesian_coordinate_cmplx_type) :: rv

            rv%x = lhs%x + rhs%x
            rv%y = lhs%y + rhs%y
            rv%z = lhs%z + rhs%z

        end function

        function lib_math_cartesian_operator_sub(lhs, rhs) result(rv)
            implicit none
            ! dummy
            type (cartesian_coordinate_cmplx_type), intent(in) :: lhs
            type (cartesian_coordinate_cmplx_type), intent(in) :: rhs
            type (cartesian_coordinate_cmplx_type) :: rv

            rv%x = lhs%x - rhs%x
            rv%y = lhs%y - rhs%y
            rv%z = lhs%z - rhs%z

        end function

        function lib_math_cartesian_operator_scalar_real_mul(lhs, rhs) result(rv)
            implicit none
            ! dummy
            real(kind=8), intent(in) :: lhs
            type (cartesian_coordinate_cmplx_type), intent(in) :: rhs
            type (cartesian_coordinate_cmplx_type) :: rv

            rv%x = lhs * rhs%x
            rv%y = lhs * rhs%y
            rv%z = lhs * rhs%z

        end function

        function lib_math_cartesian_operator_scalar_cmplx_mul(lhs, rhs) result(rv)
            implicit none
            ! dummy
            complex(kind=8), intent(in) :: lhs
            type (cartesian_coordinate_cmplx_type), intent(in) :: rhs
            type (cartesian_coordinate_cmplx_type) :: rv

            rv%x = lhs * rhs%x
            rv%y = lhs * rhs%y
            rv%z = lhs * rhs%z

        end function

        function lib_math_cartesian_operator_mul_scalar_real(lhs, rhs) result(rv)
            implicit none
            ! dummy
            type (cartesian_coordinate_cmplx_type), intent(in) :: lhs
            real(kind=8), intent(in) :: rhs
            type (cartesian_coordinate_cmplx_type) :: rv

            rv%x = lhs%x * rhs
            rv%y = lhs%y * rhs
            rv%z = lhs%z * rhs

        end function

        function lib_math_cartesian_operator_mul_scalar_cmplx(lhs, rhs) result(rv)
            implicit none
            ! dummy
            type (cartesian_coordinate_cmplx_type), intent(in) :: lhs
            complex(kind=8), intent(in) :: rhs
            type (cartesian_coordinate_cmplx_type) :: rv

            rv%x = lhs%x * rhs
            rv%y = lhs%y * rhs
            rv%z = lhs%z * rhs

        end function

        function lib_math_cartesian_operator_divide_by_scalar_real(lhs, rhs) result(rv)
            implicit none
            ! dummy
            type (cartesian_coordinate_cmplx_type), intent(in) :: lhs
            real(kind=8), intent(in) :: rhs
            type (cartesian_coordinate_cmplx_type) :: rv

            rv%x = lhs%x / rhs
            rv%y = lhs%y / rhs
            rv%z = lhs%z / rhs

        end function

        function lib_math_cartesian_operator_divide_by_scalar_cmplx(lhs, rhs) result(rv)
            implicit none
            ! dummy
            type (cartesian_coordinate_cmplx_type), intent(in) :: lhs
            complex(kind=8), intent(in) :: rhs
            type (cartesian_coordinate_cmplx_type) :: rv

            rv%x = lhs%x / rhs
            rv%y = lhs%y / rhs
            rv%z = lhs%z / rhs

        end function

        function lib_math_cartesian_operator_abs_real(rhs) result(rv)
            implicit none
            ! dummy
            type (cartesian_coordinate_real_type), intent(in) :: rhs

            real(kind=lib_math_type_kind) :: rv

            rv = sqrt(rhs%x*rhs%x + rhs%y*rhs%y + rhs%z*rhs%z)

        end function

        function lib_math_cartesian_operator_abs_cmplx(rhs) result(rv)
            implicit none
            ! dummy
            type (cartesian_coordinate_cmplx_type), intent(in) :: rhs

            complex(kind=8) :: rv

            rv = sqrt(rhs%x*rhs%x + rhs%y*rhs%y + rhs%z*rhs%z)

        end function

! ---- array cartesian coordinate ----
        function lib_math_cartesian_operator_add_array(lhs, rhs) result(rv)
            implicit none
            ! dummy
            type (cartesian_coordinate_cmplx_type), dimension(:), intent(in) :: lhs
            type (cartesian_coordinate_cmplx_type), dimension(lbound(lhs,1):ubound(lhs,1)), intent(in) :: rhs
            type (cartesian_coordinate_cmplx_type), dimension(lbound(lhs,1):ubound(lhs,1)) :: rv

            ! auxiliary
            integer :: i

            !$OMP PARALLEL DO PRIVATE(i)
            do i=lbound(lhs,1), ubound(lhs,1)
                rv(i) = lhs(i) + rhs(i)
            end do
            !$OMP END PARALLEL DO

        end function

        function lib_math_cartesian_operator_sub_array(lhs, rhs) result(rv)
            implicit none
            ! dummy
            type (cartesian_coordinate_cmplx_type), dimension(:), intent(in) :: lhs
            type (cartesian_coordinate_cmplx_type), dimension(lbound(lhs,1):ubound(lhs,1)), intent(in) :: rhs
            type (cartesian_coordinate_cmplx_type), dimension(lbound(lhs,1):ubound(lhs,1)) :: rv

            ! auxiliary
            integer :: i

            !$OMP PARALLEL DO PRIVATE(i)
            do i=lbound(lhs,1), ubound(lhs,1)
                rv(i) = lhs(i) - rhs(i)
            end do
            !$OMP END PARALLEL DO

        end function

        function lib_math_cartesian_operator_0d_add_array(lhs, rhs) result(rv)
            implicit none
            ! dummy
            type (cartesian_coordinate_cmplx_type), intent(in) :: lhs
            type (cartesian_coordinate_cmplx_type), dimension(:), intent(in) :: rhs
            type (cartesian_coordinate_cmplx_type), dimension(lbound(rhs,1):ubound(rhs,1)) :: rv

            ! auxiliary
            integer :: i

            !$OMP PARALLEL DO PRIVATE(i)
            do i=lbound(rhs,1), ubound(rhs,1)
                rv(i) = lhs + rhs(i)
            end do
            !$OMP END PARALLEL DO

        end function

        function lib_math_cartesian_operator_array_add_0d(lhs, rhs) result(rv)
            implicit none
            ! dummy
            type (cartesian_coordinate_cmplx_type), dimension(:), intent(in) :: lhs
            type (cartesian_coordinate_cmplx_type), intent(in) :: rhs
            type (cartesian_coordinate_cmplx_type), dimension(lbound(lhs,1):ubound(lhs,1)) :: rv

            ! auxiliary
            integer :: i

            !$OMP PARALLEL DO PRIVATE(i)
            do i=lbound(lhs,1), ubound(lhs,1)
                rv(i) = lhs(i) + rhs
            end do
            !$OMP END PARALLEL DO

        end function

        function lib_math_cartesian_operator_array_sub_0d(lhs, rhs) result(rv)
            implicit none
            ! dummy
            type (cartesian_coordinate_cmplx_type), dimension(:), intent(in) :: lhs
            type (cartesian_coordinate_cmplx_type), intent(in) :: rhs
            type (cartesian_coordinate_cmplx_type), dimension(lbound(lhs,1):ubound(lhs,1)) :: rv

            ! auxiliary
            integer :: i

            !$OMP PARALLEL DO PRIVATE(i)
            do i=lbound(lhs,1), ubound(lhs,1)
                rv(i) = lhs(i) - rhs
            end do
            !$OMP END PARALLEL DO

        end function

        function lib_math_cartesian_operator_array_real_mul_array(lhs, rhs) result(rv)
            implicit none
            ! dummy
            real(kind=lib_math_type_kind), dimension(:), intent(in) :: lhs
            type (cartesian_coordinate_cmplx_type), dimension(lbound(lhs,1):ubound(lhs,1)), intent(in) :: rhs
            type (cartesian_coordinate_cmplx_type), dimension(lbound(lhs,1):ubound(lhs,1)) :: rv

            ! auxiliary
            integer :: i

            !$OMP PARALLEL DO PRIVATE(i)
            do i=lbound(lhs,1), ubound(lhs,1)
                rv(i) = lhs(i) * rhs(i)
            end do
            !$OMP END PARALLEL DO

        end function

        function lib_math_cartesian_operator_array_cmplx_mul_array(lhs, rhs) result(rv)
            implicit none
            ! dummy
            complex(kind=lib_math_type_kind), dimension(:), intent(in) :: lhs
            type (cartesian_coordinate_cmplx_type), dimension(lbound(lhs,1):ubound(lhs,1)), intent(in) :: rhs
            type (cartesian_coordinate_cmplx_type), dimension(lbound(lhs,1):ubound(lhs,1)) :: rv

            ! auxiliary
            integer :: i

            !$OMP PARALLEL DO PRIVATE(i)
            do i=lbound(lhs,1), ubound(lhs,1)
                rv(i) = lhs(i) * rhs(i)
            end do
            !$OMP END PARALLEL DO

        end function

        function lib_math_cartesian_operator_0d_real_mul_array(lhs, rhs) result(rv)
            implicit none
            ! dummy
            real(kind=lib_math_type_kind), intent(in) :: lhs
            type (cartesian_coordinate_cmplx_type), dimension(:), intent(in) :: rhs
            type (cartesian_coordinate_cmplx_type), dimension(lbound(rhs,1):ubound(rhs,1)) :: rv

            ! auxiliary
            integer :: i

            !$OMP PARALLEL DO PRIVATE(i)
            do i=lbound(rhs,1), ubound(rhs,1)
                rv(i) = lhs * rhs(i)
            end do
            !$OMP END PARALLEL DO

        end function

        function lib_math_cartesian_operator_0d_cmplx_mul_array(lhs, rhs) result(rv)
            implicit none
            ! dummy
            complex(kind=lib_math_type_kind), intent(in) :: lhs
            type (cartesian_coordinate_cmplx_type), dimension(:), intent(in) :: rhs
            type (cartesian_coordinate_cmplx_type), dimension(lbound(rhs,1):ubound(rhs,1)) :: rv

            ! auxiliary
            integer :: i

            !$OMP PARALLEL DO PRIVATE(i)
            do i=lbound(rhs,1), ubound(rhs,1)
                rv(i) = lhs * rhs(i)
            end do
            !$OMP END PARALLEL DO

        end function

        function lib_math_cartesian_operator_array_divide_by_scalar_array_real(lhs, rhs) result(rv)
            implicit none
            ! dummy
            type (cartesian_coordinate_cmplx_type), dimension(:), intent(in) :: lhs
            real(kind=lib_math_type_kind), dimension(lbound(lhs,1):ubound(lhs,1)), intent(in) :: rhs
            type (cartesian_coordinate_cmplx_type), dimension(lbound(lhs,1):ubound(lhs,1)) :: rv

            ! auxiliary
            integer :: i

            !$OMP PARALLEL DO PRIVATE(i)
            do i=lbound(lhs,1), ubound(lhs,1)
                rv(i) = lhs(i) / rhs(i)
            end do
            !$OMP END PARALLEL DO

        end function

        function lib_math_cartesian_operator_array_divide_by_scalar_array_cmplx(lhs, rhs) result(rv)
            implicit none
            ! dummy
            type (cartesian_coordinate_cmplx_type), dimension(:), intent(in) :: lhs
            complex(kind=lib_math_type_kind), dimension(lbound(lhs,1):ubound(lhs,1)), intent(in) :: rhs
            type (cartesian_coordinate_cmplx_type), dimension(lbound(lhs,1):ubound(lhs,1)) :: rv

            ! auxiliary
            integer :: i

            !$OMP PARALLEL DO PRIVATE(i)
            do i=lbound(lhs,1), ubound(lhs,1)
                rv(i) = lhs(i) / rhs(i)
            end do
            !$OMP END PARALLEL DO

        end function

        function lib_math_cartesian_operator_array_divide_by_0d_real(lhs, rhs) result(rv)
            implicit none
            ! dummy
            type (cartesian_coordinate_cmplx_type), dimension(:), intent(in) :: lhs
            real(kind=lib_math_type_kind), intent(in) :: rhs
            type (cartesian_coordinate_cmplx_type), dimension(lbound(lhs,1):ubound(lhs,1)) :: rv

            ! auxiliary
            integer :: i

            !$OMP PARALLEL DO PRIVATE(i)
            do i=lbound(lhs,1), ubound(lhs,1)
                rv(i) = lhs(i) / rhs
            end do
            !$OMP END PARALLEL DO

        end function

        function lib_math_cartesian_operator_array_divide_by_0d_cmplx(lhs, rhs) result(rv)
            implicit none
            ! dummy
            type (cartesian_coordinate_cmplx_type), dimension(:), intent(in) :: lhs
            complex(kind=lib_math_type_kind), intent(in) :: rhs
            type (cartesian_coordinate_cmplx_type), dimension(lbound(lhs,1):ubound(lhs,1)) :: rv

            ! auxiliary
            integer :: i

            !$OMP PARALLEL DO PRIVATE(i)
            do i=lbound(lhs,1), ubound(lhs,1)
                rv(i) = lhs(i) / rhs
            end do
            !$OMP END PARALLEL DO

        end function

! ---- list_cartesian_coordinate ----
        function lib_math_list_cartesian_operator_add_array_cmplx(lhs, rhs) result(rv)
            implicit none
            ! dummy
            type (list_cartesian_coordinate_cmplx_type), intent(in) :: lhs
            type (list_cartesian_coordinate_cmplx_type), intent(in) :: rhs
            type (list_cartesian_coordinate_cmplx_type) :: rv

            ! auxiliary
            integer :: i

            if (lbound(lhs%coordinate,1) .eq. lbound(rhs%coordinate,1) .and. &
                ubound(lhs%coordinate,1) .eq. ubound(rhs%coordinate,1)) then

                allocate(rv%coordinate(lbound(rhs%coordinate,1):ubound(rhs%coordinate,1)))

                !$OMP PARALLEL DO PRIVATE(i)
                do i=lbound(lhs%coordinate,1), ubound(lhs%coordinate,1)
                    rv%coordinate(i) = lhs%coordinate(i) + rhs%coordinate(i)
                end do
                !$OMP END PARALLEL DO
            end if

        end function

        function lib_math_list_cartesian_operator_sub_array_cmplx(lhs, rhs) result(rv)
            implicit none
            ! dummy
            type (list_cartesian_coordinate_cmplx_type), intent(in) :: lhs
            type (list_cartesian_coordinate_cmplx_type), intent(in) :: rhs
            type (list_cartesian_coordinate_cmplx_type) :: rv

            ! auxiliary
            integer :: i

            if (lbound(lhs%coordinate,1) .eq. lbound(rhs%coordinate,1) .and. &
                ubound(lhs%coordinate,1) .eq. ubound(rhs%coordinate,1)) then

                allocate(rv%coordinate(lbound(rhs%coordinate,1):ubound(rhs%coordinate,1)))

                !$OMP PARALLEL DO PRIVATE(i)
                do i=lbound(lhs%coordinate,1), ubound(lhs%coordinate,1)
                    rv%coordinate(i) = lhs%coordinate(i) - rhs%coordinate(i)
                end do
                !$OMP END PARALLEL DO
            end if
        end function

        function lib_math_list_cartesian_operator_0d_add_array_cmplx(lhs, rhs) result(rv)
            implicit none
            ! dummy
            type (cartesian_coordinate_cmplx_type), intent(in) :: lhs
            type (list_cartesian_coordinate_cmplx_type), intent(in) :: rhs
            type (list_cartesian_coordinate_cmplx_type) :: rv

            ! auxiliary
            integer :: i

            allocate(rv%coordinate(lbound(rhs%coordinate,1):ubound(rhs%coordinate,1)))

            !$OMP PARALLEL DO PRIVATE(i)
            do i=lbound(rhs%coordinate,1), ubound(rhs%coordinate,1)
                rv%coordinate(i) = lhs + rhs%coordinate(i)
            end do
            !$OMP END PARALLEL DO

        end function

        function lib_math_list_cartesian_operator_array_add_0d_cmplx(lhs, rhs) result(rv)
            implicit none
            ! dummy
            type (list_cartesian_coordinate_cmplx_type), intent(in) :: lhs
            type (cartesian_coordinate_cmplx_type), intent(in) :: rhs
            type (list_cartesian_coordinate_cmplx_type) :: rv

            ! auxiliary
            integer :: i

            allocate(rv%coordinate(lbound(lhs%coordinate,1):ubound(lhs%coordinate,1)))

            !$OMP PARALLEL DO PRIVATE(i)
            do i=lbound(lhs%coordinate,1), ubound(lhs%coordinate,1)
                rv%coordinate(i) = lhs%coordinate(i) + rhs
            end do
            !$OMP END PARALLEL DO

        end function

        function lib_math_list_cartesian_operator_array_sub_0d_cmplx(lhs, rhs) result(rv)
            implicit none
            ! dummy
            type (list_cartesian_coordinate_cmplx_type), intent(in) :: lhs
            type (cartesian_coordinate_cmplx_type), intent(in) :: rhs
            type (list_cartesian_coordinate_cmplx_type) :: rv

            ! auxiliary
            integer :: i

            allocate(rv%coordinate(lbound(lhs%coordinate,1):ubound(lhs%coordinate,1)))

            !$OMP PARALLEL DO PRIVATE(i)
            do i=lbound(lhs%coordinate,1), ubound(lhs%coordinate,1)
                rv%coordinate(i) = lhs%coordinate(i) - rhs
            end do
            !$OMP END PARALLEL DO

        end function

        function lib_math_list_cartesian_operator_array_real_mul_array_c(lhs, rhs) result(rv)
            implicit none
            ! dummy
            real(kind=lib_math_type_kind), dimension(:), intent(in) :: lhs
            type (list_cartesian_coordinate_cmplx_type), intent(in) :: rhs
            type (list_cartesian_coordinate_cmplx_type) :: rv

            ! auxiliary
            integer :: i

            if (lbound(lhs,1) .eq. lbound(rhs%coordinate,1) .and. &
                ubound(lhs,1) .eq. ubound(rhs%coordinate,1)) then

                allocate(rv%coordinate(lbound(lhs,1):ubound(lhs,1)))

                !$OMP PARALLEL DO PRIVATE(i)
                do i=lbound(lhs,1), ubound(lhs,1)
                    rv%coordinate(i) = lhs(i) * rhs%coordinate(i)
                end do
                !$OMP END PARALLEL DO
            else
                print *, "lib_math_list_cartesian_operator_array_real_mul_array_c: ERROR"
            end if

        end function

        function lib_math_list_cartesian_operator_array_c_mul_array_cmplx(lhs, rhs) result(rv)
            implicit none
            ! dummy
            complex(kind=lib_math_type_kind), dimension(:), intent(in) :: lhs
            type (list_cartesian_coordinate_cmplx_type), intent(in) :: rhs
            type (list_cartesian_coordinate_cmplx_type) :: rv

            ! auxiliary
            integer :: i

            if (lbound(lhs,1) .eq. lbound(rhs%coordinate,1) .and. &
                ubound(lhs,1) .eq. ubound(rhs%coordinate,1)) then

                allocate(rv%coordinate(lbound(rhs%coordinate,1):ubound(rhs%coordinate,1)))

                !$OMP PARALLEL DO PRIVATE(i)
                do i=lbound(lhs,1), ubound(lhs,1)
                    rv%coordinate(i) = lhs(i) * rhs%coordinate(i)
                end do
                !$OMP END PARALLEL DO

            end if

        end function

        function lib_math_list_cartesian_operator_real_mul_array_cmplx(lhs, rhs) result(rv)
            implicit none
            ! dummy
            real(kind=lib_math_type_kind), intent(in) :: lhs
            type (list_cartesian_coordinate_cmplx_type), intent(in) :: rhs
            type (list_cartesian_coordinate_cmplx_type) :: rv

            ! auxiliary
            integer :: i

            allocate(rv%coordinate(lbound(rhs%coordinate,1):ubound(rhs%coordinate,1)))

            !$OMP PARALLEL DO PRIVATE(i)
            do i=lbound(rhs%coordinate,1), ubound(rhs%coordinate,1)
                rv%coordinate(i) = lhs * rhs%coordinate(i)
            end do
            !$OMP END PARALLEL DO

        end function

        function lib_math_list_cartesian_operator_cmplx_mul_array_cmplx(lhs, rhs) result(rv)
            implicit none
            ! dummy
            complex(kind=lib_math_type_kind), intent(in) :: lhs
            type (list_cartesian_coordinate_cmplx_type), intent(in) :: rhs
            type (list_cartesian_coordinate_cmplx_type) :: rv

            ! auxiliary
            integer :: i

            allocate(rv%coordinate(lbound(rhs%coordinate,1):ubound(rhs%coordinate,1)))

            !$OMP PARALLEL DO PRIVATE(i)
            do i=lbound(rhs%coordinate,1), ubound(rhs%coordinate,1)
                rv%coordinate(i) = lhs * rhs%coordinate(i)
            end do
            !$OMP END PARALLEL DO

        end function

        function lib_math_list_cartesian_operator_array_divide_by_real_array(lhs, rhs) result(rv)
            implicit none
            ! dummy
            type (list_cartesian_coordinate_cmplx_type), intent(in) :: lhs
            real(kind=lib_math_type_kind), dimension(:), intent(in) :: rhs
            type (list_cartesian_coordinate_cmplx_type) :: rv

            ! auxiliary
            integer :: i

            if (lbound(lhs%coordinate,1) .eq. lbound(rhs,1) .and. &
                ubound(lhs%coordinate,1) .eq. ubound(rhs,1)) then

                allocate(rv%coordinate(lbound(lhs%coordinate,1):ubound(lhs%coordinate,1)))

                !$OMP PARALLEL DO PRIVATE(i)
                do i=lbound(lhs%coordinate,1), ubound(lhs%coordinate,1)
                    rv%coordinate(i) = lhs%coordinate(i) / rhs(i)
                end do
                !$OMP END PARALLEL DO

            end if

        end function

        function lib_math_list_cartesian_operator_c_array_divide_by_c_array(lhs, rhs) result(rv)
            implicit none
            ! dummy
            type (list_cartesian_coordinate_cmplx_type), intent(in) :: lhs
            complex(kind=lib_math_type_kind), dimension(:), intent(in) :: rhs
            type (list_cartesian_coordinate_cmplx_type) :: rv

            ! auxiliary
            integer :: i

            if (lbound(lhs%coordinate,1) .eq. lbound(rhs,1) .and. &
                ubound(lhs%coordinate,1) .eq. ubound(rhs,1)) then

                allocate(rv%coordinate(lbound(lhs%coordinate,1):ubound(lhs%coordinate,1)))

                !$OMP PARALLEL DO PRIVATE(i)
                do i=lbound(lhs%coordinate,1), ubound(lhs%coordinate,1)
                    rv%coordinate(i) = lhs%coordinate(i) / rhs(i)
                end do
                !$OMP END PARALLEL DO

            end if

        end function

        function lib_math_list_cartesian_operator_array_divide_by_real(lhs, rhs) result(rv)
            implicit none
            ! dummy
            type (list_cartesian_coordinate_cmplx_type), intent(in) :: lhs
            real(kind=lib_math_type_kind), intent(in) :: rhs
            type (list_cartesian_coordinate_cmplx_type) :: rv

            ! auxiliary
            integer :: i

            allocate(rv%coordinate(lbound(lhs%coordinate,1):ubound(lhs%coordinate,1)))

            !$OMP PARALLEL DO PRIVATE(i)
            do i=lbound(lhs%coordinate,1), ubound(lhs%coordinate,1)
                rv%coordinate(i) = lhs%coordinate(i) / rhs
            end do
            !$OMP END PARALLEL DO

        end function

        function lib_math_list_cartesian_operator_array_divide_by_cmplx(lhs, rhs) result(rv)
            implicit none
            ! dummy
            type (list_cartesian_coordinate_cmplx_type), intent(in) :: lhs
            complex(kind=lib_math_type_kind), intent(in) :: rhs
            type (list_cartesian_coordinate_cmplx_type) :: rv

            ! auxiliary
            integer :: i

            allocate(rv%coordinate(lbound(lhs%coordinate,1):ubound(lhs%coordinate,1)))

            !$OMP PARALLEL DO PRIVATE(i)
            do i=lbound(lhs%coordinate,1), ubound(lhs%coordinate,1)
                rv%coordinate(i) = lhs%coordinate(i) / rhs
            end do
            !$OMP END PARALLEL DO

        end function

! ---- list_list ----

        ! Arguments
        ! ----
        !   list: type (list_list_logical)
        !       derived data type to initialize
        !   fnu: integer
        !       start index
        !   n: integer
        !       number of elements, n .GE. 1
        !
        subroutine lib_math_list_list_logical_init(list, fnu, n)
            implicit none
            ! dummy
            type (list_list_logical), intent(inout) :: list
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

        end subroutine lib_math_list_list_logical_init

        ! Arguments
        ! ----
        !   list: type (list_list_integer_sys)
        !       derived data type to initialize
        !   fnu: integer
        !       start index
        !   n: integer
        !       number of elements, n .GE. 1
        !
        subroutine lib_math_list_list_integer_sys_init(list, fnu, n)
            implicit none
            ! dummy
            type (list_list_integer_sys), intent(inout) :: list
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

        end subroutine lib_math_list_list_integer_sys_init

        ! Arguments
        ! ----
        !   list: type (list_list_integer)
        !       derived data type to initialize
        !   fnu: integer
        !       start index
        !   n: integer
        !       number of elements, n .GE. 1
        !
        subroutine lib_math_list_list_integer_init(list, fnu, n)
            implicit none
            ! dummy
            type (list_list_integer), intent(inout) :: list
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

        end subroutine lib_math_list_list_integer_init

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

        ! reference: http://mathworld.wolfram.com/SphericalCoordinates.html
        ! ISO 31-11
        subroutine lib_math_spherical_point_to_cartesian_point(lhs, rhs)
            implicit none
            ! dummy
            type(spherical_coordinate_real_type), intent(in) :: rhs

            type(cartesian_coordinate_real_type), intent(inout) :: lhs

            ! auxiliary
            real(kind=lib_math_type_kind) :: cos_theta
            real(kind=lib_math_type_kind) :: sin_theta
            real(kind=lib_math_type_kind) :: cos_phi
            real(kind=lib_math_type_kind) :: sin_phi

            real(kind=lib_math_type_kind) :: r_sin_theta

            cos_theta = cos(rhs%theta)
            sin_theta = sin(rhs%theta)
            cos_phi = cos(rhs%phi)
            sin_phi = sin(rhs%phi)

            r_sin_theta = rhs%rho * sin_theta

            lhs%x = r_sin_theta * cos_phi
            lhs%y = r_sin_theta * sin_phi
            lhs%z = rhs%rho * cos_theta

        end subroutine lib_math_spherical_point_to_cartesian_point

        ! reference: http://mathworld.wolfram.com/SphericalCoordinates.html
        ! ISO 31-11
        subroutine lib_math_cartesian_point_to_spherical_point(lhs, rhs)
            implicit none
            ! dummy
            type(cartesian_coordinate_real_type), intent(in) :: rhs

            type(spherical_coordinate_real_type), intent(inout) :: lhs

            ! auxiliary
            real(kind=lib_math_type_kind) :: r

            r = sqrt(rhs%x * rhs%x + rhs%y * rhs%y + rhs%z * rhs%z)

            lhs%rho = r
            lhs%theta = acos(rhs%z / r)
            lhs%phi = atan(rhs%y, rhs%x)

        end subroutine lib_math_cartesian_point_to_spherical_point

        ! Converts a complex vector field from a sherical to a cartesian coordinate system
        !
        ! Argument
        ! ----
        !   rhs: type(spherical_coordinate_cmplx_type)
        !       complex vector (spherical coordinate system)
        !   theta: real
        !       polar coordinate of the complex vector (rhs)
        !   phi: real
        !       azimuthal coordinate of the complex vector (rhs)
        !
        ! Returns
        ! ----
        !   lhs: type(cartesian_coordinate_cmplx_type)
        !       complex vector (cartesian coordinate system)
        !
        ! LaTeX:  $$ \left[v_{r}(\theta, \varphi) v_{\theta}(\theta, \varphi) v_{\varphi}(\theta, \varphi)\right]=\left[v_{x}(\theta, \varphi) v_{y}(\theta, \varphi) v_{z}(\theta, \varphi)\right) ] \mathbb{T}(\theta, \varphi) $$
        !         $$ \mathbb{T}(\theta, \varphi)=\left[\begin{array}{ccc}{\sin \theta \cos \varphi} & {\sin \theta \sin \varphi} & {\cos \theta} \\ {\cos \theta \cos \varphi} & {\cos \theta \sin \varphi} & {-\sin \theta} \\ {-\sin \varphi} & {\cos \varphi} & {0}\end{array}\right] $$
        !
        ! Reference: Accurate calculation of spherical and vector spherical harmonic
        !            expansions via spectral element grids, Wang^2, Xie, eq. 3.23
        function lib_math_spherical_components_to_cartesian_components_cmplx_a(rhs, theta, phi) result (lhs)
            implicit none
            ! dummy
            type(spherical_coordinate_cmplx_type), intent(in) :: rhs
            real(kind=lib_math_type_kind), intent(in) :: theta
            real(kind=lib_math_type_kind), intent(in) :: phi

            type(cartesian_coordinate_cmplx_type) :: lhs

            ! auxiliary
            real(kind=lib_math_type_kind) :: cos_theta
            real(kind=lib_math_type_kind) :: sin_theta
            real(kind=lib_math_type_kind) :: cos_phi
            real(kind=lib_math_type_kind) :: sin_phi

            cos_theta = cos(theta)
            sin_theta = sin(theta)
            cos_phi = cos(phi)
            sin_phi = sin(phi)

!            lhs%x = rhs%rho * (sin_theta * cos_phi) &
!                    + rhs%theta * (sin_theta * sin_phi) &
!                    + rhs%phi * cos_theta
!            lhs%y = rhs%rho * (cos_theta * cos_phi) &
!                      + rhs%theta * (cos_theta * sin_phi) &
!                      - rhs%phi * sin_theta
!            lhs%z = - rhs%rho * sin_phi &
!                    + rhs%theta * cos_phi

            lhs%x = rhs%rho * (sin_theta * cos_phi) &
                    + rhs%theta * (cos_theta * cos_phi) &
                    - rhs%phi * sin_phi
            lhs%y = rhs%rho * (sin_theta * sin_phi) &
                      + rhs%theta * (cos_theta * sin_phi) &
                      + rhs%phi * cos_phi
            lhs%z = rhs%rho * cos_theta &
                    - rhs%theta * sin_theta

        end function lib_math_spherical_components_to_cartesian_components_cmplx_a

        ! Converts a complex vector field from a sherical to a cartesian coordinate system
        !
        ! Argument
        ! ----
        !   rhs: type(spherical_coordinate_cmplx_type)
        !       complex vector (spherical coordinate system)
        !   coordinate: type(cartesian_coordinate_real_type)
        !       cartesian coordinate of the complex vector (rhs)
        !
        ! Returns
        ! ----
        !   lhs: type(cartesian_coordinate_cmplx_type)
        !       complex vector (cartesian coordinate system)
        !
        ! LaTeX:  $$ \left[v_{r}(\theta, \varphi) v_{\theta}(\theta, \varphi) v_{\varphi}(\theta, \varphi)\right]=\left[v_{x}(\theta, \varphi) v_{y}(\theta, \varphi) v_{z}(\theta, \varphi)\right) ] \mathbb{T}(\theta, \varphi) $$
        !         $$ \mathbb{T}(\theta, \varphi)=\left[\begin{array}{ccc}{\sin \theta \cos \varphi} & {\sin \theta \sin \varphi} & {\cos \theta} \\ {\cos \theta \cos \varphi} & {\cos \theta \sin \varphi} & {-\sin \theta} \\ {-\sin \varphi} & {\cos \varphi} & {0}\end{array}\right] $$
        !
        ! Reference: Accurate calculation of spherical and vector spherical harmonic
        !            expansions via spectral element grids, Wang^2, Xie, eq. 3.23
        function lib_math_spherical_components_to_cartesian_components_cmplx_c(rhs, coordinate) &
                                                                                      result (lhs)
            implicit none
            ! dummy
            type(spherical_coordinate_cmplx_type), intent(in) :: rhs
            type(cartesian_coordinate_real_type), intent(in) :: coordinate

            type(cartesian_coordinate_cmplx_type) :: lhs

            ! auxiliary
            type(spherical_coordinate_real_type) :: coordinate_spherical

            coordinate_spherical = coordinate

            lhs = lib_math_spherical_components_to_cartesian_components_cmplx_a(rhs, &
                                                                                coordinate_spherical%theta, &
                                                                                coordinate_spherical%phi)
        end function lib_math_spherical_components_to_cartesian_components_cmplx_c

        ! Converts a complex vector field from a sherical to a cartesian coordinate system
        !
        ! Argument
        ! ----
        !   rhs: type(spherical_coordinate_cmplx_type)
        !       complex vector (spherical coordinate system)
        !   coordinate: type(spherical_coordinate_real_type)
        !       spherical coordinate of the complex vector (rhs)
        !
        ! Returns
        ! ----
        !   lhs: type(cartesian_coordinate_cmplx_type)
        !       complex vector (cartesian coordinate system)
        !
        ! LaTeX:  $$ \left[v_{r}(\theta, \varphi) v_{\theta}(\theta, \varphi) v_{\varphi}(\theta, \varphi)\right]=\left[v_{x}(\theta, \varphi) v_{y}(\theta, \varphi) v_{z}(\theta, \varphi)\right) ] \mathbb{T}(\theta, \varphi) $$
        !         $$ \mathbb{T}(\theta, \varphi)=\left[\begin{array}{ccc}{\sin \theta \cos \varphi} & {\sin \theta \sin \varphi} & {\cos \theta} \\ {\cos \theta \cos \varphi} & {\cos \theta \sin \varphi} & {-\sin \theta} \\ {-\sin \varphi} & {\cos \varphi} & {0}\end{array}\right] $$
        !
        ! Reference: Accurate calculation of spherical and vector spherical harmonic
        !            expansions via spectral element grids, Wang^2, Xie, eq. 3.23
        function lib_math_spherical_components_to_cartesian_components_cmplx_s(rhs, coordinate) &
                                                                                      result (lhs)
            implicit none
            ! dummy
            type(spherical_coordinate_cmplx_type), intent(in) :: rhs
            type(spherical_coordinate_real_type), intent(in) :: coordinate

            type(cartesian_coordinate_cmplx_type) :: lhs

            lhs = lib_math_spherical_components_to_cartesian_components_cmplx_a(rhs, &
                                                                                    coordinate%theta, &
                                                                                    coordinate%phi)
        end function lib_math_spherical_components_to_cartesian_components_cmplx_s

        ! Converts a complex vector field from a cartesian to a spherical coordinate system
        !
        ! Argument
        ! ----
        !   rhs: type(cartesian_coordinate_cmplx_type)
        !       complex vector (cartesian coordinate system)
        !   theta: real
        !       polar coordinate of the complex vector (rhs)
        !   phi: real
        !       azimuthal coordinate of the complex vector (rhs)
        !
        ! Returns
        ! ----
        !   lhs: type(spherical_coordinate_cmplx_type)
        !       complex vector (spherical coordinate system)
        !
        ! LaTeX:  $$ \left[v_{r}(\theta, \varphi) v_{\theta}(\theta, \varphi) v_{\varphi}(\theta, \varphi)\right]=\left[v_{x}(\theta, \varphi) v_{y}(\theta, \varphi) v_{z}(\theta, \varphi)\right) ] \mathbb{T}(\theta, \varphi) $$
        !         $$ \mathbb{T}(\theta, \varphi)=\left[\begin{array}{ccc}{\sin \theta \cos \varphi} & {\sin \theta \sin \varphi} & {\cos \theta} \\ {\cos \theta \cos \varphi} & {\cos \theta \sin \varphi} & {-\sin \theta} \\ {-\sin \varphi} & {\cos \varphi} & {0}\end{array}\right] $$
        !
        ! Reference: Accurate calculation of spherical and vector spherical harmonic
        !            expansions via spectral element grids, Wang^2, Xie, eq. 3.23
        function lib_math_cartesian_components_to_spherical_components_cmplx_a(rhs, theta, phi) result (lhs)
            implicit none
            ! dummy
            type(cartesian_coordinate_cmplx_type), intent(in) :: rhs
            real(kind=lib_math_type_kind), intent(in) :: theta
            real(kind=lib_math_type_kind), intent(in) :: phi

            type(spherical_coordinate_cmplx_type) :: lhs

            ! auxiliary
            complex(kind=lib_math_type_kind) :: cos_theta
            complex(kind=lib_math_type_kind) :: sin_theta
            complex(kind=lib_math_type_kind) :: cos_phi
            complex(kind=lib_math_type_kind) :: sin_phi

            cos_theta = cos(theta)
            sin_theta = sin(theta)
            cos_phi = cos(phi)
            sin_phi = sin(phi)

!            lhs%rho = rhs%x * (sin_theta * cos_phi) &
!                    + rhs%y * (cos_theta * cos_phi) &
!                    - rhs%z * sin_phi
!            lhs%theta = rhs%x * (sin_theta * sin_phi) &
!                      + rhs%y * (cos_theta * sin_phi) &
!                      + rhs%z * cos_phi
!            lhs%phi = rhs%x * cos_theta &
!                    - rhs%y * sin_theta

            lhs%rho = rhs%x * (sin_theta *cos_phi) &
                    + rhs%y * (sin_theta * sin_phi) &
                    + rhs%z * cos_theta
            lhs%theta = rhs%x * (cos_theta * cos_phi) &
                      + rhs%y * (cos_theta * sin_phi) &
                      - rhs%z * sin_theta
            lhs%phi = - rhs%x * sin_phi &
                      + rhs%y * cos_phi

        end function lib_math_cartesian_components_to_spherical_components_cmplx_a

        ! Converts a complex vector field from a cartesian to a spherical coordinate system
        !
        ! Argument
        ! ----
        !   rhs: type(cartesian_coordinate_cmplx_type)
        !       complex vector (cartesian coordinate system)
        !   coordinate: type(cartesian_coordinate_real_type)
        !       coordinate of the complex vector (rhs)
        !
        ! Returns
        ! ----
        !   lhs: type(spherical_coordinate_cmplx_type)
        !       complex vector (spherical coordinate system)
        !
        ! LaTeX:  $$ \left[v_{r}(\theta, \varphi) v_{\theta}(\theta, \varphi) v_{\varphi}(\theta, \varphi)\right]=\left[v_{x}(\theta, \varphi) v_{y}(\theta, \varphi) v_{z}(\theta, \varphi)\right) ] \mathbb{T}(\theta, \varphi) $$
        !         $$ \mathbb{T}(\theta, \varphi)=\left[\begin{array}{ccc}{\sin \theta \cos \varphi} & {\sin \theta \sin \varphi} & {\cos \theta} \\ {\cos \theta \cos \varphi} & {\cos \theta \sin \varphi} & {-\sin \theta} \\ {-\sin \varphi} & {\cos \varphi} & {0}\end{array}\right] $$
        !
        ! Reference: Accurate calculation of spherical and vector spherical harmonic
        !            expansions via spectral element grids, Wang^2, Xie, eq. 3.23
        function lib_math_cartesian_components_to_spherical_components_cmplx_c(rhs, coordinate) &
                                                                                      result (lhs)
            implicit none
            ! dummy
            type(cartesian_coordinate_cmplx_type), intent(in) :: rhs
            type(cartesian_coordinate_real_type), intent(in) :: coordinate

            type(spherical_coordinate_cmplx_type) :: lhs

            ! auxiliary
            type(spherical_coordinate_real_type) :: coordinate_spherical

            coordinate_spherical = coordinate

            lhs = lib_math_cartesian_components_to_spherical_components_cmplx_a(rhs, &
                                                                                coordinate_spherical%theta, &
                                                                                coordinate_spherical%phi)
        end function lib_math_cartesian_components_to_spherical_components_cmplx_c

        ! Converts a complex vector field from a cartesian to a spherical coordinate system
        !
        ! Argument
        ! ----
        !   rhs: type(cartesian_coordinate_cmplx_type)
        !       complex vector (cartesian coordinate system)
        !   coordinate: type(spherical_coordinate_real_type)
        !       coordinate of the complex vector (rhs)
        !
        ! Returns
        ! ----
        !   lhs: type(spherical_coordinate_cmplx_type)
        !       complex vector (spherical coordinate system)
        !
        ! LaTeX:  $$ \left[v_{r}(\theta, \varphi) v_{\theta}(\theta, \varphi) v_{\varphi}(\theta, \varphi)\right]=\left[v_{x}(\theta, \varphi) v_{y}(\theta, \varphi) v_{z}(\theta, \varphi)\right) ] \mathbb{T}(\theta, \varphi) $$
        !         $$ \mathbb{T}(\theta, \varphi)=\left[\begin{array}{ccc}{\sin \theta \cos \varphi} & {\sin \theta \sin \varphi} & {\cos \theta} \\ {\cos \theta \cos \varphi} & {\cos \theta \sin \varphi} & {-\sin \theta} \\ {-\sin \varphi} & {\cos \varphi} & {0}\end{array}\right] $$
        !
        ! Reference: Accurate calculation of spherical and vector spherical harmonic
        !            expansions via spectral element grids, Wang^2, Xie, eq. 3.23
        function lib_math_cartesian_components_to_spherical_components_cmplx_s(rhs, coordinate) &
                                                                                      result (lhs)
            implicit none
            ! dummy
            type(cartesian_coordinate_cmplx_type), intent(in) :: rhs
            type(spherical_coordinate_real_type), intent(in) :: coordinate

            type(spherical_coordinate_cmplx_type) :: lhs

            lhs = lib_math_cartesian_components_to_spherical_components_cmplx_a(rhs, &
                                                                                coordinate%theta, &
                                                                                coordinate%phi)
        end function lib_math_cartesian_components_to_spherical_components_cmplx_s

        function lib_math_cartesian_cross_product_real (lhs, rhs) result (rv)
            implicit none
            ! dummy
            type(cartesian_coordinate_real_type), intent(in) :: lhs
            type(cartesian_coordinate_real_type), intent(in) :: rhs

            type(cartesian_coordinate_real_type) :: rv

            rv%x = lhs%y * rhs%z - lhs%z * rhs%y
            rv%y = lhs%z * rhs%x - lhs%x * rhs%z
            rv%z = lhs%x * rhs%y - lhs%y * rhs%x

        end function

        function lib_math_cartesian_cross_product_real_list (lhs, rhs) result (rv)
            implicit none
            ! dummy
            type(cartesian_coordinate_real_type), dimension(:), intent(in) :: lhs
            type(cartesian_coordinate_real_type), dimensioN(lbound(lhs,1):ubound(lhs,1)), intent(in) :: rhs

            type(cartesian_coordinate_real_type), dimension(lbound(lhs,1):ubound(lhs,1)) :: rv

            ! auxiliary
            integer :: i

            !$OMP PARALLEL DO PRIVATE(i)
            do i=lbound(lhs,1), ubound(lhs,1)
                rv(i) = lib_math_cartesian_cross_product_real(lhs(i), rhs(i))
            end do
            !$OMP END PARALLEL DO

        end function

        function lib_math_cartesian_cross_product_cmplx (lhs, rhs) result (rv)
            implicit none
            ! dummy
            type(cartesian_coordinate_cmplx_type), intent(in) :: lhs
            type(cartesian_coordinate_cmplx_type), intent(in) :: rhs

            type(cartesian_coordinate_cmplx_type) :: rv

            rv%x = lhs%y * rhs%z - lhs%z * rhs%y
            rv%y = lhs%z * rhs%x - lhs%x * rhs%z
            rv%z = lhs%x * rhs%y - lhs%y * rhs%x

        end function

        function lib_math_cartesian_cross_product_cmplx_list (lhs, rhs) result (rv)
            implicit none
            ! dummy
            type(cartesian_coordinate_cmplx_type), dimension(:), intent(in) :: lhs
            type(cartesian_coordinate_cmplx_type), dimensioN(lbound(lhs,1):ubound(lhs,1)), intent(in) :: rhs

            type(cartesian_coordinate_cmplx_type), dimension(lbound(lhs,1):ubound(lhs,1)) :: rv

            ! auxiliary
            integer :: i

            !$OMP PARALLEL DO PRIVATE(i)
            do i=lbound(lhs,1), ubound(lhs,1)
                rv(i) = lib_math_cartesian_cross_product_cmplx(lhs(i), rhs(i))
            end do
            !$OMP END PARALLEL DO

        end function

        function lib_math_cartesian_dot_product_real (lhs, rhs) result (rv)
            implicit none
            ! dummy
            type(cartesian_coordinate_real_type), intent(in) :: lhs
            type(cartesian_coordinate_real_type), intent(in) :: rhs

            real(kind=lib_math_type_kind) :: rv

            rv = lhs%x * rhs%x + lhs%y * rhs%y + lhs%z * rhs%z

        end function

        function lib_math_cartesian_dot_product_real_list (lhs, rhs) result (rv)
            implicit none
            ! dummy
            type(cartesian_coordinate_real_type), dimension(:), intent(in) :: lhs
            type(cartesian_coordinate_real_type), dimensioN(lbound(lhs,1):ubound(lhs,1)), intent(in) :: rhs

            real(kind=lib_math_type_kind), dimension(lbound(lhs,1):ubound(lhs,1)) :: rv

            ! auxiliary
            integer :: i

            !$OMP PARALLEL DO PRIVATE(i)
            do i=lbound(lhs,1), ubound(lhs,1)
                rv(i) = lib_math_cartesian_dot_product_real(lhs(i), rhs(i))
            end do
            !$OMP END PARALLEL DO

        end function

        function lib_math_cartesian_dot_product_cmplx (lhs, rhs) result (rv)
            implicit none
            ! dummy
            type(cartesian_coordinate_cmplx_type), intent(in) :: lhs
            type(cartesian_coordinate_cmplx_type), intent(in) :: rhs

            complex(kind=lib_math_type_kind) :: rv

            rv = lhs%x * rhs%x + lhs%y * rhs%y + lhs%z * rhs%z

        end function

        function lib_math_cartesian_dot_product_cmplx_list (lhs, rhs) result (rv)
            implicit none
            ! dummy
            type(cartesian_coordinate_cmplx_type), dimension(:), intent(in) :: lhs
            type(cartesian_coordinate_cmplx_type), dimensioN(lbound(lhs,1):ubound(lhs,1)), intent(in) :: rhs

            complex(kind=lib_math_type_kind), dimension(lbound(lhs,1):ubound(lhs,1)) :: rv

            ! auxiliary
            integer :: i

            !$OMP PARALLEL DO PRIVATE(i)
            do i=lbound(lhs,1), ubound(lhs,1)
                rv(i) = lib_math_cartesian_dot_product_cmplx(lhs(i), rhs(i))
            end do
            !$OMP END PARALLEL DO

        end function

        function lib_math_type_operator_test_functions() result (rv)
            implicit none
            ! dummy
            integer :: rv

            ! parameter
            double precision, parameter :: PI=4.D0*atan(1.D0)   ! maximum precision, platform independet

            rv = 0

            if (.not. test_lib_math_list_spherical_operator_add_array_cmplx()) then
                rv = rv + 1
            end if
            if (.not. test_lib_math_list_spherical_operator_sub_array_cmplx()) then
                rv = rv + 1
            end if
            if (.not. test_lib_math_list_spherical_operator_0d_add_array_cmplx()) then
                rv = rv + 1
            end if
            if (.not. test_lib_math_list_spherical_operator_array_add_0d_cmplx()) then
                rv = rv + 1
            end if
            if (.not. test_lib_math_list_spherical_operator_array_sub_0d_cmplx()) then
                rv = rv + 1
            end if
            if (.not. test_lib_math_list_spherical_operator_array_real_mul_array_c()) then
                rv = rv + 1
            end if
            if (.not. test_lib_math_list_spherical_operator_array_c_mul_array_cmplx()) then
                rv = rv + 1
            end if
            if (.not. test_lib_math_list_spherical_operator_real_mul_array_cmplx()) then
                rv = rv + 1
            end if
            if (.not. test_lib_math_list_spherical_operator_cmplx_mul_array_cmplx()) then
                rv = rv + 1
            end if
            if (.not. test_lib_math_list_spherical_operator_c_array_divide_by_r_array()) then
                rv = rv + 1
            end if
            if (.not. test_lib_math_list_spherical_operator_c_array_divide_by_c_array()) then
                rv = rv + 1
            end if
            if (.not. test_lib_math_list_spherical_operator_array_divide_by_real()) then
                rv = rv + 1
            end if
            if (.not. test_lib_math_list_spherical_operator_array_divide_by_cmplx()) then
                rv = rv + 1
            end if

            if (.not. test_lib_math_list_cartesian_operator_add_array_cmplx()) then
                rv = rv + 1
            end if
            if (.not. test_lib_math_list_cartesian_operator_sub_array_cmplx()) then
                rv = rv + 1
            end if
            if (.not. test_lib_math_list_cartesian_operator_0d_add_array_cmplx()) then
                rv = rv + 1
            end if
            if (.not. test_lib_math_list_cartesian_operator_array_add_0d_cmplx()) then
                rv = rv + 1
            end if
            if (.not. test_lib_math_list_cartesian_operator_array_sub_0d_cmplx()) then
                rv = rv + 1
            end if
            if (.not. test_lib_math_list_cartesian_operator_array_real_mul_array_c()) then
                rv = rv + 1
            end if
            if (.not. test_lib_math_list_cartesian_operator_array_c_mul_array_cmplx()) then
                rv = rv + 1
            end if
            if (.not. test_lib_math_list_cartesian_operator_real_mul_array_cmplx()) then
                rv = rv + 1
            end if
            if (.not. test_lib_math_list_cartesian_operator_cmplx_mul_array_cmplx()) then
                rv = rv + 1
            end if
            if (.not. test_lib_math_list_cartesian_operator_c_array_divide_by_r_array()) then
                rv = rv + 1
            end if
            if (.not. test_lib_math_list_cartesian_operator_c_array_divide_by_c_array()) then
                rv = rv + 1
            end if
            if (.not. test_lib_math_list_cartesian_operator_array_divide_by_real()) then
                rv = rv + 1
            end if
            if (.not. test_lib_math_list_cartesian_operator_array_divide_by_cmplx()) then
                rv = rv + 1
            end if

            if (.not. test_lib_math_cartesian_point_to_spherical_point()) then
                rv = rv + 1
            end if
            if (.not. test_lib_math_spherical_point_to_cartesian_point()) then
                rv = rv + 1
            end if
            if (.not. test_lib_math_spherical_components_to_cartesian_components_c_a()) then
                rv = rv + 1
            end if
            if (.not. test_lib_math_cartesian_components_to_spherical_components_c_a()) then
                rv = rv + 1
            end if
            if (.not. test_lib_math_cartesian_cross_product_real_list()) then
                rv = rv + 1
            end if
            if (.not. test_lib_math_cartesian_cross_product_cmplx_list()) then
                rv = rv + 1
            end if
            if (.not. test_lib_math_cartesian_dot_product_real_list()) then
                rv = rv + 1
            end if
            if (.not. test_lib_math_cartesian_dot_product_cmplx_list()) then
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

                function get_test_values_cartesian_1() result(rv)
                    implicit none
                    ! dummy
                    type (cartesian_coordinate_cmplx_type) :: rv

                    rv%y = cmplx(1,2, kind=8)
                    rv%x = cmplx(3,4, kind=8)
                    rv%z = cmplx(5,6, kind=8)
                end function

                function get_test_values_cartesian_2() result(rv)
                    implicit none
                    ! dummy
                    type (cartesian_coordinate_cmplx_type) :: rv

                    rv%y = cmplx(6,5, kind=8)
                    rv%x = cmplx(4,3, kind=8)
                    rv%z = cmplx(2,1, kind=8)
                end function

                function get_test_values_spherical_1() result(rv)
                    implicit none
                    ! dummy
                    type (spherical_coordinate_cmplx_type) :: rv

                    rv%phi = cmplx(1,2, kind=8)
                    rv%rho = cmplx(3,4, kind=8)
                    rv%theta = cmplx(5,6, kind=8)
                end function

                function get_test_values_spherical_2() result(rv)
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

                function test_lib_math_list_spherical_operator_add_array_cmplx() result (rv)
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
                    lhs%coordinate(1) = get_test_values_spherical_1()
                    lhs%coordinate(2) = get_test_values_spherical_1()

                    ! rhs
                    rhs%coordinate(1) = get_test_values_spherical_1()
                    rhs%coordinate(2) = get_test_values_spherical_2()

                    ! ground truth
                    ground_truth_value%coordinate(1)%phi = cmplx(2,4, kind=8)
                    ground_truth_value%coordinate(1)%rho = cmplx(6,8, kind=8)
                    ground_truth_value%coordinate(1)%theta = cmplx(10,12, kind=8)

                    ground_truth_value%coordinate(2)%phi = cmplx(7,7, kind=8)
                    ground_truth_value%coordinate(2)%rho = cmplx(7,7, kind=8)
                    ground_truth_value%coordinate(2)%theta = cmplx(7,7, kind=8)

                    value = lib_math_list_spherical_operator_add_array_cmplx(lhs, rhs)

                    rv = .true.
                    print *, "test_lib_math_list_spherical_operator_add_array_cmplx:"
                    do i=1, d
                        buffer = abs(value%coordinate(i)%phi - ground_truth_value%coordinate(i)%phi)
                        if (.not. evaluate(buffer, i, 'phi')) rv = .false.

                        buffer = abs(value%coordinate(i)%rho - ground_truth_value%coordinate(i)%rho)
                        if (.not. evaluate(buffer, i, "rho")) rv = .false.

                        buffer = abs(value%coordinate(i)%theta - ground_truth_value%coordinate(i)%theta)
                        if (.not. evaluate(buffer, i, "theta")) rv = .false.
                    end do

                end function test_lib_math_list_spherical_operator_add_array_cmplx

                function test_lib_math_list_spherical_operator_sub_array_cmplx() result (rv)
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
                    lhs%coordinate(1) = get_test_values_spherical_1()
                    lhs%coordinate(2) = get_test_values_spherical_1()

                    ! rhs
                    rhs%coordinate(1) = get_test_values_spherical_1()
                    rhs%coordinate(2) = get_test_values_spherical_2()

                    ! ground truth
                    ground_truth_value%coordinate(1)%phi = cmplx(0,0, kind=8)
                    ground_truth_value%coordinate(1)%rho = cmplx(0,0, kind=8)
                    ground_truth_value%coordinate(1)%theta = cmplx(0,0, kind=8)

                    ground_truth_value%coordinate(2)%phi = cmplx(-5,-3, kind=8)
                    ground_truth_value%coordinate(2)%rho = cmplx(-1,1, kind=8)
                    ground_truth_value%coordinate(2)%theta = cmplx(3,5, kind=8)

                    value = lib_math_list_spherical_operator_sub_array_cmplx(lhs, rhs)

                    rv = .true.
                    print *, "test_lib_math_list_spherical_operator_sub_array_cmplx:"
                    do i=1, d
                        buffer = abs(value%coordinate(i)%phi - ground_truth_value%coordinate(i)%phi)
                        if (.not. evaluate(buffer, i, 'phi')) rv = .false.

                        buffer = abs(value%coordinate(i)%rho - ground_truth_value%coordinate(i)%rho)
                        if (.not. evaluate(buffer, i, "rho")) rv = .false.

                        buffer = abs(value%coordinate(i)%theta - ground_truth_value%coordinate(i)%theta)
                        if (.not. evaluate(buffer, i, "theta")) rv = .false.
                    end do

                end function test_lib_math_list_spherical_operator_sub_array_cmplx

                function test_lib_math_list_spherical_operator_0d_add_array_cmplx() result (rv)
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
                    lhs = get_test_values_spherical_1()

                    ! rhs
                    rhs%coordinate(1) = get_test_values_spherical_1()
                    rhs%coordinate(2) = get_test_values_spherical_2()

                    ! ground truth
                    ground_truth_value%coordinate(1)%phi = cmplx(2,4, kind=8)
                    ground_truth_value%coordinate(1)%rho = cmplx(6,8, kind=8)
                    ground_truth_value%coordinate(1)%theta = cmplx(10,12, kind=8)

                    ground_truth_value%coordinate(2)%phi = cmplx(7,7, kind=8)
                    ground_truth_value%coordinate(2)%rho = cmplx(7,7, kind=8)
                    ground_truth_value%coordinate(2)%theta = cmplx(7,7, kind=8)

                    value = lib_math_list_spherical_operator_0d_add_array_cmplx(lhs, rhs)

                    rv = .true.
                    print *, "test_lib_math_list_spherical_operator_0d_add_array_cmplx:"
                    do i=1, d
                        buffer = abs(value%coordinate(i)%phi - ground_truth_value%coordinate(i)%phi)
                        if (.not. evaluate(buffer, i, 'phi')) rv = .false.

                        buffer = abs(value%coordinate(i)%rho - ground_truth_value%coordinate(i)%rho)
                        if (.not. evaluate(buffer, i, "rho")) rv = .false.

                        buffer = abs(value%coordinate(i)%theta - ground_truth_value%coordinate(i)%theta)
                        if (.not. evaluate(buffer, i, "theta")) rv = .false.
                    end do

                end function test_lib_math_list_spherical_operator_0d_add_array_cmplx

                function test_lib_math_list_spherical_operator_array_add_0d_cmplx() result (rv)
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
                    lhs%coordinate(1) = get_test_values_spherical_1()
                    lhs%coordinate(2) = get_test_values_spherical_2()

                    ! rhs
                    rhs = get_test_values_spherical_1()

                    ! ground truth
                    ground_truth_value%coordinate(1)%phi = cmplx(2,4, kind=8)
                    ground_truth_value%coordinate(1)%rho = cmplx(6,8, kind=8)
                    ground_truth_value%coordinate(1)%theta = cmplx(10,12, kind=8)

                    ground_truth_value%coordinate(2)%phi = cmplx(7,7, kind=8)
                    ground_truth_value%coordinate(2)%rho = cmplx(7,7, kind=8)
                    ground_truth_value%coordinate(2)%theta = cmplx(7,7, kind=8)

                    value = lib_math_list_spherical_operator_array_add_0d_cmplx(lhs, rhs)

                    rv = .true.
                    print *, "test_lib_math_list_spherical_operator_array_add_0d_cmplx:"
                    do i=1, d
                        buffer = abs(value%coordinate(i)%phi - ground_truth_value%coordinate(i)%phi)
                        if (.not. evaluate(buffer, i, 'phi')) rv = .false.

                        buffer = abs(value%coordinate(i)%rho - ground_truth_value%coordinate(i)%rho)
                        if (.not. evaluate(buffer, i, "rho")) rv = .false.

                        buffer = abs(value%coordinate(i)%theta - ground_truth_value%coordinate(i)%theta)
                        if (.not. evaluate(buffer, i, "theta")) rv = .false.
                    end do

                end function test_lib_math_list_spherical_operator_array_add_0d_cmplx

                function test_lib_math_list_spherical_operator_array_sub_0d_cmplx() result (rv)
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
                    lhs%coordinate(1) = get_test_values_spherical_1()
                    lhs%coordinate(2) = get_test_values_spherical_2()

                    ! rhs
                    rhs = get_test_values_spherical_1()

                    ! ground truth
                    ground_truth_value%coordinate(1)%phi = cmplx(0,0, kind=8)
                    ground_truth_value%coordinate(1)%rho = cmplx(0,0, kind=8)
                    ground_truth_value%coordinate(1)%theta = cmplx(0,0, kind=8)

                    ground_truth_value%coordinate(2)%phi = cmplx(5,3, kind=8)
                    ground_truth_value%coordinate(2)%rho = cmplx(1,-1, kind=8)
                    ground_truth_value%coordinate(2)%theta = cmplx(-3,-5, kind=8)

                    value = lib_math_list_spherical_operator_array_sub_0d_cmplx(lhs, rhs)

                    rv = .true.
                    print *, "test_lib_math_list_spherical_operator_array_sub_0d_cmplx:"
                    do i=1, d
                        buffer = abs(value%coordinate(i)%phi - ground_truth_value%coordinate(i)%phi)
                        if (.not. evaluate(buffer, i, 'phi')) rv = .false.

                        buffer = abs(value%coordinate(i)%rho - ground_truth_value%coordinate(i)%rho)
                        if (.not. evaluate(buffer, i, "rho")) rv = .false.

                        buffer = abs(value%coordinate(i)%theta - ground_truth_value%coordinate(i)%theta)
                        if (.not. evaluate(buffer, i, "theta")) rv = .false.
                    end do

                end function test_lib_math_list_spherical_operator_array_sub_0d_cmplx

                function test_lib_math_list_spherical_operator_array_real_mul_array_c() result (rv)
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
                    rhs%coordinate(1) = get_test_values_spherical_1()
                    rhs%coordinate(2) = get_test_values_spherical_2()

                    ! ground truth
                    ground_truth_value%coordinate(1)%phi = cmplx(2,4, kind=8)
                    ground_truth_value%coordinate(1)%rho = cmplx(6,8, kind=8)
                    ground_truth_value%coordinate(1)%theta = cmplx(10,12, kind=8)

                    ground_truth_value%coordinate(2)%phi = cmplx(30,25, kind=8)
                    ground_truth_value%coordinate(2)%rho = cmplx(20,15, kind=8)
                    ground_truth_value%coordinate(2)%theta = cmplx(10,5, kind=8)

                    value = lib_math_list_spherical_operator_array_real_mul_array_c(lhs, rhs)

                    rv = .true.
                    print *, "test_lib_math_list_spherical_operator_array_real_mul_array_c:"
                    do i=1, d
                        buffer = abs(value%coordinate(i)%phi - ground_truth_value%coordinate(i)%phi)
                        if (.not. evaluate(buffer, i, 'phi')) rv = .false.

                        buffer = abs(value%coordinate(i)%rho - ground_truth_value%coordinate(i)%rho)
                        if (.not. evaluate(buffer, i, "rho")) rv = .false.

                        buffer = abs(value%coordinate(i)%theta - ground_truth_value%coordinate(i)%theta)
                        if (.not. evaluate(buffer, i, "theta")) rv = .false.
                    end do

                end function test_lib_math_list_spherical_operator_array_real_mul_array_c

                function test_lib_math_list_spherical_operator_array_c_mul_array_cmplx() result (rv)
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
                    rhs%coordinate(1) = get_test_values_spherical_1()
                    rhs%coordinate(2) = get_test_values_spherical_2()

                    ! ground truth
                    ground_truth_value%coordinate(1)%phi = cmplx(-2,1, kind=8)
                    ground_truth_value%coordinate(1)%rho = cmplx(-4,3, kind=8)
                    ground_truth_value%coordinate(1)%theta = cmplx(-6,5, kind=8)

                    ground_truth_value%coordinate(2)%phi = cmplx(12,10, kind=8)
                    ground_truth_value%coordinate(2)%rho = cmplx(8,6, kind=8)
                    ground_truth_value%coordinate(2)%theta = cmplx(4,2, kind=8)

                    value = lib_math_list_spherical_operator_array_c_mul_array_cmplx(lhs, rhs)

                    rv = .true.
                    print *, "test_lib_math_list_spherical_operator_array_c_mul_array_cmplx:"
                    do i=1, d
                        buffer = abs(value%coordinate(i)%phi - ground_truth_value%coordinate(i)%phi)
                        if (.not. evaluate(buffer, i, 'phi')) rv = .false.

                        buffer = abs(value%coordinate(i)%rho - ground_truth_value%coordinate(i)%rho)
                        if (.not. evaluate(buffer, i, "rho")) rv = .false.

                        buffer = abs(value%coordinate(i)%theta - ground_truth_value%coordinate(i)%theta)
                        if (.not. evaluate(buffer, i, "theta")) rv = .false.
                    end do

                end function test_lib_math_list_spherical_operator_array_c_mul_array_cmplx

                function test_lib_math_list_spherical_operator_real_mul_array_cmplx() result (rv)
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
                    rhs%coordinate(1) = get_test_values_spherical_1()
                    rhs%coordinate(2) = get_test_values_spherical_2()

                    ! ground truth
                    ground_truth_value%coordinate(1)%phi = cmplx(4,8, kind=8)
                    ground_truth_value%coordinate(1)%rho = cmplx(12,16, kind=8)
                    ground_truth_value%coordinate(1)%theta = cmplx(20,24, kind=8)

                    ground_truth_value%coordinate(2)%phi = cmplx(24,20, kind=8)
                    ground_truth_value%coordinate(2)%rho = cmplx(16,12, kind=8)
                    ground_truth_value%coordinate(2)%theta = cmplx(8,4, kind=8)

                    value = lib_math_list_spherical_operator_real_mul_array_cmplx(lhs, rhs)

                    rv = .true.
                    print *, "test_lib_math_list_spherical_operator_real_mul_array_cmplx:"
                    do i=1, d
                        buffer = abs(value%coordinate(i)%phi - ground_truth_value%coordinate(i)%phi)
                        if (.not. evaluate(buffer, i, 'phi')) rv = .false.

                        buffer = abs(value%coordinate(i)%rho - ground_truth_value%coordinate(i)%rho)
                        if (.not. evaluate(buffer, i, "rho")) rv = .false.

                        buffer = abs(value%coordinate(i)%theta - ground_truth_value%coordinate(i)%theta)
                        if (.not. evaluate(buffer, i, "theta")) rv = .false.
                    end do

                end function test_lib_math_list_spherical_operator_real_mul_array_cmplx

                function test_lib_math_list_spherical_operator_cmplx_mul_array_cmplx() result (rv)
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
                    rhs%coordinate(1) = get_test_values_spherical_1()
                    rhs%coordinate(2) = get_test_values_spherical_2()

                    ! ground truth
                    ground_truth_value%coordinate(1)%phi = cmplx(-8,4, kind=8)
                    ground_truth_value%coordinate(1)%rho = cmplx(-16,12, kind=8)
                    ground_truth_value%coordinate(1)%theta = cmplx(-24,20, kind=8)

                    ground_truth_value%coordinate(2)%phi = cmplx(-20,24, kind=8)
                    ground_truth_value%coordinate(2)%rho = cmplx(-12,16, kind=8)
                    ground_truth_value%coordinate(2)%theta = cmplx(-4,8, kind=8)

                    value = lib_math_list_spherical_operator_cmplx_mul_array_cmplx(lhs, rhs)

                    rv = .true.
                    print *, "test_lib_math_list_spherical_operator_cmplx_mul_array_cmplx:"
                    do i=1, d
                        buffer = abs(value%coordinate(i)%phi - ground_truth_value%coordinate(i)%phi)
                        if (.not. evaluate(buffer, i, 'phi')) rv = .false.

                        buffer = abs(value%coordinate(i)%rho - ground_truth_value%coordinate(i)%rho)
                        if (.not. evaluate(buffer, i, "rho")) rv = .false.

                        buffer = abs(value%coordinate(i)%theta - ground_truth_value%coordinate(i)%theta)
                        if (.not. evaluate(buffer, i, "theta")) rv = .false.
                    end do

                end function test_lib_math_list_spherical_operator_cmplx_mul_array_cmplx

                function test_lib_math_list_spherical_operator_c_array_divide_by_r_array() result (rv)
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
                    lhs%coordinate(1) = get_test_values_spherical_1()
                    lhs%coordinate(2) = get_test_values_spherical_2()

                    ! rhs
                    rhs(1) = 2.0
                    rhs(2) = 4.0

                    ! ground truth
                    ground_truth_value%coordinate(1)%phi = cmplx(0.5,1, kind=8)
                    ground_truth_value%coordinate(1)%rho = cmplx(1.5,2, kind=8)
                    ground_truth_value%coordinate(1)%theta = cmplx(2.5,3, kind=8)

                    ground_truth_value%coordinate(2)%phi = cmplx(1.5,1.25, kind=8)
                    ground_truth_value%coordinate(2)%rho = cmplx(1,0.75, kind=8)
                    ground_truth_value%coordinate(2)%theta = cmplx(0.5,0.25, kind=8)

                    value = lib_math_list_spherical_operator_array_divide_by_real_array(lhs, rhs)

                    rv = .true.
                    print *, "test_lib_math_list_spherical_operator_c_array_divide_by_r_array:"
                    do i=1, d
                        buffer = abs(value%coordinate(i)%phi - ground_truth_value%coordinate(i)%phi)
                        if (.not. evaluate(buffer, i, 'phi')) rv = .false.

                        buffer = abs(value%coordinate(i)%rho - ground_truth_value%coordinate(i)%rho)
                        if (.not. evaluate(buffer, i, "rho")) rv = .false.

                        buffer = abs(value%coordinate(i)%theta - ground_truth_value%coordinate(i)%theta)
                        if (.not. evaluate(buffer, i, "theta")) rv = .false.
                    end do

                end function test_lib_math_list_spherical_operator_c_array_divide_by_r_array

                function test_lib_math_list_spherical_operator_c_array_divide_by_c_array() result (rv)
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
                    lhs%coordinate(1) = get_test_values_spherical_1()
                    lhs%coordinate(2) = get_test_values_spherical_2()

                    ! rhs
                    rhs(1) = cmplx(0, 2, kind=8)
                    rhs(2) = cmplx(0, 4, kind=8)

                    ! ground truth
                    ground_truth_value%coordinate(1)%phi = cmplx(1, -0.5, kind=8)
                    ground_truth_value%coordinate(1)%rho = cmplx(2, -1.5, kind=8)
                    ground_truth_value%coordinate(1)%theta = cmplx(3, -2.5, kind=8)

                    ground_truth_value%coordinate(2)%phi = cmplx(1.25, -1.5, kind=8)
                    ground_truth_value%coordinate(2)%rho = cmplx(0.75, -1, kind=8)
                    ground_truth_value%coordinate(2)%theta = cmplx(0.25, -0.5, kind=8)

                    value = lib_math_list_spherical_operator_c_array_divide_by_c_array(lhs, rhs)

                    rv = .true.
                    print *, "test_lib_math_list_spherical_operator_c_array_divide_by_c_array:"
                    do i=1, d
                        buffer = abs(value%coordinate(i)%phi - ground_truth_value%coordinate(i)%phi)
                        if (.not. evaluate(buffer, i, 'phi')) rv = .false.

                        buffer = abs(value%coordinate(i)%rho - ground_truth_value%coordinate(i)%rho)
                        if (.not. evaluate(buffer, i, "rho")) rv = .false.

                        buffer = abs(value%coordinate(i)%theta - ground_truth_value%coordinate(i)%theta)
                        if (.not. evaluate(buffer, i, "theta")) rv = .false.
                    end do

                end function test_lib_math_list_spherical_operator_c_array_divide_by_c_array

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
                    lhs%coordinate(1) = get_test_values_spherical_1()
                    lhs%coordinate(2) = get_test_values_spherical_2()

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
                        if (.not. evaluate(buffer, i, 'phi')) rv = .false.

                        buffer = abs(value%coordinate(i)%rho - ground_truth_value%coordinate(i)%rho)
                        if (.not. evaluate(buffer, i, "rho")) rv = .false.

                        buffer = abs(value%coordinate(i)%theta - ground_truth_value%coordinate(i)%theta)
                        if (.not. evaluate(buffer, i, "theta")) rv = .false.
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
                    lhs%coordinate(1) = get_test_values_spherical_1()
                    lhs%coordinate(2) = get_test_values_spherical_2()

                    ! rhs
                    rhs = cmplx(0, 2, kind=8)

                    ! ground truth
                    ground_truth_value%coordinate(1)%phi = cmplx(1, -0.5, kind=8)
                    ground_truth_value%coordinate(1)%rho = cmplx(2, -1.5, kind=8)
                    ground_truth_value%coordinate(1)%theta = cmplx(3, -2.5, kind=8)

                    ground_truth_value%coordinate(2)%phi = cmplx(2.5, -3, kind=8)
                    ground_truth_value%coordinate(2)%rho = cmplx(1.5, -2, kind=8)
                    ground_truth_value%coordinate(2)%theta = cmplx(0.5, -1, kind=8)

                    value = lib_math_list_spherical_operator_array_divide_by_cmplx(lhs, rhs)

                    rv = .true.
                    print *, "test_lib_math_list_spherical_operator_array_divide_by_cmplx:"
                    do i=1, d
                        buffer = abs(value%coordinate(i)%phi - ground_truth_value%coordinate(i)%phi)
                        if (.not. evaluate(buffer, i, "phi  ")) rv = .false.

                        buffer = abs(value%coordinate(i)%rho - ground_truth_value%coordinate(i)%rho)
                        if (.not. evaluate(buffer, i, "rho  ")) rv = .false.

                        buffer = abs(value%coordinate(i)%theta - ground_truth_value%coordinate(i)%theta)
                        if (.not. evaluate(buffer, i, "theta")) rv = .false.
                    end do

                end function test_lib_math_list_spherical_operator_array_divide_by_cmplx

                function test_lib_math_list_cartesian_operator_add_array_cmplx() result (rv)
                    implicit none
                    ! dummy
                    logical :: rv

                    ! parameter
                    integer(kind=4), parameter :: d = 2

                    ! auxiliary
                    integer(kind=4) :: i
                    real(kind=8) :: buffer
                    type (list_cartesian_coordinate_cmplx_type) :: lhs
                    type (list_cartesian_coordinate_cmplx_type) :: rhs
                    type (list_cartesian_coordinate_cmplx_type) :: value

                    type (list_cartesian_coordinate_cmplx_type) :: ground_truth_value

                    allocate(lhs%coordinate(d))
                    allocate(rhs%coordinate(d))
                    allocate(ground_truth_value%coordinate(d))

                    ! lhs
                    lhs%coordinate(1) = get_test_values_cartesian_1()
                    lhs%coordinate(2) = get_test_values_cartesian_1()

                    ! rhs
                    rhs%coordinate(1) = get_test_values_cartesian_1()
                    rhs%coordinate(2) = get_test_values_cartesian_2()

                    ! ground truth
                    ground_truth_value%coordinate(1)%y = cmplx(2,4, kind=8)
                    ground_truth_value%coordinate(1)%x = cmplx(6,8, kind=8)
                    ground_truth_value%coordinate(1)%z = cmplx(10,12, kind=8)

                    ground_truth_value%coordinate(2)%y = cmplx(7,7, kind=8)
                    ground_truth_value%coordinate(2)%x = cmplx(7,7, kind=8)
                    ground_truth_value%coordinate(2)%z = cmplx(7,7, kind=8)

                    value = lib_math_list_cartesian_operator_add_array_cmplx(lhs, rhs)

                    rv = .true.
                    print *, "test_lib_math_list_cartesian_operator_add_array_cmplx:"
                    do i=1, d
                        buffer = abs(value%coordinate(i)%y - ground_truth_value%coordinate(i)%y)
                        if (.not. evaluate(buffer, i, 'y')) rv = .false.

                        buffer = abs(value%coordinate(i)%x - ground_truth_value%coordinate(i)%x)
                        if (.not. evaluate(buffer, i, "x")) rv = .false.

                        buffer = abs(value%coordinate(i)%z - ground_truth_value%coordinate(i)%z)
                        if (.not. evaluate(buffer, i, "z")) rv = .false.
                    end do

                end function test_lib_math_list_cartesian_operator_add_array_cmplx

                function test_lib_math_list_cartesian_operator_sub_array_cmplx() result (rv)
                    implicit none
                    ! dummy
                    logical :: rv

                    ! parameter
                    integer(kind=4), parameter :: d = 2

                    ! auxiliary
                    integer(kind=4) :: i
                    real(kind=8) :: buffer
                    type (list_cartesian_coordinate_cmplx_type) :: lhs
                    type (list_cartesian_coordinate_cmplx_type) :: rhs
                    type (list_cartesian_coordinate_cmplx_type) :: value

                    type (list_cartesian_coordinate_cmplx_type) :: ground_truth_value

                    allocate(lhs%coordinate(d))
                    allocate(rhs%coordinate(d))
                    allocate(ground_truth_value%coordinate(d))

                    ! lhs
                    lhs%coordinate(1) = get_test_values_cartesian_1()
                    lhs%coordinate(2) = get_test_values_cartesian_1()

                    ! rhs
                    rhs%coordinate(1) = get_test_values_cartesian_1()
                    rhs%coordinate(2) = get_test_values_cartesian_2()

                    ! ground truth
                    ground_truth_value%coordinate(1)%y = cmplx(0,0, kind=8)
                    ground_truth_value%coordinate(1)%x = cmplx(0,0, kind=8)
                    ground_truth_value%coordinate(1)%z = cmplx(0,0, kind=8)

                    ground_truth_value%coordinate(2)%y = cmplx(-5,-3, kind=8)
                    ground_truth_value%coordinate(2)%x = cmplx(-1,1, kind=8)
                    ground_truth_value%coordinate(2)%z = cmplx(3,5, kind=8)

                    value = lib_math_list_cartesian_operator_sub_array_cmplx(lhs, rhs)

                    rv = .true.
                    print *, "test_lib_math_list_cartesian_operator_sub_array_cmplx:"
                    do i=1, d
                        buffer = abs(value%coordinate(i)%y - ground_truth_value%coordinate(i)%y)
                        if (.not. evaluate(buffer, i, 'y')) rv = .false.

                        buffer = abs(value%coordinate(i)%x - ground_truth_value%coordinate(i)%x)
                        if (.not. evaluate(buffer, i, "x")) rv = .false.

                        buffer = abs(value%coordinate(i)%z - ground_truth_value%coordinate(i)%z)
                        if (.not. evaluate(buffer, i, "z")) rv = .false.
                    end do

                end function test_lib_math_list_cartesian_operator_sub_array_cmplx

                function test_lib_math_list_cartesian_operator_0d_add_array_cmplx() result (rv)
                    implicit none
                    ! dummy
                    logical :: rv

                    ! parameter
                    integer(kind=4), parameter :: d = 2

                    ! auxiliary
                    integer(kind=4) :: i
                    real(kind=8) :: buffer
                    type (cartesian_coordinate_cmplx_type) :: lhs
                    type (list_cartesian_coordinate_cmplx_type) :: rhs
                    type (list_cartesian_coordinate_cmplx_type) :: value

                    type (list_cartesian_coordinate_cmplx_type) :: ground_truth_value

                    allocate(rhs%coordinate(d))
                    allocate(ground_truth_value%coordinate(d))

                    ! lhs
                    lhs = get_test_values_cartesian_1()

                    ! rhs
                    rhs%coordinate(1) = get_test_values_cartesian_1()
                    rhs%coordinate(2) = get_test_values_cartesian_2()

                    ! ground truth
                    ground_truth_value%coordinate(1)%y = cmplx(2,4, kind=8)
                    ground_truth_value%coordinate(1)%x = cmplx(6,8, kind=8)
                    ground_truth_value%coordinate(1)%z = cmplx(10,12, kind=8)

                    ground_truth_value%coordinate(2)%y = cmplx(7,7, kind=8)
                    ground_truth_value%coordinate(2)%x = cmplx(7,7, kind=8)
                    ground_truth_value%coordinate(2)%z = cmplx(7,7, kind=8)

                    value = lib_math_list_cartesian_operator_0d_add_array_cmplx(lhs, rhs)

                    rv = .true.
                    print *, "test_lib_math_list_cartesian_operator_0d_add_array_cmplx:"
                    do i=1, d
                        buffer = abs(value%coordinate(i)%y - ground_truth_value%coordinate(i)%y)
                        if (.not. evaluate(buffer, i, 'y')) rv = .false.

                        buffer = abs(value%coordinate(i)%x - ground_truth_value%coordinate(i)%x)
                        if (.not. evaluate(buffer, i, "x")) rv = .false.

                        buffer = abs(value%coordinate(i)%z - ground_truth_value%coordinate(i)%z)
                        if (.not. evaluate(buffer, i, "z")) rv = .false.
                    end do

                end function test_lib_math_list_cartesian_operator_0d_add_array_cmplx

                function test_lib_math_list_cartesian_operator_array_add_0d_cmplx() result (rv)
                    implicit none
                    ! dummy
                    logical :: rv

                    ! parameter
                    integer(kind=4), parameter :: d = 2

                    ! auxiliary
                    integer(kind=4) :: i
                    real(kind=8) :: buffer
                    type (list_cartesian_coordinate_cmplx_type) :: lhs
                    type (cartesian_coordinate_cmplx_type) :: rhs
                    type (list_cartesian_coordinate_cmplx_type) :: value

                    type (list_cartesian_coordinate_cmplx_type) :: ground_truth_value

                    allocate(lhs%coordinate(d))
                    allocate(ground_truth_value%coordinate(d))

                    ! lhs
                    lhs%coordinate(1) = get_test_values_cartesian_1()
                    lhs%coordinate(2) = get_test_values_cartesian_2()

                    ! rhs
                    rhs = get_test_values_cartesian_1()

                    ! ground truth
                    ground_truth_value%coordinate(1)%y = cmplx(2,4, kind=8)
                    ground_truth_value%coordinate(1)%x = cmplx(6,8, kind=8)
                    ground_truth_value%coordinate(1)%z = cmplx(10,12, kind=8)

                    ground_truth_value%coordinate(2)%y = cmplx(7,7, kind=8)
                    ground_truth_value%coordinate(2)%x = cmplx(7,7, kind=8)
                    ground_truth_value%coordinate(2)%z = cmplx(7,7, kind=8)

                    value = lib_math_list_cartesian_operator_array_add_0d_cmplx(lhs, rhs)

                    rv = .true.
                    print *, "test_lib_math_list_cartesian_operator_array_add_0d_cmplx:"
                    do i=1, d
                        buffer = abs(value%coordinate(i)%y - ground_truth_value%coordinate(i)%y)
                        if (.not. evaluate(buffer, i, 'y')) rv = .false.

                        buffer = abs(value%coordinate(i)%x - ground_truth_value%coordinate(i)%x)
                        if (.not. evaluate(buffer, i, "x")) rv = .false.

                        buffer = abs(value%coordinate(i)%z - ground_truth_value%coordinate(i)%z)
                        if (.not. evaluate(buffer, i, "z")) rv = .false.
                    end do

                end function test_lib_math_list_cartesian_operator_array_add_0d_cmplx

                function test_lib_math_list_cartesian_operator_array_sub_0d_cmplx() result (rv)
                    implicit none
                    ! dummy
                    logical :: rv

                    ! parameter
                    integer(kind=4), parameter :: d = 2

                    ! auxiliary
                    integer(kind=4) :: i
                    real(kind=8) :: buffer
                    type (list_cartesian_coordinate_cmplx_type) :: lhs
                    type (cartesian_coordinate_cmplx_type) :: rhs
                    type (list_cartesian_coordinate_cmplx_type) :: value

                    type (list_cartesian_coordinate_cmplx_type) :: ground_truth_value

                    allocate(lhs%coordinate(d))
                    allocate(ground_truth_value%coordinate(d))

                    ! lhs
                    lhs%coordinate(1) = get_test_values_cartesian_1()
                    lhs%coordinate(2) = get_test_values_cartesian_2()

                    ! rhs
                    rhs = get_test_values_cartesian_1()

                    ! ground truth
                    ground_truth_value%coordinate(1)%y = cmplx(0,0, kind=8)
                    ground_truth_value%coordinate(1)%x = cmplx(0,0, kind=8)
                    ground_truth_value%coordinate(1)%z = cmplx(0,0, kind=8)

                    ground_truth_value%coordinate(2)%y = cmplx(5,3, kind=8)
                    ground_truth_value%coordinate(2)%x = cmplx(1,-1, kind=8)
                    ground_truth_value%coordinate(2)%z = cmplx(-3,-5, kind=8)

                    value = lib_math_list_cartesian_operator_array_sub_0d_cmplx(lhs, rhs)

                    rv = .true.
                    print *, "test_lib_math_list_cartesian_operator_array_sub_0d_cmplx:"
                    do i=1, d
                        buffer = abs(value%coordinate(i)%y - ground_truth_value%coordinate(i)%y)
                        if (.not. evaluate(buffer, i, 'y')) rv = .false.

                        buffer = abs(value%coordinate(i)%x - ground_truth_value%coordinate(i)%x)
                        if (.not. evaluate(buffer, i, "x")) rv = .false.

                        buffer = abs(value%coordinate(i)%z - ground_truth_value%coordinate(i)%z)
                        if (.not. evaluate(buffer, i, "z")) rv = .false.
                    end do

                end function test_lib_math_list_cartesian_operator_array_sub_0d_cmplx

                function test_lib_math_list_cartesian_operator_array_real_mul_array_c() result (rv)
                    implicit none
                    ! dummy
                    logical :: rv

                    ! parameter
                    integer(kind=4), parameter :: d = 2

                    ! auxiliary
                    integer(kind=4) :: i
                    real(kind=8) :: buffer
                    real(kind=lib_math_type_kind), dimension(d) :: lhs
                    type (list_cartesian_coordinate_cmplx_type)  :: rhs
                    type (list_cartesian_coordinate_cmplx_type) :: value

                    type (list_cartesian_coordinate_cmplx_type) :: ground_truth_value

                    allocate(rhs%coordinate(d))
                    allocate(ground_truth_value%coordinate(d))

                    ! lhs
                    lhs(1) = 2.0
                    lhs(2) = 5.0

                    ! rhs
                    rhs%coordinate(1) = get_test_values_cartesian_1()
                    rhs%coordinate(2) = get_test_values_cartesian_2()

                    ! ground truth
                    ground_truth_value%coordinate(1)%y = cmplx(2,4, kind=8)
                    ground_truth_value%coordinate(1)%x = cmplx(6,8, kind=8)
                    ground_truth_value%coordinate(1)%z = cmplx(10,12, kind=8)

                    ground_truth_value%coordinate(2)%y = cmplx(30,25, kind=8)
                    ground_truth_value%coordinate(2)%x = cmplx(20,15, kind=8)
                    ground_truth_value%coordinate(2)%z = cmplx(10,5, kind=8)

                    value = lib_math_list_cartesian_operator_array_real_mul_array_c(lhs, rhs)

                    rv = .true.
                    print *, "test_lib_math_list_cartesian_operator_array_real_mul_array_c:"
                    do i=1, d
                        buffer = abs(value%coordinate(i)%y - ground_truth_value%coordinate(i)%y)
                        if (.not. evaluate(buffer, i, 'y')) rv = .false.

                        buffer = abs(value%coordinate(i)%x - ground_truth_value%coordinate(i)%x)
                        if (.not. evaluate(buffer, i, "x")) rv = .false.

                        buffer = abs(value%coordinate(i)%z - ground_truth_value%coordinate(i)%z)
                        if (.not. evaluate(buffer, i, "z")) rv = .false.
                    end do

                end function test_lib_math_list_cartesian_operator_array_real_mul_array_c

                function test_lib_math_list_cartesian_operator_array_c_mul_array_cmplx() result (rv)
                    implicit none
                    ! dummy
                    logical :: rv

                    ! parameter
                    integer(kind=4), parameter :: d = 2

                    ! auxiliary
                    integer(kind=4) :: i
                    real(kind=8) :: buffer
                    complex(kind=lib_math_type_kind), dimension(d) :: lhs
                    type (list_cartesian_coordinate_cmplx_type) :: rhs
                    type (list_cartesian_coordinate_cmplx_type) :: value

                    type (list_cartesian_coordinate_cmplx_type) :: ground_truth_value

                    allocate(rhs%coordinate(d))
                    allocate(ground_truth_value%coordinate(d))

                    ! lhs
                    lhs(1) = cmplx(0, 1, kind=8)
                    lhs(2) = cmplx(2, 0, kind=8)

                    ! rhs
                    rhs%coordinate(1) = get_test_values_cartesian_1()
                    rhs%coordinate(2) = get_test_values_cartesian_2()

                    ! ground truth
                    ground_truth_value%coordinate(1)%y = cmplx(-2,1, kind=8)
                    ground_truth_value%coordinate(1)%x = cmplx(-4,3, kind=8)
                    ground_truth_value%coordinate(1)%z = cmplx(-6,5, kind=8)

                    ground_truth_value%coordinate(2)%y = cmplx(12,10, kind=8)
                    ground_truth_value%coordinate(2)%x = cmplx(8,6, kind=8)
                    ground_truth_value%coordinate(2)%z = cmplx(4,2, kind=8)

                    value = lib_math_list_cartesian_operator_array_c_mul_array_cmplx(lhs, rhs)

                    rv = .true.
                    print *, "test_lib_math_list_cartesian_operator_array_c_mul_array_cmplx:"
                    do i=1, d
                        buffer = abs(value%coordinate(i)%y - ground_truth_value%coordinate(i)%y)
                        if (.not. evaluate(buffer, i, 'y')) rv = .false.

                        buffer = abs(value%coordinate(i)%x - ground_truth_value%coordinate(i)%x)
                        if (.not. evaluate(buffer, i, "x")) rv = .false.

                        buffer = abs(value%coordinate(i)%z - ground_truth_value%coordinate(i)%z)
                        if (.not. evaluate(buffer, i, "z")) rv = .false.
                    end do

                end function test_lib_math_list_cartesian_operator_array_c_mul_array_cmplx

                function test_lib_math_list_cartesian_operator_real_mul_array_cmplx() result (rv)
                    implicit none
                    ! dummy
                    logical :: rv

                    ! parameter
                    integer(kind=4), parameter :: d = 2

                    ! auxiliary
                    integer(kind=4) :: i
                    real(kind=8) :: buffer
                    real(kind=lib_math_type_kind) :: lhs
                    type (list_cartesian_coordinate_cmplx_type) :: rhs
                    type (list_cartesian_coordinate_cmplx_type) :: value

                    type (list_cartesian_coordinate_cmplx_type) :: ground_truth_value

                    allocate(rhs%coordinate(d))
                    allocate(ground_truth_value%coordinate(d))

                    ! lhs
                    lhs = 4

                    ! rhs
                    rhs%coordinate(1) = get_test_values_cartesian_1()
                    rhs%coordinate(2) = get_test_values_cartesian_2()

                    ! ground truth
                    ground_truth_value%coordinate(1)%y = cmplx(4,8, kind=8)
                    ground_truth_value%coordinate(1)%x = cmplx(12,16, kind=8)
                    ground_truth_value%coordinate(1)%z = cmplx(20,24, kind=8)

                    ground_truth_value%coordinate(2)%y = cmplx(24,20, kind=8)
                    ground_truth_value%coordinate(2)%x = cmplx(16,12, kind=8)
                    ground_truth_value%coordinate(2)%z = cmplx(8,4, kind=8)

                    value = lib_math_list_cartesian_operator_real_mul_array_cmplx(lhs, rhs)

                    rv = .true.
                    print *, "test_lib_math_list_cartesian_operator_real_mul_array_cmplx:"
                    do i=1, d
                        buffer = abs(value%coordinate(i)%y - ground_truth_value%coordinate(i)%y)
                        if (.not. evaluate(buffer, i, 'y')) rv = .false.

                        buffer = abs(value%coordinate(i)%x - ground_truth_value%coordinate(i)%x)
                        if (.not. evaluate(buffer, i, "x")) rv = .false.

                        buffer = abs(value%coordinate(i)%z - ground_truth_value%coordinate(i)%z)
                        if (.not. evaluate(buffer, i, "z")) rv = .false.
                    end do

                end function test_lib_math_list_cartesian_operator_real_mul_array_cmplx

                function test_lib_math_list_cartesian_operator_cmplx_mul_array_cmplx() result (rv)
                    implicit none
                    ! dummy
                    logical :: rv

                    ! parameter
                    integer(kind=4), parameter :: d = 2

                    ! auxiliary
                    integer(kind=4) :: i
                    real(kind=8) :: buffer
                    complex(kind=lib_math_type_kind) :: lhs
                    type (list_cartesian_coordinate_cmplx_type) :: rhs
                    type (list_cartesian_coordinate_cmplx_type) :: value

                    type (list_cartesian_coordinate_cmplx_type) :: ground_truth_value

                    allocate(rhs%coordinate(d))
                    allocate(ground_truth_value%coordinate(d))

                    ! lhs
                    lhs = cmplx(0, 4, kind=8)

                    ! rhs
                    rhs%coordinate(1) = get_test_values_cartesian_1()
                    rhs%coordinate(2) = get_test_values_cartesian_2()

                    ! ground truth
                    ground_truth_value%coordinate(1)%y = cmplx(-8,4, kind=8)
                    ground_truth_value%coordinate(1)%x = cmplx(-16,12, kind=8)
                    ground_truth_value%coordinate(1)%z = cmplx(-24,20, kind=8)

                    ground_truth_value%coordinate(2)%y = cmplx(-20,24, kind=8)
                    ground_truth_value%coordinate(2)%x = cmplx(-12,16, kind=8)
                    ground_truth_value%coordinate(2)%z = cmplx(-4,8, kind=8)

                    value = lib_math_list_cartesian_operator_cmplx_mul_array_cmplx(lhs, rhs)

                    rv = .true.
                    print *, "test_lib_math_list_cartesian_operator_cmplx_mul_array_cmplx:"
                    do i=1, d
                        buffer = abs(value%coordinate(i)%y - ground_truth_value%coordinate(i)%y)
                        if (.not. evaluate(buffer, i, 'y')) rv = .false.

                        buffer = abs(value%coordinate(i)%x - ground_truth_value%coordinate(i)%x)
                        if (.not. evaluate(buffer, i, "x")) rv = .false.

                        buffer = abs(value%coordinate(i)%z - ground_truth_value%coordinate(i)%z)
                        if (.not. evaluate(buffer, i, "z")) rv = .false.
                    end do

                end function test_lib_math_list_cartesian_operator_cmplx_mul_array_cmplx

                function test_lib_math_list_cartesian_operator_c_array_divide_by_r_array() result (rv)
                    implicit none
                    ! dummy
                    logical :: rv

                    ! parameter
                    integer(kind=4), parameter :: d = 2

                    ! auxiliary
                    integer(kind=4) :: i
                    real(kind=8) :: buffer
                    type (list_cartesian_coordinate_cmplx_type) :: lhs
                    real(kind=lib_math_type_kind), dimension(d) :: rhs
                    type (list_cartesian_coordinate_cmplx_type) :: value

                    type (list_cartesian_coordinate_cmplx_type) :: ground_truth_value

                    allocate(lhs%coordinate(d))
                    allocate(ground_truth_value%coordinate(d))

                    ! lhs
                    lhs%coordinate(1) = get_test_values_cartesian_1()
                    lhs%coordinate(2) = get_test_values_cartesian_2()

                    ! rhs
                    rhs(1) = 2.0
                    rhs(2) = 4.0

                    ! ground truth
                    ground_truth_value%coordinate(1)%y = cmplx(0.5,1, kind=8)
                    ground_truth_value%coordinate(1)%x = cmplx(1.5,2, kind=8)
                    ground_truth_value%coordinate(1)%z = cmplx(2.5,3, kind=8)

                    ground_truth_value%coordinate(2)%y = cmplx(1.5,1.25, kind=8)
                    ground_truth_value%coordinate(2)%x = cmplx(1,0.75, kind=8)
                    ground_truth_value%coordinate(2)%z = cmplx(0.5,0.25, kind=8)

                    value = lib_math_list_cartesian_operator_array_divide_by_real_array(lhs, rhs)

                    rv = .true.
                    print *, "test_lib_math_list_cartesian_operator_c_array_divide_by_r_array:"
                    do i=1, d
                        buffer = abs(value%coordinate(i)%y - ground_truth_value%coordinate(i)%y)
                        if (.not. evaluate(buffer, i, 'y')) rv = .false.

                        buffer = abs(value%coordinate(i)%x - ground_truth_value%coordinate(i)%x)
                        if (.not. evaluate(buffer, i, "x")) rv = .false.

                        buffer = abs(value%coordinate(i)%z - ground_truth_value%coordinate(i)%z)
                        if (.not. evaluate(buffer, i, "z")) rv = .false.
                    end do

                end function test_lib_math_list_cartesian_operator_c_array_divide_by_r_array

                function test_lib_math_list_cartesian_operator_c_array_divide_by_c_array() result (rv)
                    implicit none
                    ! dummy
                    logical :: rv

                    ! parameter
                    integer(kind=4), parameter :: d = 2

                    ! auxiliary
                    integer(kind=4) :: i
                    real(kind=8) :: buffer
                    type (list_cartesian_coordinate_cmplx_type) :: lhs
                    complex(kind=lib_math_type_kind), dimension(d) :: rhs
                    type (list_cartesian_coordinate_cmplx_type) :: value

                    type (list_cartesian_coordinate_cmplx_type) :: ground_truth_value

                    allocate(lhs%coordinate(d))
                    allocate(ground_truth_value%coordinate(d))

                    ! lhs
                    lhs%coordinate(1) = get_test_values_cartesian_1()
                    lhs%coordinate(2) = get_test_values_cartesian_2()

                    ! rhs
                    rhs(1) = cmplx(0, 2, kind=8)
                    rhs(2) = cmplx(0, 4, kind=8)

                    ! ground truth
                    ground_truth_value%coordinate(1)%y = cmplx(1, -0.5, kind=8)
                    ground_truth_value%coordinate(1)%x = cmplx(2, -1.5, kind=8)
                    ground_truth_value%coordinate(1)%z = cmplx(3, -2.5, kind=8)

                    ground_truth_value%coordinate(2)%y = cmplx(1.25, -1.5, kind=8)
                    ground_truth_value%coordinate(2)%x = cmplx(0.75, -1, kind=8)
                    ground_truth_value%coordinate(2)%z = cmplx(0.25, -0.5, kind=8)

                    value = lib_math_list_cartesian_operator_c_array_divide_by_c_array(lhs, rhs)

                    rv = .true.
                    print *, "test_lib_math_list_cartesian_operator_c_array_divide_by_c_array:"
                    do i=1, d
                        buffer = abs(value%coordinate(i)%y - ground_truth_value%coordinate(i)%y)
                        if (.not. evaluate(buffer, i, 'y')) rv = .false.

                        buffer = abs(value%coordinate(i)%x - ground_truth_value%coordinate(i)%x)
                        if (.not. evaluate(buffer, i, "x")) rv = .false.

                        buffer = abs(value%coordinate(i)%z - ground_truth_value%coordinate(i)%z)
                        if (.not. evaluate(buffer, i, "z")) rv = .false.
                    end do

                end function test_lib_math_list_cartesian_operator_c_array_divide_by_c_array

                function test_lib_math_list_cartesian_operator_array_divide_by_real() result (rv)
                    implicit none
                    ! dummy
                    logical :: rv

                    ! parameter
                    integer(kind=4), parameter :: d = 2

                    ! auxiliary
                    integer(kind=4) :: i
                    real(kind=8) :: buffer
                    type (list_cartesian_coordinate_cmplx_type) :: lhs
                    real(kind=lib_math_type_kind) :: rhs
                    type (list_cartesian_coordinate_cmplx_type) :: value

                    type (list_cartesian_coordinate_cmplx_type) :: ground_truth_value

                    allocate(lhs%coordinate(d))
                    allocate(ground_truth_value%coordinate(d))

                    ! lhs
                    lhs%coordinate(1) = get_test_values_cartesian_1()
                    lhs%coordinate(2) = get_test_values_cartesian_2()

                    ! rhs
                    rhs = 2

                    ! ground truth
                    ground_truth_value%coordinate(1)%y = cmplx(0.5, 1, kind=8)
                    ground_truth_value%coordinate(1)%x = cmplx(1.5, 2, kind=8)
                    ground_truth_value%coordinate(1)%z = cmplx(2.5, 3, kind=8)

                    ground_truth_value%coordinate(2)%y = cmplx(3, 2.5, kind=8)
                    ground_truth_value%coordinate(2)%x = cmplx(2, 1.5, kind=8)
                    ground_truth_value%coordinate(2)%z = cmplx(1, 0.5, kind=8)

                    value = lib_math_list_cartesian_operator_array_divide_by_real(lhs, rhs)

                    rv = .true.
                    print *, "test_lib_math_list_cartesian_operator_array_divide_by_real:"
                    do i=1, d
                        buffer = abs(value%coordinate(i)%y - ground_truth_value%coordinate(i)%y)
                        if (.not. evaluate(buffer, i, 'y')) rv = .false.

                        buffer = abs(value%coordinate(i)%x - ground_truth_value%coordinate(i)%x)
                        if (.not. evaluate(buffer, i, "x")) rv = .false.

                        buffer = abs(value%coordinate(i)%z - ground_truth_value%coordinate(i)%z)
                        if (.not. evaluate(buffer, i, "z")) rv = .false.
                    end do

                end function test_lib_math_list_cartesian_operator_array_divide_by_real

                function test_lib_math_list_cartesian_operator_array_divide_by_cmplx() result (rv)
                    implicit none
                    ! dummy
                    logical :: rv

                    ! parameter
                    integer(kind=4), parameter :: d = 2

                    ! auxiliary
                    integer(kind=4) :: i
                    real(kind=8) :: buffer
                    type (list_cartesian_coordinate_cmplx_type) :: lhs
                    complex(kind=lib_math_type_kind) :: rhs
                    type (list_cartesian_coordinate_cmplx_type) :: value

                    type (list_cartesian_coordinate_cmplx_type) :: ground_truth_value

                    allocate(lhs%coordinate(d))
                    allocate(ground_truth_value%coordinate(d))

                    ! lhs
                    lhs%coordinate(1) = get_test_values_cartesian_1()
                    lhs%coordinate(2) = get_test_values_cartesian_2()

                    ! rhs
                    rhs = cmplx(0, 2, kind=8)

                    ! ground truth
                    ground_truth_value%coordinate(1)%y = cmplx(1, -0.5, kind=8)
                    ground_truth_value%coordinate(1)%x = cmplx(2, -1.5, kind=8)
                    ground_truth_value%coordinate(1)%z = cmplx(3, -2.5, kind=8)

                    ground_truth_value%coordinate(2)%y = cmplx(2.5, -3, kind=8)
                    ground_truth_value%coordinate(2)%x = cmplx(1.5, -2, kind=8)
                    ground_truth_value%coordinate(2)%z = cmplx(0.5, -1, kind=8)

                    value = lib_math_list_cartesian_operator_array_divide_by_cmplx(lhs, rhs)

                    rv = .true.
                    print *, "test_lib_math_list_cartesian_operator_array_divide_by_cmplx:"
                    do i=1, d
                        buffer = abs(value%coordinate(i)%y - ground_truth_value%coordinate(i)%y)
                        if (.not. evaluate(buffer, i, "y  ")) rv = .false.

                        buffer = abs(value%coordinate(i)%x - ground_truth_value%coordinate(i)%x)
                        if (.not. evaluate(buffer, i, "x  ")) rv = .false.

                        buffer = abs(value%coordinate(i)%z - ground_truth_value%coordinate(i)%z)
                        if (.not. evaluate(buffer, i, "z")) rv = .false.
                    end do

                end function test_lib_math_list_cartesian_operator_array_divide_by_cmplx

                function test_lib_math_cartesian_point_to_spherical_point() result (rv)
                    implicit none
                    ! dummy
                    logical :: rv

                    ! parameter
                    integer, parameter :: d = 1

                    ! auxiliary
                    integer :: i
                    type(cartesian_coordinate_real_type) :: rhs
                    type(spherical_coordinate_real_type), dimension(d) :: lhs
                    type(spherical_coordinate_real_type), dimension(d) :: ground_truth_lhs

                    real(kind=lib_math_type_kind) :: buffer

                    rhs%x = 1
                    rhs%y = 1
                    rhs%z = 1

                    ground_truth_lhs(1)%rho = sqrt(3.0_8)
                    ground_truth_lhs(1)%theta = 0.9553166181245092_8
                    ground_truth_lhs(1)%phi = 0.25 * PI

                    i=1

                    call lib_math_cartesian_point_to_spherical_point(lhs(i), rhs)

                    rv = .true.
                    print *, "test_lib_math_cartesian_point_to_spherical_point:"

                    buffer = abs(lhs(i)%rho - ground_truth_lhs(i)%rho)
                    if (.not. evaluate(buffer, i, "rho  ")) rv = .false.

                    buffer = abs(lhs(i)%theta - ground_truth_lhs(i)%theta)
                    if (.not. evaluate(buffer, i, "theta")) rv = .false.

                    buffer = abs(lhs(i)%phi - ground_truth_lhs(i)%phi)
                    if (.not. evaluate(buffer, i, "phi  ")) rv = .false.

                end function test_lib_math_cartesian_point_to_spherical_point

                function test_lib_math_spherical_point_to_cartesian_point() result (rv)
                    implicit none
                    ! dummy
                    logical :: rv

                    ! parameter
                    integer, parameter :: d = 1

                    ! auxiliary
                    integer :: i
                    type(spherical_coordinate_real_type) :: rhs
                    type(cartesian_coordinate_real_type), dimension(d) :: lhs
                    type(cartesian_coordinate_real_type), dimension(d) :: ground_truth_lhs

                    real(kind=lib_math_type_kind) :: buffer

                    rhs%rho = 1
                    rhs%theta = 0.25 * PI
                    rhs%phi = 0.25 * PI

                    ground_truth_lhs(1)%x = 0.5_8
                    ground_truth_lhs(1)%y = 0.5_8
                    ground_truth_lhs(1)%z = 1.0_8 / sqrt(2.0_8)

                    i=1

                    call lib_math_spherical_point_to_cartesian_point(lhs(i), rhs)

                    rv = .true.
                    print *, "test_lib_math_spherical_point_to_cartesian_point:"

                    buffer = abs(lhs(i)%x - ground_truth_lhs(i)%x)
                    if (.not. evaluate(buffer, i, "x")) rv = .false.

                    buffer = abs(lhs(i)%y - ground_truth_lhs(i)%y)
                    if (.not. evaluate(buffer, i, "y")) rv = .false.

                    buffer = abs(lhs(i)%z - ground_truth_lhs(i)%z)
                    if (.not. evaluate(buffer, i, "z")) rv = .false.

                end function test_lib_math_spherical_point_to_cartesian_point

                function test_lib_math_spherical_components_to_cartesian_components_c_a() result (rv)
                    implicit none
                    ! dummy
                    logical :: rv

                    ! parameter
                    integer, parameter :: d = 6

                    ! auxiliary
                    integer :: i
                    type(spherical_coordinate_cmplx_type) :: rhs
                    real(kind=lib_math_type_kind), dimension(d) :: theta
                    real(kind=lib_math_type_kind), dimension(d) :: phi

                    type(cartesian_coordinate_cmplx_type), dimension(d) :: lhs
                    type(cartesian_coordinate_cmplx_type), dimension(d) :: ground_truth_lhs

                    real(kind=lib_math_type_kind) :: buffer

                    rhs%rho = cmplx(1,0)
                    rhs%theta = cmplx(0,0)
                    rhs%phi = cmplx(0,0)

                    theta(1) = 0.25 * PI
                    phi(1) = 0

                    theta(2) = 0.25 * PI
                    phi(2) = 0.25 * PI

                    theta(3) = 0.25 * PI
                    phi(3) = 0.75 * PI

                    theta(4) = 0.25 * PI
                    phi(4) = 1.25 * PI

                    theta(5) = 0.25 * PI
                    phi(5) = 1.5 * PI

                    theta(6) = 0
                    phi(6) = 0

                    ground_truth_lhs(1)%x = 1.0_8 / sqrt(2.0_8)
                    ground_truth_lhs(1)%y = cmplx(0,0)
                    ground_truth_lhs(1)%z = 1.0_8 / sqrt(2.0_8)

                    ground_truth_lhs(2)%x = 0.5_8
                    ground_truth_lhs(2)%y = 0.5_8
                    ground_truth_lhs(2)%z = 1.0_8 / sqrt(2.0_8)

                    ground_truth_lhs(3)%x = -0.5_8
                    ground_truth_lhs(3)%y = 0.5_8
                    ground_truth_lhs(3)%z = 1.0_8 / sqrt(2.0_8)

                    ground_truth_lhs(4)%x = - 0.5_8
                    ground_truth_lhs(4)%y = - 0.5_8
                    ground_truth_lhs(4)%z = 1.0_8 / sqrt(2.0_8)

                    ground_truth_lhs(5)%x = cmplx(0,0)
                    ground_truth_lhs(5)%y = - 1.0_8 / sqrt(2.0_8)
                    ground_truth_lhs(5)%z = 1.0_8 / sqrt(2.0_8)

                    ground_truth_lhs(6)%x = cmplx(0,0)
                    ground_truth_lhs(6)%y = cmplx(0,0)
                    ground_truth_lhs(6)%z = 1.0_8

                    do i=1, d
                        lhs(i) = lib_math_spherical_components_to_cartesian_components_cmplx_a(rhs, &
                                                                                               theta(i), &
                                                                                               phi(i))
                    end do

                    rv = .true.
                    print *, "test_lib_math_spherical_components_to_cartesian_components_c_a:"
                    do i=1, d
                        buffer = abs(lhs(i)%x - ground_truth_lhs(i)%x)
                        if (.not. evaluate(buffer, i, "x")) rv = .false.

                        buffer = abs(lhs(i)%y - ground_truth_lhs(i)%y)
                        if (.not. evaluate(buffer, i, "y")) rv = .false.

                        buffer = abs(lhs(i)%z - ground_truth_lhs(i)%z)
                        if (.not. evaluate(buffer, i, "z")) rv = .false.
                    end do

                end function test_lib_math_spherical_components_to_cartesian_components_c_a

                function test_lib_math_cartesian_components_to_spherical_components_c_a() result(rv)
                    implicit none
                    ! dummy
                    logical :: rv

                    ! parameter
                    integer, parameter :: d = 6

                    ! auxiliary
                    integer :: i
                    type(cartesian_coordinate_cmplx_type) :: rhs
                    real(kind=lib_math_type_kind), dimension(d) :: theta
                    real(kind=lib_math_type_kind), dimension(d) :: phi

                    type(spherical_coordinate_cmplx_type), dimension(d) :: lhs
                    type(spherical_coordinate_cmplx_type), dimension(d) :: ground_truth_lhs

                    real(kind=lib_math_type_kind) :: buffer

                    rhs%x = cmplx(0,0)
                    rhs%y = cmplx(0,0)
                    rhs%z = cmplx(1,0)

                    theta(1) = 0.25 * PI
                    phi(1) = 0

                    theta(2) = 0.25 * PI
                    phi(2) = 0.25 * PI

                    theta(3) = 0.25 * PI
                    phi(3) = 0.75 * PI

                    theta(4) = 0.25 * PI
                    phi(4) = 1.25 * PI

                    theta(5) = 0.25 * PI
                    phi(5) = 1.5 * PI

                    theta(6) = 0
                    phi(6) = 0

                    ground_truth_lhs(1)%rho = 1.0_8 / sqrt(2.0_8)
                    ground_truth_lhs(1)%theta = - 1.0_8 / sqrt(2.0_8)
                    ground_truth_lhs(1)%phi = cmplx(0,0)

                    ground_truth_lhs(2)%rho = 1.0_8 / sqrt(2.0_8)
                    ground_truth_lhs(2)%theta = - 1.0_8 / sqrt(2.0_8)
                    ground_truth_lhs(2)%phi = cmplx(0,0)

                    ground_truth_lhs(3)%rho = 1.0_8 / sqrt(2.0_8)
                    ground_truth_lhs(3)%theta = - 1.0_8 / sqrt(2.0_8)
                    ground_truth_lhs(3)%phi = cmplx(0,0)

                    ground_truth_lhs(4)%rho = 1.0_8 / sqrt(2.0_8)
                    ground_truth_lhs(4)%theta = - 1.0_8 / sqrt(2.0_8)
                    ground_truth_lhs(4)%phi = cmplx(0,0)

                    ground_truth_lhs(5)%rho = 1.0_8 / sqrt(2.0_8)
                    ground_truth_lhs(5)%theta = - 1.0_8 / sqrt(2.0_8)
                    ground_truth_lhs(5)%phi = cmplx(0,0)

                    ground_truth_lhs(6)%rho = 1.0_8
                    ground_truth_lhs(6)%theta = cmplx(0,0)
                    ground_truth_lhs(6)%phi = cmplx(0,0)

                    do i=1, d
                        lhs(i) = lib_math_cartesian_components_to_spherical_components_cmplx_a(rhs, &
                                                                                               theta(i), &
                                                                                               phi(i))
                    end do

                    rv = .true.
                    print *, "test_lib_math_cartesian_components_to_spherical_components_c_a:"
                    do i=1, d
                        buffer = abs(lhs(i)%rho - ground_truth_lhs(i)%rho)
                        if (.not. evaluate(buffer, i, "rho  ")) rv = .false.

                        buffer = abs(lhs(i)%theta - ground_truth_lhs(i)%theta)
                        if (.not. evaluate(buffer, i, "theta")) rv = .false.

                        buffer = abs(lhs(i)%phi - ground_truth_lhs(i)%phi)
                        if (.not. evaluate(buffer, i, "phi  ")) rv = .false.
                    end do

                end function test_lib_math_cartesian_components_to_spherical_components_c_a

                function test_lib_math_cartesian_cross_product_real_list() result(rv)
                    implicit none
                    ! dummy
                    logical :: rv

                    integer, parameter :: d = 4

                    ! auxiliary
                    integer :: i
                    type(cartesian_coordinate_real_type), dimension(d) :: lhs
                    type(cartesian_coordinate_real_type), dimensioN(d) :: rhs
                    type(cartesian_coordinate_real_type), dimension(d) :: res

                    type(cartesian_coordinate_real_type), dimension(d) :: ground_truth_res

                    real(kind=lib_math_type_kind) :: buffer

                    lhs(1)%x = 1
                    lhs(1)%y = 0
                    lhs(1)%z = 0
                    lhs(2)%x = 0
                    lhs(2)%y = 1
                    lhs(2)%z = 0
                    lhs(3)%x = 0
                    lhs(3)%y = 0
                    lhs(3)%z = 1
                    lhs(4)%x = 4
                    lhs(4)%y = 45.5
                    lhs(4)%z = 50

                    rhs(1)%x = 0
                    rhs(1)%y = 1
                    rhs(1)%z = 0
                    rhs(2)%x = 0
                    rhs(2)%y = 0
                    rhs(2)%z = 1
                    rhs(3)%x = 1
                    rhs(3)%y = 0
                    rhs(3)%z = 0
                    rhs(4)%x = 5
                    rhs(4)%y = 1
                    rhs(4)%z = 0

                    ground_truth_res(1)%x = 0
                    ground_truth_res(1)%y = 0
                    ground_truth_res(1)%z = 1
                    ground_truth_res(2)%x = 1
                    ground_truth_res(2)%y = 0
                    ground_truth_res(2)%z = 0
                    ground_truth_res(3)%x = 0
                    ground_truth_res(3)%y = 1
                    ground_truth_res(3)%z = 0
                    ground_truth_res(4)%x = -50
                    ground_truth_res(4)%y = 250
                    ground_truth_res(4)%z = -223.5

                    res = lib_math_cartesian_cross_product_real_list(lhs, rhs)

                    rv = .true.
                    print *, "test_lib_math_cartesian_cross_product_real_list:"
                    do i=1, d
                        buffer = abs(res(i)%x - ground_truth_res(i)%x)
                        if (.not. evaluate(buffer, i, "x")) rv = .false.

                        buffer = abs(res(i)%y - ground_truth_res(i)%y)
                        if (.not. evaluate(buffer, i, "y")) rv = .false.

                        buffer = abs(res(i)%z - ground_truth_res(i)%z)
                        if (.not. evaluate(buffer, i, "z")) rv = .false.
                    end do

                end function test_lib_math_cartesian_cross_product_real_list

                function test_lib_math_cartesian_cross_product_cmplx_list() result(rv)
                    implicit none
                    ! dummy
                    logical :: rv

                    integer, parameter :: d = 4

                    ! auxiliary
                    integer :: i
                    type(cartesian_coordinate_cmplx_type), dimension(d) :: lhs
                    type(cartesian_coordinate_cmplx_type), dimensioN(d) :: rhs
                    type(cartesian_coordinate_cmplx_type), dimension(d) :: res

                    type(cartesian_coordinate_cmplx_type), dimension(d) :: ground_truth_res

                    real(kind=lib_math_type_kind) :: buffer

                    lhs(1)%x = (1, 1)
                    lhs(1)%y = 0
                    lhs(1)%z = 0
                    lhs(2)%x = 0
                    lhs(2)%y = (1, 1)
                    lhs(2)%z = 0
                    lhs(3)%x = 0
                    lhs(3)%y = 0
                    lhs(3)%z = (1, 1)
                    lhs(4)%x = 4
                    lhs(4)%y = 45.5
                    lhs(4)%z = 50

                    rhs(1)%x = 0
                    rhs(1)%y = (1, 1)
                    rhs(1)%z = 0
                    rhs(2)%x = 0
                    rhs(2)%y = 0
                    rhs(2)%z = (1, 1)
                    rhs(3)%x = (1, 1)
                    rhs(3)%y = 0
                    rhs(3)%z = 0
                    rhs(4)%x = 5
                    rhs(4)%y = 1
                    rhs(4)%z = 0

                    ground_truth_res(1)%x = 0
                    ground_truth_res(1)%y = 0
                    ground_truth_res(1)%z = (0, 2)
                    ground_truth_res(2)%x = (0, 2)
                    ground_truth_res(2)%y = 0
                    ground_truth_res(2)%z = 0
                    ground_truth_res(3)%x = 0
                    ground_truth_res(3)%y = (0, 2)
                    ground_truth_res(3)%z = 0
                    ground_truth_res(4)%x = -50
                    ground_truth_res(4)%y = 250
                    ground_truth_res(4)%z = -223.5

                    res = lib_math_cartesian_cross_product_cmplx_list(lhs, rhs)

                    rv = .true.
                    print *, "test_lib_math_cartesian_cross_product_cmplx_list:"
                    do i=1, d
                        buffer = abs(res(i)%x - ground_truth_res(i)%x)
                        if (.not. evaluate(buffer, i, "x")) rv = .false.

                        buffer = abs(res(i)%y - ground_truth_res(i)%y)
                        if (.not. evaluate(buffer, i, "y")) rv = .false.

                        buffer = abs(res(i)%z - ground_truth_res(i)%z)
                        if (.not. evaluate(buffer, i, "z")) rv = .false.
                    end do

                end function test_lib_math_cartesian_cross_product_cmplx_list

                function test_lib_math_cartesian_dot_product_real_list() result(rv)
                    implicit none
                    ! dummy
                    logical :: rv

                    integer, parameter :: d = 4

                    ! auxiliary
                    integer :: i
                    type(cartesian_coordinate_real_type), dimension(d) :: lhs
                    type(cartesian_coordinate_real_type), dimensioN(d) :: rhs
                    real(kind=lib_math_type_kind), dimension(d) :: res

                    real(kind=lib_math_type_kind), dimension(d) :: ground_truth_res

                    real(kind=lib_math_type_kind) :: buffer

                    lhs(1)%x = 1
                    lhs(1)%y = 0
                    lhs(1)%z = 0
                    lhs(2)%x = 0
                    lhs(2)%y = 1
                    lhs(2)%z = 0
                    lhs(3)%x = 0
                    lhs(3)%y = 0
                    lhs(3)%z = 1
                    lhs(4)%x = 4
                    lhs(4)%y = 45.5
                    lhs(4)%z = 50

                    rhs(1)%x = 0
                    rhs(1)%y = 1
                    rhs(1)%z = 0
                    rhs(2)%x = 0
                    rhs(2)%y = 0
                    rhs(2)%z = 1
                    rhs(3)%x = 1
                    rhs(3)%y = 0
                    rhs(3)%z = 0
                    rhs(4)%x = 5
                    rhs(4)%y = 1
                    rhs(4)%z = 0

                    ground_truth_res(1) = 0
                    ground_truth_res(2) = 0
                    ground_truth_res(3) = 0
                    ground_truth_res(4) = 65.5

                    res = lib_math_cartesian_dot_product_real_list(lhs, rhs)

                    rv = .true.
                    print *, "test_lib_math_cartesian_dot_product_real_list:"
                    do i=1, d
                        buffer = abs(res(i) - ground_truth_res(i))
                        if (.not. evaluate(buffer, i)) rv = .false.
                    end do

                end function test_lib_math_cartesian_dot_product_real_list

                function test_lib_math_cartesian_dot_product_cmplx_list() result(rv)
                    implicit none
                    ! dummy
                    logical :: rv

                    integer, parameter :: d = 4

                    ! auxiliary
                    integer :: i
                    type(cartesian_coordinate_cmplx_type), dimension(d) :: lhs
                    type(cartesian_coordinate_cmplx_type), dimensioN(d) :: rhs
                    complex(kind=lib_math_type_kind), dimension(d) :: res

                    complex(kind=lib_math_type_kind), dimension(d) :: ground_truth_res

                    real(kind=lib_math_type_kind) :: buffer

                    lhs(1)%x = 1
                    lhs(1)%y = 0
                    lhs(1)%z = 0
                    lhs(2)%x = 0
                    lhs(2)%y = 1
                    lhs(2)%z = 0
                    lhs(3)%x = 0
                    lhs(3)%y = 0
                    lhs(3)%z = 1
                    lhs(4)%x = 4
                    lhs(4)%y = 45.5
                    lhs(4)%z = 50

                    rhs(1)%x = 0
                    rhs(1)%y = 1
                    rhs(1)%z = 0
                    rhs(2)%x = 0
                    rhs(2)%y = 0
                    rhs(2)%z = 1
                    rhs(3)%x = 1
                    rhs(3)%y = 0
                    rhs(3)%z = 0
                    rhs(4)%x = 5
                    rhs(4)%y = 1
                    rhs(4)%z = 0

                    ground_truth_res(1) = 0
                    ground_truth_res(2) = 0
                    ground_truth_res(3) = 0
                    ground_truth_res(4) = 65.5

                    res = lib_math_cartesian_dot_product_cmplx_list(lhs, rhs)

                    rv = .true.
                    print *, "test_lib_math_cartesian_dot_product_cmplx_list:"
                    do i=1, d
                        buffer = abs(res(i) - ground_truth_res(i))
                        if (.not. evaluate(buffer, i)) rv = .false.
                    end do

                end function test_lib_math_cartesian_dot_product_cmplx_list

        end function lib_math_type_operator_test_functions

end module lib_math_type_operator
