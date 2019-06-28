module lib_math_type_operator
    use lib_math_type
    implicit none

    private

    ! --- public ---
    public :: operator (+)
    public :: operator (-)
    public :: operator (*)
    public :: operator (/)

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

        function lib_math_type_operator_test_functions() result (rv)
            implicit none
            ! dummy
            integer :: rv

            rv = 0

            if (.not. test_lib_math_list_spherical_operator_add_array()) then
                rv = rv + 1
            end if

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

        end function lib_math_type_operator_test_functions

end module lib_math_type_operator
