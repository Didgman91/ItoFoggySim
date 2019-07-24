module lib_math_factorial
    implicit none

    private

    ! --- public functions ---
    public :: lib_math_factorial_get_factorial
    public :: lib_math_factorial_get_n_minus_m_divided_by_n_plus_m
    public :: lib_math_factorial_get_n_plus_m_divided_by_n_minus_m
    public :: lib_math_factorial_get_n_plus_m_divided_by_m_fac_n_minus_m

    public :: lib_math_factorial_test_functions

    ! parameter
    double precision, parameter :: ground_truth_e = 10.0_8**(-14.0_8)

    contains

        ! calculates the division of two factorials
        !
        ! formula: rv = (n-m)!  / (n+m)!
        !
        ! Arguments
        ! ----
        !   n: integer
        !       1 <= n
        !   m: integer
        !       1 <= m <= n
        !
        ! Result
        ! ----
        !   rv: double precision
        !
        function lib_math_factorial_get_n_minus_m_divided_by_n_plus_m(n,m) result (rv)
            implicit none
            ! dummy
            integer(kind=4) :: n
            integer(kind=4) :: m

            real(kind=8) :: rv

            ! auxiliary
            integer(kind=4) :: i

            if (n .eq. m) then
                ! (n - m)! / (n + m)! = 0! / (n + m)!
                ! = 1 / (n + m)!
                rv = 1.0_8 / lib_math_factorial_get_factorial(n+m)
            else
                ! 1 / ( (n-m+1) * (n-m+2) * ... * (n+m) )
                rv = 1
                do i = n-m+1, n+m
                   rv = rv * i
                end do
                rv = 1D0 / rv
            end if

        end function lib_math_factorial_get_n_minus_m_divided_by_n_plus_m

        ! calculates the division of two factorials
        !
        ! formula: rv = (n+m)!  / (n-m)!
        !
        ! Arguments
        ! ----
        !   n: integer
        !       0 <= n
        !   m: integer
        !       0 <= m <= n
        !
        ! Result
        ! ----
        !   rv: double precision
        !
        function lib_math_factorial_get_n_plus_m_divided_by_n_minus_m(n,m) result (rv)
            implicit none
            ! dummy
            integer(kind=4) :: n
            integer(kind=4) :: m

            real(kind=8) :: rv

            ! auxiliary
            integer(kind=4) :: i

            if (n .eq. m) then
                ! (n + m)! / (n - m)! = (n + m)! / 0!
                ! = (n + m)!
                rv = lib_math_factorial_get_factorial(n+m)
            else
                ! ( (n-m+1) * (n-m+2) * ... * (n+m) ) / 1
                rv = 1
                do i = n-m+1, n+m
                   rv = rv * i
                end do
            end if

        end function lib_math_factorial_get_n_plus_m_divided_by_n_minus_m

        ! calculates the division of two factorials
        !
        ! formula: rv = (n+m)!  / ( m!(n-m)! )
        !
        ! Arguments
        ! ----
        !   n: integer
        !       0 <= n
        !   m: integer
        !       0 <= m <= n
        !
        ! Result
        ! ----
        !   rv: double precision
        !
        function lib_math_factorial_get_n_plus_m_divided_by_m_fac_n_minus_m(n, m) result (rv)
            implicit none
            ! dummy
            integer(kind=4) :: n
            integer(kind=4) :: m

            real(kind=8) :: rv

            ! auxiliary

            rv = lib_math_factorial_get_n_plus_m_divided_by_n_minus_m(n, m) / lib_math_factorial_get_factorial(m)

        end function lib_math_factorial_get_n_plus_m_divided_by_m_fac_n_minus_m

        ! calculates the factorial
        !
        ! formula: rv = n!
        !
        ! Argument
        ! ----
        !   n: integer
        !       n .ge. 0
        !
        ! Retruns
        ! ----
        !   rv: real
        !       factorial of n
        function lib_math_factorial_get_factorial(n) result (rv)
            implicit none
            ! dummy
            integer(kind=4) :: n

            real(kind=8) :: rv

            ! auxiliary
            integer(kind=4) :: i

            rv = 1.0_8

            do i=1, n
                rv = real(i, kind=8) * rv
            end do

        end function lib_math_factorial_get_factorial

        function lib_math_factorial_test_functions() result (rv)
            implicit none
            ! dummy
            integer :: rv

            rv = 0

            if (.not. test_lib_math_factorial_get_factorial()) then
                rv = rv + 1
            end if
            if (.not. test_lib_math_factorial_get_n_minus_m_divided_by_n_plus_m()) then
                rv = rv + 1
            end if
            if (.not. test_lib_math_factorial_get_n_plus_m_divided_by_n_minus_m()) then
                rv = rv + 1
            end if
            if (.not. test_factorial_get_n_plus_m_divided_by_m_fac_n_minus_m()) then
                rv = rv + 1
            end if


            print *, "-------------lib_math_factorial_test_functions----------------"
            if (rv == 0) then
                print *, "lib_math_factorial_test_functions tests: OK"
            else
                print *, rv,"lib_math_factorial_test_functions test(s) FAILED"
            end if
            print *, "-----------------------------------------------------------"

            contains

                function test_lib_math_factorial_get_n_minus_m_divided_by_n_plus_m() result (rv)
                    implicit none
                    ! dummy
                    logical :: rv

                    ! parameter
                    integer(kind=1), parameter :: d = 6

                    ! auxiliary
                    integer(kind=4) :: i
                    integer(kind=4), dimension(d) :: n
                    integer(kind=4), dimension(d) :: m
                    real(kind=8), dimension(d) :: value
                    real(kind=8), dimension(d) :: ground_truth_value
                    real(kind=8) :: buffer

                    n = (/ 1, 5, 20, 30, 80, 80 /)
                    m = (/ 1, 1, 15, 30, 40, 80 /)

                    ! Values were generated with sageMath
                    !
                    ! source code:
                    ! >>> var('n')
                    ! >>> var('m')
                    ! >>> f(n,m) = factorial(n-m)/factorial(n+m)
                    ! >>>
                    ! >>> n = [1, 5, 20, 30, 80, 80]
                    ! >>> m = [1, 1, 15, 30, 40, 80]
                    ! >>>
                    ! >>> for i in range(0,6):
                    ! >>>     value = numerical_approx(f(n[i],m[i]))
                    ! >>>     print("n = {}, m = {}: {}".format(n[i], m[i], value))
                    ground_truth_value(1) = 0.500000000000000D0
                    ground_truth_value(2) = 0.0333333333333333D0
                    ground_truth_value(3) = 1.16131115503583D-38
                    ground_truth_value(4) = 1.20178049364932D-82
                    ground_truth_value(5) = 1.21969493668583D-151
                    ground_truth_value(6) = 2.12101509485313D-285

                    do i=1, d
                        value(i) = lib_math_factorial_get_n_minus_m_divided_by_n_plus_m(n(i), m(i))
                    end do

                    rv = .true.
                    print *, "test_lib_math_factorial_get_n_minus_m_divided_by_n_plus_m:"
                    do i=1, d
                        buffer = value(i) - ground_truth_value(i)
                        if (abs(buffer) .gt. ground_truth_e) then
                            print *, "  ", i , "difference: ", buffer, " : FAILED"
                            rv = .false.
                        else
                            print *, "  ", i, ": OK"
                        end if
                    end do

                end function test_lib_math_factorial_get_n_minus_m_divided_by_n_plus_m

                function test_lib_math_factorial_get_n_plus_m_divided_by_n_minus_m() result (rv)
                    implicit none
                    ! dummy
                    logical :: rv

                    ! parameter
                    integer(kind=1), parameter :: d = 7

                    ! auxiliary
                    integer(kind=4) :: i
                    integer(kind=4), dimension(d) :: n
                    integer(kind=4), dimension(d) :: m
                    real(kind=8), dimension(d) :: value
                    real(kind=8), dimension(d) :: ground_truth_value
                    real(kind=8) :: buffer

                    n = (/ 0, 1, 5, 20, 30, 80, 80 /)
                    m = (/ 0, 1, 1, 15, 30, 40, 80 /)

                    ! Values were generated with sageMath
                    !
                    ! source code:
                    ! >>> var('n')
                    ! >>> var('m')
                    ! >>> f(n,m) = factorial(n+m)/factorial(n-m)
                    ! >>>
                    ! >>> n = [0, 1, 5, 20, 30, 80, 80]
                    ! >>> m = [0, 1, 1, 15, 30, 40, 80]
                    ! >>>
                    ! >>> for i in range(0,6):
                    ! >>>     value = numerical_approx(f(n[i],m[i]))
                    ! >>>     print("n = {}, m = {}: {}".format(n[i], m[i], value))
                    ground_truth_value(1) = 1.00000000000000_8
                    ground_truth_value(2) = 2.00000000000000_8
                    ground_truth_value(3) = 30.0000000000000_8
                    ground_truth_value(4) = 8.61095663865512d37
                    ground_truth_value(5) = 8.32098711274139d81
                    ground_truth_value(6) = 8.19877142982340d150
                    ground_truth_value(7) = 4.71472363599206d284

                    do i=1, d
                        value(i) = lib_math_factorial_get_n_plus_m_divided_by_n_minus_m(n(i), m(i))
                    end do

                    rv = .true.
                    print *, "test_lib_math_factorial_get_n_plus_m_divided_by_n_minus_m:"
                    do i=1, d
                        buffer = log(value(i)) - log(ground_truth_value(i))
                        if (abs(buffer) .gt. ground_truth_e) then
                            print *, "  ", i , "difference: ", buffer, " : FAILED"
                            rv = .false.
                        else
                            print *, "  ", i, ": OK"
                        end if
                    end do

                end function test_lib_math_factorial_get_n_plus_m_divided_by_n_minus_m

                function test_factorial_get_n_plus_m_divided_by_m_fac_n_minus_m() result (rv)
                    implicit none
                    ! dummy
                    logical :: rv

                    ! parameter
                    integer(kind=1), parameter :: d = 7

                    ! auxiliary
                    integer(kind=4) :: i
                    integer(kind=4), dimension(d) :: n
                    integer(kind=4), dimension(d) :: m
                    real(kind=8), dimension(d) :: value
                    real(kind=8), dimension(d) :: ground_truth_value
                    real(kind=8) :: buffer

                    n = (/ 0, 1, 5, 20, 30, 80, 80 /)
                    m = (/ 0, 1, 1, 15, 30, 40, 80 /)

                    ! Values were generated with sageMath
                    !
                    ! source code:
                    ! >>> var('n')
                    ! >>> var('m')
                    ! >>> f(n,m) = factorial(n+m)/( factorial(m) * factorial(n-m) )
                    ! >>>
                    ! >>> n = [0, 1, 5, 20, 30, 80, 80]
                    ! >>> m = [0, 1, 1, 15, 30, 40, 80]
                    ! >>>
                    ! >>> for i in range(0,6):
                    ! >>>     value = numerical_approx(f(n[i],m[i]))
                    ! >>>     print("n = {}, m = {}: {}".format(n[i], m[i], value))
                    ground_truth_value(1) = 1.00000000000000_8
                    ground_truth_value(2) = 2.00000000000000_8
                    ground_truth_value(3) = 30.0000000000000_8
                    ground_truth_value(4) = 6.58493953033965d25
                    ground_truth_value(5) = 3.13700184745716d49
                    ground_truth_value(6) = 1.00485572438191d103
                    ground_truth_value(7) = 6.58761967824400d165

                    do i=1, d
                        value(i) = lib_math_factorial_get_n_plus_m_divided_by_m_fac_n_minus_m(n(i), m(i))
                    end do

                    rv = .true.
                    print *, "test_lib_math_factorial_get_n_plus_m_divided_by_m_fac_n_minus_m:"
                    do i=1, d
                        buffer = log(value(i)) - log(ground_truth_value(i))
                        if (abs(buffer) .gt. ground_truth_e) then
                            print *, "  ", i , "difference: ", buffer, " : FAILED"
                            rv = .false.
                        else
                            print *, "  ", i, ": OK"
                        end if
                    end do

                end function test_factorial_get_n_plus_m_divided_by_m_fac_n_minus_m

                function test_lib_math_factorial_get_factorial() result (rv)
                    implicit none
                    ! dummy
                    logical :: rv

                    ! parameter
                    integer(kind=1), parameter :: d = 5

                    ! auxiliary
                    integer(kind=4) :: i
                    integer(kind=4), dimension(d) :: n
                    real(kind=8), dimension(d) :: value
                    real(kind=8), dimension(d) :: ground_truth_value
                    real(kind=8) :: buffer

                    n = (/ 0, 1, 2, 5, 7 /)
                    ground_truth_value = (/ 1, &
                                            1, &
                                            2, &
                                            120, &
                                            5040 /)

                    do i=1, d
                        value(i) = lib_math_factorial_get_factorial(n(i))
                    end do

                    rv = .true.
                    print *, "lib_math_factorial_get_factorial:"
                    do i=1, d
                        buffer = value(i) - ground_truth_value(i)
                        if (abs(buffer) .gt. ground_truth_e) then
                            print *, "  ", i , "difference: ", buffer, " : FAILED"
                            rv = .false.
                        else
                            print *, "  ", i, ": OK"
                        end if
                    end do

                end function test_lib_math_factorial_get_factorial

        end function lib_math_factorial_test_functions

end module lib_math_factorial
