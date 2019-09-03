module lib_math_auxiliaries
    implicit none

    private

    public :: isinf

    public :: test_lib_math_auxiliaries

    interface isinf
        module procedure lib_math_real_4_isinf
        module procedure lib_math_real_8_isinf
    end interface

    contains

    function lib_math_real_4_isinf(value) result (rv)
        implicit none
        ! dummy
        real(kind=4) :: value

        logical :: rv

        ! auxiliary
        real(kind=4) :: inf

        inf = huge(inf)

        if (value .gt. inf &
            .or. value .lt. -inf) then
            rv = .true.
        else
            rv = .false.
        end if
    end function

    function lib_math_real_8_isinf(value) result (rv)
        implicit none
        ! dummy
        real(kind=8) :: value

        logical :: rv

        ! auxiliary
        real(kind=8) :: inf

        inf = huge(inf)

        if (value .gt. inf &
            .or. value .lt. -inf) then
            rv = .true.
        else
            rv = .false.
        end if
    end function

    function test_lib_math_auxiliaries() result (rv)
        implicit none
        ! dummy
        integer :: rv

        rv = 0

        if (.not. test_lib_math_real_4_isinf()) then
            rv = rv + 1
        end if
        if (.not. test_lib_math_real_8_isinf()) then
            rv = rv + 1
        end if


        print *, "-------------test_lib_math_auxiliaries----------------"
        if (rv == 0) then
            print *, "test_lib_math_auxiliaries tests: OK"
        else
            print *, rv,"test_lib_math_auxiliaries test(s) FAILED"
        end if
        print *, "------------------------------------------------------"

        contains

            function test_lib_math_real_4_isinf() result (rv)
                implicit none
                ! dummy
                logical :: rv

                ! parameter
                integer, parameter :: d = 3

                ! auxiliary
                integer :: i
                real(kind=4), dimension(d) :: value
                logical, dimension(d) :: res
                logical, dimension(d) :: ground_truth_res

                value(1) = huge(value)
                value(1) = value(1) + value(1)
                value(2) = -value(1)
                value(3) = huge(value)

                ground_truth_res = (/ .true., .true., .false. /)

                do i=1, d
                    res(i) = lib_math_real_4_isinf(value(i))
                end do

                rv = .true.
                print *, "test_lib_math_real_4_isinf:"
                do i=1, d
                    if (res(i) .neqv. ground_truth_res(i)) then
                        rv = .false.
                        print *, "  ", i, ": FAILED"
                    else
                        print *, "  ", i, ": OK"
                    end if
                end do

            end function test_lib_math_real_4_isinf

            function test_lib_math_real_8_isinf() result (rv)
                implicit none
                ! dummy
                logical :: rv

                ! parameter
                integer, parameter :: d = 3

                ! auxiliary
                integer :: i
                real(kind=8), dimension(d) :: value
                logical, dimension(d) :: res
                logical, dimension(d) :: ground_truth_res

                value(1) = huge(value)
                value(1) = value(1) + value(1)
                value(2) = -value(1)
                value(3) = huge(value)

                ground_truth_res = (/ .true., .true., .false. /)

                do i=1, d
                    res(i) = lib_math_real_8_isinf(value(i))
                end do

                rv = .true.
                print *, "test_lib_math_real_8_isinf:"
                do i=1, d
                    if (res(i) .neqv. ground_truth_res(i)) then
                        rv = .false.
                        print *, "  ", i, ": FAILED"
                    else
                        print *, "  ", i, ": OK"
                    end if
                end do

            end function test_lib_math_real_8_isinf

    end function test_lib_math_auxiliaries

end module lib_math_auxiliaries
