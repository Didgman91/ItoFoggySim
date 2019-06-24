module lib_math_legendre_polynomial
    use libmath
    implicit none

    private

    public :: lib_math_legendre_polynomial_test_functions

    contains

        ! calculates the associated Legendre polynomial
        !
        ! Argument
        ! ----
        !   x: double precision
        !       input value
        !   m: integer
        !       order of the polynomial, >= 0
        !   n: integer
        !       degree of the polynomial, >= 0
        !   condon_shortly_phase: boolean, optional(std: false)
        !       true: with Condon–Shortley phase
        !       false: without Condon–Shortley phase
        !
        ! Results
        ! ----
        !   pm: double precision, dimension(0:n)
        !       result of the associated Legendre polynomial
        !   pd: double precision, dimension(0:n)
        !       result of the deriviative of the associated Legendre polynomial
        !
        subroutine lib_math_associated_legendre_polynomial(x, m, n, pm, pd, condon_shortley_phase)
            implicit none
            ! dummy
            double precision, intent(in) :: x
            integer(kind=4), intent(in) :: m
            integer(kind=4), intent(in) :: n
            double precision, intent(inout) :: pm(0:n)
            double precision, intent(inout) :: pd(0:n)
            logical, optional :: condon_shortley_phase

            call LPMNS(m, n, x, pm, pd)

            if (present(condon_shortley_phase)) then
                if (condon_shortley_phase) then
                    if (IAND(m, 1) .eq. 1) then
                        pm = -pm
                        pd = -pd
                    end if
!                    pm = (-1.0_8)**m * pm
!                    pd = (-1.0_8)**m * pd
                end if
            end if

        end subroutine lib_math_associated_legendre_polynomial


        function lib_math_legendre_polynomial_test_functions() result (rv)
            implicit none
            ! dummy
            integer :: rv

            ! auxiliaray
            double precision, parameter :: ground_truth_e = 10.0_8**(-6.0_8)

            rv = 0

            if (.not. test_lib_math_associated_legendre_polynomial_m1()) then
                rv = rv + 1
            end if
            if (.not. test_lib_math_associated_legendre_polynomial_m1_without_phase()) then
                rv = rv + 1
            end if
            if (.not. test_lib_math_associated_legendre_polynomial_m2()) then
                rv = rv + 1
            end if
            if (.not. test_lib_math_associated_legendre_polynomial_m2_without_phase()) then
                rv = rv + 1
            end if

            contains

            function test_lib_math_associated_legendre_polynomial_m1() result (rv)
                implicit none
                ! dummy
                logical :: rv

                ! auxiliary
                integer(kind=4) :: i
                double precision :: x
                integer(kind=4), parameter :: n = 5
                integer(kind=4), parameter :: m = 1
                double precision, dimension(0:n) :: pm
                double precision, dimension(0:n) :: pd

                double precision, dimension(0:n) :: ground_truth_pm
                double precision, dimension(0:n) :: ground_truth_pd

                double precision :: buffer

                ground_truth_pm(0) = 0.000000000000000_8
                ground_truth_pm(1) = -0.979795897113271_8
                ground_truth_pm(2) = -0.587877538267963_8
                ground_truth_pm(3) = 1.17575507653593_8
                ground_truth_pm(4) = 1.33252242007405_8
                ground_truth_pm(5) = -0.870058756636585_8

                ground_truth_pd(0) = -0.000000000000000_8
                ground_truth_pd(1) = 0.204124145231932_8
                ground_truth_pd(2) = -2.81691320420066_8
                ground_truth_pd(3) = -3.18433666561813_8
                ground_truth_pd(4) = 5.01328900689624_8
                ground_truth_pd(5) = 9.23457633029258_8

                x = 0.2

                call lib_math_associated_legendre_polynomial(x, m, n, pm, pd, .true.)


                rv = .true.
                print *, "test_lib_math_associated_legendre_polynomial_m1:"
                do i=0, n
                    buffer = abs(pm(i) - ground_truth_pm(i))
                    if (buffer .gt. ground_truth_e) then
                        print *, "  ", i , "difference: ", buffer, " : FAILED"
                        rv = .false.
                    else
                        print *, "  ", i, ": OK"
                    end if
                end do

                print*, "  deriviation:"
                do i=0, n
                    buffer = abs(pd(i) - ground_truth_pd(i))
                    if (buffer .gt. ground_truth_e) then
                        print *, "  ", i , "difference: ", buffer, " : FAILED"
                        rv = .false.
                    else
                        print *, "  ", i, ": OK"
                    end if
                end do

            end function test_lib_math_associated_legendre_polynomial_m1

            function test_lib_math_associated_legendre_polynomial_m1_without_phase() result (rv)
                implicit none
                ! dummy
                logical :: rv

                ! auxiliary
                integer(kind=4) :: i
                double precision :: x
                integer(kind=4), parameter :: n = 5
                integer(kind=4), parameter :: m = 1
                double precision, dimension(0:n) :: pm
                double precision, dimension(0:n) :: pd

                double precision, dimension(0:n) :: ground_truth_pm
                double precision, dimension(0:n) :: ground_truth_pd

                double precision :: buffer

                ground_truth_pm(0) = -0.000000000000000_8
                ground_truth_pm(1) = 0.979795897113271_8
                ground_truth_pm(2) = 0.587877538267963_8
                ground_truth_pm(3) = -1.17575507653593_8
                ground_truth_pm(4) = -1.33252242007405_8
                ground_truth_pm(5) = 0.870058756636585_8

                ground_truth_pd(0) = 0.000000000000000_8
                ground_truth_pd(1) = -0.204124145231932_8
                ground_truth_pd(2) = 2.81691320420066_8
                ground_truth_pd(3) = 3.18433666561813_8
                ground_truth_pd(4) = -5.01328900689624_8
                ground_truth_pd(5) = -9.23457633029258_8

                x = 0.2

!                call lib_math_associated_legendre_polynomial(x, m, n, pm, pd, .false.)
                call lib_math_associated_legendre_polynomial(x, m, n, pm, pd)


                rv = .true.
                print *, "test_lib_math_associated_legendre_polynomial_m1_without_phase:"
                do i=0, n
                    buffer = abs(pm(i) - ground_truth_pm(i))
                    if (buffer .gt. ground_truth_e) then
                        print *, "  ", i , "difference: ", buffer, " : FAILED"
                        rv = .false.
                    else
                        print *, "  ", i, ": OK"
                    end if
                end do

                print*, "  deriviation:"
                do i=0, n
                    buffer = abs(pd(i) - ground_truth_pd(i))
                    if (buffer .gt. ground_truth_e) then
                        print *, "  ", i , "difference: ", buffer, " : FAILED"
                        rv = .false.
                    else
                        print *, "  ", i, ": OK"
                    end if
                end do

            end function test_lib_math_associated_legendre_polynomial_m1_without_phase

            function test_lib_math_associated_legendre_polynomial_m2() result (rv)
                implicit none
                ! dummy
                logical :: rv

                ! auxiliary
                integer(kind=4) :: i
                double precision :: x
                integer(kind=4), parameter :: n = 5
                integer(kind=4), parameter :: m = 2
                double precision, dimension(0:n) :: pm
                double precision, dimension(0:n) :: pd

                double precision, dimension(0:n) :: ground_truth_pm
                double precision, dimension(0:n) :: ground_truth_pd

                double precision :: buffer

                ground_truth_pm(0) = 0.000000000000000_8
                ground_truth_pm(1) = -1.06581410364015D-16
                ground_truth_pm(2) = 2.88000000000000_8
                ground_truth_pm(3) = 2.88000000000000_8
                ground_truth_pm(4) = -5.18400000000000_8
                ground_truth_pm(5) = -8.87040000000000_8

                ground_truth_pd(0) = -0.000000000000000_8
                ground_truth_pd(1) = 4.81867632215780D-16
                ground_truth_pd(2) = -1.20000000000000_8
                ground_truth_pd(3) = 13.2000000000000_8
                ground_truth_pd(4) = 22.3200000000000_8
                ground_truth_pd(5) = -28.5600000000000_8

                x = 0.2

                call lib_math_associated_legendre_polynomial(x, m, n, pm, pd, .true.)


                rv = .true.
                print *, "test_lib_math_associated_legendre_polynomial_m2:"
                do i=0, n
                    buffer = abs(pm(i) - ground_truth_pm(i))
                    if (buffer .gt. ground_truth_e) then
                        print *, "  ", i , "difference: ", buffer, " : FAILED"
                        rv = .false.
                    else
                        print *, "  ", i, ": OK"
                    end if
                end do

                print*, "  deriviation:"
                do i=0, n
                    buffer = abs(pd(i) - ground_truth_pd(i))
                    if (buffer .gt. ground_truth_e) then
                        print *, "  ", i , "difference: ", buffer, " : FAILED"
                        rv = .false.
                    else
                        print *, "  ", i, ": OK"
                    end if
                end do

            end function test_lib_math_associated_legendre_polynomial_m2

            function test_lib_math_associated_legendre_polynomial_m2_without_phase() result (rv)
                implicit none
                ! dummy
                logical :: rv

                ! auxiliary
                integer(kind=4) :: i
                double precision :: x
                integer(kind=4), parameter :: n = 5
                integer(kind=4), parameter :: m = 2
                double precision, dimension(0:n) :: pm
                double precision, dimension(0:n) :: pd

                double precision, dimension(0:n) :: ground_truth_pm
                double precision, dimension(0:n) :: ground_truth_pd

                double precision :: buffer

                ground_truth_pm(0) = 0.000000000000000_8
                ground_truth_pm(1) = -1.06581410364015D-16
                ground_truth_pm(2) = 2.88000000000000_8
                ground_truth_pm(3) = 2.88000000000000_8
                ground_truth_pm(4) = -5.18400000000000_8
                ground_truth_pm(5) = -8.87040000000000_8

                ground_truth_pd(0) = -0.000000000000000_8
                ground_truth_pd(1) = 4.81867632215780D-16
                ground_truth_pd(2) = -1.20000000000000_8
                ground_truth_pd(3) = 13.2000000000000_8
                ground_truth_pd(4) = 22.3200000000000_8
                ground_truth_pd(5) = -28.5600000000000_8

                x = 0.2

                call lib_math_associated_legendre_polynomial(x, m, n, pm, pd, .false.)


                rv = .true.
                print *, "test_lib_math_associated_legendre_polynomial_m2_without_phase:"
                do i=0, n
                    buffer = abs(pm(i) - ground_truth_pm(i))
                    if (buffer .gt. ground_truth_e) then
                        print *, "  ", i , "difference: ", buffer, " : FAILED"
                        rv = .false.
                    else
                        print *, "  ", i, ": OK"
                    end if
                end do

                print*, "  deriviation:"
                do i=0, n
                    buffer = abs(pd(i) - ground_truth_pd(i))
                    if (buffer .gt. ground_truth_e) then
                        print *, "  ", i , "difference: ", buffer, " : FAILED"
                        rv = .false.
                    else
                        print *, "  ", i, ": OK"
                    end if
                end do

            end function test_lib_math_associated_legendre_polynomial_m2_without_phase

        end function

end module lib_math_legendre_polynomial
