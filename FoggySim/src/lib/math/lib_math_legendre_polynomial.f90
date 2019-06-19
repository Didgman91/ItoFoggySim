module lib_math_legendre_polynomial
    use libmath
    implicit none

    private

    public :: lib_math_legendre_polynomial_test_functions

    contains

        subroutine lib_math_associated_legendre_polynomial(x, m, n, pm, pd)
            implicit none
            ! dummy
            double precision, intent(in) :: x
            integer(kind=4), intent(in) :: m
            integer(kind=4), intent(in) :: n
            double precision, intent(inout) :: pm(0:n)
            double precision, intent(inout) :: pd(0:n)

            call LPMNS(m, n, x, pm, pd)

            pm = -pm
            pd = -pd

        end subroutine


        function lib_math_legendre_polynomial_test_functions() result (rv)
            implicit none
            ! dummy
            integer :: rv

            ! auxiliaray
            double precision, parameter :: ground_truth_e = 10.0_8**(-6.0_8)

            rv = 0

            if (.not. test_lib_math_associated_legendre_polynomial()) then
                rv = rv + 1
            end if

            contains

            function test_lib_math_associated_legendre_polynomial() result (rv)
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

                call lib_math_associated_legendre_polynomial(x, m, n, pm, pd)


                rv = .true.
                print *, "test_lib_math_associated_legendre_polynomial:"
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

            end function test_lib_math_associated_legendre_polynomial

        end function

end module lib_math_legendre_polynomial
