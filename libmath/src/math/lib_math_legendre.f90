module lib_math_legendre
    implicit none

    private

    ! --- public ---
    public :: lib_math_associated_legendre_polynomial
    public :: lib_math_associated_legendre_polynomial_theta
    public :: lib_math_legendre_polynomial

    public :: lib_math_legendre_test_functions

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
        !   condon_shortley_phase: boolean, optional(std: false)
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

        ! calculates the associated Legendre polynomial
        !
        ! Formula
        ! ----
        !   pm = P_ml(cos(theta)) / sin(theta)
        !
        !   pd = derivation[P_ml(cos(theta)), theta ]
        !
        !
        ! Argument
        ! ----
        !   theta: double precision
        !       input value
        !   m: integer
        !       order of the polynomial, >= 0
        !   n: integer
        !       degree of the polynomial, >= 0
        !   condon_shortley_phase: boolean, optional(std: false)
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
        subroutine lib_math_associated_legendre_polynomial_theta(theta, m, n, pm, pd, condon_shortley_phase)
            implicit none
            ! dummy
            double precision, intent(in) :: theta
            integer(kind=4), intent(in) :: m
            integer(kind=4), intent(in) :: n

            double precision, intent(inout) :: pm(0:n)
            double precision, intent(inout) :: pd(0:n)
            logical, optional :: condon_shortley_phase

            ! auxiliary
            double precision :: cos_theta
            double precision :: sin_theta

            cos_theta = cos(theta)
            sin_theta = sin(theta)

            call LPMNS(m, n, cos_theta, pm, pd)

            pm = pm / sin_theta

            if (m .eq. 1) then
                pm(0) = 0.0_8
                pm(1) = 1.0_8
            end if


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

        end subroutine lib_math_associated_legendre_polynomial_theta

        ! calculates the Legendre polynomial
        !
        ! Argument
        ! ----
        !   x: double precision
        !       input value
        !   n: integer
        !       degree of the polynomial, >= 0
        !
        ! Results
        ! ----
        !   pm: double precision, dimension(0:n)
        !       result of the Legendre polynomial
        !   pd: double precision, dimension(0:n)
        !       result of the deriviative of the Legendre polynomial
        !
        subroutine lib_math_legendre_polynomial(x, n, pn, pd)
            implicit none
            ! dummy
            double precision, intent(in) :: x
            integer(kind=4), intent(in) :: n

            double precision, dimension(0:n), intent(inout) :: pn
            double precision, dimension(0:n), intent(inout) :: pd

            call lpn ( n, x, pn, pd )

        end subroutine lib_math_legendre_polynomial

        function lib_math_legendre_test_functions() result (rv)
            implicit none
            ! dummy
            integer :: rv

            ! auxiliaray
            double precision, parameter :: ground_truth_e = 10.0_8**(-6.0_8)

            rv = 0

            ! associated Legendre Polynomial
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
            if (.not. test_lib_math_associated_legendre_polynomial_m0()) then
                rv = rv + 1
            end if
            if (.not. test_lib_math_associated_legendre_polynomial_theta_m1_without_p()) then
                rv = rv + 1
            end if

            ! Legendre Polynomial
            if (.not. test_lib_math_legendre_polynomial()) then
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

                ! Values were generated with sageMath
                !
                ! source code:
                !  >>>  var('l')
                !  >>>  var('m')
                !  >>>  var('x')
                !  >>>
                !  >>>  m=1
                !  >>>
                !  >>>  P_d(l,x) = derivative(legendre_P(l,x), x, m)
                !  >>>  P_a(l,x) = (-1)**m * (1-x**2)**(m/2) * P_d(l,x)
                !  >>>  P_ad(l,x) = derivative(P_a(l,x), x).simplify_full()
                !  >>>
                !  >>>  x=0.2
                !  >>>
                !  >>>  for i in range(0,6):
                !  >>>      value = numerical_approx(P_a(i,x))
                !  >>>      print("l = {}: {}".format(i, value))
                !  >>>
                !  >>>  for i in range(0,6):
                !  >>>      value = numerical_approx(P_ad(i,x))
                !  >>>      print("l = {}: {}".format(i, value))
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

                ! Values were generated with sageMath
                !
                ! source code:
                !  >>>  var('l')
                !  >>>  var('m')
                !  >>>  var('x')
                !  >>>
                !  >>>  m=1
                !  >>>
                !  >>>  P_d(l,x) = derivative(legendre_P(l,x), x, m)
                !  >>>  P_a(l,x) = (1-x**2)**(m/2) * P_d(l,x)
                !  >>>  P_ad(l,x) = derivative(P_a(l,x), x).simplify_full()
                !  >>>
                !  >>>  x=0.2
                !  >>>
                !  >>>  for i in range(0,6):
                !  >>>      value = numerical_approx(P_a(i,x))
                !  >>>      print("l = {}: {}".format(i, value))
                !  >>>
                !  >>>  for i in range(0,6):
                !  >>>      value = numerical_approx(P_ad(i,x))
                !  >>>      print("l = {}: {}".format(i, value))
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

                ! Values were generated with sageMath
                !
                ! source code:
                !  >>>  var('l')
                !  >>>  var('m')
                !  >>>  var('x')
                !  >>>
                !  >>>  m=2
                !  >>>
                !  >>>  P_d(l,x) = derivative(legendre_P(l,x), x, m)
                !  >>>  P_a(l,x) = (-1)**m * (1-x**2)**(m/2) * P_d(l,x)
                !  >>>  P_ad(l,x) = derivative(P_a(l,x), x).simplify_full()
                !  >>>
                !  >>>  x=0.2
                !  >>>
                !  >>>  for i in range(0,6):
                !  >>>      value = numerical_approx(P_a(i,x))
                !  >>>      print("l = {}: {}".format(i, value))
                !  >>>
                !  >>>  for i in range(0,6):
                !  >>>      value = numerical_approx(P_ad(i,x))
                !  >>>      print("l = {}: {}".format(i, value))
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

                ! Values were generated with sageMath
                !
                ! source code:
                !  >>>  var('l')
                !  >>>  var('m')
                !  >>>  var('x')
                !  >>>
                !  >>>  m=2
                !  >>>
                !  >>>  P_d(l,x) = derivative(legendre_P(l,x), x, m)
                !  >>>  P_a(l,x) = (1-x**2)**(m/2) * P_d(l,x)
                !  >>>  P_ad(l,x) = derivative(P_a(l,x), x).simplify_full()
                !  >>>
                !  >>>  x=0.2
                !  >>>
                !  >>>  for i in range(0,6):
                !  >>>      value = numerical_approx(P_a(i,x))
                !  >>>      print("l = {}: {}".format(i, value))
                !  >>>
                !  >>>  for i in range(0,6):
                !  >>>      value = numerical_approx(P_ad(i,x))
                !  >>>      print("l = {}: {}".format(i, value))
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

            function test_lib_math_associated_legendre_polynomial_m0() result (rv)
                implicit none
                ! dummy
                logical :: rv

                ! auxiliary
                integer(kind=4) :: i
                double precision :: x
                integer(kind=4), parameter :: n = 5
                integer(kind=4), parameter :: m = 0
                double precision, dimension(0:n) :: pm
                double precision, dimension(0:n) :: pd

                double precision, dimension(0:n) :: ground_truth_pm
                double precision, dimension(0:n) :: ground_truth_pd

                double precision :: buffer

                ! Values were generated with sageMath
                !
                ! source code:
                !  >>>  var('l')
                !  >>>  var('m')
                !  >>>  var('x')
                !  >>>
                !  >>>  m=0
                !  >>>
                !  >>>  P_d(l,x) = derivative(legendre_P(l,x), x, m)
                !  >>>  P_a(l,x) = (-1)**m * (1-x**2)**(m/2) * P_d(l,x)
                !  >>>  P_ad(l,x) = derivative(P_a(l,x), x).simplify_full()
                !  >>>
                !  >>>  x=0.2
                !  >>>
                !  >>>  for i in range(0,6):
                !  >>>      value = numerical_approx(P_a(i,x))
                !  >>>      print("l = {}: {}".format(i, value))
                !  >>>
                !  >>>  for i in range(0,6):
                !  >>>      value = numerical_approx(P_ad(i,x))
                !  >>>      print("l = {}: {}".format(i, value))
                ground_truth_pm(0) = 1.00000000000000_8
                ground_truth_pm(1) = 0.200000000000000_8
                ground_truth_pm(2) = -0.440000000000000_8
                ground_truth_pm(3) = -0.280000000000000_8
                ground_truth_pm(4) = 0.232000000000000_8
                ground_truth_pm(5) = 0.307520000000000_8

                ground_truth_pd(0) = -0.000000000000000_8
                ground_truth_pd(1) = 1.00000000000000_8
                ground_truth_pd(2) = 0.600000000000000_8
                ground_truth_pd(3) = -1.20000000000000_8
                ground_truth_pd(4) = -1.36000000000000_8
                ground_truth_pd(5) = 0.888000000000000_8

                x = 0.2

                call lib_math_associated_legendre_polynomial(x, m, n, pm, pd, .true.)


                rv = .true.
                print *, "test_lib_math_associated_legendre_polynomial_m0:"
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

            end function test_lib_math_associated_legendre_polynomial_m0

            function test_lib_math_associated_legendre_polynomial_theta_m1_without_p() result(rv)
                implicit none
                ! dummy
                logical :: rv

                ! auxiliary
                integer(kind=4) :: i
                double precision :: theta
                integer(kind=4), parameter :: n = 5
                integer(kind=4), parameter :: m = 1
                double precision, dimension(0:n) :: pm
                double precision, dimension(0:n) :: pd

                double precision, dimension(0:n) :: ground_truth_pm
                double precision, dimension(0:n) :: ground_truth_pd

                double precision :: buffer

                ! Values were generated with sageMath
                !
                ! source code:
                !  >>>  var('l')
                !  >>>  var('m')
                !  >>>  var('x')
                !  >>>
                !  >>>  m=1
                !  >>>
                !  >>>  P_d(l,x) = derivative(legendre_P(l,x), x, m)
                !  >>>  P_a(l,x) = (1-x**2)**(m/2) * P_d(l,x)
                !  >>>  P_ad(l,x) = derivative(P_a(l,x), x).simplify_full()
                !  >>>
                !  >>>  x=cos(0.2)
                !  >>>
                !  >>>  for i in range(0,6):
                !  >>>      value = numerical_approx(P_ad(i,x))
                !  >>>      print("l = {}: {}".format(i, value))
                ground_truth_pd(0) = 0.000000000000000_8
                ground_truth_pd(1) = -4.93315487558688_8
                ground_truth_pd(2) = -13.9084526582466_8
                ground_truth_pd(3) = -25.2179729025493_8
                ground_truth_pd(4) = -36.4802655764097_8
                ground_truth_pd(5) = -44.8436276630205_8

                theta = 0.2

                call genlgp(theta, ground_truth_pm, n+1)

                call lib_math_associated_legendre_polynomial_theta(theta, m, n, pm, pd, .false.)


                rv = .true.
                print *, "test_lib_math_associated_legendre_polynomial_m0:"
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

            end function test_lib_math_associated_legendre_polynomial_theta_m1_without_p

            subroutine genlgp(theta,pnmllg,nc)
                  implicit none
            !     ........................................................
            !     .  calculate associated Legendre functions (argument   .
            !     .    cos(theta)) divided by sin(theta) for m = 1       .
            !     .  generate first two orders by formula and remaining  .
            !     .    orders by recursion                               .
            !     .                                                      .
            !     .  pnmllg = associated Legendre function/sin(theta)    .
            !     .  nc = number of orders (0 to nc-1)                   .
            !     .  the order of the associated Legendre functions is   .
            !     .    incremented by one in the pnmllg(*) array         .
            !     ........................................................
                  real(kind=8) theta, pnmllg, costh, rn
                  integer nc, n
                  dimension pnmllg(nc)
                  costh = cos(theta)
            !     ..............................
            !     .  calculate orders 0 and 1  .
            !     ..............................
                  pnmllg(1) = 0.0                                                   !eq 4.70a
                  pnmllg(2) = 1.0                                                   !eq 4.70b
            !     .................................................
            !     .  recur upward to obtain all remaining orders  .
            !     .................................................
                  do 10 n = 3,nc
                  rn = real(n-1)
                  pnmllg(n) = ((2.0*rn-1.0)*costh*pnmllg(n-1) &
                              -rn*pnmllg(n-2))/(rn-1.0)                            !eq 4.71
            10    continue
                  return
          end subroutine genlgp

            function test_lib_math_legendre_polynomial() result (rv)
                implicit none
                ! dummy
                logical :: rv

                ! auxiliary
                integer(kind=4) :: i
                double precision :: x
                integer(kind=4), parameter :: n = 5
                double precision, dimension(0:n) :: pm
                double precision, dimension(0:n) :: pd

                double precision, dimension(0:n) :: ground_truth_pm
                double precision, dimension(0:n) :: ground_truth_pd

                double precision :: buffer

                ! Values were generated with sageMath
                !
                ! source code:
                !  >>> var('l')
                !  >>> var('x')
                !  >>>
                !  >>> P(l,x) = legendre_P(l,x)
                !  >>> P_d(l,x) = derivative(legendre_P(l,x), x)
                !  >>>
                !  >>> x=0.4
                !  >>>
                !  >>> for i in range(0,6):
                !  >>>     value = numerical_approx(P(i,x))
                !  >>>     print("l = {}: {}".format(i, value))
                !  >>>
                !  >>> print("\nderivative:")
                !  >>> for i in range(0,6):
                !  >>>     value = numerical_approx(P_d(i,x))
                !  >>>     print("l = {}: {}".format(i, value))
                ground_truth_pm(0) = 1.00000000000000_8
                ground_truth_pm(1) = 0.400000000000000_8
                ground_truth_pm(2) = -0.260000000000000_8
                ground_truth_pm(3) = -0.440000000000000_8
                ground_truth_pm(4) = -0.113000000000000_8
                ground_truth_pm(5) = 0.270640000000000_8

                ground_truth_pd(0) = -0.000000000000000_8
                ground_truth_pd(1) = 1.00000000000000_8
                ground_truth_pd(2) = 1.20000000000000_8
                ground_truth_pd(3) = -0.300000000000000_8
                ground_truth_pd(4) = -1.88000000000000_8
                ground_truth_pd(5) = -1.31700000000000_8

                x = 0.4

                call lib_math_legendre_polynomial(x, n, pm, pd)


                rv = .true.
                print *, "test_lib_math_legendre_polynomial:"
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

            end function test_lib_math_legendre_polynomial

        end function lib_math_legendre_test_functions

end module lib_math_legendre
