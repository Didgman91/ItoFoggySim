module lib_math_legendre
    use lib_math_factorial
    implicit none

    private

    ! --- public ---
    public :: lib_math_associated_legendre_polynomial
    public :: lib_math_associated_legendre_polynomial_with_negative_m
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
        !   pm: double precision, dimension(2, 0:n)
        !       pm(1,:): result of the associated Legendre polynomial with a negativ m
        !       pm(2,:): result of the associated Legendre polynomial with a positiv m
        !   pd: double precision, dimension(2, 0:n)
        !       pd(1,:): result of the deriviative of the associated Legendre polynomial with a negativ m
        !       pd(2,:): result of the deriviative of the associated Legendre polynomial with a positiv m
        !
        subroutine lib_math_associated_legendre_polynomial_with_negative_m(x, m, n, pm, pd, condon_shortley_phase)
            implicit none
            ! dummy
            double precision, intent(in) :: x
            integer(kind=4), intent(in) :: m
            integer(kind=4), intent(in) :: n

            double precision, dimension(2, 0:n), intent(inout) :: pm
            double precision, dimension(2, 0:n), intent(inout) :: pd
            logical, optional :: condon_shortley_phase

            ! auxiliary
            double precision, dimension(0:n) :: buffer_pm
            double precision, dimension(0:n) :: buffer_pd
            real(kind=8) :: buffer

            call LPMNS(m, n, x, buffer_pm, buffer_pd)

            pm(2, :) = buffer_pm
            pd(2, :) = buffer_pd

            ! calculation with negative m
            buffer = lib_math_factorial_get_n_minus_m_divided_by_n_plus_m(n, m)

            if ( IAND(m, 1) .eq. 1) then
                ! m is odd
                buffer = -buffer
            end if

            pm(1, :) = buffer * pm(2, :)
            pd(1, :) = buffer * pd(2, :)

            if (present(condon_shortley_phase)) then
                if (condon_shortley_phase) then
                    if (IAND(m, 1) .eq. 1) then
                        ! m is odd
                        pm(2,:) = -pm(2,:)
                        pd(2,:) = -pd(2,:)
                    end if
!                    pm = (-1.0_8)**m * pm
!                    pd = (-1.0_8)**m * pd
                end if
            end if

        end subroutine lib_math_associated_legendre_polynomial_with_negative_m

        ! calculates the associated Legendre polynomial
        !
        ! Formula
        ! ----
        !   pm = P_1l(cos(theta)) / sin(theta)
        !
        !   pd = derivation[P_1l(cos(theta)), theta ]
        !
        !
        ! Argument
        ! ----
        !   theta: double precision
        !       input value
        !   n_max: integer
        !       maximum degree of the polynomial, > 0
        !   condon_shortley_phase: boolean, optional(std: false)
        !       true: with Condon–Shortley phase
        !       false: without Condon–Shortley phase
        !
        ! Results
        ! ----
        !   pi_mn: double precision, dimension(n, -n:n)
        !       result of the associated Legendre polynomial (cos(theta)) / sin(theta)
        !       1st-dimension: degree n of the associated Legendre function
        !       2nd-dimension: order m of the associated Legendre function
        !   tau_mn: double precision, dimension(n, -n:n)
        !       result of the deriviative of the associated Legendre polynomial
        !       1st-dimension: degree n of the deriviative of the associated Legendre function
        !       2nd-dimension: order m of the deriviative of the associated Legendre function
        !
        !
        ! Formula
        ! ----
        !   pi_mn
        !       LaTeX: $$ \pi_{m n}(\cos \theta)=\frac{m}{\sin \theta} P_{n}^{m}(\cos \theta) $$
        !
        !   tau_mn
        !       LaTeX: $$ \tau_{m n}(\cos \theta)=\frac{d}{d \theta} P_{n}^{m}(\cos \theta) $$
        !
        ! Reference: Electromagnetic scattering by an aggregate of spheres, Yu-lin Xu, Appendix A
        !
        !
        subroutine lib_math_associated_legendre_polynomial_theta(theta, n_max, pi_nm, tau_nm)
            implicit none
            ! dummy
            double precision, intent(in) :: theta
            integer(kind=4), intent(in) :: n_max

            real(kind=8), dimension(0:n_max, -n_max:n_max), intent(inout) :: pi_nm
            real(kind=8), dimension(0:n_max, -n_max:n_max), intent(inout) :: tau_nm

            ! auxiliary
            logical :: func_rv
            integer(kind=4) :: n
            integer(kind=4) :: m

            real(kind=8), dimension(0:n_max+1, -n_max:n_max) :: m_pi_nm

            real(kind=8) :: x

            x = cos(theta)

            ! Recurrence Formulas for pi_mn, tau_mn
            ! where -1 .le. x .le. 1
            ! Special Values, eq. A6
            m_pi_nm(0,0) = 0.0_8
            m_pi_nm(1,0) = 0
            m_pi_nm(0,1) = 0
            m_pi_nm(1,1) = 1

            tau_nm(0,0) = 0
            tau_nm(1,0) = -sqrt(1.0_8-x*x)
            tau_nm(0,1) = 0
            tau_nm(1,1) = x

            if (x .eq. 1.0_8) then
                ! Special Values, eq. 7
                m_pi_nm = 0
                tau_nm = 0
                do n = 1, n_max
                    m_pi_nm(n, -1) = 0.5_8
                    m_pi_nm(n, 1) = real(n*(n + 1), kind=8) / 2.0_8

                    tau_nm(n, -1) = -0.5_8
                    tau_nm(n, 1) = real(n*(n + 1), kind=8) / 2.0_8
                end do
            else if (x .eq. -1.0_8) then
                ! Special Values, eq. 7
                m_pi_nm = 0
                tau_nm = 0
                do n = 1, n_max
                    m_pi_nm(n, -1) = 0.5_8
                    tau_nm(n, -1) = -0.5_8
                    ! check if n+1 is even or odd
                    if (IAND(n+1_4, 1_4) .eq. 1) then
                        ! n+1 is odd
                        m_pi_nm(n, -1) = -m_pi_nm(n, -1)
                    else
                        ! n+1 is even
                        ! --> n is odd
                        tau_nm(n, -1) = -tau_nm(n, -1)
                    end if

                    m_pi_nm(n, 1) = real(n*(n + 1), kind=8) / 2.0_8
                    tau_nm(n, 1) = real(n*(n + 1), kind=8) / 2.0_8
                    ! check if n+1 is even or odd
                    if (IAND(n+1_4, 1_4) .eq. 1) then
                        ! is odd
                        m_pi_nm(n, -1) = -m_pi_nm(n, -1)
                    else
                        ! n+1 is even
                        ! --> n is odd
                        tau_nm(n, 1) = -tau_nm(n, 1)
                    end if
                end do
            else
                ! Recurrence Relations eq. A2, A3, A4, A5
                ! eq. A3: Electromagnetic scattering by an aggregate of spheres: errata

                ! Calculation of pi_mn
                ! ----
                !
                ! Step 1
                ! -----
                !        __|0 --> n                 n_max|
                !    n_max -------------------------------
                ! m        |  ^  |  ^  |     |     |     |
                !  ^       |--|-----|--------------------|   pi_m+1,n
                !  |       |  |  |  |  |     |     |     |
                !          |-----------------------------|
                !          |  0  |  1  |     |     |     |
                !          |-----------------------------|
                !      0 --|  0  |  0  |     |     |     |
                !          |-----------------------------|
                !          |     |     |     |     |     |
                !          |-----------------------------|
                !          |     |     |     |     |     |
                !          |-----------------------------|
                !          |     |     |     |     |     |
                !  -n_max__-------------------------------
                !
                do n=0, 1
                    do m=2, n_max-1
                        func_rv = set_pi_m_plus_1_n(x, m, n, m_pi_nm)
                    end do
                end do

                ! Step 2
                ! -----
                !        __|0 --> n                 n_max|
                !    n_max -------------------------------
                ! m        |  x  |  x  |     |     |     |
                !  ^       |-----------------------------|
                !  |       |  x  |  x  |     |     |     |
                !          |-----------------------------|
                !          |  0  |  1  |     |     |     |
                !          |-----------------------------|
                !      0 --|  0  |  0  |     |     |     |
                !          |-----------------------------|
                !          |  |  |  |  |     |     |     |
                !          |--|-----|--------------------|
                !          |  |  |  |  |     |     |     |   pi_-m,n
                !          |--|-----|--------------------|
                !          |  V  |  V  |     |     |     |
                !  -n_max__-------------------------------             x: calculated element
                !
                do n=0, 1
                    do m=1, n_max
                        func_rv = set_pi_minus_m_n(m, n, m_pi_nm)
                    end do
                end do

                ! Step 3
                ! -----
                !        __|0 --> n                 n_max|
                !    n_max -------------------------------
                ! m        |  x  |  x  |  ------------>  |
                !  ^       |-----------------------------|
                !  |       |  x  |  x  |  ------------>  |
                !          |-----------------------------|
                !          |  0  |  1  |  ------------>  |
                !          |-----------------------------|
                !      0 --|  0  |  0  |  ------------>  |
                !          |-----------------------------|
                !          |  x  |  x  |  ------------>  |
                !          |-----------------------------|
                !          |  x  |  x  |  ------------>  |   pi_m,n+1
                !          |-----------------------------|
                !          |  x  |  x  |  ------------>  |
                !  -n_max__-------------------------------             x: calculated element
                !
                do n=1, n_max
                    do m=-n_max, n_max
                        func_rv = set_pi_m_n_plus_1(x, m, n, m_pi_nm)
                    end do
                end do

                ! Calculation of tau_mn
                ! ----
                !
                ! Step 1
                ! -----
                !        __|0 --> n                 n_max|
                !    n_max -------------------------------
                ! m        |     |     |     |     |     |
                !  ^       |-----------------------------|
                !  |       |     |     |     |     |     |
                !          |-----------------------------|
                !          |  0  |  x  |     |     |     |
                !          |-----------------------------|
                !      0 --|  0  |  x  -------------->   |   tau_0,n+1
                !          |-----------------------------|
                !          |     |     |     |     |     |
                !          |-----------------------------|
                !          |     |     |     |     |     |
                !          |-----------------------------|
                !          |     |     |     |     |     |
                !  -n_max__-------------------------------             x: calculated element
                !
                do n=2, n_max-1
                    func_rv = set_tau_0_n_plus_1(x, n, tau_nm)
                end do

                ! Step 2
                ! -----
                !        __|0 --> n                 n_max|
                !    n_max -------------------------------
                ! m        |  ^  |  ^  |  ^  |  ^  |  ^  |
                !  ^       |--|-----|-----|-----|-----|--|
                !  |       |  |  |  |  |  |  |  |  |  |  |   tau_mn
                !          |--------------|-----|-----|--|
                !          |  0  |  x  |  |  |  |  |  |  |
                !          |-----------------------------|
                !      0 --|  0  |  x  |  x  |  x  |  x  |
                !          |-----------------------------|
                !          |     |     |     |     |     |
                !          |-----------------------------|
                !          |     |     |     |     |     |
                !          |-----------------------------|
                !          |     |     |     |     |     |
                !  -n_max__-------------------------------             x: calculated element
                !
                do n=0, 1
                    do m=2, n_max
                        func_rv = set_tau_m_n(x, m, n, m_pi_nm, tau_nm)
                    end do
                end do

                do n=2, n_max
                    do m=1, n_max
                        func_rv = set_tau_m_n(x, m, n, m_pi_nm, tau_nm)
                    end do
                end do

                ! Step 3
                ! -----
                !        __|0 --> n                 n_max|
                !    n_max -------------------------------
                ! m        |  x  |  x  |  x  |  x  |  x  |
                !  ^       |-----------------------------|
                !  |       |  x  |  x  |  x  |  x  |  x  |
                !          |-----------------------------|
                !          |  0  |  x  |  x  |  x  |  x  |
                !          |-----------------------------|
                !      0 --|  0  |  x  |  x  |  x  |  x  |
                !          |-----------------------------|
                !          |  |  |  |  |  |  |  |  |  |  |
                !          |--|-----|-----|-----|-----|--|
                !          |  |  |  |  |  |  |  |  |  |  |
                !          |--|-----|-----|-----|-----|--|
                !          |  V  |  V  |  V  |  V  |  V  |
                !  -n_max__-------------------------------             x: calculated element
                !
                do n=0, n_max
                    do m=1, n_max
                        func_rv = set_tau_minus_m_n(m, n, tau_nm)
                    end do
                end do

!                ! Remove the factor from the definition (eq. A1) of tau to obtain
!                ! only the derivative of the associated Legendre polynomial.
!                tau_nm = tau_nm / (-sqrt(1-x*x))
            end if

            pi_nm = m_pi_nm(0:n_max, -n_max:n_max)

            ! reference: Absorption and Scattering of Light by smal particles eq. 4.47
            !  --> m=1
!            pm(0) = 0
!            pm(1) = 1
!            do i=2, n
!                pm(i) = (2*n - 1)/(n - 1) * cos_theta * pm(i-1) - n/(n - 1) * pm(i-2)
!                pd(i) = i * cos_theta * pm(i) - (n + 1) * pm(i-1)
!            end do

!            if (present(condon_shortley_phase)) then
!                if (condon_shortley_phase) then
!                    if (IAND(m, 1) .eq. 1) then
!                        pm = -pm
!                        pd = -pd
!                    end if
!                end if
!            end if

            contains

                ! first line of eq. A1
                ! Restriction
                ! ----
                !   - pi_nm(n, m) has to be known
                !   - pi_nm(n-1, m) has to be known
                function set_pi_m_n_plus_1(x, m, n, pi_nm) result (rv)
                    implicit none
                    ! dummy
                    real(kind=8) :: x
                    integer(kind=4) :: m
                    integer(kind=4) :: n
                    real(kind=8), dimension(0:n_max+1, -n_max:n_max), intent(inout) :: pi_nm

                    logical :: rv

                    ! auxiliaray
                    real(kind=8) :: denominator

                    denominator = n - m + 1

                    pi_nm(n+1, m) = 0.0_8

                    if (pi_nm(n,m) .ne. 0.0_8) then
                        pi_nm(n+1, m) = (2*n + 1) / denominator * x * pi_nm(n,m)
                    end if
                    if (pi_nm(n-1,m) .ne. 0.0_8) then
                        pi_nm(n+1, m) = pi_nm(n+1, m) - (n + m) / denominator * pi_nm(n-1, m)
                    end if

!                    pi_nm(n+1, m) = (2*n + 1) / denominator * x * pi_nm(n,m) &
!                                    -(n + m) / denominator * pi_nm(n-1, m)

                    rv = .true.
                end function

                ! second line of eq. A1
                ! Restriction
                ! ----
                !   - m .ne. 1
                !   - pi_nm(n, m) has to be known
                !   - pi_nm(n, m-1) has to be known
                function set_pi_m_plus_1_n(x, m, n, pi_nm) result (rv)
                    implicit none
                    ! dummy
                    real(kind=8) :: x
                    integer(kind=4) :: m
                    integer(kind=4) :: n
                    real(kind=8), dimension(0:n_max+1, -n_max:n_max), intent(inout) :: pi_nm

                    logical :: rv

                    pi_nm(n, m+1) = 0.0_8

                    if (pi_nm(n,m) .ne. 0.0_8) then
                        pi_nm(n, m+1) = 2*(m + 1)*x / sqrt(1 - x*x) * pi_nm(n,m)
                    end if
                    if (pi_nm(n,m-1) .ne. 0.0_8) then
                        pi_nm(n, m+1) = pi_nm(n, m+1) - (m + 1)*(n + m)*(n - m + 1)/(m - 1) * pi_nm(n,m-1)
                    end if

!                    pi_nm(n, m+1) = 2*(m + 1)*x / sqrt(1 - x*x) * pi_nm(n,m) &
!                                    -(m - 1)*(n + m)*(n - m + 1)/(m - 1) * pi_nm(n,m-1)

                    rv = .true.
                end function

                ! eq. A2
                ! Restriction
                ! ----
                !   - n .nq. 1
                !   - pi_nm(n-1, n-1) has to be known
                function set_pi_n_n(x, n, pi_nm) result (rv)
                    implicit none
                    ! dummy
                    real(kind=8) :: x
                    integer(kind=4) :: n
                    real(kind=8), dimension(0:n_max+1, -n_max:n_max), intent(inout) :: pi_nm

                    logical :: rv

                    pi_nm(n,n) = sqrt(1 - x) * n*(2*n - 1)/(n - 1) * pi_nm(n-1,n-1)

                    rv = .true.
                end function

                ! first line of eq. A4
                !
                ! Restriction
                ! ----
                !   - pi_nm(n, m) has to be known
                function set_pi_minus_m_n(m,n,pi_nm) result (rv)
                    implicit none
                    ! dummy
                    integer(kind=4) :: m
                    integer(kind=4) :: n
                    real(kind=8), dimension(0:n_max+1, -n_max:n_max), intent(inout) :: pi_nm

                    logical :: rv

                    if (pi_nm(n, m) .eq. 0.0_8) then
                        pi_nm(n, -m) = 0.0_8
                    else
                        pi_nm(n, -m) = lib_math_factorial_get_n_minus_m_divided_by_n_plus_m(n,m) * pi_nm(n,m)
                        ! check if m+1 is even or odd
                        if (IAND(m+1_4, 1_4) .eq. 1) then
                            ! m+1 is odd
                            pi_nm(n, -m) = -pi_nm(n, -m)
                        end if
                    end if

                    rv = .true.
                end function

                ! first line eq. A3
                !
                ! Restriction
                ! ----
                !   - m .ne. 0
                !   - pi_nm(n+1, m) has to be known
                !   - pi_nm(n, m) has to be known
                function set_tau_m_n(x, m, n, pi_nm, tau_nm) result (rv)
                    implicit none
                    ! dummy
                    real(kind=8) :: x
                    integer(kind=4) :: m
                    integer(kind=4) :: n
                    real(kind=8), dimension(0:n_max+1, -n_max:n_max), intent(inout) :: pi_nm
                    real(kind=8), dimension(0:n_max, -n_max:n_max), intent(inout) :: tau_nm

                    logical :: rv

                    tau_nm(n,m) = 0.0_8

                    if (pi_nm(n+1, m) .ne. 0.0_8) then
                        tau_nm(n,m) = (n - m + 1)/m * pi_nm(n+1, m)
                    end if

                    if (pi_nm(n, m) .ne. 0.0_8) then
                         tau_nm(n,m) = tau_nm(n,m) - (n + 1)/m * x * pi_nm(n, m)
                    end if

                    rv = .true.
                end function

                ! second line eq. A3
                !
                ! eq. A3: Electromagnetic scattering by an aggregate of spheres: errata 2001
                !
                ! Restriction
                ! ----
                !   - n .ne. 0
                !   - tau_nm(n, m) has to be known
                !   - tau_nm(n-1, m) has to be known
                function set_tau_0_n_plus_1(x, n, tau_nm) result (rv)
                    implicit none
                    ! dummy
                    real(kind=8) :: x
                    integer(kind=4) :: n
                    real(kind=8), dimension(0:n_max, -n_max:n_max), intent(inout) :: tau_nm

                    logical :: rv

                    tau_nm(n+1,0) = 0.0_8

                    if (tau_nm(n, 0) .ne. 0.0_8) then
                        tau_nm(n+1,0) = (2*n + 1)/n * x * tau_nm(n, 0)
                    end if

                    if (tau_nm(n-1, 0) .ne. 0.0_8) then
                        tau_nm(n+1,0) = tau_nm(n+1,0) - (n + 1)/n * tau_nm(n-1, 0)
                    end if

                    rv = .true.
                end function

                ! first line eq. A5
                !
                ! Restriction
                ! ----
                !   - tau_nm(n, m) has to be known
                function set_tau_minus_m_n(m, n, tau_nm) result (rv)
                    implicit none
                    ! dummy
                    integer(kind=4) :: m
                    integer(kind=4) :: n
                    real(kind=8), dimension(0:n_max, -n_max:n_max), intent(inout) :: tau_nm

                    logical :: rv

                    if (tau_nm(n, m) .eq. 0.0_8) then
                        tau_nm(n, -m) = 0.0_8
                    else
                        tau_nm(n, -m) = lib_math_factorial_get_n_minus_m_divided_by_n_plus_m(n,m) * tau_nm(n,m)
                        ! check if m is even or odd
                        if (IAND(m, 1_4) .eq. 1) then
                            ! m is odd
                            tau_nm(n, -m) = -tau_nm(n, -m)
                        end if
                    end if

                    rv = .true.
                end function

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
!            if (.not. test_lib_math_associated_legendre_polynomial_theta_m1()) then
!                rv = rv + 1
!            end if
            if (.not. test_lib_math_associated_legendre_polynomial_theta_n0_3()) then
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

            function test_lib_math_associated_legendre_polynomial_theta_m1() result(rv)
                implicit none
                ! dummy
                logical :: rv

                ! parameter
                integer(kind=4), parameter :: n_max = 5

                ! auxiliary
                integer(kind=4) :: i
                double precision :: theta

                real(kind=8), dimension(0:n_max, -n_max:n_max) :: pi_nm
                real(kind=8), dimension(0:n_max, -n_max:n_max) :: tau_nm

                double precision, dimension(0:n_max) :: ground_truth_pi_n1
                double precision, dimension(0:n_max) :: ground_truth_tau_n1

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
                ground_truth_tau_n1(0) = 0.000000000000000_8
                ground_truth_tau_n1(1) = -4.93315487558688_8
                ground_truth_tau_n1(2) = -13.9084526582466_8
                ground_truth_tau_n1(3) = -25.2179729025493_8
                ground_truth_tau_n1(4) = -36.4802655764097_8
                ground_truth_tau_n1(5) = -44.8436276630205_8

                theta = 0.2_8

                call genlgp(theta, ground_truth_pi_n1, n_max+1)

                call lib_math_associated_legendre_polynomial_theta(theta, n_max, pi_nm, tau_nm)


                rv = .true.
                print *, "test_lib_math_associated_legendre_polynomial_theta_m1:"
                do i=0, n_max
                    buffer = abs(pi_nm(i, 1) - ground_truth_pi_n1(i))
                    if (buffer .gt. ground_truth_e) then
                        print *, "  ", i , "difference: ", buffer, " : FAILED"
                        rv = .false.
                    else
                        print *, "  ", i, ": OK"
                    end if
                end do

                print*, "  deriviation:"
                do i=0, n_max
                    buffer = abs(tau_nm(i, 1) - ground_truth_tau_n1(i))
                    if (buffer .gt. ground_truth_e) then
                        print *, "  ", i , "difference: ", buffer, " : FAILED"
                        rv = .false.
                    else
                        print *, "  ", i, ": OK"
                    end if
                end do

            end function test_lib_math_associated_legendre_polynomial_theta_m1

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

          function test_lib_math_associated_legendre_polynomial_theta_n0_3() result(rv)
                implicit none
                ! dummy
                logical :: rv

                ! parameter
                integer(kind=4), parameter :: n_max = 3

                ! auxiliary
                integer(kind=4) :: i
                integer(kind=4) :: ii
                double precision :: theta

                real(kind=8), dimension(0:n_max, -n_max:n_max) :: pi_nm
                real(kind=8), dimension(0:n_max, -n_max:n_max) :: tau_nm

                double precision, dimension(0:n_max, -n_max:n_max) :: ground_truth_pi_nm
                double precision, dimension(0:n_max, -n_max:n_max) :: ground_truth_tau_nm

                double precision :: buffer

                ! Values were generated with WolframAplpha
                !
                ! source code:
                !  >>>  N[m*LegendreP[n, m, cos (t)]/sin (t), 16], where n = 0, m = [-2, -1, 0, 1, 2], t = 2
                !  >>>  N[m*LegendreP[n, m, cos (t)]/sin (t), 16], where n = 1, m = [-2, -1, 0, 1, 2], t = 2
                !  >>>  N[m*LegendreP[n, m, cos (t)]/sin (t), 16], where n = 2, m = [-2, -1, 0, 1, 2], t = 2

                ground_truth_pi_nm(0, -3) = 2.077165092735346_8
                ground_truth_pi_nm(0, -2) = -2.667464736243829_8
                ground_truth_pi_nm(0, -1) = 1.712759410407380_8
                ground_truth_pi_nm(0,  0) = 0_8
                ground_truth_pi_nm(0,  1) = 0_8
                ground_truth_pi_nm(0,  2) = 0_8
                ground_truth_pi_nm(0,  3) = 0_8

                ground_truth_pi_nm(1, -3) = 1.341772398969518_8
                ground_truth_pi_nm(1, -2) = -1.408290820299577_8
                ground_truth_pi_nm(1, -1) = 0.5000000000000000_8
                ground_truth_pi_nm(1,  0) = 0_8
                ground_truth_pi_nm(1,  1) = 1.000000000000000_8
                ground_truth_pi_nm(1,  2) = 0_8
                ground_truth_pi_nm(1,  3) = 0_8

                ground_truth_pi_nm(2, -3) = 0.4958414335756772_8
                ground_truth_pi_nm(2, -2) = -0.2273243567064204_8
                ground_truth_pi_nm(2, -1) = -0.2080734182735712_8
                ground_truth_pi_nm(2,  0) = 0_8
                ground_truth_pi_nm(2,  1) = -1.248440509641427_8
                ground_truth_pi_nm(2,  2) = 5.455784560954090_8
                ground_truth_pi_nm(2,  3) = 0_8

                ground_truth_pi_nm(3, -3) = 0.05167636315198787_8
                ground_truth_pi_nm(3, -2) = 0.09460031191349103_8
                ground_truth_pi_nm(3, -1) = -0.01676363151987872_8
                ground_truth_pi_nm(3,  0) = 0_8
                ground_truth_pi_nm(3,  1) = -0.2011635782385447_8
                ground_truth_pi_nm(3,  2) = -11.35203742961892_8
                ground_truth_pi_nm(3,  3) = 37.20698146943127_8

                theta = 2.0_8

                call lib_math_associated_legendre_polynomial_theta(theta, n_max, pi_nm, tau_nm)


                rv = .true.
                print *, "test_lib_math_associated_legendre_polynomial_theta_n0_3:"
                do i=0, n_max
                    do ii=-n_max, n_max
                        buffer = abs(pi_nm(i, ii) - ground_truth_pi_nm(i, ii))
                        if (buffer .gt. ground_truth_e) then
                            print *, "  n: ", i, " m: ", ii , "difference: ", buffer, " : FAILED"
                            rv = .false.
                        else
                            print *, "  n: ", i, " m: ", ii, "OK"
                        end if
                    end do
                    print*, ""
                end do

!                print*, "  deriviation:"
!                do i=0, n_max
!                    buffer = abs(tau_nm(i, 1) - ground_truth_tau_n1(i))
!                    if (buffer .gt. ground_truth_e) then
!                        print *, "  ", i , "difference: ", buffer, " : FAILED"
!                        rv = .false.
!                    else
!                        print *, "  ", i, ": OK"
!                    end if
!                end do

            end function test_lib_math_associated_legendre_polynomial_theta_n0_3

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
