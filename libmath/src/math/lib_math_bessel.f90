#define _DEBUG_

module lib_math_bessel
    use lib_math_factorial
    implicit none

    private

    ! --- interface ---
    ! calculates the the spherical Bessel function of the first kind
    !
    ! symbol: j_n
    !
    ! Argument
    ! ----
    !   x: double precision .OR. complex
    !   fnu: integer
    !       ORDER OF INITIAL FUNCTION, FNU.GE.0
    !   n: integer
    !       NUMBER OF MEMBERS OF THE SEQUENCE, N.GE.1
    !
    ! Returns
    ! ----
    !   rv: array<double precision>
    !
    ! LaTeX: $$ j_{n}(\rho)=\sqrt{\frac{\pi}{2 \rho}} J_{n+1 / 2}(\rho) $$
    interface lib_math_bessel_spherical_first_kind
        module procedure lib_math_bessel_spherical_first_kind_real
        module procedure lib_math_bessel_spherical_first_kind_cmplx
    end interface

    ! calculates the the spherical Bessel function of the second kind
    !
    ! symbol: y_n
    !
    ! Argument
    ! ----
    !   x: double precision .OR. complex
    !   fnu: integer
    !       ORDER OF INITIAL FUNCTION, FNU.GE.0
    !   n: integer
    !       NUMBER OF MEMBERS OF THE SEQUENCE, N.GE.1
    !
    ! Returns
    ! ----
    !   rv: array<double precision>
    !
    ! LaTeX: $$ y_{n}(\rho)=\sqrt{\frac{\pi}{2 \rho}} Y_{n+1 / 2}(\rho) $$
    interface lib_math_bessel_spherical_second_kind
        module procedure lib_math_bessel_spherical_second_kind_real
        module procedure lib_math_bessel_spherical_second_kind_cmplx
    end interface

    ! calculates the the spherical Bessel function of the third kind
    !
    ! symbol: h_n^(1)
    !
    ! Argument
    ! ----
    !   x: double precision .OR. complex
    !   fnu: integer
    !       ORDER OF INITIAL FUNCTION, FNU.GE.0
    !   n: integer
    !       NUMBER OF MEMBERS OF THE SEQUENCE, N.GE.1
    !
    ! Returns
    ! ----
    !   rv: array<double precision>
    !
    ! LaTeX: $$ h_{n}^{(1)}(\rho)=j_{n}(\rho)+i y_{n}(\rho) $$
    interface lib_math_hankel_spherical_1
        module procedure lib_math_bessel_spherical_third_kind_1_real
        module procedure lib_math_bessel_spherical_third_kind_1_cmplx
    end interface

    ! calculates the the spherical Bessel function of the third kind
    !
    ! symbol: h_n^(2)
    !
    ! Argument
    ! ----
    !   x: double precision .OR. complex
    !   fnu: integer
    !       ORDER OF INITIAL FUNCTION, FNU.GE.0
    !   n: integer
    !       NUMBER OF MEMBERS OF THE SEQUENCE, N.GE.1
    !
    ! Returns
    ! ----
    !   rv: array<double precision>
    !
    !
    ! LaTeX: $$ h_{n}^{(2)}(\rho)=j_{n}(\rho) - i y_{n}(\rho) $$
    interface lib_math_hankel_spherical_2
        module procedure lib_math_bessel_spherical_third_kind_2_real
        module procedure lib_math_bessel_spherical_third_kind_2_cmplx
    end interface

    ! calculates the the Riccati Bessel function
    !
    ! symbol: S
    !
    ! Formula: S = x * j_n(x)
    !
    ! Argument
    ! ----
    !   x: double precision .OR. complex
    !   fnu: integer
    !       ORDER OF INITIAL FUNCTION, FNU.GE.0
    !   n: integer
    !       NUMBER OF MEMBERS OF THE SEQUENCE, N.GE.1
    !
    ! Returns
    ! ----
    !   rv: array<double precision>
    !
    ! LaTeX: $$ S_{n}(z) \equiv z j_{n}(z)=\sqrt{\frac{\pi z}{2}} J_{n+1 / 2}(z) $$
    interface lib_math_riccati_s
        module procedure lib_math_bessel_riccati_s_real
        module procedure lib_math_bessel_riccati_s_cmplx
    end interface

    ! calculates the the Riccat Bessel function
    !
    ! symbol: C
    !
    ! Formula: C = -x * y_n(x)
    !
    ! Argument
    ! ----
    !   x: double precision .OR. complex
    !   fnu: integer
    !       ORDER OF INITIAL FUNCTION, FNU.GE.0
    !   n: integer
    !       NUMBER OF MEMBERS OF THE SEQUENCE, N.GE.1
    !
    ! Returns
    ! ----
    !   rv: array<double precision>
    !
    ! LaTeX: $$ C_{n}(z) \equiv-z n_{n}(z)=-\sqrt{\frac{\pi z}{2}} N_{n+1 / 2}(z) $$
    interface lib_math_riccati_c
        module procedure lib_math_bessel_riccati_c_real
        module procedure lib_math_bessel_riccati_c_cmplx
    end interface

    ! calculates the the Riccat Bessel function
    !
    ! symbol: Xi
    !
    ! Formula: Xi = x * h^(1)_n(x)
    !
    ! Argument
    ! ----
    !   x: double precision .OR. complex
    !   fnu: integer
    !       ORDER OF INITIAL FUNCTION, FNU.GE.0
    !   n: integer
    !       NUMBER OF MEMBERS OF THE SEQUENCE, N.GE.1
    !
    ! Returns
    ! ----
    !   rv: array<double precision>
    !
    ! LaTeX: $$ \xi_{n}(\rho)= \sqrt{\frac{\pi \rho}{2}} J_{n+1 / 2}(\rho) + i \sqrt{\frac{\pi \rho}{2}} Y_{n+1 / 2}(\rho) $$
    interface lib_math_riccati_xi
        module procedure lib_math_bessel_riccati_xi_real
        module procedure lib_math_bessel_riccati_xi_cmplx
    end interface

    ! calculates the the Riccat Bessel function
    !
    ! symbol: Zeta
    !
    ! Formula: Zeta = x * h^(2)_n(x)
    !
    ! Argument
    ! ----
    !   x: double precision .OR.complex
    !   fnu: integer
    !       ORDER OF INITIAL FUNCTION, FNU.GE.0
    !   n: integer
    !       NUMBER OF MEMBERS OF THE SEQUENCE, N.GE.1
    !
    ! Returns
    ! ----
    !   rv: array<double precision>
    !
    ! LaTeX: $$ \zeta_{n}(x)=x h_{n}^{(2)}(x)=\sqrt{\frac{\pi x}{2}} H_{n+\frac{1}{2}}^{(2)}(x)=S_{n}(x)+i C_{n}(x) $$
    interface lib_math_riccati_zeta
        module procedure lib_math_bessel_riccati_zeta_real
        module procedure lib_math_bessel_riccati_zeta_cmplx
    end interface

    ! Calculates the derivative of the spherical Bessel function of the first kind
    !
    ! symbol: j'_n(a x)
    !
    ! Argument
    ! ----
    !   a: integer, optional
    !       coefficient
    !   x: double precision .OR. complex
    !       control variable
    !   fnu: integer
    !       ORDER OF INITIAL FUNCTION, FNU.GE.0
    !   n: integer
    !       NUMBER OF MEMBERS OF THE SEQUENCE, N.GE.1
    !
    ! Returns
    ! ----
    !   j_n: array<double precision>
    !       spherical Bessel function of the first kind [fnu..fnu+n-1]
    !   rv: array<double precision>
    !
    !
    ! LaTeX: $$ \frac{\partial}{\partial x}\left(j_{n}(x)\right)=\frac{1}{2}\left(j_{n-1}(x)-\frac{j_{n}(x)+x j_{n+1}(x)}{x}\right) $$
    interface lib_math_bessel_spherical_first_kind_derivative
        module procedure lib_math_bessel_spherical_first_kind_derivative_real
        module procedure lib_math_bessel_spherical_first_kind_derivative_cmplx
        module procedure lib_math_bessel_spherical_first_kind_derivative_coeff_real
        module procedure lib_math_bessel_spherical_first_kind_derivative_coeff_cmplx
    end interface

    ! Calculates the derivative of the spherical Bessel function of the second kind
    !
    ! symbol: y'_n
    !
    ! Argument
    ! ----
    !   x: double precision
    !   fnu: integer
    !       ORDER OF INITIAL FUNCTION, FNU.GE.0
    !   n: integer
    !       NUMBER OF MEMBERS OF THE SEQUENCE, N.GE.1
    !
    ! Returns
    ! ----
    !   y_n: array<double precision>
    !       spherical Bessel function of the second kind [fnu..fnu+n-1]
    !   rv: array<double precision>
    !
    ! LaTeX: $$ y_{v}(z) \equiv \sqrt{\frac{\pi}{2 z}} Y_{v+1 / 2}(z) $$
    interface lib_math_bessel_spherical_second_kind_derivative
        module procedure lib_math_bessel_spherical_second_kind_derivative_real
        module procedure lib_math_bessel_spherical_second_kind_derivative_cmplx
    end interface

    ! Calculates the derivative of the spherical Bessel function of the third kind 1
    !
    ! symbol: h'^(1)_n
    !
    ! Argument
    ! ----
    !   x: double precision
    !   fnu: integer
    !       ORDER OF INITIAL FUNCTION, FNU.GE.1
    !   n: integer
    !       NUMBER OF MEMBERS OF THE SEQUENCE, N.GE.1
    !
    ! Returns
    ! ----
    !   h_1_n: array<double precision>
    !       spherical Bessel function of the third kind 1 [fnu-1..fnu+n-1]
    !   rv: array<double precision>
    !
    ! LaTeX: $$ \frac{d}{d z} h_{n}^{(1)}(z)=\frac{1}{2}\left[h_{n-1}^{(1)}(z)-\frac{h_{n}^{(1)}(z)+z h_{n+1}^{(1)}(z)}{z}\right] $$
    interface lib_math_hankel_spherical_1_derivative
        module procedure lib_math_bessel_spherical_third_kind_1_derivative_real
        module procedure lib_math_bessel_spherical_third_kind_1_derivative_cmplx
    end interface

    ! Calculates the derivative of the spherical Bessel function of the third kind 1
    !
    ! symbol: h'^(2)_n
    !
    ! Argument
    ! ----
    !   x: double precision
    !   fnu: integer
    !       ORDER OF INITIAL FUNCTION, FNU.GE.1
    !   n: integer
    !       NUMBER OF MEMBERS OF THE SEQUENCE, N.GE.1
    !
    ! Returns
    ! ----
    !   h_2_n: array<double precision>
    !       spherical Bessel function of the third kind 2 [fnu-1..fnu+n-1]
    !   rv: array<double precision>
    !
    ! LaTeX: $$ \frac{d}{d z} h_{n}^{(2)}(z)=\frac{1}{2}\left[h_{n-1}^{(2)}(z)-\frac{h_{n}^{(2)}(z)+z h_{n+1}^{(2)}(z)}{z}\right] $$
    interface lib_math_hankel_spherical_2_derivative
        module procedure lib_math_bessel_spherical_third_kind_2_derivative_real
        module procedure lib_math_bessel_spherical_third_kind_2_derivative_cmplx
    end interface

    ! Calculates the derivative of the Riccati-Bessel function
    !
    ! symbol: S'
    !
    ! Formula: S' = [ x * j_n(x) ]'
    !
    ! Argument
    ! ----
    !   x: double precision
    !   fnu: integer
    !       ORDER OF INITIAL FUNCTION, FNU.GE.1
    !   n: integer
    !       NUMBER OF MEMBERS OF THE SEQUENCE, N.GE.1
    !
    ! Returns
    ! ----
    !   S_n: array<double precision>
    !       Riccati-Bessel function [fnu..fnu+n-1]
    !   rv: array<double precision>
    !
    ! LaTeX: $$  S^\prime_n = -x j_{n+1}(x)+(n+1) j_{n}(x) $$
    interface lib_math_riccati_s_derivative
        module procedure lib_math_bessel_riccati_s_derivative_real
        module procedure lib_math_bessel_riccati_s_derivative_cmplx
    end interface

    ! Calculates the derivative of the Riccati-Bessel function
    !
    ! symbol: C'
    !
    ! Formula: C' = [ -x * y_n(x) ]'
    !
    ! Argument
    ! ----
    !   x: double precision
    !   fnu: integer
    !       ORDER OF INITIAL FUNCTION, FNU.GE.1
    !   n: integer
    !       NUMBER OF MEMBERS OF THE SEQUENCE, N.GE.1
    !
    ! Returns
    ! ----
    !   S_n: array<double precision>
    !       Riccati-Bessel function [fnu..fnu+n-1]
    !   rv: array<double precision>
    !
    ! LaTeX: $$ C^\prime_n = x y_{n+1}(x)-(n+1) y_{n}(x) $$
    interface lib_math_riccati_c_derivative
        module procedure lib_math_bessel_riccati_c_derivative_real
        module procedure lib_math_bessel_riccati_c_derivative_cmplx
    end interface

    ! Calculates the derivative of the Riccati-Bessel function
    !
    ! symbol: Xi'
    !
    ! Formula: Xi' = [ x * h^(1)_n(x) ]'
    !
    ! Argument
    ! ----
    !   x: double precision .OR. compelx
    !   fnu: integer
    !       ORDER OF INITIAL FUNCTION, FNU.GE.1
    !   n: integer
    !       NUMBER OF MEMBERS OF THE SEQUENCE, N.GE.1
    !
    ! Returns
    ! ----
    !   xi_n: array<double precision>
    !       Riccati-Bessel function [fnu..fnu+n-1]
    !   rv: array<double precision>
    !
    ! LaTeX: $$ \xi^\prime_n (x) = \frac{\partial}{\partial x}\left(x h_{n}^{(1)}(x)\right)=\frac{1}{2}\left(x h_{n-1}^{(1)}(x)+h_{n}^{(1)}(x)-x h_{n+1}^{(1)}(x)\right) $$
    interface lib_math_riccati_xi_derivative
        module procedure lib_math_bessel_riccati_xi_derivative_real
        module procedure lib_math_bessel_riccati_xi_derivative_cmplx
    end interface

    ! Calculates the derivative of the Riccati-Bessel function
    !
    ! symbol: Zeta'
    !
    ! Formula: Zeta' = [ x * h^(2)_n(x) ]'
    !
    ! Argument
    ! ----
    !   x: double precision .OR. complex
    !   fnu: integer
    !       ORDER OF INITIAL FUNCTION, FNU.GE.1
    !   n: integer
    !       NUMBER OF MEMBERS OF THE SEQUENCE, N.GE.1
    !
    ! Returns
    ! ----
    !   xi_n: array<double precision>
    !       Riccati-Bessel function [fnu..fnu+n-1]
    !   rv: array<double precision>
    !
    ! LaTeX: $$ \zeta^\prime_n (x) = \frac{\partial}{\partial x}\left(x h_{n}^{(2)}(x)\right)=\frac{1}{2}\left(x h_{n-1}^{(2)}(x)+h_{n}^{(2)}(x)-x h_{n+1}^{(2)}(x)\right) $$
    interface lib_math_riccati_zeta_derivative
        module procedure lib_math_bessel_riccati_zeta_derivative_real
        module procedure lib_math_bessel_riccati_zeta_derivative_cmplx
    end interface

    ! --- public functions ---
    public :: lib_math_bessel_spherical_first_kind
    public :: lib_math_bessel_spherical_second_kind
    public :: lib_math_hankel_spherical_1
    public :: lib_math_hankel_spherical_2
    public :: lib_math_bessel_spherical_first_kind_derivative
    public :: lib_math_bessel_spherical_second_kind_derivative
    public :: lib_math_hankel_spherical_1_derivative
    public :: lib_math_hankel_spherical_2_derivative
    public :: lib_math_riccati_s_derivative
    public :: lib_math_riccati_c_derivative
    public :: lib_math_riccati_xi_derivative
    public :: lib_math_riccati_zeta_derivative

    public :: lib_math_bessel_test_functions

    ! --- paraemters ---
    double precision, parameter :: PI=4.D0*atan(1.D0)   ! maximum precision, platform independet

    contains

    ! calculates the the spherical Bessel function of the first kind
    !
    ! symbol: j_n
    !
    ! Argument
    ! ----
    !   x: double precision
    !   fnu: integer
    !       ORDER OF INITIAL FUNCTION, FNU.GE.0
    !   n: integer
    !       NUMBER OF MEMBERS OF THE SEQUENCE, N.GE.1
    !
    ! Returns
    ! ----
    !   rv: array<double precision>
    !
    ! LaTeX: $$ j_{n}(\rho)=\sqrt{\frac{\pi}{2 \rho}} J_{n+1 / 2}(\rho) $$
    !
    ! Reference: Absorption and Scattering of Light by Small Particles Absorption and Scattering of Light by Small Particles (4.9)
    !
    function lib_math_bessel_spherical_first_kind_real(x, fnu, n) result (rv)
        implicit none
        ! dummy
        double precision, intent(in) :: x
        integer(kind=4), intent(in) :: fnu
        integer(kind=4), intent(in) :: n

        double precision, dimension(n) :: rv

        ! auxiliary
        integer(kind=4) order
        integer(kind=4) nm
        real(kind=8) dj(0:fnu+n-1)
        real(kind=8) sj(0:fnu+n-1)

        order = fnu+n-1

        call SPHJ( order, x, nm, sj, dj)

        if (order .ne. nm) then
            print *, "lib_math_bessel_spherical_first_kind_real: ERROR"
            print *, "  calculated highest order / requested: ", nm, " / ", order
        end if

        rv = sj(fnu:order)

    end function lib_math_bessel_spherical_first_kind_real

    ! calculates the the spherical Bessel function of the first kind
    !
    ! symbol: j_n
    !
    ! Argument
    ! ----
    !   x: complex
    !   fnu: integer
    !       ORDER OF INITIAL FUNCTION, FNU.GE.0
    !   n: integer
    !       NUMBER OF MEMBERS OF THE SEQUENCE, N.GE.1
    !   kode: integer, optional (std_value = 1)
    !       A PARAMETER TO INDICATE THE SCALING OPTION
    !           KODE= 1  RETURNS
    !                    CY(I)=J(FNU+I-1,Z), I=1,...,N
    !               = 2  RETURNS
    !                    CY(I)=J(FNU+I-1,Z)EXP(-ABS(Y)), I=1,...,N
    !
    ! Returns
    ! ----
    !   rv: array<cmplx>
    !
    ! LaTeX: $$ j_{n}(\rho)=\sqrt{\frac{\pi}{2 \rho}} J_{n+1 / 2}(\rho) $$
    !
    ! Reference: Absorption and Scattering of Light by Small Particles Absorption and Scattering of Light by Small Particles (4.9)
    function lib_math_bessel_spherical_first_kind_cmplx(z, fnu, n) result (rv)
        implicit none
        ! dummy
        integer(kind=4), intent(in) :: fnu
        integer(kind=4), intent(in) :: n
        complex(kind=8), intent(in) :: z

        complex(kind=8), dimension(n) :: rv

        ! auxiliary
        integer(kind=4) order
        integer(kind=4) nm
        complex(kind=8) dj(0:fnu+n-1)
        complex(kind=8) sj(0:fnu+n-1)
        complex(kind=8) dy(0:fnu+n-1)
        complex(kind=8) sy(0:fnu+n-1)

        order = fnu+n-1

        call CSPHJY( order, z, nm, sj, dj, sy, dy)

        if (order .ne. nm) then
            print *, "lib_math_bessel_spherical_first_kind_cmplx: ERROR"
            print *, "  calculated highest order / requested: ", nm, " / ", order
        end if

        rv = sj(fnu:order)

    end function lib_math_bessel_spherical_first_kind_cmplx


    ! calculates the the spherical Bessel function of the second kind
    !
    ! symbol: y_n
    !
    ! Argument
    ! ----
    !   x: double precision
    !   fnu: integer
    !       ORDER OF INITIAL FUNCTION, FNU.GE.0
    !   n: integer
    !       NUMBER OF MEMBERS OF THE SEQUENCE, N.GE.1
    !
    ! Returns
    ! ----
    !   rv: array<double precision>
    !
    ! LaTeX: $$ y_{n}(\rho)=\sqrt{\frac{\pi}{2 \rho}} Y_{n+1 / 2}(\rho) $$
    !
    ! Reference: Absorption and Scattering of Light by Small Particles Absorption and Scattering of Light by Small Particles (4.10)
    !
    function lib_math_bessel_spherical_second_kind_real(x, fnu, n) result (rv)
        implicit none
        ! dummy
        double precision, intent(in) :: x
        integer(kind=4), intent(in) :: fnu
        integer(kind=4), intent(in) :: n

        double precision, dimension(n) :: rv

        ! auxiliary
        integer(kind=4) order
        integer(kind=4) nm
        real(kind=8) dy(0:fnu+n-1)
        real(kind=8) sy(0:fnu+n-1)

        order = fnu+n-1

        call SPHY( order, x, nm, sy, dy)

        if (order .ne. nm) then
            print *, "lib_math_bessel_spherical_second_kind_real: ERROR"
            print *, "  calculated highest order / requested: ", nm, " / ", order
        end if

        rv = sy(fnu:order)

    end function lib_math_bessel_spherical_second_kind_real

    ! calculates the the spherical Bessel function of the first kind
    !
    ! symbol: y_n
    !
    ! Argument
    ! ----
    !   x: complex
    !   fnu: integer
    !       ORDER OF INITIAL FUNCTION, FNU.GE.0
    !   n: integer
    !       NUMBER OF MEMBERS OF THE SEQUENCE, N.GE.1
    !   kode: integer, optional (std_value = 1)
    !       A PARAMETER TO INDICATE THE SCALING OPTION
    !           KODE= 1  RETURNS
    !                    CY(I)=Y(FNU+I-1,Z), I=1,...,N
    !               = 2  RETURNS
    !                    CY(I)=Y(FNU+I-1,Z)EXP(-ABS(Y)), I=1,...,N
    !
    ! Returns
    ! ----
    !   rv: array<cmplx>
    !
    ! LaTeX: $$ j_{n}(\rho)=\sqrt{\frac{\pi}{2 \rho}} J_{n+1 / 2}(\rho) $$
    !
    ! Reference: Absorption and Scattering of Light by Small Particles Absorption and Scattering of Light by Small Particles (4.9)
    function lib_math_bessel_spherical_second_kind_cmplx(z, fnu, n) result (rv)
    !
        implicit none
        ! dummy
        integer(kind=4), intent(in) :: fnu
        integer(kind=4), intent(in) :: n
        complex(kind=8), intent(in) :: z

        complex(kind=8), dimension(n) :: rv

        ! auxiliary
        integer(kind=4) order
        integer(kind=4) nm
        complex(kind=8) dj(0:fnu+n-1)
        complex(kind=8) sj(0:fnu+n-1)
        complex(kind=8) dy(0:fnu+n-1)
        complex(kind=8) sy(0:fnu+n-1)

        order = fnu+n-1

        call CSPHJY( order, z, nm, sj, dj, sy, dy)

        if (order .ne. nm) then
            print *, "lib_math_bessel_spherical_second_kind_cmplx: ERROR"
            print *, "  calculated highest order / requested: ", nm, " / ", order
        end if

        rv = sy(fnu:order)

    end function lib_math_bessel_spherical_second_kind_cmplx

    ! calculates the the spherical Bessel function of the third kind
    !
    ! symbol: h_n^(1)
    !
    ! Argument
    ! ----
    !   x: double precision
    !       x .gt. 0.0
    !   fnu: integer
    !       ORDER OF INITIAL FUNCTION, FNU.GE.0
    !   n: integer
    !       NUMBER OF MEMBERS OF THE SEQUENCE, N.GE.1
    !
    ! Returns
    ! ----
    !   rv: array<double precision>
    !
    !
    ! LaTeX: $$ h_{n}^{(1)}(z)=i^{-n-1} z^{-1} e^{i z} \sum_{0}^{n}\left(n+\frac{1}{2}, k\right)(-2 i z)^{-k} $$
    !        $$ \left(n+\frac{1}{2}, k\right)=\frac{(n+k) !}{k ! \Gamma(n-k+1)} $$
    !
    ! Reference: Handbook of Mathematical Functions, Abramowitz and Stegun, pp. 437-442, eq. 10.1.16
    !
    function lib_math_bessel_spherical_third_kind_1_real(x, fnu, n) result (rv)
        implicit none
        ! dummy
        double precision, intent(in) :: x
        integer(kind=4), intent(in) :: fnu
        integer(kind=4), intent(in) :: n
        complex(kind=8), dimension(n) :: rv

        ! auxiliary
!        integer(kind=4) :: i
!        integer(kind=4) :: k
!        integer(kind=4) :: order
!        real(kind=8) :: nk
!        complex(kind=8) :: prefactor
!
!        complex(kind=8), dimension(fnu:fnu+n-1) :: i_pow_n
!        complex(kind=8), dimension(0:fnu+n-1) :: x_pow_k
!        complex(kind=8) :: summation
!
!        complex(kind=8), dimension(fnu:fnu+n-1) :: h
!
!
!        ! pre-calulations
!        order = fnu + n - 1
!
!        prefactor = cmplx(cos(x), sin(x), kind=8) / x
!
!        x_pow_k(0) = cmplx(1.0, 0.0, kind=8)
!        do k=1, order
!            x_pow_k(k) = x_pow_k(k-1) * cmplx(0.0_8, -2.0_8 * x, kind=8)
!        end do
!
!        do i=fnu, order
!            i_pow_n(i) = cmplx(0.0, 1.0, kind=8)**(-i-1)
!        end do
!
!
!        ! calculate Hankel functions of degree fnu until order
!        do i=fnu, order
!            ! calculate the spherical Hankel function with degree i
!            summation = cmplx(0.0, 0.0, kind=8)
!            do k=0, i
!                nk = lib_math_factorial_get_n_plus_m_divided_by_m_fac_n_minus_m(i, k)
!                summation = summation + nk / x_pow_k(k)
!            end do
!
!            h(i) = i_pow_n(i) * prefactor * summation
!
!        end do
!
!        rv = h

        rv = cmplx(lib_math_bessel_spherical_first_kind_real(x, fnu, n), &
                   lib_math_bessel_spherical_second_kind_real(x, fnu, n), &
                   kind=8)

    end function lib_math_bessel_spherical_third_kind_1_real

    ! calculates the the spherical Bessel function of the third kind
    !
    ! symbol: h_n^(1)
    !
    ! Argument
    ! ----
    !   x: cmplx
    !       x .gt. 0.0
    !   fnu: integer
    !       ORDER OF INITIAL FUNCTION, FNU.GE.0
    !   n: integer
    !       NUMBER OF MEMBERS OF THE SEQUENCE, N.GE.1
    !
    ! Returns
    ! ----
    !   rv: array<double precision>
    !
    !
    ! LaTeX: $$ h_{n}^{(1)}(z)=i^{-n-1} z^{-1} e^{i z} \sum_{0}^{n}\left(n+\frac{1}{2}, k\right)(-2 i z)^{-k} $$
    !        $$ \left(n+\frac{1}{2}, k\right)=\frac{(n+k) !}{k ! \Gamma(n-k+1)} $$
    !
    ! Reference: Handbook of Mathematical Functions, Abramowitz and Stegun, pp. 437-442, eq. 10.1.16
    !
    function lib_math_bessel_spherical_third_kind_1_cmplx(x, fnu, n) result (rv)
        implicit none
        ! dummy
        complex(kind=8), intent(in) :: x
        integer(kind=4), intent(in) :: fnu
        integer(kind=4), intent(in) :: n
        complex(kind=8), dimension(n) :: rv

        ! auxiliary
!        integer(kind=4) :: i
!        integer(kind=4) :: k
!        integer(kind=4) :: order
!        real(kind=8) :: nk
!        complex(kind=8) :: prefactor
!
!        complex(kind=8), dimension(fnu:fnu+n-1) :: i_pow_n
!        complex(kind=8), dimension(0:fnu+n-1) :: x_pow_k
!        complex(kind=8) :: summation
!
!        complex(kind=8), dimension(fnu:fnu+n-1) :: h
!
!
!        ! pre-calulations
!        order = fnu + n - 1
!
!        prefactor = exp(cmplx(0.0, 1.0, kind=8) * x) / x
!
!        x_pow_k(0) = cmplx(1.0, 0.0, kind=8)
!        do k=1, order
!            x_pow_k(k) = x_pow_k(k-1) * (-2.0_8) * x * cmplx(0.0, 1.0, kind=8)
!        end do
!
!        do i=fnu, order
!            i_pow_n(i) = cmplx(0.0, 1.0, kind=8)**(-i-1)
!        end do
!
!
!        ! calculate Hankel functions of degree fnu until order
!        do i=fnu, order
!            ! calculate the spherical Hankel function with degree i
!            summation = cmplx(0.0, 0.0, kind=8)
!            do k=0, i
!                nk = lib_math_factorial_get_n_plus_m_divided_by_m_fac_n_minus_m(i, k)
!                summation = summation + nk / x_pow_k(k)
!            end do
!
!            h(i) = i_pow_n(i) * prefactor * summation
!
!        end do
!
!        rv = h

        rv = lib_math_bessel_spherical_first_kind_cmplx(x, fnu, n) &
             + cmplx(0,1, kind=8) * lib_math_bessel_spherical_second_kind_cmplx(x, fnu, n)

    end function lib_math_bessel_spherical_third_kind_1_cmplx

    ! calculates the the spherical Bessel function of the third kind
    !
    ! symbol: h_n^(2)
    !
    ! Argument
    ! ----
    !   x: double precision
    !       x .gt. 0.0
    !   fnu: integer
    !       ORDER OF INITIAL FUNCTION, FNU.GE.0
    !   n: integer
    !       NUMBER OF MEMBERS OF THE SEQUENCE, N.GE.1
    !
    ! Returns
    ! ----
    !   rv: array<double precision>
    !
    !
    ! LaTeX: $$ h_{\mathrm{n}}^{(2)}(z)=i^{n+1} z^{-1} e^{-i z} \sum_{0}^{n}\left(n+\frac{1}{2}, k\right)(2 i z)^{-k} $$
    !        $$ \left(n+\frac{1}{2}, k\right)=\frac{(n+k) !}{k ! \Gamma(n-k+1)} $$
    !
    ! Reference: Handbook of Mathematical Functions, Abramowitz and Stegun, pp. 437-442, eq. 10.1.17
    !
    function lib_math_bessel_spherical_third_kind_2_real(x, fnu, n) result (rv)
        implicit none
        ! dummy
        double precision, intent(in) :: x
        integer(kind=4), intent(in) :: fnu
        integer(kind=4), intent(in) :: n
        complex(kind=8), dimension(n) :: rv

        ! auxiliary
!        integer(kind=4) :: i
!        integer(kind=4) :: k
!        integer(kind=4) :: order
!        real(kind=8) :: nk
!        complex(kind=8) :: prefactor
!
!        complex(kind=8), dimension(fnu:fnu+n-1) :: i_pow_n
!        complex(kind=8), dimension(0:fnu+n-1) :: x_pow_k
!        complex(kind=8) :: summation
!
!        complex(kind=8), dimension(fnu:fnu+n-1) :: h
!
!
!        ! pre-calulations
!        order = fnu + n - 1
!
!        prefactor = cmplx(cos(x), -sin(x), kind=8) / x
!
!        x_pow_k(0) = cmplx(1.0, 0.0, kind=8)
!        do k=1, order
!            x_pow_k(k) = x_pow_k(k-1) * cmplx(0.0_8, 2.0_8 * x, kind=8)
!        end do
!
!        do i=fnu, order
!            i_pow_n(i) = cmplx(0.0, 1.0, kind=8)**(i+1)
!        end do
!
!
!        ! calculate Hankel functions of degree fnu until order
!        do i=fnu, order
!            ! calculate the spherical Hankel function with degree i
!            summation = cmplx(0.0, 0.0, kind=8)
!            do k=0, i
!                nk = lib_math_factorial_get_n_plus_m_divided_by_m_fac_n_minus_m(i, k)
!                summation = summation + nk / x_pow_k(k)
!            end do
!
!            h(i) = i_pow_n(i) * prefactor * summation
!
!        end do
!
!        rv = h

        rv = cmplx(lib_math_bessel_spherical_first_kind_real(x, fnu, n), &
                   -lib_math_bessel_spherical_second_kind_real(x, fnu, n), &
                   kind=8)

    end function lib_math_bessel_spherical_third_kind_2_real

    ! calculates the the spherical Bessel function of the third kind
    !
    ! symbol: h_n^(2)
    !
    ! Argument
    ! ----
    !   z: complex
    !       z .ne. (0.0, 0.0)
    !   fnu: integer
    !       ORDER OF INITIAL FUNCTION, FNU.GE.0
    !   n: integer
    !       NUMBER OF MEMBERS OF THE SEQUENCE, N.GE.1
    !
    ! Returns
    ! ----
    !   rv: array<double precision>
    !
    !
    ! LaTeX: $$ h_{\mathrm{n}}^{(2)}(z)=i^{n+1} z^{-1} e^{-i z} \sum_{0}^{n}\left(n+\frac{1}{2}, k\right)(2 i z)^{-k} $$
    !        $$ \left(n+\frac{1}{2}, k\right)=\frac{(n+k) !}{k ! \Gamma(n-k+1)} $$
    !
    ! Reference: Handbook of Mathematical Functions, Abramowitz and Stegun, pp. 437-442, eq. 10.1.17
    !
    function lib_math_bessel_spherical_third_kind_2_cmplx(z, fnu, n) result (rv)
        implicit none
        ! dummy
        complex(kind=8), intent(in) :: z
        integer(kind=4), intent(in) :: fnu
        integer(kind=4), intent(in) :: n
        complex(kind=8), dimension(n) :: rv

        ! auxiliary
!        integer(kind=4) :: i
!        integer(kind=4) :: k
!        integer(kind=4) :: order
!        real(kind=8) :: nk
!        complex(kind=8) :: prefactor
!
!        complex(kind=8), dimension(fnu:fnu+n-1) :: i_pow_n
!        complex(kind=8), dimension(0:fnu+n-1) :: x_pow_k
!        complex(kind=8) :: summation
!
!        complex(kind=8), dimension(fnu:fnu+n-1) :: h
!
!
!        ! pre-calulations
!        order = fnu + n - 1
!
!        prefactor = exp(-cmplx(0.0, 1.0, kind=8) * z) / z
!
!        x_pow_k(0) = cmplx(1.0, 0.0, kind=8)
!        do k=1, order
!            x_pow_k(k) = x_pow_k(k-1) * (2.0_8) * z * cmplx(0.0, 1.0, kind=8)
!        end do
!
!        do i=fnu, order
!            i_pow_n(i) = cmplx(0.0, 1.0, kind=8)**(i+1)
!        end do
!
!
!        ! calculate Hankel functions of degree fnu until order
!        do i=fnu, order
!            ! calculate the spherical Hankel function with degree i
!            summation = cmplx(0.0, 0.0, kind=8)
!            do k=0, i
!                nk = lib_math_factorial_get_n_plus_m_divided_by_m_fac_n_minus_m(i, k)
!                summation = summation + nk / x_pow_k(k)
!            end do
!
!            h(i) = i_pow_n(i) * prefactor * summation
!
!        end do
!
!        rv = h

        rv = lib_math_bessel_spherical_first_kind_cmplx(z, fnu, n) &
             - cmplx(0,1, kind=8) * lib_math_bessel_spherical_second_kind_cmplx(z, fnu, n)

    end function lib_math_bessel_spherical_third_kind_2_cmplx

    ! calculates the the Riccati Bessel function
    !
    ! symbol: S
    !
    ! Argument
    ! ----
    !   x: double precision
    !   fnu: integer
    !       ORDER OF INITIAL FUNCTION, FNU.GE.0
    !   n: integer
    !       NUMBER OF MEMBERS OF THE SEQUENCE, N.GE.1
    !
    ! Returns
    ! ----
    !   rv: array<double precision>
    !
    ! LaTeX: $$ S_{n}(z) \equiv z j_{n}(z)=\sqrt{\frac{\pi z}{2}} J_{n+1 / 2}(z) $$
    !
    ! Reference: http://mathworld.wolfram.com/Riccati-BesselFunctions.html
    !
    function lib_math_bessel_riccati_s_real(x, fnu, n) result (rv)
        implicit none
        ! dummy
        double precision, intent(in) :: x
        integer(kind=4), intent(in) :: fnu
        integer(kind=4), intent(in) :: n

        double precision, dimension(n) :: rv

        ! auxiliary
        integer(kind=4) order
        integer(kind=4) nm
        real(kind=8) rf(0:fnu+n-1)
        real(kind=8) df(0:fnu+n-1)

        order = fnu+n-1

        call rctj( order, x, nm, rf, df)

        if (order .ne. nm) then
            print *, "lib_math_bessel_riccati_s_real: ERROR"
            print *, "  calculated highest order / requested: ", nm, " / ", order
        end if

        rv = rf(fnu:order)

    end function lib_math_bessel_riccati_s_real

    ! calculates the the Riccati Bessel function
    !
    ! symbol: S
    !
    ! Argument
    ! ----
    !   x: complex
    !   fnu: integer
    !       ORDER OF INITIAL FUNCTION, FNU.GE.0
    !   n: integer
    !       NUMBER OF MEMBERS OF THE SEQUENCE, N.GE.1
    !
    ! Returns
    ! ----
    !   rv: array<double precision>
    !
    ! LaTeX: $$ S_{n}(z) \equiv z j_{n}(z)=\sqrt{\frac{\pi z}{2}} J_{n+1 / 2}(z) $$
    !
    ! Reference: http://mathworld.wolfram.com/Riccati-BesselFunctions.html
    !
    function lib_math_bessel_riccati_s_cmplx(z, fnu, n) result (rv)
        implicit none
        ! dummy
        complex(kind=8), intent(in) :: z
        integer(kind=4), intent(in) :: fnu
        integer(kind=4), intent(in) :: n

        complex(kind=8), dimension(n) :: rv

        rv = z * lib_math_bessel_spherical_first_kind_cmplx(z, fnu, n)

    end function lib_math_bessel_riccati_s_cmplx

    ! calculates the the Riccat Bessel function
    !
    ! symbol: C
    !
    ! Argument
    ! ----
    !   x: double precision
    !   fnu: integer
    !       ORDER OF INITIAL FUNCTION, FNU.GE.0
    !   n: integer
    !       NUMBER OF MEMBERS OF THE SEQUENCE, N.GE.1
    !
    ! Returns
    ! ----
    !   rv: array<double precision>
    !
    ! LaTeX: $$ C_{n}(z) \equiv-z n_{n}(z)=-\sqrt{\frac{\pi z}{2}} N_{n+1 / 2}(z) $$
    !
    ! Reference: http://mathworld.wolfram.com/Riccati-BesselFunctions.html
    !
    function lib_math_bessel_riccati_c_real(x, fnu, n) result (rv)
        implicit none
        ! dummy
        double precision, intent(in) :: x
        integer(kind=4), intent(in) :: fnu
        integer(kind=4), intent(in) :: n

        double precision, dimension(n) :: rv

        ! auxiliary
        integer(kind=4) order
        integer(kind=4) nm
        real(kind=8) rf(0:fnu+n-1)
        real(kind=8) df(0:fnu+n-1)

        order = fnu+n-1

        call rcty( order, x, nm, rf, df)

        if (order .ne. nm) then
            print *, "lib_math_bessel_riccati_c_real: ERROR"
            print *, "  calculated highest order / requested: ", nm, " / ", order
        end if

        ! tests with wolfram alpha and sageMath results yield a prefactor of "-1"
        ! wolfram alphs
        !   >>> riccatiBesselC[2, 4]
        ! sageMath
        !   >>> C(n,x) = -x*spherical_bessel_Y(n,x)
        !   >>> numerical_approx(C(2,4))
        rv = -rf(fnu:order)

    end function lib_math_bessel_riccati_c_real

    ! calculates the the Riccat Bessel function
    !
    ! symbol: C
    !
    ! Argument
    ! ----
    !   x: double precision
    !   fnu: integer
    !       ORDER OF INITIAL FUNCTION, FNU.GE.0
    !   n: integer
    !       NUMBER OF MEMBERS OF THE SEQUENCE, N.GE.1
    !
    ! Returns
    ! ----
    !   rv: array<double precision>
    !
    ! LaTeX: $$ C_{n}(z) \equiv-z n_{n}(z)=-\sqrt{\frac{\pi z}{2}} N_{n+1 / 2}(z) $$
    !
    ! Reference: http://mathworld.wolfram.com/Riccati-BesselFunctions.html
    !
    function lib_math_bessel_riccati_c_cmplx(z, fnu, n) result (rv)
        implicit none
        ! dummy
        complex(kind=8), intent(in) :: z
        integer(kind=4), intent(in) :: fnu
        integer(kind=4), intent(in) :: n

        complex(kind=8), dimension(n) :: rv

        rv = -z * lib_math_bessel_spherical_second_kind_cmplx(z, fnu, n)

    end function lib_math_bessel_riccati_c_cmplx

    ! calculates the the Riccat Bessel function
    !
    ! symbol: Xi
    !
    ! Argument
    ! ----
    !   x: double precision
    !   fnu: integer
    !       ORDER OF INITIAL FUNCTION, FNU.GE.0
    !   n: integer
    !       NUMBER OF MEMBERS OF THE SEQUENCE, N.GE.1
    !
    ! Returns
    ! ----
    !   rv: array<double precision>
    !
    ! LaTeX: $$ \xi_{n}(\rho)= \sqrt{\frac{\pi \rho}{2}} J_{n+1 / 2}(\rho) + i \sqrt{\frac{\pi \rho}{2}} Y_{n+1 / 2}(\rho) $$
    !
    ! Reference: Absorption and Scattering of Light by Small Particles Absorption and Scattering of Light by Small Particles (page 101)
    !
    function lib_math_bessel_riccati_xi_real(x, fnu, n) result (rv)
        implicit none
        ! dummy
        double precision, intent(in) :: x
        integer(kind=4), intent(in) :: fnu
        integer(kind=4), intent(in) :: n

        complex(kind=8), dimension(n) :: rv

        ! auxiliary
        double precision, dimension(n) :: j
        double precision, dimension(n) :: y
        double precision :: rv_1

        ! nz: integer
        !    number of components of Y set to zero due to
        !    underflow,
        !    NZ=0   , normal return, computation completed
        !    NZ .NE. 0, last NZ components of Y set to zero,
        !             Y(K)=0.0D0, K=N-NZ+1,...,N.
        integer :: nz

        rv_1 = sqrt(PI/2.D0 * x)

        call DBESJ (X, FNU+0.5D0, N, j, nz)

#ifdef _DEBUG_
        if (nz .ne. 0) then
            print *, "lib_math_bessel_riccati_s_real: WARNING"
            print *, "  number of components of Y set to zero due to underflow: ", nz
        end if
#endif

        call DBESY (X, FNU+0.5D0, N, y)

        rv(:) = cmplx(rv_1 * j(:), rv_1 * y(:), kind=8)

    end function lib_math_bessel_riccati_xi_real

    ! calculates the the Riccat Bessel function
    !
    ! symbol: Xi
    !
    ! Argument
    ! ----
    !   z: complex
    !   fnu: integer
    !       ORDER OF INITIAL FUNCTION, FNU.GE.0
    !   n: integer
    !       NUMBER OF MEMBERS OF THE SEQUENCE, N.GE.1
    !
    ! Returns
    ! ----
    !   rv: array<double precision>
    !
    ! LaTeX: $$ \xi_{n}(\rho)= \sqrt{\frac{\pi \rho}{2}} J_{n+1 / 2}(\rho) + i \sqrt{\frac{\pi \rho}{2}} Y_{n+1 / 2}(\rho) $$
    !
    ! Reference: Absorption and Scattering of Light by Small Particles Absorption and Scattering of Light by Small Particles (page 101)
    !
    function lib_math_bessel_riccati_xi_cmplx(z, fnu, n) result (rv)
        implicit none
        ! dummy
        complex(kind=8), intent(in) :: z
        integer(kind=4), intent(in) :: fnu
        integer(kind=4), intent(in) :: n

        complex(kind=8), dimension(n) :: rv

        rv = z * lib_math_bessel_spherical_third_kind_1_cmplx(z, fnu, n)

    end function lib_math_bessel_riccati_xi_cmplx

    ! calculates the the Riccat Bessel function
    !
    ! symbol: Zeta
    !
    ! Argument
    ! ----
    !   x: double precision
    !   fnu: integer
    !       ORDER OF INITIAL FUNCTION, FNU.GE.0
    !   n: integer
    !       NUMBER OF MEMBERS OF THE SEQUENCE, N.GE.1
    !
    ! Returns
    ! ----
    !   rv: array<double precision>
    !
    ! LaTeX: $$ \zeta_{n}(x)=x h_{n}^{(2)}(x)=\sqrt{\frac{\pi x}{2}} H_{n+\frac{1}{2}}^{(2)}(x)=S_{n}(x)+i C_{n}(x) $$
    !
    ! Reference:
    !
    function lib_math_bessel_riccati_zeta_real(x, fnu, n) result (rv)
        implicit none
        ! dummy
        double precision, intent(in) :: x
        integer(kind=4), intent(in) :: fnu
        integer(kind=4), intent(in) :: n

        complex(kind=8), dimension(n) :: rv

        ! auxiliary
        double precision, dimension(n) :: j
        double precision, dimension(n) :: y
        double precision :: rv_1

        ! nz: integer
        !    number of components of Y set to zero due to
        !    underflow,
        !    NZ=0   , normal return, computation completed
        !    NZ .NE. 0, last NZ components of Y set to zero,
        !             Y(K)=0.0D0, K=N-NZ+1,...,N.
        integer :: nz

        rv_1 = sqrt(PI * x/2.D0)

        call DBESJ (X, FNU+0.5D0, N, j, nz)

#ifdef _DEBUG_
        if (nz .ne. 0) then
            print *, "lib_math_bessel_riccati_zeta_real: WARNING"
            print *, "  number of components of Y set to zero due to underflow: ", nz
        end if
#endif

        call DBESY (X, FNU+0.5D0, N, y)

        rv(:) = cmplx(rv_1 * j(:), -rv_1 * y(:), kind=8)

    end function lib_math_bessel_riccati_zeta_real

    ! calculates the the Riccat Bessel function
    !
    ! symbol: Zeta
    !
    ! Argument
    ! ----
    !   z: complex
    !   fnu: integer
    !       ORDER OF INITIAL FUNCTION, FNU.GE.0
    !   n: integer
    !       NUMBER OF MEMBERS OF THE SEQUENCE, N.GE.1
    !
    ! Returns
    ! ----
    !   rv: array<double precision>
    !
    ! LaTeX: $$ \zeta_{n}(x)=x h_{n}^{(2)}(x)=\sqrt{\frac{\pi x}{2}} H_{n+\frac{1}{2}}^{(2)}(x)=S_{n}(x)+i C_{n}(x) $$
    !
    ! Reference:
    !
    function lib_math_bessel_riccati_zeta_cmplx(z, fnu, n) result (rv)
        implicit none
        ! dummy
        complex(kind=8), intent(in) :: z
        integer(kind=4), intent(in) :: fnu
        integer(kind=4), intent(in) :: n

        complex(kind=8), dimension(n) :: rv

        rv = z * lib_math_bessel_spherical_third_kind_2_cmplx(z, fnu, n)

    end function lib_math_bessel_riccati_zeta_cmplx

    ! Calculates the derivative of the spherical Bessel function of the first kind
    !
    ! symbol: j'_n
    !
    ! Argument
    ! ----
    !   x: double precision
    !   fnu: integer
    !       ORDER OF INITIAL FUNCTION, FNU.GE.0
    !   n: integer
    !       NUMBER OF MEMBERS OF THE SEQUENCE, N.GE.1
    !
    ! Returns
    ! ----
    !   j_n: array<double precision>
    !       spherical Bessel function of the first kind [fnu..fnu+n-1]
    !   rv: array<double precision>
    !
    !
    ! LaTeX: $$ \frac{\partial}{\partial x}\left(j_{n}(x)\right)=\frac{1}{2}\left(j_{n-1}(x)-\frac{j_{n}(x)+x j_{n+1}(x)}{x}\right) $$
    !
    ! Reference: WolframAlpha: https://www.wolframalpha.com/input/?i=derivative%5BSphericalBesselJ%5Bn,x%5D,+x%5D
    !               >>> derivative[SphericalBesselJ[n,x], x]
    function lib_math_bessel_spherical_first_kind_derivative_real(x, fnu, n, j_n) result (rv)
        implicit none
        ! dummy
        double precision, intent(in) :: x
        integer(kind=4), intent(in) :: fnu
        integer(kind=4), intent(in) :: n
        double precision, dimension(n), intent(inout) :: j_n
        double precision, dimension(n) :: rv

        ! auxiliary
        integer(kind=4) order
        integer(kind=4) nm
        real(kind=8) dj(0:fnu+n-1)
        real(kind=8) sj(0:fnu+n-1)

        order = fnu+n-1

        call SPHJ( order, x, nm, sj, dj)

        if (order .ne. nm) then
            print *, "lib_math_bessel_spherical_first_kind_derivative_real: ERROR"
            print *, "  calculated highest order / requested: ", nm, " / ", order
        end if

        j_n = sj(fnu:order)
        rv = dj(fnu:order)

    end function lib_math_bessel_spherical_first_kind_derivative_real

    ! Calculates the derivative of the spherical Bessel function of the first kind
    !
    ! symbol: j'_n(a x)
    !
    ! Argument
    ! ----
    !   a: double precision
    !       coefficient
    !   x: double precision
    !       control variable
    !   fnu: integer
    !       ORDER OF INITIAL FUNCTION, FNU.GE.0
    !   n: integer
    !       NUMBER OF MEMBERS OF THE SEQUENCE, N.GE.1
    !
    ! Returns
    ! ----
    !   j_n: array<double precision>
    !       spherical Bessel function of the first kind [fnu..fnu+n-1]
    !   rv: array<double precision>
    !
    !
    ! LaTeX: $$ \frac{\partial}{\partial x}\left(j_{n}(a x)\right) = \frac{n j_n(a x)}{x}-a j_{n+1}(a x) $$
    !
    ! Mathematica Source Code:
    !   >>> jdn2[n_, x_, t_] := FullSimplify[D[SphericalBesselJ[n, t y], y] /. y -> x];
    function lib_math_bessel_spherical_first_kind_derivative_coeff_real(a, x, fnu, n, j_n) result (rv)
        implicit none
        ! dummy
        double precision, intent(in) :: a
        double precision, intent(in) :: x
        integer(kind=4), intent(in) :: fnu
        integer(kind=4), intent(in) :: n
        double precision, dimension(n), intent(inout) :: j_n
        double precision, dimension(n) :: rv

        ! auxiliary
        integer(kind=4) order
        integer(kind=4) nm
        real(kind=8) dj(0:fnu+n)
        real(kind=8) sj(0:fnu+n)

        order = fnu+n-1

        call SPHJ( order, a * x, nm, sj, dj)

        if (order .ne. nm) then
            print *, "lib_math_bessel_spherical_first_kind_derivative_real: ERROR"
            print *, "  calculated highest order / requested: ", nm, " / ", order
        end if

        j_n = sj(fnu:order)
        rv = a * dj(fnu:order)

    end function lib_math_bessel_spherical_first_kind_derivative_coeff_real

    ! Calculates the derivative of the spherical Bessel function of the first kind
    !
    ! symbol: j'_n
    !
    ! Argument
    ! ----
    !   z: complex
    !   fnu: integer
    !       ORDER OF INITIAL FUNCTION, FNU.GE.0
    !   n: integer
    !       NUMBER OF MEMBERS OF THE SEQUENCE, N.GE.1
    !
    ! Returns
    ! ----
    !   j_n: array<double precision>
    !       spherical Bessel function of the first kind [fnu..fnu+n-1]
    !   rv: array<double precision>
    !
    !
    ! LaTeX: $$ \frac{\partial}{\partial x}\left(j_{n}(x)\right)=\frac{1}{2}\left(j_{n-1}(x)-\frac{j_{n}(x)+x j_{n+1}(x)}{x}\right) $$
    !
    ! Reference: WolframAlpha: https://www.wolframalpha.com/input/?i=derivative%5BSphericalBesselJ%5Bn,x%5D,+x%5D
    !               >>> derivative[SphericalBesselJ[n,x], x]
    function lib_math_bessel_spherical_first_kind_derivative_cmplx(z, fnu, n, j_n) result (rv)
        implicit none
        ! dummy
        complex(kind=8), intent(in) :: z
        integer(kind=4), intent(in) :: fnu
        integer(kind=4), intent(in) :: n
        complex(kind=8), dimension(n), intent(inout) :: j_n
        complex(kind=8), dimension(n) :: rv

        ! auxiliary
        integer(kind=4) order
        integer(kind=4) nm
        complex(kind=8) csj(0:fnu+n-1)
        complex(kind=8) cdj(0:fnu+n-1)
        complex(kind=8) csy(0:fnu+n-1)
        complex(kind=8) cdy(0:fnu+n-1)

        order = fnu+n-1

        call csphjy ( order, z, nm, csj, cdj, csy, cdy )

        if (order .ne. nm) then
            print *, "lib_math_bessel_spherical_first_kind_derivative_cmplx: ERROR"
            print *, "  calculated highest order / requested: ", nm, " / ", order
        end if

        j_n = csj(fnu:order)
        rv = cdj(fnu:order)

    end function lib_math_bessel_spherical_first_kind_derivative_cmplx

    ! Calculates the derivative of the spherical Bessel function of the first kind
    !
    ! symbol: j'_n
    !
    ! Argument
    ! ----
    !   a: double precision
    !       coefficient
    !   z: complex
    !   fnu: integer
    !       ORDER OF INITIAL FUNCTION, FNU.GE.0
    !   n: integer
    !       NUMBER OF MEMBERS OF THE SEQUENCE, N.GE.1
    !
    ! Returns
    ! ----
    !   j_n: array<double precision>
    !       spherical Bessel function of the first kind [fnu..fnu+n-1]
    !   rv: array<double precision>
    !
    !
    ! LaTeX: $$ \frac{\partial}{\partial x}\left(j_{n}(a x)\right)= j'_{n}(a x) (a x)'$$
    !
    ! Reference: WolframAlpha: https://www.wolframalpha.com/input/?i=derivative%5BSphericalBesselJ%5Bn,x%5D,+x%5D
    !               >>> derivative[SphericalBesselJ[n,x], x]
    function lib_math_bessel_spherical_first_kind_derivative_coeff_cmplx(a, z, fnu, n, j_n) result (rv)
        implicit none
        ! dummy
        real(kind=8), intent(in) :: a
        complex(kind=8), intent(in) :: z
        integer(kind=4), intent(in) :: fnu
        integer(kind=4), intent(in) :: n
        complex(kind=8), dimension(n), intent(inout) :: j_n
        complex(kind=8), dimension(n) :: rv

        ! auxiliary
        integer(kind=4) order
        integer(kind=4) nm
        complex(kind=8) csj(0:fnu+n-1)
        complex(kind=8) cdj(0:fnu+n-1)
        complex(kind=8) csy(0:fnu+n-1)
        complex(kind=8) cdy(0:fnu+n-1)

        order = fnu+n-1

        call csphjy ( order, a * z, nm, csj, cdj, csy, cdy )

        if (order .ne. nm) then
            print *, "lib_math_bessel_spherical_first_kind_derivative_coeff_cmplx: ERROR"
            print *, "  calculated highest order / requested: ", nm, " / ", order
        end if

        j_n = csj(fnu:order)
        rv = a * cdj(fnu:order)

    end function lib_math_bessel_spherical_first_kind_derivative_coeff_cmplx

    ! Calculates the derivative of the spherical Bessel function of the second kind
    !
    ! symbol: y'_n
    !
    ! Argument
    ! ----
    !   x: double precision
    !   fnu: integer
    !       ORDER OF INITIAL FUNCTION, FNU.GE.0
    !   n: integer
    !       NUMBER OF MEMBERS OF THE SEQUENCE, N.GE.1
    !
    ! Returns
    ! ----
    !   y_n: array<double precision>
    !       spherical Bessel function of the second kind [fnu..fnu+n-1]
    !   rv: array<double precision>
    !
    ! LaTeX: $$ y_{v}(z) \equiv \sqrt{\frac{\pi}{2 z}} Y_{v+1 / 2}(z) $$
    !
    ! Reference: http://mathworld.wolfram.com/SphericalBesselFunctionoftheSecondKind.html
    !
    function lib_math_bessel_spherical_second_kind_derivative_real(x, fnu, n, y_n) result (rv)
        implicit none
        ! dummy
        double precision, intent(in) :: x
        integer(kind=4), intent(in) :: fnu
        integer(kind=4), intent(in) :: n
        double precision, dimension(n), intent(inout) :: y_n
        double precision, dimension(n) :: rv

        ! auxiliary
        integer(kind=4) order
        integer(kind=4) nm
        real(kind=8) dy(0:fnu+n-1)
        real(kind=8) sy(0:fnu+n-1)

        order = fnu+n-1

        call SPHY( order, x, nm, sy, dy)

        if (order .ne. nm) then
            print *, "lib_math_bessel_spherical_second_kind_derivative_real: ERROR"
            print *, "  calculated highest order / requested: ", nm, " / ", order
        end if

        y_n = sy(fnu:order)
        rv = dy(fnu:order)

    end function lib_math_bessel_spherical_second_kind_derivative_real

    ! Calculates the derivative of the spherical Bessel function of the second kind
    !
    ! symbol: y'_n
    !
    ! Argument
    ! ----
    !   x: double precision
    !   fnu: integer
    !       ORDER OF INITIAL FUNCTION, FNU.GE.0
    !   n: integer
    !       NUMBER OF MEMBERS OF THE SEQUENCE, N.GE.1
    !
    ! Returns
    ! ----
    !   y_n: array<double precision>
    !       spherical Bessel function of the second kind [fnu..fnu+n-1]
    !   rv: array<double precision>
    !
    ! LaTeX: $$ y_{v}(z) \equiv \sqrt{\frac{\pi}{2 z}} Y_{v+1 / 2}(z) $$
    !
    ! Reference: http://mathworld.wolfram.com/SphericalBesselFunctionoftheSecondKind.html
    !
    function lib_math_bessel_spherical_second_kind_derivative_cmplx(z, fnu, n, y_n) result (rv)
        implicit none
        ! dummy
        complex(kind=8), intent(in) :: z
        integer(kind=4), intent(in) :: fnu
        integer(kind=4), intent(in) :: n
        complex(kind=8), dimension(n), intent(inout) :: y_n
        complex(kind=8), dimension(n) :: rv

        ! auxiliary
        integer(kind=4) :: order
        integer(kind=4) :: nm
        complex(kind=8) csj(0:fnu+n-1)
        complex(kind=8) cdj(0:fnu+n-1)
        complex(kind=8) csy(0:fnu+n-1)
        complex(kind=8) cdy(0:fnu+n-1)

        order = fnu+n-1

        call csphjy ( order, z, nm, csj, cdj, csy, cdy )

        if (order .ne. nm) then
            print *, "lib_math_bessel_spherical_second_kind_derivative_cmplx: ERROR"
            print *, "  calculated highest order / requested: ", nm, " / ", order
        end if

        y_n = csy(fnu:order)
        rv = cdy(fnu:order)

    end function lib_math_bessel_spherical_second_kind_derivative_cmplx

    ! Calculates the derivative of the spherical Bessel function of the third kind 1
    !
    ! symbol: h'^(1)_n
    !
    ! Argument
    ! ----
    !   x: double precision
    !   fnu: integer
    !       ORDER OF INITIAL FUNCTION, FNU.GE.1
    !   n: integer
    !       NUMBER OF MEMBERS OF THE SEQUENCE, N.GE.1
    !
    ! Returns
    ! ----
    !   h_1_n: array<double precision>
    !       spherical Bessel function of the third kind 1 [fnu-1..fnu+n-1]
    !   rv: array<double precision>
    !
    ! LaTeX: $$ \frac{d}{d z} h_{n}^{(1)}(z)=\frac{1}{2}\left[h_{n-1}^{(1)}(z)-\frac{h_{n}^{(1)}(z)+z h_{n+1}^{(1)}(z)}{z}\right] $$
    !
    ! Reference: http://mathworld.wolfram.com/SphericalHankelFunctionoftheFirstKind.html
    !
    function lib_math_bessel_spherical_third_kind_1_derivative_real(x, fnu, n, h_1_n) result (rv)
        implicit none
        ! dummy
        double precision, intent(in) :: x
        integer(kind=4), intent(in) :: fnu
        integer(kind=4), intent(in) :: n
        complex(kind=8), dimension(n), intent(inout) :: h_1_n
        complex(kind=8), dimension(n) :: rv

        ! auxiliary
        integer :: i
        complex(kind=8), dimension(n+2) :: buffer_h_1_n

        buffer_h_1_n = lib_math_bessel_spherical_third_kind_1_real(x, fnu-1, n+2)

        h_1_n = buffer_h_1_n(2:n+1)

        do i=1, n
             rv(i) = 0.5D0 * (buffer_h_1_n(i) - (buffer_h_1_n(i+1) + x * buffer_h_1_n(i+2))/x)
        end do

    end function lib_math_bessel_spherical_third_kind_1_derivative_real

    ! Calculates the derivative of the spherical Bessel function of the third kind 1
    !
    ! symbol: h'^(1)_n
    !
    ! Argument
    ! ----
    !   x: double precision
    !   fnu: integer
    !       ORDER OF INITIAL FUNCTION, FNU.GE.1
    !   n: integer
    !       NUMBER OF MEMBERS OF THE SEQUENCE, N.GE.1
    !
    ! Returns
    ! ----
    !   h_1_n: array<double precision>
    !       spherical Bessel function of the third kind 1 [fnu-1..fnu+n-1]
    !   rv: array<double precision>
    !
    ! LaTeX: $$ \frac{d}{d z} h_{n}^{(1)}(z)=\frac{1}{2}\left[h_{n-1}^{(1)}(z)-\frac{h_{n}^{(1)}(z)+z h_{n+1}^{(1)}(z)}{z}\right] $$
    !
    ! Reference: http://mathworld.wolfram.com/SphericalHankelFunctionoftheFirstKind.html
    !
    function lib_math_bessel_spherical_third_kind_1_derivative_cmplx(z, fnu, n, h_1_n) result (rv)
        implicit none
        ! dummy
        complex(kind=8), intent(in) :: z
        integer(kind=4), intent(in) :: fnu
        integer(kind=4), intent(in) :: n
        complex(kind=8), dimension(n), intent(inout) :: h_1_n
        complex(kind=8), dimension(n) :: rv

        ! auxiliary
        complex(kind=8), dimension(n+2) :: buffer_h_1_n
        integer :: i

        buffer_h_1_n = lib_math_bessel_spherical_third_kind_1_cmplx(z, fnu-1, n+2)
        h_1_n = buffer_h_1_n(2:n+1)

        do i=1, n
             rv(i) = 0.5D0 * (buffer_h_1_n(i) - (buffer_h_1_n(i+1) + z * buffer_h_1_n(i+2))/z)
        end do

    end function lib_math_bessel_spherical_third_kind_1_derivative_cmplx

    ! Calculates the derivative of the spherical Bessel function of the third kind 1
    !
    ! symbol: h'^(2)_n
    !
    ! Argument
    ! ----
    !   x: double precision
    !   fnu: integer
    !       ORDER OF INITIAL FUNCTION, FNU.GE.1
    !   n: integer
    !       NUMBER OF MEMBERS OF THE SEQUENCE, N.GE.1
    !
    ! Returns
    ! ----
    !   h_2_n: array<double precision>
    !       spherical Bessel function of the third kind 2 [fnu-1..fnu+n-1]
    !   rv: array<double precision>
    !
    ! LaTeX: $$ \frac{d}{d z} h_{n}^{(2)}(z)=\frac{1}{2}\left[h_{n-1}^{(2)}(z)-\frac{h_{n}^{(2)}(z)+z h_{n+1}^{(2)}(z)}{z}\right] $$
    !
    ! Reference: http://mathworld.wolfram.com/SphericalHankelFunctionoftheSecondKind.html
    !
    function lib_math_bessel_spherical_third_kind_2_derivative_real(x, fnu, n, h_2_n) result (rv)
        implicit none
        ! dummy
        double precision, intent(in) :: x
        integer(kind=4), intent(in) :: fnu
        integer(kind=4), intent(in) :: n
        complex(kind=8), dimension(n), intent(inout) :: h_2_n
        complex(kind=8), dimension(n) :: rv

        ! auxiliary
        integer :: i
        complex(kind=8), dimension(n+2) :: buffer_h_2_n

        buffer_h_2_n = lib_math_bessel_spherical_third_kind_2_real(x, fnu-1, n+2)

        h_2_n = buffer_h_2_n(2:n+1)

        do i=1, n
             rv(i) = 0.5D0 * (buffer_h_2_n(i) - (buffer_h_2_n(i+1) + x * buffer_h_2_n(i+2))/x)
        end do

    end function lib_math_bessel_spherical_third_kind_2_derivative_real

    ! Calculates the derivative of the spherical Bessel function of the third kind 1
    !
    ! symbol: h'^(2)_n
    !
    ! Argument
    ! ----
    !   x: double precision
    !   fnu: integer
    !       ORDER OF INITIAL FUNCTION, FNU.GE.1
    !   n: integer
    !       NUMBER OF MEMBERS OF THE SEQUENCE, N.GE.1
    !
    ! Returns
    ! ----
    !   h_2_n: array<double precision>
    !       spherical Bessel function of the third kind 2 [fnu-1..fnu+n-1]
    !   rv: array<double precision>
    !
    ! LaTeX: $$ \frac{d}{d z} h_{n}^{(2)}(z)=\frac{1}{2}\left[h_{n-1}^{(2)}(z)-\frac{h_{n}^{(2)}(z)+z h_{n+1}^{(2)}(z)}{z}\right] $$
    !
    ! Reference: http://mathworld.wolfram.com/SphericalHankelFunctionoftheSecondKind.html
    !
    function lib_math_bessel_spherical_third_kind_2_derivative_cmplx(z, fnu, n, h_2_n) result (rv)
        implicit none
        ! dummy
        complex(kind=8), intent(in) :: z
        integer(kind=4), intent(in) :: fnu
        integer(kind=4), intent(in) :: n
        complex(kind=8), dimension(n), intent(inout) :: h_2_n
        complex(kind=8), dimension(n) :: rv

        ! auxiliary
        complex(kind=8), dimension(n+2) :: buffer_h_2_n
        integer :: i

        buffer_h_2_n = lib_math_bessel_spherical_third_kind_2_cmplx(z, fnu-1, n+2)
        h_2_n = buffer_h_2_n(2:n+1)

        do i=1, n
             rv(i) = 0.5D0 * (buffer_h_2_n(i) - (buffer_h_2_n(i+1) + z * buffer_h_2_n(i+2))/z)
        end do

    end function lib_math_bessel_spherical_third_kind_2_derivative_cmplx

    ! Calculates the derivative of the Riccati-Bessel function
    !
    ! symbol: S'
    !
    ! Argument
    ! ----
    !   x: double precision
    !   fnu: integer
    !       ORDER OF INITIAL FUNCTION, FNU.GE.1
    !   n: integer
    !       NUMBER OF MEMBERS OF THE SEQUENCE, N.GE.1
    !
    ! Returns
    ! ----
    !   S_n: array<double precision>
    !       Riccati-Bessel function [fnu..fnu+n-1]
    !   rv: array<double precision>
    !
    ! LaTeX: $$  S^\prime_n = -x j_{n+1}(x)+(n+1) j_{n}(x) $$
    !
    ! Reference:
    !
    function lib_math_bessel_riccati_s_derivative_real(x, fnu, n, r_n) result (rv)
        implicit none
        ! dummy
        double precision, intent(in) :: x
        integer(kind=4), intent(in) :: fnu
        integer(kind=4), intent(in) :: n
        double precision, intent(inout), dimension(n) :: r_n

        double precision, dimension(n) :: rv

        ! auxiliary
        integer(kind=4) order
        integer(kind=4) nm
        real(kind=8) rf(0:fnu+n-1)
        real(kind=8) df(0:fnu+n-1)

        order = fnu+n-1

        call rctj( order, x, nm, rf, df)

        if (order .ne. nm) then
            print *, "lib_math_bessel_riccati_s_derivative_real: ERROR"
            print *, "  calculated highest order / requested: ", nm, " / ", order
        end if

        rv = df(fnu:order)
        r_n = rf(fnu:order)

    end function lib_math_bessel_riccati_s_derivative_real

    ! Calculates the derivative of the Riccati-Bessel function
    !
    ! symbol: S'
    !
    ! Argument
    ! ----
    !   x: double precision
    !   fnu: integer
    !       ORDER OF INITIAL FUNCTION, FNU.GE.1
    !   n: integer
    !       NUMBER OF MEMBERS OF THE SEQUENCE, N.GE.1
    !
    ! Returns
    ! ----
    !   S_n: array<double precision>
    !       Riccati-Bessel function [fnu..fnu+n-1]
    !   rv: array<double precision>
    !
    ! LaTeX: $$  S^\prime_n = $$ \frac{1}{2}\left(x j_{n-1}(x)+j_{n}(x)-x j_{n+1}(x)\right) $$
    !
    ! Reference: Wolfram Alpha: https://www.wolframalpha.com/input/?i=derivative%5Bx+SphericalBesselJ%5Bn,x%5D,+x%5D
    !               >>> derivative[x SphericalBesselJ[n,x], x]
    !
    function lib_math_bessel_riccati_s_derivative_cmplx(z, fnu, n, s_n) result (rv)
        implicit none
        ! dummy
        complex(kind=8), intent(in) :: z
        integer(kind=4), intent(in) :: fnu
        integer(kind=4), intent(in) :: n
        complex(kind=8), intent(inout), dimension(n) :: s_n

        complex(kind=8), dimension(n) :: rv

        ! auxiliary
        integer(kind=4) i
        integer(kind=4) order
        complex(kind=8), dimension(n+2) :: jf
        complex(kind=8), dimension(n+2) :: buffer_r_n

        order = fnu+n-1

        jf = lib_math_bessel_spherical_first_kind_cmplx(z, fnu-1, n+2)

        buffer_r_n = z * jf

        s_n = buffer_r_n(fnu+1:order+1)

        do i=2, n+1
            rv(i-1) = (buffer_r_n(i-1) + jf(i) - buffer_r_n(i+1)) / 2
        end do

    end function lib_math_bessel_riccati_s_derivative_cmplx

    ! Calculates the derivative of the Riccati-Bessel function
    !
    ! symbol: C'
    !
    ! Argument
    ! ----
    !   x: double precision
    !   fnu: integer
    !       ORDER OF INITIAL FUNCTION, FNU.GE.1
    !   n: integer
    !       NUMBER OF MEMBERS OF THE SEQUENCE, N.GE.1
    !
    ! Returns
    ! ----
    !   S_n: array<double precision>
    !       Riccati-Bessel function [fnu..fnu+n-1]
    !   rv: array<double precision>
    !
    ! LaTeX: $$ C^\prime_n = x y_{n+1}(x)-(n+1) y_{n}(x) $$
    !
    ! Reference:
    !
    function lib_math_bessel_riccati_c_derivative_real(x, fnu, n, r_n) result (rv)
        implicit none
        ! dummy
        double precision, intent(in) :: x
        integer(kind=4), intent(in) :: fnu
        integer(kind=4), intent(in) :: n
        double precision, intent(inout), dimension(n) :: r_n

        double precision, dimension(n) :: rv

        ! auxiliary
        integer(kind=4) order
        integer(kind=4) nm
        real(kind=8) rf(0:fnu+n-1)
        real(kind=8) df(0:fnu+n-1)

        order = fnu+n-1

        call rcty( order, x, nm, rf, df)

        if (order .ne. nm) then
            print *, "lib_math_bessel_riccati_c_derivative_real: ERROR"
            print *, "  calculated highest order / requested: ", nm, " / ", order
        end if

        ! tests with wolfram alpha and sageMath results yield a prefactor of "-1"
        ! wolfram alphs
        !   >>> derivative[ riccatiBesselC[n,x], x]
        !   >>> N[1/2*(-x sphericalBesselY[n-1,x] - sphericalBesselY[n,x] + x sphericalBesselY[n+1,x] ), 20] where n = 3, x=4
        ! sageMath
        !   >>> var('x')
        !   >>> var('n')
        !   >>> C(n,x) = -x*spherical_bessel_Y(n,x)
        !   >>> C_(n,x) = derivative(C(n,x), x)
        !   >>> n=5
        !   >>> x=4.0
        !   >>> for i in range(1,6):
        !   >>>     value = numerical_approx(C_(i,x))
        !   >>>     print("n = {}: {}".format(i, value))
        rv = -df(fnu:order)
        ! tests with wolfram alpha and sageMath results yield a prefactor of "-1"
        ! wolfram alphs
        !   >>> riccatiBesselC[2, 4]
        ! sageMath
        !   >>> C(n,x) = -x*spherical_bessel_Y(n,x)
        !   >>> numerical_approx(C(2,4))
        r_n = -rf(fnu:order)

    end function lib_math_bessel_riccati_c_derivative_real

    ! Calculates the derivative of the Riccati-Bessel function
    !
    ! symbol: C'
    !
    ! Argument
    ! ----
    !   x: double precision
    !   fnu: integer
    !       ORDER OF INITIAL FUNCTION, FNU.GE.1
    !   n: integer
    !       NUMBER OF MEMBERS OF THE SEQUENCE, N.GE.1
    !
    ! Returns
    ! ----
    !   c_n: array<double precision>
    !       Riccati-Bessel function [fnu..fnu+n-1]
    !   rv: array<double precision>
    !
    ! LaTeX: $$ C^\prime_n = \frac{1}{2}\left(-x y_{n-1}(x)-y_{n}(x)+x y_{n+1}(x)\right) $$
    !
    ! Reference: WolframAlpha: https://www.wolframalpha.com/input/?i=derivative%5B-x+SphericalBesselY%5Bn,x%5D,+x%5D
    !               >>> derivative[-x SphericalBesselY[n,x], x]
    !
    function lib_math_bessel_riccati_c_derivative_cmplx(z, fnu, n, c_n) result (rv)
        implicit none
        ! dummy
        complex(kind=8), intent(in) :: z
        integer(kind=4), intent(in) :: fnu
        integer(kind=4), intent(in) :: n
        complex(kind=8), intent(inout), dimension(n) :: c_n

        complex(kind=8), dimension(n) :: rv

        ! auxiliary
        integer(kind=4) i
        integer(kind=4) order
        complex(kind=8), dimension(n+2) :: yf
        complex(kind=8), dimension(n+2) :: buffer_c_n

        order = fnu+n-1

        yf = lib_math_bessel_spherical_second_kind_cmplx(z, fnu-1, n+2)

        buffer_c_n = z * yf

        c_n = -buffer_c_n(fnu+1:order+1)

        do i=2, n+1
            rv(i-1) = (-buffer_c_n(i-1) - yf(i) + buffer_c_n(i+1)) / 2
        end do

    end function lib_math_bessel_riccati_c_derivative_cmplx

    ! Calculates the derivative of the Riccati-Bessel function
    !
    ! symbol: Xi'
    !
    ! Argument
    ! ----
    !   x: double precision
    !   fnu: integer
    !       ORDER OF INITIAL FUNCTION, FNU.GE.1
    !   n: integer
    !       NUMBER OF MEMBERS OF THE SEQUENCE, N.GE.1
    !
    ! Returns
    ! ----
    !   xi_n: array<double precision>
    !       Riccati-Bessel function [fnu..fnu+n-1]
    !   rv: array<double precision>
    !
    ! LaTeX: $$ \xi^\prime_n (x) = \frac{\partial}{\partial x}\left(x h_{n}^{(1)}(x)\right)=\frac{1}{2}\left(x h_{n-1}^{(1)}(x)+h_{n}^{(1)}(x)-x h_{n+1}^{(1)}(x)\right) $$
    !
    ! Reference: WolframAlpha
    !           >>> derivative[x*SphericalHankelH1(n,x), x]
    !
    function lib_math_bessel_riccati_xi_derivative_real(x, fnu, n, xi_n) result (rv)
        implicit none
        ! dummy
        double precision, intent(in) :: x
        integer(kind=4), intent(in) :: fnu
        integer(kind=4), intent(in) :: n
        complex(kind=8), intent(inout), dimension(n) :: xi_n

        complex(kind=8), dimension(n) :: rv

        ! auxiliary
        integer(kind=4) i
        integer(kind=4) order
        complex(kind=8), dimension(n+2) :: h_n
        complex(kind=8), dimension(n+2) :: buffer_xi_n

        order = fnu+n-1

        h_n = lib_math_bessel_spherical_third_kind_1_real(x, fnu-1, n+2)

        buffer_xi_n = x * h_n

        xi_n = buffer_xi_n(fnu+1:order+1)

        do i=2, n+1
            rv(i-1) = (buffer_xi_n(i-1) + h_n(i) - buffer_xi_n(i+1)) / 2.0_8
        end do

    end function lib_math_bessel_riccati_xi_derivative_real

    ! Calculates the derivative of the Riccati-Bessel function
    !
    ! symbol: Xi'
    !
    ! Argument
    ! ----
    !   x: double precision
    !   fnu: integer
    !       ORDER OF INITIAL FUNCTION, FNU.GE.1
    !   n: integer
    !       NUMBER OF MEMBERS OF THE SEQUENCE, N.GE.1
    !
    ! Returns
    ! ----
    !   xi_n: array<double precision>
    !       Riccati-Bessel function [fnu..fnu+n-1]
    !   rv: array<double precision>
    !
    ! LaTeX: $$ \xi^\prime_n (x) = \frac{\partial}{\partial x}\left(x h_{n}^{(1)}(x)\right)=\frac{1}{2}\left(x h_{n-1}^{(1)}(x)+h_{n}^{(1)}(x)-x h_{n+1}^{(1)}(x)\right) $$
    !
    ! Reference: WolframAlpha
    !           >>> derivative[x*SphericalHankelH1(n,x), x]
    !
    function lib_math_bessel_riccati_xi_derivative_cmplx(x, fnu, n, xi_n) result (rv)
        implicit none
        ! dummy
        complex(kind=8), intent(in) :: x
        integer(kind=4), intent(in) :: fnu
        integer(kind=4), intent(in) :: n
        complex(kind=8), intent(inout), dimension(n) :: xi_n

        complex(kind=8), dimension(n) :: rv

        ! auxiliary
        integer(kind=4) i
        integer(kind=4) order
        complex(kind=8), dimension(n+2) :: h_n
        complex(kind=8), dimension(n+2) :: buffer_xi_n

        order = fnu+n-1

        h_n = lib_math_bessel_spherical_third_kind_1_cmplx(x, fnu-1, n+2)

        buffer_xi_n = x * h_n

        xi_n = buffer_xi_n(fnu+1:order+1)

        do i=2, n+1
            rv(i-1) = (buffer_xi_n(i-1) + h_n(i) - buffer_xi_n(i+1)) / 2.0_8
        end do

    end function lib_math_bessel_riccati_xi_derivative_cmplx

    ! Calculates the derivative of the Riccati-Bessel function
    !
    ! symbol: Zeta'
    !
    ! Argument
    ! ----
    !   x: double precision
    !   fnu: integer
    !       ORDER OF INITIAL FUNCTION, FNU.GE.1
    !   n: integer
    !       NUMBER OF MEMBERS OF THE SEQUENCE, N.GE.1
    !
    ! Returns
    ! ----
    !   xi_n: array<double precision>
    !       Riccati-Bessel function [fnu..fnu+n-1]
    !   rv: array<double precision>
    !
    ! LaTeX: $$ \zeta^\prime_n (x) = \frac{\partial}{\partial x}\left(x h_{n}^{(2)}(x)\right)=\frac{1}{2}\left(x h_{n-1}^{(2)}(x)+h_{n}^{(2)}(x)-x h_{n+1}^{(2)}(x)\right) $$
    !
    ! Reference: WolframAlpha
    !           >>> derivative[x*SphericalHankelH1(n,x), x]
    !
    function lib_math_bessel_riccati_zeta_derivative_real(x, fnu, n, zeta_n) result (rv)
        implicit none
        ! dummy
        double precision, intent(in) :: x
        integer(kind=4), intent(in) :: fnu
        integer(kind=4), intent(in) :: n
        complex(kind=8), intent(inout), dimension(n) :: zeta_n

        complex(kind=8), dimension(n) :: rv

        ! auxiliary
        integer(kind=4) i
        integer(kind=4) order
        complex(kind=8), dimension(n+2) :: h_n
        complex(kind=8), dimension(n+2) :: buffer_zeta_n

        order = fnu+n-1

        h_n = lib_math_bessel_spherical_third_kind_2_real(x, fnu-1, n+2)

        buffer_zeta_n = x * h_n

        zeta_n = buffer_zeta_n(fnu+1:order+1)

        do i=2, n+1
            rv(i-1) = (buffer_zeta_n(i-1) + h_n(i) - buffer_zeta_n(i+1)) / 2.0_8
        end do

    end function lib_math_bessel_riccati_zeta_derivative_real

    ! Calculates the derivative of the Riccati-Bessel function
    !
    ! symbol: Zeta'
    !
    ! Argument
    ! ----
    !   z: complex
    !   fnu: integer
    !       ORDER OF INITIAL FUNCTION, FNU.GE.1
    !   n: integer
    !       NUMBER OF MEMBERS OF THE SEQUENCE, N.GE.1
    !
    ! Returns
    ! ----
    !   xi_n: array<double precision>
    !       Riccati-Bessel function [fnu..fnu+n-1]
    !   rv: array<double precision>
    !
    ! LaTeX: $$ \zeta^\prime_n (x) = \frac{\partial}{\partial x}\left(x h_{n}^{(2)}(x)\right)=\frac{1}{2}\left(x h_{n-1}^{(2)}(x)+h_{n}^{(2)}(x)-x h_{n+1}^{(2)}(x)\right) $$
    !
    ! Reference: WolframAlpha
    !           >>> derivative[x*SphericalHankelH1(n,x), x]
    !
    function lib_math_bessel_riccati_zeta_derivative_cmplx(z, fnu, n, zeta_n) result (rv)
        implicit none
        ! dummy
        complex(kind=8), intent(in) :: z
        integer(kind=4), intent(in) :: fnu
        integer(kind=4), intent(in) :: n
        complex(kind=8), intent(inout), dimension(n) :: zeta_n

        complex(kind=8), dimension(n) :: rv

        ! auxiliary
        integer(kind=4) i
        integer(kind=4) order
        complex(kind=8), dimension(n+2) :: h_n
        complex(kind=8), dimension(n+2) :: buffer_zeta_n

        order = fnu+n-1

        h_n = lib_math_bessel_spherical_third_kind_2_cmplx(z, fnu-1, n+2)

        buffer_zeta_n = z * h_n

        zeta_n = buffer_zeta_n(fnu+1:order+1)

        do i=2, n+1
            rv(i-1) = (buffer_zeta_n(i-1) + h_n(i) - buffer_zeta_n(i+1)) / 2.0_8
        end do

    end function lib_math_bessel_riccati_zeta_derivative_cmplx

    function lib_math_bessel_test_functions() result (rv)
        implicit none
        ! dummy
        integer :: rv

        rv = 0

        ! test: Bessel functions
        if (.not. test_lib_math_bessel_spherical_first_kind_real()) then
            rv = rv + 1
        end if
        if (.not. test_lib_math_bessel_spherical_second_kind_real()) then
            rv = rv + 1
        end if
        if (.not. test_lib_math_bessel_spherical_third_kind_1_real()) then
            rv = rv + 1
        end if
        if (.not. test_lib_math_bessel_spherical_third_kind_2_real()) then
            rv = rv + 1
        end if
        if (.not. test_lib_math_bessel_riccati_s_real()) then
            rv = rv + 1
        end if
        if (.not. test_lib_math_bessel_riccati_c_real()) then
            rv = rv + 1
        end if
        if (.not. test_lib_math_bessel_riccati_xi_real()) then
            rv = rv + 1
        end if
        if (.not. test_lib_math_bessel_riccati_zeta_real()) then
            rv = rv + 1
        end if

        ! test: derivatives of the Bessel functions
        if (.not. test_lib_math_bessel_spherical_first_kind_derivative_real()) then
            rv = rv + 1
        end if
        if (.not. test_lib_math_bessel_spherical_first_kind_derivative_coeff_real()) then
            rv = rv + 1
        end if
        if (.not. test_lib_math_bessel_spherical_second_kind_derivative_real()) then
            rv = rv + 1
        end if
        if (.not. test_lib_math_bessel_spherical_third_kind_1_derivative_real()) then
            rv = rv + 1
        end if
        if (.not. test_lib_math_bessel_spherical_third_kind_2_derivative_real()) then
            rv = rv + 1
        end if
        if (.not. test_lib_math_bessel_riccati_s_derivative_real()) then
            rv = rv + 1
        end if
        if (.not. test_lib_math_bessel_riccati_c_derivative_real()) then
            rv = rv + 1
        end if
        if (.not. test_lib_math_bessel_riccati_xi_derivative_real()) then
            rv = rv + 1
        end if
        if (.not. test_lib_math_bessel_riccati_zeta_derivative_real()) then
            rv = rv + 1
        end if


        ! test: Bessel functions with a complex argument
        if (.not. test_lib_math_bessel_spherical_first_kind_cmplx()) then
            rv = rv + 1
        end if
        if (.not. test_lib_math_bessel_spherical_second_kind_cmplx()) then
            rv = rv + 1
        end if
        if (.not. test_lib_math_bessel_spherical_third_kind_1_cmplx()) then
            rv = rv + 1
        end if
        if (.not. test_lib_math_bessel_spherical_third_kind_2_cmplx()) then
            rv = rv + 1
        end if
        if (.not. test_lib_math_bessel_riccati_s_cmplx()) then
            rv = rv + 1
        end if
        if (.not. test_lib_math_bessel_riccati_c_cmplx()) then
            rv = rv + 1
        end if
        if (.not. test_lib_math_bessel_riccati_xi_cmplx()) then
            rv = rv + 1
        end if
        if (.not. test_lib_math_bessel_riccati_zeta_cmplx()) then
            rv = rv + 1
        end if

        ! test: derivatives of Bessel functions with complex argument
        if (.not. test_lib_math_bessel_spherical_first_kind_derivative_cmplx()) then
            rv = rv + 1
        end if
        if (.not. test_lib_math_bessel_spherical_first_kind_derivative_coeff_c()) then
            rv = rv + 1
        end if
        if (.not. test_lib_math_bessel_spherical_second_kind_derivative_cmplx()) then
            rv = rv + 1
        end if
        if (.not. test_lib_math_bessel_spherical_third_kind_1_derivative_cmplx()) then
            rv = rv + 1
        end if
        if (.not. test_lib_math_bessel_spherical_third_kind_2_derivative_cmplx()) then
            rv = rv + 1
        end if
        if (.not. test_lib_math_bessel_riccati_s_derivative_cmplx()) then
            rv = rv + 1
        end if
        if (.not. test_lib_math_bessel_riccati_c_derivative_cmplx()) then
            rv = rv + 1
        end if
        if (.not. test_lib_math_bessel_riccati_xi_derivative_cmplx()) then
            rv = rv + 1
        end if
        if (.not. test_lib_math_bessel_riccati_zeta_derivative_cmplx()) then
            rv = rv + 1
        end if


        print *, "-------------lib_math_bessel_test_functions----------------"
        if (rv == 0) then
            print *, "lib_math_bessel_test_functions tests: OK"
        else
            print *, rv,"lib_math_bessel_test_functions test(s) FAILED"
        end if
        print *, "-----------------------------------------------------------"

        contains

        function test_lib_math_bessel_spherical_first_kind_real() result (rv)
            implicit none
            ! dummy
            logical :: rv

            ! parameter
            integer, parameter :: n = 5
            double precision, parameter :: ground_truth_e = 10.0_8**(-13.0_8)

            ! auxiliary
            double precision :: x = 5
            double precision, dimension(n) :: j
            integer :: fnu = 0

            double precision, dimension(n) :: ground_truth_j

            integer :: i
            double precision :: buffer

            ! Values were generated with sageMath
            !
            ! source code:
            !  >>> x=5.0 +I*0.0
            !  >>> for n in range(0,5):
            !  >>>     value = numerical_approx(spherical_bessel_J(n, x))
            !  >>>     print("n = {}: {}".format(n, value))
            data ground_truth_j / -0.191784854932628_8, -0.0950894080791708_8, 0.134731210085125_8,&
                                  0.229820618164296_8, 0.187017655344889_8/

            j = lib_math_bessel_spherical_first_kind_real(x, fnu, n)

            rv = .true.
            print *, "test_lib_math_bessel_spherical_first_kind_real:"
            do i=1, n
                buffer = j(i) - ground_truth_j(i)
                if (abs(buffer) .gt. ground_truth_e) then
                    print *, "  ", i , "difference: ", buffer, " : FAILED"
                    rv = .false.
                else
                    print *, "  ", i, ": OK"
                end if
            end do

        end function

        function test_lib_math_bessel_spherical_first_kind_cmplx() result (rv)
            implicit none
            ! dummy
            logical :: rv

            ! parameter
            integer, parameter :: n = 5
            double precision, parameter :: ground_truth_e = 10.0_8**(-13.0_8)

            ! auxiliary
            complex(kind=8) :: x = cmplx(2, 2, kind=8)
            complex(kind=8), dimension(n) :: y
            integer :: fnu = 0

            complex(kind=8), dimension(n) :: ground_truth_y

            integer :: i
            double precision :: buffer

            ! Values were generated with sageMath
            !
            ! source code:
            !  >>> x=5.0 +I*0.0
            !  >>> for n in range(0,5):
            !  >>>     value = numerical_approx(spherical_bessel_J(n, x))
            !  >>>     print("n = {}: {}".format(n, value))
            ground_truth_y(1) = cmplx(0.477912093948349_8, -1.23256533661016_8, kind=8)
            ground_truth_y(2) = cmplx(1.02721685724129_8, +0.00544789260924664_8, kind=8)
            ground_truth_y(3) = cmplx(0.296586468439555_8, +0.466238613136122_8, kind=8)
            ground_truth_y(4) = cmplx(-0.0736855052716959_8, +0.206617288261462_8, kind=8)
            ground_truth_y(5) = cmplx(-0.0639558482074643_8, +0.0242912755469043_8, kind=8)

            y = lib_math_bessel_spherical_first_kind_cmplx(x, fnu, n)

            rv = .true.
            print *, "test_lib_math_bessel_spherical_first_kind_cmplx:"
            do i=1, n
                buffer = abs(y(i) - ground_truth_y(i))
                if (buffer .gt. ground_truth_e) then
                    print *, "  ", i , "difference: ", buffer, " : FAILED"
                    rv = .false.
                else
                    print *, "  ", i, ": OK"
                end if
            end do

        end function test_lib_math_bessel_spherical_first_kind_cmplx

        function test_lib_math_bessel_spherical_second_kind_real() result (rv)
            implicit none
            ! dummy
            logical :: rv

            ! parameter
            integer, parameter :: n = 5

            ! auxiliary
            double precision :: x = 5
            double precision, dimension(n) :: y
            integer :: fnu = 0

            double precision, dimension(n) :: ground_truth_y
            double precision :: ground_truth_e = 10.0_8**(-13.0_8)

            integer :: i
            double precision :: buffer

            ! Values were generated with sageMath
            !
            ! source code:
            !  >>> x=5.0 +I*0.0
            !  >>> for n in range(0,5):
            !  >>>     value = numerical_approx(spherical_bessel_Y(n, x))
            !  >>>     print("n = {}: {}".format(n, value))
            data ground_truth_y / -0.0567324370926452_8, 0.180438367514099_8, 0.164995457601104_8,&
                                  -0.0154429099129942_8, -0.186615531479296_8/

            y = lib_math_bessel_spherical_second_kind_real(x, fnu, n)

            rv = .true.
            print *, "test_lib_math_bessel_spherical_second_kind_real:"
            do i=1, n
                buffer = y(i) - ground_truth_y(i)
                if (abs(buffer) .gt. ground_truth_e) then
                    print *, "  ", i , "difference: ", buffer, " : FAILED"
                    rv = .false.
                else
                    print *, "  ", i, ": OK"
                end if
            end do

        end function test_lib_math_bessel_spherical_second_kind_real

        function test_lib_math_bessel_spherical_second_kind_cmplx() result (rv)
            implicit none
            ! dummy
            logical :: rv

            ! parameter
            integer, parameter :: n = 5

            ! auxiliary
            complex(kind=8) :: x = cmplx(5, 2, kind = 8)
            complex(kind=8), dimension(n) :: y
            integer :: fnu = 0

            complex(kind=8), dimension(n) :: ground_truth_y
            double precision :: ground_truth_e = 10.0_8**(-13.0_8)

            integer :: i
            double precision :: buffer

            ! Values were generated with sageMath
            !
            ! source code:
            !  >>> x=5.0 +I*2.0
            !  >>> for n in range(0,5):
            !  >>>     value = numerical_approx(spherical_bessel_Y(n, x))
            !  >>>     print("n = {}: {}".format(n, value))
            ground_truth_y(1) = cmplx(-0.423852835557376_8, -0.526035762956881_8, kind = 8)
            ground_truth_y(2) = cmplx(0.441702133380871_8, -0.487648635875712_8, kind = 8)
            ground_truth_y(3) = cmplx(0.551426635055955_8, +0.182417061632022_8, kind = 8)
            ground_truth_y(4) = cmplx(0.0965680905059940_8, +0.454757607952988_8, kind = 8)
            ground_truth_y(5) = cmplx(-0.215340784054175_8, +0.319809249101448_8, kind = 8)

            y = lib_math_bessel_spherical_second_kind_cmplx(x, fnu, n)

            rv = .true.
            print *, "test_lib_math_bessel_spherical_second_kind_cmplx:"
            do i=1, n
                buffer = abs(y(i) - ground_truth_y(i))
                if (buffer .gt. ground_truth_e) then
                    print *, "  ", i , "difference: ", buffer, " : FAILED"
                    rv = .false.
                else
                    print *, "  ", i, ": OK"
                end if
            end do

        end function test_lib_math_bessel_spherical_second_kind_cmplx

        function test_lib_math_bessel_spherical_third_kind_1_real() result (rv)
            implicit none
            ! dummy
            logical :: rv

            ! parameter
            integer, parameter :: n = 41

            ! auxiliary
            double precision :: x = 2
            complex(kind=8), dimension(n) :: y
            integer :: fnu = 0

            complex(kind=8), dimension(n) :: ground_truth_y
            double precision :: ground_truth_e = 10.0_8**(-14.0_8)

            integer :: i
            double precision :: buffer

            ! Values were generated with sageMath
            !
            ! source code:
            !  >>> x=5.0 +I*0.0
            !  >>> for n in range(0,5):
            !  >>>     value = numerical_approx(spherical_hankel1(n, x))
            !  >>>     print("n = {}: {}".format(n, value))
            ground_truth_y(1)  = cmplx(0.454648713412841_8, + 0.208073418273571_8, kind=8)
            ground_truth_y(2)  = cmplx(0.435397774979992_8, - 0.350612004276055_8, kind=8)
            ground_truth_y(3)  = cmplx(0.198447949057147_8, - 0.733991424687654_8, kind=8)
            ground_truth_y(4)  = cmplx(0.0607220976628748_8, - 1.48436655744308_8, kind=8)
            ground_truth_y(5)  = cmplx(0.0140793927629153_8, - 4.46129152636313_8, kind=8)
            ground_truth_y(6)  = cmplx(0.00263516977024412_8, - 18.5914453111910_8, kind=8)
            ground_truth_y(7)  = cmplx(0.000414040973427324_8, - 97.7916576851873_8, kind=8)
            ground_truth_y(8)  = cmplx(0.0000560965570334895_8, - 617.054329642526_8, kind=8)
            ground_truth_y(9)  = cmplx(6.68320432384702d-6, - 4530.11581463376_8, kind=8)
            ground_truth_y(10) = cmplx(7.10679719210186d-7, - 37888.9300947444_8, kind=8)
            ground_truth_y(11) = cmplx(6.82530086497472d-8, - 355414.720085438_8, kind=8)
            ground_truth_y(12) = cmplx(5.97687161216011d-9, - 3.69396563080236d6, kind=8)
            ground_truth_y(13) = cmplx(4.81014890094075d-10, - 4.21251900341417d7, kind=8)
            ground_truth_y(14) = cmplx(3.58145140158186d-11, - 5.22870909795969d8, kind=8)
            ground_truth_y(15) = cmplx(2.48104911947672d-12, - 7.01663209221144d9, kind=8)
            ground_truth_y(16) = cmplx(1.60698216593841d-13, - 1.01218294427270d11, kind=8)
            ground_truth_y(17) = cmplx(9.77323772781462d-15, - 1.56186693153047d12, kind=8)
            ground_truth_y(18) = cmplx(5.60205915100117d-16, - 2.56695860758255d13, kind=8)
            ground_truth_y(19) = cmplx(3.03657864374240d-17, - 4.47655889395416d14, kind=8)
            ground_truth_y(20) = cmplx(1.56113399222733d-18, - 8.25596436773937d15, kind=8)
            ground_truth_y(21) = cmplx(7.63264110088761d-20, - 1.60543649281522d17, kind=8)
            ground_truth_y(22) = cmplx(3.55743345463290d-21, - 3.28288884590347d18, kind=8)
            ground_truth_y(23) = cmplx(1.58408265731184d-22, - 7.04215665376430d19, kind=8)
            ground_truth_y(24) = cmplx(6.75252431874067d-24, - 1.58120235825106d21, kind=8)
            ground_truth_y(25) = cmplx(2.76055759221868d-25, - 3.70878338523624d22, kind=8)
            ground_truth_y(26) = cmplx(1.08417821950868d-26, - 9.07070727024627d23, kind=8)
            ground_truth_y(27) = cmplx(4.09686752846144d-28, - 2.30932157052756d25, kind=8)
            ground_truth_y(28) = cmplx(1.49167553359982d-29, - 6.11063145462780d26, kind=8)
            ground_truth_y(29) = cmplx(5.24018893806175d-31, - 1.67811432845212d28, kind=8)
            ground_truth_y(30) = cmplx(1.77831374777941d-32, - 4.77651520463390d29, kind=8)
            ground_truth_y(31) = cmplx(5.83661788752249d-34, - 1.40739387103855d31, kind=8)
            ground_truth_y(32) = cmplx(1.85470791494552d-35, - 4.28777479146294d32, kind=8)
            ground_truth_y(33) = cmplx(5.71204455589985d-37, - 1.34924166543979d34, kind=8)
            ground_truth_y(34) = cmplx(1.70656572193294d-38, - 4.38074763788785d35, kind=8)
            ground_truth_y(35) = cmplx(4.95061257549874d-40, - 1.46620121702699d37, kind=8)
            ground_truth_y(36) = cmplx(1.39561661412576d-41, - 5.05401345110523d38, kind=8)
            ground_truth_y(37) = cmplx(3.82640464769116d-43, - 1.79270857392533d40, kind=8)
            ground_truth_y(38) = cmplx(1.02108228151588d-44, - 6.53833228137634d41, kind=8)
            ground_truth_y(39) = cmplx(2.65390799339793d-46, - 2.45008189694220d43, kind=8)
            ground_truth_y(40) = cmplx(6.72295942323183d-48, - 9.42627697094610d44, kind=8)
            ground_truth_y(41) = cmplx(1.66097877863811d-49, - 3.72092932162677d46, kind=8)

            y = lib_math_bessel_spherical_third_kind_1_real(x, fnu, n)

            rv = .true.
            print *, "test_lib_math_bessel_spherical_third_kind_1_real:"
            do i=1, n
                buffer = abs(log(real(y(i))) - log(real(ground_truth_y(i))))
                if (buffer .gt. ground_truth_e) then
                    print *, "  ", i , " (Re) difference: ", buffer, " : FAILED"
                    rv = .false.
                else
                    print *, "  ", i, ": (Re) OK"
                end if

                buffer = abs(log(aimag(y(i))) - log(aimag(ground_truth_y(i))))
                if (buffer .gt. ground_truth_e) then
                    print *, "  ", i , " (Im) difference: ", buffer, " : FAILED"
                    rv = .false.
                else
                    print *, "  ", i, ": (Im) OK"
                end if
            end do

        end function test_lib_math_bessel_spherical_third_kind_1_real

        function test_lib_math_bessel_spherical_third_kind_1_cmplx() result (rv)
            implicit none
            ! dummy
            logical :: rv

            ! parameter
            integer, parameter :: n = 41
            double precision, parameter :: ground_truth_e = 10.0_8**(-13.0_8)

            ! auxiliary
            complex(kind=8) :: x = cmplx(5, 2, kind=8)

            complex(kind=8), dimension(n) :: y
            integer :: fnu = 0

            complex(kind=8), dimension(n) :: ground_truth_y

            integer :: i
            double precision :: buffer

            ! Values were generated with sageMath
            !
            ! source code:
            !  >>> x=5.0 +I*2.0
            !  >>> for n in range(0,5):
            !  >>>     value = numerical_approx(spherical_hankel1(n, x))
            !  >>>     print("n = {}: {}".format(n, value))
            ground_truth_y(1)  = cmplx(-0.0250227739998745_8, + 0.00233120915731335_8, kind=8)
            ground_truth_y(2)  = cmplx(-0.00182228917664341_8, + 0.0271504151649199_8, kind=8)
            ground_truth_y(3)  = cmplx(0.0296975379081458_8, + 0.0120891343783301_8, kind=8)
            ground_truth_y(4)  = cmplx(0.0315922819865381_8, - 0.0269692779105477_8, kind=8)
            ground_truth_y(5)  = cmplx(-0.00458857312258838_8, - 0.0598897093673198_8, kind=8)
            ground_truth_y(6)  = cmplx(-0.0758854047150979_8, - 0.0631149498592040_8, kind=8)
            ground_truth_y(7)  = cmplx(-0.187212328816476_8, - 0.00224281954730325_8, kind=8)
            ground_truth_y(8)  = cmplx(-0.345739239467690_8, + 0.225933476709331_8, kind=8)
            ground_truth_y(9)  = cmplx(-0.473216590452379_8, + 0.944214817382837_8, kind=8)
            ground_truth_y(10) = cmplx(0.0657355705905944_8, + 3.09639836994315_8, kind=8)
            ground_truth_y(11) = cmplx(4.74590684093329_8, + 9.11302288820878_8, kind=8)
            ground_truth_y(12) = cmplx(30.3158913122633_8, + 23.0256470073922_8, kind=8)
            ground_truth_y(13) = cmplx(151.995722926319_8, + 34.1083014388942_8, kind=8)
            ground_truth_y(14) = cmplx(683.645503092378_8, - 138.068697574778_8, kind=8)
            ground_truth_y(15) = cmplx(2773.39852701965_8, - 1949.83696832108_8, kind=8)
            ground_truth_y(16) = cmplx(9283.67319536370_8, - 15157.9131980699_8, kind=8)
            ground_truth_y(17) = cmplx(14439.6610247403_8, - 98914.4486804165_8, kind=8)
            ground_truth_y(18) = cmplx(-152242.622982445_8, - 580493.868177847_8, kind=8)
            ground_truth_y(19) = cmplx(-2.33433724014119d6, - 3.03658359691793d6, kind=8)
            ground_truth_y(20) = cmplx(-2.24877427424674d7, - 1.28342306028349d7, kind=8)
            ground_truth_y(21) = cmplx(-1.83396001442696d8, - 2.27782796320593d7, kind=8)
            ground_truth_y(22) = cmplx(-1.33833981676034d9, + 3.70383361421073d8, kind=8)
            ground_truth_y(23) = cmplx(-8.64040025791115d9, + 6.73759369159482d9, kind=8)
            ground_truth_y(24) = cmplx(-4.47894748655319d10, + 7.87190857358492d10, kind=8)
            ground_truth_y(25) = cmplx(-9.91503767845082d10, + 7.76338123042356d11, kind=8)
            ground_truth_y(26) = cmplx(1.83062719024299d12, + 6.81505943392825d12, kind=8)
            ground_truth_y(27) = cmplx(4.01662536792895d13, + 5.27104267820238d13, kind=8)
            ground_truth_y(28) = cmplx(5.57871526737558d14, + 3.28034602884404d14, kind=8)
            ground_truth_y(29) = cmplx(6.49426395908323d15, + 9.41898119841410d14, kind=8)
            ground_truth_y(30) = cmplx(6.69677013698362d16, - 1.66006251953219d16, kind=8)
            ground_truth_y(31) = cmplx(6.07181533663459d17, - 4.42299594473656d17, kind=8)
            ground_truth_y(32) = cmplx(4.45819841040840d18, - 7.18952087209458d18, kind=8)
            ground_truth_y(33) = cmplx(1.65808484454651d19, - 9.70208753855694d19, kind=8)
            ground_truth_y(34) = cmplx(-2.53559510663783d20, - 1.15444133423896d21, kind=8)
            ground_truth_y(35) = cmplx(-8.27994549880366d21, - 1.20671471432732d22, kind=8)
            ground_truth_y(36) = cmplx(-1.55672216449990d23, - 1.03001878858670d23, kind=8)
            ground_truth_y(37) = cmplx(-2.40171673166249d24, - 4.86560861785316d23, kind=8)
            ground_truth_y(38) = cmplx(-3.25224137448419d25, + 6.07044768130982d24, kind=8)
            ground_truth_y(39) = cmplx(-3.86747869548312d26, + 2.47203110593422d26, kind=8)
            ground_truth_y(40) = cmplx(-3.78915519814181d27, + 5.32952850021206d27, kind=8)
            ground_truth_y(41) = cmplx(-2.21874176557107d28, + 9.29890134028607d28, kind=8)

            y = lib_math_bessel_spherical_third_kind_1_cmplx(x, fnu, n)

            rv = .true.
            print *, "test_lib_math_bessel_spherical_third_kind_1_cmplx:"
            do i=1, n
                buffer = abs(log(real(y(i))) - log(real(ground_truth_y(i))))
                if (buffer .gt. ground_truth_e) then
                    print *, "  ", i , " (Re) difference: ", buffer, " : FAILED"
                    rv = .false.
                else
                    print *, "  ", i, ": (Re) OK"
                end if

                buffer = abs(log(aimag(y(i))) - log(aimag(ground_truth_y(i))))
                if (buffer .gt. ground_truth_e) then
                    print *, "  ", i , " (Im) difference: ", buffer, " : FAILED"
                    rv = .false.
                else
                    print *, "  ", i, ": (Im) OK"
                end if
            end do

        end function test_lib_math_bessel_spherical_third_kind_1_cmplx

        function test_lib_math_bessel_spherical_third_kind_2_real() result (rv)
            implicit none
            ! dummy
            logical :: rv

            ! parameter
            integer, parameter :: n = 41

            ! auxiliary
            double precision :: x = 2
            complex(kind=8), dimension(n) :: y
            integer :: fnu = 0

            complex(kind=8), dimension(n) :: ground_truth_y
            double precision :: ground_truth_e = 10.0_8**(-13.0_8)

            integer :: i
            double precision :: buffer

            ! Values were generated with sageMath
            !
            ! source code:
            !  >>> x=2.0 +I*0.0
            !  >>> for n in range(0,5):
            !  >>>     value = numerical_approx(spherical_hankel1(n, x))
            !  >>>     print("n = {}: {}".format(n, value))
            ground_truth_y(1)  = cmplx(0.454648713412841_8, - 0.208073418273571_8, kind=8)
            ground_truth_y(2)  = cmplx(0.435397774979992_8, + 0.350612004276055_8, kind=8)
            ground_truth_y(3)  = cmplx(0.198447949057147_8, + 0.733991424687654_8, kind=8)
            ground_truth_y(4)  = cmplx(0.0607220976628748_8, + 1.48436655744308_8, kind=8)
            ground_truth_y(5)  = cmplx(0.0140793927629153_8, + 4.46129152636313_8, kind=8)
            ground_truth_y(6)  = cmplx(0.00263516977024412_8, + 18.5914453111910_8, kind=8)
            ground_truth_y(7)  = cmplx(0.000414040973427324_8, + 97.7916576851873_8, kind=8)
            ground_truth_y(8)  = cmplx(0.0000560965570334895_8, + 617.054329642526_8, kind=8)
            ground_truth_y(9)  = cmplx(6.68320432384702d-6, + 4530.11581463376_8, kind=8)
            ground_truth_y(10) = cmplx(7.10679719210186d-7, + 37888.9300947444_8, kind=8)
            ground_truth_y(11) = cmplx(6.82530086497472d-8, + 355414.720085438_8, kind=8)
            ground_truth_y(12) = cmplx(5.97687161216011d-9, + 3.69396563080236d6, kind=8)
            ground_truth_y(13) = cmplx(4.81014890094075d-10, + 4.21251900341417d7, kind=8)
            ground_truth_y(14) = cmplx(3.58145140158186d-11, + 5.22870909795969d8, kind=8)
            ground_truth_y(15) = cmplx(2.48104911947672d-12, + 7.01663209221144d9, kind=8)
            ground_truth_y(16) = cmplx(1.60698216593841d-13, + 1.01218294427270d11, kind=8)
            ground_truth_y(17) = cmplx(9.77323772781462d-15, + 1.56186693153047d12, kind=8)
            ground_truth_y(18) = cmplx(5.60205915100117d-16, + 2.56695860758255d13, kind=8)
            ground_truth_y(19) = cmplx(3.03657864374240d-17, + 4.47655889395416d14, kind=8)
            ground_truth_y(20) = cmplx(1.56113399222733d-18, + 8.25596436773937d15, kind=8)
            ground_truth_y(21) = cmplx(7.63264110088761d-20, + 1.60543649281522d17, kind=8)
            ground_truth_y(22) = cmplx(3.55743345463290d-21, + 3.28288884590347d18, kind=8)
            ground_truth_y(23) = cmplx(1.58408265731184d-22, + 7.04215665376430d19, kind=8)
            ground_truth_y(24) = cmplx(6.75252431874067d-24, + 1.58120235825106d21, kind=8)
            ground_truth_y(25) = cmplx(2.76055759221868d-25, + 3.70878338523624d22, kind=8)
            ground_truth_y(26) = cmplx(1.08417821950868d-26, + 9.07070727024627d23, kind=8)
            ground_truth_y(27) = cmplx(4.09686752846144d-28, + 2.30932157052756d25, kind=8)
            ground_truth_y(28) = cmplx(1.49167553359982d-29, + 6.11063145462780d26, kind=8)
            ground_truth_y(29) = cmplx(5.24018893806175d-31, + 1.67811432845212d28, kind=8)
            ground_truth_y(30) = cmplx(1.77831374777941d-32, + 4.77651520463390d29, kind=8)
            ground_truth_y(31) = cmplx(5.83661788752249d-34, + 1.40739387103855d31, kind=8)
            ground_truth_y(32) = cmplx(1.85470791494552d-35, + 4.28777479146294d32, kind=8)
            ground_truth_y(33) = cmplx(5.71204455589985d-37, + 1.34924166543979d34, kind=8)
            ground_truth_y(34) = cmplx(1.70656572193294d-38, + 4.38074763788785d35, kind=8)
            ground_truth_y(35) = cmplx(4.95061257549874d-40, + 1.46620121702699d37, kind=8)
            ground_truth_y(36) = cmplx(1.39561661412576d-41, + 5.05401345110523d38, kind=8)
            ground_truth_y(37) = cmplx(3.82640464769116d-43, + 1.79270857392533d40, kind=8)
            ground_truth_y(38) = cmplx(1.02108228151588d-44, + 6.53833228137634d41, kind=8)
            ground_truth_y(39) = cmplx(2.65390799339793d-46, + 2.45008189694220d43, kind=8)
            ground_truth_y(40) = cmplx(6.72295942323183d-48, + 9.42627697094610d44, kind=8)
            ground_truth_y(41) = cmplx(1.66097877863811d-49, + 3.72092932162677d46, kind=8)

            y = lib_math_bessel_spherical_third_kind_2_real(x, fnu, n)

            rv = .true.
            print *, "test_lib_math_bessel_spherical_third_kind_2_real:"
            do i=1, n
                buffer = abs(log(real(y(i))) - log(real(ground_truth_y(i))))
                if (buffer .gt. ground_truth_e) then
                    print *, "  ", i , " (Re) difference: ", buffer, " : FAILED"
                    rv = .false.
                else
                    print *, "  ", i, ": (Re) OK"
                end if

                buffer = abs(log(aimag(y(i))) - log(aimag(ground_truth_y(i))))
                if (buffer .gt. ground_truth_e) then
                    print *, "  ", i , " (Im) difference: ", buffer, " : FAILED"
                    rv = .false.
                else
                    print *, "  ", i, ": (Im) OK"
                end if
            end do

        end function test_lib_math_bessel_spherical_third_kind_2_real

        function test_lib_math_bessel_spherical_third_kind_2_cmplx() result (rv)
            implicit none
            ! dummy
            logical :: rv

            ! parameter
            integer, parameter :: n = 41

            ! auxiliary
            complex(kind=8) :: x = cmplx(2_8, 2_8, kind=8)
            complex(kind=8), dimension(n) :: y
            integer :: fnu = 0

            complex(kind=8), dimension(n) :: ground_truth_y
            double precision :: ground_truth_e = 10.0_8**(-13.0_8)

            integer :: i
            double precision :: buffer

            ! Values were generated with sageMath
            !
            ! source code:
            !  >>> x=2.0 +I*2.0
            !  >>> for n in range(0,5):
            !  >>>     value = numerical_approx(spherical_hankel2(n, x))
            !  >>>     print("n = {}: {}".format(n, value))
            ground_truth_y(1)  = cmplx(0.910979344197223_8, - 2.44844550451690_8, kind=8)
            ground_truth_y(2)  = cmplx(2.06407896443698_8, + 0.0711231320186915_8, kind=8)
            ground_truth_y(3)  = cmplx(0.690422228144533_8, + 0.953728630203184_8, kind=8)
            ground_truth_y(4)  = cmplx(-0.00889039150233644_8, + 0.258009870554623_8, kind=8)
            ground_truth_y(5)  = cmplx(-0.254463139803031_8, - 0.486653171603505_8, kind=8)
            ground_truth_y(6)  = cmplx(-1.65862130916237_8, - 0.780437442105689_8, kind=8)
            ground_truth_y(7)  = cmplx(-6.45294842618414_8, + 2.90165880600938_8, kind=8)
            ground_truth_y(8)  = cmplx(-9.88306995640557_8, + 31.1829109467346_8, kind=8)
            ground_truth_y(9)  = cmplx(86.3273521399181_8, + 151.095769580766_8, kind=8)
            ground_truth_y(10) = cmplx(1018.93133726931_8, + 244.082863176871_8, kind=8)
            ground_truth_y(11) = cmplx(5912.99009997946_8, - 3831.62602151987_8, kind=8)
            ground_truth_y(12) = cmplx(9908.23007464352_8, - 51403.3175010484_8, kind=8)
            ground_truth_y(13) = cmplx(-244509.742801807_8, - 348709.772538709_8, kind=8)
            ground_truth_y(14) = cmplx(-3.71753020095287d6, - 599846.868354584_8, kind=8)
            ground_truth_y(15) = cmplx(-2.88977854750235d7, + 2.13930722675771d7, kind=8)
            ground_truth_y(16) = cmplx(-5.06916405530332d7, + 3.65208565502209d8, kind=8)
            ground_truth_y(17) = cmplx(2.46640395383114d9, + 3.20183352466055d9, kind=8)
            ground_truth_y(18) = cmplx(4.68136508381095d10, + 5.70208539384046d9, kind=8)
            ground_truth_y(19) = cmplx(4.57046288075731d11, - 3.62928031162014d11, kind=8)
            ground_truth_y(20) = cmplx(8.23780225613767d11, - 7.59046453834298d12, kind=8)
            ground_truth_y(21) = cmplx(-6.64322183371856d13, - 8.16759584174163d13, kind=8)
            ground_truth_y(22) = cmplx(-1.51893259196028d15, - 1.48657871284022d14, kind=8)
            ground_truth_y(23) = cmplx(-1.78601652615391d16, + 1.48121292056872d16, kind=8)
            ground_truth_y(24) = cmplx(-3.27714730363732d16, + 3.67711970627580d17, kind=8)
            ground_truth_y(25) = cmplx(3.95341101195822d18, + 4.69086833384576d18, kind=8)
            ground_truth_y(26) = cmplx(1.05925193459135d20, + 8.66614022249484d18, kind=8)
            ground_truth_y(27) = cmplx(1.45708609342882d21, - 1.24474379710101d21, kind=8)
            ground_truth_y(28) = cmplx(2.70761023288441d21, - 3.58079121897428d22, kind=8)
            ground_truth_y(29) = cmplx(-4.56586238000232d23, - 5.28343689514023d23, kind=8)
            ground_truth_y(30) = cmplx(-1.40379590773110d25, - 9.86735771881784d23, kind=8)
            ground_truth_y(31) = cmplx(-2.21157662787594d26, + 1.93033887444595d26, kind=8)
            ground_truth_y(32) = cmplx(-4.14849614903414d26, + 6.31740787681276d27, kind=8)
            ground_truth_y(33) = cmplx(9.31864502878598d28, + 1.05840021607085d29, kind=8)
            ground_truth_y(34) = cmplx(3.23459501790776d30, + 1.99303126060599d29, kind=8)
            ground_truth_y(35) = cmplx(5.74246074611821d31, - 5.09469792100470d31, kind=8)
            ground_truth_y(36) = cmplx(1.08504492314173d32, - 1.86960917320476d33, kind=8)
            ground_truth_y(37) = cmplx(-3.13170326932691d34, - 3.50605705837511d34, kind=8)
            ground_truth_y(38) = cmplx(-1.21149976429793d36, - 6.64499573280903d34, kind=8)
            ground_truth_y(39) = cmplx(-2.39302402477947d37, + 2.15047444512683d37, kind=8)
            ground_truth_y(40) = cmplx(-4.54792943188346d37, + 8.74689905414290d38, kind=8)
            ground_truth_y(41) = cmplx(1.64008398093830d40, + 1.81518369502779d40, kind=8)

            y = lib_math_bessel_spherical_third_kind_2_cmplx(x, fnu, n)

            rv = .true.
            print *, "test_lib_math_bessel_spherical_third_kind_2_cmplx:"
            do i=1, n
                buffer = abs(log(real(y(i))) - log(real(ground_truth_y(i))))
                if (buffer .gt. ground_truth_e) then
                    print *, "  ", i , " (Re) difference: ", buffer, " : FAILED"
                    rv = .false.
                else
                    print *, "  ", i, ": (Re) OK"
                end if

                buffer = abs(log(aimag(y(i))) - log(aimag(ground_truth_y(i))))
                if (buffer .gt. ground_truth_e) then
                    print *, "  ", i , " (Im) difference: ", buffer, " : FAILED"
                    rv = .false.
                else
                    print *, "  ", i, ": (Im) OK"
                end if
            end do

        end function test_lib_math_bessel_spherical_third_kind_2_cmplx

        function test_lib_math_bessel_riccati_s_real() result(rv)
            implicit none
            ! dummy
            logical :: rv

            ! parameter
            integer, parameter :: n = 5

            ! auxiliary
            double precision :: x = 4
            double precision, dimension(n) :: s
            integer :: fnu = 1

            double precision, dimension(n) :: ground_truth_s
            double precision :: ground_truth_e = 10.0_8**(-13.0_8)

            integer :: i
            double precision :: buffer

            ! Values were generated with sageMath
            !
            ! source code:
            !  >>> var('x')
            !  >>> var('n')
            !  >>> S(n,x) = x*spherical_bessel_J(n,x)
            !  >>>
            !  >>> from IPython.display import Math
            !  >>> t = latex(S)
            !  >>> display(Math("""{}""".format(t)))
            !  >>>
            !  >>> x=4.0
            !  >>> for i in range(1,6):
            !  >>>     value = numerical_approx(S(i,x))
            !  >>>     print("n = {}: {}".format(i, value))

            ground_truth_s(1) = 0.464442997036630_8
            ground_truth_s(2) = 1.10513474308540_8
            ground_truth_s(3) = 0.916975431820121_8
            ground_truth_s(4) = 0.499572262599811_8
            ground_truth_s(5) = 0.207062159029454_8

            s = lib_math_bessel_riccati_s_real(x, fnu, n)

            rv = .true.
            print *, "test_lib_math_bessel_riccati_s_real:"
            do i=1, n
                buffer = abs(s(i) - ground_truth_s(i))
                if (buffer .gt. ground_truth_e) then
                    print *, "  ", i , "difference: ", buffer, " : FAILED"
                    rv = .false.
                else
                    print *, "  ", i, ": OK"
                end if
            end do
        end function test_lib_math_bessel_riccati_s_real

        function test_lib_math_bessel_riccati_s_cmplx() result(rv)
            implicit none
            ! dummy
            logical :: rv

            ! parameter
            integer, parameter :: n = 5

            ! auxiliary
            complex(kind=8) :: x = cmplx(4, 2, kind = 8)
            complex(kind=8), dimension(n) :: s
            integer :: fnu = 1

            complex(kind=8), dimension(n) :: ground_truth_s
            double precision :: ground_truth_e = 10.0_8**(-13.0_8)

            integer :: i
            double precision :: buffer

            ! Values were generated with sageMath
            !
            ! source code:
            !  >>> var('x')
            !  >>> var('n')
            !  >>> S(n,x) = x*spherical_bessel_J(n,x)
            !  >>>
            !  >>> from IPython.display import Math
            !  >>> t = latex(S)
            !  >>> display(Math("""{}""".format(t)))
            !  >>>
            !  >>> x=4.0 + I*2
            !  >>> for i in range(1,6):
            !  >>>     value = numerical_approx(S(i,x))
            !  >>>     print("n = {}: {}".format(i, value))

            ground_truth_s(1) = cmplx(1.65261997961242_8, -2.93422793197767_8, kind = 8)
            ground_truth_s(2) = cmplx(2.95854269502298_8, +0.114351416281673_8, kind = 8)
            ground_truth_s(3) = cmplx(1.36309842355140_8, +1.56930800074786_8, kind = 8)
            ground_truth_s(4) = cmplx(0.0483106984724765_8, +1.12851088827935_8, kind = 8)
            ground_truth_s(5) = cmplx(-0.260479366849524_8, +0.418531969529742_8, kind = 8)

            s = lib_math_bessel_riccati_s_cmplx(x, fnu, n)

            rv = .true.
            print *, "test_lib_math_bessel_riccati_s_cmplx:"
            do i=1, n
                buffer = abs(s(i) - ground_truth_s(i))
                if (buffer .gt. ground_truth_e) then
                    print *, "  ", i , "difference: ", buffer, " : FAILED"
                    rv = .false.
                else
                    print *, "  ", i, ": OK"
                end if
            end do
        end function test_lib_math_bessel_riccati_s_cmplx

        function test_lib_math_bessel_riccati_c_real() result(rv)
            implicit none
            ! dummy
            logical :: rv

            ! parameter
            integer, parameter :: n = 5

            ! auxiliary
            double precision :: x = 4
            double precision, dimension(n) :: c
            integer :: fnu = 1

            double precision, dimension(n) :: ground_truth_c
            double precision :: ground_truth_e = 10.0_8**(-13.0_8)

            integer :: i
            double precision :: buffer

            ! Values were generated with sageMath
            !
            ! source code:
            !  >>> var('x')
            !  >>> var('n')
            !  >>> C(n,x) = -x*spherical_bessel_Y(n,x)
            !  >>>
            !  >>> from IPython.display import Math
            !  >>> t = latex(C)
            !  >>> display(Math("""{}""".format(t)))
            !  >>>
            !  >>> x=4.0
            !  >>> for i in range(1,6):
            !  >>>     value = numerical_approx(C(i,x))
            !  >>>     print("n = {}: {}".format(i, value))

            ground_truth_c(1) = -0.920213400523831_8
            ground_truth_c(2) = -0.0365164295292615_8
            ground_truth_c(3) = 0.874567863612254_8
            ground_truth_c(4) = 1.56701019085071_8
            ground_truth_c(5) = 2.65120506580184_8

            c = lib_math_bessel_riccati_c_real(x, fnu, n)

            rv = .true.
            print *, "test_lib_math_bessel_riccati_c_real:"
            do i=1, n
                buffer = abs(c(i) - ground_truth_c(i))
                if (buffer .gt. ground_truth_e) then
                    print *, "  ", i , "difference: ", buffer, " : FAILED"
                    rv = .false.
                else
                    print *, "  ", i, ": OK"
                end if
            end do
        end function test_lib_math_bessel_riccati_c_real

        function test_lib_math_bessel_riccati_c_cmplx() result(rv)
            implicit none
            ! dummy
            logical :: rv

            ! parameter
            integer, parameter :: n = 5

            ! auxiliary
            complex(kind=8) :: x = cmplx(4, 2, kind = 8)
            complex(kind=8), dimension(n) :: c
            integer :: fnu = 1

            complex(kind=8), dimension(n) :: ground_truth_c
            double precision :: ground_truth_e = 10.0_8**(-13.0_8)

            integer :: i
            double precision :: buffer

            ! Values were generated with sageMath
            !
            ! source code:
            !  >>> var('x')
            !  >>> var('n')
            !  >>> C(n,x) = -x*spherical_bessel_Y(n,x)
            !  >>>
            !  >>> from IPython.display import Math
            !  >>> t = latex(C)
            !  >>> display(Math("""{}""".format(t)))
            !  >>>
            !  >>> x=4.0 + I*2.0
            !  >>> for i in range(1,6):
            !  >>>     value = numerical_approx(C(i,x))
            !  >>>     print("n = {}: {}".format(i, value))

            ground_truth_c(1) = cmplx(-3.06458442895309_8, -1.57579724660183_8, kind = 8)
            ground_truth_c(2) = cmplx(0.147645382564980_8, -2.77092002606733_8, kind = 8)
            ground_truth_c(3) = cmplx(1.82676979848441_8, -1.26894547074798_8, kind = 8)
            ground_truth_c(4) = cmplx(1.52157050578960_8, -0.284342491918936_8, kind = 8)
            ground_truth_c(5) = cmplx(0.656148869209830_8, -0.612284469916739_8, kind = 8)

            c = lib_math_bessel_riccati_c_cmplx(x, fnu, n)

            rv = .true.
            print *, "test_lib_math_bessel_riccati_c_cmplx:"
            do i=1, n
                buffer = abs(c(i) - ground_truth_c(i))
                if (buffer .gt. ground_truth_e) then
                    print *, "  ", i , "difference: ", buffer, " : FAILED"
                    rv = .false.
                else
                    print *, "  ", i, ": OK"
                end if
            end do
        end function test_lib_math_bessel_riccati_c_cmplx

        function test_lib_math_bessel_riccati_xi_real() result(rv)
            implicit none
            ! dummy
            logical :: rv

            ! parameter
            integer, parameter :: n = 5

            ! auxiliary
            double precision :: x = 4
            complex(kind=8), dimension(n) :: xi
            integer :: fnu = 1

            complex(kind=8), dimension(n) :: ground_truth_xi
            double precision :: ground_truth_e = 10.0_8**(-13.0_8)

            integer :: i
            double precision :: buffer

            ! Values were generated with sageMath
            !
            ! source code:
            !  >>> var('x')
            !  >>> var('n')
            !  >>> S(n,x) = x*spherical_bessel_J(n,x)
            !  >>>
            !  >>> from IPython.display import Math
            !  >>> t = latex(S)
            !  >>> display(Math("""{}""".format(t)))
            !  >>>
            !  >>> x=4.0
            !  >>> for i in range(1,6):
            !  >>>     value = numerical_approx(S(i,x))
            !  >>>     print("n = {}: {}".format(i, value))

            ground_truth_xi(1) = cmplx(0.464442997036630_8, +0.920213400523831_8, kind=8)
            ground_truth_xi(2) = cmplx(1.10513474308540_8, +0.0365164295292615_8, kind=8)
            ground_truth_xi(3) = cmplx(0.916975431820121_8, -0.874567863612254_8, kind=8)
            ground_truth_xi(4) = cmplx(0.499572262599811_8, -1.56701019085071_8, kind=8)
            ground_truth_xi(5) = cmplx(0.207062159029454_8, -2.65120506580184_8, kind=8)

            xi = lib_math_bessel_riccati_xi_real(x, fnu, n)

            rv = .true.
            print *, "test_lib_math_bessel_riccati_xi_real:"
            do i=1, n
                buffer = abs(xi(i) - ground_truth_xi(i))
                if (buffer .gt. ground_truth_e) then
                    print *, "  ", i , "difference: ", buffer, " : FAILED"
                    rv = .false.
                else
                    print *, "  ", i, ": OK"
                end if
            end do
        end function test_lib_math_bessel_riccati_xi_real

        function test_lib_math_bessel_riccati_xi_cmplx() result(rv)
            implicit none
            ! dummy
            logical :: rv

            ! parameter
            integer, parameter :: n = 5

            ! auxiliary
            complex(kind=8) :: x = cmplx(4, 2)
            complex(kind=8), dimension(n) :: xi
            integer :: fnu = 1

            complex(kind=8), dimension(n) :: ground_truth_xi
            double precision :: ground_truth_e = 10.0_8**(-13.0_8)

            integer :: i
            double precision :: buffer

            ! Values were generated with sageMath
            !
            ! source code:
            !  >>> var('x')
            !  >>> var('n')
            !  >>> S(n,x) = x*spherical_bessel_J(n,x)
            !  >>>
            !  >>> from IPython.display import Math
            !  >>> t = latex(S)
            !  >>> display(Math("""{}""".format(t)))
            !  >>>
            !  >>> x=4.0+2.0*I
            !  >>> for i in range(1,6):
            !  >>>     value = numerical_approx(S(i,x))
            !  >>>     print("n = {}: {}".format(i, value))

            ground_truth_xi(1) = cmplx(0.0768227330105855_8, +0.130356496975417_8, kind=8)
            ground_truth_xi(2) = cmplx(0.187622668955650_8, -0.0332939662833071_8, kind=8)
            ground_truth_xi(3) = cmplx(0.0941529528034112_8, -0.257461797736550_8, kind=8)
            ground_truth_xi(4) = cmplx(-0.236031793446459_8, -0.393059617510250_8, kind=8)
            ground_truth_xi(5) = cmplx(-0.872763836766263_8, -0.237616899680087_8, kind=8)

            xi = lib_math_bessel_riccati_xi_cmplx(x, fnu, n)

            rv = .true.
            print *, "test_lib_math_bessel_riccati_xi_cmplx:"
            do i=1, n
                buffer = abs(xi(i) - ground_truth_xi(i))
                if (buffer .gt. ground_truth_e) then
                    print *, "  ", i , "difference: ", buffer, " : FAILED"
                    rv = .false.
                else
                    print *, "  ", i, ": OK"
                end if
            end do
        end function test_lib_math_bessel_riccati_xi_cmplx

        function test_lib_math_bessel_riccati_zeta_real() result(rv)
            implicit none
            ! dummy
            logical :: rv

            ! parameter
            integer, parameter :: n = 5

            ! auxiliary
            double precision :: x = 4
            complex(kind=8), dimension(n) :: zeta
            integer :: fnu = 1

            complex(kind=8), dimension(n) :: ground_truth_zeta
            double precision :: ground_truth_e = 10.0_8**(-13.0_8)

            integer :: i
            double precision :: buffer

            ! Values were generated with sageMath
            !
            ! source code:
            !  >>> var('x')
            !  >>> var('n')
            !  >>> S(n,x) = x*spherical_bessel_J(n,x)
            !  >>>
            !  >>> from IPython.display import Math
            !  >>> t = latex(S)
            !  >>> display(Math("""{}""".format(t)))
            !  >>>
            !  >>> x=4.0
            !  >>> for i in range(1,6):
            !  >>>     value = numerical_approx(S(i,x))
            !  >>>     print("n = {}: {}".format(i, value))

            ground_truth_zeta(1) = cmplx(0.464442997036630_8, -0.920213400523831_8, kind=8)
            ground_truth_zeta(2) = cmplx(1.10513474308540_8, -0.0365164295292615_8, kind=8)
            ground_truth_zeta(3) = cmplx(0.916975431820121_8, +0.874567863612254_8, kind=8)
            ground_truth_zeta(4) = cmplx(0.499572262599811_8, +1.56701019085071_8, kind=8)
            ground_truth_zeta(5) = cmplx(0.207062159029454_8, +2.65120506580184_8, kind=8)

            zeta = lib_math_bessel_riccati_zeta_real(x, fnu, n)

            rv = .true.
            print *, "test_lib_math_bessel_riccati_zeta_real:"
            do i=1, n
                buffer = abs(zeta(i) - ground_truth_zeta(i))
                if (buffer .gt. ground_truth_e) then
                    print *, "  ", i , "difference: ", buffer, " : FAILED"
                    rv = .false.
                else
                    print *, "  ", i, ": OK"
                end if
            end do
        end function test_lib_math_bessel_riccati_zeta_real

        function test_lib_math_bessel_riccati_zeta_cmplx() result(rv)
            implicit none
            ! dummy
            logical :: rv

            ! parameter
            integer, parameter :: n = 5

            ! auxiliary
            complex(kind=8) :: z = cmplx(4, 2)
            complex(kind=8), dimension(n) :: zeta
            integer :: fnu = 1

            complex(kind=8), dimension(n) :: ground_truth_zeta
            double precision :: ground_truth_e = 10.0_8**(-13.0_8)

            integer :: i
            double precision :: buffer

            ! Values were generated with sageMath
            !
            ! source code:
            !  >>> var('x')
            !  >>> var('n')
            !  >>> S(n,x) = x*spherical_bessel_J(n,x)
            !  >>>
            !  >>> from IPython.display import Math
            !  >>> t = latex(S)
            !  >>> display(Math("""{}""".format(t)))
            !  >>>
            !  >>> x=4.0+2.0*I
            !  >>> for i in range(1,6):
            !  >>>     value = numerical_approx(S(i,x))
            !  >>>     print("n = {}: {}".format(i, value))

            ground_truth_zeta(1) = cmplx(3.22841722621425_8, -5.99881236093076_8, kind=8)
            ground_truth_zeta(2) = cmplx(5.72946272109030_8, +0.261996798846654_8, kind=8)
            ground_truth_zeta(3) = cmplx(2.63204389429938_8, +3.39607779923226_8, kind=8)
            ground_truth_zeta(4) = cmplx(0.332653190391413_8, +2.65008139406895_8, kind=8)
            ground_truth_zeta(5) = cmplx(0.351805103067215_8, +1.07468083873957_8, kind=8)

            zeta = lib_math_bessel_riccati_zeta_cmplx(z, fnu, n)

            rv = .true.
            print *, "test_lib_math_bessel_riccati_zeta_cmplx:"
            do i=1, n
                buffer = abs(zeta(i) - ground_truth_zeta(i))
                if (buffer .gt. ground_truth_e) then
                    print *, "  ", i , "difference: ", buffer, " : FAILED"
                    rv = .false.
                else
                    print *, "  ", i, ": OK"
                end if
            end do
        end function test_lib_math_bessel_riccati_zeta_cmplx

        function test_lib_math_bessel_spherical_first_kind_derivative_real() result (rv)
            implicit none
            ! dummy
            logical :: rv

            ! parameter
            integer, parameter :: n = 5
            double precision, parameter :: ground_truth_e = 10.0_8**(-13.0_8)

            ! auxiliary
            double precision :: x = 4
            double precision, dimension(n) :: y
            double precision, dimension(n+1) :: y_n
            integer :: fnu = 1

            double precision, dimension(n) :: ground_truth_y
            double precision, dimension(n) :: ground_truth_y_n

            integer :: i
            double precision :: buffer
            double precision :: buffer_y_n

            ! Values were generated with sageMath
            !
            ! source code:
            !  >>> for n in range(0,5):
            !  >>>     var('x')
            !  >>>     f(x) = derivative(spherical_bessel_J(n, x))
            !  >>>     x=4.0
            !  >>>     value = numerical_approx(f(x))
            !  >>>     print("n = {}: {}".format(n, value))
            data ground_truth_y / -0.247255998456561_8, &
                                  -0.0911020150693551_8, &
                                  0.0470398278163199_8, &
                                  0.0731275258925893_8, &
                                  0.0472447560139076_8/

            ground_truth_y_n = lib_math_bessel_spherical_first_kind_real(x, fnu, n)

            y = lib_math_bessel_spherical_first_kind_derivative_real(x, fnu, n, y_n)

            rv = .true.
            print *, "test_lib_math_bessel_spherical_first_kind_derivative_real:"
            do i=1, n
                buffer = y(i) - ground_truth_y(i)
                if (abs(buffer) .gt. ground_truth_e) then
                    print *, "  ", i , "difference: ", buffer, " : FAILED"
                    rv = .false.
                else
                    print *, "  ", i, ": OK"
                end if

                buffer_y_n = abs(y_n(i) - ground_truth_y_n(i))
                if (buffer_y_n .gt. ground_truth_e) then
                    print *, "  ", i, "antiderivative: difference: ", buffer_y_n, " : FAILED"
                    rv = .false.
                end if
            end do

        end function test_lib_math_bessel_spherical_first_kind_derivative_real

        function test_lib_math_bessel_spherical_first_kind_derivative_coeff_real() result (rv)
            implicit none
            ! dummy
            logical :: rv

            ! parameter
            integer, parameter :: n = 5
            double precision, parameter :: ground_truth_e = 10.0_8**(-6.0_8)

            ! auxiliary
            double precision :: a = 1.5
            double precision :: x = 4
            double precision, dimension(n) :: y_d
            double precision, dimension(n+1) :: y_n
            integer :: fnu = 1

            double precision, dimension(n) :: ground_truth_y_d
            double precision, dimension(n) :: ground_truth_y_n

            integer :: i
            double precision :: buffer
            double precision :: buffer_y_n

            ! Values were generated with mathematica
            !
            ! source code:
            !  >>> dn2[n1_, x1_, t_] := FullSimplify[D[SphericalBesselJ[n1, t y], y] /. y -> x1];
            !  >>> n = {1, 2, 3, 4, 5};
            !  >>> a = 1.5;
            !  >>> x = 4;
            !  >>> N[jdn2[n, x, a], 16]
            ground_truth_y_d = (/ 0.0140411_8, &
                             -0.223691_8, &
                             -0.192674_8, &
                             -0.0409619_8, &
                             0.0574339_8/)

            ground_truth_y_n = lib_math_bessel_spherical_first_kind_real(a*x, fnu, n)

            y_d = lib_math_bessel_spherical_first_kind_derivative_coeff_real(a, x, fnu, n, y_n)

            rv = .true.
            print *, "test_lib_math_bessel_spherical_first_kind_derivative_coeff_real:"
            do i=1, n
                buffer = y_d(i) - ground_truth_y_d(i)
                if (abs(buffer) .gt. ground_truth_e) then
                    print *, "  ", i , "difference: ", buffer, " : FAILED"
                    rv = .false.
                else
                    print *, "  ", i, ": OK"
                end if

                buffer_y_n = abs(y_n(i) - ground_truth_y_n(i))
                if (buffer_y_n .gt. ground_truth_e) then
                    print *, "  ", i, "antiderivative: difference: ", buffer_y_n, " : FAILED"
                    rv = .false.
                end if
            end do

        end function test_lib_math_bessel_spherical_first_kind_derivative_coeff_real

        function test_lib_math_bessel_spherical_first_kind_derivative_cmplx() result (rv)
            implicit none
            ! dummy
            logical :: rv

            ! parameter
            integer, parameter :: n = 5
            double precision, parameter :: ground_truth_e = 10.0_8**(-13.0_8)

            ! auxiliary
            complex(kind=8) :: x = cmplx(4, 2, kind = 8)
            complex(kind=8), dimension(n) :: y
            complex(kind=8), dimension(n) :: y_n
            integer :: fnu = 1

            complex(kind=8), dimension(n) :: ground_truth_y
            complex(kind=8), dimension(n) :: ground_truth_y_n

            integer :: i
            double precision :: buffer
            double precision :: buffer_y_n

            ! Values were generated with sageMath
            !
            ! source code:
            !  >>> for n in range(1,6):
            !  >>>     var('x')
            !  >>>     f(x) = derivative(spherical_bessel_J(n, x))
            !  >>>     x=4.0 + I*2
            !  >>>     value = numerical_approx(f(x))
            !  >>>     print("n = {}: {}".format(n, value))
            ground_truth_y(1) = cmplx(-0.670934198523497_8, +0.118852349102136_8, kind=8)
            ground_truth_y(2) = cmplx(-0.242889809781152_8, -0.407374088419369_8, kind=8)
            ground_truth_y(3) = cmplx(0.188482589686938_8, -0.243205198567482_8, kind=8)
            ground_truth_y(4) = cmplx(0.196601702358324_8, +0.0179372642470247_8, kind=8)
            ground_truth_y(5) = cmplx(0.0689518418682064_8, +0.0830203052493826_8, kind=8)

            y = lib_math_bessel_spherical_first_kind_derivative_cmplx(x, fnu, n, y_n)
            ground_truth_y_n = lib_math_bessel_spherical_first_kind_cmplx(x, fnu, n)

            rv = .true.
            print *, "test_lib_math_bessel_spherical_first_kind_derivative_cmplx:"
            do i=1, n
                buffer = abs(y(i) - ground_truth_y(i))
                if (buffer .gt. ground_truth_e) then
                    print *, "  ", i , "difference: ", buffer, " : FAILED"
                    rv = .false.
                else
                    print *, "  ", i, ": OK"
                end if

                buffer_y_n = abs(y_n(i) - ground_truth_y_n(i))
                if (buffer_y_n .gt. ground_truth_e) then
                    print *, "  ", i, "antiderivative: difference: ", buffer_y_n, " : FAILED"
                    rv = .false.
                end if
            end do

        end function test_lib_math_bessel_spherical_first_kind_derivative_cmplx

        function test_lib_math_bessel_spherical_first_kind_derivative_coeff_c() result (rv)
            implicit none
            ! dummy
            logical :: rv

            ! parameter
            integer, parameter :: n = 5
            double precision, parameter :: ground_truth_e = 10.0_8**(-5.0_8)

            ! auxiliary
            real(kind=8) :: a = 1.5
            complex(kind=8) :: x = cmplx(4, 2, kind = 8)
            complex(kind=8), dimension(n) :: y_d
            complex(kind=8), dimension(n) :: y_n
            integer :: fnu = 1

            complex(kind=8), dimension(n) :: ground_truth_y_d
            complex(kind=8), dimension(n) :: ground_truth_y_n

            integer :: i
            double precision :: buffer
            double precision :: buffer_y_n

            ! Values were generated with sageMath
            !
            ! source code:
            !  >>> for n in range(1,6):
            !  >>>     var('x')
            !  >>>     f(x) = derivative(spherical_bessel_J(n, x))
            !  >>>     x=4.0 + I*2
            !  >>>     value = numerical_approx(f(x))
            !  >>>     print("n = {}: {}".format(n, value))
            ground_truth_y_d(1) = cmplx(0.846179_8, +1.75688_8, kind=8)
            ground_truth_y_d(2) = cmplx(-1.2786_8, +1.07342_8, kind=8)
            ground_truth_y_d(3) = cmplx(-1.13322_8, -0.607966_8, kind=8)
            ground_truth_y_d(4) = cmplx(0.0302722_8, -0.888809_8, kind=8)
            ground_truth_y_d(5) = cmplx(0.508252_8, -0.277546_8, kind=8)

            y_d = lib_math_bessel_spherical_first_kind_derivative_coeff_cmplx(a, x, fnu, n, y_n)
            ground_truth_y_n = lib_math_bessel_spherical_first_kind_cmplx(a * x, fnu, n)

            rv = .true.
            print *, "test_lib_math_bessel_spherical_first_kind_derivative_coeff_cmplx:"
            do i=1, n
                buffer = abs(y_d(i) - ground_truth_y_d(i))
                if (buffer .gt. ground_truth_e) then
                    print *, "  ", i , "difference: ", buffer, " : FAILED"
                    rv = .false.
                else
                    print *, "  ", i, ": OK"
                end if

                buffer_y_n = abs(y_n(i) - ground_truth_y_n(i))
                if (buffer_y_n .gt. ground_truth_e) then
                    print *, "  ", i, "antiderivative: difference: ", buffer_y_n, " : FAILED"
                    rv = .false.
                end if
            end do

        end function test_lib_math_bessel_spherical_first_kind_derivative_coeff_c

        function test_lib_math_bessel_spherical_second_kind_derivative_real() result (rv)
            implicit none
            ! dummy
            logical :: rv

            ! parameter
            integer, parameter :: n = 5
            double precision, parameter :: ground_truth_e = 10.0_8**(-13.0_8)

            ! auxiliary
            double precision :: x = 4
            double precision, dimension(n) :: y
            double precision, dimension(n+1) :: y_n
            integer :: fnu = 1

            double precision, dimension(n) :: ground_truth_y
            double precision, dimension(n) :: ground_truth_y_n

            integer :: i
            double precision :: buffer
            double precision :: buffer_y_n

            ! Values were generated with sageMath
            !
            ! source code:
            !  >>> for n in range(1,6):
            !  >>>     var('x')
            !  >>>     f(x) = derivative(spherical_bessel_Y(n, x))
            !  >>>     x=4.0
            !  >>>     value = numerical_approx(f(x))
            !  >>>     print("n = {}: {}".format(n, value))
            data ground_truth_y / 0.0483842301504241_8, &
                                  0.223206519594221_8, &
                                  0.227771073285379_8, &
                                  0.271048718737782_8, &
                                  0.602449351963012_8 /

            y = lib_math_bessel_spherical_second_kind_derivative_real(x, fnu, n, y_n)
            ground_truth_y_n = lib_math_bessel_spherical_second_kind_real(x, fnu, n)

            rv = .true.
            print *, "test_lib_math_bessel_spherical_second_kind_derivative_real:"
            do i=1, n
                buffer = y(i) - ground_truth_y(i)
                if (abs(buffer) .gt. ground_truth_e) then
                    print *, "  ", i , "difference: ", buffer, " : FAILED"
                    rv = .false.
                else
                    print *, "  ", i, ": OK"
                end if

                buffer_y_n = abs(y_n(i) - ground_truth_y_n(i))
                if (buffer_y_n .gt. ground_truth_e) then
                    print *, "  ", i, "antiderivative: difference: ", buffer_y_n, " : FAILED"
                    rv = .false.
                end if
            end do

        end function test_lib_math_bessel_spherical_second_kind_derivative_real

        function test_lib_math_bessel_spherical_second_kind_derivative_cmplx() result (rv)
            implicit none
            ! dummy
            logical :: rv

            ! parameter
            integer, parameter :: n = 5
            double precision, parameter :: ground_truth_e = 10.0_8**(-13.0_8)

            ! auxiliary
            complex(kind=8) :: x = cmplx(4, 2, kind=8)
            complex(kind=8), dimension(n) :: y
            complex(kind=8), dimension(n+1) :: y_n
            integer :: fnu = 1

            complex(kind=8), dimension(n) :: ground_truth_y
            complex(kind=8), dimension(n) :: ground_truth_y_n

            integer :: i
            double precision :: buffer
            double precision :: buffer_y_n

            ! Values were generated with sageMath
            !
            ! source code:
            !  >>> for n in range(1,6):
            !  >>>     var('x')
            !  >>>     f(x) = derivative(spherical_bessel_Y(n, x))
            !  >>>     x=4.0 + I*2
            !  >>>     value = numerical_approx(f(x))
            !  >>>     print("n = {}: {}".format(n, value))
            ground_truth_y(1) = cmplx(-0.0925935033610706_8, -0.644258003230032_8, kind=8)
            ground_truth_y(2) = cmplx(0.451274291753570_8, -0.258399241828799_8, kind=8)
            ground_truth_y(3) = cmplx(0.263744026592188_8, +0.124391919222700_8, kind=8)
            ground_truth_y(4) = cmplx(-0.0670923351374301_8, +0.0895005990522773_8, kind=8)
            ground_truth_y(5) = cmplx(-0.304721328288274_8, -0.0586613842326251_8, kind=8)

            y = lib_math_bessel_spherical_second_kind_derivative_cmplx(x, fnu, n, y_n)
            ground_truth_y_n = lib_math_bessel_spherical_second_kind_cmplx(x, fnu, n)

            rv = .true.
            print *, "test_lib_math_bessel_spherical_second_kind_derivative_cmplx:"
            do i=1, n
                buffer = abs(y(i) - ground_truth_y(i))
                if (buffer .gt. ground_truth_e) then
                    print *, "  ", i , "difference: ", buffer, " : FAILED"
                    rv = .false.
                else
                    print *, "  ", i, ": OK"
                end if

                buffer_y_n = abs(y_n(i) - ground_truth_y_n(i))
                if (buffer_y_n .gt. ground_truth_e) then
                    print *, "  ", i, "antiderivative: difference: ", buffer_y_n, " : FAILED"
                    rv = .false.
                end if
            end do

        end function test_lib_math_bessel_spherical_second_kind_derivative_cmplx

        function test_lib_math_bessel_spherical_third_kind_1_derivative_real() result (rv)
            implicit none
            ! dummy
            logical :: rv

            ! parameter
            integer, parameter :: n = 5
            double precision, parameter :: ground_truth_e = 10.0_8**(-13.0_8)

            ! auxiliary
            double precision :: x = 4
            complex(kind=8), dimension(n) :: h_1
            complex(kind=8), dimension(n) :: h_1_n
            integer :: fnu = 1

            complex(kind=8), dimension(n) :: ground_truth_h_1
            complex(kind=8), dimension(n) :: ground_truth_h_1_n

            integer :: i
            double precision :: buffer
            double precision :: buffer_y_n

            ! Values were generated with sageMath
            !
            ! source code:
            !  >>> for n in range(1,6):
            !  >>>     var('x')
            !  >>>     f(x) = derivative(spherical_hankel1(n, x))
            !  >>>     x=4.0
            !  >>>     value = numerical_approx(f(x))
            !  >>>     print("n = {}: {}".format(n, value))
            ground_truth_h_1(1) = cmplx(-0.247255998456561_8, +0.0483842301504241_8, kind=8)
            ground_truth_h_1(2) = cmplx(-0.0911020150693551_8, +0.223206519594221_8, kind=8)
            ground_truth_h_1(3) = cmplx(0.0470398278163199_8, +0.227771073285379_8, kind=8)
            ground_truth_h_1(4) = cmplx(0.0731275258925893_8, +0.271048718737782_8, kind=8)
            ground_truth_h_1(5) = cmplx(0.0472447560139076_8, +0.602449351963012_8, kind=8)

            h_1 = lib_math_bessel_spherical_third_kind_1_derivative_real(x, fnu, n, h_1_n)
            ground_truth_h_1_n = lib_math_bessel_spherical_third_kind_1_real(x, fnu, n)

            rv = .true.
            print *, "test_lib_math_bessel_spherical_third_kind_1_derivative_real:"
            do i=1, n
                buffer = abs(h_1(i) - ground_truth_h_1(i))
                if (buffer .gt. ground_truth_e) then
                    print *, "  ", i , "difference: ", buffer, " : FAILED"
                    rv = .false.
                else
                    print *, "  ", i, ": OK"
                end if

                buffer_y_n = abs(h_1_n(i) - ground_truth_h_1_n(i))
                if (buffer_y_n .gt. ground_truth_e) then
                    print *, "  ", i, "antiderivative: difference: ", buffer_y_n, " : FAILED"
                    rv = .false.
                end if
            end do

        end function test_lib_math_bessel_spherical_third_kind_1_derivative_real

        function test_lib_math_bessel_spherical_third_kind_1_derivative_cmplx() result (rv)
            implicit none
            ! dummy
            logical :: rv

            ! parameter
            integer, parameter :: n = 5
            double precision, parameter :: ground_truth_e = 10.0_8**(-13.0_8)

            ! auxiliary
            complex(kind=8) :: x = cmplx(4, 2, kind=8)
            complex(kind=8), dimension(n) :: h_1
            complex(kind=8), dimension(n) :: h_1_n
            integer :: fnu = 1

            complex(kind=8), dimension(n) :: ground_truth_h_1
            complex(kind=8), dimension(n) :: ground_truth_h_1_n

            integer :: i
            double precision :: buffer
            double precision :: buffer_y_n

            ! Values were generated with sageMath
            !
            ! source code:
            !  >>> f(n,x) = derivative(spherical_hankel1(n, x), x)
            !  >>>
            !  >>> x=4.0 + I*2
            !  >>> diff = 0
            !  >>> for n in range(1,6):
            !  >>>     value = numerical_approx(f(n,x))
            !  >>>     print("n = {}: {}".format(n, value))
            ground_truth_h_1(1) = cmplx(-0.0266761952934651_8, +0.0262588457410656_8, kind=8)
            ground_truth_h_1(2) = cmplx(0.0155094320476472_8, +0.0439002033342006_8, kind=8)
            ground_truth_h_1(3) = cmplx(0.0640906704642379_8, +0.0205388280247053_8, kind=8)
            ground_truth_h_1(4) = cmplx(0.107101103306046_8, -0.0491550708904054_8, kind=8)
            ground_truth_h_1(5) = cmplx(0.127613226100831_8, -0.221701023038892_8, kind=8)

            h_1 = lib_math_bessel_spherical_third_kind_1_derivative_cmplx(x, fnu, n, h_1_n)
            ground_truth_h_1_n = lib_math_bessel_spherical_third_kind_1_cmplx(x, fnu, n)

            rv = .true.
            print *, "test_lib_math_bessel_spherical_third_kind_1_derivative_cmplx:"
            do i=1, n
                buffer = abs(h_1(i) - ground_truth_h_1(i))
                if (buffer .gt. ground_truth_e) then
                    print *, "  ", i , "difference: ", buffer, " : FAILED"
                    rv = .false.
                else
                    print *, "  ", i, ": OK"
                end if

                buffer_y_n = abs(h_1_n(i) - ground_truth_h_1_n(i))
                if (buffer_y_n .gt. ground_truth_e) then
                    print *, "  ", i, "antiderivative: difference: ", buffer_y_n, " : FAILED"
                    rv = .false.
                end if
            end do

        end function test_lib_math_bessel_spherical_third_kind_1_derivative_cmplx

        function test_lib_math_bessel_spherical_third_kind_2_derivative_real() result (rv)
            implicit none
            ! dummy
            logical :: rv

            ! parameter
            integer, parameter :: n = 5
            double precision, parameter :: ground_truth_e = 10.0_8**(-13.0_8)

            ! auxiliary
            double precision :: x = 4
            complex(kind=8), dimension(n) :: h_2
            complex(kind=8), dimension(n) :: h_2_n
            integer :: fnu = 1

            complex(kind=8), dimension(n) :: ground_truth_h_2
            complex(kind=8), dimension(n) :: ground_truth_h_2_n

            integer :: i
            double precision :: buffer
            double precision :: buffer_y_n

            ! Values were generated with sageMath
            !
            ! source code:
            !  >>> for n in range(1,6):
            !  >>>     var('x')
            !  >>>     f(x) = derivative(spherical_hankel2(n, x))
            !  >>>     x=4.0
            !  >>>     value = numerical_approx(f(x))
            !  >>>     print("n = {}: {}".format(n, value))
            ground_truth_h_2(1) = cmplx(-0.247255998456561_8, -0.0483842301504241_8, kind=8)
            ground_truth_h_2(2) = cmplx(-0.0911020150693552_8, -0.223206519594221_8, kind=8)
            ground_truth_h_2(3) = cmplx(0.0470398278163199_8, -0.227771073285379_8, kind=8)
            ground_truth_h_2(4) = cmplx(0.0731275258925893_8, -0.271048718737782_8, kind=8)
            ground_truth_h_2(5) = cmplx(0.0472447560139076_8, -0.602449351963012_8, kind=8)

            h_2 = lib_math_bessel_spherical_third_kind_2_derivative_real(x, fnu, n, h_2_n)
            ground_truth_h_2_n = lib_math_bessel_spherical_third_kind_2_real(x, fnu, n)

            rv = .true.
            print *, "test_lib_math_bessel_spherical_third_kind_2_derivative_real:"
            do i=1, n
                buffer = abs(h_2(i) - ground_truth_h_2(i))
                if (buffer .gt. ground_truth_e) then
                    print *, "  ", i , "difference: ", buffer, " : FAILED"
                    rv = .false.
                else
                    print *, "  ", i, ": OK"
                end if

                buffer_y_n = abs(h_2_n(i) - ground_truth_h_2_n(i))
                if (buffer_y_n .gt. ground_truth_e) then
                    print *, "  ", i, "antiderivative: difference: ", buffer_y_n, " : FAILED"
                    rv = .false.
                end if
            end do

        end function test_lib_math_bessel_spherical_third_kind_2_derivative_real

        function test_lib_math_bessel_spherical_third_kind_2_derivative_cmplx() result (rv)
            implicit none
            ! dummy
            logical :: rv

            ! parameter
            integer, parameter :: n = 5
            double precision, parameter :: ground_truth_e = 10.0_8**(-13.0_8)

            ! auxiliary
            complex(kind=8) :: x = cmplx(4, 2, kind=8)
            complex(kind=8), dimension(n) :: h_2
            complex(kind=8), dimension(n) :: h_2_n
            integer :: fnu = 1

            complex(kind=8), dimension(n) :: ground_truth_h_2
            complex(kind=8), dimension(n) :: ground_truth_h_2_n

            integer :: i
            double precision :: buffer
            double precision :: buffer_y_n

            ! Values were generated with sageMath
            !
            ! source code:
            !  >>> for n in range(1,6):
            !  >>>     var('x')
            !  >>>     f(x) = derivative(spherical_hankel2(n, x))
            !  >>>     x=4.0 + I*2
            !  >>>     value = numerical_approx(f(x))
            !  >>>     print("n = {}: {}".format(n, value))
            ground_truth_h_2(1) = cmplx(-1.31519220175353_8, +0.211445852463207_8, kind=8)
            ground_truth_h_2(2) = cmplx(-0.501289051609952_8, -0.858648380172939_8, kind=8)
            ground_truth_h_2(3) = cmplx(0.312874508909638_8, -0.506949225159670_8, kind=8)
            ground_truth_h_2(4) = cmplx(0.286102301410601_8, +0.0850295993844548_8, kind=8)
            ground_truth_h_2(5) = cmplx(0.0102904576355813_8, +0.387741633537657_8, kind=8)

            h_2 = lib_math_bessel_spherical_third_kind_2_derivative_cmplx(x, fnu, n, h_2_n)
            ground_truth_h_2_n = lib_math_bessel_spherical_third_kind_2_cmplx(x, fnu, n)

            rv = .true.
            print *, "test_lib_math_bessel_spherical_third_kind_2_derivative_cmplx:"
            do i=1, n
                buffer = abs(h_2(i) - ground_truth_h_2(i))
                if (buffer .gt. ground_truth_e) then
                    print *, "  ", i , "difference: ", buffer, " : FAILED"
                    rv = .false.
                else
                    print *, "  ", i, ": OK"
                end if

                buffer_y_n = abs(h_2_n(i) - ground_truth_h_2_n(i))
                if (buffer_y_n .gt. ground_truth_e) then
                    print *, "  ", i, "antiderivative: difference: ", buffer_y_n, " : FAILED"
                    rv = .false.
                end if
            end do

        end function test_lib_math_bessel_spherical_third_kind_2_derivative_cmplx

        function test_lib_math_bessel_riccati_s_derivative_real() result(rv)
            implicit none
            ! dummy
            logical :: rv

            ! parameter
            integer, parameter :: n = 5

            ! auxiliary
            double precision :: x = 4
            double precision, dimension(n) :: s
            double precision, dimension(n) :: s_n

            integer :: fnu = 1

            double precision, dimension(n) :: ground_truth_s
            double precision, dimension(n) :: ground_truth_s_n
            double precision :: ground_truth_e = 10.0_8**(-13.0_8)

            integer :: i
            double precision :: buffer
            double precision :: buffer_y_n

            ! Values were generated with sageMath
            !
            ! source code:
            !  >>> var('x')
            !  >>> var('n')
            !  >>> S(n,x) = x*spherical_bessel_J(n,x)
            !  >>>
            !  >>> from IPython.display import Math
            !  >>> t = latex(S)
            !  >>> display(Math("""{}""".format(t)))
            !  >>>
            !  >>> x=4.0
            !  >>> for i in range(1,6):
            !  >>>     value = numerical_approx(S(i,x))
            !  >>>     print("n = {}: {}".format(i, value))

            ground_truth_s(1) = -0.872913244567086_8
            ground_truth_s(2) = -0.0881243745060704_8
            ground_truth_s(3) = 0.417403169220310_8
            ground_truth_s(4) = 0.417403169220310_8
            ground_truth_s(5) = 0.240744563812994_8

            s = lib_math_bessel_riccati_s_derivative_real(x, fnu, n, s_n)
            ground_truth_s_n = lib_math_bessel_riccati_s_real(x, fnu, n)

            rv = .true.
            print *, "test_lib_math_bessel_riccati_s_derivative_real:"
            do i=1, n
                buffer = abs(s(i) - ground_truth_s(i))
                if (buffer .gt. ground_truth_e) then
                    print *, "  ", i , "difference: ", buffer, " : FAILED"
                    rv = .false.
                else
                    print *, "  ", i, ": OK"
                end if

                buffer_y_n = abs(s_n(i) - ground_truth_s_n(i))
                if (buffer_y_n .gt. ground_truth_e) then
                    print *, "  ", i, "antiderivative: difference: ", buffer_y_n, " : FAILED"
                    rv = .false.
                end if
            end do
        end function test_lib_math_bessel_riccati_s_derivative_real

        function test_lib_math_bessel_riccati_s_derivative_cmplx() result(rv)
            implicit none
            ! dummy
            logical :: rv

            ! parameter
            integer, parameter :: n = 5

            ! auxiliary
            complex(kind=8) :: x = cmplx(4, 2, kind=8)
            complex(kind=8), dimension(n) :: s
            complex(kind=8), dimension(n) :: s_n

            integer :: fnu = 1

            complex(kind=8), dimension(n) :: ground_truth_s
            complex(kind=8), dimension(n) :: ground_truth_s_n
            double precision :: ground_truth_e = 10.0_8**(-13.0_8)

            integer :: i
            double precision :: buffer
            double precision :: buffer_y_n

            ! Values were generated with sageMath
            !
            ! source code:
            !  >>> var('x')
            !  >>> var('n')
            !  >>> S(n,x) = x*spherical_bessel_J(n,x)
            !  >>> S_(n,x) = derivative(S(n,x), x)
            !  >>> # WolframAlpha: S_w = x SphericalBesselJ[-1 + n, x] + SphericalBesselJ[n, x] - x SphericalBesselJ[1 + n, x])/2
            !  >>> S_w(n,x) = 1/2 * (x*spherical_bessel_J(n-1, x) + spherical_bessel_J(n,x) - x*spherical_bessel_J(n+1,x))
            !  >>>
            !  >>> from IPython.display import Math
            !  >>> t = latex(S_w)
            !  >>> display(Math("""{}""".format(t)))
            !  >>>
            !  >>> x=4.0+I*2.0
            !  >>> print("WolframAlpha:")
            !  >>> for i in range(1,6):
            !  >>>     value = numerical_approx(S_w(i,x))
            !  >>>     print("n = {}: {}".format(i, value))
            ground_truth_s(1) = cmplx(-2.88434028957354_8, -1.61856658499523_8, kind=8)
            ground_truth_s(2) = cmplx(0.446332618346892_8, -2.38825995948575_8, kind=8)
            ground_truth_s(3) = cmplx(1.66989124066778_8, -0.418303857101622_8, kind=8)
            ground_truth_s(4) = cmplx(0.873045509461675_8, +0.685823569513368_8, kind=8)
            ground_truth_s(5) = cmplx(0.0995240805571297_8, +0.579739235324844_8, kind=8)

            s = lib_math_bessel_riccati_s_derivative_cmplx(x, fnu, n, s_n)
            ground_truth_s_n = lib_math_bessel_riccati_s_cmplx(x, fnu, n)

            rv = .true.
            print *, "test_lib_math_bessel_riccati_s_derivative_cmplx:"
            do i=1, n
                buffer = abs(s(i) - ground_truth_s(i))
                if (buffer .gt. ground_truth_e) then
                    print *, "  ", i , "difference: ", buffer, " : FAILED"
                    rv = .false.
                else
                    print *, "  ", i, ": OK"
                end if

                buffer_y_n = abs(s_n(i) - ground_truth_s_n(i))
                if (buffer_y_n .gt. ground_truth_e) then
                    print *, "  ", i, "antiderivative: difference: ", buffer_y_n, " : FAILED"
                    rv = .false.
                end if
            end do
        end function test_lib_math_bessel_riccati_s_derivative_cmplx

        function test_lib_math_bessel_riccati_c_derivative_real() result(rv)
            implicit none
            ! dummy
            logical :: rv

            ! parameter
            integer, parameter :: n = 5

            ! auxiliary
            double precision :: x = 4
            double precision, dimension(n) :: c
            double precision, dimension(n) :: c_n
            integer :: fnu = 1

            double precision, dimension(n) :: ground_truth_c
            double precision, dimension(n) :: ground_truth_c_n
            double precision :: ground_truth_e = 10.0_8**(-13.0_8)

            integer :: i
            double precision :: buffer
            double precision :: buffer_y_n

            ! Values were generated with sageMath
            !
            ! source code:
            !  >>> var('x')
            !  >>> var('n')
            !  >>> C(n,x) = -x*spherical_bessel_Y(n,x)
            !  >>>
            !  >>> n=5
            !  >>>
            !  >>> from IPython.display import Math
            !  >>> t = latex(C_.simplify_full())
            !  >>> display(Math("""{}""".format(t)))
            !  >>>
            !  >>> x=4.0
            !  >>> for i in range(1,6):
            !  >>>     value = numerical_approx(C_(i,x))
            !  >>>     print("n = {}: {}".format(i, value))
            ground_truth_c(1) = -0.423590270732654_8
            ground_truth_c(2) = -0.901955185759200_8
            ground_truth_c(3) = -0.692442327238452_8
            ground_truth_c(4) = -0.692442327238452_8
            ground_truth_c(5) = -1.74699614140159_8

            c = lib_math_bessel_riccati_c_derivative_real(x, fnu, n, c_n)
            ground_truth_c_n = lib_math_bessel_riccati_c_real(x, fnu, n)

            rv = .true.
            print *, "test_lib_math_bessel_riccati_c_derivative_real:"
            do i=1, n
                buffer = abs(c(i) - ground_truth_c(i))
                if (buffer .gt. ground_truth_e) then
                    print *, "  ", i , "difference: ", buffer, " : FAILED"
                    rv = .false.
                else
                    print *, "  ", i, ": OK"
                end if

                buffer_y_n = abs(c_n(i) - ground_truth_c_n(i))
                if (buffer_y_n .gt. ground_truth_e) then
                    print *, "  ", i, "antiderivative: difference: ", buffer_y_n, " : FAILED"
                    rv = .false.
                end if
            end do
        end function test_lib_math_bessel_riccati_c_derivative_real

        function test_lib_math_bessel_riccati_c_derivative_cmplx() result(rv)
            implicit none
            ! dummy
            logical :: rv

            ! parameter
            integer, parameter :: n = 5

            ! auxiliary
            complex(kind=8) :: x = cmplx(4, 2, kind=8)
            complex(kind=8), dimension(n) :: c
            complex(kind=8), dimension(n) :: c_n
            integer :: fnu = 1

            complex(kind=8), dimension(n) :: ground_truth_c
            complex(kind=8), dimension(n) :: ground_truth_c_n
            double precision :: ground_truth_e = 10.0_8**(-13.0_8)

            integer :: i
            double precision :: buffer
            double precision :: buffer_y_n

            ! Values were generated with sageMath
            !
            ! source code:
            !  >>> var('x')
            !  >>> var('n')
            !  >>> C(n,x) = -x*spherical_bessel_Y(n,x)
            !  >>> C_(n,x) = derivative(C(n,x), x)
            !  >>>
            !  >>> n=5
            !  >>>
            !  >>> from IPython.display import Math
            !  >>> t = latex(C_.simplify_full())
            !  >>> display(Math("""{}""".format(t)))
            !  >>>
            !  >>> x=4.0+I*2
            !  >>> for i in range(1,6):
            !  >>>     value = numerical_approx(C_(i,x))
            !  >>>     print("n = {}: {}".format(i, value))
            ground_truth_c(1) = cmplx(-1.68863860346658_8, +2.75351801321721_8, kind=8)
            ground_truth_c(2) = cmplx(-2.56945857676562_8, -0.437900159661906_8, kind=8)
            ground_truth_c(3) = cmplx(-0.567732855301268_8, -1.46152180407321_8, kind=8)
            ground_truth_c(4) = cmplx(0.723250390620301_8, -0.432843274896996_8, kind=8)
            ground_truth_c(5) = cmplx(1.17156387153814_8, +0.656016412602718_8, kind=8)

            c = lib_math_bessel_riccati_c_derivative_cmplx(x, fnu, n, c_n)
            ground_truth_c_n = lib_math_bessel_riccati_c_cmplx(x, fnu, n)

            rv = .true.
            print *, "test_lib_math_bessel_riccati_c_derivative_cmplx:"
            do i=1, n
                buffer = abs(c(i) - ground_truth_c(i))
                if (buffer .gt. ground_truth_e) then
                    print *, "  ", i , "difference: ", buffer, " : FAILED"
                    rv = .false.
                else
                    print *, "  ", i, ": OK"
                end if

                buffer_y_n = abs(c_n(i) - ground_truth_c_n(i))
                if (buffer_y_n .gt. ground_truth_e) then
                    print *, "  ", i, "antiderivative: difference: ", buffer_y_n, " : FAILED"
                    rv = .false.
                end if
            end do
        end function test_lib_math_bessel_riccati_c_derivative_cmplx

        function test_lib_math_bessel_riccati_xi_derivative_real() result(rv)
            implicit none
            ! dummy
            logical :: rv

            ! parameter
            integer, parameter :: n = 5

            ! auxiliary
            double precision :: x = 4
            complex(kind=8), dimension(n) :: xi
            complex(kind=8), dimension(n) :: xi_n
            integer :: fnu = 1

            complex(kind=8), dimension(n) :: ground_truth_xi
            complex(kind=8), dimension(n) :: ground_truth_xi_n
            double precision :: ground_truth_e = 10.0_8**(-13.0_8)

            integer :: i
            double precision :: buffer
            double precision :: buffer_y_n

            ! Values were generated with sageMath
            !
            ! source code:
            !  >>> var('x')
            !  >>> var('n')
            !  >>> Xi(n,x) = x*spherical_hankel1(n,x)
            !  >>> Xi_(n,x) = derivative(Xi(n,x), x)
            !  >>>
            !  >>> from IPython.display import Math
            !  >>> t = latex(Xi_.simplify_full())
            !  >>> display(Math("""{}""".format(t)))
            !  >>>
            !  >>> x=4.0
            !  >>> for i in range(1,6):
            !  >>>     value = numerical_approx(Xi_(i,x))
            !  >>>     print("n = {}: {}".format(i, value))
            !  >>>
            ground_truth_xi(1) = cmplx(-0.872913244567086_8, +0.423590270732654_8, kind=8)
            ground_truth_xi(2) = cmplx(-0.0881243745060705_8, +0.901955185759200_8, kind=8)
            ground_truth_xi(3) = cmplx(0.417403169220310_8, +0.692442327238452_8, kind=8)
            ground_truth_xi(4) = cmplx(0.417403169220310_8, +0.692442327238452_8, kind=8)
            ground_truth_xi(5) = cmplx(0.240744563812994_8, +1.74699614140159_8, kind=8)

            xi = lib_math_bessel_riccati_xi_derivative_real(x, fnu, n, xi_n)
            ground_truth_xi_n = lib_math_bessel_riccati_xi_real(x, fnu, n)

            rv = .true.
            print *, "test_lib_math_bessel_riccati_xi_derivative_real:"
            do i=1, n
                buffer = abs(xi(i) - ground_truth_xi(i))
                if (buffer .gt. ground_truth_e) then
                    print *, "  ", i , "difference: ", buffer, " : FAILED"
                    rv = .false.
                else
                    print *, "  ", i, ": OK"
                end if

                buffer_y_n = abs(xi_n(i) - ground_truth_xi_n(i))
                if (buffer_y_n .gt. ground_truth_e) then
                    print *, "  ", i, "antiderivative: difference: ", buffer_y_n, " : FAILED"
                    rv = .false.
                end if
            end do
        end function test_lib_math_bessel_riccati_xi_derivative_real

        function test_lib_math_bessel_riccati_xi_derivative_cmplx() result(rv)
            implicit none
            ! dummy
            logical :: rv

            ! parameter
            integer, parameter :: n = 5

            ! auxiliary
            complex(kind=8) :: x = cmplx(4, 2)
            complex(kind=8), dimension(n) :: xi
            complex(kind=8), dimension(n) :: xi_n
            integer :: fnu = 1

            complex(kind=8), dimension(n) :: ground_truth_xi
            complex(kind=8), dimension(n) :: ground_truth_xi_n
            double precision :: ground_truth_e = 10.0_8**(-13.0_8)

            integer :: i
            double precision :: buffer
            double precision :: buffer_y_n

            ! Values were generated with sageMath
            !
            ! source code:
            !  >>> var('x')
            !  >>> var('n')
            !  >>> Xi(n,x) = x*spherical_hankel1(n,x)
            !  >>> Xi_(n,x) = derivative(Xi(n,x), x)
            !  >>>
            !  >>> from IPython.display import Math
            !  >>> t = latex(Xi_.simplify_full())
            !  >>> display(Math("""{}""".format(t)))
            !  >>>
            !  >>> x=4.0+2.0*I
            !  >>> for i in range(1,6):
            !  >>>     value = numerical_approx(Xi_(i,x))
            !  >>>     print("n = {}: {}".format(i, value))
            !  >>>
            ground_truth_xi(1) = cmplx(-0.130822276356333_8, +0.0700720184713570_8, kind=8)
            ground_truth_xi(2) = cmplx(0.00843245868498678_8, +0.181198617279870_8, kind=8)
            ground_truth_xi(3) = cmplx(0.208369436594568_8, +0.149428998199646_8, kind=8)
            ground_truth_xi(4) = cmplx(0.440202234564679_8, -0.0374268211069332_8, kind=8)
            ground_truth_xi(5) = cmplx(0.755540493159848_8, -0.591824636213295_8, kind=8)

            xi = lib_math_bessel_riccati_xi_derivative_cmplx(x, fnu, n, xi_n)
            ground_truth_xi_n = lib_math_bessel_riccati_xi_cmplx(x, fnu, n)

            rv = .true.
            print *, "test_lib_math_bessel_riccati_xi_derivative_cmplx:"
            do i=1, n
                buffer = abs(xi(i) - ground_truth_xi(i))
                if (buffer .gt. ground_truth_e) then
                    print *, "  ", i , "difference: ", buffer, " : FAILED"
                    rv = .false.
                else
                    print *, "  ", i, ": OK"
                end if

                buffer_y_n = abs(xi_n(i) - ground_truth_xi_n(i))
                if (buffer_y_n .gt. ground_truth_e) then
                    print *, "  ", i, "antiderivative: difference: ", buffer_y_n, " : FAILED"
                    rv = .false.
                end if
            end do
        end function test_lib_math_bessel_riccati_xi_derivative_cmplx

        function test_lib_math_bessel_riccati_zeta_derivative_real() result(rv)
            implicit none
            ! dummy
            logical :: rv

            ! parameter
            integer, parameter :: n = 5

            ! auxiliary
            double precision :: x = 4
            complex(kind=8), dimension(n) :: zeta
            complex(kind=8), dimension(n) :: zeta_n
            integer :: fnu = 1

            complex(kind=8), dimension(n) :: ground_truth_zeta
            complex(kind=8), dimension(n) :: ground_truth_zeta_n
            double precision :: ground_truth_e = 10.0_8**(-13.0_8)

            integer :: i
            double precision :: buffer
            double precision :: buffer_y_n

            ! Values were generated with sageMath
            !
            ! source code:
            !  >>> var('x')
            !  >>> var('n')
            !  >>> Zeta(n,x) = x*spherical_hankel2(n,x)
            !  >>> Zeta_(n,x) = derivative(Zeta(n,x), x)
            !  >>>
            !  >>> from IPython.display import Math
            !  >>> t = latex(Zeta_.simplify_full())
            !  >>> display(Math("""{}""".format(t)))
            !  >>>
            !  >>> x=4.0
            !  >>> for i in range(1,6):
            !  >>>     value = numerical_approx(Zeta_(i,x))
            !  >>>     print("n = {}: {}".format(i, value))
            !  >>>
            ground_truth_zeta(1) = cmplx(-0.872913244567086_8, -0.423590270732654_8, kind=8)
            ground_truth_zeta(2) = cmplx(-0.0881243745060704_8, -0.901955185759200_8, kind=8)
            ground_truth_zeta(3) = cmplx(0.417403169220310_8, -0.692442327238452_8, kind=8)
            ground_truth_zeta(4) = cmplx(0.417403169220310_8, -0.692442327238452_8, kind=8)
            ground_truth_zeta(5) = cmplx(0.240744563812994_8, -1.74699614140159_8, kind=8)

            zeta = lib_math_bessel_riccati_zeta_derivative_real(x, fnu, n, zeta_n)
            ground_truth_zeta_n = lib_math_bessel_riccati_zeta_real(x, fnu, n)

            rv = .true.
            print *, "test_lib_math_bessel_riccati_zeta_derivative_real:"
            do i=1, n
                buffer = abs(zeta(i) - ground_truth_zeta(i))
                if (buffer .gt. ground_truth_e) then
                    print *, "  ", i , "difference: ", buffer, " : FAILED"
                    rv = .false.
                else
                    print *, "  ", i, ": OK"
                end if

                buffer_y_n = abs(zeta_n(i) - ground_truth_zeta_n(i))
                if (buffer_y_n .gt. ground_truth_e) then
                    print *, "  ", i, "antiderivative: difference: ", buffer_y_n, " : FAILED"
                    rv = .false.
                end if
            end do
        end function test_lib_math_bessel_riccati_zeta_derivative_real

        function test_lib_math_bessel_riccati_zeta_derivative_cmplx() result(rv)
            implicit none
            ! dummy
            logical :: rv

            ! parameter
            integer, parameter :: n = 5

            ! auxiliary
            complex(kind=8) :: x = cmplx(4,2)
            complex(kind=8), dimension(n) :: zeta
            complex(kind=8), dimension(n) :: zeta_n
            integer :: fnu = 1

            complex(kind=8), dimension(n) :: ground_truth_zeta
            complex(kind=8), dimension(n) :: ground_truth_zeta_n
            double precision :: ground_truth_e = 10.0_8**(-13.0_8)

            integer :: i
            double precision :: buffer
            double precision :: buffer_y_n

            ! Values were generated with sageMath
            !
            ! source code:
            !  >>> var('x')
            !  >>> var('n')
            !  >>> Zeta(n,x) = x*spherical_hankel2(n,x)
            !  >>> Zeta_(n,x) = derivative(Zeta(n,x), x)
            !  >>>
            !  >>> from IPython.display import Math
            !  >>> t = latex(Zeta_.simplify_full())
            !  >>> display(Math("""{}""".format(t)))
            !  >>>
            !  >>> x=4.0+I*2.0
            !  >>> for i in range(1,6):
            !  >>>     value = numerical_approx(Zeta_(i,x))
            !  >>>     print("n = {}: {}".format(i, value))
            !  >>>
            ground_truth_zeta(1) = cmplx(-5.63785830279076_8, -3.30720518846181_8, kind=8)
            ground_truth_zeta(2) = cmplx(0.884232778008798_8, -4.95771853625136_8, kind=8)
            ground_truth_zeta(3) = cmplx(3.13141304474100_8, -0.986036712402890_8, kind=8)
            ground_truth_zeta(4) = cmplx(1.30588878435867_8, +1.40907396013367_8, kind=8)
            ground_truth_zeta(5) = cmplx(-0.556492332045588_8, +1.75130310686298_8, kind=8)

            zeta = lib_math_bessel_riccati_zeta_derivative_cmplx(x, fnu, n, zeta_n)
            ground_truth_zeta_n = lib_math_bessel_riccati_zeta_cmplx(x, fnu, n)

            rv = .true.
            print *, "test_lib_math_bessel_riccati_zeta_derivative_cmplx:"
            do i=1, n
                buffer = abs(zeta(i) - ground_truth_zeta(i))
                if (buffer .gt. ground_truth_e) then
                    print *, "  ", i , "difference: ", buffer, " : FAILED"
                    rv = .false.
                else
                    print *, "  ", i, ": OK"
                end if

                buffer_y_n = abs(zeta_n(i) - ground_truth_zeta_n(i))
                if (buffer_y_n .gt. ground_truth_e) then
                    print *, "  ", i, "antiderivative: difference: ", buffer_y_n, " : FAILED"
                    rv = .false.
                end if
            end do
        end function test_lib_math_bessel_riccati_zeta_derivative_cmplx

    end function lib_math_bessel_test_functions

end module lib_math_bessel
