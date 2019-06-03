module lib_math_bessel
    implicit none

    private

    interface lib_math_bessel_spherical_first_kind
        module procedure lib_math_bessel_spherical_first_kind_real
        !module procedure lib_math_bessel_spherical_first_kind_cmplx    ! <--- todo: implement the complex version
    end interface

    interface lib_math_bessel_spherical_second_kind
        module procedure lib_math_bessel_spherical_second_kind_real
        !module procedure lib_math_bessel_spherical_second_kind_cmplx   ! <--- todo: implement the complex version
    end interface

    interface lib_math_bessel_spherical_third_kind_1
        module procedure lib_math_bessel_spherical_third_kind_1_real
    end interface

    interface lib_math_bessel_spherical_third_kind_2
        module procedure lib_math_bessel_spherical_third_kind_2_real
    end interface

    double precision, parameter :: PI=4.D0*atan(1.D0)  ! maximum precision

    contains

    ! calculates the the spherical Bessel function of the first kind
    !
    ! symbol: j_n
    !
    ! Argument
    ! ----
    !   n: integer
    !       order
    !   x: double precision
    !
    ! Returns
    ! ----
    !   rv: double precision
    !
    ! LaTeX: $$ j_{n}(\rho)=\sqrt{\frac{\pi}{2 \rho}} J_{n+1 / 2}(\rho) $$
    !
    ! Reference: Absorption and Scattering of Light by Small Particles Absorption and Scattering of Light by Small Particles (4.9)
    !
    function lib_math_bessel_spherical_first_kind_real(n, x) result (rv)
        implicit none
        ! dummy
        integer, intent(in) :: n
        double precision, intent(in) :: x
        double precision :: rv

!        rv = sqrt(PI/(2.D0*x)) * bessel_jn(real(n)+0.5,x)

    end function lib_math_bessel_spherical_first_kind_real


    ! calculates the the spherical Bessel function of the second kind
    !
    ! symbol: y_n
    !
    ! Argument
    ! ----
    !   n: integer
    !       order
    !   x: double precision
    !
    ! Returns
    ! ----
    !   rv: double precision
    !
    ! LaTeX: $$ y_{n}(\rho)=\sqrt{\frac{\pi}{2 \rho}} Y_{n+1 / 2}(\rho) $$
    !
    ! Reference: Absorption and Scattering of Light by Small Particles Absorption and Scattering of Light by Small Particles (4.10)
    !
    function lib_math_bessel_spherical_second_kind_real(n, x) result (rv)
        implicit none
        ! dummy
        integer, intent(in) :: n
        double precision, intent(in) :: x
        double precision :: rv

!        rv = sqrt(PI/(2.D0*x)) * bessel_yn(real(n)+0.5,x)

    end function lib_math_bessel_spherical_second_kind_real

    ! calculates the the spherical Bessel function of the third kind
    !
    ! symbol: h_n^(1)
    !
    ! Argument
    ! ----
    !   n: integer
    !       order
    !   x: double precision
    !
    ! Returns
    ! ----
    !   rv: double precision
    !
    !
    ! LaTeX: $$ h_{n}^{(1)}(\rho)=j_{n}(\rho)+i y_{n}(\rho) $$
    !
    ! Reference: Absorption and Scattering of Light by Small Particles Absorption and Scattering of Light by Small Particles (4.13)
    !
    function lib_math_bessel_spherical_third_kind_1_real(n, x) result (rv)
        implicit none
        ! dummy
        integer, intent(in) :: n
        double precision, intent(in) :: x
        complex(kind=8) :: rv

!        rv = cmplx(lib_math_bessel_spherical_first_kind_real, lib_math_bessel_spherical_second_kind_real)

    end function lib_math_bessel_spherical_third_kind_1_real

    ! calculates the the spherical Bessel function of the third kind
    !
    ! symbol: h_n^(2)
    !
    ! Argument
    ! ----
    !   n: integer
    !       order
    !   x: double precision
    !
    ! Returns
    ! ----
    !   rv: double precision
    !
    !
    ! LaTeX: $$ h_{n}^{(2)}(\rho)=j_{n}(\rho) - i y_{n}(\rho) $$
    !
    ! Reference: Absorption and Scattering of Light by Small Particles Absorption and Scattering of Light by Small Particles (4.14)
    !
    function lib_math_bessel_spherical_third_kind_2_real(n, x) result (rv)
        implicit none
        ! dummy
        integer, intent(in) :: n
        double precision, intent(in) :: x
        complex(kind=8) :: rv

!        rv = cmplx(lib_math_bessel_spherical_first_kind_real, -lib_math_bessel_spherical_second_kind_real)

    end function lib_math_bessel_spherical_third_kind_2_real

    ! Calculates the derivative of the spherical Bessel function of the first kind
    !
    ! symbol: j'_n
    !
    ! Argument
    ! ----
    !   n: integer
    !       order
    !   x: double precision
    !
    ! Returns
    ! ----
    !   rv: double precision
    !
    !
    ! LaTeX: $$ j_{n}^{\prime}=\left(\frac{n}{x}\right) j_{n}(x)-j_{n+1}(x) $$
    !
    ! Reference: https://www.boost.org/doc/libs/1_57_0/libs/math/doc/html/math_toolkit/bessel/bessel_derivatives.html
    !
    function lib_math_bessel_spherical_first_kind_real_derivative(n, x) result (rv)
        implicit none
        ! dummy
        integer, intent(in) :: n
        double precision, intent(in) :: x
        double precision :: rv

        rv = real(n, 8) / x * lib_math_bessel_spherical_first_kind_real(n,x) - lib_math_bessel_spherical_first_kind_real(n+1, x)

    end function lib_math_bessel_spherical_first_kind_real_derivative

    ! Calculates the derivative of the spherical Bessel function of the second kind
    !
    ! symbol: y'_n
    !
    ! Argument
    ! ----
    !   n: integer
    !       order
    !   x: double precision
    !
    ! Returns
    ! ----
    !   rv: double precision
    !
    !
    ! LaTeX: $$ y_{n}^{\prime}=\left(\frac{n}{x}\right) y_{n}(x)-y_{n+1}(x) $$
    !
    ! Reference: https://www.boost.org/doc/libs/1_57_0/libs/math/doc/html/math_toolkit/bessel/bessel_derivatives.html
    !
    function lib_math_bessel_spherical_second_kind_real_derivative(n, x) result (rv)
        implicit none
        ! dummy
        integer, intent(in) :: n
        double precision, intent(in) :: x
        double precision :: rv

        rv = real(n, 8) / x * lib_math_bessel_spherical_second_kind_real(n,x) - lib_math_bessel_spherical_second_kind_real(n+1, x)

    end function lib_math_bessel_spherical_second_kind_real_derivative

end module lib_math_bessel
