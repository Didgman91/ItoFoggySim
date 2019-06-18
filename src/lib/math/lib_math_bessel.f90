#define _DEBUG_

module lib_math_bessel
    use libmath
    implicit none

    private

    ! --- interface ---
    interface lib_math_bessel_spherical_first_kind
        module procedure lib_math_bessel_spherical_first_kind_real
        module procedure lib_math_bessel_spherical_first_kind_cmplx
    end interface

    interface lib_math_bessel_spherical_second_kind
        module procedure lib_math_bessel_spherical_second_kind_real
        module procedure lib_math_bessel_spherical_second_kind_cmplx
    end interface

    interface lib_math_hankel_spherical_1
        module procedure lib_math_bessel_spherical_third_kind_1_real
        module procedure lib_math_bessel_spherical_third_kind_1_cmplx
    end interface

    interface lib_math_hankel_spherical_2
        module procedure lib_math_bessel_spherical_third_kind_2_real
    end interface

    interface lib_math_bessel_spherical_first_kind_derivative
        module procedure lib_math_bessel_spherical_first_kind_derivative_real
    end interface

    interface lib_math_bessel_spherical_second_kind_derivative
        module procedure lib_math_bessel_spherical_second_kind_derivative_real
    end interface

    interface lib_math_hankel_1_derivative
        module procedure lib_math_bessel_spherical_third_kind_1_derivative_real
    end interface

    interface lib_math_hankel_2_derivative
        module procedure lib_math_bessel_spherical_third_kind_2_derivative_real
    end interface

    ! --- public functions ---
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
    !       ORDER OF INITIAL J FUNCTION, FNU.GE.0
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
    !       ORDER OF INITIAL J FUNCTION, FNU.GE.0
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
    !       ORDER OF INITIAL J FUNCTION, FNU.GE.0
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
    !       ORDER OF INITIAL J FUNCTION, FNU.GE.0
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
    !   fnu: integer
    !       ORDER OF INITIAL J FUNCTION, FNU.GE.0
    !   n: integer
    !       NUMBER OF MEMBERS OF THE SEQUENCE, N.GE.1
    !
    ! Returns
    ! ----
    !   rv: array<double precision>
    !
    !
    ! LaTeX: $$ h_{n}^{(1)}(\rho)=j_{n}(\rho)+i y_{n}(\rho) $$
    !
    ! Reference: Absorption and Scattering of Light by Small Particles Absorption and Scattering of Light by Small Particles (4.13)
    !
    function lib_math_bessel_spherical_third_kind_1_real(x, fnu, n) result (rv)
        implicit none
        ! dummy
        double precision, intent(in) :: x
        integer(kind=4), intent(in) :: fnu
        integer(kind=4), intent(in) :: n
        complex(kind=8), dimension(n) :: rv

        ! auxiliary
!        integer(kind=4) order
!        integer(kind=4) nm
!        complex(kind=8) chf1(0:fnu+n-1)
!        complex(kind=8) chd1(0:fnu+n-1)
!        complex(kind=8) chf2(0:fnu+n-1)
!        complex(kind=8) chd2(0:fnu+n-1)

!        order = fnu+n-1

        rv = cmplx(lib_math_bessel_spherical_first_kind_real(x, fnu, n), &
                   lib_math_bessel_spherical_second_kind_real(x, fnu, n), &
                   kind=8)

!        call CH12N( order, cmplx(x, 0, kind=8), nm, chf1, chd1, chf2, chd2)
!
!        if (order .ne. nm) then
!            print *, "lib_math_bessel_spherical_third_kind_1_real: ERROR"
!            print *, "  calculated highest order / requested: ", nm, " / ", order
!        end if
!
!        rv = chf1(fnu:order)

    end function lib_math_bessel_spherical_third_kind_1_real

    ! calculates the the spherical Bessel function of the third kind
    !
    ! symbol: h_n^(1)
    !
    ! Argument
    ! ----
    !   x: double precision
    !   fnu: integer
    !       ORDER OF INITIAL J FUNCTION, FNU.GE.0
    !   n: integer
    !       NUMBER OF MEMBERS OF THE SEQUENCE, N.GE.1
    !
    ! Returns
    ! ----
    !   rv: array<double precision>
    !
    !
    ! LaTeX: $$ h_{n}^{(1)}(\rho)=j_{n}(\rho)+i y_{n}(\rho) $$
    !
    ! Reference: http://mathworld.wolfram.com/SphericalHankelFunctionoftheFirstKind.html
    !
    function lib_math_bessel_spherical_third_kind_1_cmplx(z, fnu, n) result (rv)
        implicit none
        ! dummy
        complex(kind=8), intent(in) :: z
        integer(kind=4), intent(in) :: fnu
        integer(kind=4), intent(in) :: n
        complex(kind=8), dimension(n) :: rv

        rv = lib_math_bessel_spherical_first_kind_cmplx(z, fnu, n) &
             + cmplx(0,1, kind=8) * lib_math_bessel_spherical_second_kind_cmplx(z, fnu, n)

    end function lib_math_bessel_spherical_third_kind_1_cmplx

    ! calculates the the spherical Bessel function of the third kind
    !
    ! symbol: h_n^(2)
    !
    ! Argument
    ! ----
    !   x: double precision
    !   fnu: integer
    !       ORDER OF INITIAL J FUNCTION, FNU.GE.0
    !   n: integer
    !       NUMBER OF MEMBERS OF THE SEQUENCE, N.GE.1
    !
    ! Returns
    ! ----
    !   rv: array<double precision>
    !
    !
    ! LaTeX: $$ h_{n}^{(2)}(\rho)=j_{n}(\rho) - i y_{n}(\rho) $$
    !
    ! Reference: Absorption and Scattering of Light by Small Particles Absorption and Scattering of Light by Small Particles (4.14)
    !
    function lib_math_bessel_spherical_third_kind_2_real(x, fnu, n) result (rv)
        implicit none
        ! dummy
        double precision, intent(in) :: x
        integer(kind=4), intent(in) :: fnu
        integer(kind=4), intent(in) :: n
        complex(kind=8), dimension(n) :: rv

        ! auxiliary
!        integer(kind=4) order
!        integer(kind=4) nm
!        complex(kind=8) chf1(0:fnu+n-1)
!        complex(kind=8) chd1(0:fnu+n-1)
!        complex(kind=8) chf2(0:fnu+n-1)
!        complex(kind=8) chd2(0:fnu+n-1)
!
!        order = fnu+n-1

        rv = cmplx(lib_math_bessel_spherical_first_kind_real(x, fnu, n), &
                   - lib_math_bessel_spherical_second_kind_real(x, fnu, n), &
                   kind=8)

!        call CH12N( order, cmplx(x, 0, kind=8), nm, chf1, chd1, chf2, chd2)
!
!        if (order .ne. nm) then
!            print *, "lib_math_bessel_spherical_third_kind_2_real: ERROR"
!            print *, "  calculated highest order / requested: ", nm, " / ", order
!        end if
!
!        rv = chf2(fnu:order)

    end function lib_math_bessel_spherical_third_kind_2_real

    ! calculates the the spherical Bessel function of the third kind
    !
    ! symbol: h_n^(2)
    !
    ! Argument
    ! ----
    !   x: double precision
    !   fnu: integer
    !       ORDER OF INITIAL J FUNCTION, FNU.GE.0
    !   n: integer
    !       NUMBER OF MEMBERS OF THE SEQUENCE, N.GE.1
    !
    ! Returns
    ! ----
    !   rv: array<double precision>
    !
    !
    ! LaTeX: $$ h_{n}^{(2)}(\rho)=j_{n}(\rho) - i y_{n}(\rho) $$
    !
    ! Reference:
    !
    function lib_math_bessel_spherical_third_kind_2_cmplx(z, fnu, n) result (rv)
        implicit none
        ! dummy
        complex(kind=8), intent(in) :: z
        integer(kind=4), intent(in) :: fnu
        integer(kind=4), intent(in) :: n
        complex(kind=8), dimension(n) :: rv

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
    !       ORDER OF INITIAL J FUNCTION, FNU.GE.0
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
    !       ORDER OF INITIAL J FUNCTION, FNU.GE.0
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
    !       ORDER OF INITIAL J FUNCTION, FNU.GE.0
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
    !       ORDER OF INITIAL J FUNCTION, FNU.GE.0
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
    !       ORDER OF INITIAL J FUNCTION, FNU.GE.0
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
    ! symbol: Zeta
    !
    ! Argument
    ! ----
    !   x: double precision
    !   fnu: integer
    !       ORDER OF INITIAL J FUNCTION, FNU.GE.0
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
    function lib_math_bessel_riccati_zetta_real(x, fnu, n) result (rv)
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
            print *, "lib_math_bessel_riccati_s_real: WARNING"
            print *, "  number of components of Y set to zero due to underflow: ", nz
        end if
#endif

        call DBESY (X, FNU+0.5D0, N, y)

        rv(:) = cmplx(rv_1 * j(:), -rv_1 * y(:), kind=8)

    end function lib_math_bessel_riccati_zetta_real

    ! Calculates the derivative of the spherical Bessel function of the first kind
    !
    ! symbol: j'_n
    !
    ! Argument
    ! ----
    !   x: double precision
    !   fnu: integer
    !       ORDER OF INITIAL J FUNCTION, FNU.GE.0
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
    ! LaTeX:
    !
    ! Reference:
    !
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
    ! symbol: j'_n
    !
    ! Argument
    ! ----
    !   z: complex
    !   fnu: integer
    !       ORDER OF INITIAL J FUNCTION, FNU.GE.0
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
    ! LaTeX:
    !
    ! Reference:
    !
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

    ! Calculates the derivative of the spherical Bessel function of the second kind
    !
    ! symbol: y'_n
    !
    ! Argument
    ! ----
    !   x: double precision
    !   fnu: integer
    !       ORDER OF INITIAL J FUNCTION, FNU.GE.0
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
    !       ORDER OF INITIAL J FUNCTION, FNU.GE.0
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
    !       ORDER OF INITIAL J FUNCTION, FNU.GE.1
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
        complex(kind=8), dimension(n+2), intent(inout) :: h_1_n
        complex(kind=8), dimension(n) :: rv

        ! auxiliary
        integer :: i

        h_1_n = lib_math_bessel_spherical_third_kind_1_real(x, fnu-1, n+2)

        do i=1, n
             rv(i) = 0.5D0 * (h_1_n(i) - (h_1_n(i+1) + x * h_1_n(i+2))/x)
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
    !       ORDER OF INITIAL J FUNCTION, FNU.GE.1
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
        complex(kind=8), dimension(n+2), intent(inout) :: h_1_n
        complex(kind=8), dimension(n) :: rv

        ! auxiliary
        integer :: i

        h_1_n = lib_math_bessel_spherical_third_kind_1_cmplx(z, fnu-1, n+2)

        do i=1, n
             rv(i) = 0.5D0 * (h_1_n(i) - (h_1_n(i+1) + z * h_1_n(i+2))/z)
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
    !       ORDER OF INITIAL J FUNCTION, FNU.GE.1
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
        complex(kind=8), dimension(n+2), intent(inout) :: h_2_n
        complex(kind=8), dimension(n) :: rv

        ! auxiliary
        integer :: i

        h_2_n = lib_math_bessel_spherical_third_kind_2_real(x, fnu-1, n+2)

        do i=1, n
             rv(i) = 0.5D0 * (h_2_n(i) - (h_2_n(i+1) + x * h_2_n(i+2))/x)
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
    !       ORDER OF INITIAL J FUNCTION, FNU.GE.1
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
        complex(kind=8), dimension(n+2), intent(inout) :: h_2_n
        complex(kind=8), dimension(n) :: rv

        ! auxiliary
        integer :: i

        h_2_n = lib_math_bessel_spherical_third_kind_2_cmplx(z, fnu-1, n+2)

        do i=1, n
             rv(i) = 0.5D0 * (h_2_n(i) - (h_2_n(i+1) + z * h_2_n(i+2))/z)
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
    !       ORDER OF INITIAL J FUNCTION, FNU.GE.1
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
    !       ORDER OF INITIAL J FUNCTION, FNU.GE.1
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
    function lib_math_bessel_riccati_s_derivative_cmplx(z, fnu, n, r_n) result (rv)
        implicit none
        ! dummy
        complex(kind=8), intent(in) :: z
        integer(kind=4), intent(in) :: fnu
        integer(kind=4), intent(in) :: n
        complex(kind=8), intent(inout), dimension(n) :: r_n

        complex(kind=8), dimension(n) :: rv

        ! auxiliary
        integer(kind=4) i
        integer(kind=4) order
        complex(kind=8), dimension(fnu-1:fnu+n) :: jf

        order = fnu+n-1

        jf = lib_math_bessel_spherical_first_kind_cmplx(z, fnu-1, n+1)

        r_n = z * jf(fnu:order)

        do i=2, n+1
            rv(i-1) = r_n(i-1) - i*jf(i)
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
    !       ORDER OF INITIAL J FUNCTION, FNU.GE.1
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
        if (.not. test_lib_math_bessel_riccati_zetta_real()) then
            rv = rv + 1
        end if

        ! test: derivatives of the Bessel functions
        if (.not. test_lib_math_bessel_spherical_first_kind_derivative_real()) then
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
!        if (.not. test_lib_math_bessel_riccati_xi_derivative_real()) then
!            rv = rv + 1
!        end if
!        if (.not. test_lib_math_bessel_riccati_zetta_derivative_real()) then
!            rv = rv + 1
!        end if


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

        ! test: derivatives of Bessel functions with complex argument
        if (.not. test_lib_math_bessel_spherical_first_kind_derivative_cmplx()) then
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
            integer, parameter :: n = 5

            ! auxiliary
            double precision :: x = 5
            complex(kind=8), dimension(n) :: y
            integer :: fnu = 0

            complex(kind=8), dimension(n) :: ground_truth_y
            double precision :: ground_truth_e = 10.0_8**(-13.0_8)

            integer :: i
            double precision :: buffer

            ! Values were generated with sageMath
            !
            ! source code:
            !  >>> x=5.0 +I*0.0
            !  >>> for n in range(0,5):
            !  >>>     value = numerical_approx(spherical_hankel1(n, x))
            !  >>>     print("n = {}: {}".format(n, value))
            ground_truth_y(1) = cmplx(-0.191784854932628_8, -0.0567324370926452_8, kind=8)
            ground_truth_y(2) = cmplx(-0.0950894080791708_8, +0.180438367514099_8, kind=8)
            ground_truth_y(3) = cmplx(0.134731210085125_8, +0.164995457601104_8, kind=8)
            ground_truth_y(4) = cmplx(0.229820618164296_8, -0.0154429099129942_8, kind=8)
            ground_truth_y(5) = cmplx(0.187017655344889_8, -0.186615531479296_8, kind=8)

            y = lib_math_bessel_spherical_third_kind_1_real(x, fnu, n)

            rv = .true.
            print *, "test_lib_math_bessel_spherical_third_kind_1_real:"
            do i=1, n
                buffer = abs(y(i) - ground_truth_y(i))
                if (buffer .gt. ground_truth_e) then
                    print *, "  ", i , "difference: ", buffer, " : FAILED"
                    rv = .false.
                else
                    print *, "  ", i, ": OK"
                end if
            end do

        end function test_lib_math_bessel_spherical_third_kind_1_real

        function test_lib_math_bessel_spherical_third_kind_1_cmplx() result (rv)
            implicit none
            ! dummy
            logical :: rv

            ! parameter
            integer, parameter :: n = 5
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
            ground_truth_y(1) = cmplx(-0.0250227739998745_8, +0.00233120915731335_8, kind = 8)
            ground_truth_y(2) = cmplx(-0.00182228917664341_8, +0.0271504151649199_8, kind = 8)
            ground_truth_y(3) = cmplx(0.0296975379081458_8, +0.0120891343783301_8, kind = 8)
            ground_truth_y(4) = cmplx(0.0315922819865381_8, -0.0269692779105477_8, kind = 8)
            ground_truth_y(5) = cmplx(-0.00458857312258838_8, -0.0598897093673198_8, kind = 8)

            y = lib_math_bessel_spherical_third_kind_1_cmplx(x, fnu, n)

            rv = .true.
            print *, "test_lib_math_bessel_spherical_third_kind_1_cmplx:"
            do i=1, n
                buffer = abs(y(i) - ground_truth_y(i))
                if (buffer .gt. ground_truth_e) then
                    print *, "  ", i , "difference: ", buffer, " : FAILED"
                    rv = .false.
                else
                    print *, "  ", i, ": OK"
                end if
            end do

        end function test_lib_math_bessel_spherical_third_kind_1_cmplx

        function test_lib_math_bessel_spherical_third_kind_2_real() result (rv)
            implicit none
            ! dummy
            logical :: rv

            ! parameter
            integer, parameter :: n = 5

            ! auxiliary
            double precision :: x = 5
            complex(kind=8), dimension(n) :: y
            integer :: fnu = 0

            complex(kind=8), dimension(n) :: ground_truth_y
            double precision :: ground_truth_e = 10.0_8**(-13.0_8)

            integer :: i
            double precision :: buffer

            ! Values were generated with sageMath
            !
            ! source code:
            !  >>> x=5.0 +I*0.0
            !  >>> for n in range(0,5):
            !  >>>     value = numerical_approx(spherical_hankel1(n, x))
            !  >>>     print("n = {}: {}".format(n, value))
            ground_truth_y(1) = cmplx(-0.191784854932628_8, +0.0567324370926452_8, kind=8)
            ground_truth_y(2) = cmplx(-0.0950894080791708_8, -0.180438367514099_8, kind=8)
            ground_truth_y(3) = cmplx(0.134731210085125_8, -0.164995457601104_8, kind=8)
            ground_truth_y(4) = cmplx(0.229820618164296_8, +0.0154429099129942_8, kind=8)
            ground_truth_y(5) = cmplx(0.187017655344889_8, +0.186615531479296_8, kind=8)

            y = lib_math_bessel_spherical_third_kind_2_real(x, fnu, n)

            rv = .true.
            print *, "test_lib_math_bessel_spherical_third_kind_2_real:"
            do i=1, n
                buffer = abs(y(i) - ground_truth_y(i))
                if (buffer .gt. ground_truth_e) then
                    print *, "  ", i , "difference: ", buffer, " : FAILED"
                    rv = .false.
                else
                    print *, "  ", i, ": OK"
                end if
            end do

        end function test_lib_math_bessel_spherical_third_kind_2_real

        function test_lib_math_bessel_spherical_third_kind_2_cmplx() result (rv)
            implicit none
            ! dummy
            logical :: rv

            ! parameter
            integer, parameter :: n = 5

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
            ground_truth_y(1) = cmplx(0.910979344197223_8, -2.44844550451690_8, kind=8)
            ground_truth_y(2) = cmplx(2.06407896443698_8, +0.0711231320186915_8, kind=8)
            ground_truth_y(3) = cmplx(0.690422228144533_8, +0.953728630203184_8, kind=8)
            ground_truth_y(4) = cmplx(-0.00889039150233644_8, +0.258009870554623_8, kind=8)
            ground_truth_y(5) = cmplx(-0.254463139803031_8, -0.486653171603505_8, kind=8)

            y = lib_math_bessel_spherical_third_kind_2_cmplx(x, fnu, n)

            rv = .true.
            print *, "test_lib_math_bessel_spherical_third_kind_2_cmplx:"
            do i=1, n
                buffer = abs(y(i) - ground_truth_y(i))
                if (buffer .gt. ground_truth_e) then
                    print *, "  ", i , "difference: ", buffer, " : FAILED"
                    rv = .false.
                else
                    print *, "  ", i, ": OK"
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

        function test_lib_math_bessel_riccati_zetta_real() result(rv)
            implicit none
            ! dummy
            logical :: rv

            ! parameter
            integer, parameter :: n = 5

            ! auxiliary
            double precision :: x = 4
            complex(kind=8), dimension(n) :: zetta
            integer :: fnu = 1

            complex(kind=8), dimension(n) :: ground_truth_zetta
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

            ground_truth_zetta(1) = cmplx(0.464442997036630_8, -0.920213400523831_8, kind=8)
            ground_truth_zetta(2) = cmplx(1.10513474308540_8, -0.0365164295292615_8, kind=8)
            ground_truth_zetta(3) = cmplx(0.916975431820121_8, +0.874567863612254_8, kind=8)
            ground_truth_zetta(4) = cmplx(0.499572262599811_8, +1.56701019085071_8, kind=8)
            ground_truth_zetta(5) = cmplx(0.207062159029454_8, +2.65120506580184_8, kind=8)

            zetta = lib_math_bessel_riccati_zetta_real(x, fnu, n)

            rv = .true.
            print *, "test_lib_math_bessel_riccati_zetta_real:"
            do i=1, n
                buffer = abs(zetta(i) - ground_truth_zetta(i))
                if (buffer .gt. ground_truth_e) then
                    print *, "  ", i , "difference: ", buffer, " : FAILED"
                    rv = .false.
                else
                    print *, "  ", i, ": OK"
                end if
            end do
        end function test_lib_math_bessel_riccati_zetta_real

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

            integer :: i
            double precision :: buffer

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
            end do

        end function test_lib_math_bessel_spherical_first_kind_derivative_real

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

            integer :: i
            double precision :: buffer

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
            end do

        end function test_lib_math_bessel_spherical_first_kind_derivative_cmplx

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

            integer :: i
            double precision :: buffer

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

            integer :: i
            double precision :: buffer

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
            complex(kind=8), dimension(n+2) :: h_1_n
            integer :: fnu = 1

            complex(kind=8), dimension(n) :: ground_truth_h_1

            integer :: i
            double precision :: buffer

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
            complex(kind=8), dimension(n+2) :: h_1_n
            integer :: fnu = 1

            complex(kind=8), dimension(n) :: ground_truth_h_1

            integer :: i
            double precision :: buffer

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
            complex(kind=8), dimension(n+2) :: h_2_n
            integer :: fnu = 1

            complex(kind=8), dimension(n) :: ground_truth_h_2

            integer :: i
            double precision :: buffer

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
            complex(kind=8), dimension(n+2) :: h_2_n
            integer :: fnu = 1

            complex(kind=8), dimension(n) :: ground_truth_h_2

            integer :: i
            double precision :: buffer

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

            ground_truth_s(1) = -0.872913244567086_8
            ground_truth_s(2) = -0.0881243745060704_8
            ground_truth_s(3) = 0.417403169220310_8
            ground_truth_s(4) = 0.417403169220310_8
            ground_truth_s(5) = 0.240744563812994_8

            s = lib_math_bessel_riccati_s_derivative_real(x, fnu, n, s_n)

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
            !  >>> x=4.0 + I*2.0
            !  >>> for i in range(1,6):
            !  >>>     value = numerical_approx(S(i,x))
            !  >>>     print("n = {}: {}".format(i, value))

            ground_truth_s(1) = cmplx(-2.88434028957354_8, -1.61856658499523_8, kind=8)
            ground_truth_s(2) = cmplx(0.446332618346892_8, -2.38825995948575_8, kind=8)
            ground_truth_s(3) = cmplx(1.66989124066778_8, -0.418303857101622_8, kind=8)
            ground_truth_s(4) = cmplx(0.873045509461675_8, +0.685823569513368_8, kind=8)
            ground_truth_s(5) = cmplx(0.0995240805571297_8, +0.579739235324844_8, kind=8)

            s = lib_math_bessel_riccati_s_derivative_cmplx(x, fnu, n, s_n)

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
            end do
        end function test_lib_math_bessel_riccati_c_derivative_real

!        function test_lib_math_bessel_riccati_xi_derivative_real() result(rv)
!            implicit none
!            ! dummy
!            logical :: rv
!
!            ! parameter
!            integer, parameter :: n = 5
!
!            ! auxiliary
!            double precision :: x = 4
!            complex(kind=8), dimension(n) :: xi
!            integer :: fnu = 1
!
!            complex(kind=8), dimension(n) :: ground_truth_xi
!            double precision :: ground_truth_e = 10.0_8**(-13.0_8)
!
!            integer :: i
!            double precision :: buffer
!
!            ! Values were generated with sageMath
!            !
!            ! source code:
!            !  >>> var('x')
!            !  >>> var('n')
!            !  >>> S(n,x) = x*spherical_bessel_J(n,x)
!            !  >>>
!            !  >>> from IPython.display import Math
!            !  >>> t = latex(S)
!            !  >>> display(Math("""{}""".format(t)))
!            !  >>>
!            !  >>> x=4.0
!            !  >>> for i in range(1,6):
!            !  >>>     value = numerical_approx(S(i,x))
!            !  >>>     print("n = {}: {}".format(i, value))
!
!            ground_truth_xi(1) = cmplx(0.464442997036630_8, +0.920213400523831_8, kind=8)
!            ground_truth_xi(2) = cmplx(1.10513474308540_8, +0.0365164295292615_8, kind=8)
!            ground_truth_xi(3) = cmplx(0.916975431820121_8, -0.874567863612254_8, kind=8)
!            ground_truth_xi(4) = cmplx(0.499572262599811_8, -1.56701019085071_8, kind=8)
!            ground_truth_xi(5) = cmplx(0.207062159029454_8, -2.65120506580184_8, kind=8)
!
!            xi = lib_math_bessel_riccati_xi_real(x, fnu, n)
!
!            rv = .true.
!            print *, "test_lib_math_bessel_riccati_xi_derivative_real:"
!            do i=1, n
!                buffer = abs(xi(i) - ground_truth_xi(i))
!                if (buffer .gt. ground_truth_e) then
!                    print *, "  ", i , "difference: ", buffer, " : FAILED"
!                    rv = .false.
!                else
!                    print *, "  ", i, ": OK"
!                end if
!            end do
!        end function test_lib_math_bessel_riccati_xi_derivative_real
!
!        function test_lib_math_bessel_riccati_zetta_derivative_real() result(rv)
!            implicit none
!            ! dummy
!            logical :: rv
!
!            ! parameter
!            integer, parameter :: n = 5
!
!            ! auxiliary
!            double precision :: x = 4
!            complex(kind=8), dimension(n) :: zetta
!            integer :: fnu = 1
!
!            complex(kind=8), dimension(n) :: ground_truth_zetta
!            double precision :: ground_truth_e = 10.0_8**(-13.0_8)
!
!            integer :: i
!            double precision :: buffer
!
!            ! Values were generated with sageMath
!            !
!            ! source code:
!            !  >>> var('x')
!            !  >>> var('n')
!            !  >>> S(n,x) = x*spherical_bessel_J(n,x)
!            !  >>>
!            !  >>> from IPython.display import Math
!            !  >>> t = latex(S)
!            !  >>> display(Math("""{}""".format(t)))
!            !  >>>
!            !  >>> x=4.0
!            !  >>> for i in range(1,6):
!            !  >>>     value = numerical_approx(S(i,x))
!            !  >>>     print("n = {}: {}".format(i, value))
!
!            ground_truth_zetta(1) = cmplx(0.464442997036630_8, -0.920213400523831_8, kind=8)
!            ground_truth_zetta(2) = cmplx(1.10513474308540_8, -0.0365164295292615_8, kind=8)
!            ground_truth_zetta(3) = cmplx(0.916975431820121_8, +0.874567863612254_8, kind=8)
!            ground_truth_zetta(4) = cmplx(0.499572262599811_8, +1.56701019085071_8, kind=8)
!            ground_truth_zetta(5) = cmplx(0.207062159029454_8, +2.65120506580184_8, kind=8)
!
!            zetta = lib_math_bessel_riccati_zetta_real(x, fnu, n)
!
!            rv = .true.
!            print *, "test_lib_math_bessel_riccati_zetta_derivative_real:"
!            do i=1, n
!                buffer = abs(zetta(i) - ground_truth_zetta(i))
!                if (buffer .gt. ground_truth_e) then
!                    print *, "  ", i , "difference: ", buffer, " : FAILED"
!                    rv = .false.
!                else
!                    print *, "  ", i, ": OK"
!                end if
!            end do
!        end function test_lib_math_bessel_riccati_zetta_derivative_real

    end function lib_math_bessel_test_functions

end module lib_math_bessel
