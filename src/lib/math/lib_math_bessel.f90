module lib_math_bessel
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

    interface lib_math_bessel_spherical_third_kind_1
        module procedure lib_math_bessel_spherical_third_kind_1_real
        module procedure lib_math_bessel_spherical_third_kind_1_cmplx
    end interface

    interface lib_math_bessel_spherical_third_kind_2
        module procedure lib_math_bessel_spherical_third_kind_2_real
    end interface

    interface lib_math_bessel_spherical_first_kind_derivative
        module procedure lib_math_bessel_spherical_first_kind_real_derivative
    end interface

    interface lib_math_bessel_spherical_second_kind_derivative
        module procedure lib_math_bessel_spherical_second_kind_real_derivative
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
        integer, intent(in) :: fnu
        integer, intent(in) :: n

        double precision, dimension(n) :: rv

        ! auxiliary
        double precision :: rv_1

        ! nz: integer
        !    number of components of Y set to zero due to
        !    underflow,
        !    NZ=0   , normal return, computation completed
        !    NZ .NE. 0, last NZ components of Y set to zero,
        !             Y(K)=0.0D0, K=N-NZ+1,...,N.
        integer :: nz

!        rv = sqrt(PI/(2.D0*x)) * bessel_jn(real(n)+0.5,x)
        rv_1 = sqrt(PI/(2.D0*x))

        call DBESJ (X, FNU+0.5D0, N, rv, nz)

#ifdef _DEBUG_
        if (nz .ne. 0) then
            print *, "lib_math_bessel_spherical_first_kind_real: WARNING"
            print *, "  number of components of Y set to zero due to underflow: ", nz
        end if
#endif

        rv(:) = rv_1 * rv(:)

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
    function lib_math_bessel_spherical_first_kind_cmplx(x, fnu, n, kode) result (rv)
    !
        implicit none
        ! dummy
        integer, intent(in) :: fnu
        integer, intent(in) :: n
        complex(kind=8), intent(in) :: x
        integer, optional, intent(in) :: kode

        complex(kind=8), dimension(n) :: rv

        ! auxiliary
        integer :: m_kode = 1
        double precision, dimension(n) :: cyr
        double precision, dimension(n) :: cyi
        integer :: nz
        integer :: ierr

        complex(kind=8) :: rv_1

        if (present(kode)) then
            m_kode = kode
        end if

        ! sqrt(PI/(2.D0*x)) * bessel_jn(real(n)+0.5,x)
        rv_1 = sqrt(PI/(2.D0*x))

        call ZBESJ(REALPART(x), imagpart(x), real(fnu, 8)+0.5D0, m_kode, n, cyr, cyi, nz, ierr)

        if (ierr .ne. 0) then

            print *, "lib_math_bessel_spherical_first_kind_cmplx: ERROR"

            select case (ierr)
                case(0)
                print *, "    IERR=0, NORMAL RETURN - COMPUTATION COMPLETED"
                case(1)
                print *, "    IERR=1, INPUT ERROR   - NO COMPUTATION"
                case(2)
                print *, "    IERR=2, OVERFLOW      - NO COMPUTATION, AIMAG(Z)"
                print *, "            TOO LARGE ON KODE=1"
                case(3)
                print *, "    IERR=3, CABS(Z) OR FNU+N-1 LARGE - COMPUTATION DONE"
                print *, "            BUT LOSSES OF SIGNIFCANCE BY ARGUMENT"
                print *, "            REDUCTION PRODUCE LESS THAN HALF OF MACHINE"
                print *, "            ACCURACY"
                case(4)
                print *, "    IERR=4, CABS(Z) OR FNU+N-1 TOO LARGE - NO COMPUTA-"
                print *, "            TION BECAUSE OF COMPLETE LOSSES OF SIGNIFI-"
                print *, "            CANCE BY ARGUMENT REDUCTION"
                case(5)
                print *, "    IERR=5, ERROR              - NO COMPUTATION,"
                print *, "            ALGORITHM TERMINATION CONDITION NOT MET"
            end select
        end if

#ifdef _DEBUG_
        if (nz .ne. 0) then
            print *, "lib_math_bessel_spherical_first_kind_real: WARNING"
            print *, "  number of components set to zero due to underflow: ", nz
        end if
#endif

        rv(:) = cmplx(cyr(:), cyi(:), kind=8)

        rv(:) = rv_1 * rv(:)

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
        integer, intent(in) :: fnu
        integer, intent(in) :: n

        double precision, dimension(n) :: rv

        ! auxiliary
        double precision :: rv_1

!        rv = sqrt(PI/(2.D0*x)) * bessel_yn(real(n)+0.5,x)
        rv_1 = sqrt(PI/(2.D0*x))

        call DBESY (X, FNU+0.5D0, N, rv)

        rv(:) = rv_1 * rv(:)

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
    function lib_math_bessel_spherical_second_kind_cmplx(x, fnu, n, kode) result (rv)
    !
        implicit none
        ! dummy
        integer, intent(in) :: fnu
        integer, intent(in) :: n
        complex(kind=8), intent(in) :: x
        integer, optional, intent(in) :: kode

        complex(kind=8), dimension(n) :: rv

        ! auxiliary
        integer :: m_kode = 1
        double precision, dimension(n) :: cyr
        double precision, dimension(n) :: cyi
        integer :: nz
        integer :: ierr
        complex(kind=8) :: CWRKR(n)
        complex(kind=8) :: CWRKI(n)

        complex(kind=8) :: rv_1

        if (present(kode)) then
            m_kode = kode
        end if

        ! sqrt(PI/(2.D0*x)) * bessel_jn(real(n)+0.5,x)
        rv_1 = sqrt(PI/(2.D0*x))

        call ZBESY(REALPART(x), imagpart(x), real(fnu, 8)+0.5D0, m_kode, n, cyr, cyi, nz, &
                   CWRKR, CWRKI, ierr)

        if (ierr .ne. 0) then

            print *, "lib_math_bessel_spherical_second_kind_cmplx: ERROR"

            select case (ierr)
                case(0)
                print *, "    IERR=0, NORMAL RETURN - COMPUTATION COMPLETED"
                case(1)
                print *, "    IERR=1, INPUT ERROR   - NO COMPUTATION"
                case(2)
                print *, "    IERR=2, OVERFLOW      - NO COMPUTATION, AIMAG(Z)"
                print *, "            TOO LARGE ON KODE=1"
                case(3)
                print *, "    IERR=3, CABS(Z) OR FNU+N-1 LARGE - COMPUTATION DONE"
                print *, "            BUT LOSSES OF SIGNIFCANCE BY ARGUMENT"
                print *, "            REDUCTION PRODUCE LESS THAN HALF OF MACHINE"
                print *, "            ACCURACY"
                case(4)
                print *, "    IERR=4, CABS(Z) OR FNU+N-1 TOO LARGE - NO COMPUTA-"
                print *, "            TION BECAUSE OF COMPLETE LOSSES OF SIGNIFI-"
                print *, "            CANCE BY ARGUMENT REDUCTION"
                case(5)
                print *, "    IERR=5, ERROR              - NO COMPUTATION,"
                print *, "            ALGORITHM TERMINATION CONDITION NOT MET"
            end select
        end if

#ifdef _DEBUG_
        if (nz .ne. 0) then
            print *, "lib_math_bessel_spherical_first_kind_real: WARNING"
            print *, "  number of components set to zero due to underflow: ", nz
        end if
#endif

        rv(:) = cmplx(cyr(:), cyi(:), kind=8)

        rv(:) = rv_1 * rv(:)

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
        integer, intent(in) :: fnu
        integer, intent(in) :: n
        complex(kind=8), dimension(n) :: rv

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
    function lib_math_bessel_spherical_third_kind_1_cmplx(x, fnu, n) result (rv)
        implicit none
        ! dummy
        complex(kind=8), intent(in) :: x
        integer, intent(in) :: fnu
        integer, intent(in) :: n
        complex(kind=8), dimension(n) :: rv

        rv = lib_math_bessel_spherical_first_kind_cmplx(x, fnu, n) &
             +lib_math_bessel_spherical_second_kind_cmplx(x, fnu, n)

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
        integer, intent(in) :: fnu
        integer, intent(in) :: n
        complex(kind=8), dimension(n) :: rv

        rv = cmplx(lib_math_bessel_spherical_first_kind_real(x, fnu, n), &
                   -lib_math_bessel_spherical_second_kind_real(x, fnu, n), &
                   kind=8)

    end function lib_math_bessel_spherical_third_kind_2_real

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
    ! LaTeX: $$ j_{n}^{\prime}=\left(\frac{n}{x}\right) j_{n}(x)-j_{n+1}(x) $$
    !
    ! Reference: https://www.boost.org/doc/libs/1_57_0/libs/math/doc/html/math_toolkit/bessel/bessel_derivatives.html
    !
    function lib_math_bessel_spherical_first_kind_real_derivative(x, fnu, n, j_n) result (rv)
        implicit none
        ! dummy
        double precision, intent(in) :: x
        integer, intent(in) :: fnu
        integer, intent(in) :: n
        double precision, dimension(n+1), intent(inout) :: j_n
        double precision, dimension(n) :: rv

        ! auxiliary
        double precision :: rv_1
        integer :: i

        j_n = lib_math_bessel_spherical_first_kind_real(x, fnu, n+1)
        rv_1 = real(n, 8) / x

        do i=1, n
            rv(i) = rv_1 * j_n(i) - j_n(i+1)
        end do

    end function lib_math_bessel_spherical_first_kind_real_derivative

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
    ! LaTeX: $$ y_{n}^{\prime}=\left(\frac{n}{x}\right) y_{n}(x)-y_{n+1}(x) $$
    !
    ! Reference: https://www.boost.org/doc/libs/1_57_0/libs/math/doc/html/math_toolkit/bessel/bessel_derivatives.html
    !
    function lib_math_bessel_spherical_second_kind_real_derivative(x, fnu, n, y_n) result (rv)
        implicit none
        ! dummy
        double precision, intent(in) :: x
        integer, intent(in) :: fnu
        integer, intent(in) :: n
        double precision, dimension(n+1), intent(inout) :: y_n
        double precision, dimension(n) :: rv

        ! auxiliary
        double precision :: rv_1
        integer :: i

        y_n = lib_math_bessel_spherical_second_kind_real(x, fnu, n+1)
        rv_1 = real(n, 8) / x

        do i=1, n
            rv(i) = rv_1 * y_n(i) - y_n(i+1)
        end do

    end function lib_math_bessel_spherical_second_kind_real_derivative

    function lib_math_bessel_test_functions() result (rv)
        implicit none
        ! dummy
        integer :: rv

        rv = 0

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
        if (.not. test_lib_math_bessel_spherical_first_kind_cmplx()) then
            rv = rv + 1
        end if
        if (.not. test_lib_math_bessel_spherical_third_kind_1_cmplx()) then
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
            double precision, dimension(n) :: y
            integer :: fnu = 0

            double precision, dimension(n) :: ground_truth_y

            integer :: i
            double precision :: buffer

            ! Values were generated with sageMath
            !
            ! source code:
            !  >>> x=5.0 +I*0.0
            !  >>> for n in range(0,5):
            !  >>>     value = numerical_approx(spherical_bessel_J(n, x))
            !  >>>     print("n = {}: {}".format(n, value))
            data ground_truth_y / -0.191784854932628_8, -0.0950894080791708_8, 0.134731210085125_8,&
                                  0.229820618164296_8, 0.187017655344889_8/

            y = lib_math_bessel_spherical_first_kind_real(x, fnu, n)

            rv = .true.
            print *, "test_lib_math_bessel_spherical_first_kind_real:"
            do i=1, n
                buffer = y(i) - ground_truth_y(i)
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

        end function

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

        end function

        function test_lib_math_bessel_spherical_third_kind_1_cmplx() result (rv)
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
            ground_truth_y(1) = cmplx(0.0448448436994762_8, -0.0166851687034122_8, kind=8)
            ground_truth_y(2) = cmplx(-0.00964524995439623_8, -0.0602273468001982_8, kind=8)
            ground_truth_y(3) = cmplx(-0.0972492912654220_8, -0.0212514039309393_8, kind=8)
            ground_truth_y(4) = cmplx(-0.138480619041055_8, +0.155224705968302_8, kind=8)
            ground_truth_y(5) = cmplx(0.126551443388103_8, +0.535235722697314_8, kind=8)

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

        end function

        function test_lib_math_bessel_spherical_first_kind_real_derivative() result (rv)
            implicit none
            ! dummy
            logical :: rv

            ! parameter
            integer, parameter :: n = 5
            double precision, parameter :: ground_truth_e = 10.0_8**(-13.0_8)

            ! auxiliary
            double precision :: x = 5
            double precision, dimension(n) :: y
            double precision, dimension(n+1) :: y_n
            integer :: fnu = 0

            double precision, dimension(n) :: ground_truth_y

            integer :: i
            double precision :: buffer

            ! Values were generated with sageMath
            !
            ! source code:
            !  >>> x=5.0 +I*0.0
            !  >>> for n in range(0,5):
            !  >>>     value = numerical_approx(spherical_bessel_J(n, x))
            !  >>>     print("n = {}: {}".format(n, value))
            data ground_truth_y / -0.191784854932628_8, -0.0950894080791708_8, 0.134731210085125_8,&
                                  0.229820618164296_8, 0.187017655344889_8/

            y = lib_math_bessel_spherical_first_kind_real_derivative(x, fnu, n, y_n)

            rv = .true.
            print *, "test_lib_math_bessel_spherical_first_kind_real_derivative:"
            do i=1, n
                buffer = y(i) - ground_truth_y(i)
                if (abs(buffer) .gt. ground_truth_e) then
                    print *, "  ", i , "difference: ", buffer, " : FAILED"
                    rv = .false.
                else
                    print *, "  ", i, ": OK"
                end if
            end do

        end function

    end function lib_math_bessel_test_functions

end module lib_math_bessel
