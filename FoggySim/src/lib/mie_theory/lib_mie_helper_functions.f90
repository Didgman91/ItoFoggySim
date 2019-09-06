module lib_mie_helper_functions
    use libmath
    use lib_mie_type
    use lib_mie_vector_spherical_harmonics
    implicit none

    private

    public :: lib_mie_hf_get_An
    public :: lib_mie_hf_get_n_c
    public :: lib_mie_hf_calc_triple_sum

    interface lib_mie_hf_get_An
        module procedure get_An_real
        module procedure get_An_cmplx
    end interface lib_mie_hf_get_An

    contains

        ! Calculates the logarithmic derivative An
        !
        ! Formula: A_n = [mx j_n(mx)]' / (mx j_n(mx))
        !
        ! Argument
        ! ----
        !   x: real
        !       size parameter: x = k*r = 2 * PI * N * r / lambda
        !       k: wavenumber
        !       r: distance
        !       N: refractive index of the medium
        !       lambda: wave length
        !   m: real
        !       relative refractive index: m = N_1 / N
        !       N_1: refractive index of the particle
        !       N: refractive index of the medium
        !   fnu: integer
        !       order of initial function, fnu.GE.0
        !   n: integer
        !       number of members of the sequence, n.GE.1
        !
        ! Returns
        ! ----
        !   An: complex(fnu:fnu + n - 1)
        !       value of the logarithmic derivative of the last interation
        !
        ! Reference: Light Scattering by Particles, Barber, Hill, eq. (4.19)
        function get_An_real(x, m, fnu, n) result (An)
            implicit none
            ! dummy
            real(kind=8) :: x
            real(kind=8) :: m
            integer :: fnu
            integer :: n

            real(kind=8), dimension(fnu:fnu+n-1) :: An

            ! parameter

            ! auxiliary
            integer :: i
            integer :: m_n_max
            integer :: m_n_mx
            real(kind=8) :: m_mx
            real(kind=8) :: m_buffer


            m_mx = m * x
            m_n_max = fnu + n - 1

            ! eq. (4.20)
            m_n_mx = max(lib_mie_hf_get_n_c(abs(x)), int(abs(m_mx))) + 15

            m_buffer = cmplx(0.0, 0.0, kind=8)
            do i=m_n_mx, fnu+1, -1
                m_buffer = get_An_minus_1_real(i, m_mx, m_buffer)

                if (i .le. m_n_max+1) then
                    An(i-1) = m_buffer
                end if
            end do

        end function get_An_real

        ! Calculates the logarithmic derivative An
        !
        ! Formula: A_n = [mx j_n(mx)]' / (mx j_n(mx))
        !
        ! Argument
        ! ----
        !   x: real
        !       size parameter: x = k*r = 2 * PI * r / lambda
        !       k: wavenumber
        !       r: distance
        !       N: refractive index of the medium
        !       lambda: wave length
        !   m: complex
        !       relative refractive index: m = N_1 / N
        !       N_1: refractive index of the particle
        !       N: refractive index of the medium
        !   fnu: integer
        !       order of initial function, fnu.GE.0
        !   n: integer
        !       number of members of the sequence, n.GE.1
        !
        ! Returns
        ! ----
        !   An: complex(fnu:fnu + n - 1)
        !       value of the logarithmic derivative of the last interation
        !
        ! Reference: Light Scattering by Particles, Barber, Hill, eq. (4.19)
        function get_An_cmplx(x, m, fnu, n) result (An)
            implicit none
            ! dummy
            double precision :: x
            complex(kind=8) :: m
            integer :: fnu
            integer :: n

            complex(kind=8), dimension(fnu:fnu+n-1) :: An

            ! parameter

            ! auxiliary
            integer :: i
            integer :: m_n_max
            integer :: m_n_mx
            complex(kind=8) :: m_mx
            complex(kind=8) :: m_buffer


            m_mx = m * x
            m_n_max = fnu + n - 1

            ! eq. (4.20)
            m_n_mx = max(lib_mie_hf_get_n_c(abs(x)), int(abs(m_mx))) + 15

            m_buffer = cmplx(0.0, 0.0, kind=8)
            do i=m_n_mx, fnu+1, -1
                m_buffer = get_An_minus_1_cmplx(i, m_mx, m_buffer)

                if (i .le. m_n_max+1) then
                    An(i-1) = m_buffer
                end if
            end do

        end function get_An_cmplx

        ! Argument
        ! ----
        !   n: integer
        !       degree of the logarithmic derivative
        !   mx: complex
        !       x: real
        !           size parameter: x = k*r = 2 * PI * N * r / lambda
        !           k: wavenumber
        !           r: distance
        !           N: refractive index of the medium
        !           lambda: wave length
        !       m: real
        !           relative refractive index: m = N_1 / N
        !           N_1: refractive index of the particle
        !           N: refractive index of the medium
        !   An: complex
        !       value of the logarithmic derivative of the last interation
        !
        ! Reference: Light Scattering by Particles, Barber, Hill, eq. (4.19)
        function get_An_minus_1_real(n, mx, An) result (rv)
            implicit none
            ! dummy
            integer :: n
            real(kind=8) :: mx
            real(kind=8) :: An

            real(kind=8) :: rv

            ! auxiliary
            real(kind=8) :: n_div_mx

            n_div_mx = real(n, kind=8) / mx

            rv = n_div_mx - 1.0_8/(An + n_div_mx)

        end function get_An_minus_1_real

        ! Argument
        ! ----
        !   n: integer
        !       degree of the logarithmic derivative
        !   mx: complex
        !       x: complex
        !           size parameter: x = k*r = 2 * PI * N * r / lambda
        !           k: wavenumber
        !           r: distance
        !           N: refractive index of the medium
        !           lambda: wave length
        !       m: complex
        !           relative refractive index: m = N_1 / N
        !           N_1: refractive index of the particle
        !           N: refractive index of the medium
        !   An: complex
        !       value of the logarithmic derivative of the last interation
        !
        ! Reference: Light Scattering by Particles, Barber, Hill, eq. (4.19)
        function get_An_minus_1_cmplx(n, mx, An) result (rv)
            implicit none
            ! dummy
            integer :: n
            complex(kind=8) :: mx
            complex(kind=8) :: An

            complex(kind=8) :: rv

            ! auxiliary
            complex(kind=8) :: n_div_mx

            n_div_mx = real(n, kind=8) / mx

            rv = n_div_mx - 1.0_8/(An + n_div_mx)

        end function get_An_minus_1_cmplx

        ! calculates highest degree n for a convergent algorithm
        !
        ! Argument
        ! ----
        !   x: double precision
        !       size parameter
        !       x = r / lambda
        !       r: radius
        !       lambda: wavelength
        !
        ! Reference: Light Scattering by Particles: Computational Methods, PW Barber, S C Hill, eq. 4.16
        function lib_mie_hf_get_n_c(x) result (rv)
            implicit none
            ! dummy
            double precision :: x

            integer :: rv

            ! auxiliary
            double precision :: dummy

            dummy = x + 4.05_8 * x**(1.0_8/3.0_8) + 2.0_8

            rv = int(ceiling(dummy))

        end function lib_mie_hf_get_n_c

        subroutine lib_mie_hf_calc_triple_sum(lambda, n_medium, &
                                              sphere, sphere_parameter, sphere_j, &
                                              z_selector, &
                                              a, b, &
                                              a_old, b_old, &
                                              f)
            implicit none
            ! dummy
            double precision, intent(in) :: lambda
            double precision, intent(in) :: n_medium
            type(sphere_type), dimension(:), intent(in) :: sphere
            type(sphere_parameter_type), dimension(:), intent(in) :: sphere_parameter
            integer :: sphere_j
            integer(kind=1) :: z_selector

            type(list_list_cmplx), intent(inout) :: a
            type(list_list_cmplx), intent(inout) :: b

            type(list_list_cmplx), intent(in), optional :: a_old
            type(list_list_cmplx), intent(in), optional :: b_old
            double precision, intent(in), optional :: f

            ! auxiliary
            integer :: m
            integer :: n

            integer :: l
            integer :: i

            integer, dimension(2) :: n_range
            integer, dimension(2) :: nu_range

            type(list_list_cmplx) :: p
            type(list_list_cmplx) :: q
            type(list_cmplx) :: a_j_n
            type(list_cmplx) :: b_j_n
            type(list_cmplx) :: a_l_nu
            type(list_cmplx) :: b_l_nu
            type(list_list_cmplx) :: a_l_numu
            type(list_list_cmplx) :: b_l_numu
            type(list_4_cmplx) :: A_mnmunu
            type(list_4_cmplx) :: B_mnmunu
            type(cartesian_coordinate_real_type) :: d_0_j
            type(cartesian_coordinate_real_type) :: d_0_l
            type(cartesian_coordinate_real_type) :: x

            double complex, dimension(:), allocatable :: buffer_sum_a
            double complex, dimension(:), allocatable :: buffer_sum_b

            double complex :: buffer_cmplx

            i = sphere(sphere_j)%sphere_parameter_index
            n_range = sphere_parameter(i)%n_range
            d_0_j = sphere(sphere_j)%d_0_j

            allocate(buffer_sum_a(lbound(sphere, 1):ubound(sphere, 1)))
            allocate(buffer_sum_b(lbound(sphere, 1):ubound(sphere, 1)))

            ! calculate the triple sum
            do l=ubound(sphere, 1), lbound(sphere, 1)
                if (l .ne. sphere_j) then
                    i = sphere(l)%sphere_parameter_index
                    nu_range = sphere_parameter(i)%n_range
                    a_l_nu = sphere_parameter(i)%a_n
                    b_l_nu = sphere_parameter(i)%b_n
                    d_0_l = sphere(l)%d_0_j

                    x = d_0_l - d_0_j

!                    do nu=nu_range(1), nu_range(2)
!                        a_l_numu%item(nu)%item(:) = a_l_nu%item(nu) *
!                    end do
!
!                    call lib_mie_vector_spherical_harmonics_translation_coefficient(x, &
!                                                                                    n_range, nu_range, z_selector,&
!                                                                                    A_mnmunu, B_mnmunu)
!
!
!                    call calc_inner_sum(m ,n, &
!                         a_l_munu, b_L_munu, a_munumn, b_munumn, &
!                         nu_range, &
!                         buffer_sum_a(l), buffer_sum_b(l))

                else
                    buffer_sum_a(l) = cmplx(0,0)
                    buffer_sum_b(l) = cmplx(0,0)
                end if
            end do
        end subroutine

        ! Calculates the two inner sums of eq. 30 (interactive scattering coefficients)
        !
        ! Argument
        ! ----
        !   a_l_munu: type(list_list_cmplx)
        !       Mie coefficient of the l-th sphere
        !       restriction: same dimension as b_l_munu
        !   b_l_munu: type(list_list_cmplx)
        !       Mie coefficient of the l-th sphere
        !       restriction: same dimension as a_l_munu
        !   a_munumn: type(list_4_cmplx)
        !       translation coefficient
        !       restriction: same dimension as b_munumn
        !   b_munumn: type(list_4_cmplx)
        !       translation coefficient
        !       restriction: same dimension as a_munumn
        !   nu_range: integer, dimension(2), optional(std: bounds of a_munumn)
        !       first element: start index
        !       second element: last index
        !       CONDITION:
        !           - first element .le. second element
        !           - 0 <= nu
        !
        ! Returns
        ! ----
        !   a_2sum: double complex
        !
        !
        ! Reference:Electromagnetic scattering by an aggregate of spheres, Yu-lin Xu, eq. 30
        subroutine calc_inner_sum(m, n, &
                                             a_l_munu, b_L_munu, a_munumn, b_munumn, &
                                             nu_range, &
                                             a_2sum, b_2sum)
            implicit none
            ! dummy
            integer, intent(in) :: m
            integer, intent(in) :: n
            type(list_list_cmplx), intent(in) :: a_l_munu
            type(list_list_cmplx), intent(in) :: b_l_munu
            type(list_4_cmplx), intent(in) :: a_munumn
            type(list_4_cmplx), intent(in) :: b_munumn

            integer, dimension(2), intent(in), optional :: nu_range

            double complex, intent(inout) :: a_2sum
            double complex, intent(inout) :: b_2sum

            ! auxiliary
            integer :: mu
            integer :: nu

            integer, dimension(2) :: m_nu_range

            double complex :: buffer_cmplx

            ! set boundaries
            if (present(nu_range)) then
                m_nu_range = nu_range
            else
                m_nu_range(1) = max(lbound(a_munumn%item(n)%item(m)%item, 1), 1)
                m_nu_range(2) = ubound(a_munumn%item(n)%item(m)%item, 1)
            end if

            ! calculate the summand
            do nu=nu_range(1), nu_range(2)
                do mu=-nu, nu
                    ! summand of a_j_mn
                    ! first line eq. 30
                    buffer_cmplx = a_l_munu%item(nu)%item(mu) * a_munumn%item(n)%item(m)%item(nu)%item(mu)
                    buffer_cmplx = buffer_cmplx &
                                   + b_l_munu%item(nu)%item(mu) * b_munumn%item(n)%item(m)%item(nu)%item(mu)
                    a_2sum = buffer_cmplx

                    ! summand of b_j_mn
                    ! second line eq. 30
                    buffer_cmplx = a_l_munu%item(nu)%item(mu) * b_munumn%item(n)%item(m)%item(nu)%item(mu)
                    buffer_cmplx = buffer_cmplx &
                                   + b_l_munu%item(nu)%item(mu) * a_munumn%item(n)%item(m)%item(nu)%item(mu)
                    b_2sum = buffer_cmplx
                end do
            end do

        end subroutine calc_inner_sum

end module lib_mie_helper_functions
