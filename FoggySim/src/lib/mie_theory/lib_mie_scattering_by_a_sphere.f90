! Scattering by a sphere
! ----
!
! - incident plane wave along the z-axis
!
#define _DEBUG_

module lib_mie_scattering_by_a_sphere
    use libmath
    use lib_mie_vector_spherical_harmonics
    implicit none

    private

    ! --- public ---
    public:: lib_mie_scattering_by_a_sphere_test_functions

    contains

!        ! calculates the scatterd electical field of a sphere
!        !
!        ! Setup
!        ! ----
!        !
!        !        ^ z
!        !        |
!        !        o -> x
!        !    __________
!        !    __________
!        !
!        !    o: sphere
!        !    _: plane wave
!        !       - propagation direction: z
!        !
!        ! Argument
!        ! ----
!        !   theta: double precision
!        !       polar angle
!        !   phi: double precision
!        !       azimuthal angle
!        !   rho: double precision
!        !       dimensionless varibale rho = k*a
!        !       k: wavenumber
!        !       a: distance
!        !   rho_particle: double precision
!        !       dimensionless varibale rho = k*r
!        !       k: wavenumber
!        !       r: radius of the sphere
!        !   e_field_0: double precision
!        !       amplitude of the incident wave
!        !   n_particle: double precision
!        !       refractive index of the sphere
!        !   n_medium: double precision
!        !       refractive index of the medium
!        !
!        ! Reference: Absorption and Scattering of Light by Small Particles, Bohren & Huffman
!        !
!        function get_e_field_scattered_bh(theta, phi, rho, e_field_0, rho_particle, n_particle, n_medium) result (e_field_s)
!            implicit none
!            ! dummy
!            double precision, intent(in) :: theta
!            double precision, intent(in) :: phi
!            double precision, intent(in) :: rho
!            double precision, intent(in) :: rho_particle
!            double precision, intent(in) :: e_field_0
!            double precision, intent(in) :: n_particle
!            double precision, intent(in) :: n_medium
!
!            type(spherical_coordinate_cmplx_type) :: e_field_s
!
!            ! parameter
!            integer(kind=4), dimension(2), parameter :: n = (/1, 10/)
!
!            ! auxiliary
!            double precision :: buffer_real
!            integer(kind=4) :: i
!            double precision :: mu
!            double precision :: mu1
!
!            complex(kind=8), dimension(n(2)-n(1)+1) :: a_n
!            complex(kind=8), dimension(n(2)-n(1)+1) :: b_n
!
!
!            integer(kind=4), dimension(2) :: m
!            integer(kind=1) :: z_selector
!
!            type(list_spherical_coordinate_cmplx_type), dimension(:), allocatable :: M_emn
!            type(list_spherical_coordinate_cmplx_type), dimension(:), allocatable :: M_omn
!            type(list_spherical_coordinate_cmplx_type), dimension(:), allocatable :: N_emn
!            type(list_spherical_coordinate_cmplx_type), dimension(:), allocatable :: N_omn
!
!            complex(kind=8), dimension(n(2)-n(1)+1) :: e_field_n
!            type(spherical_coordinate_cmplx_type), dimension(n(2)-n(1)+1) :: e_field_n_s
!
!            m = (/1,1/)
!
!            mu = 1
!            mu1 = 1
!
!            call get_coefficients_a_b_real(rho_particle, n_particle/n_medium, mu, mu1, n, a_n, b_n)
!
!            z_selector = 3
!
!            call lib_mie_vector_spherical_harmonics_components(theta, phi, rho, m, n, z_selector, &
!                                                               M_emn, M_omn, N_emn, N_omn, &
!                                                               not_calc_Memn=.true., &
!                                                               not_calc_Momn=.false., &
!                                                               not_calc_Nemn=.false., &
!                                                               not_calc_Nomn=.true.)
!            ! page 93
!            do i=n(1), n(2)
!                buffer_real = e_field_0 * (2*i+1)/(i*(i+1))
!                e_field_n(i) = cmplx(buffer_real,0, kind=8)
!                e_field_n(i) = e_field_n(i) * (cmplx(0,1, kind=8)**i)
!            end do
!
!            ! eq. 4.45
!            do i=n(1), n(2)
!                e_field_n_s(i) = e_field_n(i) * (cmplx(0,1, kind=8) * a_n(i) * N_emn(m(1))%coordinate(i)&
!                                                 - b_n(i) * M_omn(m(1))%coordinate(i))
!            end do
!
!
!            e_field_s%theta = cmplx(0,0,kind=8)
!            e_field_s%phi = cmplx(0,0,kind=8)
!            e_field_s%rho = cmplx(0,0,kind=8)
!            do i=n(1), n(2)
!                e_field_s = e_field_s + e_field_n_s(i)
!            end do
!
!            print*, "done"
!
!        end function

        ! calculates the scatterd electical field of a sphere
        !
        ! Setup
        ! ----
        !
        !        ^ z
        !        |
        !        o -> x
        !    __________
        !    __________
        !
        !    o: sphere
        !    _: plane wave
        !       - propagation direction: z
        !
        ! Argument
        ! ----
        !   theta: double precision
        !       polar angle [radian]
        !   phi: double precision
        !       azimuthal angle [radian]
        !   r: double precision
        !       distance
        !   e_field_0: double precision
        !       amplitude of the incident wave
        !   lambda: double precision
        !       wave length
        !   n_medium: double precision
        !       refractive index of the medium
        !   r_particle: double precision
        !       radius of the sphere
        !   n_particle: double precision
        !       refractive index of the sphere
        !   n_range: integer, dimension(2)
        !       first and last index (degree) of the sum to calculate the electical field
        !       e.g. n = (/ 1, 80 /) <-- size parameter
        !
        ! Reference: Electromagnetic scattering by an aggregate of spheres, Yu-lin Xu
        !
        function get_e_field_scattered_xu(theta, phi, r, e_field_0, lambda, n_medium, r_particle, n_particle, n_range) &
                                         result (e_field_s)
            implicit none
            ! dummy
            double precision, intent(in) :: theta
            double precision, intent(in) :: phi
            double precision, intent(in) :: r
            double precision, intent(in) :: e_field_0
            double precision, intent(in) :: lambda
            double precision, intent(in) :: n_medium
            double precision, intent(in) :: r_particle
            double precision, intent(in) :: n_particle
            integer(kind=4), dimension(2) :: n_range

            type(spherical_coordinate_cmplx_type) :: e_field_s

            ! auxiliary
            double precision :: k0 ! wave number = 2 pi / lambda
            double precision :: rho_particle ! r_particle * k0 * n_particle

            double precision :: buffer_real
            complex(kind=8) :: buffer_cmplx
            integer(kind=4) :: i
            integer(kind=4) :: ii
            integer(kind=4) :: m
            integer(kind=4) :: n
            double precision :: mu
            double precision :: mu1

            complex(kind=8), dimension(n_range(2)-n_range(1)+1) :: a_n
            complex(kind=8), dimension(n_range(2)-n_range(1)+1) :: b_n


            integer(kind=1) :: z_selector

            type(list_spherical_coordinate_cmplx_type), dimension(:), allocatable :: M_nm
            type(list_spherical_coordinate_cmplx_type), dimension(:), allocatable :: N_nm

            type(list_list_cmplx) :: e_field_nm
            type(spherical_coordinate_cmplx_type), dimension(n_range(2)-n_range(1)+1) :: e_field_n_s

            k0 = 2.0_8 * PI / lambda
            rho_particle = r_particle * k0 * n_particle

            mu = 1
            mu1 = 1.5

            call init_list(e_field_nm, n_range(1), n_range(2)-n_range(1)+1)

            call get_coefficients_a_b_real_barberh(rho_particle, n_particle/n_medium, n_range, a_n, b_n)
!            call get_coefficients_a_b_real(rho_particle, n_particle/n_medium, mu, mu1, n_range, a_n, b_n)

            z_selector = 3

            ! todo:
            call lib_mie_vector_spherical_harmonics_components(theta, phi, r, k0 * n_medium, n_range, z_selector, &
                                                               M_nm, N_nm)
            ! eq. (5)
            do n=n_range(1), n_range(2)
                do m=-n, n
                    buffer_real = abs(e_field_0) * real((2*n+1), kind=8) * lib_math_factorial_get_n_minus_m_divided_by_n_plus_m(n,m)

#ifdef _DEBUG_
                    if (isnan(buffer_real)) then
                        print *, "get_e_field_scattered_xu: ERROR"
                        print *, "  buffer_real is NaN"
                        print * , "  n = ", n
                        print * , "  m = ", m
                    end if
#endif
                    e_field_nm%item(n)%item(m) = buffer_real * cmplx(0,1, kind=8)**n
                end do
            end do

            ! first line eq. (4)
            do n= n_range(1), n_range(2)
                i = n - n_range(1) + 1
                do m=-n, n
                    buffer_cmplx = cmplx(0, 1, kind=8) * e_field_nm%item(n)%item(m)
                    e_field_n_s(i) = buffer_cmplx * (a_n(i)*N_nm(n)%coordinate(m) + b_n(i)*M_nm(n)%coordinate(m))
#ifdef _DEBUG_
                    if (isnan(real(e_field_n_s(i)%rho)) .or. isnan(aimag(e_field_n_s(i)%rho)) &
                        .or. isnan(real(e_field_n_s(i)%phi)) .or. isnan(aimag(e_field_n_s(i)%phi)) &
                        .or. isnan(real(e_field_n_s(i)%theta)) .or. isnan(aimag(e_field_n_s(i)%theta)) ) then
                        print *, "get_e_field_scattered_xu: ERROR"
                        print *, "  e_field_n_s is NaN"
                        print * , "  n = ", n
                        print * , "  m = ", m
                    end if
#endif
                end do
            end do
            e_field_s%theta = cmplx(0,0,kind=8)
            e_field_s%phi = cmplx(0,0,kind=8)
            e_field_s%rho = cmplx(0,0,kind=8)
            do i=1, n_range(2)-n_range(1)+1
                e_field_s = e_field_s + e_field_n_s(i)
            end do

        end function get_e_field_scattered_xu

        ! calculates the scattering coefficients
        !
        ! HINT: The refractive index has only a real value.
        !
        ! Arguments
        ! ----
        !   x: double precision
        !       size parameter: x = k*r = 2 * PI * N * r / lambda
        !       k: wavenumber
        !       r: distance
        !       N: refractive index of the medium
        !       lambda: wave length
        !   m: double precision
        !       relative refractive index: m = N_1 / N
        !       N_1: refractive index of the particle
        !       N: refractive index of the medium
        !   mu: double precision
        !       permeability of the medium
        !   mu1: double precision
        !       permeability of the particle
        !   n: integer, dimension(2)
        !       first element: start index
        !       second element: last index
        !       HINT: first element .le. second element
        !
        ! Reference: Absorption and Scattering of Light by Small Particles, eq. (4.53)
        subroutine get_coefficients_a_b_real(x, m, mu, mu1, n, a_n, b_n)
            implicit none
            ! dummy
            double precision :: x
            double precision :: m
            double precision :: mu
            double precision :: mu1
            integer(kind=4), dimension(2) :: n

            complex(kind=8), dimension(n(2)-n(1)+1) :: a_n
            complex(kind=8), dimension(n(2)-n(1)+1) :: b_n

            ! auxiliary
            integer :: i
            complex(kind=8), dimension(n(2)-n(1)+1) :: numerator
            complex(kind=8), dimension(n(2)-n(1)+1) :: denominator

            double precision, dimension(n(2)-n(1)+1) :: j_n_x
            double precision, dimension(n(2)-n(1)+1) :: j_n_mx
            double precision, dimension(n(2)-n(1)+1) :: s_n_x
            double precision, dimension(n(2)-n(1)+1) :: s_dn_x
            double precision, dimension(n(2)-n(1)+1) :: s_n_mx
            double precision, dimension(n(2)-n(1)+1) :: s_dn_mx

            complex(kind=8), dimension(n(2)-n(1)+1) :: h_n_x
            complex(kind=8), dimension(n(2)-n(1)+1) :: xi_n_x
            complex(kind=8), dimension(n(2)-n(1)+1) :: xi_dn_x

            double precision :: mx

            integer(kind=4) :: number_of_n

            number_of_n = n(2) - n(1) + 1

            mx = m*x

            s_dn_x = lib_math_riccati_s_derivative(x, n(1), number_of_n, s_n_x)
            j_n_x = s_n_x / x

            s_dn_mx = lib_math_riccati_s_derivative(mx, n(1), number_of_n, s_n_mx) * m
            j_n_mx = s_n_mx / mx

            xi_dn_x = lib_math_riccati_xi_derivative(x, n(1), number_of_n, xi_n_x)
            h_n_x = xi_n_x / x

            numerator = cmplx(mu * m*m * j_n_mx * s_dn_x - mu1 * j_n_x * s_dn_mx, 0.0, kind=8)
            denominator = cmplx(mu * m*m * j_n_mx, 0.0, kind=8) * xi_dn_x - mu1 * h_n_x * s_dn_mx

            a_n = numerator / denominator

            numerator = cmplx(mu1 * j_n_mx * s_dn_x - mu * j_n_x * s_dn_mx, 0.0, kind=8)
            denominator = cmplx(mu1 * j_n_mx, 0.0, kind=8) * xi_dn_x - mu * h_n_x * s_dn_mx

            b_n = numerator / denominator

        end subroutine get_coefficients_a_b_real

        ! calculates the scattering coefficients
        !
        ! Arguments
        ! ----
        !   x: complex
        !       size parameter: x = k*r = 2 * PI * N * r / lambda
        !       k: wavenumber
        !       r: distance
        !       N: refractive index of the medium
        !       lambda: wave length
        !   m: complex
        !       relative refractive index: m = N_1 / N
        !       N_1: refractive index of the particle
        !       N: refractive index of the medium
        !   mu: double precision
        !       permeability of the medium
        !   mu1: double precision
        !       permeability of the particle
        !   n: integer, dimension(2)
        !       first element: start index
        !       second element: last index
        !       HINT: first element .le. second element
        !
        ! Reference: Absorption and Scattering of Light by Small Particles, eq. (4.53)
        subroutine get_coefficients_a_b_cmplx(x, m, mu, mu1, n, a_n, b_n)
            implicit none
            ! dummy
            complex(kind=8) :: x
            complex(kind=8) :: m
            double precision :: mu
            double precision :: mu1
            integer(kind=4), dimension(2) :: n

            complex(kind=8), dimension(n(2)-n(1)+1) :: a_n
            complex(kind=8), dimension(n(2)-n(1)+1) :: b_n

            ! auxiliary
            complex(kind=8), dimension(n(2)-n(1)+1) :: numerator
            complex(kind=8), dimension(n(2)-n(1)+1) :: denominator

            complex(kind=8), dimension(n(2)-n(1)+1) :: j_n_x
            complex(kind=8), dimension(n(2)-n(1)+1) :: j_n_mx
            complex(kind=8), dimension(n(2)-n(1)+1) :: s_n_x
            complex(kind=8), dimension(n(2)-n(1)+1) :: s_dn_x
            complex(kind=8), dimension(n(2)-n(1)+1) :: s_n_mx
            complex(kind=8), dimension(n(2)-n(1)+1) :: s_dn_mx

            complex(kind=8), dimension(n(2)-n(1)+1) :: h_n_x
            complex(kind=8), dimension(n(2)-n(1)+1) :: xi_n_x
            complex(kind=8), dimension(n(2)-n(1)+1) :: xi_dn_x

            complex(kind=8) :: mx

            integer(kind=4) :: number_of_n

            number_of_n = n(2) - n(1) + 1

            mx = m*x

            s_dn_x = lib_math_riccati_s_derivative(x, n(1), number_of_n, s_n_x)
            j_n_x = s_n_x / x

            s_dn_mx = lib_math_riccati_s_derivative(mx, n(1), number_of_n, s_n_mx) * m
            j_n_mx = s_n_mx / mx

            xi_dn_x = lib_math_riccati_xi_derivative(x, n(1), number_of_n, xi_n_x)
            h_n_x = xi_n_x / x

            ! --- test ---

            numerator = mu * m*m * j_n_mx * s_dn_x
            numerator = - mu1 * j_n_x * s_dn_mx

            denominator = mu * m*m * j_n_mx * xi_dn_x
            denominator = - mu1 * h_n_x * s_dn_mx

            ! ~~~ test ~~~
            numerator = mu * m*m * j_n_mx * s_dn_x - mu1 * j_n_x * s_dn_mx
            denominator = mu * m*m * j_n_mx * xi_dn_x - mu1 * h_n_x * s_dn_mx

            a_n = numerator / denominator

            numerator = mu1 * j_n_mx * s_dn_x - mu * j_n_x * s_dn_mx
            denominator = mu1 * j_n_mx * xi_dn_x - mu * h_n_x * s_dn_mx

            b_n = numerator / denominator

        end subroutine get_coefficients_a_b_cmplx

        ! calculates the scattering coefficients
        !
        ! Arguments
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
        !   n: integer, dimension(2)
        !       first element: start index
        !       second element: last index
        !       HINT: first element .le. second element
        !
        ! Reference: Light Scattering by Particles, Barber, Hill, eq. (4.18)
        subroutine get_coefficients_a_b_real_barberh(x, m, n, a_n, b_n)
            implicit none
            ! dummy
            real(kind=8) :: x
            real(kind=8) :: m
            integer(kind=4), dimension(2) :: n

            complex(kind=8), dimension(n(2)-n(1)+1) :: a_n
            complex(kind=8), dimension(n(2)-n(1)+1) :: b_n

            ! auxiliary
            integer :: i
            complex(kind=8) :: numerator
            complex(kind=8) :: denominator

            complex(kind=8), dimension(n(2)-n(1)+2) :: j_n_x
            complex(kind=8), dimension(n(2)-n(1)+2) :: h_n_x

            real(kind=8), dimension(n(2)-n(1)+1) :: An

            real(kind=8) :: mx
            real(kind=8) :: n_div_x
            real(kind=8) :: mAn
            real(kind=8) :: An_div_m
            real(kind=8) :: buffer

            integer(kind=4) :: number_of_n

            number_of_n = n(2) - n(1) + 1

            mx = m*x

            j_n_x = lib_math_bessel_spherical_first_kind(x, n(1)-1, number_of_n+1)
            h_n_x = lib_math_hankel_spherical_1(x, n(1)-1, number_of_n+1)

            An = get_An_real(x, m, n(1), number_of_n)

            do i=1, number_of_n
                n_div_x = real(n(1)+i-1, kind=8) / x
                mAn = m * An(i)
                An_div_m = An(i) / m

                buffer = An_div_m + n_div_x
                numerator = buffer * j_n_x(i+1) - j_n_x(i)
                denominator = buffer * h_n_x(i+1) - h_n_x(i)

                a_n(i) = numerator / denominator

                buffer = mAn + n_div_x
                numerator = buffer * j_n_x(i+1) - j_n_x(i)
                denominator = buffer * h_n_x(i+1) - h_n_x(i)

                b_n(i) = numerator / denominator
            end do

        end subroutine get_coefficients_a_b_real_barberh

        ! calculates the scattering coefficients
        !
        ! Arguments
        ! ----
        !   x: complex
        !       size parameter: x = k*r = 2 * PI * N * r / lambda
        !       k: wavenumber
        !       r: distance
        !       N: refractive index of the medium
        !       lambda: wave length
        !   m: complex
        !       relative refractive index: m = N_1 / N
        !       N_1: refractive index of the particle
        !       N: refractive index of the medium
        !   n: integer, dimension(2)
        !       first element: start index
        !       second element: last index
        !       HINT: first element .le. second element
        !
        ! Reference: Light Scattering by Particles, Barber, Hill, eq. (4.18)
        subroutine get_coefficients_a_b_cmplx_barberh(x, m, n, a_n, b_n)
            implicit none
            ! dummy
            complex(kind=8) :: x
            complex(kind=8) :: m
            integer(kind=4), dimension(2) :: n

            complex(kind=8), dimension(n(2)-n(1)+1) :: a_n
            complex(kind=8), dimension(n(2)-n(1)+1) :: b_n

            ! auxiliary
            integer :: i
            complex(kind=8) :: numerator
            complex(kind=8) :: denominator

            complex(kind=8), dimension(n(2)-n(1)+2) :: j_n_x
            complex(kind=8), dimension(n(2)-n(1)+2) :: h_n_x

            complex(kind=8), dimension(n(2)-n(1)+1) :: An

            complex(kind=8) :: mx
            complex(kind=8) :: n_div_x
            complex(kind=8) :: mAn
            complex(kind=8) :: An_div_m
            complex(kind=8) :: buffer

            integer(kind=4) :: number_of_n

            number_of_n = n(2) - n(1) + 1

            mx = m*x

            j_n_x = lib_math_bessel_spherical_first_kind(x, n(1)-1, number_of_n+1)
            h_n_x = lib_math_hankel_spherical_1(x, n(1)-1, number_of_n+1)

            An = get_An_cmplx(x, m, n(1), number_of_n)

            do i=1, number_of_n
                n_div_x = real(n(1)+i-1, kind=8) / x
                mAn = m * An(i)
                An_div_m = An(i) / m

                buffer = An_div_m + n_div_x
                numerator = buffer * j_n_x(i+1) - j_n_x(i)
                denominator = buffer * h_n_x(i+1) - h_n_x(i)

                a_n(i) = numerator / denominator

                buffer = mAn + n_div_x
                numerator = buffer * j_n_x(i+1) - j_n_x(i)
                denominator = buffer * h_n_x(i+1) - h_n_x(i)

                b_n(i) = numerator / denominator
            end do

        end subroutine get_coefficients_a_b_cmplx_barberh

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
            m_n_mx = max(get_n_c(abs(x)), int(abs(m_mx))) + 15

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
        !   x: complex
        !       size parameter: x = k*r = 2 * PI * N * r / lambda
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
            complex(kind=8) :: x
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
            m_n_mx = max(get_n_c(abs(x)), int(abs(m_mx))) + 15

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
        function get_n_c(x) result (rv)
            implicit none
            ! dummy
            double precision :: x

            integer :: rv

            ! auxiliary
            double precision :: dummy

            dummy = x + 4.05_8 * x**(1.0_8/3.0_8) + 2.0_8

            rv = int(ceiling(dummy))

        end function get_n_c

        function lib_mie_scattering_by_a_sphere_test_functions() result (rv)
            implicit none
            ! dummy
            integer :: rv

            ! auxiliaray
            double precision, parameter :: ground_truth_e = 10.0_8**(-6.0_8)

            rv = 0

            if (.not. test_get_coefficients_a_b_real()) then
                rv = rv + 1
            end if
!            if (.not. test_get_coefficients_a_b_cmplx()) then
!                rv = rv + 1
!            end if
            if (.not. test_get_coefficients_a_b_real_barberh()) then
                rv = rv + 1
            end if
            if (.not. test_get_coefficients_a_b_cmplx_barberh()) then
                rv = rv + 1
            end if
            if (.not. test_get_e_field_scattered()) then
                rv = rv + 1
            end if

            print *, "----lib_mie_scattering_by_a_sphere_test_functions----"
            if (rv == 0) then
                print *, "lib_mie_scattering_by_a_sphere_test_functions tests: OK"
            else
                print *, rv,"lib_mie_scattering_by_a_sphere_test_functions test(s) FAILED"
            end if
            print *, "---------------------------------------------------------"

            contains

                function test_get_coefficients_a_b_real() result (rv)
                    implicit none
                    ! dummy
                    logical :: rv

                    ! auxiliary
                    integer(kind=4) :: i
                    real(kind=8) :: buffer

                    real(kind=8) :: x
                    real(kind=8) :: m
                    double precision :: mu
                    double precision :: mu1
                    integer(kind=4), dimension(2), parameter :: n = (/1, 4/)

                    complex(kind=8), dimension(n(2)-n(1)+1) :: a_n
                    complex(kind=8), dimension(n(2)-n(1)+1) :: b_n

                    complex(kind=8), dimension(n(2)-n(1)+1) :: ground_truth_a_n
                    complex(kind=8), dimension(n(2)-n(1)+1) :: ground_truth_b_n

                    mu = 1
                    mu1 = 1

                    m = 1.5
                    x= 10

                    ground_truth_a_n(1) = cmplx(0.938111_8, +0.240954_8, kind=8)
                    ground_truth_a_n(2) = cmplx(0.962707_8, +0.189478_8, kind=8)
                    ground_truth_a_n(3) = cmplx(0.994996_8, +0.0705625_8, kind=8)
                    ground_truth_a_n(4) = cmplx(0.99737_8, -0.0512159_8, kind=8)

                    ground_truth_b_n(1) = cmplx(0.980385_8, -0.138672_8, kind=8)
                    ground_truth_b_n(2) = cmplx(0.805329_8, +0.395947_8, kind=8)
                    ground_truth_b_n(3) = cmplx(0.940931_8, -0.235754_8, kind=8)
                    ground_truth_b_n(4) = cmplx(0.999006_8, -0.0315159_8, kind=8)

                    call get_coefficients_a_b_real(x, m, mu, mu1, n, a_n, b_n)

                    rv = .true.
                    print *, "test_get_coefficients_a_b_real:"
                    do i=n(1), n(2)
                        buffer = abs(a_n(i) - ground_truth_a_n(i))
                        if (buffer .gt. ground_truth_e) then
                            print *, "  a_n: ", i , "difference: ", buffer, " : FAILED"
                            rv = .false.
                        else
                            print *, "  a_n: ", i, ": OK"
                        end if

                        buffer = abs(b_n(i) - ground_truth_b_n(i))
                        if (buffer .gt. ground_truth_e) then
                            print *, "  b_n: ", i , "difference: ", buffer, " : FAILED"
                            rv = .false.
                        else
                            print *, "  b_n: ", i, ": OK"
                        end if
                    end do

                end function test_get_coefficients_a_b_real

                function test_get_coefficients_a_b_cmplx() result (rv)
                    implicit none
                    ! dummy
                    logical :: rv

                    ! auxiliary
                    integer(kind=4) :: i
                    real(kind=8) :: buffer

                    complex(kind=8) :: x
                    complex(kind=8) :: m
                    double precision :: mu
                    double precision :: mu1
                    integer(kind=4), dimension(2), parameter :: n = (/1, 4/)

                    complex(kind=8), dimension(n(2)-n(1)+1) :: a_n
                    complex(kind=8), dimension(n(2)-n(1)+1) :: b_n

                    complex(kind=8), dimension(n(2)-n(1)+1) :: ground_truth_a_n
                    complex(kind=8), dimension(n(2)-n(1)+1) :: ground_truth_b_n

                    ! Reference: Electromagnetic scattering on spherical polydispersions,  D.Deirmendjian, p. 27
                    m = cmplx(1.28, -1.37, kind=8)
                    x= cmplx(20, 0, kind=8)
!                    ground_truth_a_n(1) = cmplx(-0.22686_8+0.5_8, -0.12863_8, kind=8)
!                    ground_truth_b_n(1) = cmplx(0.22864_8+0.5_8, 0.13377_8, kind=8)

                    ground_truth_a_n(1) = cmplx(-181.13_8, -327.306_8, kind=8)
                    ground_truth_a_n(2) = cmplx(81.2324_8, +94.3237_8, kind=8)
                    ground_truth_a_n(3) = cmplx(-51.6918_8, -32.7439, kind=8)
                    ground_truth_a_n(4) = cmplx(36.652_8, +5.67418_8, kind=8)

                    ground_truth_b_n(1) = cmplx(0.367587_8, -0.463775_8, kind=8)
                    ground_truth_b_n(2) = cmplx(0.722992_8, +0.427339_8, kind=8)
                    ground_truth_b_n(3) = cmplx(0.159304_8, -0.340386_8, kind=8)
                    ground_truth_b_n(4) = cmplx(0.947238_8, +0.177162_8, kind=8)


                    call get_coefficients_a_b_cmplx(x, m, mu, mu1, n, a_n, b_n)

                    rv = .true.
                    print *, "test_get_coefficients_a_b_cmplx:"
                    do i=n(1), n(2)
                        buffer = abs(a_n(i) - ground_truth_a_n(i))
                        if (buffer .gt. ground_truth_e) then
                            print *, "  a_n: ", i , "difference: ", buffer, " : FAILED"
                            rv = .false.
                        else
                            print *, "  a_n: ", i, ": OK"
                        end if

                        buffer = abs(b_n(i) - ground_truth_b_n(i))
                        if (buffer .gt. ground_truth_e) then
                            print *, "  b_n: ", i , "difference: ", buffer, " : FAILED"
                            rv = .false.
                        else
                            print *, "  b_n: ", i, ": OK"
                        end if
                    end do

                end function test_get_coefficients_a_b_cmplx

                function test_get_coefficients_a_b_real_barberh() result (rv)
                    implicit none
                    ! dummy
                    logical :: rv

                    ! auxiliary
                    integer(kind=4) :: i
                    real(kind=8) :: buffer

                    real(kind=8) :: x
                    real(kind=8) :: m
                    double precision :: mu
                    double precision :: mu1
                    integer(kind=4), dimension(2), parameter :: n = (/1, 4/)

                    complex(kind=8), dimension(n(2)-n(1)+1) :: a_n
                    complex(kind=8), dimension(n(2)-n(1)+1) :: b_n

                    complex(kind=8), dimension(n(2)-n(1)+1) :: ground_truth_a_n
                    complex(kind=8), dimension(n(2)-n(1)+1) :: ground_truth_b_n


                    m = 1.5_8
                    x = 10.0_8
                    ! Reference: Light Scattering by Particles: Computational Methods, Barber Hill, program S2
                    ! WARNING: same algorithm
                    ground_truth_a_n(1) = cmplx(0.825333297_8, 0.379681736_8, kind=8)
                    ground_truth_a_n(2) = cmplx(0.999948084_8, 0.0072028893_8, kind=8)
                    ground_truth_a_n(3) = cmplx(0.970794678_8, 0.168381661_8, kind=8)
                    ground_truth_a_n(4) = cmplx(0.995235085_8, -0.068863928_8, kind=8)

                    ground_truth_b_n(1) = cmplx(0.997406423_8,0.0508609638_8, kind=8)
                    ground_truth_b_n(2) = cmplx(0.885268927_8,0.318697095_8, kind=8)
                    ground_truth_b_n(3) = cmplx(0.995330393_8,-0.0681746677_8, kind=8)
                    ground_truth_b_n(4) = cmplx(0.998444855_8,-0.0394046269_8, kind=8)

                    call get_coefficients_a_b_real_barberh(x, m, n, a_n, b_n)

                    rv = .true.
                    print *, "test_get_coefficients_a_b_real_barberh:"
                    do i=n(1), n(2)
                        buffer = abs(a_n(i) - ground_truth_a_n(i))
                        if (buffer .gt. ground_truth_e) then
                            print *, "  a_n: ", i , "difference: ", buffer, " : FAILED"
                            rv = .false.
                        else
                            print *, "  a_n: ", i, ": OK"
                        end if

                        buffer = abs(b_n(i) - ground_truth_b_n(i))
                        if (buffer .gt. ground_truth_e) then
                            print *, "  b_n: ", i , "difference: ", buffer, " : FAILED"
                            rv = .false.
                        else
                            print *, "  b_n: ", i, ": OK"
                        end if
                    end do

                end function test_get_coefficients_a_b_real_barberh

                function test_get_coefficients_a_b_cmplx_barberh() result (rv)
                    implicit none
                    ! dummy
                    logical :: rv

                    ! auxiliary
                    integer(kind=4) :: i
                    real(kind=8) :: buffer

                    complex(kind=8) :: x
                    complex(kind=8) :: m
                    double precision :: mu
                    double precision :: mu1
                    integer(kind=4), dimension(2), parameter :: n = (/1, 4/)

                    complex(kind=8), dimension(n(2)-n(1)+1) :: a_n
                    complex(kind=8), dimension(n(2)-n(1)+1) :: b_n

                    complex(kind=8), dimension(n(2)-n(1)+1) :: ground_truth_a_n
                    complex(kind=8), dimension(n(2)-n(1)+1) :: ground_truth_b_n

                    m = cmplx(1.28, -1.37, kind=8)
                    x= cmplx(10, 0, kind=8)

                    ! Reference: Light Scattering by Particles: Computational Methods, Barber Hill, program S2
                    ! WARNING: same algorithm
                    ground_truth_a_n(1) = cmplx(-0.333899468_8, 0.472839594_8, kind=8)
                    ground_truth_a_n(2) = cmplx(1.1025393_8, -0.765782595_8, kind=8)
                    ground_truth_a_n(3) = cmplx(0.419178337_8, 0.996963501_8, kind=8)
                    ground_truth_a_n(4) = cmplx(-0.174542636_8, -0.790816486_8, kind=8)

                    ground_truth_b_n(1) = cmplx(1.31455839_8, -0.476593971_8, kind=8)
                    ground_truth_b_n(2) = cmplx(-0.0436286479_8, 0.753336608_8, kind=8)
                    ground_truth_b_n(3) = cmplx(0.495494395_8, -0.906942308_8, kind=8)
                    ground_truth_b_n(4) = cmplx(1.16279888_8, 0.575285077_8, kind=8)

                    call get_coefficients_a_b_cmplx_barberh(x, m, n, a_n, b_n)

                    rv = .true.
                    print *, "test_get_coefficients_a_b_cmplx_barberh:"
                    do i=n(1), n(2)
                        buffer = abs(a_n(i) - ground_truth_a_n(i))
                        if (buffer .gt. ground_truth_e) then
                            print *, "  a_n: ", i , "difference: ", buffer, " : FAILED"
                            rv = .false.
                        else
                            print *, "  a_n: ", i, ": OK"
                        end if

                        buffer = abs(b_n(i) - ground_truth_b_n(i))
                        if (buffer .gt. ground_truth_e) then
                            print *, "  b_n: ", i , "difference: ", buffer, " : FAILED"
                            rv = .false.
                        else
                            print *, "  b_n: ", i, ": OK"
                        end if
                    end do

                end function test_get_coefficients_a_b_cmplx_barberh

                function test_get_e_field_scattered() result (rv)
                    use file_io
                    implicit none
                    ! dummy
                    logical :: rv

                    ! parameter
                    double precision, parameter :: start_angle = -90
                    double precision, parameter :: stop_angle = 90
                    integer(kind=8), parameter :: number_of_values = 720

                    ! auxiliary
                    integer(kind=4) :: i
                    integer :: u
                    complex(kind=8) :: buffer_cmplx
                    character(len=25), dimension(2) :: header
                    double precision :: theta
                    double precision :: phi
                    double precision :: r
                    double precision :: e_field_0
                    double precision :: lambda
                    double precision :: r_particle
                    double precision :: n_particle
                    double precision :: n_medium

                    double precision :: k0

                    integer(kind=4), dimension(2) :: n_range

                    double precision, dimension(number_of_values) :: degree_list
                    type(spherical_coordinate_cmplx_type), dimension(number_of_values) :: e_field_s
                    type(cartesian_coordinate_cmplx_type), dimension(number_of_values) :: e_field_s_cart
                    real(kind=8), dimension(number_of_values) :: i_field_s

                    phi = 0.0
                    theta = Pi/2.0_8
                    r = 6 * unit_mu
                    r_particle = 5 * unit_mu

                    e_field_0 = 1
                    lambda = 1 * unit_mu

                    n_particle = 1.5
                    n_medium = 1

                    k0 = 2 * PI / lambda

                    n_range(1) = 1
                    n_range(2) = get_n_c(r_particle * k0 * n_particle)

                    do i=1, number_of_values
                        degree_list(i) = start_angle + (i-1) * (stop_angle - start_angle) / number_of_values
                        theta = degree_list(i) * PI / 180.0_8
                        e_field_s(i) = get_e_field_scattered_xu(theta, phi, r, e_field_0, lambda, n_medium, &
                                                                r_particle, n_particle, n_range)
                        e_field_s_cart(i) = make_cartesian(e_field_s(i), theta, phi)

                        ! calculate the intensities
                        buffer_cmplx = abs(spherical_abs(e_field_s(i)))
                        i_field_s(i) = real(buffer_cmplx * buffer_cmplx)
                    end do

                    ! write to csv
                    header(1) = "degree"
                    header(2) = "i_field_s"
                    u = 99
                    open(unit=u, file="i_field_s.csv", status='unknown')
                    rv = write_csv(u, header, degree_list, i_field_s)
                    close(u)

                    ! --- spherical components ---
                    header(1) = "degree"
                    header(2) = "e_field_s_rho"
                    open(unit=u, file="e_field_s_rho.csv", status='unknown')
                    rv = write_csv(u, header, degree_list &
                                            , real(e_field_s(:)%rho))
                    close(u)

                    header(2) = "e_field_s_theta"
                    open(unit=u, file="e_field_s_theta.csv", status='unknown')
                    rv = write_csv(u, header, degree_list &
                                            , real(e_field_s(:)%theta))
                    close(u)

                    header(2) = "e_field_s_phi"
                    open(unit=u, file="e_field_s_phi.csv", status='unknown')
                    rv = write_csv(u, header, degree_list &
                                            , real(e_field_s(:)%phi))
                    close(u)

                    ! --- cartesian components ---
                    header(1) = "degree"
                    header(2) = "e_field_s_x"
                    open(unit=u, file="e_field_s_x.csv", status='unknown')
                    rv = write_csv(u, header, degree_list &
                                            , e_field_s_cart(:)%x)
                    close(u)

                    header(2) = "e_field_s_y"
                    open(unit=u, file="e_field_s_y.csv", status='unknown')
                    rv = write_csv(u, header, degree_list &
                                            , e_field_s_cart(:)%y)
                    close(u)

                    header(2) = "e_field_s_z"
                    open(unit=u, file="e_field_s_z.csv", status='unknown')
                    rv = write_csv(u, header, degree_list &
                                            , e_field_s_cart(:)%z)
                    close(u)


                end function

                function test_get_e_field_scattered_plane_section() result (rv)
                    use file_io
                    implicit none
                    ! dummy
                    logical :: rv

                    ! parameter
                    double precision, dimension(2), parameter :: x_range = (/ -5, 5 /)
                    double precision, dimension(2), parameter :: z_range = (/ 0, 5 /)
                    real(kind=8), parameter :: step_size = 0.5

                    ! auxiliary
                    integer(kind=4) :: i
                    integer :: u
                    complex(kind=8) :: buffer_cmplx
                    character(len=25), dimension(2) :: header
                    double precision :: theta
                    double precision :: phi
                    double precision :: rho
                    double precision :: e_field_0
                    double precision :: rho_particle
                    double precision :: n_particle
                    double precision :: n_medium

                    integer(kind=4), dimension(2) :: n

!                    double precision, dimension(number_of_values) :: degree_list
!                    type(spherical_coordinate_cmplx_type), dimension(number_of_values) :: e_field_s
!                    real(kind=8), dimension(number_of_values) :: i_field_s
!
!                    phi = 0.0
!                    theta = Pi/2.0_8
!                    rho = 10
!                    rho_particle = 10
!
!                    e_field_0 = 1
!                    n_particle = 1.5
!                    n_medium = 1
!
!                    n(1) = 1
!                    n(2) = get_n_c(rho_particle)
!
!                    do i=1, number_of_values
!                        degree_list(i) = start_angle + (i-1) * (stop_angle - start_angle) / number_of_values
!!                        theta = degree_list(i) * PI / 180.0_8
!                        rho = degree_list(i)
!                        e_field_s(i) = get_e_field_scattered_xu(theta, phi, rho, e_field_0, rho_particle, n_particle, n_medium, n)
!
!                        ! calculate the intensities
!                        buffer_cmplx = abs(spherical_abs(e_field_s(i)))
!                        i_field_s(i) = real(buffer_cmplx * buffer_cmplx)
!                    end do
!
!
!
!                    ! write to csv
!                    header(1) = "degree"
!                    header(2) = "i_field_s"
!                    u = 99
!                    open(unit=u, file="i_field_s.csv", status='unknown')
!                    rv = write_csv(u, header, degree_list, i_field_s)
!                    close(u)
!
!                    header(1) = "degree"
!                    header(2) = "e_field_s_rho"
!                    open(unit=u, file="e_field_s_rho.csv", status='unknown')
!                    rv = write_csv(u, header, degree_list &
!                                            , real(e_field_s(:)%rho))
!                    close(u)
!
!                    header(2) = "e_field_s_theta"
!                    open(unit=u, file="e_field_s_theta.csv", status='unknown')
!                    rv = write_csv(u, header, degree_list &
!                                            , real(e_field_s(:)%theta))
!                    close(u)
!
!                    header(2) = "e_field_s_phi"
!                    open(unit=u, file="e_field_s_phi.csv", status='unknown')
!                    rv = write_csv(u, header, degree_list &
!                                            , real(e_field_s(:)%phi))
!                    close(u)


                end function

        end function lib_mie_scattering_by_a_sphere_test_functions

end module lib_mie_scattering_by_a_sphere

