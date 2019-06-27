! Scattering by a sphere
! ----
!
! - incident plane wave along the z-axis
!


module lib_mie_scattering_by_a_sphere
    use lib_math_bessel
    use lib_math_type
    use lib_mie_vector_spherical_harmonics
    implicit none

    private

    contains

        function get_e_field_scattered(theta, phi, rho, e_field_0, n_particle, n_medium) result (e_field_s)
            implicit none
            ! dummy
            double precision, intent(in) :: theta
            double precision, intent(in) :: phi
            double precision, intent(in) :: rho
            double precision, intent(in) :: e_field_0
            double precision, intent(in) :: n_particle
            double precision, intent(in) :: n_medium

            type(list_spherical_coordinate_cmplx_type) :: e_field_s

            ! parameter
            integer(kind=4), dimension(2), parameter :: n = (/1, 10/)

            ! auxiliary
            integer(kind=4) :: i
            double precision :: mu
            double precision :: mu1

            complex(kind=8), dimension(n(2)-n(1)+1) :: a_n
            complex(kind=8), dimension(n(2)-n(1)+1) :: b_n


            integer(kind=4), dimension(2) :: m
            integer(kind=1) :: z_selector

            type(list_spherical_coordinate_cmplx_type), dimension(:), allocatable :: M_emn
            type(list_spherical_coordinate_cmplx_type), dimension(:), allocatable :: M_omn
            type(list_spherical_coordinate_cmplx_type), dimension(:), allocatable :: N_emn
            type(list_spherical_coordinate_cmplx_type), dimension(:), allocatable :: N_omn

            complex(kind=8), dimension(n(2)-n(1)+1) :: e_field_n
            type(list_spherical_coordinate_cmplx_type), dimension(n(2)-n(1)+1) :: e_field_n_s

            m = (/1,1/)

            mu = 1
            mu1 = 1

            call get_coefficients_a_b_real(rho, n_particle/n_medium, mu, mu1, n, a_n, b_n)

            z_selector = 3

            call lib_mie_vector_spherical_harmonics_components(theta, phi, rho, m, n, z_selector, &
                                                               M_emn, M_omn, N_emn, N_omn, &
                                                               not_calc_Memn=.true., &
                                                               not_calc_Momn=.false., &
                                                               not_calc_Nemn=.false., &
                                                               not_calc_Nomn=.true.)
!            do i=n(1), n(2)
!                e_field_n(i) = cmplx(e_field_0 * (2*i+1)/(n*(n+1)), 0, kind=8)
!                e_field_n(i) = e_field_n(i) * (cmplx(0,1, kind=8)**i)
!            end do
!
!            do i=n(1), n(2)
!                e_field_n_s(i) = e_field_n(i) * (cmplx(0,1, kind=8) * a_n(i) * N_emn(m(1),i) - b_n(i) * M_omn(m(1), i))
!            end do
!
!            e_field_s = sum(e_field_n_s)

        end function

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
        ! Reference: Absorption and Scattering of Light by Small Particles, eq. 4.52, 4.53
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

            s_dn_mx = lib_math_riccati_s_derivative(mx, n(1), number_of_n, s_n_mx)
            j_n_mx = s_n_mx / mx

            xi_dn_x = lib_math_riccati_xi_derivative(x, n(1), number_of_n, xi_n_x)
            h_n_x = xi_n_x / x


            numerator = cmplx(mu * m*m * j_n_mx * s_dn_x - mu1 * j_n_x * s_dn_mx, 0, kind=8)
            denominator = mu * m*m * j_n_mx * xi_dn_x - mu1 * h_n_x * s_dn_mx

            a_n = numerator / denominator

            numerator = cmplx(mu1 * j_n_mx * s_dn_x - mu * j_n_x * s_dn_mx, 0, kind=8)
            denominator = mu1 * j_n_mx * xi_dn_x - mu * h_n_x * s_dn_mx

            b_n = numerator / denominator

        end subroutine get_coefficients_a_b_real

end module lib_mie_scattering_by_a_sphere
