! Scattering by a sphere
! ----
!
! - incident plane wave along the z-axis
!


module lib_mie_scattering_by_a_sphere
    use lib_math_bessel
    use lib_mie_vector_spherical_harmonics
    implicit none

    private

    contains

        function get_E_scattered() result (rv)
            implicit none
            ! dummy
            logical :: rv
            ! auxiliary

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
        function get_a_n_real(x, m, mu, mu1, n) result (a_n)
            implicit none
            ! dummy
            double precision :: x
            double precision :: m
            double precision :: mu
            double precision :: mu1
            integer(kind=4), dimension(2) :: n

            complex(kind=8), dimension(n(2)-n(1)+1) :: a_n

            ! auxiliary
            complex(kind=8), dimension(n(2)-n(1)+1) :: m_numerator
            complex(kind=8), dimension(n(2)-n(1)+1) :: m_denominator

            double precision, dimension(n(2)-n(1)+1) :: j_n_x
            double precision, dimension(n(2)-n(1)+1) :: s_n_x
            double precision, dimension(n(2)-n(1)+1) :: s_dn_x
            double precision, dimension(n(2)-n(1)+1) :: s_n_mx
            double precision, dimension(n(2)-n(1)+1) :: s_dn_mx

            double precision, dimension(n(2)-n(1)+1) :: h_n_x


        end function get_a_n_real

end module lib_mie_scattering_by_a_sphere
