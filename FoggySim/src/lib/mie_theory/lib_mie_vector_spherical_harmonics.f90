module lib_mie_vector_spherical_harmonics
    use lib_math_bessel
    use lib_math_legendre_polynomial
    use lib_math_types
    implicit none

    private

    integer(kind=1), parameter :: VECTOR_SPHERICAL_HARMONICS_COMPONENT_NUMBER_KIND = 4

    contains

        ! calculation of the vector spherical harmonic M_emn
        !
        ! Argument
        ! ----
        !   theta: double precision
        !       polar angle
        !   phi: double precision
        !       azimuthal angle
        !   rho: double precision
        !       dimensionless varibale rho = k*r
        !       k: wavenumber
        !       r: distance
        !   m: intger, dimension(2)
        !       first element: start index
        !       second element: last index
        !       HINT: first element .le. second element
        !   n: integer, dimension(2)
        !       first element: start index
        !       second element: last index
        !       HINT: first element .le. second element
        !   z_selector: integer
        !       1: spherical Bessel function first kind   j_n
        !       2: spherical Bessel function second kind  y_n
        !       3: spherical Hankel function first kind   h^(1)_n
        !       4: spherical Hankel function second kind  h^(2)_n
        !
        ! Results
        ! ----
        !   rv: complex, dimension(3)
        !       values of the spherical coordinates (rho, theta, phi)
        !
        ! LaTeX: $$ \begin{aligned} \mathbf{M}_{e m n}=& \frac{-m}{\sin \theta} \sin m \phi P_{n}^{m}(\cos \theta) z_{n}(\rho) \hat{\mathbf{e}}_{\theta} \\ &-\cos m \phi \frac{d P_{n}^{m}(\cos \theta)}{d \theta} z_{n}(\rho) \hat{\mathbf{e}}_{\phi} \end{aligned} $$
        !
        ! Reference: Absorption and Scattering of Light by Small Particles, eq. 4.17
        subroutine lib_mie_vector_spherical_harmonics_components(theta, phi, rho, m, n, z_selector, &
                                                      M_emn, M_omn, N_emn, N_omn)
            implicit none
            ! dummy
            double precision, intent(in) :: theta
            double precision, intent(in) :: phi
            double precision, intent(in) :: rho
            integer(kind=VECTOR_SPHERICAL_HARMONICS_COMPONENT_NUMBER_KIND), dimension(2) :: m
            integer(kind=VECTOR_SPHERICAL_HARMONICS_COMPONENT_NUMBER_KIND), dimension(2) :: n
            integer(kind=1) :: z_selector

            type(list_spherical_coordinate_type), dimension(m(1):m(2)), intent(inout) :: M_emn
            type(list_spherical_coordinate_type), dimension(m(1):m(2)), intent(inout) :: M_omn
            type(list_spherical_coordinate_type), dimension(m(1):m(2)), intent(inout) :: N_emn
            type(list_spherical_coordinate_type), dimension(m(1):m(2)), intent(inout) :: N_omn

!            complex(kind=8), dimension(m(1):m(2), n(1):n(2), 3), intent(inout) :: M_emn
!            complex(kind=8), dimension(m(1):m(2), n(1):n(2), 3), intent(inout) :: M_omn
!            complex(kind=8), dimension(m(1):m(2), n(1):n(2), 3), intent(inout) :: N_emn
!            complex(kind=8), dimension(m(1):m(2), n(1):n(2), 3), intent(inout) :: N_omn

            ! auxiliary
            integer(kind=VECTOR_SPHERICAL_HARMONICS_COMPONENT_NUMBER_KIND) :: i
            double precision, dimension(n(1):n(2)) :: z_real
            complex(kind=8), dimension(n(1):n(2)) :: z_cmplx
            double precision, dimension(n(1):n(2)) :: z_d_real ! deriviative
            complex(kind=8), dimension(n(1):n(2)) :: z_d_cmplx ! deriviative
            integer(kind=VECTOR_SPHERICAL_HARMONICS_COMPONENT_NUMBER_KIND) :: number_of_members_m
            integer(kind=VECTOR_SPHERICAL_HARMONICS_COMPONENT_NUMBER_KIND) :: number_of_members_n
            double precision :: cos_theta
            double precision :: sin_theta
            double precision, dimension(m(1):m(2)) :: cos_m_phi
            double precision, dimension(m(1):m(2)) :: sin_m_phi

            number_of_members_m = m(2) - m(1) + 1
            number_of_members_n = n(2) - n(1) + 1

!            allocate (M_emn(:)%coordinate(number_of_members_n))
!            allocate (M_omn(:)%coordinate(number_of_members_n))
!            allocate (N_emn(:)%coordinate(number_of_members_n))
!            allocate (N_omn(:)%coordinate(number_of_members_n))

            ! --- pre-calculation ---
            cos_theta = cos(theta)
            sin_theta = sin(theta)

            do i=m(1), m(2)
                cos_m_phi(i) = cos(i*phi)
                sin_m_phi(i) = sin(i*phi)
            end do

            select case (z_selector)
                case(1)
                    ! spherical Bessel function first kind j_n
                    z_real(n(1):n(2)) = lib_math_bessel_spherical_first_kind(rho, n(1), number_of_members_n)
                case(2)
                    ! spherical Bessel function second kind y_n
                    z_real(n(1):n(2)) = lib_math_bessel_spherical_second_kind(rho, n(1), number_of_members_n)
                case(3)

                case(4)

                case default
                    z_real = 0
                    z_cmplx = cmplx(0,0)
                    print*, "lib_mie_vector_spherical_harmonics_M_emn: ERROR"
                    print*, "  undefined z_selector value: ", z_selector
            end select

        end subroutine lib_mie_vector_spherical_harmonics_components


end module lib_mie_vector_spherical_harmonics
