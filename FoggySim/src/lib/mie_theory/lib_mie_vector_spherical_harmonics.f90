module lib_mie_vector_spherical_harmonics
    use lib_math_bessel
    use lib_math_legendre
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
        subroutine lib_mie_vector_spherical_harmonics_components_real(theta, phi, rho, m, n, z_selector, &
                                                      M_emn, M_omn, N_emn, N_omn)
            implicit none
            ! dummy
            double precision, intent(in) :: theta
            double precision, intent(in) :: phi
            double precision, intent(in) :: rho
            integer(kind=VECTOR_SPHERICAL_HARMONICS_COMPONENT_NUMBER_KIND), dimension(2) :: m
            integer(kind=VECTOR_SPHERICAL_HARMONICS_COMPONENT_NUMBER_KIND), dimension(2) :: n
            integer(kind=1) :: z_selector

            type(list_spherical_coordinate_cmplx_type), dimension(m(2)-m(1)+1), intent(inout) :: M_emn
            type(list_spherical_coordinate_cmplx_type), dimension(m(2)-m(1)+1), intent(inout) :: M_omn
            type(list_spherical_coordinate_cmplx_type), dimension(m(2)-m(1)+1), intent(inout) :: N_emn
            type(list_spherical_coordinate_cmplx_type), dimension(m(2)-m(1)+1), intent(inout) :: N_omn

            ! auxiliary
            integer(kind=VECTOR_SPHERICAL_HARMONICS_COMPONENT_NUMBER_KIND) :: i
            integer(kind=VECTOR_SPHERICAL_HARMONICS_COMPONENT_NUMBER_KIND) :: number_of_members_m
            integer(kind=VECTOR_SPHERICAL_HARMONICS_COMPONENT_NUMBER_KIND) :: number_of_members_n

            double precision, dimension(m(2)-m(1)+1, n(2)-n(1)+1) :: p_mn
            double precision, dimension(m(2)-m(1)+1, n(2)-n(1)+1) :: p_dmn
            double precision, dimension(n(2)-n(1)+1) :: buffer_p_n
            double precision, dimension(n(2)-n(1)+1) :: buffer_p_dn

            double precision, dimension(n(2)-n(1)+1) :: z_real
            complex(kind=8), dimension(n(2)-n(1)+1) :: z_cmplx
            double precision, dimension(n(2)-n(1)+1) :: z_d_real ! deriviative
            complex(kind=8), dimension(n(2)-n(1)+1) :: z_d_cmplx ! deriviative

            double precision :: cos_theta
            double precision :: sin_theta
            double precision, dimension(n(2)-n(1)+1) :: cos_m_phi
            double precision, dimension(n(2)-n(1)+1) :: sin_m_phi


            number_of_members_m = m(2) - m(1) + 1
            number_of_members_n = n(2) - n(1) + 1

            ! --- init ---
            do i=1, number_of_members_m
                allocate (M_emn(i)%coordinate(number_of_members_n))
                allocate (M_omn(i)%coordinate(number_of_members_n))
                allocate (N_emn(i)%coordinate(number_of_members_n))
                allocate (N_omn(i)%coordinate(number_of_members_n))
            end do

            ! --- pre-calculation ---
            cos_theta = cos(theta)
            sin_theta = sin(theta)

            do i=1, number_of_members_m
                cos_m_phi(i) = cos(i*phi)
                sin_m_phi(i) = sin(i*phi)
            end do

            select case (z_selector)
                case(1)
                    ! spherical Bessel function first kind j_n
                    ! internal: calculation with Riccati-Bessel functions: S_n
                    z_d_real = lib_math_riccati_s_derivative(rho, n(1), number_of_members_n, z_real)
                case(2)
                    ! spherical Bessel function second kind y_n
                    ! internal: calculation with Riccati-Bessel functions: C_n
                    z_d_real = lib_math_riccati_c_derivative(rho, n(1), number_of_members_n, z_real)
                case(3)
                    ! spherical Hankel function first kind   h^(1)_n
                    ! internal: calculation with Riccati-Bessel functions: Xi_n
                    z_d_cmplx = lib_math_riccati_xi_derivative(rho, n(1), number_of_members_n, z_cmplx)
                case(4)
                    ! spherical Hankel function first kind   h^(2)_n
                    ! internal: calculation with Riccati-Bessel functions: Zeta_n
                    z_d_cmplx = lib_math_riccati_zeta_derivative(rho, n(1), number_of_members_n, z_cmplx)
                case default
                    z_real = 0
                    z_cmplx = cmplx(0,0)
                    z_d_real = 0
                    z_d_cmplx = cmplx(0,0)
                    print*, "lib_mie_vector_spherical_harmonics_M_emn: ERROR"
                    print*, "  undefined z_selector value: ", z_selector
            end select

            do i=1, number_of_members_m
                call lib_math_associated_legendre_polynomial(cos_theta, m, n, buffer_p_n,  buffer_p_dn)
                p_mn(i, :) = buffer_p_n
                p_dnm(i, :) = buffer_pdn
            end do


        end subroutine lib_mie_vector_spherical_harmonics_components_real


end module lib_mie_vector_spherical_harmonics
