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
                                                      M_emn, M_omn, N_emn, N_omn, &
                                                      not_calc_Memn, not_calc_Momn, not_calc_Nemn, not_calc_Nomn)
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

            logical, optional :: not_calc_Memn
            logical, optional :: not_calc_Momn
            logical, optional :: not_calc_Nemn
            logical, optional :: not_calc_Nomn

            ! auxiliary
            integer(kind=VECTOR_SPHERICAL_HARMONICS_COMPONENT_NUMBER_KIND) :: i
            integer(kind=VECTOR_SPHERICAL_HARMONICS_COMPONENT_NUMBER_KIND) :: ii
            integer(kind=VECTOR_SPHERICAL_HARMONICS_COMPONENT_NUMBER_KIND) :: number_of_members_m
            integer(kind=VECTOR_SPHERICAL_HARMONICS_COMPONENT_NUMBER_KIND) :: number_of_members_n

            double precision, dimension(m(2)-m(1)+1, n(2)-n(1)+1) :: p_mn
            double precision, dimension(m(2)-m(1)+1, n(2)-n(1)+1) :: p_dmn
            double precision, dimension(n(2)-n(1)+1) :: buffer_p_n
            double precision, dimension(n(2)-n(1)+1) :: buffer_p_dn

            double precision, dimension(n(2)-n(1)+1) :: z_n_real
            complex(kind=8), dimension(n(2)-n(1)+1) :: z_n_cmplx
            double precision, dimension(n(2)-n(1)+1) :: z_d_real ! deriviative
            complex(kind=8), dimension(n(2)-n(1)+1) :: z_d_cmplx ! deriviative

            double precision, dimension(n(2)-n(1)+1) :: z_divided_by_rho_real
            complex(kind=8), dimension(n(2)-n(1)+1) :: z_divided_by_rho_cmplx

            ! Riccati-Bessel
            double precision, dimension(n(2)-n(1)+1) :: r_real
            complex(kind=8), dimension(n(2)-n(1)+1) :: r_cmplx
            double precision, dimension(n(2)-n(1)+1) :: r_d_real ! deriviative
            complex(kind=8), dimension(n(2)-n(1)+1) :: r_d_cmplx ! deriviative

            double precision :: cos_theta
            double precision :: sin_theta
            double precision :: minus_sin_theta
            double precision, dimension(m(2)-m(1)+1) :: m_divided_by_sin_theta

            double precision, dimension(n(2)-n(1)+1) :: cos_m_phi
            double precision, dimension(n(2)-n(1)+1) :: sin_m_phi


            logical :: m_not_calc_Memn
            logical :: m_not_calc_Momn
            logical :: m_not_calc_Nemn
            logical :: m_not_calc_Nomn

            double precision :: buffer_real
            complex(kind=8) :: buffer_cmplx


            number_of_members_m = m(2) - m(1) + 1
            number_of_members_n = n(2) - n(1) + 1

            ! --- standard value ---
            m_not_calc_Memn = .false.
            m_not_calc_Momn = .false.
            m_not_calc_Nemn = .false.
            m_not_calc_Nomn = .false.

            if (present(not_calc_Memn)) then
                m_not_calc_Memn = not_calc_Memn
            end if

            if (present(not_calc_Momn)) then
                m_not_calc_Momn = not_calc_Momn
            end if

            if (present(not_calc_Nemn)) then
                m_not_calc_Nemn = not_calc_Nemn
            end if
            if (present(not_calc_Nomn)) then
                m_not_calc_Nomn = not_calc_Nomn
            end if

            ! --- init ---
            if (.not. m_not_calc_Memn) then
                do i=1, number_of_members_m
                    allocate (M_emn(i)%coordinate(number_of_members_n))
                end do
            end if

            if (.not. m_not_calc_Momn) then
                do i=1, number_of_members_m
                    allocate (M_omn(i)%coordinate(number_of_members_n))
                end do
            end if

            if (.not. m_not_calc_Nemn) then
                do i=1, number_of_members_m
                    allocate (N_emn(i)%coordinate(number_of_members_n))
                end do
            end if

            if (.not. m_not_calc_Nemn) then
                do i=1, number_of_members_m
                    allocate (N_omn(i)%coordinate(number_of_members_n))
                end do
            end if


            ! --- pre-calculation ---
            cos_theta = cos(theta)
            sin_theta = sin(theta)
            minus_sin_theta = -sin(theta)


            do i=1, number_of_members_m
                cos_m_phi(i) = cos(i*phi)
                sin_m_phi(i) = sin(i*phi)
                m_divided_by_sin_theta(i) = i / sin_theta
            end do

            select case (z_selector)
                case(1)
                    ! spherical Bessel function first kind j_n
                    ! internal: calculation with Riccati-Bessel functions: S_n
                    r_d_real = lib_math_riccati_s_derivative(rho, n(1), number_of_members_n, r_real)
                    z_n_real = r_real / rho
                case(2)
                    ! spherical Bessel function second kind y_n
                    ! internal: calculation with Riccati-Bessel functions: C_n
                    r_d_real = lib_math_riccati_c_derivative(rho, n(1), number_of_members_n, r_real)
                    z_n_real = r_real / rho
                case(3)
                    ! spherical Hankel function first kind   h^(1)_n
                    ! internal: calculation with Riccati-Bessel functions: Xi_n
                    r_d_cmplx = lib_math_riccati_xi_derivative(rho, n(1), number_of_members_n, r_cmplx)
                    z_n_cmplx = r_cmplx / rho
                case(4)
                    ! spherical Hankel function first kind   h^(2)_n
                    ! internal: calculation with Riccati-Bessel functions: Zeta_n
                    r_d_cmplx = lib_math_riccati_zeta_derivative(rho, n(1), number_of_members_n, r_cmplx)
                    z_n_cmplx = r_cmplx / rho
                case default
                    z_n_real = 0
                    z_n_cmplx = cmplx(0,0)
                    z_d_real = 0
                    z_d_cmplx = cmplx(0,0)

                    r_real = 0
                    r_cmplx = cmplx(0,0)
                    r_d_real = 0
                    r_d_cmplx = cmplx(0,0)
                    print*, "lib_mie_vector_spherical_harmonics_M_emn: ERROR"
                    print*, "  undefined z_selector value: ", z_selector
                    return
            end select

            do i=1, number_of_members_m
                call lib_math_associated_legendre_polynomial(cos_theta, i, n(2), buffer_p_n,  buffer_p_dn, .false.)
                p_mn(i, :) = buffer_p_n(n(1):n(2))
                ! final step of the calculation of >>> derivation[LegendreP[l, m, cos (x)], x]
                p_dmn(i, :) = minus_sin_theta * buffer_p_dn(n(1):n(2))
            end do


            ! --- calculations of the components M and N ---
            ! M_emn
            if (.not. m_not_calc_Memn) then
                select case (z_selector)
                    case (1,2)
                        ! z = [j_n, y_n] ==> z: real
                        do i=1, number_of_members_m
                            do ii=1, number_of_members_n
                                buffer_real = - m_divided_by_sin_theta(i) * sin_m_phi(i) * p_mn(i, ii) * z_n_real(ii)
                                M_emn(i)%coordinate(ii)%theta = cmplx(buffer_real, 0, kind=8)

                                buffer_real = - cos_m_phi(i) * p_dmn(i, ii) * z_n_real(ii)
                                M_emn(i)%coordinate(ii)%phi = cmplx(buffer_real, 0, kind=8)

                                M_emn(i)%coordinate(ii)%rho = cmplx(0,0, kind=8)
                            end do
                        end do
                    case (3,4)
                        ! z = [h^(1)_n, h^(2)_n] ==> z: complex
                        do i=1, number_of_members_m
                            do ii=1, number_of_members_n
                                buffer_real = - m_divided_by_sin_theta(i) * sin_m_phi(i) * p_mn(i, ii)
                                buffer_cmplx = buffer_real * z_n_cmplx(ii)
                                M_emn(i)%coordinate(ii)%theta = buffer_cmplx

                                buffer_real = - cos_m_phi(i) * p_dmn(i, ii)
                                buffer_cmplx = buffer_real * z_n_cmplx(ii)
                                M_emn(i)%coordinate(ii)%phi = buffer_cmplx

                                M_emn(i)%coordinate(ii)%rho = cmplx(0,0, kind=8)
                            end do
                        end do
                end select
            end if

            ! M_omn
            if (.not. m_not_calc_Momn) then
                select case (z_selector)
                    case (1,2)
                        ! z = [j_n, y_n] ==> z: real
                        do i=1, number_of_members_m
                            do ii=1, number_of_members_n
                                buffer_real = m_divided_by_sin_theta(i) * cos_m_phi(i) * p_mn(i, ii) * z_n_real(ii)
                                M_omn(i)%coordinate(ii)%theta = cmplx(buffer_real, 0, kind=8)

                                buffer_real = - sin_m_phi(i) * p_dmn(i, ii) * z_n_real(ii)
                                M_omn(i)%coordinate(ii)%phi = cmplx(buffer_real, 0, kind=8)

                                M_omn(i)%coordinate(ii)%rho = cmplx(0,0, kind=8)
                            end do
                        end do
                    case (3,4)
                        ! z = [h^(1)_n, h^(2)_n] ==> z: complex
                        do i=1, number_of_members_m
                            do ii=1, number_of_members_n
                                buffer_real = m_divided_by_sin_theta(i) * cos_m_phi(i) * p_mn(i, ii)
                                buffer_cmplx = buffer_real * z_n_cmplx(ii)
                                M_omn(i)%coordinate(ii)%theta = buffer_cmplx

                                buffer_real = - sin_m_phi(i) * p_dmn(i, ii)
                                buffer_cmplx = buffer_real * z_n_cmplx(ii)
                                M_omn(i)%coordinate(ii)%phi = buffer_cmplx

                                M_omn(i)%coordinate(ii)%rho = cmplx(0,0, kind=8)
                            end do
                        end do
                end select
            end if

            ! N_emn
            if (.not. m_not_calc_Nemn) then
                select case (z_selector)
                    case (1,2)
                        ! z = [j_n, y_n] ==> z: real
                        do i=1, number_of_members_m
                            do ii=1, number_of_members_n
                                buffer_real = z_divided_by_rho_real(ii) * cos_m_phi(i) * ii*(ii+1) * p_mn(i, ii)
                                N_emn(i)%coordinate(ii)%rho = cmplx(buffer_real, 0, kind=8)

                                buffer_real = cos_m_phi(i) * p_dmn(i, ii) / rho * r_d_real(ii)
                                N_emn(i)%coordinate(ii)%theta = cmplx(buffer_real, 0, kind=8)

                                buffer_real = - i * sin_m_phi(i) * p_mn(i, ii) / (sin_theta * rho) * r_d_real(ii)
                                N_emn(i)%coordinate(ii)%phi = cmplx(buffer_real,0, kind=8)
                            end do
                        end do
                    case (3,4)
                        ! z = [h^(1)_n, h^(2)_n] ==> z: complex
                        do i=1, number_of_members_m
                            do ii=1, number_of_members_n
                                buffer_real = cos_m_phi(i) * ii*(ii+1) * p_mn(i, ii)
                                buffer_cmplx = z_divided_by_rho_cmplx(ii) * buffer_real
                                N_emn(i)%coordinate(ii)%rho = cmplx(buffer_real, 0, kind=8)

                                buffer_real = cos_m_phi(i) * p_dmn(i, ii) / rho
                                buffer_cmplx = buffer_real * r_d_cmplx(ii)
                                N_emn(i)%coordinate(ii)%theta = buffer_cmplx

                                buffer_real = - i * sin_m_phi(i) * p_mn(i, ii) / (sin_theta * rho)
                                buffer_cmplx = buffer_real * r_d_cmplx(ii)
                                N_emn(i)%coordinate(ii)%phi = cmplx(buffer_real,0, kind=8)
                            end do
                        end do
                end select
            end if

            ! N_omn
            if (.not. m_not_calc_Nomn) then
                select case (z_selector)
                    case (1,2)
                        ! z = [j_n, y_n] ==> z: real
                        do i=1, number_of_members_m
                            do ii=1, number_of_members_n
                                buffer_real = z_divided_by_rho_real(ii) * sin_m_phi(i) * ii*(ii+1) * p_mn(i, ii)
                                N_omn(i)%coordinate(ii)%rho = cmplx(buffer_real, 0, kind=8)

                                buffer_real = sin_m_phi(i) * p_dmn(i, ii) / rho * r_d_real(ii)
                                N_omn(i)%coordinate(ii)%theta = cmplx(buffer_real, 0, kind=8)

                                buffer_real = i * cos_m_phi(i) * p_mn(i, ii) / (sin_theta * rho) * r_d_real(ii)
                                N_omn(i)%coordinate(ii)%phi = cmplx(buffer_real,0, kind=8)
                            end do
                        end do
                    case (3,4)
                        ! z = [h^(1)_n, h^(2)_n] ==> z: complex
                        do i=1, number_of_members_m
                            do ii=1, number_of_members_n
                                buffer_real = sin_m_phi(i) * ii*(ii+1) * p_mn(i, ii)
                                buffer_cmplx = z_divided_by_rho_cmplx(ii) * buffer_real
                                N_omn(i)%coordinate(ii)%rho = cmplx(buffer_real, 0, kind=8)

                                buffer_real = sin_m_phi(i) * p_dmn(i, ii) / rho
                                buffer_cmplx = buffer_real * r_d_cmplx(ii)
                                N_omn(i)%coordinate(ii)%theta = buffer_cmplx

                                buffer_real = i * cos_m_phi(i) * p_mn(i, ii) / (sin_theta * rho)
                                buffer_cmplx = buffer_real * r_d_cmplx(ii)
                                N_omn(i)%coordinate(ii)%phi = buffer_cmplx
                            end do
                        end do
                end select
            end if

        end subroutine lib_mie_vector_spherical_harmonics_components_real


end module lib_mie_vector_spherical_harmonics
