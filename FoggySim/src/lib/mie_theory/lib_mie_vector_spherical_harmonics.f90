module lib_mie_vector_spherical_harmonics
    use libmath
    implicit none

    private

    ! --- public ---
    public :: lib_mie_vector_spherical_harmonics_components

    public :: lib_mie_vector_spherical_harmonics_components_real_xu

    interface lib_mie_vector_spherical_harmonics_components
        module procedure lib_mie_vector_spherical_harmonics_components_real
        module procedure lib_mie_vector_spherical_harmonics_components_cmplx
    end interface

    public :: lib_mie_vector_spherical_harmonics_test_functions

    ! --- parameter ---
    integer(kind=1), parameter :: VECTOR_SPHERICAL_HARMONICS_COMPONENT_NUMBER_KIND = 4

    contains

        ! calculation of the components of the vector spherical harmonic
        !
        ! Argument
        ! ----
        !   theta: double precision
        !       polar angle
        !   phi: double precision
        !       azimuthal angle
        !   rho: double precision
        !       dimensionless varibale rho = k*r > 0
        !       k: wavenumber
        !       r: distance
        !   m: integer, dimension(2)
        !       first element: start index
        !       second element: last index
        !       CONDITION:
        !           - first element .le. second element
        !           - 0 <= m <= n
        !   n: integer, dimension(2)
        !       first element: start index
        !       second element: last index
        !       CONDITION:
        !           - first element .le. second element
        !           - 0 <= n
        !
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
        ! LaTeX: $$ \mathbf{M}_{m n}^{(J)}=\left[\mathbf{i}_{\theta} i \pi_{m n}(\cos \theta)-\mathbf{i}_{\phi} \tau_{m n}(\cos \theta)\right] z_{n}^{(J)} )(k r) \exp (i m \phi) $$
        !        $$ \mathbf{N}_{m n}^{(J)}=\mathbf{i}_r n\left(n+1 ) P_{n}^{m} (\cos \theta\right) \frac{z_{n}^{(J)}(k r)}{k r} \exp(i m \phi)
        !           + \left[\mathbf{i}_{\theta} \tau_{m n}(\cos \theta)+\mathbf{i}_{\phi} i \pi_{m n}(\cos \theta)\right]
        !           \times \frac{1}{k r} \frac{d}{d r}\left[r z_{n}^{(J)}(k r)\right] \exp (i m \phi) $$
        !
        ! Reference: Electromagnetic scattering by an aggregate of spheres, Yu-lin Xu, eq. 2
        subroutine lib_mie_vector_spherical_harmonics_components_real_xu(theta, phi, rho, m, n, z_selector, &
                                                      M_mn, N_mn)
            implicit none
            ! dummy
            double precision, intent(in) :: theta
            double precision, intent(in) :: phi
            double precision, intent(in) :: rho
            integer(kind=VECTOR_SPHERICAL_HARMONICS_COMPONENT_NUMBER_KIND), dimension(2) :: m
            integer(kind=VECTOR_SPHERICAL_HARMONICS_COMPONENT_NUMBER_KIND), dimension(2) :: n
            integer(kind=1) :: z_selector

            type(list_spherical_coordinate_cmplx_type), dimension(:), allocatable, intent(inout) :: M_mn
            type(list_spherical_coordinate_cmplx_type), dimension(:), allocatable, intent(inout) :: N_mn

            ! auxiliary
            integer(kind=VECTOR_SPHERICAL_HARMONICS_COMPONENT_NUMBER_KIND) :: i
            integer(kind=VECTOR_SPHERICAL_HARMONICS_COMPONENT_NUMBER_KIND) :: ii
            integer(kind=VECTOR_SPHERICAL_HARMONICS_COMPONENT_NUMBER_KIND) :: number_of_members_m
            integer(kind=VECTOR_SPHERICAL_HARMONICS_COMPONENT_NUMBER_KIND) :: number_of_members_n

            double precision, dimension(0:n(2), -n(2):n(2)) :: pi_nm
            double precision, dimension(0:n(2), -n(2):n(2)) :: tau_nm

            double precision, dimension(0:n(2), -n(2):n(2)) :: p_nm
            double precision, dimension(0:n(2), -n(2):n(2)) :: p_d_nm

            double precision, dimension(n(2)-n(1)+1) :: buffer_p_n
            double precision, dimension(n(2)-n(1)+1) :: buffer_p_d_n

            double precision, dimension(2, n(2)-n(1)+1) :: buffer_p_n_m_neg
            double precision, dimension(2, n(2)-n(1)+1) :: buffer_p_d_n_m_neg

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

            double precision :: buffer_real
            complex(kind=8) :: buffer_cmplx

            complex(kind=8), dimension(m(2)-m(1)+1) :: exp_i_m_phi
            double precision :: sin_theta


            number_of_members_m = m(2) - m(1) + 1
            number_of_members_n = n(2) - n(1) + 1


            ! --- init ---
            allocate(M_mn(m(2)-m(1)+1))
            do i=1, number_of_members_m
                allocate (M_mn(i)%coordinate(number_of_members_n))
            end do

            allocate(N_mn(m(2)-m(1)+1))
            do i=1, number_of_members_m
                allocate (N_mn(i)%coordinate(number_of_members_n))
            end do


            ! --- pre-calculation ---
            sin_theta = sin(theta)

            do i=1, number_of_members_m
                buffer_real = i+m(1)-1
                exp_i_m_phi(i)= cmplx(cos(buffer_real * phi), sin(buffer_real * phi), kind=8)
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


            call lib_math_associated_legendre_polynomial_theta(theta, n(2), pi_nm, tau_nm)

            call lib_math_associated_legendre_polynomial(sin_theta, 0, n(1), number_of_members_n, buffer_p_n, buffer_p_d_n, .false.)
            p_nm(:, 0) = buffer_p_n
            do i=1, m(2)
                call lib_math_associated_legendre_polynomial_with_negative_m(sin_theta, i, n(1), number_of_members_n, &
                                                                             buffer_p_n_m_neg, buffer_p_d_n_m_neg, .false.)
                p_nm(:, i) = buffer_p_n_m_neg(2, :)
                p_nm(:,-i) = buffer_p_n_m_neg(1, :)
            end do

            ! --- calculations of the components M and N ---
            ! M_mn
            ! first line eq. (2)
            select case (z_selector)
                case (1,2)
                    ! z = [j_n, y_n] ==> z: real
                    do i=1, number_of_members_m
                        do ii=1, number_of_members_n
                            buffer_real = pi_nm(n(1)+ii-1, m(1)+i-1)
                            M_mn(i)%coordinate(ii)%theta = cmplx(0, buffer_real, kind=8)

                            buffer_real = -tau_nm(n(1)+ii-1, m(1)+i-1) * z_n_real(ii)
                            buffer_cmplx = buffer_real * exp_i_m_phi(i)
                            M_mn(i)%coordinate(ii)%phi = buffer_cmplx

                            M_mn(i)%coordinate(ii)%rho = cmplx(0,0, kind=8)
                        end do
                    end do
                case (3,4)
                    ! z = [h^(1)_n, h^(2)_n] ==> z: complex
                    do i=1, number_of_members_m
                        do ii=1, number_of_members_n
                            buffer_real = pi_nm(n(1)+ii-1, m(1)+i-1)
                            M_mn(i)%coordinate(ii)%theta = cmplx(0, buffer_real, kind=8)

                            buffer_real = -tau_nm(n(1)+ii-1, m(1)+i-1)
                            buffer_cmplx = buffer_real * z_n_cmplx(ii) * exp_i_m_phi(i)
                            M_mn(i)%coordinate(ii)%phi = buffer_cmplx

                            M_mn(i)%coordinate(ii)%rho = cmplx(0,0, kind=8)
                        end do
                    end do
            end select

            ! N_mn
            ! second line eq. (2)
            select case (z_selector)
                case (1,2)
                    ! z = [j_n, y_n] ==> z: real
                    do i=1, number_of_members_m
                        do ii=1, number_of_members_n
                            ! n(n+1)
                            buffer_real = n(1)+i-1
                            buffer_real = buffer_real * (buffer_real + 1)

                            buffer_real = buffer_real * p_nm(n(1)+ii-1, m(1)+i-1) * z_n_real(ii) / rho
                            buffer_cmplx = buffer_real * exp_i_m_phi(i)
                            N_mn(i)%coordinate(ii)%rho = buffer_cmplx

                            buffer_cmplx = z_d_real(ii) / rho * exp_i_m_phi(i)

                            N_mn(i)%coordinate(ii)%theta = tau_nm(n(1)+ii-1, m(1)+i-1) * buffer_cmplx

                            N_mn(i)%coordinate(ii)%phi = cmplx(0, pi_nm(n(1)+ii-1, m(1)+i-1), kind=8) * buffer_cmplx
                        end do
                    end do
                case (3,4)
                    ! z = [h^(1)_n, h^(2)_n] ==> z: complex
                    do i=1, number_of_members_m
                        do ii=1, number_of_members_n
                            ! n(n+1)
                            buffer_real = n(1)+i-1
                            buffer_real = buffer_real * (buffer_real + 1)

                            buffer_real = buffer_real * p_nm(n(1)+ii-1, m(1)+i-1) / rho
                            buffer_cmplx = buffer_real * exp_i_m_phi(i) * z_n_cmplx(ii)
                            N_mn(i)%coordinate(ii)%rho = buffer_cmplx

                            buffer_cmplx = z_d_cmplx(ii) / rho * exp_i_m_phi(i)

                            N_mn(i)%coordinate(ii)%theta = tau_nm(n(1)+ii-1, m(1)+i-1) * buffer_cmplx

                            N_mn(i)%coordinate(ii)%phi = cmplx(0, pi_nm(n(1)+ii-1, m(1)+i-1), kind=8) * buffer_cmplx
                        end do
                    end do
            end select

        end subroutine lib_mie_vector_spherical_harmonics_components_real_xu

        ! calculation of the components of the vector spherical harmonic
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
        !   m: integer, dimension(2)
        !       first element: start index
        !       second element: last index
        !       CONDITION: first element .le. second element
        !   n: integer, dimension(2)
        !       first element: start index
        !       second element: last index
        !       CONDITION: first element .le. second element
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
        ! Reference: Absorption and Scattering of Light by Small Particles, eq. 4.17, 4.18, 4.19, 4.20
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

            type(list_spherical_coordinate_cmplx_type), dimension(:), allocatable, intent(inout) :: M_emn
            type(list_spherical_coordinate_cmplx_type), dimension(:), allocatable, intent(inout) :: M_omn
            type(list_spherical_coordinate_cmplx_type), dimension(:), allocatable, intent(inout) :: N_emn
            type(list_spherical_coordinate_cmplx_type), dimension(:), allocatable, intent(inout) :: N_omn

            logical, optional :: not_calc_Memn
            logical, optional :: not_calc_Momn
            logical, optional :: not_calc_Nemn
            logical, optional :: not_calc_Nomn

            ! auxiliary
            integer(kind=VECTOR_SPHERICAL_HARMONICS_COMPONENT_NUMBER_KIND) :: i
            integer(kind=VECTOR_SPHERICAL_HARMONICS_COMPONENT_NUMBER_KIND) :: ii
            integer(kind=VECTOR_SPHERICAL_HARMONICS_COMPONENT_NUMBER_KIND) :: number_of_members_m
            integer(kind=VECTOR_SPHERICAL_HARMONICS_COMPONENT_NUMBER_KIND) :: number_of_members_n

            double precision, dimension(0:n(2), -n(2):n(2)) :: p_nm
            double precision, dimension(0:n(2), -n(2):n(2)) :: p_dnm
            double precision, dimension(0:n(2), -n(2):n(2)) :: p_nm_divided_by_sin_theta

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
                allocate(M_emn(m(2)-m(1)+1))
                do i=1, number_of_members_m
                    allocate (M_emn(i)%coordinate(number_of_members_n))
                end do
            end if

            if (.not. m_not_calc_Momn) then
                allocate(M_omn(m(2)-m(1)+1))
                do i=1, number_of_members_m
                    allocate (M_omn(i)%coordinate(number_of_members_n))
                end do
            end if

            if (.not. m_not_calc_Nemn) then
                allocate(N_emn(m(2)-m(1)+1))
                do i=1, number_of_members_m
                    allocate (N_emn(i)%coordinate(number_of_members_n))
                end do
            end if

            if (.not. m_not_calc_Nemn) then
                allocate(N_omn(m(2)-m(1)+1))
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


            call lib_math_associated_legendre_polynomial_theta(theta, n(2), p_nm_divided_by_sin_theta,  p_dnm)
            p_nm = p_nm_divided_by_sin_theta * sin_theta ! pontential of an error

            ! --- calculations of the components M and N ---
            ! M_emn
            ! eq. (4.17)
            if (.not. m_not_calc_Memn) then
                select case (z_selector)
                    case (1,2)
                        ! z = [j_n, y_n] ==> z: real
                        do i=1, number_of_members_m
                            do ii=1, number_of_members_n
                                buffer_real = -i * sin_m_phi(i) * p_nm_divided_by_sin_theta(n(1)+ii-1, m(1)+i-1) * z_n_real(ii)
                                M_emn(i)%coordinate(ii)%theta = cmplx(buffer_real, 0, kind=8)

                                buffer_real = -cos_m_phi(i) * p_dnm(n(1)+ii-1, m(1)+i-1) * z_n_real(ii)
                                M_emn(i)%coordinate(ii)%phi = cmplx(buffer_real, 0, kind=8)

                                M_emn(i)%coordinate(ii)%rho = cmplx(0,0, kind=8)
                            end do
                        end do
                    case (3,4)
                        ! z = [h^(1)_n, h^(2)_n] ==> z: complex
                        do i=1, number_of_members_m
                            do ii=1, number_of_members_n
                                buffer_real = -i * sin_m_phi(i) * p_nm_divided_by_sin_theta(n(1)+ii-1, m(1)+i-1)
                                buffer_cmplx = buffer_real * z_n_cmplx(ii)
                                M_emn(i)%coordinate(ii)%theta = buffer_cmplx

                                buffer_real = - cos_m_phi(i) * p_dnm(n(1)+ii-1, m(1)+i-1)
                                buffer_cmplx = buffer_real * z_n_cmplx(ii)
                                M_emn(i)%coordinate(ii)%phi = buffer_cmplx

                                M_emn(i)%coordinate(ii)%rho = cmplx(0,0, kind=8)
                            end do
                        end do
                end select
            end if

            ! M_omn
            ! eq. (4.18)
            if (.not. m_not_calc_Momn) then
                select case (z_selector)
                    case (1,2)
                        ! z = [j_n, y_n] ==> z: real
                        do i=1, number_of_members_m
                            do ii=1, number_of_members_n
                                buffer_real = i * cos_m_phi(i) * p_nm_divided_by_sin_theta(n(1)+ii-1, m(1)+i-1) * z_n_real(ii)
                                M_omn(i)%coordinate(ii)%theta = cmplx(buffer_real, 0, kind=8)

                                buffer_real = - sin_m_phi(i) * p_dnm(n(1)+ii-1, m(1)+i-1) * z_n_real(ii)
                                M_omn(i)%coordinate(ii)%phi = cmplx(buffer_real, 0, kind=8)

                                M_omn(i)%coordinate(ii)%rho = cmplx(0,0, kind=8)
                            end do
                        end do
                    case (3,4)
                        ! z = [h^(1)_n, h^(2)_n] ==> z: complex
                        do i=1, number_of_members_m
                            do ii=1, number_of_members_n
                                buffer_real = i * cos_m_phi(i) * p_nm_divided_by_sin_theta(n(1)+ii-1, m(1)+i-1)
                                buffer_cmplx = buffer_real * z_n_cmplx(ii)
                                M_omn(i)%coordinate(ii)%theta = buffer_cmplx

                                buffer_real = - sin_m_phi(i) * p_dnm(n(1)+ii-1, m(1)+i-1)
                                buffer_cmplx = buffer_real * z_n_cmplx(ii)
                                M_omn(i)%coordinate(ii)%phi = buffer_cmplx

                                M_omn(i)%coordinate(ii)%rho = cmplx(0,0, kind=8)
                            end do
                        end do
                end select
            end if

            ! N_emn
            ! eq. (4.19)
            if (.not. m_not_calc_Nemn) then
                select case (z_selector)
                    case (1,2)
                        ! z = [j_n, y_n] ==> z: real
                        do i=1, number_of_members_m
                            do ii=1, number_of_members_n
                                buffer_real = z_divided_by_rho_real(ii) * cos_m_phi(i) * ii*(ii+1) * p_nm(n(1)+ii-1, m(1)+i-1)
                                N_emn(i)%coordinate(ii)%rho = cmplx(buffer_real, 0, kind=8)

                                buffer_real = cos_m_phi(i) * p_dnm(n(1)+ii-1, m(1)+i-1) / rho * r_d_real(ii)
                                N_emn(i)%coordinate(ii)%theta = cmplx(buffer_real, 0, kind=8)

                                buffer_real = - i * sin_m_phi(i) * p_nm_divided_by_sin_theta(n(1)+ii-1, m(1)+i-1) / rho &
                                              * r_d_real(ii)
                                N_emn(i)%coordinate(ii)%phi = cmplx(buffer_real,0, kind=8)
                            end do
                        end do
                    case (3,4)
                        ! z = [h^(1)_n, h^(2)_n] ==> z: complex
                        do i=1, number_of_members_m
                            do ii=1, number_of_members_n
                                buffer_real = cos_m_phi(i) * ii*(ii+1) * p_nm(n(1)+ii-1, m(1)+i-1)
                                buffer_cmplx = z_divided_by_rho_cmplx(ii) * buffer_real
                                N_emn(i)%coordinate(ii)%rho = buffer_cmplx

                                buffer_real = cos_m_phi(i) * p_dnm(n(1)+ii-1, m(1)+i-1) / rho
                                buffer_cmplx = buffer_real * r_d_cmplx(ii)
                                N_emn(i)%coordinate(ii)%theta = buffer_cmplx

                                buffer_real = - i * sin_m_phi(i) * p_nm_divided_by_sin_theta(n(1)+ii-1, m(1)+i-1) / rho
                                buffer_cmplx = buffer_real * r_d_cmplx(ii)
                                N_emn(i)%coordinate(ii)%phi = buffer_cmplx
                            end do
                        end do
                end select
            end if

           ! N_omn
           ! eq. (4.20)
            if (.not. m_not_calc_Nomn) then
                select case (z_selector)
                    case (1,2)
                        ! z = [j_n, y_n] ==> z: real
                        do i=1, number_of_members_m
                            do ii=1, number_of_members_n
                                buffer_real = z_divided_by_rho_real(ii) * sin_m_phi(i) * ii*(ii+1) * p_nm(n(1)+ii-1, m(1)+i-1)
                                N_omn(i)%coordinate(ii)%rho = cmplx(buffer_real, 0, kind=8)

                                buffer_real = sin_m_phi(i) * p_dnm(n(1)+ii-1, m(1)+i-1) / rho * r_d_real(ii)
                                N_omn(i)%coordinate(ii)%theta = cmplx(buffer_real, 0, kind=8)

                                buffer_real = i * cos_m_phi(i) * p_nm_divided_by_sin_theta(n(1)+ii-1, m(1)+i-1) / rho * r_d_real(ii)
                                N_omn(i)%coordinate(ii)%phi = cmplx(buffer_real,0, kind=8)
                            end do
                        end do
                    case (3,4)
                        ! z = [h^(1)_n, h^(2)_n] ==> z: complex
                        do i=1, number_of_members_m
                            do ii=1, number_of_members_n
                                buffer_real = sin_m_phi(i) * ii*(ii+1) * p_nm(n(1)+ii-1, m(1)+i-1)
                                buffer_cmplx = z_divided_by_rho_cmplx(ii) * buffer_real
                                N_omn(i)%coordinate(ii)%rho = buffer_cmplx

                                buffer_real = sin_m_phi(i) * p_dnm(n(1)+ii-1, m(1)+i-1) / rho
                                buffer_cmplx = buffer_real * r_d_cmplx(ii)
                                N_omn(i)%coordinate(ii)%theta = buffer_cmplx

                                buffer_real = i * cos_m_phi(i) * p_nm_divided_by_sin_theta(n(1)+ii-1, m(1)+i-1) / rho
                                buffer_cmplx = buffer_real * r_d_cmplx(ii)
                                N_omn(i)%coordinate(ii)%phi = buffer_cmplx
                            end do
                        end do
                end select
            end if


        end subroutine lib_mie_vector_spherical_harmonics_components_real

        ! calculation of the components of the vector spherical harmonic
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
        subroutine lib_mie_vector_spherical_harmonics_components_cmplx(theta, phi, rho, m, n, z_selector, &
                                                      M_emn, M_omn, N_emn, N_omn, &
                                                      not_calc_Memn, not_calc_Momn, not_calc_Nemn, not_calc_Nomn)
            implicit none
            ! dummy
            double precision, intent(in) :: theta
            double precision, intent(in) :: phi
            complex(kind=8), intent(in) :: rho
            integer(kind=VECTOR_SPHERICAL_HARMONICS_COMPONENT_NUMBER_KIND), dimension(2) :: m
            integer(kind=VECTOR_SPHERICAL_HARMONICS_COMPONENT_NUMBER_KIND), dimension(2) :: n
            integer(kind=1) :: z_selector

            type(list_spherical_coordinate_cmplx_type), dimension(:), allocatable, intent(inout) :: M_emn
            type(list_spherical_coordinate_cmplx_type), dimension(:), allocatable, intent(inout) :: M_omn
            type(list_spherical_coordinate_cmplx_type), dimension(:), allocatable, intent(inout) :: N_emn
            type(list_spherical_coordinate_cmplx_type), dimension(:), allocatable, intent(inout) :: N_omn

            logical, optional :: not_calc_Memn
            logical, optional :: not_calc_Momn
            logical, optional :: not_calc_Nemn
            logical, optional :: not_calc_Nomn

            ! auxiliary
            integer(kind=VECTOR_SPHERICAL_HARMONICS_COMPONENT_NUMBER_KIND) :: i
            integer(kind=VECTOR_SPHERICAL_HARMONICS_COMPONENT_NUMBER_KIND) :: ii
            integer(kind=VECTOR_SPHERICAL_HARMONICS_COMPONENT_NUMBER_KIND) :: number_of_members_m
            integer(kind=VECTOR_SPHERICAL_HARMONICS_COMPONENT_NUMBER_KIND) :: number_of_members_n

            double precision, dimension(0:n(2), -n(2):n(2)) :: p_nm
            double precision, dimension(0:n(2), -n(2):n(2)) :: p_dnm
            double precision, dimension(0:n(2), -n(2):n(2)) :: p_nm_divided_by_sin_theta

            complex(kind=8), dimension(n(2)-n(1)+1) :: z_n_cmplx
            complex(kind=8), dimension(n(2)-n(1)+1) :: z_d_cmplx ! deriviative

            complex(kind=8), dimension(n(2)-n(1)+1) :: z_divided_by_rho_cmplx

            ! Riccati-Bessel
            complex(kind=8), dimension(n(2)-n(1)+1) :: r_cmplx
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
                allocate(M_emn(m(2)-m(1)+1))
                do i=1, number_of_members_m
                    allocate (M_emn(i)%coordinate(number_of_members_n))
                end do
            end if

            if (.not. m_not_calc_Momn) then
                allocate(M_omn(m(2)-m(1)+1))
                do i=1, number_of_members_m
                    allocate (M_omn(i)%coordinate(number_of_members_n))
                end do
            end if

            if (.not. m_not_calc_Nemn) then
                allocate(N_emn(m(2)-m(1)+1))
                do i=1, number_of_members_m
                    allocate (N_emn(i)%coordinate(number_of_members_n))
                end do
            end if

            if (.not. m_not_calc_Nemn) then
                allocate(N_omn(m(2)-m(1)+1))
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
                    r_d_cmplx = lib_math_riccati_s_derivative(rho, n(1), number_of_members_n, r_cmplx)
                    z_n_cmplx = r_cmplx / rho
                case(2)
                    ! spherical Bessel function second kind y_n
                    ! internal: calculation with Riccati-Bessel functions: C_n
                    r_d_cmplx = lib_math_riccati_c_derivative(rho, n(1), number_of_members_n, r_cmplx)
                    z_n_cmplx = r_cmplx / rho
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
                    z_n_cmplx = cmplx(0,0)
                    z_d_cmplx = cmplx(0,0)

                    r_cmplx = cmplx(0,0)
                    r_d_cmplx = cmplx(0,0)
                    print*, "lib_mie_vector_spherical_harmonics_M_emn: ERROR"
                    print*, "  undefined z_selector value: ", z_selector
                    return
            end select

            call lib_math_associated_legendre_polynomial_theta(theta, n(2), p_nm_divided_by_sin_theta,  p_dnm)
            p_nm = p_nm_divided_by_sin_theta * sin_theta

            ! --- calculations of the components M and N ---
            ! M_emn
            ! eq. (4.17)
            if (.not. m_not_calc_Memn) then
                do i=1, number_of_members_m
                    do ii=1, number_of_members_n
                        buffer_real = -i * sin_m_phi(i) * p_nm_divided_by_sin_theta(n(1)+ii-1, m(1)+i-1)
                        buffer_cmplx = buffer_real * z_n_cmplx(ii)
                        M_emn(i)%coordinate(ii)%theta = buffer_cmplx

                        buffer_real = - cos_m_phi(i) * p_dnm(n(1)+ii-1, m(1)+i-1)
                        buffer_cmplx = buffer_real * z_n_cmplx(ii)
                        M_emn(i)%coordinate(ii)%phi = buffer_cmplx

                        M_emn(i)%coordinate(ii)%rho = cmplx(0,0, kind=8)
                    end do
                end do
            end if

            ! M_omn
            ! eq. (4.18)
            if (.not. m_not_calc_Momn) then
                do i=1, number_of_members_m
                    do ii=1, number_of_members_n
                        buffer_real = i * cos_m_phi(i) * p_nm_divided_by_sin_theta(n(1)+ii-1, m(1)+i-1)
                        buffer_cmplx = buffer_real * z_n_cmplx(ii)
                        M_omn(i)%coordinate(ii)%theta = buffer_cmplx

                        buffer_real = - sin_m_phi(i) * p_dnm(n(1)+ii-1, m(1)+i-1)
                        buffer_cmplx = buffer_real * z_n_cmplx(ii)
                        M_omn(i)%coordinate(ii)%phi = buffer_cmplx

                        M_omn(i)%coordinate(ii)%rho = cmplx(0,0, kind=8)
                    end do
                end do
            end if

            ! N_emn
            ! eq. (4.19)
            if (.not. m_not_calc_Nemn) then
                do i=1, number_of_members_m
                    do ii=1, number_of_members_n
                        buffer_real = cos_m_phi(i) * ii*(ii+1) * p_nm(n(1)+ii-1, m(1)+i-1)
                        buffer_cmplx = z_divided_by_rho_cmplx(ii) * buffer_real
                        N_emn(i)%coordinate(ii)%rho = buffer_cmplx

                        buffer_real = cos_m_phi(i) * p_dnm(n(1)+ii-1, m(1)+i-1)
                        buffer_cmplx = buffer_real * r_d_cmplx(ii) / rho
                        N_emn(i)%coordinate(ii)%theta = buffer_cmplx

                        buffer_real = - i * sin_m_phi(i) * p_nm_divided_by_sin_theta(n(1)+ii-1, m(1)+i-1)
                        buffer_cmplx = buffer_real * r_d_cmplx(ii) / rho
                        N_emn(i)%coordinate(ii)%phi = buffer_cmplx
                    end do
                end do
            end if

           ! N_omn
           ! eq. (4.20)
            if (.not. m_not_calc_Nomn) then
                do i=1, number_of_members_m
                    do ii=1, number_of_members_n
                        buffer_real = sin_m_phi(i) * ii*(ii+1) * p_nm(n(1)+ii-1, m(1)+i-1)
                        buffer_cmplx = z_divided_by_rho_cmplx(ii) * buffer_real
                        N_omn(i)%coordinate(ii)%rho = buffer_cmplx

                        buffer_real = sin_m_phi(i) * p_dnm(n(1)+ii-1, m(1)+i-1)
                        buffer_cmplx = buffer_real * r_d_cmplx(ii) / rho
                        N_omn(i)%coordinate(ii)%theta = buffer_cmplx

                        buffer_real = i * cos_m_phi(i) * p_nm_divided_by_sin_theta(n(1)+ii-1, m(1)+i-1)
                        buffer_cmplx = buffer_real * r_d_cmplx(ii) / rho
                        N_omn(i)%coordinate(ii)%phi = buffer_cmplx
                    end do
                end do
            end if


        end subroutine lib_mie_vector_spherical_harmonics_components_cmplx


        function lib_mie_vector_spherical_harmonics_test_functions() result (rv)
            implicit none
            ! dummy
            integer :: rv

            ! auxiliaray
            double precision, parameter :: ground_truth_e = 10.0_8**(-6.0_8)

            rv = 0

            if (.not. test_lib_mie_vector_spherical_harmonics_components_real_xu()) then
                rv = rv + 1
            end if

            print *, "----lib_mie_vector_spherical_harmonics_test_functions----"
            if (rv == 0) then
                print *, "lib_mie_vector_spherical_harmonics_test_functions tests: OK"
            else
                print *, rv,"lib_mie_vector_spherical_harmonics_test_functions test(s) FAILED"
            end if
            print *, "---------------------------------------------------------"

            contains

            function test_lib_mie_vector_spherical_harmonics_components_real() result (rv)
                implicit none
                ! dummy
                logical :: rv

                ! auxiliary
                double precision :: theta
                double precision :: phi
                double precision :: rho
                integer(kind=VECTOR_SPHERICAL_HARMONICS_COMPONENT_NUMBER_KIND), dimension(2) :: m
                integer(kind=VECTOR_SPHERICAL_HARMONICS_COMPONENT_NUMBER_KIND), dimension(2) :: n
                integer(kind=1) :: z_selector

                type(list_spherical_coordinate_cmplx_type), dimension(:), allocatable :: M_emn
                type(list_spherical_coordinate_cmplx_type), dimension(:), allocatable :: M_omn
                type(list_spherical_coordinate_cmplx_type), dimension(:), allocatable :: N_emn
                type(list_spherical_coordinate_cmplx_type), dimension(:), allocatable :: N_omn

!                logical :: not_calc_Memn
!                logical :: not_calc_Momn
!                logical :: not_calc_Nemn
!                logical :: not_calc_Nomn

                theta = 0
                phi = 0
                rho = 100

                z_selector = 3

                m = (/ 1, 1 /)
                n = (/ 1, 3 /)

                call lib_mie_vector_spherical_harmonics_components_real(theta, phi, rho, m, n, z_selector, &
                                                                        M_emn, M_omn, N_emn, N_omn)


            end function

            function test_lib_mie_vector_spherical_harmonics_components_real_xu() result (rv)
                implicit none
                ! dummy
                logical :: rv

                ! auxiliary
                double precision :: theta
                double precision :: phi
                double precision :: rho
                integer(kind=VECTOR_SPHERICAL_HARMONICS_COMPONENT_NUMBER_KIND), dimension(2) :: m
                integer(kind=VECTOR_SPHERICAL_HARMONICS_COMPONENT_NUMBER_KIND), dimension(2) :: n
                integer(kind=1) :: z_selector

                type(list_spherical_coordinate_cmplx_type), dimension(:), allocatable :: M_mn
                type(list_spherical_coordinate_cmplx_type), dimension(:), allocatable :: N_mn

!                logical :: not_calc_Memn
!                logical :: not_calc_Momn
!                logical :: not_calc_Nemn
!                logical :: not_calc_Nomn

                theta = 0
                phi = 0
                rho = 100

                z_selector = 3

                m = (/ 1, 1 /)
                n = (/ 1, 3 /)

                call lib_mie_vector_spherical_harmonics_components_real_xu(theta, phi, rho, m, n, z_selector, &
                                                                           M_mn, N_mn)


            end function

        end function lib_mie_vector_spherical_harmonics_test_functions

end module lib_mie_vector_spherical_harmonics
