#define _DEBUG_


module lib_mie_vector_spherical_harmonics
    use libmath
    implicit none

    private

    ! --- public ---
    public :: lib_mie_vector_spherical_harmonics_components

    interface lib_mie_vector_spherical_harmonics_components
        module procedure lib_mie_vector_spherical_harmonics_components_real_xu
        module procedure lib_mie_vector_spherical_harmonics_components_cmplx_xu
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
        !   r: double precision
        !       distance [m]
        !   k: double precision
        !       wavenumber [1/m]
        !   n_range: integer, dimension(2)
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
        !   M_nm: type(list_spherical_coordinate_cmplx_type)
        !       M component of the vector spherical harmonic
        !   N_nm: type(list_spherical_coordinate_cmplx_type)
        !       N component of the vector spherical harmonic
        !
        ! LaTeX: $$ \mathbf{M}_{m n}^{(J)}=\left[\mathbf{i}_{\theta} i \pi_{m n}(\cos \theta)-\mathbf{i}_{\phi} \tau_{m n}(\cos \theta)\right] z_{n}^{(J)} )(k r) \exp (i m \phi) $$
        !        $$ \mathbf{N}_{m n}^{(J)}=\mathbf{i}_r n\left(n+1 ) P_{n}^{m} (\cos \theta\right) \frac{z_{n}^{(J)}(k r)}{k r} \exp(i m \phi)
        !           + \left[\mathbf{i}_{\theta} \tau_{m n}(\cos \theta)+\mathbf{i}_{\phi} i \pi_{m n}(\cos \theta)\right]
        !           \times \frac{1}{k r} \frac{d}{d r}\left[r z_{n}^{(J)}(k r)\right] \exp (i m \phi) $$
        !
        ! Reference: Electromagnetic scattering by an aggregate of spheres, Yu-lin Xu, eq. 2
        subroutine lib_mie_vector_spherical_harmonics_components_real_xu(theta, phi, r, k, n_range, z_selector, &
                                                      M_nm, N_nm)
            implicit none
            
            ! dummy
            double precision, intent(in) :: theta
            double precision, intent(in) :: phi
            double precision, intent(in) :: r
            double precision, intent(in) :: k
            integer(kind=4), dimension(2) :: n_range
            integer(kind=1) :: z_selector

            type(list_spherical_coordinate_cmplx_type), dimension(:), allocatable, intent(inout) :: M_nm
            type(list_spherical_coordinate_cmplx_type), dimension(:), allocatable, intent(inout) :: N_nm

            ! auxiliary
            integer :: n
            integer :: m
            integer :: i
            integer :: ii

            double precision :: rho

            integer(kind=4) :: number_of_members_n

            type(list_list_real) :: pi_nm
            type(list_list_real) :: tau_nm

            type(list_list_real) :: p_nm
            type(list_list_real) :: p_d_nm

            double precision, dimension(n_range(2)-n_range(1)+1) :: buffer_p_n
            double precision, dimension(n_range(2)-n_range(1)+1) :: buffer_p_d_n

            double precision, dimension(2, n_range(2)-n_range(1)+1) :: buffer_p_n_m_neg
            double precision, dimension(2, n_range(2)-n_range(1)+1) :: buffer_p_d_n_m_neg

            double precision, dimension(n_range(2)-n_range(1)+1) :: z_n_real
            complex(kind=8), dimension(n_range(2)-n_range(1)+1) :: z_n_cmplx
            double precision, dimension(n_range(2)-n_range(1)+1) :: z_d_real ! deriviative
            complex(kind=8), dimension(n_range(2)-n_range(1)+1) :: z_d_cmplx ! deriviative

            double precision, dimension(n_range(2)-n_range(1)+1) :: z_divided_by_rho_real
            complex(kind=8), dimension(n_range(2)-n_range(1)+1) :: z_divided_by_rho_cmplx

            ! Riccati-Bessel
            double precision, dimension(n_range(2)-n_range(1)+1) :: r_real
            complex(kind=8), dimension(n_range(2)-n_range(1)+1) :: r_cmplx
            double precision, dimension(n_range(2)-n_range(1)+1) :: r_d_real ! deriviative
            complex(kind=8), dimension(n_range(2)-n_range(1)+1) :: r_d_cmplx ! deriviative

            double precision :: buffer_real
            complex(kind=8) :: buffer_cmplx

            complex(kind=8), dimension(-n_range(2):n_range(2)) :: exp_i_m_phi
            double precision :: cos_theta

            number_of_members_n = n_range(2) - n_range(1) + 1

            ! --- init ---
            call init_list(p_nm, n_range(1), number_of_members_n)
            call init_list(p_d_nm, n_range(1), number_of_members_n)

            allocate( M_nm(n_range(1):n_range(2)) )
            do i=n_range(1), n_range(2)
                allocate (M_nm(i)%coordinate(-i:i))
            end do

            allocate( N_nm(n_range(1):n_range(2)) )
            do i=n_range(1), n_range(2)
                allocate (N_nm(i)%coordinate(-i:i))
            end do


            ! --- pre-calculation ---
            rho = k * r

            cos_theta = cos(theta)

            do i=-n_range(2), n_range(2)
                exp_i_m_phi(i)= cmplx(cos(i * phi), sin(i * phi), kind=8)
            end do

            select case (z_selector)
                case(1)
                    ! spherical Bessel function first kind j_n
                    ! internal: calculation with Riccati-Bessel functions: S_n
                    r_d_real = lib_math_riccati_s_derivative(rho, n_range(1), number_of_members_n, r_real)
                    z_n_real = r_real / rho
                case(2)
                    ! spherical Bessel function second kind y_n
                    ! internal: calculation with Riccati-Bessel functions: C_n
                    r_d_real = lib_math_riccati_c_derivative(rho, n_range(1), number_of_members_n, r_real)
                    z_n_real = r_real / rho
                case(3)
                    ! spherical Hankel function first kind   h^(1)_n
                    ! internal: calculation with Riccati-Bessel functions: Xi_n
                    r_d_cmplx = lib_math_riccati_xi_derivative(rho, n_range(1), number_of_members_n, r_cmplx)
                    z_n_cmplx = r_cmplx / rho
                case(4)
                    ! spherical Hankel function second kind   h^(2)_n
                    ! internal: calculation with Riccati-Bessel functions: Zeta_n
                    r_d_cmplx = lib_math_riccati_zeta_derivative(rho, n_range(1), number_of_members_n, r_cmplx)
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
                    print*, "lib_mie_vector_spherical_harmonics_components_real_xu: ERROR"
                    print*, "  undefined z_selector value: ", z_selector
                    return
            end select


            call lib_math_associated_legendre_polynomial_theta(theta, n_range(2), pi_nm, tau_nm, .false.)


            ! set p_nm(m .eq. 0)
            call lib_math_associated_legendre_polynomial(cos_theta, 0, n_range(1), number_of_members_n , &
                                                         buffer_p_n, buffer_p_d_n, .false.)
            do n=n_range(1), n_range(2)
                i = n - n_range(1) + 1
                p_nm%item(n)%item(0) = buffer_p_n(i)
            end do

            ! set p_nm(m .ne. 0)
            !       ->n
            !          +
            !         ++
            !    ^   +++  <-  m
            !   m|  ++++
            !        +++  <- -m
            !         ++
            !          +
            do m=1, n_range(2)
                call lib_math_associated_legendre_polynomial_with_negative_m(cos_theta, m, m, n_range(2) - m + 1, &
                                                                         buffer_p_n_m_neg, buffer_p_d_n_m_neg, .false.)
                do n=m, n_range(2)
                    i = n - m + 1
                    p_nm%item(n)%item( m) = buffer_p_n_m_neg(2, i)
                    p_nm%item(n)%item(-m) = buffer_p_n_m_neg(1, i)
                end do
            end do

            ! --- calculations of the components M and N ---
            ! M_mn
            ! first line eq. (2)
            select case (z_selector)
                case (1,2)
                    ! z = [j_n, y_n] ==> z: real
                    do n=n_range(1), n_range(2)
                        i = n - n_range(1) + 1
                        do m=-n, n
                            buffer_real = pi_nm%item(n)%item(m) * z_n_real(i)
                            M_nm(n)%coordinate(m)%theta = cmplx(0, buffer_real, kind=8) * exp_i_m_phi(m)

                            buffer_real = -tau_nm%item(n)%item(m) * z_n_real(i)
                            buffer_cmplx = cmplx(buffer_real, 0, kind=8) * exp_i_m_phi(m)
                            M_nm(n)%coordinate(m)%phi = buffer_cmplx

                            M_nm(n)%coordinate(m)%rho = cmplx(0,0, kind=8)
#ifdef _DEBUG_
                    if (isnan(real(M_nm(n)%coordinate(m)%rho))   .or. isnan(aimag(M_nm(n)%coordinate(m)%rho)) .or. &
                        isnan(real(M_nm(n)%coordinate(m)%phi))   .or. isnan(aimag(M_nm(n)%coordinate(m)%phi)) .or. &
                        isnan(real(M_nm(n)%coordinate(m)%theta)) .or. isnan(aimag(M_nm(n)%coordinate(m)%theta)) ) then
                        print *, "lib_mie_vector_spherical_harmonics_components_real_xu: ERROR"
                        print *, "  M_nm(n)%coordinate(m) is NaN"
                        print *, "  n = ", n
                        print *, "  m = ", m
                        print *, "  z_selector: ", z_selector
                    end if
#endif
                        end do
                    end do
                case (3,4)
                    ! z = [h^(1)_n, h^(2)_n] ==> z: complex
                    do n=n_range(1), n_range(2)
                        i = n - n_range(1) + 1
                        do m=-n, n
                            buffer_real = pi_nm%item(n)%item(m)
                            M_nm(n)%coordinate(m)%theta = cmplx(0, buffer_real, kind=8) * z_n_cmplx(i) * exp_i_m_phi(m)

                            buffer_real = -tau_nm%item(n)%item(m)
                            buffer_cmplx = cmplx(buffer_real, 0, kind=8) * z_n_cmplx(i) * exp_i_m_phi(m)
                            M_nm(n)%coordinate(m)%phi = buffer_cmplx

                            M_nm(n)%coordinate(m)%rho = cmplx(0,0, kind=8)
#ifdef _DEBUG_
                    if (isnan(real(M_nm(n)%coordinate(m)%rho))   .or. isnan(aimag(M_nm(n)%coordinate(m)%rho)) .or. &
                        isnan(real(M_nm(n)%coordinate(m)%phi))   .or. isnan(aimag(M_nm(n)%coordinate(m)%phi)) .or. &
                        isnan(real(M_nm(n)%coordinate(m)%theta)) .or. isnan(aimag(M_nm(n)%coordinate(m)%theta)) ) then
                        print *, "lib_mie_vector_spherical_harmonics_components_real_xu: ERROR"
                        print *, "  M_nm(n)%coordinate(m) is NaN"
                        print *, "  n = ", n
                        print *, "  m = ", m
                        print *, "  z_selector: ", z_selector
                    end if
#endif
                        end do
                    end do
            end select

            ! N_mn
            ! second line eq. (2)
            select case (z_selector)
                case (1,2)
                    ! z = [j_n, y_n] ==> z: real
                    do n=n_range(1), n_range(2)
                        i = n - n_range(1) + 1
                        do m=-n, n
                            buffer_real = real(n*(n + 1), kind=8) * p_nm%item(n)%item(m) * z_n_real(i) / rho
                            buffer_cmplx = cmplx(buffer_real, 0, kind=8) * exp_i_m_phi(m)
                            N_nm(n)%coordinate(m)%rho = buffer_cmplx

                            buffer_cmplx = cmplx(r_d_real(i) / (k * rho), 0, kind=8) * exp_i_m_phi(m)

                            N_nm(n)%coordinate(m)%theta = tau_nm%item(n)%item(m) * buffer_cmplx

                            N_nm(n)%coordinate(m)%phi = cmplx(0, pi_nm%item(n)%item(m), kind=8) * buffer_cmplx
#ifdef _DEBUG_
                    if (isnan(real(N_nm(n)%coordinate(m)%rho))   .or. isnan(aimag(N_nm(n)%coordinate(m)%rho)) .or. &
                        isnan(real(N_nm(n)%coordinate(m)%phi))   .or. isnan(aimag(N_nm(n)%coordinate(m)%phi)) .or. &
                        isnan(real(N_nm(n)%coordinate(m)%theta)) .or. isnan(aimag(N_nm(n)%coordinate(m)%theta)) ) then
                        print *, "lib_mie_vector_spherical_harmonics_components_real_xu: ERROR"
                        print *, "  N_nm(n)%coordinate(m) is NaN"
                        print *, "  n = ", n
                        print *, "  m = ", m
                        print *, "  z_selector: ", z_selector
                    end if
#endif
                        end do
                    end do
                case (3,4)
                    ! z = [h^(1)_n, h^(2)_n] ==> z: complex
                    do n=n_range(1), n_range(2)
                        i = n - n_range(1) + 1
                        do m=-n, n
                            buffer_real = real(n*(n + 1), kind=8) * p_nm%item(n)%item(m) / rho
                            buffer_cmplx = cmplx(buffer_real, 0, kind=8) * exp_i_m_phi(m) * z_n_cmplx(i)
                            N_nm(n)%coordinate(m)%rho = buffer_cmplx

                            buffer_cmplx = r_d_cmplx(i) * exp_i_m_phi(m) / (k * rho)

                            N_nm(n)%coordinate(m)%theta = tau_nm%item(n)%item(m) * buffer_cmplx

                            N_nm(n)%coordinate(m)%phi = cmplx(0, pi_nm%item(n)%item(m), kind=8) * buffer_cmplx
#ifdef _DEBUG_
                    if (isnan(real(N_nm(n)%coordinate(m)%rho))   .or. isnan(aimag(N_nm(n)%coordinate(m)%rho)) .or. &
                        isnan(real(N_nm(n)%coordinate(m)%phi))   .or. isnan(aimag(N_nm(n)%coordinate(m)%phi)) .or. &
                        isnan(real(N_nm(n)%coordinate(m)%theta)) .or. isnan(aimag(N_nm(n)%coordinate(m)%theta)) ) then
                        print *, "lib_mie_vector_spherical_harmonics_components_real_xu: ERROR"
                        print *, "  N_nm(n)%coordinate(m) is NaN"
                        print *, "  n = ", n
                        print *, "  m = ", m
                        print *, "  z_selector: ", z_selector
                    end if
#endif
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
        !   k: complex
        !       wavenumber [1/m]
        !   r: double precision
        !       distance [m]
        !   n_range: integer, dimension(2)
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
        subroutine lib_mie_vector_spherical_harmonics_components_cmplx_xu(theta, phi, k, r, n_range, z_selector, &
                                                      M_nm, N_nm)
            implicit none

            ! dummy
            double precision, intent(in) :: theta
            double precision, intent(in) :: phi
            complex(kind=8), intent(in) :: k
            double precision, intent(in) :: r
            integer(kind=4), dimension(2) :: n_range
            integer(kind=1) :: z_selector

            type(list_spherical_coordinate_cmplx_type), dimension(:), allocatable, intent(inout) :: M_nm
            type(list_spherical_coordinate_cmplx_type), dimension(:), allocatable, intent(inout) :: N_nm

            ! auxiliary
            integer :: n
            integer :: m
            integer :: i

            complex(kind=8) :: rho

            integer(kind=4) :: number_of_members_n

            type(list_list_real) :: pi_nm
            type(list_list_real) :: tau_nm

            type(list_list_real) :: p_nm
            type(list_list_real) :: p_d_nm

            double precision, dimension(n_range(2)-n_range(1)+1) :: buffer_p_n
            double precision, dimension(n_range(2)-n_range(1)+1) :: buffer_p_d_n

            double precision, dimension(2, n_range(2)-n_range(1)+1) :: buffer_p_n_m_neg
            double precision, dimension(2, n_range(2)-n_range(1)+1) :: buffer_p_d_n_m_neg

            complex(kind=8), dimension(n_range(2)-n_range(1)+1) :: z_n_cmplx
            complex(kind=8), dimension(n_range(2)-n_range(1)+1) :: z_d_cmplx ! deriviative

            ! Riccati-Bessel
            complex(kind=8), dimension(n_range(2)-n_range(1)+1) :: r_cmplx
            complex(kind=8), dimension(n_range(2)-n_range(1)+1) :: r_d_cmplx ! deriviative

            double precision :: buffer_real
            complex(kind=8) :: buffer_cmplx

            complex(kind=8), dimension(-n_range(2):n_range(2)) :: exp_i_m_phi
            double precision :: cos_theta

            number_of_members_n = n_range(2) - n_range(1) + 1

            ! --- init ---
            call init_list(p_nm, n_range(1), number_of_members_n)
            call init_list(p_d_nm, n_range(1), number_of_members_n)

            allocate( M_nm(n_range(1):n_range(2)) )
            do i=n_range(1), n_range(2)
                allocate (M_nm(i)%coordinate(-i:i))
            end do

            allocate( N_nm(n_range(1):n_range(2)) )
            do i=n_range(1), n_range(2)
                allocate (N_nm(i)%coordinate(-i:i))
            end do


            ! --- pre-calculation ---
            rho = k * r

            cos_theta = cos(theta)

            do i=-n_range(2), n_range(2)
                exp_i_m_phi(i)= cmplx(cos(i * phi), sin(i * phi), kind=8)
            end do

            select case (z_selector)
                case(1)
                    ! spherical Bessel function first kind j_n
                    ! internal: calculation with Riccati-Bessel functions: S_n
                    r_d_cmplx = lib_math_riccati_s_derivative(rho, n_range(1), number_of_members_n, r_cmplx)
                    z_n_cmplx = r_cmplx / rho
                case(2)
                    ! spherical Bessel function second kind y_n
                    ! internal: calculation with Riccati-Bessel functions: C_n
                    r_d_cmplx = lib_math_riccati_c_derivative(rho, n_range(1), number_of_members_n, r_cmplx)
                    z_n_cmplx = r_cmplx / rho
                case(3)
                    ! spherical Hankel function first kind   h^(1)_n
                    ! internal: calculation with Riccati-Bessel functions: Xi_n
                    r_d_cmplx = lib_math_riccati_xi_derivative(rho, n_range(1), number_of_members_n, r_cmplx)
                    z_n_cmplx = r_cmplx / rho
                case(4)
                    ! spherical Hankel function first kind   h^(2)_n
                    ! internal: calculation with Riccati-Bessel functions: Zeta_n
                    r_d_cmplx = lib_math_riccati_zeta_derivative(rho, n_range(1), number_of_members_n, r_cmplx)
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


            call lib_math_associated_legendre_polynomial_theta(theta, n_range(2), pi_nm, tau_nm, .false.)


            ! set p_nm(m .eq. 0)
            call lib_math_associated_legendre_polynomial(cos_theta, 0, n_range(1), number_of_members_n , &
                                                         buffer_p_n, buffer_p_d_n, .false.)
            do n=n_range(1), n_range(2)
                i = n - n_range(1) + 1
                p_nm%item(n)%item(0) = buffer_p_n(i)
            end do

            ! set p_nm(m .ne. 0)
            !       ->n
            !          +
            !         ++
            !    ^   +++  <-  m
            !   m|  ++++
            !        +++  <- -m
            !         ++
            !          +
            do m=1, n_range(2)
                call lib_math_associated_legendre_polynomial_with_negative_m(cos_theta, m, m, n_range(2) - m + 1, &
                                                                         buffer_p_n_m_neg, buffer_p_d_n_m_neg, .false.)
                do n=m, n_range(2)
                    i = n - m + 1
                    p_nm%item(n)%item( m) = buffer_p_n_m_neg(2, i)
                    p_nm%item(n)%item(-m) = buffer_p_n_m_neg(1, i)
                end do
            end do

            ! --- calculations of the components M and N ---
            ! M_mn
            ! first line eq. (2)
            do n=n_range(1), n_range(2)
                i = n - n_range(1) + 1
                do m=-n, n
                    buffer_real = pi_nm%item(n)%item(m)
                    M_nm(n)%coordinate(m)%theta = cmplx(0, buffer_real, kind=8) * z_n_cmplx(i) * exp_i_m_phi(m)

                    buffer_real = -tau_nm%item(n)%item(m)
                    buffer_cmplx = cmplx(buffer_real, 0, kind=8) * z_n_cmplx(i) * exp_i_m_phi(m)
                    M_nm(n)%coordinate(m)%phi = buffer_cmplx

                    M_nm(n)%coordinate(m)%rho = cmplx(0,0, kind=8)
#ifdef _DEBUG_
            if (isnan(real(M_nm(n)%coordinate(m)%rho))   .or. isnan(aimag(M_nm(n)%coordinate(m)%rho)) .or. &
                isnan(real(M_nm(n)%coordinate(m)%phi))   .or. isnan(aimag(M_nm(n)%coordinate(m)%phi)) .or. &
                isnan(real(M_nm(n)%coordinate(m)%theta)) .or. isnan(aimag(M_nm(n)%coordinate(m)%theta)) ) then
                print *, "lib_mie_vector_spherical_harmonics_components_real_xu: ERROR"
                print *, "  M_nm(n)%coordinate(m) is NaN"
                print *, "  n = ", n
                print *, "  m = ", m
                print *, "  z_selector: ", z_selector
            end if
#endif
                end do
            end do

            ! N_mn
            ! second line eq. (2)
            do n=n_range(1), n_range(2)
                i = n - n_range(1) + 1
                do m=-n, n
                    buffer_real = real(n*(n + 1), kind=8) * p_nm%item(n)%item(m)
                    buffer_cmplx = cmplx(buffer_real, 0, kind=8) / rho * exp_i_m_phi(m) * z_n_cmplx(i)
                    N_nm(n)%coordinate(m)%rho = buffer_cmplx

                    buffer_cmplx = r_d_cmplx(i) * exp_i_m_phi(m) / (k * rho)

                    N_nm(n)%coordinate(m)%theta = tau_nm%item(n)%item(m) * buffer_cmplx

                    N_nm(n)%coordinate(m)%phi = cmplx(0, pi_nm%item(n)%item(m), kind=8) * buffer_cmplx
#ifdef _DEBUG_
            if (isnan(real(N_nm(n)%coordinate(m)%rho))   .or. isnan(aimag(N_nm(n)%coordinate(m)%rho)) .or. &
                isnan(real(N_nm(n)%coordinate(m)%phi))   .or. isnan(aimag(N_nm(n)%coordinate(m)%phi)) .or. &
                isnan(real(N_nm(n)%coordinate(m)%theta)) .or. isnan(aimag(N_nm(n)%coordinate(m)%theta)) ) then
                print *, "lib_mie_vector_spherical_harmonics_components_real_xu: ERROR"
                print *, "  N_nm(n)%coordinate(m) is NaN"
                print *, "  n = ", n
                print *, "  m = ", m
                print *, "  z_selector: ", z_selector
            end if
#endif
                end do
            end do

        end subroutine lib_mie_vector_spherical_harmonics_components_cmplx_xu

!        ! calculation of the components of the vector spherical harmonic
!        !
!        ! Argument
!        ! ----
!        !   theta: double precision
!        !       polar angle
!        !   phi: double precision
!        !       azimuthal angle
!        !   rho: double precision
!        !       dimensionless varibale rho = k*r
!        !       k: wavenumber
!        !       r: distance
!        !   m: integer, dimension(2)
!        !       first element: start index
!        !       second element: last index
!        !       CONDITION: first element .le. second element
!        !   n: integer, dimension(2)
!        !       first element: start index
!        !       second element: last index
!        !       CONDITION: first element .le. second element
!        !   z_selector: integer
!        !       1: spherical Bessel function first kind   j_n
!        !       2: spherical Bessel function second kind  y_n
!        !       3: spherical Hankel function first kind   h^(1)_n
!        !       4: spherical Hankel function second kind  h^(2)_n
!        !
!        ! Results
!        ! ----
!        !   rv: complex, dimension(3)
!        !       values of the spherical coordinates (rho, theta, phi)
!        !
!        ! LaTeX: $$ \begin{aligned} \mathbf{M}_{e m n}=& \frac{-m}{\sin \theta} \sin m \phi P_{n}^{m}(\cos \theta) z_{n}(\rho) \hat{\mathbf{e}}_{\theta} \\ &-\cos m \phi \frac{d P_{n}^{m}(\cos \theta)}{d \theta} z_{n}(\rho) \hat{\mathbf{e}}_{\phi} \end{aligned} $$
!        !
!        ! Reference: Absorption and Scattering of Light by Small Particles, eq. 4.17, 4.18, 4.19, 4.20
!        subroutine lib_mie_vector_spherical_harmonics_components_real(theta, phi, rho, m, n, z_selector, &
!                                                      M_emn, M_omn, N_emn, N_omn, &
!                                                      not_calc_Memn, not_calc_Momn, not_calc_Nemn, not_calc_Nomn)
!            implicit none
!            ! dummy
!            double precision, intent(in) :: theta
!            double precision, intent(in) :: phi
!            double precision, intent(in) :: rho
!            integer(kind=VECTOR_SPHERICAL_HARMONICS_COMPONENT_NUMBER_KIND), dimension(2) :: m
!            integer(kind=VECTOR_SPHERICAL_HARMONICS_COMPONENT_NUMBER_KIND), dimension(2) :: n
!            integer(kind=1) :: z_selector
!
!            type(list_spherical_coordinate_cmplx_type), dimension(:), allocatable, intent(inout) :: M_emn
!            type(list_spherical_coordinate_cmplx_type), dimension(:), allocatable, intent(inout) :: M_omn
!            type(list_spherical_coordinate_cmplx_type), dimension(:), allocatable, intent(inout) :: N_emn
!            type(list_spherical_coordinate_cmplx_type), dimension(:), allocatable, intent(inout) :: N_omn
!
!            logical, optional :: not_calc_Memn
!            logical, optional :: not_calc_Momn
!            logical, optional :: not_calc_Nemn
!            logical, optional :: not_calc_Nomn
!
!            ! auxiliary
!            integer(kind=VECTOR_SPHERICAL_HARMONICS_COMPONENT_NUMBER_KIND) :: i
!            integer(kind=VECTOR_SPHERICAL_HARMONICS_COMPONENT_NUMBER_KIND) :: ii
!            integer(kind=VECTOR_SPHERICAL_HARMONICS_COMPONENT_NUMBER_KIND) :: number_of_members_m
!            integer(kind=VECTOR_SPHERICAL_HARMONICS_COMPONENT_NUMBER_KIND) :: number_of_members_n
!
!            double precision, dimension(0:n(2), -n(2):n(2)) :: p_nm
!            double precision, dimension(0:n(2), -n(2):n(2)) :: p_dnm
!            double precision, dimension(0:n(2), -n(2):n(2)) :: p_nm_divided_by_sin_theta
!
!            double precision, dimension(n(2)-n(1)+1) :: z_n_real
!            complex(kind=8), dimension(n(2)-n(1)+1) :: z_n_cmplx
!            double precision, dimension(n(2)-n(1)+1) :: z_d_real ! deriviative
!            complex(kind=8), dimension(n(2)-n(1)+1) :: z_d_cmplx ! deriviative
!
!            double precision, dimension(n(2)-n(1)+1) :: z_divided_by_rho_real
!            complex(kind=8), dimension(n(2)-n(1)+1) :: z_divided_by_rho_cmplx
!
!            ! Riccati-Bessel
!            double precision, dimension(n(2)-n(1)+1) :: r_real
!            complex(kind=8), dimension(n(2)-n(1)+1) :: r_cmplx
!            double precision, dimension(n(2)-n(1)+1) :: r_d_real ! deriviative
!            complex(kind=8), dimension(n(2)-n(1)+1) :: r_d_cmplx ! deriviative
!
!            double precision :: cos_theta
!            double precision :: sin_theta
!            double precision :: minus_sin_theta
!            double precision, dimension(m(2)-m(1)+1) :: m_divided_by_sin_theta
!
!            double precision, dimension(n(2)-n(1)+1) :: cos_m_phi
!            double precision, dimension(n(2)-n(1)+1) :: sin_m_phi
!
!
!            logical :: m_not_calc_Memn
!            logical :: m_not_calc_Momn
!            logical :: m_not_calc_Nemn
!            logical :: m_not_calc_Nomn
!
!            double precision :: buffer_real
!            complex(kind=8) :: buffer_cmplx
!
!
!            number_of_members_m = m(2) - m(1) + 1
!            number_of_members_n = n(2) - n(1) + 1
!
!            ! --- standard value ---
!            m_not_calc_Memn = .false.
!            m_not_calc_Momn = .false.
!            m_not_calc_Nemn = .false.
!            m_not_calc_Nomn = .false.
!
!            if (present(not_calc_Memn)) then
!                m_not_calc_Memn = not_calc_Memn
!            end if
!
!            if (present(not_calc_Momn)) then
!                m_not_calc_Momn = not_calc_Momn
!            end if
!
!            if (present(not_calc_Nemn)) then
!                m_not_calc_Nemn = not_calc_Nemn
!            end if
!            if (present(not_calc_Nomn)) then
!                m_not_calc_Nomn = not_calc_Nomn
!            end if
!
!            ! --- init ---
!            if (.not. m_not_calc_Memn) then
!                allocate(M_emn(m(2)-m(1)+1))
!                do i=1, number_of_members_m
!                    allocate (M_emn(i)%coordinate(number_of_members_n))
!                end do
!            end if
!
!            if (.not. m_not_calc_Momn) then
!                allocate(M_omn(m(2)-m(1)+1))
!                do i=1, number_of_members_m
!                    allocate (M_omn(i)%coordinate(number_of_members_n))
!                end do
!            end if
!
!            if (.not. m_not_calc_Nemn) then
!                allocate(N_emn(m(2)-m(1)+1))
!                do i=1, number_of_members_m
!                    allocate (N_emn(i)%coordinate(number_of_members_n))
!                end do
!            end if
!
!            if (.not. m_not_calc_Nemn) then
!                allocate(N_omn(m(2)-m(1)+1))
!                do i=1, number_of_members_m
!                    allocate (N_omn(i)%coordinate(number_of_members_n))
!                end do
!            end if
!
!
!            ! --- pre-calculation ---
!            cos_theta = cos(theta)
!            sin_theta = sin(theta)
!            minus_sin_theta = -sin(theta)
!
!
!            do i=1, number_of_members_m
!                cos_m_phi(i) = cos(i*phi)
!                sin_m_phi(i) = sin(i*phi)
!                m_divided_by_sin_theta(i) = i / sin_theta
!            end do
!
!            select case (z_selector)
!                case(1)
!                    ! spherical Bessel function first kind j_n
!                    ! internal: calculation with Riccati-Bessel functions: S_n
!                    r_d_real = lib_math_riccati_s_derivative(rho, n(1), number_of_members_n, r_real)
!                    z_n_real = r_real / rho
!                case(2)
!                    ! spherical Bessel function second kind y_n
!                    ! internal: calculation with Riccati-Bessel functions: C_n
!                    r_d_real = lib_math_riccati_c_derivative(rho, n(1), number_of_members_n, r_real)
!                    z_n_real = r_real / rho
!                case(3)
!                    ! spherical Hankel function first kind   h^(1)_n
!                    ! internal: calculation with Riccati-Bessel functions: Xi_n
!                    r_d_cmplx = lib_math_riccati_xi_derivative(rho, n(1), number_of_members_n, r_cmplx)
!                    z_n_cmplx = r_cmplx / rho
!                case(4)
!                    ! spherical Hankel function first kind   h^(2)_n
!                    ! internal: calculation with Riccati-Bessel functions: Zeta_n
!                    r_d_cmplx = lib_math_riccati_zeta_derivative(rho, n(1), number_of_members_n, r_cmplx)
!                    z_n_cmplx = r_cmplx / rho
!                case default
!                    z_n_real = 0
!                    z_n_cmplx = cmplx(0,0)
!                    z_d_real = 0
!                    z_d_cmplx = cmplx(0,0)
!
!                    r_real = 0
!                    r_cmplx = cmplx(0,0)
!                    r_d_real = 0
!                    r_d_cmplx = cmplx(0,0)
!                    print*, "lib_mie_vector_spherical_harmonics_M_emn: ERROR"
!                    print*, "  undefined z_selector value: ", z_selector
!                    return
!            end select
!
!
!            call lib_math_associated_legendre_polynomial_theta(theta, n(2), p_nm_divided_by_sin_theta,  p_dnm)
!            p_nm = p_nm_divided_by_sin_theta * sin_theta ! pontential of an error
!
!            ! --- calculations of the components M and N ---
!            ! M_emn
!            ! eq. (4.17)
!            if (.not. m_not_calc_Memn) then
!                select case (z_selector)
!                    case (1,2)
!                        ! z = [j_n, y_n] ==> z: real
!                        do i=1, number_of_members_m
!                            do ii=1, number_of_members_n
!                                buffer_real = -i * sin_m_phi(i) * p_nm_divided_by_sin_theta(n(1)+ii-1, m(1)+i-1) * z_n_real(ii)
!                                M_emn(i)%coordinate(ii)%theta = cmplx(buffer_real, 0, kind=8)
!
!                                buffer_real = -cos_m_phi(i) * p_dnm(n(1)+ii-1, m(1)+i-1) * z_n_real(ii)
!                                M_emn(i)%coordinate(ii)%phi = cmplx(buffer_real, 0, kind=8)
!
!                                M_emn(i)%coordinate(ii)%rho = cmplx(0,0, kind=8)
!                            end do
!                        end do
!                    case (3,4)
!                        ! z = [h^(1)_n, h^(2)_n] ==> z: complex
!                        do i=1, number_of_members_m
!                            do ii=1, number_of_members_n
!                                buffer_real = -i * sin_m_phi(i) * p_nm_divided_by_sin_theta(n(1)+ii-1, m(1)+i-1)
!                                buffer_cmplx = buffer_real * z_n_cmplx(ii)
!                                M_emn(i)%coordinate(ii)%theta = buffer_cmplx
!
!                                buffer_real = - cos_m_phi(i) * p_dnm(n(1)+ii-1, m(1)+i-1)
!                                buffer_cmplx = buffer_real * z_n_cmplx(ii)
!                                M_emn(i)%coordinate(ii)%phi = buffer_cmplx
!
!                                M_emn(i)%coordinate(ii)%rho = cmplx(0,0, kind=8)
!                            end do
!                        end do
!                end select
!            end if
!
!            ! M_omn
!            ! eq. (4.18)
!            if (.not. m_not_calc_Momn) then
!                select case (z_selector)
!                    case (1,2)
!                        ! z = [j_n, y_n] ==> z: real
!                        do i=1, number_of_members_m
!                            do ii=1, number_of_members_n
!                                buffer_real = i * cos_m_phi(i) * p_nm_divided_by_sin_theta(n(1)+ii-1, m(1)+i-1) * z_n_real(ii)
!                                M_omn(i)%coordinate(ii)%theta = cmplx(buffer_real, 0, kind=8)
!
!                                buffer_real = - sin_m_phi(i) * p_dnm(n(1)+ii-1, m(1)+i-1) * z_n_real(ii)
!                                M_omn(i)%coordinate(ii)%phi = cmplx(buffer_real, 0, kind=8)
!
!                                M_omn(i)%coordinate(ii)%rho = cmplx(0,0, kind=8)
!                            end do
!                        end do
!                    case (3,4)
!                        ! z = [h^(1)_n, h^(2)_n] ==> z: complex
!                        do i=1, number_of_members_m
!                            do ii=1, number_of_members_n
!                                buffer_real = i * cos_m_phi(i) * p_nm_divided_by_sin_theta(n(1)+ii-1, m(1)+i-1)
!                                buffer_cmplx = buffer_real * z_n_cmplx(ii)
!                                M_omn(i)%coordinate(ii)%theta = buffer_cmplx
!
!                                buffer_real = - sin_m_phi(i) * p_dnm(n(1)+ii-1, m(1)+i-1)
!                                buffer_cmplx = buffer_real * z_n_cmplx(ii)
!                                M_omn(i)%coordinate(ii)%phi = buffer_cmplx
!
!                                M_omn(i)%coordinate(ii)%rho = cmplx(0,0, kind=8)
!                            end do
!                        end do
!                end select
!            end if
!
!            ! N_emn
!            ! eq. (4.19)
!            if (.not. m_not_calc_Nemn) then
!                select case (z_selector)
!                    case (1,2)
!                        ! z = [j_n, y_n] ==> z: real
!                        do i=1, number_of_members_m
!                            do ii=1, number_of_members_n
!                                buffer_real = z_divided_by_rho_real(ii) * cos_m_phi(i) * ii*(ii+1) * p_nm(n(1)+ii-1, m(1)+i-1)
!                                N_emn(i)%coordinate(ii)%rho = cmplx(buffer_real, 0, kind=8)
!
!                                buffer_real = cos_m_phi(i) * p_dnm(n(1)+ii-1, m(1)+i-1) / rho * r_d_real(ii)
!                                N_emn(i)%coordinate(ii)%theta = cmplx(buffer_real, 0, kind=8)
!
!                                buffer_real = - i * sin_m_phi(i) * p_nm_divided_by_sin_theta(n(1)+ii-1, m(1)+i-1) / rho &
!                                              * r_d_real(ii)
!                                N_emn(i)%coordinate(ii)%phi = cmplx(buffer_real,0, kind=8)
!                            end do
!                        end do
!                    case (3,4)
!                        ! z = [h^(1)_n, h^(2)_n] ==> z: complex
!                        do i=1, number_of_members_m
!                            do ii=1, number_of_members_n
!                                buffer_real = cos_m_phi(i) * ii*(ii+1) * p_nm(n(1)+ii-1, m(1)+i-1)
!                                buffer_cmplx = z_divided_by_rho_cmplx(ii) * buffer_real
!                                N_emn(i)%coordinate(ii)%rho = buffer_cmplx
!
!                                buffer_real = cos_m_phi(i) * p_dnm(n(1)+ii-1, m(1)+i-1) / rho
!                                buffer_cmplx = buffer_real * r_d_cmplx(ii)
!                                N_emn(i)%coordinate(ii)%theta = buffer_cmplx
!
!                                buffer_real = - i * sin_m_phi(i) * p_nm_divided_by_sin_theta(n(1)+ii-1, m(1)+i-1) / rho
!                                buffer_cmplx = buffer_real * r_d_cmplx(ii)
!                                N_emn(i)%coordinate(ii)%phi = buffer_cmplx
!                            end do
!                        end do
!                end select
!            end if
!
!           ! N_omn
!           ! eq. (4.20)
!            if (.not. m_not_calc_Nomn) then
!                select case (z_selector)
!                    case (1,2)
!                        ! z = [j_n, y_n] ==> z: real
!                        do i=1, number_of_members_m
!                            do ii=1, number_of_members_n
!                                buffer_real = z_divided_by_rho_real(ii) * sin_m_phi(i) * ii*(ii+1) * p_nm(n(1)+ii-1, m(1)+i-1)
!                                N_omn(i)%coordinate(ii)%rho = cmplx(buffer_real, 0, kind=8)
!
!                                buffer_real = sin_m_phi(i) * p_dnm(n(1)+ii-1, m(1)+i-1) / rho * r_d_real(ii)
!                                N_omn(i)%coordinate(ii)%theta = cmplx(buffer_real, 0, kind=8)
!
!                                buffer_real = i * cos_m_phi(i) * p_nm_divided_by_sin_theta(n(1)+ii-1, m(1)+i-1) / rho * r_d_real(ii)
!                                N_omn(i)%coordinate(ii)%phi = cmplx(buffer_real,0, kind=8)
!                            end do
!                        end do
!                    case (3,4)
!                        ! z = [h^(1)_n, h^(2)_n] ==> z: complex
!                        do i=1, number_of_members_m
!                            do ii=1, number_of_members_n
!                                buffer_real = sin_m_phi(i) * ii*(ii+1) * p_nm(n(1)+ii-1, m(1)+i-1)
!                                buffer_cmplx = z_divided_by_rho_cmplx(ii) * buffer_real
!                                N_omn(i)%coordinate(ii)%rho = buffer_cmplx
!
!                                buffer_real = sin_m_phi(i) * p_dnm(n(1)+ii-1, m(1)+i-1) / rho
!                                buffer_cmplx = buffer_real * r_d_cmplx(ii)
!                                N_omn(i)%coordinate(ii)%theta = buffer_cmplx
!
!                                buffer_real = i * cos_m_phi(i) * p_nm_divided_by_sin_theta(n(1)+ii-1, m(1)+i-1) / rho
!                                buffer_cmplx = buffer_real * r_d_cmplx(ii)
!                                N_omn(i)%coordinate(ii)%phi = buffer_cmplx
!                            end do
!                        end do
!                end select
!            end if
!
!
!        end subroutine lib_mie_vector_spherical_harmonics_components_real

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

!            call lib_math_associated_legendre_polynomial_theta(theta, n(2), p_nm_divided_by_sin_theta,  p_dnm)
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

        ! Calculation of the tranlation transformation coefficients
        ! from the l-th coordinate system to the j-th coordinate system
        !
        ! Argument
        ! ----
        !   x: double precision
        !       normalized distance: x = k * d_lj
        !       k: wave number
        !       d_lj: distance from origin l to origin j
        !   theta: double precision
        !       polar coordinate [rad]
        !   phi: double precision
        !       azimuthal coordinate [rad]
        !   n_range: integer, dimension(2)
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
        ! Returns
        ! ----
        !   A:
        !
        ! Reference: Experimental and theoretical results of light scattering by aggregates of spheres, Yu-lin Xu and Bo Å. S. Gustafson
        subroutine lib_mie_vector_spherical_harmonics_tranlation_coefficient_real(x, theta, phi, n_range, z_selector,&
                                                                                  A_mnkl, B_mnkl)
            implicit none
            ! dummy
            double precision, intent(in) :: x
            double precision, intent(in) :: theta
            double precision, intent(in) :: phi
            integer(kind=4), dimension(2) :: n_range
            integer(kind=1) :: z_selector

            type(list_spherical_coordinate_cmplx_type), dimension(:), allocatable, intent(inout) :: A_mnkl
            type(list_spherical_coordinate_cmplx_type), dimension(:), allocatable, intent(inout) :: B_mnkl

            ! auxiliary

            print*, "test"

        end subroutine lib_mie_vector_spherical_harmonics_tranlation_coefficient_real


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
            if (.not. test_lib_mie_vector_spherical_harmonics_components_cmplx_xu()) then
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

!            function test_lib_mie_vector_spherical_harmonics_components_real() result (rv)
!                implicit none
!                ! dummy
!                logical :: rv
!
!                ! auxiliary
!                double precision :: theta
!                double precision :: phi
!                double precision :: rho
!                integer(kind=VECTOR_SPHERICAL_HARMONICS_COMPONENT_NUMBER_KIND), dimension(2) :: m
!                integer(kind=VECTOR_SPHERICAL_HARMONICS_COMPONENT_NUMBER_KIND), dimension(2) :: n
!                integer(kind=1) :: z_selector
!
!                type(list_spherical_coordinate_cmplx_type), dimension(:), allocatable :: M_emn
!                type(list_spherical_coordinate_cmplx_type), dimension(:), allocatable :: M_omn
!                type(list_spherical_coordinate_cmplx_type), dimension(:), allocatable :: N_emn
!                type(list_spherical_coordinate_cmplx_type), dimension(:), allocatable :: N_omn
!
!!                logical :: not_calc_Memn
!!                logical :: not_calc_Momn
!!                logical :: not_calc_Nemn
!!                logical :: not_calc_Nomn
!
!                theta = 0
!                phi = 0
!                rho = 100
!
!                z_selector = 3
!
!                m = (/ 1, 1 /)
!                n = (/ 1, 3 /)
!
!                call lib_mie_vector_spherical_harmonics_components_real(theta, phi, rho, m, n, z_selector, &
!                                                                        M_emn, M_omn, N_emn, N_omn)
!
!
!            end function

            function test_lib_mie_vector_spherical_harmonics_components_real_xu() result (rv)
                use file_io
                implicit none
                ! dummy
                logical :: rv

                ! parameter
                character(len=*), parameter :: file_name_M_mn = &
                     "src/lib/mie_theory/lib_mie_vector_spherical_harmonics/ground_truth_M_mn.csv"
                 character(len=*), parameter :: file_name_N_mn = &
                     "src/lib/mie_theory/lib_mie_vector_spherical_harmonics/ground_truth_N_mn.csv"

                ! auxiliary
                integer :: i
                integer :: ii
                integer :: n_value
                integer :: m_value

                double precision :: buffer
                complex(kind=8) :: buffer_cmplx

                double precision :: theta
                double precision :: phi
                double precision :: k
                double precision :: r
                integer(kind=VECTOR_SPHERICAL_HARMONICS_COMPONENT_NUMBER_KIND), dimension(2) :: n
                integer(kind=1) :: z_selector

                type(list_spherical_coordinate_cmplx_type), dimension(:), allocatable :: M_mn
                type(list_spherical_coordinate_cmplx_type), dimension(:), allocatable :: N_mn

                type(list_spherical_coordinate_cmplx_type), dimension(:), allocatable :: ground_truth_M_mn
                type(list_spherical_coordinate_cmplx_type), dimension(:), allocatable :: ground_truth_N_mn

                double precision, dimension(:,:), allocatable :: csv_data
                integer :: csv_columns

                theta = 0.2
                phi = 0
                ! n = 1
                ! lam = 10**-6 m
                ! k = n * 2 Pi / lam
                k = 2.0_8 * PI * 10.0_8**(6.0_8)
                ! r = 50 * 10**-6 m
                r = 50.0_8 * 10.0_8**(-6.0_8)

                z_selector = 3

                n = (/ 1, 4 /)

                rv = .false.

                call init_list(ground_truth_M_mn, n(1), n(2)-n(1)+1)
                call init_list(ground_truth_N_mn, n(1), n(2)-n(1)+1)

                ! load ground truth M_mn
                if (file_exists(file_name_M_mn)) then
                    csv_columns = 8
                    call read_csv(file_name_M_mn, csv_columns, csv_data)

                    do i=lbound(csv_data, 1), ubound(csv_data, 1)
                        n_value = int(csv_data(i, 1))
                        m_value = int(csv_data(i, 2))

                        buffer_cmplx = cmplx(csv_data(i, 3), csv_data(i, 4), kind=8)
                        ground_truth_M_mn(n_value)%coordinate(m_value)%rho = buffer_cmplx

                        buffer_cmplx = cmplx(csv_data(i, 5), csv_data(i, 6), kind=8)
                        ground_truth_M_mn(n_value)%coordinate(m_value)%theta = buffer_cmplx

                        buffer_cmplx = cmplx(csv_data(i, 7), csv_data(i, 8), kind=8)
                        ground_truth_M_mn(n_value)%coordinate(m_value)%phi = buffer_cmplx
                    end do
                else
                    print *, "test_lib_mie_vector_spherical_harmonics_components_real_xu: ERROR"
                    print *, "  file does not exist"
                    print *, "  file_name: ", file_name_M_mn

                    return
                end if

                ! load ground truth M_mn
                if (file_exists(file_name_N_mn)) then
                    csv_columns = 8
                    call read_csv(file_name_N_mn, csv_columns, csv_data)

                    do i=lbound(csv_data, 1), ubound(csv_data, 1)
                        n_value = int(csv_data(i, 1))
                        m_value = int(csv_data(i, 2))

                        buffer_cmplx = cmplx(csv_data(i, 3), csv_data(i, 4), kind=8)
                        ground_truth_N_mn(n_value)%coordinate(m_value)%rho = buffer_cmplx

                        buffer_cmplx = cmplx(csv_data(i, 5), csv_data(i, 6), kind=8)
                        ground_truth_N_mn(n_value)%coordinate(m_value)%theta = buffer_cmplx

                        buffer_cmplx = cmplx(csv_data(i, 7), csv_data(i, 8), kind=8)
                        ground_truth_N_mn(n_value)%coordinate(m_value)%phi = buffer_cmplx
                    end do
                else
                    print *, "test_lib_mie_vector_spherical_harmonics_components_real_xu: ERROR"
                    print *, "  file does not exist"
                    print *, "  file_name: ", file_name_N_mn

                    return
                end if

                ! calculate M_mn, N_mn
                call lib_mie_vector_spherical_harmonics_components_real_xu(theta, phi, r, k, n, z_selector, &
                                                                           M_mn, N_mn)

                ! evaluate
                rv = .true.
                print *, "test_lib_mie_vector_spherical_harmonics_components_real_xu:"
                print *, "  M_mn:"
                do i=n(1), n(2)
                    do ii=-i, i
                        buffer = abs(spherical_abs(M_mn(i)%coordinate(ii) - ground_truth_M_mn(i)%coordinate(ii)))
                        if (buffer .gt. ground_truth_e) then
                            print *, "  n: ", i ," m: ", ii, "difference: ", buffer, " : FAILED"
                            rv = .false.
                        else
                            print *, "  n: ", i ," m: ", ii, ": OK"
                        end if
                    end do
                end do

                print *, "  N_mn:"
                do i=n(1), n(2)
                    do ii=-i, i
                        buffer = abs(spherical_abs(N_mn(i)%coordinate(ii) - ground_truth_N_mn(i)%coordinate(ii)))
                        if (buffer .gt. ground_truth_e) then
                            print *, "  n: ", i ," m: ", ii, "difference: ", buffer, " : FAILED"
                            rv = .false.
                        else
                            print *, "  n: ", i ," m: ", ii, ": OK"
                        end if
                    end do
                end do

            end function test_lib_mie_vector_spherical_harmonics_components_real_xu

            function test_lib_mie_vector_spherical_harmonics_components_cmplx_xu() result (rv)
                use file_io
                implicit none
                ! dummy
                logical :: rv

                ! parameter
                character(len=*), parameter :: file_name_M_mn = &
                     "src/lib/mie_theory/lib_mie_vector_spherical_harmonics/ground_truth_M_mn_cmplx.csv"
                 character(len=*), parameter :: file_name_N_mn = &
                     "src/lib/mie_theory/lib_mie_vector_spherical_harmonics/ground_truth_N_mn_cmplx.csv"

                ! auxiliary
                integer :: i
                integer :: ii
                integer :: n_value
                integer :: m_value

                double precision :: buffer
                complex(kind=8) :: buffer_cmplx

                double precision :: theta
                double precision :: phi
                complex(kind=8) :: k
                double precision :: r
                integer(kind=VECTOR_SPHERICAL_HARMONICS_COMPONENT_NUMBER_KIND), dimension(2) :: n
                integer(kind=1) :: z_selector

                type(list_spherical_coordinate_cmplx_type), dimension(:), allocatable :: M_mn
                type(list_spherical_coordinate_cmplx_type), dimension(:), allocatable :: N_mn

                type(list_spherical_coordinate_cmplx_type), dimension(:), allocatable :: ground_truth_M_mn
                type(list_spherical_coordinate_cmplx_type), dimension(:), allocatable :: ground_truth_N_mn

                double precision, dimension(:,:), allocatable :: csv_data
                integer :: csv_columns

                theta = 0.2
                phi = 0
                k = cmplx(1.0, 0.4, kind=8) * 2.0 * PI / 10.0_8**(-6.0_8)
                r = 50.0_8 * 10.0_8**(-6.0_8)

                z_selector = 3

                n = (/ 1, 4 /)

                rv = .false.

                call init_list(ground_truth_M_mn, n(1), n(2)-n(1)+1)
                call init_list(ground_truth_N_mn, n(1), n(2)-n(1)+1)

                ! load ground truth M_mn
                if (file_exists(file_name_M_mn)) then
                    csv_columns = 8
                    call read_csv(file_name_M_mn, csv_columns, csv_data)

                    do i=lbound(csv_data, 1), ubound(csv_data, 1)
                        n_value = int(csv_data(i, 1))
                        m_value = int(csv_data(i, 2))

                        buffer_cmplx = cmplx(csv_data(i, 3), csv_data(i, 4), kind=8)
                        ground_truth_M_mn(n_value)%coordinate(m_value)%rho = buffer_cmplx

                        buffer_cmplx = cmplx(csv_data(i, 5), csv_data(i, 6), kind=8)
                        ground_truth_M_mn(n_value)%coordinate(m_value)%theta = buffer_cmplx

                        buffer_cmplx = cmplx(csv_data(i, 7), csv_data(i, 8), kind=8)
                        ground_truth_M_mn(n_value)%coordinate(m_value)%phi = buffer_cmplx
                    end do
                else
                    print *, "test_lib_mie_vector_spherical_harmonics_components_cmplx_xu: ERROR"
                    print *, "  file does not exist"
                    print *, "  file_name: ", file_name_M_mn

                    return
                end if

                ! load ground truth M_mn
                if (file_exists(file_name_N_mn)) then
                    csv_columns = 8
                    call read_csv(file_name_N_mn, csv_columns, csv_data)

                    do i=lbound(csv_data, 1), ubound(csv_data, 1)
                        n_value = int(csv_data(i, 1))
                        m_value = int(csv_data(i, 2))

                        buffer_cmplx = cmplx(csv_data(i, 3), csv_data(i, 4), kind=8)
                        ground_truth_N_mn(n_value)%coordinate(m_value)%rho = buffer_cmplx

                        buffer_cmplx = cmplx(csv_data(i, 5), csv_data(i, 6), kind=8)
                        ground_truth_N_mn(n_value)%coordinate(m_value)%theta = buffer_cmplx

                        buffer_cmplx = cmplx(csv_data(i, 7), csv_data(i, 8), kind=8)
                        ground_truth_N_mn(n_value)%coordinate(m_value)%phi = buffer_cmplx
                    end do
                else
                    print *, "test_lib_mie_vector_spherical_harmonics_components_cmplx_xu: ERROR"
                    print *, "  file does not exist"
                    print *, "  file_name: ", file_name_N_mn

                    return
                end if

                ! calculate M_mn, N_mn
                call lib_mie_vector_spherical_harmonics_components_cmplx_xu(theta, phi, k, r, n, z_selector, &
                                                                           M_mn, N_mn)

                ! evaluate
                rv = .true.
                print *, "test_lib_mie_vector_spherical_harmonics_components_cmplx_xu:"
                print *, "  M_mn:"
                do i=n(1), n(2)
                    do ii=-i, i
                        buffer = abs(spherical_abs(M_mn(i)%coordinate(ii) - ground_truth_M_mn(i)%coordinate(ii)))
                        if (buffer .gt. ground_truth_e) then
                            print *, "  n: ", i ," m: ", ii, "difference: ", buffer, " : FAILED"
                            rv = .false.
                        else
                            print *, "  n: ", i ," m: ", ii, ": OK"
                        end if
                    end do
                end do

                print *, "  N_mn:"
                do i=n(1), n(2)
                    do ii=-i, i
                        buffer = abs(spherical_abs(N_mn(i)%coordinate(ii) - ground_truth_N_mn(i)%coordinate(ii)))
                        if (buffer .gt. ground_truth_e) then
                            print *, "  n: ", i ," m: ", ii, "difference: ", buffer, " : FAILED"
                            rv = .false.
                        else
                            print *, "  n: ", i ," m: ", ii, ": OK"
                        end if
                    end do
                end do

            end function test_lib_mie_vector_spherical_harmonics_components_cmplx_xu

        end function lib_mie_vector_spherical_harmonics_test_functions

end module lib_mie_vector_spherical_harmonics
