! Scattering by a sphere
! ----
!
! - incident plane wave along the z-axis
!
#define _DEBUG_

module lib_mie_scattering_by_a_sphere
    !$  use omp_lib
    use libmath
    use lib_constants
    use lib_mie_vector_spherical_harmonics
    use lib_mie_type
    implicit none

    interface lib_mie_scattering_by_a_sphere_init
        module procedure lib_mie_scattering_by_a_sphere_init_coeff_ab_real
        module procedure lib_mie_scattering_by_a_sphere_init_coeff_ab_cmplx
    end interface

    ! caching
    type cache_coefficients_a_b_real_barberh_type
        double precision :: x
        double precision :: m
        integer :: n_max
        complex(kind=8), dimension(:), allocatable :: a_n
        complex(kind=8), dimension(:), allocatable :: b_n
    end type

    type cache_coefficients_a_b_cmplx_barberh_type
        double precision :: x
        complex(kind=8) :: m
        integer :: n_max
        complex(kind=8), dimension(:), allocatable :: a_n
        complex(kind=8), dimension(:), allocatable :: b_n
    end type

    type cache_coefficients_p_0_q_0_type
        double precision :: alpha
        double precision :: beta
        integer :: n_max
        type(list_list_cmplx), allocatable :: p_0
        type(list_list_cmplx), allocatable :: q_0
    end type

    type(cache_coefficients_a_b_real_barberh_type), dimension(:), allocatable :: cache_coefficients_a_b_real_barberh
    logical :: cache_coefficients_a_b_real_barberh_enabled = .false.

    type(cache_coefficients_a_b_cmplx_barberh_type), dimension(:), allocatable :: cache_coefficients_a_b_cmplx_barberh
    logical :: cache_coefficients_a_b_cmplx_barberh_enabled = .false.

    type(cache_coefficients_p_0_q_0_type), dimension(:), allocatable :: cache_coefficients_p_0_q_0
    logical :: cache_coefficients_p_0_q_0_enabled = .false.

    private

    ! --- public ---
    public:: lib_mie_scattering_by_a_sphere_test_functions

    contains

        ! Argument
        ! ----
        !   x: double precision
        !       size parameter: x = k*r = 2 * PI * N * r / lambda
        !       k: wavenumber
        !       r: distance
        !       N: refractive index of the medium
        !       lambda: wave length
        !   m: real
        !       relative refractive index: m = N_1 / N
        !       N_1: refractive index of the particle
        !       N: refractive index of the medium
        subroutine lib_mie_scattering_by_a_sphere_init_coeff_ab_real(x, m, n_max)
            implicit none
            ! dummy
            double precision, dimension(:), intent(in) :: x
            double precision, dimension(size(x)), intent(in) :: m
            integer, dimension(size(x)), intent(in) :: n_max

            ! auxiliary
            integer :: i

            complex(kind=8), dimension(:), allocatable :: a_n
            complex(kind=8), dimension(:), allocatable :: b_n

            ! --- init: cache_coefficients_a_b_real_barberh_x ---
            if (allocated(cache_coefficients_a_b_real_barberh)) then
                deallocate(cache_coefficients_a_b_real_barberh)
            end if

            allocate(cache_coefficients_a_b_real_barberh(size(x)))
            cache_coefficients_a_b_real_barberh_enabled = .false.

            do i=1, size(x)
                cache_coefficients_a_b_real_barberh%x = x(i)
                cache_coefficients_a_b_real_barberh%m = m(i)
                cache_coefficients_a_b_real_barberh%n_max = n_max(i)

                allocate(a_n(n_max(i)))
                allocate(b_n(n_max(i)))

                call get_coefficients_a_b_real_barberh(x(i), m(i), (/ 1, n_max(i) /), a_n, b_n)

                cache_coefficients_a_b_real_barberh(i)%a_n = a_n
                cache_coefficients_a_b_real_barberh(i)%b_n = b_n

                deallocate(a_n)
                deallocate(b_n)
            end do

            cache_coefficients_a_b_real_barberh_enabled = .true.

        end subroutine lib_mie_scattering_by_a_sphere_init_coeff_ab_real

        ! Argument
        ! ----
        !   x: double precision
        !       size parameter: x = k*r = 2 * PI * N * r / lambda
        !       k: wavenumber
        !       r: distance
        !       N: refractive index of the medium
        !       lambda: wave length
        !   m: complex
        !       relative refractive index: m = N_1 / N
        !       N_1: refractive index of the particle
        !       N: refractive index of the medium
        subroutine lib_mie_scattering_by_a_sphere_init_coeff_ab_cmplx(x, m, n_max)
            implicit none
            ! dummy
            double precision, dimension(:), intent(in) :: x
            complex(kind=8), dimension(size(x)), intent(in) :: m
            integer, dimension(size(x)), intent(in) :: n_max

            ! auxiliary
            integer :: i

            complex(kind=8), dimension(:), allocatable :: a_n
            complex(kind=8), dimension(:), allocatable :: b_n

            ! --- init: cache_coefficients_a_b_cmplx_barberh_x ---
            if (allocated(cache_coefficients_a_b_cmplx_barberh)) then
                deallocate(cache_coefficients_a_b_cmplx_barberh)
            end if

            allocate(cache_coefficients_a_b_cmplx_barberh(size(x)))
            cache_coefficients_a_b_cmplx_barberh_enabled = .false.

            do i=1, size(x)
                cache_coefficients_a_b_cmplx_barberh%x = x(i)
                cache_coefficients_a_b_cmplx_barberh%m = m(i)
                cache_coefficients_a_b_cmplx_barberh%n_max = n_max(i)

                allocate(a_n(n_max(i)))
                allocate(b_n(n_max(i)))

                call get_coefficients_a_b_cmplx_barberh(x(i), m(i), (/ 1, n_max(i) /), a_n, b_n)

                cache_coefficients_a_b_cmplx_barberh(i)%a_n = a_n
                cache_coefficients_a_b_cmplx_barberh(i)%b_n = b_n

                deallocate(a_n)
                deallocate(b_n)
            end do

            cache_coefficients_a_b_cmplx_barberh_enabled = .true.

        end subroutine lib_mie_scattering_by_a_sphere_init_coeff_ab_cmplx

        subroutine lib_mie_scattering_by_a_sphere_init_coeff_p0_q0(alpha, beta, n_max)
            implicit none
            ! dummy
            double precision, dimension(:), intent(in) :: alpha
            double precision, dimension(lbound(alpha, 1):ubound(alpha, 1)), intent(in) :: beta
            integer(kind=4), dimension(lbound(alpha, 1):ubound(alpha, 1)) :: n_max

            ! auxiliary
            integer :: i

            ! --- init: cache_coefficients_a_b_cmplx_barberh_x ---
            if (allocated(cache_coefficients_p_0_q_0)) then
                deallocate(cache_coefficients_p_0_q_0)
            end if

            allocate(cache_coefficients_p_0_q_0(size(alpha)))
            cache_coefficients_p_0_q_0_enabled = .false.

            do i=1, size(alpha)
                cache_coefficients_p_0_q_0%alpha = alpha(i)
                cache_coefficients_p_0_q_0%beta = beta(i)
                cache_coefficients_p_0_q_0%n_max = n_max(i)

                call lib_mie_get_p_q_j_j_core(alpha(i), beta(i), (/1, n_max(i)/), &
                                              cache_coefficients_p_0_q_0(i)%p_0, &
                                              cache_coefficients_p_0_q_0(i)%q_0)
            end do

            cache_coefficients_p_0_q_0_enabled = .true.

        end subroutine

        subroutine destructor
            implicit none

            ! --- deallocate caching ---
            if (allocated(cache_coefficients_a_b_real_barberh)) then
                deallocate(cache_coefficients_a_b_real_barberh)
            end if
            cache_coefficients_a_b_real_barberh_enabled = .false.

            if (allocated(cache_coefficients_a_b_cmplx_barberh)) then
                deallocate(cache_coefficients_a_b_cmplx_barberh)
            end if
            cache_coefficients_a_b_cmplx_barberh_enabled = .false.

            if (allocated(cache_coefficients_p_0_q_0)) then
                deallocate(cache_coefficients_p_0_q_0)
            end if
            cache_coefficients_p_0_q_0_enabled = .false.

        end subroutine destructor

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
        !       e.g. n = (/ 1, 45 /) <-- size parameter
        !   alpha: double precision, optional(std: 0)
        !       incident angle with respect to the z-axis
        !       codomain: 0..Pi
        !   beta: double precision, optional(std: 0)
        !       angle between the x axis and the projection of the wave vector on the x-y plane
        !       codomain: 0..2Pi
        !
        ! Results
        ! ----
        !   field_s: type(spherical_coordinate_cmplx_type), dimension(2)
        !       1: scattered e field
        !       2: scattered h field
        !
        ! Reference: Electromagnetic scattering by an aggregate of spheres, Yu-lin Xu
        !
        function get_field_scattered_xu_real(theta, phi, r, &
                                             e_field_0, lambda, &
                                             n_medium, r_particle, n_particle, &
                                             n_range, &
                                             alpha, beta) &
                                           result (field_s)
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
            integer(kind=4), dimension(2),intent(in) :: n_range
            double precision, intent(in), optional :: alpha
            double precision, intent(in), optional :: beta

            type(spherical_coordinate_cmplx_type), dimension(2) :: field_s

            ! auxiliary
            type(spherical_coordinate_cmplx_type) :: e_field_s
            type(spherical_coordinate_cmplx_type) :: h_field_s

            double precision :: k0 ! wave number = 2 pi / lambda
            double precision :: size_parameter ! r_particle * k0

            double precision :: buffer_real
            complex(kind=8) :: buffer_cmplx
            integer(kind=4) :: i
!            integer(kind=4) :: ii
            integer(kind=4) :: m
            integer(kind=4) :: n
!            double precision :: mu
!            double precision :: mu1

            complex(kind=8), dimension(n_range(2)-n_range(1)+1) :: a_n
            complex(kind=8), dimension(n_range(2)-n_range(1)+1) :: b_n

            integer(kind=1) :: z_selector

            type(list_spherical_coordinate_cmplx_type), dimension(:), allocatable :: M_nm
            type(list_spherical_coordinate_cmplx_type), dimension(:), allocatable :: N_nm

            type(list_list_cmplx) :: e_field_nm
            type(spherical_coordinate_cmplx_type), dimension(n_range(2)-n_range(1)+1) :: e_field_n_s
            type(spherical_coordinate_cmplx_type) :: buffer_e_field_n_s
            type(spherical_coordinate_cmplx_type), dimension(n_range(2)-n_range(1)+1) :: h_field_n_s
            type(spherical_coordinate_cmplx_type) :: buffer_h_field_n_s

            type(list_list_cmplx) :: p0_nm
            type(list_list_cmplx) :: q0_nm
            type(list_list_real) :: pi_nm
            type(list_list_real) :: tau_nm

            double precision :: cos_beta
            double precision :: sin_beta

            double precision :: m_alpha
            double precision :: m_beta

            type(list_list_logical) :: calc_order_m

            ! --- init ---
            ! standard values
            if (present(alpha)) then
                m_alpha = alpha
            else
                m_alpha = 0
            end if
            if (present(beta)) then
                m_beta = beta
            else
                m_beta = 0
            end if

            k0 = 2.0_8 * PI / lambda
            size_parameter = r_particle * k0 * n_medium ! eq. 11

            cos_beta = cos(m_beta)
            sin_beta = sin(m_beta)

!            mu = 1
!            mu1 = 1

            call init_list(e_field_nm, n_range(1), n_range(2)-n_range(1)+1)
            call init_list(p0_nm, n_range(1), n_range(2)-n_range(1)+1)
            call init_list(q0_nm, n_range(1), n_range(2)-n_range(1)+1)
            call init_list(calc_order_m, n_range(1), n_range(2)-n_range(1)+1)

            ! --- pre-calc ---
            call get_coefficients_a_b_real_barberh(size_parameter, n_particle/n_medium, n_range, a_n, b_n)
!            call get_coefficients_a_b_real(size_parameter, n_particle/n_medium, mu, mu1, n_range, a_n, b_n)

            ! errata eq. (1) Equations (21) on p. 4577
            call lib_math_associated_legendre_polynomial_theta(m_alpha, n_range(2), pi_nm, tau_nm)
            !$OMP PARALLEL DO PRIVATE(n, m, buffer_real, buffer_cmplx)
            do n=n_range(1), n_range(2)
                do m=-n, n
                    buffer_real = -m * m_beta
                    buffer_cmplx = cmplx(cos(buffer_real), sin(buffer_real), kind=8)

                    buffer_real = 1.0_8 / real(n * (n + 1), kind=8)
                    buffer_cmplx = buffer_real * buffer_cmplx

                    p0_nm%item(n)%item(m) = tau_nm%item(n)%item(m) * buffer_cmplx
                    q0_nm%item(n)%item(m) = pi_nm%item(n)%item(m) * buffer_cmplx

                    if(p0_nm%item(n)%item(m) .eq. cmplx(0,0) &
                       .and. q0_nm%item(n)%item(m) .eq. cmplx(0,0)) then
                        calc_order_m%item(n)%item(m) = .false.
                    end if
                end do
            end do
            !$OMP END PARALLEL DO


            if ( r .ge. r_particle) then
                ! outside the sphere
                z_selector = 3

                call lib_mie_vector_spherical_harmonics_components(theta, phi, r, k0 * n_medium, n_range, z_selector, &
                                                                   M_nm, N_nm)
                ! eq. (5)
                !$OMP PARALLEL DO PRIVATE(n, m, buffer_real)
                do n=n_range(1), n_range(2)
                    do m=-n, n
                        if (calc_order_m%item(n)%item(m)) then
                            buffer_real = abs(e_field_0) * real((2*n+1), kind=8) &
                                          * lib_math_factorial_get_n_minus_m_divided_by_n_plus_m(n,m)

#ifdef _DEBUG_
                            if (isnan(buffer_real)) then
                                print *, "get_e_field_scattered_xu: ERROR"
                                print *, "  buffer_real is NaN"
                                print * , "  n = ", n
                                print * , "  m = ", m
                            end if
#endif
                            e_field_nm%item(n)%item(m) = buffer_real * cmplx(0,1, kind=8)**n

                            if(buffer_real .eq. 0.0) then
                                calc_order_m%item(n)%item(m) = .false.
                            end if
                        end if
                    end do
                end do
                !$OMP END PARALLEL DO

                ! first line eq. (4)
                !$OMP PARALLEL DO PRIVATE(n, m, i, buffer_cmplx)
                do n= n_range(1), n_range(2)
                    i = n - n_range(1) + 1

                    e_field_n_s(i)%rho = cmplx(0,0)
                    e_field_n_s(i)%theta = cmplx(0,0)
                    e_field_n_s(i)%phi = cmplx(0,0)
                    buffer_e_field_n_s%rho = cmplx(0,0)
                    buffer_e_field_n_s%theta = cmplx(0,0)
                    buffer_e_field_n_s%phi = cmplx(0,0)

                    h_field_n_s(i)%rho = cmplx(0,0)
                    h_field_n_s(i)%theta = cmplx(0,0)
                    h_field_n_s(i)%phi = cmplx(0,0)
                    buffer_h_field_n_s%rho = cmplx(0,0)
                    buffer_h_field_n_s%theta = cmplx(0,0)
                    buffer_h_field_n_s%phi = cmplx(0,0)

                    do m=-n, n
                        if (calc_order_m%item(n)%item(m)) then
                            buffer_cmplx = cmplx(0, 1, kind=8) * e_field_nm%item(n)%item(m)
                            buffer_e_field_n_s = buffer_cmplx * (a_n(i)*N_nm(n)%coordinate(m) * p0_nm%item(n)%item(m) &
                                                                 +b_n(i)*M_nm(n)%coordinate(m) * q0_nm%item(n)%item(m))

                            buffer_h_field_n_s = e_field_nm%item(n)%item(m) &
                                                 * (b_n(i)*N_nm(n)%coordinate(m) * p0_nm%item(n)%item(m) &
                                                    +a_n(i)*M_nm(n)%coordinate(m) * q0_nm%item(n)%item(m))
#ifdef _DEBUG_
                            if (isnan(real(buffer_e_field_n_s%rho)) .or. isnan(aimag(buffer_e_field_n_s%rho)) &
                                .or. isnan(real(buffer_e_field_n_s%phi)) .or. isnan(aimag(buffer_e_field_n_s%phi)) &
                                .or. isnan(real(buffer_e_field_n_s%theta)) .or. isnan(aimag(buffer_e_field_n_s%theta)) ) then
                                print *, "get_e_field_scattered_xu: ERROR"
                                print *, "  e_field_n_s is NaN"
                                print * , "  n = ", n
                                print * , "  m = ", m
                            end if
#endif
                            e_field_n_s(i) = e_field_n_s(i) + buffer_e_field_n_s
                            h_field_n_s(i) = h_field_n_s(i) + buffer_h_field_n_s
                        end if
                    end do
                end do
                !$OMP END PARALLEL DO

                e_field_s%theta = cmplx(0,0,kind=8)
                e_field_s%phi = cmplx(0,0,kind=8)
                e_field_s%rho = cmplx(0,0,kind=8)

                h_field_s%theta = cmplx(0,0,kind=8)
                h_field_s%phi = cmplx(0,0,kind=8)
                h_field_s%rho = cmplx(0,0,kind=8)

                do i=n_range(2)-n_range(1)+1, 1, -1
                    e_field_s = e_field_s + e_field_n_s(i)
                    h_field_s = h_field_s + h_field_n_s(i)
                end do
            else
                ! inside the sphere
                e_field_s%theta = cmplx(0,0,kind=8)
                e_field_s%phi = cmplx(0,0,kind=8)
                e_field_s%rho = cmplx(0,0,kind=8)

                h_field_s%theta = cmplx(0,0,kind=8)
                h_field_s%phi = cmplx(0,0,kind=8)
                h_field_s%rho = cmplx(0,0,kind=8)
            end if

            field_s(1) = e_field_s

            ! omega = k * v
            ! v = c0 / n_medium
            ! mu = 1 <-- definition
            !
            !   k / (omega * mu) = k / ( k * v * 1 )
            ! = k / ( k * c0 / n_medium )
            ! = n_medium / c0
            field_s(2) = h_field_s * n_medium / real(const_c0, kind=8)

        end function get_field_scattered_xu_real

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
        !       e.g. n = (/ 1, 45 /) <-- size parameter
        !   alpha: double precision, optional(std: 0)
        !       incident angle with respect to the z-axis
        !   beta: double precision, optional(std: 0)
        !       angle between the x axis and the projection of the wave vector on the x-y plane
        !
        ! Reference: Electromagnetic scattering by an aggregate of spheres, Yu-lin Xu
        !
        function get_field_scattered_xu_cmplx(theta, phi, r, &
                                              e_field_0, lambda, &
                                              n_medium, r_particle, n_particle, &
                                              n_range, &
                                              alpha, beta) &
                                            result (field_s)
            implicit none
            ! dummy
            double precision, intent(in) :: theta
            double precision, intent(in) :: phi
            double precision, intent(in) :: r
            double precision, intent(in) :: e_field_0
            double precision, intent(in) :: lambda
            double precision, intent(in) :: n_medium
            double precision, intent(in) :: r_particle
            complex(kind=8), intent(in) :: n_particle
            integer(kind=4), dimension(2) :: n_range
            double precision, intent(in), optional :: alpha
            double precision, intent(in), optional :: beta

            type(spherical_coordinate_cmplx_type), dimension(2) :: field_s

            ! auxiliary
            type(spherical_coordinate_cmplx_type) :: e_field_s
            type(spherical_coordinate_cmplx_type) :: h_field_s

            double precision :: k0 ! wave number = 2 pi / lambda
            double precision :: size_parameter !  = r_particle * k0 * n_medium ! eq. 11

            double precision :: buffer_real
            complex(kind=8) :: buffer_cmplx
            integer(kind=4) :: i
            integer(kind=4) :: m
            integer(kind=4) :: n

            complex(kind=8), dimension(n_range(2)-n_range(1)+1) :: a_n
            complex(kind=8), dimension(n_range(2)-n_range(1)+1) :: b_n


            integer(kind=1) :: z_selector

            type(list_spherical_coordinate_cmplx_type), dimension(:), allocatable :: M_nm
            type(list_spherical_coordinate_cmplx_type), dimension(:), allocatable :: N_nm

            type(list_list_cmplx) :: e_field_nm
            type(spherical_coordinate_cmplx_type), dimension(n_range(2)-n_range(1)+1) :: e_field_n_s
            type(spherical_coordinate_cmplx_type) :: buffer_e_field_n_s
            type(spherical_coordinate_cmplx_type), dimension(n_range(2)-n_range(1)+1) :: h_field_n_s
            type(spherical_coordinate_cmplx_type) :: buffer_h_field_n_s

            type(list_list_cmplx) :: p0_nm
            type(list_list_cmplx) :: q0_nm
            type(list_list_real) :: pi_nm
            type(list_list_real) :: tau_nm

            double precision :: cos_beta
            double precision :: sin_beta

            double precision :: m_alpha
            double precision :: m_beta

            type(list_list_logical) :: calc_order_m


            ! --- init ---
            ! standard values
            if (present(alpha)) then
                m_alpha = alpha
            else
                m_alpha = 0
            end if
            if (present(beta)) then
                m_beta = beta
            else
                m_beta = 0
            end if


            k0 = 2.0_8 * PI / lambda
            size_parameter = r_particle * k0 * n_medium ! eq. 11

            cos_beta = cos(m_beta)
            sin_beta = sin(m_beta)


            call init_list(e_field_nm, n_range(1), n_range(2)-n_range(1)+1)
            call init_list(p0_nm, n_range(1), n_range(2)-n_range(1)+1)
            call init_list(q0_nm, n_range(1), n_range(2)-n_range(1)+1)
            call init_list(calc_order_m, n_range(1), n_range(2)-n_range(1)+1)

            ! --- pre-calc ---
            call get_coefficients_a_b_cmplx_barberh(size_parameter, n_particle/n_medium, n_range, a_n, b_n)

            ! errata eq. (1) Equations (21) on p. 4577
            call lib_math_associated_legendre_polynomial_theta(m_alpha, n_range(2), pi_nm, tau_nm)
            !$OMP PARALLEL DO PRIVATE(n, m, buffer_real, buffer_cmplx)
            do n=n_range(1), n_range(2)
                do m=-n, n
                    buffer_real = -m * m_beta
                    buffer_cmplx = cmplx(cos(buffer_real), sin(buffer_real), kind=8)

                    buffer_real = 1.0_8 / real(n * (n + 1), kind=8)
                    buffer_cmplx = buffer_real * buffer_cmplx

                    p0_nm%item(n)%item(m) = tau_nm%item(n)%item(m) * buffer_cmplx
                    q0_nm%item(n)%item(m) = pi_nm%item(n)%item(m) * buffer_cmplx

                    if(p0_nm%item(n)%item(m) .eq. cmplx(0,0) &
                       .and. q0_nm%item(n)%item(m) .eq. cmplx(0,0)) then
                        calc_order_m%item(n)%item(m) = .false.
                    end if
                end do
            end do
            !$OMP END PARALLEL DO

            if ( r .ge. r_particle) then
                z_selector = 3

                call lib_mie_vector_spherical_harmonics_components(theta, phi, r, k0 * n_medium, n_range, z_selector, &
                                                                   M_nm, N_nm)
                ! eq. (5)
                !$OMP PARALLEL DO PRIVATE(n, buffer_real)
                do n=n_range(1), n_range(2)
                    do m=-n, n
                        if (calc_order_m%item(n)%item(m)) then
                            buffer_real = abs(e_field_0) * real((2*n+1), kind=8) &
                                          * lib_math_factorial_get_n_minus_m_divided_by_n_plus_m(n,m)

#ifdef _DEBUG_
                            if (isnan(buffer_real)) then
                                print *, "get_e_field_scattered_xu: ERROR"
                                print *, "  buffer_real is NaN"
                                print * , "  n = ", n
                                print * , "  m = ", m
                            end if
#endif
                            e_field_nm%item(n)%item(m) = buffer_real * cmplx(0,1, kind=8)**n

                            if(buffer_real .eq. 0.0) then
                                calc_order_m%item(n)%item(m) = .false.
                            end if
                        end if
                    end do
                end do
                !$OMP END PARALLEL DO

                ! first line eq. (4)
                !$OMP PARALLEL DO PRIVATE(n, i, buffer_cmplx)
                do n= n_range(1), n_range(2)
                    i = n - n_range(1) + 1

                    e_field_n_s(i)%rho = cmplx(0,0)
                    e_field_n_s(i)%theta = cmplx(0,0)
                    e_field_n_s(i)%phi = cmplx(0,0)
                    buffer_e_field_n_s%rho = cmplx(0,0)
                    buffer_e_field_n_s%theta = cmplx(0,0)
                    buffer_e_field_n_s%phi = cmplx(0,0)

                    h_field_n_s(i)%rho = cmplx(0,0)
                    h_field_n_s(i)%theta = cmplx(0,0)
                    h_field_n_s(i)%phi = cmplx(0,0)
                    buffer_h_field_n_s%rho = cmplx(0,0)
                    buffer_h_field_n_s%theta = cmplx(0,0)
                    buffer_h_field_n_s%phi = cmplx(0,0)

                    do m=-n, n
                        if (calc_order_m%item(n)%item(m)) then
                            buffer_cmplx = cmplx(0, 1, kind=8) * e_field_nm%item(n)%item(m)
                            buffer_e_field_n_s = buffer_cmplx * (a_n(i)*N_nm(n)%coordinate(m) * p0_nm%item(n)%item(m) &
                                                             +b_n(i)*M_nm(n)%coordinate(m) * q0_nm%item(n)%item(m))

                            buffer_h_field_n_s = e_field_nm%item(n)%item(m) &
                                                 * (b_n(i)*N_nm(n)%coordinate(m) * p0_nm%item(n)%item(m) &
                                                    +a_n(i)*M_nm(n)%coordinate(m) * q0_nm%item(n)%item(m))

#ifdef _DEBUG_
                            if (isnan(real(e_field_n_s(i)%rho)) .or. isnan(aimag(e_field_n_s(i)%rho)) &
                                .or. isnan(real(e_field_n_s(i)%phi)) .or. isnan(aimag(e_field_n_s(i)%phi)) &
                                .or. isnan(real(e_field_n_s(i)%theta)) .or. isnan(aimag(e_field_n_s(i)%theta)) ) then
                                print *, "get_e_field_scattered_xu: ERROR"
                                print *, "  e_field_n_s is NaN"
                                print * , "  n = ", n
                                print * , "  m = ", m
                            end if

!                            if (isinf(real(e_field_n_s(i)%rho)) .or. isinf(aimag(e_field_n_s(i)%rho)) &
!                                .or. isinf(real(e_field_n_s(i)%phi)) .or. isinf(aimag(e_field_n_s(i)%phi)) &
!                                .or. isinf(real(e_field_n_s(i)%theta)) .or. isinf(aimag(e_field_n_s(i)%theta)) ) then
!                                print *, "get_e_field_scattered_xu: ERROR"
!                                print *, "  e_field_n_s is infinity"
!                                print * , "  n = ", n
!                                print * , "  m = ", m
!                            end if
#endif
                            e_field_n_s(i) = e_field_n_s(i) + buffer_e_field_n_s
                            h_field_n_s(i) = h_field_n_s(i) + buffer_h_field_n_s
                        end if
                    end do
                end do
                !$OMP END PARALLEL DO

                e_field_s%theta = cmplx(0,0,kind=8)
                e_field_s%phi = cmplx(0,0,kind=8)
                e_field_s%rho = cmplx(0,0,kind=8)

                h_field_s%theta = cmplx(0,0,kind=8)
                h_field_s%phi = cmplx(0,0,kind=8)
                h_field_s%rho = cmplx(0,0,kind=8)

!                do i=1, n_range(2)-n_range(1)+1
                do i=n_range(2)-n_range(1)+1, 1, -1
                    e_field_s = e_field_s + e_field_n_s(i)
                    h_field_s = h_field_s + h_field_n_s(i)
                end do
            else
                e_field_s%theta = cmplx(0,0,kind=8)
                e_field_s%phi = cmplx(0,0,kind=8)
                e_field_s%rho = cmplx(0,0,kind=8)

                h_field_s%theta = cmplx(0,0,kind=8)
                h_field_s%phi = cmplx(0,0,kind=8)
                h_field_s%rho = cmplx(0,0,kind=8)
            end if

            field_s(1) = e_field_s

            ! omega = k * v
            ! v = c0 / n_medium
            ! mu = 1 <-- definition
            !
            !   k / (omega * mu) = k / ( k * v * 1 )
            ! = k / ( k * c0 / n_medium )
            ! = n_medium / c0
            field_s(2) = h_field_s * n_medium / real(const_c0, kind=8)

        end function get_field_scattered_xu_cmplx

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
        ! Reference: Absorption and Scattering of Light by Small Particles,Bohren and Huffman,  eq. (4.53)
        subroutine get_coefficients_a_b_real_bohrenh(x, m, mu, mu1, n, a_n, b_n)
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

        end subroutine get_coefficients_a_b_real_bohrenh

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
        subroutine get_coefficients_a_b_cmplx_bohrenh(x, m, mu, mu1, n, a_n, b_n)
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

        end subroutine get_coefficients_a_b_cmplx_bohrenh

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
        !   caching: logical (std: .true.)
        !
        ! Reference: Light Scattering by Particles, Barber, Hill, eq. (4.18)
        subroutine get_coefficients_a_b_real_barberh(x, m, n, a_n, b_n, caching)
            implicit none
            ! dummy
            real(kind=8), intent(in) :: x
            real(kind=8), intent(in) :: m
            integer(kind=4), dimension(2), intent(in) :: n

            complex(kind=8), dimension(n(2)-n(1)+1), intent(inout) :: a_n
            complex(kind=8), dimension(n(2)-n(1)+1), intent(inout) :: b_n

            logical, intent(in), optional :: caching

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

            integer :: cache_no

            ! 0: calculate
            ! >0: array number
            ! -1: search
            cache_no = 0

            if (cache_coefficients_a_b_real_barberh_enabled) then
                if (present(caching)) then
                    if (caching) then
                        cache_no = -1
                    end if
                else
                    cache_no = -1
                end if
            end if

            ! search
            if (cache_no .eq. -1) then
                do i=1, size(cache_coefficients_a_b_real_barberh)
                    if ((cache_coefficients_a_b_real_barberh(i)%x .eq. x) &
                        .and. (cache_coefficients_a_b_real_barberh(i)%m .eq. m) &
                        .and. (cache_coefficients_a_b_real_barberh(i)%n_max .ge. n(2)) ) then
                        cache_no = i
                        exit
                    end if
                end do
            end if

            if (cache_no .gt. 0) then
                a_n = cache_coefficients_a_b_real_barberh(cache_no)%a_n(n(1):n(2))
                b_n = cache_coefficients_a_b_real_barberh(cache_no)%b_n(n(1):n(2))
            else

                number_of_n = n(2) - n(1) + 1

                mx = m*x

                j_n_x = lib_math_bessel_spherical_first_kind(x, n(1)-1, number_of_n+1)
                h_n_x = lib_math_hankel_spherical_1(x, n(1)-1, number_of_n+1)

                An = get_An_real(x, m, n(1), number_of_n)

                !$OMP PARALLEL DO PRIVATE(i, mAn, n_div_x, buffer, numerator, denominator)
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
                !$OMP END PARALLEL DO

            end if

        end subroutine get_coefficients_a_b_real_barberh

        ! calculates the scattering coefficients
        !
        ! Arguments
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
        !   n: integer, dimension(2)
        !       first element: start index
        !       second element: last index
        !       HINT: first element .le. second element
        !
        ! Reference: Light Scattering by Particles, Barber, Hill, eq. (4.18)
        subroutine get_coefficients_a_b_cmplx_barberh(x, m, n, a_n, b_n, caching)
            implicit none
            ! dummy
            double precision :: x
            complex(kind=8) :: m
            integer(kind=4), dimension(2) :: n

            complex(kind=8), dimension(n(2)-n(1)+1) :: a_n
            complex(kind=8), dimension(n(2)-n(1)+1) :: b_n

            logical, intent(in), optional :: caching

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

            integer :: cache_no

            ! 0: calculate
            ! >0: array number
            ! -1: search
            cache_no = 0

            if (cache_coefficients_a_b_cmplx_barberh_enabled) then
                if (present(caching)) then
                    if (caching) then
                        cache_no = -1
                    end if
                else
                    cache_no = -1
                end if
            end if

            ! search
            if (cache_no .eq. -1) then
                do i=1, size(cache_coefficients_a_b_cmplx_barberh)
                    if ((cache_coefficients_a_b_cmplx_barberh(i)%x .eq. x) &
                        .and. (cache_coefficients_a_b_cmplx_barberh(i)%m .eq. m) &
                        .and. (cache_coefficients_a_b_cmplx_barberh(i)%n_max .ge. n(2)) ) then
                        cache_no = i
                        exit
                    end if
                end do
            end if

            if (cache_no .gt. 0) then
                a_n = cache_coefficients_a_b_cmplx_barberh(cache_no)%a_n(n(1):n(2))
                b_n = cache_coefficients_a_b_cmplx_barberh(cache_no)%b_n(n(1):n(2))
            else

                number_of_n = n(2) - n(1) + 1

                mx = m*x

                j_n_x = lib_math_bessel_spherical_first_kind(x, n(1)-1, number_of_n+1)
                h_n_x = lib_math_hankel_spherical_1(x, n(1)-1, number_of_n+1)

                An = get_An_cmplx(x, m, n(1), number_of_n)

                !$OMP PARALLEL DO PRIVATE(i, mAn, n_div_x, An_div_m, buffer, numerator, denominator)
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
                !$OMP END PARALLEL DO
            end if

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

        ! Calculates the coefficients (of vector spherical components) of a plane incident wave
        ! - transverse magnetic mode (TM)
        !
        !                 z
        !                 ^
        !             K_j |
        !                 --> x
        !                ^
        !               /
        !           z  /d_j_0
        !           ^ /
        !       K_0 |/
        !           --> x
        !      ____________
        !      ____k^______
        !      _____|______
        !      ____________ plane wave
        !
        ! Argument
        ! ----
        !   k: type(cartesian_coordinate_real_type)
        !       wave vector
        !   d_0_j: type(cartesian_coordinate_real_type)
        !       vector from the
        !   n_range: integer, dimension(2)
        !       first and last index (degree) of the sum to calculate the electical field
        !       e.g. n = (/ 1, 45 /) <-- size parameter
        !
        ! Returns
        ! ----
        !   p: type(list_list_cmplx)
        !       coefficient of vector spherical componets
        !   q: type(list_list_cmplx)
        !       coefficient of vector spherical componets
        !
        ! Reference: Electromagnetic scattering by an aggregate of spheres Yu-lin Xu, eq. 20
        subroutine lib_mie_get_p_q_j_j(k, d_0_j, n_range, p, q)
            implicit none
            ! dummy
            type(cartesian_coordinate_real_type), intent(in) :: k
            type(cartesian_coordinate_real_type), intent(in) :: d_0_j
            integer(kind=4), dimension(2),intent(in) :: n_range

            type(list_list_cmplx), intent(inout) :: p
            type(list_list_cmplx), intent(inout) :: q

            ! auxiliary
            double precision :: alpha
            double precision :: beta

            type(spherical_coordinate_real_type) :: k_spherical

            double precision :: dot_k_d

            complex(kind=8) :: exp_k_d

            ! --- pre-calc ---
            k_spherical = k

            alpha = k_spherical%theta
            beta = k_spherical%phi

            if (abs(d_0_j) .gt. 0.0) then
                ! j-th coordinate system
                dot_k_d = dot_product(k, d_0_j)
                exp_k_d = cmplx(cos(dot_k_d), sin(dot_k_d), kind=8)
                call lib_mie_get_p_q_j_j_core(alpha, beta, n_range, p, q, exp_k_d)
            else
                ! 0-th coordinate system
                call lib_mie_get_p_q_j_j_core(alpha, beta, n_range, p, q)
            end if

        end subroutine

        ! Calculates the coefficients (of vector spherical components) of a plane incident wave
        ! - transverse magnetic mode (TM)
        !
        !                 z
        !                 ^
        !             K_j |
        !                 --> x
        !                ^
        !               /
        !           z  /d_j_0
        !           ^ /
        !       K_0 |/
        !           --> x
        !      ____________
        !      ____k^______
        !      _____|______
        !      ____________ plane wave
        !
        ! Argument
        ! ----
        !   alpha: double precision
        !       polar angle [0..Pi)
        !       if 0: propagation direction along the z-axis
        !   beta: double precision
        !       azimuthal angle [0..2Pi)
        !       if 0 and alpha = 0: The E-field oscillates along the x-axis.
        !   n_range: integer, dimension(2)
        !       first and last index (degree) of the sum to calculate the electical field
        !       e.g. n = (/ 1, 45 /) <-- size parameter
        !   exp_k_d: complex, optional (std: 1.0, 0-th coordinate system)
        !       formula: exp(i dot(k, d_0_j))
        !           k: wave vector
        !           d_0_j: vector from 0-th to j-th coordinate system
        !   caching: logical (std: .true.)
        !
        ! Returns
        ! ----
        !   p: type(list_list_cmplx)
        !
        !   q: type(list_list_cmplx)
        !
        ! Reference: Electromagnetic scattering by an aggregate of spheres Yu-lin Xu, eq. 20
        subroutine lib_mie_get_p_q_j_j_core(alpha, beta, n_range, p, q, exp_k_d, caching)
            implicit none
            ! dummy
            double precision, intent(in) :: alpha
            double precision, intent(in) :: beta
            integer(kind=4), dimension(2),intent(in) :: n_range

            type(list_list_cmplx), intent(inout) :: p
            type(list_list_cmplx), intent(inout) :: q

            complex(kind=8), intent(in), optional :: exp_k_d

            logical, optional :: caching

            ! auxiliary
            integer :: i
            double precision :: buffer_real_n
            integer(kind=4) :: m
            integer(kind=4) :: n

            type(list_list_real) :: pi_nm
            type(list_list_real) :: tau_nm

            double precision :: sin_beta
            double precision :: cos_beta

            logical :: m_caching
            ! 0: p_0 and q_0 are not pre calculated
            ! >0: element number of the cache array
            integer :: cache_no

            if (cache_coefficients_p_0_q_0_enabled) then
                if (present(caching)) then
                    m_caching = caching
                else
                    m_caching = .true.
                end if
            else
                m_caching = .false.
            end if

            do i=1, size(cache_coefficients_p_0_q_0)
                if ((cache_coefficients_p_0_q_0(i)%alpha .eq. alpha) &
                    .and. (cache_coefficients_p_0_q_0(i)%beta .eq. beta) &
                    .and. (cache_coefficients_p_0_q_0(i)%n_max .le. n_range(2))) then
                    cache_no = i
                else
                    cache_no = 0
                end if
            end do

            ! --- init ---
            if (alpha .eq. 0.d0) then
                allocate(p%item(n_range(1):n_range(2)))
                allocate(q%item(n_range(1):n_range(2)))

                do i=n_range(1), n_range(2)
                    allocate(p%item(i)%item(-1:1))
                    allocate(q%item(i)%item(-1:1))
                end do
            else
                call init_list(p, n_range(1), n_range(2)-n_range(1)+1)
                call init_list(q, n_range(1), n_range(2)-n_range(1)+1)
            end if

            ! --- pre-calc ---
            if (cache_no .gt. 0) then
                sin_beta = sin(beta)
                cos_beta = cos(beta)

                call lib_math_associated_legendre_polynomial_theta(alpha, n_range(2), pi_nm, tau_nm)
            end if

            ! errata eq. (1) Equations (21) on p. 4577
            do n=n_range(1), n_range(2)
                buffer_real_n = 1.0d0 / real(n*(n+1), kind=8)
                if (alpha .eq. 0.d0) then
                    p%item(n)%item(:) = 0.d0
                    q%item(n)%item(:) = 0.d0
                    call get_value(-1_4, n, p%item(n)%item(-1), q%item(n)%item(-1))
                    call get_value(1_4, n, p%item(n)%item(1), q%item(n)%item(1))
                else
                    do m=-n, n
                        call get_value(m, n, p%item(n)%item(m), q%item(n)%item(m))
                    end do
                end if
            end do

            contains
                ! errata eq. (1) Equations (21) on p. 4577
                subroutine get_value(m, n, p, q)
                    implicit none
                    ! dummy
                    integer(kind=4), intent(in) :: m
                    integer(kind=4), intent(in) :: n

                    complex(kind=lib_math_type_kind), intent(inout) :: p
                    complex(kind=lib_math_type_kind), intent(inout) :: q

                    ! auxiliaray
                    double precision :: buffer_real
                    double complex :: buffer_cmplx

                    if (cache_no .gt. 0) then
                        p = cache_coefficients_p_0_q_0(cache_no)%p_0%item(m)%item(n)
                        q = cache_coefficients_p_0_q_0(cache_no)%q_0%item(n)%item(n)
                    else
                        buffer_real = -m * beta
                        buffer_cmplx = cmplx(cos(buffer_real), sin(buffer_real), kind=8)

                        buffer_cmplx = buffer_real_n * buffer_cmplx

                        p = tau_nm%item(n)%item(m) * buffer_cmplx
                        q = pi_nm%item(n)%item(m) * buffer_cmplx
                    end if

                    if (present(exp_k_d)) then
                        p = exp_k_d * p
                        q = exp_k_d * q
                    end if

                end subroutine

        end subroutine lib_mie_get_p_q_j_j_core

        ! Reference: , eq. 30
        subroutine lib_mie_get_interactive_scattering_coefficients(k, lambda, n_medium, &
                                                                   sphere, sphere_parameter, sphere_j, &
                                                                   a, b)
            implicit none
            ! dummy
            type(spherical_coordinate_real_type), intent(in) :: k
            double precision, intent(in) :: lambda
            double precision, intent(in) :: n_medium
            type(sphere_type), dimension(:), intent(in) :: sphere
            type(sphere_parameter_type), dimension(:), intent(in) :: sphere_parameter
            integer :: sphere_j
            type(list_list_cmplx), intent(inout) :: a
            type(list_list_cmplx), intent(inout) :: b

            ! auxiliary
            integer :: m
            integer :: n
            integer :: mu
            integer :: nu
            integer, dimension(2) :: n_range
            integer, dimension(2) :: nu_range

            integer :: l
            integer :: i

            type(list_list_cmplx) :: p
            type(list_list_cmplx) :: q
            type(list_cmplx) :: a_j_n
            type(list_cmplx) :: b_j_n
            type(list_cmplx) :: a_l_n
            type(list_cmplx) :: b_l_n

            i = sphere(sphere_j)%sphere_parameter_index
            n_range = sphere_parameter(i)%n_range

            do l=ubound(sphere, 1), lbound(sphere, 1)
                if (l .ne. sphere_j) then
                    i = sphere(l)%sphere_parameter_index
                    nu_range = sphere_parameter(i)%n_range
                    a_l_n = sphere_parameter(i)%a_n
                    b_l_n = sphere_parameter(i)%b_n

                    if (n .le. nu_range(2)) then
                        do n=n_range(1), n_range(2)
                            do m=-n, n


                            end do
                        end do
                    end if
                end if
            end do

!            contains
!                subroutine calc_summand(a_munu, a_munumn, b_munu, b_munumn, )
!                    implicit none
!
!                end subroutine calc_summand

        end subroutine lib_mie_get_interactive_scattering_coefficients

        subroutine lib_mie_set_sphere_parameter(lambda, n_medium, &
                                                r_particle, n_particle,&
                                                n_range, sphere_parameter)
            implicit none
            ! dummy
            double precision, intent(in) :: lambda
            double precision, intent(in) :: n_medium
            double precision, intent(in) :: r_particle
            double complex, intent(in) :: n_particle
            integer, dimension(2) :: n_range

            type(sphere_parameter_type), intent(inout) :: sphere_parameter

            ! auxiliary
            double precision :: size_parameter

            size_parameter = 2 * PI * n_medium * r_particle / lambda

            sphere_parameter%n_range = n_range
            sphere_parameter%size_parameter = size_parameter
            sphere_parameter%radius = r_particle
            sphere_parameter%refractive_index = n_particle

            if (aimag(n_particle) .eq. 0d0) then
                call get_coefficients_a_b_real_barberh(size_parameter, real(n_particle)/n_medium, n_range, &
                                                       sphere_parameter%a_n%item, sphere_parameter%b_n%item)
            else
                call get_coefficients_a_b_cmplx_barberh(size_parameter, n_particle/n_medium, n_range, &
                                                           sphere_parameter%a_n%item, sphere_parameter%b_n%item)
            end if
        end subroutine

        function lib_mie_scattering_by_a_sphere_test_functions() result (rv)
            implicit none
            ! dummy
            integer :: rv

            ! auxiliaray
            double precision, parameter :: ground_truth_e = 10.0_8**(-6.0_8)
            ! CPU-time
            real :: test_start, test_finish
            ! WALL-time
            INTEGER :: test_count_start, test_count_finish, test_count_rate

            rv = 0

            call system_clock(test_count_start, test_count_rate)
            call cpu_time(test_start)

            if (.not. test_get_coefficients_a_b_real_bohrenh()) then
                rv = rv + 1
            end if
!            if (.not. test_get_coefficients_a_b_cmplx_bohrenh()) then
!                rv = rv + 1
!            end if
            if (.not. test_get_coefficients_a_b_real_barberh()) then
                rv = rv + 1
            end if
            if (.not. test_get_coefficients_a_b_cmplx_barberh()) then
                rv = rv + 1
            end if
!            if (.not. test_get_e_field_scattered_real()) then
!                rv = rv + 1
!            end if
            if (.not. test_get_field_scattered_plane_section_real()) then
                rv = rv + 1
            end if
            if (.not. test_get_field_scattered_plane_section_cmplx()) then
                rv = rv + 1
            end if

            call cpu_time(test_finish)
            call system_clock(test_count_finish, test_count_rate)

            print *, ""
            print *, "------lib_mie_scattering_by_a_sphere_test_functions------"
            print '("  CPU-Time = ",f10.3," seconds.")',test_finish-test_start
            print '("  WALL-Time = ",f10.3," seconds.")',(test_count_finish-test_count_start) / real(test_count_rate)
            print *, ""
            if (rv == 0) then
                print *, "lib_mie_scattering_by_a_sphere_test_functions tests: OK"
            else
                print *, rv,"lib_mie_scattering_by_a_sphere_test_functions test(s) FAILED"
            end if
            print *, "---------------------------------------------------------"
            print *, ""

            contains

                function test_get_coefficients_a_b_real_bohrenh() result (rv)
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

                    call get_coefficients_a_b_real_bohrenh(x, m, mu, mu1, n, a_n, b_n)

                    rv = .true.
                    print *, "test_get_coefficients_a_b_real_bohrenh:"
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

                end function test_get_coefficients_a_b_real_bohrenh

                function test_get_coefficients_a_b_cmplx_bohrenh() result (rv)
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


                    call get_coefficients_a_b_cmplx_bohrenh(x, m, mu, mu1, n, a_n, b_n)

                    rv = .true.
                    print *, "test_get_coefficients_a_b_cmplx_bohrenh:"
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

                end function test_get_coefficients_a_b_cmplx_bohrenh

                function test_get_coefficients_a_b_real_barberh() result (rv)
                    implicit none
                    ! dummy
                    logical :: rv

                    ! auxiliary
                    integer(kind=4) :: i
                    real(kind=8) :: buffer

                    real(kind=8) :: x
                    real(kind=8) :: m
!                    double precision :: mu
!                    double precision :: mu1
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

                    double precision :: x
                    complex(kind=8) :: m
!                    double precision :: mu
!                    double precision :: mu1
                    integer(kind=4), dimension(2), parameter :: n = (/1, 4/)

                    complex(kind=8), dimension(n(2)-n(1)+1) :: a_n
                    complex(kind=8), dimension(n(2)-n(1)+1) :: b_n

                    complex(kind=8), dimension(n(2)-n(1)+1) :: ground_truth_a_n
                    complex(kind=8), dimension(n(2)-n(1)+1) :: ground_truth_b_n

                    m = cmplx(1.28, -1.37, kind=8)
                    x= 10.0_8

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

                function test_get_e_field_scattered_real() result (rv)
                    use file_io
                    implicit none
                    ! dummy
                    logical :: rv

                    ! parameter
                    double precision, parameter :: start_angle = 0
                    double precision, parameter :: stop_angle = 180
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
                    type(spherical_coordinate_cmplx_type), dimension(number_of_values, 2) :: field_s
                    type(cartesian_coordinate_cmplx_type), dimension(number_of_values) :: e_field_s_cart
                    type(cartesian_coordinate_cmplx_type), dimension(number_of_values) :: h_field_s_cart
                    type(cartesian_coordinate_cmplx_type), dimension(number_of_values) :: poynting_s

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
                        field_s(i, 1:2) = get_field_scattered_xu_real(theta, phi, r, e_field_0, lambda, n_medium, &
                                                                r_particle, n_particle, n_range)
                        e_field_s_cart(i) = make_cartesian(field_s(i, 1), theta, phi)
                        h_field_s_cart(i) = make_cartesian(field_s(i, 2), theta, phi)

                        ! calculate the Poynting vector: S = E x H
                        poynting_s(i) = cross_product(e_field_s_cart(i), h_field_s_cart(i))
                    end do

                    ! write to csv
                    header(1) = "degree"
                    header(2) = "i_field_s"
                    u = 99
                    open(unit=u, file="poynting_s.csv", status='unknown')
                    rv = write_csv(u, header, degree_list, real(poynting_s%x), &
                                                           real(poynting_s%y), &
                                                           real(poynting_s%z))
                    close(u)

!                    ! --- spherical components ---
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

                end function test_get_e_field_scattered_real

                function test_get_field_scattered_plane_section_real() result (rv)
                    use file_io
                    implicit none
                    ! dummy
                    logical :: rv

                    ! parameter

                    ! auxiliary
                    integer :: i
                    integer :: ii
                    double precision :: x
                    double precision :: y
                    double precision :: z
                    integer :: u

                    double precision :: lambda
                    double precision :: k0
                    double precision :: e_field_0
                    double precision :: r_particle
                    double precision :: n_particle
                    double precision :: n_medium

                    integer(kind=4), dimension(2) :: n_range

                    double precision, dimension(2) :: x_range
                    double precision, dimension(2) :: z_range
                    real(kind=8) :: step_size

                    integer :: no_x_values
                    integer :: no_z_values

                    type(spherical_coordinate_cmplx_type), dimension(2) :: buffer_field
                    type(cartesian_coordinate_real_type) :: point_cartesian
                    type(spherical_coordinate_real_type) :: point_spherical

                    type(cartesian_coordinate_cmplx_type), dimension(:, :), allocatable :: e_field_s
                    double precision, dimension(:, :), allocatable :: e_field_s_real_x
                    double precision, dimension(:, :), allocatable :: e_field_s_real_y
                    double precision, dimension(:, :), allocatable :: e_field_s_real_z

                    type(cartesian_coordinate_cmplx_type), dimension(:, :), allocatable :: h_field_s
                    double precision, dimension(:, :), allocatable :: h_field_s_real_x
                    double precision, dimension(:, :), allocatable :: h_field_s_real_y
                    double precision, dimension(:, :), allocatable :: h_field_s_real_z

                    type(cartesian_coordinate_real_type), dimension(:, :), allocatable :: poynting_s
                    type(cartesian_coordinate_cmplx_type) :: buffer_cartesian_cmplx
                    double precision, dimension(:, :), allocatable :: poynting_s_abs

                    ! CPU-time
                    real :: start, finish
                    ! WALL-time
                    INTEGER :: count_start, count_finish, count_rate

                    x_range = (/ -5.0_8 * unit_mu, 5.0_8 * unit_mu /)
                    z_range = (/ -5_8 * unit_mu, 10.0_8 * unit_mu /)
!                    step_size = 0.02_8 * unit_mu
                    step_size = 0.05_8 * unit_mu

                    no_x_values = abs(int(floor((x_range(2)-x_range(1))/step_size)))
                    no_z_values = abs(int(floor((z_range(2)-z_range(1))/step_size)))

                    allocate(e_field_s(no_x_values, no_z_values))
                    allocate(e_field_s_real_x(no_x_values, no_z_values))
                    allocate(e_field_s_real_y(no_x_values, no_z_values))
                    allocate(e_field_s_real_z(no_x_values, no_z_values))

                    allocate(h_field_s(no_x_values, no_z_values))
                    allocate(h_field_s_real_x(no_x_values, no_z_values))
                    allocate(h_field_s_real_y(no_x_values, no_z_values))
                    allocate(h_field_s_real_z(no_x_values, no_z_values))

                    allocate(poynting_s(no_x_values, no_z_values))
                    allocate(poynting_s_abs(no_x_values, no_z_values))

                    x = 0
                    y = 0
                    z = 0

                    r_particle = 2 * unit_mu

                    e_field_0 = 1
                    lambda = 0.7 * unit_mu

                    n_particle = 1.5
                    n_medium = 1

                    k0 = 2 * PI / lambda

                    n_range(1) = 1
                    n_range(2) = get_n_c(r_particle * k0)! * n_particle)
                    if (n_range(2) .gt. 45) then
                        print *, "WARNING: max degree (45) reached: ", n_range(2)
                        n_range(2) = 45
                    else
                        print *, "NOTE: max degree = ", n_range(2)
                    end if

                    call lib_math_factorial_initialise_caching(n_range(2))

                    call lib_mie_scattering_by_a_sphere_init((/ r_particle * k0 * n_medium /), &
                                                             (/ n_particle/n_medium /), &
                                                             (/ n_range(2) /))

                    call system_clock(count_start, count_rate)
                    call cpu_time(start)
                    !$OMP PARALLEL DO PRIVATE(i, ii, point_cartesian, point_spherical, buffer_field) &
                    !$OMP  FIRSTPRIVATE(x, y, z)
                    do i=1, no_x_values
                        x = x_range(1) + (i-1) * step_size
                        do ii=1, no_z_values
                            z = z_range(1) + (ii-1) * step_size

                            point_cartesian%x = x
                            point_cartesian%y = y
                            point_cartesian%z = z

                            point_spherical = point_cartesian

                            buffer_field = get_field_scattered_xu_real(point_spherical%theta, point_spherical%phi, &
                                                                       point_spherical%rho, &
                                                                       e_field_0, lambda, n_medium, &
                                                                       r_particle, n_particle, n_range)!, &
!                                                                       beta=PI/2.0_8)
                            e_field_s(i, ii) = make_cartesian(buffer_field(1), point_spherical%theta, point_spherical%phi)
                            h_field_s(i, ii) = make_cartesian(buffer_field(2), point_spherical%theta, point_spherical%phi)

                            e_field_s_real_x(i, ii) = real(e_field_s(i, ii)%x)
                            e_field_s_real_y(i, ii) = real(e_field_s(i, ii)%y)
                            e_field_s_real_z(i, ii) = real(e_field_s(i, ii)%z)

                            h_field_s_real_x(i, ii) = real(h_field_s(i, ii)%x)
                            h_field_s_real_y(i, ii) = real(h_field_s(i, ii)%y)
                            h_field_s_real_z(i, ii) = real(h_field_s(i, ii)%z)

                            ! calculate the Poynting vector: S = E x H*
                            ! eq. 43
                            h_field_s(i, ii)%x = conjg(h_field_s(i, ii)%x)
                            h_field_s(i, ii)%y = conjg(h_field_s(i, ii)%y)
                            h_field_s(i, ii)%z = conjg(h_field_s(i, ii)%z)

                            buffer_cartesian_cmplx = cross_product(e_field_s(i, ii), h_field_s(i, ii))

                            poynting_s(i, ii)%x = real(buffer_cartesian_cmplx%x) / 2.0_8
                            poynting_s(i, ii)%y = real(buffer_cartesian_cmplx%y) / 2.0_8
                            poynting_s(i, ii)%z = real(buffer_cartesian_cmplx%z) / 2.0_8

                            poynting_s_abs(i, ii) = abs(poynting_s(i, ii))
                        end do
                        print *, "  x-Value: ", x
                    end do
                    !$OMP END PARALLEL DO
                    call cpu_time(finish)
                    call system_clock(count_finish, count_rate)
                    print *, "test_get_e_field_scattered_plane_section_real"
                    print '("  CPU-Time = ",f10.3," seconds.")',finish-start

                    print '("  WALL-Time = ",f10.3," seconds.")',(count_finish-count_start) / real(count_rate)

                    ! --- wirte to PPM ---
                    ! e field
                    u = 99
                    open(unit=u, file="temp/real/e_field_s_x.ppm", status='unknown')
                    rv = write_ppm_p3(u, e_field_s_real_x)
                    close(u)

                    open(unit=u, file="temp/real/e_field_s_y.ppm", status='unknown')
                    rv = write_ppm_p3(u, e_field_s_real_y)
                    close(u)

                    open(unit=u, file="temp/real/e_field_s_z.ppm", status='unknown')
                    rv = write_ppm_p3(u, e_field_s_real_z)
                    close(u)

                    u = 99
                    open(unit=u, file="temp/real/e_field_s_x_log.ppm", status='unknown')
                    rv = write_ppm_p3(u, e_field_s_real_x, logarithmic = .true.)
                    close(u)

                    open(unit=u, file="temp/real/e_field_s_y_log.ppm", status='unknown')
                    rv = write_ppm_p3(u, e_field_s_real_y, logarithmic = .true.)
                    close(u)

                    open(unit=u, file="temp/real/e_field_s_z_log.ppm", status='unknown')
                    rv = write_ppm_p3(u, e_field_s_real_z, logarithmic = .true.)
                    close(u)

                    ! h field
                    u = 99
                    open(unit=u, file="temp/real/h_field_s_x.ppm", status='unknown')
                    rv = write_ppm_p3(u, h_field_s_real_x)
                    close(u)

                    open(unit=u, file="temp/real/h_field_s_y.ppm", status='unknown')
                    rv = write_ppm_p3(u, h_field_s_real_y)
                    close(u)

                    open(unit=u, file="temp/real/h_field_s_z.ppm", status='unknown')
                    rv = write_ppm_p3(u, h_field_s_real_z)
                    close(u)

                    u = 99
                    open(unit=u, file="temp/real/h_field_s_x_log.ppm", status='unknown')
                    rv = write_ppm_p3(u, h_field_s_real_x, logarithmic = .true.)
                    close(u)

                    open(unit=u, file="temp/real/h_field_s_y_log.ppm", status='unknown')
                    rv = write_ppm_p3(u, h_field_s_real_y, logarithmic = .true.)
                    close(u)

                    open(unit=u, file="temp/real/h_field_s_z_log.ppm", status='unknown')
                    rv = write_ppm_p3(u, h_field_s_real_z, logarithmic = .true.)
                    close(u)

                    ! Poynting
                    u = 99
                    open(unit=u, file="temp/real/poynting_s_x.ppm", status='unknown')
                    rv = write_ppm_p3(u, poynting_s(:,:)%x)
                    close(u)

                    open(unit=u, file="temp/real/poynting_s_y.ppm", status='unknown')
                    rv = write_ppm_p3(u, poynting_s(:,:)%y)
                    close(u)

                    open(unit=u, file="temp/real/poynting_s_z.ppm", status='unknown')
                    rv = write_ppm_p3(u, poynting_s(:,:)%z)
                    close(u)

                    u = 99
                    open(unit=u, file="temp/real/poynting_s_x_log.ppm", status='unknown')
                    rv = write_ppm_p3(u, poynting_s(:,:)%x, logarithmic = .true.)
                    close(u)

                    open(unit=u, file="temp/real/poynting_s_y_log.ppm", status='unknown')
                    rv = write_ppm_p3(u, poynting_s(:,:)%y, logarithmic = .true.)
                    close(u)

                    open(unit=u, file="temp/real/poynting_s_z_log.ppm", status='unknown')
                    rv = write_ppm_p3(u, poynting_s(:,:)%z, logarithmic = .true.)
                    close(u)

                    ! Poynting abs
                    open(unit=u, file="temp/real/poynting_s_abs.ppm", status='unknown')
                    rv = write_ppm_p3(u, poynting_s_abs)
                    close(u)
                    open(unit=u, file="temp/real/poynting_s_abs_log.ppm", status='unknown')
                    rv = write_ppm_p3(u, poynting_s_abs, logarithmic = .true.)
                    close(u)

                end function test_get_field_scattered_plane_section_real

                function test_get_field_scattered_plane_section_cmplx() result (rv)
                    use file_io
                    implicit none
                    ! dummy
                    logical :: rv

                    ! parameter

                    ! auxiliary
                    integer :: i
                    integer :: ii
                    double precision :: x
                    double precision :: y
                    double precision :: z
                    integer :: u

                    double precision :: lambda
                    double precision :: k0
                    double precision :: e_field_0
                    double precision :: r_particle
                    complex(kind=8) :: n_particle
                    double precision :: n_medium

                    integer(kind=4), dimension(2) :: n_range

                    double precision, dimension(2) :: x_range
                    double precision, dimension(2) :: z_range
                    real(kind=8) :: step_size

                    integer :: no_x_values
                    integer :: no_z_values

                    type(spherical_coordinate_cmplx_type), dimension(2) :: buffer
                    type(cartesian_coordinate_real_type) :: point_cartesian
                    type(spherical_coordinate_real_type) :: point_spherical

                    type(cartesian_coordinate_cmplx_type), dimension(:, :), allocatable :: e_field_s
                    double precision, dimension(:, :), allocatable :: e_field_s_real_x
                    double precision, dimension(:, :), allocatable :: e_field_s_real_y
                    double precision, dimension(:, :), allocatable :: e_field_s_real_z

                    type(cartesian_coordinate_cmplx_type), dimension(:, :), allocatable :: h_field_s
                    double precision, dimension(:, :), allocatable :: h_field_s_real_x
                    double precision, dimension(:, :), allocatable :: h_field_s_real_y
                    double precision, dimension(:, :), allocatable :: h_field_s_real_z

                    type(cartesian_coordinate_real_type), dimension(:, :), allocatable :: poynting_s
                    type(cartesian_coordinate_cmplx_type) :: buffer_cartesian_cmplx
                    double precision, dimension(:, :), allocatable :: poynting_s_abs

                    ! CPU-time
                    real :: start, finish
                    ! WALL-time
                    INTEGER :: count_start, count_finish, count_rate

                    x_range = (/ -5.0_8 * unit_mu, 5.0_8 * unit_mu /)
                    z_range = (/ -5.0_8 * unit_mu, 10.0_8 * unit_mu /)
!                    step_size = 0.02_8 * unit_mu
                    step_size = 0.05_8 * unit_mu

                    no_x_values = abs(int(floor((x_range(2)-x_range(1))/step_size)))
                    no_z_values = abs(int(floor((z_range(2)-z_range(1))/step_size)))

                    allocate(e_field_s(no_x_values, no_z_values))
                    allocate(e_field_s_real_x(no_x_values, no_z_values))
                    allocate(e_field_s_real_y(no_x_values, no_z_values))
                    allocate(e_field_s_real_z(no_x_values, no_z_values))

                    allocate(h_field_s(no_x_values, no_z_values))
                    allocate(h_field_s_real_x(no_x_values, no_z_values))
                    allocate(h_field_s_real_y(no_x_values, no_z_values))
                    allocate(h_field_s_real_z(no_x_values, no_z_values))

                    allocate(poynting_s(no_x_values, no_z_values))
                    allocate(poynting_s_abs(no_x_values, no_z_values))

                    x = 0
                    y = 0
                    z = 0

                    r_particle = 2 * unit_mu

                    e_field_0 = 1
                    lambda = 0.7 * unit_mu

                    ! https://refractiveindex.info/?shelf=main&book=Ag&page=Johnson
                    n_particle = cmplx(0.040000, 7.1155, kind=8)
!                    n_particle = cmplx(1.5, 0, kind=8)
                    n_medium = 1

                    k0 = 2 * PI / lambda

                    n_range(1) = 1
!                    n_range(2) = min(45, get_n_c(r_particle * k0 * abs(n_particle))) ! todo: abs??
                    n_range(2) = get_n_c(r_particle * k0)
                    if (n_range(2) .gt. 45) then
                        print *, "WARNING: max degree (45) reached: ", n_range(2)
                        n_range(2) = 45
                    else
                        print *, "NOTE: max degree = ", n_range(2)
                    end if

                    call system_clock(count_start, count_rate)
                    call cpu_time(start)
                    !$OMP PARALLEL DO PRIVATE(i, ii, point_cartesian, point_spherical, buffer) &
                    !$OMP  FIRSTPRIVATE(x, y, z)
                    do i=1, no_x_values
                        x = x_range(1) + (i-1) * step_size
                        do ii=1, no_z_values
                            z = z_range(1) + (ii-1) * step_size

                            point_cartesian%x = x
                            point_cartesian%y = y
                            point_cartesian%z = z

                            point_spherical = point_cartesian

                            buffer = get_field_scattered_xu_cmplx(point_spherical%theta, point_spherical%phi, &
                                                                  point_spherical%rho, &
                                                                  e_field_0, lambda, n_medium, &
                                                                  r_particle, n_particle, n_range)

                            e_field_s(i, ii) = make_cartesian(buffer(1), point_spherical%theta, point_spherical%phi)
                            h_field_s(i, ii) = make_cartesian(buffer(2), point_spherical%theta, point_spherical%phi)

                            e_field_s_real_x(i, ii) = real(e_field_s(i, ii)%x)
                            e_field_s_real_y(i, ii) = real(e_field_s(i, ii)%y)
                            e_field_s_real_z(i, ii) = real(e_field_s(i, ii)%z)

                            h_field_s_real_x(i, ii) = real(h_field_s(i, ii)%x)
                            h_field_s_real_y(i, ii) = real(h_field_s(i, ii)%y)
                            h_field_s_real_z(i, ii) = real(h_field_s(i, ii)%z)

                            ! calculate the Poynting vector: S = E x H
                            ! eq. 43
                            h_field_s(i, ii)%x = conjg(h_field_s(i, ii)%x)
                            h_field_s(i, ii)%y = conjg(h_field_s(i, ii)%y)
                            h_field_s(i, ii)%z = conjg(h_field_s(i, ii)%z)

                            buffer_cartesian_cmplx = cross_product(e_field_s(i, ii), h_field_s(i, ii))

                            poynting_s(i, ii)%x = real(buffer_cartesian_cmplx%x) / 2.0_8
                            poynting_s(i, ii)%y = real(buffer_cartesian_cmplx%y) / 2.0_8
                            poynting_s(i, ii)%z = real(buffer_cartesian_cmplx%z) / 2.0_8

                            poynting_s_abs(i, ii) = abs(poynting_s(i, ii))
                        end do
                        print *, "  x-Value: ", x
                    end do
                    !$OMP END PARALLEL DO
                    call cpu_time(finish)
                    call system_clock(count_finish, count_rate)
                    print *, "test_get_e_field_scattered_plane_section_real"
                    print '("  CPU-Time = ",f10.3," seconds.")',finish-start

                    print '("  WALL-Time = ",f10.3," seconds.")',(count_finish-count_start) / real(count_rate)

                    ! --- wirte to PPM ---
                    ! e field
                    u = 99
                    open(unit=u, file="temp/cmplx/e_field_s_x.ppm", status='unknown')
                    rv = write_ppm_p3(u, e_field_s_real_x)
                    close(u)

                    open(unit=u, file="temp/cmplx/e_field_s_y.ppm", status='unknown')
                    rv = write_ppm_p3(u, e_field_s_real_y)
                    close(u)

                    open(unit=u, file="temp/cmplx/e_field_s_z.ppm", status='unknown')
                    rv = write_ppm_p3(u, e_field_s_real_z)
                    close(u)

                    u = 99
                    open(unit=u, file="temp/cmplx/e_field_s_x_log.ppm", status='unknown')
                    rv = write_ppm_p3(u, e_field_s_real_x, logarithmic = .true.)
                    close(u)

                    open(unit=u, file="temp/cmplx/e_field_s_y_log.ppm", status='unknown')
                    rv = write_ppm_p3(u, e_field_s_real_y, logarithmic = .true.)
                    close(u)

                    open(unit=u, file="temp/cmplx/e_field_s_z_log.ppm", status='unknown')
                    rv = write_ppm_p3(u, e_field_s_real_z, logarithmic = .true.)
                    close(u)

                    ! h field
                    u = 99
                    open(unit=u, file="temp/cmplx/h_field_s_x.ppm", status='unknown')
                    rv = write_ppm_p3(u, h_field_s_real_x)
                    close(u)

                    open(unit=u, file="temp/cmplx/h_field_s_y.ppm", status='unknown')
                    rv = write_ppm_p3(u, h_field_s_real_y)
                    close(u)

                    open(unit=u, file="temp/cmplx/h_field_s_z.ppm", status='unknown')
                    rv = write_ppm_p3(u, h_field_s_real_z)
                    close(u)

                    u = 99
                    open(unit=u, file="temp/cmplx/h_field_s_x_log.ppm", status='unknown')
                    rv = write_ppm_p3(u, h_field_s_real_x, logarithmic = .true.)
                    close(u)

                    open(unit=u, file="temp/cmplx/h_field_s_y_log.ppm", status='unknown')
                    rv = write_ppm_p3(u, h_field_s_real_y, logarithmic = .true.)
                    close(u)

                    open(unit=u, file="temp/cmplx/h_field_s_z_log.ppm", status='unknown')
                    rv = write_ppm_p3(u, h_field_s_real_z, logarithmic = .true.)
                    close(u)

                    ! Poynting
                    u = 99
                    open(unit=u, file="temp/cmplx/poynting_s_x.ppm", status='unknown')
                    rv = write_ppm_p3(u, poynting_s(:,:)%x)
                    close(u)

                    open(unit=u, file="temp/cmplx/poynting_s_y.ppm", status='unknown')
                    rv = write_ppm_p3(u, poynting_s(:,:)%y)
                    close(u)

                    open(unit=u, file="temp/cmplx/poynting_s_z.ppm", status='unknown')
                    rv = write_ppm_p3(u, poynting_s(:,:)%z)
                    close(u)

                    u = 99
                    open(unit=u, file="temp/cmplx/poynting_s_x_log.ppm", status='unknown')
                    rv = write_ppm_p3(u, poynting_s(:,:)%x, logarithmic = .true.)
                    close(u)

                    open(unit=u, file="temp/cmplx/poynting_s_y_log.ppm", status='unknown')
                    rv = write_ppm_p3(u, poynting_s(:,:)%y, logarithmic = .true.)
                    close(u)

                    open(unit=u, file="temp/cmplx/poynting_s_z_log.ppm", status='unknown')
                    rv = write_ppm_p3(u, poynting_s(:,:)%z, logarithmic = .true.)
                    close(u)

                    ! Poynting abs
                    open(unit=u, file="temp/cmplx/poynting_s_abs.ppm", status='unknown')
                    rv = write_ppm_p3(u, poynting_s_abs)
                    close(u)
                    open(unit=u, file="temp/cmplx/poynting_s_abs_log.ppm", status='unknown')
                    rv = write_ppm_p3(u, poynting_s_abs, logarithmic = .true.)
                    close(u)

                end function test_get_field_scattered_plane_section_cmplx

        end function lib_mie_scattering_by_a_sphere_test_functions

end module lib_mie_scattering_by_a_sphere

