! Scattering by a sphere
! ----
!
! - incident plane wave along the z-axis
!
#define _DEBUG_

!#define _ONLY_M_ 1

module lib_mie_single_sphere
    !$  use omp_lib
    use libmath
    use lib_constants
    use lib_mie_vector_spherical_harmonics
    use lib_mie_type
    use lib_mie_ss_helper_functions
    implicit none

    private

    ! --- public ---
    public :: lib_mie_ss_destructor

    public:: lib_mie_ss_test_functions

    contains

        subroutine lib_mie_ss_destructor
            implicit none

            call lib_mie_ss_hf_destructor

        end subroutine lib_mie_ss_destructor

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
        ! Reference: Electromagnetic scattering by an aggregate of spheres, Yu-lin Xu, eq. 17
        !
        function get_field_initial_incident_xu_real(theta, phi, r, &
                                                    e_field_0, lambda, &
                                                    n_medium, &
                                                    n_range, &
                                                    alpha, beta) &
                                                result (field_0)
            implicit none
            ! dummy
            double precision, intent(in) :: theta
            double precision, intent(in) :: phi
            double precision, intent(in) :: r
            double precision, intent(in) :: e_field_0
            double precision, intent(in) :: lambda
            double precision, intent(in) :: n_medium
            integer(kind=4), dimension(2),intent(in) :: n_range
            double precision, intent(in), optional :: alpha
            double precision, intent(in), optional :: beta

            type(spherical_coordinate_cmplx_type), dimension(2) :: field_0

            ! auxiliary
            integer(kind=4) :: i
            type(spherical_coordinate_cmplx_type) :: e_field_incident_0
            type(spherical_coordinate_cmplx_type) :: h_field_incident_0

            double precision :: k0 ! wave number = 2 pi / lambda

            double precision :: buffer_real
            complex(kind=8) :: buffer_cmplx

            integer(kind=4) :: m
            integer(kind=4) :: n

            integer(kind=1) :: z_selector

            type(list_spherical_coordinate_cmplx_type), dimension(:), allocatable :: M_nm
            type(list_spherical_coordinate_cmplx_type), dimension(:), allocatable :: N_nm

            type(list_list_cmplx) :: e_field_nm
            type(spherical_coordinate_cmplx_type), dimension(n_range(2)-n_range(1)+1) :: e_field_n_incident_0
            type(spherical_coordinate_cmplx_type) :: buffer_e_field_n_incident_0
            type(spherical_coordinate_cmplx_type), dimension(n_range(2)-n_range(1)+1) :: h_field_n_incident_0
            type(spherical_coordinate_cmplx_type) :: buffer_h_field_n_incident_0

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

            cos_beta = cos(m_beta)
            sin_beta = sin(m_beta)

!            mu = 1
!            mu1 = 1

            call init_list(e_field_nm, n_range(1), n_range(2)-n_range(1)+1, dcmplx(0,0))
            call init_list(p0_nm, n_range(1), n_range(2)-n_range(1)+1, dcmplx(0,0))
            call init_list(q0_nm, n_range(1), n_range(2)-n_range(1)+1, dcmplx(0,0))
            call init_list(calc_order_m, n_range(1), n_range(2)-n_range(1)+1, .true.)

            ! --- pre-calc ---
            ! errata eq. (1) Equations (21) on p. 4577
            call lib_math_associated_legendre_polynomial_theta(m_alpha, n_range(2), pi_nm, tau_nm)
            !$OMP PARALLEL DO PRIVATE(n, m, buffer_real, buffer_cmplx)
            do n=n_range(1), n_range(2)
#ifdef _ONLY_M_
                m=_ONLY_M_
#else
                do m=-n, n
#endif
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
#ifndef _ONLY_M_
                end do
#endif
            end do
            !$OMP END PARALLEL DO


            z_selector = 1

            call lib_mie_vector_spherical_harmonics_components(theta, phi, r, k0 * n_medium, n_range, z_selector, &
                                                               M_nm, N_nm)
            ! eq. (5)
            !$OMP PARALLEL DO PRIVATE(n, m, buffer_real)
            do n=n_range(1), n_range(2)
#ifdef _ONLY_M_
                m=_ONLY_M_
#else
                do m=-n, n
#endif
                    if (calc_order_m%item(n)%item(m)) then
                        buffer_real = abs(e_field_0) * real((2*n+1), kind=8) &
                                      * lib_math_factorial_get_n_minus_m_divided_by_n_plus_m(n,m)

#ifdef _DEBUG_
                        if (isnan(buffer_real)) then
                            print *, "get_field_initial_incident_xu_real: ERROR"
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
#ifndef _ONLY_M_
                end do
#endif
            end do
            !$OMP END PARALLEL DO

            ! eq. (17)
            !$OMP PARALLEL DO PRIVATE(n, m, i, buffer_cmplx, buffer_e_field_n_incident_0, buffer_h_field_n_incident_0)
            do n= n_range(1), n_range(2)
                i = n - n_range(1) + 1

                e_field_n_incident_0(i)%rho = cmplx(0,0)
                e_field_n_incident_0(i)%theta = cmplx(0,0)
                e_field_n_incident_0(i)%phi = cmplx(0,0)
                buffer_e_field_n_incident_0%rho = cmplx(0,0)
                buffer_e_field_n_incident_0%theta = cmplx(0,0)
                buffer_e_field_n_incident_0%phi = cmplx(0,0)

                h_field_n_incident_0(i)%rho = cmplx(0,0)
                h_field_n_incident_0(i)%theta = cmplx(0,0)
                h_field_n_incident_0(i)%phi = cmplx(0,0)
                buffer_h_field_n_incident_0%rho = cmplx(0,0)
                buffer_h_field_n_incident_0%theta = cmplx(0,0)
                buffer_h_field_n_incident_0%phi = cmplx(0,0)

#ifdef _ONLY_M_
                m=_ONLY_M_
#else
                do m=-n, n
#endif
                    if (calc_order_m%item(n)%item(m)) then
                        buffer_cmplx = cmplx(0, 1, kind=8) * e_field_nm%item(n)%item(m)
                        buffer_e_field_n_incident_0 = buffer_cmplx &
                                                      * (N_nm(n)%coordinate(m) * p0_nm%item(n)%item(m) &
                                                         +M_nm(n)%coordinate(m) * q0_nm%item(n)%item(m))

                        buffer_h_field_n_incident_0 = e_field_nm%item(n)%item(m) &
                                             * (M_nm(n)%coordinate(m) * p0_nm%item(n)%item(m) &
                                                +N_nm(n)%coordinate(m) * q0_nm%item(n)%item(m))

!                        buffer_cmplx = cmplx(0, 1, kind=8)
!                        buffer_e_field_n_incident_0 = buffer_cmplx &
!                                                      * (N_nm(n)%coordinate(m) &
!                                                         +M_nm(n)%coordinate(m) )
!
!                        buffer_h_field_n_incident_0 = e_field_nm%item(n)%item(m) &
!                                             * (N_nm(n)%coordinate(m) &
!                                                +M_nm(n)%coordinate(m) )
#ifdef _DEBUG_
                        if (isnan(real(buffer_e_field_n_incident_0%rho)) &
                            .or. isnan(aimag(buffer_e_field_n_incident_0%rho)) &
                            .or. isnan(real(buffer_e_field_n_incident_0%phi)) &
                            .or. isnan(aimag(buffer_e_field_n_incident_0%phi)) &
                            .or. isnan(real(buffer_e_field_n_incident_0%theta)) &
                            .or. isnan(aimag(buffer_e_field_n_incident_0%theta)) ) then
                            print *, "get_field_initial_incident_xu_real: ERROR"
                            print *, "  buffer_e_field_n_incident_0 is NaN"
                            print * , "  n = ", n
                            print * , "  m = ", m
                        end if
#endif
                        e_field_n_incident_0(i) = e_field_n_incident_0(i) + buffer_e_field_n_incident_0
                        h_field_n_incident_0(i) = h_field_n_incident_0(i) + buffer_h_field_n_incident_0
                    end if
#ifndef _ONLY_M_
                end do
#endif
            end do
            !$OMP END PARALLEL DO

            e_field_incident_0%theta = cmplx(0,0,kind=8)
            e_field_incident_0%phi = cmplx(0,0,kind=8)
            e_field_incident_0%rho = cmplx(0,0,kind=8)

            h_field_incident_0%theta = cmplx(0,0,kind=8)
            h_field_incident_0%phi = cmplx(0,0,kind=8)
            h_field_incident_0%rho = cmplx(0,0,kind=8)

            do i=n_range(2)-n_range(1)+1, 1, -1
                e_field_incident_0 = e_field_incident_0 + e_field_n_incident_0(i)
                h_field_incident_0 = h_field_incident_0 + h_field_n_incident_0(i)
            end do

            field_0(1) = (-1D0) * e_field_incident_0

            ! omega = k * v
            ! v = c0 / n_medium
            ! mu = 1 <-- definition
            !
            !   k / (omega * mu) = k / ( k * v * 1 )
            ! = k / ( k * c0 / n_medium )
            ! = n_medium / c0
            field_0(2) = (-1D0) * h_field_incident_0 * n_medium / real(const_c0, kind=8)

        end function get_field_initial_incident_xu_real

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
        ! Reference: Electromagnetic scattering by an aggregate of spheres, Yu-lin Xu, eq. 17
        !
        function get_field_initial_incident_multi_wave_real(theta, phi, r, &
                                                    e_field_0, &
                                                    n_medium, &
                                                    n_range, &
                                                    illumination) &
                                                result (field_0)
            implicit none
            ! dummy
            double precision, intent(in) :: theta
            double precision, intent(in) :: phi
            double precision, intent(in) :: r
            double precision, intent(in) :: e_field_0
            double precision, intent(in) :: n_medium
            integer(kind=4), dimension(2), intent(in) :: n_range
            type(lib_mie_illumination_parameter), dimension(:), intent(in) :: illumination

            type(spherical_coordinate_cmplx_type), dimension(2) :: field_0

            ! auxiliary
            integer(kind=4) :: i
            type(spherical_coordinate_cmplx_type) :: e_field_incident_0
            type(spherical_coordinate_cmplx_type) :: h_field_incident_0

            double precision :: lambda
            double precision :: k0 ! wave number = 2 pi / lambda

            double precision :: buffer_real
            complex(kind=8) :: buffer_cmplx

            integer(kind=4) :: m
            integer(kind=4) :: n

            integer(kind=1) :: z_selector

            type(list_spherical_coordinate_cmplx_type), dimension(:), allocatable :: M_nm
            type(list_spherical_coordinate_cmplx_type), dimension(:), allocatable :: N_nm

            type(list_list_cmplx) :: e_field_nm
            type(spherical_coordinate_cmplx_type), dimension(n_range(2)-n_range(1)+1) :: e_field_n_incident_0
            type(spherical_coordinate_cmplx_type) :: buffer_e_field_n_incident_0
            type(spherical_coordinate_cmplx_type), dimension(n_range(2)-n_range(1)+1) :: h_field_n_incident_0
            type(spherical_coordinate_cmplx_type) :: buffer_h_field_n_incident_0

            type(list_list_cmplx) :: p0_nm
            type(list_list_cmplx) :: q0_nm

            type(cartesian_coordinate_real_type) :: d_0_j

            type(list_list_logical) :: calc_order_m

            ! --- init ---
            ! standard values

            lambda = illumination(lbound(illumination, 1))%lambda_0 / n_medium

            k0 = 2.0_8 * PI / lambda

!            mu = 1
!            mu1 = 1

            call init_list(e_field_nm, n_range(1), n_range(2)-n_range(1)+1, dcmplx(0,0))
            call init_list(p0_nm, n_range(1), n_range(2)-n_range(1)+1, dcmplx(0,0))
            call init_list(q0_nm, n_range(1), n_range(2)-n_range(1)+1, dcmplx(0,0))
            call init_list(calc_order_m, n_range(1), n_range(2)-n_range(1)+1, .true.)

            ! --- pre-calc ---
            d_0_j%x = 0
            d_0_j%y = 0
            d_0_j%z = 0

            call lib_mie_ss_hf_get_p_q_j_j(illumination, n_medium, d_0_j, n_range, p0_nm, q0_nm)

            z_selector = 1

            call lib_mie_vector_spherical_harmonics_components(theta, phi, r, k0 * n_medium, n_range, z_selector, &
                                                               M_nm, N_nm)
            ! eq. (5)
            !$OMP PARALLEL DO PRIVATE(n, m, buffer_real)
            do n=n_range(1), n_range(2)
#ifdef _ONLY_M_
                m=_ONLY_M_
#else
                do m=-n, n
#endif
                    if (calc_order_m%item(n)%item(m)) then
                        buffer_real = abs(e_field_0) * real((2*n+1), kind=8) &
                                      * lib_math_factorial_get_n_minus_m_divided_by_n_plus_m(n,m)

#ifdef _DEBUG_
                        if (isnan(buffer_real)) then
                            print *, "get_field_initial_incident_xu_real: ERROR"
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
#ifndef _ONLY_M_
                end do
#endif
            end do
            !$OMP END PARALLEL DO

            ! eq. (17)
            !$OMP PARALLEL DO PRIVATE(n, m, i, buffer_cmplx, buffer_e_field_n_incident_0, buffer_h_field_n_incident_0)
            do n= n_range(1), n_range(2)
                i = n - n_range(1) + 1

                e_field_n_incident_0(i)%rho = cmplx(0,0)
                e_field_n_incident_0(i)%theta = cmplx(0,0)
                e_field_n_incident_0(i)%phi = cmplx(0,0)
                buffer_e_field_n_incident_0%rho = cmplx(0,0)
                buffer_e_field_n_incident_0%theta = cmplx(0,0)
                buffer_e_field_n_incident_0%phi = cmplx(0,0)

                h_field_n_incident_0(i)%rho = cmplx(0,0)
                h_field_n_incident_0(i)%theta = cmplx(0,0)
                h_field_n_incident_0(i)%phi = cmplx(0,0)
                buffer_h_field_n_incident_0%rho = cmplx(0,0)
                buffer_h_field_n_incident_0%theta = cmplx(0,0)
                buffer_h_field_n_incident_0%phi = cmplx(0,0)

#ifdef _ONLY_M_
                m=_ONLY_M_
#else
                do m=-n, n
#endif
                    if (calc_order_m%item(n)%item(m)) then
                        buffer_cmplx = cmplx(0, 1, kind=8) * e_field_nm%item(n)%item(m)
                        buffer_e_field_n_incident_0 = buffer_cmplx &
                                                      * (N_nm(n)%coordinate(m) * p0_nm%item(n)%item(m) &
                                                         +M_nm(n)%coordinate(m) * q0_nm%item(n)%item(m))

                        buffer_h_field_n_incident_0 = e_field_nm%item(n)%item(m) &
                                             * (M_nm(n)%coordinate(m) * p0_nm%item(n)%item(m) &
                                                +N_nm(n)%coordinate(m) * q0_nm%item(n)%item(m))

!                        buffer_cmplx = cmplx(0, 1, kind=8)
!                        buffer_e_field_n_incident_0 = buffer_cmplx &
!                                                      * (N_nm(n)%coordinate(m) &
!                                                         +M_nm(n)%coordinate(m) )
!
!                        buffer_h_field_n_incident_0 = e_field_nm%item(n)%item(m) &
!                                             * (N_nm(n)%coordinate(m) &
!                                                +M_nm(n)%coordinate(m) )
#ifdef _DEBUG_
                        if (isnan(real(buffer_e_field_n_incident_0%rho)) &
                            .or. isnan(aimag(buffer_e_field_n_incident_0%rho)) &
                            .or. isnan(real(buffer_e_field_n_incident_0%phi)) &
                            .or. isnan(aimag(buffer_e_field_n_incident_0%phi)) &
                            .or. isnan(real(buffer_e_field_n_incident_0%theta)) &
                            .or. isnan(aimag(buffer_e_field_n_incident_0%theta)) ) then
                            print *, "get_field_initial_incident_xu_real: ERROR"
                            print *, "  buffer_e_field_n_incident_0 is NaN"
                            print * , "  n = ", n
                            print * , "  m = ", m
                        end if
#endif
                        e_field_n_incident_0(i) = e_field_n_incident_0(i) + buffer_e_field_n_incident_0
                        h_field_n_incident_0(i) = h_field_n_incident_0(i) + buffer_h_field_n_incident_0
                    end if
#ifndef _ONLY_M_
                end do
#endif
            end do
            !$OMP END PARALLEL DO

            e_field_incident_0%theta = cmplx(0,0,kind=8)
            e_field_incident_0%phi = cmplx(0,0,kind=8)
            e_field_incident_0%rho = cmplx(0,0,kind=8)

            h_field_incident_0%theta = cmplx(0,0,kind=8)
            h_field_incident_0%phi = cmplx(0,0,kind=8)
            h_field_incident_0%rho = cmplx(0,0,kind=8)

            do i=n_range(2)-n_range(1)+1, 1, -1
                e_field_incident_0 = e_field_incident_0 + e_field_n_incident_0(i)
                h_field_incident_0 = h_field_incident_0 + h_field_n_incident_0(i)
            end do

            field_0(1) = (-1D0) * e_field_incident_0

            ! omega = k * v
            ! v = c0 / n_medium
            ! mu = 1 <-- definition
            !
            !   k / (omega * mu) = k / ( k * v * 1 )
            ! = k / ( k * c0 / n_medium )
            ! = n_medium / c0
            field_0(2) = (-1D0) * h_field_incident_0 * n_medium / real(const_c0, kind=8)

        end function get_field_initial_incident_multi_wave_real

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
            call init_list(calc_order_m, n_range(1), n_range(2)-n_range(1)+1, .true.)

            ! --- pre-calc ---
            call lib_mie_ss_hf_get_coefficients_a_n_b_n(size_parameter, n_particle/n_medium, n_range, a_n, b_n)
!            call get_coefficients_a_b_real(size_parameter, n_particle/n_medium, mu, mu1, n_range, a_n, b_n)

            ! errata eq. (1) Equations (21) on p. 4577
            call lib_math_associated_legendre_polynomial_theta(m_alpha, n_range(2), pi_nm, tau_nm)
            !$OMP PARALLEL DO PRIVATE(n, m, buffer_real, buffer_cmplx)
            do n=n_range(1), n_range(2)
#ifdef _ONLY_M_
                m=_ONLY_M_
#else
                do m=-n, n
#endif
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
#ifndef _ONLY_M_
                end do
#endif
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
#ifdef _ONLY_M_
                    m=_ONLY_M_
#else
                    do m=-n, n
#endif
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
#ifndef _ONLY_M_
                    end do
#endif
                end do
                !$OMP END PARALLEL DO

                ! first line eq. (4)
                !$OMP PARALLEL DO PRIVATE(n, m, i, buffer_cmplx, buffer_e_field_n_s, buffer_h_field_n_s)
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

#ifdef _ONLY_M_
                    m=_ONLY_M_
#else
                    do m=-n, n
#endif
                        if (calc_order_m%item(n)%item(m)) then
                            buffer_cmplx = cmplx(0, 1, kind=8) * e_field_nm%item(n)%item(m)
                            buffer_e_field_n_s = buffer_cmplx * (a_n(i)*N_nm(n)%coordinate(m) * p0_nm%item(n)%item(m) &
                                                                 +b_n(i)*M_nm(n)%coordinate(m) * q0_nm%item(n)%item(m))

                            buffer_h_field_n_s = e_field_nm%item(n)%item(m) &
                                                 * (b_n(i)*N_nm(n)%coordinate(m) * q0_nm%item(n)%item(m) &
                                                    +a_n(i)*M_nm(n)%coordinate(m) * p0_nm%item(n)%item(m))
#ifdef _DEBUG_
                            if (isnan(real(buffer_e_field_n_s%rho)) .or. isnan(aimag(buffer_e_field_n_s%rho)) &
                                .or. isnan(real(buffer_e_field_n_s%phi)) .or. isnan(aimag(buffer_e_field_n_s%phi)) &
                                .or. isnan(real(buffer_e_field_n_s%theta)) .or. isnan(aimag(buffer_e_field_n_s%theta)) ) then
                                print *, "get_field_scattered_xu_real: ERROR"
                                print *, "  e_field_n_s is NaN"
                                print * , "  n = ", n
                                print * , "  m = ", m
                            end if
#endif
                            e_field_n_s(i) = e_field_n_s(i) + buffer_e_field_n_s
                            h_field_n_s(i) = h_field_n_s(i) + buffer_h_field_n_s
                        end if
#ifndef _ONLY_M_
                    end do
#endif
                end do
                !$OMP END PARALLEL DO

                e_field_s%theta = cmplx(0,0,kind=8)
                e_field_s%phi = cmplx(0,0,kind=8)
                e_field_s%rho = cmplx(0,0,kind=8)

                h_field_s%theta = cmplx(0,0,kind=8)
                h_field_s%phi = cmplx(0,0,kind=8)
                h_field_s%rho = cmplx(0,0,kind=8)

!                do i=n_range(2)-n_range(1)+1, 1, -1
                do i=1, n_range(2)-n_range(1)+1
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
            call init_list(calc_order_m, n_range(1), n_range(2)-n_range(1)+1, .true.)

            ! --- pre-calc ---
            call lib_mie_ss_hf_get_coefficients_a_n_b_n(size_parameter, n_particle/n_medium, n_range, a_n, b_n)

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
                !$OMP PARALLEL DO PRIVATE(n, m, i, buffer_cmplx, buffer_e_field_n_s, buffer_h_field_n_s)
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
                                                 * (b_n(i)*N_nm(n)%coordinate(m) * q0_nm%item(n)%item(m) &
                                                    +a_n(i)*M_nm(n)%coordinate(m) * p0_nm%item(n)%item(m))

#ifdef _DEBUG_
                            if (isnan(real(e_field_n_s(i)%rho)) .or. isnan(aimag(e_field_n_s(i)%rho)) &
                                .or. isnan(real(e_field_n_s(i)%phi)) .or. isnan(aimag(e_field_n_s(i)%phi)) &
                                .or. isnan(real(e_field_n_s(i)%theta)) .or. isnan(aimag(e_field_n_s(i)%theta)) ) then
                                print *, "get_field_scattered_xu_cmplx: ERROR"
                                print *, "  e_field_n_s is NaN"
                                print * , "  n = ", n
                                print * , "  m = ", m
                            end if

                            if (isinf(real(e_field_n_s(i)%rho)) .or. isinf(aimag(e_field_n_s(i)%rho)) &
                                .or. isinf(real(e_field_n_s(i)%phi)) .or. isinf(aimag(e_field_n_s(i)%phi)) &
                                .or. isinf(real(e_field_n_s(i)%theta)) .or. isinf(aimag(e_field_n_s(i)%theta)) ) then
                                print *, "get_field_scattered_xu_cmplx: ERROR"
                                print *, "  e_field_n_s is infinity"
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

        ! Argument
        ! ----
        !   lambda: double precsision
        !       vacuum wave length
        !   n_medium: double precision
        !       refractive index of the medium
        !   r_particle: double precision
        !       radius of the sphere
        !   n_particle: double complex
        !       refractive index of the medium
        !       aimag(n_particle) > 0: absorbing medium
        !   n_range: integer, dimension(2)
        !       first and last index (degree) of the sum to calculate the electical field
        !       e.g. n = (/ 1, 45 /) <-- size parameter
        !   sphere_parameter
        !
        !
        !
        !

        function lib_mie_ss_test_functions() result (rv)
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

!            if (.not. test_get_e_field_scattered_real()) then
!                rv = rv + 1
!            end if
!            if (.not. test_get_field_initial_incident_xu_real()) then
!                rv = rv + 1
!            end if
            if (.not. test_get_field_initial_incident_multi_wave_real()) then
                rv = rv + 1
            end if
            if (.not. test_get_field_scattered_plane_section_real()) then
                rv = rv + 1
            end if
!            if (.not. test_get_field_scattered_plane_section_cmplx()) then
!                rv = rv + 1
!            end if

            call cpu_time(test_finish)
            call system_clock(test_count_finish, test_count_rate)

            print *, ""
            print *, "------lib_mie_ss_scattering_by_a_sphere_test_functions------"
            print '("  CPU-Time = ",f10.3," seconds.")',test_finish-test_start
            print '("  WALL-Time = ",f10.3," seconds.")',(test_count_finish-test_count_start) / real(test_count_rate)
            print *, ""
            if (rv == 0) then
                print *, "lib_mie_ss_test_functions tests: OK"
            else
                print *, rv,"lib_mie_ss_test_functions test(s) FAILED"
            end if
            print *, "------------------------------------------------------------"
            print *, ""

            contains

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
                    n_range(2) = lib_mie_ss_hf_get_n_c(r_particle * k0 * n_particle)

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

                function test_get_field_initial_incident_xu_real() result (rv)
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

                    double precision :: alpha
                    double precision :: beta

                    double precision :: lambda
                    double precision :: k0
                    double precision :: e_field_0
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

                    x_range = (/ -5_8 * unit_mu, 5.0_8 * unit_mu /)
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

                    e_field_0 = 1
                    lambda = 0.7 * unit_mu
                    alpha = PI / 4D0
                    beta = 0

                    n_medium = 1

                    k0 = 2 * PI / lambda

                    n_range(1) = 1
                    n_range(2) = N_MAX
                    print *, "NOTE: max degree = ", n_range(2)

                    call lib_math_factorial_initialise_caching(n_range(2))

                    call system_clock(count_start, count_rate)
                    call cpu_time(start)
                    !$OMP PARALLEL DO PRIVATE(i, ii, point_cartesian, point_spherical, buffer_field) &
                    !$OMP  PRIVATE(buffer_cartesian_cmplx) &
                    !$OMP  FIRSTPRIVATE(x, y, z)
                    do i=1, no_x_values
                        x = x_range(1) + (i-1) * step_size
                        do ii=1, no_z_values
                            z = z_range(1) + (ii-1) * step_size

                            point_cartesian%x = x
                            point_cartesian%y = y
                            point_cartesian%z = z

                            point_spherical = point_cartesian

                            buffer_field = get_field_initial_incident_xu_real(point_spherical%theta, point_spherical%phi, &
                                                                       point_spherical%rho, &
                                                                       e_field_0, lambda, n_medium, &
                                                                       n_range, alpha=alpha, beta=beta)
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
                    print *, "get_field_initial_incident_xu_real"
                    print '("  CPU-Time = ",f10.3," seconds.")',finish-start

                    print '("  WALL-Time = ",f10.3," seconds.")',(count_finish-count_start) / real(count_rate)

                    ! --- wirte to PPM ---
                    ! e field
                    u = 99
                    open(unit=u, file="temp/real/e_field_0_x.ppm", status='unknown')
                    rv = write_ppm_p3(u, e_field_s_real_x)
                    close(u)

                    open(unit=u, file="temp/real/e_field_0_y.ppm", status='unknown')
                    rv = write_ppm_p3(u, e_field_s_real_y)
                    close(u)

                    open(unit=u, file="temp/real/e_field_0_z.ppm", status='unknown')
                    rv = write_ppm_p3(u, e_field_s_real_z)
                    close(u)

                    u = 99
                    open(unit=u, file="temp/real/e_field_0_x_log.ppm", status='unknown')
                    rv = write_ppm_p3(u, e_field_s_real_x, logarithmic = .true.)
                    close(u)

                    open(unit=u, file="temp/real/e_field_0_y_log.ppm", status='unknown')
                    rv = write_ppm_p3(u, e_field_s_real_y, logarithmic = .true.)
                    close(u)

                    open(unit=u, file="temp/real/e_field_0_z_log.ppm", status='unknown')
                    rv = write_ppm_p3(u, e_field_s_real_z, logarithmic = .true.)
                    close(u)

                    ! h field
                    u = 99
                    open(unit=u, file="temp/real/h_field_0_x.ppm", status='unknown')
                    rv = write_ppm_p3(u, h_field_s_real_x)
                    close(u)

                    open(unit=u, file="temp/real/h_field_0_y.ppm", status='unknown')
                    rv = write_ppm_p3(u, h_field_s_real_y)
                    close(u)

                    open(unit=u, file="temp/real/h_field_0_z.ppm", status='unknown')
                    rv = write_ppm_p3(u, h_field_s_real_z)
                    close(u)

                    u = 99
                    open(unit=u, file="temp/real/h_field_0_x_log.ppm", status='unknown')
                    rv = write_ppm_p3(u, h_field_s_real_x, logarithmic = .true.)
                    close(u)

                    open(unit=u, file="temp/real/h_field_0_y_log.ppm", status='unknown')
                    rv = write_ppm_p3(u, h_field_s_real_y, logarithmic = .true.)
                    close(u)

                    open(unit=u, file="temp/real/h_field_0_z_log.ppm", status='unknown')
                    rv = write_ppm_p3(u, h_field_s_real_z, logarithmic = .true.)
                    close(u)

                    ! Poynting
                    u = 99
                    open(unit=u, file="temp/real/poynting_0_x.ppm", status='unknown')
                    rv = write_ppm_p3(u, poynting_s(:,:)%x)
                    close(u)

                    open(unit=u, file="temp/real/poynting_0_y.ppm", status='unknown')
                    rv = write_ppm_p3(u, poynting_s(:,:)%y)
                    close(u)

                    open(unit=u, file="temp/real/poynting_0_z.ppm", status='unknown')
                    rv = write_ppm_p3(u, poynting_s(:,:)%z)
                    close(u)

                    u = 99
                    open(unit=u, file="temp/real/poynting_0_x_log.ppm", status='unknown')
                    rv = write_ppm_p3(u, poynting_s(:,:)%x, logarithmic = .true.)
                    close(u)

                    open(unit=u, file="temp/real/poynting_0_y_log.ppm", status='unknown')
                    rv = write_ppm_p3(u, poynting_s(:,:)%y, logarithmic = .true.)
                    close(u)

                    open(unit=u, file="temp/real/poynting_0_z_log.ppm", status='unknown')
                    rv = write_ppm_p3(u, poynting_s(:,:)%z, logarithmic = .true.)
                    close(u)

                    ! Poynting abs
                    open(unit=u, file="temp/real/poynting_0_abs.ppm", status='unknown')
                    rv = write_ppm_p3(u, poynting_s_abs)
                    close(u)
                    open(unit=u, file="temp/real/poynting_0_abs_log.ppm", status='unknown')
                    rv = write_ppm_p3(u, poynting_s_abs, logarithmic = .true.)
                    close(u)

                end function test_get_field_initial_incident_xu_real

                function test_get_field_initial_incident_multi_wave_real() result (rv)
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

                    double precision :: alpha
                    double precision :: beta

                    double precision :: lambda
                    double precision :: k0
                    double precision :: e_field_0
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

                    type(cartesian_coordinate_real_type) :: buffer_car
                    type(lib_mie_illumination_parameter), dimension(2) :: illumination

                    ! CPU-time
                    real :: start, finish
                    ! WALL-time
                    INTEGER :: count_start, count_finish, count_rate

                    x_range = (/ -5_8 * unit_mu, 5.0_8 * unit_mu /)
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

                    e_field_0 = 1
                    lambda = 0.7 * unit_mu
                    alpha = PI / 4D0
                    beta = 0

                    n_medium = 1

                    k0 = 2 * PI / lambda

                    illumination(:)%lambda_0 = lambda

                    buffer_car%x = 0
                    buffer_car%y = 0
                    buffer_car%z = 0
                    illumination(:)%d_0_i = buffer_car

                    buffer_car%x = 1
                    buffer_car%y = 0
                    buffer_car%z = 1
                    buffer_car = buffer_car / abs(buffer_car) / lambda
                    illumination(1)%wave_vector_0 = buffer_car

                    buffer_car%x = -1
                    buffer_car%y = 0
                    buffer_car%z = 1
                    buffer_car = buffer_car / abs(buffer_car) / lambda
                    illumination(2)%wave_vector_0 = buffer_car

!                    buffer_car%x = 0
!                    buffer_car%y = 0
!                    buffer_car%z = 1
!                    buffer_car = buffer_car / abs(buffer_car) / lambda
!                    illumination(3)%wave_vector_0 = buffer_car

                    n_range(1) = 1
                    n_range(2) = 20!N_MAX
                    print *, "NOTE: max degree = ", n_range(2)

                    call lib_math_factorial_initialise_caching(n_range(2))

                    call system_clock(count_start, count_rate)
                    call cpu_time(start)
                    !$OMP PARALLEL DO PRIVATE(i, ii, point_cartesian, point_spherical, buffer_field) &
                    !$OMP  PRIVATE(buffer_cartesian_cmplx) &
                    !$OMP  FIRSTPRIVATE(x, y, z)
                    do i=1, no_x_values
                        x = x_range(1) + (i-1) * step_size
                        do ii=1, no_z_values
                            z = z_range(1) + (ii-1) * step_size

                            point_cartesian%x = x
                            point_cartesian%y = y
                            point_cartesian%z = z

                            point_spherical = point_cartesian

                            buffer_field = get_field_initial_incident_multi_wave_real(point_spherical%theta, &
                                                                       point_spherical%phi, &
                                                                       point_spherical%rho, &
                                                                       e_field_0, n_medium, &
                                                                       n_range, illumination)
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
                    print *, "get_field_initial_incident_xu_real"
                    print '("  CPU-Time = ",f10.3," seconds.")',finish-start

                    print '("  WALL-Time = ",f10.3," seconds.")',(count_finish-count_start) / real(count_rate)

                    ! --- wirte to PPM ---
                    ! e field
                    u = 99
                    open(unit=u, file="temp/real/e_field_0_x.ppm", status='unknown')
                    rv = write_ppm_p3(u, e_field_s_real_x)
                    close(u)

                    open(unit=u, file="temp/real/e_field_0_y.ppm", status='unknown')
                    rv = write_ppm_p3(u, e_field_s_real_y)
                    close(u)

                    open(unit=u, file="temp/real/e_field_0_z.ppm", status='unknown')
                    rv = write_ppm_p3(u, e_field_s_real_z)
                    close(u)

                    u = 99
                    open(unit=u, file="temp/real/e_field_0_x_log.ppm", status='unknown')
                    rv = write_ppm_p3(u, e_field_s_real_x, logarithmic = .true.)
                    close(u)

                    open(unit=u, file="temp/real/e_field_0_y_log.ppm", status='unknown')
                    rv = write_ppm_p3(u, e_field_s_real_y, logarithmic = .true.)
                    close(u)

                    open(unit=u, file="temp/real/e_field_0_z_log.ppm", status='unknown')
                    rv = write_ppm_p3(u, e_field_s_real_z, logarithmic = .true.)
                    close(u)

                    ! h field
                    u = 99
                    open(unit=u, file="temp/real/h_field_0_x.ppm", status='unknown')
                    rv = write_ppm_p3(u, h_field_s_real_x)
                    close(u)

                    open(unit=u, file="temp/real/h_field_0_y.ppm", status='unknown')
                    rv = write_ppm_p3(u, h_field_s_real_y)
                    close(u)

                    open(unit=u, file="temp/real/h_field_0_z.ppm", status='unknown')
                    rv = write_ppm_p3(u, h_field_s_real_z)
                    close(u)

                    u = 99
                    open(unit=u, file="temp/real/h_field_0_x_log.ppm", status='unknown')
                    rv = write_ppm_p3(u, h_field_s_real_x, logarithmic = .true.)
                    close(u)

                    open(unit=u, file="temp/real/h_field_0_y_log.ppm", status='unknown')
                    rv = write_ppm_p3(u, h_field_s_real_y, logarithmic = .true.)
                    close(u)

                    open(unit=u, file="temp/real/h_field_0_z_log.ppm", status='unknown')
                    rv = write_ppm_p3(u, h_field_s_real_z, logarithmic = .true.)
                    close(u)

                    ! Poynting
                    u = 99
                    open(unit=u, file="temp/real/poynting_0_x.ppm", status='unknown')
                    rv = write_ppm_p3(u, poynting_s(:,:)%x)
                    close(u)

                    open(unit=u, file="temp/real/poynting_0_y.ppm", status='unknown')
                    rv = write_ppm_p3(u, poynting_s(:,:)%y)
                    close(u)

                    open(unit=u, file="temp/real/poynting_0_z.ppm", status='unknown')
                    rv = write_ppm_p3(u, poynting_s(:,:)%z)
                    close(u)

                    u = 99
                    open(unit=u, file="temp/real/poynting_0_x_log.ppm", status='unknown')
                    rv = write_ppm_p3(u, poynting_s(:,:)%x, logarithmic = .true.)
                    close(u)

                    open(unit=u, file="temp/real/poynting_0_y_log.ppm", status='unknown')
                    rv = write_ppm_p3(u, poynting_s(:,:)%y, logarithmic = .true.)
                    close(u)

                    open(unit=u, file="temp/real/poynting_0_z_log.ppm", status='unknown')
                    rv = write_ppm_p3(u, poynting_s(:,:)%z, logarithmic = .true.)
                    close(u)

                    ! Poynting abs
                    open(unit=u, file="temp/real/poynting_0_abs.ppm", status='unknown')
                    rv = write_ppm_p3(u, poynting_s_abs)
                    close(u)
                    open(unit=u, file="temp/real/poynting_0_abs_log.ppm", status='unknown')
                    rv = write_ppm_p3(u, poynting_s_abs, logarithmic = .true.)
                    close(u)

                end function test_get_field_initial_incident_multi_wave_real

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
                    n_range(2) = lib_mie_ss_hf_get_n_c(r_particle * k0)! * n_particle)
                    if (n_range(2) .gt. N_MAX) then
                        print *, "WARNING: max degree (", N_MAX, ") reached: ", n_range(2)
                        n_range(2) = N_MAX
                    else
                        print *, "NOTE: max degree = ", n_range(2)
                    end if

                    call lib_math_factorial_initialise_caching(n_range(2))

                    call lib_mie_ss_hf_init_coeff_a_n_b_n((/ r_particle * k0 * n_medium /), &
                                                       (/ n_particle/n_medium /), &
                                                       (/ n_range(2) /))

                    call system_clock(count_start, count_rate)
                    call cpu_time(start)
                    !$OMP PARALLEL DO PRIVATE(i, ii, point_cartesian, point_spherical, buffer_field) &
                    !$OMP  PRIVATE(buffer_cartesian_cmplx) &
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

                    double precision :: alpha
                    double precision :: beta

                    ! CPU-time
                    real :: start, finish
                    ! WALL-time
                    INTEGER :: count_start, count_finish, count_rate

                    x_range = (/ -2.5_8 * unit_mu, 2.5_8 * unit_mu /)
                    z_range = (/ -2.5_8 * unit_mu, 2.5_8 * unit_mu /)
!                    step_size = 0.02_8 * unit_mu
                    step_size = 0.02_8 * unit_mu

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

                    r_particle = 0.35 * unit_mu

                    e_field_0 = 1
                    lambda = 0.7 * unit_mu
                    alpha = 0 !PI / 4D0
                    beta = 0

                    ! https://refractiveindex.info/?shelf=main&book=Ag&page=Johnson
                    n_particle = cmplx(0.040000, 7.1155, kind=8)
!                    n_particle = cmplx(1.5, 0, kind=8)
                    n_medium = 1

                    k0 = 2 * PI / lambda

                    n_range(1) = 1
!                    n_range(2) = min(45, lib_mie_ss_hf_get_n_c(r_particle * k0 * abs(n_particle))) ! todo: abs??
                    n_range(2) = lib_mie_ss_hf_get_n_c(r_particle * k0)
                    if (n_range(2) .gt. N_MAX) then
                        print *, "WARNING: max degree (", N_MAX, ") reached: ", n_range(2)
                        n_range(2) = N_MAX
                    else
                        print *, "NOTE: max degree = ", n_range(2)
                    end if

                    call system_clock(count_start, count_rate)
                    call cpu_time(start)
                    !$OMP PARALLEL DO PRIVATE(i, ii, point_cartesian, point_spherical, buffer) &
                    !$OMP  PRIVATE(buffer_cartesian_cmplx) &
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
                                                                  r_particle, n_particle, n_range, &
                                                                  alpha, beta)

                            e_field_s(i, ii) = make_cartesian(buffer(1), point_spherical%theta, point_spherical%phi)
                            h_field_s(i, ii) = make_cartesian(buffer(2), point_spherical%theta, point_spherical%phi)

                            e_field_s_real_x(i, ii) = real(e_field_s(i, ii)%x)
!                            e_field_s_real_x(i, ii) = abs(e_field_s(i, ii)%x)**2
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

        end function lib_mie_ss_test_functions

end module lib_mie_single_sphere

