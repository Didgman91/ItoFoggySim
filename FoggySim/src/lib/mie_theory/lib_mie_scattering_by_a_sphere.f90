! Scattering by a sphere
! ----
!
! - incident plane wave along the z-axis
!


module lib_mie_scattering_by_a_sphere
    use lib_math_public
    use lib_mie_vector_spherical_harmonics
    implicit none

    private

    ! --- public ---
    public:: lib_mie_scattering_by_a_sphere_test_functions

    contains

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
        !       polar angle
        !   phi: double precision
        !       azimuthal angle
        !   rho: double precision
        !       dimensionless varibale rho = k*a
        !       k: wavenumber
        !       a: distance
        !   rho_particle: double precision
        !       dimensionless varibale rho = k*r
        !       k: wavenumber
        !       r: radius of the sphere
        !   e_field_0: double precision
        !       amplitude of the incident wave
        !   n_particle: double precision
        !       refractive index of the sphere
        !   n_medium: double precision
        !       refractive index of the medium
        !
        function get_e_field_scattered(theta, phi, rho, e_field_0, rho_particle, n_particle, n_medium) result (e_field_s)
            implicit none
            ! dummy
            double precision, intent(in) :: theta
            double precision, intent(in) :: phi
            double precision, intent(in) :: rho
            double precision, intent(in) :: rho_particle
            double precision, intent(in) :: e_field_0
            double precision, intent(in) :: n_particle
            double precision, intent(in) :: n_medium

            type(spherical_coordinate_cmplx_type) :: e_field_s

            ! parameter
            integer(kind=4), dimension(2), parameter :: n = (/1, 10/)

            ! auxiliary
            double precision :: buffer_real
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
            type(spherical_coordinate_cmplx_type), dimension(n(2)-n(1)+1) :: e_field_n_s

            m = (/1,1/)

            mu = 1
            mu1 = 1

            call get_coefficients_a_b_real(rho_particle, n_particle/n_medium, mu, mu1, n, a_n, b_n)

            z_selector = 3

            call lib_mie_vector_spherical_harmonics_components(theta, phi, rho, m, n, z_selector, &
                                                               M_emn, M_omn, N_emn, N_omn, &
                                                               not_calc_Memn=.true., &
                                                               not_calc_Momn=.false., &
                                                               not_calc_Nemn=.false., &
                                                               not_calc_Nomn=.true.)
            ! page 93
            do i=n(1), n(2)
                buffer_real = e_field_0 * (2*i+1)/(i*(i+1))
                e_field_n(i) = cmplx(buffer_real,0, kind=8)
                e_field_n(i) = e_field_n(i) * (cmplx(0,1, kind=8)**i)
            end do

            ! eq. 4.45
            do i=n(1), n(2)
                e_field_n_s(i) = e_field_n(i) * (cmplx(0,1, kind=8) * a_n(i) * N_emn(m(1))%coordinate(i)&
                                                 - b_n(i) * M_omn(m(1))%coordinate(i))
            end do


            e_field_s%theta = cmplx(0,0,kind=8)
            e_field_s%phi = cmplx(0,0,kind=8)
            e_field_s%rho = cmplx(0,0,kind=8)
            do i=n(1), n(2)
                e_field_s = e_field_s + e_field_n_s(i)
            end do

            print*, "done"

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
        subroutine get_coefficients_a_b_cmplx(x, m, mu, mu1, n, a_n, b_n)
            implicit none
            ! dummy
            double precision :: x
            complex(kind=8) :: m
            double precision :: mu
            double precision :: mu1
            integer(kind=4), dimension(2) :: n

            complex(kind=8), dimension(n(2)-n(1)+1) :: a_n
            complex(kind=8), dimension(n(2)-n(1)+1) :: b_n

            ! auxiliary
            complex(kind=8), dimension(n(2)-n(1)+1) :: numerator
            complex(kind=8), dimension(n(2)-n(1)+1) :: denominator

            double precision, dimension(n(2)-n(1)+1) :: j_n_x
            complex(kind=8), dimension(n(2)-n(1)+1) :: j_n_mx
            double precision, dimension(n(2)-n(1)+1) :: s_n_x
            double precision, dimension(n(2)-n(1)+1) :: s_dn_x
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

            s_dn_mx = lib_math_riccati_s_derivative(mx, n(1), number_of_n, s_n_mx)
            j_n_mx = s_n_mx / mx

            xi_dn_x = lib_math_riccati_xi_derivative(x, n(1), number_of_n, xi_n_x)
            h_n_x = xi_n_x / x


            numerator = mu * m*m * j_n_mx * s_dn_x - mu1 * j_n_x * s_dn_mx
            denominator = mu * m*m * j_n_mx * xi_dn_x - mu1 * h_n_x * s_dn_mx

            a_n = numerator / denominator

            numerator = mu1 * j_n_mx * s_dn_x - mu * j_n_x * s_dn_mx
            denominator = mu1 * j_n_mx * xi_dn_x - mu * h_n_x * s_dn_mx

            b_n = numerator / denominator

        end subroutine get_coefficients_a_b_cmplx

        function lib_mie_scattering_by_a_sphere_test_functions() result (rv)
            implicit none
            ! dummy
            integer :: rv

            ! auxiliaray
            double precision, parameter :: ground_truth_e = 10.0_8**(-6.0_8)

            rv = 0

            if (.not. test_get_coefficients_a_b_cmplx()) then
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

                function test_get_coefficients_a_b_cmplx() result (rv)
                    implicit none
                    ! dummy
                    logical :: rv

                    ! auxiliary
                    integer(kind=4) :: i
                    real(kind=8) :: buffer

                    double precision :: x
                    complex(kind=8) :: m
                    double precision :: mu
                    double precision :: mu1
                    integer(kind=4), dimension(2), parameter :: n = (/1, 1/)

                    complex(kind=8), dimension(n(2)-n(1)+1) :: a_n
                    complex(kind=8), dimension(n(2)-n(1)+1) :: b_n

                    complex(kind=8), dimension(n(2)-n(1)+1) :: ground_truth_a_n
                    complex(kind=8), dimension(n(2)-n(1)+1) :: ground_truth_b_n

                    ! Reference: Electromagnetic scattering on spherical polydispersions, D.Deirmendjian, p. 27
                    m = cmplx(1.28, -1.37, kind=8)
                    x=20
                    ground_truth_a_n(1) = cmplx(-0.22686_8+0.5_8, -0.12863_8, kind=8)
                    ground_truth_b_n(1) = cmplx(0.22864_8+0.5_8, 0.13377_8, kind=8)

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

                end function

                function test_get_e_field_scattered() result (rv)
                    implicit none
                    ! dummy
                    logical :: rv

                    ! auxiliary
                    double precision :: theta
                    double precision :: phi
                    double precision :: rho
                    double precision :: e_field_0
                    double precision :: rho_particle
                    double precision :: n_particle
                    double precision :: n_medium

                    type(spherical_coordinate_cmplx_type) :: e_field_s

                    theta = PI / 4.0
                    phi = 0.0
                    rho = 100
                    rho_particle = 10

                    e_field_0 = 1
                    n_particle = 1.5
                    n_medium = 1

                    e_field_s = get_e_field_scattered(theta, phi, rho, e_field_0, rho_particle, n_particle, n_medium)

                end function

        end function lib_mie_scattering_by_a_sphere_test_functions

end module lib_mie_scattering_by_a_sphere
