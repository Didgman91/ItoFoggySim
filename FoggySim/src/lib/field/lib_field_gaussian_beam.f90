module lib_field_gaussian_beam
    use libmath
    use lib_constants
    use lib_field_polarisation
    implicit none

    private

    ! public function
    public :: lib_field_gaussian_beam_hermite_get_field

    public :: lib_field_gaussian_beam_test_functions

    ! public types
    public :: lib_field_gaussian_beam_hermite_type

    type lib_field_gaussian_beam_hermite_type
        double precision :: e_field_0                   ! [V/m]
        double precision :: wave_length_0               ! [m]
        double precision :: refractive_index_medium     ! [1]
        double precision :: waist_x0                    ! [m]
        double precision :: waist_y0                    ! [m]
        integer :: tem_m                                ! Hermite-Gaussian Mode
        integer :: tem_n                                ! Hermite-Gaussian Mode
        type(jones_vector_type) :: polarisation
    end type

    contains

        ! Argument
        ! ----
        !   parameter: type(lib_field_gaussian_beam_hermite_type)
        !       parameter of the Hermite-Gaussian beam
        !   evaluation_point_x: type(cartesian_coordinate_real_type)
        !       evaluation point at the beam koordinate system
        !   polarisation: type(jones_vector_type)
        !       beam polarisation
        !   psi: double precsision, optional (std: 0)
        !       rotation about the x-axis [rad]
        !   theta: double precsision, optional (std: 0)
        !       rotation about the y-axis [rad]
        !   phi: double precision, optional (std: 0)
        !       rotation about the z-axis [rad]
        !
        ! Returns
        ! ----
        !   e_field: type(cartesian_coordinate_cmplx_type)
        !       electical field component at "evaluation_point_x"
        !   h_field: type(cartesian_coordinate_cmplx_type)
        !       magnetical field component at "evaluation_point_x"
        !
        ! Convention
        ! ----
        !   Tait-Bryan angles: z-y'-x'' (intrinsic rotations) or
        !                      x-y-z   (extrinsic rotations)
        !   The intrinsic rotations are known as: yaw, pitch and roll
        !
        subroutine lib_field_gaussian_beam_hermite_get_field(parameter, evaluation_point_x, &
                                                   e_field, h_field,&
                                                   phi, theta, psi)
            implicit none

            ! dummy
            type(lib_field_gaussian_beam_hermite_type), intent(in) :: parameter
            type(cartesian_coordinate_real_type), intent(in) :: evaluation_point_x

            double precision, intent(in), optional :: phi
            double precision, intent(in), optional :: theta
            double precision, intent(in), optional :: psi


            type(cartesian_coordinate_cmplx_type) :: e_field
            type(cartesian_coordinate_cmplx_type) :: h_field

            ! auxiliary
            double complex :: field
            double precision :: wave_impedance

            double precision :: m_phi
            double precision :: m_theta
            double precision :: m_psi
            type(cartresian_coordinate_rot_matrix_type) :: rot
            type(cartesian_coordinate_real_type) :: point_x

            m_phi = 0
            if (present(phi)) m_phi = phi

            m_theta = 0
            if (present(theta)) m_theta = theta

            m_psi = 0
            if (present(psi)) m_psi = psi

            if (m_phi .ne. 0 .or. m_theta .ne. 0 .or. m_psi .ne. 0) then
                rot  = lib_math_get_matrix_rot_x_y_z(m_phi, m_theta, m_psi)
                point_x = rot * evaluation_point_x
            else
                point_x = evaluation_point_x
            end if

            field = lib_fiel_gaussian_beam_hermite_scalar(parameter%e_field_0, &
                                                          parameter%wave_length_0, &
                                                          parameter%refractive_index_medium, &
                                                          parameter%waist_x0, parameter%waist_y0, &
                                                          point_x, &
                                                          tem_m=parameter%tem_m, &
                                                          tem_n=parameter%tem_n)
            e_field = field * parameter%polarisation;

            ! Quantum Optics: An Introduction, Mark Fox
            ! eq. 2.25 but plane wave
            wave_impedance = const_z_0 / parameter%refractive_index_medium

            rot  = lib_math_get_matrix_rot_x_y_z(0d0, 0d0, PI / 2d0)
            h_field = rot * e_field / wave_impedance

            if (m_phi .ne. 0 .or. m_theta .ne. 0 .or. m_psi .ne. 0) then
                rot  = lib_math_get_matrix_rot_x_y_z(m_phi, m_theta, m_psi)
                e_field = rot * e_field
                h_field = rot * h_field
            end if

        end subroutine

        ! Argument
        ! ----
        !   e_field_0: double precision
        !       electical field [V/m]
        !   wave_length: double precision
        !       vacuum wave length [m]
        !   n_medium: double precision
        !       refractive index of the medium
        !   waist_x0: double precision
        !       beam waist at z = 0 in along the x axis
        !   waist_y0: double precision
        !       beam waist at z = 0 in along the y axis
        !   x: type(cartesian_coordinate_real)
        !       evaluation point [m^3]
        !   tem_m, tem_n: integer, optional(std: 0)
        !       order of the Hermite-Gaussian mode
        !
        ! Returns
        ! ----
        !   field: double cmplx
        !       Gaussian field distribution with a propagation along the z-axis
        !
        !
        !          \     z       / <-- Gaussian beam
        !           \     ^     /
        !            \ K_B|    /
        !          -->|   --> |<-- 2 w_0x,y
        !            /     x,y \
        !           /           \
        !
        !       K_B: beam coordinate system
        !       z:   propagation along the z-axis
        !
        ! TeX: $$ \begin{aligned} E_{m n}(x, y, z, t)=& \sqrt{\frac{w_{x, 0}}{w_{x}(z)}} H_{m}\left(\frac{x \sqrt{2}}{w_{x}(z)}\right) e^{-\frac{x^{2}}{w_{x}^{2}(z)}} e^{-i \frac{k}{2 R_{x}(z)}} \\ &\times \sqrt{\frac{w_{y, 0}}{w_{y}(z)} H_{n}\left(\frac{y \sqrt{2}}{w_{y}(z)}\right.}{w_{y}(z)}) e^{\frac{y^{2}}{w_{y}^{2}(z)}}{e^{w_{y}^{2}(z)}}-i \frac{k}{2 R_{y}(z)} \\ & \times E_{0} \sqrt{\frac{1}{2^{m+n} m ! n ! \pi}} e^{-i(k z-\omega t)} \times \\ & \times e^{i(m+1 / 2) \arctan \left(z / z_{R, x}\right)} e^{i(n+1 / 2) \arctan \left(z / z_{R, y}\right)} \end{aligned}$$
        !
        ! Reference: Laser in der Fertigung Helmut Hügel, Thomas Graf
        !            eq. 2.28 without the time-dependent term
        function lib_fiel_gaussian_beam_hermite_scalar(e_field_0, wave_length_0, n_medium, &
                                                  waist_x0, waist_y0, &
                                                  x, tem_m, tem_n) result(field)
            implicit none
            ! dummy
            double precision, intent(in) :: e_field_0
            double precision, intent(in) :: wave_length_0
            double precision, intent(in) :: n_medium
            double precision, intent(in) :: waist_x0
            double precision, intent(in) :: waist_y0
            type(cartesian_coordinate_real_type), intent(in) :: x
            integer, intent(in), optional :: tem_m
            integer, intent(in), optional :: tem_n

            double complex :: field

            ! auxiliaray
            integer :: m_tem_m
            integer :: m_tem_n

            double precision :: w_x
            double precision :: zr_x
            double precision :: cr_x

            double precision :: w_y
            double precision :: zr_y
            double precision :: cr_y

            double precision :: k

            integer :: h_min_i
            integer :: h_max_i
            double precision, dimension(:,:), allocatable :: h

            double precision :: buffer_real

            double complex :: buffer_cmplx

            m_tem_m = 0
            if (present(tem_m)) m_tem_m = tem_m

            m_tem_n = 0
            if (present(tem_n)) m_tem_n = tem_n


            ! pre calc
            k = 2d0 * PI * n_medium / wave_length_0

            zr_x = lib_field_gaussian_beam_rayleigh_length(waist_x0, wave_length_0, n_medium)
            w_x = lib_field_gaussian_beam_waist(waist_x0, zr_x, x%z)
            cr_x = lib_field_gaussian_beam_curvature(zr_x, x%z)

            zr_y = lib_field_gaussian_beam_rayleigh_length(waist_y0, wave_length_0, n_medium)
            w_y = lib_field_gaussian_beam_waist(waist_y0, zr_y, x%z)
            cr_y = lib_field_gaussian_beam_curvature(zr_y, x%z)


            h_min_i = min(m_tem_m, m_tem_n)
            h_max_i = max(m_tem_m, m_tem_n)

            call lib_math_hermite_polynomial((/ x%x * sqrt(2d0) / w_x, &
                                               x%y * sqrt(2d0) / w_y /), &
                                            h_min_i, h_max_i - h_min_i + 1, h)

            buffer_real = sqrt(waist_x0 / w_x) * h(1, m_tem_m) * exp(- x%x**2 / w_x**2)
            buffer_cmplx = buffer_real * exp(dcmplx(0, -k * x%x**2 * cr_x / 2))

            buffer_real = sqrt(waist_y0 / w_y) * h(2, m_tem_n) * exp(- x%y**2 / w_y**2)
            buffer_cmplx = buffer_real * exp(dcmplx(0, -k * x%y**2 * cr_y / 2)) * buffer_cmplx

            buffer_real = e_field_0 * sqrt(2d0**dble(-m_tem_m - m_tem_n) &
                                           / (dble(lib_math_factorial_get_factorial(m_tem_m) &
                                              * lib_math_factorial_get_factorial(m_tem_n)) * PI))
            buffer_cmplx = buffer_real * exp(dcmplx(0, -k * x%z)) * buffer_cmplx

            buffer_real = (dble(m_tem_m) + 0.5d0) * atan2(x%z, zr_x) &
                          + (dble(m_tem_n) + 0.5d0) * atan2(x%z, zr_y)
            buffer_cmplx = exp(dcmplx(0,buffer_real)) * buffer_cmplx

            field = buffer_cmplx

        end function lib_fiel_gaussian_beam_hermite_scalar

        ! Argument
        ! ----
        !   w0: double precision
        !       beam waist at z = 0
        !   zr: double precision
        !       Rayleigh length
        !   z: double precision
        !       distance along the z axis
        !
        ! Retruns
        ! ----
        !   w: double precision
        !       beam waist at z
        !
        ! Reference: Laser in der Fertigung Helmut Hügel, Thomas Graf
        !            eq. 2.29
        function lib_field_gaussian_beam_waist(w0, zr, z) result(w)
            implicit none
            ! dummy
            double precision, intent(in) :: w0
            double precision, intent(in) :: zr
            double precision, intent(in) :: z

            double precision :: w

            w = w0 * sqrt(1 + (z/zr)**2)

        end function lib_field_gaussian_beam_waist

        ! Argument
        ! ----
        !   zr: double precision
        !       Rayleigh length
        !   z: double precision
        !       distance along the z axis
        !
        ! Retruns
        ! ----
        !   c: double precision
        !       radius of curvature of the beam's wavefront
        !
        ! Formula
        ! ----
        !   c = R^-1
        !
        ! Reference: Laser in der Fertigung Helmut Hügel, Thomas Graf
        !            inverse of eq. 2.35
        function lib_field_gaussian_beam_curvature(zr, z) result(c)
            implicit none
            ! dummy
            double precision, intent(in) :: zr
            double precision, intent(in) :: z

            double precision :: c

            ! R[zr_, z_] := zr(z/zr + zr/z);
            ! cr[zr_, z_] := R[zr, z]^-1;
            c = z / (z**2 + zr**2)

        end function lib_field_gaussian_beam_curvature

        ! Argument
        ! ----
        !   w0: double precision
        !       beam waist at z = 0
        !   z: double precision
        !       distance along the z axis
        !
        ! Retruns
        ! ----
        !   zr: double precision
        !       Rayleigh length
        !
        ! Formula
        ! ----
        !   c = R^-1
        !
        ! Reference: Laser in der Fertigung Helmut Hügel, Thomas Graf
        !            eq. 2.30
        function lib_field_gaussian_beam_rayleigh_length(w0, wave_length_0, n_medium) result(zr)
            implicit none
            ! dummy
            double precision, intent(in) :: w0
            double precision, intent(in) :: wave_length_0
            double precision, intent(in) :: n_medium

            double precision :: zr

            zr = PI * w0**2 * n_medium / wave_length_0

        end function lib_field_gaussian_beam_rayleigh_length

        function lib_field_gaussian_beam_test_functions() result(rv)
            use lib_field
            implicit none
            ! dummy
            integer :: rv

            rv = 0

            if (.not. test_lib_field_gaussian_beam_hermite_get_field()) rv = rv + 1
!            if (.not. test_lib_field_gaussian_beam_hermite_get_field_2()) rv = rv + 1

            contains

            function test_lib_field_gaussian_beam_hermite_get_field() result(rv)
                implicit none
                ! dummy
                logical :: rv

                ! auxiliary
                integer :: i
                integer :: ii

                double precision :: x
                double precision :: y
                double precision :: z
                type(lib_field_gaussian_beam_hermite_type) :: gauss_parameter

                double precision, dimension(2) :: x_range
                double precision, dimension(2) :: y_range
                real(kind=8) :: step_size

                integer :: no_x_values
                integer :: no_y_values

                type(cartesian_coordinate_real_type) :: point_cartesian

                type(cartesian_coordinate_cmplx_type) :: buffer_e_field
                type(cartesian_coordinate_cmplx_type) :: buffer_h_field

                type(cartesian_coordinate_cmplx_type), dimension(:, :), allocatable :: e_field
                type(cartesian_coordinate_cmplx_type), dimension(:, :), allocatable :: h_field

                x_range = (/ -10_8 * unit_mu, 10.0_8 * unit_mu /)
                y_range = (/ -10_8 * unit_mu, 10.0_8 * unit_mu /)
!                    step_size = 0.02_8 * unit_mu
                step_size = 0.075_8 * unit_mu


                no_x_values = abs(int(floor((x_range(2)-x_range(1))/step_size)))
                no_y_values = abs(int(floor((y_range(2)-y_range(1))/step_size)))

                allocate(e_field(no_x_values, no_y_values))
                allocate(h_field(no_x_values, no_y_values))

                x = 0
                y = 0
                z = 10 * unit_mu

                gauss_parameter%e_field_0 = 1
                gauss_parameter%refractive_index_medium = 1
                gauss_parameter%tem_m = 0
                gauss_parameter%tem_n = 0
!                gauss_parameter%polarisation = lib_field_polarisation_jones_vector_type_get_linear_rot(PI/4d0)
                gauss_parameter%polarisation = lib_field_polarisation_jones_vector_type_get_linear_h()
!                gauss_parameter%polarisation = lib_field_polarisation_jones_vector_type_get_circular_plus()
                gauss_parameter%waist_x0 = 2.5 * unit_mu
                gauss_parameter%waist_y0 = 2.5 * unit_mu
                gauss_parameter%wave_length_0 = 0.55 * unit_mu

!                !$OMP PARALLEL DO PRIVATE(i, ii) FIRSTPRIVATE x, y, z)
                do i=1, no_x_values
                    x = x_range(1) + (i-1) * step_size
                    do ii= no_y_values, 1, -1
                        y = y_range(1) + (ii-1) * step_size

                        point_cartesian%x = x
                        point_cartesian%y = y
                        point_cartesian%z = z

                        call lib_field_gaussian_beam_hermite_get_field(gauss_parameter, &
                                                                       point_cartesian, &
                                                                       buffer_e_field, buffer_h_field, &
                                                                       phi = PI / 8d0) !, &
!                                                                       psi = PI / 8d0)
                        e_field(i,ii) = buffer_e_field
                        h_field(i,ii) = buffer_h_field
                    end do
                end do
!                !$OMP END PARALLEL DO

                rv = lib_field_export(e_field, h_field, "temp/real/gauss_")

            end function test_lib_field_gaussian_beam_hermite_get_field

            function test_lib_field_gaussian_beam_hermite_get_field_2() result(rv)
                implicit none
                ! dummy
                logical :: rv

                ! auxiliary
                integer :: i
                integer :: ii

                double precision :: x
                double precision :: y
                double precision :: z
                type(lib_field_gaussian_beam_hermite_type), dimension(2) :: gauss_parameter

                double precision, dimension(2) :: x_range
                double precision, dimension(2) :: y_range
                real(kind=8) :: step_size

                integer :: no_x_values
                integer :: no_y_values

                type(cartesian_coordinate_real_type) :: point_cartesian

                type(cartesian_coordinate_cmplx_type) :: buffer_e_field
                type(cartesian_coordinate_cmplx_type) :: buffer_h_field

                type(cartesian_coordinate_cmplx_type), dimension(:, :), allocatable :: e_field
                type(cartesian_coordinate_cmplx_type), dimension(:, :), allocatable :: h_field

                x_range = (/ -10_8 * unit_mu, 10.0_8 * unit_mu /)
                y_range = (/ -10_8 * unit_mu, 10.0_8 * unit_mu /)
!                    step_size = 0.02_8 * unit_mu
                step_size = 0.075_8 * unit_mu


                no_x_values = abs(int(floor((x_range(2)-x_range(1))/step_size)))
                no_y_values = abs(int(floor((y_range(2)-y_range(1))/step_size)))

                allocate(e_field(no_x_values, no_y_values))
                allocate(h_field(no_x_values, no_y_values))

                x = 0
                y = 0
                z = 10 * unit_mu

                gauss_parameter(:)%e_field_0 = 1
                gauss_parameter(:)%refractive_index_medium = 1
                gauss_parameter(:)%tem_m = 0
                gauss_parameter(:)%tem_n = 0
!                gauss_paramete(:)r%polarisation = lib_field_polarisation_jones_vector_type_get_linear_rot(PI/4d0)
                gauss_parameter(:)%polarisation = lib_field_polarisation_jones_vector_type_get_linear_h()
!                gauss_paramete(:)r%polarisation = lib_field_polarisation_jones_vector_type_get_circular_plus()
                gauss_parameter(:)%waist_x0 = 2.5 * unit_mu
                gauss_parameter(:)%waist_y0 = 2.5 * unit_mu
                gauss_parameter(:)%wave_length_0 = 0.55 * unit_mu

!                !$OMP PARALLEL DO PRIVATE(i, ii) FIRSTPRIVATE x, y, z)
                do i=1, no_x_values
                    x = x_range(1) + (i-1) * step_size
                    do ii=1, no_y_values
                        y = y_range(1) + (ii-1) * step_size

                        point_cartesian%x = x
                        point_cartesian%y = y
                        point_cartesian%z = z

                        call lib_field_gaussian_beam_hermite_get_field(gauss_parameter(1), &
                                                                       point_cartesian, &
                                                                       buffer_e_field, buffer_h_field , &
                                                                       theta = PI / 16d0)!, &
!                                                                       psi = PI / 16d0)
                        e_field(i,ii) = buffer_e_field
                        h_field(i,ii) = buffer_h_field

                        call lib_field_gaussian_beam_hermite_get_field(gauss_parameter(2), &
                                                                       point_cartesian, &
                                                                       buffer_e_field, buffer_h_field , &
                                                                       theta = PI / 16d0)!, &
!                                                                       psi = PI / 16d0)
                        e_field(i,ii) = e_field(i,ii) + buffer_e_field
                        h_field(i,ii) = h_field(i,ii) + buffer_h_field
                    end do
                end do
!                !$OMP END PARALLEL DO

                rv = lib_field_export(e_field, h_field, "temp/real/gauss_")

            end function

        end function

end module lib_field_gaussian_beam
