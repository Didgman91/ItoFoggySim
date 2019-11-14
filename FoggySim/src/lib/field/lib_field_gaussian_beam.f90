module lib_field_gaussian_beam
    use libmath
    implicit none

    contains



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
        ! TeX: $$ \begin{aligned} E_{m n}(x, y, z, t)=& \sqrt{\frac{w_{x, 0}}{w_{x}(z)}} H_{m}\left(\frac{x \sqrt{2}}{w_{x}(z)}\right) e^{-\frac{x^{2}}{w_{x}^{2}(z)}} e^{-i \frac{k}{2 R_{x}(z)}} \\ &\times \sqrt{\frac{w_{y, 0}}{w_{y}(z)} H_{n}\left(\frac{y \sqrt{2}}{w_{y}(z)}\right.}{w_{y}(z)}) e^{\frac{y^{2}}{w_{y}^{2}(z)}}{e^{w_{y}^{2}(z)}}-i \frac{k}{2 R_{y}(z)} \\ & \times E_{0} \sqrt{\frac{1}{2^{m+n} m ! n ! \pi}} e^{-i(k z-\omega t)} \times \\ & \times e^{i(m+1 / 2) \arctan \left(z / z_{R, x}\right)} e^{i(n+1 / 2) \arctan \left(z / z_{R, y}\right)} \end{aligned}$$
        !
        ! Reference: Laser in der Fertigung Helmut Hügel, Thomas Graf
        !            eq. 2.28 without the time-dependent term
        subroutine lib_fiel_gaussian_beam_hermite_scalar(e_field_0, wave_length_0, n_medium, &
                                                  waist_x0, waist_y0, &
                                                  x, tem_m, tem_n, &
                                                  field)
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

            buffer_real = e_field_0 * sqrt(dble(2**(-m_tem_m - m_tem_n)) &
                                           / (dble(lib_math_factorial_get_factorial(m_tem_m) &
                                              * lib_math_factorial_get_factorial(m_tem_n)) * PI))
            buffer_cmplx = buffer_real * exp(dcmplx(0, -k * x%z)) * buffer_cmplx

            buffer_real = (dble(m_tem_m) + 0.5d0) * atan2(x%z, zr_x) &
                          + (dble(m_tem_n) + 0.5d0) * atan2(x%z, zr_y)
            buffer_cmplx = exp(dcmplx(0,buffer_real)) * buffer_cmplx

            field = buffer_cmplx

        end subroutine lib_fiel_gaussian_beam_hermite_scalar

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

end module lib_field_gaussian_beam
