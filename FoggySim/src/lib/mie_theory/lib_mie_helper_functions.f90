module lib_mie_helper_functions
    use libmath
    use lib_mie_type
    use lib_mie_vector_spherical_harmonics
    implicit none

    private

    public :: lib_mie_hf_init_coeff_p0_q0
    public :: lib_mie_hf_destructor
    public :: lib_mie_hf_get_An
    public :: lib_mie_hf_get_n_c
    public :: lib_mie_hf_get_p_q_j_j
    public :: lib_mie_hf_calc_triple_sum

    interface lib_mie_hf_get_An
        module procedure get_An_real
        module procedure get_An_cmplx
    end interface lib_mie_hf_get_An

    type cache_coefficients_p_0_q_0_type
        double precision :: alpha
        double precision :: beta
        integer :: n_max
        type(list_list_cmplx), allocatable :: p_0
        type(list_list_cmplx), allocatable :: q_0
    end type

    type(cache_coefficients_p_0_q_0_type), dimension(:), allocatable :: cache_coefficients_p_0_q_0
    logical :: cache_coefficients_p_0_q_0_enabled = .false.

    contains

        subroutine lib_mie_hf_init_coeff_p0_q0(alpha, beta, n_max)
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

                call get_p_q_j_j_core(alpha(i), beta(i), (/1, n_max(i)/), &
                                      cache_coefficients_p_0_q_0(i)%p_0, &
                                      cache_coefficients_p_0_q_0(i)%q_0)
            end do

            cache_coefficients_p_0_q_0_enabled = .true.

        end subroutine lib_mie_hf_init_coeff_p0_q0

        subroutine lib_mie_hf_destructor
            implicit none
            ! dummy
            if (allocated(cache_coefficients_p_0_q_0)) then
                deallocate(cache_coefficients_p_0_q_0)
            end if
            cache_coefficients_p_0_q_0_enabled = .false.
        end subroutine

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
            m_n_mx = max(lib_mie_hf_get_n_c(abs(x)), int(abs(m_mx))) + 15

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
            m_n_mx = max(lib_mie_hf_get_n_c(abs(x)), int(abs(m_mx))) + 15

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
        function lib_mie_hf_get_n_c(x) result (rv)
            implicit none
            ! dummy
            double precision :: x

            integer :: rv

            ! auxiliary
            double precision :: dummy

            dummy = x + 4.05_8 * x**(1.0_8/3.0_8) + 2.0_8

            rv = int(ceiling(dummy))

        end function lib_mie_hf_get_n_c

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
        subroutine lib_mie_hf_get_p_q_j_j(k, d_0_j, n_range, p, q)
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
                call get_p_q_j_j_core(alpha, beta, n_range, p, q, exp_k_d)
            else
                ! 0-th coordinate system
                call get_p_q_j_j_core(alpha, beta, n_range, p, q)
            end if

        end subroutine lib_mie_hf_get_p_q_j_j

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
        subroutine get_p_q_j_j_core(alpha, beta, n_range, p, q, exp_k_d, caching)
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

        end subroutine get_p_q_j_j_core

        ! Arguments
        ! ----
        !   simulation_parameter: type(simulation_parameter_type)
        !   sphere: type(sphere_type), dimension(:)
        !       list of spheres
        !   sphere_parameter: type(sphere_parameter_type), dimension(:)
        !       list of shared sphere parameters
        !   z_selector: integer
        !       parameter of the spherical harmonics
        !       values:
        !           1: spherical Bessel function first kind   j_n
        !           2: spherical Bessel function second kind  y_n
        !           3: spherical Hankel function first kind   h^(1)_n
        !           4: spherical Hankel function second kind  h^(2)_n
        !   f: double precision, optional (std: 1)
        !       numerical factor (0, 1]
        !       "In our actual calculations, some multi-
        !        sphere systems do not converge if f 5 1, but they do
        !        converge when the value of f is reduced to, say, 0.7."[1]
        !
        ! Returns
        ! ----
        !   a_j_nm: type(list_list_cmplx)
        !   b_j_nm: type(list_list_cmplx)
        !
        ! Reference: Electromagnetic scattering by an aggregate of spheres, Yu-lin Xu, eq. 30, 35
        subroutine lib_mie_hf_calc_triple_sum(simulation_parameter, &
                                              sphere, sphere_parameter, sphere_j, &
                                              z_selector, &
                                              a_j_nm, b_j_nm, &
                                              f)
            implicit none
            ! dummy
            type(simulation_parameter_type) :: simulation_parameter
            type(sphere_type), dimension(:), intent(in) :: sphere
            type(sphere_parameter_type), dimension(:), intent(in) :: sphere_parameter
            integer :: sphere_j
            integer(kind=1) :: z_selector

            type(list_list_cmplx), intent(inout) :: a_j_nm
            type(list_list_cmplx), intent(inout) :: b_j_nm

            double precision, intent(in), optional :: f

            ! auxiliary
            integer :: n
            integer :: nu

            integer :: l
            integer :: i

            double precision :: lambda
            double precision :: n_medium
            type(cartesian_coordinate_real_type) :: k

            integer, dimension(2) :: n_range
            integer, dimension(2) :: nu_range

            type(list_list_cmplx) :: p
            type(list_list_cmplx) :: q
            type(list_cmplx) :: a_j_n
            type(list_cmplx) :: b_j_n
            type(list_cmplx) :: a_l_nu
            type(list_cmplx) :: b_l_nu
            type(list_list_cmplx) :: a_l_numu
            type(list_list_cmplx) :: b_l_numu
            type(list_4_cmplx) :: A_nmnumu
            type(list_4_cmplx) :: B_nmnumu
            type(cartesian_coordinate_real_type) :: d_0_j
            type(cartesian_coordinate_real_type) :: d_0_l
            type(cartesian_coordinate_real_type) :: x

            type(list_list_cmplx) :: buffer_sum_a
            type(list_list_cmplx) :: buffer_sum_b
            type(list_list_cmplx), dimension(:), allocatable :: buffer_sum_a_l
            type(list_list_cmplx), dimension(:), allocatable :: buffer_sum_b_l

            type(list_list_cmplx) :: a_j_nm_old
            type(list_list_cmplx) :: b_j_nm_old

            i = sphere(sphere_j)%sphere_parameter_index
            n_range = sphere_parameter(i)%n_range
            d_0_j = sphere(sphere_j)%d_0_j
            a_j_n = sphere_parameter(i)%a_n
            b_j_n = sphere_parameter(i)%b_n

            if (present(f)) then
                a_j_nm_old = sphere(sphere_j)%a_nm
                b_j_nm_old = sphere(sphere_j)%b_nm
            end if

            lambda = simulation_parameter%lambda
            n_medium = simulation_parameter%refractive_index_medium
            k = simulation_parameter%wave_vector


            allocate(buffer_sum_a_l(lbound(sphere, 1):ubound(sphere, 1)))
            allocate(buffer_sum_b_l(lbound(sphere, 1):ubound(sphere, 1)))

            ! calculate the triple sum
            !$OMP PARALLEL DO PRIVATE(l, i, nu_range, a_l_nu, b_l_nu, d_0_l, &
            !$OMP&  a_nmnumu, b_nmnumu, x, p, q, nu)
            do l=ubound(sphere, 1), lbound(sphere, 1)
                if (l .ne. sphere_j) then
                    i = sphere(l)%sphere_parameter_index
                    nu_range = sphere_parameter(i)%n_range
                    a_l_numu = sphere(l)%a_nm
                    b_l_numu = sphere(l)%b_nm
                    d_0_l = sphere(l)%d_0_j

                    x = (d_0_l - d_0_j) * n_medium / lambda

                    call lib_mie_vector_spherical_harmonics_translation_coefficient(x, &
                                                                                    n_range, nu_range, z_selector,&
                                                                                    a_nmnumu, b_nmnumu)


                    call calc_inner_2sum(a_l_numu, b_l_numu, a_nmnumu, b_nmnumu, &
                                         buffer_sum_a_l(l), buffer_sum_b_l(l))

                end if
            end do
            !$OMP END PARALLEL DO

            call init_list(buffer_sum_a, n_range(1), n_range(2) - n_range(1) + 1, cmplx(0, 0, kind=lib_math_type_kind))
            call init_list(buffer_sum_b, n_range(1), n_range(2) - n_range(1) + 1, cmplx(0, 0, kind=lib_math_type_kind))
            ! final summation of the triple sum
            do l=lbound(sphere, 1), ubound(sphere, 1)
                buffer_sum_a = buffer_sum_a + buffer_sum_a_l(l)
                buffer_sum_b = buffer_sum_b + buffer_sum_b_l(l)
            end do

            deallocate(buffer_sum_a_l)
            deallocate(buffer_sum_b_l)


            call lib_mie_hf_get_p_q_j_j(k, d_0_j, n_range, p, q)

            call init_list(a_j_nm, n_range(1), n_range(2) - n_range(1) + 1)
            call init_list(b_j_nm, n_range(1), n_range(2) - n_range(1) + 1)

            a_j_nm = p - buffer_sum_a
            b_j_nm = q - buffer_sum_b

            !$OMP PARALLEL DO PRIVATE(n)
            do n=n_range(1), n_range(2)
                a_j_nm%item(n)%item(:) = a_j_n%item(n) * a_j_nm%item(n)%item(:)
                b_j_nm%item(n)%item(:) = b_j_n%item(n) * b_j_nm%item(n)%item(:)
            end do
            !$OMP END PARALLEL DO

            if (present(f)) then
                a_j_nm = (1.0d0 - f) * a_j_nm_old + f * a_j_nm
                b_j_nm = (1.0d0 - f) * b_j_nm_old + f * b_j_nm
            end if

        end subroutine

        ! Calculates the two inner sums of eq. 30 (interactive scattering coefficients)
        !
        ! Argument
        ! ----
        !   a_l_munu: type(list_list_cmplx)
        !       Mie coefficient of the l-th sphere
        !       restriction: same dimension as b_l_munu
        !   b_l_munu: type(list_list_cmplx)
        !       Mie coefficient of the l-th sphere
        !       restriction: same dimension as a_l_munu
        !   a_munumn: type(list_4_cmplx)
        !       translation coefficient
        !       restriction: same dimension as b_munumn
        !   b_munumn: type(list_4_cmplx)
        !       translation coefficient
        !       restriction: same dimension as a_munumn
        !
        ! Returns
        ! ----
        !   a_2sum: double complex
        !
        !
        ! Reference:Electromagnetic scattering by an aggregate of spheres, Yu-lin Xu, eq. 30
        subroutine calc_inner_2sum(a_l_munu, b_l_munu, a_munumn, b_munumn, &
                                   a_2sum, b_2sum)
            implicit none
            ! dummy
            type(list_list_cmplx), intent(in) :: a_l_munu
            type(list_list_cmplx), intent(in) :: b_l_munu
            type(list_4_cmplx), intent(in) :: a_munumn
            type(list_4_cmplx), intent(in) :: b_munumn

            type(list_list_cmplx), intent(inout) :: a_2sum
            type(list_list_cmplx), intent(inout) :: b_2sum

            ! auxiliary
            integer :: n
            integer :: m
            integer :: nu
            integer :: mu

            integer, dimension(2) :: n_range
            integer, dimension(2) :: nu_range

            double complex :: buffer_cmplx

            ! set boundaries
            n_range(1) = lbound(a_munumn%item, 1)
            n_range(2) = ubound(a_munumn%item, 1)
            nu_range(1) = max(lbound(a_munumn%item(n_range(1))%item(0)%item, 1), 1)
            nu_range(2) = ubound(a_munumn%item(n_range(1))%item(0)%item, 1)

            call init_list(a_2sum, n_range(1), n_range(2)-n_range(1)+1)
            call init_list(b_2sum, n_range(1), n_range(2)-n_range(1)+1)

            ! calculate the summand
            !$OMP PARALLEL DO PRIVATE(n, m)
            do n=n_range(1), n_range(2)
                do m=-n, n
                    !$OMP PARALLEL DO PRIVATE(nu, mu, buffer_cmplx)
                    do nu=nu_range(1), nu_range(2)
                        do mu=-nu, nu
                            ! summand of a_j_mn
                            ! first line eq. 30
                            buffer_cmplx = a_l_munu%item(nu)%item(mu) * a_munumn%item(n)%item(m)%item(nu)%item(mu)
                            buffer_cmplx = buffer_cmplx &
                                           + b_l_munu%item(nu)%item(mu) * b_munumn%item(n)%item(m)%item(nu)%item(mu)
                            a_2sum%item(n)%item(m) = buffer_cmplx

                            ! summand of b_j_mn
                            ! second line eq. 30
                            buffer_cmplx = a_l_munu%item(nu)%item(mu) * b_munumn%item(n)%item(m)%item(nu)%item(mu)
                            buffer_cmplx = buffer_cmplx &
                                           + b_l_munu%item(nu)%item(mu) * a_munumn%item(n)%item(m)%item(nu)%item(mu)
                            b_2sum%item(n)%item(m) = buffer_cmplx
                        end do
                    end do
                    !$OMP END PARALLEL DO
                end do
            end do
            !$OMP END PARALLEL DO

        end subroutine calc_inner_2sum

end module lib_mie_helper_functions
