module lib_mie_illumination
    use libmath
    implicit none

    private

    ! public functions
    public :: lib_mie_illumination_get_p_q_j_j
    public :: lib_mie_illumination_init_plane_wave
    public :: lib_mie_illumination_destructor

    ! public  types
    public :: lib_mie_illumination_parameter
    public :: lib_mie_illumination_plane_wave_parameter

    interface lib_mie_illumination_get_p_q_j_j
        module procedure lib_mie_illumination_get_p_q_j_j_single_plane_wave
        module procedure lib_mie_illumination_get_p_q_j_j_multi_plane_wave
    end interface

    type lib_mie_illumination_plane_wave_parameter
        double precision :: g ! ratio of the e-field of this plane wave to the "global" e-field:  e_x_field_0 / e_field_0
        type(cartesian_coordinate_real_type) :: d_0_i ! [m]
        type(cartesian_coordinate_real_type) :: wave_vector_0 ! |k| = 2 Pi / lambda, wave_vector = [1/m]
    end type lib_mie_illumination_plane_wave_parameter

    type lib_mie_illumination_elliptical_gaussian_beam_parameter
        double precision :: g ! ratio of the e-field of this Gaussian beam to the "global" e-field:  e_x_field_0 / e_field_0
        type(cartesian_coordinate_real_type) :: d_0_i ! [m]
        double precision :: wave_number_0 ! |k| = 2 Pi / lambda, wave_number = [1/m]
        double precision :: w0x ! beam waist radius along the x-axis
        double precision :: w0y ! beam waist radius along the y-axis
    end type lib_mie_illumination_elliptical_gaussian_beam_parameter

    ! Illumination Parameter
    ! ----
    !
    ! Plane Wave
    ! -----
    !             _________
    !             ___k^____
    !             ____|____
    !             _________ plane wave
    !                 z
    !                 ^
    !             K_i |
    !                 --> x
    !                ^
    !               /
    !           z  /d_0_i
    !           ^ /
    !       K_0 |/
    !           --> x
    !
    ! K_0: world coordinate system
    ! K_i: illumination coordinate system
    !
    !
    ! Gaussian Beam
    ! -----
    !          _______________
    !            ___________
    !             ___k,z___
    !               __^__   Gaussian beam
    !                _|_
    !                 --> x
    !                ^
    !               /
    !           z  /d_0_i
    !           ^ /
    !       K_0 |/
    !           --> x
    !
    ! K_0: world coordinate system
    ! K_i: illumination coordinate system
    !
    type lib_mie_illumination_parameter
        integer :: type    ! 1: plane wave
        double precision :: e_field_0 ! [V/m]
        double precision :: lambda_0 ! wave_length_vaccum
        type(lib_mie_illumination_plane_wave_parameter), dimension(:), allocatable :: plane_wave
        type(lib_mie_illumination_elliptical_gaussian_beam_parameter), dimension(:), allocatable :: elliptical_gaussian_beam
    end type lib_mie_illumination_parameter

    ! --- caching ---
    type cache_coefficients_p_0_q_0_plane_wave_type
        double precision :: alpha
        double precision :: beta
        integer :: n_max
        type(list_list_cmplx) :: p_0
        type(list_list_cmplx) :: q_0
    end type

    type(cache_coefficients_p_0_q_0_plane_wave_type), dimension(:), allocatable :: cache_coefficients_p_0_q_0_plane_wave
    logical :: cache_coefficients_p_0_q_0_plane_wave_enabled = .false.
    ! ~~~ caching ~~~

    contains

        ! Argument
        ! ----
        !   alpha: double precision, dimension(:)
        !       polar angle [0, Pi)
        !   beta: double precision, dimension(:)
        !       azimuthal angle [0, 2 Pi)
        !   n_max: integer, dimension(:)
        !       maximum degree of the polynomials
        subroutine lib_mie_illumination_init_plane_wave(alpha, beta, n_max)
            implicit none
            ! dummy
            double precision, dimension(:), intent(in) :: alpha
            double precision, dimension(lbound(alpha, 1):ubound(alpha, 1)), intent(in) :: beta
            integer(kind=4), dimension(lbound(alpha, 1):ubound(alpha, 1)) :: n_max

            ! auxiliary
            integer :: i

            ! --- init: cache_coefficients_a_b_cmplx_barberh_x ---
            if (allocated(cache_coefficients_p_0_q_0_plane_wave)) then
                deallocate(cache_coefficients_p_0_q_0_plane_wave)
            end if

            allocate(cache_coefficients_p_0_q_0_plane_wave(size(alpha)))
            cache_coefficients_p_0_q_0_plane_wave_enabled = .false.

            do i=1, size(alpha)
                cache_coefficients_p_0_q_0_plane_wave(i)%alpha = alpha(i)
                cache_coefficients_p_0_q_0_plane_wave(i)%beta = beta(i)
                cache_coefficients_p_0_q_0_plane_wave(i)%n_max = n_max(i)

                call get_p_q_j_j_plane_wave_core(alpha(i), beta(i), (/1, n_max(i)/), &
                                                 cache_coefficients_p_0_q_0_plane_wave(i)%p_0, &
                                                 cache_coefficients_p_0_q_0_plane_wave(i)%q_0, &
                                                 caching=.false.)
            end do

            cache_coefficients_p_0_q_0_plane_wave_enabled = .true.

        end subroutine lib_mie_illumination_init_plane_wave

        subroutine lib_mie_illumination_destructor
            implicit none

            if (allocated(cache_coefficients_p_0_q_0_plane_wave)) then
                deallocate(cache_coefficients_p_0_q_0_plane_wave)
            end if
            cache_coefficients_p_0_q_0_plane_wave_enabled = .false.
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
        !   k: type(cartesian_coordinate_real_type)
        !       wave vector
        !   d_0_j: type(cartesian_coordinate_real_type)
        !       vector from the
        !   n_range: integer, dimension(2)
        !       first and last index (degree) of the sum to calculate the electical field
        !       e.g. n = (/ 1, 45 /) <-- size parameter
        !   caching: logical (std: .true.)
        !
        ! Returns
        ! ----
        !   p: type(list_list_cmplx)
        !       coefficient of vector spherical componets
        !   q: type(list_list_cmplx)
        !       coefficient of vector spherical componets
        !
        ! Reference: Electromagnetic scattering by an aggregate of spheres Yu-lin Xu, eq. 20
        subroutine lib_mie_illumination_get_p_q_j_j_single_plane_wave(k, d_0_j, n_range, p, q, caching)
            implicit none
            ! dummy
            type(cartesian_coordinate_real_type), intent(in) :: k
            type(cartesian_coordinate_real_type), intent(in) :: d_0_j
            integer(kind=4), dimension(2),intent(in) :: n_range

            type(list_list_cmplx), intent(inout) :: p
            type(list_list_cmplx), intent(inout) :: q

            logical, optional :: caching

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
                call get_p_q_j_j_plane_wave_core(alpha, beta, n_range, p, q, exp_k_d, caching=caching)
            else
                ! 0-th coordinate system
                call get_p_q_j_j_plane_wave_core(alpha, beta, n_range, p, q, caching=caching)
            end if

        end subroutine lib_mie_illumination_get_p_q_j_j_single_plane_wave

        ! Calculates the coefficients (of vector spherical components) of multiple plane incident waves
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
        !   caching: logical (std: .true.)
        !
        ! Returns
        ! ----
        !   p: type(list_list_cmplx)
        !       coefficient of vector spherical componets
        !   q: type(list_list_cmplx)
        !       coefficient of vector spherical componets
        !
        ! Reference: Electromagnetic scattering by an aggregate of spheres Yu-lin Xu, eq. 20
        subroutine lib_mie_illumination_get_p_q_j_j_multi_plane_wave(illumination, n_medium, d_0_j, n_range, p_nm, q_nm, caching)
            implicit none
            ! dummy
            type(lib_mie_illumination_parameter), intent(in) :: illumination
            double precision :: n_medium
            type(cartesian_coordinate_real_type), intent(in) :: d_0_j
            integer(kind=4), dimension(2),intent(in) :: n_range

            type(list_list_cmplx), intent(inout) :: p_nm
            type(list_list_cmplx), intent(inout) :: q_nm

            logical, optional :: caching

            ! auxiliary
            integer :: i

            type(cartesian_coordinate_real_type) :: k
            type(cartesian_coordinate_real_type) :: d_i_j

            type(list_list_cmplx), dimension(size(illumination%plane_wave)) :: buffer_p_nm
            type(list_list_cmplx), dimension(size(illumination%plane_wave)) :: buffer_q_nm

            !$OMP PARALLEL DO PRIVATE(i, d_i_j, k)
            do i = 1, size(illumination%plane_wave)
                d_i_j = d_0_j - illumination%plane_wave(i)%d_0_i
                k = illumination%plane_wave(i)%wave_vector_0 * n_medium

                call lib_mie_illumination_get_p_q_j_j_single_plane_wave(k, d_i_j, n_range, &
                                                                        buffer_p_nm(i), buffer_q_nm(i),&
                                                                        caching=caching)
            end do
            !$OMP END PARALLEL DO

            call init_list(p_nm, n_range(1), n_range(2)-n_range(1)+1, dcmplx(0,0))
            call init_list(q_nm, n_range(1), n_range(2)-n_range(1)+1, dcmplx(0,0))

            do i = 1, size(illumination%plane_wave)
                p_nm = p_nm + illumination%plane_wave(i)%g * buffer_p_nm(i)
                q_nm = q_nm + illumination%plane_wave(i)%g * buffer_q_nm(i)
            end do

        end subroutine lib_mie_illumination_get_p_q_j_j_multi_plane_wave

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
        subroutine get_p_q_j_j_plane_wave_core(alpha, beta, n_range, p, q, exp_k_d, caching)
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

            if (cache_coefficients_p_0_q_0_plane_wave_enabled) then
                if (present(caching)) then
                    m_caching = caching
                else
                    m_caching = .true.
                end if
            else
                m_caching = .false.
            end if

            if (m_caching) then
                do i=1, size(cache_coefficients_p_0_q_0_plane_wave)
                    if ((cache_coefficients_p_0_q_0_plane_wave(i)%alpha .eq. alpha) &
                        .and. (cache_coefficients_p_0_q_0_plane_wave(i)%beta .eq. beta) &
                        .and. (cache_coefficients_p_0_q_0_plane_wave(i)%n_max .le. n_range(2))) then
                        cache_no = i
                    else
                        cache_no = 0
                    end if
                end do
            else
                cache_no = 0
            end if

            ! --- init ---
!            if (alpha .eq. 0.d0) then
!                allocate(p%item(n_range(1):n_range(2)))
!                allocate(q%item(n_range(1):n_range(2)))
!
!                do i=n_range(1), n_range(2)
!                    allocate(p%item(i)%item(-1:1))
!                    allocate(q%item(i)%item(-1:1))
!                end do
!            else
                call init_list(p, n_range(1), n_range(2)-n_range(1)+1)
                call init_list(q, n_range(1), n_range(2)-n_range(1)+1)
!            end if

            ! --- pre-calc ---
            if (cache_no .eq. 0) then
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
                        p = cache_coefficients_p_0_q_0_plane_wave(cache_no)%p_0%item(n)%item(m)
                        q = cache_coefficients_p_0_q_0_plane_wave(cache_no)%q_0%item(n)%item(m)
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

        end subroutine get_p_q_j_j_plane_wave_core
end module lib_mie_illumination
