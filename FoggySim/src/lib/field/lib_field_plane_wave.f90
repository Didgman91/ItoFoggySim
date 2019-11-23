module lib_field_plane_wave
    use libmath
    use lib_constants
    use lib_field_polarisation
    implicit none

    private

    ! public function
    public :: lib_field_plane_wave_get_field

    ! public types
    public :: lib_field_plane_wave_type

    type lib_field_plane_wave_type
        double precision :: e_field_0                   ! [V/m]
        double precision :: wave_length_0               ! [m]
        double precision :: refractive_index_medium     ! [1]
        double precision :: phase = 0
        type(jones_vector_type) :: polarisation         ! Jones vector
        double precision :: theta = 0                   ! propagation direction: polar angle [0, Pi] [rad]
        double precision :: phi = 0                     ! propagation direction: azimuthal angle [0, 2 Pi) [rad]
        ! 1: exp(-i(k*z - omega*t))
        ! 2: exp(i(k*z - omega*t))
        integer :: convention = 1
    end type

    contains

        ! Argument
        ! ----
        !   parameter: type(lib_field_gaussian_beam_hermite_type)
        !       parameter of the Hermsubroutnite-Gaussian beam
        !   evaluation_point_x: type(cartesian_coordinate_real_type)
        !       evaluation point at the beam koordinate system
        !
        !          z            y
        !           ^ theta      ^    k
        !           |<->/ k      |   /
        !           |  /         |  /^
        !           | /          | / | phi
        !           |/           |/  v
        !           ------>      ------>
        !                  x            x
        !
        ! Returns
        ! ----
        !   e_field: type(cartesian_coordinate_cmplx_type)
        !       electical field component at "evaluation_point_x"
        !   h_field: type(cartesian_coordinate_cmplx_type)
        !       magnetical field component at "evaluation_point_x"
        !
        subroutine lib_field_plane_wave_get_field(parameter, evaluation_point_x, &
                                                  e_field, h_field)
            implicit none
            ! dummy
            type(lib_field_plane_wave_type), intent(in) :: parameter
            type(cartesian_coordinate_real_type), intent(in) :: evaluation_point_x

            type(cartesian_coordinate_cmplx_type) :: e_field
            type(cartesian_coordinate_cmplx_type) :: h_field

            ! auxiliary
            type(cartresian_coordinate_rot_matrix_type) :: rot
            type(cartesian_coordinate_real_type) :: point_x

            double precision :: wave_impedance
            double complex :: field

            if (parameter%phi .ne. 0 .or. parameter%theta .ne. 0) then
                rot  = lib_math_get_matrix_rot(parameter%phi, parameter%theta, parameter%phi)
                point_x = rot * evaluation_point_x
            else
                point_x = evaluation_point_x
            end if

            field = lib_field_plane_wave_scalar(parameter%e_field_0, &
                                                parameter%wave_length_0, &
                                                parameter%refractive_index_medium, &
                                                point_x%z, &
                                                phase=parameter%phase, &
                                                convention=parameter%convention)

            ! Quantum Optics: An Introduction, Mark Fox
            ! eq. 2.25 but plane wave
            wave_impedance = const_z_0 / parameter%refractive_index_medium

            rot  = lib_math_get_matrix_rot(-PI / 2d0, 0d0, 0d0)
            h_field = rot * e_field / wave_impedance

            if (parameter%phi .ne. 0 .or. parameter%theta .ne. 0) then
                rot  = lib_math_get_matrix_rot(-parameter%phi, -parameter%theta, -parameter%phi)
                e_field = rot * e_field
                h_field = rot * h_field
            end if

        end subroutine lib_field_plane_wave_get_field

        ! Argument
        ! ----
        !   e_field_0: double precision
        !       electical field [V/m]
        !   wave_length: double precision
        !       vacuum wave length [m]
        !   n_medium: double precision
        !       refractive index of the medium
        !   phase: double precision
        !       additional phase of the plane wave
        !   convention: integer, optional (std: 1)
        !       1: exp(-i(k*z - omega*t))
        !       2: exp(i(k*z - omega*t))
        !
        ! Returns
        ! ----
        !   field: double cmplx
        !       field of a plane wave with a propagation along the z-axis at position "z"
        !
        !          z,k
        !           ^
        !        ___|___
        !        ___|___
        !           ---->
        !
        ! Reference: Experimentalphysik 2, Demtröder, eq. 7.9d
        function lib_field_plane_wave_scalar(e_field_0, wave_length_0, n_medium, z, phase, convention) &
                                             result (field)
            implicit none
            ! dummy
            double precision, intent(in) :: e_field_0
            double precision, intent(in) :: wave_length_0
            double precision, intent(in) :: n_medium
            double precision, intent(in) :: z
            double precision, intent(in), optional :: phase
            integer, intent(in), optional :: convention

            double complex :: field

            ! auxiliary
            double precision :: buffer_real

            double precision :: m_phase
            integer :: m_convention

            m_convention = 1
            if (present(convention)) m_convention = convention

            m_phase = 0
            if (present(phase)) m_phase = phase


            select case(m_convention)
                case (1)
                    ! 1: exp(-i(k*z - omega*t))
                    buffer_real = - 2d0 * PI * z * n_medium / wave_length_0 + m_phase
                case (2)
                    ! exp(i(k*z - omega*t))
                    buffer_real = 2d0 * PI * z * n_medium / wave_length_0 + m_phase

                case default
                    print *, "lib_field_plane_wave_scalar: ERROR"
                    print *, "  convention is not defined: ", m_convention
             end select

            field = e_field_0 * exp(dcmplx(0, buffer_real))

        end function lib_field_plane_wave_scalar
end module lib_field_plane_wave
