module lib_mie_multi_sphere
    use lib_mie_type
    use lib_mie_single_sphere
    implicit none

    private

    contains

        ! Arguments
        ! ----
        !   k: type(spherical_coordinate_real_type)
        !       ?
        !   lambda: double precision
        !       vacuum wave length
        !   n_medium: double precision
        !       refractive index of the medium
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
        !   a_old, b_old: type(list_list_cmplx), optional
        !       coefficient of the previous calculation, if f < 1
        !   f: double precision, optional (std: 1)
        !       numerical factor (0, 1]
        !       "In our actual calculations, some multi-
        !        sphere systems do not converge if f 5 1, but they do
        !        converge when the value of f is reduced to, say, 0.7."[1]
        ! Returns
        ! ----
        !
        !
        !
        !
        ! Reference: [1] Electromagnetic scattering by an aggregate of spheres, Yu-lin Xu, eq. 30
        subroutine lib_mie_ms_get_interactive_scattering_coefficients(k, lambda, n_medium, &
                                                                   sphere, sphere_parameter, sphere_j, &
                                                                   z_selector, &
                                                                   a, b)
            implicit none
            ! dummy
            type(spherical_coordinate_real_type), intent(in) :: k
            double precision, intent(in) :: lambda
            double precision, intent(in) :: n_medium
            type(lib_mie_sphere_type), dimension(:), intent(in) :: sphere
            type(lib_mie_sphere_parameter_type), dimension(:), intent(in) :: sphere_parameter
            integer :: sphere_j
            integer(kind=1) :: z_selector

            type(list_list_cmplx), intent(inout) :: a
            type(list_list_cmplx), intent(inout) :: b

            ! auxiliary
            integer :: mu
            integer :: nu


        end subroutine lib_mie_ms_get_interactive_scattering_coefficients
end module lib_mie_multi_sphere
