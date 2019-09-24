module lib_field
    use libmath
    use file_io
    implicit none

    contains

        function lib_field_export(e_field, h_field, path) result(rv)
            implicit none
            ! dummy
            type(cartesian_coordinate_cmplx_type), dimension(:,:), intent(in) :: e_field
            type(cartesian_coordinate_cmplx_type), dimension(lbound(e_field, 1):ubound(e_field, 1), &
                                                             lbound(e_field, 2):ubound(e_field, 2)),&
                                                   intent(in) :: h_field
            character(len=*), intent(in) :: path

            logical :: rv

            ! auxiliary
            integer :: i
            integer :: ii

            integer :: u

            type(cartesian_coordinate_cmplx_type), dimension(:,:), allocatable :: h_field_conjg

            double precision, dimension(:, :), allocatable :: e_field_real_x
            double precision, dimension(:, :), allocatable :: e_field_real_y
            double precision, dimension(:, :), allocatable :: e_field_real_z

            double precision, dimension(:, :), allocatable :: h_field_real_x
            double precision, dimension(:, :), allocatable :: h_field_real_y
            double precision, dimension(:, :), allocatable :: h_field_real_z

            type(cartesian_coordinate_real_type), dimension(: , :), allocatable :: poynting
            type(cartesian_coordinate_cmplx_type) :: buffer_cartesian_cmplx
            double precision, dimension(:, :), allocatable :: poynting_abs

            allocate(e_field_real_x(lbound(e_field, 1):ubound(e_field, 1), &
                                    lbound(e_field, 2):ubound(e_field, 2)))
            allocate(poynting(lbound(e_field, 1):ubound(e_field, 1), &
                              lbound(e_field, 2):ubound(e_field, 2)))

            allocate(e_field_real_y, mold=e_field_real_x)
            allocate(e_field_real_z, mold=e_field_real_x)
            allocate(h_field_real_x, mold=e_field_real_x)
            allocate(h_field_real_y, mold=e_field_real_x)
            allocate(h_field_real_z, mold=e_field_real_x)

            allocate(poynting_abs, mold=e_field_real_x)

            allocate(h_field_conjg, mold=h_field)

            do i = lbound(e_field, 1), ubound(e_field, 1)
                do ii = lbound(e_field, 2), ubound(e_field, 2)
                    e_field_real_x(i, ii) = real(e_field(i, ii)%x)
                    e_field_real_y(i, ii) = real(e_field(i, ii)%y)
                    e_field_real_z(i, ii) = real(e_field(i, ii)%z)

                    h_field_real_x(i, ii) = real(h_field(i, ii)%x)
                    h_field_real_y(i, ii) = real(h_field(i, ii)%y)
                    h_field_real_z(i, ii) = real(h_field(i, ii)%z)

                    ! calculate the Poynting vector: S = E x H*
                    ! eq. 43
                    h_field_conjg(i, ii)%x = conjg(h_field(i, ii)%x)
                    h_field_conjg(i, ii)%y = conjg(h_field(i, ii)%y)
                    h_field_conjg(i, ii)%z = conjg(h_field(i, ii)%z)

                    buffer_cartesian_cmplx = cross_product(e_field(i, ii), h_field_conjg(i, ii))

                    poynting(i, ii)%x = real(buffer_cartesian_cmplx%x) / 2D0
                    poynting(i, ii)%y = real(buffer_cartesian_cmplx%y) / 2D0
                    poynting(i, ii)%z = real(buffer_cartesian_cmplx%z) / 2D0

                    poynting_abs(i, ii) = abs(poynting(i, ii))
                end do
            end do

            ! --- wirte to PPM ---
            ! e field
            u = 99
            open(unit=u, file= trim(path) // "e_field_x.ppm", status='unknown')
            rv = write_ppm_p3(u, e_field_real_x)
            close(u)

            open(unit=u, file=trim(path) // "e_field_y.ppm", status='unknown')
            rv = write_ppm_p3(u, e_field_real_y)
            close(u)

            open(unit=u, file= trim(path) // "e_field_z.ppm", status='unknown')
            rv = write_ppm_p3(u, e_field_real_z)
            close(u)

            u = 99
            open(unit=u, file= trim(path) // "e_field_x_log.ppm", status='unknown')
            rv = write_ppm_p3(u, e_field_real_x, logarithmic = .true.)
            close(u)

            open(unit=u, file= trim(path) // "e_field_y_log.ppm", status='unknown')
            rv = write_ppm_p3(u, e_field_real_y, logarithmic = .true.)
            close(u)

            open(unit=u, file= trim(path) // "e_field_z_log.ppm", status='unknown')
            rv = write_ppm_p3(u, e_field_real_z, logarithmic = .true.)
            close(u)

            ! h field
            u = 99
            open(unit=u, file= trim(path) // "h_field_x.ppm", status='unknown')
            rv = write_ppm_p3(u, h_field_real_x)
            close(u)

            open(unit=u, file= trim(path) // "h_field_y.ppm", status='unknown')
            rv = write_ppm_p3(u, h_field_real_y)
            close(u)

            open(unit=u, file= trim(path) // "h_field_z.ppm", status='unknown')
            rv = write_ppm_p3(u, h_field_real_z)
            close(u)

            u = 99
            open(unit=u, file= trim(path) // "h_field_x_log.ppm", status='unknown')
            rv = write_ppm_p3(u, h_field_real_x, logarithmic = .true.)
            close(u)

            open(unit=u, file= trim(path) // "h_field_y_log.ppm", status='unknown')
            rv = write_ppm_p3(u, h_field_real_y, logarithmic = .true.)
            close(u)

            open(unit=u, file= trim(path) // "h_field_z_log.ppm", status='unknown')
            rv = write_ppm_p3(u, h_field_real_z, logarithmic = .true.)
            close(u)

            ! Poynting
            u = 99
            open(unit=u, file= trim(path) // "poynting_x.ppm", status='unknown')
            rv = write_ppm_p3(u, poynting(:,:)%x)
            close(u)

            open(unit=u, file= trim(path) // "poynting_y.ppm", status='unknown')
            rv = write_ppm_p3(u, poynting(:,:)%y)
            close(u)

            open(unit=u, file= trim(path) // "poynting_z.ppm", status='unknown')
            rv = write_ppm_p3(u, poynting(:,:)%z)
            close(u)

            u = 99
            open(unit=u, file= trim(path) // "poynting_x_log.ppm", status='unknown')
            rv = write_ppm_p3(u, poynting(:,:)%x, logarithmic = .true.)
            close(u)

            open(unit=u, file= trim(path) // "poynting_y_log.ppm", status='unknown')
            rv = write_ppm_p3(u, poynting(:,:)%y, logarithmic = .true.)
            close(u)

            open(unit=u, file= trim(path) // "poynting_z_log.ppm", status='unknown')
            rv = write_ppm_p3(u, poynting(:,:)%z, logarithmic = .true.)
            close(u)

            ! Poynting abs
            open(unit=u, file= trim(path) // "poynting_abs.ppm", status='unknown')
            rv = write_ppm_p3(u, poynting_abs)
            close(u)
            open(unit=u, file= trim(path) // "poynting_abs_log.ppm", status='unknown')
            rv = write_ppm_p3(u, poynting_abs, logarithmic = .true.)
            close(u)

        end function
end module lib_field
