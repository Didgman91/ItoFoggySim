module lib_refractive_index
    use toolbox
    use file_io
    implicit none

!    integer, dimension()
!
!    contains
!
!        function lib_refractive_index(medium, wave_length) result (rv)
!            implicit none
!            ! dummy
!            integer, intent(in) :: medium
!            double precision, intent(in) :: wave_length
!
!            double complex :: rv
!
!            ! auxiliary
!
!            call data_interpolation(data_refractive_index_particle, 1, lambda / unit_mu, &
!                                                    data_refractive_index_particle_interpolation)
!            rv = dcmplx(data_refractive_index_particle_interpolation(2), &
!                                data_refractive_index_particle_interpolation(3))
!
!
!        end function
end module lib_refractive_index
