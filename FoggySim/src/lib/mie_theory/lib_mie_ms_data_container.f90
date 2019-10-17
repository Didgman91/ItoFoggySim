module lib_mie_ms_data_container
    use lib_mie_type
    implicit none

    type(lib_mie_simulation_parameter_type) :: simulation_data

    contains

        subroutine lib_mie_ms_data_container_destructor()
            implicit none

            if (allocated(simulation_data%sphere_list)) deallocate(simulation_data%sphere_list)
            if (allocated(simulation_data%sphere_parameter_list)) deallocate(simulation_data%sphere_parameter_list)
            if (allocated(simulation_data%illumination%plane_wave)) deallocate(simulation_data%illumination%plane_wave)
        end subroutine

end module lib_mie_ms_data_container
