program main
    use file_io
    use libmath
    use lib_data_types
    use lib_tree
    use lib_tree_helper_functions
    use lib_hash_function
    use lib_sort

    use lib_ml_fmm
    use ml_fmm_math
    use lib_ml_fmm_type_operator

    use lib_mie_vector_spherical_harmonics

    use lib_mie_scattering_by_a_sphere

    use light_scattering
    implicit none

    integer :: error_counter

    integer(kind=4) :: l
    integer(kind=4) :: m
    double precision :: x
    double precision :: erg

    ! CPU-time
    real :: test_start, test_finish
    ! WALL-time
    INTEGER :: test_count_start, test_count_finish, test_count_rate


    l = 3
    m = 1
    x = 0.01_8

!    polarization: parallel (1) perpendicular (2)
!    size parameter: x
!    index of refraction: real,imaginary (+ for absorption)
!    npnts: number of grid points (npnts by npnts)
!
!    e.g.:
!        ip = 1        ! polarisation
!        x = 2*pi * 10 ! sphere radius: a=10um;  wave length: 1um; a / wave length = 10
!        cmr = 1.33    ! water: 1.33
!        cmi = 0       ! water: 0
!        npnts = 10
!   call S2(1, 2*3.14159265358979*10, 1.33, 0.0, 10)
    call S2(1, 10.0, 1.50, 0.0, 10)

    call test_file_io

    call system_clock(test_count_start, test_count_rate)
    call cpu_time(test_start)

    error_counter = 0
    error_counter = error_counter + test_lib_math()
    error_counter = error_counter + lib_mie_vector_spherical_harmonics_test_functions()
    error_counter = error_counter + lib_mie_scattering_by_a_sphere_test_functions()
    error_counter = error_counter + lib_sort_test_functions()
    error_counter = error_counter + lib_test_hash_function()
    error_counter = error_counter + lib_tree_hf_test_functions()
    error_counter = error_counter + lib_tree_test_functions()
    error_counter = error_counter + lib_ml_fmm_type_operator_test_functions()
    error_counter = error_counter + lib_ml_fmm_hf_test_functions()
    error_counter = error_counter + lib_ml_fmm_test_functions()

    call cpu_time(test_finish)
    call system_clock(test_count_finish, test_count_rate)

    print *, ""
    print *, "-------------MAIN------------------"
    print '("  CPU-Time = ",f10.3," seconds.")',test_finish-test_start
    print '("  WALL-Time = ",f10.3," seconds.")',(test_count_finish-test_count_start) / real(test_count_rate)
    print *, ""
    if (error_counter == 0) then
        print *, "All tests: OK"
    else
        print *, error_counter,"test(s) FAILED"
    end if
    print *, "-----------------------------------"

!    call lib_tree_benchmark()
!     call lib_tree_hf_benchmark()

    call lib_tree_destructor()
!    call lib_tree_hf_destructor()

end program main
