module libmath
    use lib_math_public
    implicit none
!    public :: libmath_hello
!
!    interface libmath_hello
!        module procedure hello
!    end interface

    contains

        function test_lib_math() result (error_counter)
            implicit none
            ! dummy
            integer :: error_counter

            error_counter = 0

            error_counter = error_counter + test_lib_math_auxiliaries()
            error_counter = error_counter + lib_math_factorial_test_functions()
            error_counter = error_counter + lib_math_type_operator_test_functions()
            error_counter = error_counter + lib_math_bessel_test_functions()
            error_counter = error_counter + lib_math_legendre_test_functions()
            error_counter = error_counter + lib_math_hermite_test_functions()
            error_counter = error_counter + lib_math_wigner_test_functions()
            error_counter = error_counter + lib_math_solver_test_functions()
            error_counter = error_counter + lib_math_vector_spherical_harmonics_test_functions()

            print *, "-------------test_lib_math----------------"
            if (error_counter == 0) then
                print *, "test_lib_math tests: OK"
            else
                print *, error_counter,"test_lib_math test(s) FAILED"
            end if
            print *, "------------------------------------------"
        end function

end module
