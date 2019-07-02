module libmath
    use lib_math_public
    implicit none

!    public :: libmath_hello
!
!    interface libmath_hello
!        module procedure hello
!    end interface

    contains

    subroutine hello
        implicit none

        print *, "hello from libmath"
    end subroutine hello

end module
