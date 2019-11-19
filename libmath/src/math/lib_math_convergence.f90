module lib_math_convergence
    implicit none

    private

    ! public functions
    public :: lib_math_convergence_root_test

    contains

    ! Argument
    ! ----
    !   a_n: double precision, dimension(:)
    !       serie to test
    !       NOTE: 1 <= lbound(a_n)
    !
    ! Results
    ! ----
    !   rv: integer
    !       -1: series may converge or diverge
    !        0: series diverges
    !       >0: series converges absolutely (first element with r < 1)
    !
    ! Reference: http://mathworld.wolfram.com/RootTest.html
    function lib_math_convergence_root_test(a_n) result(rv)
        implicit none
        ! dummy
        double precision, dimension(:), intent(in) :: a_n

        integer :: rv

        ! auxiliary
        integer :: i
        double precision, dimension(lbound(a_n, 1): ubound(a_n, 1)) :: r

        r = lib_math_convergence_root_test_core(a_n)

        if (r(ubound(a_n, 1)) .lt. 1D0) then
            do i=lbound(a_n, 1), ubound(a_n ,1)
                if (r(i) < 1) then
                    rv = i
                    exit
                end if
            end do
        else if (r(ubound(a_n, 1)) .gt. 1D0) then
            rv = 0
        else
            rv = -1
        end if

    end function

    ! Argument
    ! ----
    !   a_n: double precision, dimension(:)
    !       serie to test
    !       NOTE: 1 <= lbound(a_n)
    !
    ! Results
    ! ----
    !   r: double precision, dimension(lbound(a_n, 1): ubound(a_n, 1))
    !
    ! Reference: http://mathworld.wolfram.com/RootTest.html
    function lib_math_convergence_root_test_core(a_n) result(r)
        implicit none
        ! dummy
        double precision, dimension(:), intent(in) :: a_n

        double precision, dimension(lbound(a_n, 1): ubound(a_n, 1)) :: r

        ! auxiliary
        integer :: i

        !$OMP PARALLEL DO PRIVATE(i)
        do i=lbound(a_n, 1), ubound(a_n, 1)
            r(i) = (abs(a_n(i)))**(1D0/real(i))
        end do
        !$OMP END PARALLEL DO

    end function lib_math_convergence_root_test_core

end module lib_math_convergence
