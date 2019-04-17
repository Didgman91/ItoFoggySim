
module toolbox
    implicit none
    public :: reallocate_1d

contains
    ! reallocate a 1-dimensional array
    !
    ! Arguments
    ! ----
    !   a: 1-dimensional integer array
    !       original array, will be replaced with the resized array
    !   n: pos. integer
    !       number of additional array elements
    !
    ! change log
    ! ----
    !    - ni_new changed to relative length (n)
    !
    ! source: https://gist.github.com/ponderomotion/3527522
    !! Daniel Fletcher 2012
    !! module for increasing array sizes dynamically
    !! currently new indices must be larger than old
    SUBROUTINE reallocate_1d(a,n)
        implicit none

        INTEGER,DIMENSION(:),ALLOCATABLE,INTENT(INOUT) :: a
        INTEGER,DIMENSION(:),ALLOCATABLE :: temp
        INTEGER,INTENT(IN) :: n
        INTEGER :: ni_old

        ni_old = SIZE(a)

        ALLOCATE(temp(ni_old+n))

        temp(1:ni_old) = a

        CALL MOVE_ALLOC(temp,a)

    END SUBROUTINE reallocate_1d

end module toolbox
