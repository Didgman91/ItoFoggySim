program main
use lib_data_types
implicit none


!    type element
!      integer(kind=4)   :: start
!      integer(kind=4)   :: end
!      integer(kind=4)   :: parent_element
!    end type element

    type(octree_element), dimension(2) :: test

    test(1)%start = 34
    test(1)%end = 978
    test(1)%parent_element = 78


end program main
