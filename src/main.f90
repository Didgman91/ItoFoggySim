program main
use lib_data_types
implicit none

    type(lib_octree_element), dimension(2) :: test

    test(1)%start = 34
    test(1)%end = 978
    test(1)%parent_element = 78


end program main
