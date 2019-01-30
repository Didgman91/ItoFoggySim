program main
use lib_data_types
use lib_octree
implicit none

    type(lib_octree_element), dimension(2) :: test
    integer, dimension(5) :: universal_index
    integer(kind=1), parameter :: l = 2

    test(1)%start = 34
    test(1)%end = 978
    test(1)%parent_element = 78

    universal_index(1) = get_universal_index(0.125, l)
    universal_index(2) = get_universal_index(0.3125, l)
    universal_index(3) = get_universal_index(0.375, l)
    universal_index(4) = get_universal_index(0.625, l)
    universal_index(5) = get_universal_index(0.825, l)


end program main
