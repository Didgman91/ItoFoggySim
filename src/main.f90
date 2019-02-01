program main
use lib_data_types
use lib_octree
implicit none

    type(lib_octree_element), dimension(2) :: test
    integer, dimension(5) :: universal_index
    integer(kind=1), parameter :: l = 3

    test(1)%start = 34
    test(1)%end = 978
    test(1)%parent_element = 78

    universal_index(1) = get_universal_index(0d+0, l)
    universal_index(1) = get_universal_index(1d+0, l)

    universal_index(1) = get_universal_index(0.125d+0, l)
    universal_index(2) = get_universal_index(0.3125d+0, l)
    universal_index(3) = get_universal_index(0.375d+0, l)
    universal_index(4) = get_universal_index(0.625d+0, l)
    universal_index(5) = get_universal_index(0.825d+0, l)


end program main
