program main
use lib_data_types
use lib_octree
implicit none

    type(lib_octree_element), dimension(2) :: test
    integer(kind=4), dimension(5) :: universal_index
    integer(kind=4), dimension(5) :: parent_index
    integer(kind=4), dimension(2**1,5) :: children_all_index
    double precision, dimension(5) :: point_c
    integer(kind=4), dimension(2,5) :: neighbours

    integer(kind=1), parameter :: l = 3
    integer(kind=1), parameter :: k = 1

    integer(kind=1) :: i

    real, dimension(3) :: f
    integer(kind=16) :: coord


    f(1) = 0.5
    f(2) = 0.25
    f(3) = 0.125
    coord = lib_octree_hf_get_coordinate_binary_number_3D_float(f)

    test(1)%start = 34
    test(1)%end = 978
    test(1)%parent_element = 78

    universal_index(1) = lib_octree_hf_get_universal_index(0d+0, l)
    universal_index(1) = lib_octree_hf_get_universal_index(1d+0, l)

    universal_index(1) = lib_octree_hf_get_universal_index(0.125d+0, l)
    universal_index(2) = lib_octree_hf_get_universal_index(0.3125d+0, l)
    universal_index(3) = lib_octree_hf_get_universal_index(0.375d+0, l)
    universal_index(4) = lib_octree_hf_get_universal_index(0.625d+0, l)
    universal_index(5) = lib_octree_hf_get_universal_index(0.825d+0, l)

    do i = 1, 5
        parent_index(i) = lib_octree_hf_get_parent(universal_index(i))
        children_all_index(:,i) = lib_octree_hf_get_children_all(universal_index(i))
        point_c(i) = lib_octree_hf_get_centre_of_box(universal_index(i), l)
        neighbours(:,i) = lib_octree_hf_get_neighbour_all_1D(k, universal_index(i), l)
    end do

end program main
