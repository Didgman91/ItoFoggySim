program main
    use lib_data_types
    use lib_octree
    implicit none

    type(lib_octree_element), dimension(2) :: test
    type(lib_octree_universal_index), dimension(5) :: universal_index
    type(lib_octree_universal_index), dimension(5) :: parent_index
!    type(lib_octree_universal_index), dimension(2**fmm_dimensions,5) :: children_all_index
    double precision, dimension(5) :: point_c
    integer(kind=4), dimension(2,5) :: neighbours

    integer(kind=1), parameter :: l = 1
    integer(kind=1), parameter :: k = 1

    integer(kind=1) :: i

!    type(lib_octree_spatial_point), dimension(fmm_dimensions) :: f
    integer(kind=16), dimension(5) :: coord

    type(lib_octree_spatial_point), dimension(5) :: point
    type(lib_octree_spatial_point) :: point_buffer
    real :: f_buffer



    integer(kind=1), dimension(2) :: x
    integer(kind=1), dimension(2) :: buffer

    x(1) = 2
    x(2) = 0

    buffer = lib_octree_hf_interleave_bits_use_lut(x)

    buffer = lib_octree_hf_interleave_bits(x)

!    f(1) = 0.875
!    f(2) = 0.875
!    f(3) = 0.875
!    coord = lib_octree_hf_get_coordinate_binary_number_3D_float(f)


!    test(1)%start = 34
!    test(1)%end = 978
!    test(1)%parent_element = 78

!    point(1)%x = (/0.125, 0.125, 0.125/)
!    point(2)%x = (/0.25, 0.25, 0.25/)
!    point(3)%x = (/0.5, 0.5, 0.5/)
!    point(4)%x = (/0.75, 0.75, 0.75/)
!    point(5)%x = (/0.875, 0.875, 0.875/)

!    point(1)%x = (/0.125, 0.0, 0.0/)
!    point(2)%x = (/0.25, 0.0, 0.0/)
!    point(3)%x = (/0.5, 0.0, 0.0/)
!    point(4)%x = (/0.75, 0.0, 0.0/)
!    point(5)%x = (/0.875, 0.0, 0.0/)

!    point(1)%x = (/0.0, 0.125, 0.0/)
!    point(2)%x = (/0.0, 0.25, 0.0/)
!    point(3)%x = (/0.0, 0.5, 0.0/)
!    point(4)%x = (/0.0, 0.75, 0.0/)
!    point(5)%x = (/0.0, 0.875, 0.0/)

!    point(1)%x = (/0.0, 0.0, 0.125/)
!    point(2)%x = (/0.0, 0.0, 0.25/)
!    point(3)%x = (/0.0, 0.0, 0.5/)
!    point(4)%x = (/0.0, 0.0, 0.75/)
!    point(5)%x = (/0.0, 0.0, 0.875/)

!   !1D
!   point(1)%x = (/0.125/)
!   point(2)%x = (/0.25/)
!   point(3)%x = (/0.5/)
!   point(4)%x = (/0.75/)
!   point(5)%x = (/0.875/)

!    do i = 1, 5
!        universal_index(i) = lib_octree_hf_get_universal_index(point(i), l)
!    end do

!
!    do i = 1, 5
!        parent_index(i) = lib_octree_hf_get_parent(universal_index(i))
!        children_all_index(:,i) = lib_octree_hf_get_children_all(universal_index(i))
!        point(i) = lib_octree_hf_get_centre_of_box(universal_index(i), l)
!        neighbours(:,i) = lib_octree_hf_get_neighbour_all_1D(k, universal_index(i), l)
!    end do


end program main
