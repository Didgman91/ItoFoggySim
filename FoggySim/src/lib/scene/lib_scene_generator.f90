module lib_scene_generator
    use libmath
    use lib_scene_type
    implicit none

    private

    ! public functions
    public :: lib_scene_gerator_hcp_lattice
    public ::lib_scene_generator_hcp_lattice_fill_cuboid
    public ::lib_scene_generator_hcp_lattice_fill_sphere


    contains

        ! Argument
        ! ----
        !   i, j, k: integer
        !       sphere index starting at 0 for the x-, y- and z-coordinates
        !   sphere_radius: doubel precision
        !       radius of one sphere [m]
        !
        ! Returns
        ! ----
        !   rv: type(cartesian_coordinate_real_type)
        !       coordinate point of the lattice for a hexagonal close packing
        !
        ! Reference: http://mathworld.wolfram.com/HexagonalClosePacking.html
        !            https://en.wikipedia.org/wiki/Close-packing_of_equal_spheres <-- todo: replace
        !
        function lib_scene_gerator_hcp_lattice(i, j, k, sphere_radius) result(rv)
            implicit none
            ! dummy
            integer, intent(in) :: i
            integer, intent(in) :: j
            integer, intent(in) :: k
            double precision, intent(in) :: sphere_radius

            type(cartesian_coordinate_real_type) :: rv

            rv%x = dble(2 * i + mod(dble(j + k), 2d0))
            rv%y = sqrt(3d0) * (dble(j) + mod(dble(k), 2d0) / 3d0)
            rv%z = 2d0 * sqrt(3d0) * dble(k) / 3d0

            rv = sphere_radius * rv

        end function lib_scene_gerator_hcp_lattice

        ! Argument
        ! ----
        !   size_x, size_y, size_z: doubel precision
        !       edge length of the cuboid with a hexagonal close packing scheme
        !   sphere_radius: double precision
        !       radius of the spheres at the lattice coordinate points
        !
        ! Returns
        ! ----
        !   rv: type(lib_scene_object_hcp_cuboid)
        !       an arangement of sphers with a hexagonal cloase packing (hcp) in a cuboid
        !
        !   <---      size_x     --->
        !    . . . . . . . . . . . .  ^
        !   . . . . . . . . . . . . . |
        !    . . . . . z . . . . . .  |
        !   . . . . . . ^ . . . . . . |
        !    . . . .K_o |. . . . . .  |
        !   . . . . . . --> x . . . . | size_z
        !    . . . . . . . . . . . .  |
        !   . . . . . . . . . . . . . |
        !    . . . . . . . . . . . .  |
        !   . . . . . . . . . . . . . |
        !    . . . . . . . . . . . .  v
        !
        !
        !   K_o: object coordinate system
        !   .: spheres at the hcp_lattice_coordiantes
        !
        function lib_scene_generator_hcp_lattice_fill_cuboid(size_x, size_y, size_z, sphere_radius) result(rv)
            implicit none
            ! dummy
            double precision, intent(in) :: size_x
            double precision, intent(in) :: size_y
            double precision, intent(in) :: size_z
            double precision, intent(in) :: sphere_radius

            type(lib_scene_object_hcp_cuboid_type) :: rv

            ! auxiliary
            integer :: i
            integer :: j
            integer :: k

            integer :: max_i
            integer :: max_j
            integer :: max_k

            type(cartesian_coordinate_real_type) :: offset

            ! The sphere with the index (0,0,0) has the coordinate (0,0,0) m
            ! and the cuboid has a corner at (0,0,0) m. An additional offset brings all spheres
            ! into the cuboid.
            !
            ! Position of the spheres with the first offset:
            !      . . . . . .
            !       . . . . . .
            !    z . . . . . .
            !     ^ . . . . . .
            !     |. . . . . .
            !     --> x
            !
            offset = make_cartesian(sphere_radius / 2d0, &
                                    sphere_radius / 2d0, &
                                    sphere_radius / 2d0)

            ! move the coordinate system to the centre of cuboid
            offset = offset - make_cartesian(size_x/ 2d0, &
                                             size_y/ 2d0, &
                                             size_z/ 2d0)

            max_k = int(floor((size_z - 2d0 * sphere_radius) / sphere_radius * 3d0 / (2d0 * sqrt(6d0))))
            max_j = int(floor((size_y - 2d0 * sphere_radius) / (sphere_radius * sqrt(3d0) - mod(dble(max_k), 2d0)/3d0)))
            max_i = int(floor(((size_x - 2d0 * sphere_radius) / sphere_radius - mod(dble(max_j + max_k), 2d0) ) / 2d0))

            allocate(rv%hcp_lattice_coordiantes(0:max_i, 0:max_j, 0:max_k))

            do i = 0, max_i
                do j = 0, max_j
                    do k = 0, max_k
                        rv%hcp_lattice_coordiantes(i, j, k) = lib_scene_gerator_hcp_lattice(i, j, k, sphere_radius) &
                                                              + offset
                    end do
                end do
            end do

        end function lib_scene_generator_hcp_lattice_fill_cuboid


        ! Argument
        ! ----
        !   size_x, size_y, size_z: doubel precision
        !       edge length of the cuboid with a hexagonal close packing scheme
        !   sphere_radius: double precision
        !       radius of the sphere containing the lattice spheres (long distance order)
        !   lattice_sphere_radius: double precision
        !       radius of the spheres at the lattice coordinate points
        !   spherical_cap_point: type(cartesian_coordinate_real_type)
        !       point on the plane intersecting the sphere
        !   spherical_cap_normal: type(cartesian_coordinate_real_type)
        !       plane normal vector, in this direction the lattice spheres are removed
        !
        !
        ! Returns
        ! ----
        !   rv: type(lib_scene_object_hcp_cuboid)
        !       an arangement of sphers with a hexagonal cloase packing (hcp) in a cuboid
        !
        !   <---      size_x     --->
        !    x x x x . . . . x x x x  ^
        !    x x  . . . . . . . x x x |
        !    x . . . . z . . . . . x  |
        !   x . . . . . ^ . . . . . x |
        !    . . . .K_o |. . . . . .  |
        !   . . . . . . --> x . . . . | size_z
        !    . . . . . . . . . . . .  |
        !   x . . . . . . . . . . . x |
        !    x . . . . . . . . . . x  |
        !   x x x . . . . . . . x x x |
        !    x x x x . . . . x x x x  v
        !
        !   K_o: object coordinate system
        !   .: spheres at the hcp_lattice_coordiantes inside the sphere with the radius "sphere_radius"
        !   x: spheres at the hcp_lattice_coordiantes outside the sphere with the radius "sphere_radius"
        !
        !
        ! with spherical cap
        !
        !   <---      size_x     --->
        !    x x x x . . \ x x x x x  ^
        !    x x  . . . . \ x x x x x |
        !    x . . . . z . \ x x x x  |
        !   x . . . . . ^ . \ x x x x |
        !    . . . .K_o |. . \ x x x  |
        !   . . . . . . --> x \ x x x | size_z
        !    . . . . . . . . . \ x x  |
        !   x . . . . . . . . . \ x x |
        !    x . . . . . . . . . \ x  |
        !   x x x . . . . . . . x \ x |
        !    x x x x . . . . x x x \  v
        !
        !
        function lib_scene_generator_hcp_lattice_fill_sphere(sphere_radius, lattice_sphere_radius, &
                                                             cap_plane_point, cap_plane_normal) result(rv)
            implicit none
            ! dummy
            double precision, intent(in) :: sphere_radius
            double precision, intent(in) :: lattice_sphere_radius

            type(cartesian_coordinate_real_type), intent(in), optional :: cap_plane_point
            type(cartesian_coordinate_real_type), intent(in), optional :: cap_plane_normal

            type(lib_scene_object_hcp_sphere_type) :: rv

            ! auxiliary
            integer :: i
            integer :: j
            integer :: k

            integer, dimension(2):: i_range
            integer, dimension(2):: j_range
            integer, dimension(2):: k_range

            type(spherical_coordinate_real_type) :: spherical_coord


            logical :: use_spherical_cap
            type(cartesian_coordinate_real_type) :: m_cap_plane_point
            type(cartesian_coordinate_real_type) :: m_cap_plane_normal
            double precision :: distance

            if (present(cap_plane_point)) m_cap_plane_point = cap_plane_point
            if (present(cap_plane_normal)) m_cap_plane_normal = cap_plane_normal / abs(cap_plane_normal)

            use_spherical_cap = .false.
            if (present(cap_plane_point) .and. present(cap_plane_normal)) use_spherical_cap = .true.

            rv%hcp_cuboid = lib_scene_generator_hcp_lattice_fill_cuboid(2d0 * sphere_radius, &
                                                                     2d0 * sphere_radius, &
                                                                     2d0 * sphere_radius, &
                                                                     lattice_sphere_radius)

            i_range(1) = lbound(rv%hcp_cuboid%hcp_lattice_coordiantes, 1)
            i_range(2) = ubound(rv%hcp_cuboid%hcp_lattice_coordiantes, 1)

            j_range(1) = lbound(rv%hcp_cuboid%hcp_lattice_coordiantes, 2)
            j_range(2) = ubound(rv%hcp_cuboid%hcp_lattice_coordiantes, 2)

            k_range(1) = lbound(rv%hcp_cuboid%hcp_lattice_coordiantes, 3)
            k_range(2) = ubound(rv%hcp_cuboid%hcp_lattice_coordiantes, 3)

            allocate (rv%inside_sphere(i_range(1):i_range(2), &
                                       j_range(1):j_range(2), &
                                       k_range(1):k_range(2)))

            do i = i_range(1), i_range(2)
                do j = j_range(1), j_range(2)
                    do k = k_range(1), k_range(2)
                        spherical_coord = rv%hcp_cuboid%hcp_lattice_coordiantes(i, j, k)
                        if ( spherical_coord%rho .lt. sphere_radius - lattice_sphere_radius) then
                            if (use_spherical_cap) then
                                distance = dot_product(rv%hcp_cuboid%hcp_lattice_coordiantes(i, j, k) - cap_plane_point, &
                                                       cap_plane_normal)
                                if (distance .lt. 0d0) then
                                    rv%inside_sphere(i,j,k) = .true.
                                else
                                    rv%inside_sphere(i,j,k) = .false.
                                end if
                            else
                                rv%inside_sphere(i,j,k) = .true.
                            end if
                        else
                            rv%inside_sphere(i,j,k) = .false.
                        end if
                    end do
                end do
            end do


        end function lib_scene_generator_hcp_lattice_fill_sphere
end module lib_scene_generator
