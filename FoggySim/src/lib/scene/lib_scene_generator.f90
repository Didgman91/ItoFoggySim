module lib_scene_generator
    use libmath
    use lib_scene_type
    implicit none

    private


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

            type(lib_scene_object_hcp_cuboid) :: rv

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
        !       radius of the spheres at the lattice coordinate points
        !
        ! Returns
        ! ----
        !   rv: type(lib_scene_object_hcp_cuboid)
        !       an arangement of sphers with a hexagonal cloase packing (hcp) in a cuboid
        !
        !            <- size_x ->
        !            . . o o . .  ^
        !             . o o o . . |
        !          z . o o o o .  | size_z
        !           ^ . o o o . . |
        !       K_o |. . o o . .  v
        !           --> x
        !
        !   K_o: object coordinate system
        !   o: spheres at the hcp_lattice_coordiantes inside the sphere with the radius "sphere_radius"
        !   .: spheres at the hcp_lattice_coordiantes outside the sphere with the radius "sphere_radius"
        !
        function lib_scene_generator_hcp_fill_sphere(sphere_radius, lattice_sphere_radius, h) result(rv)
            implicit none
            ! dummy
            double precision, intent(in) :: sphere_radius
            double precision, intent(in) :: lattice_sphere_radius

            double precision, dimension(6), intent(in), optional :: h

            type(lib_scene_object_hcp_sphere) :: rv

            ! auxiliary
            double precision :: size_x
            double precision :: size_y
            double precision :: size_z
            type(lib_scene_object_hcp_cuboid) :: hcp_cuboid

            hcp_cuboid = lib_scene_generator_hcp_lattice_fill_cuboid(size_x, size_y, size_z, sphere_radius)

        end function lib_scene_generator_hcp_fill_sphere
end module lib_scene_generator
