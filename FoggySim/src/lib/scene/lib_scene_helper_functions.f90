module lib_scene_helper_functions
    use libmath
    implicit none

    contains

        ! Argument
        ! ----
        !   list: type(cartesian_coordinate_real_type), dimension(:,:,:)
        !
!        function make_coordinate_list_1d(list, use_point) result(rv)
!            implicit none
!            ! dummy
!            type(cartesian_coordinate_real_type), dimension(:,:,:), allocatable, intent(in) :: list
!            logical, dimension(lbound(list, 1):ubound(list, 1), &
!                               lbound(list, 2):ubound(list, 2), &
!                               lbound(list, 3):ubound(list, 3)), intent(in), optional :: use_point
!
!            type(cartesian_coordinate_real_type), dimension(:), allocatable :: rv
!
!            ! auxiliary
!            integer :: i
!            integer :: j
!            integer :: k
!            integer :: x
!
!            integer, dimension(2):: i_range
!            integer, dimension(2):: j_range
!            integer, dimension(2):: k_range
!
!            i_range(1) = lbound(list, 1)
!            i_range(2) = ubound(list, 1)
!
!            j_range(1) = lbound(list, 2)
!            j_range(2) = ubound(list, 2)
!
!            k_range(1) = lbound(list, 3)
!            k_range(2) = ubound(list, 3)
!
!            if (present(use_point)) then
!
!                x = count(use_point)
!                allocate(rv(x))
!
!                x = 0
!                do i = i_range(1), i_range(2)
!                    do j = j_range(1), j_range(2)
!                        do k = k_range(1), k_range(2)
!                            if (use_point(i, j, k)) then
!                                x = x + 1
!                                rv(x) = list(i, j, k)
!                            end if
!                        end do
!                    end do
!                end do
!
!            else
!                allocate(rv((i_range(2) - i_range(1) + 1) &
!                            * (j_range(2) - j_range(1) + 1) &
!                            * (k_range(2) - k_range(1) + 1)))
!
!                x = 0
!                do i = i_range(1), i_range(2)
!                    do j = j_range(1), j_range(2)
!                        do k = k_range(1), k_range(2)
!                            x = x + 1
!                            rv(x) = list(i, j, k)
!                        end do
!                    end do
!                end do
!
!            end if
!
!        end function make_coordinate_list_1d

        function make_coordinate_list_1d(list, use_point) result(rv)
            implicit none
            ! dummy
            type(cartesian_coordinate_real_type), dimension(:), intent(in) :: list
            logical, dimension(lbound(list, 1):ubound(list, 1)), intent(in), optional :: use_point

            type(cartesian_coordinate_real_type), dimension(:), allocatable :: rv

            ! auxiliary
            integer :: i
            integer :: x

            integer, dimension(2):: i_range

            i_range(1) = lbound(list, 1)
            i_range(2) = ubound(list, 1)

            if (present(use_point)) then

                x = count(use_point)
                allocate(rv(x))

                x = 0
                do i = i_range(1), i_range(2)
                    if (use_point(i)) then
                        x = x + 1
                        rv(x) = list(i)
                    end if
                end do

            else
                allocate(rv((i_range(2) - i_range(1) + 1)))

                x = 0
                do i = i_range(1), i_range(2)
                    x = x + 1
                    rv(x) = list(i)
                end do
            end if

        end function make_coordinate_list_1d

end module lib_scene_helper_functions
