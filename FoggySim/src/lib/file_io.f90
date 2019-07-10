module file_io

    interface operator( .f. )
      module procedure file_exists
    end interface

    interface write_csv
        module procedure write_csv_real
        module procedure write_csv_cmplx
        module procedure write_csv_real_cmplx
    end interface

    contains

    function file_exists(filename) result(res)
      implicit none
      character(len=*),intent(in) :: filename
      logical                     :: res

      ! Check if the file exists
      inquire( file=trim(filename), exist=res )
    end function

    ! writes up to three columns into a csv file
    !
    ! Argument
    !   u: integer
    !       >>> open(unit=u, file=trim(filename), status='new')
    !   header:
    !
    !
    function write_csv_real(u, header, c1, c2, c3, c4) result (rv)
        implicit none
        ! dummy
        integer, intent(in) :: u
        character(len=*), dimension(:) :: header
        double precision, dimension(:) :: c1
        double precision, dimension(size(c1)), optional :: c2
        double precision, dimension(size(c1)), optional :: c3
        double precision, dimension(size(c1)), optional :: c4

        logical :: rv

        ! parameter
        character, parameter :: delimiter = ','

        ! auxiliary
        integer :: i
        character(len=511) :: dummy_str

        dummy_str = trim(header(1))
        do i=2, size(header)
            dummy_str = trim(dummy_str) // delimiter // trim(header(i))
        end do
        write(u, *) trim(dummy_str)

        do i=1, size(c1)
            if (present(c2) .and. .not. present(c3) .and. .not. present(c4)) then
!                write(u, '(F2.15,A,ES2.15)') c1(i), delimiter, c2(i)
                write(u, *) c1(i), delimiter, c2(i)
            else if (present(c2) .and. present(c3) .and. .not. present(c4)) then
                write(u, '(ES2.15,A,ES2.15,A,ES2.15)') c1(i), delimiter, c2(i), delimiter, c3(i)
            else if (present(c2) .and. present(c3) .and. present(c4)) then
                write(u, '(ES2.15,A,ES2.15,A,ES2.15,A,ES2.15)') c1(i), delimiter, c2(i), delimiter, c3(i), delimiter, c4(i)
            else
                write(u, '(ES2.15)') c1(i)
            end if
        end do

        rv = .true.

    end function

    ! writes up to three columns into a csv file
    !
    ! Argument
    !   u: integer
    !       >>> open(unit=u, file=trim(filename), status='new')
    !   header:
    !
    !
    function write_csv_cmplx(u, header, c1, c2, c3) result (rv)
        implicit none
        ! dummy
        integer, intent(in) :: u
        character(len=*), dimension(:) :: header
        complex(kind=8), dimension(:) :: c1
        complex(kind=8), dimension(size(c1)), optional :: c2
        complex(kind=8), dimension(size(c1)), optional :: c3

        logical :: rv

        ! parameter
        character, parameter :: delimiter = ','

        ! auxiliary
        integer :: i
        character(len=511) :: dummy_str

        dummy_str = trim(header(1))
        do i=2, size(header)
            dummy_str = trim(dummy_str) // delimiter // trim(header(i))
        end do
        write(u, *) trim(dummy_str)

        do i=1, size(c1)
            if (present(c2) .and. .not. present(c3)) then
                write(u, '(ES2.15,A,ES2.15,A,ES2.15,A,ES2.15)') &
                                real(c1(i)), delimiter, aimag(c1(i)), delimiter, &
                                real(c2(i)), delimiter, aimag(c2(i))
            else if (present(c2) .and. present(c3)) then
                write(u, '(ES2.15,A,ES2.15,A,ES2.15,A,ES2.15,A,ES2.15,A,ES2.15)') &
                                real(c1(i)), delimiter, aimag(c1(i)), delimiter, &
                                real(c2(i)), delimiter, aimag(c2(i)), delimiter, &
                                real(c3(i)), delimiter, aimag(c3(i))
            else
                write(u, '(ES2.15,A,ES2.15)') real(c1(i)), delimiter, aimag(c2(i))
            end if
        end do

        rv = .true.

    end function

    ! writes up to three columns into a csv file
    !
    ! Argument
    !   u: integer
    !       >>> open(unit=u, file=trim(filename), status='new')
    !   header:
    !
    !
    function write_csv_real_cmplx(u, header, c1, c2, c3, c4) result (rv)
        implicit none
        ! dummy
        integer, intent(in) :: u
        character(len=*), dimension(:) :: header
        real(kind=8), dimension(:) :: c1
        complex(kind=8), dimension(size(c1)) :: c2
        complex(kind=8), dimension(size(c1)), optional :: c3
        complex(kind=8), dimension(size(c1)), optional :: c4

        logical :: rv

        ! parameter
        character, parameter :: delimiter = ','

        ! auxiliary
        integer :: i
        character(len=511) :: dummy_str

        dummy_str = trim(header(1))
        do i=2, size(header)
            dummy_str = trim(dummy_str) // delimiter // trim(header(i))
        end do
        write(u, *) trim(dummy_str)

        do i=1, size(c1)
            if (.not. present(c3) .and. .not. present(c4)) then
                write(u, '(ES2.15,A,ES2.15,A,ES2.15)') &
                                c1(i), delimiter, &
                                real(c2(i)), delimiter, aimag(c2(i))
            else if (present(c3) .and. .not. present(c4)) then
                write(u, '(ES2.15,A,ES2.15,A,ES2.15,A,ES2.15,A,ES2.15)') &
                                c1(i), delimiter, &
                                real(c2(i)), delimiter, aimag(c2(i)), delimiter, &
                                real(c3(i)), delimiter, aimag(c3(i))
            else if (present(c3) .and. present(c4)) then
                write(u, '(ES2.15,A,ES2.15,A,ES2.15,A,ES2.15,A,ES2.15,A,ES2.15,A,ES2.15)') &
                                c1(i), delimiter, &
                                real(c2(i)), delimiter, aimag(c2(i)), delimiter, &
                                real(c3(i)), delimiter, aimag(c3(i)), delimiter, &
                                real(c4(i)), delimiter, aimag(c4(i))
            end if
        end do

        rv = .true.

    end function

end module
