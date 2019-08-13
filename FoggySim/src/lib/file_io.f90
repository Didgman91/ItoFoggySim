module file_io

    private

    ! public functions
    public :: operator( .f. )
    public :: file_exists
    public :: write_csv
    public :: read_csv

    ! interface
    interface operator( .f. )
      module procedure file_exists
    end interface

    interface write_csv
        module procedure write_csv_real
        module procedure write_csv_cmplx
        module procedure write_csv_real_cmplx
    end interface

    interface read_csv
        module procedure file_io_csv_run_get_double
    end interface


    ! parameter
    integer, parameter :: MAX_LINES = 65535

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
!                write(u, '(ES2.15,A,ES2.15,A,ES2.15)') c1(i), delimiter, c2(i), delimiter, c3(i)
                write(u, *) c1(i), delimiter, c2(i), delimiter, c3(i)
            else if (present(c2) .and. present(c3) .and. present(c4)) then
!                write(u, '(ES2.15,A,ES2.15,A,ES2.15,A,ES2.15)') c1(i), delimiter, c2(i), delimiter, c3(i), delimiter, c4(i)
                write(u, *) c1(i), delimiter, c2(i), delimiter, c3(i), delimiter, c4(i)
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
!                write(u, '(ES2.15,A,ES2.15,A,ES2.15,A,ES2.15)') &
                write(u, *) real(c1(i)), delimiter, aimag(c1(i)), delimiter, &
                            real(c2(i)), delimiter, aimag(c2(i))
            else if (present(c2) .and. present(c3)) then
!                write(u, '(ES2.15,A,ES2.15,A,ES2.15,A,ES2.15,A,ES2.15,A,ES2.15)') &
                write(u, *) real(c1(i)), delimiter, aimag(c1(i)), delimiter, &
                            real(c2(i)), delimiter, aimag(c2(i)), delimiter, &
                            real(c3(i)), delimiter, aimag(c3(i))
            else
!                write(u, '(ES2.15,A,ES2.15)') real(c1(i)), delimiter, aimag(c2(i))
                write(u, *) real(c1(i)), delimiter, aimag(c2(i))
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
            dummy_str = trim(dummy_str) // delimiter // trim(header(i)) // "(real)" // delimiter // trim(header(i)) // "(cmplx)"
        end do
        write(u, *) trim(dummy_str)

        do i=1, size(c1)
            if (.not. present(c3) .and. .not. present(c4)) then
                write(u, *) c1(i), delimiter, real(c2(i)), delimiter, aimag(c2(i))
            else if (present(c3) .and. .not. present(c4)) then
!                write(u, '(ES2.15,A,ES2.15,A,ES2.15,A,ES2.15,A,ES2.15)') &
                write(u, *) c1(i), delimiter, &
                            real(c2(i)), delimiter, aimag(c2(i)), delimiter, &
                            real(c3(i)), delimiter, aimag(c3(i))
            else if (present(c3) .and. present(c4)) then
!                write(u, '(ES2.15,A,ES2.15,A,ES2.15,A,ES2.15,A,ES2.15,A,ES2.15,A,ES2.15)') &
                write(u, *) c1(i), delimiter, &
                            real(c2(i)), delimiter, aimag(c2(i)), delimiter, &
                            real(c3(i)), delimiter, aimag(c3(i)), delimiter, &
                            real(c4(i)), delimiter, aimag(c4(i))
            end if
        end do

        rv = .true.

    end function

    ! reads a csv file with double entries
    !
    ! Arguments
    ! ----
    !   file_name: character(len=*)
    !       name of the csv file
    !   header_lines: integer, optional (std: 1)
    !       number of header lines
    subroutine file_io_csv_run_get_double(file_name, columns, data, header_lines, data_lines)
        implicit none
        ! dummy
        character(len=*) :: file_name
        integer :: columns
        integer, optional :: header_lines
        integer, optional :: data_lines

        double precision, dimension(:,:), allocatable, intent(inout) :: data

        ! auxiliary
        integer :: m_nlines
        integer :: m_header_lines
        integer :: m_data_lines
        integer, dimension(2) :: m_line_range

        ! --- init ---
        m_header_lines = 1
        if (present(header_lines)) then
            m_header_lines = header_lines
        end if

        m_data_lines = -1
        if (present(data_lines)) then
            m_data_lines = data_lines
        end if


        ! get number of lines
        m_nlines = file_io_csv_get_number_of_lines(file_name)

        if (m_nlines .lt. m_header_lines) then
            print *, "file_io_csv_run_get_double: ERROR"
            print *, "  no. of lines < header_lines"
            print *, "  no. of lines: ", m_nlines
            print *, "  header lines: ", m_header_lines
            return
        end if

        if (m_data_lines .gt. 0 .and.&
            m_data_lines .le. m_nlines-m_header_lines) then
            m_line_range = (/ m_header_lines + 1, m_header_lines + m_data_lines /)
        else
            m_line_range = (/ m_header_lines + 1, m_nlines /)
        end if

        ! read data
        m_nlines = file_io_csv_read_double(file_name, m_line_range, columns, data)

    end subroutine file_io_csv_run_get_double

    ! Counts the number of lines in a text file
    !
    ! Argument
    ! ----
    !   file_name: character(len=*)
    !       file name
    !
    ! Result
    ! ----
    !   lines: integer
    !       number of lines of the file
    !       specific values:
    !           -1: open statement is NOT successful
    !           -2: the maximum number of rows has been reached (nax = 65535)
    !
    function file_io_csv_get_number_of_lines(file_name) result(lines)
        implicit none
        ! dummy
        character(len=*) :: file_name

        integer :: lines

        ! auxiliary
        integer :: io
        integer :: unit

        lines = 0

        open (newunit = unit, file = file_name)

        do
            if (lines .gt. MAX_LINES) then
                lines = -2
                exit
            end if

            read(unit,*,iostat=io)


            if (io > 0) then
                print *, "file_io_csv_get_number_of_lines:"
                print *, "  Check input.  Something was wrong"
                print *, "  line: ", lines+1
                lines = -1
                exit
            else if (io < 0) then
                ! end of file
                exit
            else
                lines = lines + 1
            end if
        end do

        close(unit)

    end function

    ! Counts the number of lines in a text file
    !
    ! Argument
    ! ----
    !   unit: integer
    !       unit number, range 9-99
    !
    ! Result
    ! ----
    !   lines: integer
    !       number of lines of the file
    !       specific values:
    !           -1: open statement is NOT successful
    !           -2: the maximum number of rows has been reached (nax = 65535)
    !
    function file_io_csv_read_double(file_name, line_range, number_of_columns, data) result(lines)
        implicit none
        ! dummy
        character(len=*) :: file_name
        integer, dimension(2) :: line_range
        integer :: number_of_columns

        double precision, dimension(:,:), allocatable, intent(inout) :: data
        integer :: lines

        ! auxiliary
        integer :: unit
        integer :: io

        if (allocated(data)) then
            deallocate(data)
        end if
        allocate(data(line_range(1):line_range(2), number_of_columns))

        lines = 0

        open (newunit = unit, file = file_name)

        do
            if (lines .gt. MAX_LINES) then
                lines = -2
                exit
            end if

            if (lines .lt. line_range(1)-1) then
                ! last line < line_range(1) - 1
                read(unit,*,iostat=io)
            else if (lines .ge. line_range(2)) then
                ! last line >= line_range(2)
                exit
            else
                read(unit,*,iostat=io) data(lines+1, :)
            end if


            if (io > 0) then
                print *, "file_io_csv_get_number_of_lines:"
                print *, "  Check input.  Something was wrong"
                print *, "  line: ", lines+1
                lines = -1
                exit
            else if (io < 0) then
                ! end of file
                exit
            else
                lines = lines + 1
            end if
        end do

        close(unit)

    end function file_io_csv_read_double

end module
