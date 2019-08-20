module lib_math_wigner
    use fwigxjpf
    implicit none

    private

    ! fwigxjpf lib functions
    public :: fwig_table_init
    public :: fwig_temp_init

    ! functions
    public :: lib_math_wigner_3j

    contains

        ! Argument
        ! ----
        !
        function lib_math_wigner_3j(j1, j2, j3, m1, m2, m3, double) result(rv)
            implicit none
            ! dummy
            integer, intent(in) :: j1
            integer, intent(in) :: j2
            integer, intent(in) :: j3
            integer, intent(in) :: m1
            integer, intent(in) :: m2
            integer, intent(in) :: m3
            logical, intent(in), optional :: double

            real(kind=8) :: rv


            if (present(double)) then
                if (double) then
                    rv = fwig3jj(j1, j2, j3, &
                                 m1, m2, m3)
                else
                    rv = fwig3jj(2_4*j1, 2_4*j2, 2_4*j3, &
                                 2_4*m1, 2_4*m2, 2_4*m3)
                end if
            else
                rv = fwig3jj(2_4*j1, 2_4*j2, 2_4*j3, &
                             2_4*m1, 2_4*m2, 2_4*m3)
            end if

        end function

end module lib_math_wigner
