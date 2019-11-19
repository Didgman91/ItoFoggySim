module lib_constants
    implicit none

    ! --- parameter ---

    ! speed of light in vacuum
    ! https://physics.nist.gov/cgi-bin/cuu/Value?c
    integer, parameter :: const_c0 =  299792458_4 ! m/s

    ! vacuum electric permittivity
    ! https://physics.nist.gov/cgi-bin/cuu/Value?ep0|search_for=permitivity
    double precision, parameter :: const_epsilon_0 = 8.8541878128d-12! [F m-1]

    ! characteristic impedance of vacuum
    ! https://physics.nist.gov/cgi-bin/cuu/Value?z0|search_for=impedance
    double precision, parameter :: const_z_0 = 376.730313668d0 ! [Ohm]

end module lib_constants
