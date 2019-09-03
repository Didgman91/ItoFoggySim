module lib_math_public
    use lib_math_auxiliaries
    use lib_math_type
    use lib_math_type_operator
    use lib_math_bessel
    use lib_math_legendre
    use lib_math_factorial
    use lib_math_wigner
    implicit none

    double precision, parameter :: PI=4.D0*atan(1.D0)   ! maximum precision, platform independet

    double precision, parameter :: unit_m = 1.0_8
    double precision, parameter :: unit_mm = 10.0_8**(-3)
    double precision, parameter :: unit_mu = 10.0_8**(-6)
    double precision, parameter :: unit_nm = 10.0_8**(-9)
end module lib_math_public
