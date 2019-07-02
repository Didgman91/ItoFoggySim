module lib_math_public
    use lib_math_type
    use lib_math_type_operator
    use lib_math_bessel
    use lib_math_legendre
    implicit none

    double precision, parameter :: PI=4.D0*atan(1.D0)   ! maximum precision, platform independet
end module lib_math_public
