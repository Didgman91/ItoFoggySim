module lib_math_public
    use lib_math_auxiliaries
    use lib_math_bessel
    use lib_math_convergence
    use lib_math_factorial
    use lib_math_legendre
    use lib_math_hermite
    use lib_math_solver
    use lib_math_type
    use lib_math_type_operator
#ifdef __GFORTRAN__
    use lib_math_wigner
    use lib_math_vector_spherical_harmonics
#endif
    use lib_math_constants
    use lib_hash_function
    implicit none

end module lib_math_public
