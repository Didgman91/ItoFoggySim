module lib_octree_helper_functions
    implicit none

    private

    public :: get_universal_index

    integer, private, parameter :: octree_integer_kind = 4
    integer, private, parameter :: fmm_dimensions = 1 ! dimensions

    contains

    function lib_octree_hf_get_neighbourhood_size_S(R_c) result (k)
    implicit none
        ! "The S-expansion (10) near the center of the nth box at level l for
        !  x_i ∈ E_1 (n,l) is valid for any y in the domain E_3 (n,l)."
        !
        ! Equation
        !   k >= 1/2 (R_c * d^(1/2) - 1)   (25)
        !
        ! Reference: Data_Structures_Optimal_Choice_of_Parameters_and_C
        !
        ! Arguments
        ! ----
        !   R_c: double precision
        !
        !
        ! Returns
        ! ----
        !   k: integer
        !       the minimum number of neighbours
        !
        !
        double precision, intent (in) :: R_c
        double precision :: buffer
        integer(kind=octree_integer_kind) :: k

        buffer = 0.5 * ( R_c * sqrt(real(fmm_dimensions)) - 1 )

        k = ceiling(buffer)

    end function lib_octree_hf_get_neighbourhood_size_S

    function lib_octree_hf_get_neighbourhood_size_R(r_c) result (k)
    implicit none
        ! "The R-expansion (8) near the center of the nth box at level l for x i ∈ E_3 (n,l) 
        ! is valid for any y from the domain E_1 (n,l)."
        !
        ! Equation
        !   k >= 1/2 (1/r_c * d^(1/2) - 1)   (26)
        !
        ! Reference: Data_Structures_Optimal_Choice_of_Parameters_and_C
        !
        ! Arguments
        ! ----
        !   R_c: double precision
        !       number of the node, with n ranging form 0 to 2^(3*l) in a three-dimensional space
        !
        ! Returns
        ! ----
        !   k: integer
        !       the minimum number of neighbours
        !
        !
        double precision, intent (in) :: r_c
        double precision :: buffer
        integer(kind=octree_integer_kind) :: k

        buffer = 0.5 * (1/r_c * sqrt(real(fmm_dimensions)) -1 )

        k = ceiling(buffer)

    end function lib_octree_hf_get_neighbourhood_size_R

    function lib_octree_hf_get_neighbourhood_size(R_c1, r_c2) result (k)
    implicit none
        ! "The R-expansion (8) near the center of the nth box at level l for x i ∈ E_3 (n,l) 
        ! is valid for any y from the domain E_1 (n,l)."
        !
        ! Equation
        !   k >= 1/2 (max(1/r_c, R_c) * d^(1/2) - 1)   (26)
        !
        ! Reference: Data_Structures_Optimal_Choice_of_Parameters_and_C
        !
        ! Arguments
        ! ----
        !   R_c: double precision
        !       number of the node, with n ranging form 0 to 2^(3*l) in a three-dimensional space
        !
        ! Returns
        ! ----
        !   k: integer
        !       the minimum number of neighbours
        !
        !
        double precision, intent (in) :: R_c1
        double precision, intent (in) :: r_c2
        double precision :: buffer
        integer(kind=octree_integer_kind) :: k

        buffer = 0.5 * (max(1/r_c2, R_c1) * sqrt(real(fmm_dimensions)) -1 )

        k = ceiling(buffer)

    end function lib_octree_hf_get_neighbourhood_size

    function get_universal_index(point_x, l) result(n)
    implicit none
        ! Calculates the universal index *n* of a given normalised floating point *point_x*.
        ! Related to the level *l*.
        !
        !   n = (2**d)**(l-1)*N_1 + (2**d)**(l-2)*N_2 + ... + (2**d)*N_l-1 + N_l      (49)

        ! Reference: Data_Structures_Optimal_Choice_of_Parameters_and_C
        !
        !
        ! Example:
        !   Point |  x_10    |  x_2    |  n(l=2)
        !   --------------------------------------------
        !      1  |  0.125   | 0.001   |  2*0 + 1*0 = 0
        !      2  |  0.3125  | 0.0101  |  2*0 + 1*1 = 1
        !      3  |  0.375   | 0.011   |  2*0 + 1*1 = 1
        !      4  |  0.625   | 0.101   |  2*1 + 1*0 = 2
        !      5  |  0.825   | 0.111   |  2*1 + 1*1 = 3
        !
        !
        ! Arguments
        ! ----
        !   point_x: double precision
        !       normalised floating point number (0.0 .. 1.0)
        !       HINT: datatype real is also possible, use *_1D_float() instead
        !
        !   l: integer(kind=1)
        !       number of layers
        !
        ! Returns
        ! ----
        !   the universal index *n*.
        !
        !   n: integer(kind=4)
        !

        ! dummy arguments
        double precision, intent (in) :: point_x
        integer(kind=1), intent (in) :: l
        integer(kind=4) :: n

        ! auxiliary
        integer :: i
        integer(kind=1) :: buffer

        integer(kind=8) :: coordinate_binary
        integer, parameter :: COORDINATE_BINARY_NUMBER_OF_BYTES = 8

!        coordinate_binary = get_coordinate_binary_number_1D_float(point_x)
        coordinate_binary = get_coordinate_binary_number_1D_double(point_x)

!        do i = 1, COORDINATE_BINARY_NUMBER_OF_BYTES*8
!            if (btest(coordinate_binary, i-1)) then ! bit number starts at 0
!                point_n(COORDINATE_BINARY_NUMBER_OF_BYTES*8-i+1) = 1
!            else
!                point_n(COORDINATE_BINARY_NUMBER_OF_BYTES*8-i+1) = 0
!            end if
!        end do

        n = 0
        do i = 1, l
            if (btest(coordinate_binary, COORDINATE_BINARY_NUMBER_OF_BYTES*8-i)) then ! bit number starts at 0
                buffer = 1
            else
                buffer = 0
            end if

            n = n +  (2**fmm_dimensions)**(l-i) * buffer

!            n = n +  2**(l-i) * point_n(i)
        end do

    end function

    function get_coordinate_binary_number_1D_float(f) result (coordinate_binary)
    implicit none
        ! Calculates the binary coordinate of a normalised floating point (0..1).
        !
        !   String(n, l)=(N_1, N_2, ..., N_l), N_j=0, ..., 2**d−1, j=1, ..., l (48)
        !
        ! Reference: Data_Structures_Optimal_Choice_of_Parameters_and_C
        !
        ! Example:
        !   Point |  x_10    |  x_2
        !   --------------------------
        !      1  |  0.125   | 0.001    -> String(n,l)=(0,0,1)
        !      2  |  0.3125  | 0.0101   -> String(n,l)=(0,1,0,1)
        !      3  |  0.375   | 0.011    -> String(n,l)=(0,1,1)
        !      4  |  0.625   | 0.101    -> String(n,l)=(1,0,1)
        !      5  |  0.825   | 0.111    -> String(n,l)=(1,1,1)
        !
        !
        !
        ! Floting point bitwise structure
        !   Bit no.: 7 6 5 4   3 2 1 0
        !   Byte 4: |S|E|E|E| |E|E|E|E|
        !   Byte 3: |E|M|M|M| |M|M|M|M|
        !   Byte 2: |M|M|M|M| |M|M|M|M|
        !   Byte 1: |M|M|M|M| |M|M|M|M|
        !
        ! legend:
        !   S: algebraic sign
        !   E: Exponent
        !   M: Mantissa
        !   Point x_a: (base a)
        !
        !
        ! Arguments
        ! ----
        !   f: float
        !       normalised single precision floating point number (0.0 .. 1.0)
        !
        ! Returns
        ! ----
        !   the binary representation of the floating point number (only the decimal place).
        !
        !   coordinate_binary: 4 bytes
        !
        !


        ! dummy arguments
        real, intent(in) :: f
        integer(kind=4) :: coordinate_binary

        ! variable for the binary access
        real :: f_buffer
        byte, dimension(4) :: f_byte
        equivalence (f_buffer, f_byte)

        byte :: f_exponent
        integer(kind=4) :: f_mantissa
        integer(kind=4) :: f_integer_buffer
        equivalence (f_integer_buffer, f_byte)

        ! auxiiary variables
        integer(kind=4) :: shift

        ! parametres
        integer(kind=1), parameter :: INTEGER_SIGNED_MIN_ABS = 127
        integer(kind=1), parameter :: BITS_SIGN = 1
        integer(kind=1), parameter :: BITS_EXPONENT = 8
        integer(kind=1), parameter :: BITS_MANTISSA = 23


        f_buffer = f

        ! --- extract the exponent from byte 4 and 3 ---
        ! Bit no.: 7 6 5 4   3 2 1 0
        ! Byte 4: |S|E|E|E| |E|E|E|E|
        ! Byte 3: |E|M|M|M| |M|M|M|M|
        !
        ! legend:
        !   S: algebraic sign
        !   E: Exponent
        !   M: Mantissa
        f_exponent = ishft(f_byte(4), BITS_SIGN)

        if (btest(f_byte(3), 7)) then
            f_exponent = ibset(f_exponent, 0)
        else
            f_exponent = ibclr(f_exponent, 0)
        end if

        ! --- extract the mantissa from byte 3, 2 and 1 ---
        ! f_integer_buffer:
        !   Bit no.: 7 6 5 4   3 2 1 0
        !   Byte 4: |S|E|E|E| |E|E|E|E|
        !   Byte 3: |E|M|M|M| |M|M|M|M|
        !   Byte 2: |M|M|M|M| |M|M|M|M|
        !   Byte 1: |M|M|M|M| |M|M|M|M|
        !
        ! f_mantissa
        !   Bit no.: 7 6 5 4   3 2 1 0
        !   Byte 4: |1|M|M|M| |M|M|M|M|   bit 7: "virtual 1 of mantissa: 1.MMMMM"
        !   Byte 3: |M|M|M|M| |M|M|M|M|
        !   Byte 2: |M|M|M|M| |M|M|M|M|
        !   Byte 2: |0|0|0|0| |0|0|0|0|

        f_mantissa = ishft(f_integer_buffer, BITS_SIGN + BITS_EXPONENT - 1)
        f_mantissa = ibset(f_mantissa, 31)

        ! --- convert exponent ---
        ! binary representation: sigend integer, but interpreted as unsigned integer
        ! e.g.: -2 (base 10) = 0111 1101 (base 2) displayed as 125 (base 10)
        !
        ! -> f_exponent = INT_MIN_ABS - f_exponent
        !               = 127 - f_exponent
        !

        ! --- generate binary coordinate (only the decimal place) ---
        shift = INTEGER_SIGNED_MIN_ABS - f_exponent - 1
        coordinate_binary = ishft(f_mantissa, -shift)

    end function get_coordinate_binary_number_1D_float

    function get_coordinate_binary_number_1D_double(f) result (coordinate_binary)
    implicit none
        ! Calculates the binary coordinate of a normalised floating point (0..1).
        !
        !   String(n, l)=(N_1, N_2, ..., N_l), N_j=0, ..., 2**d−1, j=1, ..., l (48)
        !
        ! Reference: Data_Structures_Optimal_Choice_of_Parameters_and_C
        !
        ! Example:
        !   Point |  x_10    |  x_2
        !   --------------------------
        !      1  |  0.125   | 0.001    -> String(n,l)=(0,0,1)
        !      2  |  0.3125  | 0.0101   -> String(n,l)=(0,1,0,1)
        !      3  |  0.375   | 0.011    -> String(n,l)=(0,1,1)
        !      4  |  0.625   | 0.101    -> String(n,l)=(1,0,1)
        !      5  |  0.825   | 0.111    -> String(n,l)=(1,1,1)
        !
        !
        !
        ! Floting point bitwise structure
        !   Bit no.: 7 6 5 4   3 2 1 0
        !   Byte 7: |S|E|E|E| |E|E|E|E|
        !   Byte 6: |E|E|E|E| |M|M|M|M|
        !   Byte 5: |M|M|M|M| |M|M|M|M|
        !   Byte 4: |M|M|M|M| |M|M|M|M|
        !   Byte 3: |M|M|M|M| |M|M|M|M|
        !   Byte 2: |M|M|M|M| |M|M|M|M|
        !   Byte 1: |M|M|M|M| |M|M|M|M|
        !   Byte 0: |M|M|M|M| |M|M|M|M|
        !
        ! legend:
        !   S: algebraic sign
        !   E: Exponent
        !   M: Mantissa
        !   Point x_a: (base a)
        !
        !
        ! Arguments
        ! ----
        !   f: float
        !       normalised single precision floating point number (0.0 .. 1.0)
        !
        ! Returns
        ! ----
        !   the binary representation of the floating point number (only the decimal place).
        !
        !   coordinate_binary: 4 bytes
        !
        !


        ! dummy arguments
        double precision, intent(in) :: f
        integer(kind=8) :: coordinate_binary

        ! variable for the binary access
        double precision :: f_buffer
        integer(kind=8) :: f_integer_buffer
        equivalence (f_integer_buffer, f_buffer)

        integer(kind=8) :: f_exponent
        integer(kind=8) :: f_mantissa


        ! auxiiary variables
        integer(kind=8) :: shift

        ! parametres
        integer(kind=2), parameter :: EXPONENT_CONVENTION = 1023
        integer(kind=1), parameter :: BITS_SIGN = 1
        integer(kind=1), parameter :: BITS_EXPONENT = 11
        integer(kind=1), parameter :: BITS_MANTISSA = 52


        f_buffer = f

        ! --- extract the exponent from byte 7 and 6 ---
        !  Bit no.: 7 6 5 4   3 2 1 0
        !   Byte 7: |S|E|E|E| |E|E|E|E|
        !   Byte 6: |E|E|E|E| |M|M|M|M|
        !
        ! Result: f_exponent
        !  Bit no.: 7 6 5 4   3 2 1 0
        !   Byte 1: |0|0|0|0| |0|E|E|E|
        !   Byte 0: |E|E|E|E| |E|E|E|E|
        !
        ! legend:
        !   S: algebraic sign
        !   E: Exponent
        !   M: Mantissa
        f_exponent = ishft(f_integer_buffer, -BITS_MANTISSA)    ! shift right

        ! --- extract the mantissa from byte 0-6 ---
        ! f_integer_buffer:
        !   Bit no.: 7 6 5 4   3 2 1 0
        !   Byte 7: |S|E|E|E| |E|E|E|E|
        !   Byte 6: |E|E|E|E| |M|M|M|M|
        !   Byte 5: |M|M|M|M| |M|M|M|M|
        !   Byte 4: |M|M|M|M| |M|M|M|M|
        !   Byte 3: |M|M|M|M| |M|M|M|M|
        !   Byte 2: |M|M|M|M| |M|M|M|M|
        !   Byte 1: |M|M|M|M| |M|M|M|M|
        !   Byte 0: |M|M|M|M| |M|M|M|M|
        !
        ! Result: f_mantissa
        !   Bit no.: 7 6 5 4   3 2 1 0
        !   Byte 4: |1|M|M|M| |M|M|M|M|   bit 7: "virtual 1 of mantissa: 1.MMMMM"
        !   Byte 6: |M|M|M|M| |M|M|M|M|
        !   Byte 5: |M|M|M|M| |M|M|M|M|
        !   Byte 4: |M|M|M|M| |M|M|M|M|
        !   Byte 3: |M|M|M|M| |M|M|M|M|
        !   Byte 2: |M|M|M|M| |M|M|M|M|
        !   Byte 1: |M|M|M|M| |M|0|0|0|
        !   Byte 0: |0|0|0|0| |0|0|0|0|
        f_mantissa = ishft(f_integer_buffer, BITS_SIGN + BITS_EXPONENT - 1)
        f_mantissa = ibset(f_mantissa, 63)


        ! --- convert exponent ---
        ! binary representation: sigend integer, but interpreted as unsigned integer
        ! e.g.: -2 (base 10) = 0000 0011 1111 1101 (base 2) displayed as 1021 (base 10)
        !
        ! -> f_exponent = EXPONENT_CONVENTION - f_exponent
        !               = 1023 - f_exponent
        !

        ! --- generate binary coordinate (only the decimal place) ---
        shift = EXPONENT_CONVENTION - int(f_exponent) - 1
        coordinate_binary = ishft(f_mantissa, -shift)   ! right shift

    end function get_coordinate_binary_number_1D_double

end module lib_octree_helper_functions
