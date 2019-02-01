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
        !
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
        !   point_x: real
        !       normalised single precision floating point number (0.0 .. 1.0)
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
        real, intent (in) :: point_x
        integer(kind=1), intent (in) :: l
        integer(kind=4) :: n

        ! auxiliary
        integer :: i
        integer(kind=1) :: buffer

        integer(kind=4) :: coordinate_binary
        integer, parameter :: COORDINATE_BINARY_NUMBER_OF_BYTES = 4

        coordinate_binary = get_coordinate_binary_number_1D_float(point_x)

!        ! use these lines instead to load all binary coordinates into the *point_n* variable.
!        do i = 1, COORDINATE_BINARY_NUMBER_OF_BYTES*8
!            if (btest(coordinate_binary, i-1)) then ! bit number starts at 0
!                point_n(COORDINATE_BINARY_NUMBER_OF_BYTES*8-i+1) = 1
!            else
!                point_n(COORDINATE_BINARY_NUMBER_OF_BYTES*8-i+1) = 0
!            end if
!        end do
!
!        n = 0
!        do i = 1, l
!            n = n +  2**(l-i) * point_n(i)
!        end do

        n = 0
        do i = 1, l
            if (btest(coordinate_binary, COORDINATE_BINARY_NUMBER_OF_BYTES*8-i)) then ! bit number starts at 0
                buffer = 1
            else
                buffer = 0
            end if

            n = n +  (2**fmm_dimensions)**(l-i) * buffer
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
        logical :: bit
        byte :: mask_bit_0_true
        data mask_bit_0_true/ B'00000001' /

        integer(kind=1), parameter :: INT_MIN_ABS = 127
        integer(kind=1), parameter :: SHIFT_EXPONENT = 1
        integer(kind=1), parameter :: SHIFT_MANTISSA = 9
        integer(kind=1) :: shift

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
        f_exponent = ishft(f_byte(4), SHIFT_EXPONENT)

        if (btest(f_byte(3), 7)) then
            bit = .true.
        else
            bit = .false.
        end if

        if (bit) then
            f_exponent = or(f_exponent, mask_bit_0_true)
        else
            f_exponent = and(f_exponent, not(mask_bit_0_true))
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

        f_mantissa = ishft(f_integer_buffer, SHIFT_MANTISSA-1)
        f_mantissa = ibset(f_mantissa, 31)  ! set virtual 1 of mantissa

        ! --- convert exponent ---
        ! binary representation: sigend integer, but interpreted as unsigned integer
        ! e.g.: -2 (base 10) = 0111 1110 (base 2) displayed as 125 (base 10)
        !
        ! -> f_exponent = INT_MIN_ABS - f_exponent
        !               = 127 - f_exponent
        !

        ! --- generate binary coordinate (only the decimal place) ---
        shift = INT_MIN_ABS-f_exponent-1    ! -1: to respect the virtual 1 of mantissa
        coordinate_binary = ishft(f_mantissa, -shift)


    end function get_coordinate_binary_number_1D_float


end module lib_octree_helper_functions
