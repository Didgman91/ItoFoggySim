! File: lib_octree_helper_functions.F90
!
! Created on Tue Feb  5 16:33:05 2019
!
! @author: itodaiber
!
!
! Notes
! ----
!
! Preprocessor
! -----
!   standard value: _FMM_DIMENSION_ = 3
!
! Eclipse settings
! ------
!   Project properties -> Fortran General -> Paths and Symbols -> Symbols
!

! spatial dimension, value = [2,3]
#define _FMM_DIMENSION_ 3

! 1: true, 0: false (-> spatial point is real)
#define _SPATIAL_POINT_IS_DOUBLE_ 1

! integer kind of the bit interleaving process, value = [1,2,4,8], default = 1
! a value of 8 is only possible if the spatial point variable is of type double
#define _INTERLEAVE_BITS_INTEGER_KIND_ 1

module lib_octree_helper_functions
    use file_io
    implicit none

    private

    ! parameter
    integer(kind=1), public, parameter :: OCTREE_DIMENSIONS = _FMM_DIMENSION_ ! dimensions
    integer(kind=1), private, parameter :: OCTREE_INTEGER_KIND = 4
    integer(kind=1), private, parameter :: NUMBER_OF_BITS_PER_BYTE = 8

#if(_SPATIAL_POINT_IS_DOUBLE_ == 1)
    integer(kind=1), parameter :: COORDINATE_BINARY_BYTES = 8
#elif(_SPATIAL_POINT_IS_DOUBLE_ == 0)
    integer(kind=1), parameter :: COORDINATE_BINARY_BYTES = 4
#endif
    ! ~ parameter ~

    ! type definitions
    type lib_octree_spatial_point
#if (_SPATIAL_POINT_IS_DOUBLE_ == 1)
        double precision, dimension(OCTREE_DIMENSIONS) :: x
#elif (_SPATIAL_POINT_IS_DOUBLE_ == 0)
        real, dimension(OCTREE_DIMENSIONS) :: x
#endif
    end type lib_octree_spatial_point

    type lib_octree_universal_index
        integer(kind=COORDINATE_BINARY_BYTES) :: n
        integer(kind=COORDINATE_BINARY_BYTES), dimension(OCTREE_DIMENSIONS) :: n_per_dimension
        integer(kind=4) :: l
    end type lib_octree_universal_index
    ! ~ type definitions ~

    ! module global variable
#if (_FMM_DIMENSION_ == 2)
    integer(kind=_INTERLEAVE_BITS_INTEGER_KIND_), dimension (:,:,:) &
                                                , allocatable :: lib_octree_interleave_bits_lut
    logical :: lib_octree_interleave_bits_lut_initialised = .false.

    integer(kind=_INTERLEAVE_BITS_INTEGER_KIND_), dimension (:,:,:) &
                                                , allocatable :: lib_octree_deinterleave_bits_lut
    logical :: lib_octree_deinterleave_bits_lut_initialised = .false.
#elif (_FMM_DIMENSION_ == 3)
    integer(kind=_INTERLEAVE_BITS_INTEGER_KIND_), dimension (:,:,:,:) &
                                                , allocatable :: lib_octree_interleave_bits_lut
    logical :: lib_octree_interleave_bits_lut_initialised = .false.

    integer(kind=_INTERLEAVE_BITS_INTEGER_KIND_), dimension (:,:,:,:) &
                                                , allocatable :: lib_octree_deinterleave_bits_lut
    logical :: lib_octree_deinterleave_bits_lut_initialised = .false.
#endif
    ! ~ module global variable ~

    public :: lib_octree_hf_destructor

    public :: lib_octree_spatial_point
    public :: lib_octree_universal_index

    public :: lib_octree_hf_get_universal_index
    public :: lib_octree_hf_get_parent

    public :: lib_octree_hf_get_children_all
    public :: lib_octree_hf_get_centre_of_box
    public :: lib_octree_hf_get_neighbour_all_1D

    ! test
    public :: lib_octree_hf_test_functions
    public :: lib_octree_hf_benchmark

contains

    ! Cleans up the memory
    subroutine lib_octree_hf_destructor()
        implicit none


#if (_FMM_DIMENSION_ == 2)
        if (lib_octree_interleave_bits_lut_initialised) then
           lib_octree_interleave_bits_lut_initialised = .false.
           deallocate (lib_octree_interleave_bits_lut)
        end if
        if (lib_octree_deinterleave_bits_lut_initialised) then
            lib_octree_deinterleave_bits_lut_initialised = .false.
           deallocate (lib_octree_deinterleave_bits_lut)
        end if
#elif (_FMM_DIMENSION_ == 3)
        if (lib_octree_interleave_bits_lut_initialised) then
            lib_octree_interleave_bits_lut_initialised = .false.
            deallocate (lib_octree_interleave_bits_lut)
        end if

        if (lib_octree_deinterleave_bits_lut_initialised) then
            lib_octree_deinterleave_bits_lut_initialised = .false.
            deallocate (lib_octree_deinterleave_bits_lut)
        end if
#endif

    end subroutine lib_octree_hf_destructor

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
    function lib_octree_hf_get_neighbourhood_size_S(R_c) result (k)
        implicit none

        double precision, intent (in) :: R_c
        double precision :: buffer
        integer(kind=OCTREE_INTEGER_KIND) :: k

        buffer = 0.5 * ( R_c * sqrt(real(OCTREE_DIMENSIONS)) - 1 )

        k = ceiling(buffer)

    end function lib_octree_hf_get_neighbourhood_size_S

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
    function lib_octree_hf_get_neighbourhood_size_R(r_c) result (k)
        implicit none

        double precision, intent (in) :: r_c
        double precision :: buffer
        integer(kind=OCTREE_INTEGER_KIND) :: k

        buffer = 0.5 * (1/r_c * sqrt(real(OCTREE_DIMENSIONS)) -1 )

        k = ceiling(buffer)

    end function lib_octree_hf_get_neighbourhood_size_R

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
    function lib_octree_hf_get_neighbourhood_size(R_c1, r_c2) result (k)
        implicit none

        double precision, intent (in) :: R_c1
        double precision, intent (in) :: r_c2
        double precision :: buffer
        integer(kind=OCTREE_INTEGER_KIND) :: k

        buffer = 0.5 * (max(1/r_c2, R_c1) * sqrt(real(OCTREE_DIMENSIONS)) -1 )

        k = ceiling(buffer)

    end function lib_octree_hf_get_neighbourhood_size

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
    !   point_x: type(lib_octree_spatial_point)
    !       normalised floating point number (0.0 .. 1.0)
    !       HINT: datatype real is also possible, use *_1D_float() instead
    !
    !   l: integer(kind=1)
    !       number of layers
    !
    ! Returns
    ! ----
    !   the universal index *uindex*.
    !
    !   uindex: type(lib_octree_universal_index)
    !
    function lib_octree_hf_get_universal_index(point_x, l) result(uindex)
        implicit none

        ! dummy arguments
        type(lib_octree_spatial_point), intent(in) :: point_x
        integer(kind=1), intent (in) :: l
        type(lib_octree_universal_index) :: uindex

        ! auxiliary
        integer(kind=COORDINATE_BINARY_BYTES), dimension(OCTREE_DIMENSIONS) :: coordinate_binary

        ! Example
        ! ----
        !
        ! coordinate_binary: kind=4, dimension=3                                            2**-1
        ! element   byte representation                          base(2)                     v        base(10)
        ! 1:       |---04---|---03---|---02---|---01---|   e.g.: |00000000|00000000|00000000|11000000|  0.75
        ! 2:       |---04---|---03---|---02---|---01---|   e.g.: |00000000|00000000|00000000|10000000|  0.5
        ! 3:       |---04---|---03---|---02---|---01---|   e.g.: |00000000|00000000|10000000|00000001|  0.005859375
        !
        ! cb_buffer: kind=1, dimension=3*4=12
        ! element   byte representation                          base(2)                                base(10)
        ! 1-4:     |---01---|---02---|---03---|---04---|   e.g.: |00000000|00000000|00000000|11000000|  0.75
        ! 5-8:     |---05---|---06---|---07---|---08---|   e.g.: |00000000|00000000|00000000|10000000|  0.5
        ! 9-12:    |---09---|---10---|---11---|---12---|   e.g.: |00000000|00000000|10000000|00000001|  0.005859375
        !                                     <-- *
        !                    interleave column wise
        !
        ! equivalence (coordinate_binary, cb_buffer)
        !
        ! Interleave bits
        ! -----
        !
        ! interleaved_bits: kind=1, dimension=3*4=12
        ! element   byte representation
        ! 1-4:     |04-08-12|04-08-12|04-08-12|03-07-11|   e.g.: |11010000|00000000|00000001|00100000|
        ! 5-8:     |03-07-11|03-07-11|02-06-10|02-06-10|   e.g.: |00000000|00000000|00000000|00000000|
        ! 9-12:    |02-06-10|01-05-09|01-05-09|01-05-09|   e.g.: |00000000|00000000|00000000|00000000|
        !
        !   xx-yy-zz describes the bytes which were interleaved
        !
        ! ib_buffer: kind=1, dimension=3
        ! element   byte representation
        ! 1-3:     |---03---|---02---|---01---|
        !
        !
        ! ib_buffer = interleave(cb_buffer(9), cb_buffer(5), cb_buffer(1))
        ! iterleaved_bits(1:3) = ib_buffer
        !
        ! ...
        !
        ! ib_buffer = interleave(cb_buffer(12), cb_buffer(8), cb_buffer(4))
        ! interleaved_bits(10:12) = ib_buffer
        !
        integer(kind=_INTERLEAVE_BITS_INTEGER_KIND_) &
            ,dimension(OCTREE_DIMENSIONS * COORDINATE_BINARY_BYTES/_INTERLEAVE_BITS_INTEGER_KIND_) &
            :: interleaved_bits
        !
        ! doubel precision
        !   Number of decimal digits: ca. 16
        !   bit precision: 53
        !
        ! 3d example
        ! ----
        ! smalest cube
        !   edge length: 10**-16
        !   volume: 10**-48
        !
        ! largest cube
        !   edge length: 1
        !   volume: 1
        !
        ! number of smalest cubes in the biggest cube
        !   number = 1 / 10**-48 = 10**48
        !
        ! determination of the integer kind
        ! -----
        ! 32 bit word range: −(2**31) to 2**31 − 1
        ! 64 bit word range: −(2**63) to 2**63 − 1
        !
        !
        integer(kind=COORDINATE_BINARY_BYTES) :: interleaved_bits_dimension_0

        !   OCTREE_DIMENSIONS * COORDINATE_BINARY_BYTES/_INTERLEAVE_BITS_INTEGER_KIND_ - COORDINATE_BINARY_BYTES /_INTERLEAVE_BITS_INTEGER_KIND_ + 1
        ! = COORDINATE_BINARY_BYTES/_INTERLEAVE_BITS_INTEGER_KIND_ * ( OCTREE_DIMENSIONS - 1) + 1
        equivalence (interleaved_bits(COORDINATE_BINARY_BYTES/_INTERLEAVE_BITS_INTEGER_KIND_ * ( OCTREE_DIMENSIONS - 1) + 1), &
                     interleaved_bits_dimension_0)

        integer(kind=_INTERLEAVE_BITS_INTEGER_KIND_) &
            ,dimension(OCTREE_DIMENSIONS * COORDINATE_BINARY_BYTES/_INTERLEAVE_BITS_INTEGER_KIND_) &
            :: cb_buffer

        equivalence (coordinate_binary, cb_buffer)

        integer(kind=_INTERLEAVE_BITS_INTEGER_KIND_), dimension(OCTREE_DIMENSIONS) :: ob_buffer_DIM
        integer(kind=_INTERLEAVE_BITS_INTEGER_KIND_), dimension(OCTREE_DIMENSIONS) :: ib_buffer

        integer(kind=1) :: i
        integer(kind=1) :: ii

        ! calculate the x-dimensional binary coordinate
        coordinate_binary = lib_octree_hf_get_coordinate_binary_number_xD(point_x%x)

        ! interleave bits
        do i=COORDINATE_BINARY_BYTES/_INTERLEAVE_BITS_INTEGER_KIND_, 1, -1 ! interleave column wise
            do ii=1, OCTREE_DIMENSIONS  ! get column entries
                ob_buffer_DIM(ii) = cb_buffer(i + (ii-1)*COORDINATE_BINARY_BYTES/_INTERLEAVE_BITS_INTEGER_KIND_)
            end do
            ib_buffer = lib_octree_hf_interleave_bits_use_lut(ob_buffer_DIM)

            ! e.g.    12                4       (4..1)        3
            ! ii = total length - (total columns - i + 1) * length(ib_buffer) + 1
            ! ii = OCTREE_DIMENSIONS * COORDINATE_BINARY_BYTES/_INTERLEAVE_BITS_INTEGER_KIND_ - (COORDINATE_BINARY_BYTES/_INTERLEAVE_BITS_INTEGER_KIND_ - i + 1) * OCTREE_DIMENSIONS + 1
            ! ii = OCTREE_DIMENSIONS * (COORDINATE_BINARY_BYTES/_INTERLEAVE_BITS_INTEGER_KIND_ - COORDINATE_BINARY_BYTES/_INTERLEAVE_BITS_INTEGER_KIND_ + i - 1) + 1
            ! ii = OCTREE_DIMENSIONS * (i - 1) + 1
            ii = int(OCTREE_DIMENSIONS * (i-1) + 1, 1)

            interleaved_bits(ii:ii+OCTREE_DIMENSIONS+1) = ib_buffer(:)
        end do

        ! calculate the universal index n
        uindex%n = ibits(interleaved_bits_dimension_0, &
                        COORDINATE_BINARY_BYTES*NUMBER_OF_BITS_PER_BYTE-l*OCTREE_DIMENSIONS, &
                        l*OCTREE_DIMENSIONS)

        do i = 1, OCTREE_DIMENSIONS
            uindex%n_per_dimension(i) = ibits(coordinate_binary(i), COORDINATE_BINARY_BYTES*NUMBER_OF_BITS_PER_BYTE-l, l)
        end do

        uindex%l = l

    end function lib_octree_hf_get_universal_index

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
    function lib_octree_hf_get_coordinate_binary_number_1D_float(f) result (coordinate_binary)
        implicit none

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
        f_mantissa = ibset(f_mantissa, 31)  ! set virtual 1 of mantissa

        ! --- convert exponent ---
        ! binary representation: sigend integer, but interpreted as unsigned integer
        ! e.g.: -2 (base 10) = 0111 1101 (base 2) displayed as 125 (base 10)
        !
        ! -> f_exponent = INT_MIN_ABS - f_exponent
        !               = 127 - f_exponent
        !

        ! --- generate binary coordinate (only the decimal place) ---
        shift = INTEGER_SIGNED_MIN_ABS-f_exponent-1    ! -1: to respect the virtual 1 of mantissa
        coordinate_binary = ishft(f_mantissa, -shift)

    end function lib_octree_hf_get_coordinate_binary_number_1D_float

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
    !   coordinate_binary: 8 bytes
    !
    function lib_octree_hf_get_coordinate_binary_number_1D_double(f) result (coordinate_binary)
        implicit none

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
        !   Byte 7: |1|M|M|M| |M|M|M|M|   bit 7: "virtual 1 of mantissa: 1.MMMMM"
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

    end function lib_octree_hf_get_coordinate_binary_number_1D_double

!    function lib_octree_hf_get_coordinate_binary_number_3D_float(f) result (coordinate_binary_3D)
!        implicit none
!        ! Calculates the binary coordinate of a normalised floating point vector (0..1, 0..1, 0..1).
!        !
!        !   x = (0.N_1 N_2 N_3 ... N_j...)_(2_d), N_j=(b_(1j) b_(2j)...b_(dj) )_2 , j=1, 2, ..., N_j=0, ..., 2_d−1. (77)
!        !
!        ! Reference: Data_Structures_Optimal_Choice_of_Parameters_and_C
!        !
!        ! Arguments
!        ! ----
!        !   f: float, dimension(3)
!        !       normalised single precision floating point number (0.0..1.0, 0.0..1.0, 0.0..1.0)
!        !
!        ! Returns
!        ! ----
!        !   the binary representation of the floating point vector (only the decimal place).
!        !
!        !   coordinate_binary_D: 16 bytes
!        !
!
!        ! parameters
!        integer(kind=1), parameter :: DIMENSION = 3         ! do not change, function is specialised in three-dimensional data points
!        integer(kind=1), parameter :: NUMBER_OF_BYTES_COORDINATE_3D = 16  ! space for 3 (=DIMENSION) integer of kind 4 (OCTREE_INTEGER_KIND)
!        integer(kind=1), parameter :: NUMBER_OF_BITS_COORDINATE_1D = 4 * NUMBER_OF_BITS_PER_BYTE
!        integer(kind=2), parameter :: NUMBER_OF_BITS_COORDINATE_3D = NUMBER_OF_BYTES_COORDINATE_3D * NUMBER_OF_BITS_PER_BYTE
!
!        ! dummy arguments
!        real, dimension(DIMENSION), intent(in) :: f
!        integer(kind=NUMBER_OF_BYTES_COORDINATE_3D) :: coordinate_binary_3D
!
!        ! auxiliary variables
!        real :: f_buffer
!        integer(kind=1) :: i
!        integer(kind=1) :: ii
!        integer(kind=4), dimension(DIMENSION) :: coordinate_binary_1D
!
!        do i = 1, DIMENSION
!            f_buffer = f(i)
!            coordinate_binary_1D(i) = lib_octree_hf_get_coordinate_binary_number_1D_float(f_buffer)
!        end do
!
!        ! make out of three binary coordinates one binary coordinate
!        !
!        ! Example
!        ! ----
!        !   x1: 0.100   => 0.5   (base 10)
!        !   x2: 0.010   => 0.25  (base 10)
!        !   x3: 0.001   => 0.125 (base 10)
!        !
!        !   x1: 0.1  |0  |0
!        !   x2: 0. 0 | 1 | 0
!        !   x3: 0.  0|  0|  1
!        !  ------------------
!        !  x3D: 0.100|010|001
!        !
!        ! Bit number example
!        ! ----
!        !   integer(kind=1): 1000 0100
!        !        bit number |7..4 3..0|
!        !
!
!        coordinate_binary_3D = 0 !ishft(coordinate_binary_3D, NUMBER_OF_BITS_COORDINATE_3D) ! set every bit to 0
!        do i = 1, DIMENSION
!            do ii = 0, NUMBER_OF_BITS_COORDINATE_1D - 1  ! bit operations: index starts at 0
!                if (btest(coordinate_binary_1D(i), NUMBER_OF_BITS_COORDINATE_1D - 1 - ii)) then
!                    coordinate_binary_3D = ibset(coordinate_binary_3D, NUMBER_OF_BITS_COORDINATE_3D - ii*DIMENSION - i)
!                else
!                    coordinate_binary_3D= ibclr(coordinate_binary_3D, NUMBER_OF_BITS_COORDINATE_3D - ii*DIMENSION - i)
!                end if
!            end do
!        end do
!
!
!    end function lib_octree_hf_get_coordinate_binary_number_3D_float
!
!    function lib_octree_hf_get_coordinate_binary_number_3D_double(f) result (coordinate_binary_3D)
!        implicit none
!        ! Calculates the binary coordinate of a normalised floating point vector (0..1, 0..1, 0..1).
!        !
!        !   x = (0.N_1 N_2 N_3 ... N_j...)_(2_d), N_j=(b_(1j) b_(2j)...b_(dj) )_2 , j=1, 2, ..., N_j=0, ..., 2_d−1. (77)
!        !
!        ! Reference: Data_Structures_Optimal_Choice_of_Parameters_and_C
!        !
!        ! Arguments
!        ! ----
!        !   f: double precision, dimension(3)
!        !       normalised single precision floating point number (0.0..1.0, 0.0..1.0, 0.0..1.0)
!        !
!        ! Returns
!        ! ----
!        !   the binary representation of the floating point vector (only the decimal place).
!        !
!        !   coordinate_binary_D: 16 bytes
!        !
!
!        ! parameters
!        integer(kind=1), parameter :: DIMENSION = 3         ! do not change, function is specialised in three-dimensional data points
!        integer(kind=1), parameter :: NUMBER_OF_BYTES_COORDINATE_3D = 16  ! space for 3 (=DIMENSION) integer of kind 8
!        integer(kind=1), parameter :: NUMBER_OF_BITS_COORDINATE_1D = 8 * NUMBER_OF_BITS_PER_BYTE
!        integer(kind=2), parameter :: NUMBER_OF_BITS_COORDINATE_3D = NUMBER_OF_BYTES_COORDINATE_3D * NUMBER_OF_BITS_PER_BYTE
!
!        ! dummy arguments
!        double precision, dimension(DIMENSION), intent(in) :: f
!        integer(kind=NUMBER_OF_BYTES_COORDINATE_3D) :: coordinate_binary_3D
!
!        ! auxiliary variables
!        double precision :: f_buffer
!        integer(kind=1) :: i
!        integer(kind=1) :: ii
!        integer(kind=8), dimension(DIMENSION) :: coordinate_binary_1D
!        integer(kind=1) :: bit_number_3D
!        integer(kind=1) :: bit_number_1D
!
!        do i = 1, DIMENSION
!            f_buffer = f(i)
!            coordinate_binary_1D(i) = lib_octree_hf_get_coordinate_binary_number_1D_double(f_buffer)
!        end do
!
!        ! make out of three binary coordinates one binary coordinate
!        !
!        ! Example
!        ! ----
!        !   x1: 0.100   => 0.5   (base 10)
!        !   x2: 0.010   => 0.25  (base 10)
!        !   x3: 0.001   => 0.125 (base 10)
!        !
!        !   x1: 0.1  |0  |0
!        !   x2: 0. 0 | 1 | 0
!        !   x3: 0.  0|  0|  1
!        !  ------------------
!        !  x3D: 0.100|010|001
!        !
!        ! Bit number example
!        ! ----
!        !   integer(kind=1): 1000 0100
!        !        bit number |7..4 3..0|
!
!        !
!        coordinate_binary_3D = 0 !ishft(coordinate_binary_3D, NUMBER_OF_BITS_COORDINATE_3D) ! set every bit to 0
!                bit_number_3D = NUMBER_OF_BITS_COORDINATE_3D - ii*DIMENSION - i
!        do i = 1, DIMENSION
!            do ii = 0, NUMBER_OF_BITS_COORDINATE_1D - 1  ! bit operations: index starts at 0
!                if (bit_number_3D >= 0) then
!                    bit_number_1D = NUMBER_OF_BITS_COORDINATE_1D - 1 - ii
!                    if (btest(coordinate_binary_1D(i), bit_number_1D)) then
!                        coordinate_binary_3D = ibset(coordinate_binary_3D, bit_number_3D)
!                    else
!                        coordinate_binary_3D= ibclr(coordinate_binary_3D, bit_number_3D)
!                    end if
!                end if
!            end do
!        end do
!
!
!    end function lib_octree_hf_get_coordinate_binary_number_3D_double
    
    ! Calculates the binary coordinate of a normalised floating point vector (0..1, 0..1, 0..1).
    ! For each element separately.
    !
    ! Arguments
    ! ----
    !   f: double precision, dimension(OCTREE_DIMENSIONS)
    !       normalised single precision floating point number (0.0..1.0, 0.0..1.0, 0.0..1.0)
    !
    ! Returns
    ! ----
    !   the binary representation of the floating point vector (only the decimal place).
    !
    !   coordinate_binary_D: vector<int(kind=8)>
    !
    function lib_octree_hf_get_coordinate_binary_number_xD(f) result (coordinate_binary_xD)
        implicit none

        ! dummy arguments
#if(_SPATIAL_POINT_IS_DOUBLE_ == 1)
        double precision, dimension(OCTREE_DIMENSIONS), intent(in) :: f
#elif(_SPATIAL_POINT_IS_DOUBLE_ == 0)
        real, dimension(OCTREE_DIMENSIONS), intent(in) :: f
#endif
        integer(kind=COORDINATE_BINARY_BYTES), dimension(OCTREE_DIMENSIONS) :: coordinate_binary_xD

        ! auxiliary variables
#if(_SPATIAL_POINT_IS_DOUBLE_ == 1)
        double precision :: f_buffer
#elif(_SPATIAL_POINT_IS_DOUBLE_ == 0)
        real :: f_buffer
#endif
        integer(kind=1) :: i

        do i = 1, OCTREE_DIMENSIONS
            f_buffer = f(i)
#if(_SPATIAL_POINT_IS_DOUBLE_ == 1)
            coordinate_binary_xD(i) = lib_octree_hf_get_coordinate_binary_number_1D_double(f_buffer)
#elif(_SPATIAL_POINT_IS_DOUBLE_ == 0)
            coordinate_binary_xD(i) = lib_octree_hf_get_coordinate_binary_number_1D_float(f_buffer)
#endif
        end do

    end function lib_octree_hf_get_coordinate_binary_number_xD

    ! Calculates the universal index of the parent box.
    !
    !   Parent(n) = [n/(2^d)]     (57)
    !
    ! Reference: Data_Structures_Optimal_Choice_of_Parameters_and_C
    !
    ! case: one-dimensional tree
    ! l=1
    !  |     0     |     1     |
    !  -------------------------
    ! l=2             ^
    !                 |
    !  |  0  |  1  |  2  |  3  |
    !  -------------------------
    !
    ! example:
    !   Box(n,l) = (2,2)
    !   Parent(2,2) = (1,1)
    !
    !
    ! Arguments
    ! ----
    !   uindex: type(lib_octree_universal_index)
    !       universal index of a box
    !
    ! Returns
    ! ----
    !   the universal index of the parent box.
    !
    !   parent_uindex: type(lib_octree_universal_index)
    !
    function lib_octree_hf_get_parent(uindex) result (parent_uindex)
        implicit none

        ! dummy arguments
        type(lib_octree_universal_index), intent (in) :: uindex
        type(lib_octree_universal_index) :: parent_uindex

        parent_uindex%n = uindex%n/(2**OCTREE_DIMENSIONS)
        parent_uindex%l = uindex%l

    end function

    ! Calculates the universal index of all children's boxes
    !
    !   ChildrenAll(n) = {2^d * n + j}, j=0, ..., 2^d − 1   (58)
    !
    ! Reference: Data_Structures_Optimal_Choice_of_Parameters_and_C
    !
    ! example: one-dimensional tree
    ! ----
    !    l=1
    !     |     0     |     1     |
    !     -------------------------
    !                      / \
    !    l=2              v   v
    !     |  0  |  1  |  2  |  3  |
    !     -------------------------
    !
    !      Box(n,l) = (1,1)
    !      ChildrenAll(1,1) = {(2,2), (3,2)}
    !
    ! Arguments
    ! ----
    !   n: integer(kind=4)
    !       universal index of a box
    !
    ! Returns
    ! ----
    !   the universal indexes of all children boxes of the given box.
    !
    !   children_n: Integer(kind=4), dimension(2^d)
    function lib_octree_hf_get_children_all(n) result (children_n)
        implicit none

        ! dummy arguments
        integer(kind=OCTREE_INTEGER_KIND), intent (in) :: n
        integer(kind=OCTREE_INTEGER_KIND), dimension(2**OCTREE_DIMENSIONS) :: children_n

        ! auxiliary variables
        integer(kind=1) :: j

        do j = 0, 2**OCTREE_DIMENSIONS - 1
            children_n(j+1) = 2**OCTREE_DIMENSIONS * n + j
        end  do

    end function

    ! Calculates the centre of a box(n,l).
    !
    !   x_m(n, l) = 2^(−l) ( n + 2^(−1))       (87)
    !
    ! Reference: Data_Structures_Optimal_Choice_of_Parameters_and_C
    !
    ! Arguments
    ! ----
    !   n: integer
    !       number of the node, with n ranging form 0 to 2^(3*l) in a three-dimensional space
    !   l: integer
    !       number of the level
    ! Returns
    ! ----
    !   x_m: [float, double]
    !       counting system undependet value
    !
    function lib_octree_hf_get_centre_of_box(n,l) result (point)
        implicit none

        ! dummy arguments
        integer(kind=OCTREE_INTEGER_KIND), intent (in) :: n
        integer(kind=1), intent (in) :: l
        type(lib_octree_spatial_point) :: point

        ! auxiliary variables
#if _SPATIAL_POINT_IS_DOUBLE_ == 1
        double precision :: x_m
#elif _SPATIAL_POINT_IS_DOUBLE_ == 0
        real :: x_m
#endif

!        x_m = 2d+0**(-l) * (n + 0.5d+0)
        x_m = 2**(-l) * (n + 0.5)

    end function

    ! Calculates the universal index of all k-th neigbour's boxes
    !
    !   N_(min)^(Neighbours)(d) = 2^d − 1                (19)
    !   N_(max)^(Neighbours)(d) = 3^d − 1                (19)
    !
    !   NeighborAll^(k)(n, l) = {(n−k, l), (n+k, l)}     (74)
    !
    ! Reference: Data_Structures_Optimal_Choice_of_Parameters_and_C
    !
    ! example: one-dimensional tree
    ! ----
    !    l=1
    !     |     0     |     1     |
    !     -------------------------
    !    l=2
    !     |  0  |  1  |  2  |  3  |
    !     -------------------------
    !             -k <-- * --> +k
    !
    !      Box(n,l) = (2,2)
    !      NeighbourAll(1,2,2) = {(1,2), (3,2)}
    !
    ! Arguments
    ! ----
    !   k: integer(kind=OCTREE_INTEGER_KIND)
    !       k-th neighbour
    !   n: integer(kind=OCTREE_INTEGER_KIND)
    !       universal index of a box
    !   l: integer(kind=OCTREE_INTEGER_KIND)
    !       number of the level
    !
    ! Returns
    ! ----
    !   the universal indexes of all neigbour boxes of the given box.
    !
    !   neighbour_all: Integer(kind=OCTREE_INTEGER_KIND), dimension(2^d)
    function lib_octree_hf_get_neighbour_all_1D(k,n,l) result (neighbour_all)
        implicit none

        ! dummy arguments
        integer(kind=1), intent (in) :: k
        integer(kind=OCTREE_INTEGER_KIND), intent (in) :: n
        integer(kind=1), intent (in) :: l
        !integer(kind=OCTREE_INTEGER_KIND), dimension(3**OCTREE_DIMENSIONS -1) :: neighbour_all
        integer(kind=OCTREE_INTEGER_KIND), dimension(2) :: neighbour_all

        ! auxiliary variables
        integer(kind=OCTREE_INTEGER_KIND), parameter :: ignore_entry = -1
        integer(kind=OCTREE_INTEGER_KIND), parameter :: lower_boundary = 0
        integer(kind=OCTREE_INTEGER_KIND) :: upper_boundary
        integer(kind=OCTREE_INTEGER_KIND) :: buffer_n

        upper_boundary = 2**l - 1

        ! calculate left neighbour
        buffer_n = n-k
        if (buffer_n < lower_boundary) then
            buffer_n = ignore_entry
        end if

        neighbour_all(1) = buffer_n

        ! calculate right neigbour
        buffer_n = n+k
        if (buffer_n > upper_boundary) then
            buffer_n = ignore_entry
        end if

        neighbour_all(2) = buffer_n

    end function

    ! Returns the look-up table (LUT) for interleaved bits. If necessary, the LUT is recalculated.
    !
    ! *Hint*
    !   This function has only global dependencies.
    !
    ! Dependences
    ! ----
    !   _FMM_DIMENSION_
    !       number of dimensions
    !   _INTERLEAVE_BITS_INTEGER_KIND_
    !       number of bytes of the integer
    !
    ! Returns
    ! ----
    !   rv: integer, dimension(:,:,:[,,:])
    !       the look-up tabel
    function lib_octree_hf_get_interleave_bits_lut() result(rv)
        implicit none
        ! parameter
        integer(kind=1), parameter :: integer_kind = _INTERLEAVE_BITS_INTEGER_KIND_
        integer(kind=integer_kind), parameter :: integer_range_high = huge(integer_range_high)
        integer(kind=integer_kind), parameter :: integer_range_low = -integer_range_high-1

#if (_FMM_DIMENSION_ == 2)
        ! parameter
        character (len = *), parameter :: file_lut="pre_calc/bit_interleaving_LUT_2d.dat"

        ! dummy
        integer(kind=integer_kind), dimension (integer_range_low:integer_range_high, &
                                               integer_range_low:integer_range_high, &
                                               1:2) :: rv

        ! check if LUT is already calculated
        ! if yes: load
        ! if not: calculate
        if (file_exists(file_lut)) then
            OPEN(UNIT=14, FILE=file_lut, ACTION="read", STATUS="old", &
                 FORM='unformatted')
            READ(14) rv
            CLOSE(UNIT=14)
        else
            rv = lib_octree_hf_creat_interleave_bits_lut()
            OPEN(UNIT=13, FILE=file_lut, ACTION="write", STATUS="replace", &
                 FORM="unformatted")
            WRITE(13) rv
            CLOSE(UNIT=13)
        end if
#elif (_FMM_DIMENSION_ == 3)
        ! parameter
        character (len = *), parameter :: file_lut="pre_calc/bit_interleaving_LUT_3d.dat"

        ! dummy
        integer(kind=integer_kind), dimension (integer_range_low:integer_range_high, &
                                           integer_range_low:integer_range_high, &
                                           integer_range_low:integer_range_high, &
                                           1:3) :: rv

        ! check if LUT is already calculated
        ! if yes: load
        ! if not: calculate
        if (file_exists(file_lut)) then
            OPEN(UNIT=14, FILE=file_lut, ACTION="read", STATUS="old", &
                 FORM='unformatted')
            READ(14) rv
            CLOSE(UNIT=14)
        else
            rv = lib_octree_hf_creat_interleave_bits_lut()
            OPEN(UNIT=13, FILE=file_lut, ACTION="write", STATUS="replace", &
                 FORM="unformatted")
            WRITE(13) rv
            CLOSE(UNIT=13)
        end if
#endif

    end function lib_octree_hf_get_interleave_bits_lut

    ! Creates the look-up table (LUT) for the interleaved bits.
    !
    ! *Hint*
    !   This function has only global dependencies.
    !
    ! Dependences
    ! ----
    !   _FMM_DIMENSION_
    !       number of dimensions
    !   _INTERLEAVE_BITS_INTEGER_KIND_
    !       number of bytes of the integer
    !
    ! Returns
    ! ----
    !   rv: integer, dimension(:,:,:[,,:])
    !       the look-up tabel
    !
    function lib_octree_hf_creat_interleave_bits_lut() result(rv)
        implicit none
        ! parameter
        integer(kind=1), parameter :: integer_kind = _INTERLEAVE_BITS_INTEGER_KIND_
        integer(kind=integer_kind), parameter :: integer_range_high = huge(integer_range_high)
        integer(kind=integer_kind), parameter :: integer_range_low = -integer_range_high-1

#if (_FMM_DIMENSION_ == 2)
        ! dummy
        integer(kind=integer_kind), dimension (integer_range_low:integer_range_high, &
                                               integer_range_low:integer_range_high, &
                                               1:2) :: rv

        ! auxiliary
        integer(kind=integer_kind*2) :: i
        integer(kind=integer_kind*2) :: ii

        integer(kind=integer_kind), dimension(2) :: buffer
        integer(kind=1) :: p
        integer(kind=1) :: p_old

        p_old = 0
        do i = integer_range_low, integer_range_high
            do ii = integer_range_low, integer_range_high
                buffer(1) = int(i,1)
                buffer(2) = int(ii,1)
                buffer = lib_octree_hf_interleave_bits(buffer)
                rv(i,ii, 1) = buffer(1)
                rv(i,ii, 2) = buffer(2)
            end do
            p = int(100.0*(i-integer_range_low)/(integer_range_high-integer_range_low), 1)
            if (int(p/10, 1) .ne. int(p_old/10, 1)) then
                print *, "Interleave bits: create LUT: ", p, "%"
            end if
            p_old = p
        end do

#elif (_FMM_DIMENSION_ == 3)
        ! dummy
        integer(kind=integer_kind), dimension (integer_range_low:integer_range_high, &
                                               integer_range_low:integer_range_high, &
                                               integer_range_low:integer_range_high, &
                                               1:3) :: rv

        ! auxiliary
        integer(kind=integer_kind*2) :: i
        integer(kind=integer_kind*2) :: ii
        integer(kind=integer_kind*2) :: iii

        integer(kind=integer_kind), dimension(3) :: buffer
        integer(kind=1) :: p
        integer(kind=1) :: p_old

        p_old = 0
        do i = integer_range_low, integer_range_high
            do ii = integer_range_low, integer_range_high
                do iii = integer_range_low, integer_range_high
                    buffer(1) = int(i,1)
                    buffer(2) = int(ii,1)
                    buffer(3) = int(iii,1)
                    buffer = lib_octree_hf_interleave_bits(buffer)
                    rv(i,ii, iii,1) = buffer(1)
                    rv(i,ii, iii,2) = buffer(2)
                    rv(i,ii, iii,3) = buffer(3)
                end do
            end do
            p = int(100.0*(i-integer_range_low)/(integer_range_high-integer_range_low), 1)
            if (int(p/10, 1) .ne. int(p_old/10, 1)) then
                print *, "Interleave bits: create LUT: ", p, "%"
            end if
            p_old = p
        end do
#endif
    end function lib_octree_hf_creat_interleave_bits_lut

    ! Calculates the bit interleaving of x-dimensional integers.
    ! The kind of the integer is defined with _INTERLEAVE_BITS_INTEGER_KIND_.
    ! The dafault value is 1 (1 byte);
    !
    ! Argument
    ! ----
    !   x: integer, x-dimensional
    !
    !
    ! Example
    ! ----
    !   x1: 0.100   => 0.5   (base 10)
    !   x2: 0.010   => 0.25  (base 10)
    !
    !   x1: 0. 1| 0| 0
    !   x2: 0.0 |1 |0
    !  ------------------
    !  x2D: 0.01|01|00
    function lib_octree_hf_interleave_bits(x) result(rv)
        implicit none
        ! parameter
        integer(kind=1), parameter :: x_kind = _INTERLEAVE_BITS_INTEGER_KIND_

        ! dummy
        integer(kind=x_kind), dimension(:), intent(in) :: x
        integer(kind=x_kind), dimension(size(x)) :: rv

        ! auxiliary
        integer(kind=1) :: i
        integer(kind=1) :: ii
        integer(kind=1) :: x_dimension
        integer(kind=1) :: bit_number
        integer(kind=1) :: x_element

        x_dimension = int(size(x), 1)

        do ii = 1, x_dimension
            rv(ii) = 0
        end do

        do ii = 1, x_dimension
            do i = 0, x_kind * NUMBER_OF_BITS_PER_BYTE - 1                                      ! e.g.: 16 bit =  2 byte * 8 bit / byte
                bit_number = i*x_dimension+ii-int(1, 1)                                         ! calculates the "global" bit number; two element example (1 byte / element): bit_number=9; element 2: |15 ... 8| element 1: |7 ... 0|
                x_element = int(bit_number / (x_kind * NUMBER_OF_BITS_PER_BYTE), 1) + int(1,1)  ! calculates the element; bit_number=9 => element=2
                bit_number = bit_number - (x_element-int(1,1))*x_kind*NUMBER_OF_BITS_PER_BYTE   ! calculates the "local" bit_number; bit_number=9, element=2 => bit_number=1
                if (btest(x(ii), i)) then       ! first bit number = 0
                    rv(x_element) = ibset(rv(x_element), bit_number)
                end if
            end do
        end do
    end function lib_octree_hf_interleave_bits

    ! Calculates the interleaved bits form a x-dimensional integer.
    !
    ! Arguments
    ! ----
    !   x: integer, dimension(x)
    !       integers to interleave
    !
    ! Returns
    ! ----
    !   rv: integer, dimension(:,:,:[,,:]) // 2- or 3-dimensional case
    !       interleaved integers
    !
    ! Example
    ! ----
    !   Interleaving of two one-byte integers
    !
    !   x1: 0 0 0 0| 0 0 1 0   => 2 (base 10)
    !   x2:  0 0 0 0 |0 0 0 0    => 0 (base 10)
    !      -------------------
    !       00000000|00001000   => 0|4 (base 10)
    !
    !   rv1 = 4
    !   rv2 = 0
    !
    function lib_octree_hf_interleave_bits_use_lut(x) result(rv)
        implicit none
        ! parameter
        integer(kind=1), parameter :: x_kind = _INTERLEAVE_BITS_INTEGER_KIND_

        integer(kind=x_kind), parameter :: integer_range_high = huge(integer_range_high)
        integer(kind=x_kind), parameter :: integer_range_low = -integer_range_high-1

        ! dummy
        integer(kind=x_kind), dimension(OCTREE_DIMENSIONS), intent(in) :: x
        integer(kind=x_kind), dimension(OCTREE_DIMENSIONS) :: rv

#if (_FMM_DIMENSION_ == 2)
        ! allocate memory
        if (.NOT. lib_octree_interleave_bits_lut_initialised) then
            allocate( lib_octree_interleave_bits_lut(integer_range_low:integer_range_high, &
                                                     integer_range_low:integer_range_high, &
                                                     1:2) )
            lib_octree_interleave_bits_lut = lib_octree_hf_creat_interleave_bits_lut()

            lib_octree_interleave_bits_lut_initialised = .true.
        end if
        rv(1) = lib_octree_interleave_bits_lut(x(2), x(1), 1)
        rv(2) = lib_octree_interleave_bits_lut(x(2), x(1), 2)
#elif (_FMM_DIMENSION_ == 3)
        ! allocate memory
        if (.NOT. lib_octree_interleave_bits_lut_initialised) then
            allocate( lib_octree_interleave_bits_lut(integer_range_low:integer_range_high, &
                                                 integer_range_low:integer_range_high, &
                                                 integer_range_low:integer_range_high, &
                                                 1:3) )
            lib_octree_interleave_bits_lut = lib_octree_hf_creat_interleave_bits_lut()

            lib_octree_interleave_bits_lut_initialised = .true.
        end if
        rv(1) = lib_octree_interleave_bits_lut(x(3), x(2), x(1), 1)
        rv(2) = lib_octree_interleave_bits_lut(x(3), x(2), x(1), 2)
        rv(3) = lib_octree_interleave_bits_lut(x(3), x(2), x(1), 3)
#else
        rv = lib_octree_hf_interleave_bits(x)
#endif
    end function lib_octree_hf_interleave_bits_use_lut

    ! Returns the look-up table (LUT) for deinterleaved bits. If necessary, the LUT is recalculated.
    !
    ! *Hint*
    !   This function has only global dependencies.
    !
    ! Dependences
    ! ----
    !   _FMM_DIMENSION_
    !       number of dimensions
    !   _INTERLEAVE_BITS_INTEGER_KIND_
    !       number of bytes of the integer
    !
    ! Returns
    ! ----
    !   rv: integer, dimension(:,:,:[,,:])
    !       the look-up tabel
    function lib_octree_hf_get_deinterleave_bits_lut() result(rv)
        implicit none
        ! parameter
        integer(kind=1), parameter :: integer_kind = _INTERLEAVE_BITS_INTEGER_KIND_
        integer(kind=integer_kind), parameter :: integer_range_high = huge(integer_range_high)
        integer(kind=integer_kind), parameter :: integer_range_low = -integer_range_high-1

#if (_FMM_DIMENSION_ == 2)
        ! parameter
        character (len = *), parameter :: file_lut="pre_calc/bit_deinterleaving_LUT_2d.dat"

        ! dummy
        integer(kind=integer_kind), dimension (integer_range_low:integer_range_high, &
                                               integer_range_low:integer_range_high, &
                                               1:2) :: rv

        ! check if LUT is already calculated
        ! if yes: load
        ! if not: calculate
        if (file_exists(file_lut)) then
            OPEN(UNIT=14, FILE=file_lut, ACTION="read", STATUS="old", &
                 FORM='unformatted')
            READ(14) rv
            CLOSE(UNIT=14)
        else
            rv = lib_octree_hf_creat_deinterleave_bits_lut()
            OPEN(UNIT=13, FILE=file_lut, ACTION="write", STATUS="replace", &
                 FORM="unformatted")
            WRITE(13) rv
            CLOSE(UNIT=13)
        end if
#elif (_FMM_DIMENSION_ == 3)
        ! parameter
        character (len = *), parameter :: file_lut="pre_calc/bit_deinterleaving_LUT_3d.dat"

        ! dummy
        integer(kind=integer_kind), dimension (integer_range_low:integer_range_high, &
                                           integer_range_low:integer_range_high, &
                                           integer_range_low:integer_range_high, &
                                           1:3) :: rv

        ! check if LUT is already calculated
        ! if yes: load
        ! if not: calculate
        if (file_exists(file_lut)) then
            OPEN(UNIT=14, FILE=file_lut, ACTION="read", STATUS="old", &
                 FORM='unformatted')
            READ(14) rv
            CLOSE(UNIT=14)
        else
            rv = lib_octree_hf_creat_deinterleave_bits_lut()
            OPEN(UNIT=13, FILE=file_lut, ACTION="write", STATUS="replace", &
                 FORM="unformatted")
            WRITE(13) rv
            CLOSE(UNIT=13)
        end if
#endif

    end function lib_octree_hf_get_deinterleave_bits_lut

    ! Creates the look-up table (LUT) for the deinterleaved bits.
    !
    ! *Hint*
    !   This function has only global dependencies.
    !
    ! Dependences
    ! ----
    !   _FMM_DIMENSION_
    !       number of dimensions
    !   _INTERLEAVE_BITS_INTEGER_KIND_
    !       number of bytes of the integer
    !
    ! Returns
    ! ----
    !   rv: integer, dimension(:,:,:[,,:])
    !       the look-up tabel
    !
    function lib_octree_hf_creat_deinterleave_bits_lut() result(rv)
        implicit none
        ! parameter
        integer(kind=1), parameter :: integer_kind = _INTERLEAVE_BITS_INTEGER_KIND_
        integer(kind=integer_kind), parameter :: integer_range_high = huge(integer_range_high)
        integer(kind=integer_kind), parameter :: integer_range_low = -integer_range_high-1

#if (_FMM_DIMENSION_ == 2)
        ! dummy
        integer(kind=integer_kind), dimension (integer_range_low:integer_range_high, &
                                               integer_range_low:integer_range_high, &
                                               1:2) :: rv

        ! auxiliary
        integer(kind=integer_kind*2) :: i
        integer(kind=integer_kind*2) :: ii

        integer(kind=integer_kind), dimension(2) :: buffer
        integer(kind=1) :: p
        integer(kind=1) :: p_old

        p_old = 0
        do i = integer_range_low, integer_range_high
            do ii = integer_range_low, integer_range_high
                buffer(1) = int(i,1)
                buffer(2) = int(ii,1)
                buffer = lib_octree_hf_deinterleave_bits(buffer)
                rv(i,ii, 1) = buffer(1)
                rv(i,ii, 2) = buffer(2)
            end do
            p = int(100.0*(i-integer_range_low)/(integer_range_high-integer_range_low), 1)
            if (int(p/10, 1) .ne. int(p_old/10, 1)) then
                print *, "Deinterleave bits: create LUT: ", p, "%"
            end if
            p_old = p
        end do

#elif (_FMM_DIMENSION_ == 3)
        ! dummy
        integer(kind=integer_kind), dimension (integer_range_low:integer_range_high, &
                                               integer_range_low:integer_range_high, &
                                               integer_range_low:integer_range_high, &
                                               1:3) :: rv

        ! auxiliary
        integer(kind=integer_kind*2) :: i
        integer(kind=integer_kind*2) :: ii
        integer(kind=integer_kind*2) :: iii

        integer(kind=integer_kind), dimension(2) :: buffer
        integer(kind=1) :: p
        integer(kind=1) :: p_old

        p_old = 0
        do i = integer_range_low, integer_range_high
            do ii = integer_range_low, integer_range_high
                do iii = integer_range_low, integer_range_high
                    buffer(1) = int(i,1)
                    buffer(2) = int(ii,1)
                    buffer(3) = int(iii,1)
                    buffer = lib_octree_hf_deinterleave_bits(buffer)
                    rv(i,ii, iii,1) = buffer(1)
                    rv(i,ii, iii,2) = buffer(2)
                    rv(i,ii, iii,3) = buffer(3)
                end do
            end do
            p = int(100.0*(i-integer_range_low)/(integer_range_high-integer_range_low), 1)
            if (int(p/10, 1) .ne. int(p_old/10, 1)) then
                print *, "Deinterleave bits: create LUT: ", p, "%"
            end if
            p_old = p
        end do
#endif
    end function lib_octree_hf_creat_deinterleave_bits_lut

    ! deinterleavs the
    ! The kind of the integer is defined with _INTERLEAVE_BITS_INTEGER_KIND_.
    ! The dafault value is 1 (1 byte);
    !
    ! Argument
    ! ----
    !   x:
    !       counting system undependet value
    !
    !
    ! Example
    ! ----
    !  x2D: 0.01|01|00
    !  ------------------
    !   x1: 0. 1| 0| 0
    !   x2: 0.0 |1 |0
    !
    !   x1: 0.100   => 0.5   (base 10)
    !   x2: 0.010   => 0.25  (base 10)
    !
    function lib_octree_hf_deinterleave_bits(x) result(rv)
        implicit none
        ! parameter
        integer(kind=1), parameter :: x_kind = _INTERLEAVE_BITS_INTEGER_KIND_

        ! dummy
        integer(kind=x_kind), dimension(:), intent(in) :: x
        integer(kind=x_kind), dimension(size(x)) :: rv

        ! auxiliary
        integer(kind=1) :: i
        integer(kind=1) :: ii
        integer(kind=1) :: x_dimension
        integer(kind=1) :: bit_number
        integer(kind=1) :: target_element

        x_dimension = int(size(x), 1)

        do ii = 1, x_dimension
            rv(ii) = 0
        end do

        do ii = 1, x_dimension
            do i = 0, x_kind * NUMBER_OF_BITS_PER_BYTE - 1                                      ! e.g.: 16 bit =  2 byte * 8 bit / byte
                ! global bit number = (ii-1)*x_kind * NUMBER_OF_BITS_PER_BYTE + i
                ! target_element = mod(gloabel bit number + 1, OCTREE_DIMENSIONS)
                ! target local bit number = int(global bit number / OCTREE_DIMENSIONS)
                target_element = int((ii-1)*x_kind * NUMBER_OF_BITS_PER_BYTE + i, 1)
                bit_number = int(target_element / OCTREE_DIMENSIONS, 1)
                target_element = int(mod(target_element, OCTREE_DIMENSIONS) + 1, 1)
                if (btest(x(ii), i)) then       ! first bit number = 0
                    rv(target_element) = ibset(rv(target_element), bit_number)
                end if
            end do
        end do
    end function lib_octree_hf_deinterleave_bits

    ! deinterleavs the
    ! The kind of the integer is defined with _INTERLEAVE_BITS_INTEGER_KIND_.
    ! The dafault value is 1 (1 byte);
    !
    ! Argument
    ! ----
    !   x:
    !       counting system undependet value
    !
    !
    ! Example
    ! ----
    !  x2D: 0.01|01|00
    !  ------------------
    !   x1: 0. 1| 0| 0
    !   x2: 0.0 |1 |0
    !
    !   x1: 0.100   => 0.5   (base 10)
    !   x2: 0.010   => 0.25  (base 10)
    !
    function lib_octree_hf_deinterleave_bits_use_lut(x) result(rv)
        implicit none
        ! parameter
        integer(kind=1), parameter :: x_kind = _INTERLEAVE_BITS_INTEGER_KIND_

        integer(kind=x_kind), parameter :: integer_range_high = huge(integer_range_high)
        integer(kind=x_kind), parameter :: integer_range_low = -integer_range_high-1

        ! dummy
        integer(kind=x_kind), dimension(OCTREE_DIMENSIONS), intent(in) :: x
        integer(kind=x_kind), dimension(OCTREE_DIMENSIONS) :: rv

#if (_FMM_DIMENSION_ == 2)
        ! allocate memory
        if (.NOT. lib_octree_deinterleave_bits_lut_initialised) then
            allocate( lib_octree_deinterleave_bits_lut(integer_range_low:integer_range_high, &
                                                       integer_range_low:integer_range_high, &
                                                       1:2) )
            lib_octree_deinterleave_bits_lut = lib_octree_hf_creat_deinterleave_bits_lut()

            lib_octree_deinterleave_bits_lut_initialised = .true.
        end if
        rv(1) = lib_octree_deinterleave_bits_lut(x(2), x(1), 1)
        rv(2) = lib_octree_deinterleave_bits_lut(x(2), x(1), 2)
#elif (_FMM_DIMENSION_ == 3)
        ! allocate memory
        if (.NOT. lib_octree_deinterleave_bits_lut_initialised) then
            allocate( lib_octree_deinterleave_bits_lut(integer_range_low:integer_range_high, &
                                                   integer_range_low:integer_range_high, &
                                                   integer_range_low:integer_range_high, &
                                                   1:3) )
            lib_octree_deinterleave_bits_lut = lib_octree_hf_creat_deinterleave_bits_lut()

            lib_octree_deinterleave_bits_lut_initialised = .true.
        end if
        rv(1) = lib_octree_deinterleave_bits_lut(x(3), x(2), x(1), 1)
        rv(2) = lib_octree_deinterleave_bits_lut(x(3), x(2), x(1), 2)
        rv(3) = lib_octree_deinterleave_bits_lut(x(3), x(2), x(1), 3)
#else
        rv = lib_octree_hf_deinterleave_bits(x)
#endif
    end function lib_octree_hf_deinterleave_bits_use_lut

! ----------------- test functions -----------------
    subroutine lib_octree_hf_test_functions()
        implicit none

        integer :: error_counter

        error_counter = 0

        if (.not. test_lib_octree_hf_get_universal_index()) then
            error_counter = error_counter + 1
        end if
        if (.not. test_lib_octree_hf_get_parent()) then
            error_counter = error_counter + 1
        end if
        if (.not. test_lib_octree_hf_get_children_all()) then
            error_counter = error_counter + 1
        end if
        if (.not. test_lib_octree_hf_get_centre_of_box()) then
            error_counter = error_counter + 1
        end if
        if (.not. test_lib_octree_hf_get_coordinate_binary_number_xD()) then
            error_counter = error_counter + 1
        end if
        if (.not. test_lib_octree_hf_interleave_bits()) then
            error_counter = error_counter + 1
        end if
        if (.not. test_lib_octree_hf_interleave_bits_2()) then
            error_counter = error_counter + 1
        end if
        if (.not. test_lib_octree_hf_deinterleave_bits()) then
            error_counter = error_counter + 1
        end if
        if (.not. test_lib_octree_hf_deinterleave_bits_2()) then
            error_counter = error_counter + 1
        end if


        if (error_counter == 0) then
            print *, "All tests: OK"
        else
            print *, error_counter,"test(s) FAILED"
        end if

        contains

        function test_lib_octree_hf_get_universal_index() result(rv)
            implicit none

            ! dummy
            logical :: rv

            type(lib_octree_spatial_point) :: point
            type(lib_octree_universal_index) :: universal_index

            integer(kind=1) :: l
            type(lib_octree_universal_index) :: universal_index_ground_trouth

            l = 1
            universal_index_ground_trouth%l = l

            point%x(1) = 0.75
            point%x(2) = 0.5

            universal_index_ground_trouth%n_per_dimension(1) = 1
            universal_index_ground_trouth%n_per_dimension(2) = 1
#if (_FMM_DIMENSION_ == 2)
            universal_index_ground_trouth%n = 3
#elif (_FMM_DIMENSION_ == 3)
            point%x(3) = 2.0**(-9.0) + 2.0**(-8)

            universal_index_ground_trouth%n = 6
            universal_index_ground_trouth%n_per_dimension(3) = 0
#else
            print *, "test_lib_octree_hf_get_parent: Dimension not defines: ", _FMM_DIMENSION_
#endif

            universal_index = lib_octree_hf_get_universal_index(point, l)

            if (universal_index%n == universal_index_ground_trouth%n) then
                rv = .true.
                print *, "test_lib_octree_hf_get_universal_index: ", "OK"
            else
                rv = .false.
                print *, "test_lib_octree_hf_get_universal_index: ", "FAILED"
            end if
        end function test_lib_octree_hf_get_universal_index

        function test_lib_octree_hf_get_parent() result (rv)
            implicit none

            ! dummy
            logical :: rv

            type(lib_octree_universal_index) :: universal_index
            type(lib_octree_universal_index) :: universal_index_parent

            integer(kind=1) :: l
            type(lib_octree_universal_index) :: universal_index_parent_ground_trouth

            l = 1
#if (_FMM_DIMENSION_ == 2)
            universal_index%n = 2
            universal_index_parent_ground_trouth%n = 0
#elif (_FMM_DIMENSION_ == 3)
            universal_index%n = 6
            universal_index_parent_ground_trouth%n = 0
#else
            print *, "test_lib_octree_hf_get_parent: Dimension not defines: ", _FMM_DIMENSION_
#endif

            universal_index_parent = lib_octree_hf_get_parent(universal_index)

            if (universal_index_parent%n == universal_index_parent_ground_trouth%n) then
                print *, "test_lib_octree_hf_get_parent: ", "ok"
                rv = .true.
            else
                print *, "test_lib_octree_hf_get_parent: ", "FAILED"
                rv = .false.
            end if

        end function test_lib_octree_hf_get_parent

        function test_lib_octree_hf_get_children_all() result (rv)
            implicit none

            ! dummy
            logical :: rv

            integer(kind=OCTREE_INTEGER_KIND) :: n
            integer(kind=OCTREE_INTEGER_KIND), dimension(2**OCTREE_DIMENSIONS) :: children_n
            integer(kind=OCTREE_INTEGER_KIND), dimension(2**OCTREE_DIMENSIONS) :: children_n_ground_truth

            n = 1
#if (_FMM_DIMENSION_ == 2)
            children_n_ground_truth(1) = 4
            children_n_ground_truth(2) = 5
            children_n_ground_truth(3) = 6
            children_n_ground_truth(4) = 7
#elif (_FMM_DIMENSION_ == 3)
            children_n_ground_truth(1) = 8
            children_n_ground_truth(2) = 9
            children_n_ground_truth(3) = 10
            children_n_ground_truth(4) = 11
            children_n_ground_truth(5) = 12
            children_n_ground_truth(6) = 13
            children_n_ground_truth(7) = 14
            children_n_ground_truth(8) = 15
#else
            print *, "test_lib_octree_hf_get_children_all: Dimension not defines: ", _FMM_DIMENSION_
#endif
            children_n = lib_octree_hf_get_children_all(n)

            if (sum(children_n) == sum(children_n_ground_truth)) then
                print *, "test_lib_octree_hf_get_children_all: ", "ok"
                rv = .true.
            else
                print *, "test_lib_octree_hf_get_children_all: ", "FAILED"
                rv = .false.
            end if

        end function test_lib_octree_hf_get_children_all

        function test_lib_octree_hf_get_centre_of_box() result (rv)
            implicit none

            ! dummy
            logical :: rv

            integer(kind=OCTREE_INTEGER_KIND) :: n
            integer(kind=1) :: l
            type(lib_octree_spatial_point) :: point

            type(lib_octree_spatial_point) :: point_ground_trouth

            integer :: i

            l=1
            n=1
#if (_FMM_DIMENSION_ == 2)
            point_ground_trouth%x(1) = 0.25
            point_ground_trouth%x(2) = 0.75
#elif (_FMM_DIMENSION_ == 3)
            point_ground_trouth%x(1) = 0.25
            point_ground_trouth%x(2) = 0.25
            point_ground_trouth%x(3) = 0.75
#else
            print *, "test_lib_octree_hf_get_centre_of_box: Dimension not defines: ", _FMM_DIMENSION_
#endif
            point = lib_octree_hf_get_centre_of_box(n,l)

            rv = .true.
            do i=1,OCTREE_DIMENSIONS
                if (point%x(i) == point_ground_trouth%x(i)) then
                    print *, "test_lib_octree_hf_get_centre_of_box (dim: ", i,"): ", "ok"
                else
                    print *, "test_lib_octree_hf_get_centre_of_box (dim: ", i,"): ", "FAILED"
                    rv = .false.
                end if
            end do

        end function test_lib_octree_hf_get_centre_of_box

        function test_lib_octree_hf_get_coordinate_binary_number_xD() result (rv)
            implicit none

            ! dummy
            logical :: rv

#if(_SPATIAL_POINT_IS_DOUBLE_ == 1)
            double precision, dimension(OCTREE_DIMENSIONS) :: f
#elif(_SPATIAL_POINT_IS_DOUBLE_ == 0)
            real, dimension(OCTREE_DIMENSIONS) :: f
#endif
            integer(kind=COORDINATE_BINARY_BYTES), dimension(OCTREE_DIMENSIONS) :: coordinate_binary_xD
            integer(kind=COORDINATE_BINARY_BYTES), dimension(OCTREE_DIMENSIONS) :: coordinate_binary_xD_ground_trouth

            integer :: i

#if (_FMM_DIMENSION_ == 2)
            f(1) = 0.25
            f(2) = 0.375
            if (COORDINATE_BINARY_BYTES == 4) then
                coordinate_binary_xD_ground_trouth(1) = 1073741824      ! |0100 0000|0000 0000| ... |0000 0000| byte 3-0
                coordinate_binary_xD_ground_trouth(2) = 1610612736      ! |0110 0000|0000 0000| ... |0000 0000| byte 3-0
            else if (COORDINATE_BINARY_BYTES == 8) then
                coordinate_binary_xD_ground_trouth(1) = 2**62           ! |0100 0000|0000 0000| ... |0000 0000| byte 7-0
                coordinate_binary_xD_ground_trouth(2) = 2**62 + 2**61   ! |0110 0000|0000 0000| ... |0000 0000| byte 7-0
            end if
#elif (_FMM_DIMENSION_ == 3)
            f(1) = 0.25
            f(2) = 0.375
            f(3) = 0.4375
            if (COORDINATE_BINARY_BYTES == 4) then
                coordinate_binary_xD_ground_trouth(1) = 1073741824          ! |0100 0000|0000 0000| ... |0000 0000| byte 3-0
                coordinate_binary_xD_ground_trouth(2) = 1610612736          ! |0110 0000|0000 0000| ... |0000 0000| byte 3-0
                coordinate_binary_xD_ground_trouth(3) = 1879048192          ! |0111 0000|0000 0000| ... |0000 0000| byte 3-0
            else if (COORDINATE_BINARY_BYTES == 8) then
                coordinate_binary_xD_ground_trouth(1) = 2**62               ! |0100 0000|0000 0000| ... |0000 0000| byte 7-0
                coordinate_binary_xD_ground_trouth(2) = 2**62+2**61         ! |0110 0000|0000 0000| ... |0000 0000| byte 7-0
                coordinate_binary_xD_ground_trouth(3) = 2**62+2**61+2**60   ! |0111 0000|0000 0000| ... |0000 0000| byte 7-0
            end if
#else
            print *, "test_lib_octree_hf_get_coordinate_binary_number_xD: Dimension not defines: ", _FMM_DIMENSION_
#endif

            coordinate_binary_xD = lib_octree_hf_get_coordinate_binary_number_xD(f)

            rv = .true.
            do i=1,OCTREE_DIMENSIONS
                if (coordinate_binary_xD(i) == coordinate_binary_xD_ground_trouth(i)) then
                    print *, "test_lib_octree_hf_get_coordinate_binary_number_xD (dim",i,"): ", "ok"
                else
                    print *, "test_lib_octree_hf_get_coordinate_binary_number_xD (dim",i,"): ", "FAILED"
                    rv = .false.
                end if
            end do

        end function test_lib_octree_hf_get_coordinate_binary_number_xD

        function test_lib_octree_hf_interleave_bits() result (rv)
            implicit none

            ! dummy
            logical :: rv

            ! parameter
            integer(kind=1), parameter :: x_kind = _INTERLEAVE_BITS_INTEGER_KIND_

            integer(kind=x_kind), dimension(OCTREE_DIMENSIONS) :: x
            integer(kind=x_kind), dimension(size(x)) :: interleaved_bits

            integer(kind=x_kind), dimension(size(x)) :: interleaved_bits_ground_trouth

            integer :: i

            x(1) = 4  ! 0100
            x(2) = 1  ! 0001
#if (_FMM_DIMENSION_ == 2)
            interleaved_bits_ground_trouth(1) = 2**1 + 2**4 !           |0001 0010|
            interleaved_bits_ground_trouth(2) = 0           ! |0000 0000|
#elif (_FMM_DIMENSION_ == 3)
            x(3) = 2  ! 0010
            interleaved_bits_ground_trouth(1) = 2**1 + 2**5 + 2**6 !                     |0110 0010|
            interleaved_bits_ground_trouth(2) = 0                  !           |0000 0000|
            interleaved_bits_ground_trouth(3) = 0                  ! |0000 0000|
#else
            print *, "test_lib_octree_hf_interleave_bits: Dimension not defines: ", _FMM_DIMENSION_
#endif

            interleaved_bits = lib_octree_hf_interleave_bits(x)

            rv = .true.
            do i=1, OCTREE_DIMENSIONS
                if (interleaved_bits(i) == interleaved_bits_ground_trouth(i)) then
                    print *, "test_lib_octree_hf_interleave_bits (dim",i,"): ", "ok"
                else
                    print *, "test_lib_octree_hf_interleave_bits (dim",i,"): ", "FAILED"
                    rv = .false.
                end if
            end do

        end function test_lib_octree_hf_interleave_bits

        function test_lib_octree_hf_interleave_bits_2() result (rv)
            implicit none

            ! dummy
            logical :: rv

            ! parameter
            integer(kind=1), parameter :: x_kind = _INTERLEAVE_BITS_INTEGER_KIND_

            integer(kind=x_kind), dimension(OCTREE_DIMENSIONS) :: x
            integer(kind=x_kind), dimension(size(x)) :: interleaved_bits

            integer(kind=x_kind), dimension(size(x)) :: interleaved_bits_ground_trouth

            integer :: i

            x(1) = 20  ! 0001 0100
            x(2) = 65  ! 0100 0001
#if (_FMM_DIMENSION_ == 2)
            interleaved_bits_ground_trouth(1) = 2**1 + 2**4  !           |0001 0010|
            interleaved_bits_ground_trouth(2) = 2**0 + 2**5  ! |0010 0001|
#elif (_FMM_DIMENSION_ == 3)
            x(3) = 40   ! 0010 1000
            interleaved_bits_ground_trouth(1) = 2**1 + 2**6  !                     |0100 0010|
            interleaved_bits_ground_trouth(2) = 2**3 + 2**4  !           |0001 1000|
            interleaved_bits_ground_trouth(3) = 2**1 + 2**3  ! |0000 1010|
#else
            print *, "test_lib_octree_hf_interleave_bits_2: Dimension not defines: ", _FMM_DIMENSION_
#endif

            interleaved_bits = lib_octree_hf_interleave_bits(x)

            rv = .true.
            do i=1, OCTREE_DIMENSIONS
                if (interleaved_bits(i) == interleaved_bits_ground_trouth(i)) then
                    print *, "test_lib_octree_hf_interleave_bits_2 (dim",i,"): ", "ok"
                else
                    print *, "test_lib_octree_hf_interleave_bits_2 (dim",i,"): ", "FAILED"
                    rv = .false.
                end if
            end do

        end function test_lib_octree_hf_interleave_bits_2

        function test_lib_octree_hf_deinterleave_bits() result (rv)
            implicit none

            ! dummy
            logical :: rv

            ! parameter
            integer(kind=1), parameter :: x_kind = _INTERLEAVE_BITS_INTEGER_KIND_

            ! dummy
            integer(kind=x_kind), dimension(OCTREE_DIMENSIONS) :: x
            integer(kind=x_kind), dimension(size(x)) :: deinterleaved_bits

            integer(kind=x_kind), dimension(size(x)) :: deinterleaved_bits_ground_trouth

            integer :: i

            deinterleaved_bits_ground_trouth(1) = 4  ! 0100
            deinterleaved_bits_ground_trouth(2) = 1  ! 0001

#if (_FMM_DIMENSION_ == 2)
            x(1) = 2**1 + 2**4 !           |0001 0010|
            x(2) = 0           ! |0000 0000|
#elif (_FMM_DIMENSION_ == 3)
            deinterleaved_bits_ground_trouth(3) = 2  ! 0010
            x(1) = 2**1 + 2**5 + 2**6 !                     |0110 0010|
            x(2) = 0                  !           |0000 0000|
            x(3) = 0                  ! |0000 0000|
#else
            print *, "test_lib_octree_hf_deinterleave_bits: Dimension not defines: ", _FMM_DIMENSION_
#endif

            deinterleaved_bits = lib_octree_hf_deinterleave_bits(x)

            rv = .true.
            do i=1, OCTREE_DIMENSIONS
                if (deinterleaved_bits(i) == deinterleaved_bits_ground_trouth(i)) then
                    print *, "test_lib_octree_hf_deinterleave_bits (dim",i,"): ", "ok"
                else
                    print *, "test_lib_octree_hf_deinterleave_bits (dim",i,"): ", "FAILED"
                    rv = .false.
                end if
            end do
        end function test_lib_octree_hf_deinterleave_bits

        function test_lib_octree_hf_deinterleave_bits_2() result (rv)
            implicit none

            ! dummy
            logical :: rv

            ! parameter
            integer(kind=1), parameter :: x_kind = _INTERLEAVE_BITS_INTEGER_KIND_

            ! dummy
            integer(kind=x_kind), dimension(OCTREE_DIMENSIONS) :: x
            integer(kind=x_kind), dimension(size(x)) :: deinterleaved_bits

            integer(kind=x_kind), dimension(size(x)) :: deinterleaved_bits_ground_trouth

            integer :: i

            deinterleaved_bits_ground_trouth(1) = 2**2 + 2**6  ! 0100 0100
            deinterleaved_bits_ground_trouth(2) = 2**0 + 2**5  ! 0010 0001

#if (_FMM_DIMENSION_ == 2)
            x(1) = 2**1 + 2**4 !           |0001 0010|
            x(2) = 2**3 + 2**4 ! |0001 1000|
#elif (_FMM_DIMENSION_ == 3)
            deinterleaved_bits_ground_trouth(3) = 2**1 + 2**4  ! 0001 0010
            x(1) = 2**1 + 2**5 + 2**6 !                     |0110 0010|
            x(2) = 2**6               !           |0100 0000|
            x(3) = 2**0 + 2**2        ! |0000 0101|
#else
            print *, "test_lib_octree_hf_deinterleave_bits_2: Dimension not defines: ", _FMM_DIMENSION_
#endif

            deinterleaved_bits = lib_octree_hf_deinterleave_bits(x)

            rv = .true.
            do i=1, OCTREE_DIMENSIONS
                if (deinterleaved_bits(i) == deinterleaved_bits_ground_trouth(i)) then
                    print *, "test_lib_octree_hf_deinterleave_bits_2 (dim",i,"): ", "ok"
                else
                    print *, "test_lib_octree_hf_deinterleave_bits_2 (dim",i,"): ", "FAILED"
                    rv = .false.
                end if
            end do
        end function test_lib_octree_hf_deinterleave_bits_2

    end subroutine lib_octree_hf_test_functions

    subroutine lib_octree_hf_benchmark
        implicit none

        call benchmark_lib_octree_hf_interleave_bits_use_lut()
        call benchmark_lib_octree_hf_deinterleave_bits_use_lut()

        contains

        subroutine benchmark_lib_octree_hf_interleave_bits_use_lut()
            implicit none

            integer(kind=1), dimension(OCTREE_DIMENSIONS) :: x
            integer(kind=1), dimension(OCTREE_DIMENSIONS) :: buffer

            integer :: number_of_runs = 100000000
            integer :: i
            real :: start, finish

            x(1) = 2
            x(2) = 0
#if (_FMM_DIMENSION_ == 3)
            x(3) = 0
#endif
            print *, "benchmark_lib_octree_hf_interleave_bits_use_lut"
            call cpu_time(start)
            buffer = lib_octree_hf_interleave_bits_use_lut(x)
            call cpu_time(finish)
            print *, "Interleave + LUT Time = ", finish-start, " seconds."

            call cpu_time(start)
            do i=1, number_of_runs
                buffer = lib_octree_hf_interleave_bits_use_lut(x)
            end do
            call cpu_time(finish)
            print *, "Interleave + LUT Time (second run) = ", (finish-start)/number_of_runs, " seconds."

            call cpu_time(start)
            do i=1, number_of_runs
                buffer = lib_octree_hf_interleave_bits(x)
            end do
            call cpu_time(finish)
            print *, "Interleave Time = ", (finish-start)/number_of_runs, " seconds."
        end subroutine benchmark_lib_octree_hf_interleave_bits_use_lut

        subroutine benchmark_lib_octree_hf_deinterleave_bits_use_lut()
            implicit none

            integer(kind=1), dimension(OCTREE_DIMENSIONS) :: x
            integer(kind=1), dimension(OCTREE_DIMENSIONS) :: buffer

            integer :: number_of_runs = 100000000
            integer :: i
            real :: start, finish

            x(1) = 2
            x(2) = 0
#if (_FMM_DIMENSION_ == 3)
            x(3) = 0
#endif
            print *, "benchmark_lib_octree_hf_interleave_bits_use_lut"
            call cpu_time(start)
            buffer = lib_octree_hf_deinterleave_bits_use_lut(x)
            call cpu_time(finish)
            print *, "Deinterleave + LUT Time = ", finish-start, " seconds."

            call cpu_time(start)
            do i=1, number_of_runs
                buffer = lib_octree_hf_deinterleave_bits_use_lut(x)
            end do
            call cpu_time(finish)
            print *, "Deinterleave + LUT Time (second run) = ", (finish-start)/number_of_runs, " seconds."

            call cpu_time(start)
            do i=1, number_of_runs
                buffer = lib_octree_hf_deinterleave_bits(x)
            end do
            call cpu_time(finish)
            print *, "Deinterleave Time = ", (finish-start)/number_of_runs, " seconds."
        end subroutine benchmark_lib_octree_hf_deinterleave_bits_use_lut

    end subroutine

end module lib_octree_helper_functions
