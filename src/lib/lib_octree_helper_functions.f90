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

! integer kind of the bit interleaving process, value = [1,2,3,4,8]
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
#elif (_FMM_DIMENSION_ == 3)
    integer(kind=_INTERLEAVE_BITS_INTEGER_KIND_), dimension (:,:,:,:) &
                                                , allocatable :: lib_octree_interleave_bits_lut
    logical :: lib_octree_interleave_bits_lut_initialised = .false.
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
    public :: lib_octree_hf_interleave_bits_use_lut
    public :: lib_octree_hf_interleave_bits

contains

    ! Cleans up the memory
    subroutine lib_octree_hf_destructor()
        implicit none

#if (_FMM_DIMENSION_ == 2)
        deallocate (lib_octree_interleave_bits_lut)
#elif (_FMM_DIMENSION_ == 3)
        deallocate (lib_octree_interleave_bits_lut)
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

!    function lib_octree_hf_get_universal_index(point_x, l) result(n)
!        implicit none
!        ! Calculates the universal index *n* of a given normalised floating point *point_x*.
!        ! Related to the level *l*.
!        !
!        !   n = (2**d)**(l-1)*N_1 + (2**d)**(l-2)*N_2 + ... + (2**d)*N_l-1 + N_l      (49)
!        !
!        ! Reference: Data_Structures_Optimal_Choice_of_Parameters_and_C
!        !
!        !
!        ! Example:
!        !   Point |  x_10    |  x_2    |  n(l=2)
!        !   --------------------------------------------
!        !      1  |  0.125   | 0.001   |  2*0 + 1*0 = 0
!        !      2  |  0.3125  | 0.0101  |  2*0 + 1*1 = 1
!        !      3  |  0.375   | 0.011   |  2*0 + 1*1 = 1
!        !      4  |  0.625   | 0.101   |  2*1 + 1*0 = 2
!        !      5  |  0.825   | 0.111   |  2*1 + 1*1 = 3
!        !
!        !
!        ! Arguments
!        ! ----
!        !   point_x: double precision
!        !       normalised floating point number (0.0 .. 1.0)
!        !       HINT: datatype real is also possible, use *_1D_float() instead
!        !
!        !   l: integer(kind=1)
!        !       number of layers
!        !
!        ! Returns
!        ! ----
!        !   the universal index *n*.
!        !
!        !   n: integer(kind=4)
!        !
!
!        ! dummy arguments
!        type(lib_octree_spatial_point), intent (in) :: point_x
!        integer(kind=1), intent (in) :: l
!        integer(kind=OCTREE_INTEGER_KIND) :: n
!
!        ! auxiliary
!        integer :: i
!        integer(kind=1) :: buffer
!#if(_FMM_DIMENSION_ == 1)
!#if(_SPATIAL_POINT_IS_DOUBLE_ == 1)
!        integer, parameter :: COORDINATE_BINARY_NUMBER_OF_BYTES = 8
!
!#elif(_SPATIAL_POINT_IS_DOUBLE_ == 0)
!        integer, parameter :: COORDINATE_BINARY_NUMBER_OF_BYTES = 4
!#endif
!#elif(_FMM_DIMENSION_ == 3)
!#if(_SPATIAL_POINT_IS_DOUBLE_ == 1)
!        integer, parameter :: COORDINATE_BINARY_NUMBER_OF_BYTES = 16
!#elif(_SPATIAL_POINT_IS_DOUBLE_ == 0)
!        integer, parameter :: COORDINATE_BINARY_NUMBER_OF_BYTES = 16
!#endif
!#endif
!        integer(kind=COORDINATE_BINARY_NUMBER_OF_BYTES) :: coordinate_binary
!        integer(kind=COORDINATE_BINARY_NUMBER_OF_BYTES) :: n_buffer
!
!        !        coordinate_binary = get_coordinate_binary_number_1D_float(point_x)
!#if(_FMM_DIMENSION_ == 1)
!#if(_SPATIAL_POINT_IS_DOUBLE_ == 1)
!        coordinate_binary = lib_octree_hf_get_coordinate_binary_number_1D_double(point_x%x(1))
!#elif(_SPATIAL_POINT_IS_DOUBLE_ == 0)
!        coordinate_binary = lib_octree_hf_get_coordinate_binary_number_1D_float(point_x%x(1))
!#endif
!#elif(_FMM_DIMENSION_ == 3)
!#if(_SPATIAL_POINT_IS_DOUBLE_ == 1)
!        coordinate_binary = lib_octree_hf_get_coordinate_binary_number_3D_double(point_x%x)
!#elif(_SPATIAL_POINT_IS_DOUBLE_ == 0)
!        coordinate_binary = lib_octree_hf_get_coordinate_binary_number_3D_float(point_x%x)
!#endif
!#endif
!
!
!        !        ! use these lines instead to load all binary coordinates into the *point_n* variable.
!        !        do i = 1, COORDINATE_BINARY_NUMBER_OF_BYTES*8
!        !            if (btest(coordinate_binary, i-1)) then ! bit number starts at 0
!        !                point_n(COORDINATE_BINARY_NUMBER_OF_BYTES*8-i+1) = 1
!        !            else
!        !                point_n(COORDINATE_BINARY_NUMBER_OF_BYTES*8-i+1) = 0
!        !            end if
!        !        end do
!        !
!        !        n = 0
!        !        do i = 1, l
!        !            n = n +  2**(l-i) * point_n(i)
!        !        end do
!
!        i = l*OCTREE_DIMENSIONS
!        n_buffer = ibits(coordinate_binary, COORDINATE_BINARY_NUMBER_OF_BYTES*NUMBER_OF_BITS_PER_BYTE-i, i)
!
!        n = n_buffer
!
!!        n = 0
!!        do i = 1, l*OCTREE_DIMENSIONS
!!            if (btest(coordinate_binary, COORDINATE_BINARY_NUMBER_OF_BYTES*NUMBER_OF_BITS_PER_BYTE-i)) then ! bit number starts at 0
!!                buffer = 1
!!            else
!!                buffer = 0
!!            end if
!!
!!            n = n +  (2**OCTREE_DIMENSIONS)**(l-i) * buffer
!!        end do
!
!    end function lib_octree_hf_get_universal_index

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

#if (_FMM_DIMENSION_ == 2)
#elif (_FMM_DIMENSION_ == 3)
        ! Example
        ! ----
        !
        ! coordinate_binary: kind=4, dimension=3
        ! element   byte representation
        ! 1:       |---01---|---02---|---03---|---04---|
        ! 2:       |---01---|---02---|---03---|---04---|
        ! 3:       |---01---|---02---|---03---|---04---|
        !
        ! cb_buffer: kind=1, dimension=3*4=12
        ! element   byte representation
        ! 1-4:     |---01---|---02---|---03---|---04---|
        ! 5-8:     |---05---|---06---|---07---|---08---|
        ! 9-12:    |---09---|---10---|---11---|---12---|
        !              *-->
        !              interleave column wise
        !
        ! equivalence (coordinate_binary, cb_buffer)
        !
        ! Interleave bits
        ! -----
        !
        ! interleaved_bits: kind=1, dimension=3*4=12
        ! element   byte representation
        ! 1-4:     |01-05-09|01-05-09|01-05-09|02-06-10|   <-- xx-yy-zz describes the bytes which were interleaved
        ! 5-8:     |02-06-10|02-06-10|03-07-11|03-07-11|
        ! 9-12:    |03-07-11|04-08-12|04-08-12|04-08-12|
        !
        ! ib_buffer: kind=1, dimension=3
        ! element   byte representation
        ! 1-3:     |---01---|---02---|---03---|
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
        integer(kind=_INTERLEAVE_BITS_INTEGER_KIND_) &
            ,dimension(OCTREE_DIMENSIONS * COORDINATE_BINARY_BYTES/_INTERLEAVE_BITS_INTEGER_KIND_) &
            :: cb_buffer

        equivalence (coordinate_binary, cb_buffer)

        integer(kind=_INTERLEAVE_BITS_INTEGER_KIND_), dimension(OCTREE_DIMENSIONS) :: ib_buffer

        integer(kind=1) :: i
        integer(kind=1) :: ii

        ! calculate the x-dimensional binary coordinate
        coordinate_binary = lib_octree_hf_get_coordinate_binary_number_xD(point_x%x)

        ! interleave bits
        do i=1, COORDINATE_BINARY_BYTES/_INTERLEAVE_BITS_INTEGER_KIND_
            do ii=1, OCTREE_DIMENSIONS
                ib_buffer(ii) = cb_buffer((i-1) * COORDINATE_BINARY_BYTES/_INTERLEAVE_BITS_INTEGER_KIND_ + ii)
            end do
            ib_buffer = lib_octree_hf_interleave_bits_use_lut(ib_buffer)

            ii = int((i-1)*COORDINATE_BINARY_BYTES/_INTERLEAVE_BITS_INTEGER_KIND_ + 1, 1)
            interleaved_bits(ii:ii+OCTREE_DIMENSIONS+1) = ib_buffer(:)
        end do

        ! calculte the universal index n
!        do i=1, int(_INTERLEAVE_BITS_INTEGER_KIND_ * NUMBER_OF_BITS_PER_BYTE / l)+1  ! maximum number of bytes to respect by a given l
!
!        end do
#else

        coordinate_binary = lib_octree_hf_get_coordinate_binary_number_xD(point_x%x)

        do i = 1, OCTREE_DIMENSIONS
            uindex%n(i) = ibits(coordinate_binary(i), COORDINATE_BINARY_BYTES*NUMBER_OF_BITS_PER_BYTE-l, l)
        end do
#endif

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
    function lib_octree_hf_get_lut() result(rv)
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
            rv = lib_octree_hf_creat_lut()
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
            rv = lib_octree_hf_creat_lut()
            OPEN(UNIT=13, FILE=file_lut, ACTION="write", STATUS="replace", &
                 FORM="unformatted")
            WRITE(13) rv
            CLOSE(UNIT=13)
        end if
#endif

    end function lib_octree_hf_get_lut

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
    function lib_octree_hf_creat_lut() result(rv)
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

        do i = integer_range_low, integer_range_high
            do ii = integer_range_low, integer_range_high
                buffer(1) = int(i,1)
                buffer(2) = int(ii,1)
                buffer = lib_octree_hf_interleave_bits(buffer)
                rv(i,ii, 1) = buffer(1)
                rv(i,ii, 2) = buffer(2)
            end do
            write(1, *) '(f3.2)',100.0*(i-integer_range_low)/(integer_range_high-integer_range_low)
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
            print *, "create LUT: ", 100.0*(i-integer_range_low)/(integer_range_high-integer_range_low), "%"
        end do
#endif
    end function lib_octree_hf_creat_lut

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
    !   x1:  0 0 0 0| 0 0 1 0   => 2 (base 10)
    !   x2: 0 0 0 0 |0 0 0 0    => 0 (base 10)
    !      -------------------
    !       00000000|00000100   => 0|4 (base 10)
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
            lib_octree_interleave_bits_lut = lib_octree_hf_creat_lut()

            lib_octree_interleave_bits_lut_initialised = .true.
        end if
        rv(1) = lib_octree_interleave_bits_lut(x(1), x(2), 1)
        rv(2) = lib_octree_interleave_bits_lut(x(1), x(2), 2)
#elif (_FMM_DIMENSION_ == 3)
        ! allocate memory
        if (.NOT. lib_octree_interleave_bits_lut_initialised) then
            allocate( lib_octree_interleave_bits_lut(integer_range_low:integer_range_high, &
                                                 integer_range_low:integer_range_high, &
                                                 integer_range_low:integer_range_high, &
                                                 1:3) )
            lib_octree_interleave_bits_lut = lib_octree_hf_creat_lut()

            lib_octree_interleave_bits_lut_initialised = .true.
        end if
        rv(1) = lib_octree_interleave_bits_lut(x(1), x(2), x(3), 1)
        rv(2) = lib_octree_interleave_bits_lut(x(1), x(2), x(3), 2)
        rv(3) = lib_octree_interleave_bits_lut(x(1), x(2), x(3), 3)
#else
        rv = lib_octree_hf_interleave_bits(x)
#endif


    end function lib_octree_hf_interleave_bits_use_lut

end module lib_octree_helper_functions
