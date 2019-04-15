module lib_hash_function
    implicit none
    private

    public :: hash_fnv1a
    public :: hash_fnv1a_8_byte
    public :: hashpp_kf
    public :: hash_kf
    public :: hashpp_kf_8_byte
    public :: hash_kf_8_byte
    public :: hashpp_kf_16_byte
    public :: hash_kf_16_byte

    ! test functions
    public :: lib_test_hash_function

    contains

    ! FNV 1a 32 bit
    !
    ! Reference: http://www.isthe.com/chongo/tech/comp/fnv/index.html#FNV-1a
    ! hash = offset_basis
    ! for each octet_of_data to be hashed
    !         hash = hash xor octet_of_data
    !         hash = hash * FNV_prime
    ! return hash
    !
    ! Source code source:
    !   http://www.isthe.com/chongo/tech/comp/fnv/fnv32.f
    !
    !    *#######################################################################
    !    * 32 bit Fowler/Noll/Vo FNV-1a hash code
    !    * Public Domain
    !    *#######################################################################
    !    * fixed 32 -> 8 bit conversion bug            2013/10/12 AJA
    !    * translated to FORTRAN                       2013/03/11 Andy Allinger
    !    *                                             andy_a@users.sourceforge.net
    !    *
    !    * Revision: 5.1                               2009/06/30 09:13:32
    !    *#######################################################################
    !    * Fowler/Noll/Vo hash
    !    *
    !    * The basis of this hash algorithm was taken from an idea sent
    !    * as reviewer comments to the IEEE POSIX P1003.2 committee by:
    !    *
    !    *      Phong Vo (http://www.research.att.com/info/kpv/)
    !    *      Glenn Fowler (http://www.research.att.com/~gsf/)
    !    *
    !    * In a subsequent ballot round:
    !    *
    !    *      Landon Curt Noll (http://www.isthe.com/chongo/)
    !    *
    !    * improved on their algorithm.  Some people tried this hash
    !    * and found that it worked rather well.  In an EMail message
    !    * to Landon, they named it the ``Fowler/Noll/Vo'' or FNV hash.
    !    *
    !    * FNV hashes are designed to be fast while maintaining a low
    !    * collision rate. The FNV speed allows one to quickly hash lots
    !    * of data while maintaining a reasonable collision rate.  See:
    !    *
    !    *      http://www.isthe.com/chongo/tech/comp/fnv/index.html
    !    *
    !    * for more details as well as other forms of the FNV hash.
    !    *
    !    *#######################################################################
    !    * Please do not copyright this code.  This code is in the public domain.
    !    *
    !    * LANDON CURT NOLL DISCLAIMS ALL WARRANTIES WITH REGARD TO THIS SOFTWARE,
    !    * INCLUDING ALL IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS. IN NO
    !    * EVENT SHALL LANDON CURT NOLL BE LIABLE FOR ANY SPECIAL, INDIRECT OR
    !    * CONSEQUENTIAL DAMAGES OR ANY DAMAGES WHATSOEVER RESULTING FROM LOSS OF
    !    * USE, DATA OR PROFITS, WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE OR
    !    * OTHER TORTIOUS ACTION, ARISING OUT OF OR IN CONNECTION WITH THE USE OR
    !    * PERFORMANCE OF THIS SOFTWARE.
    !    *
    !    * By:
    !    *   chongo <Landon Curt Noll> /\oo/\
    !    *      http://www.isthe.com/chongo/
    !    *
    !    * Share and Enjoy!  :-)
    !    *
    !    *#######################################################################
    !    * 32 bit FNV-1 and FNV-1a non-zero initial basis
    !    *
    !    * The FNV-1 initial basis is the FNV-0 hash of the following 32 octets:
    !    *
    !    *              chongo <Landon Curt Noll> /\../\
    !    *
    !    * NOTE: The \'s above are not back-slashing escape characters.
    !    * They are literal ASCII  backslash 0x5c characters.
    !    *
    !    * NOTE: The FNV-1a initial basis is the same value as FNV-1 by definition.
    !    *
    !    *#######################################################################
    !    * perform a 32 bit Fowler/Noll/Vo FNV-1a hash on a buffer
    !    *
    !    * input:
    !    *   BUFFER  - start of buffer to hash
    !    *   LENGTH  - length of buffer in octets [bytes]
    !    *   HASH    - previous hash value or [standard initial value]
    !    *
    !    * assumption:  INTEGER's are at least 32 bits
    !    *
    !    * NOTE:  To use the recommended 32 bit FNV-1a hash, use the initial value:
    !    *               HASH = X'811C9DC5'
    !    *        as the argument on the first call to FNV32
    !    *
    !    * !NOTE:  Only pass the previous hash value if you always hash the same
    !    *         data in the same order (or else you will get different answers!)
    !    *
    !    * returns:
    !    *   32 bit HASH as default INTEGER
    !    *
    !    *#######################################################################
          SUBROUTINE FNV32 (BUFFER, LENGTH, HASH)
           IMPLICIT NONE
           INTEGER LENGTH, HASH
           INTEGER*1 BUFFER
           DIMENSION BUFFER(LENGTH)

           INTEGER PRIME ; PARAMETER (PRIME = 16777619)
           INTEGER I, J, K
           INTEGER*1 B

    !    *#######################################################################
    !    *                begin
    !    *#######################################################################
    !    *          FNV-1a hash each octet in the buffer
           DO 90 J = 1, LENGTH
             B = BUFFER(J)
             K = 0
             DO 80 I = 0, 7           ! copy each bit from B to K
               IF (BTEST(B, I)) K = IBSET(K, I)
      80     CONTINUE ! next i

    !    *          xor the bottom with the current octet
             HASH = IEOR(HASH, K)

    !    *          multiply by the 32 bit FNV magic prime mod 2^32
             HASH = HASH * PRIME
             HASH = IAND(HASH, X'FFFFFFFF')      ! discard > 32 bits
      90   CONTINUE ! next j

          END !############## of file fnv32.f ##############################

    function hash_fnv1a(buffer) result(hash)
        ! dummy
        integer(kind=4), intent(in) :: buffer
        integer(kind=4) :: hash

        ! auxiliary
        integer(kind=4) :: buffer_buffer
        integer(kind=1), dimension(4) :: buffer_list

        equivalence (buffer_buffer, buffer_list)

        buffer_buffer = buffer

        hash = -2128831035
        call FNV32(buffer_list, size(buffer_list), hash)

    end function hash_fnv1a

    function hash_fnv1a_8_byte(buffer) result(hash)
        ! dummy
        integer(kind=8), intent(in) :: buffer
        integer(kind=4) :: hash

        ! auxiliary
        integer(kind=4) :: buffer_buffer
        integer(kind=1), dimension(8) :: buffer_list

        equivalence (buffer_buffer, buffer_list)

        buffer_buffer = buffer

        hash = -2128831035
        call FNV32(buffer_list, size(buffer_list), hash)

    end function hash_fnv1a_8_byte

!    ! ********************************************************* hash
!    !berechnet ndigit-Hash wert der Koordinaten Matrix a,b,max,n ,versuch i
!    ! Source: K. Frenner: fem_sparse.F90
!    integer(kind=8) function hash(a,b,max,i,idum)
!    implicit none
!    integer,intent(in)::a,b,max,i
!    integer(kind=8),intent(inout):: idum
!      integer i2
!      idum=int(a+int(b,kind=8)*int(max,kind=8),kind=8) !als start
!      do i2=1,i+2,1  !2=Offset zum einschwingen !
!        idum=int(modulo(real(16807.0D0*idum,kind=16),2147483647.0D0),kind=8)
!      end do
!      hash=int(real(idum,kind=8)/2147483647.0D0*real(max-1,kind=8),kind=8)
!    end function

!    ! ********************************************************* hash
!    !berechnet ndigit-Hash wert der Koordinaten Matrix a,b,max,n ,versuch i
!    ! Source: K. Frenner: fem_sparse.F90
!    integer(kind=8) function hash_kf(a,b,max,i)
!    implicit none
!    integer(kind=8),intent(in)::a,b,max
!    integer, intent(in)::i
!    integer(kind=8):: idum
!      integer i2
!      idum=int(a+int(b,kind=8)*int(max,kind=8),kind=8) !als start
!      do i2=1,i+2,1  !2=Offset zum einschwingen !
!        idum=int(modulo(real(16807.0D0*idum,kind=16),2147483647.0D0),kind=8)
!      end do
!      hash_kf=int(real(idum,kind=8)/2147483647.0D0*real(max-1,kind=8),kind=8)
!    end function
!
!
!    ! ********************************************************* hashpp
!    !beschleunigte Fkt
!    ! Source: K. Frenner: fem_sparse.F90
!    integer(kind=8) function hashpp_kf(max,idum)
!    implicit none
!    integer(kind=8),intent(inout):: idum
!    integer(kind=8),intent(in)::max
!      idum=int(modulo(real(16807.0D0*idum,kind=16),2147483647.0D0),kind=8)
!      hashpp_kf=int(real(idum,kind=8)/2147483647.0D0*real(max-1,kind=8),kind=8)
!    end function

    ! ********************************************************* hash
    !berechnet ndigit-Hash wert der Koordinaten Matrix a,b,max,n ,versuch i
    integer(kind=8) function hash_kf(a,b,max,i,idum) result(hash)
    implicit none
    integer(kind=8),intent(in)::a,b,max,i
    integer(kind=8),intent(inout):: idum
      integer i2
!      idum=int(a+int(b,kind=8)*int(max,kind=8),kind=8) !als start
      idum=int(a, 8)
      do i2=1,i+2,1  !2=Offset zum einschwingen !
        idum=int(modulo(real(16807.0D0*idum,kind=16),2147483647.0D0),kind=8)
      end do
      hash=int(real(idum,kind=8)/2147483647.0D0*real(max-1,kind=8),kind=8)
    end function

    ! ********************************************************* hashpp
    !beschleunigte Fkt
    integer(kind=8) function hashpp_kf(max,idum) result (hashpp)
    implicit none
    integer(kind=8),intent(inout):: idum
    integer(kind=8),intent(in)::max
      idum=int(modulo(real(16807.0D0*idum,kind=16),2147483647.0D0),kind=8)
      hashpp=int(real(idum,kind=8)/2147483647.0D0*real(max-1,kind=8),kind=8)
    end function

    ! ********************************************************* hash
    !berechnet ndigit-Hash wert der Koordinaten Matrix a,b,max,n ,versuch i
    integer(kind=8) function hash_kf_8_byte(a,b,max,i,idum) result(hash)
    implicit none
        integer(kind=8), intent(in) :: a
        integer(kind=4), intent(in)::b,max,i
        integer(kind=8), dimension(2), intent(inout):: idum
        integer(kind=4) i2

        ! auxiliary
        integer(kind=8) :: aa
        integer(kind=4), dimension(2) :: buffer_a
        real(kind=8) :: buffer

        equivalence(aa, buffer_a)

        aa = a

        idum(1)=int(buffer_a(1)+int(b,kind=8)*int(max,kind=8),kind=8) !als start
        idum(2)=int(buffer_a(2)+int(b,kind=8)*int(max,kind=8),kind=8) !als start
        do i2=1,i+2,1  !2=Offset zum einschwingen !
            idum(1)=int(modulo(real(16807.0D0*idum(1),kind=16),2147483647.0D0),kind=8)
            idum(2)=int(modulo(real(16807.0D0*idum(2),kind=16),2147483647.0D0),kind=8)
        end do
        buffer = real(idum(1),kind=8)/2147483647.0D0
        buffer = buffer * real(idum(2),kind=8)/2147483647.0D0

        hash=int(buffer*real(max-1,kind=8),kind=8)
    end function

    ! ********************************************************* hashpp
    !beschleunigte Fkt
    integer(kind=8) function hashpp_kf_8_byte(max,idum) result (hashpp)
    implicit none
        integer(kind=8), dimension(2), intent(inout):: idum
        integer(kind=4), intent(in)::max

        !auxiliary
        real(kind=8) :: buffer

        idum(1)=int(modulo(real(16807.0D0*idum(1),kind=16),2147483647.0D0),kind=8)
        idum(2)=int(modulo(real(16807.0D0*idum(2),kind=16),2147483647.0D0),kind=8)

        buffer = real(idum(1),kind=8)/2147483647.0D0
        buffer = buffer * real(idum(2),kind=8)/2147483647.0D0

        hashpp=int(buffer*real(max-1,kind=8),kind=8)
    end function

    ! ********************************************************* hash
    !berechnet ndigit-Hash wert der Koordinaten Matrix a,b,max,n ,versuch i
    integer(kind=8) function hash_kf_16_byte(a,b,max,i,idum) result(hash)
    implicit none
        integer(kind=16), intent(in) :: a
        integer(kind=8),intent(in)::b,max,i
        integer(kind=8), dimension(2), intent(inout):: idum
        integer(kind=8) i2

        ! auxiliary
        integer(kind=16) :: aa
        integer(kind=8), dimension(2) :: buffer_a
        real(kind=8) :: buffer

        equivalence(aa, buffer_a)

        aa = a

        idum(1)=int(buffer_a(1)+int(b,kind=8)*int(max,kind=8),kind=8) !als start
        idum(2)=int(buffer_a(2)+int(b,kind=8)*int(max,kind=8),kind=8) !als start
        do i2=1,i+2,1  !2=Offset zum einschwingen !
            idum(1)=int(modulo(real(16807.0D0*idum(1),kind=16),2147483647.0D0),kind=8)
            idum(2)=int(modulo(real(16807.0D0*idum(2),kind=16),2147483647.0D0),kind=8)
        end do
        buffer = real(idum(1),kind=8)/2147483647.0D0
        buffer = buffer * real(idum(2),kind=8)/2147483647.0D0
        hash=int(buffer*real(max-1,kind=8),kind=8)
    end function

    ! ********************************************************* hashpp
    !beschleunigte Fkt
    integer(kind=8) function hashpp_kf_16_byte(max,idum) result (hashpp)
    implicit none
        integer(kind=8), dimension(2),intent(inout):: idum
        integer(kind=8), intent(in)::max

        !auxiliary
        real(kind=8) :: buffer
        idum(1)=int(modulo(real(16807.0D0*idum(1),kind=16),2147483647.0D0),kind=8)
        idum(2)=int(modulo(real(16807.0D0*idum(2),kind=16),2147483647.0D0),kind=8)

        buffer = real(idum(1),kind=8)/2147483647.0D0
        buffer = buffer * real(idum(2),kind=8)/2147483647.0D0
        hashpp=int(buffer*real(max-1,kind=8),kind=8)
    end function


    ! ----------- test functions -----------
    subroutine lib_test_hash_function
        implicit none

!        if (.not. test_hash()) then
!            print *, "test_hash: error"
!        end if
!
!        if (.not. test_hashpp()) then
!            print *, "test_hashpp: error"
!        end if

        contains

!        function test_hash() result (rv)
!            implicit none
!            ! dummy
!            logical :: rv
!
!            integer ::i
!            integer(kind=8) :: a, start, max, rv_hash
!
!            a = 40
!            max = 10**7
!            i=2
!            start = 1
!
!            rv_hash = hash_kf(start, a, max, i)
!            rv_hash = hash_kf(start, a, max, i)
!            rv_hash = hash_kf(start, a, max, i)
!            rv_hash = hash_kf(start, rv_hash, max, i)
!            rv_hash = hash_kf(start, rv_hash, max, i)
!            rv_hash = hash_kf(start, rv_hash, max, i)
!
!            rv = .true.
!
!        end function

!        function test_hashpp() result (rv)
!            implicit none
!
!            ! dummy
!            logical :: rv
!
!            integer(kind=8) :: idum
!            integer(kind=8) :: max
!            integer(kind=8) :: hash
!
!            max = 2
!            idum = 1
!
!            hash = hashpp_kf(max, idum)
!            hash = hashpp_kf(max, idum)
!            hash = hashpp_kf(max, idum)
!            hash = hashpp_kf(max, idum)
!            hash = hashpp_kf(max, idum)
!            hash = hashpp_kf(max, idum)
!            hash = hashpp_kf(max, idum)
!            hash = hashpp_kf(max, idum)
!
!            rv = .true.
!
!        end function test_hashpp

    end subroutine lib_test_hash_function
end module lib_hash_function
