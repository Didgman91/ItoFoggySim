module lib_hash_function
    implicit none
    private

    public :: hash
    public :: hashpp

    ! test functions
    public :: lib_test_hash_function

    contains


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

    ! ********************************************************* hash
    !berechnet ndigit-Hash wert der Koordinaten Matrix a,b,max,n ,versuch i
    ! Source: K. Frenner: fem_sparse.F90
    integer(kind=8) function hash(a,b,max,i)
    implicit none
    integer(kind=8),intent(in)::a,b,max
    integer, intent(in)::i
    integer(kind=8):: idum
      integer i2
      idum=int(a+int(b,kind=8)*int(max,kind=8),kind=8) !als start
      do i2=1,i+2,1  !2=Offset zum einschwingen !
        idum=int(modulo(real(16807.0D0*idum,kind=16),2147483647.0D0),kind=8)
      end do
      hash=int(real(idum,kind=8)/2147483647.0D0*real(max-1,kind=8),kind=8)
    end function


    ! ********************************************************* hashpp
    !beschleunigte Fkt
    ! Source: K. Frenner: fem_sparse.F90
    integer(kind=8) function hashpp(max,idum)
    implicit none
    integer(kind=8),intent(inout):: idum
    integer(kind=8),intent(in)::max
      idum=int(modulo(real(16807.0D0*idum,kind=16),2147483647.0D0),kind=8)
      hashpp=int(real(idum,kind=8)/2147483647.0D0*real(max-1,kind=8),kind=8)
    end function

    ! ----------- test functions -----------
    subroutine lib_test_hash_function
        implicit none

        if (.not. test_hash()) then
            print *, "test_hash: error"
        end if

        if (.not. test_hashpp()) then
            print *, "test_hashpp: error"
        end if

        contains

        function test_hash() result (rv)
            implicit none
            ! dummy
            logical :: rv

            integer ::i
            integer(kind=8) :: a, start, max, rv_hash

            a = 40
            max = 10**7
            i=2
            start = 1

            rv_hash = hash(start, a, max, i)
            rv_hash = hash(start, a, max, i)
            rv_hash = hash(start, a, max, i)
            rv_hash = hash(start, rv_hash, max, i)
            rv_hash = hash(start, rv_hash, max, i)
            rv_hash = hash(start, rv_hash, max, i)

            rv = .true.

        end function

        function test_hashpp() result (rv)
            implicit none

            ! dummy
            logical :: rv

            integer(kind=8) :: idum
            integer(kind=8) :: max
            integer(kind=8) :: hash

            max = 2
            idum = 1

            hash = hashpp(max, idum)
            hash = hashpp(max, idum)
            hash = hashpp(max, idum)
            hash = hashpp(max, idum)
            hash = hashpp(max, idum)
            hash = hashpp(max, idum)
            hash = hashpp(max, idum)
            hash = hashpp(max, idum)

            rv = .true.

        end function test_hashpp

    end subroutine lib_test_hash_function
end module lib_hash_function
