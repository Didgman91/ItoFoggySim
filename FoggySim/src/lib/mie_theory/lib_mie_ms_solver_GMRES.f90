module lib_mie_ms_solver_GMRES
    use libmath
    use lib_mie_type
    use lib_mie_type_functions
    use lib_mie_ms_solver_interface
    implicit none

    ! parameter

    integer, parameter :: GMRES_ORTHOGONALIZATION_SCHEME_MGS = 0
    integer, parameter :: GMRES_ORTHOGONALIZATION_SCHEME_IMGS = 1
    integer, parameter :: GMRES_ORTHOGONALIZATION_SCHEME_CGS = 2
    integer, parameter :: GMRES_ORTHOGONALIZATION_SCHEME_ICGS = 3

    type solver_gmres_parameter_type
        integer :: max_iterations
        integer :: restart
        logical :: use_initial_guess
        double precision :: convergence_tolerance
        integer :: orthogonalization_scheme
        logical :: use_recurence_formula_at_restart
        logical :: residual_calc_explicitly

        double precision :: backward_error
    end type solver_gmres_parameter_type

    contains

        ! Returns
        ! ----
        !   rv: type(solver_gmres_parameter_type)
        !       default values of the GMRES solver
        function lib_mie_ms_solver_gmres_get_parameter_std_values() result(rv)
            implicit none
            ! dummy
            type(solver_gmres_parameter_type) :: rv

            rv%max_iterations = 1000
            rv%restart = 10
            rv%use_initial_guess = .false.
            rv%convergence_tolerance = 1D-15
            rv%orthogonalization_scheme = GMRES_ORTHOGONALIZATION_SCHEME_MGS
            rv%use_recurence_formula_at_restart = .false.
            rv%residual_calc_explicitly = .true.

        end function lib_mie_ms_solver_gmres_get_parameter_std_values

        ! Argument
        ! ----
        !   parameter: type(solver_gmres_parameter_type)
        !       parameter of the GMRES solver
        !   simulation_parameter: type(lib_mie_simulation_parameter_type)
        !       dataset of the simulation
        !   save_solution: logical, optional (std: .true.)
        !       true: save solution x into simulation_parameter
        subroutine lib_mie_ms_solver_gmres_run_without_ml_fmm(gmres_parameter, simulation_parameter, save_solution)
            implicit none
            ! dummy
            type(solver_gmres_parameter_type), intent(inout) :: gmres_parameter
            type(lib_mie_simulation_parameter_type), intent(inout) :: simulation_parameter
            logical, intent(in), optional :: save_solution

            ! parameter
            integer, parameter :: matvec = 1
            integer, parameter :: precondLeft = 2
            integer, parameter :: precondRight = 3
            integer, parameter :: dotProd = 4

            double complex, parameter :: ZERO = dcmplx(0,0)
            double complex, parameter :: ONE = dcmplx(1,0)

            ! auxiliaray
            integer :: n
            integer :: m
            integer :: lda
            integer :: lwork
            integer :: ldstrt
            integer revcom, colx, coly, colz, nbscal
            integer, dimension(5) :: irc
            integer, dimension(8) :: icntl
            integer, dimension(3) :: info

            integer nout

            double complex, dimension(:), allocatable :: work
            double precision  cntl(5), rinfo(2)

            ! auxiliary
            integer, dimension(2) :: n_range
            integer :: no_of_spheres

            double complex, dimension(:), allocatable :: vector

            integer :: m_use_initial_guess
            logical :: m_save_solution
            integer :: m_residual_calc


            ldstrt = gmres_parameter%restart

            if (gmres_parameter%use_initial_guess) then
                m_use_initial_guess = 1
            else
                m_use_initial_guess = 0
            end if

            if (gmres_parameter%residual_calc_explicitly) then
                m_residual_calc = 1
            else
                m_residual_calc = 0
            end if

            m_save_solution = .true.
            if (present(save_solution)) m_save_solution = save_solution

            nout = 6

            ! calculate size(x) = size(b)
            n_range = simulation_parameter%spherical_harmonics%n_range
            no_of_spheres = size(simulation_parameter%sphere_list)

            lda = no_of_spheres * 2 * ( (1 + n_range(2))**2 - n_range(1) )

            ! Reference: A Set of GMRES Routines for Real and Complex Arithmetics on High Performance Computers,
            !            Valerie Frayss LucGiraud Serge Gratton Julien Langou, p. 9
            if (gmres_parameter%orthogonalization_scheme .eq. GMRES_ORTHOGONALIZATION_SCHEME_MGS &
                .or. gmres_parameter%orthogonalization_scheme .eq. GMRES_ORTHOGONALIZATION_SCHEME_IMGS) then
                if (gmres_parameter%residual_calc_explicitly) then
                    lwork = ldstrt**2 + ldstrt*(lda+5) + 5*lda + 2
                else
                    lwork = ldstrt**2 + ldstrt*(lda+5) + 6*lda + 2
                end if
            else
                if (gmres_parameter%residual_calc_explicitly) then
                    lwork = ldstrt**2 + ldstrt*(lda+5) + 5*lda + ldstrt + 1
                else
                    lwork = ldstrt**2 + ldstrt*(lda+5) + 6*lda + ldstrt + 1
                end if
            end if

            allocate(work(lwork))

            ! ----------------------------------------------------
            !  Initialize the control parameters to default value
            ! ----------------------------------------------------
            call init_zgmres(icntl,cntl)

            ! -----------------------
            !  Tune some parameters
            ! -----------------------

            ! Tolerance
            cntl(1) = gmres_parameter%convergence_tolerance
            ! Save the convergence history on standard output
            icntl(3) = 6
            ! Maximum number of iterations
            icntl(7) = gmres_parameter%max_iterations

            ! preconditioner location
            icntl(4) = 0
            ! orthogonalization scheme
            icntl(5)=gmres_parameter%orthogonalization_scheme
            ! initial guess
            icntl(6) = m_use_initial_guess
            ! residual calculation strategy at restart
            icntl(8) = m_residual_calc

            ! set initial guess x_0
            call lib_mie_ms_solver_get_vector_x(simulation_parameter, vector)
            work(1:lda) = vector

            ! Initialise the right hand side b
            call lib_mie_ms_solver_get_vector_b(simulation_parameter, vector)
            work(lda+1:2*lda) = vector

            n = lda
            m = ldstrt

            ! ---------------------------------------
            !  Reverse communication implementation
            ! ---------------------------------------
            do
                call drive_zgmres(n,n,m,lwork,work, &
                                  irc,icntl,cntl,info,rinfo)
                revcom = irc(1)
                colx   = irc(2)
                coly   = irc(3)
                colz   = irc(4)
                nbscal = irc(5)

                if (revcom.eq.matvec) then
                    ! perform the matrix vector product
                    !        work(colz) <-- A * work(colx)
!                    call zgemv('N',n,n,ONE,a,lda,work(colx),1, &
!                               ZERO,work(colz),1)
                    call lib_mie_ms_solver_calculate_vector_b(simulation_parameter, work(colx:colx+lda-1), vector)
                    work(colz:colz+lda-1) = vector

                    cycle

                else if (revcom.eq.precondLeft) then
                    ! perform the left preconditioning
                    !         work(colz) <-- M^{-1} * work(colx)
                    print *, "lib_mie_ms_solver_gmres_run_without_ml_fmm: ERROR"
                    print *, "  no preconditioning (left), but required"

                    exit

                else if (revcom.eq.precondRight) then
                    ! perform the right preconditioning
                    print *, "lib_mie_ms_solver_gmres_run_without_ml_fmm: ERROR"
                    print *, "  no preconditioning (right), but required"

                    exit

                else if (revcom.eq.dotProd) then
                    ! perform the scalar product
                    ! work(colz) <-- work(colx) work(coly)

                   call zgemv('C',n,nbscal,ONE,work(colx),n, &
                              work(coly),1,ZERO,work(colz),1)

                   cycle
                endif

                ! -----------------------------
                ! dump the solution on a file
                ! -----------------------------
                !
                if (icntl(5).eq.0) then
                  write(nout,*) 'Orthogonalisation : MGS'
                elseif (icntl(5).eq.1) then
                  write(nout,*) 'Orthogonalisation : IMGS'
                elseif (icntl(5).eq.2) then
                  write(nout,*) 'Orthogonalisation : CGS'
                elseif (icntl(5).eq.3) then
                  write(nout,*) 'Orthogonalisation : ICGS'
                endif
                write(nout,*) 'Restart : ', m
                write(nout,*) 'info(1) = ',info(1),'  info(2) = ',info(2)
                write(nout,*) 'rinfo(1) = ',rinfo(1),'  rinfo(2) = ',rinfo(2)
                write(nout,*) 'Optimal workspace = ', info(3)
                write(nout,*) ' **********************************************!  '
                write(nout,*)
                write(*,*) ' **********************************************!  '
                write(*,*)
                write(*,*) '    '

                gmres_parameter%backward_error = rinfo(2)

                exit

            end do

            if (m_save_solution) then
                call lib_mie_ms_solver_set_sphere_parameter_ab_nm(work(1:lda), simulation_parameter)
            end if

        end subroutine lib_mie_ms_solver_gmres_run_without_ml_fmm

        ! Argument
        ! ----
        !   parameter: type(solver_gmres_parameter_type)
        !       parameter of the GMRES solver
        !   simulation_parameter: type(lib_mie_simulation_parameter_type)
        !       dataset of the simulation
        !   save_solution: logical, optional (std: .true.)
        !       true: save solution x into simulation_parameter
        subroutine lib_mie_ms_solver_gmres_run_with_ml_fmm(gmres_parameter, simulation_parameter, save_solution)
            implicit none
            ! dummy
            type(solver_gmres_parameter_type), intent(inout) :: gmres_parameter
            type(lib_mie_simulation_parameter_type), intent(inout) :: simulation_parameter
            logical, intent(in), optional :: save_solution

            ! parameter
            integer, parameter :: matvec = 1
            integer, parameter :: precondLeft = 2
            integer, parameter :: precondRight = 3
            integer, parameter :: dotProd = 4

            double complex, parameter :: ZERO = dcmplx(0,0)
            double complex, parameter :: ONE = dcmplx(1,0)

            ! auxiliaray
            integer :: n
            integer :: m
            integer :: lda
            integer :: lwork
            integer :: ldstrt
            integer revcom, colx, coly, colz, nbscal
            integer, dimension(5) :: irc
            integer, dimension(8) :: icntl
            integer, dimension(3) :: info

            integer nout

            double complex, dimension(:), allocatable :: work
            double precision  cntl(5), rinfo(2)

            ! auxiliary
            integer :: i
            integer, dimension(2) :: n_range
            integer :: no_of_spheres

            double complex, dimension(:), allocatable :: vector

            integer, dimension(:), allocatable :: selection

            integer :: m_use_initial_guess
            logical :: m_save_solution
            integer :: m_residual_calc


            ldstrt = gmres_parameter%restart

            if (gmres_parameter%use_initial_guess) then
                m_use_initial_guess = 1
            else
                m_use_initial_guess = 0
            end if

            if (gmres_parameter%residual_calc_explicitly) then
                m_residual_calc = 1
            else
                m_residual_calc = 0
            end if

            m_save_solution = .true.
            if (present(save_solution)) m_save_solution = save_solution

            nout = 6

            ! calculate size(x) = size(b)
            n_range = simulation_parameter%spherical_harmonics%n_range
            no_of_spheres = size(simulation_parameter%sphere_list)

            lda = no_of_spheres * 2 * ( (1 + n_range(2))**2 - n_range(1) )

            ! Reference: A Set of GMRES Routines for Real and Complex Arithmetics on High Performance Computers,
            !            Valerie Frayss LucGiraud Serge Gratton Julien Langou, p. 9
            if (gmres_parameter%orthogonalization_scheme .eq. GMRES_ORTHOGONALIZATION_SCHEME_MGS &
                .or. gmres_parameter%orthogonalization_scheme .eq. GMRES_ORTHOGONALIZATION_SCHEME_IMGS) then
                if (gmres_parameter%residual_calc_explicitly) then
                    lwork = ldstrt**2 + ldstrt*(lda+5) + 5*lda + 2
                else
                    lwork = ldstrt**2 + ldstrt*(lda+5) + 6*lda + 2
                end if
            else
                if (gmres_parameter%residual_calc_explicitly) then
                    lwork = ldstrt**2 + ldstrt*(lda+5) + 5*lda + ldstrt + 1
                else
                    lwork = ldstrt**2 + ldstrt*(lda+5) + 6*lda + ldstrt + 1
                end if
            end if

            allocate(work(lwork))

            ! ----------------------------------------------------
            !  Initialize the control parameters to default value
            ! ----------------------------------------------------
            call init_zgmres(icntl,cntl)

            ! -----------------------
            !  Tune some parameters
            ! -----------------------

            ! Tolerance
            cntl(1) = gmres_parameter%convergence_tolerance
            ! Save the convergence history on standard output
            icntl(3) = 6
            ! Maximum number of iterations
            icntl(7) = gmres_parameter%max_iterations

            ! preconditioner location
            icntl(4) = 0
            ! orthogonalization scheme
            icntl(5)=gmres_parameter%orthogonalization_scheme
            ! initial guess
            icntl(6) = m_use_initial_guess
            ! residual calculation strategy at restart
            icntl(8) = m_residual_calc

            ! set initial guess x_0
            call lib_mie_ms_solver_get_vector_x(simulation_parameter, vector)
            work(1:lda) = vector

            ! Initialise the right hand side b
            call lib_mie_ms_solver_get_vector_b(simulation_parameter, vector)
            work(lda+1:2*lda) = vector

            n = lda
            m = ldstrt

            ! ---------------------------------------
            !  Reverse communication implementation
            ! ---------------------------------------
            do
                call drive_zgmres(n,n,m,lwork,work, &
                                  irc,icntl,cntl,info,rinfo)
                revcom = irc(1)
                colx   = irc(2)
                coly   = irc(3)
                colz   = irc(4)
                nbscal = irc(5)

                if (revcom.eq.matvec) then
                    ! perform the matrix vector product
                    !        work(colz) <-- A * work(colx)
!                    call zgemv('N',n,n,ONE,a,lda,work(colx),1, &
!                               ZERO,work(colz),1)
                    if (allocated(selection)) deallocate(selection)
                    allocate(selection(size(simulation_parameter%sphere_list)))

                    do i = 1, size(simulation_parameter%sphere_list)
                        selection(i) = i
                    end do

                    call lib_mie_ms_solver_calculate_vector_b(simulation_parameter, selection, &
                                                              work(colx:colx+lda-1), vector)
                    work(colz:colz+lda-1) = vector
                    cycle

                else if (revcom.eq.precondLeft) then
                    ! perform the left preconditioning
                    !         work(colz) <-- M^{-1} * work(colx)
                    print *, "lib_mie_ms_solver_gmres_run_without_ml_fmm: ERROR"
                    print *, "  no preconditioning (left), but required"

                    exit

                else if (revcom.eq.precondRight) then
                    ! perform the right preconditioning
                    print *, "lib_mie_ms_solver_gmres_run_without_ml_fmm: ERROR"
                    print *, "  no preconditioning (right), but required"

                    exit

                else if (revcom.eq.dotProd) then
                    ! perform the scalar product
                    ! work(colz) <-- work(colx) work(coly)

                   call zgemv('C',n,nbscal,ONE,work(colx),n, &
                              work(coly),1,ZERO,work(colz),1)

                   cycle
                endif

                ! -----------------------------
                ! dump the solution on a file
                ! -----------------------------
                !
                if (icntl(5).eq.0) then
                  write(nout,*) 'Orthogonalisation : MGS'
                elseif (icntl(5).eq.1) then
                  write(nout,*) 'Orthogonalisation : IMGS'
                elseif (icntl(5).eq.2) then
                  write(nout,*) 'Orthogonalisation : CGS'
                elseif (icntl(5).eq.3) then
                  write(nout,*) 'Orthogonalisation : ICGS'
                endif
                write(nout,*) 'Restart : ', m
                write(nout,*) 'info(1) = ',info(1),'  info(2) = ',info(2)
                write(nout,*) 'rinfo(1) = ',rinfo(1),'  rinfo(2) = ',rinfo(2)
                write(nout,*) 'Optimal workspace = ', info(3)
                write(nout,*) ' **********************************************!  '
                write(nout,*)
                write(*,*) ' **********************************************!  '
                write(*,*)
                write(*,*) '    '

                gmres_parameter%backward_error = rinfo(2)

                exit

            end do

            if (m_save_solution) then
                call lib_mie_ms_solver_set_sphere_parameter_ab_nm(work(1:lda), simulation_parameter)
            end if

        end subroutine lib_mie_ms_solver_gmres_run_with_ml_fmm

        subroutine GMRES_test()
            implicit none
            ! dummy

            ! parameter
            integer, parameter :: lda = 1000
            integer, parameter :: ldstrt = 60

            integer, parameter :: matvec = 1
            integer, parameter :: precondLeft = 2
            integer, parameter :: precondRight = 3
            integer, parameter :: dotProd = 4

            double complex, parameter :: ZERO = dcmplx(0,0)
            double complex, parameter :: ONE = dcmplx(1,0)


            ! auxiliary
            integer i, j, n, m
            integer :: lwork
            integer revcom, colx, coly, colz, nbscal
            integer irc(5), icntl(8), info(3)

            integer nout

            double complex, dimension(:,:), allocatable :: a
            double complex, dimension(:), allocatable :: work
            double precision  cntl(5), rinfo(2)

            ! ------------------------------------------------------------
            !  Generate the test matrix a and set the right-hand side
            !  in positions (n+1) to 2n of the array work.
            !  The right-hand side is chosen such that the exact solution
            !  is the vector of all ones.
            ! ------------------------------------------------------------

            write(*,*) '***********************************************'
            write(*,*) 'This code is an example of use of GMRES'
            write(*,*) 'in double precision complex arithmetic'
            write(*,*) 'Results are written standard output'
            write(*,*) '***********************************************'
            write(*,*)
!            write(*,*) 'Matrix size < ', lda
            n = 100
            write(*,*) 'Matrix size = ', n
            if (n.gt.lda) then
              write(*,*) 'You are asking for a too large matrix'
              return
            end if

            lwork = ldstrt**2 + ldstrt*(lda+5) + 5*lda + 1

            allocate(a(lda, lda))
            allocate(work(lwork))

            write(*,*) "test"

            do j = 1,n
              do i = 1,n
                a(i,j) = ZERO
              end do
            end do

            do i = 1,n
              a(i,i) = dcmplx(4.d0, 0.d0)
            end do
            do i = 1,n-1
              a(i,i+1) = dcmplx(-2.d0, 1.d0)
              a(i+1,i) = dcmplx(-2.d0, 1.d0)
            end do
            nout = 1
            ! -------------------------------
            !  Choose the restart parameter
            ! -------------------------------

!            write(*,*) 'Restart  <', ldstrt
            m = ldstrt
            write(*,*) 'Restart  =', m

            ! ----------------------------------------------------
            !  Initialize the control parameters to default value
            ! ----------------------------------------------------

            call init_zgmres(icntl,cntl)

            ! -----------------------
            !  Tune some parameters
            ! -----------------------

            ! Tolerance
            cntl(1) = 1.d-15
            ! Save the convergence history on standard output
            icntl(3) = 6
            ! Maximum number of iterations
            icntl(7) = 1000

            ! preconditioner location
            icntl(4) = 1
            ! orthogonalization scheme
            icntl(5)=0
            ! initial guess
            icntl(6) = 0
            ! residual calculation strategy at restart
            icntl(8) = 0
            ! Initialise the right hand side
            do j = 1,n
                work(j) = ONE
            end do
            call ZGEMV('N',n,n,ONE,A,lda,work(1),1,ZERO,work(n+1),1)
            do j = 1,n
                work(j) = ONE/2.0
            end do

            ! ---------------------------------------
            !  Reverse communication implementation
            ! ---------------------------------------

            do
                call drive_zgmres(n,n,m,lwork,work, &
                                  irc,icntl,cntl,info,rinfo)
                revcom = irc(1)
                colx   = irc(2)
                coly   = irc(3)
                colz   = irc(4)
                nbscal = irc(5)

                if (revcom.eq.matvec) then
                    ! perform the matrix vector product
                    !        work(colz) <-- A * work(colx)
                    call zgemv('N',n,n,ONE,a,lda,work(colx),1, &
                               ZERO,work(colz),1)
                    cycle

                else if (revcom.eq.precondLeft) then
                    ! perform the left preconditioning
                    !         work(colz) <-- M^{-1} * work(colx)
                    call zcopy(n,work(colx),1,work(colz),1)
                    call ztrsm('L','L','N','N',n,1,ONE,A,lda,work(colz),n)

                    cycle

                else if (revcom.eq.precondRight) then
                    ! perform the right preconditioning
                    call zcopy(n,work(colx),1,work(colz),1)
                    call ztrsm('L','U','N','N',n,1,ONE,A,lda,work(colz),n)

                    cycle

                else if (revcom.eq.dotProd) then
                    ! perform the scalar product
                    ! work(colz) <-- work(colx) work(coly)

                   call zgemv('C',n,nbscal,ONE,work(colx),n, &
                              work(coly),1,ZERO,work(colz),1)

                   cycle
                endif

                ! -----------------------------
                ! dump the solution on a file
                ! -----------------------------
                !
                if (icntl(5).eq.0) then
                  write(nout,*) 'Orthogonalisation : MGS'
                elseif (icntl(5).eq.1) then
                  write(nout,*) 'Orthogonalisation : IMGS'
                elseif (icntl(5).eq.2) then
                  write(nout,*) 'Orthogonalisation : CGS'
                elseif (icntl(5).eq.3) then
                  write(nout,*) 'Orthogonalisation : ICGS'
                endif
                write(nout,*) 'Restart : ', m
                write(nout,*) 'info(1) = ',info(1),'  info(2) = ',info(2)
                write(nout,*) 'rinfo(1) = ',rinfo(1),'  rinfo(2) = ',rinfo(2)
                write(nout,*) 'Optimal workspace = ', info(3)
                write(nout,*) ' **********************************************!  '
                write(nout,*)
                write(*,*) ' **********************************************!  '
                write(*,*)
                write(*,*) '    '

                exit

            end do
    end subroutine

end module lib_mie_ms_solver_GMRES
