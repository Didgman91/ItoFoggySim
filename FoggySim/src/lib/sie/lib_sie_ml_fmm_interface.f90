module lib_sie_ml_fmm_interface
    use lib_tree
    use lib_tree_type
    use lib_ml_fmm
    use lib_ml_fmm_type_operator
    use lib_ml_fmm_helper_functions
    use ml_fmm_type
    use ml_fmm_math

    implicit none

    private

    public :: lib_sie_constructor
    public :: lib_sie_destructor

    public :: lib_sie_ml_fmm_calculate_vector_b

    contains

        ! Argument
        ! ----
        !   data_elements: type(lib_ml_fmm_data)
        !       positions of the data elements of the X, Y, and XY hierarchy
        !   s_opt: integer
        !       maximum number of elements at box (n, l_max)
        !
        subroutine lib_sie_constructor(data_elements, s_opt)
            implicit none
            ! dummy
            type(lib_ml_fmm_data), intent(inout) :: data_elements
            integer, intent(in) :: s_opt

            ! auxiliary

            type(ml_fmm_type_operator_procedures) :: ml_fmm_operator_procedures
            type(lib_ml_fmm_procedure_handles) :: ml_fmm_procedures


            ml_fmm_operator_procedures = ml_fmm_type_operator_get_procedures(1)
            ml_fmm_procedures = lib_sie_ml_fmm_get_procedures()

            call lib_ml_fmm_constructor(data_elements, &
                                        ml_fmm_operator_procedures, &
                                        ml_fmm_procedures, &
                                        tree_s_opt=s_opt, &
                                        final_sum_calc_y_hierarchy = .false., &
                                        final_sum_calc_xy_hierarchy = .true. ,&
                                        use_own_sum=.true.)

        end subroutine

        subroutine lib_sie_destructor()
            implicit none

            call lib_ml_fmm_destructor()

        end subroutine

        ! Argument
        ! ----
        !   vector_x: double complex, dimension(:)
        !       vector x: M x = b
        !   vector_b: double complex, dimension(:)
        !       vector x:
        subroutine lib_sie_ml_fmm_calculate_vector_b(vector_x, vector_b)
!            use lib_sie_data_container
            implicit none
            ! dummy
            double complex, dimension(:), intent(in) :: vector_x

            double complex, dimension(:), allocatable, intent(inout) :: vector_b

            ! auxiliary
            type(lib_ml_fmm_v), dimension(:), allocatable :: b
            type(lib_ml_fmm_v), dimension(:), allocatable :: x

            ! todo: pre-process: x

            call lib_ml_fmm_run(x, b)

            ! todo: post-process: b

        end subroutine


        function lib_sie_ml_fmm_get_procedures() result(handle)
            implicit none
            ! dummy
            type(lib_ml_fmm_procedure_handles) :: handle

            ! Upward Pass
            ! Step 1
            handle%get_c => null() ! todo: replace null() with own function
            ! Step 2
            handle%get_translation_SS => null() ! todo: replace null() with own function

            ! Downward Pass
            ! Step 1
            handle%get_translation_SR => null() ! todo: replace null() with own function
            ! Step 2
            handle%get_translation_RR => null() ! todo: replace null() with own function

            ! Final Summation
            handle%get_v_y_j => null() ! todo: replace null() with own function

        end function

        ! todo: define procedures
        ! - get_c
        ! - get_translation_SS
        ! - get_translation_SR
        ! - get_translation_RR
        ! - get_v_y_j
end module lib_sie_ml_fmm_interface
