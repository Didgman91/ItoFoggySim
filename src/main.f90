program main
    use lib_data_types
    use lib_tree
    use lib_tree_helper_functions
!    use lib_hash_function
    implicit none

    integer :: error_counter
!    call lib_test_hash_function()

    error_counter = 0
    error_counter = error_counter + lib_tree_hf_test_functions()
    error_counter = error_counter + lib_tree_test_functions()

    call lib_tree_benchmark()
!     call lib_tree_hf_benchmark()

    call lib_tree_destructor()
    call lib_tree_hf_destructor()

    print *, "-------------MAIN------------------"
    if (error_counter == 0) then
        print *, "All tests: OK"
    else
        print *, error_counter,"test(s) FAILED"
    end if
    print *, "-----------------------------------"

end program main
