program main
    use lib_data_types
    use lib_tree
    use lib_hash_function
    implicit none

!    call lib_test_hash_function()

    call lib_tree_test_functions()

    call lib_tree_destructor()

!     call lib_tree_hf_benchmark()
!     call lib_tree_hf_test_functions()
!     call lib_tree_hf_benchmark()
!     call lib_tree_hf_destructor()


end program main
