module lib_ml_fmm_data_container
    use  ml_fmm_type
    implicit none

    ! e.g. matrix vector product v = u*phi
    ! order restiction:
    !   1. elements of the X-hierarchy
    !   2. elements of the XY-hierarchy
    type(lib_ml_fmm_v), dimension(:), allocatable :: m_ml_fmm_u
    ! e.g. matrix vector product v = u*phi
    ! order restiction:
    !   1. elements of the Y-hierarchy
    !   2. elements of the XY-hierarchy
    type(lib_ml_fmm_v), dimension(:), allocatable :: m_ml_fmm_v
end module lib_ml_fmm_data_container
