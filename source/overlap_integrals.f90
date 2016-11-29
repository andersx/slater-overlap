module overlap_integrals

contains


! Nice wrapper to get overlap integrals between atoms a and b
! Just a wrapper for crap legacy code from who know where! :D
subroutine get_overlap_block(a, za, lmaxa, b, zb, lmaxb, r, overlaps)

    implicit none

    ! Id numbers for atoms a and b
    integer, intent(in) :: a, b

    ! Basis set exponents (zeta) for a and b
    double precision, intent(in) :: za, zb

    ! Maximal angular momentum in shells of atoms a and b
    integer, intent(in) :: lmaxa, lmaxb

    ! Distance vector from a to b
    double precision, dimension(3), intent(in) :: r

    ! Overlap integrals (ss, spx, spy, spz, sdxx, sdyy, sdzz, etc)
    double precision, allocatable, dimension(:,:), intent(out) :: overlaps

    if (.not. allocated(overlaps)) allocate(overlaps(9,9))

    


    

end subroutine get_overlap_block

end module overlap_integrals

