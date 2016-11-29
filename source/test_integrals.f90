program jamess

use SlaterOverlap, only : GetSlaterOverlap

implicit none

double precision, parameter :: version = 0.01
double precision :: overlap

write (*,*) "Running JAMESS"

! Demonstration of Amber wrapper.
! Overlap between
! GetSlaterOverlap(na, la, nb, lb, mm, zeta_a, zeta_b, rab)
overlap = GetSlaterOverlap(2, 1, 1, 0, 0, 1.55555d0, 1.200d0, 1.1d0)

write (*,*) "Overlap", overlap

end program jamess
