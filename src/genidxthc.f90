
! Copyright (C) 2024 J. K. Dewhurst and S. Sharma.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine genidxthc
use modmain
use modtdhfc
implicit none
! local variables
integer ik,jk,nst,ist
! determine the maximum number of second-variational states within the energy
! window over all k-points
nmaxthc=0
do ik=1,nkpt
  nst=0
  do ist=1,nstsv
    if (abs(evalsv(ist,ik)-efermi) < ecutthc) nst=nst+1
  end do
  nmaxthc=max(nmaxthc,nst)
end do
if (nmaxthc == 0) then
  write(*,*)
  write(*,'("Error(genidxthc): no states within energy window ecutthc")')
  write(*,*)
  stop
end if
! allocate global arrays
if (allocated(nstthc)) deallocate(nstthc)
allocate(nstthc(nkptnr))
if (allocated(idxthc)) deallocate(idxthc)
allocate(idxthc(nmaxthc,nkptnr))
if (allocated(ithc)) deallocate(ithc)
allocate(ithc(nmaxthc,nkptnr))
! determine the number of second-variational states within the energy window for
! each k-point as well as the index to and total number of TDHFC states
nthc=0
do ik=1,nkptnr
  jk=ivkik(ivk(1,ik),ivk(2,ik),ivk(3,ik))
  nst=0
  do ist=1,nstsv
    if (abs(evalsv(ist,jk)-efermi) < ecutthc) then
      nst=nst+1
      idxthc(nst,ik)=ist
      nthc=nthc+1
      ithc(nst,ik)=nthc
    end if
  end do
  nstthc(ik)=nst
end do
end subroutine

