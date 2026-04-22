
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine eveqnit(nmatp,ngp,igpig,vpl,vgpl,vgpc,apwalm,evalfv,evecfv)
use modmain
use modomp
implicit none
! arguments
integer, intent(in) :: nmatp,ngp,igpig(ngkmax)
real(8), intent(in) :: vpl(3),vgpl(3,ngkmax),vgpc(3,ngkmax)
complex(8), intent(in) :: apwalm(ngkmax,apwordmax,lmmaxapw,natmtot)
real(8), intent(out) :: evalfv(nstfv)
complex(8), intent(out) :: evecfv(nmatmax,nstfv)
! local variables
integer ns,ist,it,i,j
integer is,ias,nthd
real(8) t1
real(8) ts1,ts0
! allocatable arrays
integer, allocatable :: idx(:)
real(8), allocatable :: w(:)
complex(8), allocatable :: h(:,:),o(:,:),hv(:,:),ov(:,:)
complex(8), allocatable :: u(:,:),hu(:,:),ou(:,:)
complex(8), allocatable :: hs(:,:),os(:,:),vs(:,:)
! external functions
real(8), external :: ddot
! compute Hamiltonian and overlap matrices
call timesec(ts0)
allocate(h(nmatp,nmatp),o(nmatp,nmatp))
call holdthd(2,nthd)
!$OMP PARALLEL SECTIONS DEFAULT(SHARED) &
!$OMP PRIVATE(j,ias,is) &
!$OMP NUM_THREADS(nthd)
!$OMP SECTION
! Hamiltonian
call hmlistl(ngp,igpig,vgpc,nmatp,h)
do j=ngp+1,nmatp
  h(1:j,j)=0.d0
end do
do ias=1,natmtot
  is=idxis(ias)
  call hmlaa(tefvr,is,ias,ngp,apwalm(:,:,:,ias),nmatp,h)
  call hmlalo(is,ias,ngp,apwalm(:,:,:,ias),nmatp,h)
  call hmllolo(is,ias,ngp,nmatp,h)
end do
!$OMP SECTION
! overlap
call olpistl(ngp,igpig,nmatp,o)
do j=ngp+1,nmatp
  o(1:j,j)=0.d0
end do
do ias=1,natmtot
  is=idxis(ias)
  call olpaa(tefvr,is,ngp,apwalm(:,:,:,ias),nmatp,o)
  call olpalo(is,ias,ngp,apwalm(:,:,:,ias),nmatp,o)
  call olplolo(is,ias,ngp,nmatp,o)
end do
!$OMP END PARALLEL SECTIONS
call freethd(nthd)
call timesec(ts1)
!$OMP ATOMIC
timemat=timemat+ts1-ts0
if ((iscl >= 2).or.(nefvit < 0)) then
! read in the eigenvectors from file if available
  call getevecfv(filext,0,vpl,vgpl,evecfv)
else
! initialise the eigenvalues/vectors using the diagonal elements of H
  allocate(idx(nmatp),w(nmatp))
  do i=1,nmatp
    w(i)=h(i,i)
  end do
  call sortidx(nmatp,w,idx)
  evecfv(1:nmatp,1:nstfv)=0.d0
  do ist=1,nstfv
    evecfv(idx(ist),ist)=1.d0
  end do
  deallocate(idx,w)
end if
call timesec(ts0)
ns=2*nstfv
allocate(hv(nmatp,nstfv),ov(nmatp,nstfv))
allocate(u(nmatp,nstfv),hu(nmatp,nstfv),ou(nmatp,nstfv))
allocate(hs(ns,ns),os(ns,ns),vs(ns,nstfv))
! iteration loop
do it=1,abs(nefvit)
  call holdthd(nstfv,nthd)
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(ist,t1) &
!$OMP NUM_THREADS(nthd)
!$OMP DO SCHEDULE(DYNAMIC)
  do ist=1,nstfv
! operate with O on the current eigenvector
    call zhemv('U',nmatp,zone,o,nmatp,evecfv(:,ist),1,zzero,ov(:,ist),1)
! normalise the eigenvector
    t1=ddot(2*nmatp,evecfv(:,ist),1,ov(:,ist),1)
    if (t1 > 0.d0) then
      t1=1.d0/sqrt(t1)
      call zdscal(nmatp,t1,evecfv(:,ist),1)
      call zdscal(nmatp,t1,ov(:,ist),1)
    end if
! operate with H on the current eigenvector
    call zhemv('U',nmatp,zone,h,nmatp,evecfv(:,ist),1,zzero,hv(:,ist),1)
! estimate the eigenvalue
    evalfv(ist)=ddot(2*nmatp,evecfv(:,ist),1,hv(:,ist),1)
! compute the residual |u⟩ = (H - eO)|v⟩
    t1=-evalfv(ist)
    u(1:nmatp,ist)=hv(1:nmatp,ist)+t1*ov(1:nmatp,ist)
! apply the overlap matrix to the residual
    call zhemv('U',nmatp,zone,o,nmatp,u(:,ist),1,zzero,ou(:,ist),1)
! apply the Hamiltonian matrix to the residual
    call zhemv('U',nmatp,zone,h,nmatp,u(:,ist),1,zzero,hu(:,ist),1)
  end do
!$OMP END DO
! compute the Hamiltonian and overlap matrices in the subspace formed by the
! eigenvectors and their residuals
!$OMP DO SCHEDULE(DYNAMIC)
  do ist=1,nstfv
    call zgemv('C',nmatp,nstfv,zone,evecfv,nmatmax,hv(:,ist),1,zzero, &
     hs(1,ist),1)
    call zgemv('C',nmatp,nstfv,zone,evecfv,nmatmax,hu(:,ist),1,zzero, &
     hs(1,nstfv+ist),1)
    call zgemv('C',nmatp,nstfv,zone,u,nmatp,hu(:,ist),1,zzero, &
     hs(nstfv+1,nstfv+ist),1)
  end do
!$OMP END DO NOWAIT
!$OMP DO SCHEDULE(DYNAMIC)
  do ist=1,nstfv
    call zgemv('C',nmatp,nstfv,zone,evecfv,nmatmax,ov(:,ist),1,zzero, &
     os(1,ist),1)
    call zgemv('C',nmatp,nstfv,zone,evecfv,nmatmax,ou(:,ist),1,zzero, &
     os(1,nstfv+ist),1)
    call zgemv('C',nmatp,nstfv,zone,u,nmatp,ou(:,ist),1,zzero, &
     os(nstfv+1,nstfv+ist),1)
  end do
!$OMP END DO
!$OMP END PARALLEL
  call freethd(nthd)
! solve the generalised eigenvalue problem in the subspace
  call eveqnzhg(ns,nstfv,ns,hs,os,evalfv,ns,vs)
! construct the new eigenvectors
  call holdthd(nstfv,nthd)
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(ist) &
!$OMP NUM_THREADS(nthd)
!$OMP DO SCHEDULE(DYNAMIC)
  do ist=1,nstfv
    call zgemv('N',nmatp,nstfv,zone,evecfv,nmatmax,vs(1,ist),1,zzero, &
     ov(:,ist),1)
    call zgemv('N',nmatp,nstfv,zone,u,nmatp,vs(nstfv+1,ist),1,zone,ov(:,ist),1)
  end do
!$OMP END DO
!$OMP DO SCHEDULE(DYNAMIC)
  do ist=1,nstfv
    call zcopy(nmatp,ov(:,ist),1,evecfv(:,ist),1)
  end do
!$OMP END DO
!$OMP END PARALLEL
  call freethd(nthd)
! end iteration loop
end do
deallocate(h,o,hv,ov,u,hu,ou,hs,os,vs)
call timesec(ts1)
!$OMP ATOMIC
timefv=timefv+ts1-ts0
end subroutine

