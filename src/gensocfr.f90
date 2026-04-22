
! Copyright (C) 2009 J. K. Dewhurst, S. Sharma and E. K. U. Gross
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.

subroutine gensocfr
use modmain
use modomp
implicit none
! local variables
integer is,ias,nthd
integer nr,nri,ir,irc
real(8) cso,rm
! automatic arrays
real(8) vr(nrmtmax),dvr(nrmtmax)
if (.not.spinorb) return
! coefficient of spin-orbit coupling
cso=y00*socscf/(4.d0*solsc**2)
call holdthd(natmtot,nthd)
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE(vr,dvr,is,nr,nri) &
!$OMP PRIVATE(ir,irc,rm) &
!$OMP SCHEDULE(DYNAMIC) &
!$OMP NUM_THREADS(nthd)
do ias=1,natmtot
  is=idxis(ias)
  nr=nrmt(is)
  nri=nrmti(is)
! radial derivative of the spherical part of the Kohn-Sham potential
  call rfmtlm(1,nr,nri,vsmt(:,ias),vr)
  call splined(nr,wcrmt(:,:,is),vr,dvr)
  do ir=1,nr,lradstp
    irc=(ir-1)/lradstp+1
    rm=1.d0-2.d0*cso*vr(ir)
    socfr(irc,ias)=cso*dvr(ir)/(rsp(ir,is)*rm**2)
  end do
end do
!$OMP END PARALLEL DO
call freethd(nthd)

contains

pure subroutine splined(n,wc,f,df)
implicit none
! arguments
integer, intent(in) :: n
real(8), intent(in) :: wc(12,n),f(n)
real(8), intent(out) :: df(n)
! local variables
integer i
real(8) f1,f2,f3,f4
f1=f(1); f2=f(2); f3=f(3); f4=f(4)
df(1)=wc(1,1)*f1+wc(2,1)*f2+wc(3,1)*f3+wc(4,1)*f4
df(2)=wc(1,2)*f1+wc(2,2)*f2+wc(3,2)*f3+wc(4,2)*f4
!$OMP SIMD LASTPRIVATE(f1,f2,f3,f4) SIMDLEN(8)
do i=3,n-2
  f1=f(i-1); f2=f(i); f3=f(i+1); f4=f(i+2)
  df(i)=wc(1,i)*f1+wc(2,i)*f2+wc(3,i)*f3+wc(4,i)*f4
end do
i=n-1
df(i)=wc(1,i)*f1+wc(2,i)*f2+wc(3,i)*f3+wc(4,i)*f4
df(n)=wc(1,n)*f1+wc(2,n)*f2+wc(3,n)*f3+wc(4,n)*f4
end subroutine

end subroutine

