Program breathingsphere

  use params, only: han, bes, dhan, dbes, dynte,dyntm,SCSbrsph,dyntmat
  implicit none

!-------------------PURPSOSE -----------------------------------!    
! Breahingsphere calculates the T-matrix scattering and         !
! absorption cross section of a homogeneous sphere              !
! (epssph, musph) with a time dependent radious RAD = RAD(t)    !
! in an environment (epsenv,muenv)                              !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!-------------------PARAMETERS ---------------------------------!    
! Below are the parameters that can be set from dyninput file   !
!                                                               !
! EPSSPH: The dielectric constant of the sphere                 !
! MUSPH: The magnetic permeability of the sphere                !
! EPSENV: The dielectric constant of the host                   !
! MUSPH: The magnetic permeability of the host                  !
! RAD: The radius of the sphere                                 !
! LMAX: Maximim L for the spherical wave expansion              !
! N0:  2*N0+1 are the number of calculated beams                !
! NFFT: Points for the fourier expansion                        !
! NSPLIT: Used for dampened oscillations. How many oscillations !
! until amplitude of oscillation resets                         !
! NLIGHT: How many light frequencies to calculate               !
! START1: First light frequency to calculate                    !
! END1: Last light frequency to calculate                       ! 
! G0: amplitute of the radius oscillation divided by radius     !
! NVIBRA: How many frequencies of oscillation to calculate      ! 
! STARTVIB: First oscillation frequency to calculate            !
! ENDVIB: Last oscillation frequency to calculate               !
! GAMMA : Coefficient of the dampening. See purpsose            !
! OPTOMECHANICAL: If set true vibration frequency is set to     !
! abs(Omega_resonance-omegalight)                               !
! OMEGA_RESONANCE: used for optomechanical if set TRUE          !
! SCALEFACTOR: If set to 1 then units are dimensionless         !
! zval ==== wmega*[R]/c. Set scalefactor= 0.1973269804 to use   !
! eV for energy and microns for lengths                         !
! SCSFILE: The file used for the output. Default is breath.scs  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



  integer,parameter :: dp=kind(1.0d0)
  integer :: lmax, N0, i,ii, ntime, niter, nkeep, ndyn, iter,idet,nc,ntest
  real(dp) :: g0, radious, em, mm, tau, gamma, split, correction,pi, taun, autocorrection
  real(dp) :: phase,scalefactor,domega,period,thz2ev,start1,end1,zval,onethz,omega_resonance
  complex(dp):: epssph,musph,epsenv,muenv,rad,rap,www0,abscs
  real(dp):: startvib,endvib,Resonance_vib,goff,gres
  integer :: nsplit,na,nb,lmax1d,np,nlight,ilight,nvibra,ivibra,ih,id,info,ngamma,l,i1,j1
  complex(kind=dp):: omegavib, w0,fval(-5:1000,50),wold,wnew,ww,derr,dh,f0,f1,f2,gamma1,gamma0
  complex(dp):: detomega,totalabs,totalcros,absorbn0
  complex(kind=dp), allocatable ::  wn(:), SCSno0(:) , te(:),th(:)
  integer,allocatable:: isort(:),ipivt(:)
  integer:: ipos,nfft
  complex(dp),allocatable:: ework(:),comvec(:,:),rc(:),dummyl(:,:),aux(:,:),dyntesq(:), &
       & dyntmsq(:), phasete(:),phasetm(:)
  real(dp),allocatable:: erwork(:),rrc(:)
  logical :: OPTOMECHANICAL
  character*20:: scsfile
  namelist /dyninput/ epssph,musph,epsenv,muenv,rad,lmax,n0,nfft,nsplit,&
       &nlight,start1,end1,g0,nvibra,startvib,endvib,gamma,OPTOMECHANICAL,&
       &omega_resonance,scalefactor,scsfile
900 format (E30.15, E30.15 ,E30.15, E30.15)
  open(99,file='adiabatic.txt')
  open(98,file='static_lscs.txt')
  pi = 4.d0*atan(1.d0)
  thz2ev = 0.004135667695680778d0
  write(98,*) '# r**2/2 norm : w1, scs,abs,ext, backsc, gtot, 8...lmax are scat, and after extinction '

  !READING INPUT

  open(16,file='dyninput')
  read(16,nml=dyninput)
  write(6,nml=dyninput)
  OPEN(432,FILE=trim(scsfile),ACCESS='SEQUENTIAL')
  write(432, "('# freq     SCScorr     SCS     SCSno0corr  ')")
  open(462,file='phasetcomponents')
  close(16)
  rad = rad / scalefactor
  goff=g0

  !END OF INPUT

  write (6,"('Calculate Light freq from :',D24.14,' to ',D24.14,' in eV,',D24.14,' with ',I6,' points' )") start1,end1 ,end1-start1,nlight
  write (6,"('Calculate Oscilation from :',F10.5,' to ',F10.5,' in eV, with ',I6,' points' )") &
       &  startvib,endvib ,nvibra
  write(6,"('Lmax',I4,' and NFF :',I4)") lmax,2*n0+1

  lmax1d = lmax + 1 
  allocate (te(lmax1d),th(lmax1d))

  na = 2*N0 + 1
  nb = 4*N0 + 1

  if (.not.allocated(dynte)) then
     allocate (dynte(na,na,lmax),dyntm(na,na,lmax))
  end if

  allocate(SCSbrsph(2*N0+1), SCSno0(2*N0+1))
  allocate( dyntesq(2*N0+1),  dyntmsq(2*N0+1), phasete(2*N0+1),phasetm(2*N0+1)) !manman
  allocate(wn(2*N0+1)) 
  
  !LOOP in frequency of light

  do ilight= 1,nlight
     zval = start1
     if (nlight>1) zval = start1 + (end1 -start1)/float(nlight-1)*float(ilight-1)
     rap = zval*rad/2.d0/pi        
     call TMTRX(lmax, lmax1d,rap, zval, epssph, epsenv, muenv, musph, te, th  )
     w0    =  zval

     ! LOOP in OMEGAS of oscillation

     do ivibra = 1, nvibra
        omegavib = startvib
        if (nvibra>1)   omegavib = startvib + (endvib -startvib)/float(nvibra-1)*float(ivibra-1)
        g0 = goff
        if (OPTOMECHANICAL) then
           omegavib = abs(omega_resonance - w0)
           write(6,"('Optomechanical mode w, vib, g0: ',5D20.12)") real(w0),real(omegavib),g0
        end if
        domega = real(omegavib)/float(nsplit)
        PERIOD = 2.d0*pi/domega  

        ! Start Adiabatic calculation

        ntime = 2*n0 +1           
        call adiabatic(omegavib,g0,gamma,ntime,period,rad,lmax,lmax1d,zval,  &
             &         epssph, epsenv, muenv, musph)

        write(6,"('Adiabatic finished for light freq no.:',I6,' , of ',I6,'  and vibration freq: ',I5,' , of ',I5)")&
             &ilight,nlight,ivibra,nvibra

        !Adiabatic finished

        do i=1,2*N0+1
           wn(i)= w0 - float(-N0-1+i)*domega 
        end do

        na = 2*N0 + 1
        nb = 4*N0 + 1

        call brsphere1(N0,NA,NB,nfft,lmax,period,w0,omegavib,gamma,phase,rad,                  &
             &           g0,abscs,epssph,musph,epsenv,muenv)
        print*,'brsphere finished!'

        !Calculate Dynamic SCS 

        SCSno0(1:2*N0+1)=SCSbrsph(1:2*N0+1)/rad/rad/pi
        SCSbrsph(1:2*N0+1)=SCSbrsph(1:2*N0+1)/rad/rad/pi
        totalcros = sum(SCSno0(1:2*N0+1))
        totalabs = abscs/rad/rad/pi
        dyntesq(1:2*N0+1)= dynte(1:2*N0+1,N0+1,1)*dconjg(dynte(1:2*N0+1,N0+1,1))
        dyntmsq(1:2*N0+1)= dyntm(1:2*N0+1,N0+1,1)*dconjg(dyntm(1:2*N0+1,N0+1,1))
        phasete(1:2*N0+1)=atan2(real(dynte(1:2*N0+1,N0+1,1)),aimag(dynte(1:2*N0+1,N0+1,1)))
        phasetm(1:2*N0+1)=atan2(real(dyntm(1:2*N0+1,N0+1,1)),aimag(dyntm(1:2*N0+1,N0+1,1)))
 
       write(462,909) real(w0), (real(dyntesq(i+n0+1)),i=-2,2),(real(dyntmsq(i+n0+1)),i=-2,2), &
             &(real(phasete(i+n0+1)),i=-2,2),(real(phasetm(i+n0+1)),i=-2,2) 

        write(432,901) real(w0), real(w0)/thz2ev, (real(SCSbrsph(i+n0+1)),i=-2,2),real(totalcros), &
             &real(abscs)
     end do 
  end do  

  do i=-20,20
     write(433,931) real(wn(i+n0+1)), real(SCSbrsph(i+n0+1))
  end do



901 format(10(E22.14,','))
931 format(2(E22.14,','))
909 format(21(E22.14,','))

  print*, "End of calculation "
end program breathingsphere



subroutine brsphere1(N0,NA,NB,nfft,lmax,period,w0,omegavib,gamma,phase,rad,g0,abscs,  &
     &                epssph,musph,epsenv,muenv)

  use params, only: dp,han, bes, dhan,dbes,dynte,dyntm,SCSbrsph
  implicit none
  integer,intent(in):: N0,NA,NB,lmax
  complex(dp),intent(in) :: epssph,musph,epsenv,muenv,w0,omegavib
  real(dp),intent(in):: phase,rad
  real(dp):: domega,period,deltat,prefl
  complex(dp)::absorbn0
  complex(dp),intent(out):: abscs
  integer :: i,nd,ii,iii, ne,INFO,NRHS, jj, nfft, l1, j, nfft0, nstart,nend
  integer, allocatable :: IPIVE(:),IPIVM(:)
  real(dp) :: scalefactor,g0, f0, e,m,em,mm,e0,m0, pi,gamma
  complex(kind=dp),allocatable:: amat(:,:,:),bmat(:,:,:),cmat(:,:,:),dmat(:,:,:)
  complex(kind=dp),allocatable :: emat(:,:,:),fmat(:,:,:)
  complex(dp),allocatable:: qn(:),wn(:),kn(:),qe(:,:), qm(:,:), imat(:),imatm(:)
  complex(dp),allocatable:: xm(:), xe(:), SCS(:)
  complex(dp):: ci, zeta, zeta1,zetam,absorb
  character*1 :: TRANS
  pi =4.0_dp*atan(1.0_dp)
  ci=cmplx(0.d0,1.d0, kind=dp)
  NRHS=1
  allocate (amat(NB,NA,lmax+1))
  allocate (bmat(NB,NA,lmax+1))
  allocate (cmat(NB,NA,lmax+1))
  allocate (dmat(NB,NA,lmax+1),emat(NB,NA,lmax+1))
  allocate (fmat(NB,NA,lmax+1))

  allocate(kn(NA),qn(NA),wn(NA))

  domega = 2.d0*pi/period
  do i=1,NA
     wn(i)= w0 - float(-N0-1+i)*domega 
  end do
  kn(1:NA)=wn(1:NA)*sqrt(epsenv*muenv) 
  qn(1:NA)=wn(1:NA)*sqrt(epssph*musph)                   
  if (.not.allocated(bes)) then
     allocate (bes(nfft,lmax+1),han(nfft,lmax+1),dbes(nfft,lmax+1),dhan(Nfft,lmax+1))
  end if

  !CONSTRUCT A, B, C D E F MATRICES 

  do i=1,NA
     call besselfourier1(lmax,NFFT,period,domega,phase,kn(i),rad,g0,omegavib,gamma)
     if (mod(nfft,2)==0) THEN
        NFFT0 = (NFFT+2)/2
     else
        NFFT0 = (NFFT+1)/2
     end if
     NSTART = NFFT0 - (NB-1)/2
     NEND   = NFFT0 + (NB-1)/2
     amat(1:NB, i, 1:lmax+1)=bes(NSTART:NEND,1:lmax+1)   
     bmat(1:NB, i, 1:lmax+1)=han(NSTART:NEND,1:lmax+1)
     fmat(1:NB, i, 1:lmax+1)=dbes(NSTART:NEND,1:lmax+1)
     emat(1:NB, i, 1:lmax+1)=dhan(NSTART:NEND,1:lmax+1)
     call besselfourier1(lmax,Nfft,period,domega,phase,qn(i),rad,g0,omegavib,gamma)
     cmat(1:nb, i, 1:lmax+1)=bes(NSTART:NEND, 1:lmax+1)
     dmat(1:nb, i, 1:lmax+1)=dbes(NSTART:NEND, 1:lmax+1)
  end do

  !CONSTRUCT Q MATRICES

  allocate(qe(4*N0+2,4*N0+2),qm(4*N0+2,4*N0+2))
  allocate(imat(4*N0+2),imatm(4*n0+2))
  allocate(IPIVE(4*N0+2),IPIVM(4*N0+2))
  allocate(xe(2*N0+1),xm(2*N0+1))
  allocate(SCS(2*N0+1))
  SCS(1:2*N0+1)=cmplx(0.0d0,0.0d0,kind=dp)
  absorbn0 = cmplx(0.0d0,0.0d0,kind=dp)
  abscs = 0.d0
  do l1 = 1,lmax
     zeta =  sqrt(epssph*muenv/epsenv/musph)
     zeta1 = sqrt(epsenv*muenv/epssph/musph)
     zetam = muenv/musph
     do i=1,Na
        do ii=1,Na
           Qe(i   , ii   )=     -bmat(NA-ii+i ,ii ,l1)
           Qe(i   , Na+ii)= zeta*cmat(NA-ii+i ,ii ,l1)
           Qe(Na+i, ii   )=     -emat(NA-ii+i ,ii ,l1)/kn(ii)
           Qe(Na+i, Na+ii)= dmat(NA-ii+i ,ii ,l1)/qn(ii) 
        end do
     end do
     do i=1,Na
        do ii=1,Na
           Qm(i,ii)=                  -bmat(NA-ii+i, ii, l1)
           Qm(i,Na+ii)=                cmat(NA-ii+i, ii, l1)
           Qm(Na+i,ii)=               -emat(NA-ii+i, ii, l1)/kn(ii)
           Qm(Na+i,Na+ii)=       zeta*dmat(NA-ii+i, ii, l1)/qn(ii)
        end do
     end do

     !CONSTRUCT I MATRIX

     TRANS='N'
     call zgetrf(4*N0+2,4*N0+2,QE,4*N0+2,IPIVE,INFO)
     call zgetrf(4*N0+2,4*N0+2,QM,4*N0+2,IPIVM,INFO)
     do jj=1,na
        do i=1,Na
           Imat(i )       = amat(2*N0 + 1 + i - jj, jj, l1)  
           Imat(2*N0+1+i) = fmat(2*N0 + 1 + i - jj, jj, l1)/kn(jj)  
        end do
        imatm(1:4*n0+2) =imat(1:4*n0+2)  
        call zgetrs(TRANS,4*N0+2,NRHS,QE,4*N0+2,IPIVE,imat,4*N0+2,INFO)
        xe(1:2*N0+1)=imat(1:2*N0+1)
        call zgetrs(TRANS,4*N0+2,NRHS,QM,4*N0+2,IPIVM,imatm,4*N0+2,INFO)
        xm(1:2*N0+1)=imatm(1:2*N0+1) 
        dynte(1:2*N0+1,jj,l1) = xe(1:2*N0+1)  
        dyntm(1:2*N0+1,jj,l1) = xm(1:2*N0+1)  
        prefl = float(2*l1+1)*2.d0*pi
        if ((jj==N0+1) .and. (abs(imag(w0))<1.d-14) ) then
           do i =1,2*N0+1
              SCS(i) = SCS(i) + prefl/kn(i)**2 * (xe(i)*dconjg(xe(i)) +xm(i)*dconjg(xm(i))) 
           end do
           absorbN0 = absorbN0 + prefl/kn(n0+1)**2* real(xe(N0+1) +xm(N0+1),dp)
        end if
     end do
  end do  
  abscs = -sum(SCS) - absorbN0
  SCSbrsph(1:2*n0+1)=SCS(1:2*n0+1)
  deallocate (amat,bmat,cmat,dmat,emat,fmat)
  deallocate (Qe, Qm, IPIVE,IPIVM, imat, imatm, xe, xm, scs)
  deallocate (kn,qn,wn)
end subroutine brsphere1


subroutine besselfourier1(lmax,nf,period,domega,phase,kvec,rad0,g0,omegavib,gamma)
  
  ! This sub calculates the fourier transform of
  !        - f_l ( kvec * R(t) ),   f are spherical bessel and hankel 
  !        - d/dR ( R * f( kvec*R(t) ) 
  !       Output is exported through a module  : bes,dbes,han,dhan 
  !       deltat : time discretization = 2*pi/domega/NF
  !       nf     : number of deltat. time = [0 ,.... (nf-1)*deltat]   nf points. 


  use complex_bessel
  use params, only: bes, dbes, han, dhan
  implicit none
  integer,parameter :: dp=kind(1.d0)
  integer,intent(in):: nf,lmax
  integer     :: n,i,lmax1,ll,ii, l1
  complex(dp) :: kvec

  real(dp),intent(in)    :: rad0,g0,gamma,phase
  complex(dp),intent(in) :: omegavib               
  real(dp),intent(in)   :: domega,period

  COMPLEX(DP),allocatable::  j0(:),y0(:),h0(:)
  integer:: nb
  real(dp):: T,pi,tend,tst,split,deltat
  real(dp),allocatable :: t0(:),tempreal(:), tempimag(:)
  complex(dp),allocatable :: arg(:)
  real(dp),allocatable::deltarad(:)
  integer :: lensav,lenwrk,ier
  real(8),allocatable:: wsave(:),work(:)
  complex(dp),allocatable :: cfun(:)

  pi = 4.0_dp*atan(1.0_dp)
  lensav = 3*NF + INT(LOG(REAL(NF))) + 4   
  lenwrk = 2*NF                                  
  allocate (t0(nf),arg(nf),deltarad(nf))
  allocate(j0(lmax+2),y0(lmax+2),h0(lmax+2),cfun(Nf))
  allocate (wsave(lensav), work(lenwrk) )     
  deltat = period/float(Nf) 
  do i=1,Nf
     t0(i) = phase + float(i-1) * deltat
  end do
  call functionradius(nf,t0,omegavib,g0,gamma,deltarad)
  arg(1:nf)=kvec*rad0*deltarad(1:nf)
  lmax1 = lmax + 1
  do i=1,nf  
     CALL BES4(J0,Y0,H0,ARG(i),lmax+2,lmax+1,.TRUE.,.TRUE.,.true.)
     bes(i,1:lmax+1) = j0(2:lmax+2)
     han(i,1:lmax+1) = h0(2:lmax+2)    
     do l1=1,lmax
        dbes(i,l1) = float(l1+1)*bes(i,l1) - arg(i)*bes(i,l1+1)   ! d [r*f_l(k*r)] / dr
        dhan(i,l1) = float(l1+1)*han(i,l1) - arg(i)*han(i,l1+1)
     end do
  end do
  
  ! Initialize wsave array
  
  call cfft1i(NF, wsave, lensav, ier )
  if (ier.ne.0) write(6,*) 'Error in FFT init, code:',ier
  
  ! FFT, look at https://people.sc.fsu.edu/~jburkardt/f_src/fftpack5.1/fftpack5.1.html
  
  do ll=1,lmax    
     cfun(1:Nf)=bes(1:Nf,ll)  
     call cfft1f( nf, 1, cfun, nf, wsave, lensav, work, lenwrk, ier )
     if (ier.ne.0) write(6,*) 'Error in FFT1, code:',ier
     call fftshift(cfun,nf)
     bes(1:Nf,ll)=cfun(1:Nf)  
     cfun(1:Nf)=han(1:Nf,ll)
     call cfft1f( nf, 1, cfun, nf, wsave, lensav, work, lenwrk, ier )  
     if (ier.ne.0) write(6,*) 'Error in FFT2, code:',ier
     call fftshift(cfun,nf)
     han(1:Nf,ll)=cfun(1:Nf)
     cfun(1:Nf)=dhan(1:Nf,ll)
     call cfft1f( nf, 1, cfun, nf, wsave, lensav, work, lenwrk, ier )  
     if (ier.ne.0) write(6,*) 'Error in FFT2, code:',ier
     call fftshift(cfun,nf)
     dhan(1:Nf,ll)=cfun(1:Nf)
     cfun(1:Nf)=dbes(1:Nf,ll)
     call cfft1f( nf, 1, cfun, nf, wsave, lensav, work, lenwrk, ier )  
     if (ier.ne.0) write(6,*) 'Error in FFT2, code:',ier
     call fftshift(cfun,nf)
     dbes(1:Nf,ll)=cfun(1:Nf)
  end do
       
  deallocate (t0,arg)
  deallocate (cfun )
  deallocate (wsave, work)

end subroutine besselfourier1

subroutine fftshift(a0,n)
  implicit none
  integer,parameter:: dp = kind(1.d0)  
  
  INTEGER, INTENT (IN) :: N
  complex(dp), INTENT (INOUT), DIMENSION (N) :: A0
  complex(dp),dimension (n) :: aux
  integer:: n1
  if (mod(n,2)==0) then
    aux(1:N/2) = a0(N/2+1:N)
    aux(N/2+1:N) = a0(1:N/2) 
    a0 = aux
 else
    n1 = n - 1    
    aux(1:n1/2) = a0(N1/2+2 : N)
    aux(N1/2+1:N) = a0(1:N1/2+1)
    a0 = aux
  end if   
  end subroutine fftshift



  subroutine functionradius(nf,t0,omega,g0,gamma,deltarad)
    implicit none
    integer,parameter :: dp=kind(1.d0)
    integer,intent(in) :: nf
    real(dp),intent(in):: t0(nf)
    real(dp),intent(in)    :: g0,omega, gamma   
    real(dp),intent(out) :: deltarad(nf)

    if (abs(gamma)<1.0d-20) then
       deltarad(1:nf)=1.d0 + g0*sin(omega*t0(1:nf))
    else   
       deltarad(1:nf)=1.d0 + g0*sin(omega*t0(1:nf))*exp(-t0(1:nf)*gamma)
    end if
  end subroutine functionradius


  SUBROUTINE TMTRX(lmax, lmax1d,rap, zval, epssph, epsmed, mumed, musph, te, th  )       
    use params, ONLY:DP
    !  use constants
    use Complex_Bessel
    IMPLICIT NONE
    !!     ------------------------------------------------------------------  
    !!> @brief This subroutine calculates the T-matrix for a homogeneous sphere
    !! @details
    !!     THIS SUBROUTINE  CALCULATES  THE  T-MATRIX FOR THE SCATTERING  
    !!     OF ELECTROMAGNETIC  FIELD  OF  WAVE-LENGHT LAMDA  BY A SINGLE  
    !!  @param[in]   SPHERE OF RADIUS S.  (RAP=S/LAMDA).  
    !!  @param[in]   EPSSPH : COMPLEX RELATIVE DIELECTRIC CONSTANT OF THE SPHERE.  
    !!  @param[in]   EPSMED : COMPLEX RELATIVE DIELECTRIC CONSTANT OF THE MEDIUM.  
    !!  @param[in]   LMAX   : MAXIMUM ANGULAR MOMENTUM 
    !!  @result  t-matrix goes in derived variable te and th
    !                                                                       c20.9.21      
    !     ------------------------------------------------------------------  
    !  
    ! ..  PARAMETER STATEMENTS  ..  
    !  
    !      INTEGER, PARAMETER:: LMAX1D=LMAXD+1  
    !  
    ! ..  SCALAR ARGUMENTS  ..  
    !  
    INTEGER, intent(in)         :: LMAX,lmax1d  
    COMPLEX(DP), intent(in):: EPSSPH,EPSMED,MUSPH,MUMED,RAP
    REAL(DP),intent(in)    :: ZVAL
    !  
    ! ..  ARRAY ARGUMENTS  ..  
    !  
    COMPLEX(DP),intent(out):: TE(LMAX1D),TH(LMAX1D)
    !      complex(dp),intent(out):: q0(2*lmax1d),qp(2*lmax1d)
    !  
    ! ..  LOCAL SCALARS  ..  
    !  
    INTEGER    ::   L1,LMAX1,L,M,LMAXI  
    COMPLEX(DP)::   C1,C2,C3,C4,C5,C6,AN,AJ,BN,BJ,ARG,ARGM,XISQ,XISQM
    COMPLEX(DP)::   AR,CI
    !   LOGICAL    ::   LCALL  
    !  
    ! ..  LOCAL ARRAYS  ..  
    !  
    COMPLEX(DP)::  J(LMAX1D+1),Y(LMAX1D+1),H(LMAX1D+1)  
    COMPLEX(DP)::  JM(LMAX1D+1),YM(LMAX1D+1),HM(LMAX1D+1)

    REAL(DP)   :: CRSSCA,DDET,CRSEXT,CRSABS,CRS0,pi,CRSSCA0(lmax1d),CRSEXT0(lmax1d)
    COMPLEX(DP):: DETE,DETH,backsc
    complex(DP):: cj(lmax1d),ch(lmax1d),cy(lmax1d)
    real(DP)   ::  fnu,gl,gtot,thz2ev
    integer         ::  nzh,nzj,ierrh,ierrj,kode,i,ii,ierry,nzy

    thz2ev = 0.004135665538538099d0  
    pi = 4.d0*atan(1.d0)
    ci = (0.d0,1.d0)
    LMAX1 = LMAX+1   
    XISQ = SQRT(EPSMED*MUMED)  
    if (Real(epsmed,dp)<0.and. Real(mumed,dp)<0) xisq = -xisq
    XISQM = SQRT(EPSSPH*MUSPH)
    if (Real(epssph,dp)<0.and. Real(musph,dp)<0) xisqm = -xisqm  
    AR=2.D0*PI*RAP  
    ARG = XISQ*AR  
    ARGM = XISQM*AR
    IF(LMAX1.GT.LMAX1D)  GO  TO   10 
    CALL BES4(J,Y,H,ARG,LMAX1D+1,LMAX1,.TRUE.,.TRUE.,.true.)
    CALL BES4(JM,YM,HM,ARGM,LMAX1D+1,LMAX1,.TRUE.,.TRUE.,.true.)
    C1 = EPSSPH - EPSMED  
    C2 = EPSMED*ARGM  
    C3 = -EPSSPH*ARG  
    C4 = MUSPH - MUMED  
    C5 = MUMED*ARGM  
    C6 = -MUSPH*ARG  
    do L1=2,LMAX1  
       AN = C1*L1*JM(L1)*Y(L1) + C2*JM(L1+1)*Y(L1) + C3*JM(L1)*Y(L1+1)  
       AJ = C1*L1*JM(L1)*J(L1) + C2*JM(L1+1)*J(L1) + C3*JM(L1)*J(L1+1)  
       BN = C4*L1*JM(L1)*Y(L1) + C5*JM(L1+1)*Y(L1) + C6*JM(L1)*Y(L1+1)  
       BJ = C4*L1*JM(L1)*J(L1) + C5*JM(L1+1)*J(L1) + C6*JM(L1)*J(L1+1)  
       TE(L1-1) = -AJ / (AJ + CI*AN)  
       TH(L1-1) = -BJ / (BJ + CI*BN)
    end do

    ! Cross Section

    CRSSCA=0.D0
    CRSEXT=0.D0
    backsc = 0.d0
    gtot = 0.d0
    DO L=1,LMAX
       crssca0(l) = float(L+L+1)*(TE(L)*DCONJG(TE(L))                  &
            &                  +TH(L)*DCONJG(TH(L)))
       CRSSCA=CRSSCA + float(L+L+1)*(TE(L)*DCONJG(TE(L))                  &
            &                  +TH(L)*DCONJG(TH(L)))
       crsext0(l) = float(L+L+1)*DREAL(TE(L)+TH(L))
       CRSEXT=CRSEXT + float(L+L+1)*DREAL(TE(L)+TH(L))

       ! Back scattering 
       BACkSC = backsc  + float(L+L+1)*(-1)**L*(th(l) - te(l))
    END DO 
    do l=1,lmax -1
       gl = l*(l+2.d0)/float(l+1)*real(te(l)*dconjg(te(l+1)) + th(l)*dconjg(th(l+1)))       &
            &                 + (2.0d0*l+1)/float(l*(l+1)) * real(te(l)*dconjg(th(l)))
       gtot = gtot + gl
    end do
    backsc = abs(backsc)
    CRS0=1.D0/(2.D0*PI*PI*RAP*RAP)
    gtot = gtot/CRSSCA/(PI*PI*RAP*RAP)/CRS0/CRSSCA
    CRSSCA=CRS0*CRSSCA
    CRSEXT=-CRS0*CRSEXT
    CRSABS=(CRSEXT-CRSSCA)
    write(98,120) zval,zval/thz2ev,CRSSCA,CRSABS,CRSEXT,backsc,gtot, crssca0(1:lmax),crsext0(1:lmax)
    RETURN  
10  WRITE(6,100) LMAX1,LMAX1D  
    STOP
111 FORMAT(10E16.8)     
112 format(I4,6D16.8)
120 format(200E22.14)
100 FORMAT(//10X,'FROM SUBROUTINE TMTRX :'/                                              &
         &         10X,'LMAX+1 =',I3,'  IS GREATER THAN DIMENSIONED:',I3)  
  END subroutine TMTRX


  subroutine adiabatic(omega,g0,gamma,ntime,period,rad,lmax,lmax1d,zval,epssph, epsenv, muenv, musph)
    ! PRB Gantzounis et al PRB 84, 104303 (2011)
    use params,only:dp
    implicit none
    real(dp),intent(in) :: g0,gamma,period,zval,rad
    integer,intent(in):: ntime,lmax,lmax1d
    complex(dp),intent(in) :: omega,epssph, epsenv, muenv, musph
    complex(dp)::kn,kn0
    integer:: lensav,lenwrk,itime,l,ier,i,nn,n0
    real(dp),allocatable:: deltarad(:),wsave(:), work(:),t0(:),adiascs(:)
    real(dp):: rap,pi,crs0,sum1,thz2ev,cs0
    complex(dp),allocatable:: te(:),th(:),cfun(:),tte(:,:),tth(:,:)
    pi=4.d0*atan(1.d0)
    thz2ev = 0.004135665538538099d0 
    if (ntime==0) return
    allocate(deltarad(ntime),t0(ntime),te(lmax1d),th(lmax1d),cfun(ntime))
    allocate(tte(ntime,lmax1d),tth(ntime,lmax1d),adiascs(ntime))
    lensav = 3*Ntime + INT(LOG(REAL(Ntime))) + 4   
    lenwrk = 2*Ntime                               
    allocate (wsave(lensav), work(lenwrk) )  
    if (ntime<2) stop 'ADIABATIC increase ntime > 1'
    do itime =1 ,ntime
       t0(itime) = float(itime-1)*period/float(ntime)
    end do
    call functionradius(ntime,t0,real(omega),g0,gamma,deltarad)
    cfun = 0.d0
    adiascs = 0.0d0
    do itime =1,ntime
       rap = zval*rad*deltarad(itime)/2.d0/pi   
       call TMTRX(lmax, lmax1d,rap, zval, epssph, epsenv, muenv, musph, te, th  )
       do l=1,lmax
          tte(itime,l) = TE(L)
          tth(itime,l) = Th(l)
       end do
    end do
    adiascs = 0.D0 
    call cfft1i(Ntime, wsave, lensav, ier )
    do l=1,lmax
       cfun(1:ntime)= tte(1:ntime,l)
       call cfft1f( ntime, 1, cfun, ntime, wsave, lensav, work, lenwrk, ier )  
       if (ier.ne.0) write(6,*) 'Error in FFT adiabatic,   code:',ier        
       call fftshift(cfun,ntime)
       tte(1:ntime,l) = cfun(1:ntime)
       cfun(1:ntime)= tth(1:ntime,l)
       call cfft1f( ntime, 1, cfun, ntime, wsave, lensav, work, lenwrk, ier )  
       if (ier.ne.0) write(6,*) 'Error in FFT adiabatic,   code:',ier
       call fftshift(cfun,ntime)
       tth(1:ntime,l) = cfun(1:ntime)
       do itime = 1,ntime
          adiascs(itime) = adiascs(itime) + float(2*l+1)*            &
               & (tte(itime,l)*conjg(tte(itime,l)) + tth(itime,l)*conjg(tth(itime,l)) )
       end do
    end do
    nn = 2
    n0 = (ntime-1)/2
    do i= n0 + 1, n0 + 1 + nn 
       kn = (zval + float( i - n0 - 1 )*(2.0*pi/period))*sqrt(epsenv*muenv)
       kn0 = zval*sqrt(epsenv*muenv)
       cs0 = 2.d0/(zval*zval*rad*rad)
       write(99,"(20E22.14)") zval,adiascs(i)*cs0
    end do
    if (allocated(wsave)) deallocate (wsave, work)   
    if (allocated(tth))deallocate (tth,tte,adiascs)
    if (allocated(deltarad)) deallocate(deltarad,t0,te,th,cfun)
  end subroutine adiabatic

  subroutine buildfullt(lmax,n0,nkeep)   !
    ! This subroutine builds the full dynamical t-matrix ( dyntmat)  from the blocks 
    ! dynte,dynth The new matrix is truncated to -nkeep , + nkeep beams 
    !  Dimensions are NDF x NDF     NDF = (2*nkeep+1)*2*lmax*(lmax+2) 
    use params, only: dyntmat,dynte,dyntm
    integer,intent(in) :: lmax,n0,nkeep
    integer :: lma,ib,i,j,jb,i0,j0,l,m,il
    LMA = lmax*(lmax+2)
    dyntmat = 0.d0
    write(6,*) 'build lmax, n0, nkeep: ',lmax,n0,nkeep
    do ib = -nkeep, nkeep
       i = ib + nkeep +1
       do jb= - nkeep, nkeep
          j= jb + nkeep + 1
          i0 = (i-1)*2*LMA
          j0 = (j-1)*2*LMA
          il = 0
          do l=1,lmax
             do m=-l,l
                il = il + 1
                dyntmat(i0+il,j0+il)         = dynte( N0 + 1 + ib , N0 + 1 + jb, l)
                dyntmat(i0+il+LMA,j0+il+LMA) = dyntm( N0 + 1 + ib,N0 + 1 + jb, l)
             end do
          end do
       end do
    end do
  end subroutine buildfullt

  SUBROUTINE rSORT0 (W,IND,MAX,POS,igd)
    ! ************************************************************************
    !     W   is the original array returned unchanged
    !     IND is an array that holds the new positions 
    !     max number of elements to be sorted
    !     pos the position where the first element is found
    ! ------------------------------------------------------------------------
    use params,only:DP
    implicit none
    INTEGER :: MAX,POS,igd
    real(kind=DP)::  W(igd)
    INTEGER :: IND(igd)
    INTEGER ::I,II,J,JJ,K
    real(kind=DP):: BOUND =1.0D-12
    real(kind=DP):: diff
    DO  I = 1,MAX
       IND(I) = I
    END DO
    J = MAX
    if (max.gt.5) J = 1
    DO WHILE (J.LT.MAX/3)
       J = 3*J+1
    END DO
    DO WHILE (J.GT.1)
       J = J/3
       JJ = 1
       DO WHILE (JJ.EQ.1)
          JJ = 0
          DO K=1,MAX-J
             DIFF = ABS( W(IND(K)) - W(IND(K+J)) )
             IF ( W(IND(K)) .GT. W(IND(K+J)) .AND.                               &
                  &           DIFF.GT.BOUND ) THEN
                II       = IND(K)
                IND(K)   = IND(K+J)
                IND(K+J) = II
                JJ = 1
             END IF
          END DO                    
       END DO                      
    END DO
    print*,'rsort3'
    DO  I=1,MAX
       IF (IND(I) .EQ. 1) POS=I
    END DO
    RETURN
  END SUBROUTINE rSORT0

