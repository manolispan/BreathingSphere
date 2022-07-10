  module besfun1
    use params,only:dp
  contains 
    FUNCTION ENVJ(N,X)
        implicit none
        real(DP):: ENVJ
        real(DP):: X
        integer:: N
        ENVJ=0.5D0*DLOG10(6.28D0*N)-N*DLOG10(1.36D0*X/N)
        RETURN
        END
 end Module besfun1

  Module BesselFunctions
  use params, only:dp
  contains
     
       
COMPLEX*16 FUNCTION CONFRA( N, ZINV )
    !C
    !C ROUTINE CODED BY W. J. WISCOMBE
    !C
    !C COMPUTE BESSEL FUNCTION RATIO CAPITAL-A-SUB-N FROM ITS
    !C CONTINUED FRACTION USING LENTZ METHOD ( REF. 1, PP. 17-20 )
    !C
    !C ZINV = RECIPROCAL OF ARGUMENT OF CAPITAL-A
    !C
    !C I N T E R N A L    V A R I A B L E S
    !C ------------------------------------
    !C
    !C CAK TERM IN CONTINUED FRACTION EXPANSION OF CAPITAL-A
    !C ( REF. 1, EQ. 25 )
    !C CAPT FACTOR USED IN LENTZ ITERATION FOR CAPITAL-A
    !C ( REF. 1, EQ. 27 )

    !C CDENOM DENOMINATOR IN -CAPT- ( REF. 1, EQ. 28B )
    !C CNUMER NUMERATOR IN -CAPT- ( REF. 1, EQ. 28A )
    !C CDTD PRODUCT OF TWO SUCCESSIVE DENOMINATORS OF -CAPTC
    !FACTORS ( REF. 1, EQ. 34C )
    !C CNTN PRODUCT OF TWO SUCCESSIVE NUMERATORS OF -CAPTC
    !FACTORS ( REF. 1, EQ. 34B )
    !C EPS1 ILL-CONDITIONING CRITERION
    !C EPS2 CONVERGENCE CRITERION
    !C KK SUBSCRIPT K OF -CAK- ( REF. 1, EQ. 25B )
    !C KOUNT ITERATION COUNTER ( USED ONLY TO PREVENT RUNAWAY )
    !C MAXIT MAX. ALLOWED NO. OF ITERATIONS
    !C  MM + 1 AND - 1, ALTERNATELY
    !C
    implicit double precision (A-H,O-Z)
    INTEGER N
    COMPLEX*16 ZINV
    COMPLEX*16 CAK, CAPT, CDENOM, CDTD, CNUMER, CNTN, CDON
    integer maxit,mm,kk,kount   
    DATA EPS1 / 1.D-2 /, EPS2 / 4.4408920985006262D-16 /  
    DATA MAXIT / 10000 /
    ZABS(CDON)=(DMAX1(DABS(DREAL(CDON)),DABS(DIMAG(CDON))))*         &
        &            DSQRT(((DMAX1(DABS(DREAL(CDON)),DABS(DIMAG(CDON))))/           &
        &           (DMAX1(DABS(DREAL(CDON)),DABS(DIMAG(CDON)),1.0D0)))**2+        &
        &          ((DMIN1(DABS(DREAL(CDON)),DABS(DIMAG(CDON))))/                 &
        &          (DMAX1(DABS(DREAL(CDON)),DABS(DIMAG(CDON)),1.0D0)))**2)

    !C *** REF. 1, EQS. 25A, 27

    CONFRA = ( N + 1 ) * ZINV
    MM = - 1
    KK = 2 * N + 3
    CAK = ( MM * KK ) * ZINV
    CDENOM = CAK
    CNUMER = CDENOM + 1.0 / CONFRA
    KOUNT = 1
20  KOUNT = KOUNT + 1

    !C ***                                            REF. 2, EQ. 25B

    MM = - MM
    KK = KK + 2
    CAK = ( MM * KK ) * ZINV
    !C                                            *** REF. 2, EQ. 32
    IF ( ZABS( CNUMER/CAK ).LE.EPS1.OR. ZABS( CDENOM/CAK ).LE.EPS1 ) THEN
        !
        !C
        !C                                         ** ILL-CONDITIONED CASE -- STRIDE
        !C                                         ** TWO TERMS INSTEAD OF ONE
        !C
        !C                                                 *** REF. 2, EQS. 34
        CNTN = CAK * CNUMER + 1.0
        CDTD = CAK * CDENOM + 1.0
        CONFRA = ( CNTN / CDTD ) * CONFRA
        !C                                                 *** REF. 2, EQ. 25B
        MM = - MM
        KK = KK + 2
        CAK = ( MM * KK ) * ZINV
        !C                                                 *** REF. 2, EQS. 35
        CNUMER = CAK + CNUMER / CNTN
        CDENOM = CAK + CDENOM / CDTD
        KOUNT = KOUNT + 1
        GO TO 20
    !C
    ELSE
        !C                                                 ** WELL-CONDITIONED CASE
        !C
        !C                                                  *** REF. 2, EQS. 26, 27
        CAPT = CNUMER / CDENOM
        CONFRA = CAPT * CONFRA
        !C                                                  ** CHECK FOR CONVERGENCE
        !C                                                  ** ( REF. 2, EQ. 31 )
        !C
        IF ( DABS( DREAL (CAPT) - 1.0 ).GE.EPS2 .OR. DABS( DIMAG(CAPT) ) .GE.EPS2 ) THEN
            !C
            !C                                                 *** REF. 2, EQS. 30A-B
            CNUMER = CAK + 1.0 / CNUMER
            CDENOM = CAK + 1.0 / CDENOM
            GO TO 20
        END IF
    END IF
    !C
    RETURN
!C
END



end Module BesselFunctions


Module Complex_Bessel

!      REMARK ON ALGORITHM 644, COLLECTED ALGORITHMS FROM ACM.
!      THIS WORK PUBLISHED IN TRANSACTIONS ON MATHEMATICAL SOFTWARE,
!      VOL. 21, NO. 4, December, 1995, P.  388--393.

! Code converted using TO_F90 by Alan Miller
! Date: 2002-02-08  Time: 17:53:05
! Latest revision - 16 April 2002
!
!    This is a general driver for spherical Bessel functions
!    There is a problem with the dimensions since usually BJ(1) == j_0
!    but in the EM problems almost all arrays start from l=1. 
!    i.e. the j_0 y_0 functions are not required.
!    There is a confucion since normally to calculate bessel up to lmax
!    the arrays need to be dimensioned lmax+1 to deal with the extra l=0
!    terms.
!
!    Explicitly dimension (0:LMAXD) would solve the problem
!    The solution with dimension LMAXD+1 is somewhat problematic 
!    but this is how the old BESSEL from kkr works.
!    
!    All subroutines in this work with (0:LMAXD).
!      
!
IMPLICIT NONE
INTEGER, PARAMETER, PUBLIC  :: dpl = SELECTED_REAL_KIND(12, 60)    
INTEGER, PARAMETER   :: dp = kind(1.d0)

PRIVATE
PUBLIC  :: bes4


CONTAINS


subroutine BES4(besj0,besY0,besH0,ARG,LMAXD,LMAX,lj,ly,lh)
    use BesselFunctions,only: confra
    implicit none
    
    integer,parameter:: dp=kind(0.d0)
    integer,parameter:: nmaxd=7000
    integer,parameter:: lcrit = 50 ! continous fraction is used for lmax>lcrit
    !
    !      SUBROUTINE SPHBESS (arg,besj,besy,besh,lmax,J,JR,Y,H,MAXN)
    !C*****CALCULATION OF SPHERICAL BESSEL FUNCTIONS!!!
    !C JN IS CALCULATED FOR COMPLEX ARGUMENT USING THE LOGARITHMIC DERIVATIVE
    !C AS FIRST DONE BY INFELD (CODING BY J. V. DAVE)
    !C YN IS CALCULATED BY UPWARD RECURRENCE
    !C HN(1) IS THEN CALCULATED AS JN + SQRT(-1) * YN
    !C ADAPTED (WITH EXTENSIONS) IN MARCH 1988 FROM ROUTINE DBMIE OF
    !C J. V. DAVE (1968)
    !C IN PARTICULAR, DAVE USED EXP(-IKR)/R FOR OUTGOING WAVES, SO HE
    !C USED THE HANKEL FUNCTION OF THE SECOND KIND AND A COMPLEX
    !C REFRACTIVE INDEX WITH NEGATIVE IMAGINARY PART. THIS ROUTINE
    !C USES EXP(+IKR)/R FOR OUTGOING WAVES, HENCE HANKEL FUNCTIONS
    !C OF THE FIRST KIND AND A COMPLEX REFRACTIVE INDEX WITH
    !C POSITIVE REAL PART.
    !C*****COPYRIGHT C. D. CANTRELL, 1988

    integer :: lmaxd,lmax
    complex*16 YM1,YM2,Y(0:lmaxd)
    real*8 argr,argi,T
    complex*16  arg,Z
    COMPLEX*16 RRFX
    complex*16 ACAP(0:nmaxd),acap1(0:nmaxd)
    COMPLEX*16 J(0:lmaxd),JNM1
    complex*16 besj0(lmaxd),besy0(lmaxd),besh0(lmaxd),ci
    !complex*16 confra
    integer maxn,maxnm1,nmx1,nmx1p1,n,i,nn,k
    logical lj,ly,lh
    !C*****COMPLEX INDEX OF REFRACTION & ITS RECIPROCAL
    !write(6,*) 'arg',arg
    !write(6,*) 'lmaxd ,lmax',lmaxd,lmax
    besj0 = 0.d0
    besy0 = 0.d0
    besh0 = 0.d0
    if ((.not.lj).and.(.not.ly).and.(.not.lh)) return
    maxn= lmax+1
    ci = cmplx(0.d0,1.d0,kind=dp)
    RRFX = 1.d0/arg 
    argr = real(arg,dp)
    argi = imag(arg)
    !      RFRX=RFR*X
    !      RFIX=RFI*X
    !C*****CALCULATE UPPER LIMIT OF ORDER FOR DOWNWARD RECURRENCE
    T = arg*dconjg(arg) !(X**2)*(RFR**2 + RFI**2)
    T = SQRT(T)
    !write(6,*) 'T ',T
    MAXNM1 = nmaxd - 1
    !NMX2 = T
    NMX1 = max(1.1D0*T,float(5*lmax))

    if ( NMX1.le.150 ) nmx1 = 150

    if ( NMX1.gt.MAXNM1 ) then
        WRITE (6,"('error dim stop from bes4')")
        RETURN
    end if
    !   write(6,*) 'nmx1 ',nmx1
    ACAP(NMX1 + 1) = ( 0.D0, 0.D0 )
    NMX1P1 = NMX1 + 1
    if (nmx1p1 .gt. nmaxd) then
        write(6,*) 'Increase nmaxd, this is rare....'
        stop
    end if
    !write(6,*) 'nmx1 ',nmx1
    do N = 1,NMX1P1
        NN = NMX1 - N + 1
        ACAP(NN) = ( NN + 1 )*RRFX - 1.D0 / ((NN+1)*RRFX + ACAP(NN+1))
    end do
     !C WRITE (9,513) (ACAP(NN),NN,NN=0,NMX2)
    if (lj) then
        !C*****CALCULATION OF JN(RX*X) BY UPWARD RECURRENCE USING THE LOG. DERIVATIVE
        J(0)=DCMPLX (DSIN(argr)*DCOSH(argi) , DCOS(argr)*DSINH(argi))*rrfx
        JNM1=J(0)
        !Jnik(0) = JNM1
        !jnm1nik=jnik(0)
        if (lmax < lcrit ) then

            do NN=1,lmax
                J(NN)=JNM1/(ACAP(NN)+NN*RRFX)
                !Jnik(NN) = JNM1nik/(ACAP1(NN)+NN*RRFX)
                JNM1=J(NN)
                !jnm1nik = Jnik(nn)
            end do
        else
        ! Use continous fractions to obtain logarithmic derivative
        ! Intented for higher lmax and arg (WGM).
            do I = 1,nmx1
                ACAP1(i) = CONFRA ( I, rrfx )
                !WRITE (9,509) I, ACAP1(i),acap(i)
            end do
            !509        FORMAT(I5,4(1PD30.15))
            do nn=1,lmax
                J(NN) = JNM1/(ACAP1(NN)+NN*RRFX)
                jnm1 = J(nn)
            end do
        end if
        ! write(9,*) 'Bessel old'
        ! WRITE (9,513) (J(NN),NN,NN=0,lmax)
        ! write(9,*) 'Bessel new'
        ! WRITE (9,513) (Jnik(NN),NN,NN=0,lmax)
        ! Bessel j_l finished
        !
        besj0(1:lmax+1) = J(0:lmax)
    end if  ! Bessel calculated
513 FORMAT(2(1PD30.15),I10)
    if (ly .or. lh) then
        !
        !  *****CALCULATION OF NEUMANN FUNCTIONS BY UPWARD RECURRENCE
        !
        YM2 =  SIN(arg) * rrfx
        YM1 = -COS(arg) * rrfx
        Y(0) = YM1
        do NN=1,lmax
            Y(NN) = float(2*NN - 1) * rrfx * YM1 - YM2
            YM2=YM1
            YM1=Y(NN)
        !     write(6,*) nn,y(nn)
        end do
        !
        ! alternative calculate in terms of jl  check!!!!!!!!
           !
        !        z = arg
        !        Y(0)=-COS(Z)/Z
        !        Y(1)=(Y(0)-SIN(Z))/Z
        !       CDY(0)=(CDSIN(Z)+CDCOS(Z)/Z)/Z
        !       CDY(1)=(2.0D0*CDY(0)-CDCOS(Z))/Z
        !        DO 30 K=2,lmax
        !           IF (ZABS(J(K-1)).GT.ZABS(J(K-2))) THEN
        !              Y(K)=(J(K)*Y(K-1)-1.0D0/(Z*Z))/J(K-1)
        !           ELSE
        !              Y(K)=(J(K)*Y(K-2)-(2.0D0*K-1.0D0)/Z/Z/Z)/J(K-2)
        !           ENDIF
        !30      CONTINUE
        !
        ! end of alternative

        besy0(1:lmax+1) = Y(0:lmax)
    end if
      
    if(lh) then
         
        do i=1,LMAX+1
            !      besj0(i) = J(i)
            !      besy0(i) = Y(i)
            besh0(i)=J(i-1)+CI*Y(i-1)
        end do
    end if
    ! write(6,*) 'exiting ....'
    RETURN
END subroutine
!C***********************************************************************
 


END Module Complex_Bessel
