!=======================================================================
      SUBROUTINE PHOT_SEC (EK,Ephot,GI,nopen) 
!=======================================================================
!
!      CALCULATES PHOTOIONIZATION CROSS SECTIONS AT GIVEN ENERGY
!      (for given initial state and final partial wave)
!
!      INPUT:
!
!      EK           - elctron energy (Ry)
!      EPHOT        - photon energy (Ry)
!      GI           - statistical weight for initial state
!      DK           - DIPOLE VECTORS, 1=LENGTH, 2=VELOCITY.
!
!      OUTPUT:
!
!      DL,DV        - channels DIPOLE VECTORS (within -i factor)
!      SL,SV        - channels CROSS SECTIONS
!      SLP,SVP      - partial  CROSS SECTIONS
!
!-----------------------------------------------------------------------
 
      Use bsr_phot, EKK => EK

      Implicit real(8) (A-H,O-Z)

      Integer,Intent(in) :: nopen 
      Real(8),Intent(in) :: Ephot, EK, GI 

      Real(8),Parameter :: zero = 0.d0, one = 1.d0 

      Real(8),Parameter :: CONST = 2.689095 ! 4/3 * PI^2 * a0^2 * alfa, Mb

!-----------------------------------------------------------------------
!     DETERMINE THE ASYMPTOTIC SOLUTIONS WITH INGOING WAVE 
!     BOUNDARY CONDITIONS
!
!     F = S + C*K;    F- = -iF / (1 - iK); 
!
!     F-  -->  -i / sqrt(PI*K)  [ sin  + cos * K ] / (1+K^2) * (1+iK)  
!
!     Really, we need  R^-1 x F:   
!     from  F = aRF' - bRF  ->   R^-1 x F = a F' - b F
!-----------------------------------------------------------------------

      delta_p = 1.d8
      Do i=1,khm; delta_p=min(delta_p,abs(EK-eps(i))); End do
      delta_r = 1.d8
      Do i=1,nhm; delta_r=min(delta_r,abs(EK-value(i))); End do
      met = 0
      if(delta_p.lt.0.01.and.delta_p.lt.delta_r)  met = 1

write(pri,'(f15.8,2i5,2e15.5)') EK, met,khm,delta_p,delta_r

! ... define F and F' functions, and then R^-1 x F:

      FF(1:nch,1:nopen) = MATMUL(G(1:nch,1:nch),KMAT(1:nch,1:nopen))
      FF(1:nch,1:nopen) = F(1:nch,1:nopen) + FF(1:nch,1:nopen)

      if(met.eq.1) then
      FFP(1:nch,1:nopen) = MATMUL(GP(1:nch,1:nch),KMAT(1:nch,1:nopen))
      FFP(1:nch,1:nopen) = FP(1:nch,1:nopen) + FFP(1:nch,1:nopen)
      FF(1:nch,1:nopen) = ra*FFP(1:nch,1:nopen) - rb*FF(1:nch,1:nopen) 
      end if
       
! ... define  1 / (1+K^2),  where K - open-open part

      AA(1:nopen,1:nopen) = MATMUL(KMAT(1:nopen,1:nopen),KMAT(1:nopen,1:nopen))
      Do i = 1,nopen
        AA(i,i) = AA(i,i) + one
      End do
      Call INV(nopen,nch,AA)

! ... F- --> F / (1+K^2)

      BB(1:nch,1:nopen) = MATMUL(FF(1:nch,1:nopen),AA(1:nopen,1:nopen))              

! ... F- --> F / sqrt(pi)

      P = ACOS(-one)
      P = one / sqrt(P)
      FF(1:nch,1:nopen) = P * BB(1:nch,1:nopen)

!-----------------------------------------------------------------------
!     dipole matrix elements for given solution can be obtained from the
!     dipole matrix elements between initial state and R-matrix basis 
!     states (array DK) by weighting them with expansion coefficients 
!
!         A(k) = (1/2a) (E(k) - E) SUM(i) [ w(i,k) (R^-1) F(i) ]
!
!     for given solution j.
!
!-----------------------------------------------------------------------

! ... AA -->  R^(-1) *  F-

      if(met.eq.0) then
       RMATI=RMAT; Call INV(nch,nch,RMATI)
       AA(1:nch,1:nopen) = MATMUL(RMATI(1:nch,1:nch),FF(1:nch,1:nopen))            
      else
       AA(1:nch,1:nopen) = FF(1:nch,1:nopen)
      end if

! ... dipole m.e.

      DLr = zero; DLr = zero; DVr = zero; DVi = zero;

      RI = one/RA
 
      DO K = 1,NHM
        CK = RI / (VALUE(K)-EK)
        CL = CK * DKL(k)
        CV = CK * DKV(k)
        v(1:nopen) = MATMUL(WMAT(1:nch,k),AA(1:nch,1:nopen))
        DLr(1:nopen) = DLr(1:nopen) + CL*v(1:nopen)
        DVr(1:nopen) = DVr(1:nopen) + CV*v(1:nopen)
      END DO


! ... define the imagine part from ( 1 - i K) :

      DLi(1:nopen) = - MATMUL(DLr(1:nopen),KMAT(1:nopen,1:nopen))
      DVi(1:nopen) = - MATMUL(DVr(1:nopen),KMAT(1:nopen,1:nopen))      

! ... correction (should signs be important for angular stuff?):

!      Do i=1,nopen
!       s1 = DLr(i); s2=DLi(i); DLr(i)=-s2; DLi(i)=-s1      
!       s1 = DVr(i); s2=DVi(i); DVr(i)=-s2; DVi(i)=-s1      
!      End do
!-----------------------------------------------------------------------
!
!     sigma = (4/3 pi^2 ao^2 alfa) (w/2L+1) (dip.m.e)^2
!
!     where w - photon energy, 2L+1 (or 2J+1)  define the stat.weight of
!     initial state, w --> 4/w for velocity form
!
!     TOTAL CROSS SECTION IS IN BARN (10^-18 CM^2)
!-----------------------------------------------------------------------
 

      CL = CONST*Ephot/GI;  CV = 4.d0*CONST/(Ephot*GI)

      SL  = zero; SV  = zero; SLP = zero; SVP = zero
 
      DO J = 1,NOPEN
       SL(J) = CL * (DLr(J)**2 + DLi(J)**2)
       SV(J) = CV * (DVr(J)**2 + DVi(J)**2)
       SLP = SLP + SL(J)
       SVP = SVP + SV(J)
      END DO
      
      END SUBROUTINE PHOT_SEC
 
