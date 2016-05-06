!=======================================================================
      Subroutine AK_COEF(nopen,EK)
!=======================================================================
!     calculates weights of different components in internal region
!
!     INPUT:
!
!     EK       -  ELECTRON ENERGY (Ry)
!     AA       -  RMAT^(-1) * F-   (from module bsr_phot)
!     WT       -  weights of the R_matrix basis components
!
!     RETURNS:
!
!     AK(i,j)  -  weight of ith channel in solution j (it is c**2 !)
!
!-----------------------------------------------------------------------
      Use bsr_phot

      Implicit real(8) (A-H,O-Z)
      Integer, intent(in) :: nopen
      Real(8), intent(in) :: EK

      AK = 0.d0; CC = 0.d0; WTch = 0.d0; if(nwt.le.nopen) Return
!-----------------------------------------------------------------------
! ... all weights:

      RI = 1.d0/RA
      Do K = 1,NHM; kk = nhm-k+1
        CK = RI / (VALUE(kk)-EK)
        Do J = 1,nopen
          CJ = CK * SUM(WMAT(1:nch,kk)*AA(1:nch,J))
          AK(1:nwt,J) = AK(1:nwt,J) + WT(1:nwt,k) * CJ*CJ
        End do
      End do

      BB(1:nch,1:nopen)=MATMUL(AA(1:nch,1:nopen),KMAT(1:nopen,1:nopen))

      Do K = 1,NHM; kk = nhm-k+1
        CK = RI / (VALUE(kk)-EK)
        Do J = 1,nopen
          CJ = CK * SUM(WMAT(1:nch,kk)*BB(1:nch,J))
          AK(1:nwt,J) = AK(1:nwt,J) + WT(1:nwt,k) * CJ*CJ
        End do
      End do
!-----------------------------------------------------------------------
!     WTch(1:nopen) - contribution of all closed channels
!                     in the solution j
!     WTch(nopen+1:nwt) - relative contribution of ith closed channal
!                         in all solutions
!-----------------------------------------------------------------------
      Do j=1,nopen;  WTch(j) = SUM(AK(nopen+1:nwt,j));  End do
      CC = SUM(WTch(1:nopen))

      Do i = nopen+1,nwt;  WTch(i) = SUM(AK(i,1:nopen));  End do
      WTch(nopen+1:nwt) = WTch(nopen+1:nwt) / CC

      End Subroutine AK_COEF


