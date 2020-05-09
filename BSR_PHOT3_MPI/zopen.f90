!=======================================================================
      SUBROUTINE ZOPEN(ETOT,NAST,ENAT,NCONAT,nch,nopen,EN)
!=======================================================================
!
!      DETERMINE NUMBER OF OPEN CHANNELS AND CHANNELS ENERGIES
!
!      INPUT:
!
!      ETOT       = ELECTRON ENERGY IN RYDS
!      NAST       = NUMBER OF TARGET STATES
!      ENAT       = TARGET ENERGIES IN RYDS
!      NCONAT     = NUMBER OF CHANNELS COUPLED TO EACH TARGET STATE
!      NCH        = NUMBER OF CHANNELS
!
!      RETURNS:
!
!      NOPEN      = NUMBER OF OPEN CHANNELS
!      EN         = CHANNEL ENERGIES IN RYDS.
!-----------------------------------------------------------------------
 
      IMPLICIT NONE
 
      INTEGER(4), INTENT(in) :: NAST, nch
      INTEGER(4), INTENT(out) :: nopen
      INTEGER(4), INTENT(in), DIMENSION(NAST) :: NCONAT
      
      REAL(8), INTENT(in) :: ETOT
      REAL(8), INTENT(in), DIMENSION(NAST) :: ENAT
      REAL(8), INTENT(out), DIMENSION(nch) :: EN
      
      Integer(4) :: K,N,NC
      Real(8) :: ECH
       
      NOPEN = 0
      K = 0
      DO N = 1,NAST
        IF (NCONAT(N).EQ.0) Cycle
        ECH = ETOT - ENAT(N)
        IF (ECH.GT.0.d0) NOPEN = NOPEN + NCONAT(N)
        DO NC = 1,NCONAT(N)
          K = K + 1
          EN(K) = ECH
        END DO
      END DO
 
      END SUBROUTINE ZOPEN
 
