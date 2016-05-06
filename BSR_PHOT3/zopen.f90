!=======================================================================
      Subroutine ZOPEN(ETOT,NAST,ENAT,NCONAT,nch,nopen,EN)
!=======================================================================
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
      Implicit none
      Integer, intent(in)  :: NAST, nch
      Integer, intent(out) :: nopen
      Integer, intent(in)  :: NCONAT(NAST)
      Real(8), intent(in)  :: ETOT, ENAT(NAST)
      Real(8), intent(out) :: EN(nch)
      Integer :: K,N,NC
      Real(8) :: ECH
       
      NOPEN = 0
      K = 0
      Do N = 1,NAST
        IF (NCONAT(N).EQ.0) Cycle
        ECH = ETOT - ENAT(N)
        if (ECH.GT.0.d0) NOPEN = NOPEN + NCONAT(N)
        Do NC = 1,NCONAT(N)
          K = K + 1
          EN(K) = ECH
        End do
      End do
 
      End Subroutine ZOPEN
