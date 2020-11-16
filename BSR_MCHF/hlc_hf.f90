!======================================================================
      Real(8) Function HLC_hf(I,J)           result(S)
!======================================================================
!     COMPUTES HL(I,J) MODIFIED BY THE INTERACTIONS WITH CLOSED SHELLS
!----------------------------------------------------------------------
      Use spline_orbitals, l => lbs
      Use spline_atomic, nclosd => kclosd

      Implicit none
      Integer, intent(in) :: i,j
      Integer :: ip, k, kmin, kmax
      Real(8) :: CA, CB
      Real(8), external :: ZCB, RK, BHL

      S = BHL(I,J)

! ... direct interaction:

      Do IP = 1,NCLOSD
       CA = -2.d0*(4*L(IP)+2)
       S = S + CA*RK(I,IP,J,IP,0)
      End do

! ... exchange contribution:

      kmin = 1000;  kmax = 0
      Do IP = 1,NCLOSD
       k = iabs(L(I)-L(IP));  if(k.lt.kmin) kmin = k
       k =      L(I)+L(IP) ;  if(k.gt.kmax) kmax = k
      End do

      Do k = kmin,kmax
       Do IP = 1,NCLOSD
        CB =  ZCB(L(I),K,L(IP)) * (4*L(IP)+2)
        if(CB.eq.0.d0) Cycle
        S = S + CB*RK(I,IP,IP,J,k)
       End do
      End do

      End Function HLC_hf
