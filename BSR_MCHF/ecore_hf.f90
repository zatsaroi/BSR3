!======================================================================
      Real(8) Function ECORE_hf()
!======================================================================
!     COMPUTES THE ENERGY OF THE COMMON CLOSED SHELLS  (CORE)
!----------------------------------------------------------------------
      Use bsr_mchf

      Implicit none

      Integer :: I,J,K
      Real(8) :: C, CC, S
      Real(8), external ::  BHL, RK

      S = 0.d0
      Do i = 1,ncore
       S = S - (4*lbs(i)+2)*BHL(i,i) / 2.d0
      End do

      Do k = 0,kmax
       Do i=1,ncore; Do j=1,i
        if(i.eq.j) c = (4*lbs(i)+2)*(4*lbs(j)+1)/2
        if(i.ne.j) c = (4*lbs(i)+2)*(4*lbs(j)+2)
        cc = c*coef(max(i,j),min(i,j),k)
        if(cc.ne.0.d0) S = S + cc*RK(i,j,i,j,k)
        cc = c*coef(min(i,j),max(i,j),k)
        if(cc.ne.0.d0.and.i.ne.j) S = S + cc*RK(i,j,j,i,k)
       End do; End do
      End do

      Ecore_hf = S

      End Function Ecore_hf
