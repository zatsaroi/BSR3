!======================================================================
      Real(8) Function I_int(i,j)
!======================================================================
!     <P_i| L |P_j> with inclusion of rel.shift if rel = .true.
!----------------------------------------------------------------------
      Use bsr_hf
      Use hf_orbitals 

      Implicit none
      Integer, intent(in) :: i,j
      Real(8), external :: BVMV

      if (iabs(lbs(i)-lbs(j)) .NE. 0) Stop ' I_int:  LI <> LJ'

      Call HLM (lbs(i))

      I_int = -0.5d0*BVMV(ns,ks,hl,'s',p(1,i),p(1,j))
 
      End Function I_int

