!======================================================================
      Real(8) Function BHL_hf(i,j)
!======================================================================
!     <P_i| L |P_j>  with inclusion of rel.shift if rel = .true.
!----------------------------------------------------------------------
      Use bsr_hf
      Use hf_orbitals

      Implicit none
      Integer, intent(in) :: i,j
      Real(8), external :: BVMV

      if (iabs(lbs(i)-lbs(j)) .NE. 0) Stop ' HL:  LI <> LJ'

      Call HLM (lbs(i))

      BHL_hf = BVMV(ns,ks,hl,'a',p(1,i),p(1,j))

      End Function BHL_hf

