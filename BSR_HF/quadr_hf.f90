!======================================================================
      Real(8) Function QUADR_hf(i,j,m)
!======================================================================
!     Evaluates   <P_i | r^m | P_j>     with respect to r
!----------------------------------------------------------------------
      Use bsr_hf
      Use hf_orbitals

      Implicit none
      Integer, intent(in) :: i,j,m
      Real(8), external :: BVMV
      Real(8) :: rm(ns,ks)

      if     ( m .eq. 1 ) then
        quadr_hf = BVMV (ns,ks, r1,'s',p(1,i),p(1,j))
      elseif ( m .eq. 0 ) then
        quadr_hf = BVMV (ns,ks, sb,'s',p(1,i),p(1,j))
      elseif ( m .eq.-1 ) then
        quadr_hf = BVMV (ns,ks,rm1,'s',p(1,i),p(1,j))
      elseif ( m .eq.-2 ) then
        quadr_hf = BVMV (ns,ks,rm2,'s',p(1,i),p(1,j))
      else
        Call MRM(m,rm)
        quadr_hf = BVMV(ns,ks,rm,'s',p(1,i),p(1,j))
      end if

      End Function  QUADR_hf
