!======================================================================
      Subroutine get_V (ii,jj,v)
!======================================================================
!     provides <.|nl>  for <kl|nl>; ii -> kl; jj -> nl;
!     if <kl|n'l'> = 0  then  nl -> nl - <nl|n'l'> n'l' 
!----------------------------------------------------------------------

      Use bsr_mat
      Use spline_param
      Use spline_orbitals
      
      Implicit none

      Integer, Intent(in) :: ii,jj
      Integer :: i,io
      Real(8) :: C,v(ns)

      if(IBORT(ii,jj).eq.0) Stop 'get_V: <kl|nl> = 0 ?'
      v(:) = qbs(:,jj)
      io = 0; C = 0.d0
      Do i=1,nbf
       if(i.eq.jj) Cycle
       if(iech(i).ne.0) Cycle
       if(lbs(i).ne.lbs(ii)) Cycle
       if(abs(OBS(i,jj)).lt.eps_ovl) Cycle
       if(IBORT(i,ii).ne.0) Cycle
       io = i; C = OBS(i,jj)
       v(:) = v(:) - C*qbs(:,i)
      End do

      End Subroutine get_V

