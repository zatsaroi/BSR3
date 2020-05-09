!======================================================================
      Subroutine get_V (kl,nl,v)
!======================================================================
!     provides <.|nl>  for <kl|nl>
!     if additionally <kl|n'l'> = 0  then  nl -> nl - <nl|n'l'> n'l' 
!----------------------------------------------------------------------
      Use bsr_mat

      Implicit none
      Integer, intent(in) :: kl,nl
      Real(8) :: v(ns), S
      Integer :: i,j, i1,i2, ich
      Real(8), external :: OBS
      Integer, external :: IBORT   

      if(IBORT(kl,nl).eq.0) then
       write(*,'(2a6,10i5)') &
        ebs(kl),ebs(nl),IBORT(kl,nl),IBORT(nl,kl),kl,nl,nbf
       write(*,'(10a6)') ebs(1:nbf)
       Call Stop_mpi(0,0,'get_V: <kl|nl> = 0 ?')
      end if

      v = qbs(:,nl)

      ich = iech(kl)
      i1 = ip_orth_chan(ich-1)+1 
      i2 = ip_orth_chan(ich)
      if(i1.gt.i2) Return

      Do j=i1,i2; i = ip_orth_orb(j)
       S = OBS(i,nl); if(abs(S).lt.eps_ovl) Cycle
       v(:) = v(:) - S*qbs(:,i)
      End do

      End Subroutine get_V


