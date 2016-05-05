!======================================================================
      Integer Function Ifind_pert(ilsp,ic)
!======================================================================
!     find perturber for given configuration 'ic'
!----------------------------------------------------------------------
      Use channels

      Implicit none
      Integer, intent(in) :: ilsp,ic
      Integer :: i, ii, jc

      Ifind_pert=ic;    if(ilsp.le.0) Return 

      Ifind_pert=0;     if(ic.le.ipconf(ilsp,nch(ilsp))) Return

      jc=ic-ipconf(ilsp,nch(ilsp)); ii = ipert(ilsp)

      Ifind_pert = -1 
      Do i = 1,npert(ilsp)
       if(jc.gt.ippert(ii+i)) Cycle
       Ifind_pert = i; Exit
      End do 

      if(Ifind_pert.eq.-1) Stop 'Ifind_pert: ic is to big'

      End Function Ifind_pert
