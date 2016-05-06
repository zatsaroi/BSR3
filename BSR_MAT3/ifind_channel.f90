!======================================================================
      Integer Function Ifind_channel(ic)
!======================================================================
!     find channel (or perturber) for given configuration 'ic'
!----------------------------------------------------------------------
      Use target; Use channel

      Implicit none
      Integer, intent(in) :: ic
      Integer :: ich

      if(ic.le.0) Stop 'Ifind_channel: ic <= 0'
      Ifind_channel = 1
      Do ich = 1,nch
       if(ic.gt.ipconf(ich)) Ifind_channel=ich+1
      End do 

      if(Ifind_channel.le.nch) Return

      Do ich = 1,npert
       if(ic.gt.ippert(ich)) Ifind_channel=nch+ich+1   
      End do 

      if(Ifind_channel.gt.nch+npert) then
       write(*,*) 'npert=',npert
       write(*,*)  ippert(0:npert)
       Stop 'Ifind_channel: configuration index, ic, is to big'
      end if

      End Function Ifind_channel
