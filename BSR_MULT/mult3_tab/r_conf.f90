!======================================================================
      Subroutine R_conf(in1,in2,nub)
!======================================================================
!     Read the configuration list from c-files (units 'in1,in2'),
!     define there angular symmetries and compare with ones in
!     data-bnk if any (unit nub)
!     Prepare the angular arrays.
!----------------------------------------------------------------------

      USE symc_list_LS
      USE symt_list_LS
      USE conf_LS

      Implicit none

      Integer, Intent(in) :: in1,in2,nub
      Integer :: nsymc_bnk, nsymt_bnk, ij, it,jt, i,k

! ... initialize arrays:

      Call Alloc_cfg_LS(0)
      Call Alloc_symc_LS(0)
      Call Alloc_symt_LS(0)
      Call Alloc_orb_LS(0)

      Call Read_symc_LS(nub); nsymc_bnk=nsymc
      Call Read_symt_LS(nub); nsymt_bnk=nsymt

! ... define symmetries in c-files:

      Call Add_conf_LS(in1,0);  ncfg1=ncfg
      Call Add_conf_LS(in2,0);  ncfg2=ncfg

      if(nsymt.gt.nsymt_bnk.or.nsymc.gt.nsymc_bnk) &
       Stop 'mult_bnk is incomplete - run mult first'

!----------------------------------------------------------------------
! ... sorting comfigurations according term symmetries:
! ... (define IT_state1 and IT_state2)

      if(allocated(IP_stat  )) Deallocate(IP_stat  )
      if(allocated(IT_state1)) Deallocate(IT_state1)
      if(allocated(IT_state2)) Deallocate(IT_state2)
      Allocate(IP_stat(ncfg),IT_state1(nsymt),IT_state2(nsymt))
      IT_state1=0; IT_state2=0

      Call SORTI (ncfg,IC_term,IP_stat )

      Do i=1,ncfg
       it=IC_term(IP_stat(i))
       if(IT_state1(it).eq.0) IT_state1(it)=i; IT_state2(it)=i 
      End do

!----------------------------------------------------------------------
! ... define if we have all needed angular coefficients:

      if(allocated(IT_done)) Deallocate(IT_done)
      Allocate(IT_done(nsymt*(nsymt+1)/2))

      Call Read_done_LS(nub)

      k = 0
      Do it = 1,nsymt;  if(IT_state1(it).eq.0) Cycle
       Do jt = 1,it;    if(IT_state1(jt).eq.0) Cycle
        ij = (it-1)*it/2+jt
        if(IT_done(ij).gt.0) Cycle
        k=1; Exit
       End do
        if(k.eq.1) Exit
      End do
      if(k.eq.1) Stop 'mult_bnk is incomplete - run mult first'

      End Subroutine R_conf

