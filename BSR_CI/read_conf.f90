!======================================================================
      Subroutine Read_conf
!======================================================================
!     Read the configuration list from c-file (unit 'nuc'),
!     define the corresponding angular and configurational symmetries
!     and compare with information from 'INT.BNK' (unit 'nub')
!----------------------------------------------------------------------
      Use bsr_ci
      Use symc_list_LS
      Use symt_list_LS
      Use term_LS

      Implicit none

      Integer :: i, k, nsymc0, nsymt0, it,jt

! ... number of core shells:

      Call R_closed(nuc); kclosd=nclosd    

! ... define symmetries in int_bnk file:

      Call Read_symc_LS(nub); nsymc0=nsymc
      Call Read_symt_LS(nub); nsymt0=nsymt

! ... define configurations in c-file:

      Call Add_conf_LS(nuc,0)

      if(nsymt.gt.nsymt0.or.nsymc.gt.nsymc0) &
       Stop 'angular-coefficent-bank-file is not complete'

!----------------------------------------------------------------------
! ... define IT_state1 and IT_state2:

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
! ... we check only e-e interaction !!! 

      if(allocated(IT_oper)) Deallocate(IT_oper)
      Allocate(IT_oper(1:noper,nsymt*(nsymt+1)/2))

      Call Read_oper_LS(nub)

      k = 0
      Do it = 1,nsymt;  if(IT_state1(it).eq.0) Cycle
       Do jt = 1,it;    if(IT_state1(jt).eq.0) Cycle
        ij = (it-1)*it/2+jt
        if(IT_oper(3,ij).gt.0) Cycle
        k=1; Exit
       End do
        if(k.eq.1) Exit
      End do
      if(k.eq.1) Stop 'angular-bank-file is not complete'

!----------------------------------------------------------------
! ... terms, parity, ne

      k=0; Do i=1,no;  k=k+ln(i)*iq(i);  End do
      parity=(-1)**k

      ne=0; Do i=1,no;  ne=ne+iq(i);  End do

      Call R_label(nuc,1,0)

      Call R_term(nuc)

!----------------------------------------------------------------
! ... print runnung information:

      if(meiv.le.0) meiv = NCFG

      write(iwrite,'(/a)')   'Parameters of c-list:'

      write(iwrite,'(/a,i4,a)') 'NCFG   = ',NCFG, &
                              '  - number of configurations'
      write(iwrite,'( a,i4,a)') 'NTERMS = ',NTERMS, &
                              '  - number of terms'
      write(iwrite,'( a,i4,a)') 'NCLOSD = ',NCLOSD, &
                              '  - number of common closed shells'
      write(iwrite,'( a,i4,a)') 'NE     = ',NE, &
                              '  - number of electrons'
      write(iwrite,'(/a,i4/)') 'one-electron orbitals: nwf    = ',nwf
      write(iwrite,'(15a5)') (ELF(i),i=1,nwf)


      End Subroutine Read_conf

