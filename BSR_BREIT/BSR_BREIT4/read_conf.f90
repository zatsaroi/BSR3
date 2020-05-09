!======================================================================
      Subroutine Read_conf
!======================================================================
!     Read the configuration list from c-file (unit 'nuc'),
!     define their angular symmetries and compare with existing ones
!     from int_inf file (unit nub).
!     Prepare the angular arrays.
!----------------------------------------------------------------------
      Use bsr_breit
      Use symc_list_LS
      Use symt_list_LS
      Use conf_LS

      Implicit none
      Integer ::  ic
      Integer, external :: Idef_ncfg, Iadd_symc_LS, Iadd_symt_LS
      Character(100) :: AS

! ... remove data from previous calculations:

      Call Alloc_symc_LS(0)
      Call Alloc_symt_LS(0)

! ... define ncfg:

      ncfg = Idef_ncfg(nuc)
      if(ncfg.eq.0) then; icalc=0; Return; end if
      if(allocated(WC)) Deallocate(WC)
      Allocate(WC(ncfg));  WC = 0.d0
      if(allocated(IC_term)) Deallocate(IC_term)
      Allocate(IC_term(ncfg))

! ... read old information, if any:

      if(new.eq.0) then
       Call Read_symc_LS(nub)
       Call Read_symt_LS(nub)
      end if

!---------------------------------------------------------------------
! ... define new symmetries from c-file:

      rewind(nuc)
      ic = 0;  parity=0
      Do 
       read(nuc,'(a)') AS
       read(AS,'(a)') CONFIG
       if(CONFIG(1:1).eq.'*') Exit
       if(CONFIG(5:5).ne.'(') Cycle
       ic = ic + 1
       if(LEN_TRIM(AS).gt.64) read(AS(65:),*) WC(ic)
       read(nuc,'(a)') COUPLE
       Call Decode_c;  Call TEST_C

! ... check angular symmetry:

       iconf = Iadd_symc_LS(Ltotal,Stotal,no,iq,ln)
       iterm = Iadd_symt_LS(iconf,no,LS)
       IT_stat(iterm) = 1  
       IC_term(ic) = iterm

      End do  ! on ic

      WC = abs(WC)

      write(pri,'(/a/)')  'c-fail (cfg.inp) contains: '
      write(pri,'(a,i5)') 'number of atomic states   = ',ncfg
      write(pri,'(a,i5)') 'number of configurations  = ',nsymc
      write(pri,'(a,i5)') 'number of ang. symmetries = ',nsymt

!----------------------------------------------------------------------
! ... sorting the terms according to configurations:

      Call Sort_terms

! ... IT_oper defines operatots and terms already done:

      Call Alloc_it_oper_LS(1)
      if(new.eq.0) Call Load_oper_LS(nub)

      Call Prepare_it_oper(eps_soo)

      icalc=0
      Do ic=1,nsymc; if(IC_need(ic).eq.0) Cycle; icalc=1; Exit; End do

      if(icalc.eq.1) &
       write(*,'(/a,a4/)') 'Need of additional calculations --> yes '
      if(icalc.eq.1) &
       write(pri,'(/a,a4/)') 'Need of additional calculations --> yes '
      if(icalc.eq.0) &
       write(*,'(/a,a4/)') 'Need of additional calculations --> no '
      if(icalc.eq.0) &
       write(pri,'(/a,a4/)') 'Need of additional calculations --> no '


      End Subroutine Read_conf

