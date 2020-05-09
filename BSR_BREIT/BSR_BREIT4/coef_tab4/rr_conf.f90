      Subroutine RR_conf(nuc,nub,pri)
!======================================================================
!     Read the configuration list from c-file (unit 'nuc'),
!     define there angular symmetries and compare with existing ones.
!     Prepare the angular arrays.
!
!     nub, nui - file units for old and new results, respectively
!----------------------------------------------------------------------
      Use param_br 
      Use symc_list_LS
      Use symt_list_LS
      Use conf_LS
      Use orb_LS

      Implicit none
      Integer, intent(in) :: nuc,nub,pri
      Integer :: nc, i, ii, ic,jc,ijc, it,jt, ik,jk
      Integer, external :: Idef_ncfg, Iadd_cfg_LS, DEF_ij
      Integer(8), external :: DEF_ij8

! ... remove data from previous calculations:

      Call Alloc_symc_LS(0)
      Call Alloc_symt_LS(0)

! ... define ncfg:

      nc=Idef_ncfg(nuc)
      if(nc.eq.0) Stop ' ncfg = 0, nothing to do '

! ... read old information, if any:

      Call Read_symc_LS(nub)
      Call Read_symt_LS(nub)

!---------------------------------------------------------------------
! ... define new symmetries from c-file:

      rewind(nuc); parity=0
      Do 
       read(nuc,'(a)') CONFIG
       if(CONFIG(1:1).eq.'*') Exit
       if(CONFIG(5:5).ne.'(') Cycle
       read(nuc,'(a)') COUPLE
       Call Decode_c;  Call TEST_C

       Call Decode_c;  Call TEST_C; i=Iadd_cfg_LS() 

      End do  ! on ic

      write(pri,'(/a/)')  'c-fail (cfg.inp) contains: '
      write(pri,'(a,i5)') 'number of atomic states   = ',nc
      write(pri,'(a,i5)') 'number of configurations  = ',nsymc
      write(pri,'(a,i5)') 'number of ang. symmetries = ',nsymt
      write(pri,*)

!----------------------------------------------------------------------
! ... sorting the terms according to configurations:

      if(allocated(IT_sort )) Deallocate(IT_sort)
                              Allocate(IT_sort(nsymt))
      if(allocated(IC_term1)) Deallocate(IC_term1)
                              Allocate(IC_term1(nsymc))
      if(allocated(IC_term2)) Deallocate(IC_term2)
                              Allocate(IC_term2(nsymc))

      Call SortI(nsymt,IT_conf,IT_sort)

      IC_term1=0
      Do i=1,nsymt
       ic=IT_conf(IT_sort(i))
       if(IC_term1(ic).eq.0) IC_term1(ic)=i
                             IC_term2(ic)=i
      End do

! ... IT_oper defines operatots and terms already done:

      Call alloc_it_oper_LS(1)
      Call Load_oper_LS(nub)
      
      End Subroutine RR_conf

