!======================================================================
      Subroutine RR_conf(nuc,nub,pri)
!======================================================================
!     Read the configuration list from c-file (unit 'nuc'),
!     define there angular symmetries and compare with existing ones.
!     Prepare the angular arrays.
!
!     nub, nui - file units for old and new results, respectively
!----------------------------------------------------------------------

      USE param_br 

      Use symc_list_LS
      Use symt_list_LS
      USE conf_LS
 
      Implicit none

      Integer, intent(in) :: nuc,nub,pri
      Integer :: nc, i, ii, ij, ijc, it,jt, ik,jk, nsymt0, ic,jc, icalc
      Integer, External :: Idef_ncfg, Iadd_cfg_LS, DEF_ij

! ... remove data from previous calculations:

      Call Alloc_cfg_LS(0)
      Call Alloc_symc_LS(0)
      Call Alloc_symt_LS(0)

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
       Call Decode_c;  Call TEST_C; i=Iadd_cfg_LS() 
      End do  

      write(pri,'(/a/)')   ' c-fail (cfg.inp) contains: '
      write(pri,'(a,i5)') ' number of atomic states   = ',NCFG
      write(pri,'(a,i5)') ' number of configurations  = ',NSYMC
      write(pri,'(a,i5)') ' number of ang. symmetries = ',NSYMT

!----------------------------------------------------------------------
! ... IT_stat defines if input contains given term:
 
      if(Allocated(IT_stat)) Deallocate(IT_stat)
                             Allocate(IT_stat(nsymt))
      IT_stat=0; Do i=1,ncfg; IT_stat(IC_term(i))=1; End do

! ... sorting the terms according to configurations:

      if(Allocated(IT_sort )) Deallocate(IT_sort)
                              Allocate(IT_sort(nsymt))
      if(Allocated(IC_term1)) Deallocate(IC_term1)
                              Allocate(IC_term1(nsymc))
      if(Allocated(IC_term2)) Deallocate(IC_term2)
                              Allocate(IC_term2(nsymc))

      Call SortI(nsymt,IT_conf,IT_sort)

      IC_term1=0
      Do i=1,nsymt
       ic=IT_conf(IT_sort(i))
       if(IC_term1(ic).eq.0) IC_term1(ic)=i
                             IC_term2(ic)=i
      End do

! ... IT_oper defines operatots and terms already done:

      if(allocated(IT_oper)) Deallocate(IT_oper)
      Allocate(IT_oper(1:noper,nsymt*(nsymt+1)/2))

      IT_oper = 0
      read(nub) ii
      if(ii.gt.nsymt*(nsymt+1)/2) Stop 'rr_conf: problems with nsymt'

      Do i=1,ii; read(nub) IT_oper(:,i); End do

! ... define strong orthogonality (repeated at every call?):      

       Do ic = 1,nsymc
        Call Get_symc_LS(ic,Ltotal1,Stotal1,no1,nn1,ln1,iq1,kn1) 
       Do jc = 1,ic
        Call Get_symc_LS(jc,Ltotal2,Stotal2,no2,nn2,ln2,iq2,kn2) 
       
        Call Jort_conf(joper)

        Do ik=IC_term1(ic),IC_term2(ic);  it=IT_sort(ik)
        Do jk=IC_term1(jc),IC_term2(jc);  jt=IT_sort(jk)
          ij = DEF_ij(it,jt)
          Do i=1,noper; if(joper(i).eq.1) IT_oper(i,ij)=1; End do
      	End do
        End do

        End do
        End do
			 
      
      End Subroutine RR_conf

