!======================================================================
      Subroutine R_conf
!======================================================================
!     Read the configuration list from c-file (unit 'nuc'),
!     define there angular symmetries and compare with existing ones.
!     Prepare the angular arrays.
!
!     nub, nui - file units for old and new results, respectively
!----------------------------------------------------------------------
      Use bsr_breit,       only: nuc,nub,pri,new,icalc,ioper,joper  
      Use symc_list_LS
      Use symt_list_LS
      Use conf_LS

      Implicit none
      Integer :: nc, i, ii, ij, ic,jc,ijc, it,jt, ik,jk
      Integer, external :: Idef_ncfg, Iadd_symc_LS, Iadd_symt_LS, DEF_ij

! ... remove data from previous calculations:

      Call Alloc_symc_LS(0)
      Call Alloc_symt_LS(0)

! ... define ncfg:

      nc=Idef_ncfg(nuc)
      if(nc.eq.0) Stop ' ncfg = 0, nothing to do '

! ... read old information, if any:

      if(new.eq.0) then
       Call Read_symc_LS(nub)
       Call Read_symt_LS(nub)
      end if

!---------------------------------------------------------------------
! ... define new symmetries from c-file:

      rewind(nuc); parity=0
      Do 
       read(nuc,'(a)') CONFIG
       if(CONFIG(1:1).eq.'*') Exit
       if(CONFIG(5:5).ne.'(') Cycle
       read(nuc,'(a)') COUPLE
       Call Decode_c;  Call TEST_C

! ... check angular symmetry:

       iconf = Iadd_symc_LS(Ltotal,Stotal,no,iq,ln)
       iterm = Iadd_symt_LS(iconf,no,LS)
       IT_stat(iterm)=1  

      End do  ! on ic

      write(pri,'(/a/)')   ' c-fail (cfg.inp) contains: '
      write(pri,'(a,i5)') ' number of atomic states   = ',NC
      write(pri,'(a,i5)') ' number of configurations  = ',NSYMC
      write(pri,'(a,i5)') ' number of ang. symmetries = ',NSYMT

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

      if(allocated(IT_oper)) Deallocate(IT_oper)
      Allocate(IT_oper(1:noper,nsymt*(nsymt+1)/2))

      IT_oper = 0;   if(new.eq.0) Call Read_oper_LS(nub)

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

		 
! ... define if we need additional calculations

      if(allocated(IT_need)) Deallocate(IT_need)
                             Allocate(IT_need(nsymt))
      if(allocated(JT_need)) Deallocate(JT_need)
                             Allocate(JT_need(nsymt*(nsymt+1)/2))

       IT_need=0; JT_need=0
       Do it=1,nsymt; Do jt=1,it; ij=it*(it-1)/2+jt
        ii = IT_stat(it)*IT_stat(jt)        
	      Do i=1,noper
         if(IT_oper(i,ij).eq.1) Cycle
         if(ioper(i).eq.0.or.ii.eq.0) then
		        IT_oper(i,ij)=-1; Cycle
		       end if
         JT_need(ij)=1; IT_need(it)=1; IT_need(jt)=1
        End do
       End do; End do

      if(allocated(IC_need)) Deallocate(IC_need)
                             Allocate(IC_need(nsymc))
      if(allocated(JC_need)) Deallocate(JC_need)
                             Allocate(JC_need(nsymc*(nsymc+1)/2))

      IC_need = 0; JC_need = 0
      Do ic = 1,nsymc; Do jc = 1,ic; ijc=ic*(ic-1)/2+jc
       Do ik=IC_term1(ic),IC_term2(ic);  it=IT_sort(ik)
       Do jk=IC_term1(jc),IC_term2(jc);  jt=IT_sort(jk)
        ij = DEF_ij(it,jt)
        if(JT_need(ij).eq.0) Cycle
		      JC_need(ijc)=1; IC_need(ic)=1; IC_need(jc)=1
       End do
       End do
      End do; End do

      icalc=0
      Do ic=1,nsymc; if(IC_need(ic).eq.0) Cycle; icalc=1; Exit; End do

      if(icalc.eq.1) &
       write(*,'(/a,a4/)')   ' Need of additional calculations --> yes '
      if(icalc.eq.0) &
       write(*,'(/a,a4/)')   ' Need of additional calculations --> no '

      End Subroutine R_conf

