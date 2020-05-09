!======================================================================
      Subroutine R_conf
!======================================================================
!     Read the configuration list from c-files (units 'in1,in2'),
!     define there angular symmetries and compare with ones in
!     data-bnk if any (unit nub)
!     Prepare the angular arrays.
!----------------------------------------------------------------------
      Use mult_par  

      Use symc_list_LS
      Use symt_list_LS
      USE conf_LS

      Implicit none

      Integer, external :: Idef_ncfg, Iadd_symc_LS, Iadd_symt_LS,  &
                           Iort_dipol, Kort_conf
      Integer :: i, ii, ic,jc,it,jt, ik,jk, nuc, k, nsymc_bnk, nsymt_bnk

      Integer(8) :: ijc
      Integer(8), external :: DEF_ij8

! ... initialize arrays:

      Call Alloc_symc_LS(0)
      Call Alloc_symt_LS(0)

! ... define ncfg:

      Call Check_file(AF1);  open(in1,file=AF1)
      ncfg1=Idef_ncfg(in1)

      ncfg2=ncfg1
      if(AF1.ne.AF2) then
      Call Check_file(AF2);  open(in2,file=AF2)
      ncfg2=Idef_ncfg(in2)
      end if

! ... read old information, if any:

      if(new.eq.0) then
       Call Read_symc_LS(nub)
       Call Read_symt_LS(nub)
      end if
      ii = nsymt*(nsymt+1)/2
      nsymc_bnk=nsymc
      nsymt_bnk=nsymt

!---------------------------------------------------------------------
! ... define new symmetries from c-files:

      nuc=in1
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
      End do  

      write(*,'(a,a,a,3i12)') 'First set : ',trim(AF1),'  ncfg=',ncfg1, nsymc, nsymt


      if(AF1.ne.AF2) then
      nuc=in2
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
       IT_stat(iterm)=2  
      End do  
      end if

      write(*,'(a,a,a,3i12)') 'Second set: ',trim(AF2),'  ncfg=',ncfg2, nsymc, nsymt

!----------------------------------------------------------------------
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

!----------------------------------------------------------------------
! ... IT_done defines terms which are already done:

      if(allocated(IT_done)) Deallocate(IT_done)
      ij_done = DEF_ij8(nsymt,nsymt)
      Allocate(IT_done(ij_done))

      IT_done = 0
      if(new.eq.0) Call Load_done_LS(nub)

! ... define strong orthogonality:      

       Do ic = 1,nsymc
        Call Get_symc_LS(ic,Ltotal1,Stotal1,no1,nn1,ln1,iq1,kn1) 
        Call Conf_parity(1)
       Do jc = 1,ic
        Call Get_symc_LS(jc,Ltotal2,Stotal2,no2,nn2,ln2,iq2,kn2) 
        Call Conf_parity(2)
        k = Kort_conf()
        Do ik=IC_term1(ic),IC_term2(ic);  it=IT_sort(ik)
        Do jk=IC_term1(jc),IC_term2(jc);  jt=IT_sort(jk)
         ij = DEF_ij8(it,jt)
         if(IT_done(ij).eq.1) Cycle
         if(k.gt.1) then; IT_done(ij)=1; Cycle; end if
         ! ... angular orthogonality:

         IT_done(ij) = Iort_dipol(kpol,ktype,LTOTAL1,STOTAL1,PTOTAL1, &
                                             LTOTAL2,STOTAL2,PTOTAL2)
      	End do
        End do
       End do
       End do

! ... define if we need additional calculations

      if(allocated(IT_need)) Deallocate(IT_need)
                             Allocate(IT_need(nsymt))
      if(allocated(JT_need)) Deallocate(JT_need)
                             Allocate(JT_need(ij_done))

       IT_need=0; JT_need=0

       Do it=1,nsymt; Do jt=1,it; ij = DEF_ij8(it,jt) 
        if(IT_done(ij).eq.1) Cycle
        ii = IT_stat(it)*IT_stat(jt)        
        if(ii.eq.0) then; IT_done(ij)=-1; Cycle; end if
        JT_need(ij)=1; IT_need(it)=1; IT_need(jt)=1
       End do; End do

      write(*,'(a)') 'IT_need'

      if(allocated(IC_need)) Deallocate(IC_need)
                             Allocate(IC_need(nsymc))
      if(allocated(JC_need)) Deallocate(JC_need)
                             ijc=DEF_ij8(nsymc,nsymc)
                             Allocate(JC_need(ijc))

      IC_need = 0; JC_need = 0
      Do ic = 1,nsymc; Do jc = 1,ic; ijc=DEF_ij8(ic,jc)
        Do ik=IC_term1(ic),IC_term2(ic);  it=IT_sort(ik)
        Do jk=IC_term1(jc),IC_term2(jc);  jt=IT_sort(jk)
         ij = DEF_ij8(it,jt)
         if(JT_need(ij).eq.0) Cycle
         JC_need(ijc)=1; IC_need(ic)=1; IC_need(jc)=1
        End do
        End do
      End do; End do

      icalc=0
      Do ic=1,nsymc; if(IC_need(ic).eq.0) Cycle; icalc=1; Exit; End do
      if(nsymt.gt.nsymt_bnk) icalc=1

      if(icalc.eq.1) then
       write(*,*)
       write(*,'(a,a,a,i12)') 'First set : ',trim(AF1),'  ncfg=',ncfg1
       write(*,'(a,a,a,i12)') 'Second set: ',trim(AF2),'  ncfg=',ncfg2
       write(*,*)
       write(*,'(a,i8)') 'Number of conf. symmetries = ',NSYMC
       write(*,'(a,i8)') 'Number of ang.  symmetries = ',NSYMT
       write(*,*)
      else
       write(*,'(/a,a/)') trim(AF_b),' is complete'
      end if 

      End Subroutine R_conf

