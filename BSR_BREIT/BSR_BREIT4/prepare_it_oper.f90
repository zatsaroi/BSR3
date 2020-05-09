!======================================================================
      Subroutine Prepare_it_oper(eps_soo)
!----------------------------------------------------------------------
      Use param_LS
      Use symc_list_LS
      Use symt_list_LS
      Use conf_LS

      Implicit none
      Real(8) :: eps_soo
      Real(8), allocatable :: C_term(:)
      Integer :: i,j, ik,jk, ic,jc, it,jt , ii
      Integer(8) :: ijc
      Integer(8), external :: Def_ij8

!-------------------------------------------------------------------------		 
! ... additional restrictions for spi_orbit interaction:

      if(eps_soo.gt.0.d0) then
       Allocate(C_term(nsymt));  C_term = 0.d0
       Do ic = 1,ncfg; it = IC_term(ic) 
         C_term(it) = max(C_term(it),WC(ic))
       End do
      end if

!-------------------------------------------------------------------------		 
! ... define strong orthogonality:      

      Do ic = 1,nsymc
       Call Get_symc_LS(ic,Ltotal1,Stotal1,no1,nn1,ln1,iq1,kn1) 
      Do jc = 1,ic
       Call Get_symc_LS(jc,Ltotal2,Stotal2,no2,nn2,ln2,iq2,kn2) 
      
       Call Jort_conf(joper)

       Do ik=IC_term1(ic),IC_term2(ic);  it=IT_sort(ik)
       Do jk=IC_term1(jc),IC_term2(jc);  jt=IT_sort(jk)
        ij = DEF_ij8(it,jt)
        Do i=1,noper; if(joper(i).eq.1) IT_oper(i,ij)=1; End do
        if(eps_soo.gt.0.d0) then 
         if(C_term(it)*C_term(jt).lt.eps_soo) then
          Do i = 4,7; if(IT_oper(i,ij).eq.0) IT_oper(i,ij)=-1; End do
         end if
        end if
       End do
       End do
       
      End do
      End do
      
      if(allocated(C_term)) Deallocate(C_term)
      if(allocated(WC)) Deallocate(WC)
      if(allocated(IC_term)) Deallocate(IC_term)
      
! ... define if we need additional calculations

      if(allocated(IT_need)) Deallocate(IT_need)
                             Allocate(IT_need(nsymt))
      if(allocated(JT_need)) Deallocate(JT_need)
                             Allocate(JT_need(ij_oper))
             
      IT_need=0; JT_need=0
      Do it=1,nsymt; Do jt=1,it; ij=DEF_ij8(it,jt)
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
                             ijc = DEF_ij8(nsymc,nsymc)
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
      
      Deallocate(JT_need)
      
      End Subroutine Prepare_it_oper
      
             