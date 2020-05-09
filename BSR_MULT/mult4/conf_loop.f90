!=======================================================================
      Subroutine Conf_loop
!=======================================================================
!     run loop over configurations 
!-----------------------------------------------------------------------

      USE mult_par,      only: nui,nud,ic,jc, kpol,ktype, &
                               JT_oper,CT_oper, qpol,mpol,spol, CNA,CNB
      USE spin_orbitals, only: Lsym1,Msym1,Ssym1,NNsym1, &
                               Lsym2,Msym2,Ssym2,NNsym2
      USE term_exp,      only: kt1,kt2, IP_kt1,IP_kt2, &
                               kd1,kd2, kdt1, kdt2, C_det1, C_det2, &
                               IM_det1,IM_det2, IS_det1,IS_det2, &
                               ILT1,ILT2, IST1,IST2, ic_case

      USE conf_LS,      only: ne
      USE symc_list_LS, only: JC_need, IC_need
      USE symt_list_LS, only: IT_done, ij
      USE coef_list,    only: ntrm
      USE zoef_list,    only: nzoef

      Implicit none 

      Integer :: k1,k2,is,js, it,jt, MLT1,MST1,MLT2,MST2, m,k
      Integer(8), external :: DEF_ij8

      Real(8) :: t0,t1,t2, zero=0.d0,one=1.d0  
      Real(8), external :: Z_3j, Clebsh

      Character(80) :: conf

      Call CPU_time(t0)

!----------------------------------------------------------------------
!                                          cycle 1 over configurations:
      rewind(nud)
      Do is=1,ic_case

       Read(nud) ic,kt1,kdt1,ILT1,IST1,MLT1,MST1

       if(Allocated(IP_kt1)) Deallocate(IP_kt1)
       Allocate(IP_kt1(kt1)); Read(nud) IP_kt1

       if(Allocated(C_det1)) Deallocate(C_det1)
       Allocate(C_det1(kt1,kdt1)); Read(nud) C_det1

       if(Allocated(IM_det1)) Deallocate(IM_det1)
       Allocate(IM_det1(ne,kdt1)); Read(nud) IM_det1

       if(Allocated(IS_det1)) Deallocate(IS_det1)
       Allocate(IS_det1(ne,kdt1)); Read(nud) IS_det1

       read(nud) NNsym1(1:ne)
       read(nud) Lsym1(1:ne)

       if(IC_need(ic).eq.0) Cycle

       Call CPU_time(t1)

!----------------------------------------------------------------------
!                                          cycle 2 over configurations:
      rewind(nud)
      Do js=1,is

       Read(nud) jc,kt2,kdt2,ILT2,IST2,MLT2,MST2

       if(Allocated(IP_kt2)) Deallocate(IP_kt2)
       Allocate(IP_kt2(kt2)); Read(nud) IP_kt2

       if(Allocated(C_det2)) Deallocate(C_det2)
       Allocate(C_det2(kt2,kdt2)); Read(nud) C_det2

       if(Allocated(IM_det2)) Deallocate(IM_det2)
       Allocate(IM_det2(ne,kdt2)); Read(nud) IM_det2

       if(Allocated(IS_det2)) Deallocate(IS_det2)
       Allocate(IS_det2(ne,kdt2)); Read(nud) IS_det2

       read(nud) NNsym2(1:ne)
       read(nud) Lsym2(1:ne)

       if(JC_need(DEF_ij8(ic,jc)).eq.0) Cycle      

!----------------------------------------------------------------------
! ...  define number of terms:

       ntrm = 0
       Do k1=1,kt1; it=IP_kt1(k1) 
       Do k2=1,kt2; jt=IP_kt2(k2)  
        if(ic.eq.jc.and.it.gt.jt) Cycle;  ntrm = ntrm + 1
       End do; End do 

!----------------------------------------------------------------------
! ...  joper and JT_oper:

       if(allocated(JT_oper)) Deallocate(JT_oper,CT_oper)
       Allocate(JT_oper(ntrm),CT_oper(ntrm))

       k = 0; m = 0; JT_oper=0
       Do k1=1,kt1; it=IP_kt1(k1) 
       Do k2=1,kt2; jt=IP_kt2(k2)  
        if(ic.eq.jc.and.it.gt.jt) Cycle;  k=k+1
        ij=DEF_ij8(it,jt) 
        if(IT_done(ij).ne.0) Cycle 
        JT_oper(k) = 1
        m = m + 1
       End do; End do 

       if(m.eq.0) Cycle

!----------------------------------------------------------------------
! ... define the normalization constants for different operators:

      qpol = (MLT1-MLT2)/2; spol = (MST1-MST2)/2; mpol = qpol + spol

      CNA = Z_3j(ILT1,-MLT1+2,2*kpol+1,MLT1-MLT2+1,ILT2,MLT2) &
            * (-1)**((ILT1-MLT1)/2)

      if(IST1.ne.IST2) CNA = 0.d0

      if(ktype.eq.'E') then
       CNB = zero
      else  
       CNB =  Z_3j(ILT1,-MLT1+2,2*kpol-1,MLT1-MLT2+1,ILT2,MLT2) &
              * (-1)**((ILT1-MLT1)/2) &
              * Z_3j(IST1,-MST1+2,3,MST1-MST2+1,IST2,MST2) &
              * (-1)**((IST1-MST1)/2) &
              * CLEBSH(2*kpol-1,MLT1-MLT2+1,3,MST1-MST2+1, &
                2*kpol+1,mpol+mpol+1)
      end if

      if(abs(CNA)+abs(CNB).eq.zero) then; Call DEF_IC; Cycle; end if

      if(CNA.ne.0.d0) CNA = one/CNA  
      if(CNB.ne.0.d0) CNB = one/CNB  

!----------------------------------------------------------------------
! ...  initial allocations for coefficients:

       Call Alloc_coef(-1)

!----------------------------------------------------------------------
! ...  calculations:

       Do kd1 = 1,kdt1

        Msym1(1:ne)=IM_det1(1:ne,kd1)
        Ssym1(1:ne)=IS_det1(1:ne,kd1)

        Call Det_mult1
 
       Do kd2 = 1,kdt2

        k = 0; m = 0; nzoef = 0; CT_oper = 0.d0
        Do k1=1,kt1; it=IP_kt1(k1) 
        Do k2=1,kt2; jt=IP_kt2(k2)  
         if(ic.eq.jc.and.it.gt.jt) Cycle
         k = k + 1
         CT_oper(k) = JT_oper(k)*C_det1(k1,kd1)*C_det2(k2,kd2)
         if(CT_oper(k).ne.0.d0) m=1
        End do; End do 

        if(m.eq.0) Cycle
 
        Msym2(1:ne)=IM_det2(1:ne,kd2)
        Ssym2(1:ne)=IS_det2(1:ne,kd2)

        Call Det_mult2; Call Term_loop 

       End do 
       End do

! ...  store results for given config.s:

       Call Add_res(nui); Call DEF_IC

      End do    ! over jc

      Call CPU_time(t2)

      Call Symc_conf(ic,conf)
      write(*  ,'(a,i6,a,i6,a,i6,a,i6,i10,2F10.2,a,3x,a)') &
        ' ic=',ic,'/',ic_case,'  nterm=',kt1,'  ndet=', kdt1, kt1*kdt1, &
          t2-t1,t2-t0,' sec.',trim(conf)
      End do    ! over ic

      End Subroutine Conf_loop

