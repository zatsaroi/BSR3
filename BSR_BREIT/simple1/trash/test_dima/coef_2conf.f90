!=======================================================================
      Subroutine coef_2conf
!=======================================================================
!     run loop over configurations 
!-----------------------------------------------------------------------
      USE bsr_breit,     only: pri,nui,nud, &
                               noper,joper,coper,JT_oper,CT_oper
      USE spin_orbitals, only: Lsym1,Msym1,Ssym1,NNsym1, &
                               Lsym2,Msym2,Ssym2,NNsym2
      USE term_exp,      only: kt1,kt2, IP_kt1,IP_kt2, &
                               kd1,kd2, kdt1, kdt2, C_det1, C_det2, &
                               IM_det1,IM_det2, IS_det1,IS_det2, &
                               ILT1,ILT2, IST1,IST2, MLT,MST, ic_case

      USE conf_LS,      only: ne
      USE symc_list_LS, only: JC_need, IC_need, nsymc
      USE coef_list,    only: ntrm,ctrm, ncoef
      USE zoef_list,    only: nzoef

      Implicit none 

      Integer :: k1,k2, it,jt, ij, MLT2,MST2, i,m,k, is,js,ic,jc
      Integer, External :: IDEF_cme, DEF_ij

      Real(8) :: t1,t2, C_ee,C_so,C_ss, zero=0.d0, one=1.d0  
      Real(8), External :: RRTC, Z_3j

      ncoef = 0

! ... checking symmetries:

      if(Ltotal1.ne.Ltotal2)  Return
      if(Stotal1.ne.Stotal2)  Return
      if(MLT2.ne.MLT1.or.MST2.ne.MST1) Return

! ... initial allocations:

      Call Alloc_coef(-1)
      Call Alloc_boef(-1)
      Call Alloc_blk (-1)

! ... calculations:

      Do kd1 = 1,kdt1

       Msym1(1:ne)=MLdet1(1:ne,kd1)
       Ssym1(1:ne)=MSdet1(1:ne,kd1)

       Call Det_breit1
 
      Do kd2 = 1,kdt2

        k = 0; m = 0; nzoef = 0
        Do k1=1,kt1; it=IP_kt1(k1) 
        Do k2=1,kt2; jt=IP_kt2(k2)  
         if(is.eq.js.and.it.gt.jt) Cycle
         k = k + 1;  ctrm(k) = C_det1(k1,kd1)*C_det2(k2,kd2)
         if(ctrm(k).ne.0.d0) m=1
        End do; End do 
        if(m.eq.0) Cycle
 
        CT_oper = 0.d0; m=0
        Do i = 1,noper
         if(joper(i).eq.0) Cycle
         Do k = 1,ntrm
          CT_oper(k,i) = Coper(i)*JT_oper(k,i)*ctrm(k)        
          if(CT_oper(k,i).ne.0.d0) m=1
         End do
        End do
        if(m.eq.0) Cycle

        Msym2(1:ne)=IM_det2(1:ne,kd2)
        Ssym2(1:ne)=IS_det2(1:ne,kd2)

        Call Det_breit2; Call Term_loop 

       End do 
       End do

! ...  store results for given config.s:

       Call Add_res(nui,is,js)

      End Subroutine coef_2conf

