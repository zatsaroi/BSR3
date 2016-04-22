!=======================================================================
      Subroutine Conf_calc
!=======================================================================
!     calculations for given <ic|H|jc> and waiting to other case if any 
!-----------------------------------------------------------------------

      USE bsr_breit

      USE spin_orbitals
      USE term_exp
      USE conf_LS,   only: ne
      USE zoef_list, only: nzoef
      USE coef_list, only: ntrm, ctrm

      Implicit none 

      Integer :: i,m,k,k1,k2,it,jt,ic,jc

      Real(8) :: t1,t2, C_ee,C_so,C_ss, zero=0.d0, one=1.d0  

      Real(8), External :: Z_3j, RRTC

! ... get the job:

    1 Call Get_det_exp(ic,jc)  

      if(ic.le.0) Return

! ...  initial allocations:        

       t1=RRTC()

       Call Alloc_coef(-1)
       Call Alloc_boef(-1)
       Call Alloc_blk (-1)

! ... define the normalization constants for different operators:

       C_so = zero
       if(joper(4)+joper(5).ne.0) &
        C_so = Z_3j(ILT1,-MLt+2,3,1,ILT2,MLt)* &
               Z_3j(IST1,-MSt+2,3,1,IST2,MSt)* &
               (-1)**((ILT1-MLt+IST1-MSt)/2)
       if(C_so.ne.zero) C_so=one/C_so

       C_ss = zero
       if(joper(6).ne.0) &
        C_ss = Z_3j(ILT1,-MLt+2,5,1,ILT2,MLt)* &
               Z_3j(IST1,-MSt+2,5,1,IST2,MSt)* &
               (-1)**((ILT1-MLt+IST1-MSt)/2)
       if(C_ss.ne.zero) C_ss=one/C_ss

       C_ee = one; if(ILT1.ne.ILT2.or.IST1.ne.IST2) C_ee = zero

       if(abs(C_ee) + abs(C_so) + abs(C_ss) .eq. zero) then
        Call Send_res(ic,jc)
        go to 1
       end if

       coper(1) = zero; if(joper(1).gt.0) coper(1) = C_ee
       coper(2) = zero; if(joper(2).gt.0) coper(2) = C_ee
       coper(3) = zero; if(joper(3).gt.0) coper(3) = C_ee
       coper(4) = zero; if(joper(4).gt.0) coper(4) = C_so
       coper(5) = zero; if(joper(5).gt.0) coper(5) = C_so
       coper(6) = zero; if(joper(6).gt.0) coper(6) = C_ss
       coper(7) = zero; if(joper(7).gt.0) coper(7) = C_ee

! ...  calculations:

       Do kd1 = 1,kdt1
 
        Msym1(1:ne)=IM_det1(1:ne,kd1)
        Ssym1(1:ne)=IS_det1(1:ne,kd1)

        Call Det_breit1

       Do kd2 = 1,kdt2
 
        k = 0; m = 0; nzoef = 0
        Do k1=1,kt1; it=IP_kt1(k1) 
        Do k2=1,kt2; jt=IP_kt2(k2)  
         if(ic.eq.jc.and.it.gt.jt) Cycle
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

       End do;  End do

! ... send the results:

      Call Send_res(ic,jc)

      t2=RRTC()
      if(pri.gt.0) write(pri,'(a,2i8,f10.2,a)') &
                   'send coef.s for ic,jc:', ic,jc, (t2-t1)/60, '  min'
      go to 1


      End Subroutine Conf_calc


