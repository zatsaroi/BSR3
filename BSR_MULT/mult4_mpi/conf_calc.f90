!=======================================================================
      Subroutine Conf_calc
!=======================================================================
!     calculations for given <ic|H|jc> and waiting to other case if any 
!-----------------------------------------------------------------------
      Use mult_par
      Use spin_orbitals
      Use term_exp
      Use conf_LS,   only: ne
      Use zoef_list, only: nzoef
      Use coef_list, only: ntrm, ctrm

      Implicit none 
      Integer :: i,m,k,k1,k2,it,jt,ic,jc
      Real(8) :: C_ee,C_so,C_ss, zero=0.d0, one=1.d0  
      Real(8), external ::  Z_3j, Clebsh

! ... get the job:

    1 Call Get_det_exp(ic,jc)  

      Call Alloc_coef(-1)

      if(ic.le.0) Return

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

      if(abs(CNA)+abs(CNB).eq.zero) then
        Call Send_res(ic,jc)
        go to 1
      end if

      if(CNA.ne.0.d0) CNA = one/CNA  
      if(CNB.ne.0.d0) CNB = one/CNB  

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

       End do;  End do

! ... send the results:

      Call Send_res(ic,jc)

      go to 1

      End Subroutine Conf_calc


