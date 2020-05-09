!=====================================================================
      Subroutine DEF_IC(ic,jc) 
!=====================================================================
!     define the operators already calculated for given configurations
!---------------------------------------------------------------------
      Use mult_par
      Use conf_LS
      Use symt_list_LS
      Use term_exp

      Implicit none
      Integer :: k1,k2,it,jt,ic,jc 
      Integer(8), external :: DEF_ij8

      Do k1=1,jt1; it=JP_kt1(k1) 
      Do k2=1,jt2; jt=JP_kt2(k2)  
        if(ic.eq.jc.and.it.gt.jt) Cycle
        ij=DEF_ij8(it,jt) 
        if(IT_done(ij).eq.0) IT_done(ij)=1 
      End do; End do 

      End Subroutine DEF_IC

