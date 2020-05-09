!=====================================================================
      Subroutine DEF_IC 
!=====================================================================
!     define the operators already calculated for given configurations
!---------------------------------------------------------------------
      Use mult_par
      Use conf_LS
      Use symt_list_LS
      Use term_exp

      Implicit none
      Integer :: k1,k2,it,jt 
      Integer, external :: DEF_ij

      Do k1=1,kt1; it=IP_kt1(k1) 
      Do k2=1,kt2; jt=IP_kt2(k2)  
        if(ic.eq.jc.and.it.gt.jt) Cycle
        ij=DEF_ij(it,jt) 
        if(IT_done(ij).eq.0) IT_done(ij)=1 
      End do; End do 

      End Subroutine DEF_IC

