!======================================================================
      MODULE inter
!======================================================================
!
!     define the interaction matrix elements under consideration
!
!     noper     - number of different operators
!     ioper(:)  - pointer on the required operators
!     joper(:)  - pointer on the operator under current consideration
!
!     IT_oper(:,:) - pointer on the done calculation for specific
!                    operators and given terms
!     JT_oper(:,:) - pointer on the required operators for given
!                    subset of term between two configurations
!
!     IC_need(:)   - define the need of calc. for two given config.s
!     JC_need(:)   - define the need of calc. for the given config.
! 
!-----------------------------------------------------------------------

      USE configs, only: nsh

      IMPLICIT NONE
      SAVE

      Integer(4), parameter :: noper = 7 

      Integer(1) ioper(noper)/1,1,1,0,0,0,0/, joper(noper)
      Real(8) :: coper(noper)

!     Operator(1)   -   overlaps
!     Operator(2)   -   kinatic energy
!     Operator(3)   -   two-electron electrostatic
!     Operator(4)   -   spin-orbit
!     Operator(5)   -   spin-other-orbit
!     Operator(6)   -   spin-spin
!     Operator(7)   -   orbit-orbit


      Integer(1), Allocatable, Dimension(:,:) :: IT_oper
      Integer(1), Allocatable, Dimension(:,:) :: JT_oper

      Integer(4), Allocatable, Dimension(:) :: IC_need
      Integer(4), Allocatable, Dimension(:) :: JC_need

! ... configurations under consideration:

      Integer(4) :: ic,jc     

! ... number of different (CONFIG or ML,MS) cases :

      Integer(4) :: ic_case   
    
      Real(8), Allocatable, Dimension(:,:) :: CT_oper

! ... list of total LS:

      Integer(4), Allocatable, Dimension(:) :: ILT_ic, IST_ic

      Integer(4) :: nzero=0  ! not used at the moment


      End MODULE inter


!=====================================================================
!  following procedures control of the done calculations and the need 
!  of additional calculations for specific operators and configuration:
!
!  Idef_cme        -  defines required operators for config.s ic,jc
!  IT_CALC (it,jt) -  defines need of calc. for terms it,jt 
!  IC_CALC (ic,jc) -  defines need of calc. for config.s ic,jc
!  DEF_IC (it,jt)  -  record what has been done for config.s ic,jc
!
!===================================================================== 



!=====================================================================
      Integer(4) Function Idef_cme()
!=====================================================================
!
!     define the operators to be considered for given configurations
!     from module 'term_exp': ioper + IT_oper -> joper + JT_oper
!    
!---------------------------------------------------------------------

      USE inter; USE term_exp; USE coef_list, Only: ntrm

      Implicit none
      Integer(4) :: it,jt,ik,jk, i,j,ij, k

       k = 0
       Do ik=1,kt1;  it=IP_kt1(ik)
       Do jk=1,kt2;  jt=IP_kt2(jk)
 
        i=max(it,jt); j=min(it,jt); ij = i*(i-1)/2+j

        Do i=1,noper
         joper(i) = 0
         if(ioper(i).gt.IT_oper(i,ij)) joper(i)=ioper(i)
        End do
    
        if(ic.eq.jc.and.it.gt.jt) Cycle;  k = k + 1

        JT_oper(k,:) = joper(:)

       End do
       End do

       joper = 0
       Do i = 1,noper
        Do it=1,ntrm
         if(JT_oper(it,i).eq.0) Cycle; joper(i) = 1; Exit 
        End do
       End do

       IDEF_cme = 0
       Do i=1,noper
        if(joper(i).eq.0) Cycle; IDEF_cme=1; Exit
       End do

      End function IDEF_cme 


!=====================================================================
      Integer(4) Function IT_calc (it,jt)
!=====================================================================
!
!     determine the need of additional calculation between two terms
!     (then IT_calc=1)
!
!---------------------------------------------------------------------

      USE inter;   USE configs, ONLY: IT_stat

      Implicit none
      Integer(4), Intent(in) :: it,jt
      Integer(4) :: i,j,ij

      IT_calc = 0

! ... are there such states ?                                           
 
      if(IT_stat(it)*IT_stat(jt).eq.0) Return

! ... already done ?

      i=max(it,jt); j=min(it,jt); ij = i*(i-1)/2+j
      Do i=1,noper
       if(ioper(i).gt.IT_oper(i,ij)) IT_calc = 1
      End do                                      

      End Function IT_calc 


!=====================================================================
      Integer(4) Function IC_calc (is,js)
!=====================================================================
!
!     determine the need (then IC_calc = 1) of additional calculation
!     between two configurations is and js
!
!--------------------------------------------------------------------
      
      USE configs, ONLY: IC_term, IP_term, ibc

      Implicit none
      Integer(4), Intent(in) :: is,js
      Integer(4) :: it,jt,ik,jk,it1,it2,jt1,jt2
      Integer(4), External :: IT_calc

      IC_calc=0
 
      it=IC_term(is); it1=mod(it,ibc); it2=it/ibc
      Do ik=it1,it2;  it=IP_term(ik)
 
       jt=IC_term(js); jt1=mod(jt,ibc); jt2=jt/ibc
       Do jk=jt1,jt2;  jt=IP_term(jk)
 
         IC_calc = IT_calc (it,jt); if(IC_calc.eq.1) Return
 
        End do
       End do
 
      End Function IC_calc



!=====================================================================
       Subroutine DEF_IC (is,js)
!=====================================================================
!
!      record what has been done for given config.s is,js 
!
!---------------------------------------------------------------------

       USE inter; Use configs

       Implicit none
       Integer(4), Intent(in) :: is,js
       Integer(4) :: i,j,ij,it,jt,ik,jk,it1,it2,jt1,jt2

       it=IC_term(is); it1=mod(it,ibc); it2=it/ibc
       Do ik=it1,it2;  it=IP_term(ik)
 
        if(IT_stat(it).eq.0) Cycle

       jt=IC_term(js); jt1=mod(jt,ibc); jt2=jt/ibc
       Do jk=jt1,jt2;  jt=IP_term(jk)
 
        if(IT_stat(jt).eq.0) Cycle

        i=max(it,jt); j=min(it,jt); ij = i*(i-1)/2+j
        Do i = 1,noper
         if(ioper(i).gt.IT_oper(i,ij)) IT_oper(i,ij)=ioper(i)
        End do 

       End do; End do

       End Subroutine DEF_IC

