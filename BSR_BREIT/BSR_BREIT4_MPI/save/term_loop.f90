!======================================================================
      Subroutine TERM_loop
!======================================================================
!     Add to the common list (coef_list) the coeff.s between two 
!     determinants (zoef_list) weighted with term-dependent factors.
!     Check also if the specific coefficient is needed to be added
!     to the bank.
!----------------------------------------------------------------------

      USE bsr_breit, only: ct_oper,eps_c,joper
      USE zoef_list, only: nzoef, zoef, iz_int, iz_df 
      USE coef_list, only: ntrm,ctrm,int,idf

      Implicit none
      Integer :: i,m
      Real(8) :: C

      if(nzoef.le.0) Return
       
      Do i=1,nzoef
       C = Zoef(i); if(abs(C).lt.EPS_c) Cycle
       int = IZ_int(i); idf = IZ_df(i);  Call Decode_met(m,int)
       Select case (m)
        Case(3,4); if(joper(7).eq.0) Cycle; Ctrm(1:ntrm)=C*CT_oper(:,7)
        Case(5);   if(joper(3).eq.0) Cycle; Ctrm(1:ntrm)=C*CT_oper(:,3)
        Case(6);   if(joper(2).eq.0) Cycle; Ctrm(1:ntrm)=C*CT_oper(:,2)
        Case(7);   if(joper(4).eq.0) Cycle; Ctrm(1:ntrm)=C*CT_oper(:,4)
        Case(8,9); if(joper(5).eq.0) Cycle; Ctrm(1:ntrm)=C*CT_oper(:,5)
        Case(10);  if(joper(6).eq.0) Cycle; Ctrm(1:ntrm)=C*CT_oper(:,6)
        Case(11);  if(joper(1).eq.0) Cycle; Ctrm(1:ntrm)=C*CT_oper(:,1)
        Case Default;  Stop ' Term_loop: unknown integral '
       End select
       Call Add_coef 
      End do

      nzoef = 0

      End  Subroutine TERM_loop

