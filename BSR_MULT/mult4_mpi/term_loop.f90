!======================================================================
      Subroutine TERM_loop
!======================================================================
!     Add  the coeff.s between two determinants (zoef_list), 
!     weighted with term-dependent factors, to the common list (coef_list).
!----------------------------------------------------------------------
      USE mult_par
      USE zoef_list
      USE coef_list,  only: int,idf,ntrm,ctrm

      Implicit none
      Integer :: i,k1,k2,jcase
      Real(8) :: C

      if(nzoef.eq.0) Return

! ... add final coefficients:

      Do i=1,nzoef
       C = Zoef(i); if(abs(C).lt.EPS_c) Cycle
       int = IZ_int(i); idf = IZ_df(i)
       Call Decode_mult(jcase,k1,k2,int)
       Select case (jcase)
        Case(1);  C = C*CNA
        Case(2);  C = C*CNA
        Case(3);  C = C*CNB
       End select

       if(abs(C).lt.EPS_c) Cycle

       Ctrm(1:ntrm) = C*CT_oper(1:ntrm)

       Call Add_coef 
      End do

! ... nulify ZOEF list:

      nzoef = 0

      End Subroutine TERM_loop

