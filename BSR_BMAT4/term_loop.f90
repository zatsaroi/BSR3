!======================================================================
      Subroutine TERM_loop(is,js,ic,jc)
!======================================================================
!     Add to the common list (coef_list) the coeff.s between two 
!     determinants (zoef_list) weighted with term-dependent factors.
!     Check also if the specific coefficient is needed to be added
!     to the bank.
!----------------------------------------------------------------------
      Use bsr_breit
      Use zoef_list, only: nzoef, zoef, iz_int, iz_df 
      Use coef_list, only: ntrm,ctrm,int,idf
      Use term_exp

      Implicit none
      Integer :: i,m,k,k1,k2,it,jt,is,js,ic,jc
      Real(8) :: C

      if(nzoef.le.0) Return

        k = 0
        Do k1=1,kt1; it=IP_kt1(k1) 
        Do k2=1,kt2; jt=IP_kt2(k2)  

         if(is.eq.js.and.it.gt.jt) Cycle
         k = k + 1;  ctrm(k) = C_det1(k1,kd1)*C_det2(k2,kd2)

         if(is.eq.js) then
           if(kd2.ne.kd1) ctrm(k) = ctrm(k) + C_det1(k2,kd1)*C_det2(k1,kd2)        
         elseif(ic.eq.jc) then
           if(it.eq.jt)  ctrm(k) = ctrm(k) + C_det1(k1,kd1)*C_det2(k2,kd2)
         end if

        End do; End do 
 
        CT_oper = 0.d0
        Do i = 1,noper
         if(joper(i).eq.0) Cycle
         Do k = 1,ntrm
          CT_oper(k,i) = Coper(i)*JT_oper(k,i)*ctrm(k)        
         End do
        End do

       
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

      End Subroutine TERM_loop

