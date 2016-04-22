!=====================================================================
      Subroutine Add_it_oper(is,js)
!=====================================================================
! ... record what has been done for given case: (is,js) 
!---------------------------------------------------------------------

      USE bsr_breit,    only: JD_oper, noper
      Use symt_list_LS, only: IT_oper 
      USE term_exp,     only: jt1,jt2, JP_kt1,JP_kt2

      Implicit none
      Integer :: it,jt,ik,jk, i,ij, k, is,js
      Integer, external :: DEF_ij

       k = 0
       Do ik=1,jt1;  it=JP_kt1(ik)
       Do jk=1,jt2;  jt=JP_kt2(jk)
        if(is.eq.js.and.it.gt.jt) Cycle;  k = k + 1
        ij=DEF_ij(it,jt)
        Do i=1,noper; if(JD_oper(k,i).gt.0) IT_oper(i,ij)=1; End do
       End do; End do

      End Subroutine Add_it_oper
