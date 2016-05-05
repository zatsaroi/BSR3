!=====================================================================
      Subroutine Add_it_oper(is,js)
!=====================================================================
! ... record what has been done for given case: is,js 
!---------------------------------------------------------------------
      Use bsr_breit,    only: JT_oper, noper
      Use symt_list_LS, only: IT_oper 
      Use term_exp,     only: kt1,kt2, IP_kt1,IP_kt2

      Implicit none
      Integer :: is,js, it,jt,ik,jk, i,ij, k
      Integer, external :: DEF_ij

       k = 0
       Do ik=1,kt1;  it=IP_kt1(ik)
       Do jk=1,kt2;  jt=IP_kt2(jk)
        if(is.eq.js.and.it.gt.jt) Cycle;  k = k + 1
        ij=DEF_ij(it,jt)
        Do i=1,noper; if(JT_oper(k,i).gt.0) IT_oper(i,ij)=1; End do
       End do; End do

      End Subroutine Add_it_oper
