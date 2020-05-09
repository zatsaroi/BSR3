!=====================================================================
      Integer Function Idef_cme(is,js)
!=====================================================================
!     define the operators to be considered for given symmetries:
!            ioper + IT_oper -> joper + JT_oper
!---------------------------------------------------------------------

      USE bsr_breit,    only: noper, ioper, joper, JT_oper
      USE term_exp,     only: kt1, kt2, IP_kt1, IP_kt2
      USE symt_list_LS, only: IT_oper 

      Implicit none
      Integer :: is,js, it,jt,ik,jk, i,ij, k
      Integer, external :: DEF_ij

       k = 0; JT_oper = 0
       Do ik=1,kt1;  it=IP_kt1(ik)
       Do jk=1,kt2;  jt=IP_kt2(jk)
        if(is.eq.js.and.it.gt.jt) Cycle;  k = k + 1
        ij=DEF_ij(it,jt)
        Do i=1,noper
         if(ioper(i).eq.1.and.IT_oper(i,ij).eq.0) JT_oper(k,i)=1
        End do
       End do
       End do

       Do i=1,noper; joper(i)=maxval(JT_oper(:,i)); End do
       
       IDEF_cme=maxval(joper(:))

      End function IDEF_cme 
