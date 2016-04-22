!======================================================================
      Subroutine Add_res(nu)
!======================================================================
!
!     records results from 'coef_list' to unit 'nu'
!
!----------------------------------------------------------------------

      USE param_br;  USE inter;  USE coef_list;  USE term_exp
      USE ndef_list, ONLY: IPF
      USE configs, ONLY: ibc

      Implicit none
      Integer(4), Intent(in) :: nu
      Integer(4) :: i,j, it,jt, k,k1,k2, nc,nd
 
! ... convert det.factors from ndef_list to common list: 

      nd = 0; Do i=1,ncoef; nd=nd+idfc(i); End do 
      if(nd.gt.0) then

      Call Ndet_Idet

      Do i=1,ncoef
       if(idfc(i).eq.0) Cycle
       j = idfc(i); idfc(i) = IPF(j)
      End do

      End if

! ... define the term poiners: 

      k = 0
      Do k1=1,kt1; it=IP_kt1(k1) 
      Do k2=1,kt2; jt=IP_kt2(k2)  
       if(ic.eq.jc.and.it.gt.jt) Cycle
       k = k + 1;  ijhm(k) = it*ibc+jt
      End do; End do 
      
      if(k.ne.ntrm) Stop 'Add_res: ij <> ntrm' 

! ... record the coef.s:

      nc = 0
      Do i = 1,ntrm 
       Do j = 1,ncoef
        if(abs(coef(i,j)).lt.Eps_C) Cycle
        write(nu) coef(i,j),ijhm(i),intc(j),idfc(j)
        nc=nc+1
       End do
      End do
      ntotc = ntotc + nc

      End Subroutine Add_res
