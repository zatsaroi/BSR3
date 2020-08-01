!======================================================================
      Subroutine Add_res(nu,is,js)
!======================================================================
!     records results from 'coef_list' to unit 'nu'
!----------------------------------------------------------------------
      USE bsr_breit
      USE coef_list
      USE term_exp,  ONLY: jt1,jt2, JP_kt1,JP_kt2
      USE ndef_list, ONLY: IPF

      Implicit none
      Integer, intent(in) :: nu
      Integer :: i,j, it,jt, is,js, k,k1,k2, nd, iihm(ktrm)

      if(ncoef.le.0) Return 
       
! ... convert det.factors from ndef_list to common list: 

      nd = 0; Do i=1,ncoef; nd=nd+idfc(i); End do 
      if(nd.gt.0) then
       Call Ndet_Idet
       Do i=1,ncoef
        if(idfc(i).eq.0) Cycle
        j = idfc(i); idfc(i) = IPF(j)
       End do
      end if

! ... define the term poiners: 

      k = 0
      Do k1=1,jt1; it=JP_kt1(k1) 
      Do k2=1,jt2; jt=JP_kt2(k2)  
       if(is.eq.js.and.it.gt.jt) Cycle
       k = k + 1;  iihm(k) = it*ibc+jt
      End do; End do 
      
      if(k.ne.ktrm) Stop 'Add_res: ij <> ktrm' 

! ... record the coef.s:

      Do j = 1,ncoef
       Do i = 1,ktrm 
        if(abs(coef(i,j)).lt.Eps_C) Cycle
        write(nu) coef(i,j),iihm(i),intc(j),idfc(j)
        nc_new = nc_new + 1
       End do
      End do

      End Subroutine Add_res
