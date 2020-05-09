!======================================================================
      Subroutine Add_res(nu,ic,jc)
!======================================================================
!     records results from 'coef_list' to unit 'nu'
!----------------------------------------------------------------------
      Use mult_par
      Use coef_list
      Use term_exp,  only: jt1,jt2, JP_kt1,JP_kt2
      Use ndef_list, ONLY: IPF 

      Implicit none
      Integer, intent(in) :: nu, ic,jc
      Integer :: i,j, it,jt, k,k1,k2, nd, ihm(ktrm),jhm(ktrm)

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
       if(ic.eq.jc.and.it.gt.jt) Cycle
       k = k + 1;  ihm(k) = it; jhm(k) = jt
      End do; End do 
      
      if(k.ne.ktrm) Stop 'Add_res: ij <> ntrm' 

! ... record the coef.s:

      Do j = 1,ncoef
       Do i = 1,ktrm 
        if(abs(coef(i,j)).lt.Eps_C) Cycle
        write(nu) coef(i,j),ihm(i),jhm(i),intc(j),idfc(j)
       End do
      End do

      End Subroutine Add_res
