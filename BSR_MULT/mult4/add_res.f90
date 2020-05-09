!======================================================================
      Subroutine Add_res(nu)
!======================================================================
!     records results from 'coef_list' to unit 'nu'
!----------------------------------------------------------------------
      USE mult_par
      USE coef_list
      USE term_exp,  ONLY: kt1,kt2, IP_kt1,IP_kt2
      USE ndef_list, ONLY: IPF 

      Implicit none
      Integer, Intent(in) :: nu
      Integer :: i,j, it,jt, k,k1,k2, nd, ihm(ntrm),jhm(ntrm)

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
      Do k1=1,kt1; it=IP_kt1(k1) 
      Do k2=1,kt2; jt=IP_kt2(k2)  
       if(ic.eq.jc.and.it.gt.jt) Cycle
       k = k + 1;  ihm(k) = it; jhm(k) = jt
      End do; End do 
      
      if(k.ne.ntrm) Stop 'Add_res: ij <> ntrm' 

! ... record the coef.s:

      Do j = 1,ncoef
       Do i = 1,ntrm 
        if(abs(coef(i,j)).lt.Eps_C) Cycle
        write(nu) coef(i,j),ihm(i),jhm(i),intc(j),idfc(j)
       End do
      End do

      End Subroutine Add_res
