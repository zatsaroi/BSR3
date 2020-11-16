!====================================================================
      Subroutine Boundary_conditions
!====================================================================
! ... apply zero conditions at r=0 and on boundary r=a
! ... iprm(i) = 0 means deletion of B-spline "i"
!--------------------------------------------------------------------
      Use bsr_mchf

      Implicit none
      Integer :: i,j

      iprm=1
      Do i=1,nbf
       j=lbs(i)+1; if(j.gt.ks-1) j=1
       if(ilzero.gt.0) j=ilzero
       iprm(1:j,i)=0 
       j=ns-ibzero+1; iprm(j:ns,i)=0 
      End do

      End Subroutine Boundary_conditions


       