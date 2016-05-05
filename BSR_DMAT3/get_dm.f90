!======================================================================
      Subroutine  Get_dm(i,j,cl,cv,xl,xv)
!======================================================================
!     provide d-matrix 
!
!     kpol=1:  <i|dv|j> = INT{ P_j(r) [d/dr + fl/r] P_i(r); dr}
!                         fl = (l_j*(l_j+1) - l_i*(l_i+1))/2
!
!     kpol=2:  <i|dv|j> = INT{ P_j(r) [r d/dr + fl] P_i(r); dr}
!                         fl = l_i + 2 for lj-li=2
!                         fl = 1/2     for lj-li=0
!                         fl =-l_i + 1 for lj-li=-2
!
!     <B_i|d/dr|B_j>  -  anti-symmetric
!     <B_i|rd/dr|B_j> -  asymmetric
!----------------------------------------------------------------------
      Use spline_param,    only: ns,ks
      Use spline_orbitals, only: lbs
      Use bsr_dmat,        only: ktype,kpol,rrbs,ddbs,ttbs

      Implicit none
      Integer, intent(in) :: i,j
      Real(8), intent(in) :: cl,cv
      Real(8), intent(out) :: xl(ns,ks+ks-1),xv(ns,ks+ks-1)
      Real(8) :: fl

      xl = cl * rrbs

      if(ktype.eq.'M') then;  xv = cv * rrbs; Return; end if
        
      Call FL_kpol(kpol,lbs(i),lbs(j),fl); fl=fl*cv

      xv = cv*ddbs + fl*ttbs

      End Subroutine Get_dm
