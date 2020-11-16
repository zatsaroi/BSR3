!======================================================================
      Subroutine UPDATE_HS (ns,ks,hm,d,sym)
!======================================================================
!     update channel block (ns,ns)  
!     Possible array representation for input d-matrix:
!     sym = 's'  -->  symmetric banded upper-column storage mode
!     sym = 'n'  -->  non-symmetric band matrix  
!     sym = 'x'  -->  non-symmetric full matrix
!----------------------------------------------------------------------
      Implicit none
      Integer, intent(in) :: ns,ks
      Real(8), intent(in) :: d(ns,*)
      Real(8), intent(inout) :: hm(ns,ns)
      Character(*), intent(in) :: sym
      Integer :: i,j, jp, imin,imax
      Real(8) :: y(ns,ns), yy(ns,ns)
      
      Select case(sym)

      Case('s')

       y = 0.d0
       do jp=1,ks;  do i=1,ns-jp+1;  j=i+jp-1
        y(i,j) = d(i,jp)
        y(j,i) = d(i,jp)
       end do; end do

      Case('n')

       y = 0.d0
       do jp = 1,ks+ks-1
         imin=max0( 1, 1 + ks-jp)
         imax=min0(ns,ns + ks-jp)
         do i = imin,imax;  j=i+jp-ks
          y(i,j) = d(i,jp)
         end do
       end do

      Case('x')

       y(1:ns,1:ns) = d(1:ns,1:ns)

      End Select

      y = y/2.d0; yy = TRANSPOSE(y)
      hm = hm + y
      hm = hm + yy
     
      End Subroutine UPDATE_HS
