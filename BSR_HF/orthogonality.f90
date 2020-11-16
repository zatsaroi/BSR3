!==================================================================
      Subroutine orthogonality(hfm,v)
!==================================================================
! ... Redefine the matrix to project out the vector v so that
! ... eigenvalues of  (H' - lamda S)x = 0 are orthogonal to v.
!              H' =>  (1 - Bvv^T ) H (1 - vv^TB) 
!------------------------------------------------------------------
      Use bsr_hf

      Implicit none
      Real(8), intent(inout) :: hfm(ns,ns)
      Real(8), intent(in) :: v(ns)
      Real(8) :: x(ns,ns), y(ns,ns), zz(ns,ns), xx(ns,ns), yy(ns,ns) 
      Integer :: i,j

      Do i=1,ns; Do j=1,ns
        zz(i,j)=v(i)*v(j)
      End do; End do        

      x = matmul (BB,zz)
      y = matmul (zz,BB)

      xx  = matmul (x,hfm)
      hfm = hfm - xx
      yy  = matmul (hfm,y)
      hfm = hfm - yy

      End Subroutine orthogonality
