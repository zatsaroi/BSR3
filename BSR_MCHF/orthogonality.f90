!==================================================================
      Subroutine apply_orthogonality(hfm,v)
!==================================================================
! ... Redefine the matrix to project out the vector v so that
! ... eigenvalues of  (A' - lamda S)x = 0 are orthogonal to v.
!                H -->  (1 - Bvv') H (1 - vv'B) 
!                H -->  (1 -  wv') H (1 -  vw') 
!                H -->  (1 -    c) H (1 -   c') 
!------------------------------------------------------------------
      Use bsr_mchf

      Implicit none
      Real(8), intent(inout) :: hfm(ns,ns)
      Real(8), intent(in) :: v(ns)
      Real(8) :: w(ns), c(ns,ns)
      Integer :: i,j

      w = matmul(BB,v)
       
      Do i=1,ns
       Do j=1,ns
        c(i,j)=-w(i)*v(j); if(i.eq.j) c(i,j)=c(i,j) + 1.d0
       End do
      End do        

      hfm = MATMUL(c,hfm);  c = TRANSPOSE(c);  hfm = MATMUL(hfm,c)

      End Subroutine apply_orthogonality
