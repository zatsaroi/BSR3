!======================================================================
      Subroutine Sym_mat1 (n,m,R,dmn,dmx,dav)
!======================================================================
!     evaluate symmetry of the matrix R(n,n), and then symmetrize it
!----------------------------------------------------------------------
      Implicit none

      Integer, intent(in) :: n,m
      Real(8), intent(out) :: dmn,dmx,dav
      Real(8) :: R(m,m)
      Integer :: i,j
      Real(8) :: dij

      dmx=0.d0; dmn=0.d0; dav=0.d0; if(n.le.1) Return

      Do i = 2,n
       Do j = 1,i-1
        dij = R(i,j)+R(j,i); if(dij.eq.0.d0) Cycle
        dij = dabs ((R(i,j)-R(j,i))/dij)
        dav = dav + dij
        if(dmx.eq.0.d0.or.dij.gt.dmx) dmx = dij
        if(dmn.eq.0.d0.or.dij.lt.dmn) dmn = dij
       End do
      End do

      if(dav.gt.1.d-15) then
       dmx = LOG10(dmx/2.d0)
       dmn = LOG10(dmn/2.d0)
       dav = LOG10(dav/n/(n-1))
      else
       dmx=-15; dmn=-15; dav=-15 
      end if

      Do i = 1,n-1
       Do j = i+1,n
        R(i,j)=(R(i,j)+R(j,i))/2.d0;  R(j,i) = R(i,j)
       End do
      End do

      End Subroutine Sym_mat1


