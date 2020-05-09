!======================================================================
      Subroutine trm_out
!======================================================================
! ... output the R-matrix solutions and 
! ... find the surface amplitudes
!----------------------------------------------------------------------
      Use bsr_hd     
      Use spline_param, only: ns,ks
      Use spline_grid

      Implicit none

      Integer, parameter :: nx = 17
      Real(8), parameter :: hx = 0.08

      Integer :: i,j,i1,i2,j1,j2, is,ich 
      Real(8), allocatable :: WX(:,:,:)
      Real(8) :: x(nx)
      Real(8), external :: bvalu2
      
      i = INDEX(AF_hx,'.'); AF = AF_hx(1:i)//ALSP
      Open(nux,file=AF,form='UNFORMATTED')

      if(allocated(v)) deallocate(v);  Allocate(v(1:nhm))
      if(allocated(WX)) deallocate(WX);  Allocate(WX(kch,nx,khm))

      x(nx) =  RA
      Do i=nx-1,1,-1;  x(i) = x(i+1) - hx;  End do

      Do is=1,khm   ! ,1,-1
       v = 0.d0
       Do ich = 1,kch; i1=(ich-1)*ns+1; i2=ich*ns
        j1 = ipsol(ich-1)+1; j2=ipsol(ich)
        Do j=j1,j2
         v(i1:i2) = v(i1:i2) + a(j,is)*bb(1:ns,j)     
        End do
        WX(ich,nx,is) = v(i2)
        WMAT(ich,is) = v(i2)
        Do i=1,nx-1
         WX(ich,i,is) = bvalu2 (t, v(i1), ns, ks, x(i), 0)
        End do  

if(is.ne.1) Cycle
write(pri,*) 'ich,is = ', ich, is, eval(is)
write(pri,'(2E15.5)') (x(i),WX(ich,i,is),i=1,nx)

       End do
      End do
      
      write(nux) kch,khm
      write(nux) nx, x
      write(nux) WX

      End Subroutine trm_out
