!======================================================================
      Subroutine rsol_out
!======================================================================
! ... output the R-matrix solutions and 
! ... find the surface amplitudes
!----------------------------------------------------------------------
      Use bsr_hd     
      Use spline_param, only: ns

      Implicit none
      Integer :: i,j,i1,i2,j1,j2, is,ich 

      if(itype.eq.1) then
       i = INDEX(AF_rsol,'.'); AF = AF_rsol(1:i)//ALSP
       Open(nur,file=AF,form='UNFORMATTED')
       write(nur) nhm,khm,kch,kcp,ns
       write(nur) (eval(i),i=1,khm)
      end if

      if(allocated(v)) deallocate(v);  Allocate(v(1:nhm))
      if(allocated(WMAT)) deallocate(WMAT);  Allocate(WMAT(kch,khm))

      Do is=1,khm
      
       v = 0.d0
       Do ich = 1,kch; i1=(ich-1)*ns+1; i2=ich*ns
        j1 = ipsol(ich-1)+1; j2=ipsol(ich)
        Do j=j1,j2
         v(i1:i2) = v(i1:i2) + a(j,is)*bb(1:ns,j)     
        End do
        WMAT(ich,is) = v(i2)
       End do
       if(kcp.gt.0) v(kch*ns+1:nhm)=a(ksol+1:khm,is)
       if(itype.eq.1) write(nur) (v(i),i=1,nhm)

      End do
      
      if(itype.eq.1) close(nur)

      End Subroutine rsol_out
