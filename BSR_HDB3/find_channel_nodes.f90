!======================================================================
      Subroutine Find_channel_nodes(ich,is,nodes)
!======================================================================

      Use bsr_hd
      USE spline_param, only: ns

      Implicit none
      Integer :: ich,is,nodes, j,j1,j2
      Real(8) :: S,vv(ns)

! ... find solution in original B-spline basis

      vv = 0.d0
      Do j=ipsol(ich-1)+1,ipsol(ich)
       vv(1:ns) = vv(1:ns) + a(j,is)*bb(1:ns,j)     
      End do

! ... find nodes:

      nodes=0
      S = maxval(abs(vv)) / 100.d0
      j1=2
      Do j=2,ns; if(abs(vv(j)).lt.S) Cycle; j1=j; Exit; End do
      j2=ns-1
      Do j=ns-1,1,-1; if(abs(vv(j)).lt.S) Cycle; j2=j; Exit; End do
      Do j=j1,j2-1
       if(vv(j)*vv(j+1).lt.0.d0) nodes=nodes+1
      End do

      End Subroutine Find_channel_nodes
