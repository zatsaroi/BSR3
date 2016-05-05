!======================================================================
      Subroutine Read_bsw (nu,kset)
!======================================================================
!     read radial orbitals in B-spline representation from file 'nu'
!     with shift of set index on kset value
!     if read orbital is not in the list - drop it!
!----------------------------------------------------------------------
      USE spline_atomic;  USE spline_param;  USE spline_orbitals
      
      Implicit none
      Integer, intent(in) :: nu,kset
      Character(4) :: el
      Real(8) :: zw,hw,hmw,rmw
      Integer :: ksw,nsw,mw,n,l,k,i
      Real(8), allocatable :: v(:)
      Integer, external :: Ifind_bsorb

      if(.not.allocated(v)) Allocate(v(ns))

      rewind(nu)
    1 read(nu,end=2) el,zw,hw,hmw,rmw,ksw,nsw,mw
      read(nu) v(1:mw)

      if(zw.ne.z) Stop ' Read_bwfn:  z <> zw'
      if(abs(hw-h).gt.1.d-12) Stop ' Read_bwfn:  h <> hw'
      if(abs(hmw-hmax).gt.1.d-12) Stop ' Read_bwfn:  hmw <> hmax'
      if(abs(rmw-rmax).gt.1.d-12) Stop ' Read_bwfn:  rmw <> rmax'
      if(ksw.ne.ks) Stop ' Read_bwfn:  ksw <> ks'
      if(nsw.ne.ns) Stop ' Read_bwfn:  nsw <> ns'

      Call EL4_nlk(el,n,l,k); k=k+kset; i = Ifind_bsorb(n,l,k)
      if(i.eq.0) go to 1

      pbs(1:mw,i)=v(1:mw); mbs(i)=mw; if(mw.lt.ns) pbs(mw+1:ns,i)=0.d0 

      go to 1
    2 Continue

      End Subroutine Read_bsw


