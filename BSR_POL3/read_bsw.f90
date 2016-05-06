!======================================================================
      Subroutine Read_bsw(nu)
!======================================================================
!     read B-spline w.f. from bsw-file (unit nu) only those orbitals
!     which are in the list "spline_orbitals"
!----------------------------------------------------------------------
      USE spline_param
      USE spline_atomic
      USE spline_orbitals

      Implicit real(8) (A-H,O-Z)
      Character elw*4

      rewind(nu)
    1 READ(nu,end=2) elw,zw,hw,hmw,rmw,ksw,nsw,mw
      if(zw.ne.z) Stop ' Read_bsw:  z <> zw'
      if(hw.ne.h) Stop ' Read_bsw:  h <> hw'
      if(abs(hmw-hmax).gt.1.d-15) Stop ' Read_bsw:  hmw <> hmax'
      if(rmw.ne.rmax) Stop ' Read_bsw:  rmw <> rmax'
      if(ksw.ne.ks) Stop ' Read_bsw:  ksw <> ks'
      if(nsw.ne.ns) Stop ' Read_bsw:  nsw <> ns'
      Call EL4_nlk(elw,nw,lw,kw)
      ii = Ifind_nlk(nw,lw,kw,0)

      if(ii.gt.0) then
       mbs(ii)  = mw
       nbs(ii)  = nw
       lbs(ii)  = lw
       kbs(ii)  = kw
       ebs(ii)  = elw
       read(nu) pbs(1:mw,ii)
       if(mw.lt.ns) pbs(mw+1:ns,ii) = 0.d0
      else
       read(nu) (x,i=1,mw)
      end if

      go to 1
    2 Close(nu)

      End Subroutine Read_bsw
