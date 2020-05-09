!======================================================================
      Subroutine read_diag
!======================================================================
!     Read diagonal channel blocks and diagonalize them.
!     Transform the overlap and interaction matrixes
!     to new basis.
!
!     Main results are placed in arrays:
!
!     bb(ns,ns*kch) - new bases 
!     bval (ns*kch) - diagonal energies
!     ipsol(0:kch ) - pointer on channel blocks in new basis
!----------------------------------------------------------------------
      Use bsr_hd
      Use channel,      only: lch, elc
      Use spline_param, only: ns,ks

      Implicit none
      Real(8) :: ad(ns,ns,kch),bd(ns,ns,kch),aa(ns,ns),dd(ns,ns),dval(ns)
      Real(8) :: S
      Integer :: i, i1,i2,i3, ic,jc, ich, idim,jdim, l, nsol,info, ii

!----------------------------------------------------------------------
! ... read diagonal blocks:

      rewind(nui)
      read(nui) i1,i2,i3

! ... diagonal overlap blocs:

! ... overlap matrix:

      diag_ovl = 1
      Do 
       read(nui) ic,jc
       if(ic.le.0) Exit
       if(ic.le.kch.and.jc.le.kch.and.ic.eq.jc) then
         read(nui) bd(:,:,ic)
       else
        idim=1; if(ic.le.kch) idim=ns
        jdim=1; if(jc.le.kch) jdim=ns
        read(nui) (S,i=1,idim*jdim)
        diag_ovl = 0
       end if
      End do

! ... diagonal interaction blocks:

      Do 
       read(nui) ic,jc
       if(ic.le.0) Exit
       if(ic.le.kch.and.jc.le.kch.and.ic.eq.jc) then
         read(nui) ad(:,:,ic)
       else
        idim=1; if(ic.le.kch) idim=ns
        jdim=1; if(jc.le.kch) jdim=ns
        read(nui) (S,i=1,idim*jdim)
       end if
      End do

      if(kcp.gt.0) diag_ovl=0

      if(diag_ovl.eq.0) & 
       write(pri,'(/a/)') 'It is a GENERALIZED eigenvalue problem' 
      if(diag_ovl.ne.0) & 
       write(pri,'(/a/)') 'It is reduced to a STANDARD eigenvalue problem' 

!----------------------------------------------------------------------
! ... find eigensolutions for each channel (new basis):

      ipsol=0; ksol=0; bb = 0.d0
      if(debug.gt.0) write(pri,'(/a/)') 'Channel eigenvalues:'
      Do ich = 1,kch
       l=lch(ich); if(l.gt.ks-2) l=ks-2; if(ilzero.eq.0) l=0; i1=l+2
       i2=ns; if(itype.eq.-1) i2=ns-ibzero
       nsol=i2-i1+1 
       aa(1:nsol,1:nsol) = ad(i1:i2,i1:i2,ich)
       dd(1:nsol,1:nsol) = bd(i1:i2,i1:i2,ich)

       Call LAP_DSYGV ('V','U',nsol,ns,aa,dd,dval,info) 

       if(info.ne.0) then
        write(pri,*) 'channel diagonalization failed, channel = ',ich
        write(*  ,*) 'channel diagonalization failed, channel = ',ich
        fail=1; Return
       end if
       if(debug.gt.0) write(pri,'(/a,a,i5,a,i5/)') elc(ich), &
         '   ich =',ich,' nsol =',nsol  
       if(debug.gt.0) write(pri,'(5f15.7)') dval(1:nsol)
       Do i=1,nsol
        if(Edmin.ne.0.d0.and.dval(i).lt.Edmin) Cycle
        if(Edmax.ne.0.d0.and.dval(i).gt.Edmax) Cycle
        if(Egap.gt.0.d0.and.abs(dval(i)).lt.Egap) Cycle
        ksol=ksol+1
        bval(ksol) = dval(i)
        bb(i1:i2,ksol) = aa(1:nsol,i)
       End do
       ipsol(ich) = ksol
      End do
      khm = ksol + kcp 
      write(pri,'(/a,i6/)') 'New basis: khm =', khm

      End Subroutine read_diag
