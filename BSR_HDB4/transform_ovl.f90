!======================================================================
      Subroutine transform_ovl
!======================================================================
!     transform overlap matrix to a new basis
!----------------------------------------------------------------------
      Use blacs
      Use bsr_hd
      Use target,       only: BFT
      Use channel,      only: elc,iptar
      Use spline_param, only: ns

      Implicit none
      Real(8) :: S, SS
      Integer :: i,j, i1,i2,j1,j2, ic,jc, ii,jj, idim,jdim, info

      Call CPU_time(t0)

      if(io_processor) then

       rewind(nui)
       read(nui) i1,i2
       if(diag_ovl.eq.1) then
        Do 
         read(nui) ic,jc
         if(ic.le.0) Exit
         idim=1; if(ic.le.kch) idim=ns
         jdim=1; if(jc.le.kch) jdim=ns
         read(nui) (S,i=1,idim*jdim)
        End do
       end if
 
      end if

      Call BLACS_BARRIER (ctxt, 'all')

      if(diag_ovl.eq.1) Return

!----------------------------------------------------------------------
! ... initialize array descriptor for the overlap matrix:
! ... np, nq are matrix size on current processor from module bsr_hd

      call descinit(descb,khm,khm,nblock,nblock,rsrc,csrc,ctxt,ld,info)
      call p_error(info,'descb descriptor error in transform_ovl')
      if(allocated(b)) deallocate(b); allocate (b(np,nq));  b=zero

! ... put unit diagonal blocks:

      add = zero; Do i=1,ns; add(i,i)=one; End do
      Do i = 1,kch
       i1=ipsol(i-1)+1; idim=ipsol(i)-ipsol(i-1)
       Call pdgeadd ('notrans', idim,idim,       &
                      one, add,  1, 1, descadd, &
                      zero,  b, i1, i1, descb    )
      Call BLACS_BARRIER (ctxt, 'All')
      End do

      if(kcp.gt.0) adp = zero

!---------------------------------------------------------------------------

      Do 
       if(io_processor) read(nui) ic,jc
       Call br_ipar(ic)
       Call br_ipar(jc)
       Call BLACS_BARRIER (ctxt, 'All')
       if(ic.le.0) Exit

      if(ic.gt.kch.and.jc.gt.kch) then   !  pert-pert

       if(io_processor) then
        read(nui) S
        ic=ic-kch; jc=jc-kch       
        adp(ic,jc)=S 
        if(ic.eq.jc.and.abs(S-one).gt.eps_d) & 
         write(pri,'(a,f10.6,2i5)') &
         'Warning: perturber overlap =',S, ic,jc 
        if(ic.ne.jc.and.abs(S).gt.eps_o) & 
         write(pri,'(a,f10.6,2i5)') &
         'Warning: perturber overlap =',S, ic,jc 
       end if

      elseif(ic.gt.kch) then             !  ch-pert

       if(io_processor) then
        read(nui) w(1:ns)
        jdim = ipsol(jc)-ipsol(jc-1); jj = ipsol(jc-1)
        add = zero; S=zero
        Do j=1,jdim
         add(1,j)=SUM(w(:)*bb(:,j+jj))
         SS = abs(add(1,j)); if(SS.gt.S) S=SS
        End do 
        if(S.gt.eps_o) & 
         write(pri,'(a,f10.3,a,i5,a,a,2x,a)') &
         'Warning: perturber-channel overlap =',S, &
         '  pertuber',ic-kch, &
         '  channel ',elc(jc),BFT(iptar(jc)) 
       end if

       i1=ic-kch+ipsol(kch); j1=ipsol(jc-1)+1
       jdim = ipsol(jc)-ipsol(jc-1)
       Call pdgeadd ('notrans', 1,jdim, one, add,  1, 1, descadd, &
                                       zero,   b, i1,j1, descb  )

      else                              !  ch-ch

       if(io_processor) then

        read(nui) cc(1:ns,1:ns)

        if(ic.eq.jc) Cycle
        if(ic.gt.ktarg.or.jc.gt.ktarg) Cycle

        i1=ipsol(ic-1)+1; i2=ipsol(ic); ii = ipsol(ic-1) 
        j1=ipsol(jc-1)+1; j2=ipsol(jc); jj = ipsol(jc-1) 
        Do i=i1,i2
         Do j=1,ns; w(j)=SUM(bb(:,i)*cc(:,j)); End do 
         Do j=j1,j2; add(i-ii,j-jj)=SUM(w(:)*bb(:,j)); End do 
        End do       

        S = zero 
        idim=ipsol(ic)-ipsol(ic-1) 
        jdim=ipsol(jc)-ipsol(jc-1) 
        Do i=1,idim;  Do j=1,jdim
         SS=abs(add(i,j)); if(SS.gt.S) S=SS
        End do; End do 
        if(S.gt.eps_o) & 
         write(pri,'(a,f10.3,a,2x,a,5x,a,2x,a)') &
         'Warning: channel-channel overlap =',S, &
          elc(ic),BFT(iptar(ic)), & 
          elc(jc),BFT(iptar(jc)) 

       end if

       if(ic.eq.jc) Cycle
       if(ic.gt.ktarg.or.jc.gt.ktarg) Cycle

       i1=ipsol(ic-1)+1; idim=ipsol(ic)-ipsol(ic-1) 
       j1=ipsol(jc-1)+1; jdim=ipsol(jc)-ipsol(jc-1) 
       Call pdgeadd ('notrans',idim,jdim, one,add,   1, 1, descadd,&
                                          zero,  b, i1,j1, descb   )
      end if

       Call BLACS_BARRIER (ctxt, 'All')
      End do

      if(kcp.gt.0) then
       i=ipsol(kch)+1
       Call pdgeadd('notrans',kcp,kcp, one, adp, 1,1, descadp, &
                                       zero,  b, i,i, descb )
      end if

      if(io_processor) then           
       Call CPU_time(t1)
       write (pri,'(/a,T30,f10.2,a)') 'Transform_ovl:,', (t1-t0)/60, ' min.'
       write (*  ,'(/a,T30,f10.2,a)') 'Transform_ovl:,', (t1-t0)/60, ' min.'
      end if

      End Subroutine transform_ovl
