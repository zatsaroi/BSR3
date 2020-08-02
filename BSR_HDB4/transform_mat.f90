!======================================================================
      Subroutine transform_mat
!======================================================================
!     transform interaction matrix to a new basis
!----------------------------------------------------------------------
      Use blacs
      Use bsr_hd
      Use spline_param, only: ns,ks

      Implicit none
      Real(8) :: S
      Integer :: i,j, i1,i2,j1,j2, ic,jc, ii,jj, idim,jdim, ich, info

      Call CPU_time(t0)

!--------------------------------------------------------------------------
! ... initialize array descriptor for the interaction matrices:
! ... np, nq are matrix size on current processor from modul bsr_hd

      call descinit(desca,khm,khm,nblock,nblock,rsrc,csrc,ctxt,ld,info)
      call p_error(info,'desca descriptor error in transform_mat')
      if(allocated(a)) deallocate(a); allocate (a(np,nq));  a=zero

      Do ich = 1,kch
       i1=ipsol(ich-1)+1; idim=ipsol(ich)-ipsol(ich-1)
       if(io_processor) then
        ii=ipsol(ich-1)
        add = zero
        Do i=1,idim; add(i,i)=bval(i+ii); End do
       end if
       Call BLACS_BARRIER (ctxt, 'All')
       Call pdgeadd ('notrans', idim,idim,      &
                      one,  add,  1,  1, descadd, &
                      zero,   a, i1, i1, desca  )
      End do

      if(kcp.gt.0) adp = 0.d0

!---------------------------------------------------------------------------

      Do 
       if(io_processor) read(nui) ic,jc
       Call br_ipar(ic)
       Call br_ipar(jc)
       Call BLACS_BARRIER (ctxt, 'All')
       if(ic.le.0) Exit

      if(ic.gt.kch.and.jc.gt.kch) then  !  pert-pert
       if(io_processor) then
        read(nui) S
        adp(ic-kch,jc-kch)=S 
       end if

      elseif(ic.gt.kch) then            !  ch-pert
       if(io_processor) then
        read(nui) w(1:ns)
        jdim = ipsol(jc)-ipsol(jc-1); jj = ipsol(jc-1)
        add = zero
        Do j=1,jdim; add(1,j)=SUM(w(:)*bb(:,j+jj)); End do 
       end if
       i1=ic-kch+ipsol(kch); j1=ipsol(jc-1)+1
       jdim = ipsol(jc)-ipsol(jc-1)
       Call pdgeadd ('notrans', 1,jdim, one, add,  1, 1, descadd, &
                                       zero,   a, i1,j1, desca  )

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
       end if

       if(ic.eq.jc) Cycle
       if(ic.gt.ktarg.or.jc.gt.ktarg) Cycle

       i1=ipsol(ic-1)+1; idim=ipsol(ic)-ipsol(ic-1) 
       j1=ipsol(jc-1)+1; jdim=ipsol(jc)-ipsol(jc-1) 
       Call pdgeadd ('notrans',idim,jdim, one, add,  1, 1, descadd, &
                                         zero,   a, i1,j1, desca  )
       end if

       Call BLACS_BARRIER (ctxt, 'All')
      End do

      if(kcp.gt.0) then
       i=ipsol(kch)+1
       Call pdgeadd('notrans',kcp,kcp, one, adp, 1,1, descadp, &
                                       zero,  a, i,i, desca  )
      end if

      if(io_processor) then           
       Call CPU_time(t1)
       write (pri,'(/a,T30,f10.2,a)') 'Transform_mat:,', (t1-t0)/60, ' min.'
       write (*  ,'(/a,T30,f10.2,a)') 'Transform_mat:,', (t1-t0)/60, ' min.'
      end if

      End Subroutine transform_mat
