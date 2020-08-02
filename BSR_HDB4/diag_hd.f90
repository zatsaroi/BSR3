!======================================================================
      Subroutine diag_hd
!======================================================================
!     Diagonalization procedure:
!
!     1. Read diagonal channel blocks and diagonalize them.
!     2. Transform the overlap and interaction matrixes
!        to new basis.
!     3. Add mass-velocity correction if any.
!     4. Form a Cholesky factorization of the overlap matrix.
!     5. Transform problem to standard eigenvalue problem.
!     6. Transform to experimental thresholds energies if any.
!     7. Solve standard eigenvalue problem. 
!     8. Call w_out for weights if required.
!     9, Backtransform eigenvectors to the original problem.
!        Solutions are still in new basis!!!
!----------------------------------------------------------------------
      Use bsr_hd
      Use blacs
      Use spline_param, only: ns,ks

      Implicit none
      Real(8) :: scale=one
      Real(8), allocatable :: work(:)
      Integer :: info,ibtype, lwork, status
      Integer, external :: numroc
      Character(1) :: uplo = 'L', job = 'V', trans = 'T', range = 'A'

!----------------------------------------------------------------------
! ... read diagonal blocks:

      Call CPU_time(t0)

      if(io_processor) Call read_diag

      Call br_ipar(fail); if(fail.ne.0) Return

      if(io_processor) then           
       Call CPU_time(t1)
       write (pri,'(/a,T30,f10.2,a)') 'read_diag:,', (t1-t0)/60, ' min.'
       write (*  ,'(/a,T30,f10.2,a)') 'read_diag:,', (t1-t0)/60, ' min.'
      end if

!----------------------------------------------------------------------
! ... broadcast main parameters:

      Call BLACS_BARRIER (ctxt, 'All')

      Call br_ipar(khm)
      Call br_ipar(kch)
      Call br_ipar(kcp)
      Call br_ipar(ksol)
      Call br_ipar(diag_ovl)
      Call br_ipar(ktarg)

      if(io_processor) then
       Call igebs2d (ictxt, 'all', ' ', kch+1, 1, ipsol, kch+1)
      else
       if(allocated(ipsol)) deallocate(ipsol); allocate(ipsol(0:kch) )
       Call igebr2d (ictxt, 'all', ' ', kch+1, 1, ipsol, kch+1, rsrc, csrc)
      end if    

!--------------------------------------------------------------------
! ... allocate working arrays for distributions:

      Call BLACS_BARRIER (ctxt, 'All')

      call descinit(descadd,ns,ns,ns,ns,rsrc,csrc,ctxt,ns,info)

      call p_error(info,'descadd descriptor error in diag_hd')

      if(allocated(add)) deallocate(add); allocate(add(ns,ns))

      if(kcp.gt.0) then

       call descinit(descadp,kcp,kcp,kcp,kcp,rsrc,csrc,ctxt,kcp,info)
       call p_error(info,'descadp descriptor error in diag_hd')
       if(allocated(adp)) deallocate(adp); allocate(adp(kcp,kcp))

       call descinit(descadv,ns,1,ns,1,rsrc,csrc,ctxt,ns,info)
       call p_error(info,'descadv descriptor error in diag_hd')
       if(allocated(adv)) deallocate(adv); allocate(adv(ns))

      end if

      Call descinit(descv,khm,1,khm,1,rsrc,csrc,ctxt,khm,info)

      call p_error(info,'descv descriptor error in diag_hd')

      if(allocated(v)) deallocate(v); allocate(v(khm))     

! ... local dimensions:

      np = numroc (khm, nblock, myrow, rsrc, p)
      nq = numroc (khm, nblock, mycol, csrc, q)
      ld = MAX(1, np)

!----------------------------------------------------------------------
! ... transform the overlap and interaction matrixes to new basis:

      Call Transform_ovl
        
      Call BLACS_BARRIER (ctxt, 'all')

      Call Transform_mat

!----------------------------------------------------------------------
! ... add mass-velocity corrections:

      if(jmvc.ne.0) Call Add_mvc

!----------------------------------------------------------------------
! ... diagonalize using the SCALAPACK routines:
!----------------------------------------------------------------------
! ... default parameters:
 
      ibtype =  1    ! A * x = lambda * B * x generalized eigenproblem
      range  = 'A'   ! compute all eigenvalues
      job    = 'V'   ! compute both eigenvalues and eigenvectors
      uplo   = 'L'   ! use lower triangles of A and B matrices
      trans  = 'T'   ! 

!----------------------------------------------------------------------
! ... Form a Cholesky factorization of the overlap matrix:
! ... The factorization has the form
! ...   A = U^T * U,    if UPLO = 'U'
! ...   A = L   * L^T,  if UPLO = 'L'
! ... (other part of matrix is not touched)

      if(diag_ovl.eq.0) then

      Call CPU_time(t0)

      call PDPOTRF (uplo, khm, b, 1, 1, descb, info)

      Call BLACS_BARRIER (ctxt, 'all')

      call p_error (info, 'pdpotrf: Cholesky error')

!----------------------------------------------------------------------
! ... Transform problem to standard eigenvalue problem:
! ...    A ->  inv(L)*A*inv(L^T),  uplo = 'L'
! ....   A ->  inv(U^T)*A*inv(U),  uplo = 'U'

      lwork = 2 
      allocate (work(lwork), stat=status)
      call p_error (status, 'error memory allocation of work array')

      lwork = -1

      Call PDSYNGST (ibtype, uplo,  khm,a,1,1,desca, b,1,1,descb, &
                   &  scale, work, lwork, info)

      lwork=work(1)+1; deallocate (work); allocate(work(lwork))

      Call PDSYNGST (ibtype, uplo,  khm,a,1,1,desca, b,1,1,descb, &
                   &  scale, work, lwork, info)

      call p_error (info, 'pdsyngst error')

      if(io_processor) then           
       Call CPU_time(t1)
       write (pri,'(/a,T30,f10.2,a)') 'Cholesky factorization:,', (t1-t0)/60, ' min.'
       write (*  ,'(/a,T30,f10.2,a)') 'Cholesky factorization:,', (t1-t0)/60, ' min.'
      end if

      end if  ! over diag_ovl

!----------------------------------------------------------------------
! ... introduction of experimental energies:

      if(iexp.ne.0) Call Add_exp

!-----------------------------------------------------------------------
! ... Solve standard eigenvalue problem 
! ... (note:  divide and concer algorith requires much more space)

      Call CPU_time(t0)

      call descinit(descz,khm,khm,nblock,nblock,rsrc,csrc,ctxt,ld,info)
      call p_error(info,'descz descriptor error')
      if(allocated(z)) deallocate(z); allocate (z(np,nq));  z=zero

      if(allocated(eval)) Deallocate(eval); Allocate(eval(khm))

      if(.not.allocated(work)) Allocate(work(10))

      lwork = -1
      Call PDSYEV (job, uplo, khm,a,1,1,desca, eval, z,1,1,descz, &
                    work, lwork,  info)

      lwork = work(1) + 1;  Deallocate(work); Allocate(work(lwork))

      if(io_processor) &
      write(*,'(/a,i6,a,i6,a,i10)') &
        'diagonalization:  nhm =',nhm,'  khm =',khm,'  lwork =',lwork

      if(io_processor) &
      write(pri,'(/a,i6,a,i6,a,i10)') &
        'diagonalization:  nhm =',nhm,'  khm =',khm,'  lwork =',lwork

      call blacs_barrier (ctxt, 'All')

      Call PDSYEV (job, uplo, khm,a,1,1,desca, eval, z,1,1,descz, &
                    work, lwork, info)
      call p_error (info, 'pdsyev error')

      if(io_processor) then           
       Call CPU_time(t1)
       write (pri,'(/a,T30,f10.2,a)') 'diagonalization:,', (t1-t0)/60, ' min.'
       write (*  ,'(/a,T30,f10.2,a)') 'diagonalization:,', (t1-t0)/60, ' min.'
      end if

      deallocate (work)

      if(io_processor) &
       write(*  ,'(/a,5f15.5)') 'Eval(1:5) =',eval(1:5)
      if(io_processor) &
       write(pri,'(/a,5f15.5)') 'Eval(1:5) =',eval(1:5)

!-----------------------------------------------------------------------
!...  define weights:

      if(cwt.gt.0.d0)  Call W_out 

      if(fail.ne.0) Return

!-----------------------------------------------------------------------
! ... Backtransform eigenvectors to the original problem.
! ... For A*x=(lambda)*B*x backtransform eigenvectors:
! ... x = inv(L)'*y or inv(U)*y

      if(diag_ovl.eq.0) then

       Call CPU_time(t0)

       Call PDTRSM ('Left', uplo, trans, 'Non-unit', khm, khm, one, &
                    b, 1,1, descb, z, 1,1, descz)
       if (scale /= one) call DSCAL(khm, scale, eval, 1)
 
       if(io_processor) then           
        Call CPU_time(t1)
        write (pri,'(/a,T30,f10.2,a)') 'Back_transform:,', (t1-t0)/60, ' min.'
        write (*  ,'(/a,T30,f10.2,a)') 'Back_transform:,', (t1-t0)/60, ' min.'
       end if

      end if ! over diag_ovl

      Call BLACS_BARRIER (ctxt, 'all')

      END SUBROUTINE Diag_hd





