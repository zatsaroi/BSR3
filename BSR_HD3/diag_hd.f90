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
!     9. Backtransform eigenvectors to the original problem.
!        Solutions are still in new basis!!!
!----------------------------------------------------------------------
      Use bsr_hd
      Use channel,      only: iptar, lch, elc
      Use target,       only: Etarg, BFT
      Use spline_param, only: ns,ks

      Implicit none
      Integer :: i1,i2,i3, diag_ovl, diag_hm, i,j,ic,jc,j1,j2, ich,jch, &
                 idim,jdim, l,ll,jj,it,jt,nsol,info,idel,n,neig,ibtype, &
                 ipseudo
      Real(8) :: w(ns), bval(ns*kch), aa(ns,ns), cc(ns,ns), vc(ns,ks), &
                 ad(ns,ns,kch), cd(ns,ns,kch), zero=0.d0, one=1.d0,    &
                 S,SS, E_ch
      Real(8), external :: BVMV

      Character(1) :: uplo = 'L', job = 'V', trans = 'T'

!----------------------------------------------------------------------
! ... read diagonal blocks:

      rewind(nui)
      read(nui) i1,i2,i3

! ... overlap matrix:

      Do ich=1,kch; read(nui) cd(:,:,ich); End do 

      Do 
       read(nui) ic,jc
       if(ic.le.0) Exit
       idim=1; if(ic.le.kch) idim=ns
       jdim=1; if(jc.le.kch) jdim=ns
       read(nui) (S,i=1,idim*jdim)
      End do
      diag_ovl = jc

! ... interaction matrix:

      Do ich=1,kch; read(nui) ad(:,:,ich); End do 

      Do 
       read(nui) ic,jc
       if(ic.le.0) Exit
       idim=1; if(ic.le.kch) idim=ns
       jdim=1; if(jc.le.kch) jdim=ns
       read(nui) (S,i=1,idim*jdim)
      End do
      diag_hm = jc

      if(kcp.gt.0) diag_ovl=0

      if(diag_ovl.eq.0) & 
       write(pri,'(/a/)') 'Problem is a GENERALIZED eigenvalue problem' 
      if(diag_ovl.ne.0) & 
       write(pri,'(/a/)') 'Problem is reduced to a STANDARD eigenvalue problem' 

!----------------------------------------------------------------------
! ... find eigensolutions for each channel (new basis):

      if(allocated(bb)   ) deallocate(bb   ); allocate(bb(ns,ns*kch))
      if(allocated(ipsol)) deallocate(ipsol); allocate(ipsol(0:kch) )
      if(allocated(isol) ) deallocate(isol ); allocate(isol(nhm)    )

      ipsol=0; ksol=0; isol=0; bb = 0.d0
      if(debug.gt.0) write(pri,'(/a/)') 'Channel eigenvalues:'
      Do ich = 1,kch;  it=iptar(ich); E_ch = Etarg(it) 
       l=lch(ich); if(l.gt.ks-2) l=ks-2; if(ilzero.eq.0) l=0; i1=l+2
       i2=ns; if(itype.eq.-1) i2=ns-ibzero
       nsol=i2-i1+1 
       aa(1:nsol,1:nsol) = ad(i1:i2,i1:i2,ich)
       cc(1:nsol,1:nsol) = cd(i1:i2,i1:i2,ich)
       Call LAP_DSYGV ('V','U',nsol,ns,aa,cc,w,info) 
       if(info.ne.0) then
        write(pri,*) 'channel diagonalization failed, channel = ',ich
        Stop 'BSR_HD: channel diagonalization failed'
       end if
       if(debug.gt.0) then
        write(pri,'(/a,a,i5,a,i5/)') elc(ich), &
         '   ich =',ich,' nsol =',nsol  
        write(pri,'(5f15.7)') w(1:nsol)
       end if
       jj = (ich-1)*ns + i1 - 1 
       ipseudo = 0
       Do i=1,nsol; j=jj+i
        if(w(i).gt.E_ch) ipseudo = ipseudo + 1
!        if(ipseudo.gt.mpseudo) Exit
        if(Edmin.ne.0.d0.and.w(i).lt.Edmin) Cycle
        if(Edmax.ne.0.d0.and.w(i).gt.Edmax) Cycle
        if(Egap.gt.0.d0.and.abs(w(i)).lt.Egap) Cycle
        ksol=ksol+1
        bval(ksol) = w(i)
        bb(i1:i2,ksol) = aa(1:nsol,i)
       End do
       ipsol(ich) = ksol
      End do
      khm = ksol + kcp;  idel = kch-ksol

      write(pri,'(/a,i6/)') 'new basis: khm =', khm

!----------------------------------------------------------------------
! ... transform the overlap matrix to new basis:

      rewind(nui)
      read(nui) i1,i2,i3
      Do ich=1,kch; read(nui) (S,i=1,ns*ns); End do 

      if(diag_ovl.eq.0) then

      if(allocated(c)) deallocate(c); allocate(c(khm,khm)); c=0.d0
      Do i=1,khm; c(i,i) = 1.d0; End do

      Do 
       read(nui) ic,jc;  if(ic.le.0) Exit

!---------------------------------------------------------------------
       if(ic.gt.kch.and.jc.gt.kch) then  !  pert-pert

        read(nui) S
        c(ic-idel,jc-idel) = S 
         
        ic=ic-kch; jc=jc-kch       
        if(ic.eq.jc.and.abs(S-one).gt.eps_d) & 
         write(pri,'(a,f10.6,2i5)') &
         'Warning: perturber overlap =',S, ic,jc 
        if(ic.ne.jc.and.abs(S).gt.eps_o) & 
         write(pri,'(a,f10.6,2i5)') &
         'Warning: perturber overlap =',S, ic,jc 

!---------------------------------------------------------------------
       elseif(ic.gt.kch) then            !  ch-pert

        read(nui) w(1:ns)
        j1=ipsol(jc-1)+1; j2=ipsol(jc); i=ic-idel
        Do j=j1,j2; c(i,j)=SUM(w(:)*bb(:,j));  End do 

        S = zero 
        Do j=j1,j2; SS=abs(c(i,j)); if(SS.gt.S) S=SS;  End do 
        if(S.gt.eps_o) & 
         write(pri,'(a,f10.3,a,i5,a,a,2x,a)') &
         'Warning: perturber-channel overlap =',S, &
         '  pertuber',ic-kch, &
         '  channel ',elc(jc),BFT(iptar(jc)) 

!---------------------------------------------------------------------
       else                              !  ch-ch

        read(nui) cc(1:ns,1:ns)
        i1=ipsol(ic-1)+1; i2=ipsol(ic)
        j1=ipsol(jc-1)+1; j2=ipsol(jc)
        Do i=i1,i2
         Do j=1,ns; w(j)=SUM(bb(:,i)*cc(:,j)); End do 
         Do j=j1,j2; c(i,j)=SUM(w(:)*bb(:,j)); End do 
        End do       

        S = zero 
        Do i=i1,i2;  Do j=j1,j2
         SS=abs(c(i,j)); if(SS.gt.S) S=SS
        End do; End do 
        if(S.gt.eps_o) & 
         write(pri,'(a,f10.3,a,2x,a,5x,a,2x,a)') &
         'Warning: channel-channel overlap =',S, &
          elc(ic),BFT(iptar(ic)), & 
          elc(jc),BFT(iptar(jc)) 

       end if

       End do

      else

       read(nui) ic,jc

      end if   ! over  diag_ovl
!----------------------------------------------------------------------
! ... transform the interaction matrix to new basis:

      Do ich=1,kch; read(nui) (S,i=1,ns*ns); End do 

      if(allocated(a)) deallocate(a); allocate(a(khm,khm)); a=0.d0
      Do i=1,ksol; a(i,i) = bval(i); End do

      Do 
       read(nui) ic,jc;  if(ic.le.0) Exit

       if(ic.gt.kch.and.jc.gt.kch) then  !  pert-pert

        read(nui) S
        a(ic-idel,jc-idel) = S 

       elseif(ic.gt.kch) then            !  ch-pert

        read(nui) w(1:ns)
        j1=ipsol(jc-1)+1; j2=ipsol(jc); i=ic-idel
        Do j=j1,j2; a(i,j)=SUM(w(:)*bb(:,j));  End do 

       else                              !  ch-ch

        read(nui) cc(1:ns,1:ns)
        i1=ipsol(ic-1)+1; i2=ipsol(ic)
        j1=ipsol(jc-1)+1; j2=ipsol(jc)
        Do i=i1,i2
         Do j=1,ns; w(j)=SUM(bb(:,i)*cc(:,j)); End do 
         Do j=j1,j2; a(i,j)=SUM(w(:)*bb(:,j)); End do 
        End do       

       end if

      End do
!----------------------------------------------------------------------
! ... add mass-velocity corrections:

      if(jmvc.ne.0) then

      Do ich = 1,kch

       it=iptar(ich);  l=lch(ich);  Call mvcv(l,vc)
       i1=ipsol(ich-1)+1; i2=ipsol(ich)

      Do jch = 1,ich

       jt=iptar(jch);  ll=lch(jch);  if(l.ne.ll) Cycle
       j1=ipsol(jch-1)+1; j2=ipsol(jch)

       if(diag_ovl.ne.0.and.ich.ne.jch) Cycle

       if(diag_ovl.eq.0) then
        S = SUM(abs(c(i1:i2,j1:j2)))
        if(S.eq.0.d0) Cycle
       end if
 
       if(debug.gt.1) &
       write(pri,'(/a,a,a/)') &
        'mass-velocity corrections for channels,  ', ELC(ich), ELC(jch)

       Do i=i1,i2; if(i-i1+1.gt.iabs(jmvc)) Cycle; !if(bval(i).gt.Etarg(it)) Cycle 
       Do j=j1,j2; if(j-j1+1.gt.iabs(jmvc)) Cycle; !if(bval(j).gt.Etarg(it)) Cycle 
         if(jmvc.lt.0.and.i.ne.j) Cycle
         S = - 0.5*BVMV(ns,ks,vc,'s',bb(:,i),bb(:,j))
         a(i,j) = a(i,j) + S
         if(debug.gt.1) write(pri,'(2i5,F10.5)') i,j,S
       End do; End do
 
       End do   ! jch
       End do   ! ich

       end if ! over jmvc

!----------------------------------------------------------------------
! ... diagonalize using the LAPACK routines:
!----------------------------------------------------------------------

!----------------------------------------------------------------------
! ... Form a Cholesky factorization of the overlap matrix:
! ... The factorization has the form
! ...   A = U^T * U,  if UPLO = 'U', or
! ...   A = L  * L^T,  if UPLO = 'L'  (other part of matrix is not touched)

      if(diag_ovl.eq.0) then

      CALL DPOTRF(uplo, khm, C, khm, INFO )

      if(info.ne.0) then
       write(pri,*) 'Cholesky factorization failed'
       Stop 'BSR_HD: Cholesky factorization failed'
      end if

! ... Transform problem to standard eigenvalue problem:
! ...    A ->  inv(L)*A*inv(L^T),  uplo = 'L'
! ....   A ->  inv(U^T)*A*inv(U),  uplo = 'U'

      ibtype = 1 !  we consider only  A x = labda B x problem

      CALL DSYGST( ibtype, uplo, khm, A, khm, C, khm, INFO )

      end if  ! over diag_ovl

!----------------------------------------------------------------------
! ... introduction of experimental energies:

      if(iexp.ne.0) then

       write(pri,'(/a/)') 'Experimental energies:'  
       Do ich = 1,kch; it=iptar(ich); S=E_exp(it)-Etarg(it)
        Do i=ipsol(ich-1)+1,ipsol(ich); a(i,i)=a(i,i)+S; End do
        write(pri,'(2i6,2f16.8,f12.5)') ich,it, &
              Etarg(it),E_exp(it),(E_exp(it)-Etarg(it))*27.2112
       End do

      end if
!-----------------------------------------------------------------------
! ... Solve standard eigenvalue problem 
! ... (note:  divide and concer algorith requires much more space)

      write(*,'(/a,i3,a,i5,a,i5)') &
            'BSR_HD:  klsp =',klsp,'   nhm =',nhm,'   khm =',khm

      if(allocated(eval)) deallocate(eval); allocate(eval(khm))

      Call LAP_DSYEV(job,uplo,khm,khm,A,eval,INFO)   

      write(*,'(/a,5f15.5)') 'Eval(1:5) =',eval(1:5)
      write(pri,'(/a,5f15.5)') 'Eval(1:5) =',eval(1:5)
!-----------------------------------------------------------------------
!...  define weights:

      if(itype.ne.0.or.cwt.gt.0.d0)  Call W_out 

!-----------------------------------------------------------------------
! ... Backtransform eigenvectors to the original problem.
! ... For A*x=(lambda)*B*x backtransform eigenvectors:
! ... x = inv(L)'*y or inv(U)*y

      if(diag_ovl.eq.0) then

      IF( uplo.eq.'U' ) trans = 'N'
      IF( uplo.eq.'L' ) trans = 'T'

      n = khm
      neig = khm

      CALL DTRSM('Left', UPLO, TRANS, 'Non-unit', N, NEIG, ONE, &
                  C, khm, A, khm )

      end if ! over diag_ovl

      if(allocated(c)) Deallocate(c)

      End Subroutine Diag_hd


