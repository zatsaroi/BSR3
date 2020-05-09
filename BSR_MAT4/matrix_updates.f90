!======================================================================
      Subroutine allocate_matrix(m)
!======================================================================
!     allocate arrays connected with the Hamiltonian matrix
!----------------------------------------------------------------------
      Use bsr_mat

      Implicit none
      Integer :: ich,jch, i,k,m, mycase

      m = 0  ! memory 
      if(allocated(hcc)) Deallocate(hcc,ACF,htarg,otarg,x,icc) 
      if(allocated(hcb)) Deallocate(hcb,hbb,icb,ibb) 

      if(nch.eq.0.or.ns.eq.0) Return

      Allocate(icc(nch,nch));  icc = 0; iicc = 0
      if(allocated(my_channel)) Deallocate(my_channel)
      Allocate(my_channel(nch+npert));  my_channel = 0
      m = m + nch*nch + nch + npert

      i=0; k=0
      Do ich=1,nch
      Do jch=1,ich
       if(nprocs.gt.1) then
        k=k+1; if(k.gt.nprocs-1) k=1  
       end if
!       if(myid.eq.0) then; icc(ich,jch) = k; icc(jch,ich) = k; end if  ???
       if(myid.ne.k) Cycle
       i=i+1; icc(ich,jch) = i; icc(jch,ich) = i
       my_channel(ich) = 1
       my_channel(jch) = 1
      End do; End do
      iicc=i
      ich = (nch+1)*nch/2
      Allocate(hcc(ns,ns,iicc),acf(nch,nch,0:mk), &
               htarg(ich),otarg(ich),x(ns,ns))
      m = m + 2*(ns*ns*iicc + nch*nch*(mk+1) + 2*ich + ns*ns)
      hcc=0.d0; acf=0.d0; htarg=0.d0; otarg=0.d0

      if(npert.gt.0) then  

      Allocate(icb(nch,npert));  icb = 0
      m = m + 2 * nch * nch
      i=0; k=0
      Do ich=1,nch
      Do jch=1,npert
       if(nprocs.gt.1) then
        k=k+1; if(k.gt.nprocs-1) k=1
       end if
!       if(myid.eq.0) then; icb(ich,jch) = k; end if ???
       if(myid.ne.k) Cycle 
       i=i+1; icb(ich,jch) = i
       my_channel(ich) = 1
       my_channel(nch+jch) = 1
      End do; End do
      iicb = i
      Allocate(hcb(ns,iicb));  hcb=0.d0
      m = m + 2 * ns * iicb

      Allocate(ibb(npert,npert));  ibb = 0
      m = m +  npert * npert
      i=0; k=0
      Do ich=1,npert
      Do jch=1,ich
       if(nprocs.gt.1) then
        k=k+1; if(k.gt.nprocs-1) k=1
       end if
!       if(myid.eq.0) then; ibb(ich,jch) = k; ibb(jch,ich) = k; end if
       if(myid.ne.k) Cycle 
       i=i+1; ibb(ich,jch) = i; ibb(jch,ich) = i
       my_channel(nch+ich) = 1
       my_channel(nch+jch) = 1
      End do; End do
      iibb = i
      Allocate(hbb(iibb));  hbb=0.d0
      m = m + iibb

      end if   ! over npert > 0

      if(allocated(imycase)) Deallocate(imycase)
      Allocate( imycase(1:nch+npert,1:nch+npert) )
      m = m + (nch+npert)*(nch+npert)
      imycase = 0

      Do ich=1,nch+npert; Do jch=1,nch+npert
       imycase(ich,jch) = mycase(ich,jch)
      End do; End do


      End Subroutine allocate_matrix


!======================================================================
     Integer Function mycase(ich,jch)
!======================================================================
! ... Check if (ich,jch) block belong to the current processor
!----------------------------------------------------------------------
     Use bsr_mat

     Integer, intent(in) :: ich,jch
     mycase = 0
     if(ich.le.nch.and.jch.le.nch) then
      mycase = icc(ich,jch)
     elseif(ich.le.nch.and.jch.gt.nch) then
      mycase = icb(ich,jch-nch)
     elseif(ich.gt.nch.and.jch.le.nch) then
      mycase = icb(jch,ich-nch)
     elseif(ich.gt.nch.and.jch.gt.nch) then
      mycase = ibb(ich-nch,jch-nch)
     else
      mycase = 0
     end if

     End Function mycase


!----------------------------------------------------------------------
!     next routines updates the matrixes;
!     we have coefficients only for the half of the matrix,  
!     so then we need   h -->  h + h*   for diagonal blocks
!----------------------------------------------------------------------

!======================================================================
      Subroutine UPDATE_HX(ich,jch,nsb,ksb,d,sym)
!======================================================================
!     update channel block 
!
!     sym = 's'  -->  symmetric banded upper-column storage mode
!     sym = 'n'  -->  non-symmetric band matrix  
!     sym = 'x'  -->  non-symmetric full matrix
!----------------------------------------------------------------------
      Use bsr_mat

      Implicit none
      Integer, intent(in) :: ich,jch,nsb,ksb
      Real(8), intent(in) :: d(nsb,*)
      Character(1), intent(in) :: sym
      Integer :: i,j,ij, ith, imin,imax

      if(nsb.ne.ns) Stop 'UPDATE_HX: ksb /= ks'
      if(ksb.ne.ks) Stop 'UPDATE_HX: ksb /= ks'

      ij = imycase(ich,jch);  if(ij.le.0) Return

      x = 0.d0
      Select case(sym)

      Case('s')

       Do ith=1,ks;  Do i=1,ns-ith+1;  j=i+ith-1
        x(i,j) = d(i,ith)
        x(j,i) = d(i,ith)
       End do; End do

      Case('n')

       Do ith = 1,ks+ks-1
        imin=max0( 1, 1 + ks-ith)
        imax=min0(ns,ns + ks-ith)
        Do i = imin,imax
         j = i + ith - ks
         x(i,j) = d(i,ith)
        End do
       End do

      Case('x')

       x(1:ns,1:ns) = d(1:ns,1:ns)

      End Select

      if(ich.ge.jch) then
       hcc(:,:,ij) = hcc(:,:,ij) + x
      else
       hcc(:,:,ij) = hcc(:,:,ij) + TRANSPOSE(x) 
      end if

      End Subroutine UPDATE_HX


!======================================================================
      Subroutine UPDATE_HL(ich,jch,nsb,ksb,d,c)
!======================================================================
!     update symmetric banded matrix in lower-column storage 
!----------------------------------------------------------------------
      Use bsr_mat

      Implicit none
      Integer, intent(in) :: ich,jch,nsb,ksb
      Real(8), intent(in) :: c
      Real(8), intent(in) :: d(nsb,ksb)
      Integer :: i,j,ij,ith
      Real(8) :: S

      ij = imycase(ich,jch)
      if(ij.le.0) Return

       Do ith = 1,ks
        Do i = ks+1-ith,ns
         j = i + ith - ks
         S = c*d(i,ith)
         hcc(i,j,ij) = hcc(i,j,ij) + S
         if(i.ne.j) hcc(j,i,ij) = hcc(j,i,ij) + S
        End do
       End do

      End Subroutine UPDATE_HL


!======================================================================
     Subroutine UPDATE_HB(ic,jc,c)
!======================================================================
!    update scalar in the bound-bound block
!----------------------------------------------------------------------
     Use bsr_mat

     Implicit none
     Integer, intent(in) :: ic,jc
     Real(8), intent(in) :: c
     Integer :: ij
      
      ij = ibb(ic,jc)
      if(ij.le.0) Return
      hbb(ij) = hbb(ij) + c

     End Subroutine UPDATE_HB


!======================================================================
     Subroutine UPDATE_HV(ich,ic,nsb,v,c)
!======================================================================
!    update vector in the channel-bound block
!----------------------------------------------------------------------
     Use bsr_mat

     Implicit none
     Integer, intent(in) :: ich,ic,nsb
     Real(8), intent(in) :: c
     Real(8), intent(in) :: v(nsb)
     Integer :: ij

     ij = icb(ich,ic);  if(ij.le.0) Return

     hcb(1:nsb,ij) = hcb(1:nsb,ij) + C*v(1:nsb)

     End Subroutine UPDATE_HV


!======================================================================
     Subroutine UPDATE_HW(ich,jch,nsb,v,w)
!======================================================================
!    update channel block with by v*w
!----------------------------------------------------------------------
     Use bsr_mat

     Implicit none
     Integer, intent(in) :: ich,jch,nsb
     Real(8), intent(in) :: v(*),w(*)
     Integer :: i,j,ij

     if(nsb.ne.ns) Stop 'UPDATE_HW: nsb /= ns'

     ij = imycase(ich,jch);  if(ij.le.0) Return
 
     Do i=1,ns;  Do j=1,ns;  x(i,j)=v(i)*w(j);  End do;  End do

     if(ich.ge.jch) then
      hcc(:,:,ij) = hcc(:,:,ij) + x
     else
      hcc(:,:,ij) = hcc(:,:,ij) + TRANSPOSE(x) 
     end if

     End Subroutine UPDATE_HW


!=======================================================================
     Subroutine UPDATE_CF(k,ich,jch,C)
!=======================================================================
!    Update asymptotic coefficients
!-----------------------------------------------------------------------
     Use bsr_mat

     Implicit none
     Integer, intent(in) :: k, ich,jch
     Real(8), intent(in) :: C

     if(ich.lt.1.or.ich.gt.nch) &
       Call Stop_mpi(0,ich,'UPDATE_CF: ich index out of range')
     if(jch.lt.1.or.jch.gt.nch) &
       Call Stop_mpi(0,jch,'UPDATE_CF: jch index out of range')
     if(k.lt.0.or.k.gt.mk) &
       Call Stop_mpi(0,k,'UPDATE_CF: multipole index out of range')

     acf(ich,jch,k) = acf(ich,jch,k) + C

     End Subroutine UPDATE_CF


!======================================================================
     Subroutine Target_h(ich,jch,C,CC)
!=====================================================================
!    update target interaction and overlap matrixes, by considering 
!    the terms with structure <kl|k'l> <target|H|target'>
!
!    It is not pure target states, but the basis states before |kl>,
!    i.e. the basis states may repeat, if one target state can
!    couple to several kl.
!---------------------------------------------------------------------
     Use bsr_mat

     Implicit none
     Integer, intent(in) :: ich,jch
     Real(8), intent(in) :: c, cc

     Integer :: i,j,ij

!     if(check_target.eq.0) Return

     if(ich.lt.1.or.ich.gt.nch) &
      Call Stop_mpi(0,ich,'UPDATE_CF: channel ich out of range')
     if(jch.lt.1.or.jch.gt.nch) &
      Call Stop_mpi(0,jch,'UPDATE_CF: channel jch out of range')

     i=max0(ich,jch); j=min0(ich,jch); ij=(i-1)*i/2+j
     htarg(ij) = htarg(ij) + C
     otarg(ij) = otarg(ij) + CC

     End Subroutine Target_h


!=====================================================================
     Subroutine Target_print(iout,eps)
!=====================================================================
     Use bsr_mat

     Implicit none
     Integer, intent(in) :: iout
     Real(8), intent(in) :: eps
     Integer :: i,j,ij,it,jt,ilen,k
     Real(8) :: c
     Character(20) :: targ1,targ2

     if(check_target.eq.0) Return

     write(iout,'(/a,1Pe9.1/)') 'Target hamiltonian errors if > eps_tar =',eps
     ilen = 0
     Do it=1,ntarg; i=Len_trim(BFT(it)); if(i.gt.ilen) ilen=i; End do

     k = 0
     Do i = 1,nch; it=iptar(i); targ1=BFT(it)
      Do j = 1,i;  jt=iptar(j); targ2=BFT(jt)
       ij=(i-1)*i/2+j
       c = abs(htarg(ij))
       if(i.eq.j) then
        htarg(ij) = htarg(ij) + EC
        c = htarg(ij) - Etarg(iptar(i))
       end if
        if(abs(C).lt.eps) Cycle
        if(iitar.ne.0) Cycle
        if(i.ne.j) write(iout,'(f14.6,5x,a,a6,5x,a,a6,4i6)') &
          C, targ1(1:ilen),ELC(i),targ2(1:ilen),ELC(j), it,i,jt,j
        if(i.eq.j) write(iout,'(f14.6,5x,a,a6,2i6)')  &
          C, targ1(1:ilen),ELC(i),it,i
        k=k+1
      End do
     End do

     write(iout,'(a,1Pe9.1/)') 'Target overlaps errors if > eps_tar =',eps

     k = 0
     Do i = 1,nch; it=iptar(i); targ1=BFT(it)
      Do j = 1,i;  jt=iptar(j); targ2=BFT(jt)
       ij=(i-1)*i/2+j
       c = otarg(ij)
       if(i.eq.j) c = c - 1.d0
       if(abs(C).lt.eps) Cycle
        if(i.ne.j) write(iout,'(f14.6,5x,a,a6,5x,a,a6,4i6)') &
          C, targ1(1:ilen),ELC(i),targ2(1:ilen),ELC(j), it,i,jt,j
        if(i.eq.j) write(iout,'(f14.6,5x,a,a6,2i6)')  &
          C, targ1(1:ilen),ELC(i),it,i
        k=k+1
       End do
     End do

     End Subroutine TARGET_print


!======================================================================
      Subroutine Target_new
!======================================================================
!     new target energies 
!----------------------------------------------------------------------
      Use bsr_mat,    only: htarg
      Use channel,    only: nch,iptar,ELC
      Use target,     only: ntarg,Etarg,BFT
      Use bsr_mat,    only: nun, AF_new, AF,pri, eps_tar, klsp,debug

      Implicit none
      Integer :: i,j,ij,it,jt,ilen,k,info
      Real(8) :: c
      Character(40) :: targ1,targ2
      Integer, allocatable :: itarget(:,:)
      Real(8), allocatable :: htarget(:,:), eval(:)

      ilen = 0
      Do it=1,ntarg; i=Len_trim(BFT(it)); if(i.gt.ilen) ilen=i; End do

      Allocate(itarget(ntarg,ntarg), htarget(ntarg,ntarg), eval(ntarg))

      k=0; itarget=0; htarget=0.d0 
      Do it=1,ntarg; htarget(it,it)=Etarg(it); End do

! ... find target matrix:

      Do i = 1,nch; it=iptar(i); targ1=BFT(it)
      Do j = 1,i;   jt=iptar(j); targ2=BFT(jt)
       if(it.eq.jt) Cycle
       ij=(i-1)*i/2+j
       c = abs(htarg(ij))
       if(abs(C).lt.eps_tar) Cycle

       if(k.eq.0) &
       write(pri,'(/a,1Pe9.1/)') &
         'Target hamiltonian errors if > eps_tar =',eps_tar
       if(i.ne.j) &
        write(pri,'(f14.6,5x,a,a6,5x,a,a6)') &
          C, targ1(1:ilen),ELC(i),targ2(1:ilen),ELC(j)
       if(i.eq.j) write(pri,'(f14.6,5x,a,a6)')  C, targ1(1:ilen),ELC(i)
       k=k+1

       itarget(it,jt)=itarget(it,jt)+1
       htarget(it,jt)=htarget(it,jt)+c

      End do
      End do

      if(k.eq.0) Return

      Do it=1,ntarg; Do jt=1,it
       if(itarget(it,jt).ne.0) htarget(it,jt)=htarget(it,jt)/itarget(it,jt)
       htarget(jt,it)=htarget(it,jt)
      End do; End do

      if(debug.gt.2) then
       write(pri,'(/a/)') 'Target hamiltonian:'
       Do i=1,ntarg
        write(pri,'(5f16.8)') htarget(i,1:i)
       End do
      end if

      Call LAP_DSYEV('N','L',ntarg,ntarg,htarget,eval,info)

      if(info.ne.0) Call Stop_mpi(0,info,'DSYEV failed in Target_new')
      
      ilen=LEN_TRIM(AF_new)+1
      AF=AF_new; write(AF(ilen:),'(a,i3.3)') '.',klsp

      write(pri,'(/a/)') 'new target energies:'      
      open(nun,file=AF)
      Do i=1,ntarg
       write(nun,'(3F20.8)') eval(i),Etarg(i),eval(i)-Etarg(i)       
       write(pri,'(3F20.8)') eval(i),Etarg(i),eval(i)-Etarg(i)       
      End do

      Deallocate(eval,itarget,htarget)

      End Subroutine TARGET_new
