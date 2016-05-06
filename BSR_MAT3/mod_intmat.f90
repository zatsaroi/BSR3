!======================================================================
      MODULE bsr_matrix
!======================================================================
!
!     contains the basic matrices:
!
!     hcc  -  interaction (or overlap) chanel-channel matrix
!     hcb  -  interaction (or overlap) chanel-bound matrix
!     hbb  -  interaction (or overlap) bound-bound matrix
!     ACF  -  array of asymptotic coefficients
!     htarg  -  interaction matrix for target states
!     otarg  -  overlap matrix for target states
!
!----------------------------------------------------------------------

      Implicit none

      Integer :: kch      !   number of channels
      Integer :: kns      !   size of one channel block (=ns)
      Integer :: kcp      !   number of perturber configurations
      Integer :: kmk      !   maximmum multipole
      Integer :: kcfg     !   number of configurations in channel blocks

      Real(8), allocatable :: hcc(:,:,:),hcb(:,:,:),hbb(:)
      Real(8), allocatable :: ACF(:,:,:)
      Real(8), allocatable :: htarg(:),otarg(:)
      Real(8), allocatable :: x(:,:)

      END MODULE bsr_matrix


!======================================================================
      Subroutine allocate_matrix(nch,ns,ncp,mk,ncfg)
!======================================================================
!     allocate arrays in module bsr_matrix
!----------------------------------------------------------------------
      Use bsr_matrix

      Implicit none
      Integer, intent(in) :: nch,ns,ncp,mk,ncfg
      Integer :: ich

      kch=nch; kns=ns; kcp=ncp; kmk=mk; kcfg=ncfg

      if(allocated(hcc)) Deallocate(hcc,ACF,htarg,otarg,x) 
      if(allocated(hcb)) Deallocate(hcb,hbb) 
      if(kch.eq.0.or.kns.eq.0) Return

      ich = nch*(nch+1)/2
      Allocate(hcc(kns,kns,ich),acf(kch,kch,0:kmk), &
               htarg(ich),otarg(ich),x(kns,kns))
      hcc=0.d0; acf=0.d0; htarg=0.d0; otarg=0.d0

      if(kcp.gt.0) then  
       Allocate(hbb(kcp*(kcp+1)/2),hcb(kns,kch,kcp))
       hbb=0.d0; hcb=0.d0
      end if

      End Subroutine allocate_matrix

!----------------------------------------------------------------------
!     next routines updates the matrixes;
!     we have coefficients only for the half of the matrix,  
!     so then we should h -->  h + h*   for diagonal blocks
!----------------------------------------------------------------------


!======================================================================
      Subroutine UPDATE_HX(ich,jch,ns,ks,d,sym)
!======================================================================
!     update channel block 
!
!     sym = 's'  -->  symmetric banded upper-column storage mode
!     sym = 'n'  -->  non-symmetric band matrix  
!     sym = 'x'  -->  non-symmetric full matrix
!----------------------------------------------------------------------
      Use bsr_matrix

      Implicit none
      Integer, intent(in) :: ich,jch,ns,ks
      Real(8), intent(in) :: d(ns,*)
      Character(1), intent(in) :: sym
      Integer ::  i,j,ij, ith, imin,imax

      if(ich.lt.1.or.ich.gt.kch) &
       Stop 'UPDATE_HX: channel index out of range'
      if(jch.lt.1.or.jch.gt.kch) &
       Stop 'UPDATE_HX: channel index out of range'
      if(ns.ne.kns) &
       Stop 'UPDATE_HL: different B-spline number'

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

      i=max(ich,jch); j=min(ich,jch); ij=i*(i-1)/2+j

      if(ich.ge.jch) then
       hcc(:,:,ij) = hcc(:,:,ij) + x
      else
       hcc(:,:,ij) = hcc(:,:,ij) + TRANSPOSE(x) 
      end if

      End Subroutine UPDATE_HX


!======================================================================
      Subroutine UPDATE_HL(ich,jch,ns,ks,d,c)
!======================================================================
!     update symmetric banded matrix in lower-column storage mode !
!----------------------------------------------------------------------
      Use bsr_matrix

      Implicit none
      Integer, intent(in) :: ich,jch,ns,ks
      Real(8), intent(in) :: c
      Real(8), intent(in) :: d(ns,ks)
      Integer :: i,j,ij,ith
      Real(8) :: S

      if(ich.lt.1.or.ich.gt.kch) &
       Stop 'UPDATE_HL: channel index out of range'
      if(jch.lt.1.or.jch.gt.kch) &
       Stop 'UPDATE_HL: channel index out of range'
      if(ns.ne.kns) &
       Stop 'UPDATE_HL: different B-spline number'

       i=max(ich,jch); j=min(ich,jch); ij=i*(i-1)/2+j

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
!     update scalar in the bound-bound block
!----------------------------------------------------------------------
      Use bsr_matrix

      Implicit none
      Integer, intent(in) :: ic,jc
      Real(8), intent(in) :: c
      Integer :: i,j,ij
      
      if(ic.gt.kcp.or.jc.gt.kcp.or.ic.le.0.or.jc.le.0) &
       Stop 'UPDATE_HB: indeces out of range'

      i=max(ic,jc); j=min(ic,jc); ij=i*(i-1)/2+j
      hbb(ij) = hbb(ij) + c

      End Subroutine UPDATE_HB


!======================================================================
      Subroutine UPDATE_HV(ich,ic,ns,v,c)
!======================================================================
!     update vector in the channel-bound block
!----------------------------------------------------------------------
      Use bsr_matrix

      Implicit none
      Integer, intent(in) :: ich,ic,ns
      Real(8), intent(in) :: c
      Real(8), intent(in) :: v(*)

      if(ich.lt.1.or.ich.gt.kch) &
       Stop 'UPDATE_HV: channel index out of range'
      if(ic.lt.1.or.ic.gt.kcp) &
       Stop 'UPDATE_HV: bound index out of range'

      hcb(1:ns,ich,ic) = hcb(1:ns,ich,ic) + C*v(1:ns)

      End Subroutine UPDATE_HV


!======================================================================
      Subroutine UPDATE_HW(ich,jch,ns,v,w)
!======================================================================
!     update channel block with by v*w
!-----------------------------------------------------------------------
      Use bsr_matrix

      Implicit none
      Integer, intent(in) :: ich,jch,ns
      Real(8), intent(in) :: v(*),w(*)
      Integer ::   i,j,ij

      if(ich.lt.1.or.ich.gt.kch) &
       Stop 'UPDATE_HW: channel index out of range'
      if(jch.lt.1.or.jch.gt.kch) &
       Stop 'UPDATE_HW: channel index out of range'
      if(ns.lt.1.or.ns.gt.kns) &
       Stop 'UPDATE_HW: B-spline index out of range'
 
      Do i = 1,ns
       Do j = 1,ns
        x(i,j) = v(i)*w(j)
       End do
      End do

      i=max(ich,jch); j=min(ich,jch); ij=i*(i-1)/2+j

      if(ich.ge.jch) then
       hcc(:,:,ij) = hcc(:,:,ij) + x
      else
       hcc(:,:,ij) = hcc(:,:,ij) + TRANSPOSE(x) 
      end if

      End Subroutine UPDATE_HW


!=======================================================================
      Subroutine UPDATE_CF(k,ich,jch,C)
!=======================================================================
!     Update asymptotic coefficients
!-----------------------------------------------------------------------
      Use bsr_matrix

      Implicit none
      Integer, intent(in) :: k, ich,jch
      Real(8), intent(in) :: C

      if(ich.lt.1.or.ich.gt.kch) &
       Stop 'UPDATE_CF: channel index out of range'
      if(jch.lt.1.or.jch.gt.kch) &
       Stop 'UPDATE_CF: channel index out of range'
      if(k.lt.0.or.k.gt.kmk) &
       Stop 'UPDATE_CF: multipole index out of range'

      acf(ich,jch,k) = acf(ich,jch,k) + C

      End Subroutine UPDATE_CF


!======================================================================
      Subroutine Target_h(ich,jch,C,CC)
!======================================================================
!     update target interaction and overlap matrixes, by considering 
!     the terms with structure <kl|k'l> <target|H|target'>
!
!     It is not pure target states, but the basis states before |kl>,
!     i.e. the basis states may repeat, if one target state can
!     couple to several kl.
!----------------------------------------------------------------------
      Use bsr_matrix

      Implicit none
      Integer, intent(in) :: ich,jch
      Real(8), intent(in) :: c, cc
      Integer :: i,j,ij

      if(ich.lt.1.or.ich.gt.kch) &
       Stop 'UPDATE_CF: channel index out of range'
      if(jch.lt.1.or.jch.gt.kch) &
       Stop 'UPDATE_CF: channel index out of range'

      i=max0(ich,jch); j=min0(ich,jch); ij=(i-1)*i/2+j
      htarg(ij) = htarg(ij) + C
      otarg(ij) = otarg(ij) + CC

      End subroutine Target_h


!======================================================================
      Subroutine Target_print(iout,EC,eps)
!======================================================================
      Use bsr_matrix
      Use channel
      Use target
      Use bsr_mat,  only: iitar

      Implicit none
      Integer, intent(in) :: iout
      Real(8), intent(in) :: EC,eps
      Integer :: i,j,ij,it,jt,ilen,k
      Real(8) :: c
      Character(20) :: targ1,targ2

      ilen = 0
      Do it=1,ntarg; i=Len_trim(BFT(it)); if(i.gt.ilen) ilen=i; End do

      k = 0
      Do i = 1,kch; it=iptar(i); targ1=BFT(it)
       Do j = 1,i;  jt=iptar(j); targ2=BFT(jt)
        ij=(i-1)*i/2+j
        c = abs(htarg(ij))
        if(i.eq.j) then
         htarg(ij) = htarg(ij) + EC
         c = htarg(ij) - Etarg(iptar(i))
        end if
         if(abs(C).lt.eps) Cycle
         if(iitar.ne.0) Cycle
         if(k.eq.0) &
         write(iout,'(/a,1Pe9.1/)') 'Target hamiltonian errors if > eps_tar =',eps
         if(i.ne.j) &
          write(iout,'(f14.6,5x,a,a6,5x,a,a6)') &
            C, targ1(1:ilen),ELC(i),targ2(1:ilen),ELC(j)
         if(i.eq.j) write(iout,'(f14.6,5x,a,a6)')  C, targ1(1:ilen),ELC(i)
         k=k+1
       End do
      End do

      k = 0
      Do i = 1,kch; it=iptar(i); targ1=BFT(it)
       Do j = 1,i;  jt=iptar(j); targ2=BFT(jt)
        ij=(i-1)*i/2+j
        c = otarg(ij)
        if(i.eq.j) c = c - 1.d0
        if(abs(C).lt.eps) Cycle
        if(k.eq.0) &
        write(iout,'(/a,1Pe9.1/)') 'Target overlaps errors if > eps_tar =',eps
        if(i.ne.j) &
         write(iout,'(f14.6,5x,a,a6,5x,a,a6)') &
           C, targ1(1:ilen),ELC(i),targ2(1:ilen),ELC(j)
         if(i.eq.j) write(iout,'(f14.6,5x,a,a6)')  C, targ1(1:ilen),ELC(i)
         k=k+1
        End do
      End do

      End Subroutine TARGET_print


!======================================================================
      Subroutine Target_new
!======================================================================
!     new target energies 
!----------------------------------------------------------------------
      Use bsr_matrix, only: htarg
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

      if(info.ne.0) Stop 'DSYEV failed in Target_new'
      
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
