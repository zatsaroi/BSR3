!======================================================================
      MODULE bsr_matrix
!======================================================================
!     contains the basic matrices
!
!     hcc  -  interaction (or overlap) chanel-channel matrix
!     hcb  -  interaction (or overlap) chanel-bound matrix
!     hbb  -  interaction (or overlap) bound-bound matrix
!     ACF  -  array of asymptotic coefficients
!     htarg  -  interaction matrix for target states
!     otarg  -  overlap matrix for target states
!----------------------------------------------------------------------

      Implicit none

      Integer :: kch      !   number of channels
      Integer :: kns      !   size of one channel block (=ns)
      Integer :: kcp      !   number of perturber configurations
      Integer :: kmk      !   maximmum multipole
      Integer :: kcfg     !   number of configurations in channel blocks
      Integer :: kprocs   !   number of processors

      Real(8), allocatable :: hcc(:,:,:),hcb(:,:),hbb(:)
      Real(8), allocatable :: ACF(:,:,:)
      Real(8), allocatable :: htarg(:),otarg(:)
      Real(8), allocatable :: x(:,:)
      Integer, allocatable :: icc(:,:), icb(:,:), ibb(:,:), imycase(:,:)

      Integer :: iicc, iicb, iibb  !  dimensions


      END MODULE bsr_matrix


!======================================================================
      Subroutine allocate_matrix(nch,ns,ncp,mk,ncfg,m)
!======================================================================
!     allocate basic matrices
!----------------------------------------------------------------------
      Use mpi
      Use bsr_matrix

      Implicit none

      Integer, intent(in) :: nch,ns,ncp,mk,ncfg
      Integer :: ich,jch, i,k, myid,ierr, m
      Integer, external :: mycase

      kch=nch; kns=ns; kcp=ncp; kmk=mk; kcfg=ncfg

      m = 0
      if(Allocated(hcc)) Deallocate(hcc,ACF,htarg,otarg,x,icc) 
      if(Allocated(hcb)) Deallocate(hcb,hbb,icb,ibb) 
      if(kch.eq.0.or.kns.eq.0) Return

      Call MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierr)
      Call MPI_COMM_SIZE(MPI_COMM_WORLD, kprocs, ierr)
      kprocs=kprocs-1
   
      Allocate(icc(kch,kch));  icc = 0; iicc = 0
      m=m+kch*kch
      k = 0; i=0
      Do ich=1,kch; Do jch=1,ich
       k=k+1; if(k.gt.kprocs) k=1
       if(myid.ne.k) Cycle       
       i=i+1; icc(ich,jch) = i; icc(jch,ich) = i
      End do; End do
      iicc=i
      ich = (kch+1)*kch/2
      Allocate(hcc(kns,kns,iicc),acf(kch,kch,0:kmk), &
               htarg(ich),otarg(ich),x(kns,kns))
      m = m + 2*(kns*kns*iicc + kch*kch*(kmk+1) + 2*ich + kns*kns)
      hcc=0.d0; acf=0.d0; htarg=0.d0; otarg=0.d0

      if(kcp.gt.0) then  

      Allocate(icb(kch,kcp));  icb = 0
      m = m + 2 * kch * kch
      k = 0; i=0
      Do ich=1,kch; Do jch=1,kcp
       k=k+1; if(k.gt.kprocs) k=1
       if(myid.ne.k) Cycle       
       i=i+1; icb(ich,jch) = i
      End do; End do
      iicb = i
      Allocate(hcb(kns,iicb));  hcb=0.d0
      m = m + 2 * kns * iicb
      Allocate(ibb(kcp,kcp));  ibb = 0
      m = m +  kcp * kcp
      k = 0; i=0
      Do ich=1,kcp; Do jch=1,ich
       k=k+1; if(k.gt.kprocs) k=1
       if(myid.ne.k) Cycle       
       i=i+1; ibb(ich,jch) = i; ibb(jch,ich) = i
      End do; End do
      iibb = i
      Allocate(hbb(iibb));  hbb=0.d0
      m = m + iibb

      end if

      if(allocated(imycase)) Deallocate(imycase)
      Allocate( imycase(1:kch+kcp,1:kch+kcp) )
      m = m + (kch+kcp)*(kch+kcp)
      imycase = 0

      Do ich=1,kch+kcp; Do jch=1,kch+kcp
       imycase(ich,jch) = mycase(ich,jch)
      End do; End do

      END subroutine allocate_matrix


!======================================================================
     Integer Function mycase(ich,jch)
!======================================================================

     Use bsr_mat, only: myid
     Use bsr_matrix

     Integer, intent(in) :: ich,jch

     if(ich.le.kch.and.jch.le.kch) then
      mycase = icc(ich,jch)
     elseif(ich.le.kch.and.jch.gt.kch) then
      mycase = icb(ich,jch-kch)
     elseif(ich.gt.kch.and.jch.le.kch) then
      mycase = icb(jch,ich-kch)
     elseif(ich.gt.kch.and.jch.gt.kch) then
      mycase = ibb(ich-kch,jch-kch)
     else
      mycase = 0
     end if

     End Function mycase

!----------------------------------------------------------------------
!    next routines for updating the matreces;
!    we have coefficients only for the half of matrix,  
!    so then we should h -->  h + h*   for diagonal blocks
!----------------------------------------------------------------------


!======================================================================
     Subroutine UPDATE_HX(ich,jch,ns,ks,d,sym)
!======================================================================
!    update channel block 
!
!    sym = 's'  -->  symmetric banded upper-column storage mode
!    sym = 'n'  -->  non-symmetric band matrix  
!    sym = 'x'  -->  non-symmetric full matrix
!----------------------------------------------------------------------

     Use bsr_matrix

     Implicit none
     Integer, intent(in) :: ich,jch,ns,ks
     Real(8), intent(in) :: d(ns,*)
     Character(1), intent(in) :: sym

     Integer ::  i,j,ij, ith, imin,imax

     if(ich.lt.1.or.ich.gt.kch) &
        Call Stop_mpi(0,ich,'UPDATE_HX: ich index out of range')
     if(jch.lt.1.or.jch.gt.kch) &
        Call Stop_mpi(0,jch,'UPDATE_HX: jch index out of range')
     if(ns.ne.kns) &
        Call Stop_mpi(0,ns,'UPDATE_HL: different B-spline number')

     x = 0.d0
     Select case(sym)

     Case('s')

       do ith=1,ks;  do i=1,ns-ith+1;  j=i+ith-1
         x(i,j) = d(i,ith)
         x(j,i) = d(i,ith)
       end do; end do

     Case('n')

       do ith = 1,ks+ks-1
         imin=max0( 1, 1 + ks-ith)
         imax=min0(ns,ns + ks-ith)
         do i = imin,imax
           j = i + ith - ks
           x(i,j) = d(i,ith)
          end do
        end do

      Case('x')

        x(1:ns,1:ns) = d(1:ns,1:ns)

      End Select

      ij = icc(ich,jch)
      if(ij.le.0.or.ij.gt.iicc) &
       Call Stop_mpi(0,ij,'UPDATE_HX: ij index out of range')

      if(ich.ge.jch) then
       hcc(:,:,ij) = hcc(:,:,ij) + x
      else
       hcc(:,:,ij) = hcc(:,:,ij) + TRANSPOSE(x) 
      end if

     END subroutine UPDATE_HX


!======================================================================
     Subroutine UPDATE_HL(ich,jch,ns,ks,d,c)
!======================================================================
!    update symmetric banded matrix in lower-column storage mode !
!----------------------------------------------------------------------

     Use bsr_matrix

     Implicit none
     Integer(4), intent(in) :: ich,jch,ns,ks
     Real(8), intent(in) :: c
     Real(8), intent(in) :: d(ns,ks)

     Integer(4) :: i,j,ij,ith
     Real(8) :: S

      if(ich.lt.1.or.ich.gt.kch) &
       Call Stop_mpi(0,ich,'UPDATE_HL: ich index out of range')
      if(jch.lt.1.or.jch.gt.kch) &
       Call Stop_mpi(0,jch,'UPDATE_HL: jch index out of range')
      if(ns.ne.kns) &
       Call Stop_mpi(0,ns,'UPDATE_HL: different B-spline number')

      ij = icc(ich,jch)
      if(ij.le.0.or.ij.gt.iicc) Stop 'UPDATE_HL: ij index out of range'

       Do ith = 1,ks
        Do i = ks+1-ith,ns
         j = i + ith - ks
         S = c*d(i,ith)
         hcc(i,j,ij) = hcc(i,j,ij) + S
         if(i.ne.j) hcc(j,i,ij) = hcc(j,i,ij) + S
        End do
       End do

     END subroutine UPDATE_HL


!======================================================================
     Subroutine UPDATE_HB(ichm,jchm,c)
!======================================================================
!    update scalar in the bound-bound block
!----------------------------------------------------------------------

     Use bsr_matrix

     Implicit none
     Integer, intent(in) :: ichm,jchm
     Real(8), intent(in) :: c
     Integer :: ic,jc,ij
      
     ic = ichm; jc = jchm

     if(ic.gt.kcp.or.jc.gt.kcp.or.ic.le.0.or.jc.le.0) &
       Call Stop_mpi(0,ic,'UPDATE_HB: indeces out of range')

     ij = ibb(ic,jc)
     if(ij.le.0.or.ij.gt.iibb) &
      Call Stop_mpi(0,ij,'UPDATE_HB: ij index out of range')

     hbb(ij) = hbb(ij) + c

     END Subroutine UPDATE_HB


!======================================================================
     Subroutine UPDATE_HV(ich,ichm,ns,v,c)
!======================================================================
!    update vector in the channel-bound block
!----------------------------------------------------------------------

     Use bsr_matrix

     Implicit none
     Integer, intent(in) :: ich,ichm,ns
     Real(8), intent(in) :: c
     Real(8), intent(in) :: v(*)
     Integer :: ic, ij

     if(ich.lt.1.or.ich.gt.kch) &
       Call Stop_mpi(0,ich,'UPDATE_HV: channel index out of range')
     ic = ichm
     if(ic.lt.1.or.ic.gt.kcp) &
       Call Stop_mpi(0,ic,'UPDATE_HV: bound index out of range')

     ij = icb(ich,ic)
     if(ij.le.0.or.ij.gt.iicb) &
       Call Stop_mpi(0,ij,'UPDATE_HV: ij index out of range')

     hcb(1:ns,ij) = hcb(1:ns,ij) + C*v(1:ns)

     END subroutine UPDATE_HV


!======================================================================
     Subroutine UPDATE_HW(ich,jch,ns,v,w)
!======================================================================
!    update channel block with by v*w
!----------------------------------------------------------------------

     Use bsr_matrix

     Implicit none
     Integer, intent(in) :: ich,jch,ns
     Real(8), intent(in) :: v(*),w(*)
     Integer ::   i,j,ij

     if(ich.lt.1.or.ich.gt.kch) then
       write(*,*) 'UPDATE_HW: ich index out of range',ich
       Return
       Call Stop_mpi(0,ich,'UPDATE_HW: ich index out of range')
     end if      

     if(jch.lt.1.or.jch.gt.kch) then
       write(*,*) 'UPDATE_HW: jch index out of range',jch
       Return
       Call Stop_mpi(0,jch,'UPDATE_HW: jch index out of range')
     end if      

     if(ns.lt.1.or.ns.gt.kns) &
       Call Stop_mpi(0,ns,'UPDATE_HW: B-spline index out of range')
 
      Do i=1,ns;  Do j=1,ns;  x(i,j)=v(i)*w(j);  End do;  End do

      ij = icc(ich,jch)
      if(ij.le.0.or.ij.gt.iicc) then
        write(*,*) 'ich,jch,ij', ich,jch,ij
        Call Stop_mpi(0,ij, 'UPDATE_HW: ij index out of range')
      end if

      if(ich.ge.jch) then
       hcc(:,:,ij) = hcc(:,:,ij) + x
      else
       hcc(:,:,ij) = hcc(:,:,ij) + TRANSPOSE(x) 
      end if


     END subroutine UPDATE_HW


!=======================================================================
     Subroutine UPDATE_CF(k,ich,jch,C)
!=======================================================================
!    Update asymptotic coefficients
!-----------------------------------------------------------------------

     Use bsr_matrix

     Implicit none
     Integer, intent(in) :: k, ich,jch
     Real(8), intent(in) :: C

     if(ich.lt.1.or.ich.gt.kch) &
       Call Stop_mpi(0,ich,'UPDATE_CF: ich index out of range')
     if(jch.lt.1.or.jch.gt.kch) &
       Call Stop_mpi(0,jch,'UPDATE_CF: jch index out of range')
     if(k.lt.0.or.k.gt.kmk) &
       Call Stop_mpi(0,k,'UPDATE_CF: multipole index out of range')

     acf(ich,jch,k) = acf(ich,jch,k) + C

     END SUBROUTINE UPDATE_CF


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
       Call Stop_mpi(0,ich,'UPDATE_CF: channel ich out of range')
      if(jch.lt.1.or.jch.gt.kch) &
       Call Stop_mpi(0,jch,'UPDATE_CF: channel jch out of range')

       i=max0(ich,jch); j=min0(ich,jch); ij=(i-1)*i/2+j
       htarg(ij) = htarg(ij) + C
       otarg(ij) = otarg(ij) + CC

      End subroutine Target_h


!======================================================================
      Subroutine Target_print(iout,EC,eps)
!======================================================================

      Use bsr_matrix
      Use channel, only: iptar
      Use target, only: Etarg

      Implicit none
      Integer, intent(in) :: iout
      Real(8), intent(in) :: EC,eps
      Integer :: i,j,ij
      Real(8) :: c

      write(iout,'(/a,e12.3/)') 'Target hamiltonian errors if > ',eps
      Do i = 1,kch
       Do j = 1,i
        ij=(i-1)*i/2+j
        if(i.eq.j) then
         htarg(ij) = htarg(ij) + EC
         c = htarg(ij) - Etarg(iptar(i))
         if(abs(C).gt.eps)  &
          write(iout,'(4i5,2f16.6)') i,j,iptar(i),iptar(j),htarg(ij),C
         else
          c = abs(htarg(ij))
          if(abs(C).gt.eps) &
          write(iout,'(2i5,f16.6)') i,j,htarg(ij)
         end if
       End do
      End do

      write(iout,'(/a,e12.3/)') 'Target overlaps errors if > ',eps
      Do i = 1,kch
       Do j = 1,i
        ij=(i-1)*i/2+j
        c = otarg(ij)
        if(i.eq.j) c = c - 1.d0
        if(abs(C).gt.eps) &
         write(iout,'(4i5,f16.6)') i,j,iptar(i),iptar(j),c
       End do
      End do

     End SUBROUTINE TARGET_print



!======================================================================
      Subroutine Collect_ACF
!======================================================================

      Use MPI

      Use bsr_mat
      Use bsr_matrix
      Use conf_LS, only: nclosd 
      Use orb_LS, only: LEF
      Use target, only: nelc

      Implicit none

      Integer :: status(MPI_STATUS_SIZE)
      Integer  :: i,j,k, i1,i2,ij
      Real(8) :: c
      Character(200) :: line

! ... Summerize the ACF - matrix:

      Call MPI_BARRIER(MPI_COMM_WORLD, ierr)

      Do i = 1,kch
       Do j = 1,kch
        if(icc(i,j).ne.0) then
         Call MPI_SEND(ACF(i,j,:),kmk+1,MPI_DOUBLE_PRECISION, &
                       0, 0, MPI_COMM_WORLD, ierr)
        end if
        if(myid.eq.0) then
         Call MPI_RECV(ACF(i,j,:),kmk+1,MPI_DOUBLE_PRECISION, MPI_ANY_SOURCE, &
                      MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr)
        end if

       Call MPI_BARRIER(MPI_COMM_WORLD, ierr)

       End do
      End do

      if(myid.ne.0) Return  

! ... Symmetrize the ACF - matrix:

      Do k = 0,kmk
       Do i = 1,kch
        Do j = 1,i
         C = ACF(i,j,k) + ACF(j,i,k); if(i.ne.j) C=C*2.d0
         if(abs(C).lt.0.00001) C=0.d0
         ACF(i,j,k) = C; ACF(j,i,k) = C
        End do
       End do
      End do

! ... add corrections from the core screening:

      j=0; Do i=1,nclosd; j=j+2*(4*LEF(i)+2); End do
      Do i=1,kch; ACF(i,i,0)=ACF(i,i,0)+j;  End do

! ... print the asymptotic coefficients:

      write(pri,'(/a,i3/)') 'Asymptotic coefficients: mk = ',mk
      write(pri,'(a,i3/)') 'For k=0 should be equal to 2*nelc = ', &
                                                        2*nelc
      write(pri,'(a,i3/)') 'Derivations (if > 0.00001):'
      Do i = 1,kch
       if(abs(ACF(i,i,0)-2*nelc).lt.0.00001) Cycle
       write(pri,'(i5,2F15.6)') i,ACF(i,i,0)-2*nelc
      End do

      if(debug.gt.1) then
      write(pri,'(/a/)') 'Asymptotic coefficients: i,j, ACF(i,j,k)'
      Do k=0,mk
       if(SUM(acf(:,:,k)).eq.0) Cycle
       write(pri,'(a,i2)') 'k = ',k
       ij = 0
       Do i=1,kch; Do j = 1,i      
        if(abs(acf(i,j,k)).lt.eps_acf) Cycle
        i1=ij*20+1; i2=i1+19 
        write(line(i1:i2),'(2i4,E12.3)') j,i,acf(i,j,k)
        ij=ij+1
        if(ij.lt.5) Cycle
        write(pri,'(a)') line; ij=0
       End do; End do
       if(ij.eq.0) Cycle
       i1=1; i2=ij*20
       write(pri,'(a)') line(i1:i2)
      End do
      end if

      End Subroutine Collect_ACF


!======================================================================
      Subroutine Collect_otarg
!======================================================================

      Use MPI
      Use bsr_matrix

      Implicit none
      Integer :: status(MPI_STATUS_SIZE)
      Integer  :: i,j,ij,myid,ierr

      Call MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierr)
      Call MPI_BARRIER(MPI_COMM_WORLD, ierr)

      Do i = 1,kch
       Do j = 1,i
        ij=(i-1)*i/2+j
        if(icc(i,j).ne.0) then
         Call MPI_SEND(otarg(ij),1,MPI_DOUBLE_PRECISION, &
                       0, i, MPI_COMM_WORLD, ierr)
        end if
        if(myid.eq.0) then
         Call MPI_RECV(otarg(ij),1,MPI_DOUBLE_PRECISION, MPI_ANY_SOURCE, &
                      MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr)
        end if
       Call MPI_BARRIER(MPI_COMM_WORLD, ierr)

       End do
      End do

      End Subroutine Collect_otarg  

!======================================================================
      Subroutine Collect_htarg
!======================================================================

      Use MPI
      Use bsr_matrix

      Implicit none

      Integer :: status(MPI_STATUS_SIZE)
      Integer  :: i,j,ij,myid,ierr

      Call MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierr)
      Call MPI_BARRIER(MPI_COMM_WORLD, ierr)

      Do i = 1,kch
       Do j = 1,i
        ij=(i-1)*i/2+j
        if(icc(i,j).ne.0) then
         Call MPI_SEND(htarg(ij),1,MPI_DOUBLE_PRECISION, &
                       0, i, MPI_COMM_WORLD, ierr)
        end if
        if(myid.eq.0) then
         Call MPI_RECV(htarg(ij),1,MPI_DOUBLE_PRECISION, MPI_ANY_SOURCE, &
                      MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr)
        end if
       Call MPI_BARRIER(MPI_COMM_WORLD, ierr)

       End do
      End do

      End Subroutine Collect_htarg  