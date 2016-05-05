!======================================================================
      Module dmatrix
!======================================================================
!     contains the dipole matrices and routines for their update
!     DL  -  dipole transition matrix (L-form)
!     DV  -  dipole transition matrix (V-form)
!----------------------------------------------------------------------
      Implicit none

      Integer :: kns          !  size of one channel block (=ns)
      Integer :: kks          !  B-spline order (=ks)
      Integer :: kch1,kch2    !  number of channels
      Integer :: knc1,knc2    !  size of channel block (=kch*kns)
      Integer :: kcp1,kcp2    !  number of perturbers (not configurations !)
      Integer :: kcfg1,kcfg2  !  number of configurations
      Integer :: kdm1,kdm2    !  full size of interaction matrix

      Real(8), allocatable :: DL(:,:),DV(:,:)

      End Module dmatrix


!======================================================================
      Subroutine allocate_matrix(ns,ks,nch1,nch2,ncp1,ncp2,ncfg1,ncfg2)
!======================================================================
!     allocate arrays in module "dmatrix"
!----------------------------------------------------------------------    
      Use dmatrix

      Implicit none
      Integer, intent(in) :: ns,ks,nch1,nch2,ncp1,ncp2,ncfg1,ncfg2

      kns=ns; kks=ks
      kch1=nch1; kch2=nch2; kcp1=ncp1; kcp2=ncp2
      kcfg1=ncfg1; kcfg2=ncfg2 

      knc1 = kns * kch1; kdm1 = knc1 + kcp1
      knc2 = kns * kch2; kdm2 = knc2 + kcp2

      if(kdm1*kdm2.le.0) Stop 'DMATRIX: kdm1*kdm2 <= 0'

      if(allocated(DL)) Deallocate(DL) 
      if(allocated(DV)) Deallocate(DV) 
     
      Allocate( DL(kdm1,kdm2), DV(kdm1,kdm2) )
      
      DL = 0.d0; DV = 0.d0

      End Subroutine allocate_matrix


!======================================================================
      Subroutine UPDATE_DX(ich,jch,ns,ks,d1,d2)
!======================================================================
!     update channel block with nonsymmetric band matrix 
!----------------------------------------------------------------------
      Use dmatrix

      Implicit none
      Integer, intent(in) :: ich,jch,ns,ks
      Real(8), intent(in) :: d1(ns,ks+ks-1),d2(ns,ks+ks-1)

      Integer :: ishft, jshft, i,j, jp, ihm, jhm, imin, imax

      ishft = (ich-1)*ns
      jshft = (jch-1)*ns

      Do jp = 1,ks+ks-1
       imin=max0( 1, 1 + ks-jp)
       imax=min0(ns,ns + ks-jp)
       Do i = imin,imax;  j=i+jp-ks
        ihm = ishft + i; jhm = jshft + j
        DL(ihm,jhm) = DL(ihm,jhm) + d1(i,jp)
        DV(ihm,jhm) = DV(ihm,jhm) + d2(i,jp)
       End do
      End do

      End Subroutine UPDATE_DX


!======================================================================
      Subroutine UPDATE_DW(ich,jch,ns,v1,w1,v2,w2)
!======================================================================
!     update full matrix with v*w
!-----------------------------------------------------------------------
      Use dmatrix

      Implicit none
      Integer, intent(in) :: ich,jch,ns
      Real(8), intent(in) :: v1(ns),w1(ns),v2(ns),w2(ns)
      Integer ::  i,j, ii,jj

      j = (jch-1)*ns
      Do jj = 1,ns;   j = j + 1
       i = (ich-1)*ns
      Do ii = 1,ns;   i = i + 1
       DL(i,j) = DL(i,j) + v1(ii)*w1(jj)
       DV(i,j) = DV(i,j) + v2(ii)*w2(jj)
      End do;  End do
     
      End Subroutine UPDATE_DW


!======================================================================
      Subroutine UPDATE_DV(ich,jch,ic,jc,ns,v1,v2)
!======================================================================
!     update full matrix with vector
!----------------------------------------------------------------------
      Use dmatrix

      Implicit none
      Integer, intent(in) :: ich,jch,ic,jc,ns
      Real(8), intent(in) :: v1(ns),v2(ns)
      Integer :: i,j

      if(ich.gt.0.and.jch.eq.0) then
       i = (ich-1)*ns
       j = knc2 + jc; if(j.gt.kdm2) Stop 'UPDATE_DV: jc too big'
       DL(i+1:i+ns,j) = DL(i+1:i+ns,j) + v1(1:ns)
       DV(i+1:i+ns,j) = DV(i+1:i+ns,j) + v2(1:ns)
      elseif(jch.gt.0.and.ich.eq.0) then
       j = (jch-1)*ns
       i = knc1 + ic; if(i.gt.kdm1) Stop 'UPDATE_DV: ic too big'
       DL(i,j+1:j+ns) = DL(i,j+1:j+ns) + v1(1:ns)
       DV(i,j+1:j+ns) = DV(i,j+1:j+ns) + v2(1:ns)
      else
       Stop ' UPDATE_DV: non-allowed combination of ich,jch' 
      end if

      End Subroutine UPDATE_DV


!======================================================================
      Subroutine UPDATE_DB(ic,jc,CL,CV)
!======================================================================
!     update full matrix with scalar
!----------------------------------------------------------------------
      Use dmatrix

      Implicit none
      Integer, intent(in) :: ic,jc
      Real(8), intent(in) :: CL,CV
      Integer ::  i, j

      i = knc1 + ic; if(i.gt.kdm1) Stop 'UPDATE_DB: ic too big'
      j = knc2 + jc; if(j.gt.kdm2) Stop 'UPDATE_DB: jc too big'
      DL(i,j) = DL(i,j) + CL
      DV(i,j) = DV(i,j) + CV

      End Subroutine UPDATE_DB

