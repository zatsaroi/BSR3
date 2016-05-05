!======================================================================
      Module dets
!======================================================================
!     contains information about overlap determinants {<i|j>}
!----------------------------------------------------------------------
!
! KPD(i)  - size of i-th determinant
!
! IPD(i)  - pointer (ip) of i-th determinants in the NPD array
!
! NPD(ip+1:ip+kd) - contains information about orbitals in overlap
!                   determinant as  ii * idet + jj, 
!                   where ii - pointer on row orbitals
!                         jj - pointer on column orbitals
!
!----------------------------------------------------------------------     
!
! KPF(i)  - number (kd) of determinants in i-th overlap factor
!
! IPF(i)  - pointer (ip) on the i-th overlap factor in the NPF array
!
! NPF(ip+1:ip+kd) - i-th overlap factors as a list of pointers 
!                   on individual determinants (ip) and its power (iext): 
!                   npf(i) = ip * idef  +  iext
!----------------------------------------------------------------------
! parameters idet,idef come from program breit_bsr
!----------------------------------------------------------------------      
      Implicit none
     
      Integer :: ndet  =   0   !  number of determinants
      Integer :: kdet  =   0   !  sum of all det. dimensions      
      Integer, parameter :: idet  = 2**15 !  pack basis 
	
      Integer, allocatable :: KPD(:), IPD(:), NPD(:), JPD(:)

      Integer :: ndef   =   0  ! number of overlap factors 
      Integer :: kdef   =   0  ! sum of all overlap factor dimensions  
      Integer, parameter :: idef   =  16  ! pack basis 
      
      Integer, allocatable :: KPF(:), IPF(:), NPF(:), JPF(:)

      END MODULE dets


!=======================================================================
      Subroutine Read_dets(nub)
!=======================================================================
!     read the overlap determinants from INT_BNK (unit nub)  
!-----------------------------------------------------------------------
      Use dets

      read(nub) ndet,ldet,jdet

      if(ndet.eq.0) Return
      if(allocated(KPD)) Deallocate (KPD,IPD,NPD,JPD)
      mdet = ndet; kdet = ldet
      Allocate(KPD(mdet),IPD(mdet),NPD(kdet),JPD(mdet))
      Do i = 1,ndet
       read(nub) kpd(i),ipd(i),jpd(i),NPD(ipd(i)+1:ipd(i)+kpd(i))
      End do

      read(nub) ndef,ldef,jdef

      if(ndef.eq.0) Return
      if(allocated(kpf)) Deallocate (KPF,IPF,NPF,JPF)
      mdef = ndef; kdef = ldef
      Allocate(KPF(mdef),IPF(mdef),NPF(kdef),JPF(mdef))
      Do i = 1,ndef
       read(nub) kpf(i),ipf(i),jpf(i),npf(ipf(i)+1:ipf(i)+kpf(i))
      End do

      End Subroutine Read_dets

