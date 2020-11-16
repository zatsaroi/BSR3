!=======================================================================
      Module hf_orbitals
!=======================================================================
!     contains description of atomic orbitals 
!-----------------------------------------------------------------------
      Implicit none

      Integer :: nbf = 0                  ! number of orbitals
      Integer, allocatable :: nbs(:)      ! n-values
      Integer, allocatable :: lbs(:)      ! l-value
      Integer, allocatable :: ibs(:)      ! set numbers, not used here
      Integer, allocatable :: mbs(:)      ! number of splines
      Integer, allocatable :: iord(:)     ! order
      Real(8), allocatable :: qsum(:)     ! ocupation
      Real(8), allocatable :: dpm(:)      ! dicrepancy
      Logical, allocatable :: clsd(:)     ! closed shells
      Real(8), allocatable :: e(:,:)      ! orbital orthogonality
      Real(8), allocatable :: p(:,:)      ! radial functions
      Integer, allocatable :: iprm(:,:)   ! boundary conditions
      Character(4), allocatable :: ebs(:) ! spectroscopic notation

      End Module hf_orbitals


!=======================================================================
      Subroutine alloc_hf_orbitals(n)
!=======================================================================
!     allocates, deallocates or reallocate arrays in module hf_orbitals 
!-----------------------------------------------------------------------
      Use hf_orbitals

      Implicit none
      Integer, intent(in) :: n

      if(allocated(nbs)) &
        Deallocate(nbs,lbs,ibs,mbs,ebs,e,iord,qsum,dpm,clsd)
        nbf=0
      if(n.le.0) return
      nbf=n
      Allocate(nbs(nbf),ibs(nbf),mbs(nbf),lbs(nbf),ebs(nbf), &
               e(nbf,nbf),iord(nbf),qsum(nbf),dpm(nbf),clsd(nbf))
      nbs=0; ibs=0; mbs=0; lbs=0; iord=0
      qsum=0.d0; dpm=0.d0; clsd = .false.
 
      End Subroutine alloc_hf_orbitals

!=======================================================================
      Subroutine alloc_hf_radial(ns)
!=======================================================================
!     allocates, deallocates or reallocate arrays in module hf_orbitals 
!-----------------------------------------------------------------------
      Use hf_orbitals

      Implicit none
      Integer, intent(in) :: ns

      if(allocated(p))  Deallocate(p,iprm)
      if(ns.le.0) return
      Allocate(p(ns,nbf),iprm(ns,nbf))
      p = 0.d0; iprm=1
 
      End Subroutine alloc_hf_radial

!=======================================================================
      Integer Function Ifind_orb(n,l,iset)
!=======================================================================
!     find orbital (n,l,iset) in the list hf_orbitals
!----------------------------------------------------------------------      
      Use hf_orbitals
      Implicit none
      Integer, intent(in) :: n,l,iset
      Integer :: i
      Ifind_orb=0

      Do i=1,nbf
       if(n.ne.nbs(i)) Cycle
       if(l.ne.lbs(i)) Cycle
       if(iset.ne.ibs(i)) Cycle
       Ifind_orb = i
       Return
      End do

      End Function Ifind_orb

