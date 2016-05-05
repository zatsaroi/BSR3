!======================================================================
      Module spin_orbitals
!======================================================================
!     contains rather fancy but effective description of two determinant
!     wave function under consideration 
!----------------------------------------------------------------------
      Implicit none

! ... possible ml and ms values in chosen order: 

      Integer, allocatable :: ml_orb(:), ms_orb(:)

! ... common list of different orbital symmetries for two determinants:
!
!     NSYM  - number of symmetries:  l,ml,ms ->  Lsym,Msym,Ssym 
!     IPSYM - pointer on the last given symmetry in the list
!     KSYM  - number of orbitals with given symmetry
!     nnsym - principal quantum numbers 
!     Isym  - pointer on the original position in configuration

      Integer :: nsym, nsym1
      Integer, allocatable ::  Lsym (:),Msym(:),Ssym(:),               &
      Lsym1(:),Msym1(:),Ssym1(:),IPsym1(:),Ksym1(:),nnsym1(:),Isym1(:),&
      Lsym2(:),Msym2(:),Ssym2(:),IPsym2(:),Ksym2(:),nnsym2(:),Isym2(:) 

! ... shell values:

!     md(i)  - the max.number of det.'s for the i-th subshell
!     in(i)  - pointers on the orbitals of given shell

      Integer, allocatable :: md(:),in(:)

      Integer :: kz1,kz2     ! number of perturbations

      End Module spin_orbitals


!======================================================================
      Subroutine Alloc_spin_orbitals(ne)
!======================================================================
!     allocate arrays in module spin-orbitals
!----------------------------------------------------------------------
      Use spin_orbitals

      Implicit none
      Integer, intent(in) :: ne
      Integer :: i,l,m,no
      Integer, external :: ML_id, MS_id
      
! ... define possible mls orbitals:

      Call Def_maxl(l); no=4*l+2
      if(Allocated(ml_orb)) Deallocate(ml_orb,ms_orb)
      Allocate(ml_orb(no), ms_orb(no))

      Do i = 1,no
       ml_orb(i) = (ML_id(i)-1)/2
       ms_orb(i) =  MS_id(i)
      End do

! ... other allocations:

      if(.not.allocated(Lsym)) then
       m = 2*ne
       Allocate(Lsym(m),Lsym1(m),Lsym2(m),  &
                Msym(m),Msym1(m),Msym2(m),  &
                Ssym(m),Ssym1(m),Ssym2(m),  &
                IPsym1(m),IPsym2(m),Ksym1(m),Ksym2(m), &
                nnsym1(m),nnsym2(m),Isym1(m),Isym2(m), &
                md(ne),in(ne) )
      end if      

      End Subroutine Alloc_spin_orbitals
