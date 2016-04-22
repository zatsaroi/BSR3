!======================================================================
      MODULE spin_orbitals
!======================================================================
!
!     contains rather fancy but effective description of two determinant
!     wave function under consideration 
!
!----------------------------------------------------------------------

      USE configs, only: nsh ! defines the max. number of shells

      IMPLICIT NONE
      SAVE

! ... possible ml and ms values in chosen order: 

      Integer(4) :: mls_max
      Integer(4), Dimension(:), Allocatable :: ml_orb, ms_orb

! ... common list of different orbital symmetries for two determinants:
!
!     NSYM  - number of symmetries:  l,ml,ms ->  Lsym,Msym,Ssym 
!     IPSYM - pointer on the last given symmetry in the list
!     KSYM  - number of orbitals with given symmetry
!     nnsym - principal quantum numbers 
!     Isym  - pointer on the original position in configuration

      Integer(4) :: NSYM
      Integer(4), Allocatable, Dimension(:) :: Lsym,Lsym1,Lsym2
      Integer(4), Allocatable, Dimension(:) :: Msym,Msym1,Msym2
      Integer(4), Allocatable, Dimension(:) :: Ssym,Ssym1,Ssym2
      Integer(4), Allocatable, Dimension(:) :: IPsym1,IPsym2 
      Integer(4), Allocatable, Dimension(:) :: Ksym1 ,Ksym2 
      Integer(4), Allocatable, Dimension(:) :: nnsym1,nnsym2 
      Integer(4), Allocatable, Dimension(:) :: Isym1,Isym2 

! ... shell values:

!     md(i)  - the max.number of det.'s for the i-th subshell
!     nd(i)  - determinant under consideration
!     in(i)  - pointers on the orbitals of given shell
!     MS,ML  - shell MS,ML
!     MSp,MLp - intermediate values MS,ML

      Integer(4), Dimension(nsh) :: md,nd,in
      Integer(4), Dimension(nsh) :: MS,ML, MSp,MLp

      Integer(4) :: kz1,kz2     ! number of perturbations

! ... pointer of orbitals for given shell

      Integer(4), Allocatable, Dimension(:) :: Idet,Jdet

! ... auxiliary arrays

      Integer(4), Allocatable, Dimension(:) :: N1,N2,N3,N4, NP  

      End module spin_orbitals


