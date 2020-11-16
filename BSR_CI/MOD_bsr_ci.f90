!======================================================================
      Module bsr_ci
!======================================================================
!     contains the common variables and allocatable arrays for ci_bnk
!----------------------------------------------------------------------
      Use conf_ls
      Use orb_ls
      Use spline_atomic
      Use spline_param
      Use spline_galerkin
      Use spline_grid
      Use spline_orbitals
      Use spline_hl
      Use spline_integrals
      Use spline_moments
      Use spline_slater

      Implicit none

! ... name of case:

      Integer, parameter :: ma = 80
      Character(ma) :: name = ' '   

! ... unit numbers:

      Integer :: iread = 55; Character(ma) :: AF_inp='name.inp'
      Integer :: iwrite= 66; Character(ma) :: AF_log='name.log'

      Integer :: nuc = 1;   Character(ma) :: AF_c='name.c'
      Integer :: nuw = 2;   Character(ma) :: AF_w='name.bsw'
      Integer :: nub = 3;   Character(ma) :: AF_b='name.bnk'
                            Character(ma) :: BF_b='int_bnk'
      Integer :: nuk = 4;   Character(ma) :: AF_k='name.knot'
                            Character(ma) :: BF_k='knot.dat'
      Integer :: nud = 10;  Character(ma) :: AF_d='name.d'
      Integer :: nul = 11;  Character(ma) :: AF_l='name.l'
      Integer :: nuj = 12;  Character(ma) :: AF_j='name.j'

      Integer :: nuh = 21   !  scratch file for LS-matrix
      Integer :: nuo = 22   !  scratch file for overlap
      Integer :: nur = 23   !  scratch file for J-part

      Character(ma) :: AF

      Character(2) :: atom

! ... matrices:

      Real(8), allocatable :: HM(:,:)    ! interaction matrix
      Real(8), allocatable :: HS(:,:)    ! overlap matrix
      Real(8), allocatable :: HD(:)      ! if NZERO > 0

      Real(8), allocatable :: PIP(:,:)  ! one-electron overlaps

      Integer :: pri_mat = 0            ! output matrix

! ... relativistic operators included: 

!     Integer :: irel =  0    !  pointer on rel.calculations -> module spline_atomic
      Integer :: iso  =  0    !  spin-orbit
      Integer :: isoo =  0    !  spin-other-orbit
      Integer :: iss  =  0    !  spin-spin
!     Integer :: ioo  =  0    !  orbit-orbit      -> module spline_atomic

! ... min/max 2J:

      Integer :: jmin = -1,  jmax = -1         

! ... other parameters:

      Integer :: nort = 0
      Integer :: meiv = 0            !  number of desirable eigenvalues
      Real(8) :: Emax = 100.d0       !  maximum physical energy

      Real(8) :: eps_ovl = 1.d-8     !  tollerance for overlaps
      Real(8) :: eps_o   = 0.3       !  tollerance for total overlaps
      Real(8) :: eps_d   = 0.001     !  tollerance for diagonal m.e.
      Real(8) :: eps_c   = 1.d-8     !  tollerance for diagonal m.e.
      Real(8) :: eps_det = 1.d-8     !  tollerance for diagonal m.e.

!---------------------------------------------------------------------
! ... parameter from INT_BNK:

!      ! packing base for overlaps in INT_BNK:

!      Integer :: ibf = 16,  ibd = 2**15     

      ! size of block for reading of INT_BNK:

      Integer :: maxnc = 100000

      Integer :: ibi = 2**15

      Integer :: nblock  =   1000         ! number of blocks in c_data       
      Integer :: mblock  =   2000         ! size of blocks
      Integer :: kblock  =   1000         ! max.nb for given case

      Integer :: km = 7

!---------------------------------------------------------------------
! ... correction (fitting) factor for spin-orbital one-electron
! ... electron z-integrals (to be in consistent with the same factor
! ... in bsr_mat calculations):

      Real(8) :: zcorr = 1.d0    

      Integer :: mdiag = 0        

      Integer :: msol = 0

      End Module bsr_ci

