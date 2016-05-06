!=====================================================================
      Module bsr_pol     
!=====================================================================
!     common variable and arrays for the program "bsr_pol"
!---------------------------------------------------------------------
      Implicit none

! ... files:

      Integer, parameter :: ma=80;  Character(ma) :: AF
      Integer :: pri = 1;  Character(ma) :: AF_log  = 'bsr_pol.nnn'
      Integer :: nup = 7;  Character(ma) :: AF_par  = 'bsr_par'
      Integer :: nut = 8;  Character(ma) :: AF_tar  = 'target'
      Integer :: nui = 9;  Character(ma) :: AF_int  = 'bsr_mat.nnn'
      Integer :: nub = 11; Character(ma) :: AF_bnd  = 'bound.nnn'
      Integer :: nuw = 12; Character(ma) :: AF_bsw  = 'target.bsw'
      Integer :: nuc = 15; Character(ma) :: AF_cfg  = 'cfg.nnn'
      Integer :: nud = 16; Character(ma) :: AF_dip  = 'dv.nnn'
      Integer :: nur = 18; Character(ma) :: AF_pol  = 'pol.nnn'
      Integer :: nua = 28    ! scratch file for int.matrix
      Integer :: nus = 22    ! scratch file for ovl.matrix
      Integer :: nuq = 23    ! scratch file for orth.states

      Integer :: klsp = 1    ! index of partial wave
      Character(3) :: ALSP ='nnn'

! ... main parameters of calculations:

      Integer :: nhm               ! size of interaction matrix
      Integer :: khm               ! number of solutions
      Integer :: mhm               ! full size of matrix

      Integer :: nort = 0          ! number of orthogonal constraints
      Integer :: nortb = 0
      Integer, allocatable :: iortb(:)

      Real(8), allocatable :: a(:,:)      ! interaction matrix
      Real(8), allocatable :: c(:,:)      ! overlap matrix
      Real(8), allocatable :: d(:)        ! dipole matrix

! ... initial state:

      Integer :: jot0              ! (2L+1) or (2J+1) 
      Real(8) :: E0
      Integer :: IP0
      Character(64) :: Label0
      
! ... tarnsition:

      Character(1) :: ktype = 'E'  ! transition type
      Integer :: kpol              ! multipole index

! ... boundary conditions on solution:

      Integer :: ilzero = 1        ! at r=0
      Integer :: ibzero = 1        ! at r=a

      Real(8) :: zero = 0.d0, one = 1.d0

      End module bsr_pol

