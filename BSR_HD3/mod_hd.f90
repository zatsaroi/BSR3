!=====================================================================
      Module bsr_hd      
!=====================================================================
!     Contains common variable and arrays for program bsr_hd
!---------------------------------------------------------------------
      Implicit none

! ... files:

      Integer, parameter :: ma = 40;  Character(ma) :: AF

      Integer :: pri =66;  Character(ma) :: AF_log  = 'bsr_hd.nnn'
      Integer :: nup = 7;  Character(ma) :: AF_par  = 'bsr_par'
      Integer :: nut = 8;  Character(ma) :: AF_tar  = 'target'
      Integer :: nui = 9;  Character(ma) :: AF_int  = 'bsr_mat.nnn'
      Integer :: nub = 11; Character(ma) :: AF_b    = 'bound.nnn'
      Integer :: nuu = 12; Character(ma) :: AF_ub   = 'ubound.nnn'
      Integer :: nuh = 13; Character(ma) :: AF_h    = 'h.nnn'
      Integer :: nuw = 14; Character(ma) :: AF_w    = 'w.nnn'
      Integer :: nuc = 15; Character(ma) :: AF_cfg  = 'cfg.nnn'
      Integer :: nue = 16; Character(ma) :: AF_exp  = 'thresholds'
      Integer :: nur = 17; Character(ma) :: AF_rsol = 'rsol.nnn'

! ... arguments:

      Integer :: itype  = 0        ! type of calculations
      Integer :: IT_max = 0        ! max. threshold
      Integer :: msol = 0          ! max.number of solutions
      Integer :: ksol = 0
      Real(8) :: Emax = 0.d0       ! max energy
      Real(8) :: Emin = 0.d0       ! max energy
      Real(8) :: Edmax = 0.d0      ! max. energy for diagonalization
      Real(8) :: Edmin = 0.d0      ! min. energy for diagonalization
      Real(8) :: Egap = 0.001      ! 
      Integer :: ilzero = 1        ! more restrictions
      Integer :: ibzero = 1        ! more restrictions

      Character(3) :: ALSP ='000'

! ... main parameters of calculations:

      Integer :: klsp=0            ! # of partial wave
      Integer :: klsp1=1,klsp2=1   ! range of partial waves
      Integer :: mhm               ! dimension for allocation
      Integer :: nhm               ! size of interaction matrix
      Integer :: khm               ! number of solutions
      Integer :: kch               ! number of channels
      Integer :: kcp               ! perturber


      Real(8), allocatable :: a(:,:)  ! interaction matrix
      Real(8), allocatable :: c(:,:)  ! overlap matrices

      Real(8), allocatable :: eval(:) ! eigenvalues
      Real(8), allocatable :: v(:)    ! solution vector
      Real(8), allocatable :: bb(:,:) ! new basis

      Integer, allocatable :: ipsol(:) ! pointer on new basis 
      Integer, allocatable :: isol(:)  ! pointer on main configuration
    
      Real(8), allocatable :: WMAT(:,:)  ! surface amplitudes

      Real(8), allocatable :: CF(:,:,:)  ! asymptotic coefficients
      Integer :: lamax

! ... target energies and overlaps:

      Real(8), allocatable :: th(:), to(:)

      Real(8) :: RA  ! R-matrix radius

! ... exp.energies: 
    
      Logical :: EXP 
      Integer :: iexp = 0, iiexp = 0
      Real(8), allocatable :: E_exp(:)
      Integer, allocatable :: ip_exp(:)
      Real(8) :: au_eV, au_cm
      Character(2) :: unit = 'au'

! ... pointer on the inclusion of mass-velocity term 
! ... in first order

      Integer :: jmvc = 0

! ... additonal output: 

      Integer :: iwt = -1    ! print channel weights
      Real(8) :: cwt = -0.01 ! cut of for weights' printing

      Real(8) :: eps_o = 0.5   ! tolerance for overlaps
      Real(8) :: eps_d = 1.d-5 ! tolerance for diag.elemets

      Integer :: debug = 0   ! print debug information

      End Module bsr_hd

