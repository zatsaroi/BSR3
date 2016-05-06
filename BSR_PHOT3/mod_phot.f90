!=====================================================================
      Module bsr_phot      
!=====================================================================
!     Contains common variable and arrays 
!---------------------------------------------------------------------
      Implicit none

! ... files:

      Integer, parameter :: ma = 40
      Integer :: pri = 66; Character(ma) :: AF_log  = 'bsr_phot.log'
      Integer :: inp =  7; Character(ma) :: AF_inp  = 'bsr_phot.inp'
      Integer :: nud = 11; Character(ma) :: AF_d    = 'd.nnn'
      Integer :: nuh = 12; Character(ma) :: AF_h    = 'h.nnn'
      Integer :: nuw = 13; Character(ma) :: AF_w    = 'w.nnn'
      Integer :: nub = 14; Character(ma) :: AF_b    = 'bound.nnn'
      Integer :: ics = 16; Character(ma) :: AF_out  = 'bsr_phot.nnn'
      Integer :: iph = 17; Character(ma) :: AF_ph   = 'photo.nnn'

      Character(ma) :: AF  ! space for names to keep the same size 

! ... arguments:

      Integer, parameter :: nrange = 10  ! number of energy intervals

      Real(8) :: ELOW(nrange),ESTEP(nrange),EHIGH(nrange)

      Real(8) :: e_exp = 0.d0      ! energy shift for initial state
      Real(8) :: AWT = 0.d0        ! nuclear weight

      Integer :: ibug = 0          ! debug level

! ... ASYPCK parameters:

      Integer :: iauto = 1         ! automatic choice of method
      Integer :: mfgi  = 300       ! number of point in outer region
      Real(8) :: AC = 0.001        ! accuracy, < 0.1, but < 0.0001
      Real(8) :: DR = 0.1          ! r1 - r2, < 0.2 

! ... target data:

      Integer :: nelc              ! number of electrons
      Integer :: nz                ! nuclear charge
      Integer :: ion               ! residual charge
      Integer :: ntarg             ! number of target states
      Real(8), allocatable :: Etarg(:) 
      Integer, allocatable :: Ltarg(:),IStarg(:),IPtarg(:) 

! ... scattering channels:

      Integer :: klsp = 0             
      Character(3) :: ALSP = 'nnn'
      Integer :: IL,IS,IP
      Real(8) :: E1

      Integer :: nch               ! number of channels
      Integer :: ncp               ! bound configurations
      Integer :: nhm               ! size of R-matrix basis
      Integer :: khm               ! number of pseudo-states
      Integer :: km                ! max.multipole

      Real(8) :: athreshold=0.0001d0, bthreshold=0.0001d0

      Integer, allocatable :: NCONAT(:),LCH(:)
      Real(8), allocatable :: VALUE(:), ECH(:)
      Real(8), allocatable :: WMAT(:,:) 
      Real(8), allocatable :: CF(:,:,:) 
      Real(8), allocatable :: RMAT(:,:), RMATI(:,:) 
      Real(8), allocatable :: KMAT(:,:), uk(:,:)
      Real(8), allocatable :: ui(:)
 
      Real(8), allocatable :: F(:,:),G(:,:),FP(:,:),GP(:,:) 

! ... weights:

      Integer :: ikm = 0                        
      Integer :: nwt = -1                         
      Real(8), allocatable :: WT(:,:) 
      Real(8), allocatable :: AK(:,:) 
      Integer, allocatable :: ipak(:,:)
      Real(8), allocatable :: WTch(:)
      Real(8) :: CC = 0.d0, CE = 0.d0

! ... RM radius:

      Real(8) :: RA,RB

! ... initial state:

      Integer :: ILI,ISI,IPI
      Real(8) :: EI
      Integer :: ndm 

! ... radiative data:

      Real(8), allocatable :: DKL(:),DKV(:) 
      Real(8), allocatable :: DLr(:),DLi(:),DVr(:),DVi(:)
      Real(8), allocatable :: SL(:),SV(:)
      Real(8) :: SLP,SVP

! ... work arrays:

      Real(8), allocatable :: AA(:,:),BB(:,:),FF(:,:),FFF(:,:),FFP(:,:) 
      Real(8), allocatable :: v(:), eps(:) 

      End Module bsr_phot

