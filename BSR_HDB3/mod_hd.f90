!=====================================================================
      Module bsr_hd      
!=====================================================================
!     Contains common variable and arrays 
!---------------------------------------------------------------------
      use blacs, only: dlen_

      Implicit none

! ... files:

      Integer, parameter :: ma = 40;  Character(ma) :: AF

      Integer :: pri = 1;  Character(ma) :: AF_nnn  = 'bsr_hd.nnn'
      Integer :: nup = 7;  Character(ma) :: AF_par  = 'bsr_par'
      Integer :: nut = 8;  Character(ma) :: AF_tar  = 'target'
      Integer :: nui = 9;  Character(ma) :: AF_int  = 'bsr_mat.nnn'
      Integer :: nub = 11; Character(ma) :: AF_bound= 'bound.nnn'
      Integer :: nuh = 13; Character(ma) :: AF_h    = 'h.nnn'
      Integer :: nuw = 14; Character(ma) :: AF_w    = 'w.nnn'
      Integer :: nuc = 15; Character(ma) :: AF_cfg  = 'cfg.nnn'
      Integer :: nue = 16; Character(ma) :: AF_exp  = 'thresholds'
      Integer :: nur = 17; Character(ma) :: AF_rsol = 'rsol.nnn'

      Integer :: nus = 99; Character(ma) :: AF_knot = 'knot.dat'

! ... arguments:

      Integer :: itype  = 0        ! type of calculations
      Integer :: it_max = 0        ! max. threshold
      Integer :: msol   =10        ! max.number of solutions
      Integer :: ilzero = 1        ! # B-soplines deleted at r=0
      Integer :: ibzero = 1        ! # B-soplines deleted at r=a
      Real(8) :: Emax = 0.d0       ! max. energy
      Real(8) :: Emin = 0.d0       ! min. energy
      Real(8) :: Edmax = 0.d0      ! max. energy for diagonalization
      Real(8) :: Edmin = 0.d0      ! min. energy for diagonalization
      Real(8) :: Egap = 0.001d0    ! tolerance for zero-energy solutions

      Integer :: klsp=0            ! # of partial wave
      Integer :: klsp1=1,klsp2=1   ! range of partial waves
      Character(3) :: ALSP ='nnn'

! ... main parameters of calculations:

      Integer :: nhm               ! size of interaction matrix
      Integer :: khm               ! number of solutions
      Integer :: kch               ! number of channels
      Integer :: kcp               ! # configurations in perturber
      Integer :: ksol = 0          ! # channel eigenvalues 
      Integer :: diag_ovl          ! flag for nonzero ovl.blocks

! ... main global arrays:

      Real(8), allocatable :: a(:,:)      ! interaction matrix
      Integer              :: desca(dlen_) 
      Real(8), allocatable :: b(:,:)      ! overlap matrices
      Integer              :: descb(dlen_) 
      Real(8), allocatable :: z(:,:)      ! solutions        
      Integer              :: descz(dlen_) 
      Real(8), allocatable :: v(:)        ! solution vector
      Integer              :: descv(dlen_) 

      Real(8), allocatable :: eval(:)     ! eigenvalues

      Real(8), allocatable :: bb(:,:)     ! new basis
      Real(8), allocatable :: bval(:)     ! basis eigenvalues
      Integer, allocatable :: ipsol(:)    ! pointer on channel blocks 
                                          ! in new basis 
      Integer, allocatable :: isol(:)     ! pointer on main configuration
    
      Real(8), allocatable :: WMAT(:,:)   ! surface amplitudes

      Real(8), allocatable :: CF(:,:,:)   ! asymptotic coefficients
      Integer              :: lamax

      Real(8) :: RA  ! R-matrix radius

! ... working arrays

      Real(8), allocatable :: add(:,:)
      Integer              :: descadd(dlen_) 
      Real(8), allocatable :: adp(:,:)
      Integer              :: descadp(dlen_) 
      Real(8), allocatable :: adv(:)
      Integer              :: descadv(dlen_) 

      Real(8), allocatable :: cc(:,:),w(:)
      Integer              :: parms(5)

! ... exp.energies: 
    
      Logical :: EXP 
      Integer :: iexp = 0, iiexp = 0
      Real(8), Allocatable :: E_exp(:)
      Integer, Allocatable :: ip_exp(:)
      Real(8) :: au_eV, au_cm
      Character(2) :: unit = 'au'

! ... pointer on the inclusion of mass-velocity term 
! ... in first order

      Integer :: jmvc = 0

! ... additonal output: 

      Integer :: iwt = -1       ! print channel weights
      Real(8) :: cwt = -0.01d0  ! cut of for weights' printing

! ... miscellaneous:

      Real(8), parameter :: zero = 0.d0, one = 1.d0
      Integer :: fail = 0
      Real(8) :: eps_o = 0.25d0, eps_d = 1.d-5

! ... debug and timing:

      Integer :: debug = 0   
      Real(8) :: t0, t1, ctime, etime
      Integer :: c0, c1, crate

      End module bsr_hd


!======================================================================
      Subroutine Get_t0
!======================================================================
      Use bsr_hd
      Use blacs

      call blacs_barrier (ctxt, 'All')
      call cpu_time (t0)

      End Subroutine Get_t0

!======================================================================
      Subroutine Get_t1(name)
!======================================================================
      Use bsr_hd
      Use blacs
      Character(*) :: name

      call blacs_barrier (ctxt, 'All')
      call cpu_time (t1)

      if(io_processor) then
       ctime = (t1-t0)/60
       write (pri,'(/a,T20,a,f10.2,a)') trim(name),' CPU     = ', ctime, ' min.'
       write (*  ,'(/a,T20,a,f10.2,a)') trim(name),' CPU     = ', ctime, ' min.'
      end if

      End Subroutine Get_t1

