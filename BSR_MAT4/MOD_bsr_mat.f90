!----------------------------------------------------------------------
      Module bsr_mat
!----------------------------------------------------------------------
!     main parameters for the BSR_MAT program
!----------------------------------------------------------------------
      Use zconst, only: c_au, time_au  
      Use target
      Use channel

      Use spline_param
      Use spline_atomic
      Use spline_grid
      Use spline_galerkin
      Use spline_orbitals

      Implicit none

! ... input/output files and units:

      Integer, parameter :: ma=40    ! limit for file-name length
      Character(ma) :: AF, BF     

      Integer :: prj =  3; Character(ma) :: AF_prj = 'bsr_mat.log'                            
      Integer :: pri = 66; Character(ma) :: AF_pri = 'mat_log.nnn'

      Integer :: nut = 10; Character(ma) :: AF_tar = 'target'
      Integer :: nup = 11; Character(ma) :: AF_par = 'bsr_par'
      Integer :: nub = 12; Character(ma) :: AF_inf = 'int_inf.nnn'
                           Character(ma) :: AF_int = 'int_int.nnn'
                           Character(ma) :: AF_lst = 'int_list.nnn'
      Integer :: nuc = 13; Character(ma) :: AF_cfg = 'cfg.nnn'
      Integer :: nuw = 14; Character(ma) :: AF_bsw = 'target.bsw'
      Integer :: nui = 15; Character(ma) :: AF_mat = 'bsr_mat.nnn'
      Integer :: nuj = 16; 
      Integer :: nud = 17; Character(ma) :: AF_deb = 'int_mat.nnn'
      Integer :: nuo = 18; Character(ma) :: AF_orb = 'target_orb'
      Integer :: nun = 19; Character(ma) :: AF_new = 'target_new'
      Integer :: nua = 20  ! just to quick use

! ... range of partial waves:

      Integer :: klsp=0, klsp1=1, klsp2=1      
      Character(3) :: ALSP

! ... treat the target states as othonormalized w.f.: 

      Integer :: iitar  =  0      

! ... relativistic corrections:

      Integer :: mrel  =   0      
      Integer :: mso   =   0
      Integer :: msoo  =   0
      Integer :: mss   =   0
      Integer :: moo   =   0

      Integer :: mlso  =   5  ! max. l for so-interaction

      Integer :: imvc  =  -2  ! mode for mass-velocity correction
      Integer :: nmvc  =   0  ! n-value for mass-velocity correction

      Real(8) :: eps_soo  = 1.0D-05  !  tolerance for two-electron
                                     !  relativistic integrals

! ... near-threshold correction for spin-orbit interaction:

      Integer :: izcorr = 0    ! cut-off for small r
      Real(8) :: zcorr = 1.d0  ! correction for l=1

! ... dimension limits:

      Integer :: mk       =      7   ! max. multipole index
      Integer :: nblock   =   2000   ! number of blocks in c_data       
      Integer :: mblock   =   2000   ! size of blocks
      Integer :: kblock   =   100    ! max.nb for the given type of integrals
	                               
! ... tolerence parameters:

      Real(8) :: eps_c    = 1.0D-10  !  tolerance for coefficients
      Real(8) :: eps_det  = 1.0D-10  !  tolerance for determinants
      Real(8) :: eps_ovl  = 1.0D-8   !  tolerance for overlaps

      Real(8) :: S_ovl    = 0.75     !  channel overlap limit 
      Real(8) :: S_pert   = 0.5      !  perturber overlap limit

      Real(8) :: Eps_acf  = 1.0D-5   !  tolerance for asympt. coeff.s
      Real(8) :: Eps_tar  = 1.0D-6   !  tolerance for target energies and overlaps

! ... packing basis for orbitals in the bsr_mat (as i*ibo+j):

      Integer, parameter :: ibo = 2**15

! ... buffer for the integral coefficients:

      Integer :: maxnc = 10000000, ncbuf
      Real(8), allocatable :: CBUF(:)
      Integer, allocatable :: itb(:),jtb(:),intb(:), idfb(:)

! ... debug level:

      Integer :: debug = 0
      Integer :: pri_f = 0
      Integer :: pri_ac= 0

      Integer :: ilcorr = 0

! ... Calculations in parts:

      Integer :: myid = 0, ierr = 0, nprocs = 1

      Integer, allocatable :: ip_channel(:)

      Real(8), allocatable :: hcc(:,:,:),hcb(:,:),hbb(:)
      Real(8), allocatable :: ACF(:,:,:)
      Real(8), allocatable :: htarg(:),otarg(:)
      Real(8), allocatable :: x(:,:)
      Integer, allocatable :: icc(:,:), icb(:,:), ibb(:,:),imycase(:,:)

      Integer :: iicc, iicb, iibb  !  dimensions

! ... structure of data:

      Integer, parameter :: ibi = 2**15 ! packing basis for orbitals

      Integer, parameter :: ncase = 11 ! number of different cases
      Integer, parameter :: mtype =  9 ! number of different structures

      Integer :: icase =  0     !  current case    

! ... IJCASE  -  correspondence between icase and itype

      Integer IJCASE(mtype)  /2,1,2,1,2,1,2,3,4/

! ... AINT   -  designation for different type of integrals

      Character AINT(ncase)/' ',' ','T','M','R','L','Z','N','V','N','O'/

      Integer :: korth = 0
      Integer, allocatable :: ip_orth_chan(:), ip_orth_orb(:),jp_orth_orb(:)
      ! ip_orth_orb  - list of orbitals paticipated in the orth. conditions
      ! ip_orth_chan - pointer to the last orbital for given channel
      ! jp_orth_orb  - just used for reallocations

      Integer :: check_target = 1

      Integer :: ipert_ch = 1

      Integer :: exch_mode = 0

      Integer :: bp_mode = 0

      Integer :: mode = 0

      Real(8) :: time0, time_delay = -1 
      Integer :: interrupt = 0, intercase = 0, max_nbuf = 0

      End Module bsr_mat

