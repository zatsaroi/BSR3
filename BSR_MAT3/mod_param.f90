!----------------------------------------------------------------------
      Module bsr_mat
!----------------------------------------------------------------------
!     main parameters for the BSR_MAT program
!----------------------------------------------------------------------
      Implicit none

! ... input/output files and units:

      Integer, parameter :: ma=20    ! limit for file-name length
      Character(ma) :: AF     

      Integer :: prj =  3; Character(ma) :: AF_prj = 'bsr_mat.log'                            
      Integer :: pri = 66; Character(ma) :: AF_pri = 'mat_log.nnn'

      Integer :: nut = 10; Character(ma) :: AF_tar = 'target'
      Integer :: nup = 11; Character(ma) :: AF_par = 'bsr_par'
      Integer :: nub = 12; Character(ma) :: AF_bnk = 'int_bnk.nnn'
      Integer :: nuc = 13; Character(ma) :: AF_cfg = 'cfg.nnn'
      Integer :: nuw = 14; Character(ma) :: AF_bsw = 'target.bsw'
      Integer :: nui = 15; Character(ma) :: AF_mat = 'bsr_mat.nnn'
      Integer :: nud = 16; Character(ma) :: AF_deb = 'int_mat.nnn'
      Integer :: nuo = 17; Character(ma) :: AF_orb = 'target_orb'
      Integer :: nun = 18; Character(ma) :: AF_new = 'target_new'

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

      Real(8) :: Eps_soo  = 1.0D-03  !  tolerance for two-electron
                                     !  relativistic integrals

! ... near-threshold correction for spin-orbit interaction:

      Integer :: izcorr = 0    ! cut-off for small r
      Real(8) :: zcorr = 1.d0  ! correction for l=1

! ... dimension limits:

      Integer :: mk    =      7   ! max. multipole index
      Integer :: nb    =   2000   ! number of blocks in cmdata       
      Integer :: mb    =   5000   ! size of blocks
      Integer :: kb    =    100   ! max.nb for given case
	                               
! ... tolerence parameters:

      Real(8) :: Eps_c    = 1.0D-10  !  tolerance for coefficients
      Real(8) :: Eps_det  = 1.0D-10  !  tolerance for determinants
      Real(8) :: Eps_ovl  = 1.0D-8   !  tolerance for overlaps

      Real(8) :: S_ovl    = 0.75     !  channel overlap limit 
      Real(8) :: S_pert   = 0.5      !  perturber overlap limit

      Real(8) :: Eps_acf  = 1.0D-5   !  tolerance for asympt. coeff.s
      Real(8) :: Eps_tar  = 1.0D-6   !  tolerance for target energies and overlaps

! ... packing basis for orbitals in the bsr_mat (as i*ibo+j):

      Integer, parameter :: ibo = 2**15

! ... buffer for the integral coefficients:

      Integer :: maxnc = 1000000
      Real(8), allocatable :: CBUF(:)
      Integer, allocatable :: ijtb(:), intb(:), idfb(:)

! ... debug level:

      Integer :: debug = 0

      End Module bsr_mat

