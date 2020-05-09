!======================================================================
      Module bsr_prep
!======================================================================
!     contains main parameters and arrays used in program bsr_prep
!----------------------------------------------------------------------
      Use spline_param; Use spline_grid; Use spline_atomic               
      Use spline_orbitals, p => pbs 
      Use conf_LS; Use orb_LS;  Use target                         
      Use channels, npert1 => npert, ipert1 => ipert, ippert1 => ippert       
      Use channel, only: npert, ipert, ippert

      Implicit none

! ... files:

      Integer, parameter :: ma=80
      Character(ma) :: AF

      Integer :: nut=1;   Character(ma) :: AF_targ = 'target'
      Integer :: nup=2;   Character(ma) :: AF_par  = 'bsr_par'
      Integer :: pri=66;  Character(ma) :: AF_log  = 'bsr_prep.log'

      Integer :: nuc=11;  Character(ma) :: AFC   ! input c-file
      Integer :: nuw=12;  Character(ma) :: AFW   ! input bsw-file
      Integer :: muc=13;  Character(ma) :: BFC   ! output c-file
      Integer :: muw=14;  Character(ma) :: BFW   ! output bsw-file

      Integer :: nua=15;  Character(ma) :: AF_sub = 'target_sub.bsw'
      Integer :: nuo=16;  Character(ma) :: AF_orb = 'target_orb'
 
! ... default value for parameter: 

      Real(8) :: eps_ovl  = 1.d-6
      Real(8) :: eps_phys = 0.25d0
      Real(8) :: eps_sub  = 0.5d0
      Real(8) :: eps_targ = 2.d-8
      Integer :: ii_sub   = 0

      Integer :: LT_min = 0
      Integer :: LT_max = 25
      Integer :: IS_min = -1
      Integer :: IS_max = -1
      Integer :: JJ_min = -1
      Integer :: JJ_max = -1

! ... auxiliary arrays:

      Character(200) :: title, AS
      Character(80), allocatable :: AFK(:)
      Integer, allocatable :: klsp(:), ipt(:) 

      Integer :: kshift = 0

      End Module bsr_prep
