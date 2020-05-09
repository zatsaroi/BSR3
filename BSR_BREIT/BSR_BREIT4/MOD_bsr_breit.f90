!======================================================================
      Module bsr_breit
!======================================================================
!     main parameters and files:
!----------------------------------------------------------------------
      Use param_LS

      Implicit none

! ... files: 

!     AF  -  standard (default) names
!     BF  -  names with indication of partial wave number

      Integer, parameter :: ma=80

      Integer :: iname=0; Character(ma) :: name = ' '
     
      Integer :: pri=10;  Character(ma) :: AF_p = 'bsr_breit.log'
      Integer :: nuc=11;  Character(ma) :: AF_c = 'cfg.inp'
                          Character(ma) :: BF_c = 'cfg.nnn'
      Integer :: nub=12;  Character(ma) :: AF_b = 'int_inf'
                          Character(ma) :: BF_b = 'int_inf.nnn'
      Integer :: nur=13;  Character(ma) :: AF_r = 'int_int'
                          Character(ma) :: BF_r = 'int_int.nnn'
      Integer :: nud=14;  Character(ma) :: AF_d = 'det_exp'
                          Character(ma) :: BF_d = 'det_exp.nnn'

!----------------------------------------------------------------------
! ... main parameters:
!----------------------------------------------------------------------

      Real(8) :: eps_c = 1.d-8          ! tolerence for coefficients   
      Real(8) :: eps_soo = 1.d-8         ! tolerence for coefficients   

      Integer :: mk = 7                 ! maximum multipole index

      Integer :: new        ! pointer on the previous calculation  
      Integer :: icalc      ! pointer to the need of new calculations

!      Integer :: mc = 10000000

! ... range of partial waves:

      Integer :: klsp=0, klsp1=0, klsp2=0

! ... restrictions for det_expansion:

      Integer :: mkt = 10000
      Integer :: mkdt = 1000000

! ... information data:

      Integer(8) :: nc_new = 0
      Real(8) :: adet,adef       
      Real(8) :: t0,t1,t2,t3,t4, time_limit = 0.d0

! ... bufer:

      Integer :: mbuf=1000000, nbuf
      Integer, allocatable :: ibuf1(:),ibuf2(:),ibuf3(:),ibuf4(:)
      Real(8), allocatable :: Cbuf(:)
      
! ... MPI insert:

      Integer :: nprocs, myid, ierr
      Integer, allocatable :: ip_proc(:)
      Integer :: debug=0

      End Module bsr_breit


