!======================================================================
      Module bsr_breit
!======================================================================
!     main parameters and files:
!----------------------------------------------------------------------
      Use param_LS

      Implicit none

! ... main files: 

!     AF  -  standard (default) names
!     BF  -  names with indication of partial wave number

      Integer, parameter :: ma=80

      Integer :: iname=0; Character(ma) :: name = ' '

      Integer :: pri=66;  Character(ma) :: AF_p = 'bsr_bmat.log'
      Integer :: nuc= 1;  Character(ma) :: AF_c = 'cfg.inp'
                          Character(ma) :: BF_c = 'cfg.nnn'
                          Character(ma) :: CF_c = 'name.c'
      Integer :: nub= 2;  Character(ma) :: AF_b = 'int_list'
                          Character(ma) :: BF_b = 'int_list.nnn'

      Integer :: nud=11;  Character(ma) :: AF_d = 'det_exp'
                          Character(ma) :: BF_d = 'det_exp.nnn'
      Integer :: nur=12;  Character(ma) :: AF_r = 'det_done'
                          Character(ma) :: BF_r = 'det_done.nnn'

      Integer :: nuf=25;  Character(ma) :: AF_f = 'fail_conf'
      Integer :: nux=26;  Character(ma) :: AF_a = 'det_log'

!----------------------------------------------------------------------
! ... main parameters:
!----------------------------------------------------------------------


      Real(8) :: eps_c  = 1.d-10         ! tolerence for coefficients   
      Real(8) :: eps_so = 0.1            ! tolerence for so-interaction   

      Integer :: mk = 7                  ! maximum multipole index

      Integer :: new        ! pointer on the previous calculation  
      Integer :: icalc      ! pointer to the need of new calculations
      Integer :: fail = 0   

! ... range of partial waves:

      Integer :: klsp=0, klsp1=0, klsp2=0

!----------------------------------------------------------------------
!     the matrix elements under consideration:
!----------------------------------------------------------------------
!
!     noper     - number of different operators
!     ioper(:)  - pointer to required operators
!     joper(:)  - pointer to required operator for given configurations
!
!     IT_oper(:,:) - pointer on the done calculation for specific
!                    operators and given terms
!     JT_oper(:,:) - pointer on the required operators for given
!                    subset of term between two configurations
!
!-----------------------------------------------------------------------
!      Integer, parameter :: noper=7      
!      Integer ioper(noper)/1,1,1,0,0,0,0/, joper(noper)
!      Integer koper(noper)     !  MPI copy
!      Real(8) :: coper(noper)

!     Operator(1)   -   overlaps
!     Operator(2)   -   kinatic energy
!     Operator(3)   -   two-electron electrostatic
!     Operator(4)   -   spin-orbit
!     Operator(5)   -   spin-other-orbit
!     Operator(6)   -   spin-spin
!     Operator(7)   -   orbit-orbit

!      Integer, allocatable :: JT_oper(:,:)
!      Real(8), allocatable :: CT_oper(:,:)
!      Integer, allocatable :: JD_oper(:,:)  ! MPI copy

! ... MPI insert:

      Integer :: nprocs=0, myid=0, ierr=0
      Integer, allocatable :: ip_proc(:)
      Integer :: debug=0

! ... restrictions:

      Integer :: mkt = 100, mkdt = 10000

      Real(8) :: t0,t1,t2,t3,t4,tt,time

      Real(8) :: time_limit = -1.d0

      End Module bsr_breit


