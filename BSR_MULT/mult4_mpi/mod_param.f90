!======================================================================
      MODULE  mult_par
!======================================================================
! ... contains main parmeters and arrays for the mult progtm
!----------------------------------------------------------------------    
      Implicit none

! ... main files:

      Integer, parameter :: ma=80
      Integer :: pri=66;   Character(ma) :: AF_p = 'mult.log'
      Integer :: in1=1;    Character(ma) :: AF1 
      Integer :: in2=2;    Character(ma) :: AF2 
      Integer :: nub=3;    Character(ma) :: AF_b = 'mult_bnk'
      Integer :: nur=4;    Character(ma) :: AF_r = 'mult_res'

! ... scratch files:

      Integer :: nui = 10  ! for intermediate results
      Integer :: nus = 11  ! for re-allocations   
      Integer :: nud = 12  ! for det.expansions
      Integer :: nua = 13  ! for accumulation of data

! ... flags:

      Integer :: new     ! flag for new calculations
      Integer :: icalc   ! flag for need of calculations

! ... tolerence for coefficients:

      Real(8) :: eps_c = 1.d-7      

! ... initial (supposed) number of coef.s:

      Integer, Parameter :: iszoef = 1000
      Integer, Parameter :: iscoef = 10000

      Integer, parameter :: me = 200  ! max.number of electrons

! ... switch to MLTPOL representation:

      Integer :: mltpol = 0 

!      Character(5) :: move = 'move '   !  for WINDOWS
      Character(5) :: move = 'mv   '   !  for INIX

! ...  multipole transition under consideration:

      Integer :: kpol = 1, qpol, mpol, spol
      Character(1) :: ktype = 'E'

! ... normalization constants:

      Real(8) :: CA, CB, CNA, CNB

      Real(8), allocatable :: CT_oper(:)
      Integer, allocatable :: JT_oper(:)
      Integer, allocatable :: JD_oper(:)

! ... list of total LS:

      Integer, allocatable :: ILT_ic(:), IST_ic(:)

! ... MPI insert:

      Integer :: nprocs, myid, ierr
      Integer, allocatable :: ip_proc(:)
      Integer :: debug=0

! ... restrictions:

      Integer :: mktkdt=1000000

      Real(8) :: t0,t1,t2,t3,t4,tt, time_limit=0.d0

      End Module mult_par

