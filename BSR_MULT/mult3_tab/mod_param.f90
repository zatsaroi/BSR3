!======================================================================
      MODULE  mult_tab
!======================================================================

      Implicit none

! ... main files:

      Integer, parameter :: ma=80
      Integer :: in1=11;    Character(ma) :: AF1 = 'initial.c' 
      Integer :: in2=12;    Character(ma) :: AF2 = 'final.c'
      Integer :: nub=13;    Character(ma) :: AF_bnk = 'mult_bnk.E1'
      Integer :: out=14;    Character(ma) :: AF_tab = 'mult_tab'

      Integer :: ic=0, jc=0, ips,jps, ipc,jpc

      Character(1)  ::  ktype
      Integer :: kpol


! ... tolerence for coefficients:

      Real(8) :: Eps_c = 1.d-7      

      Integer, parameter :: me = 200  ! max.number of electrons

      End MODULE  mult_tab

