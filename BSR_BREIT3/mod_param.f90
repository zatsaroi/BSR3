!======================================================================
      Module bsr_breit
!======================================================================
!     main parameters and files:
!----------------------------------------------------------------------

      Implicit none

! ... files: 

!     AF  -  standard (default) names
!     BF  -  names with indication of partial wave number

      Integer, parameter :: ma=80
      Integer :: pri=66; Character(ma) :: AF_p = 'bsr_breit.log'
      Integer :: nuc=1;  Character(ma) :: AF_c = 'cfg.inp'
                         Character(ma) :: BF_c = 'cfg.nnn'
      Integer :: nub=2;  Character(ma) :: AF_b = 'int_bnk'
                         Character(ma) :: BF_b = 'int_bnk.nnn'
      Integer :: nur=3;  Character(ma) :: AF_r = 'int_res'
                         Character(ma) :: BF_r = 'int_res.nnn'
! ... scratch files:

      Integer :: nui=11  ! intermediate results
      Integer :: nud=12  ! for det. expansions
      Integer :: nua=13  ! for accumulation of data

!----------------------------------------------------------------------
! ... main parameters:
!----------------------------------------------------------------------

! ... tolerence for coefficients:

      Real(8) :: Eps_c = 1.d-8      

! ... packing basis for configurations:

      Integer, parameter :: ibc = 2**15 

! ... maximum multipole index:

      Integer :: mk = 7 

      Integer :: new     ! pointer on the previous calculation  
      Integer :: icalc   ! pointer for need of new calculations

! ... range of partial waves:

      Integer :: klsp1=0, klsp2=0

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

      Integer, parameter :: noper=7      
      Integer ioper(noper)/1,1,1,0,0,0,0/, joper(noper)
      Real(8) :: coper(noper)

!     Operator(1)   -   overlaps
!     Operator(2)   -   kinatic energy
!     Operator(3)   -   two-electron electrostatic
!     Operator(4)   -   spin-orbit
!     Operator(5)   -   spin-other-orbit
!     Operator(6)   -   spin-spin
!     Operator(7)   -   orbit-orbit

      Integer, allocatable :: JT_oper(:,:)
      Real(8), allocatable :: CT_oper(:,:)

! ... restrictions:

      Integer :: mktkdt=200000

      End MODULE bsr_breit


