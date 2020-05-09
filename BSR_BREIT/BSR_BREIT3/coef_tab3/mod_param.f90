!======================================================================
      MODULE  param_br
!======================================================================

      Use symt_list_LS, only: noper, it_oper 

      Implicit none

!      Integer, parameter :: noper = 7 

      Character(7) :: oper = '1110000'
      Integer(1) ioper(noper)/1,1,1,0,0,0,0/, joper(noper)

!     Operator(1)   -   overlaps
!     Operator(2)   -   kinatic energy
!     Operator(3)   -   two-electron electrostatic
!     Operator(4)   -   spin-orbit
!     Operator(5)   -   spin-other-orbit
!     Operator(6)   -   spin-spin
!     Operator(7)   -   orbit-orbit

!      Integer(1), Allocatable :: IT_oper(:,:)

! ... tolerence for coefficients:

      Real(8) :: Eps_c = 1.d-5      

      Integer :: is_soo = 0  !  Vk -> V'k

      Integer, parameter :: me = 200  ! max.number of electrons

      Integer :: nus = 99

      End MODULE  param_br

