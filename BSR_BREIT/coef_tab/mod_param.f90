!======================================================================
      MODULE  param_br
!======================================================================

      Implicit none
      Save

! ... tolerence for coefficients:

      Real(8) :: Eps_c = 1.d-5      

      Integer(4) :: is_soo = 0  !  Vk -> V'k

! ... parking basis in the data bank (jb should be > nsh):

      Integer(4), Parameter :: jb=10, jb4=jb**4, jb8=jb**8, mk=jb4-1

! ... other parameters:

      Integer(4), Parameter :: isd = 10000  ! supposed number of det.factor
      Integer(4), Parameter :: jsd = 3      ! avarage size of overlap det.
      Integer(4), Parameter :: ibd = 2**15  ! parking basis for det.s

      Integer(4), Parameter :: isf = 100000 ! supposed number of det.factor
      Integer(4), Parameter :: jsf = 3      ! avarage number of dets. in def.
      Integer(4), Parameter :: ibf = 16     ! parking basis for det.factors

      Integer(4), Parameter :: ikt = 1000

!     Integer(4) :: ibc = 2**15       ! parking basis for conf.

      Integer(4), Parameter :: isboef = 100  ! supposed number coef in list boef
      Integer(4), Parameter :: iszoef = 1000
      

! ... space for new coefficients in module COEF_list:

      Integer(4) :: nblock = 10          ! number of blocks
      Integer(4) :: mblock = 200         ! size of block
      Integer(4) :: mtermc = 5000       

! ... memory = nblock*mblock*mtermc*10 ~ 100 Mb

      Integer(4) :: nus = 99          ! file for reallocations

      Integer(4), parameter :: me = 200  ! max.number of electrons


      End MODULE  param_br

