!=====================================================================
      MODULE conf_LS_1
!=====================================================================
!     contains description of two configurations
!---------------------------------------------------------------------
      Implicit none

      Integer, parameter :: msh = 25  ! max. number of shells behind core

      Integer :: ne = 0               ! number of electrons

! ... description of 1 configuration state:  

      Integer :: no, Ltotal, Stotal, Jtotal, Ptotal, &
                 nn(msh),ln(msh),iq(msh),kn(msh),LS(msh,5)

! ... Storing configurations as character strings:

      Character(200) :: CONFIG, COUPLE

! ... closed shells core 

      Real(8) :: Ecore = 0.d0 
      Integer :: ncore = 0
      Character(200) :: core=' '

! ... determinant expansions:

      Integer :: MLT, MST, kdt
      Real(8), allocatable :: Cdet(:)
      Integer, allocatable :: MLdet(:,:), MSdet(:,:) 

! ... inter-electron interaction coefficients

      Integer :: kmax 
      Real(8), allocatable :: coefs(:,:,:)

      END MODULE conf_LS_1

