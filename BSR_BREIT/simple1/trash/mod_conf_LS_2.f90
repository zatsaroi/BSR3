!=====================================================================
      MODULE conf_LS_2
!=====================================================================
!     contains description of two configurations
!---------------------------------------------------------------------
      Implicit none

      Integer, parameter :: msh = 25  ! max. number of shells behind core

      Integer :: ne = 0               !  number of electrons

! ... description of 1 configuration state:  

      Integer :: no, Ltotal, Stotal, Jtotal, Ptotal, &
                 nn(msh),ln(msh),iq(msh),LS(msh,5)

      Integer :: no1, Ltotal1, Stotal1, Jtotal1, Ptotal1, &
                 nn1(msh),ln1(msh),iq1(msh),LS1(msh,5)

      Integer :: no2, Ltotal2, Stotal2, Jtotal2, Ptotal2, &
                 nn2(msh),ln2(msh),iq2(msh),LS2(msh,5)

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

      Integer :: MLT1, MST1, kdt1
      Real(8), allocatable :: Cdet1(:)
      Integer, allocatable :: MLdet1(:,:), MSdet1(:,:) 

      Integer :: MLT2, MST2, kdt2
      Real(8), allocatable :: Cdet2(:)
      Integer, allocatable :: MLdet2(:,:), MSdet2(:,:) 

! ... inter-electron interaction coefficients

      Integer :: ncoef = 0, mcoef = 0
      Real(8), allocatable :: coef(:)
      Integer, allocatable :: kk(:),i1(:),i2(:),i3(:),i4(:),ie(:)

      END MODULE conf_LS_2

