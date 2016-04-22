!======================================================================
      MODULE term_exp
!======================================================================
!
!     Containes the term-dependent coefficients of det.expansion
!     of two given conf.symmetries under consideration.
!
!----------------------------------------------------------------------

      Implicit none
      Save

      Integer(4) :: ILT1,ILT2   ! total L
      Integer(4) :: IST1,IST2   ! total S
      Integer(4) :: MLT,MST     ! total ML and MS
      Integer(4) :: kd1,kd2     ! det. under consideration
      
! ... lists of ang.symmetries (1:kt)

      Integer(4) :: kt1,kt2  
      Integer(4), Allocatable, Dimension(:) :: IP_kt1,IP_kt2   
 
! ... pointer on non-zero determinants (1:ne,1:kdt)

      Integer(4) :: kdt1,kdt2  
      Integer(4), Allocatable, Dimension(:,:) :: IP_det1,IP_det2   

! ... term-dependent det. expension coefficients (1:kt,1:kdt)

      Real(8),  Allocatable, Dimension(:,:) :: C_det1, C_det2

      End MODULE term_exp



