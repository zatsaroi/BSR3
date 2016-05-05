!======================================================================
      Module term_exp
!======================================================================
!     Containes the term-dependent coefficients of det.expansion
!     of two given conf.symmetries under consideration.
!----------------------------------------------------------------------
      Implicit none

      Integer :: ILT1,ILT2   ! total L
      Integer :: IST1,IST2   ! total S
      Integer :: MLT,MST     ! total ML and MS
      Integer :: kd1,kd2     ! det. under consideration
      
! ... lists of ang.symmetries (1:kt)

      Integer :: kt1,kt2  
      Integer, allocatable :: IP_kt1(:),IP_kt2(:)   
 
! ... pointer on non-zero determinants (1:ne,1:kdt)

      Integer :: kdt1,kdt2  
      Integer, allocatable :: IM_det1(:,:),IM_det2(:,:)   
      Integer, allocatable :: IS_det1(:,:),IS_det2(:,:)   

! ... term-dependent det. expension coefficients (1:kt,1:kdt)

      Real(8),  allocatable :: C_det1(:,:), C_det2(:,:) 

! ... number of different (CONFIG or ML,MS) cases :

      Integer :: ic_case   

      End Module term_exp



