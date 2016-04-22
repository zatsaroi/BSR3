!======================================================================
      MODULE detf
!======================================================================
!
!     contains information about overlap factors
!
!----------------------------------------------------------------------
!
! KPF(i)  - number (kd) of determinants in i-th overlap factor
!
! IPF(i)  - pointer (ip) on the i-th overlap factor in the NPF array
!
! NPF(ip+1:ip+kd) - i-th overlap factors as a list of pointers 
!                   on individual determinants (np) and its power (ne): 
!                   npf(i) = np * ibf  +  ne
!----------------------------------------------------------------------      

      IMPLICIT NONE
      SAVE
     
      INTEGER(4) :: mdf =   0  ! maximum number of overlap factors 
      INTEGER(4) :: ndf =   0  ! current number of overlap factors 
      INTEGER(4) :: kdf =   0  ! sum of all overlap factor dimensions  
      INTEGER(4) :: mkdf=   0  ! maximum kdf
      INTEGER(4) :: ibf = 128  ! pack basis 
      
      INTEGER(4), DIMENSION(:), ALLOCATABLE :: KPF,IPF,NPF

      END MODULE detf



!===================================================================
      Subroutine Allocate_detf(m,mm,mbf)
!===================================================================

      USE detf
      
      IMPLICIT NONE
      Integer(4),INTENT(in) :: m,mm,mbf
      
      if(m.eq.0) then
       Deallocate(KPF,IPF,NPF)
       mdf = 0; mkdf = 0; ndf = 0; kdf = 0
      else          
       if(allocated(KPF)) Deallocate(KPF,IPF,NPF) 
       mdf = m; mkdf = mm 
       Allocate(KPF(mdf),IPF(mdf),NPF(mkdf))
       if(mbf.ne.0) ibf = mbf
      end if

      END Subroutine Allocate_detf


