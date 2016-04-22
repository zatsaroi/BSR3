!===================================================================
      MODULE DETS
!===================================================================
!
!     contains information about overlap determinants {<i|j>}
!
!-------------------------------------------------------------------
!
! KPD(i)  - kd*ibd+iext, where kd - size of i-th determinant
!                                   iext - its power
!
! IPD(i)  - pointer (ip) of i-th determinants in the NPD array
!
! NPD(ip+1:ip+kd) - contains information about orbitals in overlap
!                   determinant as  ii * ibd + jj, 
!                   where ii - pointer on row orbitals
!                         jj - pointer on column orbitals
!
!-------------------------------------------------------------------     

      IMPLICIT NONE
      SAVE
     
      INTEGER(4) :: mdet  =   0   !  max. number of determinants
      INTEGER(4) :: ndet  =   0   !  current number of determinants
      INTEGER(4) :: kdet  =   0   !  sum of all det. dimensions      
      INTEGER(4) :: mkdet =   0   !  maximum kdet
      INTEGER(4) :: ibd   = 128   !  pack basis 
	
      INTEGER(4), DIMENSION(:), ALLOCATABLE :: KPD, IPD, NPD

      END MODULE dets



!===================================================================
      Subroutine Allocate_dets(m,mm,mbn)
!===================================================================

      USE dets

      IMPLICIT NONE
      Integer(4), INTENT(in) :: m,mm,mbn
	
      if(m.eq.0) then
       Deallocate(KPD,IPD,NPD)
       mdet = 0; mkdet = 0; ndet=0; kdet =0
      else          
       if(allocated(KPD)) Deallocate(KPD,IPD,NPD) 
       mdet = m; mkdet = mm
       Allocate(KPD(mdet),IPD(mdet),NPD(mkdet))
       if(mbn.ne.0) ibd = mbn
      end if

      END Subroutine Allocate_dets
