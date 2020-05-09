!====================================================================    
    MODULE spline_integrals
!====================================================================
!
!   contains the B-spline representation of two-electron integral
!   rkb(i,j;i',j') in symmetric or non-symmetric column storage mode:
!
!            rkb(1:ns, 1:ns, 1:2*ks-1, 1:2*ks-1) 
!
!   itype - character (rk, rk1, rk2, ...) which indicates the type
!           of integral and method of calculation for integral storing
!           in the rkb array  at the moment
!
!   krk   - multipole index for the integral
!
!--------------------------------------------------------------------


    IMPLICIT NONE

    INTEGER(4) :: krk = -100

    CHARACTER(3) :: itype='aaa'
    
    REAL(8), ALLOCATABLE, DIMENSION(:,:,:,:) :: rkb
    

    END MODULE spline_integrals


!====================================================================    
      SUBROUTINE allocate_integrals
!====================================================================    
!
! ... allocates space for spline integrals
!
!--------------------------------------------------------------------

      USE spline_param
      USE spline_integrals

      INTEGER :: ierr

      if(allocated(rkb)) Deallocate (rkb)

      ALLOCATE( rkb(ns,ns,2*ks-1,2*ks-1),STAT=ierr)

      rkb = 0.d0

      itype='bbb'
      krk=-100

      END SUBROUTINE allocate_integrals


!====================================================================    
      SUBROUTINE dealloc_integrals
!====================================================================    
!
! ... deallocates arrays in module "spline_integrals"
!
!--------------------------------------------------------------------

      USE spline_integrals

      if(allocated(rkb)) DEALLOCATE(rkb)

      itype='aaa'
      krk = -100

      END SUBROUTINE dealloc_integrals


