!====================================================================
    MODULE coeffs
!====================================================================
!
!   containes the description and type of possible coefficients
!
!--------------------------------------------------------------------

    IMPLICIT NONE 
    SAVE
    
    INTEGER(4) :: ncase = 11     !  number of different cases
    INTEGER(4) :: icase =  0     !  case under consideration    
    
    ! ... number of coefficients for given case

    INTEGER(4),DIMENSION(:),ALLOCATABLE :: NCOEF,KCOEF   

    INTEGER(4) :: ntype =  9     !  number of different structures
    INTEGER(4) :: itype =  0     !  structure under consideration    
    INTEGER(4) :: kpol  =  0     !  current multipole index                   
    
    ! ... number of coefficients for given multipole and case

    INTEGER(4) :: mkpol =100     !  max. multipole                   
    INTEGER(4),DIMENSION(:,:),ALLOCATABLE :: n_kpol  

    ! ... number of coefficients for given structure

    INTEGER(4),DIMENSION(:),ALLOCATABLE :: nc_type  
    
    End MODULE coeffs


!====================================================================
     Subroutine Allocate_coeffs
!====================================================================

     USE coeffs

     if(allocated(ncoef)) Deallocate(ncoef, kcoef, nc_type, n_kpol)

     Allocate(NCOEF(ncase),KCOEF(ncase), nc_type(ntype), &
              n_kpol(ncase,0:mkpol))
    
     END Subroutine Allocate_coeffs
