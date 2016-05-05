!====================================================================
      Module coeffs
!====================================================================
!     containes the description and type of possible coefficients
!--------------------------------------------------------------------
      Implicit none 
    
      Integer :: ncase = 11     !  number of different cases
      Integer :: icase =  0     !  case under consideration    
    
      ! ... number of coefficients for given case

      Integer, allocatable :: NCOEF(:),KCOEF(:)   

      Integer :: ntype =  9     !  number of different structures
      Integer :: itype =  0     !  structure under consideration    
      Integer :: kpol  =  0     !  current multipole index                   
    
      ! ... number of coefficients for given multipole and case

      Integer :: mkpol =100     !  max. multipole                   
      Integer, allocatable :: n_kpol(:,:)  

      ! ... number of coefficients for given structure

      Integer, allocatable :: nc_type(:)  
    
      End Module coeffs


!====================================================================
      Subroutine Allocate_coeffs
!====================================================================
!     allocate arrays in module "coeffs"
!--------------------------------------------------------------------
      Use coeffs

      if(allocated(ncoef)) Deallocate(ncoef, kcoef, nc_type, n_kpol)

      Allocate(NCOEF(ncase),KCOEF(ncase), nc_type(ntype), &
               n_kpol(ncase,0:mkpol))
    
      End Subroutine Allocate_coeffs
