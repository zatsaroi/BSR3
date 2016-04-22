!====================================================================
    MODULE cmdata
!====================================================================
!
!  contains a set of coefficients with three identifiers (k1,k2,k3)
!  and one pointer (ipt)
!   
!--------------------------------------------------------------------

    IMPLICIT NONE 
    SAVE
    
    INTEGER(4) :: mcdata =  0     !  max. number of coefficients
    INTEGER(4) :: ncdata =  0     !  current number of coefficients

    REAL(8),DIMENSION(:),ALLOCATABLE :: CDATA      ! coefficients
    
    INTEGER(4),DIMENSION(:),ALLOCATABLE :: K1,K2,K3,IPT

    End MODULE cmdata


!====================================================================
     Subroutine Allocate_cmdata(m)
!====================================================================

     USE cmdata

     INTEGER(4), INTENT(in) :: m
     
     if(m.eq.0) then
       if(Allocated(CDATA)) Deallocate(CDATA,K1,K2,K3,IPT)
       mcdata = 0
       ncdata = 0
     elseif(.not.Allocated(CDATA)) then
       Allocate(CDATA(m),K1(m),K2(m),K3(m),IPT(m))
       mcdata = m       
       ncdata = 0
     elseif(m.gt.mcdata) then
       if(Allocated(CDATA)) Deallocate(CDATA,K1,K2,K3,IPT)
       Allocate(CDATA(m),K1(m),K2(m),K3(m),IPT(m))
       mcdata = m       
     end if
	        
     END Subroutine Allocate_cmdata
