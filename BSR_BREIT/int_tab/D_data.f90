!====================================================================
      MODULE D_data
!====================================================================
!
!     contains a set of D-integrals coefficients 
!   
!--------------------------------------------------------------------

      IMPLICIT NONE 
      SAVE
    
      INTEGER(4) :: ndk = 0       ! current number of coefficients
      INTEGER(4) :: mdk = 0       ! maximum dimension
      INTEGER(4) :: idk = 1000    ! initial dimension

! ... coefficients

      REAL(8),ALLOCATABLE :: cdk(:)   
    
! ... their attributes:

      INTEGER(4),ALLOCATABLE :: kd1(:),kd2(:)

! ... structure of data:

      Integer(4) ::  ibd = 2**15

! ... kd1 = (ic-1)*ic/2+jc 
! ... kd2 = i1*ibd+i2

 
      End MODULE D_data


!======================================================================
      Subroutine alloc_D_data(m)
!======================================================================

      Use constants, ONLY: nus

      Use D_data

      Implicit none
      Integer(4) :: m,i

      if(m.le.0) then
       if(allocated(cdk)) Deallocate (cdk,kd1,kd2)
       ndk = 0; mdk = 0 
      elseif(.not.allocated(cdk)) then
       mdk = m; ndk = 0
       Allocate(cdk(mdk),kd1(mdk),kd2(mdk))
      elseif(m.le.mdk) then
       Return
      elseif(ndk.eq.0) then
       Deallocate (cdk,kd1,kd2)
       mdk = m
       Allocate(cdk(mdk),kd1(mdk),kd2(mdk))
      else
       Open(nus,status='SCRATCH',form='UNFORMATTED'); rewind(nus)
       Do i = 1,ndk
        write(nus) cdk(i),kd1(i),kd2(i)
       End do
       Deallocate (cdk,kd1,kd2)
       mdk = m
       Allocate(cdk(mdk),kd1(mdk),kd2(mdk))
       rewind(nus)
       Do i = 1,ndk
        read(nus) cdk(i),kd1(i),kd2(i)
       End do
       Close(nus)
      end if

      End Subroutine alloc_D_data



!======================================================================
      Subroutine Add_D_data(C,ic,jc,i1,i2)
!======================================================================
!
!     add new data to the list (ip:jp) in module cmdata
!
!----------------------------------------------------------------------

      USE D_data

      Implicit none

      Integer(4), intent(in) :: ic,jc,i1,i2
      Real(8),intent(in) :: C
      Integer(4) :: k1,k2, i,k,l,m

      if(mdk.eq.0) Call alloc_D_data(idk)

!      k1=ic*ibd+jc; k2=i1*ibd+i2 
      k1=ic*(ic-1)/2+jc; k2=i1*ibd+i2 

! ... search position (k) for new integral

      k=1; l=ndk
    1 if(k.gt.l) go to 2              
      m=(k+l)/2
      if    (k1.lt.kd1(m)) then;       l = m - 1
      elseif(k1.gt.kd1(m)) then;       k = m + 1
      else
       if    (k2.lt.kd2(m)) then;      l = m - 1
       elseif(k2.gt.kd2(m)) then;      k = m + 1
       else
        cdk(m)=cdk(m)+C; Return ! the same integral
       end if
      end if
      go to 1
    2 Continue 

! ... shift the rest data up:

      Do i=ndk,k,-1
       m = i + 1
       cdk(m)=cdk(i); kd1(m)=kd1(i); kd2(m)=kd2(i)
      End do

! ... add new integral:

      cdk(k)=C; kd1(k)=k1; kd2(k)=k2; ndk=ndk+1

      if(ndk.eq.mdk) Call alloc_D_data(mdk+idk) 

      END Subroutine Add_D_data
