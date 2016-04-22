!====================================================================
      MODULE rk_data
!====================================================================
!
!     contains a set of RK coefficients 
!   
!--------------------------------------------------------------------

      IMPLICIT NONE 
      SAVE
    
      INTEGER(4) :: nrk = 0       ! current number of coefficients
      INTEGER(4) :: mrk = 0       ! maximum dimension
      INTEGER(4) :: irk = 10000   ! initial dimension

! ... coefficients:

      REAL(8),ALLOCATABLE :: crk(:)   
    
! ... their attributes:

      INTEGER(4),ALLOCATABLE :: kr1(:),kr2(:),kr3(:),kr4(:)

! ... structure of data:

      Integer(4) ::  ibr = 2**15

! ... kr1 = ic*ibr+jc 
! ... kr2 = k
! ... kr3 = i1*ibr+i2
! ... kr4 = i3*ibr+i4

 
      End MODULE rk_data


!======================================================================
      Subroutine alloc_rk_data(m)
!======================================================================

      Use constants, ONLY: nus

      Use rk_data

      Implicit none
      Integer(4) :: m,i

      if(m.le.0) then
       if(allocated(crk)) Deallocate (crk,kr1,kr2,kr3,kr4)
       nrk = 0; mrk = 0 
      elseif(.not.allocated(crk)) then
       mrk = m; nrk = 0
       Allocate(crk(mrk),kr1(mrk),kr2(mrk),kr3(mrk),kr4(mrk))
      elseif(m.le.mrk) then
       Return
      elseif(nrk.eq.0) then
       Deallocate (crk,kr1,kr2,kr3,kr4)
       mrk = m
       Allocate(crk(mrk),kr1(mrk),kr2(mrk),kr3(mrk),kr4(mrk))
      else
       Open(nus,status='SCRATCH',form='UNFORMATTED'); rewind(nus)
       Do i = 1,nrk
        write(nus) crk(i),kr1(i),kr2(i),kr3(i),kr4(i)
       End do
       Deallocate (crk,kr1,kr2,kr3,kr4)
       mrk = m
       Allocate(crk(mrk),kr1(mrk),kr2(mrk),kr3(mrk),kr4(mrk))
       rewind(nus)
       Do i = 1,nrk
        read(nus) crk(i),kr1(i),kr2(i),kr3(i),kr4(i)
       End do
       Close(nus)
      end if

      End Subroutine alloc_rk_data



!======================================================================
      Subroutine Add_rk_data(C,ic,jc,kr,i1,i2,i3,i4)
!======================================================================
!
!     add new data to the list 
!
!----------------------------------------------------------------------

      USE rk_data

      Implicit none

      Integer(4), intent(in) :: ic,jc,kr,i1,i2,i3,i4
      Real(8),intent(in) :: C
      Integer(4) :: k1,k2,k3,k4, i,k,l,m

      if(mrk.eq.0) Call alloc_rk_data(irk)

      if(ic.ge.jc) then
!       k1=ic*ibr+jc;
       k1=ic*(ic-1)/2+jc; k2=kr; k3=i1*ibr+i2; k4=i3*ibr+i4 
      else
!       k1=jc*ibr+ic; 
       k1=jc*(jc-1)/2+ic; k2=kr; k3=i3*ibr+i4; k4=i1*ibr+i2 
      end if

! ... search position (k) for new integral

      k=1; l=nrk
    1 if(k.gt.l) go to 2              
      m=(k+l)/2
      if    (k1.lt.kr1(m)) then;       l = m - 1
      elseif(k1.gt.kr1(m)) then;       k = m + 1
      else
       if    (k2.lt.kr2(m)) then;      l = m - 1
       elseif(k2.gt.kr2(m)) then;      k = m + 1
       else
        if    (k3.lt.kr3(m)) then;     l = m - 1
        elseif(k3.gt.kr3(m)) then;     k = m + 1
        else
         if    (k4.lt.kr4(m)) then;    l = m - 1
         elseif(k4.gt.kr4(m)) then;    k = m + 1
         else
          crk(m)=crk(m)+C; Return ! the same integral
         end if
        end if
       end if
      end if
      go to 1
    2 Continue 

! ... shift the rest data up:

      Do i=nrk,k,-1
       m = i + 1
       crk(m)=crk(i)
       kr1(m)=kr1(i); kr2(m)=kr2(i); kr3(m)=kr3(i); kr4(m)=kr4(i)
      End do

! ... add new integral:

      crk(k)=C; kr1(k)=k1; kr2(k)=k2; kr3(k)=k3; kr4(k)=k4; nrk=nrk+1

      if(nrk.eq.mrk) Call alloc_rk_data(mrk+irk) 

      END Subroutine Add_rk_data


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

!====================================================================
      MODULE L_data
!====================================================================
!
!     contains a set of L-integral coefficients 
!   
!--------------------------------------------------------------------

      IMPLICIT NONE 
      SAVE
    
      INTEGER(4) :: nlk = 0       ! current number of coefficients
      INTEGER(4) :: mlk = 0       ! maximum dimension
      INTEGER(4) :: ilk = 1000    ! initial dimension

! ... coefficients

      REAL(8),ALLOCATABLE :: clk(:)   
    
! ... their attributes:

      INTEGER(4),ALLOCATABLE :: kl1(:),kl2(:)

! ... structure of data:

      Integer(4) ::  ibl = 2**15

! ... kl1 = (ic-1)*ic/2+jc 
! ... kl2 = i1*ibl+i2

 
      End MODULE L_data


!======================================================================
      Subroutine alloc_L_data(m)
!======================================================================

      Use constants, ONLY: nus

      Use L_data

      Implicit none
      Integer(4) :: m,i

      if(m.le.0) then
       if(allocated(clk)) Deallocate (clk,kl1,kl2)
       nlk = 0; mlk = 0 
      elseif(.not.allocated(clk)) then
       mlk = m; nlk = 0
       Allocate(clk(mlk),kl1(mlk),kl2(mlk))
      elseif(m.le.mlk) then
       Return
      elseif(nlk.eq.0) then
       Deallocate (clk,kl1,kl2)
       mlk = m
       Allocate(clk(mlk),kl1(mlk),kl2(mlk))
      else
       Open(nus,status='SCRATCH',form='UNFORMATTED'); rewind(nus)
       Do i = 1,nlk
        write(nus) clk(i),kl1(i),kl2(i)
       End do
       Deallocate (clk,kl1,kl2)
       mlk = m
       Allocate(clk(mlk),kl1(mlk),kl2(mlk))
       rewind(nus)
       Do i = 1,nlk
        read(nus) clk(i),kl1(i),kl2(i)
       End do
       Close(nus)
      end if

      End Subroutine alloc_L_data



!======================================================================
      Subroutine Add_L_data(C,ic,jc,i1,i2)
!======================================================================
!
!     add new data to the list 
!
!----------------------------------------------------------------------

      USE L_data

      Implicit none

      Integer(4), intent(in) :: ic,jc,i1,i2
      Real(8),intent(in) :: C
      Integer(4) :: k1,k2, i,k,l,m

      if(mlk.eq.0) Call alloc_L_data(ilk)

      if(ic.ge.jc) then
!       k1=ic*ibl+jc; 
       k1=ic*(ic-1)/2+jc; k2=i1*ibl+i2 
      else
!       k1=jc*ibl+ic; 
       k1=jc*(jc-1)/2+ic;  k2=i2*ibl+i1 
      end if

! ... search position (k) for new integral

      k=1; l=nlk
    1 if(k.gt.l) go to 2              
      m=(k+l)/2
      if    (k1.lt.kl1(m)) then;       l = m - 1
      elseif(k1.gt.kl1(m)) then;       k = m + 1
      else
       if    (k2.lt.kl2(m)) then;      l = m - 1
       elseif(k2.gt.kl2(m)) then;      k = m + 1
       else
        clk(m)=clk(m)+C; Return ! the same integral
       end if
      end if
      go to 1
    2 Continue 

! ... shift the rest data up:

      Do i=nlk,k,-1
       m = i + 1
       clk(m)=clk(i); kl1(m)=kl1(i); kl2(m)=kl2(i)
      End do

! ... add new integral:

      clk(k)=C; kl1(k)=k1; kl2(k)=k2; nlk=nlk+1

      if(nlk.eq.mlk) Call alloc_L_data(mlk+ilk) 

      END Subroutine Add_L_data
