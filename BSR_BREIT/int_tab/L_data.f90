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
