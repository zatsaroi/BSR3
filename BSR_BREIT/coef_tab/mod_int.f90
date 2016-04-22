!====================================================================
      MODULE int_list
!====================================================================
!
!     contains a set of integrals 
!
!--------------------------------------------------------------------

      IMPLICIT NONE 
      SAVE
    
      INTEGER(4) :: mint =  0      !  max. number of coefficients
      INTEGER(4) :: nint =  0      !  current number of coefficients
      INTEGER(4) :: kint = 2**10   !  initial space  

      INTEGER(4),DIMENSION(:),ALLOCATABLE :: Jcase,Jpol,J1,J2,J3,J4

      End MODULE int_list



!======================================================================
      Subroutine alloc_INT (m)
!======================================================================

      Use int_list; Use param_br, Only: nus

      Implicit none
      Integer(4), Intent(in) :: m

      if(m.le.0) then
       if(allocated(Jcase)) Deallocate (Jcase,Jpol,J1,J2,J3,J4)
       mint = 0; nint = 0
      elseif(.not.allocated(Jcase)) then
       mint = m; Allocate(Jcase(m),Jpol(m),J1(m),J2(m),J3(m),J4(m))
      elseif(m.le.mint) then
       Return
      elseif(nint.eq.0) then
       Deallocate (Jcase,Jpol,J1,J2,J3,J4);  mint = m
       Allocate(Jcase(m),Jpol(m),J1(m),J2(m),J3(m),J4(m))
      else
       Open(nus,status='SCRATCH',form='UNFORMATTED'); rewind(nus)
       write(nus) Jcase(1:nint); write(nus) Jpol(1:nint)
       write(nus) J1(1:nint); write(nus) J2(1:nint);
       write(nus) J3(1:nint); write(nus) J4(1:nint);
       Deallocate (Jcase,Jpol,J1,J2,J3,J4)
       mint = m; Allocate(Jcase(m),Jpol(m),J1(m),J2(m),J3(m),J4(m))
       rewind(nus)
       read(nus) Jcase(1:nint); read(nus) Jpol(1:nint)
       read(nus) J1(1:nint); read(nus) J2(1:nint);
       read(nus) J3(1:nint); read(nus) J4(1:nint);
       Close(nus)
       write(*,*) ' Realloc_int: new dimension = ', mint
      end if

      End Subroutine alloc_INT


!=======================================================================
      Integer(4) Function Iadd_int(icase,kpol,i1,i2,i3,i4)
!=======================================================================
!
!     add the integral to the list 'int_list'
!
!-----------------------------------------------------------------------

      Use int_list

      Implicit none
      Integer(4), Intent(in) :: icase,kpol,i1,i2,i3,i4
      Integer(4) :: i

      if(mint.eq.0) Call Alloc_int(kint)

! ... check if the same integral is already in list:

      Do i=1,nint
       if(icase.ne.Jcase(i)) Cycle
       if(kpol.ne.Jpol(i)) Cycle
       if(i1.ne.J1(i)) Cycle
       if(i2.ne.J2(i)) Cycle
       if(i3.ne.J3(i)) Cycle
       if(i4.ne.J4(i)) Cycle
       Iadd_int=i; Return
      End do

! ... add new integral:

      if(nint.eq.mint) Call Alloc_int(mint+kint)
      nint=nint+1; Jcase(nint)=icase; Jpol(nint)=kpol 
      J1(nint)=i1; J2(nint)=i2; J3(nint)=i3; J4(nint)=i4
      Iadd_int = nint

      End Function Iadd_int


