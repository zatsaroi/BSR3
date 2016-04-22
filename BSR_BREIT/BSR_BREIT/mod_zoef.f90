!======================================================================
      MODULE ZOEF_list
!======================================================================
!
!     Containes the coefficients for two configuration symmetries.
!     Each coefficient has two identifiers: iz_int and ia_df,
!     the pointers for integral and overlap factor, respectively. 
!
!----------------------------------------------------------------------

      Use param_br, ONLY: iszoef

      Implicit none
      Save

      Integer(4) :: nzoef = 0       ! number of coefficients
      Integer(4) :: mzoef = 0       ! current dimentsion of the list
      Integer(4) :: izoef = iszoef  ! initial dimension
    
      Integer(4), Allocatable, Dimension(:) :: IZ_int,IZ_df   
      Real(8), Allocatable, Dimension(:) :: Zoef

      End MODULE ZOEF_list


!======================================================================
      Subroutine alloc_zoef(m)
!======================================================================

      Use inout_br, ONLY: nus, pri
      Use zoef_list

      Implicit none
      Integer(4), Intent(in) :: m

      if(m.le.0) then
       if(allocated(IZ_int)) Deallocate(IZ_int,IZ_DF,Zoef);  mzoef=0
      elseif(.not.allocated(IZ_int)) then
       Allocate(IZ_int(m), IZ_DF(m), Zoef(m));  mzoef=m
      elseif(m.le.mzoef) then
       Return
      elseif(nzoef.eq.0) then
       Deallocate(IZ_int,IZ_DF,Zoef)
       Allocate(IZ_int(m), IZ_DF(m), Zoef(m));  mzoef=m
      else
       Open(nus,status='SCRATCH',form='UNFORMATTED'); rewind(nus)
       write(nus) IZ_int(1:nzoef)
       write(nus) IZ_DF(1:nzoef)
       write(nus) Zoef(1:nzoef)
       Deallocate(IZ_int,IZ_DF,Zoef)
       Allocate(IZ_int(m), IZ_DF(m), Zoef(m));  mzoef = m
       rewind(nus)
       read(nus) IZ_int(1:nzoef)
       read(nus) IZ_DF(1:nzoef)
       read(nus) Zoef(1:nzoef)
       Close(nus)
       write(*,*) ' realloc_zoef: m = ', m
       write(pri,*) ' realloc_zoef: m = ', m
      end if

      End Subroutine alloc_zoef



!=======================================================================
      Subroutine Iadd_zoef(C,int,idf)
!=======================================================================
!
!     add new integral to the list 'zoef'
!
!-----------------------------------------------------------------------

      Use zoef_list

      Implicit none
      Integer(4), Intent(in) :: int,idf
      Real(8), Intent(in) :: C
      Integer(4) :: i

      if(mzoef.eq.0) Call Alloc_zoef(izoef)

! ... check if the same integral is already in list:

      Do i=1,nzoef
       if(int.ne.IZ_int(i)) Cycle
       if(idf.ne.IZ_DF(i)) Cycle
       Zoef(i)=Zoef(i)+C; Return
      End do

! ... add new integral:

      if(nzoef.eq.mzoef) Call Alloc_zoef(mzoef+izoef)
      nzoef = nzoef+1 
      Zoef(nzoef)=C; IZ_int(nzoef)=int; IZ_DF(nzoef)=idf

      End Subroutine Iadd_zoef

