!======================================================================
      Module zoef_list
!======================================================================
!     Containes the coefficients for two configuration symmetries.
!     Each coefficient has two identifiers: iz_int and iz_df,
!     the pointers for integral and overlap factor, respectively. 
!----------------------------------------------------------------------
      Implicit none

      Integer :: nzoef = 0       ! number of coefficients
      Integer :: mzoef = 0       ! current dimentsion of the list
      Integer :: izoef = 2**10   ! initial dimension

      Integer, allocatable :: IZ_int(:),IZ_df(:)   
      Real(8), allocatable :: Zoef(:)

      End Module zoef_list


!======================================================================
      Subroutine alloc_zoef(m)
!======================================================================
!     allocate arrays in module zoef
!----------------------------------------------------------------------
      Use zoef_list

      Implicit none
      Integer, intent(in) :: m
      Integer, allocatable :: ia(:)
      Real(8), allocatable :: ra(:)

      if(m.le.0) then
       mzoef=0; nzoef=0
       if(allocated(IZ_int)) Deallocate(IZ_int,IZ_DF,Zoef)
       if(m.lt.0) then
        mzoef=izoef
        Allocate(IZ_int(mzoef), IZ_DF(mzoef), Zoef(mzoef))
       end if
      elseif(.not.allocated(IZ_int)) then
       Allocate(IZ_int(m), IZ_DF(m), Zoef(m));  mzoef=m
      elseif(m.le.mzoef) then
       Return
      elseif(nzoef.eq.0) then
       Deallocate(IZ_int,IZ_DF,Zoef)
       Allocate(IZ_int(m), IZ_DF(m), Zoef(m));  mzoef=m
      else
       mzoef = m
	      Allocate(ia(nzoef))
       ia=IZ_int(1:nzoef); Deallocate(IZ_int)
	      Allocate(IZ_int(mzoef)); IZ_int(1:nzoef)=ia 
       ia=IZ_df(1:nzoef); Deallocate(IZ_df)
	      Allocate(IZ_df(mzoef)); IZ_df(1:nzoef)=ia 
       Deallocate(ia)
	      Allocate(ra(nzoef))
       ra=Zoef(1:nzoef); Deallocate(Zoef)
	      Allocate(Zoef(mzoef)); Zoef(1:nzoef)=ra 
       Deallocate(ra)
       write(*,*) ' realloc_zoef: m = ', m
      end if

      End Subroutine alloc_zoef


!=======================================================================
      Subroutine Iadd_zoef(C,int,idf)
!=======================================================================
!     add new integral to the list 'zoef'
!-----------------------------------------------------------------------
      Use zoef_list

      Implicit none
      Integer, intent(in) :: int,idf
      Real(8), intent(in) :: C
      Integer :: i

      if(mzoef.eq.0) Call Alloc_zoef(izoef)

! ... check if the same integral is already in list:

      Do i=1,nzoef
       if(int.ne.IZ_int(i)) Cycle
       if(idf.ne.IZ_DF (i)) Cycle
       Zoef(i)=Zoef(i)+C; Return
      End do

! ... add new integral:

      if(nzoef.eq.mzoef) Call Alloc_zoef(mzoef+izoef)
      nzoef = nzoef+1 
      Zoef(nzoef)=C; IZ_int(nzoef)=int; IZ_DF(nzoef)=idf

      End Subroutine Iadd_zoef

