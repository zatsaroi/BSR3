!=================================================================
      MODULE new_defs
!=================================================================
!     contains the desription of overlap deferminats which are
!     obtaned by expansion of total overlap deferminants on rows
!     and columns containing the continuum orbitals.
!     The deferminants with known bound orbitals are estimated
!     in line, so, only there values are recoded 
!-----------------------------------------------------------------
      Implicit none
      Integer :: mndef = 0      !  max.number of new overlaps
      Integer :: nndef = 0      !  curent number of new overlaps
      Integer :: indef = 10000  !  initial suggestion for mndef
      Integer :: jndef = 10000  !  increment for mndef

!     values of the bound-orbitals overlap deferminants:

      Real(8), allocatable :: Adef(:)

!     pointer on the 1- or 2-electron overlaps with continuum
!     orbitals:

      Integer, allocatable :: iof(:),jof(:)

      END MODULE new_defs

!================================================================
      Subroutine Allocate_ndefs(m)
!================================================================
!     allocate, deallocate or relocate arrays in module new_defs
!     input:  m - new dimensions
!----------------------------------------------------------------
      USE new_defs

      Implicit none
      Integer, intent(in) :: m
      Integer, allocatable :: iarr(:)
      Real(8), allocatable :: rarr(:)

      if(m.le.0) then
       if(allocated(Adef)) Deallocate(Adef,iof,jof) 
    	  mndef = 0; nndef = 0
       if(m.eq.0) Return
       mndef = indef
       Allocate(Adef(mndef),iof(mndef),jof(mndef))
      elseif(m.gt.mndef) then
       if(nndef.le.0) then
        if(Allocated(Adef)) Deallocate(Adef,iof,jof) 
        mndef = m; nndef = 0
        Allocate(Adef(mndef),iof(mndef),jof(mndef))
       else 
        Allocate(iarr(nndef))
        iarr(1:nndef)=iof(1:nndef); Deallocate(iof); Allocate(iof(m))
        iof=0; iof(1:nndef)=iarr(1:nndef)
        iarr(1:nndef)=jof(1:nndef); Deallocate(jof); Allocate(jof(m))
        jof=0; jof(1:nndef)=iarr(1:nndef)
        Deallocate(iarr)
        Allocate(rarr(nndef))
        rarr(1:nndef)=Adef(1:nndef); Deallocate(Adef); Allocate(Adef(m))
        Adef=0.d0; Adef(1:nndef)=rarr(1:nndef)
        Deallocate(rarr)
	 mndef = m
       end if
      end if

      END Subroutine Allocate_ndefs


!================================================================
      Subroutine Iadd_ndefs(io,jo,S)
!================================================================
!     add new overlap determinant to the list:
!----------------------------------------------------------------
      USE bsr_mat, ONLY: Eps_det
      USE new_defs

      Implicit none
      Integer, intent(in) :: io,jo
      Real(8), intent(in) :: S
      Integer :: i,m

      if(abs(S).lt.Eps_det) Return

      m = 0
      Do i=1,nndef
       if(io.ne.iof(i)) Cycle
       if(jo.ne.jof(i)) Cycle
       Adef(i) = Adef(i) + S
       m = 1
       Exit
      End do
      if(m.ne.0) Return

      if(nndef+1.gt.mndef) Call Allocate_ndefs(mndef+jndef)
      nndef = nndef + 1
      iof(nndef) = io
      jof(nndef) = jo
      Adef(nndef) = S

      END Subroutine Iadd_ndefs


