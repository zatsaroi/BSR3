!======================================================================
      MODULE new_defs
!======================================================================
!     contains the desription of overlap deferminats which are
!     obtaned by expansion of total overlap deferminants on rows
!     and columns containing the continuum orbitals.
!     The deferminants with known bound orbitals are estimated
!     in line, so, only there values are recoded 
!----------------------------------------------------------------------
      Implicit none
   
      Integer :: mndef = 0      !  max.number of new overlaps
      Integer :: nndef = 0      !  curent number of new overlaps
      Integer :: indef = 1000   !  initial suggestion for mndef
      Integer :: jndef = 1000   !  increment for mndef

!     values of the bound-orbitals overlap deferminants:

      Real(8), allocatable :: Adef(:)

!     pointer on the rest 1- or 2-electron overlaps with continuum
!     orbitals:

      Integer, allocatable :: iof(:),jof(:)

      End Module new_defs


!======================================================================
      Subroutine Allocate_ndefs(m)
!======================================================================
!     allocate arrays in module "ndefs"
!----------------------------------------------------------------------
      Use new_defs

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
       if(nndef.eq.0) then
        if(allocated(Adef)) Deallocate(Adef,iof,jof) 
        mndef = m; nndef = 0
        Allocate(Adef(mndef),iof(mndef),jof(mndef))
	      else 
        Allocate(iarr(nndef))
        iarr(1:nndef)=iof(1:nndef); Deallocate(iof); Allocate(iof(m))
        iof(1:nndef)=iarr(1:nndef)
        iarr(1:nndef)=jof(1:nndef); Deallocate(jof); Allocate(jof(m))
        jof(1:nndef)=iarr(1:nndef)
        Deallocate(iarr)
        Allocate(rarr(nndef))
        rarr(1:nndef)=Adef(1:nndef); Deallocate(Adef); Allocate(Adef(m))
        Adef(1:nndef)=rarr(1:nndef)
        Deallocate(rarr)
	       mndef = m
       end if
      end if

      End Subroutine Allocate_ndefs


!======================================================================
      Subroutine Iadd_ndefs(io,jo,S)
!======================================================================
!     add new entry in module "ndefs"
!----------------------------------------------------------------------
      Use bsr_dmat, ONLY: eps_ndet
      Use new_defs

      Implicit none
      Integer, intent(in) :: io,jo
      Real(8), intent(in) :: S
      Integer :: i,m

      if(abs(S).lt.Eps_ndet) Return

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

      End Subroutine Iadd_ndefs


