!=================================================================
      Module new_dets
!=================================================================
!     contains the desription of overlap determinats which are
!     obtaned by expansion of total overlap determinants on rows
!     and columns containing the continuum orbitals.
!     The determinants with known bound orbitals are estimated
!     in line, so only their values are recoded 
!-----------------------------------------------------------------
      Implicit none
  
      Integer :: mndet = 0      !  max.number of new overlaps
      Integer :: nndet = 0      !  curent number of new overlaps
      Integer :: indet = 1000   !  initial suggestion for mndet
      Integer :: jndet = 1000   !  increment for mndet

!     values of the bound-orbitals overlap determinants:

      Real(8), allocatable :: Adet(:)

!     pointer on the 1- or 2-electron overlaps with continuum
!     orbitals:

      Integer, allocatable :: IZOD(:),JZOD(:)

      End Module new_dets

!======================================================================
      Subroutine Allocate_ndets(m)
!======================================================================
!     allocate arrays in module "ndets"
!----------------------------------------------------------------------
      Use new_dets

      Integer, intent(in) :: m
      Integer, allocatable :: iarr(:)
      Real(8), allocatable :: rarr(:)

      if(m.le.0) then
       if(allocated(ADET)) Deallocate(ADET,IZOD,JZOD) 
	      mndet = 0; nndet = 0
       if(m.eq.0) Return 
       mndet=indet
       Allocate(ADET(mndet),IZOD(mndet),JZOD(mndet))
      elseif(m.gt.mndet) then
       if(nndet.eq.0) then
        if(allocated(ADET)) Deallocate(ADET,IZOD,JZOD) 
        mndet = m
        Allocate(ADET(mndet),IZOD(mndet),JZOD(mndet))
	      else 
        Allocate(iarr(nndet))
        iarr(1:nndet)=IZOD(1:nndet); Deallocate(IZOD); Allocate(IZOD(m))
        IZOD(1:nndet)=iarr(1:nndet)
        iarr(1:nndet)=JZOD(1:nndet); Deallocate(JZOD); Allocate(JZOD(m))
        JZOD(1:nndet)=iarr(1:nndet)
        Deallocate(iarr)
        Allocate(rarr(nndet))
        rarr(1:nndet)=ADET(1:nndet); Deallocate(ADET); Allocate(ADET(m))
        ADET(1:nndet)=rarr(1:nndet)
        Deallocate(rarr)
	       mndef = m
       end if
      end if

      End Subroutine Allocate_ndets


!======================================================================
      Subroutine Iadd_ndets(io,jo,S)
!======================================================================
!     add new entry in the module "ndets"
!----------------------------------------------------------------------
      Use bsr_dmat, only: Eps_ndet
      Use new_dets

      Implicit none
      Integer, intent(in) :: io,jo
      Real(8), intent(in) :: S

      if(abs(S).lt.Eps_ndet) Return

      if(nndet+1.gt.mndet) Call Allocate_ndets(mndet+jndet)
      nndet = nndet + 1
      IZOD(nndet) = io
      JZOD(nndet) = jo
      ADET(nndet) = S

      End Subroutine Iadd_ndets


