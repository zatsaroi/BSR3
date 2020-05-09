!=================================================================
      MODULE new_dets
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
      Integer :: indet = 10000  !  initial suggestion for mndet
      Integer :: jndet = 10000  !  increment for mndet

!     values of the bound-orbitals overlap determinants:

      Real(8), allocatable :: ADET(:)

!     pointer on the 1- or 2-electron overlaps with continuum
!     orbitals:

      Integer, allocatable :: IZOD(:),JZOD(:)

      END MODULE new_dets

!================================================================
      Subroutine Allocate_ndets(m)
!================================================================
!     allocate, deallocate or relocate arrays in module new_defs
!     input:  m - new dimensions
!----------------------------------------------------------------

      USE new_dets

      Implicit none
      Integer, intent(in) :: m
      Integer, allocatable :: iarr(:)
      Real(8), allocatable :: rarr(:)

      if(m.le.0) then
       if(Allocated(ADET)) Deallocate(ADET,IZOD,JZOD) 
	mndet = 0; nndet = 0
       if(m.eq.0) Return 
       mndet=indet
       Allocate(ADET(mndet),IZOD(mndet),JZOD(mndet))
      elseif(m.gt.mndet) then
       if(nndet.le.0) then
        if(Allocated(ADET)) Deallocate(ADET,IZOD,JZOD) 
        mndet = m
        Allocate(ADET(mndet),IZOD(mndet),JZOD(mndet))
       else 
        Allocate(iarr(nndet))
        iarr(1:nndet)=IZOD(1:nndet); Deallocate(IZOD); Allocate(IZOD(m))
        IZOD = 0; IZOD(1:nndet)=iarr(1:nndet)
        iarr(1:nndet)=JZOD(1:nndet); Deallocate(JZOD); Allocate(JZOD(m))
        JZOD =0; JZOD(1:nndet)=iarr(1:nndet)
        Deallocate(iarr)
        Allocate(rarr(nndet))
        rarr(1:nndet)=ADET(1:nndet); Deallocate(ADET); Allocate(ADET(m))
        ADET=0.d0; ADET(1:nndet)=rarr(1:nndet)
        Deallocate(rarr)
        mndet = m
       end if
      end if

      END Subroutine Allocate_ndets


!================================================================
      Subroutine Iadd_ndets(io,jo,S)
!================================================================

      USE bsr_mat, ONLY: Eps_det
      USE new_dets

      Implicit none

      Integer, intent(in) :: io,jo
      Real(8), intent(in) :: S

      if(abs(S).lt.Eps_det) Return

      if(nndet+1.gt.mndet) Call Allocate_ndets(mndet+jndet)
      nndet = nndet + 1
      IZOD(nndet) = io
      JZOD(nndet) = jo
      ADET(nndet) = S

      END Subroutine Iadd_ndets


