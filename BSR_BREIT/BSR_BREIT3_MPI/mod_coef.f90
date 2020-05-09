!====================================================================
      MODULE coef_list
!====================================================================
!     Contains a set of coefficients with two identifiers (intc,idfc)
!     The list of ordered according the pointer 'ipcoef'
!--------------------------------------------------------------------

      Implicit none 
    
! ... number of coefficients:

      Integer :: ncoef = 0  
      Integer :: mcoef = 0  
      Integer :: icoef = 1000
      Integer :: kcoef = 0            ! copy for MPI

! ... number of terms

      Integer :: ntrm = 1
      Integer :: ktrm = 1             ! copy for MPI

! ... coefficients (1:ntrm,1:mcoef):

      REAL(8), Allocatable :: coef(:,:)     
      REAL(8), Allocatable :: ctrm(:)  
    
! ... their attributes:

      Integer, Allocatable :: intc(:),idfc(:),ijhm(:)

! ... ordering pointer: 

      Integer, Allocatable :: ipcoef(:)   

! ... current integral under consideration:

      Integer :: int, idf

      End MODULE coef_list


!======================================================================
      Subroutine Alloc_coef(mc)
!======================================================================

      USE coef_list

      Implicit none
      Integer, Intent(in)  :: mc
      Integer, allocatable :: ia(:)
      Real(8), allocatable :: rb(:,:)

      if(mc.le.0) then
       ncoef=0; mcoef=0
       if(allocated(coef)) &
        Deallocate(coef,intc,idfc,ipcoef,ctrm,ijhm) 
       if(mc.lt.0) then
        mcoef=icoef
        Allocate(coef(ntrm,mcoef),ctrm(ntrm),ijhm(ntrm), &
                 intc(mcoef),idfc(mcoef),ipcoef(mcoef))
       end if
      elseif(.not.Allocated(coef)) then
        mcoef=mc; ncoef=0
        Allocate(coef(ntrm,mcoef),ctrm(ntrm),ijhm(ntrm), &
                 intc(mcoef),idfc(mcoef),ipcoef(mcoef))
      elseif(ncoef.le.0) then
        Deallocate(coef,intc,idfc,ipcoef,ctrm,ijhm)
        mcoef=mc; ncoef=0
        Allocate(coef(ntrm,mcoef),ctrm(ntrm),ijhm(ntrm), &
                 intc(mcoef),idfc(mcoef),ipcoef(mcoef))
      elseif(mc.gt.mcoef) then
       mcoef=mc
   	   Allocate(ia(ncoef))
       ia=intc(1:ncoef); Deallocate(intc)
	      Allocate(intc(mcoef)); intc(1:ncoef)=ia 
       ia=idfc(1:ncoef); Deallocate(idfc)
	      Allocate(idfc(mcoef)); idfc(1:ncoef)=ia 
       ia=ipcoef(1:ncoef); Deallocate(ipcoef)
	      Allocate(ipcoef(mcoef)); ipcoef(1:ncoef)=ia 
       Deallocate(ia)
	      Allocate(rb(ntrm,ncoef))
       rb=coef(1:ntrm,1:ncoef); Deallocate(coef)
	      Allocate(coef(ntrm,mcoef)); coef(1:ntrm,1:ncoef)=rb 
       Deallocate(rb)
       write(*,*) ' realloc_coef: ntrm,mcoef,icoef = ', ntrm,mcoef,icoef
      end if
	        
      END Subroutine Alloc_coef


!======================================================================
      Subroutine Add_coef
!======================================================================
!     add new coefficient to the list 
!----------------------------------------------------------------------

      USE coef_list

      Implicit none

      Integer :: k,l,m,ipm

      if(mcoef.eq.0) Call Alloc_coef(-1)

! ... look for the same integral in the list

      k=1; l=ncoef 

    1 if(k.gt.l) go to 2              
      m=(k+l)/2; ipm=ipcoef(m)
      if(int.lt.intc(ipm)) then;      l = m - 1
      elseif(int.gt.intc(ipm)) then;  k = m + 1
      else
       if(idf.lt.idfc(ipm)) then;     l = m - 1
       elseif(idf.gt.idfc(ipm)) then; k = m + 1
       else
        coef(1:ntrm,ipm) = coef(1:ntrm,ipm) + ctrm(1:ntrm)
        Return
       end if
      end if
      go to 1
    2 Continue 

! ... new coefficient:

      ncoef = ncoef + 1  
      coef(:,ncoef)=ctrm(:); intc(ncoef)=int; idfc(ncoef)=idf 
      
      if(k.eq.ncoef) then
       ipcoef(k)=ncoef
      else
       Do m=ncoef,k+1,-1; ipcoef(m)=ipcoef(m-1); End do
       ipcoef(k)=ncoef
      end if        

! ... it is time for relocation:

      if(ncoef.eq.mcoef) Call Alloc_coef(mcoef+icoef)
 
      END Subroutine Add_coef


