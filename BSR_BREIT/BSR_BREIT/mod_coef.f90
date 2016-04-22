!====================================================================
      MODULE coef_list
!====================================================================
!
!     Contains a set of coefficients with two identifiers (intc,idfc)
!     The list of ordered according the pointer 'ipcoef'
!   
!--------------------------------------------------------------------

      USE param_br, ONLY: iscoef

      IMPLICIT NONE 
      SAVE
    
! ... number of coefficients:

      INTEGER(4) :: ncoef = 0  
      INTEGER(4) :: mcoef = 0  
      INTEGER(4) :: icoef = iscoef
      INTEGER(4) :: kcoef = 0
      INTEGER(4) :: ntotc = 0
      INTEGER(4) :: ntrm  = 0         !  number of terms
      INTEGER(4) :: mtrm  = 0        

! ... coefficients (1:ntrm,1:mcoef):

      REAL(8),DIMENSION(:,:),ALLOCATABLE :: coef     
      REAL(8),DIMENSION(:),ALLOCATABLE :: ctrm  
    
! ... their attributes:

      INTEGER(4),DIMENSION(:),ALLOCATABLE :: intc,idfc,ijhm

! ... ordering pointer: 

      INTEGER(4), Allocatable, Dimension(:) :: ipcoef   

! ... current integral under consideration:

      INTEGER(4) :: int, idf, jcase

      End MODULE coef_list


!======================================================================
     Subroutine Alloc_coef(nt,mc)
!======================================================================

     USE coef_list
     USE inout_br, ONLY: nus, pri

     Implicit none
     Integer(4), Intent(in) :: nt,mc
     Integer(4) :: i

     if(nt.le.0.and.mc.gt.0) Stop ' Alloc_coef: nt = 0 '

     if(mc.le.0) then
      if(Allocated(coef)) &
       Deallocate(coef,intc,idfc,ipcoef,ctrm,ijhm)
       ncoef=0; mcoef=0; mtrm=0
     elseif(.not.Allocated(coef)) then
      Allocate(coef(nt,mc),intc(mc),idfc(mc),ipcoef(mc),&
               ctrm(nt),ijhm(nt)); mcoef=mc; mtrm=nt
     elseif(ncoef.le.0) then
      Deallocate(coef,intc,idfc,ipcoef,ctrm,ijhm)
      Allocate(coef(nt,mc),intc(mc),idfc(mc),ipcoef(mc),&
               ctrm(nt),ijhm(nt)); mcoef=mc; mtrm=nt
     elseif(ncoef.gt.0.and.(mc.gt.mcoef.or.nt.gt.mtrm)) then
      Open(nus,status='scratch',form='UNFORMATTED') 
      rewind(nus) 
      Do i = 1,mcoef
       write(nus) coef(1:mtrm,i)
       write(nus) intc(i)
       write(nus) idfc(i)
       write(nus) ipcoef(i)
      End do
      write(nus) ctrm(1:mtrm)
      write(nus) ijhm(1:mtrm)
      Deallocate(coef,intc,idfc,ipcoef,ctrm,ijhm)
      Allocate(coef(nt,mc),intc(mc),idfc(mc),ipcoef(mc),&
               ctrm(nt),ijhm(nt))
      rewind(nus) 
      Do i = 1,mcoef
       read(nus) coef(1:mtrm,i)
       read(nus) intc(i)
       read(nus) idfc(i)
       read(nus) ipcoef(i)
      End do
      read(nus) ctrm(1:mtrm)
      read(nus) ijhm(1:mtrm)
      mcoef=mc; mtrm=nt
      write(*,*) ' realloc_coef: mcoef = ', mcoef,mtrm
      write(pri,*) ' realloc_coef: mcoef = ', mcoef,mtrm
      Close(nus)
     end if
	        
     END Subroutine Alloc_coef


!======================================================================
      Subroutine Add_coef
!======================================================================
!
!     add new coefficient to the list 
!
!----------------------------------------------------------------------

      USE coef_list

      Implicit none

      Integer(4) :: k,l,m,ipm

! ... look for the same integral in the list

      k=1; l = ncoef 

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
       Do m = ncoef,k+1,-1; ipcoef(m) = ipcoef(m-1); End do
       ipcoef(k)=ncoef
      end if        

! ... it is time for relocation:

      if(ncoef.eq.mcoef) Call Alloc_coef (ntrm,mcoef+icoef)
 
      END Subroutine Add_coef


