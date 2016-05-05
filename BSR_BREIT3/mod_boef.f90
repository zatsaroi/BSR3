!======================================================================
      Module boef_list
!======================================================================
!     Containes the two-electron integrals for matrix elements
!     in uncouple nlms-representation.
!     This list is introduced to decrease the number of calls for
!     ZNO_breit subroutine. Used only in the ZNO_2ee routine.
!     The coefficients are recorded by blocks - all integrals for all
!     operators under concideration for given < 1, 2 | O | 3, 4>
!----------------------------------------------------------------------
      Implicit none

      Integer :: nboef = 0       ! number of integrals
      Integer :: mboef = 0       ! current dimension of list
      Integer :: iboef = 50000   ! initial dimension
      Integer :: jboef = 0       ! first new element

! ... ib_int(1:mboef) - integral indentifier
! ... boef(1:mboef) - correspondent angular coefficient

      Integer, allocatable :: ib_int(:)     
      Real(8), allocatable :: boef(:)

! ... packing parameters:
 
      Integer, parameter :: jbl  = 2**7      ! in Check_BOEF
      Integer, parameter :: jbm  = 2**6      ! in Check_BOEF
 
      Integer :: nblk  = 0       ! number of blocks
      Integer :: mblk  = 0       ! current dimension of list
      Integer :: iblk  = 5000    ! initial dimentsion
      Integer :: kblk  = 0       ! block under consideration

! ... identifiers of block:

      Integer, allocatable :: indl(:),indm(:),inds(:) 

! ... current identifiers:

      Integer :: inl,inm,ins     

! ... ncblk - pointer on the last element in the block
! ... ipblk - ordering pointer

      Integer, allocatable :: ipblk(:),ncblk(:) 

      End Module BOEF_list


!======================================================================
      Subroutine alloc_boef(m)
!======================================================================
!     allocate/deallocate arrays in module boef
!---------------------------------------------------------------------- 
      Use boef_list

      Implicit none
      Integer, intent(in) :: m
      Integer, allocatable :: ia(:)
      Real(8), allocatable :: ra(:)

      if(m.le.0) then
       mboef=0; nboef=0
       if(allocated(IB_int)) Deallocate(IB_int,Boef)
       if(m.lt.0) then
        mboef=iboef
        Allocate(IB_int(mboef), Boef(mboef))
       end if
      elseif(.not.allocated(IB_int)) then
       Allocate(IB_int(m), Boef(m));  mboef=m
      elseif(m.le.mboef) then
       Return
      elseif(nboef.eq.0) then
       Deallocate(IB_int,Boef)
       Allocate(IB_int(m), Boef(m));  mboef=m
      else
       mboef = m
   	   Allocate(ia(nboef))
       ia=IB_int(1:nboef); Deallocate(IB_int)
	      Allocate(IB_int(mboef)); IB_int(1:nboef)=ia 
       Deallocate(ia)
	      Allocate(ra(nboef))
       ra=Boef(1:nboef); Deallocate(Boef)
	      Allocate(Boef(mboef)); Boef(1:nboef)=ra 
       Deallocate(ra)
       write(*,*) 'realloc_boef: m = ', m, iboef
      end if

      End Subroutine alloc_boef


!======================================================================
      Subroutine alloc_blk(m)
!======================================================================
!     allocate/deallocate 'block' arrays in module boef
!---------------------------------------------------------------------- 
      Use boef_list

      Implicit none
      Integer, intent(in) :: m
      Integer, allocatable :: ia(:)

      if(m.le.0) then
       mblk = 0; nblk = 0
       if(allocated(ipblk)) Deallocate(ipblk,ncblk,indl,indm,inds)
       if(m.lt.0) then
        mblk = iblk
        Allocate(ipblk(mblk),ncblk(mblk),indl(mblk),indm(mblk),inds(mblk))
       end if
      elseif(.not.allocated(ipblk)) then
       Allocate(ipblk(m),ncblk(m),indl(m),indm(m),inds(m)); mblk = m
      elseif(m.le.mblk) then
       Return
      elseif(nblk.eq.0) then
       Deallocate(ipblk,ncblk,indl,indm,inds)
       Allocate(ipblk(m),ncblk(m),indl(m),indm(m),inds(m)); mblk = m
      else
       mblk = m
   	   Allocate(ia(nblk))
       ia=ipblk(1:nblk); Deallocate(ipblk)
	      Allocate(ipblk(mblk)); ipblk(1:nblk)=ia 
       ia=ncblk(1:nblk); Deallocate(ncblk)
	      Allocate(ncblk(mblk)); ncblk(1:nblk)=ia 
       ia=indl(1:nblk); Deallocate(indl)
	      Allocate(indl(mblk)); indl(1:nblk)=ia 
       ia=indm(1:nblk); Deallocate(indm)
	      Allocate(indm(mblk)); indm(1:nblk)=ia 
       ia=inds(1:nblk); Deallocate(inds)
	      Allocate(inds(mblk)); inds(1:nblk)=ia 
       Deallocate(ia)
       write(*,*) 'realloc_blk: m,iblk = ', m,iblk
      end if

      End Subroutine alloc_blk


!=======================================================================
      Subroutine Iadd_boef(C,int)
!=======================================================================
!     add new integral to the list 'boef'
!     (in the range of new block: jboef-nboef)
!-----------------------------------------------------------------------
      Use boef_list

      Implicit none
      Integer, Intent(in) :: int
      Real(8), Intent(in) :: C
      Integer :: i

      if(mboef.eq.0) Call alloc_boef(iboef)

! ... check if the same itegral is already in the list:

      Do i=jboef,nboef
       if(int.ne.ib_int(i)) Cycle;  Boef(i)=Boef(i)+C; Return
      End do

! ... add new integral:

      if(nboef.eq.mboef) Call Alloc_boef(mboef+iboef)

      nboef=nboef+1; Boef(nboef)=C; IB_int(nboef)=int

      End Subroutine Iadd_boef


!=======================================================================
      Subroutine Check_boef(l1,m1,s1,l2,m2,s2,l3,m3,s3,l4,m4,s4)
!=======================================================================
!     Check if the m.e. for given orbitals is already in the list,
!     otherwise - calculate them. 
!     Procedure use incoding the orbitals parameters, and that
!     restrict the max. l to 64 (see parameter jbm, jbl). 
!     We can reduce this restriction  by using four identifiers 
!     (instead three), one for each orbital.  
!----------------------------------------------------------------------
      Use boef_list
      
      Implicit none
      Integer, intent(in) :: l1,m1,s1,l2,m2,s2,l3,m3,s3,l4,m4,s4
      Integer :: i1,i2,i3,i4, k,l,m,ipm 

! ... prepare indentifiers:

      i1 = jbm + m1; i2 = jbm + m2; i3 = jbm + m3; i4 = jbm + m4

      if(i1.ge.jbl.or.i2.ge.jbl.or.i3.ge.jbl.or.i4.ge.jbl) then
       write(*,'(a,3i5,3i10)') 'lms = ',l1,m1,s1,jbm,jbl,i1
       write(*,'(a,3i5,3i10)') 'lms = ',l2,m2,s2,jbm,jbl,i2
       write(*,'(a,3i5,3i10)') 'lms = ',l3,m3,s3,jbm,jbl,i3
       write(*,'(a,3i5,3i10)') 'lms = ',l4,m4,s4,jbm,jbl,i4
       Stop 'Check_boef: ml out of limits '
      end if

      inl = ((l1*jbl+l2)*jbl+l3)*jbl+l4    
      inm = ((i1*jbl+i2)*jbl+i3)*jbl+i4     
      ins = ((s1*jbl+s2)*jbl+s3)*jbl+s4     

! ... look for the same case in the list:

      k=1; l = nblk 
    1 if(k.gt.l) go to 2              
 
      m=(k+l)/2; ipm=ipblk(m)
 
      if(inl.lt.indl(ipm)) then;      l = m - 1
      elseif(inl.gt.indl(ipm)) then;  k = m + 1
      else

      if(inm.lt.indm(ipm)) then;      l = m - 1
      elseif(inm.gt.indm(ipm)) then;  k = m + 1
      else

      if(ins.lt.inds(ipm)) then;      l = m - 1
      elseif(ins.gt.inds(ipm)) then;  k = m + 1

      else;   kblk = ipm;   Return
      end if; end if; end if

      go to 1
    2 Continue 
    
! ... new block:    
            
      jboef = nboef + 1
      Call ZNO_breit(1,l1,m1,s1,2,l2,m2,s2,3,l3,m3,s3,4,l4,m4,s4,+1)
      Call ZNO_breit(1,l1,m1,s1,2,l2,m2,s2,4,l4,m4,s4,3,l3,m3,s3,-1)

      nblk = nblk + 1;  ncblk(nblk) = nboef; kblk = nblk   
      indl(nblk)=inl; indm(nblk)=inm; inds(nblk)=ins
      
      if(k.eq.nblk) then
       ipblk(k)=nblk
      else
       Do m = nblk,k+1,-1; ipblk(m) = ipblk(m-1); End do
       ipblk(k)=nblk
      end if        

! ... it is time for re-allocation:

      if(nblk.eq.mblk) Call Alloc_blk(mblk+iblk)

      End Subroutine Check_boef










