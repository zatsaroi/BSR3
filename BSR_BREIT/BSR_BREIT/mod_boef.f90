!======================================================================
      MODULE boef_list
!======================================================================
!
!     Containes the two-electron integrals for matrix elements
!     in uncouple nlms-representation.
!     This list is introduced to decrease the number of calls for
!     ZNO_breit subroutine. Used in the ZNO_2ee routine.
!     The coefficient in the list are recorded by blocks -
!     all integrals for all operators under concideration for
!     given < 1, 2 | O | 3, 4>
!----------------------------------------------------------------------

      Use param_br, ONLY: isboef,isblk

      Implicit none
      Save

      Integer(4) :: nboef = 0       ! number of integrals
      Integer(4) :: mboef = 0       ! current dimension of list
      Integer(4) :: iboef = isboef  ! initial dimension
      Integer(4) :: jboef = 0       ! first new element

! ... IB_int(1:mboef) - integral indentifier
! ... boef(1:mboef) - correspondent angular coefficient

      Integer(4), Allocatable, Dimension(:) :: IB_int     
      Real(8), Allocatable, Dimension(:) :: boef

      Integer(4) :: nblk  = 0       ! number of blocks
      Integer(4) :: mblk  = 0       ! current dimension of list
      Integer(4) :: iblk  = isblk   ! initial dimentsion
      Integer(4) :: kblk  = 0       ! block under consideration

! ... identifiers of block:

      Integer(4), Allocatable, Dimension(:) :: indl,indm,inds 

! ... current identifiers:

      Integer(4) :: inl,inm,ins     

! ... ncblk - pointer on the last element in the block
! ... ipblk - ordering pointer

      Integer(4), Allocatable, Dimension(:) :: ipblk,ncblk 

      End MODULE BOEF_list


!======================================================================
      Subroutine alloc_boef(m)
!======================================================================

      Use inout_br, ONLY: nus,pri
      Use boef_list

      Implicit none
      Integer(4), Intent(in) :: m

      if(m.le.0) then
       if(allocated(IB_int)) Deallocate(IB_int,Boef);  mboef=0
      elseif(.not.allocated(IB_int)) then
       Allocate(IB_int(m), Boef(m));  mboef=m
      elseif(m.le.mboef) then
       Return
      elseif(nboef.eq.0) then
       Deallocate(IB_int,Boef)
       Allocate(IB_int(m), Boef(m));  mboef=m
      else
       Open(nus,status='SCRATCH',form='UNFORMATTED'); rewind(nus)
       write(nus) IB_int(1:nboef)
       write(nus) Boef(1:nboef)
       Deallocate(IB_int,Boef)
       Allocate(IB_int(m), Boef(m));  mboef=m
       rewind(nus)
       read(nus) IB_int(1:nboef)
       read(nus) Boef(1:nboef)
       Close(nus)
       write(*,*) 'realloc_boef: m = ', m
       write(pri,*) 'realloc_boef: m = ', m
      end if

      End Subroutine alloc_boef



!======================================================================
      Subroutine alloc_blk(m)
!======================================================================

      Use boef_list;  Use inout_br, ONLY: nus,pri

      Implicit none
      Integer(4), Intent(in) :: m

      if(m.le.0) then
       if(allocated(ipblk)) Deallocate(ipblk,ncblk,indl,indm,inds)
       mblk = 0; nblk = 0
      elseif(.not.allocated(ipblk)) then
       Allocate(ipblk(m),ncblk(m),indl(m),indm(m),inds(m)); mblk = m
      elseif(m.le.mblk) then
       Return
      elseif(nblk.eq.0) then
       Deallocate(ipblk,ncblk,indl,indm,inds)
       Allocate(ipblk(m),ncblk(m),indl(m),indm(m),inds(m)); mblk = m
      else
       Open(nus,status='SCRATCH',form='UNFORMATTED'); rewind(nus)
       write(nus) ipblk(1:nblk)
       write(nus) ncblk(1:nblk)
       write(nus) indl (1:nblk)
       write(nus) indm (1:nblk)
       write(nus) inds (1:nblk)
       Deallocate(ipblk,ncblk,indl,indm,inds)
       Allocate(ipblk(m),ncblk(m),indl(m),indm(m),inds(m)); mblk = m
       rewind(nus)
       read(nus) ipblk(1:nblk)
       read(nus) ncblk(1:nblk)
       read(nus) indl (1:nblk)
       read(nus) indm (1:nblk)
       read(nus) inds (1:nblk)
       Close(nus)
       write(*,*) 'realloc_blk: m,iblk = ', m,iblk
       write(pri,*) 'realloc_blk: m,iblk = ', m,iblk
      end if

      End Subroutine alloc_blk


!=======================================================================
      Subroutine Iadd_boef(C,int)
!=======================================================================
!
!     add new integral to the list 'boef'
!     (in the range of new block: jboef-nboef)
!
!-----------------------------------------------------------------------

      Use boef_list

      Implicit none
      Integer(4), Intent(in) :: int
      Real(8), Intent(in) :: C
      Integer(4) :: i

      if(mboef.eq.0) Call alloc_boef(iboef)

! ... check if the same itegral is already in the list:

      Do i=jboef,nboef
       if(int.ne.IB_int(i)) Cycle;  Boef(i)=Boef(i)+C; Return
      End do

! ... add new integral:

      if(nboef.eq.mboef) Call Alloc_boef(mboef+iboef)

      nboef=nboef+1; Boef(nboef)=C; IB_int(nboef)=int

      End Subroutine Iadd_boef


!=======================================================================
      Subroutine Check_boef(l1,m1,s1,l2,m2,s2,l3,m3,s3,l4,m4,s4)
!=======================================================================
!
!     Check if already there is the m.e. for given orbitals,
!     otherwise - calculate them. 
!     Procedure use incoding the orbitals parameters, and that
!     restrict the max. l to 64 (see parameter jbm, jbl). 
!     We can reduce this restriction  by using four identifiers 
!     (instead three), one for each orbital.  
!
!----------------------------------------------------------------------

      USE boef_list; USE param_br, ONLY: jbl,jbm

      Implicit none

      Integer(4), Intent(in) :: l1,m1,s1,l2,m2,s2,l3,m3,s3,l4,m4,s4

      Integer(4) :: i1,i2,i3,i4, k,l,m,ipm 

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










