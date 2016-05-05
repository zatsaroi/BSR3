!====================================================================
      Module cmdata
!====================================================================
!  Contains a set of coefficients with three identifiers (k1,k2,k3)
!  and order pointer (ipt). The data are devided on different types
!  (ntype, itype) and saved in memory in blocks (nblocks, mblocks).
!--------------------------------------------------------------------
      Implicit none 
    
! ... coefficients:

      Real(8), allocatable :: CLDATA(:), CVDATA(:)
    
! ... their attributes:

      Integer, allocatable :: K1(:),K2(:),K3(:),IPT(:)

! ... number of different structures and current structure 

      Integer, parameter :: ntype =  10   
      Integer :: itype =  0                

! ... the data are divided on blocks, with number and size:

      Integer, parameter :: nblocks =  50     
      Integer, parameter :: mblocks =  20000  

! ... current block for given type:

      Integer, allocatable :: iblk(:)    

! ... pointer on first(last) element in given block 

      Integer, allocatable :: ipblk(:), ipi(:)   
      Integer, allocatable :: jpblk(:), ipj(:)  

! ... number of blocks for given type: 

      Integer, allocatable :: nblk(:)
    
! ... blocks chain pointer for given type:       

      Integer, allocatable :: kblk(:,:)  

! ... total number of ordererd coefficients:

      Integer :: ncdata = 0  
    
      End Module cmdata


!======================================================================
      Subroutine Allocate_cmdata(k)
!======================================================================
!     allocate arrays in module "cmdata"
!----------------------------------------------------------------------
      Use cmdata
   
      Implicit none
      Integer, intent(in) :: k
      Integer :: i,n,m

      if(allocated(CLDATA)) Deallocate(CLDATA,CVDATA,K1,K2,K3,IPT, &
                                       iblk,ipblk,jpblk,ipi,ipj,nblk,kblk)
      if(k.eq.0) Return

      n = nblocks; m = nblocks*mblocks
      Allocate(CLDATA(m),CVDATA(m),K1(m),K2(m),K3(m),IPT(m), &
               iblk(n),nblk(ntype),kblk(ntype,n), &
               ipblk(n),jpblk(n),ipi(n),ipj(n) )

! ... initilize the blocks::

      Do i = 1,nblocks
       ipblk(i) = (i-1)*mblocks + 1; jpblk(i) = -1
      End do

! ... assign one block to each type:

      Do i = 1,ntype
       iblk(i) = i; nblk(i) = 1; kblk(i,1) = i; jpblk(i) = 0
      End do

      End Subroutine Allocate_cmdata
