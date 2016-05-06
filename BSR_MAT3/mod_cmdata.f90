!====================================================================
      MODULE cmdata
!====================================================================
!     contains a set of coefficients with four identifiers
!     (k1,k2,k3,k4) and one pointer (ipt)
!     other arrays define the decomposition of coefficients into 
!     blocks according corresponding integral type 
!--------------------------------------------------------------------
      Use bsr_mat

      Implicit none 
      Integer :: ncdata = 0                 ! number of coefficients
      Real(8), allocatable :: CDATA(:)      ! coefficients    
    
! ... their attributes:

      Integer, allocatable :: K1(:),K2(:),K3(:),K4(:),IPT(:)

! ... the data are divided on blocks: 

      Integer :: nblocks =  0       !  number of blocks
      Integer :: mblocks =  0       !  size of blocks
      Integer :: kblocks =  0       !  number of blocks for one type

! ... pointer on first(last) element in given block 

      Integer, allocatable :: ipblk(:), ipi(:)   
      Integer, allocatable :: jpblk(:), ipj(:)  

! ... current block for given type and multipole:

      Integer, allocatable :: iblk(:,:)    

! ... number of blocks for given type and multipole: 

      Integer, allocatable :: nblk(:,:)
    
! ... ponter for chain of blocks for given type:       

      Integer, allocatable :: kblk(:,:,:)

! ... structure of data:

      Integer, parameter :: ncase = 11 ! number of different cases
      Integer, parameter :: ntype =  9 ! number of different structures
      Integer            :: npol  =  0 ! number of different multipoles

      Integer, parameter :: ibi = 2**15 ! packing basis for orbitals

      Integer :: icase =  0     !  current case    
      Integer :: itype =  0     !  current structure    
      Integer :: kpol  =  0     !  current multipole index                   
    
! ... IJCASE  -  correspondence between icase and itype

      Integer IJCASE(ntype)  /2,1,2,1,2,1,2,3,4/

! ... AINT   -  designation for different type of integrals

      Character AINT(ncase)/' ',' ','T','M','R','L','Z','N','V','N','O'/

! ... debug information:

      Real(8) :: timeg (ntype)
      Integer :: ibcoef(ntype)
      Real(8) :: Tadd=0.d0, Tcoef=0.d0, Tmerge=0.d0 
 
      End MODULE cmdata


!======================================================================
      Subroutine Allocate_cmdata
!======================================================================
!     allocate arrays in module cmdata
!----------------------------------------------------------------------
      Use cmdata
      Use bsr_mat, only: nb,mb,kb
   
      m=nb*mb;  nblocks=nb; mblocks=mb; kblocks=kb

      if(m.eq.0) then
       if(allocated(CDATA)) Deallocate(CDATA,K1,K2,K3,K4,IPT, &
                            ipblk,jpblk,ipi,ipj,iblk,nblk,kblk)
      elseif(.not.allocated(CDATA)) then
       Allocate(CDATA(m),K1(m),K2(m),K3(m),K4(m),IPT(m), &
                ipblk(nb),jpblk(nb),ipi(nb),ipj(nb), &
                iblk(-1:npol,ntype),nblk(-1:npol,ntype), &
                kblk(-1:npol,ntype,nb) )
      else
       Deallocate(CDATA,K1,K2,K3,K4,IPT, &
                  ipblk,jpblk,ipi,ipj,iblk,nblk,kblk)
       Allocate(CDATA(m),K1(m),K2(m),K3(m),K4(m),IPT(m), &
                ipblk(nb),jpblk(nb),ipi(nb),ipj(nb), &
                iblk(-1:npol,ntype),nblk(-1:npol,ntype), &
                kblk(-1:npol,ntype,nb) )
      end if
	        
      End Subroutine Allocate_cmdata


!======================================================================
      Subroutine Initilize_cmdata
!======================================================================
!     clear the arrays from previous calculations
!----------------------------------------------------------------------
      Use cmdata

      Implicit none
      Integer :: i,j,k

! ... nulify blocks:

      Do i=1,nblocks;  ipblk(i)=(i-1)*mblocks+1; jpblk(i)=-1; End do

! ... assign one block to each type:

      i = 0
      Do k = -1,npol
       Do j = 1,ntype
        i = i + 1; iblk(k,j)=i; nblk(k,j)=1; kblk(k,j,1)=i; jpblk(i)=0
       End do
      End do

      End Subroutine Initilize_cmdata
     
