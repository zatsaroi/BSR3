!====================================================================
      MODULE cmdata
!====================================================================
!  contains a set of coefficients with four identifiers (k1,k2,k3,k4)
!  and one order pointer (ipt)
!--------------------------------------------------------------------
      Implicit none 

      Integer :: ncdata = 0            ! number of coefficients

      Real(8), allocatable :: CDATA(:) ! coefficients    

! ... their attributes:

      Integer, allocatable :: K1(:),K2(:),K3(:),K4(:),IPT(:)

! ... the data are divided on blocks:  

      Integer :: nblocks =  0       !  number of blocks
      Integer :: mblocks =  0       !  size of blocks
      Integer :: kblocks =  0       !  max.number of blocks for given kpol

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

      Real(8) :: timeg(ntype)
      Integer :: ibcoef(ntype), igen = 0, jgen = 0
      Real(8) :: Tadd=0.d0, Tcoef=0.d0, Tmerge=0.d0 
 
      End MODULE cmdata


!======================================================================
     Subroutine Allocate_cmdata(nb,mb,kb,m)
!======================================================================
!    allocate or de-allocate memory for "cmdata" list
!----------------------------------------------------------------------
     USE cmdata
     Implicit none
     Integer, intent(in) :: nb,mb,kb
     Integer :: m   

     m=nb*mb;  nblocks=nb; mblocks=mb; kblocks=kb

     if(allocated(CDATA)) Deallocate(CDATA,K1,K2,K3,K4,IPT, &
                           ipblk,jpblk,ipi,ipj,iblk,nblk,kblk)
     if(m.gt.0) & 
       Allocate(CDATA(m),K1(m),K2(m),K3(m),K4(m),IPT(m), &
                ipblk(nb),jpblk(nb),ipi(nb),ipj(nb), &
                iblk(-1:npol,ntype),nblk(-1:npol,ntype), &
                kblk(-1:npol,ntype,nb) )
     m = 7*m + 4*nb + (2+nb)*ntype*(npol+2) 
	        
     END Subroutine Allocate_cmdata
