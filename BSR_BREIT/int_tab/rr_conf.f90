!=======================================================================
      Subroutine R_conf(nuc,nub)
!=======================================================================
!
!     Read the configuration list from file c-file (unit 'in'),
!     define the corresponding angular and configurational symmetries
!     and compare with information from 'INT.BNK' (unit 'nui')
!
!-----------------------------------------------------------------------
      
      USE configs

      IMPLICIT NONE
      
      Integer(4), Intent(in) :: nuc,nub

      Integer(4), parameter :: mrecl = 100000
      Integer(4), parameter :: noper = 7
      Integer(1), Allocatable, Dimension(:,:) :: IT_oper

      Character(26) ::  AI   ! config. symmetry
      Character(40) ::  BI   ! angular symmetry
      Character(26), DIMENSION(:), ALLOCATABLE ::  ASYM
      Character(40), DIMENSION(:), ALLOCATABLE ::  BSYM
      Character(80) :: AS

      INTEGER(4), DIMENSION(:), ALLOCATABLE :: IT_conf, IS_term

!     IT_conf(it) gives the configuration for given term 'it'
!     IS_term(is) defines the term for given stats 'is'

      Integer(4), External :: Idef_ncfg, Jfind_orb, Idef_ne
      Integer(4) :: i,j, ic, it, is, iconf 

!----------------------------------------------------------------------
!                                            define ncfg, nsymt, nsymc:
      ncfg=Idef_ncfg(nuc)
      if(ncfg.eq.0) Return

      if(allocated(NOCCSH)) Deallocate(NOCCSH)
      Allocate(NOCCSH(ncfg))
      if(allocated(NOCORB)) Deallocate(NOCORB)
      Allocate(NOCORB(nsh,ncfg))
      if(allocated(WC)) Deallocate(WC)
      Allocate(WC(ncfg))
      if(allocated(IT_STAT)) Deallocate(IT_STAT)
      Allocate(IT_STAT(ncfg))

      read(nub) nsymt, nsymc
      
      Allocate(ASYM(nsymc),BSYM(nsymt),IT_conf(nsymt),IS_term(ncfg))

      Call R_a (nub,mrecl,nsymc,ASYM)    ! read(nub) ASYM(1:nsymc)
      Call R_a (nub,mrecl,nsymt,BSYM)    ! read(nub) BSYM(1:nsymt)
      Call R_i4(nub,mrecl,nsymt,IT_conf) ! read(nub) IT_conf(1:nsymt)

!----------------------------------------------------------------------
!                                   define orbitals and configurations:
      ne = Idef_ne(nuc)
      rewind(nuc)
      ic = 0
    1 read(nuc,'(a)',end=2) AS
      if(AS(1:1).eq.'*') go to 2
      if(AS(5:5).ne.'(') go to 1
      ic = ic + 1
      read(AS,'(a64,f11.8,i5)') CONFIG,WC(ic),IT_STAT(ic)
      read(nuc,'(a)') COUPLE
      Call Decode_c

      noccsh(ic)=no
      Do i=1,no
       nocorb(i,ic)=Jfind_orb(nn(i),ln(i),kn(i))
      End do

       Call Incode_conf(AI,BI)

       iconf=0
       Do i=1,nsymc
        if(AI.ne.ASYM(i)) Cycle
        iconf=i
        Exit
       End do
       if(iconf.eq.0) Stop ' INT_BNK is incomplete !'
       
       IS_term(ic)=0
       Do i=1,nsymt
        if(BI.ne.BSYM(i) .or. IT_conf(i).ne.iconf) Cycle
        IS_term(ic)=i
        Exit
       End do
       if(IS_term(ic).eq.0) Stop ' INT_BNK is incomplete !'

       go to 1
 
    2 if(ic.ne.ncfg) Stop ' RZ_conf: ic <> ncfh '

!----------------------------------------------------------------------
!                                           define IP_stat and IT_stat:

      if(allocated(IP_stat)) Deallocate(IP_stat)
      if(allocated(IT_state1)) Deallocate(IT_state1)
      if(allocated(IT_state2)) Deallocate(IT_state2)

      Allocate(IP_stat(ncfg),IT_state1(nsymt),IT_state2(nsymt))

      IP_stat=0; IT_state1=0; IT_state2=0

      Call Isort1(ncfg,IS_term,IP_stat)

      Do i=1,ncfg
       it=IS_term(IP_stat(i))
       if(IT_state1(it).eq.0) IT_state1(it)=i; IT_state2(it)=i 
      End do

      Deallocate(ASYM, BSYM, IT_conf, IS_term)

      read(nub) i           
      Allocate(IT_oper(noper,i))
      Call R_i1(nub,mrecl,noper,i,i,IT_oper)
      Deallocate(IT_oper)


      End Subroutine R_conf




!======================================================================
      Subroutine RR_conf(in1,in2,nub,ncfg1,ncfg2)
!======================================================================
!
!     Read the configurations from c-files (units in1,in2),
!     define there angular symmetries and compare with existing ones
!     in mult_bnk (unit nub). Prepare some angular arrays.
!
!----------------------------------------------------------------------

      USE configs

      Implicit real(8) (A-H,O-Z)

      Character(26) :: AI  
      Character(40) :: BI  

      Character(26), Allocatable :: ASYM(:)   ! config. symmetry
      Character(40), Allocatable :: BSYM(:)   ! angular symmetry

      Character(80) :: AS
      Character(1) :: ktype

! ... IT_conf(it) gives the configuration for given term 'it'
! ... IS_term(is) gives the term for given state 'is'
! ... IT_done gives indication on done calculations

      Integer(4), Allocatable :: IT_conf(:),IS_term(:)
      Integer(1) :: iii

!---------------------------------------------------------------------

      read(nub) ktype,kkpol

      read(nub) nsymt, nsymc

      if(Allocated(ASYM)) Deallocate(ASYM)
      if(Allocated(BSYM)) Deallocate(BSYM)

      Allocate(ASYM(nsymc),BSYM(nsymt),IT_conf(nsymt))

      read(nub) ASYM(1:nsymc)
      read(nub) BSYM(1:nsymt)
      read(nub) IT_conf(1:nsymt)

      ncfg1=Idef_ncfg(in1)
      ncfg2=Idef_ncfg(in2)
      ncfg = ncfg1+ncfg2

      if(Allocated(NOCCSH)) Deallocate(NOCCSH)
      if(Allocated(NOCORB)) Deallocate(NOCORB)

      Allocate(NOCCSH(ncfg),NOCORB(nsh,ncfg),IS_term(ncfg))

      if(allocated(WC)) Deallocate(WC)
      Allocate(WC(ncfg))
      if(allocated(IT_STAT)) Deallocate(IT_STAT)
      Allocate(IT_STAT(ncfg))
      if(allocated(ILterm)) Deallocate(ILterm)
      Allocate(ILterm(ncfg))
      if(allocated(ISterm)) Deallocate(ISterm)
      Allocate(ISterm(ncfg))

!---------------------------------------------------------------------
! ... define symmetries in first set:

      rewind(in1)
      parity = 0; ic = 0
      Do 

       read(in1,'(a)') AS
       if(AS(1:1).eq.'*') Exit
       if(AS(5:5).ne.'(') Cycle
       ic = ic + 1
       read(AS,'(a64,f11.8,i5)') CONFIG,WC(ic),IT_STAT(ic)
       read(in1,'(a)') COUPLE
 
       Call Decode_c

! ... check if angular symmetry is already in bank:

       Call Incode_conf(AI,BI)
       iconf=0
       Do i=1,nsymc
        if(AI.ne.ASYM(i)) Cycle; iconf=i; Exit
       End do
       if(iconf.eq.0) Stop ' MULT_BNK is incomplete !'
       jterm = 0
       Do i=1,nsymt
        if(BI.ne.BSYM(i) .or. IT_conf(i).ne.iconf) Cycle
        jterm=i; Exit
       End do
       if(jterm.eq.0) Stop ' MULT_BNK is incomplete !'
       
! ... save some spectroscopic data:

       NOCCSH(ic) = no
	   Do i = 1,no
	    NOCORB(i,ic) = Jfind_orb(nn(i),ln(i),kn(i))
	   End do	
       IS_term(ic)=jterm

      End do  ! on ic


!---------------------------------------------------------------------
! ... define symmetries in second set:

      parity = 0
      rewind(in2)

      Do 

       read(in2,'(a)') AS
       if(AS(1:1).eq.'*') Exit
       if(AS(5:5).ne.'(') Cycle
       ic = ic + 1
       read(AS,'(a64,f11.8,i5)') CONFIG,WC(ic),IT_STAT(ic)
       read(in2,'(a)') COUPLE

       Call Decode_c

! ... check if angular symmetry is already in bank:

       Call Incode_conf(AI,BI)

       iconf=0
       Do i=1,nsymc
        if(AI.ne.ASYM(i)) Cycle; iconf=i; Exit
       End do
       if(iconf.eq.0) Stop ' MULT_BNK is incomplete !'
       jterm = 0
       Do i=1,nsymt
        if(BI.ne.BSYM(i) .or. IT_conf(i).ne.iconf) Cycle
        jterm=i; Exit
       End do
       if(jterm.eq.0) Stop ' MULT_BNK is incomplete !'

! ... save some spectroscopic data:

       NOCCSH(ic) = no
	   Do i = 1,no
	    NOCORB(i,ic) = Jfind_orb(nn(i),ln(i),kn(i))
	   End do	
       IS_term(ic)=jterm
       ILterm(ic) = LS(no,4); ISterm(ic) = LS(no,5)

      End do  ! on ic

      if(ic.ne.ncfg) Stop 'R_conf: ic <> ncfg '

!----------------------------------------------------------------------
!                                           define IP_stat and IT_stat:
     
      if(allocated(IP_stat))   Deallocate(IP_stat)
      if(allocated(IT_state1)) Deallocate(IT_state1)
      if(allocated(IT_state2)) Deallocate(IT_state2)

      Allocate(IP_stat(ncfg), IT_state1(nsymt), IT_state2(nsymt))

      IP_stat=0; IT_state1=0; IT_state2=0

      Call Isort1(ncfg,IS_term,IP_stat)

      Do i=1,ncfg
       it=IS_term(IP_stat(i))
       if(IT_state1(it).eq.0) IT_state1(it)=i
       IT_state2(it)=i 
      End do
    
      Deallocate(IT_conf, IS_term, ASYM, BSYM)

! ... skip some imformation:

      read(nub) (iii,i=1,nsymt*(nsymt+1)/2)

      End Subroutine RR_conf

