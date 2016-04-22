!=======================================================================
      Subroutine RR_conf(in,nui)
!=======================================================================
!
!     Read the configuration list from file c-file (unit 'in'),
!     define the corresponding angular and configurational symmetries
!     and compare with information from 'INT.BNK' (unit 'nui')
!
!-----------------------------------------------------------------------
      
      USE configs

      IMPLICIT NONE
      
      Integer(4), Intent(in) :: in,nui

      Integer(4), parameter :: mrecl = 100000

      Character(26) ::  AI   ! config. symmetry
      Character(40) ::  BI   ! angular symmetry
      Character(26), DIMENSION(:), ALLOCATABLE ::  ASYM
      Character(40), DIMENSION(:), ALLOCATABLE ::  BSYM
      Character(20) :: AF

      INTEGER(4), DIMENSION(:), ALLOCATABLE :: IT_conf, IS_term

!     IT_conf(it) gives the configuration for given term 'it'
!     IS_term(is) defines the term for given stats 'is'

      Integer(4), External :: Idef_ncfg, Iadd_orb
      Integer(4) :: i,j, ic, it, is, ierr, iconf 

!----------------------------------------------------------------------
!                                            define ncfg, nsymt, nsymc:
      ncfg=Idef_ncfg(in)
      if(ncfg.eq.0) Stop ' RZ_conf: ncfg = 0 '

      if(allocated(NOCCSH)) Deallocate(NOCCSH, NELCSH, NOCORB, &
                                       LSA, LSL, LSS, LPL, LPS) 
      Allocate(NOCCSH(ncfg), NELCSH(nsh,ncfg), NOCORB(nsh,ncfg), &
               LSA(nsh,ncfg), LSL(nsh,ncfg), LSS(nsh,ncfg),      &
               LPL(nsh,ncfg), LPS(nsh,ncfg))

      read(nui) nsymt, nsymc
      
      Allocate(ASYM(nsymc), BSYM(nsymt), IT_conf(nsymt),IS_term(ncfg), &
               STAT = ierr)
      if(ierr.ne.0) Stop ' Problem with allocation in RZ_conf'

      Call R_a (nui,mrecl,nsymc,ASYM)    ! read(nub) ASYM(1:nsymc)
      Call R_a (nui,mrecl,nsymt,BSYM)    ! read(nub) BSYM(1:nsymt)
      Call R_i4(nui,mrecl,nsymt,IT_conf) ! read(nub) IT_conf(1:nsymt)

!----------------------------------------------------------------------
!                                   define orbitals and configurations:
      rewind(in)
      ic = 0
    1 read(in,'(a)',end=2) CONFIG
      if(CONFIG(1:1).eq.'*') go to 2
      if(CONFIG(5:5).ne.'(') go to 1
      read(in,'(a)') COUPLE
      Call Decode_c
      ic = ic + 1

      noccsh(ic)=no
      Do i=1,no
       nocorb(i,ic)=Iadd_orb(nn(i),ln(i),kn(i))
       nelcsh(i,ic)=iq(i)
       LSA(i,ic) = LS(i,1)
       LSL(i,ic) = LS(i,2)
       LSS(i,ic) = LS(i,3)
       LPL(i,ic) = LS(i,4)
       LPS(i,ic) = LS(i,5)
      End do

       Call Incode_conf(AI,BI)

       iconf=0
       Do i=1,nsymc
        if(AI.ne.ASYM(i)) Cycle
        iconf=i
        Exit
       End do
       if(iconf.eq.0) Stop ' RZ_conf: unknown configuration'
       
       IS_term(ic)=0
       Do i=1,nsymt
        if(BI.ne.BSYM(i) .or. IT_conf(i).ne.iconf) Cycle
        IS_term(ic)=i
        Exit
       End do
       if(IS_term(ic).eq.0) Stop ' RZ_conf: unknown ang.symmetry'

       go to 1
 
    2 if(ic.ne.ncfg) Stop ' RZ_conf: ic <> ncfh '
!----------------------------------------------------------------------
!                                           define IP_term and IC_term:

      Allocate(IP_term(nsymt), IC_term(nsymc), STAT = ierr)
      if(ierr.ne.0) Stop ' Problem with allocation in RZ_conf'

      Call SortI(nsymt,IT_conf,IP_term)

      IC_term = 0

      Do i=1,nsymt
       it=IP_term(i)
       ic=IT_conf(it)
       if(IC_term(ic).eq.0) then
        IC_term(ic) = i + ibc*i
       else
        j=mod(IC_term(ic),ibc)
        IC_term(ic) = j + ibc*i
       end if
      End do
!----------------------------------------------------------------------
!                                           define IP_stat and IT_stat:

      Allocate(IP_stat(ncfg), IT_stat(nsymt), STAT = ierr)
      if(ierr.ne.0) Stop ' Problem with allocation in RZ_conf'

      Call SortI(ncfg,IS_term,IP_stat)

      IT_stat = 0

      Do i=1,ncfg
       is=IP_stat(i)
       it=IS_term(is)
       if(IT_stat(it).eq.0) then
        IT_stat(it) = i + ibc*i
       else
        j=mod(IT_stat(it),ibc)
        IT_stat(it) = j + ibc*i
       end if
      End do

      Deallocate(ASYM, BSYM, IT_conf, IS_term, STAT = ierr)
      if(ierr.ne.0) Stop ' Problem with deallocation in RZ_conf'

      End Subroutine RR_conf


!====================================================================
      Subroutine Incode_conf(ASYM,BSYM)
!====================================================================
!
!     incodes the conf. and angular symmetries (ASYM and BSYM)
!     (ASYM also includes the total term!)
!--------------------------------------------------------------------

      Use configs

      IMPLICIT NONE

      Character(26), Intent(out) :: ASYM
      Character(40), Intent(out) :: BSYM
      Character(1), EXTERNAL :: AL
 
      Integer(1) :: i,k,m

      k=1
      m=1
      Do i=1,no
       write(ASYM(k:k+2),'(a1,i2)') AL(ln(i),1),iq(i)
       k=k+3
       write(BSYM(m:m+4),'(2a1,i1,2a1)') &
                    AL(LS(i,3),2),AL(LS(i,2),6),LS(i,1), &
                    AL(LS(i,5),2),AL(LS(i,4),6)
       m=m+5
      End do

      Do i=no+1,nsh
       ASYM(k:k+2)='   '
       k=k+3
       BSYM(m:m+4)='     '
       m=m+5
      End do

      ASYM(25:25) = AL(LS(no,4),6)
      ASYM(26:26) = AL(LS(no,5),2)

      End Subroutine Incode_conf


