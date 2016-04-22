!=====================================================================
!     UTILITY   I N T L S T
!
!               C O P Y R I G H T -- 2006
!
!     Written by:   Oleg Zatsarinny,   email: oleg_zoi@yahoo.com
!                         
!======================================================================
!
!     records to the file INT_LST the integral coefficients 
!     for channels interaction based on:
!
! INPUT:
!
!     channels.lst  - channel information
!     chan_even.c   - configurations for even states
!     chan_odd.c    - configurations for odd states 
!     chan_even.bnk - angular data bank for even states
!     chan_odd.bnk  - angular data bank for odd states 
!     mult_bnk      - angular data bank for dipole interaction
!
! OUTPUT:
!
!     int_lst       - angular coefficients for given case
!     
!----------------------------------------------------------------------

!=====================================================================
      Module param_intlst      
!=====================================================================
!
!     Contains specific parameters for given problem 
!
!---------------------------------------------------------------------

      Implicit none
      Save

! ... files:

      Integer(4) :: inp = 10; Character(80) :: AF_inp  = 'channels.lst'
      Integer(4) :: nco = 11; Character(80) :: AF_odd  = 'chan_odd.c'
      Integer(4) :: nce = 12; Character(80) :: AF_even = 'chan_even.c'
      Integer(4) :: nbo = 13; Character(80) :: AF_bnko = 'chan_odd.bnk'
      Integer(4) :: nbe = 14; Character(80) :: AF_bnke = 'chan_even.bnk'
      Integer(4) :: nbm = 15; Character(80) :: AF_mult = 'mult_bnk'
      Integer(4) :: out = 16; Character(80) :: AF_out  = 'channels.int'
      
! ... tollerance for coefficients:

      Real(8) :: eps_c = 1.d-7

      Integer(4) :: nch,nch1,nch2,ncomp
      Integer(4) :: nwt, nwc

      End module param_intlst




!======================================================================
!     MAIN 
!======================================================================

      Use param_intlst
      Use configs

!----------------------------------------------------------------------
! ... open files:

      Open(inp,file=AF_inp,status='OLD')
      Open(nco,file=AF_odd,status='OLD')
      Open(nce,file=AF_even,status='OLD')
      Open(nbo,file=AF_bnko,form='UNFORMATTED',status='OLD')
      Open(nbe,file=AF_bnke,form='UNFORMATTED',status='OLD')
      Open(nbm,file=AF_mult,form='UNFORMATTED',status='OLD')
      Open(out,file=AF_out)

      Call Read_rarg('c_eps',c_eps)
!----------------------------------------------------------------------
! ... read list of one-electron orbitals:

      Call Read_ipar(inp,'nwt',nwt)
      Call Read_ipar(inp,'nwc',nwc)
      nwf = nwt + nwc

      m=nwf; Call allocate_orblist(m)

      Call Read_ipar(inp,'nwt',nwt)
      read(inp,*)
      read(inp,'(20a4)') (ELF(i),i=1,nwt)
      Call Read_ipar(inp,'nwc',nwc)
      read(inp,*)
      read(inp,'(20a4)') (ELF(i),i=nwt+1,nwf)
      Do i=1,nwf
       Call EL4_NLK(ELF(i),NEF(i),LEF(i),KEF(i))
      End do

      JORT=-1; Call Pre_iort(0)

      Call Read_ipar(inp,'nch1',nch1)
      Call Read_ipar(inp,'nch2',nch2)
      nch = nch1+nch2

!----------------------------------------------------------------------
! ... processing odd configurations:

      Call R_conf(nco,nbo)
      Call R_integrals(nbo)

! ... processing even configurations:

      Call R_conf(nce,nbe)
      Call R_integrals(nbe)

! ... dipole interaction:

      Call RR_conf(nco,nce,nbm,ncfg1,ncfg2)
      Call D_integrals(nbm)

! ... output results:

      Call Output_integrals  

      End ! program INT_TAB
 


!======================================================================
      Subroutine R_integrals(nub)
!======================================================================
!
!     extracts the data from INT.BNK for specific configurations
!
!----------------------------------------------------------------------

      Use param_intlst
      USE configs 

      IMPLICIT REAL(8) (A-H,O-Z)

      Integer(4), parameter :: idet = 2**15, idef = 2**4 
      Integer(4), Dimension(:), Allocatable :: KPD,IPD,NPD
      Integer(4), Dimension(:), Allocatable :: KPF,IPF,NPF

      Integer(4), Dimension(ne) :: IPN1,IPN2,NP1,NP2,NP

!----------------------------------------------------------------------
!                                             read determinant factors:
      read(nub) ndet,kdet
      if(allocated(KPD)) Deallocate(KPD,IPD,NPD)
      Allocate(KPD(ndet),IPD(ndet),NPD(kdet))
      read(nub) KPD; read(nub) IPD; read(nub) NPD

      read(nub) ndf,kdf
      if(allocated(KPF)) Deallocate(KPF,IPF,NPF)
      Allocate(KPF(ndf),IPF(ndf),NPF(kdf))
      read(nub) KPF; read(nub) IPF; read(nub) NPF

!----------------------------------------------------------------------
!                                                  processing the data:
    1 read(nub,end=2) C,ijt,int,idf

      Call Decode_int(icase,kpol,I1,I2,I3,I4,int)

      if(icase.ne.5.and.icase.ne.6) go to 1
      if(icase.eq.6) C = -2.d0*C

! ... determine the range of states for given coeff. 

      it = ijt/ibc;   jt = mod(ijt,ibc)
      is1 = IT_state1(it); js1 = IT_state1(jt)
      is2 = IT_state2(it); js2 = IT_state2(jt)
      if(is1.eq.0.or.js1.eq.0) go to 1

      kd = 0
      if(idf.gt.0) then
       kd=KPF(idf)
       ip=IPF(idf)
       NP(1:kd)=NPF(ip+1:ip+kd)
      end if
      
! ... loop over all relevant states: 

      Do ik=is1,is2
       is=IP_stat(ik); noi=NOCCSH(is); IPN1(1:noi)=NOCORB(1:noi,is)

      Do jk=js1,js2
       js=IP_stat(jk); noj=NOCCSH(js); IPN2(1:noj)=NOCORB(1:noj,js)

! ... consider only low-half of interaction matrix

       if(it.eq.jt.and.js.gt.is) Cycle                   

! ... define integral for specific orbitals

       j1=IPN1(i1); j2=IPN1(i2); j3=IPN2(i3); j4=IPN2(i4)

! ... define the channel index: 

       ic=IT_stat(is); jc=IT_stat(js)

! ... target energies

      if(ic.eq.jc.and.ic.le.nch.and.jc.le.nch) then
       if(j1.le.nwt.and.j2.le.nwt.and.j3.le.nwt.and.j4.le.nwt) go to 1 
      end if

! ... include the expansion coefficients 

       CC = C * WC(is) * WC(js)
       if(is.ne.js.and.ic.eq.jc) CC = CC + CC
       if(abs(CC).lt.Eps_C) Cycle

! ... find overlaps for specific orbitals with extracted continuum:  

       if(idf.gt.0) then
        kz=0; jj = 1
        Do ii = 1,kd
         id=NP(ii)/idef; kdn=KPD(id); ip=IPD(id)
         Do i=1,kdn
          k=NPD(i+ip); NP1(i)=IPN1(k/idet); NP2(i)=IPN2(mod(k,idet))
         End do
         jj=IDET_SIMP(kz,kdn,NP1,NP2)
         if(jj.eq.0) Exit
         if(jj.eq.2) Stop 'IDET_SIMP = 2'
        End do 
        if(jj.eq.0) Cycle
        CC = CC * (-1)**kz
       end if

       if(icase.eq.5) Call Add_rk_data(CC,ic,jc,kpol,j1,j2,j3,j4)
       if(icase.eq.6) Call Add_L_data(CC,ic,jc,j1,j3)

      End do    !  over js
      End do	!  over is

      go to 1
    2 Continue

      End Subroutine R_integrals


!======================================================================
      Subroutine D_integrals(nub)
!======================================================================
!
!     This routine extract the dipole matrix elements
!     from MULT_BNK (unit nub) for specific configurations
!
!----------------------------------------------------------------------

      Use param_intlst
      Use configs

      Implicit real(8) (A-H,O-Z)

      Integer(4), parameter :: idet = 2**15, idef = 2**4 
      Integer(4), Dimension(:), Allocatable :: KPD,IPD,NPD
      Integer(4), Dimension(:), Allocatable :: KPF,IPF,NPF

      Integer(4), Dimension(ne) :: NP, NP1,NP2, IPN1,IPN2

!----------------------------------------------------------------------
!                                             read determinant factors:
      read(nub) ndet,kdet
      if(allocated(KPD)) Deallocate(KPD,IPD,NPD)
      Allocate(KPD(ndet),IPD(ndet),NPD(kdet))
      read(nub) KPD; read(nub) IPD; read(nub) NPD

      read(nub) ndf,kdf
      if(allocated(KPF)) Deallocate(KPF,IPF,NPF)
      Allocate(KPF(ndf),IPF(ndf),NPF(kdf))
      read(nub) KPF; read(nub) IPF; read(nub) NPF

!----------------------------------------------------------------------
!                                                  processing the data:
    1 read(nub,end=2) C, ijt,int,idf

      it=ijt/ibc;  jt=mod(ijt,ibc)

      is1=IT_state1(it);  is2=IT_state2(it)
      js1=IT_state1(jt);  js2=IT_state2(jt)

      if(is1.eq.0.or.js1.eq.0) go to 1   ! no relevant states

      kd = 0
      if(idf.gt.0) then
       kd=KPF(idf)
       ip=IPF(idf)
       NP(1:kd)=NPF(ip+1:ip+kd)
      end if

      Call Decode_mult(icase,i1,i2,int)

!----------------------------------------------------------------------
!                                                    cycle over states:
      Do ik=is1,is2
       is=IP_stat(ik); no1=NOCCSH(is); IPN1(1:no1)=NOCORB(1:no1,is)

      Do jk=js1,js2
       js=IP_stat(jk); no2=NOCCSH(js); IPN2(1:no2)=NOCORB(1:no2,js)

       if(it.eq.jt.and.ic.gt.jc) Cycle 

! ... include the expansion coefficients 

       CC = C * WC(is) * WC(js)

       if(abs(CC).lt.Eps_C) Cycle

       isign=1
       ic=IT_stat(is); jc=IT_stat(js)
       j1 = IPN1(i1); j2 = IPN2(i2)
       if(jc.gt.ic) then
        kz = ( ILterm(is) - ILterm(js) + ISterm(is) - ISterm(js) )/2 
        isign=(-1)**kz
        ic = IT_stat(js); jc=IT_stat(is)
        j1 = IPN2(i2); j2 = IPN1(i1)
       end if

! ... check the overlap factor:  

       if(idf.gt.0) then
        kz=0; jj = 1
        Do ii = 1,kd
         id=NP(ii)/idef; kdn=KPD(id); ip=IPD(id)
         Do i=1,kdn
          k=NPD(i+ip); NP1(i)=IPN1(k/idet); NP2(i)=IPN2(mod(k,idet))
         End do
         jj=IDET_SIMP(kz,kdn,NP1,NP2)
         if(jj.eq.0) Exit
         if(jj.eq.2) Stop 'IDET_SIMP = 2'
        End do 
        if(jj.eq.0) Cycle
        CC = CC * (-1)**kz
       end if

       Call Add_D_data(CC,ic,jc,j1,j2)

      End do   ! over jc
      End do   ! over ic
      
      go to 1
    2 Continue

    
      End Subroutine D_integrals



!======================================================================
      Subroutine Output_integrals
!======================================================================

      Use configs; Use rk_data; Use L_data; Use D_data
      Use param_intlst

      IMPLICIT REAL(8) (A-H,O-Z)

! ... L-data:

      m = nlk
      Do i = 1,nlk
       if(abs(clk(i)).lt.eps_c) m=m-1
      End do

      write(out,'(i6,a)') m,'   -   number of L-integrals '       
      write(*,'(i6,a)') m,'   -   number of L-integrals '       
      Do i = 1,nlk
       if(abs(clk(i)).lt.eps_c) Cycle
!       ic = kl1(i)/ibl; jc = mod(kl1(i),ibl)        
       Call Get_ij(kl1(i),ic,jc)
       i1 = kl2(i)/ibl; i2 = mod(kl2(i),ibl)
       write(out,'(f16.8,2i6,3x,a2,a4,a1,a4,a1)') &
        clk(i),ic,jc,'L(',ELF(i1),';',ELF(i2),')'
      End do
       
! ... Rk-data:

      m = nrk
      Do i = 1,nrk
       if(abs(crk(i)).lt.eps_c) m=m-1
      End do

      write(out,'(i6,a)') m,'   -   number of Rk-integrals '       
      write(*,'(i6,a)') m,'   -   number of Rk-integrals '       
      Do i = 1,nrk
       if(abs(crk(i)).lt.eps_c) Cycle
!       ic = kr1(i)/ibr; jc = mod(kr1(i),ibr)        
       Call Get_ij(kr1(i),ic,jc)
       k  = kr2(i)
       i1 = kr3(i)/ibr; i2 = mod(kr3(i),ibr)
       i3 = kr4(i)/ibr; i4 = mod(kr4(i),ibr)
       write(out,'(f16.8,2i6,3x,a1,i2,a1,4(a4,a1))') &
        crk(i),ic,jc,'R',k,'(',ELF(i1),',',ELF(i2),';', &
                               ELF(i3),',',ELF(i4),')'
      End do

! ... D-data:

      m = ndk
      Do i = 1,ndk
       if(abs(cdk(i)).lt.eps_c) m=m-1
      End do

      write(out,'(i6,a)') m,'   -   number of D-integrals '       
      write(*,'(i6,a)') m,'   -   number of D-integrals '       
      Do i = 1,ndk
       if(abs(cdk(i)).lt.eps_c) Cycle
!       ic = kd1(i)/ibd; jc = mod(kd1(i),ibd)        
       Call Get_ij(kd1(i),ic,jc)
       i1 = kd2(i)/ibd; i2 = mod(kd2(i),ibd)
       write(out,'(f16.8,2i6,3x,a2,a4,a1,a4,a1)') &
        cdk(i),ic,jc,'d(',ELF(i1),';',ELF(i2),')'
      End do

      End Subroutine Output_integrals


!======================================================================
      Subroutine Get_ij(k,i,j)
!======================================================================

      i = sqrt(real(2*k))
      j = k - i*(i-1)/2
      if(j.gt.i) then
       i=i+1
       j = k - i*(i-1)/2
      end if

      End Subroutine Get_ij




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


!----------------------------------------------------------------------
      Integer(4) FUNCTION IDET_SIMP(kz,kn,N1,N2)
!----------------------------------------------------------------------
!     simplify the determinant kn,N1,N2
!     accoding the orthogonality conditions for radial w.f.'s
!     kz - number of needed permutations
!     IDET_SIMP = 0,1,2 with overlap determinant = 0,1 or some value  
!----------------------------------------------------------------------

      Use configs, ONLY: IORT

      Implicit none

      Integer(4) :: kz,kn
      Integer(4), Dimension(kn) :: N1(kn),N2(kn)
      Integer(4) :: i,ii,i1,i2, k,kk,k1,k2, m1,m2 

      if(kn.le.0) Stop ' IDET_SIMP: kn <= 0'

      IDET_SIMP=0

!----------------------------------------------------------------------
!                       Check for a row with only one non-zero element:
    1  Do i1=1,kn                
       k=0                      
       Do i2=1,kn
        m1=max(N1(i1),N2(i2)); m2=min(N1(i1),N2(i2)); ii=IORT(m1,m2)
        if(ii.ne.0) then
         k=k+1; kk=ii;  if(k.gt.1.or.kk.ne.1) Exit
         k1=i1; k2=i2
        end if
       End do
       if(k.eq.0) Return; if(k.eq.1.and.kk.eq.1) go to 2
      End do

!----------------------------------------------------------------------
!                   Check for a colum with only one non-zero element:

      Do i2=1,kn                
       k=0                    
       Do i1=1,kn
         m1=max(N1(i1),N2(i2)); m2=min(N1(i1),N2(i2)); ii=IORT(m1,m2)
        if(ii.ne.0) then
         k=k+1; kk=ii; if(k.gt.1.or.kk.ne.1) Exit
         k1=i1; k2=i2
        end if
       End do
       if(k.eq.0) Return; if(k.eq.1.and.kk.eq.1) go to 2
      End do

      go to 3
!-----------------------------------------------------------------------
!                                                 the case of <k1|k2>=1:
    2 kn=kn-1                       
      if(kn.eq.0) then
       IDET_SIMP=1; Return
      end if
      kz=kz+k1+k2
      Do i=k1,kn; N1(i)=N1(i+1); End do
      Do i=k2,kn; N2(i)=N2(i+1); End do
      go to 1
!-----------------------------------------------------------------------
!                                                  ordering of elements:
    3 Continue                     

      Do i1=1,kn-1
       Do i2=i1+1,kn
        if(N1(i1).gt.N1(i2)) then
         kk=N1(i1);  N1(i1)=N1(i2); N1(i2)=kk; kz=kz+1
        end if
        if(N2(i1).gt.N2(i2)) then
         kk=N2(i1);  N2(i1)=N2(i2); N2(i2)=kk; kz=kz+1
        end if
       End do
      End do

      IDET_SIMP=2

      End Function IDET_SIMP


!======================================================================
      Integer(4) Function Incode_INT (itype,k,I1,I2,I3,I4)
!======================================================================
!
!     incode the integral
!
!----------------------------------------------------------------------
 
      Use configs, ONLY: jb,jb4,jb8

      Implicit none
      Integer(4), intent(in) :: itype,k,I1,I2,I3,I4       

      if(max(i1,i2,i3,i4).ge.jb) Stop ' Incode_int: i > ibase '
      if(itype.ge.21) Stop 'Incode_int: itype too large'

      Incode_INT = ((i1*jb+i2)*jb+i3)*jb+i4 + k*jb4 + itype*jb8

      End Function Incode_INT



!======================================================================
      Subroutine Decode_INT (itype,k,I1,I2,I3,I4,int)
!======================================================================
!
!     decode the integral
!
!----------------------------------------------------------------------
 
      Use configs, ONLY: jb,jb8

      Implicit none
      Integer(4), Intent(in) :: int
      Integer(4), Intent(out) :: itype,k,I1,I2,I3,I4

      K  = int;   itype = k/jb8;  k = mod(k,jb8)
      I4 = mod(K,jb);  K  =  K/jb
      I3 = mod(K,jb);  K  =  K/jb
      I2 = mod(K,jb);  K  =  K/jb
      I1 = mod(K,jb);  K  =  K/jb

      End Subroutine Decode_INT



!======================================================================
      Integer(4) Function Incode_mult (itype,i1,i2)
!======================================================================
!
!     incode the integral
!
!----------------------------------------------------------------------
 
      Use configs, ONLY: jb

      Implicit none
      Integer(4), intent(in) :: itype,i1,i2       

      if(max(i1,i2).ge.jb) then
       write(*,'(a,2i5,a,i3)')  &
               ' Incode_mult: i1,i2 =',i1,i2,'  > base =',jb
       Stop
      end if
      if(itype.gt.3) then
       write(*,'(a,i5,a)')  &
               ' Incode_mult: itype =',itype,'  is out of range (3)'
       Stop 
      end if

      Incode_mult = (itype*jb+i1)*jb+i2

      End Function Incode_mult


!======================================================================
      Subroutine Decode_mult (itype,i1,i2,int)
!======================================================================
!
!     decode the integral
!
!----------------------------------------------------------------------
 
      Use configs, ONLY: jb

      Implicit none
      Integer(4), Intent(in) :: int
      Integer(4), Intent(out) :: itype,i1,i2 
      Integer(4) :: k

      k  = int
      i2 = mod(k,jb);  k = k/jb
      i1 = mod(k,jb);  k = k/jb
      itype = mod(k,jb)

      End Subroutine Decode_mult
