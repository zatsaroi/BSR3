!=====================================================================
!     PROGRAM   test_coef_ee_1conf                      
!
!               C O P Y R I G H T -- 2013
!
!     Written by:   Oleg Zatsarinny
!                   email: oleg_zoi@yahoo.com
!======================================================================
!
!    it is a debug program to check subroutine     "coef_1conf" 
!    to generate angular coefficient for Breit_Pauli calculations 
!    for one-configuration in case of orthogonal one-electron 
!    radial functions
!
!----------------------------------------------------------------------
!
!    INPUT ARGUMENTS:  name  or   cfg=...  tab=...  
!    
!----------------------------------------------------------------------
!
!    INPUT FILE:     AF_cfg    (default  - cfg.inp) 
!    OUTPUT FILES:   AF_tab    (default  - coef.tab)
!    
!---------------------------------------------------------------------     

      Implicit none 

      Integer :: nuc=1; Character(40) :: AF_cfg = 'cfg.inp'
      Integer :: out=2; Character(40) :: AF_tab = 'coef.tab'
      Character(40) :: name = ' '

      Integer, parameter :: msh = 8
      Integer :: no, nn(msh),ln(msh),iq(msh),kn(msh), LS(5,msh)

      Character(8*msh) :: CONFIG, COUPLE
      Real(8) :: C, eps_C = 1.d-7
      Integer :: i,j, ki,kj, k, ic, ncfg
      Integer, external :: Idef_ncfg 

! ... result coefficients:

      Real(8), allocatable :: coefs(:,:,:)
      Integer :: kmax
!----------------------------------------------------------------------
      Call Read_name(name)
      if(len_trim(name).ne.0) then
       AF_cfg=trim(name)//'.c'
       AF_tab=trim(name)//'.tab'
      else
       Call Read_aarg('cfg',AF_cfg)
       Call Read_aarg('tab',AF_tab)
      end if

      Call Check_file(AF_cfg)
      Open(nuc,file=AF_cfg)
      ncfg = Idef_ncfg(nuc)
      if(ncfg.eq.0) Stop 'ncfg=0: nothing to do'
      Open(out,file=AF_tab)

      rewind(nuc)
      Do ic = 1,ncfg

      Do 
       read(nuc,'(a)') CONFIG
       if(CONFIG(5:5).ne.'(') Cycle
       read(nuc,'(a)') COUPLE
       Exit
      End do

      Call Decode_cn (CONFIG,COUPLE,no,nn,ln,iq,kn,LS)

      kmax = 2*maxval(ln(1:no))
      if(allocated(coefs)) Deallocate(coefs)
      Allocate(coefs(no,no,0:kmax))

      Call coef_ee_1conf(no,ln,iq,LS,kmax,coefs)

      write(out,'(a)') trim(CONFIG) 
      write(out,'(a)') trim(COUPLE)
      write(out,*)

      Do i=1,no;  ki = 2 + (i-1)*8
       Do j=i,no; kj = 2 + (j-1)*8
        Do k=0,kmax
		       C = coefs(i,j,k); if(abs(C).lt.eps_c) Cycle
         write(out,'(a,i2,a,a,a,a,a,f10.5)') &
		            'F',k,'(',CONFIG(ki:ki+2),',',CONFIG(kj:kj+2),')=',C        
        End do
        if(i.eq.j) Cycle
        Do k=0,kmax
		       C = coefs(j,i,k); if(abs(C).lt.eps_c) Cycle
         write(out,'(a,i2,a,a,a,a,a,f10.5)') &
		            'G',k,'(',CONFIG(ki:ki+2),',',CONFIG(kj:kj+2),')=',C        
        End do
       End do; End do
       write(out,*)

      End do  ! over ic

      END ! Program ...

!======================================================================
       Subroutine Decode_cn(CONFIG,COUPLE,no,nn,ln,iq,kn,LS)
!======================================================================
!      decode the config. from c-file format to ZOI format
!----------------------------------------------------------------------
       Implicit none
       Character(*), intent(in) :: CONFIG,COUPLE
       Integer :: no,nn(*),ln(*),iq(*),kn(*),LS(5,*)
       Integer :: ii, k,i,j 
       Integer, external :: LA

       no=0; ii=LEN_TRIM(CONFIG); ii=ii/8;  k=1; j=2
       Do i=1,ii
        if(CONFIG(k+4:k+4).ne.'(') Exit
        no=no+1
        Call EL4_nlk(CONFIG(k:k+3),nn(i),ln(i),kn(i))
        read(CONFIG(k+5:k+6),'(i2)') iq(i)
        k=k+8
        read(COUPLE(j:j),'(i1)') LS(3,i)
        LS(2,i)=2*LA(COUPLE(j+1:j+1))+1
        read(COUPLE(j+2:j+2),'(i1)') LS(1,i)
        j=j+4
       End do

       LS(4,1)=LS(2,1)
       LS(5,1)=LS(3,1)
       Do i=2,no
        read(COUPLE(j:j),'(i1)') LS(5,i)
        LS(4,i)=2*LA(COUPLE(j+1:j+1))+1
        j=j+4
       End do

       End Subroutine Decode_cn

!======================================================================
      Integer Function Idef_ncfg(nu)
!======================================================================
!     gives the number of configuration in c-file (unit nu)
!----------------------------------------------------------------------
      Implicit none
      Integer, intent(in) :: nu
      Integer :: ncfg
      Character(5) :: AS

      rewind(nu)
      ncfg=0
    1 read(nu,'(a)',end=2) AS
      if(AS(1:1).eq.'*') go to 2
      if(AS(5:5).ne.'(') go to 1
      read(nu,'(a)') AS
      ncfg=ncfg+1
      go to 1
    2 rewind(nu)
      Idef_ncfg=ncfg

      End Function Idef_ncfg

!====================================================================
      Integer FUNCTION LA(a)
!====================================================================
!     gives the value of L from spetroscopic symbol "a"
!--------------------------------------------------------------------
      Implicit none  
      Character, Intent(in) :: a
      Character(21) :: AS, AB
      Integer :: i

      DATA AS/'spdfghiklmnoqrtuvwxyz'/
      DATA AB/'SPDFGHIKLMNOQRTUVWXYZ'/

      i = INDEX (AS,a)
      if(i.eq.0) i = INDEX (AB,a)
      if(i.eq.0) i = ICHAR(a)-101
      LA = i-1

      End FUNCTION LA

!====================================================================
      Subroutine EL4_NLK(EL,n,l,k)
!====================================================================
!     decodes the A4 specroscopic notation for orbital (n,l,k)
!
!     It is allowed following notations: 1s, 2s3, 2p30, 20p3,
!     1sh, 20sh, kp, kp1, kp11, ns, ns3, ns33, ... .
!
!     Call:  LA, ICHAR
!--------------------------------------------------------------------
      Implicit none
      Character(4), intent(in) :: EL
      Integer, intent(out) :: n,l,k    
      Integer :: i,j,ic, k1,k2  
      Integer, external :: LA

      Character(61) :: ASET = &
       '123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ'

      n=0; l=-1; k=0; i=1; j=1
    1 if(EL(i:i).eq.' ') then
       i=i+1
       if(i.le.4) go to 1
      end if
      if(i.eq.4.and.j.eq.1) j=2
      if(i.gt.4) Return

      ic=ICHAR(EL(i:i))-ICHAR('1')+1

      if(j.eq.1) then                      !  n -> ?
       if(n.eq.0.and.ic.le.9) then
        n=ic; j=1
       elseif(n.eq.0.and.ic.gt.9) then
        n=ic; j=2
       elseif(ic.gt.9) then
        l=LA(EL(i:i)); j=3
       else
        n=n*10+ic
       end if

      elseif(j.eq.2) then                   !  l -> ?
       l=LA(EL(i:i)); j=3

      elseif(j.eq.3) then                   !  k -> ?
       if(i.eq.3) then
        k1 = INDEX(ASET,EL(i:i)); k2 = INDEX(ASET,EL(i+1:i+1))
        if(k2.gt.0) k = k1*61+k2
        if(k2.eq.0) k = k1
       else
        k = INDEX(ASET,EL(i:i))
       end if
       j=4
      end if

      i=i+1
      if(i.le.4.and.j.le.3) go to 1
      if(n.ge.100.or.l.lt.0) then
       write(*,*) 'EL4_nlk is fail to decode: ',EL 
       Stop ' '
      end if

      End Subroutine EL4_NLK

!======================================================================
      Subroutine Read_name(name)
!======================================================================
!     read characer argument
!----------------------------------------------------------------------

      Implicit none

      Character(*) :: name

      Integer :: iarg,i,i1,i2,iname
      Character(80) :: AS
      Integer, External :: IARGC

      iarg = IARGC(); if(iarg.eq.0) Return 
      Do i=1,iarg
       Call GETARG(i,AS)
       if(INDEX(AS,'=').ne.0) Cycle
       name=AS
       Exit
      End do

      End Subroutine Read_name

!======================================================================
      Subroutine Read_aarg(name,avalue)
!======================================================================
!     read characer argument
!----------------------------------------------------------------------
      Implicit None

      Character(*) :: name, avalue
      Integer :: iarg,i,i1,i2,iname
      Character(80) :: AS
      Integer, External :: IARGC

      iarg = IARGC(); if(iarg.eq.0) Return 
      iname=LEN_TRIM(name)
      Do i=1,iarg
       Call GETARG(i,AS)
       i1=INDEX(AS,'=')+1; i2=LEN_TRIM(AS)
       if(AS(1:i1-2).ne.name(1:iname)) Cycle
       read(AS(i1:i2),'(a)') avalue; Exit
      End do

      End Subroutine Read_aarg

!======================================================================
      Subroutine Check_file(AF)
!======================================================================

      Character(*), Intent(in) :: AF
      Logical :: EX

      Inquire(file=AF,exist=EX)
      if(.not.EX) then
       write(*,*) ' can not find file  ',AF;  Stop
      end if

      End Subroutine Check_file


!======================================================================
      Integer Function Iargc()
!======================================================================
!     replace iargc from FORTRAN-77
!----------------------------------------------------------------------
      Implicit none

      iargc = command_argument_count()

      End  Function Iargc

