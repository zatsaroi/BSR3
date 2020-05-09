!=====================================================================
!     PROGRAM       M U L T                             version 4.0
!
!               C O P Y R I G H T -- 2016
!
!     Written by:   Oleg Zatsarinny
!                   email: oleg_zoi@yahoo.com
!======================================================================
!
!     The program evaluates MULTIPOLE operators in LS-coupling
!     including the case of non-orthogonal orbitals
!     
!     The used technique is described in:
!
!     Comp.Phys.Commun. 98 (1996) 235-254
!
!======================================================================
!
!    INPUT ARGUMENTS:
!    
!     AF1  -  c-file for initial state
!     AF2  -  c-file for final state
!     AA   -  type of calculation:  E1,M1,E2,M2,...
! 
!     INPUT FILES:
!
!     AF1  -  c-file for initial state
!     AF2  -  c-file for final state
!     mult_bnk  - data-bank for angular coefficients (optional)
!
!     OUTPUT FILES:
!
!     mult_bnk     -  new data bank for angular coefficients
!                     for given transition AA
!     mult.log     -  running information
!
!-----------------------------------------------------------------------
!     ONLY ONE TYPE OF TRANSITION IN ONE RUN !!!
!-----------------------------------------------------------------------

      Use mult_par

      Character(3) :: arg

      Call inf_mult

! ... check the arguments in command line:

      iarg = IARGC()

      if(iarg.lt.2) then
       write(*,*) 'Should be at least two arguments:'
       write(*,*)
       write(*,*) 'mult bsr klsp1=... klsp2=...             or'
       write(*,*) 'mult dc klsp1=... klsp2=...              or'
       write(*,*) 'mult name1 name2 [Ek|Mk|..] AF_bnk'
       Stop
      end if

! ... open log-file:

      Open(pri,file=AF_pri)
      write(*,'(/a)') &
      'MULT - GENERATION OF DATA BANK FOR MULTIPOLE MATRIX ELEMENTS: '

      Call GETARG(1,arg)

      if(arg.eq.'bsr') then;     Call Mult_bsr
      elseif(arg.eq.'dc') then;  Call Mult_dc
      elseif(arg.eq.'ccc') then; Call Mult_ccc
      else;                      Call Mult_cc
      end if

      End ! Program MULT


!=====================================================================
      Subroutine  MULT_cc
!======================================================================
!  INPUT ARGUMENTS:
!    
!     AF1  -  c-file for initial state
!     AF2  -  c-file for final state
!     AA   -  type of calculation, E1,E2,M1,M2...
! 
!     if aggument is absent, it has been asked from terminal
!
!  INPUT FILES:
!
!     AF1  -  c-file for initial state
!     AF2  -  c-file for final state
!     mult_bnk  - data-bank for angular coefficients (optional)
!     
!  OUTPUT FILES:
!
!     mult_bnk  -  new data bank for angular coefficients
!     mult.log  -  running information
!
!-----------------------------------------------------------------------
      Use mult_par  

      Character(2) :: AA

      iarg = IARGC()

      if(iarg.lt.2) then
       write(*,*) 'Should be at least two file-names in cc-mode:'
       write(*,*)
       Stop 'mult name1 name2 [E1|M1|..] AF_bnk'
      end if

      Call GETARG(1,AF1)
      Call GETARG(2,AF2)
      
! ... define the type of calculations:

      AA = 'E1'
      if(iarg.ge.3) Call GETARG(3,AA)
      read(AA,'(a1,i1)') ktype,kpol
      write(*,'(/a,a1,i1)') 'transition -> ',ktype,kpol
      if(ktype.eq.'M'.and.kpol.eq.0) Stop ' kpol=0 for M-type ? '
      
! ... mult_bnk - if indicated

      if(iarg.ge.4)  Call GETARG(4,AF_b)

      Call Mult_calc

      End Subroutine MULT_cc


!=====================================================================
      Subroutine  MULT_ccc
!======================================================================
!
!     ccc mode:  calculations between set of c-files                        
!     --------                                                              
!     Call as:   mult  ccc  list_c=... AA=..  AF_bnk=..                     
!                                                                        
!     list_c   - file  with  list of involved c-files  (one for line)       
!
!-----------------------------------------------------------------------
      USE mult_par 

      Character(80) :: list_c,AF
      Character(80), Allocatable :: AF_list(:)
      Character(2) :: AA = 'xx'

      Call Read_aarg('list_c',list_c)
      Call Check_file(list_c)
      nul = 20; open(nul,file=list_c)

      nfile=0
      Do 
       read(nul,'(a)',end=1) AF
       if(AF(1:1).eq.'*') Exit
       nfile=nfile+1
      End do
    1 if(nfile.eq.0) Stop ' nfile = '

      Allocate(AF_list(nfile))
      rewind(nul)     
      Do i=1,nfile; read(nul,'(a)') AF_list(i); End do

! ... define the type of calculations:

      Call Read_aarg('AA',AA)
      if(AA.ne.'xx') read(AA,'(a1,i1)') ktype,kpol
      write(*,'(/a,a1,i1)') 'transition -> ',ktype,kpol
      if(ktype.eq.'M'.and.kpol.eq.0) Stop ' kpol=0 for M-type ? '
      Call Read_aarg('AF_b',AF_b)

      Do i=1,nfile-1; Do j=i+1,nfile

       AF1 = AF_list(i)
       Call Check_file(AF1)
       open(in1,file=AF1)
       Call Def_term(in1,ILT1,IST1,IPT1)

       AF2 = AF_list(j)
       Call Check_file(AF2)
       open(in2,file=AF2)
       Call Def_term(in2,ILT2,IST2,IPT2)

       if(IST1.ne.IST2) Cycle

       if(IPT1.eq.IPT2) then
        if(ktype.eq.'E'.and.mod(kpol,2).ne.0) Cycle
        if(ktype.eq.'M'.and.mod(kpol,2).eq.0) Cycle
       else
        if(ktype.eq.'E'.and.mod(kpol,2).ne.1) Cycle
        if(ktype.eq.'M'.and.mod(kpol,2).eq.1) Cycle
       end if

       k=kpol; if(IST1.eq.0) k=k+k
       if(ITRA (ILT1,ILT2,k).eq.0) Cycle

       Call Mult_calc

      End do; End do

      End Subroutine MULT_ccc


!=====================================================================
      Subroutine  MULT_bsr
!======================================================================
!
!     bsr mode:  calculations between set of BSR cfg.nnn files              
!     --------                                                              
!     Call as:   mult  bsr  klsp1=.. klsp2=..  AA=..  AF_bnk=..             
!                                                                      
!     EXAMPLE: FOR CALLING:                                                 
!
!     mult bsr klsp1=1 klsp2=5  AF_b=mult_bnk_E1                               
!     
!-----------------------------------------------------------------------
      USE mult_par  

      Character(3) :: ALSP1,ALSP2, tpar
      Character(2) :: AA

      Integer, allocatable :: lpar(:),ispar(:),ipar(:),jpar(:)

! ... define the type of calculations:

      Call Read_aarg('AA',AA)
      if(AA.ne.'xx') read(AA,'(a1,i1)') ktype,kpol
      write(*,'(/a,a1,i1)') 'transition -> ',ktype,kpol
      if(ktype.eq.'M'.and.kpol.eq.0) Stop ' kpol=0 for M-type ? '
      Call Read_aarg('AF_b',AF_b)

! ... define partial waves:

      klsp1 = 0;  Call Read_iarg('klsp1' ,klsp1 )
      klsp2 = 0;  Call Read_iarg('klsp2' ,klsp2 )

      Call Check_file('target')
      nut = 21; Open(nut,file='target')
      Call Read_ipar(nut,'nlsp',nlsp); read(nut,*)
      if(nlsp.le.0) Stop 'R_channel: nlsp <= 0 '

      Allocate(lpar(nlsp),ispar(nlsp),ipar(nlsp),jpar(nlsp))

      Do i = 1,nlsp
       read(nut,*) tpar,lpar(i),ispar(i),ipar(i) 
       jpar(i) = 0; if(ispar(i).eq.0) jpar(i) = lpar(i) + 1
      End do


      Do ilsp = 1,nlsp-1;    if(klsp1.ne.0.and.ilsp.ne.klsp1) Cycle 
      Do jlsp = ilsp+1,nlsp; if(klsp2.ne.0.and.jlsp.ne.klsp2) Cycle

! ... check the dipole transitions rules: 

       if(ipar(ilsp).eq.ipar(jlsp)) then
        if(ktype.eq.'E'.and.mod(kpol,2).ne.0) Cycle
        if(ktype.eq.'M'.and.mod(kpol,2).eq.0) Cycle
       else
        if(ktype.eq.'E'.and.mod(kpol,2).ne.1) Cycle
        if(ktype.eq.'M'.and.mod(kpol,2).eq.1) Cycle
       end if

       if(ispar(ilsp).ne.ispar(jlsp)) Cycle

       if(ispar(ilsp).ne.0) then
        i = ITRA(lpar(ilsp),kpol,lpar(jlsp))
       else
        i = ITRI(jpar(ilsp),kpol+kpol+1,jpar(jlsp))
       end if
       if(i.eq.0) Cycle

       write(ALSP1,'(i3.3)') ilsp; AF1='cfg.'//ALSP1
       write(ALSP2,'(i3.3)') jlsp; AF2='cfg.'//ALSP2

!       AF_b = AF_b(1:ii)//'.'//ALSP1//'_'//ALSP2

       Call Mult_calc

      End do; End do   ! over partial waves

      End Subroutine MULT_bsr


!=====================================================================
      Subroutine  MULT_dc
!======================================================================
!
!  INPUT ARGUMENTS:
!    
!     klsp1  -  initial partial wave
!     klsp2  -  final partial wave
! 
!     if klsp = 0 -> all partial waves
!
!  INPUT FILES:
!
!     target_dc - partial wave information
!     AF1  -  c-file for initial state, cfg.nnn
!     AF2  -  c-file for final state, cfg.mmm
!     mult_bnk.nnn_mmm  - data-bank for angular coefficients (optional)
!     
!  OUTPUT FILES:
!
!     mult_bnk.nnn_mmm  -  new(revised) data bank for angular coefficients
!     mult.log  -  running information
!
!-----------------------------------------------------------------------

      USE mult_par  

      Integer :: nut=21; Character(80) :: AF_t = 'target_dc'

      Character(3) :: ALSP1,ALSP2, tpar

      Integer(4), Allocatable :: lpar(:),ispar(:),ipar(:)

      Call Read_iarg('klsp1' ,klsp1 )
      Call Read_iarg('klsp2' ,klsp2 )

! ... define partial waves:

      Call Check_file(AF_t)
      Open(nut,file=AF_t)
      Call Read_ipar(nut,'nlsp',nlsp); read(nut,*)
      if(nlsp.le.0) Stop 'target_dc: nlsp <= 0 '

      Allocate(lpar(nlsp),ispar(nlsp),ipar(nlsp))

      nlsp = 0
      Call Read_ipar(nut,'nlsp',nlsp); read(nut,*)
      Do i = 1,nlsp
       read(nut,*) tpar,lpar(i),ispar(i),ipar(i) 
      End do

      ii = LEN_TRIM(AF_b)

      Do ilsp = 1,nlsp-1;    if(klsp1.ne.0.and.ilsp.ne.klsp1) Cycle 
      Do jlsp = ilsp+1,nlsp; if(klsp2.ne.0.and.jlsp.ne.klsp2) Cycle

! ... check the dipole transitions rules: 

       if(ipar(ilsp).eq.ipar(jlsp)) Cycle
       if(ispar(ilsp).ne.ispar(jlsp)) Cycle
       i = ITRA(lpar(ilsp),kpol,lpar(jlsp))
       if(i.eq.0) Cycle

       write(ALSP1,'(i3.3)') ilsp; AF1='cfg.'//ALSP1
       write(ALSP2,'(i3.3)') jlsp; AF2='cfg.'//ALSP2 

       AF_b = AF_b(1:ii)//'.'//ALSP1//'_'//ALSP2

       Call Mult_calc

      End do; End do   ! over partial waves

      End Subroutine MULT_dc


!=====================================================================
      Subroutine  MULT_calc
!======================================================================
!     calculations for given mult_bnk and c-files:
!-----------------------------------------------------------------------
      Use mult_par  
      Use conf_LS,  only: ne
      Use det_list, only: ndet,ldet,jdet
      Use def_list, only: ndef,ldef,jdef

      Implicit none 
      Character AS*80, kt*1
      Integer :: nc, k
      Real(8) :: t1,t2
      Integer, external :: Icheck_file

!-----------------------------------------------------------------------
! ... check the existing data-bank:

      new = 1; if(Icheck_file(AF_b).eq.1) new=0
      if(new.eq.0) then
       Open(nub,file=AF_b,form='UNFORMATTED')
       rewind(nub)
       read(nub) kt,k
       if(kt.ne.ktype.or.k.ne.kpol) new=1
       if(new.eq.1) close(nub) 
      end if
      if(new.eq.1) &
       write(*,'(/a,a)')  'new calculation for ',AF_b
      if(new.eq.0) &
       write(*,'(/a,a)')  'continued calculation for ',AF_b

!-----------------------------------------------------------------------
! ... calculations:

      Call CPU_time(t1)

! ... input of configurations:

      Call R_conf;   if(icalc.eq.0) Return

! ... extract overlap factors:

      Call Read_dets(nub,new)

! ... define possible mls orbitals:
      
      Call Alloc_spin_orbitals(ne)

! ... prepare det.expantions for input configutarions:

      Open(nua,form='UNFORMATTED')
      Open(nud,form='UNFORMATTED')
 
      Call Pre_det_exp 

! ... calculations for new angular symmetries:

      Open(nui,form='UNFORMATTED')

      Call Conf_loop 

!-----------------------------------------------------------------------
! ... record results:

      AF_r = trim(AF_b)//'_res'
      Open(nur,file=AF_r,form='UNFORMATTED')
      rewind(nur)
      write (nur) ktype, kpol
      Call Write_symc_LS(nur)
      Call Write_symt_LS(nur)
      Call Record_done_LS(nur)
      Call Write_det(nur)
      Call Write_def(nur)
      nc = 0
      if(new.eq.0) Call RW(nub,nur,nc)
      rewind(nui); Call RW(nui,nur,nc)
      close(nur); close(nub)

! ... move new results to data bank (inr_res -> int_bnk): 
 
      write(AS,'(a,1x,a,1x,a)') move,trim(AF_r),trim(AF_b)  

      Call System(trim(AS))

      Close(in1); Close(in2); Close(nub)
      Close(nui,status='DELETE')
      Close(nur,status='DELETE') 
      Close(nua,status='DELETE') 
      Close(nud,status='DELETE') 

! ... time for one partial wave:

      write(*,'(/a/)') &
          'Results for new angular symmetry calculations:'

      if(ndet.gt.0) jdet=ldet/ndet+1 
      write(*,'(a,2i10)') &
          'number of overlap determinants =', ndet,jdet
      if(ndef.gt.0) jdef=ldef/ndef+1 
      write(*,'(a,2i10)') &
          'number of overlap factors      =', ndef,jdef 
      write(*,'(a,i10)') &
          'number of angular integrals    =', nc

      Call CPU_time(t2)
      write(*,'(a,F12.2,a)') 'time:',(t2-t1)/60,' min'

      End Subroutine MULT_calc


!======================================================================
      Subroutine Read_dets(nub,new)
!======================================================================
      Implicit none

      Integer :: nub, new

      if(new.eq.1) then 
       Call Alloc_det(-1)
       Call Alloc_def(-1)
      else
       Call Load_det(nub)
       Call Load_def(nub)
      end if

      End Subroutine Read_dets

