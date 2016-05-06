!======================================================================
!     PROGRAM       B S R _ P R E P 3                     
!
!               C O P Y R I G H T -- 2014
!
!     Written by:   Oleg Zatsarinny
!                   email: oleg_zoi@yahoo.com
!
!======================================================================
!
!  1. Target states are ordered according to their energies.
!  2. The orthogonality conditions for one-electron orbitals 
!     are analized and new c-files for target states (targ_nnn.c)
!     and pertubers (pert_nnn.c and pert_nnn.bsw) 
!     are generated with consistent set-indexes.  
!     File target.bsw contains all target orbitals.
!  3. The list of spectroscopic orbitals for target states is created
!     and recorded in target_orb file.
!  4  The orthogonal set of "substitution" orbitals is readed from
!     target_sub.bsw (or created if target_sub.bsw is absent).
!
!  In order to define the same orbitals, the following criteria are 
!  used:  
!              | <p1|p2> - 1 |               <  eps_ovl
!              | <p1|r|p1> - <p2|r|p2> |     <  eps_ovl
!              | <p1|1/r|p1> - <p2|1/r|p2> | <  eps_ovl
!
!  where the parameter eps_ovl can be read from bsr_par, 
!  or as argument: eps_ovl=value
!
!======================================================================
!
!  INPUT FILES:
!
!     bsr_par   -  parameters of calculation if any
!     target    -  list of target states and partial waves 
!     c- and bsw-files for target states and perturbers if any
!     knot.dat  -  B-spline parameters
!     target_sub.bsw  - additional orbitals for orthogonality, if any 
!
!  OUTPUT FILES:
!
!     target        -  modified list of target states
!     target.bsw    -  all target one-elelectron bsw-functions
!     target_orb    -  list of physical target orbitals 
!     targ_nnn.c    -  c-file for target nnn
!     pert_nnn.c    -  c-file for perturber nnn, if any
!     pert_nnn.bsw  -  w-file for perturber nnn, if any
!     bsr_prep.log  -  running information
!
!  ARGUMENTS:   
!
!     eps_ovl  [1.d-6] - tolerance for overlaps
!     eps_phys [0.5  ] - minimum occupation number for physical orbital 
!     eps_sub  [0.5  ] - tolerance for substitution orbitals
!
!======================================================================
     
      Use spline_param; Use spline_grid; Use spline_atomic               
      Use spline_orbitals, p => pbs 
      Use conf_LS; Use orb_LS;  Use target                         
      Use channels, npert1 => npert, ipert1 => ipert, ippert1 => ippert       
      Use channel, only: npert, ipert, ippert

      Implicit real(8) (A-H,O-Z)

      Character(200) :: title, AS
      Character(80), allocatable :: AFK(:)
      Integer, allocatable :: klsp(:) 

! ... files:

      Integer, parameter :: ma=20
      Character(ma) :: AF

      Integer :: nut=1;   Character(ma) :: AF_targ = 'target'
      Integer :: nup=2;   Character(ma) :: AF_par  = 'bsr_par'
      Integer :: pri=66;  Character(ma) :: AF_log  = 'bsr_prep.log'

      Integer :: nuc=11;  Character(ma) :: AFC
      Integer :: nuw=12;  Character(ma) :: AFW
      Integer :: muc=13;  Character(ma) :: BFC
      Integer :: muw=14;  Character(ma) :: BFW

      Integer :: nua=15;  Character(ma) :: AF_sub = 'target_sub.bsw'
      Integer :: nuo=16;  Character(ma) :: AF_orb = 'target_orb'
 
! ... default value for parameter: 

      Real(8) :: eps_ovl  = 1.d-6
      Real(8) :: eps_phys = 0.25d0
      Real(8) :: eps_sub  = 0.5d0

      Call bsr_prep_inf

!----------------------------------------------------------------------
! ... create bsr_prep.log file and check the bsr_par file:

      open(pri,file=AF_log)
      if(Icheck_file(AF_par).ne.0) then
       open(nup,file=AF_par)
      else
       open(nup,file=AF_par)
       write(nup,'(a/)') 'BSR_PREP parameters:'
       write(nup,'(a,1PE7.0,T22,a)') & 
       'eps_ovl  =',eps_ovl, '!  tolerance for overlaps'
       write(nup,'(a,f7.3,T22,a)') & 
       'eps_phys =',eps_phys,'!  occupation number for physical orbital'
       write(nup,'(a,f7.3,T22,a)') & 
       'eps_sub  =',eps_sub, '!  tolerance for substitution orbitals'
       write(*,*) 'Check bsr_par, just in case ...'
      end if

! ... read parameter if any:

      Call Read_rpar(nup,'eps_ovl' ,eps_ovl )
      Call Read_rarg(    'eps_ovl' ,eps_ovl )
      Call Read_rpar(nup,'eps_phys',eps_phys)
      Call Read_rarg(    'eps_phys',eps_phys)
      Call Read_rpar(nup,'eps_sub' ,eps_sub )
      Call Read_rarg(    'eps_sub' ,eps_sub )

      write(pri,'(a/)') 'BSR_PREP:'

      write(pri,'(a,d15.5,a/)') 'eps_ovl  = ',eps_ovl, &
       ' - tolerance for orthogonality between orbitals'
      write(pri,'(a,d15.5,a/)') 'eps_sub  = ',eps_sub, &
       ' - tolerance for substitution orbitals'
      write(pri,'(a,d15.5,a/)') 'eps_phys = ',eps_phys,& 
       ' - tolerance for physical orbitals'

      write(pri,'(/72(''-'')/)')

! ... sets up grid points and initializes the values of the spline: 
    
      CALL define_grid(z);  CALL define_spline

! ... read file "target":

      if(Icheck_file(AF_targ).gt.0) then
       Open(nut,file=AF_targ); rewind(nut)
      else
       Call Print_target
       Stop 'Check or prepare target file'
      end if

      rewind(nut);  read(nut,'(a)') title
      Call Read_apar(nut,'coupling',coupling)
      nz = 0;    Call Read_ipar(nut,'nz',nz)
      if(nz.le.0)   Stop 'bsr_prep: nz = 0  ? '
      nelc = 0;  Call Read_ipar(nut,'nelc',nelc)
      if(nelc.le.0) Stop 'bsr_prep: nelc = 0  ? '

! ... define number of target states:
      
      ntarg=0; Call Read_ipar(nut,'ntarg',ntarg)
      if(ntarg.le.0) Stop 'bsr_prep: ntarg = 0  ? '
      m=ntarg; Call Allocate_target(m)
      read(nut,*)
      Do i=1,ntarg
       read(nut,*) AFT(i)
       Call Check_BSR_name(AFT(i))
      End do
 
! ... define the target energies and check cores:

      Do it=1,ntarg
       AFC = trim(AFT(it))//'.c'
       Call Check_file(AFC)
       Open(nuc,file=AFC)
       read(nuc,'(15x,F16.8)') Etarg(it)
       Call R_CLOSED(nuc)
       if(it.eq.1) then
        core = CLOSED; ncore=nclosd
        write(pri,'(a,a)') 'CORE:  ',core
        AFW = trim(AFT(it))//'.bsw'
        Call Check_bsw_file(AFW,pri)
        Open(nuw,file=AFW,form='UNFORMATTED')
        Call R_bwfn(nuw)
       else
        if(CLOSED.ne.core)  then
         write(pri,'(/a,a,a)') 'target ',AFC,' has different core'
         Stop ' check core ! ' 
        end if
       end if
       Close(nuc)
      End do
      write(pri,'(/72(''-'')/)')

! ... sorting target states according to energy:

      Do i = 1,ntarg-1
       Do j = i+1,ntarg
        if(Etarg(j).ge.Etarg(i)) Cycle
        E=Etarg(i); Etarg(i)=Etarg(j); Etarg(j)=E
        AF=AFT(i); AFT(i)=AFT(j); AFT(j)=AF
       End do
      End do
      nct = 0; nwt = nclosd;  nbf = nclosd

! ... check substitution orbitals if any

      if(Icheck_file(AF_sub).ne.0) then

       write(pri,'(a,4x,a20/)') ' substitution orbitals from',AF_sub
       ncfg = 0; lcfg = 0;  nwf = ncore

       AFC = 'no';    AFW = AF_sub;   BFC = 'no';   BFW = 'no'  

       Call Check_file(AFW)
       Open(nuw,file=AFW,form='UNFORMATTED')
       Call Read_bsw_orb_LS(nuw) ! read orbitals into orb_LS module

       kshift=0; IEF=0;  CALL SUB_check_orb

      end if
 
      nwt = nbf;  iech = 0; iech(1:nbf) = 1  ! sub.orb. pointer
      write(pri,'(/72(''-'')/)')
     
!=====================================================================
! ... proceed target files:

      open(nuo,file=AF_orb)
      Do it=1,ntarg

       write(pri,'(a,i5,4x,a20/)') ' target ',it,AFT(it)

       AFC=trim(AFT(it))//'.c'
       Call Check_file(AFC)
       Open(nuc,file=AFC)
       nwf=ncore; ncfg=0; lcfg=0;  Call R_conf_LS(nuc,0)
       
       AFW=trim(AFT(it))//'.bsw'
       
       kshift=0; IEF=0; CALL SUB_check_orb

       if(ntarg.lt.1000) then
         write(BFC,'(a,i3.3,a)') 'targ_',it
       elseif(ntarg.lt.10000) then
         write(BFC,'(a,i4.4,a)') 'targ_',it
       elseif(ntarg.lt.100000) then 
         write(BFC,'(a,i5.5,a)') 'targ_',it
       else
         Stop 'Stop in bsr_prep:  ntarg > 100000 '
       end if
       BFT(it) = BFC;  BFC = trim(BFT(it))//'.c'

       Open(muc,file=BFC); rewind(nuc)
       read(nuc,'(a)') AS;  write(muc,'(a)') AS
       read(nuc,'(a)') AS;  write(muc,'(a)') AS
       Do ic=1,ncfg; if(WC(ic).eq.0.d0) WC(ic)=1.d-8;  End do
       Do ic=1,ncfg; Call Pri_conf (muc,ic,WC(ic));  End do
       write(muc,'(a)') '*'
       write(nuo,'(a,3x,i3.3,6x,a)') 'target',it,AFT(it)

       Call SUB_phys_orb

       Call Def_term_BSR(nuc,ltarg(it),istarg(it),iptarg(it))
       nctarg(it) = ncfg
       nwtarg(it) = nbf-nwt
       nct=nct+ncfg; nwt=nbf
       write(pri,'(/72(''-'')/)')

      End do ! over it - target states

!----------------------------------------------------------------------
! ... create the target.bsw file:

      Open(muw,file='target.bsw',form='UNFORMATTED')
      Do i = 1,nbf
       write(muw) ebs(i),z,h,hmax,rmax,ks,ns,mbs(i)
       write(muw) (p(j,i),j=1,mbs(i))
      End do
      Close(muw)

! ... final list of target one-electron orbitals: 

      write(pri,'(/a,i5/)') ' target orbitals, nbf = ',nbf
      write(pri,'(14a5)') (EBS(j),j=1,nbf)
      write(pri,'(/72(''-'')/)')

! ... substitution orbitals information:
      
      nwf_sub = 0
      Do i=1,nbf; if(iech(i).gt.0) nwf_sub=nwf_sub+1; End do
      write(pri,'(/a,i5/)') &
       ' substitution orbitals for target states: nsub = ',nwf_sub
      k = 1
      Do i=1,nbf
       if(iech(i).ne.1) Cycle
       write(AS(k:),'(a5)') ebs(i)
       k=k+5
       if(k.lt.70) Cycle
       write(pri,'(a)') AS
       k=1
      End do
      if(k.gt.1) write(pri,'(a)') AS
 
      k = 0
      Do i=1,nbf; if(iech(i).ne.1) Cycle
       Do j=1,i; if(iech(j).ne.1) Cycle
        if(lbs(i).ne.lbs(j)) Cycle
        S = abs(QUADR(i,j,0)); if(i.eq.j) S=S-1.d0
        if(S.lt.eps_ovl) Cycle
        write(pri,'(2a5,f20.8)') ebs(i),ebs(j),S
        k=1
       End do
      End do
      if(k.eq.0) write(pri,'(/a/)') &
         ' all substitution orbitals are orthogonal'
      write(pri,'(/72(''-'')/)')

!=====================================================================
! ... define number of partial waves:

      kpert=0; Call Read_ipar(nut,'kpert',kpert)
      nlsp=0;  Call Read_ipar(nut,'nlsp',nlsp)
      if(nlsp.eq.0) kpert=0

      if(nlsp.gt.0) then
       Allocate(ispar(nlsp),lpar(nlsp),ipar(nlsp),Tpar(nlsp), &
                AFP(nlsp),BFP(nlsp),ncp(nlsp),nwp(nlsp))
       read(nut,*)
       Do i = 1,nlsp
        read(nut,'(a)') AS;  ii = LEN_TRIM(AS)
        read(AS,*) Tpar(i),lpar(i),ispar(i),ipar(i)
        AFP(i) = 'no'; if(ii.gt.18) read(AS(19:),*) AFP(i)
        Call Check_BSR_name(AFP(i))
        BFP(i) = 'no'
       End do

! ... define number of additional perturbers:

       kpert=0; Call Read_ipar(nut,'kpert',kpert)

       if(kpert.gt.0) then
        Allocate(klsp(kpert),AFK(kpert))
        read(nut,*)
        Do i=1,kpert
         read(nut,*) klsp(i),AFK(i)
         Call Check_BSR_name(AFK(i))
        End do

       end if

      else               ! define range of partial waves
       LT_min = 0
       Call Read_ipar(nup,'LT_min',LT_min)
       Call Read_iarg(    'LT_min',LT_min)
       LT_max = 25
       Call Read_ipar(nup,'LT_max',LT_max)
       Call Read_iarg(    'LT_max',LT_max)
       IS_min = -1
       Call Read_ipar(nup,'IS_min',IS_min)
       Call Read_iarg(    'IS_min',IS_min)
       IS_max = -1
       Call Read_ipar(nup,'IS_max',IS_max)
       Call Read_iarg(    'IS_max',IS_max)
       JJ_min = -1
       Call Read_ipar(nup,'JJ_min',JJ_min)
       Call Read_iarg(    'JJ_min',JJ_min)
       JJ_max = -1
       Call Read_ipar(nup,'JJ_max',JJ_max)
       Call Read_iarg(    'JJ_max',JJ_max)

       write(*,*) 'LT_min =',LT_min
       write(*,*) 'LT_max =',LT_max
       write(*,*) 'IS_min =',IS_min
       write(*,*) 'IS_max =',IS_max
       write(*,*) 'JJ_min =',JJ_min
       write(*,*) 'JJ_max =',JJ_max

       if(JJ_max.ne.-1.and.JJ_min.ne.-1) then

        nlsp = JJ_max-JJ_min + 2
        Allocate(ispar(nlsp),lpar(nlsp),ipar(nlsp), Tpar(nlsp), & 
                 AFP(nlsp), BFP(nlsp),ncp(nlsp),nwp(nlsp))
        i = 0
        Do JJ = JJ_min,JJ_max,2
         Do ip = -1,1,2
          i = i + 1  
          ispar(i)=0; lpar(i)=JJ; ipar(i)=ip; AFP(i) = 'no'      
          write(Tpar(i),'(i3.3)') i
         End do

        End do

       elseif(IS_max.ne.-1.and.IS_min.ne.-1) then
        nlsp = (LT_max-LT_min + 1) * (IS_max-IS_min + 2)        
        write(*,*) 'nlsp =',nlsp

        Allocate(ispar(nlsp),lpar(nlsp),ipar(nlsp), Tpar(nlsp), & 
                 AFP(nlsp), BFP(nlsp),ncp(nlsp),nwp(nlsp))

        i = 0
        Do L = LT_min,LT_max
         Do is = IS_min,IS_max,2
          Do ip = -1,1,2
           i = i + 1  
           ispar(i)=is; lpar(i)=L; ipar(i)=ip; AFP(i) = 'no'      
           write(Tpar(i),'(i3.3)') i
          End do
         End do
        End do
       end if 

       AFP ='no';  BFP = 'no'

      end if    ! over nlsp > 0
!=====================================================================
! ... proceed perturbers if any:

      Do ilsp=1,nlsp
       nbf = nwt
       BFC = 'no'; BFW = 'no'
       ncfg = 0; lcfg = 0; nwf = ncore  

       npert = 0;  Call Allocate_pert(0)

! ... usual perturber (as list of separate configurations) :

       if(LEN_TRIM(AFP(ilsp)).gt.0.and.AFP(ilsp).ne.'no') then

        AFC = trim(AFP(ilsp))//'.c'  
        write(pri,'(/a,i4,4x,a/)') 'perturber', ilsp, trim(AFC)
        Call Check_file(AFC)
        Open(nuc,file=AFC)
        Call Add_conf_LS(nuc,0)
        WC(1:ncfg) = 1.d0
        ncp(ilsp) = ncfg

        Call Allocate_pert(ncfg+ipert); npert=ncfg
 
        Do i=1,npert; ippert(i)=i; End do
        
        AFW = trim(AFP(ilsp))//'.bsw';
        Call Check_file(AFW)

        kshift=0; IEF=0; CALL  SUB_check_orb

       end if

! ... additional perturbers:

       Do iip=1,kpert
        if(klsp(iip).ne.ilsp) Cycle
        AFC = trim(AFK(iip))//'.c'
        Call Check_file(AFC)
        write(pri,'(/a,i4,4x,a/)') 'perturber', ilsp, trim(AFK(iip))
        Open(nuc,file=AFC)

        kshift = maxval(KEF(1:nwf))

        Call Add_conf_LS(nuc,kshift)
        if(npert+1.gt.mpert) Call Allocate_pert(mpert+ipert)
        npert = npert + 1        
        ippert(npert)=ncfg
        ncp(ilsp) = ncfg
        AFW = trim(AFK(iip))//'.bsw'
        Call Check_file(AFW)

        CALL SUB_check_orb
   
       End do 

! ... collective perturber configurations:

       if(ncfg.gt.0) then

        write(BFP(ilsp),'(a,i3.3)') 'pert_',ilsp
        BFC = trim(BFP(ilsp))//'.c'
        Open(muc,file=BFC); rewind(nuc)
        read(nuc,'(a)') AS;  write(muc,'(a)') AS
        read(nuc,'(a)') AS;  write(muc,'(a)') AS
        Do ic=1,ncfg; Call Pri_conf (muc,ic,WC(ic));  End do
        write(muc,'(a)') '*'

        write(muc,'(a,i5)') 'npert =',npert
        write(muc,'(20i5)') (ippert(i),i=1,npert)

        write(nuo,'(a,3x,i3.3)') 'pertub',ilsp
        Call SUB_phys_pert

       end if

! ... if we have additional perturber orbitals:

       nwp(ilsp) = nbf - nwt
       if(nwp(ilsp).gt.0) then 
        BFW = trim(BFP(ilsp))//'.bsw'
        Open(muw,file=BFW,form='UNFORMATTED',status='UNKNOWN')
        Do i = nwt+1,nbf
         write(muw) ebs(i),z,h,hmax,rmax,ks,ns,mbs(i)
         write(muw) (p(j,i),j=1,mbs(i))
        End do
        Close(muw)
       end if

       ncp(ilsp) = ncfg
 
      End do  ! over partial waves

      write(pri,'(/72(''-'')/)')

!=====================================================================
! ... update the target file:

      rewind(nut)
      write(nut,'(a)') TITLE 
      write(nut,'(72(''-''))')
      write(nut,'(a,a2,4x,a)') &
                'coupling = ',coupling, '!   coupling scheme'
      write(nut,'(a,i5,T18,a)') &
                'nz     =',nz,   '!   nuclear charge' 
      write(nut,'(a,i5,T18,a)') &
                'nelc   =',nelc, '!   number of electrons'
      write(nut,'(72(''-''))')
      write(nut,'(a,i5,T18,a)') &
                'ntarg  =',ntarg,'!   number of target states'
      write(nut,'(72(''-''))')
      Do i=1,ntarg
       write(nut,'(a20,1x,a10,1x,3i4,f18.8,2i5)') AFT(i),BFT(i), &
        ltarg(i),istarg(i),iptarg(i),etarg(i),nctarg(i),nwtarg(i)
      End do
      write(nut,'(72(''-''))')
      write(nut,'(a,i5,T18,a)') 'nct    =',nct, &
       '!   total number of target configurations' 
      write(nut,'(a,i5,T18,a)') 'nwt    =',nwt, &
       '!   total number of target orbitals' 
      write(nut,'(a,i5,T18,a)') 'nsub   =',nwf_sub, &
       '!   number of substitution orbitals' 

      write(nut,'(72(''-''))')
      write(nut,'(a,i5,T18,a)') &
           'nlsp   =',nlsp, '!   number of partial waves' 
      write(nut,'(72(''-''))')
      if(nlsp.le.0) go to 100

      Do i = 1,nlsp
       if(AFP(i).ne.'no'.or.BFP(i).ne.'no') then
         write(nut,'(a3,3i5,3x,a20,1x,a10,2i5)') &
         Tpar(i),lpar(i),ispar(i),ipar(i),AFP(i),BFP(i),ncp(i),nwp(i)
       else 
         write(nut,'(a3,3i5)') Tpar(i),lpar(i),ispar(i),ipar(i)
       end if
      End do
      write(nut,'(72(''-''))')

      if(kpert.eq.0) go to 100
      write(nut,'(a,i5,T18,a)') &
           'kpert  =',kpert,   '!   number of additional perturbers' 
      write(nut,'(72(''-''))')
      Do i = 1,kpert
       write(nut,'(i3,3x,a)') klsp(i),trim(AFK(i))      
      End do
      write(nut,'(72(''-''))')

  100 Continue

      CONTAINS
     

!======================================================================
      Subroutine SUB_check_orb
!======================================================================
!     this routine analizes one set of orbitals in file AFW
!----------------------------------------------------------------------
      Implicit real(8) (A-H,O-Z)
      Character(4) :: elw
      Character(4), external :: ELF4
      Integer, allocatable :: IP_occ(:)
!----------------------------------------------------------------------
! ... read radial functions:

      Open(nuw,file=AFW,form='UNFORMATTED')
      rewind(nuw)

    1 m=nbf+1; if(m.ge.mbf) CALL Allocate_bsorb(mbf+jbf)

      read(nuw,end=2) elw,zw,hw,hmw,rmw,ksw,nsw,mw

      read(nuw) p(1:mw,m); if(mw.lt.ns) p(mw+1:ns,m)=0.d0

      ! check consistence of B-spline basis:
      if(zw.ne.z)     Stop ' Read_bsw:  z <> zw'
      if(hw.ne.h)     Stop ' Read_bsw:  h <> hw'
      if(abs(hmw-hmax).gt.1.D-12) then
                      write(*,*) elw,hmw,hmax, hmw-hmax
                      Stop ' Read_bsw:  hmw <> hmax'
      end if
      if(rmw.ne.rmax) Stop ' Read_bsw:  rmw <> rmax'
      if(ksw.ne.ks)   Stop ' Read_bsw:  ksw <> ks'
      if(nsw.ne.ns)   Stop ' Read_bsw:  nsw <> ns'

! ... find orbital in orbital list:

       Call EL4_nlk(elw,n,l,k);  k=k+kshift
       ii=Ifind_nlk(n,l,k,0)
       if(ii.gt.0) then                              
        mbs(m)=mw; nbs(m)=n; lbs(m)=l; kbs(m)=k; ebs(m)=elw
       else
        write(pri,'(a4,12x,a)') elw,' excessive orbital'
        go to 1
       end if

! ...  define overlaps with existing orbitals:

       OBS(:,m) = 0.d0
       Do i = 1,m
        if(lbs(i).ne.lbs(m)) Cycle
        OBS(i,m)=QUADR(i,m,0)
        OBS(m,i)=OBS(i,m)
       End do

!----------------------------------------------------------------------
! ...  compare with the existing orbitals:
       
       SM1 = QUADR(m,m, 1); SM2 = QUADR(m,m, 2)
       Do i = 1,nbf;     if(abs(OBS(i,m)).lt.eps_ovl) Cycle

! ... check orthogonality to core: 

        if(i.le.nclosd.and.ii.gt.nclosd) then
         write(pri,'(a,i3,a,a,a,f12.8)') ' target ',it, '  orbital ', &
               elw, ' does not orthogonal to core:', OBS(i,m)
         Stop 'problem with orthogonality to core'
        end if

! ... define if it is approximately the same orbital:

        S  = abs(OBS(i,m)-1.d0);      if(S .gt.eps_ovl) Cycle
        S1 = abs(QUADR(i,i, 1)-SM1);  if(S1.gt.eps_ovl) Cycle
        S2 = abs(QUADR(i,i, 2)-SM2);  if(S2.gt.eps_ovl) Cycle

        IEF(ii) = i
        write(pri,'(a,a,a,a)')  elw,' --> ',ebs(i),'    the same'
        go to 1
       
       End do

!---------------------------------------------------------------------
! ...  core orbitals should be the same:   

       if(it.gt.1.and.ii.le.nclosd) then
        write(pri,'(a,a,a)') 'file ',AFW,'  has another core orbital'
        Stop ' another core orbital? '
       end if

! ... assign set index for new orbital: 

       Call Assign_index(m); nbf=m; IEF(ii)=m 

       go to 1    ! go to next orbital
   2  Continue

! ... check if bsw-file contains all radial orbitals:

      Do i=ncore+1,nwf; ii=IEF(i)
       if(ii.eq.0) then
        write(pri,'(a,a,a,a)') &
        AFT(it), '-  orbital ',ELF(i),' was not found in the w-file'
        Stop ' unknown orbitals ! '
       else
        NEF(i)=nbs(ii); LEF(i)=lbs(ii); KEF(i)=kbs(ii); ELF(i)=ebs(ii)
       end if
      End do

      End Subroutine SUB_check_orb


!======================================================================
      Subroutine SUB_phys_orb
!======================================================================
!     define physical orbitals for given target state
!----------------------------------------------------------------------
      Implicit real(8) (A-H,O-Z)
      
      Real(8), allocatable :: S_orb(:)

      if(allocated(S_orb)) Deallocate(S_orb)
      Allocate(S_orb(nbf))

! ... find the occupation numbers for all orbitals:

      S_orb = 0.d0; SS = 0.d0
      Do ic = 1,ncfg; S = WC(ic)*WC(ic); SS = SS + S
       Call Get_cfg_LS(ic)
       ip = ip_state(ic)
       Do i=1,no; ip=ip+1; ii=IEF(IP_orb(ip))   
        S_orb(ii) = S_orb(ii) + iq(i)*S 
       End do
      End do
      S_orb = S_orb / SS

      write(muc,'(/a)') 'Physical orbitals:'
      SS = 0
      Do ii=1,nbf;  if(S_orb(ii).lt.eps_phys) Cycle

! ... choose the substitution orbital with biggest overlap:

       SM = 0.d0; jj = 0 
       Do k=ncore+1,nbf
        if(iech(k).ne.1)      Cycle
        if(lbs(k).ne.lbs(ii)) Cycle
        if(nbs(ii).lt.9.and.nbs(k).lt.nbs(ii)) Cycle    !!!  ???
        s_ovl = abs(OBS(k,ii))
        if(s_ovl.lt.SM)       Cycle
        SM = s_ovl; if(s_ovl.gt.eps_sub) jj = k
       End do

       if(jj.eq.0) then
         Call Add_sub_orbital(ii,jj) 
         SM = OBS(ii,jj)
       end if

       write(muc,'(a4,f8.1,5x,a4,f8.3)') ebs(ii),S_orb(ii),ebs(jj),SM
       write(nuo,'(a4,f8.1,5x,a4,f8.3)') ebs(ii),S_orb(ii),ebs(jj),SM

       SS = SS + S_orb(ii)

      End do          

      write(muc,'(a,20x,f8.3)') '*'     !,nelc-SS
      write(nuo,'(a,20x,f8.3)') '*'     !,nelc-SS

      End Subroutine SUB_phys_orb


!======================================================================
      Subroutine SUB_phys_pert
!======================================================================
!     define physical orbitals in perturber
!----------------------------------------------------------------------
      Implicit real(8) (A-H,O-Z)
      
      Real(8), allocatable :: S_orb(:), SS_orb(:)

      if(allocated(S_orb)) Deallocate(S_orb,SS_orb)
      Allocate(S_orb(nbf),SS_orb(nbf))

      write(muc,'(/a)') 'Physical orbitals:'

      S_orb = 0.d0
      Do jp=1,npert

       SS_orb = 0.d0; SS = 0.d0
       Do ic = ippert(jp-1)+1,ippert(jp)
        S = WC(ic)*WC(ic); SS = SS + S
        Call Get_cfg_LS(ic)
        ip = ip_state(ic)
        Do i=1,no; ip=ip+1; ii=IEF(IP_orb(ip))   
         SS_orb(ii) = SS_orb(ii) + iq(i)*S 
        End do
       End do
!      SS_orb = SS_orb / SS  ! ????

       Do i=1,nbf
        S_orb(i)=max(S_orb(i),SS_orb(i))       

       End do

      End do

! ... choose the substitution orbital with biggest overlap:

      Do ii=1,nbf;  if(S_orb(ii).lt.eps_phys) Cycle

       SM = 0.d0; jj = 0 
       Do k=ncore+1,nbf
        if(iech(k).ne.1)      Cycle
        if(lbs(k).ne.lbs(ii)) Cycle
        s_ovl = abs(OBS(k,ii))
        if(S_ovl.lt.SM) Cycle
        SM = s_ovl; if(s_ovl.gt.eps_sub) jj = k
       End do
       
       if(jj.eq.0) then
        Call Add_sub_orbital(ii,jj) 
        SM = OBS(ii,jj)
       end if

       write(muc,'(a4,f8.1,5x,a4,f8.3)') ebs(ii),S_orb(ii),ebs(jj),SM
       write(nuo,'(a4,f8.1,5x,a4,f8.3)') ebs(ii),S_orb(ii),ebs(jj),SM

      End do          

      write(nuo,'(a)') '*'
      write(muc,'(a)') '*'

      End Subroutine SUB_phys_pert


!======================================================================
      Subroutine Add_sub_orbital(ii,jj)
!======================================================================

      Implicit real(8) (A-H,O-Z)

      m=nbf+1; if(m.ge.mbf) CALL Allocate_bsorb(mbf+jbf)

      nbs(m)=nbs(ii); lbs(m)=lbs(ii); kbs(m)=kbs(ii); ebs(m)=ebs(ii)
      mbs(m)=mbs(ii); p(:,m) = p(:,ii)

       OBS(:,m) = 0.d0
       Do i = 1,m
        if(lbs(i).ne.lbs(m)) Cycle
        OBS(i,m)=QUADR(i,m,0)
        OBS(m,i)=OBS(i,m)
       End do

! ... check if we need orthogonalization:

      jj = ii
      Do i=1,nbf; if(iech(i).ne.1) Cycle
       S = OBS(i,m); if(abs(S).lt.eps_ovl) Cycle           
       p(:,m) = p(:,m) - S * p(:,i);  jj=m
      End do

      if(jj.eq.ii) then; iech(ii)=1; Return; End if

      S = QUADR(m,m,0);  S=sqrt(S);  p(:,m)=p(:,m)/S

       OBS(:,m) = 0.d0
       Do i = 1,m
        if(lbs(i).ne.lbs(m)) Cycle
        OBS(i,m)=QUADR(i,m,0)
        OBS(m,i)=OBS(i,m)
       End do

! ... assign set index for new sub. orbital: 

      Call Assign_index(m); nbf=m; iech(m)=1; jj=m

      End Subroutine Add_sub_orbital


!======================================================================
      Subroutine Assign_index(m)
!======================================================================
!     assign set index for orbital m
!----------------------------------------------------------------------
      Implicit real(8) (A-H,O-Z)
      Character(4) :: ELW
      Character(4), external :: ELF4

       kbs(m)=-1; ELW = EBS(m)
       knew = New_index(lbs(m),ksmax,nbf,lbs,kbs)

       Do k = 1,knew-1
        S=0.d0           
        Do i = 1,nbf
         if(lbs(m).ne.lbs(i).or.k.ne.kbs(i)) Cycle
         S=max(S,abs(OBS(i,m)))
        End do
        if(S.lt.eps_ovl) then; kbs(m)=k; Exit; end if
       End do  

       if(kbs(m).eq.-1) then  ! the orbital belongs to new set  

        kbs(m) = knew
        EBS(m)=ELF4(nbs(m),lbs(m),kbs(m))
        write(pri,'(a,a,a,a,f15.8)') &
              elw,' --> ',EBS(m),'    new orbitals and new set index',S

       else                   ! check the same label for diff.orbitals            

        EBS(m)=ELF4(nbs(m),lbs(m),kbs(m))
        Do i = 1,nbf
         if(EBS(m).ne.EBS(i)) Cycle
         kbs(m) = knew
         EBS(m)=ELF4(nbs(m),lbs(m),kbs(m))
         write(pri,'(a,a,a,a,f15.8)') &
              elw,' --> ',EBS(m),'    new orbitals and new set index',S
        End do
        if(kbs(m).ne.knew) write(pri,'(a,a,a,a,f15.8)') &
              elw,' --> ',EBS(m),'    new orbitals but old set index',S
       
       end if

      End Subroutine Assign_index

     End  ! program BSR_PREP


!======================================================================
      Subroutine Check_bsw_file(AFW,pri)
!======================================================================
! ... check if we have bsw-file, otherwis try to call w_bsw
!----------------------------------------------------------------------
      Character(*) :: AFW
      Character(80) :: BFW, AS
      Integer :: pri

      if(Icheck_file(AFW).eq.0) then
       ilen=LEN_TRIM(AFW); BFW=AFW(1:ilen-3)//'w'
       if(Icheck_file(BFW).eq.0) then
        write(*,*) 'Can not find ',trim(AFW),' or ',trim(BFW) 
        write(pri,'(a,a,a,a)') 'Can not find ',trim(AFW),' or ',trim(BFW) 
        Stop ' '
       end if
       AS = 'w_bsw '//trim(BFW)
       i = SYSTEM(AS)
       write(*,*) 'Calling ',trim(AS),' - check then w_bsw.log' 
       write(pri,'(a,a,a,a)') 'Calling ',trim(AS),' - check then w_bsw.log'  
      end if

      End Subroutine Check_bsw_file


!======================================================================
      Subroutine Check_BSR_name(AF)
!======================================================================
!     rename old-version file-name  as name.c just to name
!----------------------------------------------------------------------

      Character(*) :: AF

      i = LEN_TRIM(AF)

      if(i.gt.2.and.AF(i-1:i).eq.'.c') AF(i-1:i) = '  '

      End Subroutine Check_BSR_name


!=====================================================================
      Subroutine bsr_prep_inf
!=====================================================================
!     print information for bsr_prep program
!---------------------------------------------------------------------
      Implicit none
      Character ::  name = ' '
      
      Call Read_name(name)
      if(name.ne.'?') Return

      write(*,'(a)') &
'                                                                            ', &
'  bsr_prep analizes the input target states and prepare them for BSR:       ', &
'                                                                            ', &
'  INPUT FILES:                                                              ', &
'                                                                            ', &
'     bsr_par   -  parameters of calculation if any                          ', &
'     target    -  list of target states and partial waves                   ', &
'     c- and bsw-files for each target state and perturber if any            ', &
'     knot.dat  -  B-spline parameters                                       ', &
'     target_sub.bsw  - additional orbitals for orthogonality, if any        ', &
'                                                                            ', &
'  OPTIONAL ARGUMENTS:                                                       ', &
'                                                                             ', &
'     eps_ovl  [1.d-6] - tolerance for overlaps                               ', &
'     eps_phys [0.5  ] - minimum occupation number for orbital to be physical ', &
'     eps_sub  [0.5  ] - tolerance for substitution orbitals                  ', &
'                                                                             ',&
'  additional arguments LT_min,LT_max, IS_min,IS_max, or JJ_min,JJ_max       ', &
'  can be used at nlsp=0 for generation of possible partial waves            ', &
'  (with S -> 2S+1  and J -> 2J representation)                              ', &
'                                                                       '   
      Stop 
                                                                             
      End Subroutine bsr_prep_inf                                             
