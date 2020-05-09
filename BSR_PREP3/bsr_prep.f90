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
!  2. The one-electron orbitals are analized and new c-files for target
!     states (targ_nnn.c)  and pertubers (pert_nnn.c and pert_nnn.bsw) 
!     are generated with consistent set-indexes.  
!     File target.bsw contains all target orbitals.
!  3. The list of spectroscopic orbitals for target states is created
!     and recorded in the "target_orb" file.
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
!     bsr_par   -  parameters of calculation (optional)
!     target    -  list of target states and partial waves 
!     c- and bsw-files for target states and perturbers if any
!     knot.dat  -  B-spline parameters
!     target_sub.bsw  - substitution orbitals (optional) 
!
!  OUTPUT FILES:
!
!     target        -  modified list of target states
!     target.bsw    -  all one-elelectron target orbitals 
!     target_orb    -  list of physical target orbitals 
!     targ_nnn.c    -  c-file for target "nnn"
!     pert_nnn.c    -  c-file for perturber "nnn", if any
!     pert_nnn.bsw  -  w-file for perturber "nnn", if any
!     bsr_prep.log  -  running information
!
!  ARGUMENTS:   
!
!     eps_ovl  [1.d-6] - tolerance for overlaps
!     eps_phys [0.5  ] - minimum occupation number for physical orbital 
!     eps_sub  [0.5  ] - tolerance for substitution orbitals
!     ii_sub   [0    ] - if > 0, prevents the generation new substitution
!                        orbitals   
!======================================================================
     
      Use bsr_prep

      Implicit real(8) (A-H,O-Z)

! ... provide information about program:

      Call bsr_prep_inf

! ... open bsr_prep.log file and check the bsr_par file:

      open(pri,file=AF_log)

! ... read parameter if any:

      Call read_arg

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
 
! ... define the target energies:

      Do it=1,ntarg
       AFC = trim(AFT(it))//'.c'
       Call Check_file(AFC)
       Open(nuc,file=AFC)
       read(nuc,'(15x,F16.8)') Etarg(it)
       Close(nuc)
      End do

! ... sorting target states according to energy:

      Do i = 1,ntarg-1
       Do j = i+1,ntarg
        if(Etarg(j).ge.Etarg(i)) Cycle
        E=Etarg(i); Etarg(i)=Etarg(j); Etarg(j)=E
        AF=AFT(i); AFT(i)=AFT(j); AFT(j)=AF
       End do
      End do

! ... check cores:

      AFC = trim(AFT(1))//'.c'
      Open(nuc,file=AFC)
      Call R_CLOSED(nuc)
      core = CLOSED; ncore=nclosd
      write(pri,'(a,a)') 'CORE:  ',core
      Close(nuc)
      AFW = trim(AFT(1))//'.bsw'
      Call Check_bsw_file(AFW,pri)
      Open(nuw,file=AFW,form='UNFORMATTED')
      Call R_bwfn(nuw)
      Close(nuw)

      Do it=2,ntarg
       AFC = trim(AFT(it))//'.c'
       Open(nuc,file=AFC)
       Call R_CLOSED(nuc)
       if(CLOSED.ne.core)  then
         write(pri,'(/a,a,a)') 'target ',AFC,' has different core'
         Stop ' check core ! ' 
       end if
       Close(nuc)
      End do

      write(pri,'(/72(''-'')/)')

      nct = 0; nwt = nclosd;  nbf = nclosd

!-----------------------------------------------------------------------------
! ... check substitution orbitals if any

      if(Icheck_file(AF_sub).ne.0) then

       write(pri,'(a,4x,a20/)') ' substitution orbitals are read from',AF_sub
       ncfg = 0; lcfg = 0;  nwf = ncore

       AFC = 'no';    AFW = AF_sub;   BFC = 'no';   BFW = 'no'  

       Call Check_file(AFW)
       Open(nuw,file=AFW,form='UNFORMATTED')
       Call Read_bsw_orb_LS(nuw) ! read orbitals into orb_LS module

       kshift=0; IEF=0;  CALL SUB_check_orb

      end if
 
      nwt = nbf;  iech = 0; iech(1:nbf) = 1  ! sub.orb. pointer
      write(pri,'(/72(''-'')/)')
     
!============================================================================
! ... proceed target files:

      open(nuo,file=AF_orb)

      Do it=1,ntarg

       write(pri,'(a,i5,4x,a20/)') ' target ',it,AFT(it)

       AFC=trim(AFT(it))//'.c'
       Call Check_file(AFC)
       Open(nuc,file=AFC)
       nwf=ncore; ncfg=0; lcfg=0;  Call R_conf_LS(nuc,0)
       
       AFW=trim(AFT(it))//'.bsw'
       
       kshift=0; IEF=0;     CALL SUB_check_orb

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

       kcfg = ncfg
       Do ic=1,ncfg; if(abs(WC(ic)).gt.eps_targ) Cycle
        WC(ic) = 0.d0;  kcfg = kcfg - 1 
       End do

       if(allocated(ipt)) deallocate(ipt); Allocate(ipt(ncfg))
       Call SORTA(ncfg,WC,ipt)

       Do jc=1,kcfg; ic=ipt(jc);  Call Pri_conf (muc,ic,WC(ic));  End do
       write(muc,'(a)') '*'
       write(nuo,'(a,3x,i3.3,6x,a)') 'target',it,AFT(it)

       Call SUB_phys_orb

       Call Def_term_BSR(nuc,ltarg(it),istarg(it),iptarg(it))
       nctarg(it) = kcfg
       nwtarg(it) = nbf-nwt
       nct=nct+kcfg; nwt=nbf  
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

! ... re-write the substitution orbitals:
 
      Open(muw,file='target_sub.bsw',form='UNFORMATTED')
      Do i = 1,nbf;  if(iech(i).ne.1) Cycle
       write(muw) ebs(i),z,h,hmax,rmax,ks,ns,mbs(i)
       write(muw) (p(jj,i),jj=1,mbs(i))
      End do
      Close(muw)

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

! ... check perturber cores:

      Do ilsp = 1,nlsp
       if(LEN_TRIM(AFP(ilsp)).eq.0) Cycle
       if(AFP(ilsp).eq.'no')  Cycle
       AFC = trim(AFP(ilsp))//'.c'
       Open(nuc,file=AFC)
       Call R_CLOSED(nuc)
       if(CLOSED.ne.core)  then
         write(pri,'(/a,a,a)') 'perturber ',AFC,' has different core'
         Stop ' check core ! ' 
       end if
       Close(nuc)
      End do

      Do iip=1,kpert
       AFC = trim(AFK(iip))//'.c'
       Open(nuc,file=AFC)
       Call R_CLOSED(nuc)
       if(CLOSED.ne.core)  then
         write(pri,'(/a,a,a)') 'perturber ',AFC,' has different core'
         Stop ' check core ! ' 
       end if
       Close(nuc)
      End do

!=====================================================================
! ... proceed perturbers if any:

      Do ilsp=1,nlsp

       nbf = nwt
       BFC = 'no'; BFW = 'no'
       ncfg = 0; lcfg = 0; nwf = ncore  

       npert = 0;  Call Allocate_pert(0)

!---------------------------------------------------------------------
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

!---------------------------------------------------------------------
! ... additional perturbers:

       Do iip=1,kpert

        if(klsp(iip).ne.ilsp) Cycle
        AFC = trim(AFK(iip))//'.c'
        Call Check_file(AFC)
        write(pri,'(/a,i4,4x,a/)') 'perturber', ilsp, trim(AFK(iip))
        Open(nuc,file=AFC)

        kshift = 0  !maxval(KEF(1:nwf))

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
        write(muc,'(10i8)') (ippert(i),i=1,npert)

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

      if(nlsp.gt.0) then
       Do i = 1,nlsp
        if(AFP(i).ne.'no'.or.BFP(i).ne.'no') then
          write(nut,'(a3,3i5,3x,a20,1x,a10,2i5)') &
          Tpar(i),lpar(i),ispar(i),ipar(i),AFP(i),BFP(i),ncp(i),nwp(i)
        else 
          write(nut,'(a3,3i5)') Tpar(i),lpar(i),ispar(i),ipar(i)
        end if
       End do
       write(nut,'(72(''-''))')
      end if

      if(kpert.gt.0) then
       write(nut,'(a,i5,T18,a)') &
            'kpert  =',kpert,   '!   number of additional perturbers' 
       write(nut,'(72(''-''))')
       Do i = 1,kpert
        write(nut,'(i3,3x,a)') klsp(i),trim(AFK(i))      
       End do
       write(nut,'(72(''-''))')
      end if

     End  ! program BSR_PREP






