!======================================================================
!     PROGRAM       B S R _ C O N F 3                 
!
!                   C O P Y R I G H T -- 2014
!
!     Written by:   Oleg Zatsarinny
!                   email: oleg_zoi@yahoo.com
!                        
!======================================================================
!     Generates configuration lists for BSR codes
!     in case of LS(LSJ), JK or JJ couplings
!     with automatic generation of minimum orthogonal constraints
!======================================================================
!
!     INPUT FILES:
!
!     bsr_par       -  parameters of calculations
!     target        -  description of target states and partial waves                        
!     target_orb    -  physical and substitution orbitals 
!     targ_nnn.c    -  c-files for target states
!     pert_nnn.c    -  c-files for perturber if any  
!
!     OUTPUT FILES:
!
!     bsr_conf.log  -  running information
!     cfg.nnn       -  c-file for given partial wave nnn
!     cfg.nnn_comp  -  compensation configurations if any    
!
!     INPUT ARGUMENT:  see subroutine read_arg.f90
!
!     REMARKS:
!
!     LS-coupling is default 
!     LSJ case is defined by zero value of total multiplicity (2S+1)
!     for partial waves (in this case L = 2J)  
!
!     ORTHOGONAL CONSTRAINTS: 
!
!     by default all continuum orbitals are supposed to be 
!     non-orthogonal to all others (bound and continuum);
!     common for all partial waves orthogonal constraints can be given 
!     in the bsr_par file as  
!
!                  < nl1| nl2>=0
!
!     orthogonal constraints, specific for given partial wave, 
!     can be given in cfg.nn file. It is supposed secondary running  
!     the BSR_CONF to get spesific set numbers for kl orbitals.
!
!     Program analizes the target expansions to put orthogonal
!     constraints, such that the weights of compansation configuration 
!
!                  W < c_comp
!
!======================================================================
! ... positions of configurations in the conf_LS list:
!
! 1.  physical configurations for each target state,
!
!     their pointers ->  ip_phys_conf
!     limit position ->  ncfg_phys
!
! 2.  target configurations   
!
!     their pointers ->  ic_targ  (=ictarg + ncfg_phys)
!     limit position ->  ncfg_targ
!
! 3.  scattering configurations (all or phys):
!
!     limit position ->  ncfg_sct
!
! 4.  perturber configurations if any:
!
!     limit position ->  ncfg_pert
!
! 5.  compensation configurations if any:
!
!     limit position ->  ncfg_comp
!======================================================================
! ... Program flow:
!
! 1.  read target information (file target)
! 2.  read argument (file bsr_par or command line)
! 3.  read target orbitals (file target.bsw)
! 4.  define spectroscopic target configurations (target.nnn + target_orb)
! 5.  read all target configurations  (files target.nnn)
! 6.  read "channel-delete" conditions if any
! 7.  read information about included partial waves (file target); 
!     this information should be prepared by hand - not convivient, 
!     I am sorry, but you may copy these lines from previous calculations) 
!----------------------------------------------------------------------
! ... loop other partial waves:  
!----------------------------------------------------------------------
! 8.  define and record channel configurations in cfg.nnn
! 9.  record pertuber configurations if any
!10.  record imposed orthogonality constrans
!11.  record channels information in file target
!----------------------------------------------------------------------
! ... check if we need additional orth.conditions:  
!----------------------------------------------------------------------
!12.  define scattering channels with phys.target states only
!13.  add pertuber physical configurations and check double-counting
!14.  define orthogonal conditions (routine def_orth_cond) and
!     record all compensation configurations in file "cfg.nnn_comp"
!15.  record "derived" orth.constraints in cfg.nnn
!---- end of partial waves loop   
!...  finish by recording overall information in file target    
!----------------------------------------------------------------------
      Use bsr_conf
      Use target; Use channel; Use conf_LS;; Use orb_LS
      Use phys_orb_LS

      Implicit real(8) (A-H,O-Z)
      Character(80) :: AS
      Character(80), allocatable :: line(:)

      Call bsr_conf_inf
!---------------------------------------------------------------------- 
! ... read target information:

      Call Check_file(AF_tar); Open(nut,file=AF_tar)
      Call R_target(nut)

! ... read arguments:

      Call Read_arg

! ... target orbitals:

      Open(nuw,file=AF_wfn,form='UNFORMATTED')
      Call Read_bsw_orb_LS(nuw)
      Close(nuw)
      if(nwf.ne.nwt) Stop 'nwf in target.bsw <> nwt'
      write(pri,'(/a,T33,i8)') 'number of target orbital:',nwf
      max_ll_targ = maxval(LEF(1:nwf))

! ... find target physical configurations:

      Call Def_phys_targ

! ... target configurations:

      Do it=1,ntarg
       AF = trim(AFT(it))//'.c'
       Open(nuc,file=AF,status='OLD')
       if(it.eq.1) Call R_closed(nuc) 
       Call Add_conf_LS(nuc,0)
       Close(nuc)
      End do    
      Call Test_a

      write(pri,'(/a,T33,i8/)') &
       'number of target configurations:',ncfg-ncfg_phys
      if(ncfg-ncfg_phys.ne.nct) Stop 'ncfg_target <> nct from target'

      ncfg_targ=ncfg; lcfg_targ=lcfg; nwf_targ=nwf

! ... read partial waves:

      kpert=0; Call Read_ipar(nut,'kpert',kpert)
      Call Read_ipar(nut,'nlsp',nlsp)
      if(nlsp.le.0) Stop 'bsr_conf: nlsp <= 0 ? '
      write(pri,'(72(''=''))')
      write(pri,'(a,i3)') 'Number of partial waves, nlsp =',nlsp
      write(pri,'(72(''=''))')
      Allocate(line(nlsp))
      read(nut,*); Do i=1,nlsp; read(nut,'(a)') line(i); End do

      if(kpert.gt.0) then;  Do i=1,kpert+3; read(nut,*); End do; end if
      write(nut,'(72(''-''))')
      write(nut,'(a)') 'channels:'
      write(nut,'(72(''-''))')

! ... read delete-channels information:

      Call Read_del   

!======================================================================
! ... loop other partial waves:

      max_ch=0;  max_nc=0;  max_wf=0
 
      Do ilsp=1,nlsp

      AF=line(ilsp); ncp = 0; nwp = 0
      read(AF,*)  Tpar,lpar,ispar,ipar
      if(LEN_trim(AF).gt.18) read(AF(19:),*) AFP,AFP,ncp,nwp

      write(pri,'(/a,i3,a,3i4/)') &
        'Partial wave ',ilsp, '   LSP =', lpar,ispar,ipar

      Jpar=lpar+1; if(coupling.eq.'LS'.and.ispar.gt.0) Jpar=0

      ncfg=ncfg_targ; lcfg=lcfg_targ; nwf=nwf_targ 

! ... read perturber configurations to get pertuber orbitals if any
! ... (to be sure that set indexes are same as in bsr_prep)

      if(ncp.ne.0)  then
       AF = trim(AFP)//'.c'
       Open(nuc,file=AF,status='OLD')
       Call Add_conf_LS(nuc,0)
      end if
      ncfg=ncfg_targ; lcfg=lcfg_targ; nwf_pert=nwf 

! ... define channels configurations:

      ic_targ(0) = ncfg_phys
      ic_targ(1:ntarg) = ictarg(1:ntarg) + ncfg_phys

      nch=0; Call Allocate_channel(imch)
      c_norm = 0.1   !  tollerance for target states normalization
      Select Case (COUPLING)
       Case('LS');  Call SUB_LS
       Case('JK');  Call SUB_JK
       Case('JJ');  Call SUB_JJ
       Case default; Stop 'unknown coupling, possible: LS,JK,JJ'
      End Select

      ncfg_sct=ncfg; lcfg_sct=lcfg; nwf_sct=nwf
      write(pri,'(/a,T40,i8/)') &
       'Number of scattering configutarions:', ncfg_sct-ncfg_targ

      min_ll_ch=0; if(nch.gt.0) min_ll_ch = minval(lch(1:nch))

! ... perturber configurations:

      if(ncp.ne.0)  then
       AF = trim(AFP)//'.c'
       Open(nuc,file=AF,status='OLD')
       Call Add_conf_LS(nuc,0)
      end if
      ncfg_pert = ncfg
   
      if(ncfg.gt.ncfg_sct) write(pri,'(a,T40,i8)') &
        'Total number of configurations:', ncfg-ncfg_targ

! ... output cfg.nnn:

      write(AF,'(a,i3.3)') 'cfg.',ilsp
      Open (nuc,file=AF); rewind(nuc)
      write(nuc,'(12x,a3,f16.8)') Tpar,Etarg(1)
      write(nuc,'(a60)') CLOSED 
      Do ic=ncfg_targ+1,ncfg; Call Pri_conf(nuc,ic,WC(ic)); End do
      write(nuc,'(a)') '*'

! ... output the imposed orth. conditions

      Call Pre_iort(nup,0)                   
      if(kort.gt.0) Call R_orth(nuc)

      if(max_ll_targ.gt.min_ll_ch) &
       write(nuc,'(/a/)') 'Imposed orth. conditions:'

      Do ich=1,nch; i=ipch(ich)
       Do j=nclosd+1,nwf 
        if(LEF(i).ne.LEF(j)) Cycle
        if(IORT(max(i,j),min(i,j)).ne.0) Cycle
        write(nuc,'(a1,a4,a1,a4,a2,i1)') '<',ELF(i),'|',ELF(j),'>=',0
        IORT(max(i,j),min(i,j)) = -1
       End do
      End do

! ... maximum dimensions:  
 
      max_nc = max(max_nc,ncfg_pert-ncfg_targ)
      max_wf = max(max_wf,nwf_pert)
      max_ch = max(max_ch,nch)

! ... recording channel information for given LSP in file target:

      write(nut,'(i3,a,2x,a3,a,i5,a,2i8)')  ilsp,'.',Tpar,  &
             '  nch =',nch, '  nc =',ncfg-ncfg_targ,ncfg-ncfg_sct
      Do ich=1,nch
        write(nut,'(a4,3i5,3x,i8,i5)') &
         ELF(ipch(ich)),lch(ich),iptar(ich),ich,ipconf(ich),jkch(ich)
      End do
      write(nut,'(72(''-''))') 

!======================================================================
! ... orthogonal conditions:   
!======================================================================
 
      if(max_ll_targ.lt.min_ll_ch) Cycle     

      ncfg=ncfg_targ; lcfg=lcfg_targ; nwf=nwf_pert 

      ic_targ=jc_targ

! ... define scattering channels with phys.target states:

      nch_save = nch 
      nch=0; Call Allocate_channel(imch)
      c_norm = 0.5   !  tollerance for spectroscopic states normalization
      Select Case (COUPLING)
       Case('LS');  Call SUB_LS
       Case('JK');  Call SUB_JK
       Case('JJ');  Call SUB_JJ
      End Select
      if(nch.ne.nch_save) Stop 'nch <> nch_save'

      ncfg_sct=ncfg; lcfg_sct=lcfg; nwf_sct=nwf

      Do ic=ncfg_targ+1,ncfg; Call Pri_conf(nuc,ic,WC(ic)); End do
      write(nuc,'(a)') '*'


! ... add pertuber physical configurations: 

      if(ncp.gt.0) Call Check_perturbers
     
      ncfg_pert=ncfg; lcfg_pert=lcfg; nwf_pert=nwf

! ... define orthogonal conditions

      Call Def_orth_cond

! ... analize and record orthogonal conditions for sct. orbitals:

      Call Record_orth

      write(pri,'(/72(''-''))')

      End do ! over partial waves (ilsp)
!----------------------------------------------------------------------
! ... overall information:

      write(nut,'(a,i8,a)') 'max_ch  = ',max_ch, &
        '  -  max.number of channels'
      write(nut,'(a,i8,a)') 'max_nc  = ',max_nc, &
        '  -  max.number of configurations'
      write(nut,'(a,i8,a)') 'max_wf  = ',max_wf, &
        '  -  max.number of orbitals'
      write(nut,'(72(''-''))')

     End  ! program BSR_CONF


