!======================================================================
!     PROGRAM       B S R _ C O N F _ M P I       
!
!                   C O P Y R I G H T -- 2016
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
      Use MPI
      Use bsr_conf
      Use target; Use channel; Use conf_LS;; Use orb_LS
      Use phys_orb_LS
      Use internal_file

      Implicit real(8) (A-H,O-Z)
      Character(124) :: AS
      Integer, external :: Iadd_line

!---------------------------------------------------------------------- 
      Call MPI_INIT(ierr)
      Call MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierr)
      call MPI_COMM_SIZE( MPI_COMM_WORLD, nprocs, ierr )

      if(myid.eq.0) write(*,'(a,i5)') 'MPI: nprocs = ',nprocs

      if(myid.eq.0) Call bsr_conf_inf

!---------------------------------------------------------------------- 
! ... read target information:

      if(myid.eq.0) then   
       Call Check_file(AF_tar)
       Open(nut,file=AF_tar)
       Call R_target(nut)
      end if
      Call br_target

!---------------------------------------------------------------------- 
! ... read arguments:

      if(myid.eq.0) Call Read_arg
      Call br_arg

!---------------------------------------------------------------------- 
! ... debug output

      if(myid.ne.0) then
       if(debug.ge.myid) then
        write(AF_log,'(a,i5.5)') 'debug.',myid
        Open(pri,file=AF_log)
       else
        pri = 0
       end if
      end if

!---------------------------------------------------------------------- 
! ... target orbitals:

      if(myid.eq.0) then
       Open(nuw,file=AF_wfn,form='UNFORMATTED')
       Call Read_bsw_orb_LS(nuw)
       Close(nuw)
      end if

      Call br_orb
      if(nwf.ne.nwt) Call Stop_mpi(pri,0,'nwf in target.bsw <> nwt')
      max_ll_targ = maxval(LEF(1:nwf))

      if(pri.gt.0) &
      write(pri,'(a,T33,i8)') 'number of target orbital:',nwf,max_ll_targ

!---------------------------------------------------------------------- 
! ... find target physical configurations:

      if(myid.eq.0) Call Def_phys_targ

      Call br_phys_orb 
 
      Call MPI_BARRIER(MPI_COMM_WORLD,ierr)

      if(myid.ne.0)  Allocate(ic_targ(0:ntarg),jc_targ(0:ntarg))

      Call MPI_BARRIER(MPI_COMM_WORLD,ierr)

      Call MPI_BCAST(ic_targ(0:ntarg),ntarg+1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(jc_targ(0:ntarg),ntarg+1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(ncfg_phys,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(lcfg_phys,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

!---------------------------------------------------------------------- 
! ... target configurations:

      if(myid.eq.0) then
       Do it=1,ntarg
        AF = trim(AFT(it))//'.c'
        Open(nuc,file=AF,status='OLD')
        if(it.eq.1) Call R_closed(nuc) 
        Call Add_conf_LS(nuc,0)
        Close(nuc)
       End do    
       Call Test_a
      end if

      Call MPI_BARRIER(MPI_COMM_WORLD,ierr)

      Call br_core1

      Call br_symc_LS
      Call br_symt_LS
      Call br_conf_LS1

      if(pri.gt.0) write(pri,'(/a,T33,i8/)') &
       'number of target configurations:',ncfg-ncfg_phys

      if(ncfg-ncfg_phys.ne.nct) &
        Call Stop_mpi(pri,0,'ncfg_target <> nct from target')

      ncfg_targ=ncfg; lcfg_targ=lcfg; nwf_targ=nwf

!---------------------------------------------------------------------- 
! ... read partial waves:

      if(myid.eq.0) then

       kpert=0; Call Read_ipar(nut,'kpert',kpert)
       nlsp =0; Call Read_ipar(nut,'nlsp',nlsp)
       Call alloc_file(nlsp+kpert)
       read(nut,*)
       Do i=1,nlsp
        read(nut,'(a)') AS
        j = Iadd_line(AS)
       End do
       if(kpert.gt.0) then
        read(nut,*);  read(nut,*);  read(nut,*)
        Do i=1,kpert
         read(nut,'(a)') AS
         j = Iadd_line(AS)
        End do
       end if

      end if

      Call br_file

      Call MPI_BCAST(nlsp,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(kpert,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
 
      if(nlsp.le.0) Call Stop_mpi(pri,0,'bsr_conf: nlsp <= 0 ? ')

      if(pri.gt.0) then
       write(pri,'(72(''=''))')
       write(pri,'(a,i3)') 'Number of partial waves, nlsp =',nlsp
       write(pri,'(72(''=''))')
      end if
!-----------------------------------------------------------------------
! ... read delete-channels information:

      if(myid.eq.0)  Call Read_del   
      Call MPI_BCAST(ndel,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

      if(myid.ne.0.and.ndel.gt.0) &
       Allocate(dlsp(ndel),dlch(ndel),dtar(ndel),djk(ndel))

      if(ndel.gt.0) then
       Call MPI_BCAST(dlsp,ndel,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
       Call MPI_BCAST(dlch,ndel,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
       Call MPI_BCAST(dtar,ndel,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
       Call MPI_BCAST(djk ,ndel,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      end if

!======================================================================
! ... loop other partial waves:
!======================================================================

      max_ch=0;  max_nc=0;  max_wf=0
 
      klsp = 0
      Do ilsp=1,nlsp

       klsp=klsp+1; if(klsp.gt.nprocs-1) klsp=1; if(myid.ne.klsp) Cycle 

      AS=aline(ilsp); ncp = 0; nwp = 0

      read(AS,*)  Tpar,lpar,ispar,ipar
      if(LEN_trim(AS).gt.18) read(AS(19:),*) AFP,AFP,ncp,nwp

      if(pri.gt.0) write(pri,'(/a,i3,a,3i4/)') &
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
       Case default
        Call Stop_mpi(pri,0,'unknown coupling, possible: LS,JK,JJ')
      End Select

      ncfg_sct=ncfg; lcfg_sct=lcfg; nwf_sct=nwf
      if(pri.gt.0) write(pri,'(/a,T40,i8/)') &
       'Number of scattering configutarions:', ncfg_sct-ncfg_targ

      min_ll_ch=0; if(nch.gt.0) min_ll_ch = minval(lch(1:nch))

! ... perturber configurations:

      if(ncp.ne.0)  then
       AF = trim(AFP)//'.c'
       Open(nuc,file=AF,status='OLD')
       Call Add_conf_LS(nuc,0)
      end if

      ncfg_pert = ncfg
   
      if(ncfg.gt.ncfg_sct.and.pri.gt.0) write(pri,'(a,T40,i8)') &
        'Total number of configurations:', ncfg-ncfg_targ

! ... output cfg.nnn:

      write(AF,'(a,i3.3)') 'cfg.',ilsp
      Open (nuc,file=AF)

      Call Pre_iort(0,0)          !  nup  ???         

      if(kort.gt.0) Call Read_orth(nuc)

      rewind(nuc)
      write(nuc,'(12x,a3,f16.8)') Tpar,Etarg(1)
      write(nuc,'(a60)') CLOSED 
      Do ic=ncfg_targ+1,ncfg; Call Pri_conf(nuc,ic,WC(ic)); End do
      write(nuc,'(a)') '*'

! ... output the imposed orth. conditions

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

      write(AS,'(i3,a,2x,a3,a,i5,a,2i8)')  ilsp,'.',Tpar,  &
        '  nch =',nch, '  nc =',ncfg-ncfg_targ,ncfg-ncfg_sct
       j = Iadd_line(AS)
      Do ich=1,nch
        write(AS,'(a4,3i5,3x,i8,i5)') &
         ELF(ipch(ich)),lch(ich),iptar(ich),ich,ipconf(ich),jkch(ich)
        j = Iadd_line(AS)
      End do

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

      if(pri.gt.0) write(pri,'(/72(''-''))')

      End do ! over partial waves (ilsp)

!----------------------------------------------------------------------
! ... Collect the information:

      if(myid.eq.0) then
       write(nut,'(72(''-''))')
       write(nut,'(a)') 'channels:'
       write(nut,'(72(''+''))')
      end if

      Call MPI_BARRIER(MPI_COMM_WORLD,ierr)

      ip = nlsp+kpert+1
      Do ilsp=1,nlsp

       Call MPI_BARRIER(MPI_COMM_WORLD,ierr)

       if(myid.ne.0) then
        if(ip.gt.nlines) Cycle
        AS = aline(ip);  read(AS(1:3),*) klsp; if(klsp.ne.ilsp) Cycle
        Call MPI_SEND(AS,124, MPI_CHARACTER,0,i,MPI_COMM_WORLD, ierr)
       end if

       if(myid.eq.0) then
        Call MPI_RECV(AS,124, MPI_CHARACTER, MPI_ANY_SOURCE, &
                      MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr)
        write(nut,'(a)') trim(AS)
        write(nut,'(72(''-''))')      ! ???
       end if

       read(AS,*) AF,AF,AF,AF,nch

       Do ich = 1,nch
        if(myid.ne.0) then
         AS = aline(ip+ich)
         Call MPI_SEND(AS, 124, MPI_CHARACTER,0,i,MPI_COMM_WORLD, ierr)
        end if
        if(myid.eq.0) then
         Call MPI_RECV(AS, 124, MPI_CHARACTER, MPI_ANY_SOURCE, &
                       MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr)
         write(nut,'(a)') trim(AS)
        end if
       End do
        
       if(myid.eq.0) write(nut,'(72(''-''))') 

       ip = ip + 1 + nch

      End do

!----------------------------------------------------------------------
! ... overall information:

      Call MPI_BARRIER(MPI_COMM_WORLD,ierr)

      Call MPI_REDUCE(max_ch,m,1,MPI_integer,MPI_MAX,0,MPI_COMM_WORLD,ierr)
      if(myid.eq.0) write(nut,'(a,i8,a)') &
         'max_ch  = ', m, '  -  max.number of channels'
      Call MPI_REDUCE(max_nc,m,1,MPI_integer,MPI_MAX,0,MPI_COMM_WORLD,ierr)
      if(myid.eq.0) write(nut,'(a,i8,a)') &
         'max_nc  = ', m, '  -  max.number of configurations'
      Call MPI_REDUCE(max_wf,m,1,MPI_integer,MPI_MAX,0,MPI_COMM_WORLD,ierr)
      if(myid.eq.0) write(nut,'(a,i8,a)') &
         'max_wf  = ',m,  '  -  max.number of orbitals'
      if(myid.eq.0) write(nut,'(72(''-''))')

!      if(myid.eq.0) then 
!       if(debug.eq.0) Call System('rm targ_*')
!       if(debug.eq.0) Call System('rm pert_comp.*')
!      end if

      Call MPI_FINALIZE(ierr)

     End  ! program BSR_CONF_MPI






