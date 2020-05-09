!=====================================================================
!     PROGRAM   B S R _ B M A T _ M P I                 version: 4
!
!               C O P Y R I G H T -- 2017
!
!     Written by:   Oleg Zatsarinny
!                   email: oleg_zoi@yahoo.com
!======================================================================
!    generates angular coefficient in non-orthogonal mode 
!    for the given set of ACS with simultanious evaluation of all 
!    overlap factors
!----------------------------------------------------------------------
!
!    INPUT ARGUMENTS:
!    
!    klsp1,klsp2  - range of partial wave in BSR calculations,
!                   then cfg.001, cfg.002, ..., are input files
!                   (default -> 0, with input file is cfg.inp)
!   
!    oper  - character(7), where each position can be 0 or 1,
!            and indicate the operator under consideration:
!            oper(1) - OVERLAPS       
!            oper(2) - KINATIC ENERGY
!            oper(3) - TWO-ELECTRON ELECTROSTATIC
!            oper(4) - SPIN ORBIT       
!            oper(5) - SPIN-OTHER-ORBIT
!            oper(6) - SPIN-SPIN
!            oper(7) - ORBIT-ORBIT       
!            Default -> 1110000 - non-relativistic calculations
!
!    mk    - max.multipole index (default -> 9, see module param_br)
!
!----------------------------------------------------------------------
!
!    example:    1.  bsr_bmat 
!                2.  bsr_bmat klsp=5 oper=1111110  
!                3.  bsr_bmat km=5
!            
!----------------------------------------------------------------------
!
!    INPUT FILES:
!    
!    cfg.nnn      -  configuration list for partial wave nnn = klsp
!                   (cfg.inp in case klsp = 0, default)
!                  
!    int_list.nnn -  input data bank for angular coefficients
!                   (optional; int_list in case klsp = 0)
!                   
!    
!    OUTPUT FILES:
!    
!    int_list.nnn  - output data bank for angular coefficients
!                   
!---------------------------------------------------------------------     
      Use MPI

      Use bsr_mat
      Use bsr_breit
      Use conf_LS,       only: ne, ncfg
      Use term_exp,      only: ic_case,IS_need

      Implicit none 
      Integer :: l,i,j,mls_max
      Integer, external :: Ifind_channel
      Character(80) :: AAS,BBS,CCS
      Real(8) :: S,SS
!----------------------------------------------------------------------
! ... initialize MPI:

      Call MPI_INIT(ierr)
      Call MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierr)
      Call MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierr)

!----------------------------------------------------------------------
! ... log - files:

      if(myid.eq.0) then
       open(pri,file=AF_p)

       write(pri,'(/a/a/a/)') &
             '===========================================',     &
             ' B R E I T - M A T  C A L C U L A T I O N S',     &
             '==========================================='
       write(  *,'(/a,i6)') 'nprocs = ', nprocs
       write(pri,'(/a,i6)') 'nprocs = ', nprocs
       Allocate(ip_proc(nprocs)); ip_proc=0
      end if

      t0 = MPI_WTIME()
!----------------------------------------------------------------------
! ... read arguments from command line:

      if(myid.eq.0) Call Read_arg;    Call br_arg   

      if(myid.gt.0) then
       write(AF_p,'(a,i4.4)') 'debug_',myid
       if(debug.eq.0) pri=0
       if(pri.gt.0) open(pri,file=AF_p)
      end if

!----------------------------------------------------------------------
! ... prepare B-spline parameters:

      if(myid.eq.0) Call define_grid(z)
      Call br_grid
      Call define_spline        

!----------------------------------------------------------------------
! ... target information:

      if(myid.eq.0) then
       Open(nut,file=AF_tar,status='OLD')
       Call R_target(nut)
      end if
      Call br_target

!----------------------------------------------------------------------
! ... read the configuration list:

      t1 = MPI_WTIME()

      if(myid.eq.0) Call open_c_file       
      if(myid.eq.0) Call Read_conf(nuc)

       Call br_symc_LS
       Call br_symt_LS
       Call br_conf_LS
 
       Call MPI_BCAST(ne,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
 
       Call Def_maxl(l);  mls_max=4*l+2
       Call Alloc_spin_orbitals(ne,mls_max)

      t2=MPI_WTIME()
      if(pri.gt.0) &
      write(pri,'(/a,F12.2,a)') 'Read_conf done:',(t2-t1)/60,' min'

!----------------------------------------------------------------------
! ... prepare det. expantions:

       if(myid.eq.0) then
        Call open_det_exp
        Call open_det_done
        i = count(IS_need(:).gt.0); S = i 
        SS = ic_case; SS = (SS-1.d0)*SS/2.d0
        write(pri,'(/a,f5.2,a)')  'Needed matrix elements: ',S/SS*100,' %'
        if(i.eq.0) go to 100
       end if

       t3=MPI_WTIME()
       if(pri.gt.0) &
       write(pri,'(/a,F12.2,a)') 'Det_exp is done:',(t3-t2)/60,' min'
       Call MPI_BCAST(ic_case,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
       Call MPI_BARRIER(MPI_COMM_WORLD, ierr)
       if(ic_case.eq.0) go to 100  
 
! ... read target information, radial orbitals and define overlaps:

      if(myid.eq.0)  Call Read_data

      t3=MPI_WTIME()
      if(pri.gt.0) &
      write(pri,'(/a,F12.2,a)') 'Read_data done:',(t3-t2)/60,' min'

      Call MPI_BARRIER(MPI_COMM_WORLD, ierr)

      Call br_channel
      Call br_bsorb
      Call br_ovl

      t4=MPI_WTIME()
      if(pri.gt.0) &
      write(pri,'(/a,F12.2,a)') 'All broadcast  done:',(t4-t3)/60,' min'

      if(allocated(IP_channel)) Deallocate(IP_channel)
      Allocate(IP_channel(ncfg))
      Do i=1,ncfg; IP_channel(i)=Ifind_channel(i); End do

! ... define output file:

      if(myid.eq.0) Call open_int_list     ! data bank name

      i=Index(BF_b,'.'); write(BF_b(i+1:),'(i3.3)') klsp 

! ... nulify data in c_blocks module

      Call alloc_c_blocks(0,mb,nb,kb,eps_c) 

! ... calculations for new angular symmetries:

      if(myid.eq.0) open(nuf,file=AF_f)

      if(myid.eq.0)     Call Conf_loop 
      if(myid.gt.0)     Call Conf_calc
      
! ... time for one partial wave:

  100  t2=MPI_WTIME()
       if(pri.gt.0) &
       write(pri,'(/a,F12.2,a)') 'Partial wave:',(t2-t1)/60,' min'

      Call MPI_FINALIZE(l)

      End ! Program BR_MPI



