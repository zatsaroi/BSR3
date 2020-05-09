!======================================================================
!     PROGRAM   B S R _ B M A T                             version 4
!
!               C O P Y R I G H T -- 2017
!
!     Written by:   Oleg Zatsarinny 
!                   email: oleg_zoi@yahoo.com
!======================================================================
!    generates angular coefficient for specific set of states
!    with expansion coefficients (including the overlap factors)
!----------------------------------------------------------------------
!
!    INPUT ARGUMENTS:
!    
!    klsp,klsp1,klsp2  - range of partial wave in BSR calculations,
!                        then cfg.001, cfg.002, ..., are input files
!                        (default -> 0, with input file is cfg.inp)
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
!    mk    - max.multipole index (default -> 9)
!
!-----------------------------------------------------------------------
!
!    example:    1.  bsr_bmat 
!                2.  bsr_bmat klsp1=1 klsp2=5 oper=1111110  
!                3.  bsr_bmat km=5
!            
!------------------------------------------------------------------------
!
!    INPUT FILES:
!    
!    cfg.nnn      -  configuration list for partial wave nnn = klsp
!                    (cfg.inp in case klsp = 0, default)
!                  
!    int_list.nnn -  already determined set of angular coefficients
!                    (optional; int_list in case klsp = 0)
!                   
!    
!    OUTPUT FILES:
!    
!    int_list.nnn  -  output set of angular coefficients
!                     (int_list in case klsp = 0)
!                   
!-------------------------------------------------------------------------     
      Use bsr_breit
      Use bsr_mat
      Use conf_LS,      only: ne,ncfg
      Use term_exp,     only: ic_case, is_need
      Use symc_list_LS, only: nsymc
      Use symt_list_LS, only: nsymt
      Use orb_LS,       only: nwf  

      Implicit none 
      Integer :: i,l,mls_max
      Integer, external :: Ifind_channel

!-----------------------------------------------------------------------
! ... HEADER:

      open(pri,file=AF_p) 
      write(pri,'(20x,a/20x,a/20x,a)') &
             '===================================',     &
             ' BREIT-PAULI  ANGULAR COEFFICIENTS',     &
             '==================================='

! ... prepare B-spline parameters:

      Call def_bs

! ... target information:

      Open(nut,file=AF_tar,status='OLD')
      Call R_target(nut)

! ... read arguments from command line:

      Call Read_arg 

!----------------------------------------------------------------------
!                                             cycle over partial waves:
      time = 0.d0
      Do klsp = klsp1,klsp2

       write(pri,'(80(''-''))') 
       if(klsp.gt.0) write(pri,'(/a,i5)') 'Partial wave: ',klsp
       if(klsp.gt.0) write(*  ,'(/a,i5)') 'Partial wave: ',klsp
       Call CPU_TIME(t1)

! ...  read the configuration list:

       Call CPU_TIME(t3)

       Call open_c_file       ! c-file
       Call Read_conf(nuc)

       Call CPU_TIME(t4)

       write(pri,'(/a/)')   'c-file data: '
       write(pri,'(a,i10,a)')  'ncfg   = ',ncfg,&
         ' -  number of configurations'
       write(pri,'(a,i10,a)')  'nsymc  = ',nsymc,&
         ' -  number of conf. symmetries'
       write(pri,'(a,i10,a)')  'nsymt  = ',nsymt,&
         ' -  number of term  symmetries'
       write(pri,'(a,i10,a)')  'nwf    = ',nwf, &
         ' -  number of orbitals: '
       write(pri,'(/a,f10.2,a)') 'Read_conf:     ',(t4-t3)/60,' min '

! ... define possible mls orbitals:
      
       Call Def_maxl(l);  mls_max=4*l+2
       Call Alloc_spin_orbitals(ne,mls_max)

! ... check  det. expantions:

       Call CPU_TIME(t3)

       Call open_det_exp

       Call open_det_done

       i = count(IS_need.gt.0) 
       write(pri,'(/a,i10)')  'Needed matrix elements: ',i
       if(i.eq.0) Cycle

       Call CPU_TIME(t4)
       write(pri,'(/a,f10.2,a)') 'Pre_det_exp:   ',(t4-t3)/60,' min '

! ... read target information, radial orbitals and define overlaps:

       Call CPU_TIME(t3)

       Call Read_data

       if(allocated(IP_channel)) Deallocate(IP_channel)
       Allocate(IP_channel(ncfg))
       Do i=1,ncfg; IP_channel(i)=Ifind_channel(i); End do

       Call CPU_TIME(t4)
       write(pri,'(/a,f10.2,a/)') 'Read_data:     ',(t4-t3)/60,' min '

! ... check the int_list file if any: 

       Call open_int_list     ! data bank information

! ... nulify data in c_blocks module

       Call alloc_c_blocks(0,mb,nb,kb,eps_c) 

! ... calculations:

       Call Conf_loop 

! ... time for one partial wave:

       Call CPU_TIME(t2); tt=(t2-t1)/60
       write(pri,'(/a,T13,F12.2,a)') 'Partial wave:',tt,' min'
       write(*,  '(/a,T13,F12.2,a)') 'Partial wave:',tt,' min'
       time = time + tt

      End do  ! over klsp

!----------------------------------------------------------------------

! ... total time:

      write(pri,'(/a,T13,F12.2,a)') 'Total time: ',time,' min'
      write(*,  '(/a,T13,F12.2,a)') 'Total time: ',time,' min'

      End ! Program Breit_bsr


