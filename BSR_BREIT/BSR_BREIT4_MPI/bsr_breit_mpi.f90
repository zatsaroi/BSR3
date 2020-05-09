!=====================================================================
!     PROGRAM   B S R _ B R E I T _ M P I                 version: 4
!
!               C O P Y R I G H T -- 2016
!
!     Written by:   Oleg Zatsarinny
!                   email: oleg_zoi@yahoo.com
!======================================================================
!    generates angular coefficient in non-orthogonal mode 
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
!    example:    1.  bsr_breit 
!                2.  bsr_breit klsp1=1 klsp2=5 oper=1111110  
!                3.  bsr_breit km=5
!            
!----------------------------------------------------------------------
!
!    INPUT FILES:
!    
!    cfg.nnn     -  configuration list for partial wave nnn = klsp
!                   (cfg.inp in case klsp = 0, default)
!                  
!    int_bnk.nnn -  input data bank for angular coefficients
!                   (optional; int_bnk in case klsp = 0)
!                   
!    
!    OUTPUT FILES:
!    
!    int_bnk.nnn  - output data bank for angular coefficients
!                   (int_bnk in case klsp = 0)
!                   
!---------------------------------------------------------------------     
      Use MPI

      Use bsr_breit
      USE conf_LS,      only: ne
      Use symc_list_LS, only: nsymc
      Use term_exp,     only: ic_case

      Implicit none 
      Integer :: l,mls_max
 
!----------------------------------------------------------------------
! ... initialize MPI:

      Call MPI_INIT(ierr)
      Call MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierr)
      Call MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierr)

      if(myid.eq.0) then
       open(pri,file=AF_p)
       write(*,'(a,i6)') 'nprocs = ', nprocs
       write(pri,'(a,i6)') 'nprocs = ', nprocs
       Allocate(ip_proc(nprocs)); ip_proc=0
      end if

      t0 = MPI_WTIME()
!----------------------------------------------------------------------
! ... read arguments from command line:

      if(myid.eq.0) Call Read_arg;    Call br_arg   

! ... log - file:

      write(AF_p,'(a,i3.3)') 'br_',myid

      if(debug.eq.0.and.myid.gt.0) pri=0
      if(myid.gt.0.and.pri.gt.0) open(pri,file=AF_p)

      if(pri.gt.0) then
      write(pri,'(/20x,a/20x,a/20x,a/)') &
             '=======================',     &
             ' B R E I T - P A U L I ',     &
             '======================='
      write(pri,'(/a,i3)')     'Max.multipole index =',mk
      write(pri,'(/a,E10.1/)') 'Tollerance for coeff.s =',eps_c
      end if

!----------------------------------------------------------------------
!                                             cycle over partial waves:
      Do klsp = klsp1,klsp2

       t1 = MPI_WTIME()

! ... open relavent files: 

       if(myid.eq.0) then
        Call open_c_file
        Call open_int_inf
        write(pri,'(80(''-''))') 
        write(pri,'(/a,i5/)') ' Partial wave: ',klsp
        if(new.eq.1) write(pri,'(a/)') ' It is new calculations '
        if(new.eq.0) write(pri,'(a/)') ' It is continued calculations '
       end if 

! ... read the configuration list:

       if(myid.eq.0) Call Read_conf

       Call MPI_BCAST(icalc,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

       if(icalc.eq.0) Cycle

       Call MPI_BCAST(ne,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
 
       if(myid.eq.0) then; Call Def_maxl(l);  mls_max=4*l+2; end if

       Call MPI_BCAST(mls_max,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

       Call Alloc_spin_orbitals(ne,mls_max)  

! ... extract old results: 

       if(myid.eq.0)   Call Read_dets(nub,new)

! ... prepare det. expantions:

       if(myid.eq.0)   Call open_det_exp

       Call MPI_BARRIER(MPI_COMM_WORLD, ierr)

       t2=MPI_WTIME()
       if(pri.gt.0) &
       write(pri,'(/a,F12.2,a)') 'Prep_det is done:',(t2-t1)/60,' min'

       Call open_int_int

! ... calculations for new angular symmetries:

       if(myid.eq.0)     Call Conf_loop 
       if(myid.gt.0)     Call Conf_calc

! ... record results:

       if(myid.eq.0) then 

        write(pri,'(/a)') ' Record results:'

        Call Record_results

       end if

! ... time for one partial wave:

       t2=MPI_WTIME()
       if(pri.gt.0) &
       write(pri,'(/a,F12.2,a)') ' Partial wave:',(t2-t1)/60,' min'

      End do  ! over klsp

      if(debug.gt.0) Call Debug_printing

      Call MPI_FINALIZE(l)

      End ! Program BSR_BREIT_MPI



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


!======================================================================
      Subroutine Record_results
!======================================================================
      Use bsr_breit
      Use det_list
      Use def_list

      Implicit none
      Character AS*80
      Integer :: i, ii, nc
      Real(8) :: adet,adef

      close(nur)

      rewind(nub)
      Call Write_symc_LS(nub)
      Call Write_symt_LS(nub)
      Call Record_oper_LS(nub)

      adet=ldet; adet=adet/ndet
      Call Write_det(nub)
      adet=ldef; adef=adef/ndef
      Call Write_def(nub)

      close(nub)

! ... print the main dimensions:      

      write(pri,'(/a/)') &
          ' Results for new angular symmetry calculations:'
      write(pri,'(a,i10,f10.1,i10)') &
          ' number of overlap determinants =', ndet,adet,ldet
      write(pri,'(a,i10,f10.1,i10)') &
          ' number of overlap factors      =', ndef,adef,ldef 
      write(pri,'(a,i10)') &
          ' new coeff.s                    =', nc_new

            
      End Subroutine Record_results


!======================================================================
      Subroutine Debug_printing
!======================================================================
      Use MPI

      Use bsr_breit
      Use coef_list
      Use zoef_list
      Use boef_list

      Integer :: nn, mt, nreloc, mreloc
      Real(8) :: mem, meb

      if(pri.gt.0) write(pri,'(/a/)') 'debug printing:'

      Call MPI_REDUCE(mem_max_zoef,mem,1,MPI_DOUBLE_PRECISION,MPI_MAX,0,MPI_COMM_WORLD,ierr)
      Call MPI_REDUCE(max_zoef,nn,1,MPI_integer,MPI_MAX,0,MPI_COMM_WORLD,ierr)
      Call MPI_REDUCE(zoef_realloc,nreloc,1,MPI_integer,MPI_MAX,0,MPI_COMM_WORLD,ierr)

      if(pri.gt.0) then
       write(pri,'(a,f10.1)') 'mem_zoef = ', mem
       write(pri,'(a,2i10)')  'max_zoef = ', nn, izoef
       write(pri,'(a,i10/)')  'max_zoef = ', nreloc
      end if

      Call MPI_REDUCE(mem_max_coef,mem,   1,MPI_DOUBLE_PRECISION,MPI_MAX,0,MPI_COMM_WORLD,ierr)
      Call MPI_REDUCE(max_coef,    nn,    1,MPI_integer,MPI_MAX,0,MPI_COMM_WORLD,ierr)
      Call MPI_REDUCE(max_term,    nt,    1,MPI_integer,MPI_MAX,0,MPI_COMM_WORLD,ierr)
      Call MPI_REDUCE(coef_realloc,nreloc,1,MPI_integer,MPI_MAX,0,MPI_COMM_WORLD,ierr)

      if(pri.gt.0) then
       write(pri,'(a,f10.1)') 'mem_coef = ', mem
       write(pri,'(a,2i10)')  'max_coef = ', nn, icoef
       write(pri,'(a,2i10)')  'max_term = ', nt
       write(pri,'(a,i10/)')  'realloc  = ', nreloc
      end if

      Call MPI_REDUCE(mem_max_boef,mem,   1,MPI_DOUBLE_PRECISION,MPI_MAX,0,MPI_COMM_WORLD,ierr)
      Call MPI_REDUCE(mem_max_blk, meb,   1,MPI_DOUBLE_PRECISION,MPI_MAX,0,MPI_COMM_WORLD,ierr)
      Call MPI_REDUCE(max_boef,    nn,    1,MPI_integer,MPI_MAX,0,MPI_COMM_WORLD,ierr)
      Call MPI_REDUCE(max_blk,     nt,    1,MPI_integer,MPI_MAX,0,MPI_COMM_WORLD,ierr)
      Call MPI_REDUCE(boef_realloc,nreloc,1,MPI_integer,MPI_MAX,0,MPI_COMM_WORLD,ierr)
      Call MPI_REDUCE(blk_realloc, mreloc,1,MPI_integer,MPI_MAX,0,MPI_COMM_WORLD,ierr)

      if(pri.gt.0) then
       write(pri,'(a,f10.1)') 'mem_boef = ', mem
       write(pri,'(a,f10.1)') 'mem_blk  = ', meb
       write(pri,'(a,2i10)')  'max_boef = ', nn, iboef
       write(pri,'(a,2i10)')  'max_blk  = ', nt
       write(pri,'(a,i10/)')  'al_boef  = ', nreloc
       write(pri,'(a,i10/)')  'al_blk   = ', mreloc
      end if

      End Subroutine Debug_printing



