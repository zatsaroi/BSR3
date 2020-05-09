!======================================================================
!     PROGRAM       B S R _ M A T _ M P I           
!
!               C O P Y R I G H T -- 2011
!
!     Written by:   Oleg Zatsarinny 
!     email:        oleg_zoi@yahoo.com
!
!======================================================================
!     Generation of interaction matrixes in B-spline representation
!======================================================================
!
!   INPUT ARGUMENTS:
!
!     klsp1, klsp2   -  range of partial wave under consideration
!
!   INPUT FILES: 
!
!     knot.dat       -  B-spline grid
!     bsr_par        -  description of target states and channels
!     cfg.nnn        -  configuration list for given partial wave
!     bnk_int.nnn    -  angular coefficient data bank
!     target.bsw     -  target w.f.'s in B-spline basis
!     pert_nnn.bsw   -  perturb w.f., if any
!
!   OUTPUT FILES:
!
!     mat_log.nnn    -  running information
!     bsr_mat.nnn    -  resulting interaction matrix
!     int_mat.nnn    -  debuging information
!
!=====================================================================

      USE MPI

      USE bsr_mat
      USE conf_LS
      USE spline_atomic
      USE orb_LS,         only: nwf

      Implicit none

      Real(8) :: t1,t2,t3
      Integer :: i,nprocs

      Call MPI_INIT(ierr)
      Call MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierr)
      call MPI_COMM_SIZE( MPI_COMM_WORLD, nprocs, ierr )

      if(myid.eq.0) write(*,*) 'nprocs = ',nprocs

! ... general output:

      if(myid.eq.0) then; Open(prj,file=AF_prj); else; prj=0; end if

! ... define arguments:

      if(myid.eq.0) then
       Open(nup,file=AF_par,status='OLD')
       Call Read_arg(nup)
      end if
      Call br_arg

! ... prepare B-spline parameters:

      if(myid.eq.0) Call define_grid(z)
      Call br_grid
      Call define_spline        

! ... target:

      if(myid.eq.0) then
       Open(nut,file=AF_tar,status='OLD')
       Call R_target(nut)
      end if
      Call br_target

! ... find nclosd and core energy:

      if(myid.eq.0) then
       i=LEN_TRIM(AF_cfg); write(AF_cfg(i-2:i),'(i3.3)') klsp1
       Open(nuc,file=AF_cfg,status='OLD')
       Call R_closed(nuc)
       Close(nuc)
       kclosd=nclosd    
       Call Allocate_bsorb(nwf)
       Open(nuw,file=AF_bsw,STATUS='OLD',form='UNFORMATTED')
       Call Read_bsw(nuw)
       write(prj,'(/a,i5)') 'nclosd =',nclosd
       Call Bcore
       write(prj,'(/a,F16.8)') 'EC    =', EC
       write(*,'(/a,F16.8)') 'EC    =', EC
       Call Read_rpar(nup,'Ecore',EC)
       write(prj,'(/a,F16.8)') 'Ecore =', EC
      end if
      Call br_core

      t1 =  MPI_WTIME()

! ... loop over partial waves:

      Do klsp = klsp1,klsp2
 
       if(myid.eq.0)  write(*,'(/a,i3/)') 'BSR_MAT:  klsp =', klsp

       t2 =  MPI_WTIME();   Call SUB1;   t3 =  MPI_WTIME()

       if(pri.gt.0) &
        write(pri,'(a,5x,f8.2,a)') 'Total time:', (t3-t2)/60, ' min'
       if(myid.eq.0) &
        write(*  ,'(a,5x,f8.2,a)') 'Total time:', (t3-t2)/60, ' min'

      End do  ! over klsp

      t2 =  MPI_WTIME()
      if(pri.gt.0) &
       write(prj,'(/a,f10.2,a)') 'Total:  ',(t2-t1)/60,'  min'
      if(myid.eq.0) &
       write(*  ,'(/a,f10.2,a)') 'Total:    ',(t2-t1)/60,'  min'

      Call MPI_FINALIZE(ierr)

      End  ! program bsr_mat


