!======================================================================
!     PROGRAM       B S R _ M A T _ M P I            version.4
!
!               C O P Y R I G H T -- 2016
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
!     target         -  description of target states and channels
!     target.bsw     -  target w.f.'s in B-spline basis
!     target_orb     -  list of target physical orbitals
!     knot.dat       -  B-spline grid
!     bsr_par        -  input parameters
!     knot.dat       -  B-spline grid
!     bsr_par        -  description of target states and channels
!     cfg.nnn        -  configuration list for given partial wave
!     int_int.nnn    -  angular coefficient data bank (+ int_inf.nnn)
!     pert_nnn.bsw   -  perturb w.f., if any
!
!   OUTPUT FILES:
!
!     bsr_mat.log   -  general running information
!     mat_log.nnn    -  running information
!     bsr_mat.nnn    -  resulting interaction matrix
!
!=====================================================================
      Use MPI

      Use bsr_mat
      Use conf_LS

      Implicit none
      Real(8) :: t1,t2,t3
      Integer :: i

      Call MPI_INIT(ierr)
      Call MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierr)
      call MPI_COMM_SIZE( MPI_COMM_WORLD, nprocs, ierr )

      if(myid.eq.0) write(*,'(a,i5)') 'MPI: nprocs = ',nprocs

      time0 = MPI_WTIME()

! ... general output:

      if(myid.eq.0) Call bsr_mat_inf

      if(myid.eq.0) then; Open(prj,file=AF_prj); else; prj=0; end if

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

! ... define arguments:

      if(myid.eq.0) then
       Open(nup,file=AF_par,status='OLD')
       Call Read_arg(nup)
      end if
      Call br_arg

! ... find nclosd and core energy:

      if(myid.eq.0) then
       Open(nuw,file=AF_bsw,STATUS='OLD',form='UNFORMATTED')
       Call R_bwfn(nuw)
       Close(nuw)
       if(nwt.ne.nbf) &
        Call Stop_mpi (0,0,' nwt (in target) <> nbf (target.bsw)')
       i=LEN_TRIM(AF_cfg); ; write(AF_cfg(i-2:i),'(i3.3)') klsp1
       Open(nuc,file=AF_cfg,status='OLD')
       Call R_closed(nuc)
       kclosd=nclosd    
       Call Bcore
       write(prj,'(/a,i4,a)') 'nclosd  =',nclosd,' - common core shells'
       write(prj,'(/a,F15.8,a)') 'Bcore   =', EC,'  -  calculated core energy'
       Call Read_rpar(nup,'Ecore',EC)
       write(prj,'(/a,F15.8,a)') 'Ecore   =', EC,'  -  Used core energy'
      end if
      Call br_core

! ... loop over partial waves:

      Do klsp = klsp1,klsp2
 
       if(myid.eq.0)  write(*,'(/a,i3/)') 'BSR_MAT:  klsp =', klsp

       t1 =  MPI_WTIME();   Call SUB1;   t2 =  MPI_WTIME()

       if(pri.gt.0) &
        write(pri,'(a,5x,T20,f10.1,a)') 'Total time:', (t2-t1)/60, ' min'
       if(myid.eq.0) &
        write(*  ,'(a,5x,T20,f10.1,a)') 'Total time:', (t2-t1)/60, ' min'

      End do  ! over klsp

      Call MPI_FINALIZE(ierr)

      End  ! program bsr_mat


