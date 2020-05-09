!======================================================================
!     PROGRAM       B S R _ M A T              version 3
!
!               C O P Y R I G H T -- 2010
!
!     Written by:   Oleg Zatsarinny 
!     email:        oleg_zoi@yahoo.com
!
!======================================================================
!     Generate the interaction matrixes in B-spline representation
!======================================================================
!
!   INPUT FILES: 
!
!     target         -  description of target states and channels
!     target.bsw     -  target w.f.'s in B-spline basis
!     target_orb     -  list of target physical orbitals
!     knot.dat       -  B-spline grid
!     bsr_par        -  input parameters
!     cfg.nnn        -  configuration list for partial wave nnn
!     int_int.nnn    -  angular coefficient data bank (+ int_inf.nnn)
!     pert_nnn.bsw   -  perturber orbitals if any
!
!   OUTPUT FILES:
!
!     bsr_mat.log   -  general running information
!     mat_log.nnn   -  running information for partial wave nnn
!     bsr_mat.nnn   -  resulting overlap and interaction matrixes
!
!=====================================================================
      Use bsr_mat
      Use conf_LS

      Implicit none
      Real(8) :: t1,t2,t3
      Integer :: i

! ... general output:

      Call bsr_mat_inf
      Open(prj,file=AF_prj)       

! ... prepare B-spline parameters:

      Call define_grid(z)
      Call define_spline        

! ... target information:

      Open(nut,file=AF_tar,status='OLD')
      Call R_target(nut)

! ... read input arguments or parameters if any:

      Open(nup,file=AF_par,status='OLD')
      Call Read_arg(nup)

! ... find the common core energy:

      Open(nuw,file=AF_bsw,STATUS='OLD',form='UNFORMATTED')
      Call R_bwfn(nuw)
      Close(nuw)
      if(nwt.ne.nbf) &
       Call Stop_mpi (0,0,' nwt (in target) <> nbf (target.bsw)')
      write(ALSP,'(i3.3)') klsp1  
      i=LEN_TRIM(AF_cfg); ; write(AF_cfg(i-2:i),'(i3.3)') klsp1
      Open(nuc,file=AF_cfg,status='OLD')
      Call R_closed(nuc)
      kclosd=nclosd    
      Call Bcore
      write(prj,'(/a,i4,a)') 'nclosd  =',nclosd,' - common core shells'
      write(prj,'(/a,F15.8,a)') 'Bcore   =', EC,'  -  calculated core energy'
      Call Read_rpar(nup,'Ecore',EC)
      write(prj,'(/a,F15.8,a)') 'Ecore   =', EC,'  -  Used core energy'

      Call CPU_time(t1)
!----------------------------------------------------------------------
! ... loop over partial waves:

      Do klsp = klsp1,klsp2
 
       write(*,'(/a,i3/)') 'BSR_MAT:  klsp =', klsp

       Call CPU_time(t2);   Call SUB1;   Call CPU_time(t3)

       write(pri,'(/a,T20,f8.2,a)') 'Total time:', (t3-t2)/60, ' min'
       write(*  ,'( a,T20,f8.2,a)') 'Total time:', (t3-t2)/60, ' min'

      End do  ! over klsp

      Call CPU_time(t2)
      write(prj,'(/a,T20,f10.2,a)') 'bsr_mat:',(t2-t1)/60,'  min'

      End  ! program bsr_mat


