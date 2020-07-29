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
!     bnk_int.nnn    -  angular-coefficients data bank
!     pert_nnn.bsw   -  perturber orbitals if any
!
!   OUTPUT FILES:
!
!     bsr_mat.log   -  general running information
!     mat_log.nnn   -  running information for partial wave nnn
!     bsr_mat.nnn   -  resulting overlap and interaction matrixes
!
!=====================================================================
      Use bsr_mat; Use target;  Use conf_LS
      Use spline_atomic; Use spline_orbitals

      Implicit none
      Real(8) :: t1,t2,t3
      Integer :: i

      Call bsr_mat_inf
      Open(prj,file=AF_prj)       

! ... prepare B-spline parameters:

      Call define_grid(z)
      Call define_spline        

! ... target information:

      Open(nut,file=AF_tar,status='OLD')
      Call R_target(nut)

      Open(nuw,file=AF_bsw,STATUS='OLD',form='UNFORMATTED')
      Call R_bwfn(nuw)
      Close(nuw)
      if(nwt.ne.nbf) Stop ' nwt (in target) <> nbf (target.bsw)'

! ... read input arguments or parameters if any:

      Open(nup,file=AF_par,status='OLD')
      Call Read_arg(nup)

! ... find the common core energy:

      write(ALSP,'(i3.3)') klsp1  
      i=LEN_TRIM(AF_cfg); AF_cfg(i-2:i)=ALSP
      Open(nuc,file=AF_cfg,status='OLD')
      Call R_closed(nuc)
      kclosd=nclosd    
      Call Bcore

      write(prj,'(/a,i4,a)') 'nclosd  =',nclosd, &
      ' - common core shells'
      write(prj,'(/a,F15.8,a)') 'Bcore   =', EC, &
       '  -  calculated core energy'
      Call Read_rpar(nup,'Ecore',EC)
      write(prj,'(/a,F15.8,a)') 'Ecore   =', EC, &
       '  -  Used core energy'

      Call CPU_time(t1)
!----------------------------------------------------------------------
! ... loop over partial waves:

      Do klsp = klsp1,klsp2
 
       write(*,'(/a,i3/)') 'BSR_MAT3:  klsp =', klsp

       Call CPU_time(t2);  Call SUB1;  Call CPU_time(t3)

       write(pri,'(/a,5x,f8.2,a)') 'Total time:', (t3-t2)/60, ' min'
       write(*  ,'( a,5x,f8.2,a)') 'Total time:', (t3-t2)/60, ' min'

      End do  ! over klsp

      Call CPU_time(t2)
      write(prj,'(/a,f10.2,a)') 'bsr_mat3:  ',(t2-t1)/60,'  min'

      End  ! program bsr_mat


