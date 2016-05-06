!======================================================================
!     PROGRAM       B S R _ P O L                          version. 3                        
!
!               C O P Y R I G H T -- 2011
!
!     Written by:   Oleg Zatsarinny
!                   email: oleg_zoi@yahoo.com
!
!======================================================================
!
!    INPUT  ARGUMENTS:
!
!     klsp  - the indexes of partial waves under consideration  
!
!    INPUT FILES:
!
!     bsr_par       -  description of partial waves
!     cfg.nnn       -  c-file for close-coupling expansion
!     bsr_mat.nnn   -  interaction matrix
!     dv.nnn        -  dipole matrix
!
!    OUTPUT FILES:
!   
!     bsr_pol.log   -  running information 
!     pol.nnn       -  bound-like solusions for polarized pseudo-states
!
!     Above, nnn indicates index of the partial wave
!
!=====================================================================
      Use bsr_pol, only: pri
      
      Implicit none
      Real(8) :: t1,t2
      Real(8), external :: RRTC

      Call inf_bsr_pol

      t1=RRTC()

! ... read data:

      Call Read_data
      
! ... read interaction matrix:

      Call R_bsrmat

! ... read transition matrix:

      Call R_dipmat

! ... additional orthogonality

      Call Read_nortb

! ... solve the dipole equation:

      Call Solv_mat

      t2=RRTC()
      write(pri,'(/a, F10.2, a )' ) 'time =',(t2-t1)/60, ' min.'
      write(*  ,'(/a, F10.2, a )' ) 'time =',(t2-t1)/60, ' min.'

      End  ! program BSR_POL


