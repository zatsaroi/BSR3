!======================================================================
!     PROGRAM       B S R _ H D                      version 4        
!
!               C O P Y R I G H T -- 2011
!
!     Written by:   Oleg Zatsarinny
!                   email: oleg_zoi@yahoo.com
!
!======================================================================
!   Program diagonalizes the Hamiltonian in the inner region and
!   generates H.DAT file for scattering or photoionization 
!   calculations in the R-matrix approach.
!   Additional option - calculation of bound (or pseudo-) states.
!======================================================================
!
!    MAIN INPUT ARGUMENTS:
!
!     itype - type of calculations
!   
!           = -1 - bound-structure calculations
!           =  0 - scattering calculations
!           =  1 - additional output of all R-matrix solutions 
!                  required for photoionization calculations
!
!     klsp  or klsp1,klsp2  - the indexes of partial waves under consideration  
!
! 
!    INPUT FILES:
!
!     bsr_par        -  input parameters
!     target         -  description of target states and partial waves
!     knot.dat       -  description of B-spline grid
!     cfg.nnn        -  c-file for close-coupling expansion
!     bsr_mat.nnn    -  interaction matrix
!
!    OUTPUT FILES:
!   
!     bsr_hd.nnn     -  running information 
!     h.nnn          -  eigenvalues and suface amplitudes
!                       required to obtain R-matrix
!     w.nnn          -  file of channel weights (optional)
!     rsol.nnn       -  R-matrix inner-region solutions (optional)
!     bound.nnn      -  bound-like solutions (old FORMATTED file)
!     ubound.nnn     -  bound-like solutions (new UNFORMATTED file)
!
!     'nnn' indicates the index of the partial wave under consideration
!
!=====================================================================
      Use bsr_hd
      Use spline_atomic, only: z
      
      Implicit none
      Real(8) :: t1,t2

      Call bsr_hd_inf
!---------------------------------------------------------------------
! ... set up B-splines:
 
      CALL define_grid(z); Call define_spline   
      Call Conv_au (Z,0.d0,au_cm,au_eV,0)

! ... read channel information:

      Open(nut,file=AF_tar,status='OLD')
      Call R_target (nut)

! ... read arguments:

      Open(nup,file=AF_par,status='OLD')
      Call R_arg(nup)
      Close(nup)

!----------------------------------------------------------------------
! ... calculations for given partial wave ...

      Do klsp = klsp1,klsp2

       write(*,'(/a,i3,a,i3)') 'BSR_HD: calculations for partial wave: ',klsp

       Call CPU_time(t1);  Call SUB1_HD;  Call CPU_time(t2)

       write(*,'(/a, F10.2, a )' ) 'time =',(t2-t1)/60, ' min.'

       write(pri,'(/a, F10.2, a )' ) 'time =',(t2-t1)/60, ' min.'
       close(pri)

      End do

      End  ! program BSR_HD


