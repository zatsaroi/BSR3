!======================================================================
!     PROGRAM       B S R _ C I                           version 3.0
!
!               C O P Y R I G H T -- 2020
!
!     Written by:   Oleg Zatsarinny
!                   free shooter
!                   email: oleg_zoi@yahoo.com
!======================================================================
!     A CONFIGURATION INTERACTION program either non-relativistic
!     or in the BREIT-PAULI approximation  including the cases of
!     NON-ORTHOGONAL radial orbitals and GENERALIZED eigenvalue problem
!!-----------------------------------------------------------------------
!
!    INPUT:  command-line arguments (optional)
!            name.inp   - parameters of calculations (optional)
!            name.c     - list of configurations in the c-file format 
!            name.bsw   - set of radial functions
!            name.bnk   - data bank for angular integrals
!                         int_bnk may be used as alternative
!
!    OUTPUT: name.l     - expansions and energies for LS-solutions
!            name.j     - expansions and energies for LSJ-solutions
!            name.d     - debug information (optional)
!            name.log   - running information
!
!----------------------------------------------------------------------
!
!    The parameters of calculations may be given in command line
!    or (and) in name.inp file as key parameters:  param=value
!
!    First, program reads name.inp, then checks the command-line
!    arguments, so the latter are dominant.
!
!    No interactive mode is applied in the program that simplifies
!    its using in the batch(script) files.
!
!    'name' of case should be given as first argument
!
!----------------------------------------------------------------------
!
!    Examples of calling:
!
!    bsr_ci 2P   -> non-relativistic calculations 
!
!    bsr_ci 2P irel=2 jmin=1 jmax=3  -> relativistic calculations
!
!    bsr_ci 2P irel=1 meiv=5 -> rel.calculations (only 2P.l file),
!                               with restrictions on solutions
!    
!----------------------------------------------------------------------
!    List of possible parameters with default values:
!----------------------------------------------------------------------
!
!    irel=0      - indicate relativistic calculations if >1
!                  1 - only scalar operators; 
!                  2 - + spin-orbit;
!                  3 - + spin-other-orbit;
!                  4 - + spin-spin;
!                  5 - + orbit-orbit;  
!
!    the individual rel.interaction also can be regulated additionally by:              
!
!    iso=-1      - include (=1) or not (=0) spin-orbit interaction
!    isoo=-1     - include (=1) or not (=0) spin-other-orbit int.
!    iss=-1      - include (=1) or not (=0) spin-spin interaction
!    ioo=-1      - include (=1) or not (=0) spin-spin interaction
!
!    jmin=-1     - min. 2J value for rel. calculations
!    jmax=-1     - min. 2J value for rel. calculations
!
!    meiv=0      - how many eigensolution you need (0 -> all)
!    nzero=0     - zero-order dimension (0 -> all configurations)
!    nud=10      - produce debug file name.d if /= 0
!    eps_s=1.d-8 - tolerance for one-electron overlaps 
!    mass=0      - mass polarization (then =1,2)                                                    
!    zmu=0       - nuclear mass (needed to estimate the RMASS)
!
!    In name.inp, name of parameter should begin in the first 
!    column, otherwise it will be disregarded
!
!    all default values are defined in module ci_data
!
!----------------------------------------------------------------------
      Use bsr_ci                                                                                   
      USE term_LS, IL => ILterm, IS => ISterm

      Implicit none
      Integer :: jot
      Real(8) :: time

! ... read 'name.inp' or arguments for parameters of calculations

      Call Read_arg

! ... read configuration list

      Call Read_conf

! ... read orbitals and prepare work arrays

      Call Read_data

! ... calculation of the LS-matrix

      Call matrix_LS

! ... perform calculation for each J-value

      if(jmax.ge.jmin.and.jmin.gt.0) then

       open(nuj,file=AF_j)

       write (nuj,'(2X,A6,A,F5.1,A,I3,A,I7)' ) &
       ATOM,'  Z = ',Z ,'  NEL = ',  NE, '   NCFG = ',NCFG

       Call matrix_rel

       Do jot = jmin, jmax, 2
        CALL matrix_J(jot)
        CALL DIAG(jot,nuj)
       End do

       write(nuj,'(a)') '***'
       Close(nuj)

! ... diagonalization of LS - matrix:

      else 

      open(nul,file=AF_l)
      write (nul,'(2X,A6,A,F5.1,A,I3,A,I7)' ) ATOM,'  Z = ',Z , &
           '  NE = ',  NE, '   NCFG = ',NCFG
      jot=0; CALL DIAG(jot,nul)
      write(nul,'(a)') '***'
      Close(nul)

      end if

! ... timing

      Call CPU_time(time)

      write(iwrite,'(/a,T30,f10.2,a/)') 'Total time:', time/60,'  min'

      End ! Program BSR_CI



