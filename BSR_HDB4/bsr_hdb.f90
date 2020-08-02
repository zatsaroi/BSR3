!======================================================================
!    PROGRAM       B S R _ H D 4     (scalapack verson)  
!
!               C O P Y R I G H T -- 2011
!
!    Written by:   Oleg Zatsarinny
!                  email: oleg_zoi@yahoo.com
!
!======================================================================
!    Diagonalize the Hamiltonian in the innere region and
!    generates H.DAT fail prir scattering or photoionization 
!    calculations in R-matrix approach.
!    Another option - calculation of bound (or pseudo-) states.
!======================================================================
!
!    INPUT  ARGUMENTS:
!
!     itype - type of calculations
!   
!           = -1 - bound-structure calculations
!           =  0 - scattering calculations
!           =  1 - additional output of all R-matrix solutions 
!                  required prir photoionization calculations
!
!     klsp, klsp1,klsp2  - the indexes of partial waves under consideration  
!
! 
!    INPUT FILES:
!
!     bsr_par        -  input parameters
!     target         -  description of target states and partial waves
!     knot.dat       -  description of B-spline grid
!     cfg.nnn        -  c-file prir close-coupling expansion
!     bsr_mat.nnn    -  interaction matrix
!
!    OUTPUT FILES:
!   
!     bsr_hd.log     -  running inprirmation 
!     h.nnn          -  eigenvalues and suface amplitudes
!                       required to obtain R-matrix
!     w.nnn          -  file of channel weights (optional)
!     rsol.nnn       -  R-matrix inner-region solutions (optional)
!     bound.nnn      -  bound-like solusions
!
!     'nnn' indicates the index of the partial wave under consideration
!
!=====================================================================

      Use blacs
      Use bsr_hd
     
      Implicit none
      Integer :: err
      Integer, external :: Icheck_file
      Real(8) :: tt1,tt2

!----------------------------------------------------------------------
! ... define the blacs grid:

      Call DEF_BLACS_GRID

!----------------------------------------------------------------------
! ... read common parameters and data:

      err = 0
      if(io_processor) then

! ... set up B-splines:

      fail=Icheck_file(AF_knot) 
      if(fail.eq.0) then
       write(*,'(a)') 'Can not find file knot.dat';  err = err + 1
      else
       CALL def_BS   
       Call Conv_au (100.d0,0.d0,au_cm,au_eV,0)       
       write(*,*) 'au_cm = ',au_cm, '     au_eV = ', au_eV
      end if

! ... read target information:

      fail=Icheck_file(AF_tar) 
      if(fail.eq.0) then
       write(*,'(a)') 'Can not find file target'; err = err + 1
      else
       Open(nut,file=AF_tar); Call R_target (nut)
      end if

! ... read arguments:

      fail=Icheck_file(AF_par) 
      if(fail.eq.0) then
       write(*,'(a)') 'Can not find file bsr_par, defaults are used'
       Call Read_arg(0)
      else
       Open(nup,file=AF_par); Call Read_arg(nup);  Close(nup)    
      end if

      end if ! io_processor

      Call br_ipar(err)
      if(err.ne.0) then; Call BLACS_EXIT(); Stop ' '; end if

! ... broadcast common parameters: 

      Call br_arg

!----------------------------------------------------------------------
! ... separate calculations for each partial wave:

      if (myrow >= 0 .and. myrow < p .and. mycol >= 0 .and. mycol < q) then
 
       Do klsp = klsp1,klsp2

        Call cpu_time (tt1)
        fail = 0
        Call BLACS_BARRIER (ctxt, 'All')

         Call SUB1_HD

        Call BLACS_BARRIER (ctxt, 'All')

        if(io_processor) then           
         Call CPU_time(tt2)
         write (pri,'(/a,T30,f10.2,a)') 'Partial wave:,', (tt2-tt1)/60, ' min.'
         write (*  ,'(/a,T30,f10.2,a)') 'Partial wave:,', (tt2-tt1)/60, ' min.'
        end if

       End do     

      end if

      Call BLACS_BARRIER (ctxt, 'All')

      Call BLACS_GRIDEXIT (ctxt)    

      End  ! program BSR_HDB




!=====================================================================
      Subroutine Def_blacs_grid
!=====================================================================

      Use blacs
      Use bsr_hd, only: nup,AF_par
      
      Implicit none
      
      integer :: imycol, imyrow, ip, iq, par(3)
      integer, external :: Icheck_file

!----------------------------------------------------------------------
! ... define liniar blacs grid to include all processes to broadcast:

      Call BLACS_PINFO (iam, nprocs) ! find process #, total # processors
      Call BLACS_GET (-1, 0, ictxt)  ! find default context, ictxt

      Call BLACS_GRIDINIT (ictxt, 'Row-major', 1, nprocs)
      Call BLACS_GRIDINFO (ictxt, ip, iq, imyrow, imycol)

      io_processor = (imycol == 0)

      if (io_processor) then 

       write (*,*)
       write (*,'(a,i5)') 'BSR_HDB: number of processors = ', nprocs

! ... read p,q,nblocks if any:

      if(Icheck_file(AF_par).ne.0) then
       Open(nup,file=AF_par)
       Call Read_ipar(nup,'p'     ,p     )
       Call Read_ipar(nup,'q'     ,q     )
       Call Read_ipar(nup,'nblock',nblock)
       Close(nup)    
      end if
      Call Read_iarg('p'     ,p     )
      Call Read_iarg('q'     ,q     )
      Call Read_iarg('nblock',nblock)

! ... check the grid parameters:

      if(p*q.eq.0) then; p=SQRT(REAL(nprocs)); q=p; end if

      write (*,*)
      write (*,'(a,3i5)') 'Blacks dimensions, p,q,nblock = ', &
                           p,q,nblock

      end if  ! over io_processor

      if(io_processor) then
       par(1)=p; par(2)=q; par(3)=nblock
       Call igebs2d (ictxt, 'all', ' ', 1, 3, par, 1)
      else      
       Call igebr2d (ictxt, 'all', ' ', 1, 3, par, 1, 0, 0)
       p=par(1); q=par(2); nblock=par(3)
      end if

      if(p*q > nprocs) then; Call BLACS_EXIT(); Stop ' '; end if

      Call BLACS_BARRIER (ictxt, 'All')
      Call BLACS_GRIDEXIT (ictxt)           ! kill initial mesh

!----------------------------------------------------------------------
! ... create p-q BLACS grid:

      Call BLACS_GET (-1,0,ctxt)    
      Call BLACS_GRIDINIT (ctxt, 'Row-major', p, q)
      Call BLACS_GRIDINFO (ctxt, ip, iq, myrow, mycol)

      io_processor = (myrow == 0 .and. mycol == 0)

      End Subroutine Def_blacs_grid


