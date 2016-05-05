!======================================================================
      Subroutine inf_bsr_dmat
!======================================================================
!     provide screen information about bsr_dmat program
!----------------------------------------------------------------------
      Implicit none
      Character :: A
      Integer :: iarg, IARGC
      Integer :: nu = 99;  Character(80) :: AF = 'bsr_dmat_inf'

      iarg = IARGC();      if(iarg.eq.0) Return
      Call GETARG(1,A);    if(A.ne.'?') Return

      open(nu,file=AF)
      write(nu,'(a)') &
'BSR_DMAT provides different multipole-transition matrixes or          ',&
'oscillator strengths and transition matrix elements                   ',&
'                                                                      ',&
'Multipole index is taken from mult_bnk, so run program MULT first     ',&
'                                                                      ',&
'INPUT ARGUMENTS:                                                      ',&
'                                                                      ',&
'four mandatory position parameters:                                   ',&
'                                                                      ',&
'1.  name1.c or cfg.nnn - c-file for initial state                     ',&
'2.  name2.c or cfg.nnn - c-file for final state                       ',&
'3.  structure mode (ctype) for the initial state  (c,l,j,b)           ',&
'4.  structure mode (ctype) for the final state    (c,l,j,b,p,q,d)     ',&
'                                                                      ',&
'additonal optional key paramters (as par=value):                      ',&
'                                                                      ',&
'istate1 [0] - index of state from 1st set of solutions                ',&
'              to be considered, 0 means all solutions                 ',&
'istate2 [0] - index of state from 2nd set of solutions                ',&
'              to be considered, 0 means all solutions                 ',&
'mstate1 [0] - restriction on the number of 1st-set solutions          ',&
'              to be considered, 0 means all solutions                 ',&
'mstate2 [0] - restriction on the number of 2st-set solutions          ',&
'              to be considered, 0 means all solutions                 ',&
'gf      [f] - output f- of gf-values                                  ',&
'ialfa   [0] - additional output radial matrix elements and            ',&
'              polarizabilities (in zf_res file)                       ',&
'                                                                      ',&
'Some name conventions are taken from MCHF package of CFF (1995 ver.)  ',&
'                                                                      ',&
'ctype = c  -  expantion coefficients are read from c-file (mcfh)      ',&
'ctype = l  -  expantion coefficients are read from l-file (CI)        ',&
'ctype = j  -  expantion coefficients are read from j-file (CI)        ',&
'ctype = b  -  expantion coefficients are read from bound.nnn file     ',&
'              after BSR bound_type calculations                       ',&
'ctype = p  -  2nd set of solutions are supposed to be R-matrix        ',&
'              solutions in rsol.nnn file after BSR calculations       ',&
'              (this option is used in photoionization calculations    ',&
'               see program BSR_PHOT)                                  ',&
'ctype = q  -  cfg.nnn is expected for 2nd set but no expansion        ',&
'              coefficients are required (this option is used in the   ',&
'              polarized psedo-state calculation, see program BSR_POL) ',&
'ctype = d  -  cfg.nnn and cfg.mmm are expected for both 1st and       ',&
'              2nd set but no expansion coefficients are required      ',&
'              (used in the strong-field calculations)                 ',& 
'                                                                      '
      write(nu,'(a)') &
'INPUT FILES:                                                          ',&
'                                                                      ',&
'mult_bnk  - data bank of ang.coefficients for multipole operator      ',&
'            (obtained with program MULT)                              ',&
'                                                                      ',&
'knot.dat  - B-spline information                                      ',&
'                                                                      ',&
'name1.c   - c-file for initial state                                  ',&
'name1.bsw - initial state wavefunctions in B-spline representation    ',&
'                                                                      ',&
'name2.c   - c-file for final state                                    ',&
'name2.bsw - initial state wavefunctions in B-spline representation    ',&
'                                                                      ',&
'cfg.nnn   - One (or both) c-files can be replaced by cfg.nnn files,   ',&
'            where nnn is a partial wave index (for ctype=b).          ',&
'            The close-coupling expansion is supposed for this case,   ',&
'            that requires also the following files (instead name.bsw) ',&
'                                                                      ',&
'target        - information about close-coupling expansion            ',&
'                (target states, channels and perturbers)              ',&
'target.bsw    - target wavefunctions in B-spline representation       ',&
'pert_nnn.bsw  - perturber w.f, if any                                 ',&
'                                                                      ',&
'bound.nnn     - expansions coefficients if ctype=b                    ',&
'rsol.nnn      - R-matrix solutions if ctype=p                         ',&
'                                                                      ',&
'OUTPUT FILES:                                                         ',&
'                                                                      ',&
'bsr_dmat.log - running information                                    ',&
'                                                                      ',&
'zf_res    - oscilator strengths ( ctype2 = c,l,j,b ), appended(!)     ',&
'd.nnn     - dipole matrix for R-matrix solutions   (ctype2 = p)       ',&
'dv.nnn    - dipole vector for polarizability calc. (ctype2 = q)       ',&
'd.nnn_mmm - clear dipole matrix between B-splines  (ctype2 = d)       ',&
'                                                                      ',&
'EXAMPLES FOR CALLING:                                                 ',&
'                                                                      ',&
'bsr_dmat 1.c 2.c c c                                                  ',&
'bsr_dmat 1.c 2.c j j  gf=g  ialfa=1                                   ',&
'bsr_dmat 1.c cfg.001 c b  mstate2=5                                   ',&
'bsr_dmat 1.c cfg.002 c p                                              ',&
'bsr_dmat cfg.001 cfg.002 b b istate1=2 istate2=5                      ',&
'bsr_dmat 1.c cfg.002 c p                                              ',&
'bsr_dmat 1.c cfg.002 c q                                              ',&
'bsr_dmat cfg.001.c cfg.002 b d                                        ',&
'                                                                      ',&
'RESTRICTIONS:  in b->j or j->b calculations, j-file may contain       ',&
'               data only for one value of J.                          ',&
'                                                                      ',&
'For automatic calculations between many states, see also utils        ',&
'                                                                      ',&
'          zf_cc_bsr   and   zf_bb_bsr                                 ',&
'                                                                      '

      write(*,'(a)') &
'                                                                       ',&
'BSR_DMAT - calculates oscillator strengths or dipole matrix element    ',&
'                                                                       ',&
'INPUT FILES depends on the mode of operation, however, the "mult_bnk"  ',&
'            file after program MULT is always expected                 ',&
'                                                                       ',&
'EXAMPLES FOR CALLING:                                                  ',&
'                                                                       ',&
'bsr_dmat 1.c 2.c c c                                                   ',&
'bsr_dmat 1.c 2.c j j  gf=g  ialfa=1                                    ',&
'bsr_dmat 1.c cfg.001 c b  mstate2=5                                    ',&
'bsr_dmat 1.c cfg.002 c p                                               ',&
'bsr_dmat cfg.001 cfg.002 b b istate1=2 istate2=5                       ',&
'                                                                       ',&
'for more detailed description of optional input arguments,             ',&
'see file "bsr_dmat_inf", created after this call                       ',&
'                                                                       ',&
'For automatic calculations between many states, see also utils         ',&
'                                                                       ',&
'          zf_cc_bsr   and   zf_bb_bsr                                  ',&
'                                                                       '

      Stop                                                               
                                                                            
      End Subroutine inf_bsr_dmat                                          
                                                                            
                                                                            
                                                                            








