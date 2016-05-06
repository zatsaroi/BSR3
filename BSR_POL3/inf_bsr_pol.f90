!======================================================================
      Subroutine inf_bsr_pol
!======================================================================
!     provide screen information about bsr_mat program
!----------------------------------------------------------------------
      Character(80) :: A

      iarg = IARGC()
      if(iarg.eq.0) Return
      Call GETARG(1,A)      
      if(A.ne.'?') Return

      write(*,*) &
'                                                                            ',&
'BSR_POL - calculation of polarized psedo-states:                            ',&
'                                                                            ',&
'INPUT FILES:                                                                ',&
'                                                                            ',&
'knot.dat       -  B-spline grid                                             ',&
'target         -  target states and channel information                     ',&
'cfg.nnn        -  configuration list for psedo-states in partial wave nnn   ',&
'bsr_mat.nnn    -  overlap and interaction matrixes for partial wave nnn     ',&
'dv.nnn         -  dipole (multipole) vector for given polarizedpseudo-state ',&
'                  (dv.nnn is preparing by running bsr_dmat with ctype2=q)   ',&
'bound.nnn      -  bound states for orth.constraits if any                   ',&
'pert_nnn.bsw   -  perturb w.f., if any                                      ',&
'target.bsw     -  target w.f. in B-spline basis                             ',&
'                                                                            ',&
'OUTPUT FILES:                                                               ',&
'                                                                            ',&
'pol.nnn        -  pseudo-state expansion (like bound.nnn)                   ',&
'bsr_pol.nnn    -  running information (and polarizability)                  ',&
'                                                                            ',&
'MAIN PARAMETERS (with default values)                                       ',&
'                                                                            ',&
'klsp  = 1  -  partial wave index (nnn)                                      ',&
'nortb = 0  -  orth.constraints to bound states                              ',&
'iortb      -  indexes the above state in bound.nnn if any                   ',&
'              =1 - target state should be orthogonal, but may not be        ',&
'ilzero =1  -  initial zeros in B-spline expansion                           ',&
'ibzero =1  -  final (border) zeros in B-spline expansion                    ',&
'                                                                            ',&
'Example:      bsr_pol klsp=4  nortb=2 iortb=2,3                             ',&
'                                                                            '
      Stop ' '                                                              
                                                                            
      End Subroutine inf_bsr_pol                                          
                                                                            
                                                                            
                                                                            








