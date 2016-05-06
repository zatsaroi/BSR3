!======================================================================
      Subroutine inf_bsr_mat
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
'BSR_MAT - sets up of interaction matrix in B-spline representation          ',&
'                                                                            ',&
'INPUT FILES:                                                                ',&
'                                                                            ',&
'knot.dat       -  B-spline grid                                             ',&
'bsr_par        -  input parameters                                          ',&
'target         -  target states and channel information                     ',&
'cfg.nnn        -  configuration list for partial wave  nnn                  ',&
'bnk_int.nnn    -  angular coefficient data bank                             ',&
'target.bsw     -  target w.f. in B-spline basis                             ',&
'pert_nnn.bsw   -  perturb w.f., if any                                      ',&
'                                                                            ',&
'OUTPUT FILES:                                                               ',&
'                                                                            ',&
'bsr_mat.nnn    -  resulting interaction matrix                              ',&
'mat_log.nnn    -  running information                                       ',&
'int_mat.nnn    -  debug output of integrals                                 ',&
'                                                                            ',&
'MAIN PARAMETERS (with default values)                                       ',&
'                                                                            ',&
'klsp1 = 1  -  first partial wave under consideration                        ',&
'klsp2 = 1  -  last partial wave under consideration                         ',&
'klsp  = 1  -  condider this partial wave only                               ',&
'                                                                            ',&
'iitar = 0  -  target states are supposed to be orthogonomal eigenstates     ',&
'              of target Hamiltonian                                         ',&
'              =1 - target state should be orthogonal, but may not be        ',&
'                   eigenstates                                              ',&
'              =2 - target states may be non-orthpgonal (general case)       ',&
'                                                                            ',&
'mrel =  0  -  relativistic corrections:                                     ',&
'                                                                            ',&
'              mrel=1 - include only one-electron scalar rel. integrals      ',&
'              mrel=2 - plus one-electron spin-orbit interaction             ',&
'              mrel=3 - plus two-electron spin-other-orbit interaction       ',&
'              mrel=4 - plus spin-spin interaction                           ',&
'              mrel=5 - plus orbit-orbit interaction                         ',&
'                                                                            ',&
'              each rel. correction can also be controled by parameters      ',&
'              mso, msoo, mss, moo with values -1,0,+1, meaning              ',&
'              +1 - include anyway, -1 - exculed anyway, 0 - follow mrel     ',&
'                                                                            ',&
'imvc = -1  -  mode for processing the mass-velocity term                    ',&
'              = +1 - include mass-velocity term directly in the Hamiltonian ',&
'              =  0 - exclude anyway                                         ',&
'              = -1 - include mass-velocity term only for the bound orbitals ',&
'                     For channel orbitals this interaction can then be      ',&
'                     included in the BSR_HD program as a first-order        ',&
'                     correction.                                            ',&
'                                                                            ',&
'izcorr= 0  - if =1, small-radius cut-off will be applied to the spin-orbut  ',&
'             interaction. Recomended for Z > 40.                            ',&
'                                                                            ',&
'zcorr= 1.0 - semiempirical correction to spin-orbit parameters with l=1     ',&
'                                                                            ',&
'EC = ...   - core energy (if needed to change)                              ',&
'                                                                            ',&
'                                                                            '
      Stop ' '                                                              
                                                                            
      End Subroutine inf_bsr_mat                                          
                                                                            
                                                                            
                                                                            








