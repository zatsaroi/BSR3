!======================================================================
      Subroutine bsr_hd_inf
!======================================================================
!     provide information about "bsr_hd" program
!----------------------------------------------------------------------
      Character :: A
      Integer :: nu = 99;  Character(80) :: AF = 'bsr_hd_inf'
      iarg = IARGC();     if(iarg.eq.0) Return
      Call GETARG(1,A);   if(A.ne.'?') Return
      open(nu,file=AF) 


      write(nu,'(a)') &
'                                                                                         ',&
'BSR_HD - diagonalization of the interaction matixes bsr_mat.nnn                          ',&
'                                                                                         ',&
'Arguments:   call as:   bsr_hd klsp1=1 klsp2=10 itype=1 ...                              ',&
'                                                                                         ',&
'klsp, klsp1, klsp2 - partial waves under consoderation                                   ',&
'                                                                                         ',&
'itype  = -1  -  bound-state calculations (output - bound.nnn, w.out)                     ',&
'          0  -  scattering calculations  (output - h.nnn), default mode                  ',&
'          1  -  photoionization          (output - h.nnn, rsol.nnn)                      ',&
'iexp   = 0|1 -  theoretical | experimental target energies                               ',&
'                (exp. energies should be given in file "thresholds"	                  ',&
'jmvc   =  0  -  number of channel orbitals for inclusion of mass-velocity correction     ',&
'Edmin  =  0  -  all one-channel solutions with  E < Edmin, E > Edmax or abs(E) < Egap    ',&
'Edmax  =  0  -  will be neglected  (0 means no restrictions)                             ',&
'Egap   = 0.1                                                                             ',&
'Emin   =  0  -  restrictions for output of bound states                                  ',&
'Emax   =  0                                                                              ',&
'msol   =  0  -  max. number of bound states for output                                   ',&
'it_max =  0  -  max. target threshold for output                                         ',&
'iwt    = 0|1 -  output of w.nnn files with channels weights                              ',&
'cwt    =  0  -  additional output of channel composition for weights > cwt (0.01, e.g)   ',&
' '                                                                                      

      write(*,'(a)') &
'                                                                                         ',&
'BSR_HD - diagonalization of the interaction matixes bsr_mat.nnn                          ',&
'                                                                                         ',&
'Call as:        bsr_hd klsp1=1 klsp2=10 itype=1 ...                                      ',&
'                                                                                         ',&
'klsp, klsp1, klsp2 - partial waves under consoderation                                   ',&
'                                                                                         ',&
'itype  = -1  -  bound-state calculations (output - bound.nnn, w.out)                     ',&
'          0  -  scattering calculations  (output - h.nnn), default mode                  ',&
'          1  -  photoionization          (output - h.nnn, rsol.nnn)                      ',&
'iexp   = 0|1 -  theoretical | experimental target energies                               ',&
'                (exp. energies should be given in file "thresholds"	                  ',&
'jmvc   =  0  -  number of channel orbitals for inclusion of mass-velocity correction     ',&
'Edmin  =  0  -  all one-channel solutions with  E < Edmin, E > Edmax or abs(E) < Egap    ',&
'Edmax  =  0  -  will be neglected  (0 means no restrictions)                             ',&
'Egap   = 0.1                                                                             ',&
'Emin   =  0  -  restrictions for output of bound states                                  ',&
'Emax   =  0                                                                              ',&
'msol   =  0  -  max. number of bound states for output                                   ',&
'it_max =  0  -  max. target threshold for output                                         ',&
'iwt    = 0|1 -  output of w.nnn files with channels weights                              ',&
'cwt    =  0  -  additional output of channel composition for weights > cwt (0.01, e.g)   ',&
' '                                                                                      
      Stop                                                                               
                                                                                         
      End Subroutine bsr_hd_inf                                                          
                                                                                         
                                                                                         
                                                                                         
                                                                                         
                                                                                         
                                                                                         
                                                                                         
                                                                                         
                                                                                         
                                                                                         
                                                                                         
                                                                                         
                                                                                         