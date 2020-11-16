!=====================================================================
      Subroutine inform
!=====================================================================
!     print full ci_bnk information in ci_bnk.inf
!---------------------------------------------------------------------
      
      write(*,'(a)') &
'                                                                            ', &
'Call as: bsr_ci  name  par1=value1  par2=value2  ...                        ', & 
'                                                                            ', &
'Examples of calling:                                                        ', &
'                                                                            ', &
'bsr_ci 2P  -> non-relativistic calculations for the case 2P                 ', &
'bsr_ci 2P irel=2 jmin=1 jmax=3  -> relativistic calculations for 2J = 1,3   ', &
'bsr_ci 2P irel=1 msol=5 -> relativistic-only-shift calculations             ', &
'                           with 5 solutions to be output                    ', &
'                                                                        ', &
'------------------------------------------------------------------------', &
'List of parameters with default values:                                 ', &
'------------------------------------------------------------------------', &
'irel=0      - indicate relativistic calculations if >1                  ', &
'              (mrel - alias for irel, as in BSR)                        ', &
'              1 - only scalar operators;                                ', &
'              2 - + spin-orbit;                                         ', &
'              3 - + spin-other-orbit;                                   ', &
'              4 - + spin-spin;                                          ', &
'              5 - + orbit-orbit;                                        ', &
'the individual rel.interaction also can be regulated additionally by:   ', &           
'iso= -1     - include (=1) or not (=0) spin-orbit interaction           ', &
'isoo=-1     - include (=1) or not (=0) spin-other-orbit int.            ', &
'iss= -1     - include (=1) or not (=0) spin-spin interaction            ', &
'ioo= -1     - include (=1) or not (=0) spin-spin interaction            ', &
'jmin=-1     - min. 2J value for rel. calculations                       ', &
'jmax=-1     - min. 2J value for rel. calculations                       ', &
'msol =0      - how many eigensolution you need (0 -> all)               ', &
'----------------------------------------------------------------------  '
      Stop 
                                                                             
      End Subroutine inform                                             
                                                                             
                                                                             
                                                                   