!======================================================================
      Subroutine inf_mult
!======================================================================
!     provide screen information about program  "mult"
!----------------------------------------------------------------------
      Implicit none
      Character :: A = ' '
      Integer :: nu = 99;  Character(80) :: AF = 'mult_inf'

      Call GET_COMMAND_ARGUMENT(1,A);   if(A.ne.'?') Return

      open(nu,file=AF)
      write(nu,'(a)') &

'MULT - calculates  angular coefficients for different multipole       ',&
'       operators such as E1, M1, E2, M2, ...                          ',&
'                                                                      ',&
'The configuration expansions  for initial and final states are given  ',&
'in c-files with CFF spectroscopic format                              ',&
'                                                                      ',&
'INPUT ARGUMENTS depend on the mode of the case:                       ',&
'                                                                      ',&
'cc mode (default):   calculations between two given states            ',&
'-------                                                               ',&
'Call as:   mult  name1  name2  [ AA   AF_bnk ]                        ',&
'                                                                      ',&
'name1    - c-file for initial state                                   ',&
'name2    - c-file for final state                                     ',&
'AA       - multipole transition index (E1,E2,M1,..., default E1)      ',&
'AF_bnk   - file for results, default - mult_bnk                       ',&
'                                                                      ',&
'ccc mode:  calculations between set of c-files                        ',&
'--------                                                              ',&
'Call as:   mult  ccc  list_c=... AA=..  AF_bnk=..                     ',&
'                                                                      ',&
'list_c   - file  with  list of involved c-files  (one for line)       ',&
'                                                                      ',&
'bsr mode:  calculations between set of BSR cfg.nnn files              ',&
'--------                                                              ',&
'Call as:   mult  bsr  klsp1=.. klsp2=..  AA=..  AF_bnk=..             ',&
'                                                                      ',&
'EXAMPLES FOR CALLING:                                                 ',&
'                                                                      ',&
'mult 1.c 2.c E1                                                       ',&
'mult ccc list_c=name                                                  ',&
'mult bsr klsp1=1 klsp2=5  AF_b=mult_bsr                               ',&
'                                                                      ',&
'RESTRICTION: ONLY ONE TYPE OF TRANSITION IN ONE RUN (!)               ',&
'                                                                      ',&
'WARNING: program moves temporary mult_res file in the final one,      ',&
'         and user should indicate "mv" or "move" command              ',&
'         at compilation - not convenient, but I have no other solution yet ...',&
'                                                                      '

      write(*,'(a)') &
'                                                                      ',&
'MULT - calculates  angular coefficients for different multipole       ',&
'       operators such as E1, M1, E2, M2, ...                          ',&
'                                                                      ',&
'The configuration expansions  for initial and final states are given  ',&
'in c-files with CFF spectroscopic format                              ',&
'                                                                      ',&
'INPUT ARGUMENTS depend on the mode of the case:                       ',&
'                                                                      ',&
'In the main cc mode  call as:    mult  name1  name2  [ AA   AF_bnk ]  ',&
'                                                                      ',&
'name1    - c-file for initial state                                   ',&
'name2    - c-file for final state                                     ',&
'AA       - multipole transition index (E1,E2,M1,..., default E1)      ',&
'AF_bnk   - file for results, default - mult_bnk                       ',&
'                                                                      ',&
'for more detailed description of optional input arguments,            ',&
'see file "mult_inf", created after this call                          ',&
'                                                                      '
      Stop                                                               
                                                                            
      End Subroutine inf_mult                                          
                                                                            
                                                                            
                                                                            








