!=====================================================================
      Subroutine bsr_conf_inf
!=====================================================================
!     print information for bsr_prep program
!---------------------------------------------------------------------
      Implicit none
      Integer :: nu=6
      Character ::  name = ' '
      
      Call Read_name(name)
      if(name.ne.'?') Return
      write(*,'(a)') &
'                                                                          ', &
'bsr_conf prepares configuration expansions for all partial waves (cfg.nnn)', &
'and assigns orthogonality conditions if any                               ', &
'                                                                          ', &
'all input files are created by bsr_prep program, except file "target"     ', &
'which can be modified with additional list of partial waves;              ', &
'                                                                          ', &
'additional file "target_del"  can be used for deleting specific channels  ', &
'                                                                          ', &
'OPTIONAL ARGUMENTS:  (-1 value means out of play)                         ', &
'                                                                          ', &
'c_comp  [1.01]  - tolerance for compensation configuration                ', &
'                  to put orth. condition                                  ', &
'max_ll  [-1]    - restiction on max. small l                              ', &
'min_ll  [-1]    - restiction on min. small l                              ', &
'max_LT  [-1]    - restiction on total L as (2L+1)                         ', &
'max_ST  [-1]    - restiction on total S as (2S+1)                         ', &
'max_it  [-1]    - restiction on number of target states                   ', &
'kort    [-1]    - if > 0, the orthogonal conditions from cfg-files        ', &
'                  are also in play                                        '
      Stop 
                                                                             
      End Subroutine bsr_conf_inf                                             
