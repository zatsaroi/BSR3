!=====================================================================
      Subroutine bsr_prep_inf
!=====================================================================
!     print information for bsr_prep program
!---------------------------------------------------------------------
      Implicit none
      Character ::  name = ' '
      
      Call Read_name(name)
      if(name.ne.'?') Return

      write(*,'(a)') &
'                                                                            ', &
'  bsr_prep analizes the input target states and prepare them for BSR:       ', &
'                                                                            ', &
'  INPUT FILES:                                                              ', &
'                                                                            ', &
'     bsr_par   -  parameters of calculation (optional)                      ', &
'     target    -  list of target states and partial waves                   ', &
'     c- and bsw-files for each target state and perturbers if any           ', &
'     knot.dat  -  B-spline parameters                                       ', &
'     target_sub.bsw  - substitution orbitals (optional)                     ', &
'                                                                            ', &
'  OPTIONAL ARGUMENTS:                                                       ', &
'                                                                            ', &
'     eps_ovl  [1.d-6] - tolerance for overlaps                              ', &
'     eps_phys [0.5  ] - minimum occupation number for orbital to be physical', &
'     eps_sub  [0.5  ] - tolerance for substitution orbitals                 ', &
'     eps_targ [2.d-8] - tolerance for configuration coefficient             ', &
'     ii_sub   [0    ] - if > 0, prevent creation new sub. oritals           ', &
'                                                                            ', &
'  additional arguments LT_min,LT_max, IS_min,IS_max, or JJ_min,JJ_max       ', &
'  can be used at nlsp=0 for generation of possible partial waves            ', &
'  (with S -> 2S+1  and J -> 2J representation)                              ', &
'                                                                       '   
      Stop ' '
                                                                             
      End Subroutine bsr_prep_inf                                             

