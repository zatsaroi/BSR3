!=====================================================================
!     PROGRAM   B S R _ B R E I T                            version 4
!
!               C O P Y R I G H T -- 2018
!
!     Written by:   Oleg Zatsarinny 
!                   email: oleg_zoi@yahoo.com
!======================================================================
!     RESULTS: angular coefficient in the non-orthogonal mode 
!----------------------------------------------------------------------
!
!     INPUT ARGUMENTS:
!     
!     klsp,klsp1,klsp2  - range of partial wave in BSR calculations,
!                         with cfg.001, cfg.002, ..., as input files
!                         (default -> klsp=0, with cfg.inp as input file)
!     
!     oper  - character(7), where each position can be 0 or 1,
!             and indicate the operator under consideration:
!             oper(1) - OVERLAPS       
!             oper(2) - KINATIC ENERGY
!             oper(3) - TWO-ELECTRON ELECTROSTATIC
!             oper(4) - SPIN ORBIT       
!             oper(5) - SPIN-OTHER-ORBIT
!             oper(6) - SPIN-SPIN
!             oper(7) - ORBIT-ORBIT       
!             Default -> 1110000 - non-relativistic calculations
!     
!      mk    - max.multipole index (default -> 9)
!     
!----------------------------------------------------------------------
!
!     example:    1.  bsr_breit 
!                 2.  bsr_breit klsp1=1 klsp2=5 oper=1111110  
!                 3.  bsr_breit km=5
!            
!----------------------------------------------------------------------
!
!     INPUT FILES:
!     
!     cfg.nnn     -  configuration list for partial wave nnn = klsp
!                    (cfg.inp in case klsp = 0, default)
!                   
!     int_inf.nnn -  description of existing angular coefficients
!                    (optional; int_inf in case klsp = 0)
!                    
!     OUTPUT FILES:
!     
!     int_int.nnn  - output data bank for angular coefficients
!                    (int_int in case klsp = 0)
!                    
!---------------------------------------------------------------------     
      Use bsr_breit
      Use conf_LS,      only: ne
      Use det_list,     only: ndet,ldet
      Use def_list,     only: ndef,ldef

      Implicit none 
      Integer :: i,ii,l,ml 
      Real(8) :: tt1,tt2, ttt, total_time

      Call CPU_TIME(t0)

!----------------------------------------------------------------------
! ... HEADER:

      Open(pri,file=AF_p)
      write(pri,'(/20x,a/20x,a/20x,a/)') &
             '=====================================',     &
             ' CALCULATION OF ANGULAR COEFFICIENTS ',     &
             '====================================='

! ... read arguments from command line:

      Call Read_arg 

!----------------------------------------------------------------------
!                                             cycle over partial waves:
      total_time = 0.d0
      Do klsp = klsp1,klsp2

       Call CPU_TIME(tt1)

       write(pri,'(80(''-''))') 
       if(klsp.gt.0) write(pri,'(/a,i5)') 'Partial wave: ',klsp
       if(klsp.gt.0) write(*  ,'(/a,i5)') 'Partial wave: ',klsp

! ... open relevant files: 

       Call open_c_file
       Call open_int_inf

       if(new.eq.1) write(pri,'(/a)') 'It is new calculations '
       if(new.eq.0) write(pri,'(/a)') 'It is continued calculations '

! ...  read configuration list:

       Call Read_conf;    if(icalc.eq.0) Cycle

! ...  load the old det.factors: 

       Call Read_dets(nub,new)

! ... define possible mls orbitals:
      
       Call Def_maxl(l);  ml=4*l+2
       Call Alloc_spin_orbitals(ne,ml)

! ... prepare det. expantions:

       Call open_det_exp

! ... calculations for new angular symmetries:

       Call open_int_int

!     *****************
       Call Conf_loop       ! main calculations
!     *****************

! ... record results:

       Call Record_int_inf

! ... print some main dimensions:      
      
       write(pri,'(/a/)') &
          'Results for new angular symmetry calculations:'
       write(pri,'(a,i10,f10.1)') &
          'number of overlap determinants =', ndet,adet
       write(pri,'(a,i10,f10.1)') &
          'number of overlap factors      =', ndef,adef 
       write(pri,'(a,i10,f10.1)') &
          'new ang. coeff.s               =', nc_new

! ... time for one partial wave:

       Call CPU_TIME(tt2); ttt=(tt2-tt1)/60
       write(pri,'(/a,T13,F12.2,a)') 'Partial wave:',ttt,' min'
       write(*,  '(/a,T13,F12.2,a)') 'Partial wave:',ttt,' min'
       total_time = total_time + ttt
       Close(nud,status='DELETE')

      End do  ! over klsp

      if(debug.gt.0) Call Debug_printing

!----------------------------------------------------------------------
! ... total time:

      write(pri,'(/a,T13,F12.2,a)') 'Total time: ',total_time,' min'
      write(*,  '(/a,T13,F12.2,a)') 'Total time: ',total_time,' min'

      End ! Program BSR_BREIT


!======================================================================
      Subroutine Read_dets(nub,new)
!======================================================================
      Implicit none

      Integer :: nub, new

      if(new.eq.1) then 
       Call Alloc_det(-1)
       Call Alloc_def(-1)
      else
       Call Load_det(nub)
       Call Load_def(nub)
      end if

      End Subroutine Read_dets


!======================================================================
      Subroutine Record_int_inf
!======================================================================
      Use bsr_breit
      Use det_list
      Use def_list

      Implicit none

      rewind(nub)
      Call Write_symc_LS(nub)
      Call Write_symt_LS(nub)

      Call Record_oper_LS(nub)

      adet=ldet; if(ndet.gt.0) adet=adet/ndet;    Call Write_det(nub)
      adef=ldef; if(ndef.gt.0) adef=adef/ndef;    Call Write_def(nub)

      close(nub)
      close(nur)
            
      End Subroutine Record_int_inf


!======================================================================
      Subroutine Debug_printing
!======================================================================
      Use bsr_breit
      Use coef_list
      Use zoef_list
      Use boef_list

      write(pri,'(/a/)') 'debug printing:'

      write(pri,'(a,f10.1)') 'mem_zoef     = ', mem_max_zoef
      write(pri,'(a,2i10)')  'max_zoef     = ', max_zoef, izoef
      write(pri,'(a,i10/)')  'realloc_zoef = ', zoef_realloc

      write(pri,'(a,f10.1)') 'mem_coef     = ', mem_max_coef
      write(pri,'(a,2i10)')  'max_coef     = ', max_coef, icoef
      write(pri,'(a,2i10)')  'max_term     = ', max_term
      write(pri,'(a,i10/)')  'realloc_coef = ', coef_realloc

      write(pri,'(a,f10.1)') 'mem_boef     = ', mem_boef
      write(pri,'(a,f10.1)') 'mem_blk      = ', mem_blk
      write(pri,'(a,2i10)')  'max_boef     = ', max_boef, iboef
      write(pri,'(a,2i10)')  'max_blk      = ', max_blk,  iblk
      write(pri,'(a,i10 )')  'realloc_boef = ', boef_realloc
      write(pri,'(a,i10/)')  'realloc_blk  = ', blk_realloc

      End Subroutine Debug_printing

