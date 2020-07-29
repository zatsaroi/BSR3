!=====================================================================
!     PROGRAM   B S R _ B R E I T                            version 3
!
!               C O P Y R I G H T -- 2013
!
!     Written by:   Oleg Zatsarinny 
!                   email: oleg_zoi@yahoo.com
!======================================================================
!    generates angular coefficient in non-orthogonal mode 
!----------------------------------------------------------------------
!
!    INPUT ARGUMENTS:
!    
!    klsp,klsp1,klsp2  - range of partial wave in BSR calculations,
!                        then cfg.001, cfg.002, ..., are input files
!                        (default -> 0, with input file is cfg.inp)
!   
!    oper  - character(7), where each position can be 0 or 1,
!            and indicate the operator under consideration:
!            oper(1) - OVERLAPS       
!            oper(2) - KINATIC ENERGY
!            oper(3) - TWO-ELECTRON ELECTROSTATIC
!            oper(4) - SPIN ORBIT       
!            oper(5) - SPIN-OTHER-ORBIT
!            oper(6) - SPIN-SPIN
!            oper(7) - ORBIT-ORBIT       
!            Default -> 1110000 - non-relativistic calculations
!
!    mk    - max.multipole index (default -> 9)
!
!----------------------------------------------------------------------
!
!    example:    1.  bsr_breit 
!                2.  bsr_breit klsp1=1 klsp2=5 oper=1111110  
!                3.  bsr_breit km=5
!            
!----------------------------------------------------------------------
!
!    INPUT FILES:
!    
!    cfg.nnn     -  configuration list for partial wave nnn = klsp
!                   (cfg.inp in case klsp = 0, default)
!                  
!    int_bnk.nnn -  input data bank for angular coefficients
!                   (optional; int_bnk in case klsp = 0)
!                   
!    
!    OUTPUT FILES:
!    
!    int_bnk.nnn  - output data bank for angular coefficients
!                   (int_bnk in case klsp = 0)
!                   
!---------------------------------------------------------------------     
      Use bsr_breit
      Use conf_LS,      only: ne
      Use det_list,     only: ndet,ldet
      Use def_list,     only: ndef,ldef

      Use symt_list_LS

      Implicit none 
      Integer :: klsp, i, ntotc, system
      Real(8) :: time, t1,t2,tt, adet,adef

      Character(80) :: AS,BS

!----------------------------------------------------------------------
! ... HEADER:

      Call open_br(pri,0) 
      write(pri,'(/20x,a/20x,a/20x,a/)') &
             '=======================',     &
             '   B S R _ B R E I T   ',     &
             '======================='

! ... read arguments from command line:

      Call Read_arg 

!----------------------------------------------------------------------
!                                             cycle over partial waves:
      time = 0.d0
      Do klsp = klsp1,klsp2

       write(pri,'(80(''-''))') 
       write(pri,'(/a,i5/)') ' Partial wave: ',klsp
       if(klsp.gt.0) write(*,'(/a,i5/)') ' Partial wave: ',klsp
       Call CPU_time(t1)

! ... open relevant files: 

       Call open_br(nuc,klsp)       ! c-file
       Call open_br(nub,klsp)       ! data bank results, if any
       Call open_br(nur,klsp)       ! new results
       Call open_br(nui,klsp)       ! intermediate results

       if(new.eq.1) write(pri,'(a/)') ' It is new calculations '
       if(new.eq.0) write(pri,'(a/)') ' It is continued calculations '

! ...  read the configuration list:

       Call Read_conf

       if(icalc.eq.0) then; Close(nur,STATUS='DELETE'); Cycle; end if        
  
! ...  extract old results: 

       if(new.eq.1) then 
        Call Alloc_det(-1)
        Call Alloc_def(-1)
       else
        Call Read_det(nub)
        Call Read_def(nub)
       end if

!----------------------------------------------------------------------
! ... define possible mls orbitals:
      
       Call Alloc_spin_orbitals(ne)

! ... prepare det. expantions:

       Call open_br(nua,0); Call open_br(nud,0)
 
       write(pri,'(/a/)') ' Determinant expansions: '

       Call Pre_det_exp 

! ... calculations for new angular symmetries:

       write(pri,'(/a/)') ' Calculations for given symmetry: '

       Call Conf_loop 

! ...  record results:

       rewind(nur)
       Call Write_symc_LS(nur)
       Call Write_symt_LS(nur)
       Call Write_oper_LS(nur)
       Call Record_det(nur)
       Call Record_def(nur)
       ntotc = 0
       if(new.eq.0) Call RW(nub,nur,ntotc)
       rewind(nui); Call RW(nui,nur,ntotc)

! ...  print the main dimensions:      

       adet=ldet; if(ndet.gt.0) adet=adet/ndet
       adef=ldef; if(ndef.gt.0) adef=adef/ndef

       write(pri,'(/a/)') &
          ' Results for new calculations:'
       write(pri,'(a,i10,f10.1)') &
          ' number of overlap determinants =', ndet,adet
       write(pri,'(a,i10,f10.1)') &
          ' number of overlap factors      =', ndef,adef 
       write(pri,'(a,i10,f10.1)') &
          ' total number of coefficient    =', ntotc

! ...  rename new results as new data bank (int_res -> int_bnk): 
 
       close(nui); close(nur); close(nub)

       if(klsp.eq.0) then
        write(AS,'(a,a,1x,a)') 'mv ',trim(AF_r),trim(AF_b)
        write(BS,'(a,a,1x,a)') 'move ',trim(AF_r),trim(AF_b)
       else
        write(AS,'(a,a,1x,a)') 'mv ',trim(BF_r),trim(BF_b)
        write(BS,'(a,a,1x,a)') 'move ',trim(BF_r),trim(BF_b)
       end if

       i = System(trim(AS))      ! Unix 
!       i = System(trim(BS))     ! Windows

! ... time for one partial wave:

       Call CPU_time(t2); tt=(t2-t1)/60
       write(pri,'(/a,F12.2,a)') ' Partial wave:',tt,' min'
       write(*,  '( a,F12.2,a)') ' Partial wave:',tt,' min'
       time = time + tt

      End do  ! over klsp
!----------------------------------------------------------------------

! ... total time:

      write(pri,'(a,F12.2,a)')   ' Total time:  ',time,' min'
      write(*,  '(a,F12.2,a)')   ' Total time:  ',time,' min'

      End ! Program Breit_bsr


!======================================================================
      Subroutine RW(nu1,nu2,nc)
!======================================================================
!     re-write bnk-data from file 'nu1' to 'nu2' by blocks
!----------------------------------------------------------------------
      Implicit none
      Integer, intent(in) :: nu1,nu2
      Integer :: i,j,nc
      Integer, parameter :: mc = 100000
      Integer, allocatable :: K1(:),K2(:),K3(:)
      Real(8), allocatable :: C(:)

      Allocate(C(mc),K1(mc),K2(mc),K3(mc))

      i = 1
    1 read(nu1,end=2) c(i),k1(i),k2(i),k3(i)
      i = i + 1; if(i.le.mc) go to 1
    2 j = i - 1
      nc = nc + j

      Do i = 1,j; write(nu2) c(i),k1(i),k2(i),k3(i); End do

      i = 1;  if(j.eq.mc) go to 1

      Deallocate(C,K1,K2,K3)

      End Subroutine RW


