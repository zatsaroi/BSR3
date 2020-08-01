!=====================================================================
!     PROGRAM   B S R _ B R E I T _ M P I              version: 3
!
!               C O P Y R I G H T -- 2011
!
!     Written by:   Oleg Zatsarinny
!                   email: oleg_zoi@yahoo.com
!======================================================================
!    generates angular coefficient in non-orthogonal mode 
!----------------------------------------------------------------------
!
!    INPUT ARGUMENTS:
!    
!    klsp1,klsp2  - range of partial wave in BSR calculations,
!                   then cfg.001, cfg.002, ..., are input files
!                   (default -> 0, with input file is cfg.inp)
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
!    mk    - max.multipole index (default -> 9, see module param_br)
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

      USE MPI

      USE bsr_breit
      USE conf_LS,       only: ne
      Use spin_orbitals, only: mls_max

      Implicit none 

      Integer :: l
      Real(8) :: t1,t2,tt
 
!----------------------------------------------------------------------
! ... initialize MPI:

      Call MPI_INIT(ierr)
      Call MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierr)
      Call MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierr)

      if(myid.eq.0) then
       open(pri,file=AF_p)
       write(*,'(a,i6)') 'nprocs = ', nprocs
       write(pri,'(a,i6)') 'nprocs = ', nprocs
       Allocate(ip_proc(nprocs)); ip_proc=0
      end if

!----------------------------------------------------------------------
! ... read arguments from command line:

      if(myid.eq.0) Call Read_arg;    Call br_arg   

! ... log - file:

      write(AF_p,'(a,i3.3)') 'br_',myid

      if(debug.eq.0.and.myid.gt.0) pri=0
      if(myid.gt.0.and.pri.gt.0) open(pri,file=AF_p)

      if(pri.gt.0) then
      write(pri,'(/20x,a/20x,a/20x,a/)') &
             '=======================',     &
             ' B R E I T - P A U L I ',     &
             '======================='
      write(pri,'(/a,i3)')     'Max.multipole index =',mk
      write(pri,'(/a,E10.1/)') 'Tollerance for coeff.s =',eps_c
      end if

!----------------------------------------------------------------------
!                                             cycle over partial waves:
      Do klsp = klsp1,klsp2

       t1 = MPI_WTIME()

! ... open relavent files: 

       if(myid.eq.0) then
        Call open_br(nuc)       ! c-file
        Call open_br(nub)       ! data bank results, if any
        Call open_br(nur)       ! new results
        Call open_br(nui)       ! intermediate results
        write(pri,'(80(''-''))') 
        write(pri,'(/a,i5/)') ' Partial wave: ',klsp
        if(new.eq.1) write(pri,'(a/)') ' It is new calculations '
        if(new.eq.0) write(pri,'(a/)') ' It is continued calculations '
       end if 

! ... read the configuration list:

       if(myid.eq.0) Call R_conf

       Call MPI_BCAST(icalc,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

       if(icalc.eq.0) Cycle

       Call MPI_BCAST(ne,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
 
       if(myid.eq.0) then; Call Def_maxl(l);  mls_max=4*l+2; end if

       Call MPI_BCAST(mls_max,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

       Call Alloc_spin_orbitals(ne)  

! ... extract old results: 

       if(myid.eq.0)   Call Read_dets(nub,new)

! ... prepare det. expantions:

       if(myid.eq.0) then
        Call open_br(nua); Call open_br(nud);   Call Pre_det_exp 
       end if

       Call MPI_BARRIER(MPI_COMM_WORLD, ierr)

! ... calculations for new angular symmetries:

       if(myid.eq.0)     Call Conf_loop 
       if(myid.gt.0)     Call Conf_calc

! ... record results:

       if(myid.eq.0) Call Record_results

! ... time for one partial wave:

       t2=MPI_WTIME(); tt=(t2-t1)/60
       if(pri.gt.0) &
       write(pri,'(/a,F12.2,a)') ' Partial wave:',tt,' min'

      End do  ! over klsp

      Call MPI_FINALIZE(l)

      END ! Program BR_MPI



!======================================================================
      Subroutine Read_dets(nub,new)
!======================================================================
!     read overlap factors if any
!----------------------------------------------------------------------
      Implicit none
      Integer :: nub, new

      if(new.eq.1) then 
       Call Alloc_det(-1)
       Call Alloc_def(-1)
      else
       Call Read_det(nub)
       Call Read_def(nub)
      end if

      End Subroutine Read_dets


!======================================================================
      Subroutine Record_results
!======================================================================
      USE bsr_breit
      USE det_list
      USE def_list

      Implicit none
      Character AS*80
      Integer :: i, ii, nc
      Real(8) :: adet,adef

      rewind(nur)
      Call Write_symc_LS(nur)
      Call Write_symt_LS(nur)
      Call Write_oper_LS(nur)

      Call Record_det(nur)
      Call Record_def(nur)

      nc_old = 0
      if(new.eq.0) Call RW(nub,nur,nc_old); Close(nub)

      write(AS,'(a,a,a,a)') 'cat ',trim(AF_i),' >> ',trim(AF_r)
      Call System(trim(AS))

 write(*,*) AS

      Close(nui,status='DELETE')

!     rewind(nui); Call RW(nui,nur,nc); Close(nui)

      Close(nur)

! ... print the main dimensions:      

      write(pri,'(/a/)') &
          ' Results for new angular symmetry calculations:'
      adet=ldet; adet=adet/ndet
      write(pri,'(a,i10,f10.1,i10)') &
          ' number of overlap determinants =', ndet,adet,ldet
      adet=ldef; adef=adef/ndef
      write(pri,'(a,i10,f10.1,i10)') &
          ' number of overlap factors      =', ndef,adef,ldef 
      nc_total = nc_old + nc_new 
      write(pri,'(a,3i10)') &
          ' total number of coeff.s        =', nc_total, nc_old, nc_new

! ... rename new results as new data bank (int_res -> int_bnk): 
 
      if(klsp.eq.0) then
       write(AS,'(a,a,a,a)') 'mv ',trim(AF_r),' ',trim(AF_b)
      else
       write(AS,'(a,a,a,a)') 'mv ',trim(AF_r),' ',trim(BF_b)
      end if

 write(*,*) AS

      Call System(trim(AS))

      End Subroutine Record_results


!======================================================================
      Subroutine RW(nu1,nu2,nc)
!======================================================================
!     re-write bnk-data from file 'nu1' to 'nu2' by blocks
!----------------------------------------------------------------------

      Implicit none

      Integer, Intent(in) :: nu1,nu2
      Integer, Intent(out) :: nc
      Integer :: i,j

      Integer,Parameter :: mc = 1000000
      Integer,Allocatable,Dimension(:) :: K1,K2,K3
      Real(8),Allocatable,Dimension(:) :: C

      Allocate(C(mc),K1(mc),K2(mc),K3(mc))

      i = 1
    1 read(nu1,end=2) c(i),k1(i),k2(i),k3(i)
      i = i + 1; if(i.le.mc) go to 1
    2 j = i - 1
      nc = nc + j

      Do i = 1,j
       write(nu2) c(i),k1(i),k2(i),k3(i)
      End do

      i = 1;  if(j.eq.mc) go to 1

      Deallocate(C,K1,K2,K3)

      End Subroutine RW


