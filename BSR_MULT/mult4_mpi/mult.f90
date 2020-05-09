!=====================================================================
!     PROGRAM       M U L T                             version 4.0
!
!               C O P Y R I G H T -- 2017
!
!     Written by:   Oleg Zatsarinny
!                   email: oleg_zoi@yahoo.com
!======================================================================
!
!     The program evaluates MULTIPOLE operators in LS-coupling
!     including the case of non-orthogonal orbitals
!     
!     The used technique is described in:
!
!     Comp.Phys.Commun. 98 (1996) 235-254
!
!======================================================================
!
!    INPUT ARGUMENTS:
!    
!     AF1  -  c-file for initial state
!     AF2  -  c-file for final state
!     AA   -  type of calculation:  E1,M1,E2,M2,...
! 
!     INPUT FILES:
!
!     AF1  -  c-file for initial state
!     AF2  -  c-file for final state
!     mult_bnk  - data-bank for angular coefficients (optional)
!
!     OUTPUT FILES:
!
!     mult_bnk     -  new data bank for angular coefficients
!                     for given transition AA
!     mult.log     -  running information
!
!-----------------------------------------------------------------------
!     ONLY ONE TYPE OF TRANSITION IN ONE RUN !!!
!-----------------------------------------------------------------------
      Use MPI
      Use mult_par

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

      t0 = MPI_WTIME()                                                           

!----------------------------------------------------------------------
! ... read arguments from command line:

      if(myid.eq.0) Call Read_arg;    Call br_arg   

! ... open log-file:

      if(myid.ne.0) write(AF_p,'(a,i4.4)') 'mult_',myid

      if(debug.eq.0.and.myid.gt.0) pri=0

      if(pri.gt.0) then
       if(myid.ne.0) open(pri,file=AF_p)
       write(pri,'(/a)') &
	    'MULT - GENERATION OF DATA BANK FOR MULTIPOLE MATRIX ELEMENTS: '

       write(pri,'(/a,a1,i1)') 'transition -> ',ktype,kpol

       if(new.eq.1) &
        write(pri,'(/a,a)')  'new calculation for ',AF_b
       if(new.eq.0) &
        write(pri,'(/a,a)')  'continued calculation for ',AF_b

        write(pri,'(/a,E10.1/)') 'Tollerance for coeff.s =',eps_c
      end if

!-----------------------------------------------------------------------
! ... calculations:

      Call Mult_calc

      Call MPI_FINALIZE(l)

      End ! Program MULT



!=====================================================================
      Subroutine  MULT_calc
!======================================================================
!     calculations for given mult_bnk and c-files:
!-----------------------------------------------------------------------
      Use MPI
      Use mult_par
      Use conf_LS,  only: ne
      Use det_list, only: ndet,ldet,jdet
      Use def_list, only: ndef,ldef,jdef
      Use spin_orbitals, only: mls_max

      Implicit none 
      Character :: AS*80, kt*1
      Integer :: nc, k, l

      t1 = MPI_WTIME()                                                           

!------------------------------------------------------------------
! ... read the configuration list:

      if(myid.eq.0) Call R_conf

      Call MPI_BCAST(icalc,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

      t2=MPI_WTIME()
      if(pri.gt.0) &
      write(pri,'(/a,F12.2,a)') 'Read_conf is done:',(t2-t1)/60,' min'

      if(icalc.eq.0) Return

      Call MPI_BCAST(ne,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
 
      if(myid.eq.0) then; Call Def_maxl(l);  mls_max=4*l+2; end if

      Call MPI_BCAST(mls_max,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

      Call Alloc_spin_orbitals(ne)  

! ... extract old results: 

      if(myid.eq.0)   Call Read_dets(nub,new)

! ... prepare det. expantions:

      if(myid.eq.0) then
       Open(nua,form='UNFORMATTED')
       Open(nud,form='UNFORMATTED')
       Call Pre_det_exp 
       Open(nui,form='UNFORMATTED')
      end if

      Call MPI_BARRIER(MPI_COMM_WORLD, ierr)

      t3=MPI_WTIME()
      if(pri.gt.0) &
      write(pri,'(/a,F12.2,a)') 'Prep_det is done:',(t3-t2)/60,' min'

!----------------------------------------------------------------------------
! ... calculations for new angular symmetries:

      if(myid.eq.0)     Call Conf_loop 
      if(myid.gt.0)     Call Conf_calc

      Call MPI_BARRIER(MPI_COMM_WORLD, ierr)

      t4=MPI_WTIME()
      if(pri.gt.0) &
      write(pri,'(/a,F12.2,a)') 'Calculations done:',(t4-t3)/60,' min'

! ... record results:

      if(myid.eq.0) &
      write(pri,*) '-> Record_results, myid',  myid

      if(myid.eq.0)  Call Record_results

      Call MPI_BARRIER(MPI_COMM_WORLD, ierr)

      t2=MPI_wtime()
      if(pri.gt.0) & 
      write(pri,'(a,F12.2,a)') 'time:',(t2-t0)/60,' min'

      End Subroutine MULT_calc



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
      Subroutine Record_results
!======================================================================
      Use mult_par
      Use det_list
      Use def_list

      Implicit none
      Character AS*80
      Integer :: i, ii, nc
      Real(8) :: adet,adef

      Open(nur,file=AF_r,form='UNFORMATTED')
  write(pri,*) 'Record_results', nur, '   ',AF_r

      rewind(nur)
      write (nur) ktype, kpol

      Call Write_symc_LS(nur)
      Call Write_symt_LS(nur)

      Call Record_done_LS(nur)    
  
      adet=ldet; adet=adet/ndet
      Call Write_det(nur)
      adef=ldef; adef=adef/ndef
      Call Write_def(nur)

      nc = 0
      if(new.eq.0) Call RW(nub,nur,nc)
  write(pri,*) 'old nc:', nc
      rewind(nui); Call RW(nui,nur,nc)
  write(pri,*) 'new nc:', nc

      Close(in1); Close(in2); Close(nub); Close(nur)
      Close(nui,status='DELETE')
      Close(nua,status='DELETE') 
      Close(nud,status='DELETE') 
 

! ... move new results to data bank (inr_res -> int_bnk): 
 
      write(AS,'(a,1x,a,1x,a)') move,trim(AF_r),trim(AF_b)  

 write(pri,*) 'SYSTEM call: ', trim(AS)

      Call System(trim(AS))
     
  write(pri,*) 'closing files'
 
!      Close(nur,status='DELETE') 

! ... print the main dimensions:      

      write(pri,'(/a/)') &
          ' Results for new angular symmetry calculations:'
      write(pri,'(a,i10,f10.1,i10)') &
          ' number of overlap determinants =', ndet,adet,ldet
      write(pri,'(a,i10,f10.1,i10)') &
          ' number of overlap factors      =', ndef,adef,ldef 
      write(pri,'(a,i10)') &
          ' new coeff.s                    =', nc

      End Subroutine Record_results




!======================================================================
      Subroutine RW(nu1,nu2,nc)
!======================================================================
!     re-write bnk-data from file 'nu1' to 'nu2' by blocks
!----------------------------------------------------------------------
      Implicit none

      Integer, intent(in) :: nu1,nu2
      Integer :: i,j,nc

      Integer, parameter :: mc = 100000
      Integer, allocatable :: K1(:),K2(:),K3(:),K4(:)
      Real(8), allocatable :: C(:)

      Allocate(C(mc),K1(mc),K2(mc),K3(mc),K4(mc))

      i = 1
    1 read(nu1,end=2) c(i),k1(i),k2(i),k3(i),k4(i)
      i = i + 1; if(i.le.mc) go to 1
    2 j = i - 1
      nc = nc + j

      Do i = 1,j
       write(nu2) c(i),k1(i),k2(i),k3(i),k4(i)
      End do

      i = 1;  if(j.eq.mc) go to 1

      Deallocate(C,K1,K2,K3,K4)

      End Subroutine RW


