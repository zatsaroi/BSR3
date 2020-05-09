!======================================================================
      Subroutine SUB1 
!======================================================================
!     drive routine for one partial wave
!----------------------------------------------------------------------
      Use mpi
      Use bsr_mat  
      Use c_data
      Use conf_LS; Use symc_list_LS; Use symt_list_LS
      Use orb_overlaps

      Implicit none
      Real(8) :: C,t1,t2
      Integer :: i, ich, m
      Integer, external :: Ifind_channel

!-----------------------------------------------------------------------
! ... read configuration expansion and orbitals information: 

      t1 = MPI_WTIME()

      if(myid.eq.0) Call Read_data

      if(exch_mode.eq.2) Return

      if(mode.eq.7) then;  Call Sub1_mso_mpi; Return; end if

! ... debug output if any:

      if(myid.ne.0) then
       if(debug.ge.myid) then
        write(AF_pri,'(a,i5.5)') 'debug.',myid
        Open(pri,file=AF_pri)
       else
        pri = 0      
       end if
      end if

      if(pri.gt.0) then
       write(pri,'(/a/  )') 'Main dimensions in bsr_matrix module:'
       write(pri,'( a,i10,a)') 'nch    = ',nch,   '  -  number of channels '
       write(pri,'( a,i10,a)') 'npert  = ',npert, '  -  number of perturbers '
       write(pri,'( a,i10,a)') 'ns     = ',ns,    '  -  number of splines ' 
       write(pri,'(/a,T33,i8)') 'Dimension of interaction matrix:',nch*ns+npert
      end if

      t2 = MPI_WTIME()
      if(pri.gt.0) &
      write(pri,'(/a,T20,f10.1,a)') 'Read_data:',(t2-t1)/60,' min '

!----------------------------------------------------------------------
! ... broadcast the information:

      Call MPI_BARRIER(MPI_COMM_WORLD, ierr)

      t1 = MPI_WTIME()

      Call br_symc_LS
      Call br_symt_LS
      Call br_conf_LS
      Call br_channel
      Call br_bsorb
      Call br_phys_orb
      Call br_dets
      Call br_ovl

      t2 = MPI_WTIME()
      if(pri.gt.0) &
      write(pri,'(/a,T20,f10.1,a)') 'Broad_cast:  ',(t2-t1)/60,' min '

!----------------------------------------------------------------------
! ... initialize arrays and check memory requirements:

      if(allocated(IP_channel)) Deallocate(IP_channel)
      Allocate(IP_channel(ncfg))
      Do i=1,ncfg; IP_channel(i)=Ifind_channel(i); End do

      Call Allocate_ndets(-1)  
      Call Allocate_ndefs(-1)  

      Call Memory_estimations

!----------------------------------------------------------------------
!                                                        overlap matrix:
      Call SUB1_overlaps;  if(interrupt.gt.0) Return

!----------------------------------------------------------------------
! ...  Interaction matrix:
!----------------------------------------------------------------------
!                                  core-energy and target-energy shift:
      if(mode.eq.0) then

       hcc = hcc * EC       
       if(npert.gt.0) then;  hcb = hcb * EC; hbb = hbb * EC; end if

       Do ich = 1,nch
        C = Etarg(iptar(ich))-EC      ! htarg - ???
        if(icc(ich,ich).ne.0) Call UPDATE_HL(ich,ich,ns,ks,sb,C)
       End do

      else

       Call Allocate_matrix(i)
       Call Read_matrix_mpi
       if(myid.eq.0) then
        read(nui) m         
        read(nui) ACF
        read(nui) htarg
        read(nui) otarg
       end if

       m = nch*nch*(mk+1)
       Call MPI_BCAST(acf  ,m,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
       m = (nch+1)*nch/2
       Call MPI_BCAST(htarg,m,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
       Call MPI_BCAST(otarg,m,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

      end if
       
!----------------------------------------------------------------------
!                                                          L-integrals:
      Call SUB1_Lintegrals;   if(interrupt.gt.0) go to 10

!----------------------------------------------------------------------
!                                                          Z-integrals:
      Call SUB1_Zintegrals;   if(interrupt.gt.0) go to 10

!----------------------------------------------------------------------
!                                                          R-integrals:
      Call SUB1_Rkintegrals;  if(interrupt.gt.0) go to 10

!----------------------------------------------------------------------
! ... orthogonal conditions:

      t1 = MPI_WTIME();   Call BS_ORTH;    t2 = MPI_WTIME()

      if(pri.gt.0) &
      write(pri,'(a,T20,f10.1,a)') 'BS_ORTH:  ',(t2-t1)/60,' min '
      if(myid.eq.0) &
      write(*  ,'(a,T20,f10.1,a)') 'BS_ORTH:  ',(t2-t1)/60,' min '

! ... record interaction matrix: 
 
   10 t1 = MPI_WTIME();   Call Record_matrix(nuj);   t2 = MPI_WTIME()

      if(pri.gt.0) &
      write(pri,'(/a,T20,f10.1,a)') 'Record matrix:',(t2-t1)/60,' min '
      if(myid.eq.0) &
      write(*  ,'( a,T20,f10.1,a)') 'Record matrix:',(t2-t1)/60,' min '

!----------------------------------------------------------------------
! ... asymptotic coefficients:

      t1 = MPI_WTIME();    Call Collect_ACF;    t2 = MPI_WTIME()

      if(pri.gt.0) &
      write(pri,'(/a,T20,f10.1,a)') 'Collect_ACF:',(t2-t1)/60,' min '
      if(myid.eq.0) &
      write(*  ,'( a,T20,f10.1,a)') 'Collect_ACF:',(t2-t1)/60,' min '

      if(myid.eq.0) then
       write(nuj) mk;  write(nuj) ACF
       if(pri_f.gt.0) Call f_values
      end if

!----------------------------------------------------------------------
! ... target interaction matrix:

      if(interrupt.eq.0) then
       t1 = MPI_WTIME();  Call Collect_otarg;  Call Collect_htarg; t2 = MPI_WTIME() 
       if(pri.gt.0) &
        write(pri,'(/a,T20,f10.1,a)') 'Collect_targ:',(t2-t1)/60,' min '
       if(myid.eq.0) &
        write(*  ,'( a,T20,f10.1,a)') 'Collect_targ:',(t2-t1)/60,' min '
      end if

      if(myid.eq.0) then
       if(iitar.ne.0.and.interrupt.eq.0) Call Target_new
       write(nuj) htarg
       write(nuj) otarg
       write(nuj) Etarg
       write(nuj) EC
       close(nuj)
       if(interrupt.eq.0) Call Target_print(pri,Eps_tar)
      end if

      End Subroutine SUB1 



!======================================================================
      Subroutine Memory_estimations
!======================================================================
      Use MPI
      Use bsr_mat
      Use c_data
      Use conf_LS; Use symc_list_LS; Use symt_list_LS
      Use orb_overlaps, only: mem_orb_overlaps

      Implicit none
      Integer :: i,k,l, m,mm, npol
      Integer, external ::  memory_splines

      mm = 0

      if(pri.gt.0) write(pri,'(/a/)') 'Memory consuming:'

! ... the < i | j > arrays: 

      m = mem_orb_overlaps 
      if(pri.gt.0) &
      write(pri,'(a,T33,f8.1,a)') 'Bound overlaps:', m*4.0/(1024*1024),'  Mb' 
      mm = mm + m

! ... c_data arrays:

      l = maxval(lbs(1:nbf))
      npol = max(l,mk); k=(npol+2)*mtype*2
      if(nblock.lt.k) then
       nblock=k
       write(pri,'(/a,i8,a)') 'nblock = ',nblock,'  -  number of blocks re-assigned !!! '
      end if
      Call Alloc_c_data(mtype,-1,npol,mblock,nblock,kblock,eps_c,m) 
      if(pri.gt.0) &
      write(pri,'(a,T33,f8.1,a)')  'Memory for c_data:', m*4.0/(1024*1024),'  Mb'
      mm = mm + m

! ... buffer:

      if(.not.allocated(CBUF)) &
      Allocate(CBUF(maxnc),itb(maxnc),jtb(maxnc),intb(maxnc),idfb(maxnc))
      m = 6*maxnc
      if(pri.gt.0) write(pri,'(a,T33,f8.1,a)') 'Buffer memory:', m*4.0/(1024*1024),'  Mb' 
      mm = mm + m

! ... splines:

      m = memory_splines()
      if(pri.gt.0)  write(pri,'(a,T33,f8.1,a)') 'B-splines arrays:',m*4.0/(1024*1024),'  Mb'
      mm = mm + m

! ... configurations;

      m = m_symc + m_symt + m_conf_LS + ncfg
      if(pri.gt.0)  write(pri,'(a,T33,f8.1,a)') 'Configurations:',m*4.0/(1024*1024),'  Mb'
      mm = mm + m

! ... interaction matrix:

      Call Allocate_matrix(i)

      Call MPI_REDUCE(i,m,1,MPI_INTEGER,MPI_MAX,0,MPI_COMM_WORLD,ierr)
      if(myid.ne.0) m=i
      mm = mm + m 
      if(pri.gt.0)  write(pri,'(a,T33,f8.1,a)') 'Total estimations: ',mm*4.0/(1024*1024),'  Mb'

      ! ... the < . | j > values: 

      Call br_file
      Call Check_orb_overlaps
      Call Get_orth_chan

      End Subroutine Memory_estimations


!======================================================================
      Subroutine SUB1_overlaps
!======================================================================
      Use bsr_mat
      Use c_data

      Implicit none
      Real(8) :: C,t1,t2
      Integer :: i,j,ij,k, ich, flag

      Call CPU_time(t1)

      if(mode.ne.0) then
       if(myid.eq.0) rewind(nui)
       if(myid.eq.0) read(nui) i,j,k
       Call Read_matrix
      end if 

      if(interrupt.gt.0.and.intercase.ne.11) then
       if(myid.eq.0) rewind(nuj)
       if(myid.eq.0) write(nuj) ns,nch,npert
       Call Record_matrix(nuj)
       Return
      end if

      flag = 0
    1 Continue             ! the point to repeat overlaps calculations: 

      ! ... B-spline overlaps:  iitar =2  ???

      if(mode.eq.0.or.flag.gt.0) then
       Do ich=1,nch
        if(icc(ich,ich).ne.0) Call UPDATE_HL(ich,ich,ns,ks,sb,1.d0)
       End do
      end if

      if(exch_mode.eq.1) go to 2

      ! ... read data from INT_BNK:

      icase = 11;   Call State_res

      ! ... symmetrize the diagonal blokcs:

      Do ich = 1,nch
       ij = icc(ich,ich); if(ij.eq.0) Cycle
       Do i = 1,ns;  Do j = 1,i
        C=hcc(i,j,ij)+hcc(j,i,ij);  hcc(i,j,ij)=C/2.d0;  hcc(j,i,ij)=C/2.d0
       End do; End do
      End do

      ! ... check big overlaps:

      Call Check_mat(flag)

      ! ... redo overlap matrix:

      if(flag.gt.0) then
       Call br_file
       Call Check_orb_overlaps
       Call Get_orth_chan
       Call Allocate_matrix(i)
       go to 1
      end if

    2 Continue

      ! ... record overlap matrix: 

      if(myid.eq.0) then;   rewind(nuj);  write(nuj) ns,nch,npert;  end if

      Call Record_matrix(nuj)  

      Call CPU_time(t2)
      if(pri.gt.0) &
      write(pri,'(/a,T20,f10.2,a)') 'Overlaps:',(t2-t1)/60,' min '
      if(myid.eq.0) &
      write(  *,'( a,T20,f10.2,a)') 'Overlaps:',(t2-t1)/60,' min '

      End Subroutine SUB1_overlaps



!======================================================================
      Subroutine SUB1_Lintegrals
!======================================================================
      Use bsr_mat

      Implicit none
      Real(8) :: t1,t2

      if(intercase.gt.0.and.intercase.ne.6) Return

      if(exch_mode.eq.1) Return

      Call CPU_time(t1)

      Call Gen_Lval;  icase=6;  Call State_res;  Call Alloc_Lcore(0,0,0)

      Call CPU_time(t2)

      if(pri.gt.0) &
      write(pri,'(a,T20,f10.2,a)') 'L-integrals:',(t2-t1)/60,' min '
      if(myid.eq.0) &
      write(*  ,'(a,T20,f10.2,a)') 'L-integrals:',(t2-t1)/60,' min '

      End Subroutine SUB1_Lintegrals 


!======================================================================
      Subroutine SUB1_Zintegrals
!======================================================================
      Use bsr_mat

      Implicit none
      Real(8) :: t1,t2
      Integer :: l

      if(intercase.gt.0.and.intercase.ne.7) Return

      if(mso.le.0) Return    

      Call CPU_time(t1)

      l = maxval(lbs); if(l.gt.mlso) l = mlso;  Call Gen_Zval(l)

      icase=7;  Call State_res;   Call Alloc_zcore(0,0)
 
      Call CPU_time(t2)
      
      if(pri.gt.0)  &
      write(pri,'(a,T20,f10.2,a)') 'Z-integrals:',(t2-t1)/60,' min '
      if(myid.eq.0) &
      write(*  ,'(a,T20,f10.2,a)') 'Z-integrals:',(t2-t1)/60,' min '

      End Subroutine SUB1_Zintegrals


!======================================================================
      Subroutine SUB1_Rkintegrals
!======================================================================
      Use bsr_mat

      Implicit none
      Real(8) :: t1,t2

      Do icase = 3,10

       if(intercase.gt.0.and.intercase.ne.icase) Cycle

       Select case(icase)
        case(6,7); Cycle
        case(8,9); if(msoo.eq.-1.or.(mrel.lt.3.and.msoo.ne.1)) Cycle
        case( 10); if(mss .eq.-1.or.(mrel.lt.4.and.mss .ne.1)) Cycle
        case(3,4); if(moo .eq.-1.or.(mrel.lt.5.and.moo .ne.1)) Cycle
       End Select

       Call CPU_time(t1);   Call State_res;    Call CPU_time(t2)

       if(pri.gt.0) &
       write(pri,'(a,a,T20,f10.2,a)') Aint(icase),'-integrals:',(t2-t1)/60,' min '
       if(myid.eq.0) &
       write(*  ,'(a,a,T20,f10.2,a)') Aint(icase),'-integrals:',(t2-t1)/60,' min '

       if(interrupt.gt.0) Return

      End do  ! over icase

      End Subroutine SUB1_Rkintegrals



!======================================================================
      Subroutine SUB1_mso
!======================================================================
!     drive routine for one partial wave
!----------------------------------------------------------------------
      Use bsr_mat
      Use c_data
      Use conf_LS; Use symc_list_LS; Use symt_list_LS

      Implicit none
      Real(8) :: C
      Integer :: i,j,k, m,mm,mmm, ich
      Integer, external :: Ifind_channel

! ... initialize arrays and check memory requirements:

      if(allocated(IP_channel)) Deallocate(IP_channel)
      Allocate(IP_channel(ncfg))
      Do i=1,ncfg; IP_channel(i)=Ifind_channel(i); End do

      Call Allocate_ndets(-1)  
      Call Allocate_ndefs(-1)  

      Call Memory_estimations(mm)

! ... re-write overlap matrix:

      if(myid.eq.0) rewind(nui)
      if(myid.eq.0) read(nui) i,j,k
      Call Read_matrix

      if(myid.eq.0) rewind(nuj)
      if(myid.eq.0) write(nuj) ns,nch,npert
      Call Record_matrix(nuj)  

! ... Interaction matrix:

      Call Allocate_matrix(m)
      Call Read_matrix
      if(myid.eq.0) read(nui) m
      if(myid.eq.0) read(nui) ACF

! ... Z-integrals:

      Call SUB1_Zintegrals

! ... symmetrize the diagonal blocks:

      Do ich = 1,nch
       ij = icc(ich,ich); if(ij.eq.0) Cycle
       Do i = 1,ns;  Do j = 1,i
        C=hcc(i,j,ij)+hcc(j,i,ij); hcc(i,j,ij)=C/2.d0; hcc(j,i,ij)=C/2.d0
       End do; End do
      End do

! ... record interaction matrix: 

      Call Record_matrix(nuj)

! ... asymptotic coefficients:  

      if(myid.eq.0) write(nuj) mk
      if(myid.eq.0) write(nuj) ACF

      close(nuj)
      
      End Subroutine SUB1_mso
