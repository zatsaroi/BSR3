!======================================================================
      Subroutine SUB1
!======================================================================
!     drive routine for one partial wave
!----------------------------------------------------------------------
      USE mpi
      USE bsr_mat;   USE bsr_matrix;   USE cmdata
      USE target;    USE channel
      USE conf_LS;   USE symc_list_LS; Use symt_list_LS
      USE spline_atomic,   only: EC,kclosd
      USE spline_param,    only: ns,ks
      USE spline_orbitals, only: mbf,nbf,lbs      
      USE spline_galerkin, only: sb
      USE spline_grid,     only: t

      Implicit none
      Real(8) :: C,t3,t4
      Integer :: i,j,k,l,ich, m,mm
      Integer, external :: Ifind_channel, memory_splines

      t3 = MPI_WTIME()

! ... read configuration expansion and orbitals information: 

      if(myid.eq.0) Call Read_data

      Call MPI_BARRIER(MPI_COMM_WORLD, ierr)

      Call br_orb_LS
      Call br_bsorb

      Call MPI_BARRIER(MPI_COMM_WORLD, ierr)

      Call br_symc_LS
      Call br_symt_LS
      Call MPI_BARRIER(MPI_COMM_WORLD, ierr)

      Call br_conf_LS

      Call br_channel
      Call br_phys_orb

      if(myid.ne.0) then
       if(debug.gt.0) then
        write(AF_pri,'(a,i5.5)') 'debug.',myid
        Open(pri,file=AF_pri)
       else
        pri = 0
       end if
      end if 

! ... initialize arrays:

      if(allocated(IP_channel)) Deallocate(IP_channel)
      Allocate(IP_channel(ncfg))
      Do i=1,ncfg; IP_channel(i)=Ifind_channel(i); End do

      Call Allocate_ndets(-1)  
      Call Allocate_ndefs(-1)  

! ... memory requirements:

      Call Allocate_matrix(nch,ns,npert,mk,ipconf(nch),i)
      Call MPI_REDUCE(i,m,1,MPI_INTEGER,MPI_MAX,0,MPI_COMM_WORLD,ierr)
      if(myid.ne.0) m=i
      if(pri.gt.0) write(pri,'(a,f10.2,a)') &
       'Matrix memory:     ', m*4.0/(1024*1024),'  Mb' 
      mm = m

      l=maxval(lbs(1:nbf))
      npol=max(l,mk); k=(npol+2)*ntype; if(nb.lt.k) nb=k
      Call Allocate_cmdata(nb,mb,kb,m)
      if(pri.gt.0) write(pri,'(a,f10.2,a)') &
       'CMDATA memory:     ', m*4.0/(1024*1024),'  Mb' 
      mm = mm + m

      if(.not.allocated(CBUF)) &
       Allocate(CBUF(maxnc),ijtb(maxnc),intb(maxnc),idfb(maxnc))
      m = 5*maxnc
      if(pri.gt.0) write(pri,'(a,f10.2,a)') &
       'Buffer memory:     ', m*4.0/(1024*1024),'  Mb' 
      mm = mm + m

      m = 2*(4*ns*ks + + 3*ns*ns + nbf*nbf + ns*nbf + ns*ns*(l+1))
      if(pri.gt.0)  write(pri,'(a,f10.2,a)') &
       'L-integrals memory:',m*4.0/(1024*1024),'  Mb'
      mm = mm + m

      m = mbf * 6 + 3*nbf*nbf + 4*ns*nbf  
      if(pri.gt.0)  write(pri,'(a,f10.2,a)') &
       'Orbitals arrays:   ',m*4.0/(1024*1024),'  Mb'
      mm = mm + m

      m = memory_splines()
      if(pri.gt.0)  write(pri,'(a,f10.2,a)') &
       'B-splines arrays:  ',m*4.0/(1024*1024),'  Mb'
      mm = mm + m

      m = m_symc + m_symt + m_conf_LS + ncfg
      if(pri.gt.0)  write(pri,'(a,f10.2,a)') &
       'Configurations:    ',m*4.0/(1024*1024),'  Mb'
      mm = mm + m

      if(pri.gt.0)  write(pri,'(/a,f10.2,a/)') &
       'Total estimations: ',mm*4.0/(1024*1024),'  Mb'

      t4 = MPI_WTIME()
      if(pri.gt.0) &
       write(pri,'(/a,f10.2,a/)') 'Initilize:     ',(t4-t3)/60,' min '
      if(myid.eq.0) &
       write(*  ,'(a,f10.2,a)')   'Initilize:     ',(t4-t3)/60,' min '

!-----------------------------------------------------------------------
!                                                        overlap matrix:
      t3 = MPI_WTIME()

! ... B-spline overlaps:

      Do ich=1,nch
       if(icc(ich,ich).ne.0) Call UPDATE_HL(ich,ich,ns,ks,sb,1.d0)
      End do

! ... read data from INT.BNK ...

      icase = 11;  kpol = 0;  Call State_res

! ... symmetrize the diagonal blokcs;

      Do ich = 1,nch
       ij = icc(ich,ich); if(ij.eq.0) Cycle
       Do i = 1,ns;  Do j = 1,i
        C=hcc(i,j,ij)+hcc(j,i,ij);  hcc(i,j,ij)=C/2.d0;  hcc(j,i,ij)=C/2.d0
       End do; End do
      End do

! ... check big overlaps:

      k=0;  if(S_ovl.gt.0.d0) Call Check_mat(k)

! ... redo overlap matrix:

      if(k.gt.0) then
       
       if(pri.gt.0) write(pri,*) 'redo overlap matrix: k =',k
       hcc=0.d0; acf=0.d0; htarg=0.d0; otarg=0.d0
       if(kcp.gt.0) then;  hcb=0.d0;  hbb=0.d0; end if
       Do ich=1,nch
        if(icc(ich,ich).ne.0) Call UPDATE_HL(ich,ich,ns,ks,sb,1.d0)
       End do
       icase = 11;  kpol = 0;  Call State_res
       Do ich = 1,nch
        ij = icc(ich,ich); if(ij.eq.0) Cycle
        Do i = 1,ns;  Do j = 1,i
         C=hcc(i,j,ij)+hcc(j,i,ij);  hcc(i,j,ij)=C/2.d0;  hcc(j,i,ij)=C/2.d0
        End do; End do
       End do

      end if

! ... record overlap matrix: 

      if(myid.eq.0) then; rewind(nui);  write(nui) ns,nch,kcp; end if

      Call Record_matrix(nui)  

      Call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      t4 = MPI_WTIME()
      if(pri.gt.0) &
       write(pri,'(/a,f10.2,a/)') 'Overlaps:     ',(t4-t3)/60,' min '
      if(myid.eq.0) &
       write(*   ,'(a,f10.2,a)') 'Overlaps:     ',(t4-t3)/60,' min '

!----------------------------------------------------------------------
!                                                   Interaction matrix:
      hcc = hcc * EC       ! core-energy shift
      if(kcp.gt.0) then
       hcb = hcb * EC
       hbb = hbb * EC
      end if
!----------------------------------------------------------------------
!                                                          L-integrals:
      t3 = MPI_WTIME()

      Call Gen_Lval

      icase=6;   Call State_res   

      Call Alloc_Lcore(0,0,0)

      t4 = MPI_WTIME()

      if(pri.gt.0) &
       write(pri,'(/a,f10.2,a/)') 'L-integrals:  ',(t4-t3)/60,' min '
      if(myid.eq.0) &
       write(*,'(a,f10.2,a)') 'L-integrals:  ',(t4-t3)/60,' min '

!----------------------------------------------------------------------
!                                                          Z-integrals:
      if(mso.gt.0) then    

       t3 = MPI_WTIME()

       Call Gen_Zval;  icase=7;  Call State_res   

       Call Alloc_zcore(0,0,0)

       t4 = MPI_WTIME()
       
       Call MPI_BARRIER(MPI_COMM_WORLD, ierr)
       if(myid.eq.0) &
        write(*,'(a,f10.2,a)') 'Z-integrals:  ',(t4-t3)/60,' min '
       if(pri.gt.0) &
        write(pri,'(/a,f10.2,a/)') 'Z-integrals:  ',(t4-t3)/60,' min '
 
      end if
!----------------------------------------------------------------------
!                                                          R-integrals:
      Do icase = 3,10

       Select case(icase)
        case(6,7); Cycle
        case(8,9); if(msoo.eq.-1.or.(mrel.lt.3.and.msoo.ne.1)) Cycle
        case( 10); if(mss.eq.-1.or.(mrel.lt.4.and.mss.ne.1)) Cycle
        case(3,4); if(moo.eq.-1.or.(mrel.lt.5.and.moo.ne.1)) Cycle
       End Select

       t3=MPI_WTIME();    Call State_res;   t4=MPI_WTIME()

       Call MPI_BARRIER(MPI_COMM_WORLD, ierr)
       if(pri.gt.0) &
        write(pri,'(a,a,f10.2,a)') Aint(icase),'-integrals:  ',&
                                      (t4-t3)/60,' min '
       if(myid.eq.0) &
        write(*,'(a,a,f10.2,a)') Aint(icase),'-integrals:  ',&
                                    (t4-t3)/60,' min '
      End do  ! over icase

!----------------------------------------------------------------------
!                                                      target energies:
      Do ich = 1,nch
        C = Etarg(iptar(ich))-EC
!       C=htarg(ich*(ich+1)/2)
       if(icc(ich,ich).ne.0) Call UPDATE_HL(ich,ich,ns,ks,sb,C)
      End do

!----------------------------------------------------------------------
! ... symmetrize the int. matrix: 

      Do ich = 1,nch
       ij = icc(ich,ich); if(ij.eq.0) Cycle
       Do i = 1,ns;  Do j = 1,i
        C=hcc(i,j,ij)+hcc(j,i,ij); hcc(i,j,ij)=C/2.d0; hcc(j,i,ij)=C/2.d0
       End do; End do
      End do

! ... orthogonal conditions:

      t3 = MPI_WTIME();   Call BS_ORTH;  t4 = MPI_WTIME()

      if(pri.gt.0) &
       write(pri,'(a,6x,f10.2,a)') 'BS_ORTH:  ',(t4-t3)/60,' min '
      if(myid.eq.0) &
       write(*  ,'(a,4x,f10.2,a)') 'BS_ORTH:  ',(t4-t3)/60,' min '

! ... record interaction matrix: 

      Call Record_matrix(nui)

!----------------------------------------------------------------------
! ... asymptotic coefficients:

      Call Collect_ACF

      if(myid.eq.0) then
       write(nui) mk;  write(nui) ACF;  write(nui) t(ns+1),ns
       if(pri_f.gt.0) Call f_values
      end if

!----------------------------------------------------------------------
! ... target interaction matrix:

      Call Collect_otarg 
      Call Collect_htarg 

      if(myid.eq.0) then
       write(nui) htarg
       write(nui) otarg
       write(nui) Etarg
       write(nui) EC
       close(nui)
 
       Call Target_print(pri,EC,Eps_tar)

      end if

      End Subroutine SUB1 



