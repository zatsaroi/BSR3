!======================================================================
      Subroutine SUB1
!======================================================================
!     drive routine for one partial wave
!----------------------------------------------------------------------
      Use bsr_mat
      Use c_data
      Use conf_LS; Use symc_list_LS; Use symt_list_LS

      Implicit none
      Real(8) :: C,t3,t4
      Integer :: i,j, m,mm, ich
      Integer, external :: Ifind_channel

!---------------------------------------------------------------------- 
! ... read configuration expansion and orbitals information: 

      Call Read_data

      if(exch_mode.eq.2) Return     ! Read_RW is already done in Read_data

      if(mode.eq.7) then;  Call Sub1_mso; Return; end if

!----------------------------------------------------------------------
! ... initialize arrays and check memory requirements:

      if(allocated(IP_channel)) Deallocate(IP_channel)
      Allocate(IP_channel(ncfg))
      Do i=1,ncfg; IP_channel(i)=Ifind_channel(i); End do

      Call Allocate_ndets(-1)  
      Call Allocate_ndefs(-1)  

      Call Memory_estimations(mm)

!----------------------------------------------------------------------
!                                                       Overlap matrix:
      Call SUB1_overlaps;   if(interrupt.gt.0) Return

!----------------------------------------------------------------------
!                                                   Interaction matrix:
! ... core-energy and target-energy shifts:

      if(interrupt.eq.0) then

       hcc = hcc * EC       
       if(npert.gt.0) then;  hcb = hcb * EC;  hbb = hbb * EC;  end if

       Do ich = 1,nch
        C = Etarg(iptar(ich))-EC     ! C = htarg(ich*(ich+1)/2) ???
        if(icc(ich,ich).ne.0) Call UPDATE_HL(ich,ich,ns,ks,sb,C)
       End do

      else 

       Call Allocate_matrix(m)
       Call Read_matrix
       read(nui) m
       read(nui) ACF
       read(nui) htarg
       read(nui) otarg

      end if

!----------------------------------------------------------------------
!                                                          L-integrals:
      Call SUB1_Lintegrals;   if(interrupt.gt.0) go to 10

!----------------------------------------------------------------------
!                                                          Z-integrals:
      Call SUB1_Zintegrals;   if(interrupt.gt.0) go to 10

!----------------------------------------------------------------------
!                                                          R-integrals:
      Call SUB1_Rkintegrals;   if(interrupt.gt.0) go to 10

!----------------------------------------------------------------------
!                                            output interaction matrix:
! ... symmetrize the diagonal blocks:

      Do ich = 1,nch
       ij = icc(ich,ich); if(ij.eq.0) Cycle
       Do i = 1,ns;  Do j = 1,i
        C=hcc(i,j,ij)+hcc(j,i,ij); hcc(i,j,ij)=C/2.d0; hcc(j,i,ij)=C/2.d0
       End do; End do
      End do

! ... orthogonal conditions:

      Call CPU_time(t3);   Call BS_ORTH;  Call CPU_time(t4)

      if(pri.gt.0) &
      write(pri,'(/a,T20,f10.2,a)') 'BS_ORTH:',(t4-t3)/60,' min '
      if(myid.eq.0) &
      write(*  ,'( a,T20,f10.2,a)') 'BS_ORTH:',(t4-t3)/60,' min '

! ... record interaction matrix: 

   10 Call Record_matrix(nuj)

!----------------------------------------------------------------------
!                                              asymptotic coefficients:  
      if(interrupt.eq.0) Call SUB1_acf

      write(nuj) mk;  write(nuj) ACF

!----------------------------------------------------------------------
!                                            target Hamiltonian matrix:

      if(iitar.ne.0.and.interrupt.eq.0) Call Target_new

      write(nuj) htarg
      write(nuj) otarg
      write(nuj) Etarg
      write(nuj) EC
      close(nuj)

      if(interrupt.eq.0) Call Target_print(pri,eps_tar)
      
      End Subroutine SUB1




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

      rewind(nui);  read(nui) i,j,k;   Call Read_matrix

      rewind(nuj);  write(nuj) ns,nch,npert;   Call Record_matrix(nuj)  

! ... Interaction matrix:

      Call Allocate_matrix(m)
      Call Read_matrix
      read(nui) m
      read(nui) ACF

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

      write(nuj) mk;  write(nuj) ACF

      close(nuj)
      
      End Subroutine SUB1_mso



!======================================================================
      Subroutine Memory_estimations(mm)
!======================================================================
      Use bsr_mat
      Use c_data
      Use conf_LS; Use symc_list_LS; Use symt_list_LS
      Use orb_overlaps, only: mem_orb_overlaps

      Implicit none
      Integer :: k,l, m,mm, npol
      Integer, external ::  memory_splines

      mm = 0 

      if(pri.gt.0) write(pri,'(/a/)') 'Memory consuming:'

      ! ... the < i | j > values: 

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

      ! ... Buffer memory:

      if(.not.allocated(CBUF)) &
      Allocate(CBUF(maxnc),itb(maxnc),jtb(maxnc),intb(maxnc),idfb(maxnc))
      m = 6*maxnc
      if(pri.gt.0) write(pri,'(a,T33,f8.1,a)')  'Buffer memory:', m*4.0/(1024*1024),'  Mb' 
      mm = mm + m

      ! ... splines:

      m = memory_splines()
      if(pri.gt.0)  write(pri,'(a,T33,f8.1,a)')  'B-splines arrays:',m*4.0/(1024*1024),'  Mb'
      mm = mm + m

      ! ... configurations:

      m = m_symc + m_symt + m_conf_LS + ncfg
      if(pri.gt.0)  write(pri,'(a,T33,f8.1,a)')  'Configurations:',m*4.0/(1024*1024),'  Mb'
      mm = mm + m

      ! ... interaction matrix:

      Call Allocate_matrix(m);   mm = mm + m
      if(pri.gt.0)  &
      write(pri,'(a,T33,f8.1,a)')  'Total estimations:',mm*4.0/(1024*1024),'  Mb'

      ! ... the < . | j > values: 

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
      Integer :: i,j,ij,k, mm, ich, flag

      Call CPU_time(t1)

      if(mode.ne.0) then
       rewind(nui);  read(nui) i,j,k;  Call Read_matrix
      end if 

      if(interrupt.gt.0.and.intercase.ne.11) then
       rewind(nuj);  write(nuj) ns,nch,npert;  Call Record_matrix(nuj);  Return
      end if

      flag = 0
    1 Continue             ! the point to repeat overlaps calculations: 

      ! ... B-spline overlaps: what if iitar=2  ???

      if(mode.eq.0.or.flag.gt.0) then
       Do ich=1,nch;  if(icc(ich,ich).eq.0) Cycle
        Call UPDATE_HL(ich,ich,ns,ks,sb,1.d0)
       End do
      end if

      if(exch_mode.eq.1) go to 2

      ! ... read data from INT.BNK:

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
       Call Allocate_matrix(mm); go to 1
       Call Check_orb_overlaps
       Call Get_orth_chan
      end if

    2 Continue

      ! ... record overlap matrix: 

      rewind(nuj);    write(nuj) ns,nch,npert;  Call Record_matrix(nuj)  

      Call CPU_time(t2)
      write(pri,'(a,T20,f10.2,a)') 'Overlaps:',(t2-t1)/60,' min '
      write(  *,'(a,T20,f10.2,a)') 'Overlaps:',(t2-t2)/60,' min '

      End Subroutine SUB1_overlaps


!======================================================================
      Subroutine SUB1_acf
!======================================================================
! ... proceed the asymptotic coefficients
!--------------------------------------------------------------------- 
      Use bsr_mat

      Implicit none
      Real(8) :: C
      Character(100) :: line
      Integer :: i,j,k, i1,i2, ii

! ... Symmetrize the ACF - matrix and add the correction
! ... from core screening for k=0:

      eps_acf = 1.d-5
      Do k = 0,mk
       Do i = 1,nch
        Do j = i,nch
         C = ACF(i,j,k) + ACF(j,i,k); if(i.ne.j) C=C*2.d0
         if(abs(C).lt.eps_acf) C=0.d0
         ACF(i,j,k) = C; ACF(j,i,k) = C
        End do
       End do
      End do

      k=0; Do i=1,kclosd; k=k+2*(4*lbs(i)+2); End do
      Do i = 1,nch;  ACF(i,i,0)=ACF(i,i,0)+k; End do

      write(pri,'(/a,i2/)') 'Asymptotic coefficients: mk = ',mk

      k=0
      Do i = 1,nch
       if(abs(ACF(i,i,0)-2*nelc).lt.eps_acf) Cycle
       if(k.eq.0) then
       write(pri,'(/a,i3/)') &
        'For k=0, afc(i,i) should be equal to 2*nelectrons =', 2*nelc
       write(pri,'(a,1Pe10.1/)') 'k=0 deviations if > eps_acf = ', eps_acf
       end if
       k=k+1
       write(pri,'(i5,2F15.6)') i,ACF(i,i,0)-2*nelc
      End do
      
      if(pri_ac.gt.0) then
       write(pri,'(/a/)') 'Asymptotic coefficients: i,j, ACF(i,j,k)'
       Do k=0,mk
        if(SUM(acf(:,:,k)).eq.0) Cycle
        write(pri,'(a,i2)') 'k = ',k
        ii = 0
        Do i=1,nch; Do j = 1,i      
         if(abs(acf(i,j,k)).lt.eps_acf) Cycle
         i1=ii*20+1; i2=i1+19 
         write(line(i1:i2),'(2i4,E12.3)') j,i,acf(i,j,k)
         ii=ii+1
         if(ii.lt.5) Cycle
         write(pri,'(a)') line; ii=0
        End do; End do
        if(ii.eq.0) Cycle
        i1=1; i2=ii*20
        write(pri,'(a)') line(i1:i2)
       End do
      end if

      if(pri_f.gt.0) Call f_values

      End Subroutine SUB1_acf


