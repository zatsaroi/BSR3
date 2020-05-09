!======================================================================
      Subroutine SUB1
!======================================================================
!     drive routine for one partial wave
!----------------------------------------------------------------------
      Use MPI
      Use bsr_recoup

      Use channels_ion, only: nlsp_LS => nlsp, nch_LS => nch, &
                              ipar_LS => ipar, lpar_LS => lpar, ispar_LS => ispar
      
      Implicit real(8) (A-H,O-Z)

      if(myid.eq.0) then  
       Call R_channel(nut ,klsp)
       i=len_trim(AF_mat)
       write(AF_mat(i-2:i),'(i3.3)') klsp
       open(nui,file=AF_mat,form='UNFORMATTED',action='WRITE')
       write(nui) ns,nch,npert
      end if
      Call br_channel

      if(pri.gt.0) then
       write(pri,'(a,i6,a)') 'jpar =',jpar, ' - total 2J'
       write(pri,'(a,i6,a)') 'ipar =',ipar, ' - parity'
       write(pri,'(a,i6,a)') 'nch  =',nch,  ' - number of channels'
       write(pri,'(a,i6,a)') 'npert=',npert,' - number of perturbers'
       write(pri,'(a,i6,a)') 'ns   =',ns,   ' - channel block dimension'
      end if

! ... overlaps

      Call CPU_TIME(t1)

      Call Allocate_matrix(m)

      Call MPI_REDUCE(m,i,1,MPI_INTEGER,MPI_MAX,0,MPI_COMM_WORLD,ierr)

      if(myid.ne.0) m=i
      if(pri.gt.0) &
      write(pri,'(/a,T20,f8.1,a/)') 'Matrix memory:',m*4.0/(1024*1024),'  Mb'

      mui = kui

      Do ilsp = 1,nlsp_LS; kch = nch_LS(ilsp); LL = 2*lpar_LS(ilsp)+1
                                               IS = ispar_LS(ilsp) 
       if(kch.eq.0) Cycle
       if(ipar.ne.ipar_LS(ilsp)) Cycle
       if(ITRI (jpar,LL,IS).eq.0) Cycle
       Call CPU_TIME(t3)
       if(myid.eq.0) then
        i=len_trim(BF_mat)
        write(BF_mat(i-2:i),'(i3.3)') ilsp
        Call Check_file(BF_mat)
        mui = mui + 1
        open(mui,file=BF_mat,form='UNFORMATTED',action='READ')
        rewind(mui)       
        read(mui) ns1, nch1, npert1
       end if

       Do 

        if(myid.eq.0)  read(mui) ic,jc
        Call MPI_BCAST(ic,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        Call MPI_BCAST(jc,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        if(ic.le.0) Exit
        if(ic.gt.kch.or.jc.gt.kch) Stop 'not channel-channel block'

        if(myid.eq.0) read(mui) x(:,:)

        Call MPI_BCAST(x,ns*ns,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

        Call Add_block(ilsp,ic,jc,1)
         if(ic.ne.jc)  Call Add_block(ilsp,jc,ic,1) 
       End do

       Call CPU_TIME(t4)
       if(pri.gt.0) &
       write(pri,'(a,4i6,f8.2,a)') &
       'Overlap: ilsp, nch, L, 2S+1 =', ilsp,kch,lpar_LS(ilsp),IS, &
        (t4-t3)/60, ' min'
      End do
 
      Call CPU_TIME(t3)
      Call Record_matrix(nui)
      Call CPU_TIME(t4)

      if(pri.gt.0) &
      write(pri,'(/a,5x,T20,f8.2,a)') 'Record matix:', (t4-t3)/60, ' min'

      Call CPU_TIME(t2)

      if(pri.gt.0) &
      write(pri,'(/a,5x,T20,f8.2,a/)') 'Overlap matix:', (t2-t1)/60, ' min'

!----------------------------------------------------------------------
! ... interaction matrix

      Call CPU_TIME(t1)
      hcc = 0.d0

      mui = kui
      Do ilsp = 1,nlsp_LS; kch = nch_LS(ilsp); LL = 2*lpar_LS(ilsp)+1
                                               IS = ispar_LS(ilsp) 
       if(kch.eq.0) Cycle
       if(ipar.ne.ipar_LS(ilsp)) Cycle
       if(ITRI (jpar,LL,IS).eq.0) Cycle
       Call CPU_TIME(t3)

       mui = mui + 1

       Do 
        if(myid.eq.0)  read(mui) ic,jc
        Call MPI_BCAST(ic,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        Call MPI_BCAST(jc,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

        if(ic.le.0) Exit
        if(ic.gt.kch.or.jc.gt.kch) Stop 'not channel-channel block'

        if(myid.eq.0) read(mui) x(:,:)
        Call MPI_BCAST(x,ns*ns,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

        Call Add_block(ilsp,ic,jc,1)
         if(ic.ne.jc)  Call Add_block(ilsp,jc,ic,1) 
       End do

       Call CPU_TIME(t4)
       if(pri.gt.0) &
       write(pri,'(a,4i6,f8.2,a)') &
         'Interaction: ilsp, nch, L, 2S+1 =', ilsp,kch,lpar_LS(ilsp),IS, &
          (t4-t3)/60, ' min'

      End do

      Call CPU_TIME(t3)
      Call Record_matrix(nui)
      Call CPU_TIME(t4)
      if(pri.gt.0) &
      write(pri,'(/a,5x,T20,f8.2,a)') 'Record matix:', (t4-t3)/60, ' min'

      Call CPU_TIME(t2)
      if(pri.gt.0) &
      write(pri,'(/a,5x,T20,f8.2,a/)') 'Interaction matix:', (t2-t1)/60, ' min'

!--------------------------------------------------------------------------------
! ... asymptotic coefficients

      Call CPU_TIME(t1)

      Deallocate(hcc)
      if(allocated(acf)) deallocate(acf)
      Allocate(acf(nch,nch,0:mk));  acf = 0.d0
      m = nch*nch*(mk+1)
      if(pri.gt.0) &
      write(pri,'(/a,T20,f8.1,a/)') 'ACF memory:',m*8.0/(1024*1024),'  Mb'

      mui = kui
      Do ilsp = 1,nlsp_LS; kch = nch_LS(ilsp); LL = 2*lpar_LS(ilsp)+1
                                               IS = ispar_LS(ilsp) 
       if(kch.eq.0) Cycle
       if(ipar.ne.ipar_LS(ilsp)) Cycle
       if(ITRI (jpar,LL,IS).eq.0) Cycle

       Call CPU_TIME(t3)

       mui = mui + 1

       if(myid.eq.0) read(mui) km
       Call MPI_BCAST(km,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

       if(allocated(bcf)) deallocate(bcf);   Allocate(bcf(kch,kch,0:km))                       

       if(myid.eq.0) read(mui) bcf
       Call MPI_BCAST(bcf,kch*kch*(km+1),MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

       Do ic=1,kch
        Do jc=1,ic
         Call Add_block(ilsp,ic,jc,2)
         if(ic.ne.jc)  Call Add_block(ilsp,jc,ic,2) 
        End do
       End do

       Call CPU_TIME(t4)
       if(pri.gt.0) &
       write(pri,'(a,4i6,f8.2,a)') &
           'ACF: ilsp, nch, L, 2S+1 =', ilsp,kch,lpar_LS(ilsp),IS, &
            (t4-t3)/60, ' min'

      End do

      Call CPU_TIME(t2)
      if(pri.gt.0) &
      write(pri,'(/a,5x,T20,f8.2,a)') 'Asymptotic:', (t2-t1)/60, ' min'

! ... Collect ACF 

      Call CPU_TIME(t1)
      Call Collect_ACF

      if(myid.eq.0)   write(nui) mk
      if(myid.eq.0)   write(nui) acf

      Call CPU_TIME(t2)
      if(pri.gt.0) &
      write(pri,'(/a,5x,T20,f8.2,a)') 'Collect_ACF:', (t2-t1)/60, ' min'

      End Subroutine SUB1

