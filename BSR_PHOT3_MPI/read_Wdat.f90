!======================================================================
      Subroutine Read_Wdat
!----------------------------------------------------------------------
!     read W.DAT file (weight coefficient data) 
!----------------------------------------------------------------------      

      Use bsr_phot
      Implicit none
      Integer :: i,j,kch,kcp

      if(nwt.lt.0) Return

      i=INDEX(AF_w,'.',BACK=.TRUE.)
      write(AF,'(a,i3.3)') AF_w(1:i),klsp
      Call Check_file(AF)
      Open(nuw, file=AF, form='UNFORMATTED')

      read(nuw) kch, kcp, khm
 
      if(kch.ne.nch) Stop ' other nch in W.DAT'
      if(khm.ne.nhm) Stop ' other nhm in W.DAT'
      ncp=kcp
 
      nwt = nch + ncp;  Allocate (WT(nwt,nhm),AK(nwt,nch),WTch(nwt))
      Do i = 1,nhm; read(nuw) (WT(j,i),j=1,nwt); End do
 
      Close(nuw)

      End Subroutine Read_Wdat

!====================================================================== 
      Subroutine Br_Wdat 
!====================================================================== 

      USE MPI 
      USE bsr_phot 

      Implicit none 
      
      Call MPI_BCAST(ncp, 1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(nwt, 1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

      if(myid.ne.0)  Allocate (WT(nwt,nhm),AK(nwt,nch),WTch(nwt))

      Call MPI_BARRIER(MPI_COMM_WORLD, ierr)

      Call MPI_BCAST(WT,nwt*nhm,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(AK,nwt*nch,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(WTch,nwt,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

      Call MPI_BARRIER(MPI_COMM_WORLD, ierr)

      End Subroutine Br_Wdat
