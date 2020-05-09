!======================================================================
      Subroutine Read_Ddat 
!----------------------------------------------------------------------
! ... read d.nnn file - dipole transition matrix for partial wave 'nnn'  
!----------------------------------------------------------------------      

      Use bsr_phot
      Implicit none
      Integer :: i
      Real(8) :: E

      i=INDEX(AF_d,'.',BACK=.TRUE.)
      write(AF,'(a,i3.3)') AF_d(1:i),klsp
      Call Check_file(AF)
      Open(nud, file=AF, form='UNFORMATTED', status='OLD')      

      read(nud) IL,IS,IP,E,ndm    !  final state
            
      read(nud) ILI,ISI,IPI,EI    !  initial state

      Allocate (DKL(ndm),DKV(ndm))
 
      read(nud) (DKL(i),DKV(i),i=1,ndm)    ! dipole / velocity
      Close(nud)      

      write(pri,'(/a/)') 'D.dat data:'
      write(pri,'(/a,a,i3,a,i3,a,i3,a,E16.8/)')  'Initial state:', &
      '   L =', ILI,'    S =', ISI,'    parity =', IPI,'    E =',EI
      write(pri,'(/a,a,i3,a,i3,a,i3/)')  'Final state:', &
      '   L =', IL,'    S =', IS,'    parity =', IP
      write(pri,'(a,i6)') 'ndm =',ndm

      End Subroutine Read_Ddat

!====================================================================== 
      Subroutine br_Ddat 
!====================================================================== 

      USE MPI 
      USE bsr_phot 

      Implicit none 
      
      Call MPI_BCAST(IL, 1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(IS, 1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(IP, 1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(ILI,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(ISI,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(IPI,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(EI, 1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(ndm,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

      if(myid.ne.0) Allocate (DKL(ndm),DKV(ndm))

      Call MPI_BARRIER(MPI_COMM_WORLD, ierr)

      Call MPI_BCAST(DKL,ndm,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(DKV,ndm,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

      Call MPI_BARRIER(MPI_COMM_WORLD, ierr)

      End Subroutine br_Ddat