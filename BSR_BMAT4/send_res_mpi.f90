!======================================================================
      Subroutine send_res(ic,jc)
!======================================================================
!     send confirmation from the given processor
!======================================================================
      Use MPI
      
      Implicit none
      Integer :: ic,jc
      Integer :: myid,ierr
      
      Call MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierr)

      Call MPI_SEND(ic,   1,MPI_INTEGER,0,myid,MPI_COMM_WORLD,ierr)
      Call MPI_SEND(jc,   1,MPI_INTEGER,0,myid,MPI_COMM_WORLD,ierr)

      End Subroutine send_res
