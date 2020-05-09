!======================================================================
      Subroutine br_core
!======================================================================
!     broadcast the core parameters
!======================================================================
      USE MPI
      USE conf_LS,       only: closed, nclosd, Ecore
      USE spline_atomic, only: kclosd, EC

      Implicit none
      Integer :: myid,ierr
      
      Call MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierr)

      Call MPI_BCAST(nclosd,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(kclosd,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

      Call MPI_BCAST(CLOSED,200,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)

      Call MPI_BCAST(EC,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(Ecore,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

      End Subroutine br_core

