!======================================================================
      Subroutine br_arg
!======================================================================
      Use MPI
      Use bsr_conf

      Implicit none

      Call MPI_BCAST(kort  , 1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(max_ll, 1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(min_ll, 1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(max_LT, 1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(max_ST, 1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(debug , 1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

      Call MPI_BCAST(iread_targ,  1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

      Call MPI_BCAST(c_comp,  1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

      End Subroutine br_arg


