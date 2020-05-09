!======================================================================
      Subroutine br_arg 
!======================================================================
!     brodcast main arguments
!----------------------------------------------------------------------
      Use mpi
      Use mult_par

      Implicit none

      Call MPI_BCAST(ktype, 1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(kpol,  1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(debug ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(eps_c,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

      End Subroutine br_arg

