!======================================================================
      Subroutine br_arg 
!======================================================================
!     brodcast main arguments
!----------------------------------------------------------------------

      Use mpi
      Use bsr_breit 

      Implicit none

      Call MPI_BCAST(klsp,  1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(klsp1, 1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(klsp2, 1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(mk,    1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(debug ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

      Call MPI_BCAST(ioper,noper,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

      Call MPI_BCAST(eps_c,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)


      End Subroutine br_arg

