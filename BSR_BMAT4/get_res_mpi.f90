!======================================================================
      Subroutine GET_res(id,is,js)
!======================================================================
! ... get confirmation from processor 'id'
!----------------------------------------------------------------------
      Use MPI

      Implicit none
      Integer :: id, is,js, status(MPI_STATUS_SIZE), ierr

      Call MPI_RECV(is,1,MPI_INTEGER, MPI_ANY_SOURCE, &
                    MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr)
      id = status(MPI_SOURCE)

      Call MPI_RECV(js,1,MPI_INTEGER, id, &
                    MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr)

      End Subroutine get_res
