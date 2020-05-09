!=======================================================================
      Subroutine br_dets
!=======================================================================
!     read the overlap determinants from INT_BNK (unit nub)  
!-----------------------------------------------------------------------
      Use mpi
      Use dets

      Integer :: myid,ierr
      
      Call MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierr)

      Call MPI_BCAST(ndet,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(kdet,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

      Call MPI_BCAST(ndef,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(kdef,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

      Call MPI_BARRIER(MPI_COMM_WORLD, ierr)

      if(myid.ne.0) then
       if(Allocated(KPD)) Deallocate(KPD,IPD,NPD)
       Allocate(KPD(ndet),IPD(ndet),NPD(kdet))
       if(Allocated(KPF)) Deallocate(KPF,IPF,NPF)
       Allocate(KPF(ndef),IPF(ndef),NPF(kdef))
      end if

      Call MPI_BCAST(KPD,ndet,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(IPD,ndet,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(NPD,kdet,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(KPF,ndef,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(IPF,ndef,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(NPF,kdef,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

      Call MPI_BARRIER(MPI_COMM_WORLD, ierr)

      End Subroutine br_dets




