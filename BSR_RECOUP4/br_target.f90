!====================================================================== 
      Subroutine br_target
!====================================================================== 
      Use MPI 
      Use target

      Implicit none
      Integer :: myid,ierr,i

      Call MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierr)

      Call MPI_BCAST(ntarg,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(nelc, 1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(nz,   1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(nct,  1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(nwt,  1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(nphys,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(coupling,2,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)

      if(myid.ne.0) then; i=ntarg; Call Allocate_target(i); end if

      Call MPI_BCAST(ltarg ,ntarg,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(jtarg ,ntarg,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(istarg,ntarg,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(iptarg,ntarg,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(nctarg,ntarg,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(nwtarg,ntarg,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(etarg ,ntarg,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      i = ntarg * 20
      Call MPI_BCAST(AFT,i,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(BFT,i,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
      nct = ictarg(ntarg)
      Call MPI_BARRIER(MPI_COMM_WORLD, ierr)

      End Subroutine br_target


