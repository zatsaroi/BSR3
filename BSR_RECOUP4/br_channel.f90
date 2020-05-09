!====================================================================== 
      Subroutine br_channel 
!====================================================================== 
      Use MPI 
      Use channel 

      Implicit none 
      Integer :: myid,ierr,i 

      Call MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierr) 

      Call MPI_BCAST(lpar, 1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(ipar, 1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(jpar, 1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(ispar,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(ncp,  1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(npert,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(nwp,  1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(nch,  1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(AFP, 20,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)

      if(myid.ne.0) then; i=nch; Call Allocate_channel(i); end if

      Call MPI_BARRIER(MPI_COMM_WORLD, ierr)

      Call MPI_BCAST(lch   ,nch,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(iptar ,nch,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(ipch  ,nch,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(jkch  ,nch,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(ipconf,nch,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      i = nch*4
      Call MPI_BCAST(ELC,i,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)

      if(npert.gt.0) then
       if(myid.ne.0) then; i=npert; Call Allocate_pert(i); end if
       Call MPI_BCAST(ippert(0:npert),npert+1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      end if

      Call MPI_BARRIER(MPI_COMM_WORLD, ierr)

      End Subroutine br_channel




