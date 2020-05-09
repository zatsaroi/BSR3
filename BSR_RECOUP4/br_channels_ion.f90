!====================================================================== 
      Subroutine br_channels_ion 
!====================================================================== 
      Use MPI 
      Use channels_ion 

      Implicit none 
      Integer :: myid,ierr,i 

      Call MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierr) 

      Call MPI_BCAST(nlsp, 1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(mch,  1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

      if(myid.ne.0) Call Allocate_channels_ion

      Call MPI_BARRIER(MPI_COMM_WORLD, ierr)

      Call MPI_BCAST(ispar ,nlsp,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(lpar  ,nlsp,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(ipar  ,nlsp,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(jpar  ,nlsp,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(nch   ,nlsp,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

      Call MPI_BCAST(ncp   ,nlsp,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(nwp   ,nlsp,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

      i = nlsp*3
      Call MPI_BCAST(Tpar   ,i,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
      i = nlsp*20
      Call MPI_BCAST(AFP    ,i,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)

      i = nlsp*mch
      Call MPI_BCAST(iptar  ,i,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(lch    ,i,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(ipch   ,i,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(ipconf ,i,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(jkch   ,i,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

      i = i*4
      Call MPI_BCAST(ELC    ,i,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)

      Call MPI_BARRIER(MPI_COMM_WORLD, ierr)

      End Subroutine br_channels_ion





