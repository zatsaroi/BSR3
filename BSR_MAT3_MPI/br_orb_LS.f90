!======================================================================
      Subroutine br_orb_LS
!======================================================================
!     broadcast data from  module orb_LS
!======================================================================
      USE MPI
      USE orb_LS

      Implicit none
      Integer :: myid,ierr
      
      Call MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierr)

      Call MPI_BCAST(mwf,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(nwf,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    
      if(myid.ne.0) then
       if(allocated(NEF)) Deallocate(NEF,LEF,KEF,IEF,ELF)
       Allocate(NEF(mwf),LEF(mwf),KEF(mwf),IEF(mwf),ELF(mwf))
      end if

      Call MPI_BCAST(NEF,nwf,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(LEF,nwf,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(KEF,nwf,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(IEF,nwf,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(ELF,nwf*4,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)


!      Call MPI_BCAST(jort,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
!      if(myid.ne.0) then
!       if(allocated(IORT)) Deallocate(IORT)
!       Allocate(IORT(mwf,mwf))
!      end if
!      Call MPI_BCAST(IORT,mwf*mwf,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
!      if(allocated(IORT)) Deallocate(IORT)

      End Subroutine br_orb_LS

