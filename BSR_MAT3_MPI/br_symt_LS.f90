!======================================================================
      Subroutine br_symt_LS
!======================================================================
!     broadcast data from  module symt_list_LS
!======================================================================

      USE MPI
      USE symt_list_LS

      Implicit none
      Integer :: myid,ierr
      
      Call MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierr)

! ... configuration synnetrs:

      Call MPI_BCAST(msymt,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(nsymt,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(ksymt,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(lsymt,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

      if(myid.ne.0) then
       if(allocated(IT_conf)) Deallocate(IT_conf)
       if(allocated(ip_term)) Deallocate(ip_term)
       if(allocated(LS_term)) Deallocate(LS_term)
       Allocate(IT_conf(msymt),ip_term(msymt),LS_term(ksymt,5))    
      end if

      Call MPI_BCAST(IT_conf,msymt,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(ip_term,msymt,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

      Call MPI_BCAST(LS_term,ksymt*5,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

      End Subroutine br_symt_LS

