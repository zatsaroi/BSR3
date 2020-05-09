!======================================================================
      Subroutine br_symc_LS
!======================================================================
!     broadcast data from  module symc_list_LS
!======================================================================

      USE MPI
      USE symc_List_LS

      Implicit none
      Integer :: myid,ierr
      
      Call MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierr)

! ... configuration synnetrs:

      Call MPI_BCAST(msymc,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(nsymc,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(ksymc,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(lsymc,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

      if(myid.ne.0) then
       if(allocated(LT_conf)) Deallocate(LT_conf)
       if(allocated(ST_conf)) Deallocate(ST_conf)
       if(allocated(no_conf)) Deallocate(no_conf)
       if(allocated(ip_conf)) Deallocate(ip_conf)
       if(allocated(iq_conf)) Deallocate(iq_conf)
       if(allocated(ln_conf)) Deallocate(ln_conf)
       Allocate(LT_conf(msymc),ST_conf(msymc),no_conf(msymc), &
                ip_conf(msymc),iq_conf(ksymc),ln_conf(ksymc))
      end if

      Call MPI_BCAST(LT_conf,msymc,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(ST_conf,msymc,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(no_conf,msymc,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(ip_conf,msymc,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(iq_conf,ksymc,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(ln_conf,ksymc,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

      End Subroutine br_symc_LS

