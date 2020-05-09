!======================================================================
      Subroutine br_conf_LS
!======================================================================
!     broadcast data from  modules conf_LS
!======================================================================

      USE MPI
      USE conf_LS
      USE symt_list_LS, only: nsymt

      Implicit none
      Integer :: myid,ierr
      
      Call MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierr)

      Call MPI_BCAST(mcfg ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(ncfg ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(kcfg ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(lcfg ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

      if(myid.ne.0) then
       if(allocated(ip_state)) Deallocate(ip_state)
       if(allocated(IC_term )) Deallocate(IC_term )
       if(allocated(ip_orb  )) Deallocate(ip_orb  )
       if(allocated(ip_stat )) Deallocate(ip_stat )
       if(allocated(WC      )) Deallocate(WC      )
       Allocate(ip_state(mcfg),IC_term(mcfg),ip_orb(kcfg),WC(mcfg),ip_stat(ncfg))
       if(allocated(it_state1)) Deallocate(it_state1)
       if(allocated(it_state2)) Deallocate(it_state2)
       Allocate(it_state1(nsymt),it_state2(nsymt))
      end if
      Call MPI_BCAST(ip_state,mcfg,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(IC_term ,mcfg,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(ip_orb  ,kcfg,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

      Call MPI_BCAST(it_state1,nsymt,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(it_state2,nsymt,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

      Call MPI_BCAST(WC,mcfg,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

      Call MPI_BCAST(ip_stat ,ncfg,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

      End Subroutine br_conf_LS

