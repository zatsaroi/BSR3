!======================================================================
      Subroutine br_phys_orb
!======================================================================

      USE MPI
      USE phys_orb_LS
      USE target, only: ntarg

      Implicit none
      
      Integer :: myid,ierr

      Call MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierr)

      Call MPI_BCAST(nphys_orb,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      Call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      if(nphys_orb.eq.0) Return

      if(myid.ne.0) then
       if(allocated(ip_tar)) Deallocate(ip_tar,ip_phy,ip_sub,S_orb,jp_sub)
       Allocate(ip_tar(ntarg),ip_phy(nphys_orb),ip_sub(nphys_orb), &
                S_orb(nphys_orb),jp_sub(nphys_orb))       
      end if

      Call MPI_BCAST(ip_tar,ntarg,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(ip_phy,nphys_orb,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(ip_sub,nphys_orb,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(S_orb,nphys_orb,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

      Call MPI_BCAST(nphys_sub,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(jp_sub,nphys_orb,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(npert_sub,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

      Call MPI_BARRIER(MPI_COMM_WORLD, ierr)

      End Subroutine br_phys_orb


