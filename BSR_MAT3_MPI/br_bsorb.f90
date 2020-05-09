!======================================================================
      Subroutine br_bsorb
!======================================================================
!     broadcast data from  module spline_orbitals
!======================================================================
      USE MPI
      USE spline_orbitals
      USE spline_param, only: ns

      Implicit none
      Integer :: myid,ierr
      
      Call MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierr)

      Call MPI_BCAST(nbf,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      if(nbf.eq.0) Return

      if(myid.ne.0) then
       if(Allocated(nbs)) & 
          Deallocate (nbs,lbs,kbs,mbs,iech,ebs,IBORT,PBS,QBS,OBS)
       mbf = nbf
       Allocate(nbs(mbf),lbs(mbf),kbs(mbf),ebs(mbf),mbs(1:mbf), &
                iech(1:mbf),IBORT(1:mbf,1:mbf),PBS(1:ns,1:mbf), &
                QBS(1:ns,1:mbf),OBS(1:mbf,1:mbf))
      end if

      Call MPI_BCAST(nbs ,nbf,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(lbs ,nbf,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(kbs ,nbf,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(mbs ,nbf,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(iech,nbf,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

      Call MPI_BCAST(ebs,nbf*4,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)

      Call MPI_BCAST(PBS,nbf*ns,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(QBS,nbf*ns,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

      Call MPI_BCAST(OBS,nbf*nbf,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

      Call MPI_BCAST(IBORT,nbf*nbf,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

      End Subroutine br_bsorb

