!======================================================================
      Subroutine read_matrix
!======================================================================
!     Read interaction matrix from the disk
!----------------------------------------------------------------------
      Use mpi
      Use bsr_mat

      Implicit none
      Real(8) :: S, w(ns)
      Integer :: ic,jc,  ii, nn
      Integer :: status(MPI_STATUS_SIZE)

      Call MPI_BARRIER(MPI_COMM_WORLD, ierr)

      nn = ns*ns

      Do 
       if(myid.eq.0) read(nui) ic,jc
       Call MPI_BCAST(ic,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
       Call MPI_BCAST(jc,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

       if(ic.le.0) Exit

      if(ic.gt.nch.and.jc.gt.nch) then  !  pert-pert

       ic=ic-nch; jc=jc-nch; ii = ibb(ic,jc)
       if(myid.eq.0)  then
        read(nui) S
        Call MPI_SEND(S,1,MPI_DOUBLE_PRECISION,ii,0, MPI_COMM_WORLD, ierr)
       elseif(ii.ne.0) then
        Call MPI_RECV(S,1,MPI_DOUBLE_PRECISION,0, &
                      MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr)
        hbb(ii) = hbb(ii) + S 
       end if
        
      elseif(ic.gt.nch) then            !  ch-pert
       ic=ic-nch; ii = icb(jc,ic)
       if(myid.eq.0)  then
        read(nui) w
        Call MPI_SEND(w,ns,MPI_DOUBLE_PRECISION,ii,0, MPI_COMM_WORLD, ierr)
       elseif(ii.ne.0) then
        Call MPI_RECV(w,ns,MPI_DOUBLE_PRECISION,0, &
                      MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr)
        hcb(:,ii) = hcb(:,ii) + w
       end if

      else                              !  ch-ch
       ii = icc(ic,jc)
       if(myid.eq.0)  then
        read(nui) x
        Call MPI_SEND(x,nn,MPI_DOUBLE_PRECISION,ii,0, MPI_COMM_WORLD, ierr)
       elseif(ii.ne.0) then
        Call MPI_RECV(x,nn,MPI_DOUBLE_PRECISION,0, &
                      MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr)
        hcc(:,:,ii) = hcc(:,:,ii) + x
       end if

      end if 

      End do

      Call MPI_BARRIER(MPI_COMM_WORLD, ierr)

      End Subroutine read_matrix

!======================================================================
      Subroutine Skip_matrix
!======================================================================
!     skip overlap matrix from the file 'nui'
!----------------------------------------------------------------------
      Use bsr_mat

      Implicit none
      Real(8) :: S, w(ns), ww(ns*ns)
      Integer ::  ic,jc

      Do 
       read(nui) ic,jc; if(ic.le.0) Exit
       if(ic.gt.nch.and.jc.gt.nch) then  !  pert-pert
        read(nui) S
       elseif(ic.gt.nch) then            !  ch-pert
        read(nui) w
       else                              !  ch-ch
        read(nui) ww
       end if 
      End do

      End Subroutine skip_matrix

