!======================================================================
      Subroutine Record_matrix(nu)
!======================================================================

      USE mpi
      USE bsr_matrix
      USE bsr_mat, only: myid,pri

      Implicit none

      Integer, intent(in) :: nu
      Integer :: i,ib,ip,jp,ich,jch,idiag,nn,ierr
      Real(8) :: S, v(kns)
      integer :: status(MPI_STATUS_SIZE)

      Call MPI_BARRIER(MPI_COMM_WORLD, ierr)

! ... diagonal blocks:

      nn = kns*kns
      Do ib = 1, kch; i = icc(ib,ib)

       if(i.ne.0) then
        Call MPI_SEND(hcc(:,:,i),nn, MPI_DOUBLE_PRECISION, &
                      0, ib, MPI_COMM_WORLD, ierr)
       end if

       if(myid.eq.0) then
        Call MPI_RECV(x, nn, MPI_DOUBLE_PRECISION, MPI_ANY_SOURCE, &
                      MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr)
        write(nu) x
       end if

       Call MPI_BARRIER(MPI_COMM_WORLD, ierr)

      End do

      idiag = 1

! ... non-diagonal blocks:

      Do ich = 2,kch; Do jch = 1,ich-1;  i=icc(ich,jch) 

       if(i.ne.0) then
        Call MPI_SEND(hcc(:,:,i),nn, MPI_DOUBLE_PRECISION, &
                      0, ib, MPI_COMM_WORLD, ierr)
       end if

       if(myid.eq.0) then
        Call MPI_RECV(x, nn, MPI_DOUBLE_PRECISION, MPI_ANY_SOURCE, &
                      MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr)
        S = SUM(abs(x))
        if(S.ne.0.d0) then
         write(nu) ich,jch
         write(nu) x
         idiag = 0
        end if
       end if

       Call MPI_BARRIER(MPI_COMM_WORLD, ierr)

      End do; End do

! ... pertubers:

      if(kcp.gt.0) then

! ... channel-perturber rows:

       Do ich = 1,kch; Do ip = 1,kcp; i = icb(ich,ip) 

       if(i.ne.0) then
        Call MPI_SEND(hcb(:,i),kns, MPI_DOUBLE_PRECISION, &
                      0, myid, MPI_COMM_WORLD, ierr)
       end if

       if(myid.eq.0) then
        Call MPI_RECV(v, kns, MPI_DOUBLE_PRECISION, MPI_ANY_SOURCE, &
                      MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr)
        S = SUM(abs(v))
        if(S.ne.0.d0) then
        write(nu) ip+kch,ich
        write(nu) v(1:kns)
        idiag = 0
        end if
       end if

       Call MPI_BARRIER(MPI_COMM_WORLD, ierr)

       End do; End do

! ... perturter-perturber elements:

       Do ip = 1,kcp; Do jp = 1,ip; i = ibb(ip,jp)

       if(i.ne.0) then
        Call MPI_SEND(hbb(i),1, MPI_DOUBLE_PRECISION, &
                      0, myid, MPI_COMM_WORLD, ierr)
       end if

       if(myid.eq.0) then
        Call MPI_RECV(S, 1, MPI_DOUBLE_PRECISION, MPI_ANY_SOURCE, &
                      MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr)
        if(S.ne.0.d0) then
        write(nu) ip+kch,jp+kch
        write(nu) S
        if(ip.ne.jp) idiag = 0
        end if
       end if

       Call MPI_BARRIER(MPI_COMM_WORLD, ierr)

       End do; End do

      end if

! ... sign of the end

      if(myid.eq.0) write(nu) 0,idiag

      End Subroutine Record_matrix

