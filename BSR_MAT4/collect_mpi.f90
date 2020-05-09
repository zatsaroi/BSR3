!======================================================================
      Subroutine Collect_ACF
!======================================================================
      Use MPI 
      Use bsr_mat 
      Use conf_LS, only: nclosd  

      Implicit none
      Integer :: status(MPI_STATUS_SIZE)
      Integer  :: i,j,k, i1,i2,ij
      Real(8) :: c, a(mk+1)
      Character(200) :: line

! ... Summerize the ACF - matrix:

      Call MPI_BARRIER(MPI_COMM_WORLD, ierr)

      Do i = 1,nch
       Do j = 1,nch
        if(imycase(i,j).ne.0) then
         Call MPI_SEND(ACF(i,j,:),mk+1,MPI_DOUBLE_PRECISION, &
                       0, 0, MPI_COMM_WORLD, ierr)
        end if
        if(myid.eq.0) then
         Call MPI_RECV(a,mk+1,MPI_DOUBLE_PRECISION, MPI_ANY_SOURCE, &
                      MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr)
         ACF(i,j,:) = ACF(i,j,:) + a(:)
        end if

        Call MPI_BARRIER(MPI_COMM_WORLD, ierr)

       End do
      End do

      if(myid.ne.0) Return  

! ... Symmetrize the ACF - matrix:

      Do k = 0,mk
       Do i = 1,nch
        Do j = 1,i
         C = ACF(i,j,k) + ACF(j,i,k); if(i.ne.j) C=C*2.d0
         if(abs(C).lt.0.00001) C=0.d0
         ACF(i,j,k) = C; ACF(j,i,k) = C
        End do
       End do
      End do

! ... add corrections from the core screening:

       j=0; Do i=1,nclosd; j=j+2*(4*lbs(i)+2); End do
       Do i=1,nch; ACF(i,i,0)=ACF(i,i,0)+j;  End do

! ... Check coefficients for k = 0: 

       write(pri,'(/a,i3/)') 'Asymptotic coefficients: mk = ',mk
       write(pri,'(a,i3/)') 'For k=0 should be equal to 2*nelc = ', &
                                                        2*nelc
       write(pri,'(a,i3/)') 'Derivations (if > 0.00001):'
       Do i = 1,nch
        if(abs(ACF(i,i,0)-2*nelc).lt.0.00001) Cycle
        write(pri,'(i5,2F15.6)') i,ACF(i,i,0)-2*nelc
       End do

! ... print the asymptotic coefficients:

      if(pri_ac.gt.0) then
      write(pri,'(/a/)') 'Asymptotic coefficients: i,j, ACF(i,j,k)'
      line = ' '
      Do k=0,mk
       if(SUM(acf(:,:,k)).eq.0) Cycle
       write(pri,'(a,i2)') 'k = ',k
       ij = 0
       Do i=1,nch; Do j = 1,i      
        if(abs(acf(i,j,k)).lt.eps_acf) Cycle
        i1=ij*20+1; i2=i1+19 
        write(line(i1:i2),'(2i4,E12.3)') j,i,acf(i,j,k)
        ij=ij+1
        if(ij.lt.5) Cycle
        write(pri,'(a)') trim(line); ij=0
       End do; End do
       if(ij.eq.0) Cycle
       i1=1; i2=ij*20
       write(pri,'(a)') line(i1:i2)
      End do
      end if

      End Subroutine Collect_ACF


!======================================================================
      Subroutine Collect_otarg
!======================================================================
      Use MPI
      Use bsr_mat

      Implicit none
      Integer :: status(MPI_STATUS_SIZE)
      Integer :: i,j,ij
      Real(8) :: C

      Do i = 1,nch
       Do j = 1,i
        ij=(i-1)*i/2+j
        if(imycase(i,j).ne.0) then
         Call MPI_SEND(otarg(ij),1,MPI_DOUBLE_PRECISION, &
                       0, i, MPI_COMM_WORLD, ierr)
        end if
        if(myid.eq.0) then
         Call MPI_RECV(C,1,MPI_DOUBLE_PRECISION, MPI_ANY_SOURCE, &
                      MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr)
         otarg(ij) = otarg(ij) + C
        end if
       Call MPI_BARRIER(MPI_COMM_WORLD, ierr)

       End do
      End do

      End Subroutine Collect_otarg  


!======================================================================
      Subroutine Collect_htarg
!======================================================================
      Use MPI
      Use bsr_mat

      Implicit none
      Integer :: status(MPI_STATUS_SIZE)
      Integer  :: i,j,ij
      Real(8) :: C

      Do i = 1,nch
       Do j = 1,i
        ij=(i-1)*i/2+j
        if(imycase(i,j).ne.0) then
         Call MPI_SEND(htarg(ij),1,MPI_DOUBLE_PRECISION, &
                       0, i, MPI_COMM_WORLD, ierr)
        end if
        if(myid.eq.0) then
         Call MPI_RECV(C,1,MPI_DOUBLE_PRECISION, MPI_ANY_SOURCE, &
                      MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr)
         htarg(ij) = htarg(ij) + C
        end if
        Call MPI_BARRIER(MPI_COMM_WORLD, ierr)

       End do
      End do

      End Subroutine Collect_htarg  
