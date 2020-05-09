!======================================================================
      Subroutine Collect_ACF
!======================================================================
      Use MPI 
      Use bsr_recoup

      Implicit none
      Integer :: status(MPI_STATUS_SIZE)
      Integer  :: i,j,k, i1,i2,ij
      Real(8) :: c

! ... Summerize the ACF - matrix:

      Call MPI_BARRIER(MPI_COMM_WORLD, ierr)

       if(allocated(bcf)) deallocate(bcf)      
       if(myid.eq.0) then
         Allocate(bcf(nch,nch,0:mk)); bcf = 0.d0
       end if

       Call MPI_REDUCE(ACF,BCF,nch*nch*(mk+1),MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)

      if(myid.ne.0) Return  

! ... Symmetrize the ACF - matrix:

      ACF = BCF

      Do i = 1,nch
       Do j = 1,i
        ACF(j,i,:) = ACF(i,j,:)
       End do
      End do

      Do i = 1,nch
       Do j = 1,nch
        Do k = 0,mk
         if(abs(ACF(i,j,k)).lt.eps_ACF) ACF(i,j,K)=0.d0
        End do
       End do
      End do

! ... print the asymptotic coefficients:

      write(pri,'(/a,i3/)') 'Asymptotic coefficients: mk = ',mk
      write(pri,'(a,i3/)') 'For k=0 should be equal to 2*nelc = ', &
                                                        2*nelc
      write(pri,'(a,i3/)') 'Derivations (if > 0.00001):'
      Do i = 1,nch
       if(abs(ACF(i,i,0)-2*nelc).lt.0.00001) Cycle
       write(pri,'(i5,2F15.6)') i,ACF(i,i,0)-2*nelc
      End do

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


