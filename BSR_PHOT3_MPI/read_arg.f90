!======================================================================
      Subroutine Read_arg(nu)
!======================================================================
!     read arguments from file unit 'nu' and command line
!----------------------------------------------------------------------
      Use bsr_phot

      Implicit None
      Integer :: nu, i,i1,i2,ie, ne
      Real(8) :: EKK
      Real(8), allocatable :: ee(:)

      rewind(nu); Elow=0.d0; Ehigh=0.d0; Estep=0.d0 
    1 read(nu,'(a)',end=2) AS
      i=INDEX(AS,'='); if(i.lt.2) go to 1
      i1=i+1; i2=INDEX(AS,'!')-1; if(i2.lt.i1) i2=LEN_TRIM(AS)
      Select Case(AS(1:i-1))
       Case('klsp'  );  read(AS(i1:i2),*) klsp
       Case('iauto' );  read(AS(i1:i2),*) iauto
       Case('mfgi ' );  read(AS(i1:i2),*) mfgi
       Case('AC'    );  read(AS(i1:i2),*) AC
       Case('awt'   );  read(AS(i1:i2),*) AWT
       Case('e_exp' );  read(AS(i1:i2),*) e_exp
       Case('DR'    );  read(AS(i1:i2),*) DR
       Case('ibug'  );  read(AS(i1:i2),*) ibug
       Case('elow'  );  read(AS(i1:i2),*,end=1) ELOW
       Case('estep' );  read(AS(i1:i2),*,end=1) ESTEP
       Case('ehigh' );  read(AS(i1:i2),*,end=1) EHIGH
       Case('nwt'   );  read(AS(i1:i2),*) nwt
       Case('ikm'   );  read(AS(i1:i2),*) ikm
       Case('ath'   );  read(AS(i1:i2),*) athreshold
       Case('bth'   );  read(AS(i1:i2),*) bthreshold
      End Select
      go to 1
    2 Continue   

      Call Read_iarg('klsp' ,klsp )
      Call Read_iarg('iauto',iauto)
      Call Read_iarg('mfgi' ,mfgi )
      Call Read_iarg('ibug' ,ibug )
      Call Read_iarg('nwt'  ,nwt  )
      Call Read_iarg('ikm'  ,ikm  )
      Call Read_rarg('AC'   ,AC   )
      Call Read_rarg('DR'   ,DR   )
      Call Read_rarg('AWT'  ,AWT  )
      Call Read_rarg('e_exp',e_exp)
      Call Read_rarg('elow' ,elow(1) )
      Call Read_rarg('estep',estep(1))
      Call Read_rarg('ehigh',ehigh(1))

      Call Read_rarg('ath'  ,athreshold)
      Call Read_rarg('bth'  ,bthreshold)

!----------------------------------------------------------------------
! ... define energies:

      i=INDEX(AF_ph,'.',BACK=.TRUE.)
      write(AF,'(a,i3.3)') AF_ph(1:i),klsp
      open(iph,file=AF)
      ne = 0
    3 read(iph,*,end=4) ekk
      ne = ne + 1   
      go to 3
    4 Allocate(ee(0:ne)); ee = -1.0d0
      rewind(iph)
      Do i = 1,ne; read(iph,*) ee(i); End do
      rewind(iph) 
     
      me = 0 
      Do i = 1,nrange
       if(Elow(i).le.0.d0) Cycle
       if(Ehigh(i).lt.Elow(i)) Cycle
       if(Estep(i).le.0.d0) Cycle
      Do EKK = Elow(i),Ehigh(i)+0.5*Estep(i),Estep(i)
       if(minval(abs(ee - ekk)).lt.5.d-8) Cycle
       me = me + 1
      End do; End do

      write(*,'(a,i6,a)') 'me =', me,' - number of energies'
      Allocate(EK(0:me),IEK(0:me)) 
      if(me.eq.0) Return
      me = 0; ie=0
      Do i = 1,nrange
       if(Elow(i).le.0.d0) Cycle
       if(Ehigh(i).lt.Elow(i)) Cycle
       if(Estep(i).le.0.d0) Cycle
      Do EKK = Elow(i),Ehigh(i)+0.5*Estep(i),Estep(i)
       if(minval(abs(ee - ekk)).lt.5.d-8) Cycle
       me = me + 1
       EK(me) = EKK
       ie = ie + 1; if(ie.ge.nprocs) ie=0
       IEK(me) = ie
       write(*,'(f16.8,i10)') EK(me),IEK(me)
      End do; End do

      End Subroutine Read_arg

!======================================================================
      Subroutine br_arg
!======================================================================

      USE MPI
      Use bsr_phot

      Implicit none
      
      Call MPI_BCAST(klsp,  1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(iauto, 1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(mfgi,  1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(nwt,   1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(ikm,   1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(ibug,  1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      
      Call MPI_BCAST(AC,   1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(DR,   1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(AWT,  1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(e_exp,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

      Call MPI_BCAST(athreshold,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(bthreshold,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

      Call MPI_BCAST(me,   1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

      if(myid.ne.0) Allocate(EK(me),IEK(me))

      Call MPI_BARRIER(MPI_COMM_WORLD, ierr)

      Call MPI_BCAST(iek,me,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(ek,me,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

      End Subroutine br_arg


