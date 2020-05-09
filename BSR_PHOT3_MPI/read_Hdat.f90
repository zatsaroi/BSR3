!======================================================================
      Subroutine Read_Hdat
!----------------------------------------------------------------------
!     read H.DAT file (unit nu) and find the information for partial
!     wave IS,IL,IP
!----------------------------------------------------------------------      

      Use bsr_phot

      Implicit none

      Integer :: i,j,k,lm,iend,IL2,IS2,IP2,more
      Real(8) :: C

      i=INDEX(AF_h,'.',BACK=.TRUE.)
      write(AF,'(a,i3.3)')  AF_h(1:i),klsp
      Call Check_file(AF)
      Open(nuh, file=AF, form='UNFORMATTED')

      READ(nuh) nelc,nz,lm,km,ntarg,RA,RB
      ion = nz - nelc
 
      Allocate(Etarg(ntarg),Ltarg(ntarg),IStarg(ntarg),IPtarg(ntarg))

      READ (nuh) Etarg
      READ (nuh) Ltarg
      READ (nuh) IStarg  !,IPtarg
 
      E1  = Etarg(1)
      DO i = 1,ntarg;  Etarg(i)  = 2.d0 * (Etarg(i)-E1); END DO
     
      WRITE (pri,'(/a/)')   'H.DAT file information:' 
      WRITE (pri,'(a,i4)')  'nuclear charge             =',NZ
      WRITE (pri,'(a,i4)')  'number of target electrons =',NELC
      WRITE (pri,'(a,i4)')  'number of target states    =',NTARG
      WRITE (pri,'(a,f10.4)') 'Boundary radius =',RA
      WRITE (pri,'(a,f10.4)') 'Log. derivative =',RB
      WRITE (pri,'(a,i4   )') 'Max. multipole  =',KM
      WRITE (pri,'(a,i4   )') 'Max. small l    =',LM-1

! ... skip BUTTLE CORRECTION    ( don't use here! )
 
      READ (nuh) ((C,i=1,3),j=1,LM)
 
! ... look for the required partial wave

      iend = 0
    1 READ (nuh) IL2,IS2,IP2,NCH,NHM,MORE
      IP2 = (-1)**IP2
      if(IL2.eq.IL.and.IS2.eq.IS.and.IP2.eq.IP) iend=1

      if(iend.eq.0) then           ! skip information
                                                   
        READ (nuh) (j,i=1,ntarg)
        READ (nuh) (j,i=1,NCH)
        READ (nuh) (((C,i=1,NCH),j=1,NCH),k=1,KM)
        READ (nuh) (C,i=1,NHM)
        READ (nuh) ((C,i=1,NCH),j=1,NHM)
        if(MORE.eq.0) then
         write(pri,'(a)') &
          ' BST_PHOT cannot find required partial wave in H.DAT'
         Stop
        end if
        go to 1
 
      else                        ! read information
                       
       Allocate(NCONAT(ntarg),LCH(nch),CF(nch,nch,km),VALUE(nhm), &
                WMAT(nch,nhm))

       READ (nuh) (NCONAT(i),i=1,ntarg)
       READ (nuh) (LCH(i),i=1,nch)
       READ (nuh) (((CF(i,j,k),i=1,nch),j=1,nch),k=1,km)

       READ (nuh) (VALUE(i),i=1,nhm)
       READ (nuh) ((WMAT(i,j),i=1,nch),j=1,nhm)
       DO i = 1,NHM;   VALUE(i) = 2.0D0* (VALUE(i)-E1); END DO
  
      end if
 
      write(pri,'(/a,3i5)') 'LSP =',IL2,IS2,IP2
      write(pri,'(/a,i5 )') 'Number of channels =',nch
      write(pri,'(/a,i5/)') 'Hamiltonian matrix =',nhm

      End Subroutine Read_HDAT


!====================================================================== 
      Subroutine Br_Hdat 
!====================================================================== 

      USE MPI 
      USE bsr_phot 

      Implicit none 
 
      Call MPI_BCAST(nelc, 1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(nz,   1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(km,   1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(ntarg,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(ion,  1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(nch,  1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(nhm,  1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

      Call MPI_BCAST(RA, 1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(RB, 1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(E1, 1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

      if(myid.ne.0) &
      Allocate(Etarg(ntarg),Ltarg(ntarg),IStarg(ntarg),IPtarg(ntarg))
      if(myid.ne.0) &
      Allocate(NCONAT(ntarg),LCH(nch),CF(nch,nch,km),VALUE(nhm), &
               WMAT(nch,nhm))

      Call MPI_BARRIER(MPI_COMM_WORLD, ierr)

      Call MPI_BCAST(Ltarg, ntarg,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(IStarg,ntarg,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(IPtarg,ntarg,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(NCONAT,ntarg,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(LCH,   nch,  MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

      Call MPI_BCAST(Etarg,ntarg,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(CF,nch*nch*km,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(Value,nhm,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(WMAT,nch*nhm,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

      Call MPI_BARRIER(MPI_COMM_WORLD, ierr)

      End Subroutine Br_Hdat

