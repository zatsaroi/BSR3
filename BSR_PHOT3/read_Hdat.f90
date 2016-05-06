!======================================================================
      Subroutine Read_Hdat(nu,ipri)
!----------------------------------------------------------------------
!     Read H.DAT file (unit nu) and find the information for partial
!     wave IS,IL,IP
!----------------------------------------------------------------------      
      Use bsr_phot

      Implicit real(8) (A-H,O-Z)
      Integer, Intent(in) :: nu,ipri

      Read(nu) nelc,nz,lm,km,ntarg,RA,RB
      ion = nz - nelc
 
      Allocate(Etarg(ntarg),Ltarg(ntarg),IStarg(ntarg),IPtarg(ntarg))

      Read (nu) Etarg
      Read (nu) Ltarg
      Read (nu) IStarg  !,IPtarg
 
      E1  = Etarg(1)
      Do i = 1,ntarg
       Etarg(i)  = 2.d0 * (Etarg(i)-E1)
      End do
     
      if(ipri.gt.0) then 
       write (ipri,'(/a/)')   'H.DAT file information:' 
       write (ipri,'(a,i4)')  'nuclear charge             =',NZ
       write (ipri,'(a,i4)')  'number of target electrons =',NELC
       write (ipri,'(a,i4)')  'number of target states    =',NTARG
       write (ipri,'(a,f10.4)') 'Boundary radius =',RA
       write (ipri,'(a,f10.4)') 'Log. derivative =',RB
       write (ipri,'(a,i4   )') 'Max. multipole  =',KM
       write (ipri,'(a,i4   )') 'Max. small l    =',LM-1
      end if

! ... skip BUTTLE CORRECTION    ( don't use here! )
 
      read (nu) ((C,i=1,3),j=1,LM)
 
! ... look for the required partial wave

      iend = 0
    1 read (nu) IL2,IS2,IP2,NCH,NHM,MORE
      IP2 = (-1)**IP2
      if(IL2.eq.IL.and.IS2.eq.IS.and.IP2.eq.IP) iend=1

      if(iend.eq.0) then           ! skip information
                                                   
        read (nu) (j,i=1,ntarg)
        read (nu) (j,i=1,NCH)
        read (nu) (((C,i=1,NCH),j=1,NCH),k=1,KM)
        read (nu) (C,i=1,NHM)
        read (nu) ((C,i=1,NCH),j=1,NHM)
        if(MORE.eq.0) then
         write(ipri,'(a)') &
          ' BST_PHOT cannot find required partial wave in H.DAT'
         Stop
        end if
        go to 1
 
      else                        ! Read information
                       
       Allocate(NCONAT(ntarg),LCH(nch),CF(nch,nch,km),VALUE(nhm), &
                WMAT(nch,nhm))

       read (nu) (NCONAT(i),i=1,ntarg)
       read (nu) (LCH(i),i=1,nch)
       read (nu) (((CF(i,j,k),i=1,nch),j=1,nch),k=1,km)

       read (nu) (VALUE(i),i=1,nhm)
       read (nu) ((WMAT(i,j),i=1,nch),j=1,nhm)
       Do i = 1,NHM
         VALUE(i) = 2.0D0* (VALUE(i)-E1)
       End do
  
      end if
 
      if(ipri.gt.0) then
       write(ipri,'(/a,3i5)') 'LSP =',IL2,IS2,IP2
       write(ipri,'(/a,i5 )') 'Number of channels =',nch
       write(ipri,'(/a,i5/)') 'Hamiltonian matrix =',nhm
      end if

      End Subroutine Read_HDAT

