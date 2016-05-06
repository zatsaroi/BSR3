!======================================================================
      Subroutine Read_Ddat (nu,ipri)
!----------------------------------------------------------------------
!     read D.DAT file (unit nu) 
!----------------------------------------------------------------------      
      Use bsr_phot

      Implicit none
      Integer, intent(in) :: nu,ipri
      Integer :: I
      Real(8) :: E

      read(nu) IL,IS,IP,E,ndm
      read(nu) ILI,ISI,IPI,EI

      Allocate (DKL(ndm),DKV(ndm))
      read(nu) (DKL(i),DKV(i),i=1,ndm)         ! dipole / velocity
      Close(nu)      

      if(ipri.gt.0) then
      write(ipri,'(/a/)') 'D.dat data:'
      write(ipri,'(/a,a,i3,a,i3,a,i3,a,E16.8/)')  'Initial state:', &
      '   L =', ILI,'    S =', ISI,'    parity =', IPI,'    E =',EI
      write(ipri,'(/a,a,i3,a,i3,a,i3/)')  'Final state:', &
      '   L =', IL,'    S =', IS,'    parity =', IP
      write(ipri,'(a,i6)') 'ndm =',ndm
      end if

      End Subroutine Read_Ddat
