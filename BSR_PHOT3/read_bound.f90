!======================================================================
      Subroutine Read_bound 
!----------------------------------------------------------------------
!     read pseudo-states energies, if any, from bound.nnn 
!----------------------------------------------------------------------      
      Use bsr_phot
      Implicit real(8) (A-H,O-Z)

      if(allocated(eps)) Deallocate(eps)

      i = LEN_TRIM(AF_b) - 3
      write(AF_b(i+1:i+3),'(i3.3)') klsp      
      i = Icheck_file(AF_b)
      if(i.eq.0) then
       khm = 1; Allocate(eps(khm)); eps = -1.d8; Return
      end if

      open(nub,file=AF_b)
      khm=0; Call Read_ipar(nub,'khm',khm)
      if(khm.eq.0) then
       khm = 1;   Allocate(eps(khm)); eps = -1.d8; Return
      end if

      Allocate(eps(khm))
      read(nub,*) 
      Do i=1,khm; read(nub,*) eps(i), eps(i); end do
      close(nub)
      
      write(pri,'(/a,a,i8/)') AF_b,' data:  khm = ', khm

      End Subroutine Read_bound
