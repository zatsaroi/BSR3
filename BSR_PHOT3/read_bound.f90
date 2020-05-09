!======================================================================
      Subroutine Read_bound 
!----------------------------------------------------------------------
!     read pseudo-states energies, if any, from bound.nnn 
!----------------------------------------------------------------------      
      Use bsr_phot
      Implicit real(8) (A-H,O-Z)

      if(allocated(eps)) Deallocate(eps);  khm=0

      i = LEN_TRIM(AF_b) - 3
      write(AF_b(i+1:i+3),'(i3.3)') klsp      

      if(Icheck_file(AF_b).ne.0) then
       open(nub,file=AF_b)
       Call Read_ipar(nub,'khm',khm)
       if(khm.ne.0) then
        Allocate(eps(khm))
        read(nub,*) 
        Do i=1,khm; read(nub,*) eps(i), eps(i); end do
        write(pri,'(/a,a,i8/)') AF_b,' data:  khm = ', khm
       end if
       Close(nub); Return
      end if

      i = LEN_TRIM(AF_db) - 3
      write(AF_db(i+1:i+3),'(i3.3)') klsp      

      if(Icheck_file(AF_db).ne.0) then
       open(nub,file=AF_db,form='UNFORMATTED',position='APPEND')
       Backspace(nub)
       Backspace(nub) 
       read(nub) khm
       write(pri,'(/a,a,i8/)') AF_db,' data:  khm = ', khm
       if(khm.gt.0) then
        Allocate(eps(khm))
        read(nub) eps(i)
       end if
       Close(nub); Return
      end if

      if(khm.eq.0) then
       khm = 1; Allocate(eps(khm)); eps = -1.d8
      end if

      End Subroutine Read_bound
