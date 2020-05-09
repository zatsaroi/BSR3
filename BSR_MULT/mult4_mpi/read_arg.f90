!=====================================================================
      Subroutine  Read_arg
!======================================================================
!    INPUT ARGUMENTS:
!    
!     AF1  -  c-file for initial state
!     AF2  -  c-file for final state
!     AA   -  type of calculation, E1,E2,M1,M2...
! 
!-----------------------------------------------------------------------
      Use mult_par  

      Character(2) :: AA
      Integer, external :: Icheck_file

      iarg = command_argument_count()

      if(iarg.lt.2) then
       write(*,*) 'Should be at least two file-names in cc-mode:'
       write(*,*)
       Call Stop_mpi(0,0,'mult name1 name2 [E1|M1|..] AF_bnk')
      end if

      Call get_command_argument(1,AF1)
      Call get_command_argument(2,AF2)
      
! ... define the type of calculations:

      AA = 'E1'
      if(iarg.ge.3) Call get_command_argument(3,AA)
      read(AA,'(a1,i1)') ktype,kpol
      if(ktype.eq.'M'.and.kpol.eq.0) &
         Call Stop_mpi(0,0,' kpol=0 for M-type ? ')

      Call Read_iarg('debug',debug )
      Call Read_rarg('time' ,time_limit)
      Call Read_rarg('eps_c',eps_c)

!-----------------------------------------------------------------------
! ... check the existing data-bank:

      new = 1; if(Icheck_file(AF_b).eq.1) new=0
      if(new.eq.0) then
       Open(nub,file=AF_b,form='UNFORMATTED')
       rewind(nub)
       read(nub) AA,k
       if(AA.ne.ktype.or.k.ne.kpol) new=1
       if(new.eq.1) close(nub) 
      end if

      End Subroutine Read_arg
