!======================================================================
      Subroutine Read_arg
!======================================================================
!     determine the main input parameters
!----------------------------------------------------------------------
      Use bsr_dmat

      Implicit none
      Integer :: iarg,i1,i2
      Integer :: IARGC

      iarg = IARGC()
      AF = '?'; if(iarg.gt.0) Call GETARG(1,AF)       
      
      if(iarg.lt.4.or.AF.eq.'?') then
       write(*,*)
       write(*,*) 'BSR_DMAT: number of arguments should be => 4:'
       write(*,*)
       write(*,*) '1 - file name for initial-state c-file'
       write(*,*) '2 - file name for final-state c-file'
       write(*,*) '3 - structure mode for initial state: c,l,j,b'
       write(*,*) '4 - structure mode for final state: c,l,j,b,p,q,d'
       write(*,*)
       write(*,*) 'all other arguments have a key structure: arg=value'
       write(*,*)
       write(*,*) 'example:  bsr_dmat 1.c cfg.002 l b gf=f istate2=1'
       write(*,*)
       write(*,*) 'see bsr_dmat.f90 for short description'
       write(*,*)
       Stop 'stop BSR_DMAT'
      end if

      Call GETARG(1,name1)       
      Call GETARG(2,name2) 

      Call GETARG(3,ctype1)
      Call GETARG(4,ctype2)

      jout=0; if(ctype2.eq.'p') jout=1
      if(ctype2.eq.'q') then; jout=2; ctype2='p'; end if       
      if(ctype2.eq.'d') then; jout=3; ctype2='p'; end if       
	  
      i1 = LEN_TRIM(name1)
      if(name1(i1:i1).eq.'c') then; ilsp1=0
      else; ALS1=name1(i1-2:i1); read(ALS1,'(i3)') ilsp1; end if
      if(ilsp1.eq.0.and.ctype1.eq.'b') Stop 'ctype1 inconsistent'

      i2 = LEN_TRIM(name2)
      if(name2(i2:i2).eq.'c') then; ilsp2=0
      else; ALS2=name2(i2-2:i2); read(ALS2,'(i3)') ilsp2; end if
      if(ilsp2.eq.0.and.ctype2.eq.'b') Stop 'ctype2 inconsistent'
      if(ilsp2.eq.0.and.ctype2.eq.'p') Stop 'ctype2 inconsistent'

! ... read parameters from arguments if any

      Call Read_aarg('gf',GF)

      Call Read_iarg('mstate1' ,mstate1)
      Call Read_iarg('mstate2' ,mstate2)
      Call Read_iarg('istate1' ,istate1)
      Call Read_iarg('istate2' ,istate2)

      Call Read_iarg('ialfa'   ,ialfa  )

      Call Read_iarg('debug'   ,debug  )

      End Subroutine Read_arg


