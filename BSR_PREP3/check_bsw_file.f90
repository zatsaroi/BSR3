!======================================================================
      Subroutine Check_bsw_file(AFW,pri)
!======================================================================
! ... check if we have bsw-file, otherwise try to call w_bsw
!----------------------------------------------------------------------
      Character(*) :: AFW
      Character(80) :: BFW, AS
      Integer :: pri

      if(Icheck_file(AFW).eq.0) then
       ilen=LEN_TRIM(AFW); BFW=AFW(1:ilen-3)//'w'
       if(Icheck_file(BFW).eq.0) then
        write(*,*) 'Can not find ',trim(AFW),' or ',trim(BFW) 
        write(pri,'(a,a,a,a)') 'Can not find ',trim(AFW),' or ',trim(BFW) 
        Stop ' '
       end if
       AS = 'w_bsw '//trim(BFW)
       i = SYSTEM(AS)
       write(*,*) 'Calling ',trim(AS),' - check w_bsw.log' 
       write(pri,'(a,a,a,a)') 'Calling ',trim(AS),' - check w_bsw.log'  
      end if

      End Subroutine Check_bsw_file
