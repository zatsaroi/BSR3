!======================================================================
      Subroutine Open_c_file
!======================================================================
!     open c-file according klsp index
!----------------------------------------------------------------------
      Use bsr_breit

      Implicit none
      Integer :: i,iarg
      Character(ma) :: AF,BF

      AF = AF_c
      if(klsp.gt.0) then
       AF=BF_c; i=Index(AF,'.'); write(AF(i+1:),'(i3.3)') klsp
      else
       iarg = command_argument_count()
       if(iarg.gt.0) then
        Call GET_COMMAND_ARGUMENT(1,BF); i=INDEX(BF,'=')
        if(i.eq.0) then 
         i=LEN_TRIM(BF)
         AF=BF; if(BF(i-1:i).ne.'.c') AF=trim(BF)//'.c'
         i=LEN_TRIM(AF)-2; name = AF(1:i); iname=i
        end if
       end if
      end if  

      Call Check_file(AF)
      Open(nuc,file=AF)
      CF_c = AF

      End Subroutine Open_c_file


!======================================================================
      Subroutine Open_int_list
!======================================================================
!     just check name for int_list file 
!----------------------------------------------------------------------
      Use bsr_breit

      Implicit none
      Character(ma) :: AF,BF
      Integer :: i
      Integer, external :: Icheck_file

      i=Index(BF_b,'.'); write(BF_b(i+1:),'(i3.3)') klsp

      End  Subroutine Open_int_list      


!======================================================================
      Subroutine Open_det_exp
!======================================================================
!     open (close) files with det. expansions
!----------------------------------------------------------------------
      Use bsr_breit
      Use term_exp

      Implicit none
      Integer :: i,j
      Integer, external :: Icheck_file
      Character(ma) :: AF,BF
      Integer(8) :: ij
      Integer(8), external :: DEF_ij8

      AF = AF_d
      if(klsp.gt.0) then
       write(BF,'(a,a,i3.3)') trim(AF),'.',klsp
       AF = BF
      end if  

      ic_case = 0
      if(Icheck_file(AF).eq.0) then

       write(pri,'(/a)') 'Prepare det. expansions:'
       write(pri,'(/a,i10 )') 'mkt  =',mkt 
       write(pri,'( a,i10 )') 'mkdt =',mkdt 

       Open(nud,file=AF,form='UNFORMATTED')

       AF = AF_a
       if(klsp.gt.0) then
        write(BF,'(a,a,i3.3)') trim(AF),'.',klsp
        AF = BF
       end if  
       Open(nux,file=AF)

       Call Pre_det_exp(nud,mkt,mkdt,nux)

       ! remove det_done if any:
       AF = AF_r
       if(klsp.gt.0) then
        write(BF,'(a,a,i3.3)') trim(AF),'.',klsp
        AF = BF
       end if 
       Open(nur,file=AF,form='UNFORMATTED')
       Close(nur,status='DELETE')

      else

       Open(nud,file=AF,form='UNFORMATTED',position='APPEND')
       BACKSPACE(nud)
       read(nud) ic_case
       write(pri,'(/a)') 'Use old det. expansions:'
       write(pri,'(/a,i10 )') 'ic_case =',ic_case 

      end if
      rewind(nud)

      End Subroutine Open_det_exp      


!======================================================================
      Subroutine Open_det_done
!======================================================================
!     open (close) files with det. expansions
!----------------------------------------------------------------------
      Use bsr_breit

      Implicit none
      Integer :: i,j
      Integer, external :: Icheck_file
      Character(ma) :: AF,BF

      AF = AF_r
      if(klsp.gt.0) then
       write(BF,'(a,a,i3.3)') trim(AF),'.',klsp
       AF = BF
      end if 

      if(Icheck_file(AF).eq.0) then
       Open(nur,file=AF,form='UNFORMATTED')
       Call Define_is_done
      else
       Open(nur,file=AF,form='UNFORMATTED')
       Call Read_is_done
      end if
      rewind(nur)

      End Subroutine Open_det_done      
