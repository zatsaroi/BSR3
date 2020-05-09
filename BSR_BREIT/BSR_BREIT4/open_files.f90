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

      End Subroutine Open_c_file


!======================================================================
      Subroutine Open_int_inf
!======================================================================
!     open int_inf file if exists
!----------------------------------------------------------------------
      Use bsr_breit

      Implicit none
      Character(ma) :: AF,BF
      Integer :: i
      Integer, external :: Icheck_file, Icheck_int_inf

      AF = AF_b
      if(klsp.gt.0) then
       AF=BF_b; i=Index(BF_b,'.'); write(AF(i+1:),'(i3.3)') klsp
      end if  
      if(iname.gt.0) AF=trim(name)//'.int_inf'

      new=1; if(Icheck_file(AF).eq.1) new = 0
      Open(nub,file=AF,form='UNFORMATTED')

      ! ... check the content:

      if(new.eq.0) new = Icheck_int_inf(nub)

      End Subroutine Open_int_inf      


!======================================================================
      Integer Function Icheck_int_inf(nu)
!======================================================================
!     check int_inf file 
!----------------------------------------------------------------------
      Implicit none
      Integer, intent(in) :: nu
      Integer :: n,l,k
      Integer(1) :: i1
      Integer(8) :: i8
      
      Icheck_int_inf = 0
      rewind(nu)

! ... symc data

      read(nu,end=10,err=10) n,l
      read(nu,end=10,err=10) n
      read(nu,end=10,err=10) n
      read(nu,end=10,err=10) n
      read(nu,end=10,err=10) n
      read(nu,end=10,err=10) n
      read(nu,end=10,err=10) n

! ... symt data

      read(nu,end=10,err=10) n,l
      read(nu,end=10,err=10) n
      read(nu,end=10,err=10) n
      read(nu,end=10,err=10) n

! ... it_oper

      read(nu,end=10,err=10) i8
      read(nu,end=10,err=10) i1

! ... det-data

      read(nu,end=10,err=10) n,l,k
      read(nu,end=10,err=10) n
      read(nu,end=10,err=10) n
      read(nu,end=10,err=10) n
      read(nu,end=10,err=10) n

! ... detf-data

      read(nu,end=10,err=10) n,l,k
      read(nu,end=10,err=10) n
      read(nu,end=10,err=10) n
      read(nu,end=10,err=10) n
      read(nu,end=10,err=10) n

      rewind(nu)
      Return

   10 Icheck_int_inf = 1
      rewind(nu)

      End Function Icheck_int_inf



!======================================================================
      Subroutine Open_int_int
!======================================================================
!     open int_inf file if exists
!----------------------------------------------------------------------
      Use bsr_breit

      Implicit none
      Character(ma) :: AF,BF
      Integer :: i

      AF = AF_r
      if(klsp.gt.0) then
       AF=BF_r; i=Index(BF_r,'.'); write(AF(i+1:),'(i3.3)') klsp
      end if  
      if(iname.gt.0) AF=trim(name)//'.int_int'

      Open(nur,file=AF,form='UNFORMATTED',position='APPEND')
      if(new.eq.1) rewind(nur)

      End Subroutine Open_int_int      


!======================================================================
      Subroutine Open_det_exp
!======================================================================
!     open file with det. expansions
!     
!     mode with existing det-file is out for the moment ...  ???
!----------------------------------------------------------------------
      Use bsr_breit
      Use term_exp,     only: ic_case
      Use symc_list_LS, only: nsymc

      Implicit none
      Integer :: i
      Integer, external :: Icheck_file
      Character(ma) :: AF,BF

      AF = AF_d
      if(klsp.gt.0) write(AF,'(a,a,i3.3)') trim(AF_d),'.',klsp
      Open(nud,file=AF,form='UNFORMATTED')
      Call Pre_det_expn(nud,mkt,mkdt,pri) 

      End Subroutine Open_det_exp      
