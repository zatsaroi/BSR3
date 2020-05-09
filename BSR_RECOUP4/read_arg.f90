!======================================================================
      Subroutine Read_arg 
!======================================================================
!     read arguments from command line and check default settings
!======================================================================
      Use bsr_recoup

      Implicit real(8) (A-H,O-Z)

! ... read arguments in command line:

      klsp = 0
      Call Read_iarg('klsp'  ,klsp  )
      Call Read_iarg('klsp1' ,klsp1 )
      Call Read_iarg('klsp2' ,klsp2 )
      if(klsp.gt.0) then;  klsp1=klsp; klsp2=klsp; end if 
      Call Read_iarg('pri_AC',pri_AC)
      Call Read_iarg('debug' ,debug )

! ... define the block size:

      i=len_trim(BF_mat)-2
      write(BF_mat(i:i+2),'(i3.3)') 1
      Call Check_file(BF_mat)
      open(mui,file=BF_mat,form='UNFORMATTED',action='READ')
      read(mui) ns
      close(mui)

      End Subroutine Read_arg
