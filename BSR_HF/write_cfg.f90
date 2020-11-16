!======================================================================
      Subroutine Write_cfg
!======================================================================
      Use bsr_hf              
      Use hf_orbitals

      Implicit none
      Character :: line*200
      Integer :: i,k
      Character(4), external :: ELF4 

! ... check conf-file:

      nconf = 0
      if(len_trim(name).ne.0) AF_conf = trim(name)//BF_conf
      Call Read_aarg('confs',AF_conf)
      if(Icheck_file(AF_conf).eq.0) Return
      Open(nuc,file=AF_conf)

! ... define number of configurations:

      rewind(nuc)
    1 read(nuc,'(a)',end=2) line  
      if(INDEX(line,'(').eq.0) go to 1
      if(line(1:1).eq.'*') goto 2
      nconf = nconf + 1
      go to 1
    2 Continue
      if(nconf.eq.0) Return

! ... output file:

      if(len_trim(name).ne.0) AF_cfg = trim(name)//BF_cfg
      Call Read_aarg('cfg',AF_cfg)
      Open(nus,file=AF_cfg)
      ncfg=0

! ... define core:

      rewind(nuc)
      read(nuc,'(/a)') core
      write(nus,'(/a)') trim(core)

! ... define ASF's:

      rewind(nuc)
    3 read(nuc,'(a)',end=4) line  
      if(INDEX(line,'(').eq.0) go to 3
      if(line(1:1).eq.'*') goto 4
      i = INDEX(line,')',BACK=.true.); no = i/8
      Do i=1,no; k = (i-1)*8+1
       Call EL4_nlk(line(k:k+3),nn(i),ln(i),in(i))
       read(line(k+5:k+6),'(i2)') iq(i)
      End do
      Call  Sum_Term
      go to 3
    4 Continue

      write(nus,'(a)') '*'

      End Subroutine Write_cfg


!======================================================================
      Subroutine Sum_Term
!======================================================================
!     exhaustion of shell-terms
!----------------------------------------------------------------------
      Use bsr_hf, only: no,ln,iq,LS

      Implicit none 
      Integer :: mt(no),nt(no), i,i1,i2,ii,IA,IL,IS
      Integer, external :: Iterm_LS

!     mt(i) - the number of term in shell i
!     nt(i) - the term inder consideration

      i1=1                     ! i1 - low  limit of shells
      i2=no                    ! i2 - high limit of shells
      Do i=i1,i2
       mt(i)=Iterm_LS(ln(i),iq(i),-1,IA,IL,IS)
      End do
      i=i1                     ! first shell under consideration
      nt(i)=1
    1 ii=Iterm_LS(ln(i),iq(i),nt(i),LS(i,1),LS(i,2),LS(i,3))
      if(i.lt.i2) then
       i=i+1; nt(i)=1; go to 1
      else
       CALL Sum_Iterm
      end if
    2 nt(i)=nt(i)+1
      if(nt(i).gt.mt(i)) then
       if(i.eq.i1) go to 3
       i=i-1; go to 2
      end if
      go to 1
    3 Return

      End Subroutine Sum_Term


!======================================================================
      Subroutine Sum_Iterm
!======================================================================
!     exhaustion of intermediate terms
!----------------------------------------------------------------------
      Use bsr_hf, only: no,LS

      Implicit none
      Integer :: LL_min(no),LL_max(no),SS_min(no),SS_max(no)
      Integer :: i,i1,i2,j1,j2

      LS(1,4)=LS(1,2)
      LS(1,5)=LS(1,3)
      if(no.eq.1) then;  CALL Output_c; Return; end if

      i1=2                 ! i1 - low  limit
      i2=no                ! i2 - high limit in array LS(...)
      i=i1
    1 j1=i-1; j2=i

      LL_min(i)=IABS(LS(j1,4)-LS(j2,2))+1
      LL_max(i)=     LS(j1,4)+LS(j2,2) -1
      SS_min(i)=IABS(LS(j1,5)-LS(j2,3))+1
      SS_max(i)=     LS(j1,5)+LS(j2,3) -1
      LS(i,4)=LL_min(i)
      LS(i,5)=SS_min(i)

    2 if(i.lt.i2) then
       i=i+1; go to 1
      else
       CALL Output_c
      end if

    3 if(LS(i,5).lt.SS_max(i)) then
       LS(i,5)=LS(i,5)+2
       go to 2
      elseif(LS(i,4).lt.LL_max(i)) then
       LS(i,4)=LS(i,4)+2
       LS(i,5)=SS_min(i)
       go to 2
      else
       if(i.le.i1) go to 4
       i=i-1; go to 3
      end if

    4 Return

      End Subroutine Sum_Iterm


!======================================================================
      Subroutine Output_c
!======================================================================
!     record the ASF in name.cfg file
!----------------------------------------------------------------------
      Use bsr_hf

      Implicit none
      Integer :: i,k
      Character(4), external :: ELF4
      Character(1), external :: AL

      if(Stotal.gt.0.and.LS(no,5).ne.Stotal) Return
      if(Ltotal.gt.0.and.LS(no,4).ne.Ltotal) Return

! ... incode the configuration:

      CONFIG = ' ';  COUPLE = ' '
      k=0
      Do i=1,no
       CONFIG(k+1:k+4) = ELF4(nn(i),ln(i),in(i))
       write(CONFIG(k+5:k+8),'(a1,i2,a1)') '(',iq(i),')'
       k=k+8
      End do
      k=0
      Do i=1,no
       write(COUPLE(k+1:k+4),'(i2,a1,i1)')  &
                             LS(i,3),AL(LS(i,2),6),LS(i,1)
       k=k+4
      End do
      Do i=2,no
       write(COUPLE(k+1:k+4),'(i2,a1,a1)')  &
                             LS(i,5),AL(LS(i,4),6),' '
       k=k+4
      End do

! ... recording:

      i = len_trim(CONFIG) 
      if(i.le.72) then
       write(nus,'(a,T65,F12.8)') trim(CONFIG), 1.d0
      else
       write(nus,'(a,F12.8)') trim(CONFIG)
      end if
      write(nus,'(a)') trim(COUPLE)

      ncfg = ncfg + 1

      End Subroutine Output_c


