!======================================================================
!  zf_bb_bsr - utility for the f-value calculations between states in
!  bound.nnn files
!======================================================================
!
!     INPUT:    zf_bb.inp
!
!     OUTPUT:   zf_res from BSR_DMAT program
!
!     SYSTEM CALL:  MULT3, BSR_DMAT3
!
!     Notes:  1. you would better delete all mult_bnk before running
!             2. results are appended to zf_res
!
!     zf_bb.inp containes:
!
!     atype = E1    -   or  E2, M1, E3 ...
!     klsp =        -   number of partial waves
!     index, mstate -   for each partial wave
!     msol =        -   if klsp=0
!     followed by the list of ilsp(i), mstates(i), i=1,klsp
!---------------------------------------------------------------------
      Use target
      Use channels

      Implicit real(8) (A-H,O-Z)

      Integer :: inp=5; Character(20) :: AF_inp = 'zf_bb.inp'
      Integer :: nut=7; Character(20) :: AF_tar = 'target'
      Character(1)  :: blank = ' '
      Character(1)  :: gf = 'f'
      Character(2)  :: atype = 'E1'
      Character(80) :: AS,BS,BI,BJ
      Integer, allocatable :: ilsp(:), mstate(:) 

      iarg = COMMAND_ARGUMENT_COUNT()
      if(iarg.gt.0) Call GET_COMMAND_ARGUMENT(1,AS) 
     
      if(AS.eq.'?') then
        write(*,'(/a)') 'zf_bb_bsr calculates f-values between solutions in bound.nnn'
        write(*,'(/a)') 'Call as:   zf_bb_bsr  [klsp=  msol=  atype= gf= ]'
        write(*,'(/a)') 'or provide input file "zf_bb_inp" with:'
        write(*,'(/a)') 'atype = E1 |E2,M1,... - multipole index' 
        write(*,'(/a)') 'gf    = f  |g         - output f- or gf-values'
        write(*,'(/a)') 'klsp= ...             - number of bound.nnn files to consider'
        write(*,'(/a)') 'msol =                - max. number of solutions in bound.nnn'
        write(*,'(/a)') 'followed by the list of ilsp(i), mstates(i), i=1,klsp'
        write(*,'(/a)') 'Results are appended to zf_res'
        write(*,'(/a)') 'Program makes SYSTEM CALLS to MULT3 and BSR_DMAT3'
        Stop 
      end if

! ... read target and channels information:

      Call Check_file(AF_tar)
      Open(nut,file=AF_tar)
      Call R_target (nut)
      Call R_channels(nut)
      Close(nut)

! ... input data:

      atype = 'E1'
      Allocate(ilsp(nlsp),mstate(nlsp))
      isol = 999
      Do i=1,nlsp
       ilsp(i)=i; mstate(i)=isol
      End do
      klsp = nlsp

      if(Icheck_file(AF_inp).ne.0) then
       open(inp,file=AF_inp)
       Call Read_apar(inp,'atype',atype)
       Call Read_apar(inp,'gf',gf)
       Call Read_ipar(inp,'klsp',klsp)
       Call Read_ipar(inp,'msol',msol)
       if(klsp.ne.nlsp.or.msol.eq.0) then 
        Do i=1,klsp;  read(inp,*) ilsp(i),mstate(i); End do
       end if
      end if

      Call Read_aarg('atype',atype)
      Call Read_aarg('gf',gf)
      msol=0; Call Read_iarg('msol',msol)
      if(msol.ne.0) mstate = msol 

      Read(atype,'(1x,i1)') kpol

      write(*,*) 'klsp=',klsp
      Do i=1,klsp
        write(*,*) ilsp(i),mstate(i)
      End do

      write(*,*) atype,kpol

!----------------------------------------------------------------------

      Do ii=1,klsp; i=ilsp(ii); write(BI,'(a,i3.3)') 'cfg.' ,i 
      Do jj=i,klsp; j=ilsp(jj); write(BJ,'(a,i3.3)') 'cfg.' ,j 

       if(atype(1:1).eq.'E'.and.mod(kpol,2).eq.1.and. &
          ipar(i).eq.ipar(j)) Cycle
       if(atype(1:1).eq.'E'.and.mod(kpol,2).eq.0.and. &
          ipar(i).ne.ipar(j)) Cycle
       if(atype(1:1).eq.'M'.and.mod(kpol,2).eq.0.and. &
          ipar(i).eq.ipar(j)) Cycle
       if(atype(1:1).eq.'M'.and.mod(kpol,2).eq.1.and. &
          ipar(i).ne.ipar(j)) Cycle

       if(ispar(i).ne.ispar(j)) Cycle

       kk=kpol; if(ispar(I).eq.0) kk=kpol+kpol

       if(ITRA(lpar(i),kk,lpar(j)).eq.0) Cycle

       if(lpar(i)+lpar(j).lt.kk) Cycle      

       write(AS,'(10(a,1x))') 'mult3',trim(BI),trim(BJ),atype,' >> zf_bb.out'

       Call System(AS)

       write(AS,'(a,a,1x,a,a,a,a,a,i3.3,a,i3.3,a)')         &
         'bsr_dmat3 ',trim(BI),trim(BJ),' b b ','gf=',gf, &
         ' mstate1=',mstate(ii),' mstate2=',mstate(jj),' >> zf_bb.out'
       Call System(AS)

      End do; End do 
       
      End ! program zf_bb_bsr

