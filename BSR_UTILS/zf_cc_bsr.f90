!======================================================================
!     zf_cc_bsr - utility for the f-value calculations 
!                 between set of c-files
!======================================================================
!
!     INPUT:    zf_cc_bsr.inp
!
!     OUTPUT:   zf_res from BSR_DMAT3 program
!
!     SYSTEM CALL:  MULT3, BSR_DMAT3
!
!     Notes:  1. you would better delete all mult_bnk before running
!             2. results are appended to zf_res
!
!     zf_inp.inp containes:
!
!     atype = E1           or  E2, M1, ...
!     nfiles = ...
!     list of c-files with energies and exp.coefficients with one J
!
!---------------------------------------------------------------------
      Implicit real(8) (A-H,O-Z)
      Integer :: inp=5; Character(20) :: AF_inp = 'zf_cc_bsr.inp'
      Integer :: nuc=1
      Character(1) :: blank = ' '
      Character(2) :: atype = 'E1'
      Character(4) :: cc = ' c c'
      Character(4) :: gf = 'f'
      Character(80) :: AS,AI,AJ, AF, param
      Character(5) ::  AC='a.c', AW ='a.bsw'
      Character(80), Allocatable :: files(:) 
      Integer, allocatable :: ILT(:),IST(:),parity(:)

      iarg = COMMAND_ARGUMENT_COUNT()
      if(iarg.gt.0) Call GET_COMMAND_ARGUMENT(1,AF) 
     
      if(AF.eq.'?') then
        write(*,'(/a)') 'zf_cc_bsr  calculates f-values between set of (c+bsw)-files'
        write(*,'(/a)') 'Call as:   zf_cc_bsr  inp=name, default: zf_cc_bsr.inp'
        write(*,'(/a)') 'input file contains:'
        write(*,'(/a)') 'atype = E1 |E2,M1,... - multipole index' 
        write(*,'(/a)') 'gf    = f  |g         - output f- or gf-values'
        write(*,'(/a)') 'param =               - additional parameters for bsr_dmat3'
        write(*,'(/a)') 'nfiles= ...           - number of c-files'
        write(*,'(/a)') 'followed by list of c-files'
        write(*,'(/a)') 'Results are appended to zf_res'
        write(*,'(/a)') 'Program makes SYSTEM CALLS to MULT3 and BSR_DMAT3'
        Stop 
      end if

      Call Read_aarg('inp',AF_inp)
      Call Check_file(AF_inp)
      open(inp,file=AF_inp)

      Call Read_apar(inp,'atype',atype)
      Call Read_aarg('atype',atype)
      Read(atype,'(1x,i1)') kpol

      Call Read_apar(inp,'gf',gf)                                     
      Call Read_aarg('gf',gf)

      param = ' '
      Call Read_apar(inp,'param',param)
      Call Read_aarg('param',param)

      Call Read_ipar(inp,'nfiles',nfiles)
      if(nfiles.le.0) Stop 'nfiles = 0'

      Allocate(files(nfiles),ILT(nfiles),IST(nfiles),parity(nfiles))

      Do i=1,nfiles
       read(inp,*) files(i)
       Call Check_file(files(i))
       open(nuc,file=files(i))
       Call Def_term(nuc,ILT(I),IST(i),parity(i)) 
       close(nuc)
      End do
!----------------------------------------------------------------------

      Do i=1,nfiles-1; AI=files(i); ii=LEN_TRIM(AI) 
      Do j=i+1,nfiles; AJ=files(j); jj=LEN_TRIM(AJ)

       if(AI.eq.AJ) then
        AS = 'copy '//AJ(1:jj)//blank//AC 
        Call System(AS)
        AF = AJ(1:jj-1)//'bsw'
        AS = 'copy '//AF(1:jj+2)//blank//AW 
        Call System(AS)
        AJ=AC; jj=LEN_TRIM(AJ)
       end if 

       if(atype(1:1).eq.'E'.and.mod(kpol,2).eq.1.and. &
          parity(i).eq.parity(j)) Cycle
       if(atype(1:1).eq.'M'.and.mod(kpol,2).eq.0.and. &
          parity(i).eq.parity(j)) Cycle
       if(atype(1:1).eq.'E'.and.mod(kpol,2).eq.0.and. &
          parity(i).ne.parity(j)) Cycle
       if(atype(1:1).eq.'M'.and.mod(kpol,2).eq.1.and. &
          parity(i).ne.parity(j)) Cycle

       if(IST(i).ne.IST(j)) Cycle

       kk=kpol+kpol+1; if(IST(I).eq.0) kk=kpol+kpol
       if(ITRA(ILT(i),kk,ILT(j)).eq.0) Cycle
      
       AS = 'mult3 '//AI(1:ii)//blank//AJ(1:jj)//blank//atype//' >> zf_cc_bsr.out'
       Call System(AS)

       AS = 'bsr_dmat3 '//AI(1:ii)//blank//AJ(1:jj)//cc//' gf='//gf//' '//  &
                          trim(param)//' >> zf_cc_bsr.out'
       Call System(AS)

      End do
      End do 
       
      End ! program zf_cc_bsr

