!=====================================================================
!     UTIL      M U L T _  T A B                          version.3
!
!               C O P Y R I G H T -- 2011
!
!     Written by:   Oleg Zatsarinny
!                   free shooter
!                   email: oleg_zoi@yahoo.com
!                         
!======================================================================
!     print,  to file mult.tab, the angular coefficients for two 
!     selected states in input c-files according to MULT.BNK file
!----------------------------------------------------------------------
      Use mult_tab
      Use conf_LS
      Use orb_LS
      Use det_list
      Use def_list

      Implicit real(8) (A-H,O-Z)

      Integer,    allocatable :: IT_conf(:)
      Integer(1), allocatable :: IT_oper(:)

      Integer, dimension(me) :: IPN1,IPN2, MP,MP1,MP2

      Integer, external :: IDET_SIMP, Iadd_int, Nadd_det, Nadd_def 

      Character(26) ::  AI1, AI2  
      Character(40) ::  BI1, BI2 
      Character(80) ::  AS

!----------------------------------------------------------------------
! ... read arguments from command line:

      iarg = IARGC()
      
      if(iarg.lt.2) then
       write(*,*)
       write(*,*) 'mult_tab prints the angular coefficients for '
       write(*,*) 'selected transitions. Call as:'
       write(*,*)
       write(*,*) 'mul_tab c1=... c2=... [bnk=... tab=... ic=... jc=... iort=...]'
       write(*,*)
       write(*,*) 'c1  - file name for initial-state c-file'
       write(*,*) 'c2  - file name for final-state c-file'
       write(*,*) 'bnk - file name for mult data_bank [mult_bnk.E1]'
       write(*,*) 'tab - file name for output [mult_tab]'
       write(*,*) 'ic  - selected configuration in initial c-file [0]'
       write(*,*) 'jc  - selected configuration in final c-file  [0]'
       write(*,*) 'jort- orthogonality mode (-1|1|0)'
       write(*,*) '    -1  - full orthogonality'
       write(*,*) '     0  - full non-orthogonality'
       write(*,*) '    +1  - partial orthogonality, default'
       write(*,*)
       write(*,*) 'example: mult_tab 1.c 2.c - all dipole transitions'
       Stop 
      end if

      Call Read_aarg('c1',AF1);      Call Check_file(AF1)
      Call Read_aarg('c2',AF2);      Call Check_file(AF2)
      Call Read_aarg('bnk',AF_bnk);  Call Check_file(AF_bnk)
      Call Read_aarg('tab',AF_tab)
      Call Read_iarg('ic',ic)
      Call Read_iarg('jc',jc)
      Call Read_iarg('jort',jort)

      Open(nub,file=AF_bnk,form='UNFORMATTED') 
      read(nub) ktype,kpol

      Open(out,file=AF_tab) 
      Open(in1,file=AF1) 
      Open(in2,file=AF2) 

! ... define configurations:

      Call R_conf(in1,in2,nub)

! ... orthogonality conditions:

      if(JORT.eq.-1)  write(out,'(a/)') &
       'Orthogonality mode:  full orthogonality (jort = -1)'
      if(JORT.eq. 0)  write(out,'(a/)') &
       'Orthogonality mode:  full non-orthogonality (jort = 0)'
      if(JORT.eq. 1)  write(out,'(a/)') &
       'Orthogonality mode:  partial orthogonality (jort = 1)'

      Call Pre_iort(0,0)

! ... Cycle over configurations:

      Do is = 1,ncfg1
       if(ic.gt.0.and.ic.ne.is) Cycle
       it = IC_term(is)
       Call Get_cfg_LS(is)
       write(out,'(71(''-''))')
       Call Prj_conf_LS (out,0.d0)
       Call Save_cfg_LS(1)
       ips = IP_state(is)
      Do js = ncfg1+1,ncfg2 
       if(jc.gt.0.and.jc.ne.js) Cycle
       jt = IC_term(js)
       Call Get_cfg_LS(js)
       Call Prj_conf_LS (out,0.d0)
       write(out,'(71(''-''))')
       write(out,*)
       Call Save_cfg_LS(2)
       jps = IP_state(js)

! ... read the bank:

      rewind(nub)
      read(nub) ktype,kpol
      Call Read_symc_LS(nub)
      Call Read_symt_LS(nub)
      Call Read_done_LS(nub)
      Call Read_det(nub)
      Call Read_def(nub)

      Call Alloc_ndet(-1)
      Call Alloc_ndef(-1)

   10 read(nub,end=20) C,ijt,int,idf

      ik=ijt/ibc; jk = mod(ijt,ibc)

      m = 0
      if(it.eq.ik.and.jt.eq.jk) m=1
      if(it.eq.jk.and.jt.eq.ik) m=2
      if(m.eq.0) go to 10

      ipc=ips; jpc=jps; if(m.eq.2) then; ipc=jps; jpc=ips; end if

! ... define integral:

      Call Decode_mult (itype,i1,i2,int)
      j1=IP_orb(ipc+i1)
      j2=IP_orb(jpc+i2)

! ... define the permutation of indexes:

      kz = 0
      if(m.eq.2.and.itype.ne.0) then
       kz = (Ltotal1-Ltotal2+Stotal1-Stotal2)/2
      end if

! ... find determinant overlaps for specific orbitals

      if(idf.gt.0) then

      kd=KPF(idf); ip=IPF(idf); NP(1:kd)=NPF(ip+1:ip+kd); md=0

      Do ii=1,kd
       id=NP(ii)/ibf; iext=mod(NP(ii),ibf); nd=KPD(id); ip=IPD(id)
       Do i=1,nd
        k=NPD(i+ip)
        NP1(i)=IP_ORB(ipc+k/ibd); NP2(i)=IP_orb(jpc+mod(k,ibd))
       End do

       mm = JDET_SIMP(kz,nd,NP1,NP2)
       if(mm.eq.1) Cycle; if(mm.eq.0) go to 10

       MP(1:nd) = NP1(1:nd)*ibd+NP2(1:nd) 
       jd = Iadd_ndet(nd,MP)
       md = md + 1; MP1(md) = jd; MP2(md) = iext

      End do 
    
      idf = 0
      if(md.gt.0) then 
       MP(1:md) = MP1(1:md)*ibf + MP2(1:md)
       idf = Iadd_ndef(md,MP)
      end if

      end if   !  idf > 0

      C = C * (-1)**kz

      Call PRI_MULT(out,C,itype,kpol,j1,j2,idf)

      write(out,*)

      go to 10
   20 Continue

      End do;  End do  ! is,js


      End ! program MULT_TAB
 

