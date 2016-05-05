!=====================================================================
!     UTIL      C O E F _  T A B 3
!
!               C O P Y R I G H T -- 2011
!
!     Written by:   Oleg Zatsarinny,   email: oleg_zoi@yahoo.com
!                         
!======================================================================
!     prints from unformated INT_BNK file to COEF.TAB file the integral
!     coefficients for selected states in CFG.INP 
!----------------------------------------------------------------------
      Use param_br
      Use conf_LS
      Use orb_LS
      Use det_list
      Use def_list

      Implicit none

      Integer :: mo,mi,mso,mee,msoo,mss,moo
      Integer :: icase,idf,kpol,int,ijt,kd,nd,kz,kzz,ip,id,jd,iext
      Integer :: ic,ic1,ic2, it,it1,it2, jt,jt1,jt2, is,is1,is2
      Integer :: I1,I2,I3,I4, J1,J2,J3,J4, IL1,IL2
      Integer :: i,j,ii,ij,k,m,mm,mt,md
      Integer :: jc1,jc2, iarg, ip1,ip2

      Integer :: NNP(me),NNP1(me),NNP2(me), MP(me),MP1(me),MP2(me)

      Integer, external :: IDET_SIMP, Iadd_int, Iadd_ndet, Iadd_ndef, Def_ij, IARGC 

      Real(8) :: C
 
      Character(80) ::  name   = ' '
      Character(80) ::  AF_cfg = 'cfg.inp';  Integer :: nuc = 1 
      Character(80) ::  AF_bnk = 'int_bnk';  Integer :: nub = 2
      Character(80) ::  AF_tab = 'coef.tab'; Integer :: nur = 3
        
      iarg = IARGC()
      if(iarg.gt.0) Call GETARG(1,name)      

      if(iarg.lt.1.or.name.eq.'?') then
       write(*,*)
       write(*,*) 'coef_tab prints to file coef.tab the integral coefficients'
       write(*,*) 'for two selected states in c-file according to int_bnk file'
       write(*,*)
       write(*,*) 'Call as: coef_tab name [c= bnk= tab= jort= jc1= jc2= oper= ]'
       write(*,*)
       write(*,*) 'name   -  first position argument, then it is supposed we have'
       write(*,*) '          name.c, name.bnk and results will be in name.tab'
       write(*,*)
       write(*,*) 'c=...    -  c-file with configurations  [cfg.inp]'
       write(*,*) 'bnk=...  -  bnk-file with coefficients [int_bnk]'
       write(*,*) 'tab=...  -  output tables [coef.tab]'
       write(*,*) 'jort={-1,0,1}  - orbital orthogonality mode:'
       write(*,*) 'jort=-1  -  full orthogonality'  
       write(*,*) 'jort= 0  -  full non-orthogonality'  
       write(*,*) 'jort= 1  -  partial orthogonality, default'  
       write(*,*) 'jc1=.. jc2=.. -  matrix-element index; [0,0] - all'  
       write(*,*) 'oper=..  -   matrix-element index, as in the bsr_breit program'  
       Stop 
      end if

      if(LEN_TRIM(name).ne.0) then
        AF_cfg = trim(name)//'.c'  
        AF_bnk = trim(name)//'.bnk'  
        AF_tab = trim(name)//'.tab'  
      end if
      Call Read_aarg('c'  ,AF_cfg)
      Call Read_aarg('bnk',AF_bnk)
      Call Read_aarg('tab',AF_tab)
      jc1=0;   Call Read_iarg('jc1',jc1)
      jc2=0;   Call Read_iarg('jc2',jc2)
      if(jc1.lt.jc2) then; i=jc1; jc1=jc2; jc2=i; end if 
      jort=-1; Call Read_iarg('jort',jort)

      oper='1110000'                                  ! ???
      Call Read_aarg('oper'  ,oper  )
      read(oper,'(7i1)') ioper

!----------------------------------------------------------------------
!                                                                files:
! ... INT.BNK file:

      Call Check_file(AF_bnk)
      Open(nub,file=AF_bnk,form='UNFORMATTED')

! ... cfg.inp file:
      
      Call Check_file(AF_cfg)
      Open(nuc,file=AF_cfg)

! ....output coef.tab file:

      Open(nur,file=AF_tab)

!----------------------------------------------------------------------
!                                                define configurations:
      Call RR_conf(nuc,nub,nur)

      write(nur,'(/a,i8/)')'nwf    = ',nwf
      write(nur,'(10a6)') (ELF(i),i=1,nwf)

!----------------------------------------------------------------------
!                                                         determinants:
      Call read_det(nub)
      Call read_def(nub)
      write(nur,*) 
      write(nur,*) 'ndet,kdet =',ndet,kdet
      write(nur,*) 'ndef,kdef =',ndef,kdef

!----------------------------------------------------------------------
!                                             orthogonality conditions:

      if(JORT.eq.-1)  write(nur,'(/a/)') &
       'Orthogonality mode:  full orthogonality (jort = -1)'
      if(JORT.eq. 0)  write(nur,'(/a/)') &
       'Orthogonality mode:  full non-orthogonality (jort = 0)'
      if(JORT.eq. 1)  write(nur,'(/a/)') &
       'Orthogonality mode: partial orthogonality (jort = 1)'

      Call Pre_iort(nuc,0)

!----------------------------------------------------------------------
! ... read the specific data:

      Call Alloc_ndet(-1)
      Call Alloc_ndef(-1)

   10 read(nub,end=20) C,ijt,int,idf

      it=ijt/ibc; jt = mod(ijt,ibc)
      Call Decode_int(icase,kpol,I1,I2,I3,I4,int)

      Do ic1=1,ncfg; Do ic2=1,ic1
       if(jc1.gt.0.and.jc1.ne.ic1) Cycle 
       if(jc2.gt.0.and.jc2.ne.ic2) Cycle 
       it1 = IC_term(ic1)
       it2 = IC_term(ic2)
       m = 0
       if(it.eq.it1.and.jt.eq.it2) m=1
       if(it.eq.it2.and.jt.eq.it1) m=2
       if(m.eq.0) Cycle

       ip1 = ip_state(ic1)
       ip2 = ip_state(ic2)

      if(m.eq.1) then
       j1 = IP_orb(ip1+i1); j2 = IP_orb(ip1+i2)
       j3 = IP_orb(ip2+i3); j4 = IP_orb(ip2+i4)
      else
       j1 = IP_orb(ip2+i1); j2 = IP_orb(ip2+i2)
       j3 = IP_orb(ip1+i3); j4 = IP_orb(ip1+i4)
      end if	 

      int = Iadd_int(icase,kpol,j1,j2,j3,j4)

! ... define the permutation of indexes:

      kzz = 0
      if(m.eq.2.and.icase.ge.7) then
       Call Term_ic(ic1,IL1,IS1)
       Call Term_ic(ic2,IL2,IS2)
       kzz = (IL1-IL2+IS1-IS2)/2
      end if

!----------------------------------------------------------------------
! ... find determinant overlaps for specific orbitals

      kz = 0
      if(idf.gt.0) then

      kd=KPF(idf); ip=IPF(idf); NNP(1:kd)=NPF(ip+1:ip+kd); md=0

      Do ii=1,kd

       id=NNP(ii)/ibf; iext=mod(NNP(ii),ibf); nd=KPD(id); ip=IPD(id)

       Do i=1,nd
        k=NPD(i+ip)
        if(m.eq.1) then
         NNP1(i)=IP_orb(ip1+k/ibd); NNP2(i)=IP_orb(ip2+mod(k,ibd))
        else
         NNP1(i)=IP_orb(ip2+k/ibd); NNP2(i)=IP_orb(ip1+mod(k,ibd))
        end if
       End do

       mm = IDET_SIMP(kz,nd,NNP1,NNP2,mwf,IORT)

       if(mm.eq.1) Cycle; if(mm.eq.0) Exit
 
       MP(1:nd) = NNP1(1:nd)*ibd+NNP2(1:nd) 
       jd = Iadd_ndet(nd,MP)
       md = md + 1; MP1(md) = jd; MP2(md) = iext
 
      End do  ! over kd
      if(mm.eq.0) Cycle

       idf = 0
       if(md.gt.0) then 
        MP(1:md) = MP1(1:md)*ibf + MP2(1:md)
        idf = Iadd_ndef(md,MP)
       end if
      
      end if   !  idf > 0

!----------------------------------------------------------------------

      C = C * (-1)**(kz+kzz);  Call Add_rk4_data(int,idf,ic1,ic2,C)

      End do; End do ! ic1, ic2

      go to 10
   20 Continue

!----------------------------------------------------------------------
!                                                        print results:

      Do ic1=1,ncfg; Do ic2=1,ic1
       if(jc1.gt.0.and.jc1.ne.ic1) Cycle 
       if(jc2.gt.0.and.jc2.ne.ic2) Cycle 

       it1 = IC_term(ic1)
       it2 = IC_term(ic2)
       ij = Def_ij (it1,it2)
       joper = IT_oper(:,ij)
       Do i=1,noper; if(ioper(i).eq.0) joper(i)=-1; End do

       mo = joper(1); mi = joper(2); mee = joper(3); mso = joper(4);
       msoo = joper(5); mss = joper(6); moo = joper(7)

       Call PRI_COEF(nur,ic1,ic2,mo,mi,mso,mee,msoo,mss,moo)

      End do; End do  

      End ! program COEF_TAB
