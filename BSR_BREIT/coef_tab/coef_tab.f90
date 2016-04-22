!=====================================================================
!     UTIL      C O E F _  T A B
!
!               C O P Y R I G H T -- 2003
!
!     Written by:   Oleg Zatsarinny,   email: oleg_zoi@yahoo.com
!                         
!======================================================================
!
!     print (to file COEF.TAB) the integral coefficients for two 
!     selected states in CFG.INP according to INT_BNK file
!
!----------------------------------------------------------------------

      Use configs; Use param_br
      Use det_list; Use def_list; Use zoef_list

      Implicit none

      Integer(4), parameter :: mrecl = 100000

      Integer(4), Parameter :: ncase = 11
      Integer(4), Dimension(ncase) :: kcoef
      Integer(4) :: ncoef

      Integer(4), parameter :: noper = 7
      Integer(1), Allocatable, Dimension(:,:) :: IT_oper
      Integer(1), Dimension(noper) :: joper

      Integer(4) :: nui,nub,nuc,nur
      Integer(4) :: mo,mi,mso,mee,msoo,mss,moo
      Integer(4) :: icase,idf,kpol,int,ijt,kd,nd,kz,kzz,ip,id,jd,iext
      Integer(4) :: ic,ic1,ic2, it,it1,it2, jt,jt1,jt2, is,is1,is2
      Integer(4) :: I1,I2,I3,I4, J1,J2,J3,J4, IL1,IL2
      Integer(4) :: i,j,ii,ij,k,m,mm,mt,md, ldet,ldef
      Integer(4) :: jc1,jc2, iarg, jdiag
      Integer(4), external :: IARGC

      Integer(4), parameter :: mme = 100
      Integer(4), Dimension(mme) :: IPN1,IPN2, NP,NP1,NP2, MP,MP1,MP2

      Integer(4), EXternal :: IDET_SIMP, Iadd_int, Nadd_det, Nadd_def 

      Real(8) :: C
 
      Character(26) ::  AI  
      Character(40) ::  BI
      Character(80) ::  name = ' '
      Character(80) ::  AF_c   = 'cfg.inp'
      Character(80) ::  AF_bnk = 'int_bnk'
      Character(80) ::  AF_tab = 'coef.tab'
        
      iarg = IARGC()
      if(iarg.lt.1) then
       write(*,*)
       write(*,*) 'coef_tab prints to file coef.tab the integral coefficients'
       write(*,*) 'for two selected states in cfg.inp according to int_bnk file'
       write(*,*)
       write(*,*) 'arguments:'
       write(*,*)
       write(*,*) 'name   -  first argument without = , then it is supposed'
       write(*,*) '          name.c, name.bnk -> name.tab'
       write(*,*)
       write(*,*) 'c=...    -  c-file with configurations  [cfg.inp]'
       write(*,*) 'bnk=...  -  bnk-file with coefficients [int_bnk]'
       write(*,*) 'tab=...  -  output tables [coef.tab]'
       write(*,*) 'jort={-1,0,1}  - orbital orthogonality mode:'
       write(*,*) 'jort=-1  -  full orthogonality'  
       write(*,*) 'jort= 0  -  full non-orthogonality'  
       write(*,*) 'jort= 1  -  partial orthogonality, default'  
       write(*,*) 'jc1=.. jc2=.. -  matrix-element index; [0,0] - all'  
       Stop ' '
      end if

      Call Read_name(name)
      if(LEN_TRIM(name).ne.0) then
        AF_c = trim(name)//'.c'  
        AF_bnk = trim(name)//'.bnk'  
        AF_tab = trim(name)//'.tab'  
      end if
      Call Read_aarg('c',AF_c)
      Call Read_aarg('bnk',AF_bnk)
      Call Read_aarg('tab',AF_tab)
      jc1=0;   Call Read_iarg('jc1',jc1)
      jc2=0;   Call Read_iarg('jc2',jc2)
      jort=-1; Call Read_iarg('jort',jort)
      jdiag=1; Call Read_iarg('jdiag',jdiag)

!----------------------------------------------------------------------
!                                                                files:
! ... INT.BNK file:

      nub=1; Open(nub,file=AF_bnk,form='UNFORMATTED',status='OLD')

! ... cfg.inp file:
      
      nuc=2; Open(nuc,file=AF_c,status='OLD')

! ....output coef.tab file:

      nur=3; Open(nur,file=AF_tab)

!----------------------------------------------------------------------
!                                                define configurations:
      Call RR_conf(nuc,nub)

      write(nur,'(/a/)') 'CFG.INP contains:'
      write(nur,'(a,i8)')  'ncfg   = ',ncfg
      write(nur,'(a,i8)')  'nsymc  = ',nsymc
      write(nur,'(a,i8)')  'nsymt  = ',nsymt      
      write(nur,'(/a,i8/)')'nwf    = ',nwf
      write(nur,'(10a6)') (ELF(i),i=1,nwf)

!----------------------------------------------------------------------
!                        information about done calculation in int.bnk:
      read(nub) mt           
      Allocate(IT_oper(noper,mt))
      Call R_i1(nub,mrecl,noper,mt,mt,IT_oper) ! read(nub) (IT_oper(:,i),i=1,ii)

!----------------------------------------------------------------------
!                                                         determinants:
      read(nub) ndet,ldet; jdet=ldet/ndet+1
      Call Alloc_det(ndet)
      Call R_i4(nub,mrecl,ndet,KPD)  ! read(nub) KPD(1:ndet)
      Call R_i4(nub,mrecl,ndet,IPD)  ! read(nub) IPD(1:ndet)
      Call R_i4(nub,mrecl,ldet,NPD)  ! read(nub) NPD(1:ldet)

      read(nub) ndef,ldef; ; jdef=ldef/ndef+1
      Call Alloc_def(ndef)
      Call R_i4(nub,mrecl,ndef,KPF)  ! read(nub) KPF(1:ndef)
      Call R_i4(nub,mrecl,ndef,IPF)  ! read(nub) IPF(1:ndef)
      Call R_i4(nub,mrecl,ldef,NPF)  ! read(nub) NPF(1:ldef)
!----------------------------------------------------------------------
!                                                       read the bank:
      ncoef=0; kcoef=0
    1 read(nub,end=2) C,ijt,int,idf
      ncoef = ncoef + 1
      Call Decode_int(icase,kpol,I1,I2,I3,I4,int)
      kcoef(icase) = kcoef(icase) + 1
      go to 1
    2 Continue

      write(nur,'(/a/)') 'INT.BNK contains: '
      write(nur,'(a,2i8)')  'ndet     = ',ndet,kdet
      write(nur,'(a,2i8)')  'ndef     = ',ndef ,kdef
      write(nur,'(/a,i8/)') 'ncoef    = ',ncoef
      Do i = 1,ncase
       if(kcoef(i).gt.0) &
       write(nur,'(a,i2,a,i8)') 'icase(',i,') = ', kcoef(i)
      End do

!----------------------------------------------------------------------
!                                             orthogonality conditions:

      if(JORT.eq.-1)  write(nur,'(/a/)') &
       'Orthogonality mode:  full orthogonality (jort = -1)'
      if(JORT.eq. 0)  write(nur,'(/a/)') &
       'Orthogonality mode:  full non-orthogonality (jort = 0)'
      if(JORT.eq. 1)  write(nur,'(/a/)') &
       'Orthogonality mode: partial orthogonality (jort = 1)'

      Call Pre_iort(nuc)

!----------------------------------------------------------------------
!                                      input required matrix elements:

      Do ic1=1,ncfg; Do ic2=1,ncfg

       if(jc1.gt.0.and.jc1.ne.ic1) Cycle 
       if(jc2.gt.0.and.jc2.ne.ic2) Cycle 
       if(jdiag.eq.0.and.ic1.ne.ic2) Cycle

! ... find corr. it,jt

      Do it = 1,nsymt
       is1=mod(IT_stat(it),ibc); is2=IT_stat(it)/ibc
       Do is = is1,is2
        ic = IP_stat(is)
        if(ic.eq.ic1) it1 = it
        if(ic.eq.ic2) it2 = it
       End do
      End do

! ... find pointer for orbitals:
  
      no1 = NOCCSH(ic1); IPN1(1:no1) = NOCORB(1:no1,ic1)
      no2 = NOCCSH(ic2); IPN2(1:no2) = NOCORB(1:no2,ic2)

! ... find the indication on done calculation:

      i = max(it1,it2); j = min(it1,it2); ij = (i-1)*i/2 + j
      joper(:) = IT_oper(:,ij); Call Read_iarr('joper',7,joper)  

! ... skip the common data:

      rewind(nub)
      read(nub) nsymt, nsymc
      Call RR_a (nub,mrecl,nsymc,AI)     
      Call RR_a (nub,mrecl,nsymt,BI)     
      Call RR_i4(nub,mrecl,nsymt)    
      read(nub) ii
      Call RR_i1(nub,mrecl,7,ii)  
      read(nub) ndet,ldet
      Call RR_i4(nub,mrecl,ndet)    
      Call RR_i4(nub,mrecl,ndet)    
      Call RR_i4(nub,mrecl,ldet)    
      read(nub) ndef,ldef
      Call RR_i4(nub,mrecl,ndef)    
      Call RR_i4(nub,mrecl,ndef)    
      Call RR_i4(nub,mrecl,ldef)    

!----------------------------------------------------------------------
! ... read the specific data:

      nzoef=0; Call Alloc_ndet(0); Call Alloc_ndef(0)
   10 read(nub,end=20) C,ijt,int,idf
      it=ijt/ibc; jt = mod(ijt,ibc)
      m = 0
      if(it.eq.it1.and.jt.eq.it2) m=1
      if(it.eq.it2.and.jt.eq.it1) m=2
      if(m.eq.0) go to 10

! ... define integral:

      Call Decode_int(icase,kpol,I1,I2,I3,I4,int)

      if(m.eq.1) then
       j1 = IPN1(i1); j2 = IPN1(i2); j3 = IPN2(i3); j4 = IPN2(i4)
      else
       j1 = IPN2(i1); j2 = IPN2(i2); j3 = IPN1(i3); j4 = IPN1(i4)
      end if	 

!      Call Check_int(C,icase,kpol,j1,j2,j3,j4)

      int = Iadd_int(icase,kpol,j1,j2,j3,j4)

! ... define the permutation of indexes:

      kzz = 0
      if(m.eq.2.and.icase.ge.7) then
       IL1=LPL(no1,ic1); IS1=LPS(no1,ic1)
       IL2=LPL(no2,ic2); IS2=LPS(no2,ic2)
       kzz = (IL1-IL2+IS1-IS2)/2
      end if

!----------------------------------------------------------------------
! ...  find determinant overlaps for specific orbitals

      kz = 0
      if(idf.gt.0) then

      kd=KPF(idf); ip=IPF(idf); NP(1:kd)=NPF(ip+1:ip+kd); md=0
      Do ii=1,kd
       id=NP(ii)/ibf; iext=mod(NP(ii),ibf); nd=KPD(id); ip=IPD(id)
       Do i=1,nd
        k=NPD(i+ip)
        if(m.eq.1) then
         NP1(i)=IPN1(k/ibd); NP2(i)=IPN2(mod(k,ibd))
        else
         NP1(i)=IPN2(k/ibd); NP2(i)=IPN1(mod(k,ibd))
        end if
       End do

       mm = IDET_SIMP(kz,nd,NP1,NP2)

       if(mm.eq.1) Cycle; if(mm.eq.0) go to 10
       MP(1:nd) = NP1(1:nd)*ibd+NP2(1:nd) 
       jd = Nadd_det(nd,MP)
       md = md + 1; MP1(md) = jd; MP2(md) = iext
      End do 
    
       idf = 0
       if(md.gt.0) then 
        MP(1:md) = MP1(1:md)*ibf + MP2(1:md)
        idf = Nadd_def(md,MP)
       end if
      end if   !  idf > 0

      C = C * (-1)**(kz+kzz);  Call Iadd_zoef(C,int,idf)

      go to 10
   20 Continue
!----------------------------------------------------------------------
!
      mo = joper(1); mi = joper(2); mee = joper(3); mso = joper(4);
      msoo = joper(5); mss = joper(6); moo = joper(7)
      Call PRI_COEF(nur,ic1,ic2,mo,mi,mso,mee,msoo,mss,moo)

      End do; End do  ! ic1, ic2

      End ! program COEF_TAB
 

!======================================================================
      Subroutine Check_int(C,icase,k,j1,j2,j3,j4)
!======================================================================

      USE param_br

      Implicit none
      Integer(4) :: icase,k,j1,j2,j3,j4,j
      Real(8) :: C

      Select case(icase)

      Case(4,5)   ! Mk,Rk

       if(j1.gt.j2) then
        j=j1;j1=j2;j2=j; j=j3;j3=j4;j4=j
       elseif(j1.eq.j2.and.j3.gt.j4) then
        j=j3;j3=j4;j4=j
       end if

      Case(3)     ! Tk

       if(j1.gt.j2) then
        j=j1;j1=j2;j2=j; j=j3;j3=j4;j4=j
       elseif(j1.eq.j2.and.j3.gt.j4) then
        j=j3;j3=j4;j4=j
       end if

       if(j1.eq.j2.and.j2.eq.j3.and.j3.eq.j4) then
        C=-C*(k-1)*(k+2)/(k+k+1)/2 ! /2 ???
        icase=4
       end if

!----------------------------------------------------------------------
!                                                        T[k] --> M[k]:
!     Here are used the relations:
!
!     U[k](a,b;c,d) + U[k](c,d;a,b) =    A(a,b;c,d) -
!
!     - (k-1)(k+2)/(2k+1) { N[k-1](a,b;c,d) + N[k-1](b,a;d,c) }
!
!     and   U[k](a,b;c,d) = T[k](a,b;c,d) + T[k](c,b;a,d).
!
!     then a=b=c=d the integrals A(a,a,a,a) disapear because
!     the T and U integrals enter only as I[k+1] - I[k-1].
!
!     In particular:  Uk(a,a,a,a)=2Tk(a,a,a,a)
!                                =-(k-1)(k+2)/(2k+1) (1/2) M[k-1]
!---------------------------------------------------------------------

      Case(9)     ! Vk(a,b;a,b)

      if(j1.eq.j3.and.j2.eq.j4) then
       if(is_soo.eq.1) then  ! ... in case of Vk -> V'k:
        C = 0.d0
       else
!      we cannot do nothing because two integrals!
       end if
      end if
!----------------------------------------------------------------------
!                                                        V[k] --> N[k]:
!
!     elseif(met.eq.9.and.j1.eq.j3.and.j2.eq.j4) then
!
!      S=-C*(k+1)/2
!      int=0
!      Call INT_pack(8,k-1,j1,j2,j3,j4,int)
!      i=Iadd_boef(S,int)
!      S= C*(k+0)/2
!      int=0
!      Call INT_pack(8,k+1,j1,j2,j3,j4,int)
!      i=Iadd_boef(S,int)
!
!      here is used the relation:   (R.Glass, Z.Physik,A292,131(1979))
!
!      2*V[k](r,s;r,s) = -(k+2)*N[k-1](s,r;s,r) + (k+1)N[k+1](r,s;r,s)
!
!      that is consequence of more general relation
!
!      V[k](r,s;r's') + V[k](r's';r,s) =
!                 -(k+2)*N[k-1](s,r;s',r') + (k+1)N[k+1](r,s;r',s')
!
!----------------------------------------------------------------------

       End Select

       End Subroutine Check_int



 
!====================================================================
      SUBROUTINE PRI_COEF (nu,ic1,ic2,mo,mi,mso,mee,msoo,mss,moo)
!====================================================================
!
!     Prints the resulting angular integrals
!
!--------------------------------------------------------------------

      USE configs;  USE zoef_list; Use int_list

      if(nzoef.eq.0) Return

! ...  print HEADER:

      write(nu,'(/70(''=''))')
      write(nu,'(12x,''< state'',i2,'' || (H-E) || state'',i2,''>'')') &
            ic1,ic2
      write(nu,'(70(''=''))') 
      Call Pri_ic (nu,ic1,0.d0);  write(nu,'(70(''-''))')
      Call Pri_ic (nu,ic2,0.d0);  write(nu,'(70(''-''))')

      if(nzoef.eq.0) Return

!----------------------------------------------------------------------
!                                                             Overlaps:
      if(mo.eq.0) then 
       write(nu,'(/a/)') 'Overlaps:  no calculations!'
      else
       write(nu,'(/12x,'' Overlaps''/)')
       Do i=1,nzoef; Call Pri_coef1(nu,i,11); End do
      end if
!---------------------------------------------------------------------
!                                                 Coulomb interaction:
      if(mi+mee.eq.0) then 
       write(nu,'(/a/)') 'Coulomb interactions:  no calculations!'
      else
       write(nu,'(/12x,'' Coulomb interaction''/)')
       Do i=1,nzoef; Call Pri_coef1(nu,i,5); End do
       Do i=1,nzoef; Call Pri_coef1(nu,i,6); End do
      end if
!----------------------------------------------------------------------
!                                              Orbit-orbit interaction:
      if(moo.eq.0) then 
       write(nu,'(/a/)') ' Orbit-orbit:          no calculations!'
      else
       write(nu,'(/12x,''  Orbit-orbit interaction ''/)')
       Do i=1,nzoef; Call Pri_coef1(nu,i,3); End do
       Do i=1,nzoef; Call Pri_coef1(nu,i,4); End do
      end if
!----------------------------------------------------------------------
!                                               Spin-orbit interaction:
      if(mso.eq.0) then 
       write(nu,'(/a/)') ' Spin-orbit:           no calculations!'
      else
       write(nu,'(/12x,''  Spin-orbit interaction ''/)')
       Do i=1,nzoef; Call Pri_coef1(nu,i,7); End do
      end if
!----------------------------------------------------------------------
!                                         Spin-other-orbit interaction:
      if(msoo.eq.0) then 
       write(nu,'(/a/)') ' Spin-other-orbit:     no calculations!'
      else
       write(nu,'(/12x,''  Spin-other-orbit interaction ''/)')
       Do i=1,nzoef; Call Pri_coef1(nu,i,8); End do
       Do i=1,nzoef; Call Pri_coef1(nu,i,9); End do
      end if
!---------------------------------------------------------------------
!                                               Spin-spin interaction:
      if(mso.eq.0) then 
       write(nu,'(/a/)') ' Spin-spin:            no calculations!'
      else
       write(nu,'(/12x,''  Spin-spin interaction ''/)')
       Do i=1,nzoef; Call Pri_coef1(nu,i,10); End do
      end if
!----------------------------------------------------------------------
      write(nu,'(70(''-'')/)')

      End SUBROUTINE PRI_COEF
 

!======================================================================
      SUBROUTINE PRI_COEF1 (nu,ii,met)
!======================================================================
!
!     Prints one angular integral
!
!----------------------------------------------------------------------

      USE configs; Use param_br;
      Use zoef_list; Use int_list; USE ndet_list; Use ndef_list

      Implicit none

      Character(1), Dimension(10) :: AI
      Data AI /' ',' ','T','M','R','L','Z','N','V','N'/
      Character(3) :: EL,EL1,EL2,EL3,EL4
      Character(220) :: A
      Character(1) :: a1,b1,a2
      Integer(4), Intent(in) :: nu,ii,met
      Integer(4) :: i,j, i1,i2,i3,i4, k,k1,k2, m1,m2
      Integer(4) :: ia,ip,np,idf,int,kd,nd,ns,ni


      if(abs(Zoef(ii)).lt.Eps_C) Return
!----------------------------------------------------------------------

      i=IZ_int(ii); int=Jcase(i); if(int.ne.met) Return
      k = Jpol(i)
      if(int.eq.4.or.int.eq.8.or.int.eq.9.or.int.eq.10) k=k-1
      i1 = J1(i); i2 = J2(i); i3 = J3(i); i4 = J4(i)

      write(A,'(220a1)') (' ',i=1,220)

      Select case(int)

      Case(11)                            ! Overlap
      ia=0

      Case(6,7)                           ! One-electron integrals
      EL1 = ELF(I1)(2:4); EL2 = ELF(I3)(2:4)
      Write(A(1:15),'(a1,a4,2(a3,a1))') AI(int),'   (',EL1,',',EL2,')'
      ia=15

      Case(3,4,5,8,9,10)                   ! Two-electron integrals
      EL1 = ELF(I1)(2:4); EL2 = ELF(I2)(2:4);
      EL3 = ELF(I3)(2:4); EL4 = ELF(I4)(2:4)
      Write(A(1:23),'(a1,i2,a2,4(a3,a1))') &
          AI(int),k,' (',EL1,',',EL2,';',EL3,',',EL4,')'
      ia=23
 
      Case default
        write(*,*) ' int=',int
        Stop ' Pri_coef: unknown case'
      End select

! ... determinant factor:

      idf=IZ_df(ii)
      if(idf.gt.0) then
       kd=KPF(idf);  ip=IPF(idf)
       Do i=1,kd
        ns=mod(NPF(ip+i),ibf); ni=NPF(ip+i)/ibf; 
        nd=KPD(ni); np=IPD(ni)
        ia=ia+1; A(ia:ia)='<'
       Do j=1,nd
        ni=NPD(np+j)/ibd;  EL = ELF(ni)(2:4)
        Write(A(ia+1:ia+5),'(10(a3,'',''))') EL
        ia=ia+4
       End do
       A(ia:ia)='|'
       Do j=1,nd
        ni=mod(NPD(np+j),ibd); EL = ELF(ni)(2:4)
        Write(A(ia+1:ia+5),'(10(a3,'',''))') EL
        ia=ia+4
       End do
       if(ns.gt.1) then
        A(ia:ia+1)='>^'; ia=ia+2; write(A(ia:ia),'(i1)') ns; ia=ia+1
       else
        A(ia:ia+1)='> '; ia=ia+2
       end if
      End do       ! Over kd
      end if

      Call Num(Zoef(ii),k1,k2,999999,1.d-9)

      a1='['; b1=':';a2=']'
      k1 = iabs(k1); m1=sqrt(1.*k1)+0.1;m2=sqrt(1.*k2)+0.1
      if(k1.eq.m1*m1.and.k2.eq.m2*m2) then
      k1=m1;k2=m2;a1='(';a2=')';end if
      Write(nu,'(f12.6,1x,a1,i6,a1,i6,a1,2x,220a1)') &
                Zoef(ii),A1,iabs(k1),B1,k2,A2,(a(i:i),i=1,ia)


      End SUBROUTINE PRI_COEF1

!====================================================================
      Subroutine NUM(z,k1,k2,i_max,acr)
!====================================================================
!
!    finds  k1 and k2 such that abs(z) = sqrt ( k1 / k2 )
!    and sign(z) = sign(k1)
!
!    k1, k2 < i_max
!    acr - tollerance for fitting
!--------------------------------------------------------------------

      Implicit none

      Real(8), intent(in) :: z, acr
      Integer(4), intent(in) :: i_max
      Integer(4), intent(out) :: k1,k2

      Integer(4) :: i_min, i1,i2, m1,m2 
      Real(8) :: zz, s,ss, s1,s2

      k1=0
      k2=1
      zz=z*z
      if(zz.lt.1.d0/i_max) Return

      s=zz
      i_min=zz
      if(i_min.lt.1) i_min=1
      Do i1=i_min,i_max
       s1=DBLE(i1)
       m1=s1/zz
       if(m1.lt.1) m1=1
       m2=m1+1
      Do i2=m1,m2
       s2=DBLE(i2)
       ss=abs(s1/s2-zz)
       if(ss.lt.s) then
        s = ss
        k1 = i1
        if(z.lt.0.d0) k1 = -k1
        k2 = i2
        end if
       if(ss.lt.acr) return
      End do
      End do

      End Subroutine NUM


