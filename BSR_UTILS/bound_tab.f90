!--------------------------------------------------------------------
!     bound.nnn  -->  bound_tab - list of states
!--------------------------------------------------------------------

      Implicit real(8) (A-H,O-Z)
      
      Integer, allocatable :: jst(:), ilsp(:), IPT(:)
      Integer, allocatable :: ISS(:), ILS(:), IPS(:)
      
      Real(8), allocatable :: E(:)
      
      Character(64) :: label1
      Character(64), Allocatable :: Label(:)
      
      Integer :: in=2;   Character(40) :: AF
      Character(40), Allocatable :: BF(:)

      Character(3) :: ALSP
      Character(200) :: AS,BS

      Call inf_bound_tab
!--------------------------------------------------------------------
!                                                 input/output files:      
      nu = 11; Open(nu,file='bound_tab')

      Z  =0.d0;  Call Read_rpar(nu,'Z',Z)
      AWT=0.d0;  Call Read_rpar(nu,'AWT',AWT)
      E0 =0.d0;  Call Read_rpar(nu,'E0',E0)
      EM =0.d0;  Call Read_rpar(nu,'EM',EM)
      klsp1=1;   Call Read_ipar(nu,'klsp1',klsp1)
      klsp2=1;   Call Read_ipar(nu,'klsp2',klsp2)
      klsp3=1;   Call Read_ipar(nu,'klsp3',klsp3)
      mstate=0;  Call Read_ipar(nu,'mstate',mstate)

      Call Read_rarg('Z',Z)
      Call Read_rarg('AWT',AWT)
      Call Read_rarg('E0',E0)
      Call Read_rarg('EM',EM)
      Call Read_iarg('klsp1',klsp1)
      Call Read_iarg('klsp2',klsp2)
      Call Read_iarg('klsp3',klsp3)
      Call Read_iarg('mstate',mstate)

      Call Conv_au (Z,AWT,au_cm,au_eV,0)      

      klsp=0;  Call Read_iarg('klsp',klsp)
      if(klsp.ne.0) then; klsp1=klsp; klsp2=klsp; end if
      if(klsp2.lt.klsp1) klsp2=klsp1

!--------------------------------------------------------------------
!                                                         nstate - ?:
      nstate=0
      nterm=0
      Do klsp=klsp1,klsp2,klsp3
       write(AF,'(a6,i3.3)') 'bound.',klsp
       if(Icheck_file(AF).eq.0) Cycle
       Open(in,file=AF,status='OLD'); rewind(in)
       nterm=nterm+1
       read(in,*) ns,nch,ncp,khm,nbound 
       if(mstate.gt.0.and.nbound.gt.mstate) nbound=mstate
       nstate=nstate+nbound
      End do  

      if(nstate.eq.0) Stop 'nstate = 0'

      Allocate(jst(nstate), E(nstate), Label(nstate),IPT(nstate))
      Allocate(BF(nterm),ISS(nterm),ILS(nterm),IPS(nterm),ilsp(nterm))

!--------------------------------------------------------------------
!                                                      read energies:              
      ii=0; is=0
      Do klsp=klsp1,klsp2,klsp3
       write(AF,'(a6,i3.3)') 'bound.',klsp
       if(Icheck_file(AF).eq.0) Cycle
       Open(in,file=AF,status='OLD'); rewind(in)
       is=is+1;  BF(is)=AF
       read(in,*) ns,nch,ncp,khm,nbound,ISS(is),ILS(is),IPS(is)   
       if(mstate.gt.0.and.nbound.gt.mstate) nbound=mstate
       ilsp(is) = klsp
       Do i = 1,nbound
        ii = ii + 1
        read(in,'(i5,2x,a)') js,Label(ii)
        read(in,*) E(ii)
        read(in,'(5D15.8)') (S,j=1,khm)
        jst(ii) = js+10000*is
       End do
      End do
      
      Call SortR(nstate,E,ipt)

!--------------------------------------------------------------------
!                                                excitation energies:
      imax=0
      Do i = 1,nstate
       ii = LEN_TRIM(Label(i)); if(ii.gt.imax) imax=ii
      End do

      if(E0.eq.0.d0) E0 = E(IPT(1))
      
      rewind(nu) ; isf=0

      BS = ' '
      jmax=imax-6; if(jmax.lt.0) jmax=1
      write(AS,'(a,a)') '  klsp  sol   label', BS(1:jmax)
      jmax=imax+15
      write(AS(jmax:),'(a)') &
        ' L  S  P       E_Ry         E_eV            E_cm         E_au'
      write(nu,'(a/)') trim(AS)

      Do j = 1,nstate;  i=IPT(j)

       ee = E(i) - E0;     if(EM.ne.0.d0.and.E(i).gt.EM) Cycle
       e_eV = ee * au_eV
       e_cm = ee * au_cm
       
       label1 = Label(i)

       js = mod(jst(i),10000); is=jst(i)/10000
       
       write(nu,'(2i5,3x,a,3i3,f14.6,F14.5,f16.1,F16.8)') & 
        ilsp(is),js,label1(1:imax),ISS(is),ILS(is),IPS(is), &
        ee,e_eV,e_cm,E(i) 
       isf = isf + 1 
      
      End do 

      write(nu,'(a)') '*'

      write(nu,'(a,f5.0)') 'Z = ',Z
      write(nu,'(a,f7.3)') 'AWT = ',AWT
      write(nu,'(a,f18.8)') 'E0 = ',E0
      write(nu,'(a,f18.8)') 'EM = ',EM
      write(nu,'(a,i3)') 'klsp1 = ',klsp1
      write(nu,'(a,i3)') 'klsp2 = ',klsp2
      write(nu,'(a,i3)') 'klsp3 = ',klsp3

      write(nu,*)
      write(nu,'(a,f16.6)') 'au_eV= ',au_eV
      write(nu,'(a,f12.2)') 'au_cm= ',au_cm

      write(nu,*)
      write(nu,'(a,i5)') 'nstate = ',isf
      write(nu,'(a,i5)') 'mstate = ',mstate
      
      End  !  program bound_tab


!======================================================================
      Subroutine inf_bound_tab
!======================================================================
!     provide screen information about bound_tab utility
!----------------------------------------------------------------------
       
      Implicit real(8) (a-h,o-z)

      Character(80) :: A

      iarg = IARGC()
      if(iarg.eq.0) Return
      Call GETARG(1,A)      
      if(A.ne.'?') Return

      write(*,'(a)') &
'                                                                                  ',&
'     BOUND_TAB creats file "bound_tab" with list of solutions from bound.nnn files ',&
'                                                                                  ',&
'     Arguments: (can be also defined in the end of the existing bound_tab file)   ',&
'                (all arguments are optional and have default values)              ',&
'                                                                                  ',&
'     Z   [0.d0]  - nuclear charge, required for a.u. -> eV                        ',&
'     AWT [0.d0]  - nuclear weight                                                 ',&
'     E0  [0.d0]  - reference point for energies in a.u.                           ',&
'     EM  [0.d0]  - maximum energy for output solutions in a.u.                    ',&
'     klsp1  [1]  - initial partial wave                                           ',&
'     klsp2  [1]  - final partial wave                                             ',&
'     klsp3  [1]  - increment (klsp=klsp1,klsp2,klsp3                              ',&
'     mstate [0]  - maximum number of states for one partial wave                  ',&
'                   (0 means all)                                                  ',&
'                                                                                  ',&
'     Examples:  bound_tab  Z=14  mstate=20  klsp1=1  klsp2=12  E0=...             ',&
'                bound_tab  klsp=5                                                 ',&
'                bound_tab  Z=14  mstate=20  klsp1=1  klsp2=12  E0=...             ',&
'                                                                                  '
      Stop                                                                     
                                                                                   
      End Subroutine inf_bound_tab                                                 

