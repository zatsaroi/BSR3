!----------------------------------------------------------------------
!     target, thresholds --> target_states
!     provide information about target states and exp. thresholds
!     put E_exp in the target states
!----------------------------------------------------------------------

      Implicit real(8) (A-H,O-Z)
      Character(64), allocatable :: conf(:), name(:)
      Real(8), allocatable :: E_au(:), E_exp(:)
      Integer, allocatable :: J2(:), Ltarg(:), IStarg(:), IPtarg(:)
      Character(80) :: AF_inp='target'
      Character(80) :: AF_out='target_states'
      Character(200) :: AF, AF1,AF2
      Character(1) :: ap

      Call Read_aarg('inp',AF_inp)
      Call Check_file(AF_inp)
      inp=1; open(inp,file=AF_inp)

      ntarg=0; Call Read_ipar(inp,'ntarg',ntarg)
      if(ntarg.le.0) Stop 'ntarg <= 0'

      Allocate(E_au(ntarg), E_exp(ntarg), name(ntarg), conf(ntarg), J2(ntarg), &
               Ltarg(ntarg), IStarg(ntarg), IPtarg(ntarg))
      read(inp,*)
      Do i = 1, ntarg
       read(inp,*) name(i), conf(i), Ltarg(i),IStarg(i),IPtarg(i), E_au(i)
      End do
      if(IStarg(1).eq.0) J2 = Ltarg

      nu = 2
      Do i = 1, ntarg
       AF = trim(name(i))//'.c'
       open(nu, file = AF)
       rewind(nu)
       Call Idef_LabelC(nu,1,0,conf(i))
      End do

      Call Read_ipar(inp,'nz',nz);  Z=nz
      Call Conv_au (Z,AWT,au_cm,au_eV,0)      

      imax=0; jmax=0
      Do i = 1,ntarg
       ii = LEN_TRIM(name(i)); if(ii.gt.imax) imax=ii
       ii = LEN_TRIM(conf(i)); if(ii.gt.jmax) jmax=ii
      End do

      Call Read_aarg('out',AF_out)
      open(nu,file=AF_out)
      write(nu,'(a,i6)') 'ntarg =', ntarg
      write(nu,*)
      Do i = 1,ntarg
       E = E_au(i) - E_au(1)
       E_Ry = 2*E
       E_eV = E*au_eV
       E_cm = E*au_cm
       AF1 = name(i); AF2 = conf(i)
       write(nu,'(i4.4,3x,a,3x,a,2x, F20.8, F17.6, F15.5, F13.1)') &
          i, AF1(1:imax),AF2(1:jmax), E_au(i), E_Ry, E_ev, E_cm

      End do 

      End ! program

      