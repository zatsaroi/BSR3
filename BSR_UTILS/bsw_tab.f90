!======================================================================
!     UTILITY       B S W _ T A B
!
!               C O P Y R I G H T -- 2003
!
!     Written by:   Oleg Zatsarinny
!
!======================================================================
!     connverts the B-spline radial orbital wbs-files into tab-files
!     suitable for grafic display
!----------------------------------------------------------------------
!     INPUT FILES:   name.bsw
!     OUTPUT FILES:  name.tab
!     Call as:       bsw_tab name.bsw
!----------------------------------------------------------------------
      Use spline_param; Use spline_atomic; Use spline_grid
      Use spline_galerkin

      Implicit real(8) (A-H,O-Z)
      Character(1) :: ans
      Character(4) :: EL4
      Character(4), allocatable :: EL(:)
      Character(40) :: AF,BF
      Real(8), allocatable :: v(:), w(:,:), R(:)
      
      iarg = COMMAND_ARGUMENT_COUNT() 
      if(iarg.gt.0) Call GET_COMMAND_ARGUMENT(1,AF)

      if(iarg.le.0.or.AF.eq.'?') then
        write(*,'(a)') 'bsw_tab converts name.bsw  to name.tab'
        write(*,'(a)') 'bsw - unformatted  BSR file for radial functions'
        write(*,'(a)') 'tab - column representation of orbitals, ready for plotting'
        write(*,'(a)') 
        write(*,'(a)') 'Call as:  bsw_tab  name.bsw '
        write(*,'(a)') 
        write(*,'(a)') 'orbitals for output will be asking'
        Stop ' '
      end if        

!----------------------------------------------------------------------
! ... input data: 
         
      inp=1; Open(inp,file=AF,status='OLD',form='UNFORMATTED')

      i=index(AF,'.') 
      iout=2; BF=AF(1:i)//'tab';  Open(iout,file=BF)

!----------------------------------------------------------------------
! ... sets up grid points and initializes the values of the spline: 
    
      CALL define_grid(zz);  CALL define_spline

!----------------------------------------------------------------------
! ... define number of orbitals:

      Allocate(v(ns))

      nn= 0
    1 read(inp,end=2) el4,zw,hw,hmw,rmw,ksw,nsw,m;  read(inp) v(1:m)
      nn=nn+1
      go to 1
    2 write(*,*) 'number of orbitals = ',nn

!----------------------------------------------------------------------
! ... sets up grid points and initializes the values of the spline: 

      NR = nv*ks+2; Allocate(R(NR),w(NR,nn),EL(nn))
      ii=1; R(1)=0.d0
      Do i=1,nv; Do j=1,ks; ii=ii+1; R(ii) = gr(i,j); End do; End do
      ii=ii+1; R(ii) = t(ns+1)

!----------------------------------------------------------------------
!                                                          
      rewind(inp)
      mm = 0
      Do i=1,nn
       read(inp) el4,zw,hw,hmw,rmw,ksw,nsw,m
       read(inp) v(1:m); if(m.lt.ns) v(m+1:ns)=0.d0
       write(*,*) EL4, '   Do you need it? (y|n): '
       read(*,'(a)') ANS
       if(ans.eq.'n') Cycle
       mm=mm+1
       Do j=1,nr; w(j,mm) = bvalu2(t,v,ns,ks,r(j),0); End do
       EL4(1:1)='a'; if(EL4(2:2).eq.' ') EL4(2:2)='a'; EL(mm) = EL4
      End do

      write(iout,'(20A13)') 'R   ',EL(1:mm)
      Do i=1,nr
       write(iout,'(20D13.5)') R(i),(w(i,j),j=1,mm)
      End do
      
      END   !  program bsr_tab
          
   