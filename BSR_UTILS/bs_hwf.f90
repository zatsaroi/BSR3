!=========================================================================
      Program bs_hwf
!=========================================================================
!     This utility provides hydrogen-like bound and continuum pseudo-state 
!     on the given B-spline basis (with output in BSR format)
!-------------------------------------------------------------------------
      Use spline_atomic
      Use spline_param
      Use spline_grid
      Use spline_hl
      Use spline_galerkin

      Implicit real(8) (A-H,O-Z)
      Character :: AF*80, ls, AL, EL4*4, ELF4*4, term*4
      Integer, allocatable :: ip(:), np(:), kp(:)
      Real(8), allocatable :: ap(:), cp(:)
      Real(8), allocatable :: hm(:,:), cm(:,:)
      Integer :: pri = 6

!------------------------------------------------------------------------------------
      iarg = COMMAND_ARGUMENT_COUNT()
      if(iarg.gt.0) Call GET_COMMAND_ARGUMENT(1,AF)

      if(iarg.le.0.or.AF.eq.'?') then
        write(*,'(/a)') 'This utility provides hydrogen-like bound and continuum pseudo-state'
        write(*,'(/a)') 'on the given B-spline basis (with output in BSR format) '
        write(*,'(/a)') 'z-value comes from the knot.dat'
        write(*,'(/a)') 'Call as:  bs_hwf  l=.. [nsol=.. ii=.. jj=.. eps_tail=.. emax=..] '
        write(*,'(/a)') 'l     -  orbital momentum '
        write(*,'( a)') 'nsol  -  number of solutions [1]'
        write(*,'( a)') 'ii    -  zero B-splines in the beginning [l+1] '
        write(*,'( a)') 'jj    -  zero B-splines in the end [1] '
        write(*,'( a)') 'emax  -  maximum energy [100] '
        write(*,'( a)') 'eps_tail -  tolerence for the last B-spline [1d-7] '
        write(*,'(/a)') 'output files:  nl_###.bsw and nl_###.c, with ### - index of solution'
        Stop 
      end if        
!---------------------------------------------------------------------
      Open(pri,file='bs_hwf.log')

      Call define_grid (z);  Call define_spline

      Allocate(hm(ns,ns),cm(ns,ns))
      Allocate(ip(ns),np(ns),kp(ns),ap(ns),cp(ns))
      
!----------------------------------------------------------------------
! ,,, read parameters:


      l = 0;         Call Read_iarg('l',l)
      ii = -1;       Call Read_iarg('ii',ii) 
      jj =  1;       Call Read_iarg('jj',jj)
                                                                  
      eps_tail = 1d-7;  Call Read_rarg('eps_tail',eps_tail)
      emax = 100d0;     Call Read_rarg('emax',emax)

      nsol = 1;     Call Read_iarg('nsol',nsol) 
   
      write(pri,'(a,i6)')   'l   =',l
      write(pri,'(a,f6.1)') 'z   =',z
      write(pri,*)
      write(pri,'(a,i6)')  'ks  =',ks
      write(pri,'(a,i6)')  'ns  =',ns
      write(pri,*)

      write(pri,'(a)') 'Boundary conditions at r = 0:'
      if(ii.eq.-1) then
       ii = 1
       if(l+1.lt.ks)  ii = l+1
      end if
      write(pri,*)
      write(pri,'(a,i2,a)') 'ii =',ii,'  - number of first zero B-splines'

      write(pri,*)
      write(pri,'(a)') 'Boundary conditions at r = a:'
      write(pri,*)
      write(pri,'(a,i2,a)') 'jj =',jj,'  - number of last zero B-splines' 

! ... some initiations:

      nhm = ns      !  full size of interaction matrix:
      ip = 1        !  pointer on zero B-splines
      emax = emax * z * z

!----------------------------------------------------------------------
! ... Coulomb Matrix:

      Call HLM(l)
      hl = -0.5d0 * hl
      Call full_mat_sym(ns,ks,hl,hm,'l')
      Call full_mat_sym(ns,ks,sb,cm,'l')

!----------------------------------------------------------------------
! ... boundary conditions at  r=0:

      if(ii.gt.0)  ip(1:ii) = 0

! ... boundary conditions at  r=a:

      if(jj.gt.0)  ip(ns-jj+1:ns) = 0

! ... apply zero condition:

      m=0
      Do j=1,nhm
       if(ip(j).eq.0) Cycle; m=m+1
       k=0; ap = 0.d0; cp = 0.d0
       Do i=1,nhm
        if(ip(i).eq.0) Cycle
        k=k+1; ap(k)=hm(i,j); cp(k)=cm(i,j)
       End do
       hm(1:nhm,m)=ap(1:nhm); cm(1:nhm,m)=cp(1:nhm)
      End do 
      khm=m                                                 

!----------------------------------------------------------------------
! ... diagonalization of interaction matrix:

      Call LAP_DSYGV('V','L',khm,nhm,hm,cm,cp,info)         

      write(pri,'(/75a1)') ('-',i=1,75)
      write(pri,'(/a,i5,a)') 'khm =',khm,' - final full size of matrix'      
      write(pri,'(/75a1)') ('-',i=1,75)

!----------------------------------------------------------------------
! ... compare the energies with exact Coulomb values:

      write(pri,'(/a/)') 'Comparison with Coulomb eigenvalues:'
      write(pri,'(/3x,a,4x,a,3(6x,a,6x)/)') &
            '# ','n', '  E', '   EC ', '  acr '

      np = ICHAR('k')-ICHAR('1')+1
      kp = 1;  kn=1; kk=1

      Do i=1,khm
       if(cp(i).lt.0.d0) then
       S=abs(cp(i)); S = sqrt((z*z)/(2*S)); n=NINT(S)
       if(n.le.l) Stop 'n <= l'
       EE = -z*z / (2.d0*n*n)
       err = abs((cp(i)-EE)/EE)
       write(pri,'(2i5,3f16.10)')  i, n, cp(i), EE, err

       if(err.lt.1d-2.and.n.lt.10) then
        np(i) = n
        kp(i) = 1
       else
        np(i) = ICHAR('n')-ICHAR('1')+1
        kp(i) = kn; kn=kn+1
       end if

       else
        kp(i) = kk; kk=kk+1
       end if

      End do

      write(pri,'(/75a1)') ('-',i=1,75)

!----------------------------------------------------------------------
! ... restore the solutions in original B-spline net:

      Do j=1,khm;   k=0; ap=0.d0
       Do i=1,nhm
        if(ip(i).eq.0) Cycle;  k=k+1; ap(i)=hm(k,j)
       End do
       hm(:,j)=ap(:)
      End do 

!----------------------------------------------------------------------
! ... output solutions in BSR format:

      ls = AL(l,1)
      write(term,'(1x,i1,a1,i1)') 2,AL(l,2),1

      if(nsol.gt.khm) khm=nsol
      Do is = 1,nsol 

! ...  bsw - file:

       ap(:) = hm(:,is)
       m = ns       
       Do i=ns,1,-1
        if(abs(ap(i)).gt.eps_tail) Exit
        m=i-1
       End do

       write(AF,'(3a,i3.3,a)') 'n',ls,'_',is,'.bsw'
       nu=1
       open(nu,file=AF,form='UNFORMATTED')

       EL4 = ELF4(np(is),l,kp(is))

       write(nu) EL4,z,h,hmax,rmax,ks,ns,m
       write(nu) ap(1:m)
       Close(nu)

       write(pri,'(a4,f16.8,5x,a,i5)')  EL4,cp(is),trim(AF),m

! ...  c - file:

       write(AF,'(3a,i3.3,a)') 'n',ls,'_',is,'.c'
       nu=1
       open(nu,file=AF)
       write(nu,'(15x,f16.8)') cp(is)
       write(nu,*)
       write(nu,'(a4,a4,56x,f11.8)') EL4,'( 1)',1.d0
       write(nu,'(a,a)') term,term        
       write(nu,'(a)') '*'
       close(nu)

      End do

      End program bs_hwf


