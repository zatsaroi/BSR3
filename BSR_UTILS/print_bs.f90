!======================================================================
    Program print_bs
!======================================================================
!   This utilty prints B-splines (in separate files bs_###)
!   for given "knot.dat" 
!----------------------------------------------------------------------
    Use spline_atomic
    Use spline_param
    Use spline_grid    

    Implicit real(8) (A-H,O-Z)
    Character(20) :: AF
    Real(8), allocatable :: c(:)

    ! .. sets up the grid points:  array t in MOLULE spline_grid

    CALL define_grid (z)

    ! .. initializes the values of the spline and its derivatives
    ! .. and evaluates the spline arrays (operators in spline basis)
    ! .. which are defined in the MODULES spline_grid and spline_galerkin

    CALL define_spline

    ! .. print B-splines in separate files:
     
    Allocate(c(ns))

    n = 100;  Call Read_iarg('n',n)

    Do is = 1,ns
     write(AF,'(a,i3.3)') 'bs_',is
     nu = 1;  open(nu,file=AF)

     c = 0.d0;  c(is)=1.d0

     ! ... initial and final points of the spline:

     x0 = t(is);  x1 = t(is+ks)
     h = (x1-x0)/n
     y0 = 0.d0; if(is.eq.1) y0=1.d0
     write(nu,'(2e16.8)') x0, y0

     Do i = 1,n
      x = x0 + i*h
      y = bvalu2 (t, c, ns, ks, x, 0)
      write(nu,'(2e16.8)') x,y
     End do

    End do

    End program print_bs


