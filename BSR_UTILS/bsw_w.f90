!======================================================================
!     UTILITY       B S W _ W
!
!               C O P Y R I G H T -- 2003
!
!     Written by:   Oleg Zatsarinny
!======================================================================
!     connvert the B-spline results in the MCHF format
!======================================================================
!     INPUT FILES:   name.bsw
!     OUTPUT FILES:  name.w
!     Call as:  bsw_w name.bsw
!----------------------------------------------------------------------
      Use spline_param; Use spline_atomic; Use spline_grid
      Use spline_galerkin

      Implicit real(8) (A-H,O-Z)
      Character(3) :: EL3
      Character(3), external :: ELF3
      Character(4) :: EL4
      Character(6) :: Atom,Term
      Character(40) :: AF,BF
      Real(8), allocatable :: v(:), w(:), R(:), R2(:)
      
      iarg = COMMAND_ARGUMENT_COUNT() 
      if(iarg.gt.0) Call GET_COMMAND_ARGUMENT(1,AF) 
     
      if(iarg.eq.0.or.AF.eq.'?') then
        write(*,'(/a)') 'bsw_w  converts  name.bsw  to name.w'
        write(*,'(/a)') 'w   - unformatted MCHF-CFF format for radial functions'
        write(*,'(/a)') 'bsw - BSR format for radial functions'
        write(*,'(/a)') 'Call as:  bsw_w  name.bsw '
        write(*,'(/a)') 'Results:  name.w'
        write(*,'(/a)') 'Warning:  set index in bsw-files may be too big for name.w format'
        Stop 
      end if

!----------------------------------------------------------------------
! ... input data: 
         
      i=index(AF,'.') 
      if(AF(i+1:i+3).ne.'bsw') Stop 'input file should have bsw-extention'

      inp=1; Open(inp,file=AF,status='OLD',form='UNFORMATTED')

      iout=2; BF=AF(1:i)//'w';  Open(iout,file=BF,form='UNFORMATTED')

!----------------------------------------------------------------------
! ... sets up grid points and initializes the values of the spline: 
    
      read(inp) el4,z,h,hmax,rmax,ks,ns; rewind(inp)

      CALL define_grid(z);  CALL define_spline

      NR=220; RHO=-4.d0; HR=1./16.D0

      Allocate(R(NR),R2(NR),v(ns),w(NR))

      DO i=1,NR; R(I)=EXP(RHO+(I-1)*HR)/Z; R2(I)=SQRT(R(I)); END DO

      Rm = t(ns+1);  mx = 0;   Do i=1,NR; if(R(i).lt.Rm) mx=i; End do
      mx = mx + 1; 
!----------------------------------------------------------------------
!                                                          
    1 read(inp,end=2) el4,zw,hw,hmw,rmw,ksw,nsw,m
      read(inp) v(1:m); if(m.lt.ns) v(m+1:ns)=0.d0

      Call EL4_nlk(el4,n,l,k); 

      EL3 = ELF3(n,l,k)

      Do j=1,mx; w(j) = bvalu2(t,v,ns,ks,r(j),0)/r2(j); End do

      AZ = azl(z,h,ks,l+1) * v(l+2)
              
      Atom = 'atom'
      TERM = 'term'
      EI = 0.d0
      ZI = 0.d0        
        
      write(iout) Atom, Term, EL3, mx, Z,EI,ZI,AZ, w(1:mx)

      S = BVMV (ns,ks, sb,'s',v,v)
      write(*,'(3a10,f16.8)') el4,el3,'norma = ',S

      go to 1               
    2 Continue
      
      END   !  program bsr_w
          


!=====================================================================
    SUBROUTINE mkgridb
!=====================================================================

!   sets up the knots for spline

    USE spline_param;  USE spline_grid; USE spline_atomic

      IMPLICIT REAL(8) (A-H,O-Z)


!      INTEGER, INTRINSIC:: NINT
      INTEGER:: n, i, m, me1, me2
!      REAL(KIND=8), INTRINSIC:: LOG, MAX
      REAL(KIND=8):: hp1, h2, tmax, tx, h0

      ! .. determine ml, the number of distinct points from 0 to 1

       h0 = h
       ml = 1.d0/h + 0.1
       h = 1.d0/ml
       hp1 = 1.d0 + h

      ! .. determine tmax

      tmax = z*rmax
      hmax = z*hmax

      ! .. determine final point of "exponential" grid
      ! .. me: number of points from 1 to (1+h)**me
      ! .. m:  number of points from (1+h)**me to tmax

      me1 = MAX(0.0d0, LOG(hmax/h)/LOG(hp1)+1.d0)
      me2 = LOG(tmax)/LOG(hp1)+1

      IF ( me2 <= me1 ) THEN
        me = me2
        m = 0
      ELSE
        me = me1
        tx = hp1**me
        h2 = h*tx/hp1
        m = NINT((tmax-tx)/h2)
      END IF

      n = ml + me + m + ks -1
      ns = n
      nv = ns - (ks -1)
      nt = ns + ks

      ! .. establish the grid for z*r

      ALLOCATE (t(nt)); t = 0.d0

      DO i = ks+1, ks+ml
        t(i) = t(i-1) + h0
      END DO

      DO i = ks+ml+1, ks+me+ml
        t(i) = t(i-1)*hp1
      END DO

      DO i = ks+me+ml+1, n+1
        t(i) = t(i-1) + h2
      END DO
      t(n+2:nt) = t(n+1)

      ! .. scale the values to the R variable

      t = t/z
      hmax = hmax/z

    END SUBROUTINE mkgridb


   