!======================================================================
!     UTILITY   W _ B S W
!
!               C O P Y R I G H T -- 2004
!
!     Written by:   Oleg Zatsarinny
!                   free shooter
!                   email: oleg_zoi@yahoo.com
!======================================================================
!     Generation B-spline representation for orbitals 
!     given in MCHF package format, w-files (C. Froese Fisher)
!----------------------------------------------------------------------
!
!     INPUT ARGUMENTS:
!     
!     name.w              -  input w-file (or file with list of w-files)
!     Eps_end             -  tollerance for function tail (1.D-7)
!     klarg               -  index for Lagrange interpolation (20)
!
!     Should be at least first argument: if the name has extention .w
!     when it is considered as input w-file, otherwise it will be
!     considered as file with list of desirable w-files
!
!     INPUT FILES:
!
!     list_w              -  list of w-files (optional)
!     name.w              -  w-files  
!     knot.dat            -  parameters of B-splines
!
!     OUTPUT FILES:
!
!     w_bsw.log           -  running information
!     name.bsw            -  B-spline representation for orbitals
!     knot.dat            -  parameters of B-splines and knot points
!
!---------------------------------------------------------------------
      Use spline_param; Use spline_atomic; Use spline_grid
      Use spline_orbitals

      Implicit real(8) (A-H,O-Z)
      Real(8) :: Eps_end = 1.d-7
      Integer :: klagr = 20
      Character(80) :: AF, BF
      
      iarg = COMMAND_ARGUMENT_COUNT() 
      if(iarg.gt.0) Call GET_COMMAND_ARGUMENT(1,AF) 
     
      if(iarg.eq.0.or.AF.eq.'?') then
        write(*,'(/a)') 'w_bsw  converts  name.w  to name.bsw'
        write(*,'(/a)') 'w   - unformatted MCHF-CFF format for radial functions'
        write(*,'(/a)') 'bsw - BSR format for radial functions'
        write(*,'(/a)') 'Call as:  w_bsw  name.w [eps_end  klagr]'
        write(*,'(/a)') 'Results:  name.bsw'
        write(*,'(/a)') 'eps_end - tolerance for function tail [1.D-7]'
        write(*,'(/a)') 'klagr   - index for Lagrange interpolation [20]'
        Stop
      end if

!----------------------------------------------------------------------
! ... sets up grid points and initializes the values of the spline: 
    
      Call define_grid (z)
      Call define_spline
      Call Allocate_bsorb(ibf)

!----------------------------------------------------------------------
! ... read arguments if any:


      Call GETARG(1,AF) 
      if(iarg.ge.2) then; Call GETARG(2,BF); read(BF,*) Eps_end; end if
      if(iarg.ge.3) then; Call GETARG(3,BF); read(BF,*) klagr; end if

      inp = 3; i=LEN_TRIM(AF); if(AF(i-1:i).eq.'.w') inp = 0
      if(inp.gt.0) Open(inp,file=AF)

      ipri=66
      Open(ipri,file='w_bsw.log',position='APPEND')
      write(ipri,'(/a/)') 'Conversion to B-spline reprisantation:'

      Do

      if(inp.gt.0) read(inp,*,end=1) AF; if(AF(1:1).eq.'*') Exit 

      ii=INDEX(AF,'.')-1;  if(ii.lt.0) ii=LEN_TRIM(AF)
      AF = AF(1:ii)//'.w'
      BF = AF(1:ii)//'.bsw'

      nuw=1; Open(nuw,file=AF,form='UNFORMATTED',status='OLD')
      nub=2; Open(nub,file=BF,form='UNFORMATTED') 
 
      write(ipri,'(/a,a)') 'Input file: ',AF
      write(ipri,'( a,a)') 'Output file: ',BF

      write(ipri,'(/a,D12.3)') 'Eps_end tolerence: ',Eps_end
      write(ipri,'(a,i5/)') 'Lagrange interpolation index: ',klagr

      write(ipri,'(a,F12.3/)') 'Border radius: ',t(ns+1)
      
      nbf = 0
      Call CONVERT_MCHF(nuw,ipri,klagr,Eps_end)
  
      Do i=1,nbf
       write(nub) ebs(i),z,h,hmax,rmax,ks,ns,mbs(i)
       write(nub) pbs(1:mbs(i),i)
      End do

      if(inp.eq.0) Exit
      End do

    1 Continue

      End    !  program w_bsw


!======================================================================
      Subroutine Convert_mchf(nuw,pri,klagr,Eps_end)
!======================================================================
!     convert w -> bsw for one w-file, unit 'nuw'    
!----------------------------------------------------------------------
      Use spline_param; Use spline_atomic; Use spline_grid
      Use spline_orbitals

      Implicit real(8) (A-H,O-Z)
      Integer, intent(in) :: nuw,pri,klagr
      Real(8), intent(in) :: Eps_end

!----------------------------------------------------------------------
! ... MCHF radial scale:
 
      Real(8), parameter :: HR = 1.d0/16.d0
      Real(8), parameter :: RHO = -4.d0
      Integer(4), parameter :: NOD = 220
      Character(6) :: Atom, Term
      Character(3) :: EL3
      Real(8) :: R(NOD), R2(NOD), P(NOD) 

      RO = RHO
      Do I=1,NOD
       R(I)=DEXP(RO)/Z;  R2(I)=DSQRT(R(I)); RO=RO+HR
      End do 

!----------------------------------------------------------------------
    1 read(nuw,end=2) Atom,Term,EL3,m,ZZ,EI,ZI,AZ,P(1:m)
      if(zz.ne.z) Stop 'z in function  <>  z-knot'
      P(m+1:NOD) = 0.d0; P(:) = P(:) * R2(:) 
      Call EL3_nlk(EL3,n,l,k)      

! ... remove unphysical oscilations in the end:

      Do i=m,1,-1
       ii=i; if(abs(P(i)).gt.Eps_end) Exit; P(i)=0.d0
      End do
      m=ii

! ... add new orbital in bsw-list:

      i = Iadd_bsorb(n,l,k)
      if(nbf.eq.mbf) Call Allocate_bsorb(mbf+jbf)

! ... convert to B-spline basis

      Call CONVSPL (l,m,mbs(i),klagr,R,P,pbs(1,i),Eps_end,PB)

! ... check the normalization and max.deviation from original orbital:

      PN=sqrt(QUADR(i,i,0)); pbs(:,i)=pbs(:,i)/PN
 
      ydm = 0.d0
      Do j=1,m
       yb = bvalu2(t,pbs(1,i),ns,ks,r(j),0)
       yd = abs(yb-P(j)); if(yd.gt.ydm) ydm = yd
      End do

      write(pri,'(a4,2(a,F10.6))') ebs(i), &
            '   mdif =',ydm,'   border =',PB

      go to 1
    2 Close(nuw)

      End Subroutine Convert_mchf


!======================================================================
     Subroutine CONVSPL(l,nr,n,klagr,r,p,cb,Eps_end,PB)
!======================================================================
!
!    This programs converts numerical orbital p(1:nr) in radial
!    points r(1:nr) to B-spline representation cb(1:n)
!
!    the cubic-splines interpolation is used for first and last
!    b-spline interval, for others - Lagrange formular with index klagr
!
!    l - orbital momemtum of orbitals
!
!    a(1:l+1) => 0.d0 to garantee the correct behavior p ~ r^(l+1)
!                     at small r
!    a(ns-1:ns) = 0.d0 to garantee the bound character of orbital
!                      i.e. with zero value for function and derivative
!                      on the end of interval
!    Eps_end - tolerence for the last points 
!    PB - value of p(r) on the B_spline border
!
!----------------------------------------------------------------------
      Use spline_param; Use spline_grid; Use spline_galerkin
      
      Implicit real(8) (A-H,O-Z)
      Integer, intent(in)  :: l, nr, klagr
      Integer, intent(out) :: n
      Real(8), intent(in)  :: r(nr),p(nr)
      Real(8), intent(out) :: cb(ns)
      Real(8), intent(in)  :: Eps_end      
      Real(8), intent(out) :: PB      
      Real(8), allocatable :: b(:), c(:), d(:), yg(:,:)

! ... allocate arrays:

      Allocate (b(nr), c(nr), d(nr), yg(nv,ks))

! ... interpolate function p(r) to cubic splines:

      Call Splin3(nr,r,p,b,c,d)

! ... evaluate the function in gaussian points:

      RM=r(nr); rk1= r(klagr/2); rk2= r(nr-klagr/2)
      nint = ns - ks + 1

      Do i=1,nint
        Do j=1,ks
         if(gr(i,j).gt.RM) then
           yg(i,j) = 0.0
         elseif(i.eq.1.or.i.eq.nint) then
           yg(i,j) = SEVAL(nr,gr(i,j),r,p,b,c,d)
         elseif(gr(i,j).lt.rk1.or.gr(i,j).gt.rk2) then
           yg(i,j) = SEVAL(nr,gr(i,j),r,p,b,c,d)
         else
           yg(i,j) = XLAGR(KLAGR,nr,r,p,gr(i,j))
         end if
         yg(i,j) = yg(i,j)*grw(i,j)
        End do
      End do
      PB = 0.d0;  if(RM.gt.t(ns+1))  PB = SEVAL(nr,t(ns+1),r,p,b,c,d)

! ... form the vector of inner products of the radial function 
! ... and the spline basis functions:

      Call VINTY (yg,cb)

! ... apply the boundary condition at the origin  (FACSB already
! ... made the zeros in the first row, but now in the END)

      ll=l; if(l.ge.ks) ll = 0

      Call FACSBL0 (ll)

      cb(1:1+ll)=0.d0; ! cb(ns-1)=0.d0; cb(ns)=0.d0

! ... solve the system of equations  Sb x = cb  

      Call dpbtrs ('U',ns,ks-1,1,bs,ks,cb,ns,ierr)  
      if(ierr.ne.0) Stop 'convspl: dpbtrs (LAPACK) failed'

      Do i=1,ns; if(t(i).gt.RM) cb(i)=0.d0; End do

      cb(ns-1)=0.d0; cb(ns)=0.d0

! ... find tail of function and remove too small values

      Do i=ns,1,-1
       n = i
       if(abs(cb(i)).lt.Eps_end) cb(i)=0.d0
       if(abs(cb(i)).ge.Eps_end) Exit
      End do

      Deallocate (b, c, d, yg)

! ... restore initial bs for other routines:

      Call FACSB 
      
      End Subroutine CONVSPL


!====================================================================
      Subroutine facsbl0(l)
!====================================================================
!     Sets up the overlap bs which is a transpose of sb, <B_i,B_j>,
!     with the correct boundary conditions, and then factorizes bs.
!
!     SUBROUTINES called:   DPBTRF (from LAPACK)
!--------------------------------------------------------------------
      Use spline_param; Use spline_galerkin

      Implicit none
      Integer, intent(in) :: l
      Integer :: i,j, ierr

      bs = TRANSPOSE(sb)

! ... apply zero boundary condition at r=0:

      Do i = 1,l+1
       bs(ks,i) = 1.d0
       Do j = 1,ks-1
        bs(j,ks-j+i) = 0.d0
       End do
      End do
 
      Call DPBTRF('U',ns,ks-1,bs,ks,ierr)
      if (ierr .ne. 0 )  Stop 'facsbl: dpbtrf (LAPACK) failed'

      End Subroutine facsbl0

