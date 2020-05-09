!======================================================================
      Module Z_core
!======================================================================
!     contain B-splne representation for Z-integrals
!     including interaction with core
!----------------------------------------------------------------------
!     zl(ns,ks) - matrix of Z-operator in the B-splne basis
!                 (in symmetric lower-column storage)
!
!     zl_dir(ns,ks) - direct interaction with core
!                     (in symmetric upper-column storage)
!
!     zl and zl_dir have different representation !!!
!
!     lz_dir = [-1|1] - (not|yes) for done of zl_dir
!     this parameter allows us to compute the zl_dir only once 
!
!     zl and zl_dir does not depend on l (can be calculated only once) 
!
!     zl_exc(ns,ns) - exchange interaction with core
!
!     zl_core(ns,ns) - final zl with core interaction
!                      (in full representation)
!
!     mlz - max. l for so-interaction
!
!     zl_full(ns,ns,mlz) - all z-matreces
!
!     Z_int - Z-integrals <P_i| Z |P_j>
! 
!     Z_vec - vectors  < . | Z | P_i>
!
!     d,x - auxiliary arrays
!
!------------------------------------------------------------------
      Implicit none

      Integer :: lz_dir,mlz

      Real(8), allocatable :: zl(:,:), zl_dir(:,:), zl_exc(:,:)
      Real(8), allocatable :: zl_core(:,:), zl_full(:,:,:)
      Real(8), allocatable :: d(:,:),x(:,:), dd(:,:),xx(:,:),yy(:,:)

      Integer :: iz1,iz2      !  Z_int under consideration
      Real(8) :: ci 
      
      Integer :: iv           !  Z_vec under consideration
      Real(8), allocatable :: zv(:)

      End MODULE Z_core


!=====================================================================
      Subroutine Alloc_zcore(ns,ks)
!=====================================================================
!     allocate arrays in Z_core module
!---------------------------------------------------------------------
      Use Z_core

      Implicit none
      Integer, intent(in) :: ns,ks

      if(ns.gt.0) then

       if(allocated(zl)) Return

       Allocate(zl(ns,ks), zl_dir(ns,ks), &
                zl_exc(ns,ns), zl_core(ns,ns), &
                zl_full(ns,ns,mlz), &
                d(ns,ks), dd(ns,ns), x(ns,ks), &
                xx(ns,ns), yy(ns,ns), zv(ns))

       d=0.d0; dd=0.d0; x=0.d0; xx=0.d0; yy=0.d0
       lz_dir = -1
       iz1=0; iz2=0;  ci=0.d0
       iv=0; zv=0.d0 

      else

       if(.not.allocated(zl)) Return

       Deallocate(zl,zl_dir,zl_exc,zl_core,zl_full, &
                  d,dd,x,xx,yy, zv)

      end if

      End Subroutine Alloc_zcore


!======================================================================
      Subroutine Gen_Zval(m)
!======================================================================
!     Generates the Z-integrals and vectors 
!----------------------------------------------------------------------
      Use bsr_mat 
      Use Z_core 

      Implicit none
      Integer :: l,m

! ... find max. l and allocate arrays:      

      mlz=m

      Call Alloc_zcore(ns,ks)

! ... generate Z-data for given l: 

      Do l=1,mlz
       Call INTZ_core(l);   zl_full(:,:,l)=zl_core(:,:)               
      End do 

      End Subroutine Gen_Zval


!=====================================================================
      Subroutine INTZ_core(l)
!=====================================================================
!     generate B-splne representation in array zl_core 
!     for Z-integrals with core interaction for given l 
!---------------------------------------------------------------------
      Use bsr_mat, xb => x
      Use Z_core

      Implicit none
      Integer, intent(in) :: l
      Integer :: ip, jp, k, kmin, kmax, lp, i,j
      Real(8) :: C
      Real(8), External :: ZCB

!----------------------------------------------------------------------
! ... check the l-independent part - direct interaction with core:

      if(lz_dir .eq. -1) then

       Call Gen_zzeta(zl,izcorr) 
       zl_dir = 0.0d0
       if(kclosd.gt.0) then
        Call MNK_cell(0)
        Do ip = 1,kclosd
         C = -(4*lbs(ip)+2)
         Call Density (ns,ks,d,pbs(1,ip),pbs(1,ip),'s')
         Call Convol  (ns,ks,x,d,1,'s','s')               
         zl_dir = zl_dir + C*x
        End do
       end if
       lz_dir = 1       

      end if

!----------------------------------------------------------------------
! ... exchange interaction ...
      
      zl_exc = 0.d0
      if(kclosd.gt.0) then

      kmin = 1000; kmax = 0
      Do ip = 1,kclosd
       k = iabs(l-lbs(ip));  if(kmin.gt.k) kmin=k
       k =      l+lbs(ip) ;  if(kmax.lt.k) kmax=k
      End do

      Do k = kmin,kmax                     ! Mk - integrals
        Call MMK_cell(k)
        Do ip = 1,kclosd
         lp = lbs(ip)
         C = 3 * (4*lp+2) * ZCB(lp,k,l) / (8*l*(l+1)) 
         if(C.eq.0.d0) Cycle
         C = ( (l*(l+1)-lp*(lp+1)) * (l*(l+1)-lp*(lp+1)+k*(k+1)) &
           +   ((l+lp+1)**2-(k+1)**2) * ((k+1)**2-(l-lp)**2) )*C/(k+1)
         if(C.eq.0.d0) Cycle
         Call Density (ns,ks,dd,pbs(1,ip),pbs(1,ip),'x')
         Call Convol  (ns,ks,xx,dd,3,'s','s')               
         zl_exc = zl_exc + C*xx
        End do
      End do    

      Do k = kmin,kmax                    ! M(k-2) - integrals
        if(k.le.0) Cycle
        Call MMK_cell(k-2)
        Do ip = 1,kclosd
         lp = lbs(ip)
         C = 3 * (4*lp+2) * ZCB(lp,k,l) / (8*l*(l+1))
         if(C.eq.0.d0) Cycle
         C = -( (l*(l+1)-lp*(lp+1)) * (l*(l+1)-lp*(lp+1)+k*(k+1)) &
             +  ((l+lp+1)**2-k**2) * (k**2-(l-lp)**2) )*C/k
         Call Density (ns,ks,dd,pbs(1,ip),pbs(1,ip),'x')
         Call Convol  (ns,ks,xx,dd,3,'s','s')               
         zl_exc = zl_exc + C*xx       
       End do
      End do    

      Do k = kmin,kmax                     ! W(k-1) - integrals
        if(k.le.0) Cycle
        Call MWK_cell(k-1)
        Do ip = 1,kclosd
         lp = lbs(ip)
         C = 3 * (4*lp+2) * ZCB(lp,k,l) / (8*l*(l+1))
         if(C.eq.0.d0) Cycle
         C = C * (l*(l+1)-lp*(lp+1)+k*(k+1)) 
         if(C.eq.0.d0) Cycle
         Call Density (ns,ks,dd,pbs(1,ip),pbs(1,ip),'x')
         Call Convol  (ns,ks,xx,dd,3,'n','s')               
         Call Convol  (ns,ks,yy,dd,4,'n','s')               
         zl_exc = zl_exc + C*(xx-Transpose(yy))
       End do
      End do    

      end if ! over kclosd

!----------------------------------------------------------------------
! ... total zl:      

      zl_core = zl_exc

      Do jp = 1,ks
       Do i = ks+1-jp,ns
        j = jp+i-ks
        zl_core(i,j) = zl_core(i,j) + zl(i,jp) 
        zl_core(j,i) = zl_core(i,j)
       End do
      End do

      Do jp = 1,ks
       Do i = 1,ns-jp+1
        j = i+jp-1
        zl_core(i,j) = zl_core(i,j) + zl_dir(i,jp)
        zl_core(j,i) = zl_core(i,j)
       End do
      End do

      End Subroutine INTZ_core


!====================================================================
      SUBROUTINE Gen_zzeta(rm,izcorr)
!====================================================================
!     B-spline matrix for Z*fine*1/r^3, Used for spin-orbit interaction;
!     if izcorr > 0,  multiplied with correction r/(r+2*z*fine) ^ 2
!--------------------------------------------------------------------
      Use spline_param; Use spline_grid;  Use spline_atomic

      Implicit none
      Integer :: m, ith,jth, i, irow,jcol, izcorr
      Real(8) :: S,B, rm(ns,ks)

      rm = 0.d0
      do ith=1,ks
       do jth=1,ith
        jcol=jth-ith+ks
        do i=1,nv
         irow=i+ith-1
         do m = 1,ks
          B = gr(i,m)/(gr(i,m)+2*fine*Z)
          S = grw(i,m)*bsp(i,m,ith)*grm(i,m)**3*bsp(i,m,jth)
          if(izcorr.gt.0) S = S * B * B
          rm(irow,jcol)=rm(irow,jcol)  + S  
         end do
        end do
       end do
      end do
      rm = rm * Z * fine

      END SUBROUTINE Gen_zzeta


!======================================================================
      Real(8) Function Z_int(i,j)
!======================================================================
! ... define integral  <P_i|Z|P_j>  for the given zl_core
!---------------------------------------------------------------------    
      Use bsr_mat
      Use Z_core 
      
      Implicit none
      Integer, intent(in) :: i,j
      Integer :: m,k

      if(i.eq.iz1.and.j.eq.iz2) then
       Z_int = ci; Return
      end if

      ci = 0.d0 
      Do m = 1,ns;  Do k = 1,m-1
       ci = ci + zl_core(m,k)*(pbs(m,i)*pbs(k,j)+pbs(k,i)*pbs(m,j))
      End do;  End do
      Do m = 1,ns
       ci = ci + zl_core(m,m)*pbs(m,i)*pbs(m,j)
      End do

!      if(lbs(i).eq.1) ci = ci * zcorr           !  ???
!      ci = ci * zcorr
      Z_int=ci; iz1=i; iz2=j

      if(nud.gt.0) Call pri_int(nud,7,0,i,j,i,j,ci)

      End Function Z_int


!======================================================================
      Subroutine Z_vec(i,v)
!======================================================================
! ... define vector  <P_i|Z|.>  for the given zl_core
!---------------------------------------------------------------------    
      Use bsr_mat
      Use Z_core 

      Implicit none
      Integer, intent(in) :: i
      Real(8) :: v(ns)
      Integer :: m

      if(i.eq.iv) then; v = zv; Return; end if

      Do m=1,ns; v(m)=SUM(zl_core(:,m)*pbs(:,i)); End do

      iv=i; zv=v

      End Subroutine Z_vec
