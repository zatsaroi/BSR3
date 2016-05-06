!======================================================================
      MODULE Z_core
!======================================================================
!     contains B-splne representation for Z-integral
!     including the interaction with core
!----------------------------------------------------------------------
!     zl(ns,ks) - matrix of Z-operator in the B-splne basis
!                 (in symmetric lower-column mode)
!
!     zl_dir(ns,ks) - direct interaction with core
!                     (in symmetric upper-column mode)
!
!     zl and zl_dir have different representation !!!
!
!     lz_dir = [-1|1] - (not|yes) for Done of zl_dir
!     (this parameter allows us to compute the zl_dir only once because 
!!     zl_dir does not depend on l) 
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
!------------------------------------------------------------------

      Implicit none
      Integer :: lz_dir, mlz
      Real(8), allocatable :: zl(:,:), zl_dir(:,:), zl_exc(:,:), &
                              zl_core(:,:), zl_full(:,:,:),      &
                              Z_int(:,:), Z_vec(:,:),            &
                              d(:,:),x(:,:),dd(:,:),xx(:,:),yy(:,:)
      End MODULE Z_core


!=====================================================================
      Subroutine Alloc_zcore(ns,ks,nbf)
!=====================================================================
!     allocate arrays in module Z_core
!---------------------------------------------------------------------
      Use Z_core

      Implicit none
      Integer, intent(in) :: ns,ks,nbf

      if(ns.gt.0) then

       if(allocated(zl)) Return

       Allocate(zl(ns,ks), zl_dir(ns,ks), &
                zl_exc(ns,ns), zl_core(ns,ns), &
                Z_int(nbf,nbf), Z_vec(ns,nbf), &
                zl_full(ns,ns,mlz), &
                d(ns,ks), dd(ns,ns), x(ns,ks), &
                xx(ns,ns), yy(ns,ns))

       d=0.d0; dd=0.d0; x=0.d0; xx=0.d0; yy=0.d0
       Z_int = 0.d0; Z_vec = 0.d0
       lz_dir = -1

      else

       if(.not.allocated(zl)) Return
       Deallocate(zl,zl_dir,zl_exc,zl_core,zl_full,Z_int,Z_vec, &
                  d,dd,x,xx,yy)
       lz_dir = -1

      end if

      End Subroutine Alloc_zcore


!======================================================================
      SUBROUTINE Gen_Zval
!======================================================================
!     Generates the Z-integrals and vectors 
!----------------------------------------------------------------------
      Use bsr_mat
      Use Z_core 
      Use spline_param; Use spline_atomic; Use spline_orbitals

      Implicit none
      Integer :: i,j,k,l,m
      Real(8) :: C

      mlz=maxval(lbs(1:nbf));  if(mlz.gt.mlso) mlz=mlso
      Call Alloc_zcore(ns,ks,nbf)

! ... generate Z-data for given l: 

      Do l=1,mlz

       Call INTZ_core(l)

       Do i = kclosd+1,nbf;  if(lbs(i).ne.l.or.iech(i).ne.0) Cycle 
        Do j = i,nbf;        if(lbs(j).ne.l.or.iech(j).ne.0) Cycle 

         C = 0.d0 
         Do m = 1,ns;  Do k = 1,m-1
          C = C + zl_core(m,k)*(pbs(m,i)*pbs(k,j)+pbs(k,i)*pbs(m,j))
         End do;  End do
         Do m = 1,ns
          C = C + zl_core(m,m)*pbs(m,i)*pbs(m,j)
         End do

         if(lbs(i).eq.1) C = C * zcorr           !  ???
 
         Z_int(i,j) = C;  Z_int(j,i) = C
 
         if(nud.gt.0) Call pri_int(nud,7,0,i,j,i,j,C)

        End do
        Do m=1,ns; Z_vec(m,i)=SUM(zl_core(:,m)*pbs(:,i)); End do
       End do

       zl_full(:,:,l) = zl_core(:,:)               

      End do 

      End Subroutine Gen_Zval


!=====================================================================
      Subroutine INTZ_core(l)
!=====================================================================
!     generate B-splne representation in array zl_core 
!     for Z-integrals with core interaction for given l 
!---------------------------------------------------------------------
      Use bsr_mat;  Use Z_core
      Use spline_param; Use spline_atomic; Use spline_orbitals
      Use spline_grid; Use spline_galerkin

      Implicit none
      Integer, intent(in) :: l
      Integer :: ip, jp, k, kmin, kmax, lp, i,j
      Real(8) :: C
      Real(8), external :: ZCB

!----------------------------------------------------------------------
! ... check the l-independent part direct interaction with core:

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
!     if izcorr > 0,  multiplied with correction r/(r+2*z*fine)
!--------------------------------------------------------------------
      Use spline_param; Use spline_grid;  Use spline_atomic

      Implicit none
      Integer :: m, ith,jth, i, irow,jcol, izcorr
      Real(8) :: S,B, rm(ns,ks)

      rm = 0.d0

      Do ith=1,ks
       Do jth=1,ith
        jcol=jth-ith+ks
        Do i=1,nv
         irow=i+ith-1
         Do m = 1,ks
          B = gr(i,m)/(gr(i,m)+2*fine*Z)
          S = grw(i,m)*bsp(i,m,ith)*grm(i,m)**3*bsp(i,m,jth)
          if(izcorr.gt.0) S = S * B * B
          rm(irow,jcol)=rm(irow,jcol)  + S  
         End do
        End do
       End do
      End do

      rm = rm * Z * fine

      End Subroutine Gen_zzeta

