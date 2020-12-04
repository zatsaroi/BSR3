!=====================================================================
      Subroutine hf_nr (i,hfm,v,eii)
!=====================================================================
!     Improve the estimate v by the Newton_Raphson method
!     subject to orthogonality.
!---------------------------------------------------------------------
      Use bsr_hf 
      Use hf_orbitals

      Integer, intent(in) :: i
      Real(8), intent(in) :: hfm(ns,ns)
      Real(8), intent(inout) :: v(ns)
      Real(8), intent(out) :: eii

      Real(8), dimension(ns+nbf,ns+nbf) :: aa, ax
      Real(8), dimension(ns+nbf) :: res, xx,yy
      Integer, dimension(ns+nbf) :: ipiv

      Real(8), external  :: a
      Real(8) :: x(ns), y(ns), w(ns)
      Real(8) :: d(ns,ks), dd(ns,ks), hx(ns,ns)
      Real(8) :: eij
      Integer :: md, mm, k, info, j,jp

      md = ns+nbf
      v = p(:,i)

! ... compute the residuals     

      aa = 0.d0
      x = matmul(hfm,v)
      eii = dot_product(x,v)
      aa(1:ns,1:ns) = (hfm -eii*bb)
      res(1:ns) =  matmul(aa(1:ns,1:ns),v)

! ... add normalization constraint for orbital i

      mm = ns+1
      y = matmul(BB,v)
      aa(1:ns,mm) = -y
      aa(mm,1:ns) = -y
      res(mm) = 0.d0
        
! ... add orthogonality constraints with orbital j    (to all ???)

      Do j = 1,nbf 
       if(lbs(i).ne.lbs(j)) Cycle
       if(j.eq.i) Cycle
       mm = mm + 1
       w = p(:,j)
       eij = dot_product(w,x)
       y = matmul(BB,w)
       res(mm) =  0.d0
       aa(1:ns,mm) = -y
       aa(mm,1:ns) = -y
       res(1:ns) = res(1:ns) + eij*y
      End do
   
! ... Add Haa contributions
        
      dd = 0.d0
      Do k = 0, 2*lbs(i), 2
       Call density(ns,ks,d,p(1,i),p(1,i),'s')
       dd = dd + (2.d0*a(i,i,k))*d
      End do
      Call Convol(ns,ks,d,dd,1,'s','s')
      hx = 0.d0
      Call Update_hs(ns,ks,hx,d,'s')          
      aa(1:ns,1:ns) = aa(1:ns,1:ns) + hx
    
! ... Apply boundary conditions (delete the extra B-splines)

      ii=0
      Do j=1,mm
       if(j.le.ns) then; if(iprm(j,i).eq.0) Cycle; end if
       ii=ii+1
       k=0
       Do jp=1,mm
       if(jp.le.ns) then; if(iprm(jp,i).eq.0) Cycle; end if
        k=k+1; xx(k)=aa(jp,j)
       End do
       ax(1:k,ii)=xx(1:k); yy(ii)=res(j)
      End do 

! ... solve the matrx equation with LAPACK routines:

      CALL DGETRF(nn, nn, ax, md, ipiv, info)
      if (info /= 0) Stop 'hf_nr: error in factorization routine DSYTRF'
      CALL DGETRS('N', nn, 1, ax, md, ipiv, yy, md, info)
      if (info /= 0) Stop 'hf_nr: error in solve routine DSYTRS'

! ... restore the solution in orginal basis:

      res=0.d0; k=0
      Do j=1,md
       if(j.le.ns) then; if(iprm(j,i).eq.0) Cycle; end if
       k=k+1; res(j)=yy(k)
      End do 
        
! ... normalize the solution

      v = v - res(1:ns)
      y = matmul(BB,v)
      v = v/sqrt(Dot_Product(v,y))
      eii = eii - res(ns+1)
 
      End Subroutine hf_nr
