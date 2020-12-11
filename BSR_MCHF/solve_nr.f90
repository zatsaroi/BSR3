!=====================================================================
      Subroutine Solve_nr (i,hfm,v,rhs,hx)
!=====================================================================
!     Improve the estimate v by the Newton_Raphson method
!     subject to orthogonality.
!---------------------------------------------------------------------
      Use bsr_mchf
 
      Implicit none
      Integer :: i
      Real(8), intent(in)  :: hfm(ns,ns), hx(ns,ns), rhs(ns)
      Real(8) :: v(ns)

      Real(8), dimension(ns+nbf,ns+nbf) :: aa, ax
      Real(8), dimension(ns+nbf) :: res, xx,yy
      Integer, dimension(ns+nbf) :: ipiv

      Real(8), external  :: a
      Real(8) :: x(ns), y(ns), w(ns)
      Real(8) :: eii, eij
        
      Integer :: mm,md, k,  info, j, jp,ii, ipos(1)

      v = p(:,i)
      md = ns+nbf

! ... compute the residuals     

      res = 0.d0
      aa = 0.d0
      x = matmul(hfm,v) - rhs
      eii = dot_product(x,v)
      aa(1:ns,1:ns) = hfm - eii*BB
      res(1:ns) =  matmul(aa(1:ns,1:ns),v) - rhs
      if(debug.gt.0) &
       write(log,'(a,i4,e12.5)') 'hf_nr: i,eii',i,eii

! ... add normalization constraint for orbital i

      mm = ns + 1
      y = matmul(BB,v)
      aa(1:ns,mm) = -y
      aa(mm,1:ns) = -y
        
! ... add orthogonality constraints with orbital j

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
       if(debug.gt.0) &
        write(log,'(a,i4,2e12.5)') 'hf_nr: j,eij',j,eij,dot_product(v,y)

      End do

! ... add H_aa contributions
        
      aa(1:ns,1:ns) = aa(1:ns,1:ns) + hx

! ... apply boundary conditions (delete the extra B-splines)

      ii=0
      Do j=1,mm
       if(j.le.ns) then; if(iprm(j,i).eq.0) Cycle; end if; ii=ii+1
       k=0
       Do jp=1,mm
        if(jp.le.ns) then; if(iprm(jp,i).eq.0) Cycle; end if
        k=k+1; xx(k)=aa(jp,j)
       End do
       ax(1:k,ii)=xx(1:k); yy(ii)=res(j)
      End do 

! ... solve the matrx equation with LAPACK routines:

      CALL DGETRF(ii, ii, ax, md, ipiv, info)
      if (info /= 0) Stop 'hf_nr: error in factorization routine DSYTRF'
      CALL DGETRS('N', ii, 1, ax, md, ipiv, yy, md, info)
      if (info /= 0) Stop 'hf_nr: error in solve routine DSYTRS'

! ... restore the solution in orginal basis:

      res=0.d0; k=0
      Do j=1,ns
       if(iprm(j,i).eq.0) Cycle
       k=k+1; res(j)=yy(k)
      End do 
       
      v = v - res(1:ns)
      y = matmul(BB,v)
      v = v/sqrt(Dot_Product(v,y))

! ... if (v(ks) < 0.d0) v = -v

      ipos=maxloc(abs(v))
      if(v(ipos(1)).lt.0.d0) v=-v

      e(i) = (eii - res(ns+1))  

      End Subroutine  Solve_nr
