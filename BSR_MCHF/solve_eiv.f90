!==================================================================
      Subroutine Solve_eiv(i,hfm,v,rhs)
!==================================================================
!     Find the eigenvector of hfm for the m'th eigenvalue after 
!     orthogonality has been applied
!------------------------------------------------------------------
      Use bsr_mchf

      Implicit none
      Integer, intent(in)    :: i
      Real(8), intent(inout) :: hfm(ns,ns), rhs(ns)
      Real(8), intent(out)   :: v(ns)

      Real(8) :: aa(ns,ns), ss(ns,ns), w(3*ns)
      Real(8) :: eval(ns), a(ns), s(ns), rh(ns), eii,C
      Integer :: j, jp, info, k,ii,m, ipos(1)
    
      if(srhs.gt.0.d0.and.method.eq.1) then
       Call solve_direct (i,hfm,v,rhs); Return
      end if

      if(debug.eq.1) write(log,'(a,f10.5)') 'method - solve_eiv, srhs =',srhs

      eii = e(i)

! ... apply orthogonality conditions for orbitals through projection: 

      m = nbs(i)-lbs(i)
      Do j = 1, nbf; if(i.eq.j) Cycle  ! i-1  ???               
       if(i.le.ncore.and.j.ge.i) Cycle
       if(lbs(i).ne.lbs(j)) Cycle
       Call apply_orthogonality(hfm,p(1,j))
       if(m.gt.1) m = m - 1   
      End do

      if(iphys(i).eq.0) m=1

! ... apply boundary conditions (delete extra B-splines):

      ii=0
      Do j=1,ns
       if(iprm(j,i).eq.0) Cycle; ii=ii+1
       k=0
       Do jp=1,ns
        if(iprm(jp,i).eq.0) Cycle
        k=k+1; a(k)=hfm(jp,j); s(k)=BB(jp,j)  
       End do
       aa(1:k,ii)=a(1:k); ss(1:k,ii)=s(1:k); rh(ii) = rhs(j)
      End do 

! ... evaluates the eigenvalues and eigenvectors (LAPACK routine):

      Call dsygv(1,'V','L',ii,aa,ns,ss,ns,eval,w,3*ns,INFO)
      if (info /= 0) then
       WRITE(scr,'(a,i6)') 'scf_eiv: error in eigenvalue routine, dsygv', info
       Stop 
      end if

      if(debug.gt.0) write(log,'(a,5E15.5)') 'eval =',eval(1:5)

!----------------------------------------------------------------------
      if (srhs.eq.0.d0) then    ! without rhs 

      ! restore the solutions in original B-spline net:

      a(1:ns) = aa(1:ns,m);  v=0.d0; k=0
      Do j=1,ns
       if(iprm(j,i).eq.0) Cycle; k=k+1; v(j)=a(k)
      End do 

      ipos=maxloc(abs(v))
      if(v(ipos(1)).lt.0.d0) v=-v

      e(i) = eval(m)

!-----------------------------------------------------------------------
! ... use inverse iteration in case of rhs: 

      else   !  with rhs

      if(debug.eq.1) write(log,'(a)') 'method - inverse iterations'

       a = 0.d0
       Do j = 1,ii
        if(abs(eval(j)).lt.0.00001) Cycle                       !  ???
        C = Dot_Product(rh(1:ii),aa(1:ii,j))/(eval(j) - eii)
        a(1:ii) = a(1:ii) + aa(1:ii,j)*C
       End do

! ... restore the solutions in original B-spline net and normalize:

       v = 0.d0; k=0
       Do j=1,ns
        if(iprm(j,i).eq.0) Cycle; k=k+1; v(j)=a(k)
       End do 

       a = Matmul(BB,v)
       v = v/sqrt(Dot_Product(v,a))

! ... if (v(ks) < 0.d0) v = -v

       ipos=maxloc(abs(v))
       if(v(ipos(1)).lt.0.d0) v=-v

      end if
      
      End Subroutine  Solve_eiv
