!=====================================================================
      Subroutine Solve_direct (i,hfm,v,rhs)
!=====================================================================
!     Direct solution of the one-electron equation:  hfm v - e v = rhs 
!---------------------------------------------------------------------
      Use bsr_mchf
 
      Implicit none
      Integer :: i
      Real(8), intent(in)  :: hfm(ns,ns), rhs(ns)
      Real(8), intent(out) :: v(ns)

      Real(8), dimension(ns+nbf,ns+nbf) :: aa, ax
      Real(8), dimension(ns+nbf) :: res, xx, yy
      Integer, dimension(ns+nbf) :: ipiv

      Real(8) :: y(ns)
      Integer :: mm, md, k, info, j, jp, ii, ipos(1)

      if(debug.eq.1) write(log,'(a)') 'method - direct'

      md = ns+nbf
      mm = ns

      aa = 0.d0;  aa(1:ns,1:ns) = hfm - e(i)*BB
      res = 0.d0; res(1:ns) = rhs
 
! ... add orthogonality with all orbitals:

      Do j = 1,i-1                             !  ???
       if(lbs(i).ne.lbs(j)) Cycle
       mm = mm + 1
       y = matmul(BB,p(:,j))
       aa(1:ns,mm) = -y
       aa(mm,1:ns) = -y
      End do

! ... apply boundary conditions (delete the extra B-splines):

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
      if (info /= 0) Stop 'solve_direct: error in factorization routine DSYTRF'
      CALL DGETRS('N', ii, 1, ax, md, ipiv, yy, md, info)
      if (info /= 0) Stop 'solve_direct: error in solve routine DSYTRS'

! ... restore the solution in orginal basis:

      v=0.d0; k=0
      Do j=1,ns
       if(iprm(j,i).eq.0) Cycle; k=k+1; v(j)=yy(k)
      End do 

! ... normalize the solution:

      y = matmul(BB,v)
      v = v/sqrt(Dot_Product(v,y))

! ... solution sign:

      ipos=maxloc(abs(v))
      if(v(ipos(1)).lt.0.d0) v=-v

      End Subroutine  Solve_direct
