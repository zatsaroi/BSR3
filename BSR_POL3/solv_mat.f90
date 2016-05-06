!======================================================================
      Subroutine Solv_mat
!======================================================================
!     solves the dipole equation: (H - E0 C) x = d
!----------------------------------------------------------------------
      Use bsr_pol
      Use channel,         only: nch,ncp,ipar,lpar,ispar,ELC,ipch,lch
      Use spline_param,    only: ns,ks
      Use spline_orbitals, only: nbf,qbs,ebs,IBORT,lbs

      Implicit none
      Integer :: iprm(mhm),ipiv(mhm)
      Integer :: i,j,m,k,l,info,ii,jj,ich
      Real(8) :: sol(mhm), aa(mhm), cc(mhm)
      Real(8) :: S,EP, alfa,fvalue,dmat

! ... save matrixes:

      Open(nua,form='UNFORMATTED',status='SCRATCH')
      rewind(nua); Do j=1,nhm; write(nua) a(1:nhm,j); End do
      write(nua) d(1:nhm)
      Do j=1,nhm; write(nua) c(1:nhm,j); End do
            
! ... orthogonal constraints:

      write(pri,'(/a/)') 'orthogonal constraints:' 

      jj = nhm
      Do ich = 1,nch; i=ipch(ich); ii=(ich-1)*ns
       Do j=1,nbf;
        if(IBORT(i,j).ne.0) Cycle 
        if(lbs(i).ne.lbs(j)) Cycle
        write(pri,'(2a6)') ebs(i),ebs(j)
        jj=jj+1
        a(ii+1:ii+ns,jj) = qbs(1:ns,j)
        a(jj,ii+1:ii+ns) = qbs(1:ns,j)
       End do
      End do
      j=jj-nhm; if(j.ne.nort) Stop 'nort <> applied'

      if(nortb.gt.0) then
      rewind(nuq)
      Do i=1,nortb
       read(nuq) aa(1:nhm)
       jj=jj+1
       a(1:nhm,jj) = aa(1:nhm)
       a(jj,1:nhm) = aa(1:nhm)
      End do
      end if

! ... apply zero conditions at r=0 and on boundary r=a:

      iprm=1
      Do i=1,nch
       j=(i-1)*ns
       l=lch(i); if(l.gt.ks-2) l=ks-2 
       if(ilzero.eq.0) l=0
       iprm(j+1:j+l+1)=0                   
       iprm(j+ns-ibzero+1:j+ns)=0  
      End do

! ... delete the extra B-splines:

      m=0
      Do j=1,mhm
       if(iprm(j).eq.0) Cycle; m=m+1
       k=0
       Do i=1,mhm
        if(iprm(i).eq.0) Cycle
        k=k+1; aa(k)=a(i,j); cc(k)=c(i,j)
       End do
       a(1:k,m)=aa(1:k); c(1:k,m)=cc(1:k);  d(m)=d(j)
      End do 
      khm=m 

      a = a - E0 * c

! ... solution of equation::

      write(pri,'(/a,i5,a)') 'khm  = ',khm,'  - final dimension'

      sol = d
      Call DGESV( khm, 1, a, mhm, ipiv, sol, mhm, INFO )

      if(info.ne.0) Stop 'BSR_POL: solution failed'

! ... restore the solutions in original B-spline net:

      d=sol; sol=0.d0; k=0
      Do i=1,nhm
       if(iprm(i).eq.0) Cycle; k=k+1; sol(i)=d(k)
      End do 

! ... restore the interaction, overlap and dipole matrixes:

      rewind(nua); Do j=1,nhm; read(nua) a(1:nhm,j); End do
      read(nua) d(1:nhm)
      Do j=1,nhm; read(nua) c(1:nhm,j); End do

! ... normalization of solution:

      S = 0.0
      Do i = 1,nhm
      Do j = 1,nhm
       S = S + sol(i)*c(i,j)*sol(j)
      End do
      End do
      S = sqrt(S);
      write(pri,'(/a,f10.5)') 'norma = ',S
      sol = sol / S

! ... energy:

      EP = 0.d0
      Do i = 1,nhm
      Do j = 1,nhm
       EP = EP + sol(i)*a(i,j)*sol(j)
      End do
      End do

      write(pri,'(/a,f16.8)') 'EP = ',EP

      write(pri,'(/a,f10.5)') 'EP-E0 = ',(EP-E0)*27.2113

! ... alpha:

      dmat   = Sum(sol(1:nhm)*d(1:nhm))
      fvalue = 2.d0*dmat*dmat/((kpol+kpol+1)*jot0)*(EP-E0)
      alfa   = 2.d0*dmat*dmat/((EP-E0)*(kpol+kpol+1)*jot0)

      write(*,*) 'alfa = ',alfa,'  kpol =',kpol

      write(pri,*)
      write(pri,'(a,f14.8)') 'dmat   = ',dmat
      write(pri,'(a,f14.8)') 'fvalue = ',fvalue
      write(pri,'(a,f14.8)') 'alfa   = ',alfa

! ... output of solution in BSR-bound format:

      i = INDEX(AF_pol,'.',BACK=.TRUE.); AF=AF_pol(1:i)//ALSP
      Open(nur,file=AF)

      write(nur,'(5i10,3i5)') ns,nch,ncp,nhm,1,lpar,ispar,ipar
      write(nur,'(i5,2x,a)') 1,ELC(1)
      write(nur,'(2F15.8,2F15.5)') EP,dmat,alfa,fvalue
      write(nur,'(5D15.8)') sol(1:nhm)

! ... check the orthpgonality:

      write(pri,'(/a/)') 'orthogonality:' 

      jj = nhm
      Do ich = 1,nch; i=ipch(ich); ii=(ich-1)*ns
       Do j=1,nbf; if(IBORT(i,j).ne.0) Cycle 
        if(lbs(i).ne.lbs(j)) Cycle
        S = SUM(sol(ii+1:ns)*qbs(1:ns,j))
        write(pri,'(2a6,f10.5)') ebs(i),ebs(j),S
       End do
      End do

      if(nortb.gt.0) then
       rewind(nuq)
       Do i=1,nortb
        read(nuq) aa(1:nhm)
        S = SUM(sol(1:nhm)*aa(1:nhm))
        write(pri,'(a,i2,f10.5)') '  nortb = ',i,S
       End do      
      end if

      End Subroutine Solv_mat


