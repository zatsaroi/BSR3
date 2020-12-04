!======================================================================
      Subroutine write_nl
!======================================================================
!     record B-spline solutions for given nl-series
!----------------------------------------------------------------------
      Use bsr_hf
      Use hf_orbitals

      Implicit none
      Integer :: i,j,m,k,kk, n,l,iset,ip,jp, start,info, nsol_nl
      Integer :: jord(nwf)
      Real(8) :: hfm(ns,ns),aa(ns,ns),ss(ns,ns), &
                 w(3*ns),eval(ns),a(ns),s(ns),v(ns)
      Character(200) :: string = ' '
      Character(4) :: EL
      Integer, external :: Ifind_orb

      Call Read_string(inp,'nl',string)     
      Call Read_aarg('nl',string)     

      if(len_trim(string).eq.0) Return

      jord = 0
      if(string(1:3)=='ALL' .or. string(1:3)=='all') then 
       Do i=1,nwf; jord(i)=i; End do
      else
       start = 1; ip = 0
       Do  
        i = index(string(start:),',')
        if (i /= 0 .or. LEN_TRIM(string(start:)) /= 0) then
         read(string(start:),*) EL
         Call EL4_NLK(EL,n,l,iset)
         j = Ifind_orb(n,l,iset)
         if(j.gt.0) then; ip=ip+1; jord(j)=j; end if
        end if
        start = start + i 
        if(i == 0 .or.  LEN_TRIM(string(start:)) == 0) Exit
       End do
       if(ip.eq.0) Return
      end if

      Do i = 1,nwf; if(jord(i).eq.0) Cycle 

      Call hf_matrix(i,hfm)

! ... apply orthogonality conditions 

      m = nbs(i)-lbs(i)
      Do j = 1,nwf 
       if(i.eq.j) Cycle 
       if(lbs(i).ne.lbs(j)) Cycle
       if(e(i,j) < 1.d-10) Cycle     
       Call orthogonality(hfm,p(1,j))
       if(j.lt.i) m = m - 1
      End do

! ... apply boundary conditions (delete extra B-splines)

      kk=0
      Do j=1,ns
       if(iprm(j,i).eq.0) Cycle; kk=kk+1
       k=0
       Do jp=1,ns
        if(iprm(jp,i).eq.0) Cycle
        k=k+1; a(k)=hfm(jp,j); s(k)=bb(jp,j)  
       End do
       aa(1:k,kk)=a(1:k); ss(1:k,kk)=s(1:k)
      End do 

! ... evaluates the eigenvalues and eigenvectors (LAPACK routine):

      Call dsygv(1,'V','L',kk,aa,ns,ss,ns,eval,w,3*ns,INFO)
      if (info /= 0) then
       write(scr,'(a,i6)') 'Error in Eigenvalue routine, dsygv', info
       Stop ' '
      end if

! ... save all solutions:

      nsol_nl = kk  
      AF_nl = trim(name)//'_'//ebs(i)//BF_nl
      Call Clean_a(AF_nl)
      open(nuw,file=AF_nl,form='UNFORMATTED')
      write(nuw) ns,ks
      write(nuw) t(1:ns+ks)
      write(nuw) nsol_nl
      Do m=1,nsol_nl
       a(1:ns) = aa(1:ns,m);  v=0.d0; k=0
       Do j=1,ns
        if(iprm(j,i).eq.0) Cycle; k=k+1; v(j)=a(k)
       End do 
       if (v(ks) < 0.d0) v = -v
       write(nuw) eval(m)
       write(nuw) v
      End do
      Close(nuw)

      End do ! over i

      End Subroutine write_nl
