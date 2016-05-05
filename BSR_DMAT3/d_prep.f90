!======================================================================
      Subroutine D_prep
!======================================================================
!     prepare B-spline multipole arrays and matrix elements
!----------------------------------------------------------------------
      Use bsr_dmat
      Use spline_param; Use spline_grid; Use spline_orbitals

      Implicit none
      Integer :: i,j,k,m
      Real(8) :: C,fl,ygr(nv,ks)
      Integer, external :: ITRA

! ... prepare B-splines arrays in the band non-symmetric storage mode:      

      Allocate(bbbs(ns,ks+ks-1),rrbs(ns,ks+ks-1), &
               ddbs(ns,ks+ks-1),ttbs(ns,ks+ks-1))

      Call ZINTYN (ns,ks,nv,bsp,bsp,grw,bbbs)

      k = kpol; if(ktype.eq.'M') k = kpol-1
      Do i=1,nv; Do m=1,ks
       ygr(i,m) = grw(i,m)*gr(i,m)**k
      End do;    End do
      Call ZINTYN (ns,ks,nv,bsp,bsp,ygr,rrbs)

      ddbs = 0.d0;  ttbs = 0.d0

      if(ktype.eq.'E'.and.kpol.eq.1)  then
       Call ZINTYN (ns,ks,nv,bsp,bspd,grw,ddbs)
       Do i=1,nv; Do m=1,ks
        ygr(i,m) = grw(i,m)*grm(i,m)
       End do;    End do 
      Call ZINTYN (ns,ks,nv,bsp,bsp,ygr,ttbs)
      end if 
       
      if(ktype.eq.'E'.and.kpol.eq.2)  then
       Do i=1,nv; Do m=1,ks
        ygr(i,m) = grw(i,m)*gr(i,m)
       End do;    End do 
       Call ZINTYN (ns,ks,nv,bsp,bspd,ygr,ddbs)
       ttbs = bbbs
      end if 

! ...  < . | p >,  < . | r | p >,  < . | d/dr | p > ,  < . | 1/r | p >

! ...     qbs           rbs           dbsr,dbsl             tbs

      Allocate(rbs(ns,nbf),dbsr(ns,nbf),dbsl(ns,nbf),tbs(ns,nbf))
      qbs = 0.d0; rbs = 0.d0; dbsr = 0.d0; dbsl = 0.d0; tbs = 0.d0

      Do i=1,nbf
       if(iech(i).ne.0) Cycle
       Call bav(ns,ks,bbbs,pbs(1,i),qbs(1,i) ,'n','r')
       Call bav(ns,ks,rrbs,pbs(1,i),rbs(1,i) ,'n','r')
       Call bav(ns,ks,ddbs,pbs(1,i),dbsr(1,i),'n','r')   
       Call bav(ns,ks,ddbs,pbs(1,i),dbsl(1,i),'n','l')   
       Call bav(ns,ks,ttbs,pbs(1,i),tbs(1,i) ,'n','r')
      End do

! ... one-electron overlaps 

      obs = 0.d0
      Do i=1,nbf; if(iech(i).ne.0) Cycle
       Do j=i,nbf; if(iech(j).ne.0) Cycle
        if(lbs(i).ne.lbs(j)) Cycle
        C = SUM(pbs(:,i)*qbs(:,j))
        if(abs(C).lt.Eps_ovl) C=0.d0
        obs(i,j)=C; obs(j,i)=C
       End do
      End do

! ... one-electron multipole L-integrals:

      Allocate(dipL(nbf,nbf)); dipL = 0.d0

      Do i=1,nbf;  if(iech(i).ne.0) Cycle
       Do j=i,nbf; if(iech(j).ne.0) Cycle
        if(ITRA(lbs(i),k,lbs(j)).eq.0) Cycle
        C = SUM(pbs(:,i)*rbs(:,j))
        dipL(i,j)=C; dipL(j,i)=C
       End do
      End do

! ... one-electron multipole V-integrals:

      Allocate(dipV(nbf,nbf)); dipV = 0.d0

      if(ktype.eq.'M')  dipV = dipL

      if(ktype.eq.'E') then
       Do i=1,nbf;  if(iech(i).ne.0) Cycle
        Do j=1,nbf; if(iech(j).ne.0) Cycle
         if(ITRA(lbs(i),kpol,lbs(j)).eq.0) Cycle
         if(mod(lbs(i)+lbs(j)+kpol,2).ne.0) Cycle
         Call FL_kpol(kpol,lbs(i),lbs(j),fl)
         dipV(i,j) = SUM(pbs(:,i)*dbsr(:,j)) + &
                     fl*SUM(pbs(:,i)*tbs(:,j))
         End do
        End do
       end if

      End Subroutine D_prep
