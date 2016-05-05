!======================================================================
      Subroutine Merge_CDATA(nn, ip, jp, nc, EPS_c)
!======================================================================
!     merge the different blocks of data in MODULE 'cmdata'
!
!     nn    - number of blocks
!     ip(.) - pointer for begining of block .
!     jp(.) - pointer for end of block .
!     nc    - number of result coeff's
!     EPS_c - all coefficients < EPS_c are ignored
!----------------------------------------------------------------------
      Use cmdata, only: CLDATA,CVDATA, K1,K2,K3, IPT

      Implicit none
      Integer :: nn,nc
      Integer :: ip(nn),jp(nn)
      Integer :: i,ii, j,jj, m,mm
      Real(8), intent(in) :: EPS_C

      nc = 0

! ... choose the non-empty block

       mm = 0
       Do m = 1,nn
        i = IP(m); if(i.eq.0) Cycle
        mm = m; Exit
       End do
       if(mm.eq.0) Return

! ...  main loop ...

    1 Continue
                             
! ... compare integrals in different blocks and merge the coefficients
! ... in case of equal integrals

      Do ii=1,nn-1
       i=IP(ii); if(i.eq.0) Cycle
       Do jj=ii+1,nn
        j=IP(jj); if(j.eq.0) Cycle
        if(K1(i).ne.K1(j)) Cycle
        if(K2(i).ne.K2(j)) Cycle
        if(K3(i).ne.K3(j)) Cycle
        CLDATA(i) = CLDATA(i) + CLDATA(j); CLDATA(j) = 0.d0
        CVDATA(i) = CVDATA(i) + CVDATA(j); CVDATA(j) = 0.d0
        mm=jj; go to 2
       End do
      End do

! ...  choose the minimum K1, then K2, then K3 

      j=IP(mm)
      Do m=1,nn
       if(IP(m).eq.0) Cycle; i=IP(m)
       if(K1(i).lt.K1(j)) then
        mm=m; j=IP(mm)
       elseif(K1(i).gt.K1(j)) then
        Cycle
       elseif(K2(i).lt.K2(j)) then
        mm=m; j=IP(mm)
       elseif(K2(i).eq.K2(j).and.K3(i).lt.K3(j)) then
        mm=m; j=IP(mm)
       end if
      End do

! ... mark the chosen coefficient 

      i=IP(mm)
      if(abs(CLDATA(i))+abs(CVDATA(i)).gt.EPS_c) then
	      nc=nc+1; IPT(nc)=i
	     end if

! ... choose next data

    2 IP(mm) = IP(mm) + 1
      if(IP(mm).le.JP(mm)) go to 1
      if(IP(mm).gt.JP(mm)) then
       IP(mm)=0
       Do m=1,nn
        if(IP(m).gt.0) then; mm=m; go to 1; end if
       End do
      end if

      End Subroutine Merge_CDATA


