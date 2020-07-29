!======================================================================
      Subroutine Merge_CDATA(nn, ip, jp, nc, eps_cc)
!======================================================================
!     merge data in different blocks in MODULE 'cmdata'
!
!     nn     - number of blocks
!     ip(.)  - pointer for begining of block .
!     jp(.)  - pointer for end of block .
!     nc     - number of result coeff's
!     eps_cc - all coefficients < EPS_c are ignored
!----------------------------------------------------------------------
	
      USE bsr_mat, only: pri,debug
      USE cmdata 

      Implicit none

      Integer :: i,ii, j,jj, m,mm, nn,nc
      Integer :: ip(nn),jp(nn)
      Real(8) :: eps_cc, t1,t2

      Call CPU_time(t1)
      nc = 0

! ... choose the non-empty block

       mm = 0
       Do m = 1,nn
        if(JP(m).le.0) IP(m)=0
        i = IP(m); if(i.eq.0) Cycle
        mm = m; Exit
       End do
       if(mm.eq.0) Return

! ...  main loop ...

    1 Continue
                             
! ...  compare integrals in different blocks and merge the coefficients
! ...  in case of equal integrals

       Do ii=1,nn-1
        i=IP(ii); if(i.eq.0) Cycle
        Do jj=ii+1,nn
         j=IP(jj); if(j.eq.0) Cycle
         if(K1(i).ne.K1(j)) Cycle
         if(K2(i).ne.K2(j)) Cycle
         if(K3(i).ne.K3(j)) Cycle
         if(K4(i).ne.K4(j)) Cycle
         CDATA(i) = CDATA(i) + CDATA(j); CDATA(j) = 0.d0
         mm=jj; go to 2
        End do
       End do

! ...  choose the minimum K1, then K2, then K3, then K4 

       j=IP(mm)
       Do m=1,nn
        if(IP(m).eq.0) Cycle; i=IP(m)
        if    (K1(i).lt.K1(j)) then;   mm=m; j=IP(mm)
        elseif(K1(i).gt.K1(j)) then;   Cycle
        else
         if    (K2(i).lt.K2(j)) then;  mm=m; j=IP(mm)
         elseif(K2(i).gt.K2(j)) then;  Cycle
         else
          if    (K3(i).lt.K3(j)) then; mm=m; j=IP(mm)
          elseif(K3(i).gt.K3(j)) then; Cycle
          elseif(K4(i).lt.K4(j)) then; mm=m; j=IP(mm)
          end if
         end if
        end if
       End do

! ...  mark the chosen coefficient 

       i=IP(mm)
       if(abs(CDATA(i)).gt.EPS_cc) then; nc=nc+1; IPT(nc)=i; end if

! ...  choose next data

    2  IP(mm) = IP(mm) + 1
       if(IP(mm).le.JP(mm)) go to 1
       if(IP(mm).gt.JP(mm)) then
        IP(mm)=0
        Do m=1,nn
         if(IP(m).gt.0) then; mm=m; go to 1; end if
        End do
       end if
       Call CPU_time(t2); Tmerge = Tmerge + (t2-t1)
       if(debug.gt.1) &
       write(pri,'(a,2i9,f10.1)') 'merge data:', nn,nc,(t2-t1)/60

      End Subroutine Merge_CDATA


