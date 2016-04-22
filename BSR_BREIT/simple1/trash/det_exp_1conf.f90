!======================================================================
      Subroutine Det_expn_1conf  
!======================================================================
!     determined all possible determinants and their coefficients
!     for configuration given in module conf_LS_1 
!     (module also contains values for ne, MLT, MST)
!---------------------------------------------------------------------
      Use conf_LS_1

      Implicit none
      Integer :: in(ne),md(ne),nd(ne),idet(ne),it(ne)
      Integer :: ML(ne),MS(ne),MLp(ne),MSp(ne)  
      Integer :: mdt,kd,i,j,k,m,ii
      Real(8) :: C
      Integer, external :: Iglq, ML_id, MS_id, Iterm_LS
      Real(8), external :: Clebsh, DETC_sh

! ... auxiliary arrays:
! ... in(i)  -  pointer to given shell
! ... md(i)  -  number of determinant for given shell
! ... it(i)  -  pointer to shell term
! ... mdt    -  max. possible number of determinants

      k=1; mdt=1 
      Do i=1,no
       in(i)=k; k=k+iq(i); md(i)=Iglq(ln(i),iq(i));  mdt=mdt*md(i) 
	      it(i)=Iterm_LS(ln(i),iq(i),0,LS(i,1),LS(i,2),LS(i,3)) 
      End do
      if(allocated(Cdet) ) Deallocate(Cdet );  Allocate(Cdet(mdt)    )
      if(allocated(MLdet)) Deallocate(MLdet);  Allocate(MLdet(ne,mdt))
      if(allocated(MSdet)) Deallocate(MSdet);  Allocate(MSdet(ne,mdt))

! ... exhasting all determinants:

      kdt=0; i=1; nd(i)=1              
    1 kd = kd + 1
      ii = in(i)
      Call DET_sh(ln(i),iq(i),nd(i),ML(i),MS(i),Idet(ii)) 

      m = iabs(ML(i)); if(ML(i).lt.0) m=m+2
      if(m.gt.LS(i,2)) go to 2
      m = iabs(MS(i)); if(MS(i).lt.0) m=m+2
      if(m.gt.LS(i,3)) go to 2

      if(i.eq.1) then
       MLp(1)=ML(1); MSp(1)=MS(1)
      else
       MLp(i) = MLp(i-1)+ML(i)-1
       MSp(i) = MSp(i-1)+MS(i)-1
       m = iabs(MLp(i)); if(MLp(i).lt.0) m = m + 2
       if(m.gt.LS(i,4)) go to 2
       m = iabs(MSp(i)); if(MSp(i).lt.0) m = m + 2
       if(m.gt.LS(i,5)) go to 2
      end if

      if(i.lt.no) then; i = i + 1; nd(i) = 1;  go to 1; end if

      if(MLp(no).ne.MLT) go to 2
      if(MSp(no).ne.MST) go to 2

! ... coefficient calculation:

      C = 1.d0
      Do j=1,no
       C=C*DETC_sh(ln(j),iq(j),it(j),nd(j))
       if(C.eq.0.d0) Exit
      End do
      Do j=2,no
        C=C*Clebsh(LS(j-1,4),MLp(j-1), &
                   LS(j  ,2),ML (j  ), &
                   LS(j  ,4),MLp(j  ))
        if(C.eq.0.d0) Exit
        C=C*Clebsh(LS(j-1,5),MSp(j-1), &
                   LS(j  ,3),MS (j  ), &
                   LS(j  ,5),MSp(j  ))
        if(C.eq.0.d0) Exit
      End do
      if(C.eq.0.d0) go to 2

      kdt=kdt+1
      Cdet(kdt) = C
      Do j = 1,ne;  
       MLdet(j,kdt)=(ML_id(idet(j))-1)/2
       MSdet(j,kdt)= MS_id(idet(j))
      End do

! ... selecting the next case

    2 nd(i)=nd(i)+1                
      if(nd(i).gt.md(i)) then
       if(i.eq.1) go to 3          
       i=i-1; go to 2
      end if
      go to 1
    3 Continue

      End Subroutine Det_expn_1conf

