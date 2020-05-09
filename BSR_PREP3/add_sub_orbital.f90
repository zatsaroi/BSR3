!======================================================================
      Subroutine Add_sub_orbital(ii,jj)
!======================================================================
      Use bsr_prep

      Implicit real(8) (A-H,O-Z)

      m=nbf+1; if(m.ge.mbf) CALL Allocate_bsorb(mbf+jbf)

      nbs(m)=nbs(ii); lbs(m)=lbs(ii); kbs(m)=kbs(ii); ebs(m)=ebs(ii)
      mbs(m)=mbs(ii); p(:,m) = p(:,ii)

       OBS(:,m) = 0.d0
       Do i = 1,m
        if(lbs(i).ne.lbs(m)) Cycle
        OBS(i,m)=QUADR(i,m,0)
        OBS(m,i)=OBS(i,m)
       End do

! ... check if we need orthogonalization:

      jj = ii
      Do i=1,nbf; if(iech(i).ne.1) Cycle
       S = OBS(i,m); if(abs(S).lt.eps_ovl) Cycle           
       p(:,m) = p(:,m) - S * p(:,i);  jj=m
       mbs(m) = max(mbs(i),mbs(m))
      End do

      if(jj.eq.ii) then; iech(ii)=1; Return; End if

       S = QUADR(m,m,0);  S=sqrt(S);  p(:,m)=p(:,m)/S  

       OBS(:,m) = 0.d0
       Do i = 1,m
        if(lbs(i).ne.lbs(m)) Cycle
        OBS(i,m)=QUADR(i,m,0)
        OBS(m,i)=OBS(i,m)
       End do

! ... assign set index for new sub. orbital: 

      Call Assign_index(m); nbf=m; iech(m)=1; jj=m

      End Subroutine Add_sub_orbital
