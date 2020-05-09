!======================================================================
      Subroutine Record_matrix(nu,ii)
!======================================================================
!     record interaction (overlap) matrix to the unit 'nu'
!----------------------------------------------------------------------
      Use bsr_recoup

      Implicit none
      Integer, intent(in) :: nu,ii
      Integer :: i,ip,jp,ich,jch,idiag
      Real(8) :: S

! ... channel-channal blocks:

      idiag = 1
      Do ich=1,nch; Do jch=1,ich;  i=icc(ich,jch);  if(i.eq.0) Cycle
       S = SUM(abs(hcc(:,:,i)))
       if(S.eq.0.d0) Cycle
       write(nu) ich,jch,ii
       write(nu) hcc(:,:,i)
       if(ich.ne.jch) idiag = 0
      End do; End do

! ... pertubers:

      if(npert.gt.0) then

       Do ich=1,nch; Do ip=1,npert; i=icb(ich,ip); if(i.eq.0) Cycle
        S = SUM(abs(hcb(:,i)))
        if(S.eq.0.d0) Cycle
        write(nu) ip+nch,ich,ii
        write(nu) hcb(:,i)
        idiag = 0
       End do; End do
       
       Do ip=1,npert; Do jp=1,ip; i=ibb(ip,jp); if(i.eq.0) Cycle
        if(hbb(i).eq.0.d0) Cycle
        write(nu) ip+nch,jp+nch,ii
        write(nu) hbb(i)
        if(ip.ne.jp) idiag = 0
       End do; End do

      end if

! ... sign of the end

      write(nu) -1,-1,idiag

      End Subroutine Record_matrix

