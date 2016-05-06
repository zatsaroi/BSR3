!======================================================================
      Subroutine Record_matrix(nu)
!======================================================================
!     record interaction (overlap) matrix to the unit 'nu'
!----------------------------------------------------------------------
      Use bsr_mat
      Use bsr_matrix

      Implicit none
      Integer, intent(in) :: nu
      Integer :: ib,ip,jp,ich,jch,idiag
      Real(8) :: S

! ... diagonal blocks:

      Do ib = 1, kch;  write(nu) hcc(:,:,ib*(ib+1)/2); End do

      idiag = 1

! ... non-diagonal blocks:

      Do ich = 1,kch; Do jch = 1,ich; if(ich.eq.jch) Cycle
       ib = ich*(ich-1)/2+jch
       S = SUM(abs(hcc(:,:,ib)))
       if(S.eq.0.d0) Cycle
       write(nu) ich,jch
       write(nu) hcc(:,:,ib)
       idiag = 0
      End do; End do

! ... pertubers:

      if(kcp.gt.0) then

       Do ich = 1,kch; Do ip = 1,kcp
        S = SUM(abs(hcb(:,ich,ip)))
        if(S.eq.0.d0) Cycle
        write(nu) ip+kch,ich
        write(nu) hcb(:,ich,ip)
        idiag = 0
       End do; End do
       
       Do ip = 1,kcp; Do jp = 1,ip
        ib = ip*(ip-1)/2+jp
        if(hbb(ib).eq.0.d0) Cycle
        write(nu) ip+kch,jp+kch
        write(nu) hbb(ib)
        if(ip.ne.jp) idiag = 0
       End do; End do

      end if

! ... sign of the end

      write(nu) 0,idiag

      End Subroutine Record_matrix

