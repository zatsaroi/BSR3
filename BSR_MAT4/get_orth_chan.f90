!======================================================================
      Subroutine Get_orth_chan
!======================================================================
!     load the channels orthogonality conditions arrays
!----------------------------------------------------------------------
      Use bsr_mat

      Implicit none
      Integer :: ich, io,jo, k
      Integer, external :: IBORT

      if(allocated(ip_orth_chan)) Deallocate(ip_orth_chan)
      if(allocated(ip_orth_orb )) Deallocate(ip_orth_orb)

      Allocate(ip_orth_chan(0:nch),ip_orth_orb(nch))

      k = 0; ip_orth_chan=0; korth = nch; ip_orth_orb=0
      Do ich = 1,nch; io = ipch(ich)
       Do jo = 1,nbf;  if(iech(jo).ne.0) Cycle
        if(lbs(io).ne.lbs(jo)) Cycle
        if(IBORT(io,jo).ne.0) Cycle
        k = k + 1
        if(k.gt.korth) then
         Allocate(jp_orth_orb(korth))
         jp_orth_orb(1:korth) = ip_orth_orb(1:korth)
         Deallocate(ip_orth_orb)
         Allocate(ip_orth_orb(korth+nch))
         ip_orth_orb(1:korth) = jp_orth_orb(1:korth)
         Deallocate(jp_orth_orb)
         korth = korth + nch
        end if
        ip_orth_orb(k) = jo
       End do
       ip_orth_chan(ich) = k
      End do

      End Subroutine Get_orth_chan
