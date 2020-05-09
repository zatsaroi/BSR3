!======================================================================
      Subroutine allocate_matrix(m)
!======================================================================
!     allocate arrays connected with the Hamiltonian matrix
!----------------------------------------------------------------------
      Use bsr_recoup

      Implicit none
      Integer :: ich,jch, i,k,m, mycase

      m = 0
      if(allocated(hcc)) Deallocate(hcc) 
      if(allocated(icc)) Deallocate(icc) 
      if(allocated(x))   Deallocate(x) 
      if(allocated(hcb)) Deallocate(hcb,hbb,icb,ibb) 

      if(nch.eq.0.or.ns.eq.0) Return

      Allocate(icc(nch,nch));  icc = 0; iicc = 0
      if(allocated(my_channel)) Deallocate(my_channel)
      Allocate(my_channel(nch+npert));  my_channel = 0
      m=m+nch*nch

      i=0; k=0
      Do ich=1,nch 
      Do jch=1,ich 
       if(nprocs.gt.1) then
        k=k+1; if(k.gt.nprocs-1) k=1; if(myid.ne.k) Cycle 
       end if
       i=i+1; icc(ich,jch) = i; icc(jch,ich) = i
       my_channel(ich) = 1
       my_channel(jch) = 1
      End do; End do
      iicc=i
      ich = (nch+1)*nch/2
      Allocate(hcc(ns,ns,iicc),x(ns,ns))
      m = m + 2*(ns*ns*iicc)
      hcc=0.d0

!-----------------------------------------------------------------------
      if(npert.gt.0) then  

      Allocate(icb(nch,npert));  icb = 0
      m = m + 2 * nch * nch
      i=0; k=0
      Do ich=1,nch
      Do jch=1,npert
       if(nprocs.gt.1) then
        k=k+1; if(k.gt.nprocs-1) k=1; if(myid.ne.k) Cycle 
       end if
       i=i+1; icb(ich,jch) = i
       my_channel(ich) = 1
       my_channel(nch+jch) = 1
      End do; End do
      iicb = i
      Allocate(hcb(ns,iicb));  hcb=0.d0
      m = m + 2 * ns * iicb

      Allocate(ibb(npert,npert));  ibb = 0
      m = m +  npert * npert
      i=0; k=0
      Do ich=1,npert; if(ich+nch.lt.I1_channel.or.ich+nch.gt.I2_channel) Cycle
      Do jch=1,ich; if(jch+nch.lt.J1_channel.or.jch+nch.gt.J2_channel) Cycle
       if(nprocs.gt.1) then
        k=k+1; if(k.gt.nprocs-1) k=1; if(myid.ne.k) Cycle 
       end if
       i=i+1; ibb(ich,jch) = i; ibb(jch,ich) = i
       my_channel(nch+ich) = 1
       my_channel(nch+jch) = 1
      End do; End do
      iibb = i
      Allocate(hbb(iibb));  hbb=0.d0
      m = m + iibb

      end if   ! over npert > 0
!-----------------------------------------------------------------------

      if(allocated(imycase)) Deallocate(imycase)
      Allocate( imycase(1:nch+npert,1:nch+npert) )
      m = m + (nch+npert)*(nch+npert)
      imycase = 0

      Do ich=1,nch+npert; Do jch=1,nch+npert
       imycase(ich,jch) = mycase(ich,jch)
      End do; End do

      End Subroutine allocate_matrix



!======================================================================
      Integer Function mycase(ich,jch)
!======================================================================
      Use bsr_recoup

      Integer, intent(in) :: ich,jch

      if(ich.le.nch.and.jch.le.nch) then
       mycase = icc(ich,jch)
      elseif(ich.le.nch.and.jch.gt.nch) then
       mycase = icb(ich,jch-nch)
      elseif(ich.gt.nch.and.jch.le.nch) then
       mycase = icb(jch,ich-nch)
      elseif(ich.gt.nch.and.jch.gt.nch) then
       mycase = ibb(ich-nch,jch-nch)
      else
       mycase = 0
      end if

      End Function mycase
