!======================================================================
      Subroutine get_estimates
!======================================================================
!     Get initial estimates:
!     (1) read from bsw.inp
!     (2) screened hydrogenic
!----------------------------------------------------------------------
      Use bsr_hf
      Use hf_orbitals
      
      Implicit none
      Real(8) :: zz,ss,s
      Real(8), external :: QUADR_hf
      Integer :: i, j, nnl
      Integer, allocatable :: nl(:)
      Real(8), allocatable :: qnl(:), snl(:)
      Integer, external :: ipointer

      Call Alloc_hf_radial(ns)

      write(log,'(//a/)') 'Initial estimations for orbitals:'

! ... Read AF_inp for the fixed orbitals or initial estimates

      AF_inp=trim(name)//BF_inp
      Call Read_apar(inp,'inp',AF_inp)
      Call Read_aarg('inp',AF_inp)
      if(Icheck_file(AF_inp).ne.0) then
       open(nuw,file=AF_inp,form='UNFORMATTED')
       i = INDEX(AF_inp,'.',BACK=.TRUE.)
       if(AF_inp(i:).eq.'.bsw')  Call Read_bwfn (nuw)
      end if
      
      Allocate(nl(nwf),qnl(nwf),snl(nwf))

      nl=0; nnl = 0
      Do i = 1, nbf
       j = ipointer(nwf,nl,nbs(i)*1000+lbs(i))
       if(j.gt.0) then
        qnl(j) = qnl(j) + qsum(i)
       else
        nnl=nnl+1; nl(nnl) =  nbs(i)*1000+lbs(i)
                   qnl(nnl) = qsum(i)
       end if
      End do

      ss = 0.d0
      Do i = 1,nnl
       s =  max(0.d0,qnl(i)-1.d0)
       snl(i) = ss + s         
       ss = ss + qnl(i)
      End do
                                                      

      Do i = 1, nbf
       j = ipointer(nnl,nl,nbs(i)*1000+lbs(i))

       ss = snl(j); zz = z-ss

! ...  We have an initial estimate

       if(mbs(i) /= 0) then
        write(log,'(a,a,a)') ebs(i), ' - read from file ', AF_inp
        Cycle   
       end if

       Call bhwf(nbs(i),lbs(i),zz,p(1,i)) 
       mbs(i) = ns - 1
       Do j=ns-1,ns/3,-1
        if(abs(p(j,i)).gt.end_tol/1000) Exit   ! ???
        mbs(i) = j - 1
       End do
  
! ... set orthogonality constraints:

       Do j = 1,i-1
        if (lbs(i) /= lbs(j)) Cycle
        S=QUADR_hf(i,j,0)
        if(abs(S).lt.Eps_ovl) Cycle
        p(:,i) = p(:,i) - s * p(:,j)
       End do
       S=QUADR_hf(i,i,0)
       p(:,i) = p(:,i)/ sqrt(S)

       write(log,'(a,a,f5.2)') ebs(i),&
       ' - hydrogenic orbital with screening ',ss

      End do  ! over orbitals

      Call Check_tails(0)
 
      Call update_int(0)
      Call Energy   

      End Subroutine get_estimates


!======================================================================
      Subroutine Read_bwfn (nu)
!======================================================================
!     read radial orbitals in B-spline representation from file 'nu'
!----------------------------------------------------------------------
      Use bsr_hf
      Use hf_orbitals
      
      Implicit none
      Integer, intent(in) :: nu
      Character(4) :: el
      Real(8) :: zw,hw,hmw,rmw
      Integer :: ksw,nsw,mw,n,l,k,i
      Integer, external :: Ifind_orb

      rewind(nu)
    1 read(nu,end=2) el,zw,hw,hmw,rmw,ksw,nsw,mw
      if(zw.ne.z) Stop ' R_bwfn:  z <> zw'
      if(abs(hw-h).gt.1.d-12) Stop ' R_bwfn:  h <> hw'
      if(abs(hmw-hmax).gt.1.d-12) Stop ' R_bwfn:  hmw <> hmax'
      if(abs(rmw-rmax).gt.1.d-12) Stop ' R_bwfn:  rmw <> rmax'
      if(ksw.ne.ks) Stop ' R_bwfn:  ksw <> ks'
      if(nsw.ne.ns) Stop ' R_bwfn:  nsw <> ns'
      Call EL4_nlk(el,n,l,k)
      i = Ifind_orb(n,l,k)
      if(i.eq.0) go to 2
      read(nu) p(1:mw,i)
      mbs(i) = mw
      if(mw.lt.ns) p(mw+1:ns,i) = 0.d0 
      go to 1
    2 Continue

      End Subroutine Read_bwfn



