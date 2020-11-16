!======================================================================
      Subroutine get_estimates
!======================================================================
!     Get initial estimates:
!     (1) read from bsw.inp
!     (2) screened hydrogenic
!----------------------------------------------------------------------
      Use bsr_mchf
      
      Implicit none
      Integer :: i, j
      Real(8) :: s, snl(nbf)
      Real(8), external :: QUADR

      write(log,'(/a/)') 'Initial estimations for orbitals:'

! ... Read AF_inp for the fixed orbitals or initial estimates

      AF_inp=trim(name)//'.bsw'
      Call Read_apar(inp,'inp',AF_inp)
      Call Read_aarg('inp',AF_inp)
      if(Icheck_file(AF_inp).ne.0) then
       Open(nuw,file=AF_inp,form='UNFORMATTED')
       i = INDEX(AF_inp,'.',BACK=.TRUE.)
       if(AF_inp(i:).eq.'.bsw')  Call Read_bwfn (nuw)
      end if
      snl = 0.d0
      Do i = 1, nbf
       if(mbs(i).ne.0) Cycle 
       Call Screen_nl(i,snl(i))
      End do 
                                                      
      Do i = 1, nbf

! ...  We have an initial estimate

       if(mbs(i) /= 0) then
        write(log,'(a,a,a)') ebs(i), ' - read from file ', AF_inp
        Cycle   
       end if

       Call bhwf(nbs(i),lbs(i),z-snl(i),p(1,i)) 

       Call Check_tails(i)
  
! ... set orthogonality constraints:

       Do j = 1,i-1
        if (lbs(i) /= lbs(j)) Cycle
        S=QUADR(i,j,0)
        if(abs(S).lt.Eps_ovl) Cycle
        p(:,i) = p(:,i) - s * p(:,j)
       End do
       S=QUADR(i,i,0)
       p(:,i) = p(:,i)/ sqrt(S)

       write(log,'(a,a,f5.2)') ebs(i),&
       ' - hydrogenic orbital with screening ',snl(i)

      End do  ! over orbitals

      End Subroutine get_estimates


!======================================================================
      Subroutine Read_bwfn (nu)
!======================================================================
!     read radial orbitals in B-spline representation from file 'nu'
!----------------------------------------------------------------------
      Use bsr_mchf
      
      Implicit none
      Integer, intent(in) :: nu
      Character(4) :: el
      Real(8) :: zw,hw,hmw,rmw
      Integer :: ksw,nsw,mw,n,l,k,i
      Integer, external :: Ifind_bsorb

      rewind(nu)
    1 read(nu,end=2) el,zw,hw,hmw,rmw,ksw,nsw,mw
      if(zw.ne.z) Stop ' R_bwfn:  z <> zw'
      if(abs(hw-h).gt.1.d-12) Stop ' R_bwfn:  h <> hw'
      if(abs(hmw-hmax).gt.1.d-12) Stop ' R_bwfn:  hmw <> hmax'
      if(abs(rmw-rmax).gt.1.d-12) Stop ' R_bwfn:  rmw <> rmax'
      if(ksw.ne.ks) Stop ' R_bwfn:  ksw <> ks'
      if(nsw.ne.ns) Stop ' R_bwfn:  nsw <> ns'
      Call EL4_nlk(el,n,l,k)
      i = Ifind_bsorb(n,l,k)
      if(i.eq.0) go to 2
      read(nu) p(1:mw,i)
      mbs(i) = mw
      if(mw.lt.ns) p(mw+1:ns,i) = 0.d0 
      go to 1
    2 Continue

      End Subroutine Read_bwfn


!======================================================================
      Subroutine Screen_nl(io,S)
!======================================================================
! ... find screenning parameter for orbital i
!----------------------------------------------------------------------
      Use bsr_mchf

      Implicit none
      Integer, intent(in) :: io
      Real(8), intent(out) :: S
      Integer :: i, ic,jc, ip
      Real(8) :: C

! ... core orbitals:

      if(io.le.ncore) then
       S = 0.d0
       Do i = 1, io-1
        S = S + qsum(i)
       End do
       S = S + qsum(io)/2      
       Return
      end if

! ... find leading configuration for this orbital:

      S = 0.d0; jc = 0
      Do ic = 1,ncfg; C = WC(ic)*WC(ic); if(C.le.S) Cycle
       Call Get_cfg_LS(ic)     
       ip = ip_state(ic)
       Do i=1,no; ip=ip+1
        if(IP_orb(ip).ne.io) Cycle
        ip = -1; Exit  
       End do
       if(ip.gt.0) Cycle;  S = C; jc = ic
      End do

! ... find screening papameter from given configuration:

      S = 0.d0
      if(ncore.gt.0) S = SUM(qsum(1:ncore))
      Call Get_cfg_ls (jc)     
      ip = ip_state(jc)
      Do i=1,no; ip=ip+1
       if(IP_orb(ip).ne.io) then
        S = S + iq(i)
       else
        S = S + iq(i)/2 
        Exit
       end if
      End do

      End Subroutine Screen_nl


