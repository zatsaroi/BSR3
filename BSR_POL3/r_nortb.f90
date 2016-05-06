!======================================================================
      Subroutine Read_nortb
!======================================================================
! ... provide oscillator strengths for given dipole matrix
!----------------------------------------------------------------------
      Use bsr_pol
      Use spline_param, only: ns
      Use channel,      only: nch,ncp

      Implicit none
      Integer :: i,j, k,n, ip, nbound,nhmb,nchb,ncpb,nsb
      Real(8) :: S, aa(nhm),bb(nhm),cc(nhm)
      Integer, external :: Ipointer
       
      if(nortb.eq.0) Return

! ... check bound states file:

       i=INDEX(AF_bnd,'.'); AF=AF_bnd(1:i)//ALSP
       Call Check_file(AF); Open(nub,file=AF); rewind(nub)

       read(nub,'(5i10,3i5)') nsb,nchb,ncpb,nhmb,nbound

       if(nsb .ne.ns ) Stop 'bsr_pol: nsb <> ns '
       if(nchb.ne.nch) Stop 'bsr_pol: nchb --> ?'
       if(ncpb.ne.ncp) Stop 'bsr_pol: ncpb --> ?'
       if(nhmb.ne.nhm) Stop 'bsr_pol: nhmb --> ?'

       if(nbound.lt.nortb) Stop 'bsr_pol: nbound --> ?'

       Do i=1,nbound
        read(nub,*) j
        read(nub,*) S  
        read(nub,'(5D15.8)') cc
        ip = Ipointer(nortb,iortb,j)
        if(ip.eq.0) Cycle
        Do j = 1,nhm
         aa(1:nhm) = c(1:nhm,j)
         bb(j) = SUM(aa(1:nhm)*cc(1:nhm))
        End do 
        write(nuq) bb
        iortb(ip) = -iortb(ip) 
        k=0; Do n=1,nortb; if(iortb(n).lt.0) Cycle; k=1; Exit; End do
        if(k.eq.0) Exit
       End do

       k=0; Do n=1,nortb; if(iortb(n).lt.0) Cycle; k=1; Exit; End do

       if(k.gt.0) Stop 'Add_nortb: not all states found'
       iortb = - iortb

       End Subroutine read_nortb



