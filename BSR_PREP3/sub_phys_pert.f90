!======================================================================
      Subroutine SUB_phys_pert
!======================================================================
!     define physical orbitals in perturber
!----------------------------------------------------------------------
      Use bsr_prep

      Implicit real(8) (A-H,O-Z)
      
      Real(8), allocatable :: S_orb(:), SS_orb(:)

      if(allocated(S_orb)) Deallocate(S_orb,SS_orb)
      Allocate(S_orb(nbf),SS_orb(nbf))

      write(muc,'(/a)') 'Physical orbitals:'

      S_orb = 0.d0
      Do jp=1,npert

       SS_orb = 0.d0; SS = 0.d0
       Do ic = ippert(jp-1)+1,ippert(jp)
        S = WC(ic)*WC(ic); SS = SS + S
        Call Get_cfg_LS(ic)
        ip = ip_state(ic)
        Do i=1,no; ip=ip+1; ii=IEF(IP_orb(ip))   
         SS_orb(ii) = SS_orb(ii) + iq(i)*S 
        End do
       End do
!      SS_orb = SS_orb / SS  ! ????

       Do i=1,nbf
        S_orb(i)=max(S_orb(i),SS_orb(i))       

       End do

      End do

! ... choose the substitution orbital with biggest overlap:

      Do ii=1,nbf;  if(S_orb(ii).lt.eps_phys) Cycle

       SM = 0.d0; jj = 0 
       Do k=ncore+1,nbf
        if(iech(k).ne.1)      Cycle
        if(lbs(k).ne.lbs(ii)) Cycle
        s_ovl = abs(OBS(k,ii))
        if(S_ovl.lt.SM) Cycle
        SM = s_ovl; if(s_ovl.gt.eps_sub) jj = k
       End do
       
       if(ii_sub.gt.0.and.jj.eq.0) Cycle

       if(jj.eq.0) then
        Call Add_sub_orbital(ii,jj) 
        SM = OBS(ii,jj)
       end if

       write(muc,'(a4,f8.1,5x,a4,f8.3)') ebs(ii),S_orb(ii),ebs(jj),SM
       write(nuo,'(a4,f8.1,5x,a4,f8.3)') ebs(ii),S_orb(ii),ebs(jj),SM

      End do          

      write(nuo,'(a)') '*'
      write(muc,'(a)') '*'

      End Subroutine SUB_phys_pert
