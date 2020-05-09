!======================================================================
      Subroutine SUB_phys_orb
!======================================================================
!     define physical orbitals for given target state
!----------------------------------------------------------------------
      Use bsr_prep

      Implicit real(8) (A-H,O-Z)
      
      Real(8), allocatable :: S_orb(:)

      if(allocated(S_orb)) Deallocate(S_orb)
      Allocate(S_orb(nbf))

! ... find the occupation numbers for all orbitals:

      S_orb = 0.d0; SS = 0.d0
      ii = 5; if(ncgf.lt.5) ii=ncfg
      Do ic = 1,ii
       S = WC(ic)*WC(ic);  SS = SS + S        
       Call Get_cfg_LS(ic)
       ip = ip_state(ic)
       Do i=1,no; ip=ip+1; ii=IEF(IP_orb(ip))   
        S_orb(ii) = S_orb(ii) + iq(i)*S 
       End do
      End do

      S_orb = S_orb / SS;   SN = SS;   SS = 0
      write(muc,'(/a)') 'Physical orbitals:'

! ... choose the substitution orbital with biggest overlap:

      Do ii=1,nbf;  if(S_orb(ii).lt.eps_phys/SN) Cycle

       SM = 0.d0; jj = 0 
       Do k=ncore+1,nbf
        if(iech(k).ne.1)      Cycle
        if(lbs(k).ne.lbs(ii)) Cycle
        s_ovl = abs(OBS(k,ii))
        if(s_ovl.lt.SM)       Cycle
        SM = s_ovl; if(s_ovl.gt.eps_sub) jj = k
       End do

       if(ii_sub.gt.0.and.jj.eq.0) Cycle

       if(jj.eq.0) then
         Call Add_sub_orbital(ii,jj) 
         SM = OBS(ii,jj)
       end if

       write(muc,'(a4,f8.1,5x,a4,f8.3)') ebs(ii),S_orb(ii),ebs(jj),SM
       write(nuo,'(a4,f8.1,5x,a4,f8.3)') ebs(ii),S_orb(ii),ebs(jj),SM

       SS = SS + S_orb(ii)
       if(nelc-SS.lt.eps_phys) Exit
      End do          

      write(muc,'(a,20x,f8.3)') '*'    
      write(nuo,'(a,20x,f8.3)') '*'    

      End Subroutine SUB_phys_orb
