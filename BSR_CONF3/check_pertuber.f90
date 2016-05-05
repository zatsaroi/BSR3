!======================================================================
      Subroutine Check_perturbers
!======================================================================     
!     read the perturber configurations, substitute physical orbitals
!     and check the "double"  configurations if any
!     We are supposed that pertuber expansions are ordered
!----------------------------------------------------------------------
      Use bsr_conf
      Use target; Use channel; Use conf_LS; Use orb_LS; Use phys_orb_LS

      Implicit none
      Character(80) :: AS
      Integer :: i,j,ii,jj,is,ic,ip,jp  
      Real(8) :: W
      Integer, external :: Ifind_cfg_LS, Ifind_nlk, Ifind_position 

      if(ncp.eq.0) Return

! ... read substitution orbitals for pertuber:

      Call Read_sub_pert_LS(nuo,ilsp)

! ... find substitution pointers for pertuber orbitals:    

      if(nwf_pert.gt.nwf_targ) IEF(nwf_targ+1:nwf_pert)=0

      Do i=1,npert_sub; ii=np_phy(i); jj=np_sub(i)
       if(IEF(ii).eq.0) then
        IEF(ii)=jj
       else
        if(IEF(ii).ne.jj) &
        Stop 'troubles with sustitution orbitals for pertuber'
       end if
      End do

! ... substitute the phys. orbitals in pertuber:

      if(mcfg.gt.ncfg) WC(ncfg+1:mcfg) = 0.d0

      AF = trim(AFP)//'.c'
      Call Check_file(AF)
      Open(nuc,file=AF,status='OLD')
      Call R_pert(nuc)
      Do ip = 1,npert; jp=1; if(ip.gt.1) jp=ippert(ip-1)+1
       
      rewind(nuc)
      is = 0
    1 read(nuc,'(a)',end=2) AS
      if(AS(1:1).eq.'*') go to 2
      if(AS(5:5).ne.'(') go to 1
      read(AS,'(a64,f11.8)') CONFIG,W
      read(nuc,'(a)') COUPLE
      is = is + 1
      if(is.lt.jp) go to 1
      Call Decode_c

      Do i=1,no
       ii=Ifind_nlk(nn(i),ln(i),kn(i),1)
       j=IEF(ii); if(j.eq.0) Exit
       nn(i)=NEF(j); ln(i)=LEF(j); kn(i)=KEF(j)
      End do

      if(j.eq.0) then
       write(pri,'(/a,i5)') &
        'Cannot find physical configuration for pertuber',ip
       Stop 'Check the pertuber'
      end if

      ic = Ifind_cfg_LS()
      WC(ic) = WC(ic) + 1.d0

    2 Continue

      End do   !  over perturbers, ip

      i = Ifind_position(nuc,'Spectroscopic configurations:')
      if(i.eq.0) then
      Close(nuc); Open(nuc,file=AF,position='APPEND')
       write(nuc,'(/a/)') 'Spectroscopic configurations:'
      else
       read(nuc,*); read(nuc,*)
      end if

      Do ic=ncfg_sct+1,ncfg
       Call Pri_conf(nuc,ic,WC(ic))
      End do
      write(nuc,'(a)') '*'

      if(ncfg-ncfg_sct.ne.npert) then
       write(pri,'(/a,i5)') &
        'Check the double pertubers for given partial wave'
        Stop 'Check the double pertubers'
      end if       

      write(pri,'(/a,T40,i8)') 'Number of perturber configurations:', ncp
      write(pri,'( a,T40,i8)') 'Number of physical perturbers:', &
                                ncfg-ncfg_sct

      End Subroutine Check_perturbers
