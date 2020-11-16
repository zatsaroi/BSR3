!======================================================================
      Subroutine Def_orbitals
!======================================================================
!     This routine prepares the orbital-related arrays
!----------------------------------------------------------------------
      Use bsr_mchf
      Use orb_ls

      Implicit none
      Integer :: i,j,ic,ip
      Real(8) :: S
      Character(200) :: AS
  
      Call alloc_spline_orbitals(nwf)
      nbf = nwf
      Allocate(iord(nbf),qsum(nbf),dpm(nbf),clsd(nbf),e(nbf),iprm(ns,nbf))
      iord = 0; qsum = 0.d0; dpm = 0.d0; clsd = .FALSE.; e = 0.d0; iprm = 1

      Do i=1,nbf
       nbs(i)=nef(i);  lbs(i)=lef(i);  ebs(i)=ELF(i);  kbs(i)=0
      End do

      ncore = nclosd; kclosd = nclosd

! ... varried orbitals:
  
      Call Def_nit
  
! ... find physical orbitals:

      Call Def_physical

! ... define qsum:            ?  just initial estimation                     

!      Do i=1,nlevels; j=level(i)
!       WC(j) = WC(j) + weight(i)        !  * ???
!      End do 

      S = sqrt(SUM(WC(1:ncfg)*WC(1:ncfg)))
      if(S.eq.0.d0) then
       WC(1)=1.d0; if(ncfg.gt.1) WC(2:ncfg)=0.01d0
      end if
      S = sqrt(SUM(WC(1:ncfg)*WC(1:ncfg)))
      WC(1:ncfg)=WC(1:ncfg)/S

      Do ic = 1,ncfg; S = WC(ic)*WC(ic) 
       Call Get_cfg_LS(ic)     
       ip = ip_state(ic)
       Do i=1,no; ip=ip+1; j=IP_orb(ip)
        qsum(j) = qsum(j) + S*iq(i)  
       End do
      End do

      Do i=1,ncore
       qsum(i) = 4*lbs(i)+2;  clsd(i) = .TRUE.
      End do
   
write(log,'(a,20E15.5)') 'QSUM: ', qsum (1:nbf)

! ... define closed shells:   ???

      Do i=ncore+1,nbf
       clsd(i) = .FALSE.; S = 4*lbs(i)+2
       if(abs(qsum(i)-S).lt.0.1) clsd(i) = .TRUE.
      End do

      Call Boundary_conditions 

! ... define if core is fixed:

      icore = 0;  if(ncore.gt.0) icore = sum(iord(1:ncore))    

! ... print orbitals:

      if(ncore.gt.0) then
       AS = ' '; ip = 0
       Do i = 1,ncore
        write(AS(ip+1:ip+6),'(2x,a4)') ebs(i); ip=ip+6
       End do
       write(log,'(/a,T20,a)') 'core orbitals:',trim(AS)
      end if

      AS = ' '; ip = 0
      Do i = ncore+1,nbf
       write(AS(ip+1:ip+6),'(2x,a4)') ebs(i); ip=ip+6
      End do
      write(log,'(/a,T20,a)') 'valence orbitals:',trim(AS)

      AS = ' '; ip = 0
      Do i = 1,nbf; if(iphys(i).eq.0) Cycle
       write(AS(ip+1:ip+6),'(2x,a4)') ebs(i); ip=ip+6
      End do
      write(log,'(/a,T20,a)') 'physical orbitals:',trim(AS)

      AS = ' '; ip = 0
      Do i = 1,nbf; if(ivaried(i).eq.0) Cycle
       write(AS(ip+1:ip+6),'(2x,a4)') ebs(i); ip=ip+6
      End do
      write(log,'(/a,T20,a)') 'varied orbitals:',trim(AS) 

      End Subroutine Def_orbitals


!======================================================================
      Subroutine Def_nit
!======================================================================
!     This routine determine the orbitals to be optimized
!     nit      -  number of optimized oprbitals
!     iord(:)  -  =1 means optimized orbital
!----------------------------------------------------------------------
      Use bsr_mchf,  string => avaried,  nit => varied

      Implicit none
      Integer :: i,j,n,l,k,ip,ii,ios,iset
      Integer, external :: Ifind_bsorb, Ipointer
      Character(4) :: EL

      Call Clean_a(string)
      iord = 0
      if(string(1:3)=='ALL' .or. string(1:3)=='all' .or. &        ! all
       LEN_TRIM(string) == 0 ) then                              
       nit = nbf 
       Do i=1,nbf; iord(i)=i; End do
      elseif(string(1:4)=='NONE' .or. string(1:4)=='none') then   ! none
       nit = 0 
      elseif (INDEX(string,'n=') /= 0) then                       ! n=
       i = INDEX(string,'=') 
       read(string(i+1:),*) n 
       nit = 0
       Do i=1,nbf
        if(nbs(i).ne.n) Cycle
        nit = nit + 1
        iord(nit) = i
       End do                                                     
      elseif (INDEX(string,'n>') /= 0) then                       ! n>
       i = INDEX(string,'>') 
       read(string(i+1:),*) n 
       nit = 0
       Do i=1,nbf
        if(nbs(i).le.n) Cycle
        nit = nit + 1
        iord(nit) = i
       End do                                                   
      elseif (INDEX(string,'=') /= 0) then                        ! =last
       i = INDEX(string,'=') 
       read(string(i+1:),*) nit 
       k=0; Do i=nbf-nit+1,nbf; k=k+1; iord(k)=i; End do
      else

       Do nit=1,nbf                                               ! list
        read(string,*,iostat=ios) (EL,i=1,NIT)
        IF (ios /= 0) Exit
       End do
       nit = nit - 1
       ip = 0
       Do ii=1,nit
        read(string,*) (EL,j=1,ii)

        Call EL4_NLK(EL,n,l,iset)
        j = Ifind_bsorb(n,l,iset)
        if(j.gt.0) then; ip=ip+1; iord(ip)=j; end if
       End do
      end if

      Allocate(ivaried(nbf));  ivaried = 0
      Do i=1,nbf
       ivaried(i) = Ipointer(nbf,iord,i)
      End do

      if(debug.gt.0) then
       write(log,'(a,100i3)') 'iord:  ',iord 
       write(log,'(a,100i3)') 'varied:',ivaried
      end if

      End Subroutine Def_nit


!======================================================================
      Subroutine Def_physical
!======================================================================
! ... read or try to find the physical orbitals; approximately,
! ... they are orbitals in the first (leading?) configurations
! ... for each block
!----------------------------------------------------------------------
      Use bsr_mchf

      Implicit none
      Character(4) :: EL 
      Integer :: i,j,start,ip,ic, n,l,iset
      Integer, external :: Ifind_bsorb

      Allocate(iphys(1:nbf)); iphys=0              ! add ipys to df_orbitals ???
      if(ncore.gt.0) iphys(1:ncore)=1
      physical = ' '
      Call Read_string(inp,'physical',physical)    
      Call Read_aarg('physical',physical)

! ... if given - just read:

      if(len_trim(physical).gt.0) then
       start = 1; ip = 0
       Do  
        i = index(physical(start:),',')
        if (i /= 0 .or. LEN_TRIM(physical(start:)) /= 0) then
         read(physical(start:),*) EL
         Call EL4_NLK(EL,n,l,iset)
         j = Ifind_bsorb(n,l,iset)
         if(j.gt.0) iphys(j) = 1
        end if
        start = start + i 
        if(i == 0 .or.  LEN_TRIM(physical(start:)) == 0) Exit
       End do
       Call Clean_a(physical)
       Return
      end if 

! ... if not given - try to anticipate:
! ... we suupose that the main configurations are in the right place ...

      Do i = 1,nlevels; ic=level(i)
       Call Get_cfg_LS(ic)     
       ip = ip_state(ic)
       Do j=1,no; ip=ip+1; iphys(IP_orb(ip))=1; End do
      End do

      start=1                                          ! ???
      Do i = ncore+1,nbf; if(iphys(i).eq.0) Cycle
       write(physical(start:start+5),'(a5,a1)') ebs(i),','
       start = start+6
      End do

      i = Len_trim(physical); physical(i:i)=' '
      Call Clean_a(physical)

      End Subroutine Def_physical
