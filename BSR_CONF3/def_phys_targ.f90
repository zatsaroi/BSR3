!======================================================================
      Subroutine Def_phys_targ
!======================================================================     
!     define the physical configurations for target expansions;
!     we do it in two stage: define and record them in targ_nnn.c,
!     and then read all of them to avoid situation of missing targets
!     and put target cofigurations in needed order
!----------------------------------------------------------------------
      Use bsr_conf
      Use target; Use channel; Use conf_LS; Use orb_LS; Use phys_orb_LS

      Implicit none
      Character(80) :: AS
      Integer :: i,j,it,ii,jj,ic     
      Real(8) :: W, S
      Integer, external :: Jfind_cfg_LS, Iadd_cfg_LS, Ifind_nlk, &
                           Ifind_position, Ipointer

! ... read substitution orbitals:

      Call Check_file(AF_orb)
      Open(nuo,file=AF_orb)
      Call Read_sub_orb_LS(nuo,ntarg)

! ... find substitution pointers:    

      IEF=0
      Do i=1,nphys_orb; ii=ip_phy(i); jj=ip_sub(i)
       if(IEF(ii).eq.0) then
        IEF(ii)=jj
       else
        if(IEF(ii).ne.jj) write(*,*) 'troubles with sustitution orbitals'
       end if
      End do

      write(pri,'(/a,T33,i8)') 'number of substitution orbitals:',nphys_sub
      write(pri,*)
      write(pri,'(15a5)') (ELF(jp_sub(i)),i=1,nphys_sub)
!      write(pri,'(/a,T33)') &
!       'substitution orbitals are supposed to be in right "AFTER" order !!!'

      if(iread_targ.eq.1) go to 5 
!----------------------------------------------------------------------
! ... find and record main target configurations

      Allocate(ip_phys_conf(0:ntarg)); ip_phys_conf=0
      Allocate(jj_phys_conf(0:ntarg))
      Call alloc_cfg_LS(0)

      Do ic = 8,1,-1;   c_conf = 0.1*ic

      Do it=1,ntarg;     if(ip_phys_conf(it).ne.0) Cycle

      jj = jtarg(it)

      AF = trim(AFT(it))//'.c'
      Open(nuc,file=AF,status='OLD')
      rewind(nuc)
    1 read(nuc,'(a)',end=2) AS
      if(AS(1:1).eq.'*') go to 2
      if(AS(5:5).ne.'(') go to 1
      read(AS,'(a)') CONFIG
      read(AS(65:),*) W;  W = W * W
      read(nuc,'(a)') COUPLE

      if(W.lt.c_conf) go to 1   !  too small coefficients

      Call Decode_c
      Do i=1,no; np(i)=Ifind_nlk(nn(i),ln(i),kn(i),1); End do
      ! ... check if all orbitals are physical:
      Do i=1,no
       j=IEF(np(i)); if(j.eq.0) go to 1
       nn(i)=NEF(j); ln(i)=LEF(j); kn(i)=KEF(j); np(i)=j
      End do
      LS(no,5) = LS(no,5)*1000+jj

      i=Jfind_cfg_LS();  if(i.lt.0) go to 1            

      ip_phys_conf(it) = i
      WC(i) = sqrt(W) 

      i = Ifind_position(nuc,'Spectroscopic configuration')
      if(i.eq.0) then
      Close(nuc); Open(nuc,file=AF,position='APPEND')
       write(nuc,'(/a/)') 'Spectroscopic configuration:'
      else
       read(nuc,*); read(nuc,*)
      end if

      LS(no,5) = LS(no,5)/1000 
      Call Incode_c
      LS(no,5) = LS(no,5)*1000+jj
      
      S=1.d0; if(ip_phys_conf(it).eq.0) S=0.d0
      write(nuc,'(a64,2F11.8)') CONFIG,S,sqrt(w)
      write(nuc,'(a)') trim(COUPLE)

    2 Continue

      End do  !  over target states,  it

      if(Ipointer(ntarg,ip_phys_conf,0).eq.0) Exit

      End do  !  ic,  min. coef. 

!----------------------------------------------------------------------
!     record physical target configurations for inspection if needed

      open(nuc,file='target_conf') 
      Do it=1,ntarg;  ic=ip_phys_conf(it)
       if(ic.ne.0) then
        Call Get_cfg_LS(ic)
        LS(no,5) = LS(no,5)/1000 
        Call Incode_c
        write(nuc,'(a64,F11.8,2i8)') CONFIG,WC(ic),it,jtarg(ic)
        write(nuc,'(a)') trim(COUPLE)
       else      
        write(nuc,'(75x,i8,a)') it,' missing'
        write(nuc,*)
       end if 
      End do
      Close(nuc)

!-----------------------------------------------------------------------
      i = 0
      Do it = 1,ntarg; if(ip_phys_conf(it).ne.0) Cycle
       i = i + 1  
      End do
      if(i.ne.0) write(*,'(/a,i8/)') &
        'WARNING: bsr_conf failed to find phys.conf. for some targets:',i
      if(i.ne.0) write(pri,'(/a/)') &
        'WARNING: failed to find phys.conf. for targets:'
!-----------------------------------------------------------------------
! ... read dominant configurations for target states

    5 Continue
      Allocate(ic_targ(0:ntarg),jc_targ(0:ntarg))
      ic_targ = 0; jc_targ = 0      
      Call alloc_cfg_LS(0)

      Do it=1,ntarg

      AF = trim(AFT(it))//'.c'
      Open(nuc,file=AF,status='OLD')
      i = Ifind_position(nuc,'Spectroscopic configuration')
      if(i.eq.0) then
       write(pri,*) 'Cannot find spectroscopic configuration for target', it
      end if
   10 read(nuc,'(a)',end=20) AS
      if(AS(1:1).eq.'*') go to 20
      if(AS(5:5).ne.'(') go to 10
      read(AS,'(a)') CONFIG
      read(AS(65:),*) W
      read(nuc,'(a)') COUPLE
      Call Decode_c
      i=Iadd_cfg_LS()
      WC(i) = W 
   20 Continue

      jc_targ(it) = ncfg
      End do  !  over target states,  it

      ncfg_phys = ncfg;  lcfg_phys = lcfg  
      write(pri,'(/a,T33,i8)') &
       'number of phys. target config.s:',ncfg_phys

      Call Test_a

      End Subroutine Def_phys_targ


