!======================================================================
      Subroutine Def_phys_targ
!======================================================================     
!     define the physical configurations for target expansions
!     and record them in the end of targ_nnn.c.
!----------------------------------------------------------------------
      Use bsr_conf
      Use target; Use channel; Use conf_LS; Use orb_LS; Use phys_orb_LS

      Implicit none
      Character(80) :: AS
      Integer :: i,j,it,ii,jj,ic,ic1,ic2     
      Real(8) :: W, S, SS
      Integer, external :: Jfind_cfg_LS, Iadd_cfg_LS, Ifind_nlk, &
                           Ifind_position, Ipointer

! ... read substitution orbitals from target_orb:

      Call Check_file(AF_orb)
      Open(nuo,file=AF_orb)
      Call Read_sub_orb_LS(nuo,ntarg)

! ... find substitution pointers:    

      IEF=0
      Do i=1,nphys_orb; ii=ip_phy(i); jj=ip_sub(i)
       if(IEF(ii).eq.0) then
        IEF(ii)=jj
       else
        if(IEF(ii).ne.jj) write(*,*) 'BSR_CONF: troubles with sustitution orbitals'
       end if
      End do

      write(pri,'(/a,T33,i8)') 'number of substitution orbitals:',nphys_sub
      write(pri,*)
      write(pri,'(15a5)') (ELF(jp_sub(i)),i=1,nphys_sub)

!------------------------------------------------------------------------------

      if(iread_targ.eq.1) go to 5   ! skip if physical configurartion are given

!------------------------------------------------------------------------------
! ... find and record main target configurations

      if(allocated(jc_targ)) Deallocate(jc_targ) 
      Allocate(jc_targ(0:ntarg));  jc_targ = 0      
      if(allocated(ic_targ)) Deallocate(ic_targ) 
      Allocate(ic_targ(0:ntarg));  ic_targ = 0      

      Call alloc_cfg_LS(0)

      Do it=1,ntarg

      SS = 0.d0

      AF = trim(AFT(it))//'.c'
      Open(nuc,file=AF,status='OLD')
      rewind(nuc)
    1 read(nuc,'(a)',end=2) AS
      if(AS(1:1).eq.'*') go to 2
      if(AS(5:5).ne.'(') go to 1
      read(AS,'(a)') CONFIG
      read(AS(65:),*) W;  W = W * W
      read(nuc,'(a)') COUPLE
      if(W.lt.c_conf)  go to 2   !  too small coefficients

      Call Decode_c
      Do i=1,no; np(i)=Ifind_nlk(nn(i),ln(i),kn(i),1); End do

      ! ... check if all orbitals are physical:

      Do i=1,no
       j=IEF(np(i)); if(j.eq.0) go to 1
       nn(i)=NEF(j); ln(i)=LEF(j); kn(i)=KEF(j); np(i)=j
      End do

      i = Iadd_cfg_LS();  jc_targ(it) = i; WC(i) = sqrt(W)

      SS = SS + W 

      if(SS.lt.C_phys) go to 1 

    2 Continue

      i = Ifind_position(nuc,'Spectroscopic configuration')
      if(i.eq.0) then
      Close(nuc); Open(nuc,file=AF,position='APPEND')
       write(nuc,'(/a,f15.3/)') 'Spectroscopic configuration:',SS
      else
       read(nuc,*); read(nuc,*)
      end if

      ic1 = jc_targ(it-1)+1; ic2 = jc_targ(it)

      if(ic1.gt.ic2) then
       write(pri,*) 'Can not find physical configuration for target ',it
       Stop 'Can not find physical configuration'
      end if

      S = sqrt(SS)

      Do ic = ic1,ic2; Call Pri_conf(nuc,ic,WC(ic)/S); End do
      write(nuc,'(a)') '*'

      End do  !  over target states,  it

      go to 30

!-----------------------------------------------------------------------
! ... read dominant configurations for target states

    5 Continue

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
!-------------------------------------------------------------------------

   30 Continue
      ncfg_phys = ncfg;  lcfg_phys = lcfg  
      write(pri,'(/a,T33,i8)') &
       'number of phys. target config.s:',ncfg_phys

      Call Test_a

      End Subroutine Def_phys_targ


