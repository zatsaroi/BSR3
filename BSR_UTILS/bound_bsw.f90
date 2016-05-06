!======================================================================
!     UTILITY       B O U N D _ B S W                      version 3
!
!               C O P Y R I G H T -- 2016
!
!     Written by:   Oleg Zatsarinny
!                   email: oleg_zoi@yahoo.com 
!
!======================================================================
!     converts bound.nnn files to (c- + bsw-) files
!======================================================================
!
!     INPUT FILES:
!
!      bound_bsw.inp - input data
!      target  - information about channels
!      knot.dat - B-spline information
!      cfg.nnn  - list of configurations for the partial waves nnn
!      bound.nnn  - file with solutions for the partial wave nnn
!      target.bsw  - target B-spline w.f.
!      pert_nnn.bsw - perturber orbitals
!
!     OUTPUT FILES:
!
!      c- and bsw- files for the indicated states
!      bound_bsw.log - running information
!----------------------------------------------------------------------
      Use channels; Use target; Use conf_LS; Use orb_LS
      Use spline_param; Use spline_atomic; Use spline_grid
      Use spline_orbitals

      Implicit real(8) (A-H,O-Z)

! ... files:

      Integer :: inp=5;    Character(20) :: AF_inp = 'bound_bsw.inp'
      Integer :: pri=6;    Character(20) :: AF_log = 'bound_bsw.log'
      Integer :: nuw=11;   Character(20) :: AF_w   = 'target.bsw'
      Integer :: nut=12;   Character(20) :: AF_t   = 'target'
      Integer :: nuc=13;   Character(20) :: AF_c   = 'cfg.nnn'
      Integer :: nub=14;   Character(20) :: AF_b   = 'bound.nnn'

      Character(20) :: AF,BF,CF,   AFb=' ', mode=' '
      Character( 3) :: ALSP
      Character(64) :: Lab
      Character(80) :: AS
      Character( 4), external :: ELF4

      Real(8) :: eps_c = 1.d-8

      Real(8), allocatable :: A(:), WCC(:), v(:)
      Integer, allocatable :: ipt(:)

! ... short instructions:

      Call get_command_argument(1,AF)  
      if(AF.eq.'?') then
       write(*,'(/a)') 'bound_bsw extracts solution and records it as pair of c- and bsw-files'   
       write(*,'(/a,a)') 'Call as: ',  &
        'bound_bsw  klsp=  sol=... name=... [mode=...]' 
       write(*,'(/a,a,a)') 'or use bound_bsw.inp with ',  &
        'klsp  sol  name  for each state in one line, with * marks end of the list' 
       write(*,'(/a)') 'klsp - partial wave index'
       write(*,'( a)') 'sol  - solution index'
       write(*,'( a)') 'name - results will be in name.c and name.bsw files'
       write(*,'( a)') 'mode - all names have the same style: [mode]_klsp_sol'
       Stop 
      end if

      Call Read_aarg('mode',mode)
      Call Read_rarg('eps_c',eps_c)
      Call Read_aarg('AFb',AFb)

!----------------------------------------------------------------------
! ... input file:          

      jnp = Icheck_file(AF_inp)
      if(len_trim(AFb).ne.0) jnp=0
      if(jnp.ne.0) Open(inp,file=AF_inp)
      Open(pri,file=AF_log)

!----------------------------------------------------------------------
! ... read target and channel information from "target":
   
      Call Check_file(AF_t)
      Open(nut,file=AF_t)
      Call R_target(nut);
      Call R_channels(nut)
      Close(nut)
     
!----------------------------------------------------------------------
! ... sets up grid points and initializes the values of the spline
! ... according to the file "knot.dat" : 
    
      Call define_grid(z);   Call define_spline;   Allocate(v(ns))

!----------------------------------------------------------------------
! ... target orbitals:

      Call Allocate_bsorb(ibf)
      Call Check_file(AF_w)
      Open(nuw,file=AF_w,form='UNFORMATTED')
      Call Read_bsw(nuw)
      Close(nuw)

      nbf_targ=nbf
!----------------------------------------------------------------------
! ... loop for all input states:

      klsp=0; Call Read_iarg('klsp',klsp); ilsp = klsp
      is=0;   Call Read_iarg('sol' ,is)  
      BF=' '; Call Read_aarg('name',BF)   

      if(jnp.ne.0) rewind(inp)
   10 Continue

      if(jnp.ne.0.and.ilsp.eq.0) then

       read(inp,'(a)',end=20,err=10) AS              
       if(AS(1:1).eq.'*') go to 20
       if(Len_trim(AS).eq.0) go to 10
       read(AS,*,err=10) klsp,is,BF

      elseif(ilsp.eq.0) then
       Stop 'not enough input data'
      end if

      CF = BF
      if(Len_trim(mode).gt.0) then
       write(BF,'(a,a1,i3.3,a1,i3.3)') trim(mode),'_',klsp,'_',is
      end if
      if(klsp.eq.0) Stop 'klsp=0' 
      if(is.eq.0)   Stop 'sol=0'  
      if(Len_trim(BF).eq.0) Stop 'name - ?'                         

      write(ALSP,'(i3.3)') klsp
      nbf=nbf_targ

! ... pertuber orbitals:

      if(nwp(klsp).gt.0) then
       AF = trim(AFP(klsp))//'.bsw'
       Call Check_file(AF)
       Open(nuw,file=AF,form='UNFORMATTED')
       Call Read_bsw(nuw)
       Close(nuw)
      end if      

! ... allocate space for outer orbitals:

      nbt = nbf; kch = nch(klsp)
      if(mbf.lt.nbt+kch) Call Allocate_bsorb(nbt+kch)
      nbf = nbt+kch

! ... read configurations:

      write(AF_c,'(a,i3.3)') 'cfg.',klsp
      Open(nuc,file=AF_c,status='OLD')
      Call R_CLOSED(nuc)
      ncfg=0; lcfg=0; Call Add_conf_LS(nuc,0)
      Close(nuc)
  
      write(pri,'(/a,a)') 'Partial wave:  ',trim(AF_c)

      write(pri,'(/a,i6,a)') 'ncfg = ',ncfg,' - number of configurations'
      write(pri,'( a,i6,a)') 'nbt  = ',nbt, ' - number of target orbitals'
      write(pri,'( a,i6,a)') 'nbf  = ',nbf, ' - total number of orbitals'

      if(allocated(WCC)) deallocate(WCC); allocate(WCC(ncfg))
     
!----------------------------------------------------------------------
! ... read bound.nnn or ubound.nnn:

      write(AF_b,'(a,i3.3)') 'ubound.',klsp
      iform = Icheck_file(AF_b)
      if(iform.ne.0.and.len_trim(AFb).eq.0) then
       Open(nub,file=AF_b,form='UNFORMATTED')
       iform = 10
      else
       write(AF_b,'(a,i3.3)') 'bound.',klsp
       if(len_trim(AFb).ne.0) AF_b = AFb
       iform = Icheck_file(AF_b)
       if(iform.ne.0) then
        Open(nub,file=AF_b)
        iform = 20
       else
        write(*,*) 'Can not find ',trim(AF_b)
        Stop ' '
       end if 
      end if

      rewind(nub)
      if(iform.eq.10) &
      read(nub) ii,kch,kcp,nhm,nstate,ILT,IST,iparity_state
      if(iform.eq.20) &
      read(nub,*) ii,kch,kcp,nhm,nstate,ILT,IST,iparity_state

      if(ii.ne.ns)           Stop ' ns  --> ?'
      if(kch.ne.nch(klsp))   Stop ' nch --> ?'
      if(kcp.ne.npert(klsp)) Stop ' ncp --> ?'
      if(nhm.ne.ns*kch+kcp)  Stop ' nhm --> ?'

      write(pri,'(a,i6,a)') 'kch  = ',kch,' - number of channels'
      write(pri,'(a,i6,a)') 'kcp  = ',kcp,' - number of perturbers'

      if(is.lt.1.or.is.gt.nstate) Stop ' required solution is out of range '
          
      if(allocated(A)) Deallocate(A); Allocate(A(nhm))

      if(iform.eq.10) then
       Do i = 1,is
        read(nub) ii,Lab; read(nub) ET; read(nub) A
       End do 
      else
       Do i = 1,is
        read(nub,*) ii,Lab;  read(nub,*) ET;  read(nub,'(5D15.8)') A
       End do 
      end if

      Close(nub)
        
      write(pri,'(/a,i6)'   ) 'Solution:',is
      write(pri,'( a,T15,a)') 'Label: ',trim(LAB)
      write(pri,'( a,f16.8)') 'Energy:   ',ET

!----------------------------------------------------------------------
!  ... transfer the solution 'is' to the p-functions:

        ic1 = 0; ic2 = 0
        
        write(pri,'(/a/)') 'Channel decomposition and renormalization:'

        Do ich = 1,kch
         ishft=(ich-1)*ns;  v(1:ns)=A(ishft+1:ishft+ns)
         ii=nbt+ich; pbs(1:ns,ii) = v;  mbs(ii)=ns-1 
         s = QUADR(ii,ii,0); ss = sqrt(s)
         if(V(lch(klsp,ich)+2).lt.0.d0) ss = -ss
         v = v / ss; pbs(1:ns,ii) = v; mbs(ii)=ns-1
         ic1=ic2+1; ic2=ipconf(klsp,ich)
         WCC(ic1:ic2) = WC(ic1:ic2) * ss
         s = SUM(WCC(ic1:ic2)*WCC(ic1:ic2))
         write(pri,'(a5,2i5,2f16.8)') elc(klsp,ich),ic1,ic2,s,ss
        End do

!  ...  weights of perturbers:

        write(pri,'(/a/)') 'Perturber contribution and coefficient:'

        i1=ipert(klsp)+1; i2=ipert(klsp)+npert(klsp)
        if(kcp.gt.0) then
         ishft =kch*ns
         jshft =ipconf(klsp,kch)
         j = ipert(klsp)
         Do i=1,npert(klsp); 
          ss = A(ishft+i)
          ic1=jshft+ippert(j+i-1)+1
          if(i.eq.1) ic1=jshft+1
          ic2=jshft+ippert(j+i)
          WCC(ic1:ic2) = WC(ic1:ic2) * ss
          s = SUM(WCC(ic1:ic2)*WCC(ic1:ic2))
          write(pri,'(i4,1x,2i5,2f16.8)') i,ic1,ic2,s,ss
         End do         
        end if

!----------------------------------------------------------------------
!  ...  output c - file for given solution:

        AF = trim(BF)//'.c'     
        Open(nuc,file=AF)

        if(allocated(ipt)) Deallocate(ipt); Allocate(ipt(ncfg))
        Call SORTA(ncfg,WCC,ipt)

        if(IST.gt.0) then
         write(nuc,'(15x,f16.8,5x,a,a)') ET,'conf: ',trim(LAB)
        else
         write(nuc,'(15x,f16.8,a,i3,5x,a,a)') ET, '  2J =',ILT,'conf: ',trim(LAB)
        end if        
        write(nuc,'(a)') trim(CLOSED)

        Do jc = 1,ncfg; ic=IPT(jc)
         if(abs(WCC(ic)).lt.eps_c) Cycle
         Call Get_cfg_LS(ic)
         Call Prj_conf_LS (nuc,WCC(ic))
        End do
        write(nuc,'(a)') '*'

! ... insersion to cover the case with the same configurations in pertuber:

        if(ncp(klsp).gt.0) then

        ncfg=0; lcfg=0; WC = 0.d0

        Call RR_conf_LS(nuc,0)

        Call SORTA(ncfg,WC,ipt)

        rewind(nuc)
        if(IST.gt.0) then
         write(nuc,'(15x,f16.8,5x,a,a)') ET,'conf: ',trim(LAB)
        else
         write(nuc,'(15x,f16.8,a,i3,5x,a,a)') ET, '  2J =',ILT,'conf: ',trim(LAB)
        end if        
        write(nuc,'(a)') trim(CLOSED)

        Do jc = 1,ncfg; ic=IPT(jc)
         if(abs(WC(ic)).lt.eps_c) Cycle
         Call Get_cfg_LS(ic)
         Call Prj_conf_LS (nuc,WC(ic))
        End do
        write(nuc,'(a)') '*'

        end if

        Close(nuc)
                
!----------------------------------------------------------------------
! ...  output bsw-file:
        
        Do ich = 1,kch
         ebs(nbt+ich)=ELC(klsp,ich)
        End do

        AF = trim(BF)//'.bsw'     
        Open(nuw,file=AF,form='UNFORMATTED')
        Do i = 1,nbf
         write(nuw) ebs(i),z,h,hmax,rmax,ks,ns,mbs(i),t
         write(nuw) (pbs(j,i),j=1,mbs(i))
        End do

       write(pri,'(/80(''-'')/)')

       if(jnp.ne.0.and.ilsp.eq.0) go to 10  !  new state

   20 Continue

      End   !  utility BOUND_BSW
          

!======================================================================
      Subroutine Read_bsw(nu)
!======================================================================
!     read B-spline radial orbitals from unit 'nu' and omit 
!     the orbitals with nlk already in memory 
!----------------------------------------------------------------------
      USE spline_param
      USE spline_atomic
      USE spline_orbitals

      Implicit real(8) (A-H,O-Z)

      Character elw*4

      rewind(nu)

    1 READ(nu,end=2) elw,zw,hw,hmw,rmw,ksw,nsw,mw

      if(zw.ne.z) Stop ' Rb_wfn:  z <> zw'
      if(hw.ne.h) Stop ' Rb_wfn:  h <> hw'
      if(abs(hmw-hmax).gt.1.d-15) Stop ' Rb_wfn:  hmw <> hmax'
      if(rmw.ne.rmax) Stop ' Rb_wfn:  rmw <> rmax'
      if(ksw.ne.ks) Stop ' Rb_wfn:  ksw <> ks'
      if(nsw.ne.ns) Stop ' Rb_wfn:  nsw <> ns'

      Call EL4_nlk(elw,nw,lw,kw);  ii = Ifind_bsorb(nw,lw,kw)

      if(ii.eq.0) then
        m=nbf+1; if(m.eq.mbf) Call Allocate_bsorb(mbf+jbf)
        mbs(m)  = mw
        nbs(m)  = nw
        lbs(m)  = lw
        kbs(m)  = kw
        ebs(m)  = elw
        read(nu) pbs(1:mw,m)
        if(mw.lt.ns) pbs(mw+1:ns,m) = 0.d0
        nbf=m
      else
        read(nu) (x,i=1,mw)
      end if

      go to 1
    2 Close(nu)

      End Subroutine Read_bsw


!====================================================================
      Subroutine RR_conf_LS(nu,kshift)
!====================================================================
!     reads configurations from c-file (unit in) and saves them in
!     module conf_LS
!     the dublicated configurations are omitted
!--------------------------------------------------------------------
      Use conf_LS
      
      Implicit none
      Integer, intent(in) :: nu
      Integer, intent(in) :: kshift
      Integer :: ic
      Integer, External :: Ifind_cfg_LS
      Character(100) :: AS
      Real(8) :: W
      
      if(mcfg.eq.0) Call alloc_cfg_LS(icfg)

      rewind(nu)
    1 read(nu,'(a)',end=2) AS
      if(AS(1:1).eq.'*') go to 2
      if(AS(5:5).ne.'(') go to 1
      read(AS,'(a)') CONFIG
      W = 0.d0
      if(LEN_TRIM(AS).gt.64) read(AS(65:),*) W
      read(nu,'(a)') COUPLE
      Call Decode_c
      kn=kn+kshift
      ic = Ifind_cfg_LS()
      WC(ic) = WC(ic) + W
      go to 1
    2 Continue

      End Subroutine RR_conf_LS
         