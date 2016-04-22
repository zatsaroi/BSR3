!======================================================================
      Subroutine R_conf
!======================================================================
!
!     Read the configuration list from c-file (unit 'nuc'),
!     define there angular symmetries and compare with existing ones.
!     Prepare the angular arrays.
!
!     nub, nui - file units for old and new results, respectively
!
!----------------------------------------------------------------------

      USE inout_br;  USE configs; USE inter; USE term_exp

      Implicit none

      Character(26) :: AI  ! config. symmetry
      Character(40) :: BI  ! angular symmetry

      Character(26), Allocatable, Dimension(:) :: ASYM
      Character(40), Allocatable, Dimension(:) :: BSYM

! ... IT_conf(it) gives the configuration for given term 'it'

      Integer(4), Allocatable, Dimension(:) :: IT_conf
      Integer(4), External :: Iterm, Idef_ncfg, IC_calc
      Integer(4) :: i,j, ii,it, msymc,msymt,iconf,jterm,ndim1,ndim2

! ... define ncfg:

      ncfg=Idef_ncfg(nuc)
      if(ncfg.eq.0) Stop ' ncfg = 0, nothing to do '

! ... read old information, if any:

      if(new) then
       nsymt=0; nsymc=0; msymt = ncfg; msymc = ncfg
       if(Allocated(IT_conf)) Deallocate (IT_conf,ASYM,BSYM,IP_stat)
       Allocate(IT_conf(msymt),ASYM(msymc),BSYM(msymt),IP_stat(msymt))
       IT_conf = 0
      else
       read(nub) nsymt, nsymc
       write(pri,'(a/)')    ' int_bnk contains: '
       write(pri,'(a,i5))') ' number of configurations  = ',NSYMC
       write(pri,'(a,i5))') ' number of ang. symmetries = ',NSYMT
       msymt = nsymt + ncfg
       msymc = nsymc + ncfg
       if(Allocated(IT_conf)) Deallocate (IT_conf,ASYM,BSYM,IP_stat)
       Allocate(IT_conf(msymt),ASYM(msymc),BSYM(msymt),IP_stat(msymt))
       Call R_a (nub,mrecl,nsymc,ASYM)    ! read(nub) ASYM(1:nsymc)
       Call R_a (nub,mrecl,nsymt,BSYM)    ! read(nub) BSYM(1:nsymt)
       Call R_i4(nub,mrecl,nsymt,IT_conf) ! read(nub) IT_conf(1:nsymt)
      end if
      IP_stat = 0

!---------------------------------------------------------------------
! ... define new symmetries from c-file:

      rewind(nuc); parity=0
      Do 

       read(nuc,'(a)') CONFIG
       if(CONFIG(1:1).eq.'*') Exit
       if(CONFIG(5:5).ne.'(') Cycle
       read(nuc,'(a)') COUPLE
 
       Call Decode_c;  Call TEST_C

! ... check angular symmetry:

       Call Incode_conf(AI,BI)

       iconf=0
       Do i=1,nsymc
        if(AI.ne.ASYM(i)) Cycle; iconf=i; Exit
       End do

       if(iconf.eq.0) then
        nsymc=nsymc+1; iconf=nsymc; ASYM(nsymc)=AI
       end if

       jterm = 0
       Do i=1,nsymt
        if(BI.ne.BSYM(i) .or. IT_conf(i).ne.iconf) Cycle
        jterm=i; Exit
       End do

       if(jterm.eq.0) then
        nsymt=nsymt+1; BSYM(nsymt)=BI; IT_conf(nsymt)=iconf
        jterm = nsymt
       end if
       IP_stat(jterm) = 1       ! simply indicates existence

      End do  ! on ic

      write(pri,'(/a/)')   ' c-fail (cfg.inp) contains: '
      write(pri,'(a,i5))') ' number of atomic states   = ',NCFG
      write(pri,'(a,i5))') ' number of configurations  = ',NSYMC
      write(pri,'(a,i5))') ' number of ang. symmetries = ',NSYMT

      write(*,'(/a/)')   ' c-fail (cfg.inp) contains: '
      write(*,'(a,i5))') ' number of atomic states   = ',NCFG
      write(*,'(a,i5))') ' number of configurations  = ',NSYMC
      write(*,'(a,i5))') ' number of ang. symmetries = ',NSYMT

      if(nsymt.ge.ibc) Stop ' R_conf: nsymt > ibc '

!----------------------------------------------------------------------
! ... define IP_term and IC_term pointers:

      if(Allocated(IP_term)) Deallocate(IP_term,IC_term)
      Allocate(IP_term(nsymt),IC_term(nsymc))

      Call sortI(nsymt,IT_conf,IP_term)

      IC_term=0
      Do i=1,nsymt
       it=IP_term(i); ic=IT_conf(it)
       if(IC_term(ic).eq.0) then
        IC_term(ic) = i + ibc*i
       else
        j=mod(IC_term(ic),ibc); IC_term(ic) = j + ibc*i
       end if
      End do

      if(Allocated(IT_stat)) Deallocate(IT_stat)
      Allocate(IT_stat(nsymt)); IT_stat = IP_stat(1:nsymt)

!----------------------------------------------------------------------
!     define if we need additional calculations

      if(allocated(IT_oper)) Deallocate(IT_oper, IC_need, JC_need)
      Allocate(IT_oper(1:noper,nsymt*(nsymt+1)/2),IC_need(nsymc),  &
               JC_need(nsymc*(nsymc+1)/2))

      IT_oper = 0
      if(new) then
       icalc=.TRUE.; IC_need = 1; JC_need = 1
      else

       read(nub) ii; ndim1=noper; ndim2=nsymt*(nsymt+1)/2
       Call R_i1(nub,mrecl,ndim1,ndim2,ii,IT_oper) ! read(nub) (IT_oper(:,i),i=1,ii)

       IC_need=0; icalc=.FALSE.
       Do ic=1,nsymc
        Do jc=1,ic
         i = IC_calc(ic,jc); JC_need(ic*(ic-1)/2+jc) = i
         if(i.gt.0) then
		        IC_need(ic) = 1; IC_need(jc) = 1; icalc=.TRUE.
		       end if
        End do
       End do
      end if

      if(icalc) then
       write(pri,'(/a,a4/)') ' Need of additional calculations --> yes '
       write(*,'(/a,a4/)')   ' Need of additional calculations --> yes '
      else
       write(pri,'(/a,a4/)') ' Need of additional calculations --> no '
       write(*,'(/a,a4/)')   ' Need of additional calculations --> no '
      end if
      
      if(.NOT.icalc) then
       Deallocate(IT_conf, IP_stat, ASYM, BSYM)
       Return
      end if

!----------------------------------------------------------------------
! ... record new information:

      write(nur) nsymt, nsymc
      Call W_a (nur,mrecl,nsymc,ASYM)    ! write(nur) ASYM(1:nsymc)
      Call W_a (nur,mrecl,nsymt,BSYM)    ! write(nur) BSYM(1:nsymt)
      Call W_i4(nur,mrecl,nsymt,IT_conf) ! write(nur) IT_conf(1:nsymt)

!----------------------------------------------------------------------
! ... define angular arrays:

      if(Allocated(NOCCSH)) Deallocate(NOCCSH, NELCSH, NOCORB, &
               LSA, LSL, LSS, LPL, LPS, ILT_ic,IST_ic) 

      Allocate(NOCCSH(nsymt), NELCSH(nsh,nsymt), NOCORB(nsh,nsymt), &
               LSA(nsh,nsymt), LSL(nsh,nsymt), LSS(nsh,nsymt),      &
               LPL(nsh,nsymt), LPS(nsh,nsymt),                      &
               ILT_ic(nsymc),IST_ic(nsymc)  )

      Do it = 1,nsymt
     
       ic = IT_conf(it)
       Call Decode_conf(ASYM(ic),BSYM(it))

       noccsh(it)=no
       Do i=1,no
        nelcsh(i,it)= iq(i)
        nocorb(i,it)= ln(i)
        LSA(i,it) = Iterm (ln(i),iq(i),0,LS(i,1),LS(i,2),LS(i,3))
        LSL(i,it) = LS(i,2)
        LSS(i,it) = LS(i,3)
        LPL(i,it) = LS(i,4)
        LPS(i,it) = LS(i,5)
       End do

       ILT_ic(ic) = LS(no,4)
       IST_ic(ic) = LS(no,5)

      End do

      Deallocate(IT_conf, IP_stat, ASYM, BSYM)

      End Subroutine R_conf