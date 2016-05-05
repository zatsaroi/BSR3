!======================================================================      
      Subroutine SUB_LS
!======================================================================  
!     define channels orbitals in LS(J) case
!----------------------------------------------------------------------
      Use bsr_conf
      Use target; Use channel; Use conf_LS; Use orb_LS

      Implicit none
      Integer :: l,n,k, ll, lch_min,lch_max
      Integer :: IST_min,IST_max,ILT_min,ILT_max
      Integer :: i,it,jk,ip, ic,ic1,ic2
      Integer, external :: NEW_INDEX, Iadd_cfg_LS, Ifind_nlk, Icheck_del
      Real(8) :: CT

      Do it=1,ntarg
 
! ... range of total spin ( = ispar for LS case)

       IST_min = iabs(istarg(it)-2) + 1 
       IST_max = iabs(istarg(it)+2) - 1 
       Do IST = IST_min,IST_max,2
       if(max_ST.gt.0.and.IST.gt.max_ST) Cycle

! ... range of total L ( = 2*lpar+1 for LS case)

       if(Jpar.gt.0) then
        ILT_min = iabs(Jpar-IST) + 1 
        ILT_max = iabs(Jpar+IST) - 1 
       else
        if(IST.ne.ispar) Cycle  
        ILT_min = 2*lpar + 1 
        ILT_max = 2*lpar + 1 
       end if

       Do ILT = ILT_min,ILT_max,2

        if(max_LT.gt.0.and.ILT.gt.max_LT) Cycle

! ... range of small l:

        ll = (ILT-1)/2
        lch_min=iabs(ll-ltarg(it))
        lch_max=     ll+ltarg(it)
        if((-1)**lch_min*iptarg(it).ne.ipar) then
         lch_min=lch_min+1
         lch_max=lch_max-1
        end if
        if(lch_min.gt.lch_max) Cycle

! ... cycle over channels:

        Do l=lch_min,lch_max,2

         if(max_ll.gt.0.and.l.gt.max_ll) Cycle
         if(min_ll.gt.0.and.l.lt.min_ll) Cycle

         n = ICHAR('k')-ICHAR('1')+1
         k = NEW_INDEX(l,ksmax,nwf,LEF,KEF)
         ip = Ifind_nlk(n,l,k,2)
         jk = (ILT-1)*50 + IST

         if(max_it.gt.0.and.it.gt.max_it) Cycle
         if(Icheck_del(ilsp,l,it,jk).eq.1) Cycle

         if(nch.eq.mch) Call Allocate_channel(mch+jmch)

         nch = nch + 1
         lch(nch) = l
         iptar(nch)= it
         ipch(nch) = ip
         jkch(nch) = jk
         elc(nch) = ELF(ip)

! ...  define partial configurations:

         ic1=ic_targ(it-1)+1; ic2=ic_targ(it)
         CT = 0.d0
         Do ic=ic1,ic2
          Call Get_cfg_LS(ic)
          no=no+1; nn(no) = n; ln(no) = l; iq(no) = 1; kn(no) = k
          LS(no,1)=1; LS(no,2)=l+l+1; LS(no,3)=2
          LS(no,4)=ILT; LS(no,5)=IST
          i=Iadd_cfg_LS()
          WC(ncfg)=WC(ic); CT = CT + WC(ic)**2
         End do  ! over target configuration
         ipconf(nch)=ncfg-ncfg_targ
         if(abs(CT-1.d0).gt.c_norm) &
         write(pri,'(i5,2x,a4,3x,2a20,3x,a,f8.5,a)') &
          nch,ELC(nch),AFT(it),BFT(it),'norm =',CT,' - check the normalization'
        End do   ! over channels (small l)
        End do   ! over ILT (large l)
       End do    ! over IST (total spin)
      End do     ! over targets

      End Subroutine SUB_LS 

