!======================================================================
      Subroutine SUB_JK
!====================================================================== 
!     define channels orbitals in JK case
!----------------------------------------------------------------------
      Use bsr_conf
      Use target; Use channel; Use conf_LS; Use orb_LS

      Implicit none
      Integer :: kk,k_min,k_max
      Integer :: n,l,k, ll, lch_min,lch_max
      Integer :: IS, IL, ISTmin,ISTmax, ILTmin,ILTmax
      Integer :: ic,ic1,ic2, i,it,jt, ip
      Real(8) :: C,CC,CT

      Character(4), external :: ELF4
      Integer, external :: NEW_INDEX, Iadd_cfg_LS, Ifind_nlk, Icheck_del,ITRI
      Real(8), external :: T_LS_jk      

       K_min = iabs(Jpar-2) + 1 
       K_max = iabs(Jpar+2) - 1 

       Do it=1,ntarg

        JT = jtarg(it)

       Do KK = K_min,K_max,2

        lch_min=iabs(KK-JT)/2
        lch_max=iabs(KK+JT)/2-1
        if((-1)**lch_min*iptarg(it).ne.ipar)  lch_min=lch_min+1
        if((-1)**lch_max*iptarg(it).ne.ipar)  lch_max=lch_max-1
        if(lch_min.gt.lch_max) Cycle

        Do l=lch_min,lch_max,2

         if(max_ll.gt.0.and.l.gt.max_ll) Cycle
         if(min_ll.gt.0.and.l.lt.min_ll) Cycle

         if(Icheck_del(ilsp,l,it,kk).eq.1) Cycle

         if(nch.eq.mch) Call Allocate_channel(mch+jmch)

         nch = nch + 1; lch(nch) = l; iptar(nch)= it; ll = 2*l+1

         n = ICHAR('k')  ! -ICHAR('1')+1
         k = NEW_INDEX(l,ksmax,nwf,LEF,KEF)
         ip = Ifind_nlk(n,l,k,2)
         CC = 0.d0
         ipch(nch) = ip
         jkch(nch) = kk
         elc(nch)  = ELF(ip)

! ...  define partial configurations ...

         ic1=ic_targ(it-1)+1; ic2=ic_targ(it)

         CT = 0.d0
         Do ic = ic1,ic2
          Call Get_cfg_LS(ic)
          IL = LS(no,4); IS = LS(no,5); no=no+1
          CT = CT + WC(ic)**2
 
          ISTmin = iabs(IS-2)+1;  ILTmin = iabs(IL-ll)+1
          ISTmax = iabs(IS+2)-1;  ILTmax = iabs(IL+ll)-1

          Do IST = ISTmin,ISTmax,2
          Do ILT = ILTmin,ILTmax,2

!         the recoupling coefficient from LS- to jK-coupling:
!         < (l1,l2)L,(s1,s2)S;J | (((l1,s1)js,l)K,s2)J >

           if(ITRI(IST,ILT,Jpar).eq.0) Cycle
           C = WC(ic) * T_LS_jk (IL,ll,IS,2,ILT,IST,JT,KK,Jpar)
           CC = CC + C**2;  if(C.eq.0.d0) Cycle
           
           nn(no) = n; ln(no) = l; iq(no) = 1; kn(no) = k
           LS(no,1)=1; LS(no,2)=ll; LS(no,3)=2
           LS(no,4)=ILT; LS(no,5)=IST
           i=Iadd_cfg_LS(); WC(ncfg)=C
          End do   ! over ILT (large l)
          End do   ! over IST (total spin)
         End do   ! over target configuration (ic)
         ipconf(nch)=ncfg-ncfg_targ
!         if(abs(CT-1.d0).gt.c_norm) &
!         write(pri,'(i5,2x,a4,3x,2a20,3x,a,f8.5,a)') &
!          nch,ELC(nch),AFT(it),BFT(it),'norm =',CT,' - check the normalization'
        End do     ! over small l

       End do     ! over KK
       End do     ! over targets

      End Subroutine SUB_JK
