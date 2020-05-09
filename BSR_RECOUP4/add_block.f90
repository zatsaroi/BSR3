!======================================================================
      Subroutine add_block(ilsp,ic,jc,met)
!======================================================================
!     add LS channel block (ic,jc) to all relevant JK channel blocks
!----------------------------------------------------------------------
      Use bsr_recoup

      Use target_ion,   only: ltarg_LS => ltarg, istarg_LS => istarg
      Use channels_ion, only: nlsp_LS => nlsp, &
                              ipar_LS => ipar, lpar_LS => lpar, ispar_LS => ispar, &
                              lch_LS  => lch, iptar_LS => iptar
      
      Implicit real(8) (A-H,O-Z)

      ILT = lpar_LS(ilsp); IST = ispar_LS(ilsp)

      li = lch_LS(ilsp,ic);  iti = iptar_LS(ilsp,ic);  ILI = ltarg_LS(iti); ISI = istarg_LS(iti) 
      lj = lch_LS(ilsp,jc);  itj = iptar_LS(ilsp,jc);  ILJ = ltarg_LS(itj); ISJ = istarg_LS(itj)

      Do ich = 1,nch; if(lch(ich).ne.li) Cycle

       it = iptar(ich)
       ip1 = ip_expn(it-1)+1; ip2=ip_expn(it)       
       ci = 0.d0
       Do i=ip1,ip2; if(it_expn(i).ne.iti) Cycle; ci=c_expn(i); Exit; End do
       if(ci.eq.0.d0) Cycle

       ti =  T_LS_jk(ILI+ILI+1,li+li+1,ISI,2,ILT+ILT+1,IST,jtarg(it),jkch(ich),jpar)
       if(ti.eq.0.d0) Cycle

      Do jch = 1,ich; if(lch(jch).ne.lj) Cycle

       if(imycase(ich,jch).eq.0) Cycle

       jt = iptar(jch)
       jp1 = ip_expn(jt-1)+1; jp2=ip_expn(jt)       
       cj = 0.d0
       Do j=jp1,jp2; if(it_expn(j).ne.itj) Cycle; cj=c_expn(j); Exit; End do
       if(cj.eq.0.d0) Cycle
       
       tj =  T_LS_jk(ILJ+ILJ+1,lj+lj+1,ISJ,2,ILT+ILT+1,IST,jtarg(jt),jkch(jch),jpar)
       if(tj.eq.0.d0) Cycle

       C = ci*ti*cj*tj 

       if(met.eq.1) then   !  channel block

        ij = icc(ich,jch)

        if(ij.le.0.or.ij.gt.iicc) &
           Call Stop_mpi(0,ij,'UPDATE_HX: ij index out of range')
  
        hcc(:,:,ij) = hcc(:,:,ij) + c*x(:,:)

       else

        acf(ich,jch,0:km) = acf(ich,jch,0:km) + c*bcf(ic,jc,0:km)

       end if

      End do; End do

      End Subroutine add_block

