!======================================================================
      Subroutine Def_orth_cond
!======================================================================
!     define orthogonal conditions to avoid "over-compensation"
!----------------------------------------------------------------------
      Use bsr_conf
      Use target; Use channel; Use conf_LS;; Use orb_LS; Use phys_orb_LS

      Implicit none
      Integer :: i,j, ich,it,ii,ie,ic,ic1,ic2, io 

      if(debug.gt.0) &
      write(pri,'(/a/)') 'Define orth.conditions:'

! ... find orth. conditions by checking the compensation configurations:

      if(mcfg.gt.ncfg_pert) WC(ncfg_pert+1:mcfg)=0.d0

      Do ich = 1,nch;  it=iptar(ich)
        ii = ipch(ich); ii_comp=ii       ! position of channel orbital

        ic1=1; if(ich.gt.1) ic1=ipconf(ich-1)+1; ic2=ipconf(ich)
        ic1=ic1+ncfg_targ; ic2=ic2+ncfg_targ

        Do ic = ic1,ic2
         Do io=1,nphys_sub; ie=jp_sub(io); ie_comp=ie
          i = max(ii,ie); j=min(ii,ie); 

          if(IORT(i,j).le.0) Cycle

          if(debug.gt.0) &
          write(pri,'(a,i4,2x,2a6)')  'it=',it,ELF(ii),ELF(ie) 

          Call Get_cfg_LS(ic)
          WC_comp = WC(ic)*WC(ic)     

          if(debug.gt.0) Call Pri_conf (pri,ic,WC_comp)

          igen_conf = 0
          Call Gen_conf
          if(igen_conf.eq.0) IORT(i,j)=0

          if(debug.gt.0) write(pri,'(32x,a,i2)') 'IORT = ',IORT(i,j)
          
         End do
        End do ! ic, over configurations 

      End do  ! ich, over channels

      ncfg_comp = ncfg;  lcfg_comp=lcfg
      write(pri,'(/a,T40,i8)') &
       'Number of compensation configurations:', ncfg_comp-ncfg_pert

      if(ncfg_comp-ncfg_pert.le.0) Return

      write(AF,'(a,i3.3,a)') 'pert_comp.',ilsp
      Open(nuc,file=AF)
      write(nuc,'(12x,a3,f16.8)') Tpar,Etarg(1)
      write(nuc,'(a60)') CLOSED 
      Do ic=ncfg_pert+1,ncfg_comp
       Call Get_cfg_LS(ic)
       Call Incode_c
       write(nuc,'(a64,F10.3)') CONFIG,WC(ic)
       write(nuc,'(a)') trim(COUPLE)
      End do
      write(nuc,'(a)') '*' 
      Close(nuc)

      End Subroutine Def_orth_cond




















