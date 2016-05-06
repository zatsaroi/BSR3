!======================================================================
      Subroutine Read_arg(nu)
!======================================================================
!     read input parameters
!----------------------------------------------------------------------
      Use bsr_mat;  USE spline_atomic;  Use target

      Implicit none
      Integer, Intent(in) :: nu
      Integer :: mud=0

! ... read parameters from file if any

      Call Read_ipar(nu,'klsp'  ,klsp  )
      Call Read_ipar(nu,'klsp1' ,klsp1 )
      Call Read_ipar(nu,'klsp2' ,klsp2 )
      Call Read_ipar(nu,'mk'    ,mk    )
      Call Read_ipar(nu,'mb'    ,mb    )
      Call Read_ipar(nu,'nb'    ,nb    )
      Call Read_ipar(nu,'kb'    ,kb    )
      Call Read_ipar(nu,'mrel'  ,mrel  )
      Call Read_ipar(nu,'mso'   ,mso   )
      Call Read_ipar(nu,'msoo'  ,msoo  )
      Call Read_ipar(nu,'mss '  ,mss   )
      Call Read_ipar(nu,'moo '  ,moo   )
      Call Read_ipar(nu,'imvc'  ,imvc  )
      Call Read_ipar(nu,'nmvc'  ,nmvc  )
      Call Read_ipar(nu,'maxnc' ,maxnc )
      Call Read_ipar(nu,'nud '  ,mud   )
      Call Read_ipar(nu,'iitar' ,iitar )
      Call Read_ipar(nu,'izcorr',izcorr)
      Call Read_rpar(nu,'zcorr' ,zcorr )
      Call Read_ipar(nu,'debug' ,debug )

      Call Read_rpar(nu,'s_ovl' ,s_ovl )
      Call Read_rpar(nu,'s_pert',s_pert)
      Call Read_rpar(nu,'eps_c' ,eps_c )
      Call Read_rpar(nu,'eps_det',eps_det)
      Call Read_rpar(nu,'eps_soo',eps_soo)
      Call Read_rpar(nu,'eps_tar',eps_tar)
      Call Read_rpar(nu,'eps_acf',eps_acf)

! ... overwrite parameters from arguments if any

      Call Read_iarg('klsp'  ,klsp  )
      Call Read_iarg('klsp1' ,klsp1 )
      Call Read_iarg('klsp2' ,klsp2 )
      Call Read_iarg('mk'    ,mk    )
      Call Read_iarg('mb'    ,mb    )
      Call Read_iarg('nb'    ,nb    )
      Call Read_iarg('kb'    ,kb    )
      Call Read_iarg('mrel'  ,mrel  )
      Call Read_iarg('mso'   ,mso   )
      Call Read_iarg('msoo'  ,msoo  )
      Call Read_iarg('mss '  ,mss   )
      Call Read_iarg('moo '  ,moo   )
      Call Read_iarg('imvc'  ,imvc  )
      Call Read_iarg('nmvc'  ,nmvc  )
      Call Read_iarg('maxnc' ,maxnc )
      Call Read_iarg('nud '  ,mud   )
      Call Read_iarg('iitar' ,iitar )
      Call Read_iarg('izcorr',izcorr)
      Call Read_rarg('zcorr' ,zcorr )
      Call Read_iarg('debug' ,debug )

      Call Read_rarg('s_ovl' ,s_ovl )
      Call Read_rarg('s_pert',s_pert)
      Call Read_rarg('eps_c' ,eps_c )
      Call Read_rarg('eps_det',eps_det)
      Call Read_rarg('eps_soo',eps_soo)
      Call Read_rarg('eps_tar',eps_tar)
      Call Read_rarg('eps_acf',eps_acf)

      if(klsp.gt.0) then; klsp1=klsp; klsp2=klsp; end if
      if(klsp2.lt.klsp1) klsp2=klsp1
      if(mud.eq.0) nud=0

! ... print the used parameters:

      write(prj,'(a)') 'BSR_MAT parameters:'
      write(prj,*) 
      write(prj,'(a,i4,a)') 'klsp1   =',klsp1, &
       ' - fist partial wave under consideration'
      write(prj,'(a,i4,a)') 'klsp2   =',klsp2, &
       ' - last partial wave under consideration'
      write(prj,*) 
      if(iitar.eq.0) write(prj,'(a,i4,a/T14,a)') 'iitar   =',iitar, &
       ' - target states are supposed to be orthogonal', &
       '   eigenfunctions of target Hamiltonian'
      if(iitar.eq.1) write(prj,'(a,i4,a/T14,a)') 'iitar   =',iitar, &
       ' - target states are supposed to be orthogonal', &
       '   but may not diagonalize the target Hamiltonian'
      if(iitar.eq.2) write(prj,'(a,i4,a/T14,a)') 'iitar   =',iitar, &
       ' - target states may be non-orthogonal', &
       '   and may not diagonalize the target Hamiltonian'
      write(prj,*) 
      write(prj,'(a,i4,a)') 'mk      =',mk, &
       ' - maximum multipole index '
      write(prj,*) 
      write(prj,'(a,d11.2,a)') 'eps_c   =',eps_c, &
       '  -  tolerance for integral coefficients'
      write(prj,'(a,d11.2,a)') 'eps_det =',eps_det, &
       '  -  tolerance for total overlap factors'
      write(prj,'(a,d11.2,a)') 'eps_ovl =',eps_ovl, &
       '  -  tolerance for one-elecytron overlaps'
      write(prj,*) 
      write(prj,'(a,d11.2,a)') 's_ovl   =',s_ovl, &
       '  -  channel overlap limit for additional orthogonality conditions'
      write(prj,'(a,d11.2,a)') 's_pert  =',s_pert, &
       '  -  perturber overlap limit for additional orthogonality conditions'

! ... relativistic corrections:

      if(mrel.ge.2.and.mso .eq.0) mso =1
      if(mrel.ge.3.and.msoo.eq.0) msoo=1
      if(mrel.ge.4.and.mss .eq.0) mss =1
      if(mrel.ge.5.and.moo .eq.0) moo =1

! ... move rel.parameters to B-spline modules:

      irel=mrel; if(irel.gt.0) rel=.TRUE. 
      if(moo.gt.0) ioo = 1      

      write(prj,*) 
      if(mrel.gt.0) then
       write(prj,'(a,i4,a)') 'mrel    =',mrel, & 
        ' - relativistic corrections is included'
       if(mso.gt.0) &
       write(prj,'(a,i4,a)') 'mso     =',mso, &
        ' - spin-orbit interaction is included'
       if(msoo.gt.0) &
       write(prj,'(a,i4,a)') 'msoo    =',msoo, &
        ' - spin-other-orbit interaction is included'
       if(mss.gt.0) &
       write(prj,'(a,i4,a)') 'mss     =',mss, &
        ' - spin-spin interaction is included'
       if(msoo.gt.0) &
       write(prj,'(a,i4,a)') 'msoo    =',moo, &
        ' - orbit-orbit interaction is included'
       write(prj,'(a,i4,a)') 'imvc    =',imvc, &
        ' - mode for inclusion of mass-velocity term'
       if(izcorr.eq.0.and.nz.gt.40) izcorr = 1
       write(prj,'(a,i4,a)') 'izcorr  =',izcorr, &
        ' -  small-r cut-off correction for spin-orbit interaction'
       write(prj,'(a,i4,a)') 'mlso    =',mlso, &
        ' -  max. l for spin-orbit interaction'
       write(prj,*) 
       write(prj,'(a,f11.6,a)') 'eps_soo =',eps_soo, &
        '  -  tolerance for spin-orbit coefficients'
       write(prj,'(a,f11.6,a)') 'zcorr   =',zcorr, &
        '  -  semiempirical correction for spin-orbit paramters with l=1'
      else
       write(prj,'(a,i4,a)') 'mrel    =',mrel, & 
        ' - relativistic corrections is not included'
      end if

! ... other parameters:

      write(prj,*) 
      write(prj,'(a,i7,a)') 'nb      =',nb, &
       ' - number of blocks in module cmdata'
      write(prj,'(a,i7,a)') 'kb      =',kb, &
       ' - max. number of blocks for given type of integrals'
      write(prj,'(a,i7,a,f8.2,a)') 'mb      =',mb, &
       ' - size of block  => ', &
       28*mb*nb/(1024d0*1024d0),' Mb'
      write(prj,*) 
      write(prj,'(a,d11.2,a)') 'eps_acf =',eps_acf, &
       '  -  tolerance for asympt. coefficients'
      write(prj,'(a,d11.2,a)') 'eps_tar =',eps_tar, &
       '  -  tolerance for target energies and overlaps'
      write(prj,*) 
      write(prj,'(a,i4,a)') 'ma      =',ma, &
       ' - limit for file-name length'
      write(prj,*) 
      write(prj,'(a,i4,a)') 'debug   =',debug, &
       ' - debug level'
      write(prj,*) 
      write(prj,'(a,i10,a,f8.2,a)') 'maxnc   =',maxnc, &
       ' - size of buffer for coefficients => ', &
       20*maxnc/(1024d0*1024d0),' Mb'

      End Subroutine Read_arg


