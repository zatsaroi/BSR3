!======================================================================
      Subroutine Read_arg
!======================================================================
!     read input parameters
!----------------------------------------------------------------------
      Use bsr_prep

      Implicit none 

      open(nup,file=AF_par)

      Call Read_rpar(nup,'eps_ovl' ,eps_ovl )
      Call Read_rarg(    'eps_ovl' ,eps_ovl )
      Call Read_rpar(nup,'eps_phys',eps_phys)
      Call Read_rarg(    'eps_phys',eps_phys)
      Call Read_rpar(nup,'eps_sub' ,eps_sub )
      Call Read_rarg(    'eps_sub' ,eps_sub )
      Call Read_rpar(nup,'eps_targ',eps_targ)
      Call Read_rarg(    'eps_targ',eps_targ)
      Call Read_ipar(nup,'ii_sub'  ,ii_sub  )
      Call Read_iarg(    'ii_sub'  ,ii_sub  )

      Call Read_ipar(nup,'LT_min',LT_min)
      Call Read_iarg(    'LT_min',LT_min)
      Call Read_ipar(nup,'LT_max',LT_max)
      Call Read_iarg(    'LT_max',LT_max)
      Call Read_ipar(nup,'IS_min',IS_min)
      Call Read_iarg(    'IS_min',IS_min)
      Call Read_ipar(nup,'IS_max',IS_max)
      Call Read_iarg(    'IS_max',IS_max)
      Call Read_ipar(nup,'JJ_min',JJ_min)
      Call Read_iarg(    'JJ_min',JJ_min)
      Call Read_ipar(nup,'JJ_max',JJ_max)
      Call Read_iarg(    'JJ_max',JJ_max)

      write(pri,'(a/)') 'BSR_PREP:'

      write(pri,'(a,d15.5,a/)') 'eps_ovl  = ',eps_ovl, &
       ' - tolerance for orthogonality between orbitals'
      write(pri,'(a,d15.5,a/)') 'eps_sub  = ',eps_sub, &
       ' - tolerance for substitution orbitals'
      write(pri,'(a,d15.5,a/)') 'eps_phys = ',eps_phys,& 
       ' - tolerance for physical orbitals'
      write(pri,'(a,i15  ,a/)') 'ii_sub   = ',ii_sub,  & 
       ' - restriction for substitution orbitals'
      write(pri,'(a,d15.5,a/)') 'eps_targ = ',eps_targ,& 
       ' - restriction to target expansions'

      write(pri,'(/72(''-'')/)')


      End Subroutine Read_arg
