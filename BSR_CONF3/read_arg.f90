!======================================================================
      Subroutine Read_arg
!======================================================================
! ... read parameters for given case
!----------------------------------------------------------------------
      USE bsr_conf
      USE target, only: coupling
      Implicit none

      Open(pri,file=AF_log)
      Call Check_file(AF_par); Open(nup,file=AF_par)

      Call Read_rpar(nup,'c_comp',c_comp)
      Call Read_rarg(    'c_comp',c_comp)
      Call Read_ipar(nup,'kort'  ,kort  )
      Call Read_iarg(    'kort'  ,kort  )
      Call Read_ipar(nup,'max_ll',max_ll)
      Call Read_iarg(    'max_ll',max_ll)
      Call Read_ipar(nup,'min_ll',min_ll)
      Call Read_iarg(    'min_ll',min_ll)
      Call Read_ipar(nup,'max_LT',max_LT)
      Call Read_iarg(    'max_LT',max_LT)
      Call Read_ipar(nup,'max_ST',max_ST)
      Call Read_iarg(    'max_ST',max_ST)
      Call Read_ipar(nup,'max_it',max_it)
      Call Read_iarg(    'max_it',max_it)
      Call Read_ipar(nup,'debug' ,debug )
      Call Read_iarg(    'debug' ,debug )
      Call Read_ipar(nup,'iread_targ',iread_targ)
      Call Read_iarg(    'iread_targ',iread_targ)

      write(pri,'(a/  )')       'BSR_CONF parameters:'
      write(pri,'(a,a/)')       'coupling = ',COUPLING
      write(pri,'(a,f5.3,a/)')  'c_comp   = ',c_comp, &
         ' - tolerance for compensation configurations'
      write(pri,'(a,i5,a/)')    'max_ll   = ',max_ll, &
         ' - restiction on max. small l'
      write(pri,'(a,i5,a/)')    'min_ll   = ',min_ll, &
         ' - restiction on min. small l'
      write(pri,'(a,i5,a/)')    'max_LT   = ',max_LT, &
         ' - restiction on total L as (2L+1)'
      write(pri,'(a,i5,a/)')    'max_ST   = ',max_ST, &
         ' - restiction on total S as (2S+1)'
      write(pri,'(a,i5,a/)')    'max_it   = ',max_it, &
         ' - restiction on number of target states '
      write(pri,'(a,i5,a/)')    'kort     = ',kort, &
         ' - if > 0, the orth.conditions from cfg-files are in play'
      write(pri,'(72(''-''))')

      End Subroutine Read_arg
