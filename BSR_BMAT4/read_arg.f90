!======================================================================
      Subroutine Read_arg 
!======================================================================
!     read arguments from command line and check default settings
!======================================================================
      Use bsr_breit
      Use bsr_mat

      Implicit none
      Character(7) :: oper='1110000'	

! ... read arguments in command line:

      Call Read_aarg('oper'  ,oper  )
      Call Read_iarg('mk'    ,mk    )
      Call Read_iarg('debug' ,debug )
      Call Read_iarg('mkt'   ,mkt   )
      Call Read_iarg('mkdt'  ,mkdt  )

! ... define the operators under consideration:

      read(oper,'(7i1)') ioper

      write(pri,'(/a/)') 'Operators included: '

      if(ioper(1).gt.0) write(pri,'(a)') 'OVERLAPS'
      if(ioper(2).gt.0) write(pri,'(a)') 'KINATIC ENERGY'
      if(ioper(3).gt.0) write(pri,'(a)') 'TWO-ELECTRON ELECTROSTATIC'

      if(ioper(4).gt.0) write(pri,'(a)') 'SPIN ORBIT'
      if(ioper(5).gt.0) write(pri,'(a)') 'SPIN-OTHER-ORBIT'
      if(ioper(6).gt.0) write(pri,'(a)') 'SPIN-SPIN'
      if(ioper(7).gt.0) write(pri,'(a)') 'ORBIT-ORBIT'

      write(pri,'(/a,i3)') 'Max.multipole index =',mk

      klsp = 0
      Call Read_iarg('klsp'  ,klsp  )
      Call Read_iarg('klsp1' ,klsp1 )
      Call Read_iarg('klsp2' ,klsp2 )
      if(klsp.gt.0) then;  klsp1=klsp; klsp2=klsp; end if 
!-----------------------------------------------------------------
! ... bsr_mat parameters:

      Call Read_iarg('mb',mb)
      Call Read_iarg('nb',nb)
      Call Read_iarg('kb',kb)

      Call Read_rarg('eps_c' , eps_c )
      Call Read_rarg('eps_det',eps_det)
      Call Read_rarg('eps_ovl',eps_ovl)
      Call Read_rarg('eps_soo',eps_soo)
      Call Read_iarg('mlso',   mlso)

!      Call Read_iarg('check_target',check_target)
      Call Read_rarg('time',time_limit)
      Call Read_rarg('time_limit',time_limit)

! ... print the used parameters:

      write(pri,*) 
      write(pri,'(a,1PE11.1,a)') 'eps_c   =',eps_c, &
       '  -  tolerance for integral coefficients'
      write(pri,'(a,1PE11.1,a)') 'eps_det =',eps_det, &
       '  -  tolerance for total overlap factors'
      write(pri,'(a,1PE11.1,a)') 'eps_ovl =',eps_ovl, &
       '  -  tolerance for one-elecytron overlaps'
      if(ioper(4)+ioper(5)+ioper(6)+ioper(7).gt.0) &
      write(pri,'(a,1PE11.1,a)') 'eps_soo =',eps_soo, &
       '  -  tolerance for spin-orbit coefficients'
      if(ioper(4).gt.0) &
      write(pri,'(a,i4,a)') 'mlso    =',mlso, &
       ' -  max. l for spin-orbit interaction'

      write(pri,*) 
      write(pri,'(a,i7,a)')        'nb      =',nb, &
       ' - number of blocks in module c_blocks'
      write(pri,'(a,i7,a)')        'kb      =',kb, &
       ' - max. number of blocks for given type of integrals'
      write(pri,'(a,i7,a,f8.2,a)') 'mb      =',mb, &
       ' - size of block  => ', &
       28*mb*nb/(1024d0*1024d0),' Mb, overall'

      if(debug.gt.0) &
      write(pri,'(/a,i4,a)') 'debug   =',debug, &
       ' - debug level'

!      if(check_target.eq.0) &
!      write(pri,'(/a,i2,a)') 'check_target =',check_target, &
!       ' - target coefficients are not included'
!      if(check_target.eq.1) &
!      write(pri,'(/a,i2,a)') 'check_target =',check_target, &
!       ' - target coefficients are included'

      if(time_limit.gt.0.d0) &  
      write(pri,'(/a,f8.1,a)') 'time_limit   =',time_limit, &
       ' - time_limit for calculations (min)'

      End Subroutine Read_arg
