!=====================================================================
      Subroutine Read_arg
!=====================================================================
!     handle the input data and defalt values, read knot.dat 
!---------------------------------------------------------------------
      Use bsr_ci

      Implicit none

      Integer :: an
      Integer, external :: Icheck_file
      Character(80) :: aa1, aa2

!----------------------------------------------------------------------
! ... check the name of the case:

      Call Read_name(name)
      if(name.eq.'?'.or.len_trim(name).eq.0) Call inform

      AF_c=trim(name)//'.c'
      AF_w=trim(name)//'.bsw'
      AF_d=trim(name)//'.d'
      AF_l=trim(name)//'.l'
      AF_j=trim(name)//'.j'
      AF_b=trim(name)//'.bnk'
      AF_k=trim(name)//'.knot'

      AF_log=trim(name)//'.log'
      Open(iwrite,file=AF_log)
      write(iwrite,'(a)') 'BSR_CI calculations:'
      write(iwrite,'(/a,a)')  'CASE:  ', trim(name)

!----------------------------------------------------------------------
! ... read parameters from file if any:

      AF_inp = trim(name)//'.inp'
      if(icheck_file(AF_inp).ne.0) then
       Open(iread,file=AF_inp) 
       Call Read_apar(iread,'name'  ,name  )
       Call Read_ipar(iread,'irel'  ,irel  )
       Call Read_ipar(iread,'mrel'  ,irel  )
       Call Read_ipar(iread,'jmax'  ,jmax  )
       Call Read_ipar(iread,'jmin'  ,jmin  )
       Call Read_ipar(iread,'meiv'  ,meiv  )
       Call Read_ipar(iread,'mass'  ,mass  )
       Call Read_ipar(iread,'iso'   ,iso   )
       Call Read_ipar(iread,'isoo'  ,isoo  )
       Call Read_ipar(iread,'iss '  ,iss   )
       Call Read_ipar(iread,'ioo '  ,ioo   )
       Call Read_rpar(iread,'eps_ovl',eps_ovl)
       Call Read_rpar(iread,'eps_o' ,eps_o )
       Call Read_rpar(iread,'eps_d' ,eps_d )
       Call Read_rpar(iread,'Emax'  ,Emax  )
       Call Read_ipar(iread,'mdiag' ,mdiag )
       Close(iread)
      end if

!----------------------------------------------------------------------
! ... overwrite parameters from arguments if any:

      Call Read_iarg('meiv'  ,meiv  )
      Call Read_iarg('msol'  ,meiv  )
      Call Read_iarg('irel'  ,irel  )
      Call Read_iarg('mrel'  ,irel  )
      Call Read_iarg('jmax'  ,jmax  )
      Call Read_iarg('jmin'  ,jmin  )
      Call Read_iarg('iso'   ,iso   )
      Call Read_iarg('isoo'  ,isoo  )
      Call Read_iarg('iss '  ,iss   )
      Call Read_iarg('ioo '  ,ioo   )
      Call Read_rarg('eps_ovl',eps_ovl)
      Call Read_rarg('eps_o' ,eps_o )
      Call Read_rarg('eps_d' ,eps_d )
      Call Read_rarg('zcorr' ,zcorr )
      Call Read_rarg('Emax'  ,Emax  )
      Call Read_iarg('mdiag' ,mdiag )
      Call Read_iarg('pri_mat',pri_mat)

      jmin = jmin+1    !  to (2j+1) form         ! ???
      jmax = jmax+1

!----------------------------------------------------------------------
!                                                  open relevant files:
! ... c-file:

      Call Check_file(AF_c)
      Open(nuc,file=AF_c,status='OLD')

! ... int_bnk:

      if(Icheck_file(AF_b).ne.0) then
       Open(nub,file=AF_b,form='UNFORMATTED')        ! name.bnk
      else
       Call Check_file(BF_b)
       Open(nub,file=BF_b,form='UNFORMATTED')        ! int_bnk
      end if

! ... knot.dat

      if(Icheck_file(AF_k).ne.0) then
       AF_grid = trim(AF_k)
      else
       Call Check_file(BF_k)
       AF_grid = trim(BF_k)
      end if
      Call Read_grid(Z)
      Close(nuk)

      Call Define_spline

      an = NINT(Z)
      Call Def_atom_LS(an,atom,aa1,aa2)


      write(iwrite,'(/a,a)') 'ATOM = ', atom
      write(iwrite,'(/a,f7.2)') 'NUCLEAR CHARGE = ',Z

!----------------------------------------------------------------------

      write(iwrite,'(/a)') 'Parameters of calculations:'

      if(irel.gt.0)  then
       REL = .TRUE.
       if(irel.ge.2.and.iso .ne.-1) iso =1; if(iso .eq.-1) iso =0
       if(irel.ge.3.and.isoo.ne.-1) isoo=1; if(isoo.eq.-1) isoo=0
       if(irel.ge.4.and.iss .ne.-1) iss =1; if(iss .eq.-1) iss =0
       if(irel.ge.5.and.ioo .ne.-1) ioo =1; if(ioo .eq.-1) ioo =0
      end if

      if(irel+iso+isoo+iss+ioo.eq.0) then
       write(iwrite,'(/a)') 'non-relativistic calculations'
      else

       write(iwrite,'(/a)') 'relativistic corrections are included for:'

       if(REL) write(iwrite,'(/a)') 'relativistic one-electron shift'
       if(iso.gt.0) write(iwrite,'(a)') 'spin-orbit interaction'
       if(isoo.gt.0)write(iwrite,'(a)') 'spin-other-orbit interaction'
       if(iss.gt.0) write(iwrite,'(a)') 'spin-spin interaction'
       if(ioo.gt.0) write(iwrite,'(a)') 'orbit-orbit interaction'

       write(iwrite,'(/a,i6)') 'max. 2J+1 = ', jmax
       write(iwrite,'( a,i6)') 'min. 2J+1 = ', jmin

      end if

!----------------------------------------------------------------------
!
      if(meiv.gt.0) &
      write(iwrite,'(/a,i6)') 'Number of required eigenvalues: ', meiv
      if(meiv.eq.0) &
      write(iwrite,'(/a,i6)') 'Number of required eigenvalues:  ALL '

      End Subroutine Read_arg

