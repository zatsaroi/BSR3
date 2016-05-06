!=====================================================================
      Subroutine SUB1_HD    
!=====================================================================
!     calculations for given partial wave
!---------------------------------------------------------------------
      Use bsr_hd
      Use target
      Use channel
      Use spline_param, only: ns
      
      Implicit none
      Integer :: i,j, i1,i2,i3
      Integer, External :: Icheck_file      

      Call R_channel(nut,klsp)
      if(ncp.gt.0) ippert = ippert + ipconf(nch)
!----------------------------------------------------------------------
! ... open log file and print main parameters:

      write(ALSP,'(i3.3)') klsp
      i = INDEX(AF_log,'.'); AF=AF_log(1:i)//ALSP
      Open(pri,file=AF)
      
      write(pri,'(a/a)') 'B S R _ H D ', &
                         '***********'
      write(pri,*)
      write(pri,'(a,i3)') 'calculations for partial wave:  klsp =',klsp
	  
      write(pri,*)
	     if(itype.eq.-1) &
       write(pri,'(a)') 'itype =   -1  -  bound-state calculations'
      if(itype.eq. 0) &
       write(pri,'(a)') 'itype =    0  -  scattering calculations'
      if(itype.gt. 0) &
       write(pri,'(a)') 'itype =    1  -  photoionization calculations'

      kch = nch; kcp = npert; nhm = kch*ns + kcp

      write(pri,*)
      write(pri,'(a,i5,a)') 'kch   =',kch,'  -  number of channels'
      write(pri,'(a,i5,a)') 'kcp   =',kcp,'  -  number of pertubers'
      write(pri,'(a,i5,a)') 'nhm   =',nhm,'  -  full size of matrix'

  
      write(pri,*)
      if(iexp.eq.0) &
      write(pri,'(a)') 'iexp  =    0  -  theoretical target energies'
      if(iexp.ne.0) &
      write(pri,'(a)') 'iexp  =    1  -  exp.target energies'

      if(iexp.gt.0) then
      if(iiexp.eq.0) &
      write(pri,'(a)') 'iiexp  =    0  -  order is not changed'
      if(iiexp.ne.0) &
      write(pri,'(a)') 'iiexp  =    1  -  order is changed'

       write(pri,'(/a/)') 'target  ip_exp  Etarg   E_exp    Etarg-E_exp'

       Do j=1,ntarg; i=ip_exp(j)
        write(pri,'(2i6,2f16.8,f10.5)') &
        j, i, Etarg(j), E_exp(j), (E_exp(j)-Etarg(j))*27.2112   
       End do
      end if

      write(pri,*)
      write(pri,'(a,i5,a)') 'jmvc  =',jmvc, &
       '  -  number of channel orbitals for inclusion of' 
      write(pri,'(17x,a)') 'mass-velocity corrections' 
             
      write(pri,*)
      write(pri,'(a)') 'Restrictions on interaction matrix:'
      write(pri,*)
      write(pri,'(a,f10.5)') 'Emin  =',Emin
      write(pri,'(a,f10.5)') 'Emax  =',Emax
      write(pri,'(a,f10.5)') 'Egap  =',Egap
      write(pri,*)
      write(pri,'(a)') 'All one-channel solutions with &
         &E<Emin, E>Emax or abs(E)<Egap will be deleted.' 

      if(itype.eq.-1) then
       write(pri,*)
       write(pri,'(a)') 'Restrictions for output of bound states:'
       write(pri,*)
       write(pri,'(a,i5,a)') 'msol   =',msol,  '  - max. number of solutions' 
       write(pri,'(a,i5,a)') 'it_max =',it_max,'  - max. target threshold'
       write(pri,*)
       write(pri,'(a)') 'Zero value of any parameter means no restrictions.'
      end if

!----------------------------------------------------------------------
! ... check the file with interaction matrix:

      i = INDEX(AF_int,'.'); AF=AF_int(1:i)//ALSP
      i = Icheck_file(AF)
      if(i.eq.0) then
       write(pri,*) 'there is no bsr_mat.nnn file for given partial wave'
       write(*  ,*) 'there is no bsr_mat.nnn file for given partial wave'
       Return
      end if

      Open(nui,file=AF,status='OLD',form='UNFORMATTED')

! ... check the dimensions:

      read(nui) i1,i2,i3
      if(i1.ne.ns ) Stop ' BSR_HD: different ns  in BSR_MAT file'
      if(i2.ne.kch) Stop ' BSR_HD: different kch in BSR_MAT file'
      if(i3.ne.kcp) Stop ' BSR_HD: different kcp in BSR_MAT file'
!----------------------------------------------------------------------
! ... diagonalize the matrix and get inner-region solutions:

      Call Diag_hd
 
!----------------------------------------------------------------------
! ... output of solutions and find the surface amplitudes:

      if(itype.ge.0)  Call rsol_out

!----------------------------------------------------------------------
! ... output of standard H.nnn file:

      if(itype.ge.0) then

      ! read max. mutipole order and asymptotic coef.s:

      read(nui) lamax
      if(allocated(CF)) Deallocate(CF)
      Allocate (CF(kch,kch,lamax+1))
      read(nui) CF
      read(nui) RA     !  RM radius:

       if(itype.ge.0)  then
        if(iiexp.eq.0) Call H_out
        if(iiexp.ne.0) Call H_out1
       end if

      end if
!----------------------------------------------------------------------
! ... output of bound states in bound.nnn:

      if(itype.ne.0) Call B_out 

      End Subroutine SUB1_HD



