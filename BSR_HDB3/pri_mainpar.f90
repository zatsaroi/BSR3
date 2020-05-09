!=====================================================================
      Subroutine pri_mainpar
!=====================================================================
!     check and print main parameters
!---------------------------------------------------------------------

      Use bsr_hd
      Use target,       only: etarg
      Use channel,      only: nch,npert
      Use spline_param, only: ns
      
      Implicit none

      Integer :: i, i1,i2,i3
      Integer, external :: Icheck_file      

!----------------------------------------------------------------------
! ... open log file and print main parameters:

      write(ALSP,'(i3.3)') klsp
      i = LEN_TRIM(AF_nnn); AF=AF_nnn(1:i-3)//ALSP
      Open(pri,file=AF)
      
      write(pri,'(a/a)') 'B S R _ H D ', &
                         '***********'
      write(pri,*)
      write(pri,'(a,i3)') 'Calculations for partial wave:  klsp =',klsp
      write(*  ,*)
      write(*  ,'(a,i3)') 'Calculations for partial wave:  klsp =',klsp
	  
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

      write(pri,*)
      write(pri,'(a,i5,a)') 'jmvc  =',jmvc, &
       '  -  number of channel orbitals for inclusion of' 
      write(pri,'(17x,a)') 'mass-velocity corrections' 
             
      write(pri,*)
      write(pri,'(a)') 'Restrictions on interaction matrix:'
      write(pri,*)
      write(pri,'(a,f10.5)') 'Edmin =',Edmin
      write(pri,'(a,f10.5)') 'Edmax =',Edmax
      write(pri,'(a,f10.5)') 'Egap  =',Egap
      write(pri,*)
      write(pri,'(a)') 'All one-channel solutions with &
         &E<Emin, E>Emax or abs(E)<Egap will be deleted.' 

      if(itype.eq.-1) then
       write(pri,*)
       write(pri,'(a)') 'Restrictions for output of bound states:'
       write(pri,*)
       write(pri,'(a,i5,a)') 'msol   =',msol,  '  - max. number of solutions' 
       write(pri,'(a,i5,a)') 'IT_max =',IT_max,'  - max. target threshold'
       write(pri,'(a,f10.5)') 'Emin   =',Emin
       write(pri,'(a,f10.5)') 'Emax   =',Emax
       write(pri,*)
       write(pri,'(a)') 'Zero value of any parameter means no restrictions.'
      end if
      if(it_max.gt.0) then
       if(Emax.gt.Etarg(it_max)) Emax=Etarg(it_max)
      end if
!----------------------------------------------------------------------
! ... check the file with interaction matrix:

      i = LEN_TRIM(AF_int); AF=AF_int(1:i-3)//ALSP
      i = Icheck_file(AF)
      if(i.eq.0) then
       write(pri,*) 'there is no bsr_mat.nnn file for given partial wave'
       write(*  ,*) 'there is no bsr_mat.nnn file for given partial wave'
       fail=1; Return
      end if

      Open(nui,file=AF,status='OLD',form='UNFORMATTED')

! ... check the dimensions:

      read(nui) i1,i2,i3
      if(i1.ne.ns ) write(*,*) 'BSR_HD: different ns  in BSR_MAT file'
      if(i2.ne.kch) write(*,*) 'BSR_HD: different kch in BSR_MAT file'
      if(i3.ne.kcp) write(*,*) 'BSR_HD: different kcp in BSR_MAT file'
      if(i1.ne.ns.or.i2.ne.kch.or.i3.ne.kcp) then
       fail=1; Return
      end if  

! ... allocate common arrays:

      if(allocated(bb)   ) deallocate(bb   ); allocate(bb(ns,ns*kch)); bb = zero
      if(allocated(bval) ) deallocate(bval ); allocate(bval (ns*kch)); bval = zero
      if(allocated(ipsol)) deallocate(ipsol); allocate(ipsol(0:kch) ); ipsol = 0
      if(allocated(cc)   ) deallocate(cc   ); allocate(cc   (ns,ns) ); cc = zero
      if(allocated(w)    ) deallocate(w    ); allocate(w    (ns)    ); w = zero

      End Subroutine pri_mainpar
