!=======================================================================
!     merging the set of w-files with choice of orbitals and 
!     with optional changing the spectroscopic notation
! 
!        1.bsw + 2.bsww + 3.bsw + ... --> res.bsw 
!
!=======================================================================
      Use spline_param; Use spline_atomic; Use spline_grid
      Use spline_orbitals
     
      Implicit double precision (A-H,O-Z)
      Character(4) :: elw,eln
      Character(4), external :: ELF4
      Character(40) :: AF
      Integer :: nu=1

!----------------------------------------------------------------------
      iarg = COMMAND_ARGUMENT_COUNT()
      if(iarg.gt.0) Call GET_COMMAND_ARGUMENT(1,AF)

      if(AF.eq.'?') then
        write(*,'(/a)') 'bsw123 merges the set of bsw-files with optional choice of orbitals'
        write(*,'(/a)') '     1.bsw + 2.bsw + 3.bsw + ... --> res.bsw'
        write(*,'(/a)') 'program has interactive input/output' 
        write(*,'(/a)') 'Call as:  bsw123  or  bsw123  name.bsw'
        write(*,'(/a)') 
        Stop 
      end if        

      if(iarg.lt.1) then
       write(*,*) 'Enter file-name for bsw-file: '
       read(*,'(a)') AF
      end if      

      Call OpenF(nu,AF,'UNFORMATTED','OLD')
      read(nu) elw,z,h,hmax,rmax,ks,ns
      rewind(nu)

! ... sets up grid points and initializes the values of the spline: 
    
      CALL define_grid(zz);  CALL define_spline

      nbf=0;  CALL Allocate_bsorb(ibf)

      i=1; go to 1
!----------------------------------------------------------------------

   10 write(*,*) 'Enter file-name for bsw-file  or  end: '
      read(*,'(a)') AF
      ii=LEN_trim(AF)
      if(AF(1:3).eq.'end') go to 3
      if(ii.eq.0) go to 10
      Call OpenF(nu,AF,'UNFORMATTED','OLD')

    1 READ(nu,end=2) ebs(i),zw,hw,hmw,rmw,ksw,nsw,mbs(i)

      read(nu) pbs(1:mbs(i),i)

      if(zw.ne.z.or.hw.ne.h.or.abs(hmw-hmax).gt.1.d-15.or.rmw.ne.rmax &
         .or.ksw.ne.ks.or.nsw.ne.ns) then
         write(*,*)  ' file  ',AF,' has another B_spline parameters'
         go to 10
      end if

      write(*,'(a6)') ebs(i);  write(*,'(a)') ' new EL ?  (d|EL,a4): '
      read(*,'(a)') ELN 

      if(ELN(1:1).eq.'d') go to 1
      ii=LEN_TRIM(ELN)
      if(ii.gt.1) then
       Call EL4_nlk(ELN,n1,l1,k1); ebs(i)=ELF4(n1,l1,k1)
      end if
 
      Do j=1,i-1
       if(ebs(j).eq.ebs(i)) then; i=i-1; Exit; end if
      End do
      i=i+1
      if(i.gt.mbf)  Call Allocate_bsorb(mbf+jbf)
      go to 1

    2 Close(nu)
      go to 10

    3 nbf=i-1
      write(*,'(a,i2,a)')  'Now there is ', nbf, '  w.f.:'
      write(*,'(12(1x,a4))') (ebs(i),i=1,nbf)
      write(*,*) 'Enter file-name for result w-file : '
      read(*,'(a)') AF
      Call OpenF(nu,AF,'UNFORMATTED','UNKNOWN')

      DO i=1,nbf
       write(nu) ebs(i),zw,hw,hmw,rmw,ksw,nsw,mbs(i)
       write(nu) pbs(1:mbs(i),i)
      End do
      Close(nu)

      End  ! utility bsw123


!=======================================================================
      Subroutine OpenF (nu,AF,AFOR,AST)
!=======================================================================
!     safe opening of the file
!-----------------------------------------------------------------------
      Character AF*(*),AFOR*(*),AST*(*)

    1 OPEN(nu,file=AF,form=AFOR,status=AST,err=2)
      Return

    2 ii=LEN_trim(AF)
      write(*,*) 'There is no file  ',AF(1:ii)
      write(*,*) 'Try again, enter file-name or end: '
      read(*,*) AF
      ii=LEN_trim(AF)
      if(AF(1:3).eq.'end') Stop ' '
      go to 1

      End Subroutine OpenF


