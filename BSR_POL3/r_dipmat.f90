!======================================================================
      Subroutine r_dipmat
!======================================================================
!     read the data from dv.nnn file for the given partial wave
!----------------------------------------------------------------------
      Use bsr_pol
      Use channel, only: nch,ncp
      
      Implicit none
      Integer :: i, i1,i2,i3

      i=LEN_TRIM(AF_dip); AF_dip(i-2:i)=ALSP
      Call Check_file(AF_dip)
      Open(nud,file=AF_dip,form='UNFORMATTED')

      read(nud) kpol,ktype

      write(pri,*)
      write(pri,'(a,a5,i1)') 'transition:', ktype,kpol
      write(pri,*)

      read(nud) E0,jot0,IP0,Label0

      write(pri,*)
      write(pri,'(a,a)') 'initial state:', Label0
      write(pri,*)
      write(pri,'(a,e16.8)') 'energy  = ', E0
      write(pri,'(a,i5)')    '2J = ', jot0

      ! check the dimensions:

      read(nud) i1,i2,i3
      if(i1.ne.nhm) Stop ' BSR_POL: different nhm in BSR_MAT file'
      if(i2.ne.nch) Stop ' BSR_POL: different kch in BSR_MAT file'
      if(i3.ne.ncp) Stop ' BSR_POL: different kcp in BSR_MAT file'

      if(allocated(d)) Deallocate(d); Allocate(d(1:mhm))
      d = zero
      read(nud) (d(i),i=1,nhm)

      End Subroutine R_dipmat 
