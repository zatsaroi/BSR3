!=====================================================================
      Subroutine Read_data
!=====================================================================
!     read radial functions and prepare work arrays
!---------------------------------------------------------------------
      Use bsr_ci

      Implicit none

      Integer :: i,j
      Real(8) :: S
      Real(8), external :: QUADR

!----------------------------------------------------------------------
! ... read B-spline expantions for bound orbitals:
   
      Call alloc_spline_orbitals(nwf)
      nbf = nwf
      Do i = 1,nbf
       ebs(i)=ELF(i); nbs(i)=NEF(i); lbs(i)=LEF(i); kbs(i)=KEF(i); mbs(i) = 0
      End do

      Call Check_file(AF_w)
      Open(nuw,file=AF_w,form='UNFORMATTED')
      Call Read_bsw(nuw)
      Close(nuw)

! ... check the correspondence between c- and bsw-files: 

      j = 0
      Do i = 1,nwf
       if(mbs(i).eq.0) then
        write(iwrite,'(a,a)') ' Absent expansion for w.f. ',ELF(i)
        j = j + 1
       end if
      End do
      if(j.gt.0) Stop 'no correspondence between c- and w- files'
      
! ... the < p | p > values 

      Allocate(OBS(nbf,nbf)); OBS = 0.d0
      Do i=1,nbf
       Do j=1,i
        if(lbs(i).ne.lbs(j)) Cycle
        S=QUADR(i,j,0); if(abs(S).lt.Eps_ovl) S=0.d0
        if(i.eq.j) S=1.d0
        OBS(i,j)=S; OBS(j,i)=S
       End do
      End do

! ... core energy: 

      Call Bcore

      write(iwrite,'(/a,i4,a)') 'nclosd  =',nclosd, ' - common core shells'
      write(iwrite,'(/a,F15.8,a)') 'Ecore   =', EC, '  -  core energy'

      End Subroutine Read_data      