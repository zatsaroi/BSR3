!======================================================================
      Subroutine Read_data
!======================================================================
!     read data for given partial wave, klsp
!----------------------------------------------------------------------
      Use bsr_breit
      Use bsr_mat
      Use conf_LS;  Use orb_LS

      Implicit none
      Integer :: i,j,m, nc,lc,kc, ich,jch
      Integer, external :: Ifind_nlk, Iadd_line, Ifind_channel
      Character(124) :: line

! ... read physical orbitals if they not in c-file:
! ... (e.g., substitution orbitals)

      Call Check_file(AF_bsw)
      Open(nuw,file=AF_bsw,form='UNFORMATTED')
      Call Read_bsw_orb_LS(nuw)
      Close(nuw)

! ... read substitute orbitals pointers:

      Call Check_file(AF_orb)
      Open(nuo,file=AF_orb)
      Call Read_sub_orb_LS(nuo,ntarg)

! ... read channel information and find pointer (orbital --> channel):

      Call R_channel(nut,klsp)

      ief = 0
      Do j = 1,nch
       Call EL4_NLK(ELC(j),nc,lc,kc)
       i=Ifind_nlk(nc,lc,kc,1);  ief(i)=j; ipch(j)=i
      End do

      write(pri,'(/a/)')   'Partial wave LSP (JP): '
      if(jpar.lt.0) &
      write(pri,'(a,i10,a)')  'lpar   = ',lpar, '  -  total L '
      if(ispar.gt.0) &
      write(pri,'(a,i10,a)')  'ispar  = ',ispar,'  -  total 2S+1'
      write(pri,'(a,i10,a)')  'ipar   = ',ipar, '  -  parity'
      if(jpar.ge.0) &
      write(pri,'(a,i10,a)')  'jpar   = ',jpar-1, '  -  total 2J'

! ... read B-spline expantions for bound orbitals:
   
      Call Allocate_bsorb(nwf)
      nbf = nwf
      Do i = 1,nbf
       ebs(i)=ELF(i); nbs(i)=NEF(i); lbs(i)=LEF(i); kbs(i)=KEF(i)
       iech(i)=ief(i); mbs(i) = 0
      End do

! ... target radial functions:

      Open(nuw, file=AF_bsw, form='UNFORMATTED')
      Call Read_bsw(nuw)
      Close(nuw)

! ... perturber radial functions:

      if(nwp.gt.0) then
       i=LEN_TRIM(AFP); AF=AFP(1:i)//'.bsw'
       Call Check_file(AF)
       Open(nuw, file=AF, form='UNFORMATTED')
       Call Read_bsw(nuw)
       Close(nuw)
      end if
      
! ... correction of pertuber pointer:

      if(ncp.gt.0) ippert = ippert + ipconf(nch)

! ... check the correspondence between c- and bsw-files: 

      j = 0
      Do i = 1,nwf
       if(iech(i).ne.0) Cycle
       if(mbs(i).eq.0) then
        write(pri,'(a,a)') ' Absent expansion for w.f. ',ELF(i)
        j = j + 1
       end if
      End do
      if(j.gt.0) Call Stop_mpi(0,0,'no correspondence between c- and w- files')

! ... the < B | p > values (convolution with B matrix)

      Do i=1,nbf
       if(iech(i).ne.0) Cycle
       Call BXV(ks,ns,sb(1,1),pbs(1,i),qbs(1,i))
      End do

! ... read orthogonal conditions:

      Call alloc_file(0)
      rewind(nuc)
    1 read(nuc,'(a)',end=2) line
      if(line(1:1).ne.'<') go to 1 
      i = Iadd_line(line)
      go to 1
    2 Continue

! ... the < p | p > values 

      Call Alloc_orb_overlaps(nbf,lbs,iech,kclosd)

      End Subroutine Read_data
    
