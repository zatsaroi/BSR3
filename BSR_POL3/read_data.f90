!======================================================================
      Subroutine Read_data
!======================================================================
!     read data for given partial wave klsp
!----------------------------------------------------------------------
      Use bsr_pol
      Use channel
      Use orb_LS
      Use spline_param,      only: ns,ks
      Use spline_atomic,     only: z
      Use spline_galerkin,   only: sb
      Use spline_orbitals   

      Implicit none
      Integer :: i,j, nc,lc,kc, ich
      Integer, external :: Ifind_nlk

! ... set up B-splines:
 
      CALL define_grid(z); Call define_spline   

! ... read arguments:

      Call R_arg

! ... log-file:

      i=LEN_TRIM(AF_log); AF_log(i-2:i)=ALSP
      Open(pri,file=AF_log)

      write(pri,'(a,i5)')   'BSR_POL:   klsp = ',klsp
      write(pri,'(a)')      '*************************'

!----------------------------------------------------------------------
! ... find obitals:

! ... target orbitals:

      Call Check_file(AF_bsw)
      Open(nuw, file=AF_bsw, form='UNFORMATTED')
      Call Read_bsw_orb_LS(nuw)
      Close(nuw)

! ... perturber radial functions:

      if(nwp.gt.0) then
       i=Index(AFP,'.'); AFP(i+1:i+3)='bsw'
       Call Check_file(AFP)
       Open(nuw, file=AFP, form='UNFORMATTED')
       Call Read_bsw_orb_LS(nuw)
       Close(nuw)
      end if

! ... read c-file for channel orbitals

      i=LEN_TRIM(AF_cfg); AF_cfg(i-2:i)=ALSP
      Call Check_file(AF_cfg)
      Open(nuc,file=AF_cfg)

      Call R_closed(nuc)
      Call R_conf_LS(nuc,0)
      Call Pre_iort(nuc,0)

!----------------------------------------------------------------------
! ... read channel information and find pointer (orbital --> channel):

      Open(nut,file=AF_tar,status='OLD')
      Call R_target(nut)
      Call R_channel(nut,klsp)

      ief=0; ipch=0
      Do ich = 1,nch
       Call EL4_NLK(ELC(ich),nc,lc,kc)
       i=Ifind_nlk(nc,lc,kc,1)
       ief(i)=ich; ipch(ich)=i
      End do

      write(pri,*)
      write(pri,'(a)')   'channel data: '
      write(pri,*)
      write(pri,'(a,i5,a)') 'lpar  =',lpar
      write(pri,'(a,i5,a)') 'ispar =',ispar
      write(pri,'(a,i5,a)') 'ipar  =',ipar
      write(pri,'(a,i5,a)') 'jpar  =',jpar
      ncp = npert
      write(pri,*)
      write(pri,'(a,i5,a)') 'nch   =',nch,'  - number of channels'
      write(pri,'(a,i5,a)') 'ncp   =',ncp,'  - number of perturbers'
      write(pri,'(a,i5,a)') 'nwf   =',nwf,'  - number of orbitals'

!----------------------------------------------------------------------
! ... read B-spline expantions for bound orbitals:
   
      Call Allocate_bsorb(nwf)
      nbf = nwf
      Do i = 1,nbf
       ebs(i)=ELF(i); nbs(i)=NEF(i); lbs(i)=LEF(i); kbs(i)=KEF(i)
       iech(i)=ief(i); mbs(i) = 0
      End do

! ... target radial functions:

      Call Check_file(AF_bsw)
      Open(nuw, file=AF_bsw, form='UNFORMATTED')
      Call Read_bsw(nuw)
      Close(nuw)

! ... perturber radial functions:

      if(nwp.gt.0) then
       i=Index(AFP,'.'); AFP(i+1:i+3)='bsw'
       Call Check_file(AFP)
       Open(nuw, file=AFP, form='UNFORMATTED')
       Call Read_bsw(nuw)
       Close(nuw)
      end if

! ... check the correspondence between c- and bsw-files: 

      j = 0
      Do i = 1,nwf
       if(iech(i).ne.0) Cycle
       if(mbs(i).eq.0) then
        write(pri,'(a,a)') ' Absent expansion for w.f. ',ELF(i)
        j = j + 1
       end if
      End do
      if(j.gt.0) Stop 'no correspondence between c- and w- files'
      
! ... the < B | p > values (convolution with B matrix)

      Do i=1,nbf
       if(iech(i).ne.0) Cycle
       Call BXV(ks,ns,sb(1,1),pbs(1,i),qbs(1,i))
      End do

! ... orthogonality conditions:

      Do i=1,nwf; Do j=1,i
       IORT(j,i)=IORT(i,j); IBORT(i,j)=IORT(i,j); IBORT(j,i)=IORT(i,j)
      End do; End do     

      nort = 0
      Do ich=1,nch; i=ipch(ich)
       Do j=1,nwf; if(ief(i).eq.0) Cycle       
        if(lbs(i).ne.lbs(j)) Cycle
        if(IORT(i,j).ne.0) Cycle
        nort=nort+1
       End do
      End do
      write(pri,*)
      write(pri,'(a,i5,a)') 'nort  = ',nort, &
                            '  - number of orth.constraints to orbitals'
      write(pri,'(a,i5,a)') 'nortb = ',nortb,&
                            '  - number of orth.constraints to bound states'
      write(pri,*)

      End Subroutine Read_data
    

