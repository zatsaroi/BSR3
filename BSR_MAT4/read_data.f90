!======================================================================
      Subroutine Read_data
!======================================================================
!     read data for given partial wave, klsp
!----------------------------------------------------------------------
      Use bsr_mat
      Use conf_LS;  Use orb_LS
      Use det_list; Use def_list  

      Implicit none
      Integer :: i,j,m, nc,lc,kc, min_lt, max_lt
      Integer, external :: Ifind_nlk,  Iadd_line
      Character(124) :: line

!----------------------------------------------------------------------
! ... files:

      write(ALSP,'(i3.3)') klsp    ! partial-wave extension

! ... log-file:

      i=LEN_TRIM(AF_pri); AF_pri(i-2:i)=ALSP
      Open(pri,file=AF_pri)

! ... output file:

      i=LEN_TRIM(AF_mat); AF_mat(i-2:i)=ALSP

      interrupt = 0
      intercase = 0
      max_nbuf  = 0

      if(mode.eq.0) then

       Open(nui,file=AF_mat,form='UNFORMATTED');   nuj = nui

      elseif(mode.gt.0) then

       Open(nui,file=AF_mat,form='UNFORMATTED',status='OLD',position='APPEND')
       BACKSPACE(nui)
       read(nui) interrupt, intercase, max_nbuf
       rewind(nui)

       write(pri,'(a,i10)') 'interrupt =',interrupt  
       write(pri,'(a,i10)') 'intercase =',intercase
       write(pri,'(a,i10)') 'max_nbuf  =',max_nbuf  

       AF = trim(AF_mat)//'_new'
       Open(nuj,file=AF,form='UNFORMATTED')

      elseif(mode.lt.0) then

       mode = abs(mode)
       Open(nui,file=AF_mat,form='UNFORMATTED')
       AF = trim(AF_mat)//'_new'
       Open(nuj,file=AF,form='UNFORMATTED')

      end if

! ... debug file:

      if(nud.gt.0) then
       i=LEN_TRIM(AF_deb); AF_deb(i-2:i)=ALSP
       Open(nud,file=AF_deb)
      end if

!----------------------------------------------------------------------
! ... HEADER:

       write(pri,'(a,i5)')   'BSR_MAT:   klsp = ',klsp
       write(pri,'(a)')      '*************************'

!----------------------------------------------------------------------
! ... read configurations and open angular coef.s file:
           
! ... c-file:

      i=LEN_TRIM(AF_cfg); AF_cfg(i-2:i)=ALSP
      Open(nuc,file=AF_cfg,status='OLD')

      if(bp_mode.eq.0) then
 
       i=LEN_TRIM(AF_inf); AF_inf(i-2:i)=ALSP
       Open(nub,file=AF_inf,status='OLD',form='UNFORMATTED')
       
       Call Read_conf(nuc,nub)
       
       Call Load_det(nub)
       Call Load_def(nub)
       
       Close(nub)
       i=LEN_TRIM(AF_int); AF_int(i-2:i)=ALSP
       if(mode.ne.11) &
       Open(nub,file=AF_int,status='OLD',form='UNFORMATTED')
       
      else

       Call Read_conf_bp(nuc)

       i=LEN_TRIM(AF_lst); AF_lst(i-2:i)=ALSP
       if(mode.ne.11) &
       Open(nub,file=AF_lst,status='OLD',form='UNFORMATTED')
      
      end if

! ... read physical orbitals if they not in c-file:

      Call Check_file(AF_bsw)
      Open(nuw,file=AF_bsw,form='UNFORMATTED')
      Call Read_bsw_orb_LS(nuw)
      Close(nuw)

! ... read substitute orbitals:

      Call Check_file(AF_orb)
      Open(nuo,file=AF_orb)
      Call Read_sub_orb_LS(nuo,ntarg)

!----------------------------------------------------------------------
! ... read channel information and find pointer (orbital --> channel):

      Call R_channel(nut,klsp)

      ief = 0
      Do j = 1,nch
       Call EL4_NLK(ELC(j),nc,lc,kc)
       i=Ifind_nlk(nc,lc,kc,1);  ief(i)=j; ipch(j)=i
      End do

      i = nch + npert

!----------------------------------------------------------------------

      write(pri,'(/a/)')   'Partial wave LSP (JP): '
      if(jpar.lt.0) &
      write(pri,'(a,i10,a)')  'lpar   = ',lpar, '  -  total L '
      if(ispar.gt.0) &
      write(pri,'(a,i10,a)')  'ispar  = ',ispar,'  -  total 2S+1'
      write(pri,'(a,i10,a)')  'ipar   = ',ipar, '  -  parity'
      if(jpar.ge.0) &
      write(pri,'(a,i10,a)')  'jpar   = ',jpar-1, '  -  total 2J'

      write(pri,'(/a/)')   'c-file data: '
      write(pri,'(a,i10,a)')  'ncfg   = ',ncfg, &
        '  -  number of configurations'
      write(pri,'(a,i10,a)')  'nwf    = ',nwf, &
        '  -  number of orbitals: '

      if(debug.gt.0) then
       write(pri,'(14a5)')  (ELF(i),i=1,nwf)
       write(pri,'(/a,2i10)') 'Determinant overlaps: ndet,kdet = ',ndet,kdet
       write(pri,'( a,2i10)') 'Determinant factors:  ndef,kdef = ',ndef,kdef
       m = 3*ndet + kdet + 3*ndef + kdef
       write(pri,'(a,T33,F8.1,a)') 'Determinant memory: ',m*4.0/(1024*1024),'  Mb'
      end if

!----------------------------------------------------------------------
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
      
! ... remove orb_LS arrays:

!      Call Alloc_orb_LS(nclosd)           !  ????

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

!----------------------------------------------------------------------
! ... cgeck no-exchange mode:

      if(exch_mode.gt.0) then

       max_lt = 0
       min_lt = 1000
       Do i=1,nbf
        if(iech(i).eq.0) then
         if(lbs(i).gt.max_lt) max_lt=lbs(i)
        else
         if(lbs(i).lt.min_lt) min_lt=lbs(i)
        end if
       End do
       write(pri,*)
       write(pri,'(a,i10,a)')  'lt_max   = ',max_lt, &
        '  -  maximum l-value in target'           
       write(pri,'(a,i10,a)')  'lt_min   = ',min_lt, &
        '  -  minumum l-value in channels'
       if(max_lt + mk .ge. min_lt) Stop 'no-exchage mode is not working'

      end if

       if(exch_mode.eq.1) then
        Close(nui,status='DELETE')
        AF = trim(AF_mat)//'.mask' 
        Open(nui,file=AF,form='UNFORMATTED')
       end if

       if(exch_mode.eq.2) then
        Call Read_aarg('AF_mask',AF)
        Open(nud,file=AF,form='UNFORMATTED',status='OLD')
        Call RW_mask
       end if

! ...  first outout of main dimentions

       rewind(nuj);  write(nuj) ns,nch,npert

      End Subroutine Read_data
    
