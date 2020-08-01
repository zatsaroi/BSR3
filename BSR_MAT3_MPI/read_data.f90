!======================================================================
      Subroutine Read_data
!======================================================================
!     read data for given partial wave klsp
!----------------------------------------------------------------------
      Use bsr_mat
      Use target;  Use channel;  Use conf_LS;  Use orb_LS
      Use spline_param; Use spline_orbitals; Use spline_galerkin

      Implicit none
      Integer :: i,j, nc,lc,kc, ich,jch
      Real(8) :: S
      Integer, external :: Ifind_nlk
      Real(8), external :: QUADR
!----------------------------------------------------------------------
! ... files:

      write(ALSP,'(i3.3)') klsp          !  partial wave's sign

! ... log-file:

      i=LEN_TRIM(AF_pri); AF_pri(i-2:i)=ALSP
      Open(pri,file=AF_pri)

! ... c-file:

      i=LEN_TRIM(AF_cfg); AF_cfg(i-2:i)=ALSP
      Open(nuc,file=AF_cfg,status='OLD')

! ... bnk-file:
 
      i=LEN_TRIM(AF_bnk); AF_bnk(i-2:i)=ALSP
      Open(nub,file=AF_bnk,status='OLD',form='UNFORMATTED')

! ... out file:

      i=LEN_TRIM(AF_mat); AF_mat(i-2:i)=ALSP
      Open(nui,file=AF_mat,form='UNFORMATTED')

! ... debug file:

      if(nud.gt.0) then
       i=LEN_TRIM(AF_deb); AF_deb(i-2:i)=ALSP
       Open(nud,file=AF_deb)
      end if

!----------------------------------------------------------------------
! ... HEADER:

       write(pri,'(a,i3)')   'BSR_MAT:   klsp = ',klsp
       write(pri,'(a)')      '***********************'

!----------------------------------------------------------------------
! ... read configurations:
           
      Call Read_conf(nuc,nub)

! ... read physical orbitals if they not in c-file:

      Call Check_file(AF_bsw)
      Open(nuw,file=AF_bsw,form='UNFORMATTED')
      Call Read_bsw_orb_LS(nuw)
      Close(nuw)

! ... read substitute orbitals:

      Call Check_file(AF_orb)
      Open(nuo,file=AF_orb)
      Call Read_sub_orb_LS(nuo,ntarg)

      write(pri,'(/a/)')   'c-file data: '
      write(pri,'(a,i10,a)')  'ncfg   = ',ncfg, &
        '  -  number of configurations'
      write(pri,'(a,i10,a/)') 'nwf    = ',nwf, &
        '  -  number of orbitals'
      if(debug.gt.0) write(pri,'(14a5)')  (ELF(i),i=1,nwf)

!----------------------------------------------------------------------
! ... channel information and pointers for orbital <--> channel:

      Call R_channel(nut,klsp)

      ief = 0
      Do j = 1,nch
       Call EL4_NLK(ELC(j),nc,lc,kc)
       i=Ifind_nlk(nc,lc,kc,1);  ief(i)=j; ipch(j)=i
      End do

      write(pri,'(/a/)')     'Partial wave LSP (JP): '
      if(jpar.lt.0) &
      write(pri,'(a,i3,a)')  'lpar   =',lpar,   '  -  total L '
      if(ispar.gt.0) &
      write(pri,'(a,i3,a)')  'ispar  =',ispar,  '  -  total 2S+1'
      write(pri,'(a,i3,a)')  'ipar   =',ipar,   '  -  parity'
      if(jpar.ge.0) &
      write(pri,'(a,i3,a)')  'jpar   =',jpar-1, '  -  total 2J'

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
       i=LEN_TRIM(AFP); 
       if(AFP(i-1:i).eq.'.c') i=i-2
       AF=AFP(1:i)//'.bsw'
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
      if(j.gt.0) Stop 'no correspondence between c- and w- files'
      
! ... the < B | p > values (convolution with B matrix)

      Do i=1,nbf
       if(iech(i).ne.0) Cycle
       Call BXV(ks,ns,sb(1,1),pbs(1,i),qbs(1,i))
      End do

!----------------------------------------------------------------------
!                                  additional orthogonality conditions:

      JORT=1
      Do i=1,nwf; Do j=1,nwf
        if(LEF(i).ne.LEF(j)) then
         IORT(i,j)=0
        else
         if(KEF(i).eq.KEF(j)) then
          if(i.ne.j) then
           IORT(i,j)=0
          else
           IORT(i,j)=2
           if(iech(i).eq.0.and.iech(j).eq.0) IORT(i,j)=1
          end if
         elseif(KEF(i)*KEF(j).eq.0) then
          IORT(i,j)=0
         else
          IORT(i,j)=2
         end if
        end if
      End do; End do

      Call R_orth(nuc)
      Call R_orth(nup)
      Do i=1,nwf; Do j=1,i
       IORT(j,i)=IORT(i,j); IBORT(i,j)=IORT(i,j); IBORT(j,i)=IORT(i,j)
      End do; End do     

      Do i=1,nwf; ich=iech(i)
      Do j=1,nwf; jch=iech(j)
       if(LEF(i).ne.LEF(j)) Cycle
       if(    ich.eq.0.and.jch.eq.0) then
        Cycle
       elseif(ich.ne.0.and.jch.ne.0) then
        IBORT(i,j) = i*ibo+j
       elseif(ich.eq.0.and.jch.ne.0) then
        if(IBORT(i,j).eq.0) Cycle
        IBORT(i,j) = j*ibo+i
       elseif(ich.ne.0.and.jch.eq.0) then
        if(IBORT(i,j).eq.0) Cycle
        IBORT(i,j) = i*ibo+j
       end if
      End do
      End do 
       
! ... the < p | p > values 

      OBS = 0.d0
      Do i=1,nbf; if(iech(i).ne.0) Cycle
       Do j=1,i;  if(iech(j).ne.0) Cycle
        if(LEF(i).ne.LEF(j)) Cycle
        S=QUADR(i,j,0); if(abs(S).lt.Eps_ovl) S=0.d0
!        if(i.eq.j) S=1.d0                                  ???
        if(iort(i,j).eq.0) S=0.d0
        OBS(i,j)=S; OBS(j,i)=S
        if(S.eq.0.d0) IBORT(i,j) = 0
        if(S.eq.0.d0) IBORT(j,i) = 0
       End do
      End do

      End Subroutine Read_data
    
