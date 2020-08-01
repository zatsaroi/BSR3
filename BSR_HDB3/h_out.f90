!======================================================================
      Subroutine H_OUT 
!======================================================================
!     Output results in H.DAT file according to
!     the BELFAST R-matrix code format
!----------------------------------------------------------------------
      USE bsr_hd
      USE target
      USE channel,  only: lch,iptar, lpar,ispar,ipar,jpar
      
      Implicit none
      Integer :: NCONAT(ntarg)
      Integer :: i,j,k,ich,it,NCHAN,LRANG2,NAST,LRGL,NSPN,NPTY,MNP2,MORE
      Real(8) :: BSTO

      Call CPU_time(t0)

      i = INDEX(AF_h,'.'); AF = AF_h(1:i)//ALSP
      Open(nuh,file=AF,form='UNFORMATTED')

! ... basic parameters:

      nchan = kch 
      LRANG2=0
      Do ich=1,nchan; if(lch(ich).gt.LRANG2) LRANG2=lch(ich); End do
      LRANG2=LRANG2+1              ! max. (l+1)

      NAST = ntarg                 ! number of target states
      BSTO = 0.d0                  ! boundary logarithmic derivative,
                                   ! always zero in BSR
                                   
      write(nuh) NELC,NZ,LRANG2,LAMAX,NAST,RA,BSTO
      write(nuh) E_exp(1:ntarg)
      if(jtarg(1).gt.0) then
       ltarg(1:ntarg)=jtarg(1:ntarg)-1; istarg(1:ntarg)=0
      end if
      write(nuh) ltarg (1:ntarg)
      write(nuh) istarg(1:ntarg),iptarg(1:ntarg)

! ... Buttle corrections (not used in BSR):

      write(nuh) (0.d0,0.d0,0.d0,i=1,LRANG2)

! ... Partial wave parameters:

      LRGL =  lpar;  NSPN =  ispar
      if(jpar.gt.0) then
       LRGL =  jpar-1;  NSPN =  0
      end if
      if(ipar.eq.+1) NPTY=0
      if(ipar.eq.-1) NPTY=1
      MNP2  = khm; MORE  = 0
      write(nuh) LRGL,NSPN,NPTY,NCHAN,MNP2,MORE

! ... number of channels, coupled to each target state:

      NCONAT=0
      Do ich=1,nchan
       it=iptar(ich); NCONAT(it)=NCONAT(it)+1
      End do
      write(nuh) NCONAT

! ... orbital angular momentum of continuum electron:

      write(nuh) (lch(i),i=1,nchan)

! ... asymptotic coefficients:

      write(nuh) (((CF(i,j,k),i=1,nchan),j=1,nchan),k=2,lamax+1)

! ... R_matrix poles:
 
      write(nuh) (eval(i),i=khm,1,-1)

! ... surface amplitudes:

      write(nuh) ((WMAT(i,j),i=1,nchan),j=khm,1,-1)

      Close(nuh)

      Call CPU_time(t1)
      write (pri,'(/a,T30,f10.2,a)') 'H_out:,', (t1-t0)/60, ' min.'
      write (*  ,'(/a,T30,f10.2,a)') 'H_out:,', (t1-t0)/60, ' min.'

      End Subroutine H_OUT


!======================================================================
      Subroutine H_OUT1 
!======================================================================
!     Output results in H.DAT file according to
!     the BELFAST R-matrix code format
!----------------------------------------------------------------------
      USE bsr_hd
      USE target
      USE channel
      
      Implicit none
      Integer :: NCONAT(ntarg), npch(kch)
      Integer :: i,j,k,ich,it,NCHAN,LRANG2,NAST,LRGL,NSPN,NPTY,MNP2,MORE
      Real(8) :: BSTO

      Call CPU_time(t0)

      i = INDEX(AF_h,'.'); AF = AF_h(1:i)//ALSP
      Open(nuh,file=AF,form='UNFORMATTED')

! ... basic parameters:

      nchan = kch 
      LRANG2=0
      Do ich=1,nchan
       if(lch(ich).gt.LRANG2) LRANG2 = lch(ich)
      End do
      LRANG2=LRANG2+1              ! max. (l+1)

      NAST = ntarg                 ! number of target states
      BSTO = 0.d0                  ! boundary logarithmic derivative,
                                   ! always zero in BSR
                                   
      write(nuh) NELC,NZ,LRANG2,LAMAX,NAST,RA,BSTO
      write(nuh) (E_exp(ip_exp(i)),i=1,ntarg)
      if(jtarg(1).gt.0) then
       ltarg(1:ntarg)=jtarg(1:ntarg)-1; istarg(1:ntarg)=0
      end if
      write(nuh) (ltarg (ip_exp(i)),i=1,ntarg) 
      write(nuh) (istarg(ip_exp(i)),i=1,ntarg), &
                 (iptarg(ip_exp(i)),i=1,ntarg)

! ... Buttle corrections - don't used in BSR !

      write(nuh) (0.d0,0.d0,0.d0,i=1,LRANG2)

! ... Partial wave parameters:

      LRGL =  lpar;  NSPN =  ispar
      if(jpar.gt.0) then
       LRGL =  jpar-1;  NSPN =  0
      end if
      if(ipar.eq.+1) NPTY=0
      if(ipar.eq.-1) NPTY=1
      MNP2  = khm; MORE  = 0
      write(nuh) LRGL,NSPN,NPTY,NCHAN,MNP2,MORE

! ... number of channels, coupled to each target state:

      NCONAT=0
      Do i=1,ntarg; it=ip_exp(i)      
       Do ich=1,nchan; if(iptar(ich).ne.it) Cycle
        NCONAT(i) = NCONAT(i) + 1
       End do
      End do
      write(nuh) NCONAT

! ... new channel order:

      k = 0
      Do i=1,ntarg; it=ip_exp(i)      
       Do ich=1,nchan; if(iptar(ich).ne.it) Cycle
        k=k+1; npch(k)=ich
       End do
      End do
      if(k.ne.nchan) Stop 'Problems with ipc in H.OUT1'

! ... orbital angular momentum of continuum electron:

      write(nuh) (lch(npch(i)),i=1,nchan)

! ... asymptotic coefficients:

      write(nuh) (((CF(npch(i),npch(j),k),i=1,nchan),j=1,nchan),k=2,lamax+1)

! ... R_matrix poles:
 
      write(nuh) (eval(i),i=khm,1,-1)

! ... surface amplitudes:

      write(nuh) ((WMAT(npch(i),j),i=1,nchan),j=khm,1,-1)

      Close(nuh)

      Call CPU_time(t1)
      write (pri,'(/a,T30,f10.2,a)') 'H1_out:,', (t1-t0)/60, ' min.'
      write (*  ,'(/a,T30,f10.2,a)') 'H1_out:,', (t1-t0)/60, ' min.'

      End Subroutine H_OUT1










