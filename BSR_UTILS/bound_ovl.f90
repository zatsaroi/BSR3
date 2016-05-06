!--------------------------------------------------------------------
!     w.nnn  -->  bound_ovl - list of ion-target weights
!--------------------------------------------------------------------
      Use target
      Use channels
      Use target_ion,   only: nps => ntarg, Eps => etarg,     &
                              LPS => Ltarg, SPS => IStarg, PPS => IPtarg

      Implicit real(8) (A-H,O-Z)
      Real(8), allocatable :: wtarg(:,:), WT(:)
      Integer, allocatable :: ip_done(:)
      Real(8), parameter   :: eps_E = 1.d-7
      
! ... files:

      Character(20) :: targ   = 'target';       Integer :: nut = 1
      Character(20) :: pseudo = 'target_ps';  
      Character(20) :: wnnn   = 'w.nnn';        Integer :: nuw = 10
      Character(20) :: ovls   = 'bound_ovl';    Integer :: nuo = 2
      Character(20) :: AF

!----------------------------------------------------------------------
! ... target information:

      Call Check_file(targ)
      Open(nut,file=targ)
      Call R_target(nut)
      Call R_channels(nut)
      ion = nz-nelc
      Close(nut)

!----------------------------------------------------------------------
! ... target-pseudo information:

      Call Check_file(pseudo)
      Open(nut,file=pseudo)
      Call R_target_ion(nut)
      Close(nut)

!--------------------------------------------------------------------
      Allocate(wtarg(ntarg,nps))
      wtarg = 0.d0

      Allocate(ip_done(nps))
      ip_done = 0

      Do klsp=1,nlsp

       write(AF,'(a6,i3.3)') 'w.',klsp
       if(Icheck_file(AF).eq.0) Cycle
       Open(nuw,file=AF,status='OLD',form='UNFORMATTED')
       rewind(nuw) 
       read(nuw) kch, kcp, khm
       if(allocated(WT)) Deallocate(WT); Allocate(WT(kch+kcp))
       kwt = kch+kcp       

       Do is = 1,khm
        read(nuw) WT,j,E

        ip = 0
        Do i=1,nps

         if(lpar(klsp).ne.LPS(i)) Cycle         
         if(ispar(klsp).ne.SPS(i)) Cycle         
         if(ipar(klsp).ne.PPS(i)) Cycle         
         if(abs(E-Eps(i)).gt.eps_E) Cycle        
         ip = i;  Exit
        End do

        if(ip.eq.0) Cycle 

        ip_done(ip) = 1
        
        jopen = 0
        Do i = 1,nch(klsp)
         if(E.lt.etarg(iptar(klsp,i))) Exit
         jopen = jopen + 1
        End do

        S = 0.d0
        Do ich = 1,jopen; it = iptar(klsp,ich)
         wtarg(it,ip)=wtarg(it,ip) + WT(ich)
         S = S + WT(ich)
        End do
! ???
        if(jopen.gt.0.and.jopen.lt.kwt.and.S.gt.0.d0) then
         SS = SUM(WT(jopen+1:kwt))/S
         jt = iptar(klsp,jopen)
         wtarg(1:jt,ip)=wtarg(1:jt,ip)*(1 + SS)
        end if

       End do   ! over is in w.nnn

      End do    ! over klsp

!----------------------------------------------------------------------
! ... check if we found all needed information:

      Do i=1,nps
       if(ip_done(i).eq.1) Cycle         
       write(*,*) i,' -  pseudo-state not found'
      End do

!----------------------------------------------------------------------
! ... output:

      Open(nuo,file=ovls)

      write(nuo,'(a,i6)') 'ntarg_ion = ', ntarg
      write(nuo,'(a,i6)') 'nps = ', nps

      Do i=1,nps
       write(nuo,'(i8,3i5,f16.8,50E15.5)')   &
             i, LPS(i),SPS(i),PPS(i),EPS(i), wtarg(:,i), SUM(wtarg(:,i))
      End do     

      Close(nuo)


      End ! utiliy bound_ovl
